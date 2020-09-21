#include <vector>
#include <array>
#include <type_traits>
#include <unordered_map>
#include <queue>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "mpi.h"

extern "C"
{
	#include "higtree-parallel.h"
	#include "higtree-serialize.h"
	#include "term-det.h"
	#include "build-fringe.h"
	#include "Debug-c.h"
}

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

namespace {

using Block = neighbor_proc::to_send_fringe;

// Compute the number of coarse grid neighbors I must expect
// from a given dimension.
constexpr unsigned neighborhood_size(unsigned dim=DIM)
{
	return dim == 1 ? 3 : 3 * neighborhood_size(dim-1);
}
typedef std::integral_constant<unsigned, neighborhood_size()> NeighborhoodSize;

// A trinary number with digits from -1 to 1.
template <unsigned D=DIM>
class Trinary
{
public:
	Trinary()
	{
		for(auto &v: val)
			v = -1;
	}

	void inc()
	{
		for(unsigned i = 0; i < D; ++i) {
			if(++val[i] == 2) {
				val[i] = -1;
			} else {
				return;
			}
		}
	}

	bool is_zero()
	{
		for(auto v: val) {
			if(v != 0)
				return false;
		}
		return true;
	}

	signed char get(unsigned idx)
	{
		return val[idx];
	}

private:
	signed char val[D];
};

template <class F, class V> static void
grid_breadth_search(std::queue<unsigned> &&queue, const int base[DIM], const F& filter, const V& visit)
{
	// Breadth search on the elements of the queue
	while(!queue.empty()) {
		unsigned curr_idx = queue.front();
		queue.pop();

		// Get the coordinates of the cell in tree grid:
		std::array<int, DIM> pp_center;
		hig_tobase(curr_idx, base, pp_center.data());

		// Process each neighbor in tree
		Trinary<> offset;
		for(unsigned j = 0; j < NeighborhoodSize::value; ++j, offset.inc()) {
			// If offset is zero, it is the own cell, skip...
			if(offset.is_zero())
				continue;

			// Calculate the neighbor from offset
			auto pp = pp_center;
			bool valid = true;
			for(unsigned k = 0; k < DIM; ++k) {
				pp[k] += offset.get(k);

				// Check if bounds were violated
				if(pp[k] < 0 || pp[k] >= base[k]) {
					valid = false;
					break;
				}
			}
			// Poor goto replacement...
			if(!valid)
				continue;

			// Get the actual index of the neighbor
			int n = hig_toposition(base, pp.data());

			// Visit it if needed
			if(filter(curr_idx, pp_center.data(), n, pp.data())) {
				queue.push(n);
			}
		}

		visit(curr_idx, pp_center.data());
	}
}

template <class F, class V> static void
grid_breadth_search(unsigned initial, const int base[DIM], const F& filter, const V& visit)
{
	grid_breadth_search(std::queue<unsigned>({initial}), base, filter, visit);
}

struct Displacement: public std::array<coordtype, DIM>
{
	Displacement() = default;

	Displacement(CPPoint p)
	{
		std::copy(p, p+DIM, begin());
	}

	Displacement operator+(coordtype s) const
	{
		Displacement ret;
		POINT_ADD_SCALAR(ret, (*this), s);
		return ret;
	}

	template <typename Array>
	Displacement operator+(const Array& p) const
	{
		return add(*this, p);
	}

	Displacement operator-() const
	{
		Displacement ret;
		POINT_MULT_SCALAR(ret, (*this), -1.0);
		return ret;
	}

	Displacement operator-(coordtype s) const
	{
		Displacement ret;
		POINT_SUB_SCALAR(ret, (*this), s);
		return ret;
	}

	PPoint as_point() {
		return &(*this)[0];
	}

	CPPoint as_point() const {
		return &(*this)[0];
	}

	template <typename ArrayA, typename ArrayB>
	static Displacement add(const ArrayA& a, const ArrayB& b)
	{
		Displacement ret;
		POINT_ADD(ret, a, b);
		return ret;
	}
};

}

// Boost stuff, to make Displacemente recognizable as a point:
namespace boost { namespace geometry { namespace traits {

template<>
struct tag<Displacement>
{
	using type = bg::point_tag;
};

template<>
struct coordinate_type<Displacement>
{
	using type = coordtype;
};

template<>
struct coordinate_system<Displacement>
{
	using type = bg::cs::cartesian;
};

template<>
struct dimension<Displacement>:
	public std::integral_constant<int, DIM>
{};

template<size_t Dimension>
struct access<Displacement, Dimension>
{
	static inline coordtype get(const Displacement& a)
	{
		return a[Dimension];
	}

	static inline void set(Displacement& a, coordtype value)
	{
		a[Dimension] = value;
	}
};

}}}

struct _FringeBuilder {
	typedef std::vector<unsigned> TreeFringe;
	typedef std::unordered_map<unsigned, TreeFringe> DisplacementFringe;
	typedef bgi::rtree<std::pair<Displacement, size_t>,
			bgi::quadratic<4>> ProcFringe;

	std::vector<DisplacementFringe> displacement_fringes;
	std::unordered_map<int, ProcFringe> proc_fringes;
	std::unordered_map<int, unsigned> reverse_neighbor_map;
	load_balancer *partition_ctx;
	_local_connectivity *lc;
	TermDetection *td;

	struct RemNeighbor {
		Rect *isect;
		Displacement disp;
		int nb_rank;
	};
	std::unordered_map<int, std::vector<RemNeighbor>> rem_nbs;

	std::vector<_lb_output_tree> output_buf;

	~_FringeBuilder() {
		_local_connectivity_destroy(lc);
	}

	TreeFringe &get_tree_fringe(int dest_proc_rank, const Displacement& d,
		unsigned tree_idx)
	{
		typedef bg::model::box<Displacement> Box;
		ProcFringe &pf = proc_fringes[dest_proc_rank];

		Box eps(d - EPSDELTA, d + EPSDELTA);
		std::vector<ProcFringe::value_type> es;
		if(!pf.query(bgi::within(eps), std::back_inserter(es))) {
			es.emplace_back(d, displacement_fringes.size());
			displacement_fringes.emplace_back();
			pf.insert(es.back());
		}

		TreeFringe &tf = displacement_fringes[es[0].second][tree_idx];
		if(tf.empty()) {
			int tree_size = hig_get_number_of_children(
				partition_ctx->output[tree_idx].tree
			);

			tf.resize(tree_size);
		}
		return tf;
	}

	unsigned tree_group(unsigned tree_idx)
	{
		return partition_ctx->output[tree_idx].group;
	}
};

_FringeBuilder *_fringe_builder_create(struct _proc_conectivity *pcon,
	_local_portals *lps, load_balancer *ctx, TermDetection *td)
{
	_FringeBuilder *fb = new _FringeBuilder;

	// Build map from neighbor MPI rank to process index
	for(unsigned i = 0; i < pcon->num_neighbors; ++i) {
		fb->reverse_neighbor_map[pcon->neighbors[i].proc] = i;
	}
	fb->partition_ctx = ctx;

	// Finds where each local tree touches each other
	fb->lc = _local_connectivity_build(ctx, lps);

	fb->td = td;

	return fb;
}

void _fringe_builder_destroy(_FringeBuilder *fb)
{
	delete fb;
}

void _fringe_builder_add_remote_tree_neighbor(_FringeBuilder *fb, int nb_rank,
	unsigned tree_idx, CPPoint displacement, Rect *isect)
{
	fb->rem_nbs[tree_idx].push_back({
		.isect=isect,
		.disp=displacement,
		.nb_rank=nb_rank
	});
}

bool _fringe_builder_get_neighbor_idx(_FringeBuilder *fb, int mpi_rank, unsigned *idx)
{
	// reverse_neighbor_map is only relevant for geometrical contact,
	// not portal contact
	auto res = fb->reverse_neighbor_map.find(mpi_rank);
	if(res == fb->reverse_neighbor_map.end()) {
		return false;
	}
	*idx = res->second;
	return true;
}

static void
get_child_bounding_box(hig_cell *tree, unsigned child_idx, Rect &bbox)
{
	hig_cell *child = hig_get_child(tree, child_idx);
	hig_get_bounding_box(child, &bbox);
}

void _fringe_builder_search_from(_FringeBuilder *fb, Rect *isect, unsigned start,
	int dest_proc_rank, CPPoint disp, unsigned tree_idx)
{
	assert(tree_idx < fb->partition_ctx->num_out_higs);
	assert(start > 0);

	int myrank;
	MPI_Comm_rank(fb->partition_ctx->group, &myrank);

	static const Point zero = {0};
	assert(dest_proc_rank != myrank || !co_equal(disp, zero));

	hig_cell *tree = fb->partition_ctx->output[tree_idx].tree;

	// Every element has been initially set to zero. Every element found to be in fringe
	// must be set to a value > 0, which is higher the closer the element is found to be
	// to the origin.
	_FringeBuilder::TreeFringe &fringe_map = fb->get_tree_fringe(dest_proc_rank, disp, tree_idx);

	// Breadth search queue, every element must have fringe value > 1
	std::queue<unsigned> search_from;

	// Mark the cells touching the intersection to bootstrap the
	// breadth search.
	// TODO: use boost's rtree to speed this up.
	{
		for(unsigned i = 0; i < fringe_map.size(); ++i) {
			Rect bbox;
			get_child_bounding_box(tree, i, bbox);
			if(rect_intersect(isect, &bbox, NULL)
				// If cell height is already marked as closer to
				// the fringe origin, don't search from it.
				&& fringe_map[i] < start)
			{
				fringe_map[i] = start;

				// If not already at the limit of the fringe...
				if(start > 1)
					search_from.push(i);
			}
		}
	}

	// We are at the limit of the fringe, the search must not continue...
	if(start <= 1)
		return;

	unsigned tree_group = fb->tree_group(tree_idx);

	// List of local tree neighbors the search must span to
	struct LocSearchSeed {
		unsigned tree;
		unsigned start;
		Displacement disp;
		Rect isect;
	};
	std::vector<LocSearchSeed> local_neighbor_work;

	std::unordered_map<int, std::vector<_search_work_start>> remote_neighbor_work;

	auto filter = [&](unsigned curr_i, int curr_pp[DIM],
		unsigned next_i, int next_pp[DIM])
	{
		unsigned next_dist = fringe_map[curr_i] - 1;

		// Visit it if needed
		if(fringe_map[next_i] < next_dist) {
			fringe_map[next_i] = next_dist;
			if(next_dist > 1)
				return true;
		}
		return false;
	};

	auto visit = [&](unsigned e_idx, int pp[DIM])
	{
		// Don't bother checking if search must span to other trees if we are
		// already on the limit of the fringe:
		if(fringe_map[e_idx] <= 1)
			return;

		// Check if this cell interfaces with another tree
		bool on_border = false;
		for(unsigned k = 0; k < DIM; ++k) {
			if(pp[k] == 0 || pp[k] == tree->numcells[k] - 1) {
				on_border = true;
				break;
			}
		}

		// Skips if cell is not on border...
		if(!on_border)
			return;

		// Find if element interfaces with other tree (local or not).
		Rect bbox;
		get_child_bounding_box(tree, e_idx, bbox);

		unsigned num_local_neighbors;
		_local_neighbor *local_neighbors = _local_connectivity_get_neighbors(fb->lc,
			tree_idx, &num_local_neighbors);

		unsigned next_dist = fringe_map[e_idx] - 1;

		// Search in locally neighboring trees
		for(unsigned i = 0; i < num_local_neighbors; ++i) {
			// Span the search only inside the same group
			if(tree_group != fb->tree_group(local_neighbors[i].neighbor_tree))
				continue;

			// Skip if doen't intersect or intersection is
			// contained in the area that spanned this search.
			Rect &nisect = local_neighbors[i].intersection;
			Rect loc_isect;
			if(!rect_intersect(&nisect, &bbox, &loc_isect)
				|| rect_contains(isect, &loc_isect))
			{
				continue;
			}

			LocSearchSeed lsearch;
			lsearch.disp = Displacement::add(
				local_neighbors[i].displacement, disp);

			// Don't span a search inside same process if
			// displacement is zero.
			if(myrank == dest_proc_rank
				&& co_equal(zero, lsearch.disp.as_point()))
			{
				continue;
			}

			lsearch.tree = local_neighbors[i].neighbor_tree;
			lsearch.start = next_dist;
			rect_translate(local_neighbors[i].displacement,
				&loc_isect, &lsearch.isect);

			local_neighbor_work.push_back(lsearch);
		}

		bool is_displaced = false;
		for(unsigned dim = 0; dim < DIM; ++dim) {
			if(!POS_EQ(disp[dim], 0)) {
				is_displaced = true;
				break;
			}
		}

		// Search contact with interfacing processes
		// in all contacts of current tree with other processes
		for(auto &rem_nb: fb->rem_nbs[tree_idx]) {
			// Skip interfaces with originator process, if not displaced
			if(!is_displaced && rem_nb.nb_rank == dest_proc_rank)
				continue;

			Rect loc_isect;
			// If cell box intersects with the interface area, and this intersection
			// is not contained in the same area that generated current search call,
			// then we forward the work to the remote process
			if(!rect_intersect(rem_nb.isect, &bbox, &loc_isect)
				|| rect_contains(isect, &loc_isect))
			{
				continue;
			}

			auto new_disp = Displacement::add(disp, rem_nb.disp);
			if(rem_nb.nb_rank == dest_proc_rank) {
				continue;
			}

			auto &works = remote_neighbor_work[rem_nb.nb_rank];

			works.emplace_back();
			auto &w = works.back();

			w.search_from_dist = next_dist;
			w.tree_group = tree_group;
			rect_translate(rem_nb.disp.as_point(), &loc_isect, &w.isect);
			POINT_ASSIGN(w.disp, new_disp);
		}
	};

	// Breadth search starting from search_from
	grid_breadth_search(std::move(search_from), tree->numcells, filter, visit);

	// Dispatch all remote search works
	for(auto &work: remote_neighbor_work) {
		_dispatch_remote_fringe_search(fb->td, work.first, dest_proc_rank,
			work.second.size(), &work.second[0]);
	}

	// Perform all local search works
	// TODO: don't do this recursively
	for(auto &work: local_neighbor_work) {
		assert(work.start > 0);
		_fringe_builder_search_from(fb, &work.isect, work.start,
			dest_proc_rank, work.disp.as_point(), work.tree);
	}
}

static unsigned
mark_visited_and_count(const Block &b, std::vector<bool>& visited)
{
	int pp[DIM];
	POINT_ASSIGN(pp, b.lo_idx);

	// Starts from 1 to account for the root node, not included
	unsigned node_count = 1;
	do {
		int i = hig_toposition(b.local_tree->numcells, pp);
		visited[i] = true;

		hig_cell *c = hig_get_child(b.local_tree, i);

		higcit_celliterator *it = higcit_create_all_higtree(c);
		node_count += higcit_count(it);
		higcit_destroy(it);
	} while(subgrid_coord_advance(pp, b.lo_idx, b.hi_idx));

	return node_count;
}

static void
expand(Block &b, const int base[DIM], const _FringeBuilder::TreeFringe &tf,
	std::vector<bool>& visited)
{
	auto filter = [&](unsigned p_i, int p_pp[DIM], unsigned i, int pp[DIM]) {
		if(!visited[i] && tf[i] > 0) {
			visited[i] = true;
			return true;
		}
		return false;
	};

	auto visit = [&](unsigned i, int pp[DIM]) {
		for(unsigned dim = 0; dim < DIM; ++dim) {
			if(pp[dim] < b.lo_idx[dim])
				b.lo_idx[dim] = pp[dim];
			else if(pp[dim] >= b.hi_idx[dim])
				b.hi_idx[dim] = pp[dim] + 1;
		}
	};

	grid_breadth_search(hig_toposition(base, b.lo_idx), base, filter, visit);
}

/*! For each tree seached in a (process,displacement),
 * build the blocks of independent subtress.
 */
static std::vector<Block>
split_blocks(hig_cell* tree, const _FringeBuilder::TreeFringe& tf,
	unsigned &total_node_count)
{
	std::vector<bool> visited(tf.size(), false);

	std::vector<Block> ret;

	total_node_count = 0;

	// Find first elem of new block...
	for(unsigned i = 0; i < tf.size(); ++i) {
		if(visited[i])
			continue;
		visited[i] = true;
		if(tf[i] > 0) {
			ret.emplace_back();
			Block &b = ret.back();
			b.local_tree = tree;

			hig_tobase(i, tree->numcells, b.lo_idx);
			POINT_OPERATOR_SCALAR(b.hi_idx, b.lo_idx, +, 1);

			// The two following steps are done separately because
			// expand may not visit every cell contained in the
			// resulting fringe block
			expand(b, tree->numcells, tf, visited);
			total_node_count += mark_visited_and_count(b, visited);
		}
	}

	return ret;
}

static unsigned
serialize_tree_block(hig_serial_tree *buf, const Block& b, const Displacement& disp)
{
	int pp[DIM];
	POINT_ASSIGN(pp, b.lo_idx);

	// Craft the new root node as a slice of the original tree
	POINT_SUB(buf->numcells, b.hi_idx, b.lo_idx);
	{
		hig_cell *lo_hig = hig_get_child_in_grid(b.local_tree, b.lo_idx);
		hig_get_lowpoint(lo_hig, buf->lowpoint);
		POINT_SUB(buf->lowpoint, buf->lowpoint, disp);
	}
	{
		int hi_pp[DIM];
		POINT_SUB_SCALAR(hi_pp, b.hi_idx, 1);
		hig_cell *hi_hig = hig_get_child_in_grid(b.local_tree, hi_pp);
		hig_get_highpoint(hi_hig, buf->highpoint);
		POINT_SUB(buf->highpoint, buf->highpoint, disp);
	}

	// Serialize all nodes below root
	unsigned node_count = 1;
	do {
		hig_cell *c = hig_get_child_in_grid(b.local_tree, pp);
		node_count += hs_serialize_displaced(&buf[node_count], c,
			(-disp).as_point());
	} while(subgrid_coord_advance(pp, b.lo_idx, b.hi_idx));

	return node_count;
}

_lb_output_tree*
_fringe_builder_assemble_result(_FringeBuilder *fb, unsigned *fringe_tree_count,
	partition_graph *pg)
{
	MPI_Comm comm = pg_get_MPI_comm(pg);
	const int tag = 674809065;

	// I am using TermDetection here because I don't trust I can know
	// beforehand how many processes will be sending fringes to me.
	TermDetection* td = term_det_init(comm);

	// Build the to_send part of the nb, and send the tree slices
	// to neighbors
	{
		// Declared outside the loop so I can realloc only when needed.
		std::vector<hig_serial_tree> node_buf;
		std::vector<unsigned> group_buf;
		std::vector<const Displacement*> disp_buf;
		for(auto& i: fb->proc_fringes) {
			int nb_rank = i.first;
			_FringeBuilder::ProcFringe &pf = i.second;
			std::vector<std::vector<Block>> list_of_blocks;

			DECL_AND_ALLOC(neighbor_proc, nb, 1);
			g_hash_table_insert(pg->neighbors,
				GINT_TO_POINTER(nb_rank), nb);

			nb->to_send_count = 0;
			unsigned total_node_count = 0;
			for(auto k = pf.qbegin(bgi::satisfies([](_FringeBuilder::ProcFringe::value_type){ return true; })); k != pf.qend(); ++k) {
				const Displacement &disp = k->first;
				_FringeBuilder::DisplacementFringe &df = fb->displacement_fringes[k->second];

				// Split the marked fringe nodes in rectangular blocks, and
				// calculate the total number of tree nodes in the process
				for(const auto& j: df) {
					auto &out = fb->partition_ctx->output[j.first];
					hig_cell *tree = out.tree;

					const _FringeBuilder::TreeFringe &tf = j.second;

					unsigned node_count;
					list_of_blocks.emplace_back(split_blocks(tree,
						tf, node_count));
					auto &bls = list_of_blocks.back();
					nb->to_send_count += bls.size();

					assert(group_buf.size() + bls.size() == nb->to_send_count);
					group_buf.resize(nb->to_send_count, out.group);

					total_node_count += node_count;
				}

				// Store the displecement of the new blocks added.
				disp_buf.resize(nb->to_send_count, &disp);
			}

			// Allocate the send buffers with the calculated size
			if(node_buf.size() < total_node_count)
				node_buf.resize(total_node_count);

			// Create and fill the structure to be stored in the
			// partition graph, while serializing the tree nodes
			// to be sent.
			ALLOC(neighbor_proc::to_send_fringe, nb->to_send,
				nb->to_send_count);
			unsigned bl_count = 0;
			unsigned node_count = 0;
			for(const auto& bls: list_of_blocks) {
				for(const auto& bl: bls) {
					nb->to_send[bl_count] = bl;
					node_count += serialize_tree_block(
						&node_buf[node_count], bl,
						*disp_buf[bl_count]
					);
					bl_count++;
				}
			}
			assert(bl_count == nb->to_send_count);
			assert(bl_count == group_buf.size());
			assert(bl_count == disp_buf.size());
			assert(node_count == total_node_count);
			DEBUG_EXEC(
				printf("Number of fringe blocks of proc %d found: %u\n",
					nb_rank, bl_count)
			);

			// Actually send the tree chunks
			{
				void *bufs[] = {
					&nb->to_send_count,
					&group_buf[0],
					&node_buf[0]
				};
				size_t sizes[] = {
					sizeof nb->to_send_count,
					nb->to_send_count * sizeof group_buf[0],
					node_count * sizeof node_buf[0],
					0
				};

				term_det_send_asyncv(td, bufs, sizes, nb_rank);
			}

			group_buf.clear();
			disp_buf.clear();
		}
	}

	// Receive the fringe trees
	int from_rank;
	int msg_size;
	void *buf;
	while((buf = term_det_wait_for_msg(td, &from_rank, &msg_size, true))) {
		auto nb = static_cast<neighbor_proc *>(
			g_hash_table_lookup(pg->neighbors,
				GINT_TO_POINTER(from_rank))
		);

		unsigned *as_uint = static_cast<unsigned *>(buf);
		nb->to_recv_count = as_uint[0];
		ALLOC(hig_cell *, nb->to_recv_trees, nb->to_recv_count);

		unsigned *groups = as_uint + 1;
		auto rtrees = reinterpret_cast<hig_serial_tree *>(
			as_uint + 1 + nb->to_recv_count
		);

		unsigned processed_nodes = 0;
		for(unsigned i = 0; i < nb->to_recv_count; ++i) {
			hig_cell *tree;
			processed_nodes += hs_deserialize(&tree,
				&rtrees[processed_nodes]);

			fb->output_buf.push_back({.tree = tree, .group = groups[i]});
			nb->to_recv_trees[i] = tree;

			DECL_AND_ALLOC(tree_properties, props, 1);
			props->is_fringe = true;
			props->partition_group = groups[i];
			g_hash_table_insert(pg->tree_props, tree, props);
		}
		assert(processed_nodes * sizeof(hig_serial_tree)
			== msg_size - (nb->to_recv_count + 1) * sizeof(unsigned));
	}

	term_det_destroy(td);

	// Fill neighbor proc idx according to its final order
	// inside the hash table.
	{
		GHashTableIter it;
		g_hash_table_iter_init(&it, pg->neighbors);
		neighbor_proc *nb;
		unsigned count = 0;
		while(g_hash_table_iter_next(&it, NULL, (gpointer *)(&nb))) {
			nb->idx = count++;
		}
	}

	*fringe_tree_count = fb->output_buf.size();
	return &fb->output_buf[0];
}
