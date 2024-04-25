#include "Debug-c.h"

#include <string.h>
#include <assert.h>

#include <mpi.h>
#include <glib.h>

#include "mapper-syncer.h"
#include "higtree-parallel.h"

#include "pdomain.h"

enum {BCTYPE_TAG = 10, VAL_TAG = 11, SIZE_TAG = 12};


static void
_neighbor_proc_destroy(struct neighbor_proc *nb)
{
	free(nb->to_send);
	free(nb->to_recv_trees);
	free(nb);
}

partition_graph *pg_create(MPI_Comm comm)
{
	DECL_AND_ALLOC(partition_graph, pg, 1);

	pg->comm = comm;
	pg->fringesize = 1;

	pg->neighbors = g_hash_table_new_full(g_direct_hash, g_direct_equal,
		NULL, (GDestroyNotify)_neighbor_proc_destroy);

	pg->tree_props = g_hash_table_new_full(g_direct_hash, g_direct_equal,
		NULL, free);

	return pg;
}

void pg_set_fringe_size(partition_graph *pg, unsigned fs)
{
	pg->fringesize = fs;
}

unsigned pg_get_fringe_size(partition_graph *pg)
{
	return pg->fringesize;
}

int pg_get_rank(partition_graph *pg)
{
	int ret;
	MPI_Comm_rank(pg->comm, &ret);
	return ret;
}

int pg_get_numnodes(partition_graph *pg)
{
	int ret;
	MPI_Comm_size(pg->comm, &ret);
	return ret;
}

unsigned pg_get_num_neighbors(partition_graph *pg)
{
	return g_hash_table_size(pg->neighbors);
}

void pg_destroy(partition_graph *pg)
{
	g_hash_table_destroy(pg->neighbors);
	g_hash_table_destroy(pg->tree_props);
	free(pg);
}

MPI_Comm pg_get_MPI_comm(partition_graph *pg)
{
	return pg->comm;
}

/** This function enumerates globally all the leaf cells of the
 * domain, determining the global enumeration of cells.
 */
static int _calc_sizes(_dp_shared *dp_data, partition_graph *pg, long numleafs)
{
	// Find the first global ID of this process.
	// O(n) in the number of processes
	// TODO maybe: think in a better way to do it
	int rank = pg_get_rank(pg);
	int size = pg_get_numnodes(pg);
	int firstid, lastid;

	MPI_Comm comm = pg_get_MPI_comm(pg);
	MPI_Status status;

	if (rank == 0) {
		firstid = 0;
	} else {
		(void)MPI_Recv(&firstid, 1, MPI_INT, rank-1, SIZE_TAG, comm, &status);
	}
	lastid = firstid + numleafs;
	dp_data->firstid = firstid;
	dp_data->local_count = numleafs;

	if (rank < size - 1) {
		(void)MPI_Send(&lastid, 1, MPI_INT, rank+1, SIZE_TAG, comm);
	}

	// Broadcast the total domain size
	int globaldomainsize;
	if (rank == size - 1) {
		globaldomainsize = lastid;
	}
	MPI_Bcast(&globaldomainsize, 1, MPI_INT, size - 1, comm);

	return globaldomainsize;
}

int psd_get_global_id(psim_domain *psd, hig_cell *c)
{
	return psd_lid_to_gid(psd, sd_get_local_id(psd->localdomain, c));
}

inline int psd_lid_to_gid(psim_domain *psd, int localid)
{
	return psd->dp_data.gid_map[localid];
}

static int *_stn_get_gids(_dp_shared *dp_data, sim_stencil *stn)
{
	int *gid_map = dp_data->gid_map;
	int *ids = stn_get_ids(stn);
	assert(stn_get_numelems(stn) >= 0);
	unsigned size = stn_get_numelems(stn);

	for(unsigned i = 0; i < size; ++i)
	{
		ids[i] = gid_map[ids[i]];
	}

	return ids;
}

int* psd_stn_get_gids(psim_domain *psd, sim_stencil *stn)
{
	return _stn_get_gids(&psd->dp_data, stn);
}

/** This function enumerates locally all the leaf cells of the
 * domain (fringe included), giving each one a place in the local storage
 * array or distributed property.
 */
static void _psd_setmapper(psim_domain *psd)
{
	// Assign one local idx to each local domain cell
	higcit_celliterator *it = sd_get_domain_celliterator(psd_get_local_domain(psd));
	sim_domain *sd = psd_get_local_domain(psd);
	mp_mapper *m = sd_get_domain_mapper(sd);
	unsigned first_free = mp_assign_from_celliterator(m, it, 0);
	higcit_destroy(it);

	// Assign one local idx to each local fringe cell
	unsigned localfringesize = 0;
	unsigned num_fringe = sd_get_num_fringe_higtrees(sd);
	for(unsigned i = 0; i < num_fringe; ++i) {
		it = higcit_create_all_leaves(sd_get_fringe_higtree(sd, i));
		unsigned first_idx = first_free;
		first_free = mp_assign_from_celliterator(m, it, first_idx);
		localfringesize += first_free - first_idx;
		higcit_destroy(it);
	}

	psd->dp_data.total_count = psd->dp_data.local_count + localfringesize;
}

typedef size_t (*get_elem_count_fn)(void*, hig_cell*, int lo[DIM], int hi[DIM]);
typedef size_t (*get_elem_ids_fn)(void*, hig_cell*, int lo[DIM], int hi[DIM], int *eids);

static void _sync_gid(_dp_shared *dp_data, psim_domain *bpsd, void *fdata,
	get_elem_count_fn elem_count, get_elem_ids_fn get_eids)
{
	const int tag = 2734811;

	// Sync the gid for the fringe.
	// Send first.
	const guint num_nbs = g_hash_table_size(bpsd->filtered_neighbors);
	int *send_buff[num_nbs];
	MPI_Request send_reqs[num_nbs];

	MPI_Comm comm = pg_get_MPI_comm(bpsd->pg);

	GHashTableIter nb_it;
	g_hash_table_iter_init(&nb_it, bpsd->filtered_neighbors);
	gpointer key;
	struct neighbor_proc *nb;
	while(g_hash_table_iter_next(&nb_it, &key, (gpointer*)&nb)) {
		const int nb_rank = GPOINTER_TO_INT(key);

		int pp[DIM];

		// Calculate send buffer size.
		size_t size = 0;
		for(unsigned i = 0; i < nb->to_send_count; ++i) {
			struct to_send_fringe *sf = &nb->to_send[i];

			size += elem_count(fdata, sf->local_tree, sf->lo_idx, sf->hi_idx);
		}

		// Fill send buffer size.
		ALLOC_INFER(send_buff[nb->idx], size);
		size_t idx = 0;
		for(unsigned i = 0; i < nb->to_send_count; ++i) {
			struct to_send_fringe *sf = &nb->to_send[i];

			idx += get_eids(fdata, sf->local_tree, sf->lo_idx,
				sf->hi_idx, send_buff[nb->idx] + idx);
		}
		assert(idx == size);

		// Send buff have the local id,
		// change it inplace to global id
		for(unsigned i = 0; i < size; ++i) {
			send_buff[nb->idx][i] =
				dp_data->gid_map[send_buff[nb->idx][i]];
		}

		// Send the buffer
		MPI_Isend(send_buff[nb->idx], size, MPI_INT, nb_rank,
			tag, comm, &send_reqs[nb->idx]);
	}

	// Receive fringe gids
	int *recv_buf = NULL;
	size_t recv_buf_size = 0;

	int *eid_buf = NULL;
	size_t eid_buf_size = 0;

	int zero[DIM] = {0};

	for(unsigned nb_idx = 0; nb_idx < num_nbs; ++nb_idx) {
		MPI_Status s;
		MPI_Probe(MPI_ANY_SOURCE, tag, comm, &s);

		int size;
		MPI_Get_count(&s, MPI_INT, &size);

		if(recv_buf_size < size) {
			free(recv_buf);
			ALLOC_INFER(recv_buf, size);
			recv_buf_size = size;
		}

		MPI_Recv(recv_buf, size, MPI_INT, s.MPI_SOURCE, tag, comm,
			MPI_STATUS_IGNORE);

		nb = g_hash_table_lookup(bpsd->filtered_neighbors,
			GINT_TO_POINTER(s.MPI_SOURCE));

		int *biter = recv_buf;
		for(unsigned k = 0; k < nb->to_recv_count; ++k) {
			hig_cell *root = nb->to_recv_trees[k];
			size_t e_count = elem_count(fdata, root, zero, root->numcells);

			if(eid_buf_size < e_count) {
				free(eid_buf);
				ALLOC_INFER(eid_buf, e_count);
				eid_buf_size = e_count;
			}

			get_eids(fdata, root, zero, root->numcells, eid_buf);
			for(size_t j = 0; j < e_count; ++j) {
				assert(eid_buf[j] >= dp_data->local_count);
				dp_data->gid_map[eid_buf[j]] = *(biter++);
			}
		}
	}
	free(eid_buf);
	free(recv_buf);

	// Cleanup send
	MPI_Waitall(num_nbs, send_reqs, MPI_STATUSES_IGNORE);
	for(unsigned i = 0; i < num_nbs; ++i) {
		free(send_buff[i]);
	}
	MPI_Barrier(comm);
}

static size_t _cell_sub_counter(void *ign, hig_cell* tree, int lo[DIM], int hi[DIM])
{
	size_t size = 0;
	higcit_celliterator *it;

	int pp[DIM];
	POINT_ASSIGN(pp, lo);
	do {
		hig_cell *c = hig_get_child_in_grid(tree, pp);
		it = higcit_create_all_leaves(c);
		size += higcit_count(it);
		higcit_destroy(it);
	} while(subgrid_coord_advance(pp, lo, hi));

	return size;
}

static size_t _cell_get_ids(void *fdata, hig_cell* tree,
	int lo[DIM], int hi[DIM], int *eids)
{
	psim_domain *psd = fdata;
	higcit_celliterator *it;
	mp_mapper *mp = sd_get_domain_mapper(psd->localdomain);

	size_t count = 0;
	int pp[DIM];
	POINT_ASSIGN(pp, lo);
	do {
		hig_cell *c = hig_get_child_in_grid(tree, pp);
		it = higcit_create_all_leaves(c);

		for(hig_cell *t = higcit_getcell(it);
			!higcit_isfinished(it);
			t = higcit_nextcell(it))
		{
			eids[count++] = mp_lookup(mp, hig_get_cid(t));
		}
		higcit_destroy(it);
	} while(subgrid_coord_advance(pp, lo, hi));

	return count;
}

static size_t _facet_sub_counter(void *fdata, hig_cell *tree, int *lo, int *hi)
{
	psim_facet_domain *psfd = fdata;
	size_t size = 0;
	higfit_facetiterator *fit;

	sim_facet_domain *sfd = psfd_get_local_domain(psfd);

	// Get tree index inside the sim_domain.
	int tree_idx = GPOINTER_TO_INT(g_hash_table_lookup(psfd->psd->tree_map, tree));

	int pp[DIM];
	POINT_ASSIGN(pp, lo);
	do {
		int cell = hig_toposition(tree->numcells, pp);
		fit = sfd_get_facetiterator_for_block(sfd, tree_idx, cell);
		size += higfit_count(fit);
		higfit_destroy(fit);
	} while(subgrid_coord_advance(pp, lo, hi));

	return size;
}

static size_t _facet_get_ids(void *fdata, hig_cell *tree, int *lo, int *hi, int *ids)
{
	psim_facet_domain *psfd = fdata;
	higcit_celliterator *it;
	higfit_facetiterator *fit;

	sim_facet_domain *sfd = psfd_get_local_domain(psfd);
	mp_mapper *mp = sfd_get_domain_mapper(sfd);

	// Get tree index inside the sim_domain.
	int tree_idx = GPOINTER_TO_INT(g_hash_table_lookup(psfd->psd->tree_map, tree));

	size_t count = 0;
	int pp[DIM];
	POINT_ASSIGN(pp, lo);
	do {
		int cell = hig_toposition(tree->numcells, pp);
		fit = sfd_get_facetiterator_for_block(sfd, tree_idx, cell);

		for(; !higfit_isfinished(fit); higfit_nextfacet(fit))
		{
			hig_facet *f = higfit_getfacet(fit);
			ids[count++] = mp_lookup(mp, hig_get_fid(f));
		}
		higfit_destroy(fit);
	} while(subgrid_coord_advance(pp, lo, hi));

	return count;
}

int psd_set_joint_gid(size_t psd_count, psim_domain *psd[],
		size_t psfd_count, psim_facet_domain *psfd[])
{
	const int tag = 72342741;

	/* Count number of local elements. */
	int local_count = 0;
	for(size_t i = 0; i < psd_count; ++i) {
		local_count += psd[i]->dp_data.local_count;
	}
	for(size_t i = 0; i < psfd_count; ++i) {
		local_count += psfd[i]->dp_data.local_count;
	}

	// Sum counts for all process up to this process.
	int total_up_to_proc;
	MPI_Scan(&local_count, &total_up_to_proc, 1,
		MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	// Assign one gid for each local element.
	const int first_proc_gid = total_up_to_proc - local_count;
	int gid = first_proc_gid;
	for(size_t i = 0; i < psd_count; ++i) {
		psd[i]->dp_data.firstid = gid;

		int dcount = psd[i]->dp_data.local_count;
		int* gid_map = psd[i]->dp_data.gid_map;
		for(int j = 0; j < dcount; ++j) {
			gid_map[j] = gid++;
		}

		// Synchronize the just assigned gids.
		_sync_gid(&psd[i]->dp_data, psd[i], psd[i],
			_cell_sub_counter, _cell_get_ids);
	}
	for(size_t i = 0; i < psfd_count; ++i) {
		psfd[i]->dp_data.firstid = gid;

		int dcount = psfd[i]->dp_data.local_count;
		int* gid_map = psfd[i]->dp_data.gid_map;
		for(int j = 0; j < dcount; ++j) {
			gid_map[j] = gid++;
		}

		// Synchronize the just assigned gids.
		_sync_gid(&psfd[i]->dp_data, psfd[i]->psd, psfd[i],
			_facet_sub_counter, _facet_get_ids);
	}
	assert(total_up_to_proc == gid);

	return first_proc_gid;
}

static void _enumerate_gid(_dp_shared *dp_data)
{
	// Allocate the mapping between fringe local ids to global ids
	ALLOC_INFER(dp_data->gid_map, dp_data->total_count);
	// Set the gid of the local elements
	for(int i = 0; i < dp_data->local_count; ++i) {
		dp_data->gid_map[i] = i + dp_data->firstid;
	}
}

void psd_synced_mapper(psim_domain *psd)
{
	higcit_celliterator *it = sd_get_domain_celliterator(psd_get_local_domain(psd));
	long numleafs = higcit_count(it);
	higcit_destroy(it);
	psd->globaldomainsize = _calc_sizes(&psd->dp_data, psd->pg, numleafs);

	_psd_setmapper(psd);
	_enumerate_gid(&psd->dp_data);
	_sync_gid(&psd->dp_data, psd, psd,
		_cell_sub_counter, _cell_get_ids);

	/*mapper_syncer_info ms = {
		.elem_counter = _cell_counter,
		.send_fill = _cell_send_fill,
		.recv_parse = _cell_recv_parse,
		.psd = psd,
		.data = sd_get_domain_mapper(psd_get_local_domain(psd))
	};*/

	/* Mapper should be local now. */
	//ms_syncmapper(&ms);
}

partition_graph *psd_get_partition_graph(psim_domain *psd) {
	return psd->pg;
}

sim_domain *psd_get_local_domain(psim_domain *psd) {
	return psd->localdomain;
}

static inline gboolean hash_table_contains(GHashTable *hash_table, gconstpointer key)
{
#if GLIB_CHECK_VERSION(2,32,0)
	return g_hash_table_contains(hash_table, key);
#else
	return g_hash_table_lookup_extended(hash_table, key, NULL, NULL);
#endif
}

/*! Takes a full partition_graph, and creates a new neighbor map removing
 * all trees not mentioned in the sim_domain.
 */
static GHashTable *
_create_filtered_neighbors(partition_graph *pg, sim_domain *sd, GHashTable *trees_map)
{
	static const int TAG = 56406;
	GHashTable *res = g_hash_table_new_full(g_direct_hash, g_direct_equal,
		NULL, (GDestroyNotify)_neighbor_proc_destroy);

	// For each neighbor, filter what trees to use...
	size_t num_nbs = g_hash_table_size(pg->neighbors);
	MPI_Request sent_reqs[num_nbs];
	unsigned *sent_trees_idx[num_nbs];

	GHashTableIter iter;
	unsigned nb_idx = 0;
	gpointer key;
	struct neighbor_proc *nb;

	g_hash_table_iter_init (&iter, pg->neighbors);
	while(g_hash_table_iter_next(&iter, &key, (gpointer)&nb))
	{
		int nb_rank = GPOINTER_TO_INT(key);
		ALLOC_INFER(sent_trees_idx[nb_idx],
			nb->to_send_count + nb->to_recv_count + 1);
		unsigned *idx_buf = sent_trees_idx[nb_idx];

		// Starts from 1 instead of 0 because first element will contain
		// the size of filtered-out outgoing elements. The rest of the
		// array will contain the filtered-out incoming elements.
		unsigned count = 1;

		// Identify missing trees on the outgoing list
		for(unsigned i = 0; i < nb->to_send_count; ++i) {
			if(!hash_table_contains(trees_map, nb->to_send[i].local_tree)) {
				idx_buf[count++] = i;
			}
		}
		// First element tells how many elements we have on the outgoing side
		idx_buf[0] = count - 1;

		// Identify missing trees on the incoming list
		for(unsigned i = 0; i < nb->to_recv_count; ++i) {
			if(!hash_table_contains(trees_map, nb->to_recv_trees[i]))
				idx_buf[count++] = i;
		}

		MPI_Isend(idx_buf, count, MPI_UNSIGNED, nb_rank, TAG,
			pg->comm, &sent_reqs[nb_idx]);
		++nb_idx;
	}
	assert(nb_idx == num_nbs);

	// Receive the remote exclusions to build the new neighbor map.
	unsigned rsize = 0;
	unsigned *rbuf = NULL;
	for(nb_idx = 0; nb_idx < num_nbs; ++nb_idx) {
		MPI_Status s;
		MPI_Probe(MPI_ANY_SOURCE, TAG, pg->comm, &s);

		int count;
		MPI_Get_count(&s, MPI_UNSIGNED, &count);
		if(count > rsize) {
			free(rbuf);
			ALLOC_INFER(rbuf, count);
		}
		MPI_Recv(rbuf, count, MPI_UNSIGNED, s.MPI_SOURCE, TAG, pg->comm,
			MPI_STATUS_IGNORE);

		// Find the original neighbor proc
		nb = g_hash_table_lookup(pg->neighbors, GINT_TO_POINTER(s.MPI_SOURCE));

		// Create the new neighbor proc
		DECL_AND_ALLOC(struct neighbor_proc, new_nb, 1);
		ALLOC_INFER(new_nb->to_recv_trees, nb->to_recv_count);

		// Fill new_nb proc with remaining elements.
		// Start from the incoming trees.
		unsigned new_iter = 0;
		unsigned riter = 1;
		unsigned rcount = rbuf[0];
		for(unsigned i = 0; i < nb->to_recv_count; ++i) {
			if(riter <= rcount) {
				if(i == rbuf[riter]) {
					// Not present in remote proc, so must
					// not be included in this one.
					++riter;
					continue;
				} else {
					assert(i < rbuf[riter]);
				}
			}
			if(hash_table_contains(trees_map, nb->to_recv_trees[i])) {
				new_nb->to_recv_trees[new_iter++] =
					nb->to_recv_trees[i];
			}
		}
		new_nb->to_recv_count = new_iter;
		REALLOC_INFER(new_nb->to_recv_trees, new_iter);

		// Finish with outgoing trees
		assert(riter == rcount + 1);
		ALLOC_INFER(new_nb->to_send, nb->to_send_count);

		new_iter = 0;
		for(unsigned i = 0; i < nb->to_send_count; ++i) {
			if(riter < count) {
				if(i == rbuf[riter]) {
					// Not present in remote proc, so must
					// not be included in this one.
					++riter;
					continue;
				} else {
					assert(i < rbuf[riter]);
				}
			}
			if(hash_table_contains(trees_map, nb->to_send[i].local_tree)) {
				memcpy(&new_nb->to_send[new_iter++], &nb->to_send[i],
					sizeof nb->to_send[i]);
			}
		}

		// Check if there is any tree left.
		if(new_nb->to_recv_count > 0 || new_iter > 0) {
			new_nb->to_send_count = new_iter;
			REALLOC_INFER(new_nb->to_send, new_iter);

			// Add to new neighbors map:
			g_hash_table_insert(res, GINT_TO_POINTER(s.MPI_SOURCE),
				new_nb);
		} else {
			_neighbor_proc_destroy(new_nb);
		}
	}
	free(rbuf);

	//! Set the neighbor idx according to its position in the hash table.
	g_hash_table_iter_init(&iter, res);
	for(nb_idx = 0; g_hash_table_iter_next(&iter, NULL, (gpointer *)&nb); ++nb_idx)
	{
		nb->idx = nb_idx;
	}

	MPI_Waitall(num_nbs, sent_reqs, MPI_STATUSES_IGNORE);
	for(nb_idx = 0; nb_idx < num_nbs; ++nb_idx) {
		free(sent_trees_idx[nb_idx]);
	}

	return res;
}

void psd_set_local_domain(psim_domain *psd, sim_domain *d)
{
	assert(psd->localdomain == NULL);
	psd->localdomain = d;

	// Set what is or isn't fringe in this local domain.
	{
		unsigned idx_buf[sd_get_num_higtrees(d)];
		unsigned max, count;

		// First we find what is fringe and should not be:
		count = 0;
		max = sd_get_num_fringe_higtrees(d);
		for(unsigned i = 0; i < max; ++i) {
			hig_cell *c = sd_get_fringe_higtree(d, i);
			struct tree_properties *p =
				g_hash_table_lookup(psd->pg->tree_props, c);
			if(p && !p->is_fringe)
				idx_buf[count++] = i;
		}
		sd_set_trees_as_local_domain(d, idx_buf, count);

		// Then we find what is not fringe, but should be:
		count = 0;
		max = sd_get_num_local_higtrees(d);
		for(unsigned i = 0; i < max; ++i) {
			hig_cell *c = sd_get_local_higtree(d, i);
			struct tree_properties *p =
				g_hash_table_lookup(psd->pg->tree_props, c);
			if(p && p->is_fringe)
				idx_buf[count++] = i;
		}
		sd_set_trees_as_fringe(d, idx_buf, count);
	}

	// Index all trees in the domain:
	psd->tree_map = g_hash_table_new(g_direct_hash, g_direct_equal);
	{
		unsigned size = sd_get_num_higtrees(d);
		for(unsigned i = 0; i < size; ++i) {
			hig_cell *c = sd_get_higtree(d, i);
			g_hash_table_insert(psd->tree_map, c, GUINT_TO_POINTER(i));
		}
	}

	// Filter out unused trees in partition graph.
	psd->filtered_neighbors = _create_filtered_neighbors(psd->pg, d, psd->tree_map);
}

psim_domain *psd_create(sim_domain *subd, partition_graph *pg) {
	DECL_AND_ALLOC(psim_domain, psd, 1);
	memset(psd, 0, sizeof *psd);

	psd->pg = pg;

	if (subd != NULL) {
		psd_set_local_domain(psd, subd);
	}

	return psd;
}

static void _dp_shared_destroy(_dp_shared *dp_data, GHashTable *fn)
{
	if(dp_data->psync) {
		property_sync_datatype *psync = dp_data->psync;
		guint s = g_hash_table_size(fn);
		for(guint i = 0; i < s; ++i) {
			MPI_Type_free(&psync[i].send);
			MPI_Type_free(&psync[i].recv);
		}
		free(psync);
	}

	free(dp_data->gid_map);
}

void psd_destroy(psim_domain *psd)
{
	_dp_shared_destroy(&psd->dp_data, psd->filtered_neighbors);

	if(psd->filtered_neighbors)
		g_hash_table_destroy(psd->filtered_neighbors);
	if(psd->tree_map)
		g_hash_table_destroy(psd->tree_map);

	free(psd);
}

void psfd_destroy(psim_facet_domain *psfd) {
	_dp_shared_destroy(&psfd->dp_data, psfd->psd->filtered_neighbors);

	if(psfd->psd_managed) {
		psd_destroy(psfd->psd);
	}
	free(psfd);
}

partition_graph *psfd_get_partition_graph(psim_facet_domain *psfd) {
	return psfd->psd->pg;
}

sim_facet_domain *psfd_get_local_domain(psim_facet_domain *psfd) {
	return psfd->localdomain;
}

static psim_facet_domain *
_psfd_create_base(sim_facet_domain *subd)
{
	DECL_AND_ALLOC(psim_facet_domain, psfd, 1);
	memset(psfd, 0, sizeof *psfd);
	psfd->localdomain = subd;

	return psfd;
}

psim_facet_domain *psfd_create(sim_facet_domain *subd, psim_domain *psd)
{
	psim_facet_domain *psfd = _psfd_create_base(subd);

	psfd->psd = psd;
	psfd->psd_managed = false;

	return psfd;
}

psim_facet_domain *
psfd_create_from_pg(sim_facet_domain *subd, partition_graph *pg)
{
	psim_facet_domain *psfd = _psfd_create_base(subd);

	psfd->psd = psd_create(subd->cdom, pg);
	psfd->psd_managed = true;

	return psfd;
}

static void
_psfd_setmapper(psim_facet_domain *psfd) {
	sim_facet_domain *sfd = psfd_get_local_domain(psfd);
	mp_mapper *m = sfd_get_domain_mapper(sfd);

	// Assign one local idx to each local domain facet
	higfit_facetiterator *fit = sfd_get_domain_facetiterator(sfd);
	unsigned first_free = mp_assign_from_facetiterator(m, fit, 0);
	higfit_destroy(fit);

	// Assign one local idx to each local fringe facet
	fit = sfd_get_fringe_facetiterator(sfd);
	unsigned localfringesize = mp_assign_from_facetiterator(m, fit, first_free)
		- first_free;
	higfit_destroy(fit);

	psfd->dp_data.total_count = psfd->dp_data.local_count + localfringesize;
}

void psfd_synced_mapper(psim_facet_domain *psfd) {
	sim_facet_domain *sfd = psfd_get_local_domain(psfd);
	higfit_facetiterator *fit = sfd_get_domain_facetiterator(sfd);
	long numleafs = higfit_count(fit);
	higfit_destroy(fit);
	psfd->globaldomainsize =
		_calc_sizes(&psfd->dp_data, psfd_get_partition_graph(psfd), numleafs);

	_psfd_setmapper(psfd);
	_enumerate_gid(&psfd->dp_data);
	_sync_gid(&psfd->dp_data, psfd->psd, psfd,
		_facet_sub_counter, _facet_get_ids);

	/*sim_facet_domain *sfd = psfd_get_local_domain(psfd);
	mapper_syncer_info ms = {
		.elem_counter = _facet_counter,
		.send_fill = _facet_send_fill,
		.recv_parse = _facet_recv_parse,
		.psd = psfd->psd,
		.data = sfd
	};

	ms_syncmapper(&ms);*/
}

int psd_get_first_id(psim_domain *psd) {
	return psd->dp_data.firstid;
}

int psd_get_local_domain_size(psim_domain *psd) {
	return psd->dp_data.local_count;
}

int psd_get_global_domain_size(psim_domain *psd) {
	return psd->globaldomainsize;
}

int psfd_get_first_id(psim_facet_domain *psfd) {
	return psfd->dp_data.firstid;
}

int psfd_get_local_domain_size(psim_facet_domain *psfd) {
	return psfd->dp_data.local_count;
}

int psfd_get_global_domain_size(psim_facet_domain *psfd) {
	return psfd->globaldomainsize;
}

real dp_get_value(distributed_property *dp, int id) {
	assert(id >= 0);
	assert(id < dp->pdata->total_count);
	return dp->values[id];
}

int* psfd_stn_get_gids(psim_facet_domain *psfd, sim_stencil *stn)
{
	return _stn_get_gids(&psfd->dp_data, stn);
}

void dp_set_value(distributed_property *dp, int id, real val) {
	assert(id >= 0);
	assert(id < dp->pdata->total_count);
	dp->values[id] = val;
}

void dp_add_value(distributed_property *dp, int id, real val) {
	dp->values[id] += val;
}

int dp_local_size(distributed_property *dp)
{
	return dp->pdata->local_count;
}

void dp_copy_values(distributed_property *to, distributed_property *from)
{
	const unsigned count = to->pdata->total_count;
	assert(count == from->pdata->total_count);
	memcpy(to->values, from->values, count * sizeof *to->values);
}

void dp_set_all_values(distributed_property *dp, real val)
{
	const real *end = dp->values + dp->pdata->total_count;

	for(real *i = dp->values; i < end; ++i) {
		*i = val;
	}
}

void dp_slv_load_from_solver(distributed_property *dp, solver *s)
{
	// TODO: compare performance of loading fringe directly from solver,
	// or syncing it later.
	slv_get_x_scatter(s, dp->pdata->local_count /* total_count */,
		dp->pdata->gid_map, dp->values);
	dp_sync(dp);
}

static distributed_property *dp_create(_dp_shared *pdata)
{
	DECL_AND_ALLOC(distributed_property, dp, 1);

	dp->pdata = pdata;
	ALLOC(real, dp->values, pdata->total_count);
	memset(dp->values, 0, pdata->total_count * sizeof *dp->values);

	return dp;
}

static void
_add_elem_to_datatype(int *current, int *sizes, int *starts, int eid)
{
	assert(eid >= 0);
	if(*current != -1 && eid == starts[*current] + sizes[*current]) {
		++sizes[*current];
	} else {
		++(*current);
		starts[*current] = eid;
		sizes[*current] = 1;
	}
}

static void
_create_sync_datatypes(_dp_shared *dp_data, GHashTable *fn, void* fdata,
	get_elem_count_fn elem_count, get_elem_ids_fn get_eids)
{
	unsigned num_nbs = g_hash_table_size(fn);
	//mp_mapper *mp = sd_get_domain_mapper(psd->localdomain);

	/* Allocate the struct array, one element per neighbor. */
	ALLOC_INFER(dp_data->psync, num_nbs);

	/* Alloc arrays with the worst case size,
	 * the real size is given by count. */
	DECL_AND_ALLOC(int, sizes, dp_data->total_count);
	DECL_AND_ALLOC(int, starts, dp_data->total_count);

	int *eids = NULL;
	size_t eids_size = 0;

	int zero[DIM] = {0};

	/* Create the data types for each neighbor. */
	GHashTableIter nbiter;
	g_hash_table_iter_init(&nbiter, fn);
	struct neighbor_proc *nb;
	while(g_hash_table_iter_next(&nbiter, NULL, (gpointer)&nb)) {
		/* Create datatype used to send. */
		int count = -1;
		sizes[0] = 0;
		for(unsigned i = 0; i < nb->to_send_count; ++i) {
			struct to_send_fringe *ts = &nb->to_send[i];

			hig_cell *tree = ts->local_tree;

			/* Alloc and get element ids. */
			size_t nids = elem_count(fdata, tree, ts->lo_idx, ts->hi_idx);
			if(nids > eids_size) {
				free(eids);
				eids_size = nids;
				ALLOC_INFER(eids, eids_size);
			}
			get_eids(fdata, tree, ts->lo_idx, ts->hi_idx, eids);

			for(size_t i = 0; i < nids; ++i) {
				_add_elem_to_datatype(&count, sizes, starts, eids[i]);
			}
		}
		++count;

		MPI_Type_indexed(count, sizes, starts, MPI_HIGREAL,
			&dp_data->psync[nb->idx].send);
		MPI_Type_commit(&dp_data->psync[nb->idx].send);

		/* Create datatype used to receive. */
		count = -1;
		sizes[0] = 0;
		for(unsigned i = 0; i < nb->to_recv_count; ++i) {
			hig_cell *tree = nb->to_recv_trees[i];

			/* Alloc and get element ids. */
			size_t nids = elem_count(fdata, tree, zero, tree->numcells);
			if(nids > eids_size) {
				free(eids);
				eids_size = nids;
				ALLOC_INFER(eids, eids_size);
			}
			get_eids(fdata, tree, zero, tree->numcells, eids);


			for(size_t i = 0; i < nids; ++i) {
				_add_elem_to_datatype(&count, sizes, starts, eids[i]);
			}
		}
		++count;

		MPI_Type_indexed(count, sizes, starts, MPI_HIGREAL,
			&dp_data->psync[nb->idx].recv);
		MPI_Type_commit(&dp_data->psync[nb->idx].recv);
	}

	free(eids);

	free(sizes);
	free(starts);
}

distributed_property *psd_create_property(psim_domain *psd) {
	if(!psd->dp_data.psync) {
		_create_sync_datatypes(&psd->dp_data, psd->filtered_neighbors,
			psd, _cell_sub_counter, _cell_get_ids);
		psd->dp_data.bpsd = psd;
	}
	return dp_create(&psd->dp_data);
}

distributed_property *psfd_create_property(psim_facet_domain *psfd)
{
	if(!psfd->dp_data.psync) {
		_create_sync_datatypes(&psfd->dp_data, psfd->psd->filtered_neighbors,
			psfd, _facet_sub_counter, _facet_get_ids);
		psfd->dp_data.bpsd = psfd->psd;
	}
	return dp_create(&psfd->dp_data);
}

real dp_interpolate_from_stencil(distributed_property *dp, sim_stencil *stn)
{
	/* Computing property at a point using the Moving Least Square method */

	/* Evaluating the property at the point using the stencil using the center of the cells */
	real res = 0.0;

	/* Geting value from boundary (Diriclet boundary condition) */
	res += -stn_get_rhs(stn);

	/* Getting the number of elements of internal centers cells */
	int numelems = stn_get_numelems(stn);

	/* Getting the center and weights (MLS calculation) from neighborhood */
	for (int i = 0; i < numelems; i++) {
		/* Getting the identifier */
		int id = stn_get_id(stn, i);

		/* Getting the weight value */
		real w = stn_get_val(stn, i);

		/* Getting the property value at the point with identifier id */
		real value = dp_get_value(dp, id);

		/* Computing the value at the given point */
		res += value * w;
	}
	return res;
}

real sd_dp_interpolate(sim_domain *sd, distributed_property *dp,
	const Point center, const Point p, sim_stencil *stn)
{
	/* Computing property at a point using the Moving Least Square method */

	/* Getting the neighborhood of the point using the center of the grid cells */

	/* Reseting the stencil */
	stn_reset(stn);

	/* Getting the stencil */
	sd_get_stencil(sd, center, p, 1.0, stn);

	/* Geting the value of the property using the Moving Least Square Method */
	return dp_interpolate_from_stencil(dp, stn);
}

real sfd_dp_interpolate(sim_facet_domain *sfd, distributed_property *dp,
	const Point center, const Point p, sim_stencil *stn)
{
	/* Computing property at a point using the Moving Least Square method */

	/* Getting the neighborhood of the point using the center of the facets
	 * of the grid cells */

	/* Reseting the stencil */
	stn_reset(stn);

	/* Getting the stencil */
	sfd_get_stencil(sfd, center, p, 1.0, stn);

	/* Geting the value of the property using the Moving Least Square Method */
	return dp_interpolate_from_stencil(dp, stn);
}

void dp_destroy(distributed_property *dp) {
	free(dp->values);
	free(dp);
}

void dp_print(const char *propname, int rank, distributed_property *dp) {
	const int lastid = dp->pdata->local_count;
	const int *gid_map = dp->pdata->gid_map;
	for(int i = 0; i < lastid; i++) {
		printf("%s - %03d (%03d) - %d - %f\n", propname, gid_map[i],
			i, rank, dp->values[i]);
	}
}

enum {SEND = 0, RECEIVE = 1};

void dp_sync(distributed_property *dp)
{
	psim_domain *bpsd = dp->pdata->bpsd;
	property_sync_datatype *psync = dp->pdata->psync;

	const int tag = 629697;

	MPI_Comm comm = pg_get_MPI_comm(bpsd->pg);

	guint nb_count = g_hash_table_size(bpsd->filtered_neighbors);
	MPI_Request reqs[nb_count * 2];

	GHashTableIter nbiter;
	g_hash_table_iter_init(&nbiter, bpsd->filtered_neighbors);
	struct neighbor_proc *nb;
	gpointer key;
	while(g_hash_table_iter_next(&nbiter, &key, (gpointer)&nb)) {
		int nb_rank = GPOINTER_TO_INT(key);

		MPI_Irecv(dp->values, 1, psync[nb->idx].recv, nb_rank,
			tag, comm, &reqs[nb->idx]);

		MPI_Isend(dp->values, 1, psync[nb->idx].send, nb_rank,
			tag, comm, &reqs[nb->idx + nb_count]);
	}

	MPI_Waitall(nb_count * 2, reqs, MPI_STATUSES_IGNORE);
}

int psfd_get_global_id(psim_facet_domain *psfd, hig_facet *f)
{
	return psfd_lid_to_gid(psfd, sfd_get_local_id(psfd->localdomain, f));
}

inline int psfd_lid_to_gid(psim_facet_domain *psfd, int localid)
{
	return psfd->dp_data.gid_map[localid];
}

static void _compute_block_facet_center(sim_facet_domain *sfd, int hig, int pp[DIM], int dir, Point p) {
	hig_cell *root = sfd_get_higtree(sfd, hig);
	hig_cell *c = hig_get_child_in_grid(root, pp);

	hig_get_center(c, p);

	int dim = sfd_get_dim(sfd);
	if (dir == 0) {
		Point lp;
		hig_get_lowpoint(c, lp);
		p[dim] = lp[dim];
	} else {
		Point hp;
		hig_get_highpoint(c, hp);
		p[dim] = hp[dim];
	}
}

static bool _is_in_neumann_bc(sim_facet_domain *sfd, int hig, int pp[DIM], int dir) {
	Point p;
	_compute_block_facet_center(sfd, hig, pp, dir, p);

	sim_domain *d = sfd->cdom;
	int numbcs = sd_get_num_bcs(d, NEUMANN);
	// TODO: This search could be better than linear...
	for(int i = 0; i < numbcs; i++) {
		sim_boundary *bc = sd_get_bc(d, NEUMANN, i);
		hig_cell *bchig = sb_get_higtree(bc);
		hig_cell *c = hig_get_cell_with_point(bchig, p);
		if (c != NULL) {
			return true;
		}
	}
	return false;
}

static int _has_flow(sim_facet_domain *sfd, int hig, int pp[DIM], int dir) {
	Point p;
	_compute_block_facet_center(sfd, hig, pp, dir, p);
	sim_domain *d = sfd->cdom;
	int numhigs = sfd_get_num_higtrees(sfd);
	// TODO: This search could be better than linear...
	for(int hig2 = 0; hig2 < numhigs; hig2++) {
		if (hig2 != hig) {
			hig_cell *root = sfd_get_higtree(sfd, hig2);
			hig_cell *c = hig_get_cell_with_point(root, p);
			if (c != NULL) {
				return 1;
			}
		}
	}
	return 0;
}

static sim_facet_block_info *
_psfd_compute_local_sfbi(psim_facet_domain *psfd, int hig) {
	sim_facet_domain *sfd = psfd_get_local_domain(psfd);
	int dim = sfd_get_dim(sfd);

	hig_cell *root = sfd_get_local_higtree(sfd, hig);
	unsigned numchildren = hig_get_number_of_children(root);
	DECL_AND_ALLOC(sim_facet_block_info, sfbi, numchildren);
	for(unsigned i = 0; i < numchildren; i++) {
		int pp[DIM];
		hig_tobase(i, root->numcells, pp);

		sfbi[i].l = _is_in_neumann_bc(sfd, hig, pp, 0);

		sfbi[i].h =
			pp[dim] < (root->numcells[dim] - 1) ||
			_is_in_neumann_bc(sfd, hig, pp, 1) ||
			_has_flow(sfd, hig, pp, 1);
	}

	return sfbi;
}

void psfd_compute_sfbi(psim_facet_domain *psfd) {
	const int tag = 84012185;
	sim_facet_domain *sfd = psfd_get_local_domain(psfd);

	/* Calculate psfd form local higs */
	unsigned numhigs = sfd_get_num_local_higtrees(sfd);
	sim_facet_block_info *sfbis[numhigs];

	for(unsigned hig = 0; hig < numhigs; ++hig) {
		sfbis[hig] = _psfd_compute_local_sfbi(psfd, hig);
		sfd_set_sfbi(sfd, hig, sfbis[hig]);
	}

	/* Prepare to receive */
	unsigned num_nbs = g_hash_table_size(psfd->psd->filtered_neighbors);
	unsigned to_recv_total = 0;

	/* Type to store one sim_facet_block_info */
	MPI_Datatype sfbi_type;
	MPI_Type_contiguous(sizeof(sim_facet_block_info), MPI_BYTE, &sfbi_type);
	MPI_Type_commit(&sfbi_type);

	/* Send the relevant blocks according to processes conectivity */
	MPI_Request *send_reqs[num_nbs];

	MPI_Comm comm = pg_get_MPI_comm(psfd->psd->pg);

	GHashTableIter iter;
	g_hash_table_iter_init(&iter, psfd->psd->filtered_neighbors);
	gpointer key;
	struct neighbor_proc *nb;
	for(unsigned i = 0; g_hash_table_iter_next(&iter, &key, (gpointer *)&nb);
		++i)
	{
		int rank = GPOINTER_TO_INT(key);
		assert(i == nb->idx);

		ALLOC_INFER(send_reqs[i], nb->to_send_count);
		for(unsigned j = 0; j < nb->to_send_count; ++j) {
			struct to_send_fringe *fgs = &nb->to_send[j];
			hig_cell *tree = fgs->local_tree;

			int ssize[DIM];
			POINT_SUB(ssize, fgs->hi_idx, fgs->lo_idx);

			MPI_Datatype subarr;

			MPI_Type_create_subarray(DIM, tree->numcells, ssize,
				fgs->lo_idx, MPI_ORDER_C, sfbi_type, &subarr);
			MPI_Type_commit(&subarr);

			unsigned idx = GPOINTER_TO_UINT(
				g_hash_table_lookup(psfd->psd->tree_map, tree)
			);

			MPI_Isend(sfbis[idx], 1, subarr, rank, tag, comm, &send_reqs[i][j]);

			MPI_Type_free(&subarr);
		}

		to_recv_total += nb->to_recv_count;
	}

	// Receive the trees
	psim_domain *psd = psfd->psd;
	unsigned recvd[num_nbs];
	memset(recvd, 0, sizeof recvd);
	for(unsigned i = 0; i < to_recv_total; ++i) {
		MPI_Status s;
		MPI_Probe(MPI_ANY_SOURCE, tag, comm, &s);

		int count;
		MPI_Get_count(&s, sfbi_type, &count);

		DECL_AND_ALLOC(sim_facet_block_info, sfbi, count);
		MPI_Recv(sfbi, count, sfbi_type, s.MPI_SOURCE, tag, comm,
			MPI_STATUS_IGNORE);

		struct neighbor_proc *nb =
			g_hash_table_lookup(psd->filtered_neighbors,
				GINT_TO_POINTER(s.MPI_SOURCE));
		unsigned nb_tree_idx = recvd[nb->idx]++;

		hig_cell *tree = nb->to_recv_trees[nb_tree_idx];
		assert(hig_get_number_of_children(tree) == count);

		unsigned local_tree_idx = GPOINTER_TO_UINT(
			g_hash_table_lookup(psd->tree_map, tree)
		);
		sfd_set_sfbi(sfd, local_tree_idx, sfbi);
	}

	MPI_Type_free(&sfbi_type);

	// Wait and free send requests
	g_hash_table_iter_init(&iter, psfd->psd->filtered_neighbors);
	for(unsigned i = 0; g_hash_table_iter_next(&iter, NULL, (gpointer *)&nb); ++i)
	{
		MPI_Waitall(nb->to_send_count, send_reqs[i], MPI_STATUSES_IGNORE);
		free(send_reqs[i]);
	}
}

// Restart simulation ---> 04_04_23

static int _save_fdp_single(psim_facet_domain *psfdu, distributed_property *dpu, MPI_File *fd, int myrank, int ntasks, int proc_offset) {
    higfit_facetiterator *fit;
    sim_facet_domain *sfdu = psfd_get_local_domain(psfdu);
    mp_mapper *mu = sfd_get_domain_mapper(sfdu);
	int buff_count = DIM + 1;
	real buff[buff_count];
	int buff_size = buff_count * sizeof(real);

	int facet_count = 0;
    for(fit = sfd_get_domain_facetiterator(sfdu); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
        hig_facet *f = higfit_getfacet(fit);
        int flid = mp_lookup(mu, hig_get_fid(f));
        Point fcenter;
        hig_get_facet_center(f, fcenter);
        real val = dp_get_value(dpu, flid);

		POINT_ASSIGN(buff, fcenter);
		buff[DIM] = val;

		MPI_File_write_at(*fd, (proc_offset + facet_count)*buff_size, buff, buff_count, MPI_HIGREAL, MPI_STATUS_IGNORE);

		assert(facet_count < dpu->pdata->local_count);
		facet_count++;
    }
    higfit_destroy(fit);

	return facet_count;
}


static int _save_dp_single(psim_domain *psdu, distributed_property *dpu, MPI_File *fd, int myrank, int ntasks, int proc_offset) {
    higcit_celliterator *cit;
    sim_domain *sdu = psd_get_local_domain(psdu);
    mp_mapper *mu = sd_get_domain_mapper(sdu);
	int buff_count = DIM + 1;
	real buff[buff_count];
	int buff_size = buff_count * sizeof(real);

	int cell_count = 0;
    for(cit = sd_get_domain_celliterator(sdu); !higcit_isfinished(cit); higcit_nextcell(cit)) {
        hig_cell *c = higcit_getcell(cit);
        int clid = mp_lookup(mu, hig_get_cid(c));
        Point ccenter;
        hig_get_center(c, ccenter);
        real val = dp_get_value(dpu, clid);

		POINT_ASSIGN(buff, ccenter);
		buff[DIM] = val;

		MPI_File_write_at(*fd, (proc_offset + cell_count)*buff_size, buff, buff_count, MPI_HIGREAL, MPI_STATUS_IGNORE);

		assert(cell_count < dpu->pdata->local_count);
		cell_count++;
    }
    higcit_destroy(cit);

	return cell_count;
}

static int _load_fdp_single(psim_facet_domain *psfdu, distributed_property *dpu, FILE *fd, int myrank, int ntasks, int global_size) {
    sim_facet_domain *sfdu = psfd_get_local_domain(psfdu);
    mp_mapper *mu = sfd_get_domain_mapper(sfdu);
	int buff_count = DIM + 1;
	real buff[buff_count];
	int buff_size = buff_count * sizeof(real);
	
	hig_facet *f, f_dummy;
	f = &f_dummy;
	int found_count_local = 0;
	for(int i = 0; i < global_size; i++) {
		fread(buff, buff_size, 1, fd);
		if(sfd_get_facet_with_point(sfdu, buff, f) != 0) {
			int flid = mp_lookup(mu, hig_get_fid(f));
			// even if the cell of the facet belongs to the local domain, the facet may not
			if(flid >= 0) { // NOT the case in which it is to the left of the first non-fringe facet
				real val = buff[DIM];
				dp_set_value(dpu, flid, val);
				found_count_local++;
			}
		}
		else continue;
	}

	return found_count_local;
}


static int _load_dp_single(psim_domain *psdu, distributed_property *dpu, FILE *fd, int myrank, int ntasks, int global_size) {
    sim_domain *sdu = psd_get_local_domain(psdu);
    mp_mapper *mu = sd_get_domain_mapper(sdu);
	int buff_count = DIM + 1;
	real buff[buff_count];
	int buff_size = buff_count * sizeof(real);
	
	hig_cell *c;
	int found_count_local = 0;
	for(int i = 0; i < global_size; i++) {
		fread(buff, buff_size, 1, fd);
		c = sd_get_cell_with_point(sdu, buff);
		if(c != NULL) {
			int clid = mp_lookup(mu, hig_get_cid(c));
			real val = buff[DIM];
			dp_set_value(dpu, clid, val);
			found_count_local++;
		}
		else continue;
	}

	return found_count_local;
}


// Save Properties
void save_dp_scalar(psim_domain *psdu, distributed_property *dpu, char *filename_base, int myrank, int ntasks) {
	dp_sync(dpu);

	char filename[1024];

	sprintf(filename, "%s", filename_base);

	MPI_File fd;
    int openerr = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fd);
    if(openerr != MPI_SUCCESS) {
        fprintf(stderr, "Error opening file %s to save property on process %d\n", filename, myrank);
		MPI_Abort(MPI_COMM_WORLD, 1);
    }

	int local_count = _save_dp_single(psdu, dpu, &fd, myrank, ntasks, dpu->pdata->firstid);
	
	MPI_File_close(&fd);

	int global_count;
	MPI_Allreduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	int global_size = psd_get_global_domain_size(psdu);
	assert(global_count == global_size); // check if all cells were saved
}

void save_fdp_vec(psim_facet_domain *psfdu[DIM], distributed_property *dpu[DIM], char *filename_base, int myrank, int ntasks) {
    for(int dim = 0; dim < DIM; dim++) {
		dp_sync(dpu[dim]);

		char filename[1024];

		sprintf(filename, "%s%d", filename_base, dim);

		MPI_File fd;
		int openerr = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fd);
		if(openerr != MPI_SUCCESS) {
			fprintf(stderr, "Error opening file %s to save property on process %d\n", filename, myrank);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

        int local_count = _save_fdp_single(psfdu[dim], dpu[dim], &fd, myrank, ntasks, dpu[dim]->pdata->firstid);

		MPI_File_close(&fd);

		int global_count;
		MPI_Allreduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		int global_size = psfd_get_global_domain_size(psfdu[dim]);
		assert(global_count == global_size); // check if all facets were saved
    }
}

void save_dp_tensor(psim_domain *psdu, distributed_property *dpu[DIM][DIM], char *filename_base, int myrank, int ntasks) {
    for (int dim = 0; dim < DIM; dim++) {
        for (int dim2 = 0; dim2 < DIM; dim2++) {
			char filename[1024];
			sprintf(filename, "%s%d%d", filename_base, dim, dim2);

			save_dp_scalar(psdu, dpu[dim][dim2], filename, myrank, ntasks);
        }
    }
}


// Load Properties
void load_dp_scalar(psim_domain *psdu, distributed_property *dpu, char *filename_base, int myrank, int ntasks) {
	FILE *fd;
	char filename[1024];

	sprintf(filename, "%s", filename_base);

	fd = fopen(filename, "rb");
	if (fd == NULL) {
		fprintf(stderr, "Error opening file %s to load property on process %d\n", filename, myrank);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	int global_size = psd_get_global_domain_size(psdu); // size without fringe

    fseek(fd, 0, SEEK_END);
    int file_size = ftell(fd);
    fseek(fd, 0, SEEK_SET);
	int buff_count = DIM + 1;
	int buff_size = buff_count * sizeof(real);
	assert(file_size == global_size * buff_size); // Check if the file size is correct

	int found_count_local = _load_dp_single(psdu, dpu, fd, myrank, ntasks, global_size);

	fclose(fd);

	int found_count;
	MPI_Allreduce(&found_count_local, &found_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	int total_count_proc = dpu->pdata->total_count;
	int total_count;
	MPI_Allreduce(&total_count_proc, &total_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	assert(found_count == total_count); // check if all cells were loaded - includes fringe

	dp_sync(dpu);
}


void load_fdp_vec(psim_facet_domain *psfdu[DIM], distributed_property *dpu[DIM], char *filename_base, int myrank, int ntasks) {
	for(int dim = 0; dim < DIM; dim++) {
		FILE *fd;
		char filename[1024];

		sprintf(filename, "%s%d", filename_base, dim);

		fd = fopen(filename, "rb");
		if (fd == NULL) {
			fprintf(stderr, "Error opening file %s to load property on process %d\n", filename, myrank);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		int global_size = psfd_get_global_domain_size(psfdu[dim]);  // size without fringe

		fseek(fd, 0, SEEK_END);
		int file_size = ftell(fd);
		fseek(fd, 0, SEEK_SET);
		int buff_count = DIM + 1;
		int buff_size = buff_count * sizeof(real);
		assert(file_size == global_size * buff_size); // Check if the file size is correct

		int found_count_local = _load_fdp_single(psfdu[dim], dpu[dim], fd, myrank, ntasks, global_size);
		
		fclose(fd);

		int found_count;
		MPI_Allreduce(&found_count_local, &found_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		int total_count_proc = dpu[dim]->pdata->total_count;
		int total_count;
		MPI_Allreduce(&total_count_proc, &total_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		assert(found_count == total_count); // check if all facets were loaded - includes fringe

		dp_sync(dpu[dim]);
    }
}


void load_dp_tensor(psim_domain *psdu, distributed_property *dpu[DIM][DIM], char *filename_base, int myrank, int ntasks) {
    for (int dim = 0; dim < DIM; dim++) {
        for (int dim2 = 0; dim2 < DIM; dim2++) {
			char filename[1024];
			sprintf(filename, "%s%d%d", filename_base, dim, dim2);

			load_dp_scalar(psdu, dpu[dim][dim2], filename, myrank, ntasks);
        }
    }
}