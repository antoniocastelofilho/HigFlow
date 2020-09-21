#include <string.h>

#include "domain.h"
#include "utils.h"
#include "wls.h"
#include "Debug-c.h"

sim_boundary *sb_create(hig_cell *bc, bc_type t, mp_mapper *m) {
	if(m == NULL) {
		/* Creating the mapper for the boundary */
		m = mp_create();

		/* Creating the iterator for the boundary condition */
		higcit_celliterator *it;
		it = higcit_create_all_leaves(bc);

		/* Assign the mapper for the iterator */
		mp_assign_from_celliterator(m, it, 0);

		/* Destroying the iterator */
		higcit_destroy(it);
	}

	DECL_AND_ALLOC(sim_boundary, sb, 1);
	sb->bc = bc;
	sb->is_own_hig = 1;
	sb->type = t;
	sb->m = m;
	higcit_celliterator *it = higcit_create_all_leaves(bc);
	sb->numvals = higcit_count(it);
	higcit_destroy(it);
	ALLOC(real, sb->vals, sb->numvals);
	return sb;
}

void sb_set_values(sim_boundary *sb, real *vals) {
	for(int i = 0; i < sb->numvals; i++) {
		sb->vals[i] = vals[i];
	}
}

void sb_set_value(sim_boundary *sb, int i, real val) {
	assert(i < sb->numvals);
	assert(i >= 0);
	sb->vals[i] = val;
}

real sb_get_value(sim_boundary *sb, int i) {
	return sb->vals[i];
}

void sb_set_userid(sim_boundary *sb, int userid) {
        sb->userid = userid;
}

void sb_set_valuetype(sim_boundary *sb, bc_valuetype valuetype) {
        sb->valuetype = valuetype;
}

int sb_get_userid(sim_boundary *sb) {
        return sb->userid;
}

bc_valuetype sb_get_valuetype(sim_boundary *sb) {
        return sb->valuetype;
}

hig_cell *sb_get_higtree(sim_boundary *sb) {
	return sb->bc;
}

void sb_destroy(sim_boundary *sb) {
	if (sb->is_own_hig) {
		hig_destroy(sb->bc);
	}
	mp_destroy(sb->m);
	free(sb->vals);
	free(sb);
}

mp_mapper *sb_get_mapper(sim_boundary *sb) {
	return sb->m;
}

higcit_celliterator *sb_get_celliterator(sim_boundary *bc) {
	return higcit_create_all_leaves(bc->bc);
}

static void cache_wls_choice_destroy(void *ptr)
{
	_cache_wls_choice *cached = (_cache_wls_choice *)ptr;
	assert(cached);

	free(cached->with_bc);
	if(cached->with_bc != cached->without_bc) {
		free(cached->without_bc);
	}
	free(cached);
}

static unsigned
_calc_num_points(unsigned dim, int order)
{
        if (dim == 2) {
                if (order == 3) {
                        return 60;
                } else if (order == 2) {
                        return 45;
                }
        } else if (dim == 3) {
                if (order == 3) {
                        return 150;
                } else if (order == 2) {
                        return 120; // increased to ensure 2nd order
                }
        }
        return 2 * dim * wls_num_min_points(dim, order);
}

static void interpolator_set_order(struct interpolator *inter, unsigned dim, int order)
{
	assert(inter->wls == NULL);
	inter->dim = dim;
	inter->order = order;
	inter->maxpts = _calc_num_points(dim, order);
	inter->numminpoints = wls_num_min_points(dim, order);
}

sim_domain *sd_create(mp_mapper *m) {
	DECL_AND_ALLOC(sim_domain, sd, 1);
	sd->max_numhigtrees = MAXHIGTREESPERDOMAIN;
	sd->numhigtrees = 0;
	ALLOC(hig_cell *, sd->higtrees, sd->max_numhigtrees);

	sd->first_fringe = sd->numhigtrees;

	sd->is_own_hig = 1;
	if (m == NULL) {
		m = mp_create();
	}
	sd->m = m;
	sd->numdirichlet_bcs = 0;
	sd->numneumann_bcs = 0;
	sd->inter.wls = NULL;
	ALLOC(sim_boundary *, sd->dirichlet_bcs, MAXBCSPERDOMAIN);
	ALLOC(sim_boundary *, sd->neumann_bcs, MAXBCSPERDOMAIN);
	sd->cwls = ptm_create(cache_wls_choice_destroy);
	sd->use_cache = 1;

	sd->inter.wls = NULL;
	sd_set_interpolator_order(sd, 3);

	// This is used on Neumann boundary condition, which is, at this
	// moment, limited to error of order 2 (degree 1), which means
	// there is no point in using a higher order interpolation for
	// the boundary condition itself.
	sd->bc_inter.wls = NULL;
	interpolator_set_order(&sd->bc_inter, DIM-1, 1);
	return sd;
}

/** Expands internal storage domain trees storage if needed. */
static void _sd_make_room_for_tree(sim_domain *d)
{
	assert(d->numhigtrees <= d->max_numhigtrees);
	assert(d->first_fringe <= d->numhigtrees);

	/* Expand storage if doesn't fit. */
	if(d->numhigtrees >= d->max_numhigtrees)
	{
		d->max_numhigtrees = 2 * d->max_numhigtrees + 1;
		REALLOC_INFER(d->higtrees, d->max_numhigtrees);
	}
}

int sd_add_higtree(sim_domain *d, hig_cell *c) {
	_sd_make_room_for_tree(d);

	/* Sets the added tree is not fringe. */
	/* Place first fringe tree on the newly added last position. */
	d->higtrees[d->numhigtrees] = d->higtrees[d->first_fringe];

	/* Place new tree in space previously occupied by first fringe tree. */
	d->higtrees[d->first_fringe++] = c;

	return d->numhigtrees++;
}

int sd_add_fringe_higtree(sim_domain *d, hig_cell *c)
{
	_sd_make_room_for_tree(d);

	d->higtrees[d->numhigtrees] = c;
	return d->numhigtrees++;
}

void sd_set_trees_as_fringe(sim_domain *d, unsigned* trees_idx,
	unsigned trees_idx_size)
{
	// Sort indexes, to guarantee no unprocessed idx will be
	// swapped with a processed one.
	qsort(trees_idx, trees_idx_size, sizeof *trees_idx,
                  (int (*)(const void *, const void *))int_cmp);

	unsigned ff = d->first_fringe;

	// Process the indices in reversed order, so that the handling
	// of higher indices will not affect the lower indices.
	for(int i = trees_idx_size - 1; i >= 0; --i) {
		unsigned tidx = trees_idx[i];
		if(tidx >= ff)
			// Tree is already fringe.
			continue;

		// For sure ff > 0.
		--ff;

		// Swaps the non-fringe tree on the new position
		// of first fringe with the newly made fringe tree.
		// The swapped tree for sure is not in fringe,
		// otherwhise it would have been processed in
		// previous iterations, since we sorted the array.
		hig_cell *tmp = d->higtrees[ff];
		d->higtrees[ff] = d->higtrees[tidx];
		d->higtrees[tidx] = tmp;
	}
	d->first_fringe = ff;
}

void sd_set_trees_as_local_domain(sim_domain *d, unsigned* trees_idx,
	unsigned trees_idx_size)
{
	// Sort indexes, to guarantee no unprocessed idx will be
	// swapped with a processed one.
	qsort(trees_idx, trees_idx_size, sizeof *trees_idx,
                  (int (*)(const void *, const void *))int_cmp);

	unsigned ff = d->first_fringe;

	for(unsigned i = 0; i < trees_idx_size; ++i) {
		unsigned tidx = trees_idx[i];
		if(tidx < ff)
			// Tree is already not in fringe.
			continue;

		// Swaps the first fringe tree with the tree to
		// be removed from fringe.
		hig_cell *tmp = d->higtrees[ff];
		d->higtrees[ff] = d->higtrees[tidx];
		d->higtrees[tidx] = tmp;

		// Mark the just swapped tree as out of fringe
		++ff;
	}
	d->first_fringe = ff;
}

void sd_set_higtrees_ownership(sim_domain *d, int is_own_hig) {
	d->is_own_hig = is_own_hig;
}

void sd_get_domain_bounding_box(sim_domain *d, Rect *bbox)
{
	assert(d->first_fringe > 0);
	hig_get_bounding_box(d->higtrees[0], bbox);
	for(unsigned i = 1; i < d->first_fringe; ++i) {
		Rect extra;
		hig_get_bounding_box(d->higtrees[i], &extra);
		rect_expand(bbox, &extra);
	}
}

higcit_celliterator *sd_get_domain_celliterator(sim_domain *d) {
	higcit_celliterator *its[d->first_fringe];
	for(int i = 0; i < d->first_fringe; i++) {
		its[i] = higcit_create_all_leaves(d->higtrees[i]);
	}
	return higcit_create_concat(its, d->first_fringe);
}

int sd_get_num_higtrees(sim_domain *d) {
	return d->numhigtrees;
}

unsigned sd_get_num_fringe_higtrees(sim_domain *d) {
	return d->numhigtrees - d->first_fringe;
}

unsigned sd_get_num_local_higtrees(sim_domain *d) {
	return d->first_fringe;
}

higcit_celliterator *sd_get_bcs_celliterator(sim_domain *d) {
	int numbcs = d->numdirichlet_bcs+d->numneumann_bcs;
	higcit_celliterator *its[numbcs];
	int j = 0;
	for(int i = 0; i < d->numdirichlet_bcs; i++, j++) {
		its[j] = higcit_create_all_leaves(d->dirichlet_bcs[i]->bc);
	}
	for(int i = 0; i < d->numneumann_bcs; i++, j++) {
		its[j] = higcit_create_all_leaves(d->neumann_bcs[i]->bc);
	}
	return higcit_create_concat(its, numbcs);
}

void sd_add_boundary(sim_domain *d, sim_boundary *bc) {
	if (bc->type == DIRICHLET) {
		assert(d->numdirichlet_bcs < MAXBCSPERDOMAIN - 1);
		d->dirichlet_bcs[d->numdirichlet_bcs++] = bc;
	} else if (bc->type == NEUMANN) {
		assert(d->numneumann_bcs < MAXBCSPERDOMAIN - 1);
		d->neumann_bcs[d->numneumann_bcs++] = bc;
	}
}

int sd_get_num_bcs(sim_domain *d, int type) {
	if (type == DIRICHLET) {
		return d->numdirichlet_bcs;
	} else {
		return d->numneumann_bcs;
	}
}

sim_boundary *sd_get_bc(sim_domain *d, int type, int i) {
	if (type == DIRICHLET) {
		assert(i < d->numdirichlet_bcs);
		return d->dirichlet_bcs[i];
	} else {
		assert(i < d->numneumann_bcs);
		return d->neumann_bcs[i];
	}
}

static wls_interpolator *
interpolator_get_wls(struct interpolator *inter) {
	if (inter->wls == NULL) {
		inter->wls = wls_create(inter->dim, inter->order, inter->maxpts);
	}
	return inter->wls;
}

void sd_set_interpolator_order(sim_domain *d, int order)
{
	interpolator_set_order(&d->inter, DIM, order);
}

void sd_use_cache(sim_domain *d, int v) {
	d->use_cache = v;
}

void sfd_use_cache(sim_facet_domain *d, int v) {
	d->cdom->use_cache = v;
}

void sfd_set_interpolator_order(sim_facet_domain *d, int order) {
	sd_set_interpolator_order(d->cdom, order);
}

static _cache_wls_item *
_create_cache_wls_item(int numpts, int numbcpts) {
	DECL_AND_ALLOC(_cache_wls_item, cwls, 1);
	cwls->numpts = numpts;
	cwls->numbcpts = numbcpts;
	ALLOC(int, cwls->ids, numpts + numbcpts);
	ALLOC(real, cwls->w, numpts + numbcpts);
	ALLOC(sim_boundary *, cwls->bcpoints, numbcpts);
	return cwls;
}

static int
_wls_compar(const wls_item *a, const wls_item *b)
{
	return a->dist < b->dist ? -1 : 1;
}

typedef struct _wls_item_list
{
	wls_item *ptr;
	unsigned total_size;
	unsigned numpts;
	real max_dist;
} _wls_item_list;

inline static wls_item*
_get_new_wls_item(unsigned maxpts, _wls_item_list *items, real fdist)
{
	// If we already have enough points and next is farther,
	// we can ignore this point.
	if(fdist >= items->max_dist) {
		if(items->numpts >= maxpts) {
			return NULL;
		}

		items->max_dist = fdist;
	}

	// If our array is full, realloc it.
	if(items->numpts >= items->total_size) {
		assert(items->numpts == items->total_size);

		items->total_size += maxpts;
		REALLOC_INFER(items->ptr, items->total_size);
	}

	wls_item *new = &items->ptr[items->numpts++];
	new->dist = fdist;

	return new;
}

static void
search_cells_in_tree_box(hig_cell *root, mp_mapper *m, CPPoint x, Rect *box,
	unsigned maxpts, _wls_item_list *items,
	struct sim_boundary* sb, unsigned char type)
{
	higcit_celliterator *it;
	for (it = higcit_create_bounding_box(root, box->lo, box->hi);
		!higcit_isfinished(it);
		higcit_nextcell(it))
	{
		hig_cell *c = higcit_getcell(it);
		int cid = mp_lookup(m, hig_get_cid(c));

		if (cid >= 0) {
			Point center;
			hig_get_center(c, center);
			assert(hig_intersect_bounding_box(c, box->lo, box->hi));

			// A new point found to be included in the WLS
			real dist = co_distance(x, center);

			wls_item *new =
				_get_new_wls_item(maxpts, items, dist);
			if(new) {
				new->id = cid;
				new->bcpoint = sb;
				new->type = type;
				POINT_ASSIGN(new->x, center);
			}
		}
	}
	higcit_destroy(it);
}

static void
calc_weight_from_points(struct interpolator *inter, CPPoint x, _wls_item_list *items, real *w)
{
	// Calculate the weights:
	wls_interpolator *wls = interpolator_get_wls(inter);

	int numusedpts = inter->numminpoints;
	bool lastiteration = false;
	while (numusedpts <= items->numpts) {
		wls_set_samples_and_calc(wls, numusedpts, items->ptr, x, w);
		if(wls_is_good_enough(inter->wls, numusedpts, items->ptr, x, w)
			|| lastiteration)
		{
			items->numpts = numusedpts;
			break;
		} else {
			numusedpts += inter->numminpoints;
			if (numusedpts > items->numpts) {
				numusedpts = items->numpts;
				lastiteration = true;
			}
		}
	}
}

// Functions used by the stencil interpolation search, implemented for
// both domain types: sim_domain and sim_facet_domain.
typedef struct stencil_search_funcs {
	void (*stencil_search_box) (void *specific, CPPoint x,
		Rect *box, _wls_item_list *items);
	bool (*find_in_center) (void *specific, CPPoint x, real alpha,
		bool *in_domain, sim_stencil *stn);
} stencil_search_funcs;

static inline void
get_stencil_interpolate(sim_domain *d, const Point x, real delta, real alpha,
	sim_stencil *stn, bool use_dirichlet, _wls_item_list *items,
	const stencil_search_funcs *funcs, void *specific)
{
	real w[d->inter.maxpts];

	bool cache_found = false;
	_cache_wls_choice *cached = NULL;
	if(d->use_cache) {
		cache_found = ptm_lookup(d->cwls, x, (void **)&cached);
	}

	if (cache_found &&
		(!use_dirichlet || cached->with_bc) &&
		(use_dirichlet || cached->without_bc))
	{
		_cache_wls_item *cwls = use_dirichlet ? cached->with_bc : cached->without_bc;

		for(int i = 0; i < cwls->numpts; i++) {
			items->ptr[i].id = cwls->ids[i];
			w[i] = cwls->w[i];
			items->ptr[i].type = 0;
		}
		for(int i = 0; i < cwls->numbcpts; i++) {
			items->ptr[i+cwls->numpts].id = cwls->ids[i+cwls->numpts];
			w[i+cwls->numpts] = cwls->w[i+cwls->numpts];
			items->ptr[i+cwls->numpts].bcpoint = cwls->bcpoints[i];
			items->ptr[i+cwls->numpts].type = 1;
		}
		items->numpts = cwls->numpts + cwls->numbcpts;
	} else {
		items->numpts = 0;
		items->max_dist = 0.0;

		Rect box;
		POINT_SUB_SCALAR(box.lo, x, 5.0 * delta);
		POINT_ADD_SCALAR(box.hi, x, 5.0 * delta);

		// Find the possible points:
		// - Inside domain:
		funcs->stencil_search_box(specific, x, &box, items);

		// - On dirichlet boundary condition:
		if (use_dirichlet) {
			for(int i = 0; i < d->numdirichlet_bcs; i++) {
				sim_boundary *sb = d->dirichlet_bcs[i];
				hig_cell *root = sb->bc;
				mp_mapper *m = sb->m;

				search_cells_in_tree_box(root, m, x, &box,
					d->inter.maxpts, items, sb, 1);
			}
		}

		// Sort the elements
		qsort(items->ptr, items->numpts, sizeof *items->ptr,
			(int (*)(const void *, const void *))_wls_compar);

		// We must clip the elements beyond maxpts
		if(items->numpts > d->inter.maxpts) {
			items->numpts = d->inter.maxpts;
		}

		calc_weight_from_points(&d->inter, x, items, w);

		// Cache just calculated stencil
		if (d->use_cache) {
			int numbcpts = 0;
			int numnonbcpts = 0;
			for(int i = 0; i < items->numpts; i++) {
				if(FLT_EQ(w[i], 0.0)) {
					continue;
				}
				if (items->ptr[i].type == 1) { // Dirichlet
					numbcpts++;
				} else {
					numnonbcpts++;
				}
			}
			_cache_wls_item *cwls =
				_create_cache_wls_item(numnonbcpts, numbcpts);

			numnonbcpts = 0;
			numbcpts = 0;
			for(int i = 0; i < items->numpts; i++) {
				if(FLT_EQ(w[i], 0.0)) {
					continue;
				}
				if (items->ptr[i].type == 0) {
					cwls->ids[numnonbcpts] = items->ptr[i].id;
					cwls->w[numnonbcpts] = w[i];
					numnonbcpts++;
				} else if (items->ptr[i].type == 1) { // Dirichlet
					cwls->ids[cwls->numpts + numbcpts] = items->ptr[i].id;
					cwls->w[cwls->numpts + numbcpts] = w[i];
					cwls->bcpoints[numbcpts] = items->ptr[i].bcpoint;
					numbcpts++;
				}
			}

			if(cached) {
				// If cached is not NULL, then, as per outer if
				// condition, either cached->without_bc or cached->with_bc
				// is null, depending whether or not using dirichlet.
				if(use_dirichlet) {
					assert(!cached->with_bc);
					assert(cached->without_bc);

					// TODO: if numbcpts == 0, maybe
					// reuse cached->without_bc...
					cached->with_bc = cwls;
				} else {
					assert(cached->with_bc);
					assert(!cached->without_bc);

					cached->without_bc = cwls;
				}
			} else {
				ALLOC_INFER(cached, 1);
				if(use_dirichlet) {
					cached->with_bc = cwls;
					if(numbcpts == 0) {
						cached->without_bc = cwls;
					} else {
						cached->without_bc = NULL;
					}
				} else {
					cached->with_bc = NULL;
					cached->without_bc = cwls;
				}

				ptm_insert(d->cwls, x, cached);
			}
		}
	}

	for(int i = 0; i < items->numpts; i++) {
		assert(items->ptr[i].id >= 0);

		if (items->ptr[i].type == 1) { // Dirichlet
			stn_add_to_rhs(stn, -alpha * w[i] * sb_get_value(items->ptr[i].bcpoint, items->ptr[i].id));
		} else {
			stn_add_to_element(stn, items->ptr[i].id, alpha * w[i]);
		}
	}
}

static void get_stencil(sim_domain *d, real delta, const Point x, real alpha,
	sim_stencil *stn, bool use_dirichlet, bool use_neumann,
	const stencil_search_funcs *funcs, void *specific);

// Must only be called if point is outside the domain.
// Returns true if stencil has been filled.
// TODO: generalize to any order...
static inline bool
get_stencil_neumann(sim_domain *d, const Point x, real alpha,
	sim_stencil *stn, const stencil_search_funcs *funcs,
	_wls_item_list *items, void *specific)
{
	// TODO: implement cache for Neumann condition

	// Check the distance of the point to each Neumann BC,
	// also check if the point projects to the BC, along its
	// normal direction.
	struct bc_data {
		Point proj_x;
		sim_boundary *nbc;
		hig_cell *nbc_cell;
		real dist;
		int proj_dir;
	} best = {
		.nbc = NULL,
		.dist = DBL_MAX
	};

	for (int i = 0; i < d->numneumann_bcs; i++) {
		struct bc_data curr;

		curr.nbc = d->neumann_bcs[i];
		hig_cell *tree = curr.nbc->bc;

		// The direction of the projection is the
		// degenerated one.
		curr.proj_dir = hig_get_narrowest_dim(tree);

		Point center;
		hig_get_center(tree, center);

		POINT_ASSIGN(curr.proj_x, x);
		curr.proj_x[curr.proj_dir] = center[curr.proj_dir];

		curr.nbc_cell = hig_get_cell_with_point(tree, curr.proj_x);

		// If no cell containing the projected point was found,
		// this Neumann BC is unrelated to the point, skip it:
		if(!curr.nbc_cell) {
			continue;
		}

		curr.dist = x[curr.proj_dir] - curr.proj_x[curr.proj_dir];

		if(fabs(curr.dist) < fabs(best.dist)) {
			best = curr;
		}
	}

	// No suitable Neumann BC found, we are done here.
	if(!best.nbc) {
		return false;
	}

	const real delta = fabs(best.dist);

	// We have a candidate Neumann BC to use, but if there
	// is a Dirichlet BC closer than the Neumann BC, use
	// it instead:
	for (int i = 0; i < d->numdirichlet_bcs; i++) {
		sim_boundary *nbc = d->dirichlet_bcs[i];
		hig_cell *tree = nbc->bc;

		Rect bbox;
		hig_get_bounding_box(tree, &bbox);
		const real dist = rect_distance_to_point(&bbox, x);
		if(dist < delta) {
			return false;
		}
	}

	// Finally, we get to work with the chosen Neumann BC:

	// Take a point best_dist inside the domain (to guarantee 2nd order):
	//const real minus_h = -2.0 * best.dist;

	// Take a point half Δx into the domain:
	real minus_h;
	{
		real delta = 1.0;
		Point deltas;
		hig_get_delta(best.nbc_cell, deltas);
		for(unsigned dim = 0; dim < DIM; ++dim) {
			if(dim != best.proj_dir) {
				delta *= deltas[dim];
			}
		}
		delta = pow(delta, 1.0 / (DIM-1));

		// Sets a Δx/2 with sign to point innards:
		minus_h = copysign(delta, best.dist) * -0.5;
	}

	Point inside_p;
	POINT_ASSIGN(inside_p, best.proj_x);
	inside_p[best.proj_dir] += minus_h;

	// Take the stencil for the point inside the domain:
	get_stencil(d, delta, inside_p, alpha, stn, true, false,
		funcs, specific);

	// Approximating the function as a line, get the line slope:
	Point proj_center;
	hig_get_center(best.nbc_cell, proj_center);
	real slope;
	if(Pointequal(proj_center, best.proj_x)) {
		// Get the exact value if projection point is over
		// boundary cell center:
		const int nbc_cid = mp_lookup(best.nbc->m, hig_get_cid(best.nbc_cell));
		slope = sb_get_value(best.nbc, nbc_cid);
	} else {
		// Projection is not over a cell center, interpolate the
		// BC value at the point.
		items->numpts = 0;
		items->max_dist = 0.0;

		Rect box;
		POINT_SUB_SCALAR(box.lo, x, 5.0 * delta);
		POINT_ADD_SCALAR(box.hi, x, 5.0 * delta);
		search_cells_in_tree_box(best.nbc->bc, best.nbc->m,
			best.proj_x, &box, d->bc_inter.maxpts, items, NULL, 2);

		// Sort the elements
		qsort(items->ptr, items->numpts, sizeof *items->ptr,
			(int (*)(const void *, const void *))_wls_compar);

		// We must clip the elements beyond maxpts
		if(items->numpts > d->bc_inter.maxpts) {
			items->numpts = d->bc_inter.maxpts;
		}

		// We must remove the projection dimension from the points:
		for(unsigned i = 0; i < items->numpts; ++i) {
			PPoint x = items->ptr[i].x;
			for(unsigned dim = best.proj_dir+1; dim < DIM; ++dim) {
				x[dim-1] = x[dim];
			}
		}

		real w[d->bc_inter.maxpts];
		calc_weight_from_points(&d->bc_inter, best.proj_x, items, w);

		// Calculate the slope with the interpolation:
		slope = 0.0;
		for(unsigned i = 0; i < items->numpts; i++) {
			const int id = items->ptr[i].id;
			slope += w[i] * sb_get_value(best.nbc, id);
		}
	}

	// Add the linear variation (negative because at the right hand size):
	stn_add_to_rhs(stn, alpha * slope * minus_h);

	return true;
}

static bool dirichlet_find_in_center(sim_domain *d, CPPoint x, real alpha,
	sim_stencil *stn)
{
	for(int i = 0; i < d->numdirichlet_bcs; i++) {
		sim_boundary *sb = d->dirichlet_bcs[i];
		hig_cell *root = sb->bc;
		mp_mapper *m = sb->m;
		hig_cell *cx = hig_get_cell_with_point(root, x);
		if (cx != NULL) {
			Point cxcenter;
			hig_get_center(cx, cxcenter);
			if (Pointequal(cxcenter, x)) {
				int cxgid = mp_lookup(m, hig_get_cid(cx));
				real bcval = sb_get_value(sb, cxgid);
				stn_add_to_rhs(stn, - alpha * bcval);
				return true;
			}
		}
	}

	return false;
}

static void get_stencil(sim_domain *d, real delta, const Point x, real alpha,
	sim_stencil *stn, bool use_dirichlet, bool use_neumann,
	const stencil_search_funcs *funcs, void *specific)
{
	static _wls_item_list items;

	// TODO: Make cache work at this level.

	// Find if point is inside the domain and if it is at a cell
	// center:
	bool in_domain;
	if(funcs->find_in_center(specific, x, alpha, &in_domain, stn)) {
		return;
	}

	if(in_domain) {
		// If in domain, there is a chance it is at the center of a
		// Dirichlet boundary condition.
		if(dirichlet_find_in_center(d, x, alpha, stn)) {
			return;
		}
	} else if(use_neumann) {
		// If outside the domain, the point value may be determined by a
		// Neumann boundary condition.
		if(get_stencil_neumann(d, x, alpha, stn, funcs, &items, specific)) {
			return;
		}
	}

	// If not at center, nor given by a Neumann BC, use the general case:
	get_stencil_interpolate(d, x, delta, alpha, stn, use_dirichlet,
		&items, funcs, specific);
}

static void
cell_search_box(void *param, CPPoint x, Rect *box, _wls_item_list *items)
{
	sim_domain *d = param;
	mp_mapper *m = sd_get_domain_mapper(d);

	for(int k = 0; k < d->numhigtrees; k++) {
		hig_cell *root = d->higtrees[k];
		search_cells_in_tree_box(root, m, x, box,
			d->inter.maxpts, items, NULL, 0);
	}
}

static bool
cell_find_in_center(void *param, CPPoint x, real alpha, bool *in_domain,
	sim_stencil *stn)
{
	sim_domain *d = param;

	*in_domain = false;
	for(int i = 0; i < d->numhigtrees; i++) {
		hig_cell *root = sd_get_higtree(d, i);

		Rect bbox;
		hig_get_bounding_box(root, &bbox);
		if(rect_contains_point(&bbox, x)) {
			*in_domain = true;

			hig_cell *cx = hig_get_cell_with_point(root, x);
			if (cx != NULL) {
				Point cxcenter;
				hig_get_center(cx, cxcenter);
				if (Pointequal(cxcenter, x)) {
					mp_mapper *m = sd_get_domain_mapper(d);
					int cxgid = mp_lookup(m, hig_get_cid(cx));
					stn_add_to_element(stn, cxgid, alpha);
					return true;
				}
			}

			// Assume domains do not overlap,
			// there is no point in keeping searching...
			break;
		}
	}
	return false;
}

static const stencil_search_funcs cell_funcs = {
	.stencil_search_box = cell_search_box,
	.find_in_center = cell_find_in_center
};

void sd_get_stencil(sim_domain *d, const Point center, const Point x, real alpha,
	sim_stencil *stn)
{
	real delta = co_distance(center, x);
	if (FLT_EQ(delta, 0.0)) {
		int numhigtrees = sd_get_num_higtrees(d);
		for(int i = 0; i < numhigtrees; i++) {
			hig_cell *root = sd_get_higtree(d, i);
			hig_cell *c = hig_get_cell_with_point(root, x);
			if (c != NULL) {
				Point delt;
				hig_get_delta(c, delt);
				delta = delt[0];
				break;
			}
		}
	}

	get_stencil(d, delta, x, alpha, stn, true, true, &cell_funcs, d);
}

void sd_get_stencil_without_bc(sim_domain *d, Point x, real delta, real alpha,
	sim_stencil *stn)
{
	get_stencil(d, delta, x, alpha, stn, false, false, &cell_funcs, d);

}

hig_cell *sd_get_cell_with_point(sim_domain *sd, Point x) {
	higcit_celliterator *it;
	int numhigs = sd_get_num_higtrees(sd);
	for(int i = 0; i < numhigs; i++) {
		hig_cell *root = sd_get_higtree(sd, i);
		hig_cell *c = hig_get_cell_with_point(root, x);
		if (c != NULL) {
			return c;
		}
	}
	return NULL;
}

void sfd_adjust_facet_ids(sim_facet_domain *sfd) {
	sim_domain *sd = sfd->cdom;
	int dim = sfd_get_dim(sfd);
	higcit_celliterator *it;
	int numhigs = sd_get_num_higtrees(sd);
	for(int i = 0; i < numhigs; i++) {
		hig_cell *root = sd_get_higtree(sd, i);
		for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
			hig_cell *c = higcit_getcell(it);
			Point cl, ch, ccenter;
			hig_get_lowpoint(c, cl);
			hig_get_highpoint(c, ch);
			hig_get_center(c, ccenter);
			Point ofcenter;
			POINT_ASSIGN(ofcenter, ccenter);
			for(int dir = 0; dir < 2; dir++) {
				if(dir == 0) {
					ofcenter[dim] = cl[dim]-EPSDELTA;
				} else {
					ofcenter[dim] = ch[dim]+EPSDELTA;
				}
				int isoppositefacetdifferent = 0;
				int isoppositefacetsmaller = 0;
				hig_cell *ofc = hig_get_cell_with_point(root, ofcenter);
				for(int j = 0; j < numhigs && ofc == NULL; j++) {
					if (j != i) {
						hig_cell *oroot = sd_get_higtree(sd, j);
						ofc = hig_get_cell_with_point(oroot, ofcenter);
					}
				}
				if (ofc != NULL) {
					Point ofl, ofh;
					hig_get_lowpoint(ofc, ofl);
					hig_get_highpoint(ofc, ofh);
					for(int k = 0; k < DIM; k++) {
						if (k != dim) {
							if (FLT_NE(cl[k], ofl[k]) || FLT_NE(ch[k], ofh[k])) {
								isoppositefacetdifferent = 1;
								if (FLT_GT((ch[k] - cl[k]), (ofh[k] - ofl[k]))) {
									isoppositefacetsmaller = 1;
								}
							}
						}
					}
				}
				if (isoppositefacetdifferent) {
					if (isoppositefacetsmaller) {
						hig_set_fid_with_dim_dir(c, dim, dir, 0);
					}
				} else if (ofc != NULL) {
					if (dir == 0) {
						hig_set_fid_with_dim_dir(c, dim, dir, - hig_get_fid_with_dim_dir(ofc, dim, 1-dir));
					}
				}
			}

		}
		higcit_destroy(it);
	}
}

static void
facet_search_box(void *param, CPPoint x, Rect *box, _wls_item_list *items)
{
	sim_facet_domain* sfd = param;
	higfit_facetiterator *fit;

	sim_domain *d = sfd->cdom;
	mp_mapper *m = sfd->fm;
	unsigned maxpts = d->inter.maxpts;
	int numfidneg = 0;
	for(int k = 0; k < d->numhigtrees; k++) {
		hig_cell *root = d->higtrees[k];
		for (fit = higfit_create_bounding_box_facets(root,
				sfd->dimofinterest, box->lo, box->hi);
			!higfit_isfinished(fit);
			higfit_nextfacet(fit))
		{
			hig_facet *f = higfit_getfacet(fit);
			int fid = mp_lookup(m, hig_get_fid(f));
			Point fcenter;

			if (fid >= 0) {
				hig_get_facet_center(f, fcenter);
				assert(hig_facet_intersect_bounding_box(f,
					box->lo, box->hi));

				// A new point found to be included in the WLS
				real fdist = co_distance(x, fcenter);

				wls_item *new =
					_get_new_wls_item(maxpts, items, fdist);
				if(new) {
					new->id = fid;
					new->type = 0;
					POINT_ASSIGN(new->x, fcenter);
				}
			}
		}
		higfit_destroy(fit);
	}
}

static bool
facet_find_in_center(void *param, CPPoint x, real alpha,
	bool *in_domain, sim_stencil *stn)
{
	sim_facet_domain *sfd = param;
	sim_domain *d = sfd->cdom;

	*in_domain = false;
	for(int i = 0; i < d->numhigtrees; i++) {
		hig_cell *root = sd_get_higtree(d, i);

		Rect bbox;
		hig_get_bounding_box(root, &bbox);
		if(rect_contains_point(&bbox, x)) {
			*in_domain = true;

			hig_facet fx;
			bool hasf = hig_get_facet_with_point(root, sfd->dim, x, &fx);
			if (hasf) {
				Point fxcenter;
				hig_get_facet_center(&fx, fxcenter);
				if (Pointequal(fxcenter, x)) {
					mp_mapper *m = sfd->fm;
					int fxgid = mp_lookup(m, hig_get_fid(&fx));
					if (fxgid >= 0) {
						stn_add_to_element(stn, fxgid, alpha);
						return true;
					}
				}
			}
		}
	}
	return false;
}

static const stencil_search_funcs facet_funcs = {
	.stencil_search_box = facet_search_box,
	.find_in_center = facet_find_in_center
};

void sfd_get_stencil(sim_facet_domain *sfd, const Point center, const Point x,
	real alpha, sim_stencil *stn)
{
	real delta = co_distance(center, x);
	if (FLT_EQ(delta, 0.0)) {
		int numhigtrees = sfd_get_num_higtrees(sfd);
		for(int i = 0; i < numhigtrees; i++) {
			hig_cell *root = sfd_get_higtree(sfd, i);
			hig_cell *c = hig_get_cell_with_point(root, x);
			if (c != NULL) {
				Point delt;
				hig_get_delta(c, delt);
				delta = delt[0];
				break;
			}
		}
	}

	get_stencil(sfd->cdom, delta, x, alpha, stn, true, true, &facet_funcs, sfd);
}

void sfd_get_stencil_without_bc(sim_facet_domain *sfd, Point x, real delta,
	real alpha, sim_stencil *stn)
{
	get_stencil(sfd->cdom, delta, x, alpha, stn, false, false, &facet_funcs, sfd);
}

int sfd_get_dim(sim_facet_domain *sfd) {
	return sfd->dim;
}

void sfd_set_sfbi(sim_facet_domain *sfd, int i, sim_facet_block_info *sfbi) {
	sfd->sfbi[i] = sfbi;
}

void sd_destroy(sim_domain *d) {
	if(d->inter.wls) {
		wls_destroy(d->inter.wls);
	}
	for(int i = 0; i < d->numdirichlet_bcs; i++) {
		sb_destroy(d->dirichlet_bcs[i]);
	}
	for(int i = 0; i < d->numneumann_bcs; i++) {
		sb_destroy(d->neumann_bcs[i]);
	}
	if (d->is_own_hig) {
		for(int i = 0; i < d->numhigtrees; i++) {
			hig_destroy(d->higtrees[i]);
		}
	}
	mp_destroy(d->m);
	free(d->dirichlet_bcs);
	free(d->neumann_bcs);
	free(d->higtrees);
	ptm_destroy(d->cwls);
	free(d);
}

void sd_create_boundary(hig_cell *root, sim_domain *sd, int type) {
	Point epsdelta;
	int nc[DIM];
	hig_get_cells_per_dim(root, nc);
	POINT_ASSIGN_SCALAR(epsdelta, EPSDELTA);
	Point l, h;

	hig_get_lowpoint(root, l);
	hig_get_highpoint(root, h);
	higcit_celliterator *it;
	for(int d = 0; d < DIM;	d++) {
		Point bl, bh;
		POINT_ASSIGN(bl, l);
		POINT_ASSIGN(bh, h);
		int bnc[DIM];
		POINT_ASSIGN(bnc, nc);
		for (int dir = 0; dir < 2; dir++) {
			POINT_ASSIGN(bl, l);
			POINT_ASSIGN(bh, h);
			bnc[d] = 1;
			if (dir == 0) {
				bl[d] = h[d];
			} else {
				bh[d] = l[d];
			}
			POINT_SUB(bl, bl, epsdelta);
			POINT_ADD(bh, bh, epsdelta);
			hig_cell *bcg = hig_create_root(bl, bh);
			hig_refine_uniform(bcg, bnc);
			int refinenc[DIM];
			POINT_ASSIGN_SCALAR(refinenc, 2);
			refinenc[d] = 1;
			for(it = higcit_create_bounding_box(root, bl, bh); !higcit_isfinished(it); higcit_nextcell(it)) {
				hig_cell *c = higcit_getcell(it);
				Point ccenter;
				Point cdelta;
				hig_get_center(c, ccenter);
				hig_cell *bcell = hig_get_cell_with_point(bcg, ccenter);
				if (bcell != NULL) {
					hig_get_delta(c, cdelta);
					Point bcdelta;
					hig_get_delta(bcell, bcdelta);
					int toRefine = 0;
					for (int i = 0; i < DIM; i++) {
						if (i != d) {
							if (FLT_LT(cdelta[i], bcdelta[i])) {
								toRefine = 1;
								break;
							}
						}
					}
					if (toRefine) {
						hig_refine_uniform(bcell, refinenc);
					}
				}
			}
			higcit_destroy(it);
			mp_mapper *bm = mp_create();
			it = higcit_create_all_leaves(bcg);
			mp_assign_from_celliterator(bm, it, 0);
			higcit_destroy(it);
			sim_boundary *bc = sb_create(bcg, type, bm);
			sd_add_boundary(sd, bc);
		}
	}
}


hig_cell * sd_get_higtree(sim_domain *d, int i) {
	assert(i < d->numhigtrees);
	return d->higtrees[i];
}

hig_cell *sd_get_fringe_higtree(sim_domain *d, unsigned i)
{
	register unsigned idx = i + d->first_fringe;
	assert(i + d->first_fringe < d->numhigtrees);
	return d->higtrees[idx];
}

hig_cell *sd_get_local_higtree(sim_domain *d, unsigned i)
{
	assert(i < d->first_fringe);
	return sd_get_higtree(d, i);
}

hig_cell *sd_get_local_higtree(sim_domain *d, unsigned i);


mp_mapper * sd_get_domain_mapper(sim_domain *d) {
	return d->m;
}

int sd_get_local_id(sim_domain *sd, hig_cell *c)
{
	const int cid = hig_get_cid(c);
	assert(cid >= 0);
	return mp_lookup(sd->m, cid);
}

sim_facet_domain *sfd_create(mp_mapper *fm, int dim) {
	DECL_AND_ALLOC(sim_facet_domain, sfd, 1);
	sfd->__cm = mp_create();
	sfd->cdom = sd_create(sfd->__cm);
	sfd->dim = dim;
	if (fm == NULL) {
		fm = mp_create();
	}
	sfd->fm = fm;
	POINT_ASSIGN_SCALAR(sfd->dimofinterest, 0);
	sfd->dimofinterest[dim] = 1;
	memset(sfd->sfbi, 0, MAXHIGTREESPERDOMAIN * sizeof *sfd->sfbi);
	return sfd;
}

int sfd_add_higtree(sim_facet_domain *d, hig_cell *c) {
	return sd_add_higtree(d->cdom, c);
}

void sfd_set_higtrees_ownership(sim_facet_domain *d, int is_own_hig) {
	sd_set_higtrees_ownership(d->cdom, is_own_hig);
}

void sfd_add_boundary(sim_facet_domain *d, sim_boundary *bc) {
	sd_add_boundary(d->cdom, bc);
}

void sfd_destroy(sim_facet_domain *sfd) {
	mp_destroy(sfd->fm);
	sd_destroy(sfd->cdom);
	for(unsigned i = 0; i < MAXHIGTREESPERDOMAIN; ++i) {
		free(sfd->sfbi[i]);
	}
	free(sfd);
}

void sfd_copy_higtrees_from_center_domain(sim_facet_domain *sfd, sim_domain *sd) {
	for(int i = 0; i < sd->numhigtrees; i++) {
		sfd_add_higtree(sfd, sd_get_higtree(sd, i));
	}
	sfd_set_higtrees_ownership(sfd, 0);

	// Copy fringe information
	sfd->cdom->first_fringe = sd_get_num_local_higtrees(sd);
}

int sfd_get_num_higtrees(sim_facet_domain *d) {
	return sd_get_num_higtrees(d->cdom);
}

unsigned sfd_get_num_fringe_higtrees(sim_facet_domain *d)
{
	return sd_get_num_fringe_higtrees(d->cdom);
}

unsigned sfd_get_num_local_higtrees(sim_facet_domain *d)
{
	return sd_get_num_local_higtrees(d->cdom);
}

int sfd_get_num_bcs(sim_facet_domain *d, int type) {
	return sd_get_num_bcs(d->cdom, type);
}

void sfd_create_boundary(hig_cell *root, sim_facet_domain *sd, int type) {
	sd_create_boundary(root, sd->cdom, type);
}

sim_boundary *sfd_get_bc(sim_facet_domain *d, int type, int i) {
	return sd_get_bc(d->cdom, type, i);
}

static higfit_facetiterator *_sfd_get_facetiterator(sim_facet_domain *sfd,
	unsigned base_hig_idx, unsigned numhigs)
{
	int numblocks = 0;
	for(int hig = 0; hig < numhigs; hig++) {
		hig_cell *domaingrid = sd_get_higtree(sfd->cdom,
			base_hig_idx + hig);
		numblocks += hig_get_number_of_children(domaingrid);
	}

	int dim = sfd_get_dim(sfd);
	DECL_AND_ALLOC(higfit_facetiterator *, its, numblocks);
	int numits = 0;
	for(int hig = 0; hig < numhigs; hig++) {
		hig_cell *domaingrid = sd_get_higtree(sfd->cdom,
			base_hig_idx + hig);
		int numchildren = hig_get_number_of_children(domaingrid);
		for(int i = 0; i < numchildren; i++) {
			its[numits] = sfd_get_facetiterator_for_block(sfd,
				base_hig_idx + hig, i);
			numits++;
		}
	}
	higfit_facetiterator *it = higfit_create_concat(its, numits);

	free(its);

	return it;
}

higfit_facetiterator *sfd_get_domain_facetiterator(sim_facet_domain *sfd) {
	int numhigs = sfd_get_num_local_higtrees(sfd);
	return _sfd_get_facetiterator(sfd, 0, numhigs);
}

higfit_facetiterator *sfd_get_fringe_facetiterator(sim_facet_domain *sfd)
{
	int numhigs = sfd_get_num_fringe_higtrees(sfd);
	return _sfd_get_facetiterator(sfd,
		sd_get_num_local_higtrees(sfd->cdom), numhigs);
}

higfit_facetiterator *sfd_get_facetiterator_for_block(sim_facet_domain *sfd, int hig, int pos) {
	hig_cell *domaingrid = sd_get_higtree(sfd->cdom, hig);
	int dim = sfd_get_dim(sfd);
	hig_cell *c = hig_get_child(domaingrid, pos);
	Point l, h;
	hig_get_lowpoint(c, l);
	hig_get_highpoint(c, h);
	POINT_ADD_SCALAR(l, l, EPSDELTA);
	POINT_SUB_SCALAR(h, h, EPSDELTA);
	if (sfd->sfbi[hig][pos].l) {
		l[dim] -= 2.0*EPSDELTA;
	}
	if (sfd->sfbi[hig][pos].h) {
		h[dim] += 2.0*EPSDELTA;
	}
	higfit_facetiterator *it = higfit_create_bounding_box_facets(domaingrid, sfd->dimofinterest, l, h);
	//return higfit_create_bounding_box_facets(domaingrid, sfd->dimofinterest, sfd->cdom->lp[0], sfd->cdom->hp[0]);
	return it;
}

higcit_celliterator *sfd_get_bcs_celliterator(sim_facet_domain *d) {
	return sd_get_bcs_celliterator(d->cdom);
}

hig_cell *sfd_get_higtree(sim_facet_domain *d, int i) {
	return sd_get_higtree(d->cdom, i);
}

hig_cell *sfd_get_fringe_higtree(sim_facet_domain *d, unsigned i)
{
	return sd_get_fringe_higtree(d->cdom, i);
}

hig_cell *sfd_get_local_higtree(sim_facet_domain *d, unsigned i)
{
	return sd_get_local_higtree(d->cdom, i);
}

int sfd_get_local_id(sim_facet_domain *sfd, hig_facet *f)
{
	const int fid = hig_get_fid(f);
	if(fid < 0) {
		return -1;
	}
	return mp_lookup(sfd_get_domain_mapper(sfd), fid);
}

mp_mapper *sfd_get_domain_mapper(sim_facet_domain *d) {
	return d->fm;
}

hig_cell *sfd_get_cell_with_point(sim_facet_domain *sfd, Point x) {
	return sd_get_cell_with_point(sfd->cdom, x);
}

int sfd_get_facet_with_point(sim_facet_domain *sfd, Point x, hig_facet *facet) {
	sim_domain *sd = sfd->cdom;
	int numhigs = sfd_get_num_higtrees(sfd);
	int dim = sfd_get_dim(sfd);
	for(int hig = 0; hig < numhigs; hig++) {
		hig_cell *root = sfd_get_higtree(sfd, hig);
		if (hig_get_facet_with_point(root, dim, x, facet)) {
			return 1;
		}
	}
	return 0;
}

sim_stencil *stn_create() {
	DECL_AND_ALLOC(sim_stencil, stn, 1);
	stn->maxelems = 100;
	ALLOC(real, stn->vals, stn->maxelems);
	ALLOC(int, stn->ids, stn->maxelems);
	stn->rhs = 0.0;
	stn->numelems = 0;
	return stn;
}

void stn_reset(sim_stencil *stn) {
	stn->numelems = 0;
	stn->rhs = 0.0;
}

void stn_add_to_rhs(sim_stencil *stn, real v) {
    //DEBUG_PASS;
	// if (stn != NULL)
    stn->rhs += v; // Adicionei o teste pra ver se eh != NULL pois o codigo estava quebrando aqui. Depois disso passou! *********************
    //DEBUG_PASS;
}

void stn_set_rhs(sim_stencil *stn, real v) {
	stn->rhs = v;
}

static void
_stn_insert_element(sim_stencil *stn, int id, real v) {
	if (stn->numelems >= stn->maxelems) {
		int maxelems = 2*stn->maxelems;
		DECL_AND_ALLOC(real, vals, maxelems);
		DECL_AND_ALLOC(int, ids, maxelems);
		for(int i = 0; i < stn->numelems; i++) {
			vals[i] = stn->vals[i];
			ids[i] = stn->ids[i];
		}
		free(stn->vals);
		free(stn->ids);
		stn->vals = vals;
		stn->ids = ids;
		stn->maxelems = maxelems;
	}
	stn->ids[stn->numelems] = id;
	stn->vals[stn->numelems] = v;
	stn->numelems++;
}

void stn_add_to_element(sim_stencil *stn, int id, real v) {

    //DEBUG_PASS;
    //DEBUG_INSPECT(stn->numelems, %d);

    for(int i = 0; i < stn->numelems; i++) {

        //DEBUG_INSPECT(i, %d);
        //DEBUG_INSPECT(stn->ids[i], %d);
		if (id == stn->ids[i]) {
			stn->vals[i] += v;
			return;
		}
	}
    //DEBUG_PASS;
    _stn_insert_element(stn, id, v);
    //DEBUG_PASS;
}

void stn_add_stencil(sim_stencil *stnto, sim_stencil *stnfrom) {
	for(int i = 0; i < stnfrom->numelems; i++) {
		stn_add_to_element(stnto, stnfrom->ids[i], stnfrom->vals[i]);
	}
}

void stn_mult_scalar_add(sim_stencil *stnto, sim_stencil *stnfrom, real s)
{
	for(int i = 0; i < stnfrom->numelems; i++) {
		stn_add_to_element(stnto, stnfrom->ids[i], s*stnfrom->vals[i]);
	}
	stn_add_to_rhs(stnto, s * stn_get_rhs(stnfrom));
}

void stn_set_element(sim_stencil *stn, int id, real v) {
	for(int i = 0; i < stn->numelems; i++) {
		if (id == stn->ids[i]) {
			stn->vals[i] = v;
			return;
		}
	}
	_stn_insert_element(stn, id, v);
}

int stn_get_id(sim_stencil *stn, int i) {
	return stn->ids[i];
}

real stn_get_val(sim_stencil *stn, int i) {
	return stn->vals[i];
}

int *stn_get_ids(sim_stencil *stn) {
	return stn->ids;
}

real *stn_get_vals(sim_stencil *stn) {
	return stn->vals;
}

real stn_get_rhs(sim_stencil *stn) {
	return stn->rhs;
}

int stn_get_numelems(sim_stencil *stn) {
	return stn->numelems;
}

real stn_mult_vector(sim_stencil *stn, real *v) {
	real res = 0.0;
	for(int i = 0; i < stn->numelems; i++) {
		res += stn->vals[i]*v[stn->ids[i]];
	}
	return res;
}

void stn_destroy(sim_stencil *stn) {
	free(stn->ids);
	free(stn->vals);
	free(stn);
}
