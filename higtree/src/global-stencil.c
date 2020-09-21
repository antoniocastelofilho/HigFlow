#include <string.h>
#include "global-stencil.h"

global_stencil *gs_create()
{
	DECL_AND_ALLOC(global_stencil, gs, 1);
	memset(gs, 0, sizeof * gs);

	return gs;
}

void gs_destroy(global_stencil *gs)
{
	free(gs->gids);
	free(gs->vals);
	free(gs);
}

static void resize(global_stencil *gs, size_t new_size)
{
	REALLOC_INFER(gs->gids, new_size);
	REALLOC_INFER(gs->vals, new_size);
	gs->allocsize = new_size;
}

void gs_reset(global_stencil *gs)
{
	gs->numelems = 0;
	gs->rhs = 0.0;
}

void gs_add_to_element(global_stencil *gs, int gid, real value)
{
	for(size_t j = 0; j < gs->numelems; ++j) {
		if(gs->gids[j] == gid) {
			gs->vals[j] += value;
			return;
		}
	}

	// GID not found. Make room and insert it.
	if(gs->allocsize <= gs->numelems) {
		assert(gs->allocsize == gs->numelems);
		resize(gs, 2*gs->allocsize + 1);
	}

	assert(gs->numelems < gs->allocsize);
	gs->gids[gs->numelems] = gid;
	gs->vals[gs->numelems] = value;
	++gs->numelems;
}

typedef int (*fn_gid_to_lid)(void *param, int localid);

static inline void add_stencil(global_stencil *gs, sim_stencil *stn,
	real scale, fn_gid_to_lid lid2gid, void *param)
{
	const size_t numelems = stn_get_numelems(stn);
	const real *vals = stn_get_vals(stn);
	const int *lids = stn_get_ids(stn);

	if(gs->allocsize < numelems) {
		resize(gs, gs->allocsize + numelems);
	}

	for(size_t i = 0; i < numelems; ++i) {
		int src_gid = lid2gid(param, lids[i]);

		gs_add_to_element(gs, src_gid, scale * vals[i]);
	}
}

void gs_add_cell_stencil(global_stencil *gs, psim_domain *psd,
	sim_stencil *stn, real scale)
{
	add_stencil(gs, stn, scale, (fn_gid_to_lid)psd_lid_to_gid, psd);
}

void gs_add_facet_stencil(global_stencil *gs, psim_facet_domain *psfd,
	sim_stencil *stn, real scale)
{
	add_stencil(gs, stn, scale, (fn_gid_to_lid)psfd_lid_to_gid, psfd);
}

void gs_add_to_rhs(global_stencil *gs, real value)
{
	gs->rhs += value;
}

size_t gs_get_numelems(global_stencil *gs)
{
	return gs->numelems;
}

int *gs_get_ids(global_stencil *gs)
{
	return gs->gids;
}

real *gs_get_vals(global_stencil *gs)
{
	return gs->vals;
}

real gs_get_rhs(global_stencil *gs)
{
	return gs->rhs;
}
