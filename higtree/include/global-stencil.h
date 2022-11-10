#ifndef GLOBAL_STENCIL_H
#define GLOBAL_STENCIL_H

#include "pdomain.h"

struct global_stencil
{
	size_t numelems;
	size_t allocsize;

	int *gids;
	real *vals;

	real rhs;
};
typedef struct global_stencil global_stencil;

global_stencil *gs_create();

void gs_destroy(global_stencil *gs);

void gs_reset(global_stencil *gs);

void gs_add_cell_stencil(global_stencil *gs, psim_domain *psd,
	sim_stencil *stn, real scale);

void gs_add_facet_stencil(global_stencil *gs, psim_facet_domain *psfd,
	sim_stencil *stn, real scale);

void gs_add_to_element(global_stencil *gs, int gid, real value);

void gs_add_to_rhs(global_stencil *gs, real value);

size_t gs_get_numelems(global_stencil *gs);

int *gs_get_ids(global_stencil *gs);

real *gs_get_vals(global_stencil *gs);

real gs_get_rhs(global_stencil *gs);

#endif
