
#include<stdio.h>
#include<stdlib.h>

#include "utils.h"
#include "higtree.h"
#include "higtree-iterator.h"
#include "mapper.h"
#include "wls.h"
#include "glib.h"
#include "domain.h"



int main(int argc, char *argv[]) {
	Point l, h;
	POINT_ASSIGN_REALS(l, 0.0, 0.0);
	POINT_ASSIGN_REALS(h, 1.0, 1.0);
	hig_cell *root = hig_create_root(l, h);
	int nc[DIM];
	POINT_ASSIGN_INTS(nc, 5, 5);
	hig_refine_uniform(root, nc);
	mp_mapper *m = mp_create();
	higcit_celliterator *it;
	it = higcit_create_all_leaves(root);
	mp_assign_from_celliterator(m, it, 0);
	higcit_destroy(it);

	sim_domain *sd = sd_create(m);
	sd_add_higtree(sd, root);

	Point ls[4], hs[4];
	int ncs[4][DIM];
	hig_cell *bcs[4];
	mp_mapper *ms[4];
	sim_boundary *sbs[4];
	bc_type types[4];

	POINT_ASSIGN_REALS(ls[0], 0.0 - EPSDELTA, 0.0);
	POINT_ASSIGN_REALS(hs[0], 0.0 + EPSDELTA, 1.0);
	POINT_ASSIGN_INTS(ncs[0], 1, 5);
	types[0] = DIRICHLET;

	POINT_ASSIGN_REALS(ls[1], 1.0 - EPSDELTA, 0.0);
	POINT_ASSIGN_REALS(hs[1], 1.0 + EPSDELTA, 1.0);
	POINT_ASSIGN_INTS(ncs[1], 1, 5);
	types[1] = DIRICHLET;

	POINT_ASSIGN_REALS(ls[2], 0.0, 0.0 - EPSDELTA);
	POINT_ASSIGN_REALS(hs[2], 1.0, 0.0 + EPSDELTA);
	POINT_ASSIGN_INTS(ncs[2], 5, 1);
	types[2] = DIRICHLET;

	POINT_ASSIGN_REALS(ls[3], 0.0, 1.0 - EPSDELTA);
	POINT_ASSIGN_REALS(hs[3], 1.0, 1.0 + EPSDELTA);
	POINT_ASSIGN_INTS(ncs[3], 5, 1);
	types[3] = NEUMANN;

	for(int i = 0; i < 4; i++) {
		bcs[i] = hig_create_root(ls[i], hs[i]);
		hig_refine_uniform(bcs[i], ncs[i]);
		ms[i] = mp_create();
		it = higcit_create_all_leaves(bcs[i]);
		mp_assign_from_celliterator(ms[i], it, 0);
		higcit_destroy(it);
		sbs[i] = sb_create(bcs[i], types[i], ms[i]);
		sd_add_boundary(sd, sbs[i]);
	}

	sim_stencil *stn = stn_create();
	Point origin, x;
	int *ids;
	real *vals;
	int numelems;

	POINT_ASSIGN_REALS(origin, 0.5, 0.5);
	POINT_ASSIGN_REALS(x, 0.3, 0.5);
	stn_reset(stn);
	sd_set_interpolator_order(sd, 3);
	sd_get_stencil(sd, origin, x, 1.0, stn);
	numelems = stn_get_numelems(stn);
	ids = stn_get_ids(stn);
	vals = stn_get_vals(stn);
	for(int j = 0; j < numelems; j++) {
		printf("[%d] = %g\n", ids[j], vals[j]);
	}
	printf("rhs = %g\n\n", stn_get_rhs(stn));

	POINT_ASSIGN_REALS(origin, 0.1, 0.5);
	POINT_ASSIGN_REALS(x, -0.1, 0.5);
	stn_reset(stn);
	real bcvals[] = {1.0, 10.0, 100.0, 1000.0, 10000.0};
	sb_set_values(sbs[0], bcvals);
	sb_set_values(sbs[3], bcvals);
	sd_get_stencil(sd, origin, x, 1.0, stn);
	numelems = stn_get_numelems(stn);
	ids = stn_get_ids(stn);
	vals = stn_get_vals(stn);
	for(int j = 0; j < numelems; j++) {
		printf("[%d] = %g\n", ids[j], vals[j]);
	}
	printf("rhs = %g\n\n", stn_get_rhs(stn));

	POINT_ASSIGN_REALS(origin, 0.5, 0.9);
	POINT_ASSIGN_REALS(x, 0.5, 1.1);
	stn_reset(stn);
	sd_get_stencil(sd, origin, x, 1.0, stn);
	numelems = stn_get_numelems(stn);
	ids = stn_get_ids(stn);
	vals = stn_get_vals(stn);
	for(int j = 0; j < numelems; j++) {
		printf("[%d] = %g\n", ids[j], vals[j]);
	}
	printf("rhs = %g\n\n", stn_get_rhs(stn));

	POINT_ASSIGN_REALS(origin, 0.9, 0.1);
	POINT_ASSIGN_REALS(x, 1.1, 0.1);
	stn_reset(stn);
	sd_get_stencil(sd, origin, x, 1.0, stn);
	numelems = stn_get_numelems(stn);
	ids = stn_get_ids(stn);
	vals = stn_get_vals(stn);
	for(int j = 0; j < numelems; j++) {
		printf("[%d] = %g\n", ids[j], vals[j]);
	}
	printf("rhs = %g\n\n", stn_get_rhs(stn));
	sd_destroy(sd);
	stn_destroy(stn);

	return 0;
}
