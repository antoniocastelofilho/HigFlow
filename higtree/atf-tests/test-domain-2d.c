#include <stdlib.h>
#include <string.h>
#include <atf-c.h>
#include "domain.h"
#include "utils.h"
#define DEBUG
#include "Debug-c.h"

static void __refine_grid(hig_cell *root, int nc[DIM]) {
	int numChildren = hig_get_number_of_children(root);
	for(int i = 0; i < numChildren; i++) {
		int pp[DIM];
		hig_tobase(i, root->numcells, pp);
		if ((pp[0] + pp[1])%2==0) {
			hig_refine_uniform(hig_get_child(root, i), nc);
		}
	}
}

ATF_TC(domain_boundary_condition);
ATF_TC_HEAD(domain_boundary_condition, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(domain_boundary_condition, tc)
{
	Point l, h;
	POINT_ASSIGN_REALS(l, -1.0, -1.0);
	POINT_ASSIGN_REALS(h, 1.0, 1.0);
	hig_cell *root = hig_create_root(l, h);
	int nc[2];
	POINT_ASSIGN_INTS(nc, 4, 1);
	hig_refine_uniform(root, nc);
	mp_mapper *m = mp_create();
	higcit_celliterator *it = higcit_create_all_leaves(root);
	int size = mp_assign_from_celliterator(m, it, 0);
	higcit_destroy(it);

	sim_boundary *sb = sb_create(root, DIRICHLET, m);

	sb_set_value(sb, 0, 5.0);
	ATF_CHECK(FLT_EQ(sb_get_value(sb, 0), 5.0));
	real vals[4] = {1.0, 2.0, 3.0, 4.0};
	sb_set_values(sb, vals);
	ATF_CHECK(FLT_EQ(sb_get_value(sb, 0), 1.0));
	ATF_CHECK(FLT_EQ(sb_get_value(sb, 1), 2.0));
	ATF_CHECK(FLT_EQ(sb_get_value(sb, 2), 3.0));
	ATF_CHECK(FLT_EQ(sb_get_value(sb, 3), 4.0));
	
	sb_destroy(sb);
}

ATF_TC(domain);
ATF_TC_HEAD(domain, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(domain, tc)
{
	Point l, h;
	POINT_ASSIGN_REALS(l, -1.0, -1.0);
	POINT_ASSIGN_REALS(h, 1.0, 1.0);
	hig_cell *root = hig_create_root(l, h);
	int nc[2];
	int side = 5;
	POINT_ASSIGN_INTS(nc, side, side);
	hig_refine_uniform(root, nc);
	mp_mapper *m = mp_create();
	higcit_celliterator *it = higcit_create_all_leaves(root);
	int size = mp_assign_from_celliterator(m, it, 0);
	higcit_destroy(it);

	sim_domain *sd = sd_create(m);
	sd_add_higtree(sd, root);
	ATF_CHECK(sd_get_higtree(sd, 0) == root);
	ATF_CHECK(sd_get_domain_mapper(sd) == m);
	it = sd_get_domain_celliterator(sd);
	ATF_CHECK(higcit_count(it) == side*side);
	higcit_destroy(it);

	Point ls[4], hs[4];
	int ncs[4][DIM];
	hig_cell *bcs[4];
	mp_mapper *ms[4];
	sim_boundary *sbs[4];
	bc_type types[4];

	POINT_ASSIGN_REALS(ls[0], l[0] - EPSDELTA, l[1]);
	POINT_ASSIGN_REALS(hs[0], l[0] + EPSDELTA, h[1]);
	POINT_ASSIGN_INTS(ncs[0], 1, side);
	types[0] = NEUMANN;

	POINT_ASSIGN_REALS(ls[1], h[0] - EPSDELTA, l[1]);
	POINT_ASSIGN_REALS(hs[1], h[0] + EPSDELTA, h[1]);
	POINT_ASSIGN_INTS(ncs[1], 1, side);
	types[1] = DIRICHLET;

	POINT_ASSIGN_REALS(ls[2], l[0], l[1] - EPSDELTA);
	POINT_ASSIGN_REALS(hs[2], h[0], l[1] + EPSDELTA);
	POINT_ASSIGN_INTS(ncs[2], side, 1);
	types[2] = NEUMANN;

	POINT_ASSIGN_REALS(ls[3], l[0], h[1] - EPSDELTA);
	POINT_ASSIGN_REALS(hs[3], h[0], h[1] + EPSDELTA);
	POINT_ASSIGN_INTS(ncs[3], side, 1);
	types[3] = DIRICHLET;

	for(int i = 0; i < 4; i++) {
		bcs[i] = hig_create_root(ls[i], hs[i]);
		hig_refine_uniform(bcs[i], ncs[i]);
		ms[i] = mp_create();
		it = higcit_create_all_leaves(bcs[i]);
		mp_assign_from_celliterator(ms[i], it, 0);
		higcit_destroy(it);
		sbs[i] = sb_create(bcs[i], types[i], ms[i]);
		for(it = higcit_create_all_leaves(bcs[i]); !higcit_isfinished(it); higcit_nextcell(it)) {
			hig_cell *c = higcit_getcell(it);
			int cgid = mp_lookup(ms[i], hig_get_id(c, 0));
			sb_set_value(sbs[i], cgid, 10.0 + i);
		}
		higcit_destroy(it);
		sd_add_boundary(sd, sbs[i]);
	}

	sim_stencil *stn = stn_create();

	sd_set_interpolator_order(sd, 3);
	Point origin, x;
	POINT_ASSIGN_REALS(origin, 0.4, 0.0);
	POINT_ASSIGN_REALS(x, 0.0, 0.0);

	sd_set_interpolator_center(sd, origin);
	sd_get_stencil(sd, x, 1.0, stn);

	ATF_CHECK(stn_get_numelems(stn) == 1);
	ATF_CHECK(FLT_EQ(stn_get_rhs(stn), 0.0));
	stn_reset(stn);

	POINT_ASSIGN_REALS(origin, 0.8, 0.0);
	POINT_ASSIGN_REALS(x, 1.2, 0.0);

	sd_set_interpolator_center(sd, origin);
	sd_get_stencil(sd, x, 1.0, stn);
	ATF_CHECK(stn_get_numelems(stn) >= 3);
	ATF_CHECK(FLT_NE(stn_get_rhs(stn), 0.0));
	stn_reset(stn);

	POINT_ASSIGN_REALS(origin, -0.8, 0.0);
	POINT_ASSIGN_REALS(x, -1.2, 0.0);

	sd_set_interpolator_center(sd, origin);
	sd_get_stencil(sd, x, 1.0, stn);
	ATF_CHECK(stn_get_numelems(stn) >= 1);
	ATF_CHECK(FLT_NE(stn_get_rhs(stn), 0.0));
	stn_reset(stn);

	POINT_ASSIGN_REALS(origin, 0.0, 0.8);
	POINT_ASSIGN_REALS(x, 0.0, 1.2);

	sd_set_interpolator_center(sd, origin);
	sd_get_stencil(sd, x, 1.0, stn);
	ATF_CHECK(stn_get_numelems(stn) >= 3);
	ATF_CHECK(FLT_NE(stn_get_rhs(stn), 0.0));
	stn_reset(stn);

	POINT_ASSIGN_REALS(origin, 0.0, -0.8);
	POINT_ASSIGN_REALS(x, 0.0, -1.2);

	sd_set_interpolator_center(sd, origin);
	sd_get_stencil(sd, x, 1.0, stn);
	ATF_CHECK(stn_get_numelems(stn) >= 1);
	ATF_CHECK(FLT_NE(stn_get_rhs(stn), 0.0));
	stn_reset(stn);

	sd_set_interpolator_center(sd, origin);
	sd_get_stencil_without_bc(sd, x, 1.0, 1.0, stn);
	ATF_CHECK(stn_get_numelems(stn) >= 3);
	ATF_CHECK(FLT_EQ(stn_get_rhs(stn), 0.0));
	stn_reset(stn);

	it = sd_get_bcs_celliterator(sd);
	ATF_CHECK(higcit_count(it) == 4*side);
	higcit_destroy(it);

	stn_destroy(stn);

	sd_destroy(sd);
}

ATF_TC(domain_stencil);
ATF_TC_HEAD(domain_stencil, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(domain_stencil, tc)
{
	sim_stencil *stn = stn_create();
	
	stn_set_rhs(stn, 5.0);
	ATF_CHECK(FLT_EQ(stn_get_rhs(stn), 5.0));
	stn_add_to_rhs(stn, 5.0);
	ATF_CHECK(FLT_EQ(stn_get_rhs(stn), 10.0));

	stn_add_to_element(stn, 5, 5.0);
	ATF_CHECK(stn_get_numelems(stn) == 1);
	ATF_CHECK(stn_get_id(stn, 0) == 5);
	ATF_CHECK(FLT_EQ(stn_get_val(stn, 0), 5.0));
	stn_add_to_element(stn, 5, 3.0);
	ATF_CHECK(FLT_EQ(stn_get_val(stn, 0), 8.0));
	stn_add_to_element(stn, 10, 3.0);
	ATF_CHECK(FLT_EQ(stn_get_val(stn, 0), 8.0));
	ATF_CHECK(FLT_EQ(stn_get_val(stn, 1), 3.0));

	stn_set_element(stn, 10, 10.0);
	stn_set_element(stn, 11, 15.0);
	ATF_CHECK(stn_get_numelems(stn) == 3);
	int *ids = stn_get_ids(stn);
	real *vals = stn_get_vals(stn);
	ATF_CHECK(ids[0] == 5);
	ATF_CHECK(ids[1] == 10);
	ATF_CHECK(ids[2] == 11);
	ATF_CHECK(FLT_EQ(vals[0], 8.0));
	ATF_CHECK(FLT_EQ(vals[1], 10.0));
	ATF_CHECK(FLT_EQ(vals[2], 15.0));

	stn_reset(stn);
	ATF_CHECK(stn_get_numelems(stn) == 0);
	ATF_CHECK(FLT_EQ(stn_get_rhs(stn), 0.0));

	for(int i = 0; i < 101; i++) {
		stn_add_to_element(stn, i, i);
	}
	ATF_CHECK(stn_get_numelems(stn) == 101);
	
	sim_stencil *stn2 = stn_create();
	for(int i = 95; i < 110; i++) {
		stn_add_to_element(stn2, i, i);
	}
	stn_add_stencil(stn, stn2);
	stn_destroy(stn2);
	ATF_CHECK(stn_get_numelems(stn) == 110);
	ATF_CHECK(FLT_EQ(stn_get_val(stn, 1), 1.0));
	ATF_CHECK(FLT_EQ(stn_get_val(stn, 100), 200.0));

	stn_reset(stn);
	stn_set_element(stn, 0, 10.0);
	stn_set_element(stn, 1, 11.0);
	stn_set_element(stn, 2, 12.0);
	real v[3] = {1.0, -1.0, 0.5};
	ATF_CHECK(FLT_EQ(stn_mult_vector(stn, v), 5.0));

	stn_destroy(stn);
}

ATF_TP_ADD_TCS(tp)
{
	ATF_TP_ADD_TC(tp, domain);
	ATF_TP_ADD_TC(tp, domain_stencil);
	ATF_TP_ADD_TC(tp, domain_boundary_condition);
}
