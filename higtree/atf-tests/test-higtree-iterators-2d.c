#include <atf-c.h>
#include "higtree.h"
#include "higtree-iterator.h"
#include "utils.h"
#define DEBUG
#include "Debug-c.h"

static void __refine_grid_two_levels(hig_cell *root, int nc[DIM]) {
	int numChildren = hig_get_number_of_children(root);
	for(int i = 0; i < numChildren; i++) {
		int pp[DIM];
		hig_tobase(i, root->numcells, pp);
		if ((pp[0] + pp[1])%2==0) {
		  hig_refine_uniform(hig_get_child(root, i), nc);
		}
	}
}

void __check_cell(int *index, hig_cell *root, higcit_celliterator *it) {
	while (*index != -2) {
		hig_cell *c = root;
		while (*index != -1) {
			int i = *index;
			c = hig_get_child(c, i);
			index++;
		}
		hig_cell *p = c;
		DEBUG_PASS;
		while (p->posinparent != -1) {
			DEBUG_INSPECT(p->posinparent, %d);
			p = p->parent;
		}
		DEBUG_INSPECT(p->posinparent, %d);
		ATF_CHECK(c == higcit_getcell(it));
		higcit_nextcell(it);
		index++;
	}
	ATF_CHECK(higcit_isfinished(it));
}

ATF_TC(all_cells_one_level);
ATF_TC_HEAD(all_cells_one_level, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(all_cells_one_level, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 5);
	hig_refine_uniform(root, nc);
	higcit_celliterator *it;
	int cnt = 0;
	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		ATF_CHECK(c != NULL);
		ATF_CHECK(c == hig_get_child(root, cnt));
		cnt++;
	}
	higcit_destroy(it);
	hig_destroy(root);
}

ATF_TC(all_cells_two_levels);
ATF_TC_HEAD(all_cells_two_levels, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(all_cells_two_levels, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 3);
	hig_refine_uniform(root, nc);
	hig_refine_uniform(hig_get_child(root, 0), nc);
	hig_refine_uniform(hig_get_child(root, 2), nc);
	hig_refine_uniform(hig_get_child(root, 4), nc);
	higcit_celliterator *it;
	it = higcit_create_all_leaves(root);
	int pos[] = {
		0,0,-1,
		0,1,-1,
		0,2,-1,
		0,3,-1,
		0,4,-1,
		0,5,-1,
		0,6,-1,
		0,7,-1,
		0,8,-1,
		1,-1,
		2,0,-1,
		2,1,-1,
		2,2,-1,
		2,3,-1,
		2,4,-1,
		2,5,-1,
		2,6,-1,
		2,7,-1,
		2,8,-1,
		3,-1,
		4,0,-1,
		4,1,-1,
		4,2,-1,
		4,3,-1,
		4,4,-1,
		4,5,-1,
		4,6,-1,
		4,7,-1,
		4,8,-1,
		5,-1,
		6,-1,
		7,-1,
		8,-1,
		-2,
	};
	__check_cell(pos, root, it);
	higcit_destroy(it);
	hig_destroy(root);
}

ATF_TC(all_cells_two_levels_clone);
ATF_TC_HEAD(all_cells_two_levels_clone, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(all_cells_two_levels_clone, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 3);
	hig_refine_uniform(root, nc);
	hig_refine_uniform(hig_get_child(root, 0), nc);
	hig_refine_uniform(hig_get_child(root, 2), nc);
	hig_refine_uniform(hig_get_child(root, 4), nc);
	higcit_celliterator *it;
	it = higcit_create_all_leaves(root);
	higcit_celliterator *itclone1 = higcit_clone(it);
	for(int i = 0; i < 6; i++) {
		higcit_nextcell(it);
	}
	higcit_celliterator *itclone2 = higcit_clone(it);
	for(int i = 0; i < 6; i++) {
		higcit_nextcell(it);
	}
	higcit_celliterator *itclone3 = higcit_clone(it);
	int pos[] = {
		0,0,-1,
		0,1,-1,
		0,2,-1,
		0,3,-1,
		0,4,-1,
		0,5,-1,
		0,6,-1,
		0,7,-1,
		0,8,-1,
		1,-1,
		2,0,-1,
		2,1,-1,
		2,2,-1,
		2,3,-1,
		2,4,-1,
		2,5,-1,
		2,6,-1,
		2,7,-1,
		2,8,-1,
		3,-1,
		4,0,-1,
		4,1,-1,
		4,2,-1,
		4,3,-1,
		4,4,-1,
		4,5,-1,
		4,6,-1,
		4,7,-1,
		4,8,-1,
		5,-1,
		6,-1,
		7,-1,
		8,-1,
		-2,
	};
	__check_cell(&pos[0], root, itclone1);
	__check_cell(&pos[18], root, itclone2);
	__check_cell(&pos[35], root, itclone3);
	higcit_destroy(it);
	higcit_destroy(itclone1);
	higcit_destroy(itclone2);
	higcit_destroy(itclone3);
	hig_destroy(root);
}

ATF_TC(bounding_box_two_levels1);
ATF_TC_HEAD(bounding_box_two_levels1, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(bounding_box_two_levels1, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid_two_levels(root, nc);
	lp[0] = -1.2;
	lp[1] = -1.2;
	hp[0] = -0.9;
	hp[1] = 0.1;
	higcit_celliterator *it;
	it = higcit_create_bounding_box(root, lp, hp);
	int pos[] = {
		0,0,-1,
		0,1,-1,
		1,-1,
		2,0,-1,
		-2,
	};
	__check_cell(pos, root, it);
	higcit_destroy(it);
	hig_destroy(root);
}

ATF_TC(bounding_box_two_levels2);
ATF_TC_HEAD(bounding_box_two_levels2, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(bounding_box_two_levels2, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid_two_levels(root, nc);
	lp[0] = -1.2;
	lp[1] = -1.2;
	hp[0] = -0.9;
	hp[1] = 0.1;
	higcit_celliterator *it;
	it = higcit_create_bounding_box(root, lp, hp);
	int pos[] = {
		0,0,-1,
		0,1,-1,
		1,-1,
		2,0,-1,
		-2
	};
	__check_cell(pos, root, it);
	higcit_destroy(it);
	hig_destroy(root);
}

ATF_TC(bounding_box_two_levels_outside);
ATF_TC_HEAD(bounding_box_two_levels_outside, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(bounding_box_two_levels_outside, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid_two_levels(root, nc);
	POINT_ASSIGN_SCALAR(lp, -1.4);
	POINT_ASSIGN_SCALAR(hp, -1.1);
	higcit_celliterator *it;
	it = higcit_create_bounding_box(root, lp, hp);
	ATF_CHECK(higcit_isfinished(it));
	higcit_destroy(it);
	hig_destroy(root);
}

ATF_TC(bounding_box_two_levels3);
ATF_TC_HEAD(bounding_box_two_levels3, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(bounding_box_two_levels3, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid_two_levels(root, nc);
	lp[0] = 0.4;
	lp[1] = -0.8;
	hp[0] = 0.6;
	hp[1] = 0.6;
	higcit_celliterator *it;
	hig_cell *c;
	it = higcit_create_bounding_box(root, lp, hp);
	int pos[] = {
		8,2,-1,
		8,3,-1,
		9,-1,
		10,2,-1,
		10,3,-1,
		11,-1,
		12,-1,
		13,0,-1,
		13,1,-1,
		14,-1,
		15,0,-1,
		-2
	};
	__check_cell(pos, root, it);
	higcit_destroy(it);
	hig_destroy(root);
}

ATF_TC(bounding_box_two_levels_clone);
ATF_TC_HEAD(bounding_box_two_levels_clone, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(bounding_box_two_levels_clone, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid_two_levels(root, nc);
	lp[0] = 0.4;
	lp[1] = -0.8;
	hp[0] = 0.6;
	hp[1] = 0.6;
	higcit_celliterator *it;
	hig_cell *c;
	it = higcit_create_bounding_box(root, lp, hp);
	higcit_celliterator *itclone1 = higcit_clone(it);
	for(int i = 0; i < 6; i++) {
		higcit_nextcell(it);
	}
	higcit_celliterator *itclone2 = higcit_clone(it);
	for(int i = 0; i < 3; i++) {
		higcit_nextcell(it);
	}
	higcit_celliterator *itclone3 = higcit_clone(it);
	int pos[] = {
		8,2,-1,
		8,3,-1,
		9,-1,
		10,2,-1,
		10,3,-1,
		11,-1,
		12,-1,
		13,0,-1,
		13,1,-1,
		14,-1,
		15,0,-1,
		-2
	};
	__check_cell(&pos[0], root, itclone1);
	__check_cell(&pos[16], root, itclone2);
	__check_cell(&pos[24], root, itclone3);
	higcit_destroy(it);
	higcit_destroy(itclone1);
	higcit_destroy(itclone2);
	higcit_destroy(itclone3);
	hig_destroy(root);
}


ATF_TC(bounding_box_two_levels_one_cell_with_border);
ATF_TC_HEAD(bounding_box_two_levels_one_cell_with_border, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(bounding_box_two_levels_one_cell_with_border, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid_two_levels(root, nc);
	lp[0] = 0.0;
	lp[1] = 0.0;
	hp[0] = 0.25;
	hp[1] = 0.25;
	higcit_celliterator *it;
	hig_cell *c;
	it = higcit_create_bounding_box(root, lp, hp);
	int pos[] = {
		5,3,-1,
		6,-1,
		9,-1,
		10,0,-1,
		-2
	};
	__check_cell(pos, root, it);
	higcit_destroy(it);
	hig_destroy(root);
}

ATF_TC(bounding_box_two_levels_one_point);
ATF_TC_HEAD(bounding_box_two_levels_one_point, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(bounding_box_two_levels_one_point, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid_two_levels(root, nc);
	higcit_celliterator *it;
	hig_cell *c;
	lp[0] = 0.1;
	lp[1] = 0.1;
	hp[0] = 0.1;
	hp[1] = 0.1;
	it = higcit_create_bounding_box(root, lp, hp);
	c = higcit_getcell(it);
	ATF_CHECK(c == hig_get_child(hig_get_child(root, 10), 0));
	higcit_nextcell(it);
	ATF_CHECK(higcit_isfinished(it));
	higcit_destroy(it);

	lp[0] = 0.0;
	lp[1] = 0.0;
	hp[0] = 0.0;
	hp[1] = 0.0;
	it = higcit_create_bounding_box(root, lp, hp);
	c = higcit_getcell(it);
	ATF_CHECK(c == hig_get_child(hig_get_child(root, 5), 3));
	higcit_nextcell(it);
	ATF_CHECK(higcit_isfinished(it));
	higcit_destroy(it);

	lp[0] = 1.0;
	lp[1] = 1.0;
	hp[0] = 1.0;
	hp[1] = 1.0;
	it = higcit_create_bounding_box(root, lp, hp);
	c = higcit_getcell(it);
	ATF_CHECK(c == hig_get_child(hig_get_child(root, 15), 3));
	higcit_nextcell(it);
	ATF_CHECK(higcit_isfinished(it));
	higcit_destroy(it);

	lp[0] = -1.0;
	lp[1] = -1.0;
	hp[0] = -1.0;
	hp[1] = -1.0;
	it = higcit_create_bounding_box(root, lp, hp);
	c = higcit_getcell(it);
	ATF_CHECK(c == hig_get_child(hig_get_child(root, 0), 0));
	higcit_nextcell(it);
	ATF_CHECK(higcit_isfinished(it));
	higcit_destroy(it);

	lp[0] = -1.0;
	lp[1] = 1.0;
	hp[0] = -1.0;
	hp[1] = 1.0;
	it = higcit_create_bounding_box(root, lp, hp);
	c = higcit_getcell(it);
	ATF_CHECK(c == hig_get_child(root, 3));
	higcit_nextcell(it);
	ATF_CHECK(higcit_isfinished(it));
	higcit_destroy(it);

	hig_destroy(root);
}

ATF_TC(neighbours_two_levels);
ATF_TC_HEAD(neighbours_two_levels, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(neighbours_two_levels, tc)
{
	higcit_celliterator *it;
	hig_cell *c;
	hig_cell *n;
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);

	it = higcit_create_neighbours(root);
	int pos6[] = {
		-2,
	};
	__check_cell(pos6, root, it);
	higcit_destroy(it);

	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid_two_levels(root, nc);
	c = hig_get_child(hig_get_child(root, 5), 3);
	it = higcit_create_neighbours(c);
	int pos[] = {
		5,0,-1,
		5,1,-1,
		5,2,-1,
		6,-1,
		9,-1,
		10,0,-1,
		-2
	};
	__check_cell(pos, root, it);
	higcit_destroy(it);

	c = hig_get_child(root, 9);
	it = higcit_create_neighbours(c);
	int pos2[] = {
		4,-1,
		5,2,-1,
		5,3,-1,
		6,-1,
		8,1,-1,
		8,3,-1,
		10,0,-1,
		10,2,-1,
		12,-1,
		13,0,-1,
		13,1,-1,
		14,-1,
		-2
	};
	__check_cell(pos2, root, it);
	higcit_destroy(it);

	c = hig_get_child(hig_get_child(root, 15), 3);
	it = higcit_create_neighbours(c);
	int pos3[] = {
		15,0,-1,
		15,1,-1,
		15,2,-1,
		-2
	};
	__check_cell(pos3, root, it);
	higcit_destroy(it);

	c = hig_get_child(hig_get_child(root, 0), 0);
	it = higcit_create_neighbours(c);
	int pos4[] = {
		0,1,-1,
		0,2,-1,
		0,3,-1,
		-2
	};
	__check_cell(pos4, root, it);
	higcit_destroy(it);

	c = hig_get_child(root, 3);
	it = higcit_create_neighbours(c);
	int pos5[] = {
		2,1,-1,
		2,3,-1,
		6,-1,
		7,0,-1,
		7,1,-1,
		-2
	};
	__check_cell(pos5, root, it);
	higcit_destroy(it);

	hig_destroy(root);
}

ATF_TC(neighbours_two_levels_clone);
ATF_TC_HEAD(neighbours_two_levels_clone, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(neighbours_two_levels_clone, tc)
{
	hig_cell *c;
	hig_cell *n;
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);

	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid_two_levels(root, nc);

	c = hig_get_child(root, 9);
	higcit_celliterator *it = higcit_create_neighbours(c);
	higcit_celliterator *itclone1 = higcit_clone(it);
	for(int i = 0; i < 3; i++) {
		higcit_nextcell(it);
	}
	higcit_celliterator *itclone2 = higcit_clone(it);
	for(int i = 0; i < 3; i++) {
		higcit_nextcell(it);
	}
	higcit_celliterator *itclone3 = higcit_clone(it);
	int pos[] = {
		4,-1,
		5,2,-1,
		5,3,-1,
		6,-1,
		8,1,-1,
		8,3,-1,
		10,0,-1,
		10,2,-1,
		12,-1,
		13,0,-1,
		13,1,-1,
		14,-1,
		-2
	};
	__check_cell(&pos[0], root, itclone1);
	__check_cell(&pos[8], root, itclone2);
	__check_cell(&pos[16], root, itclone3);
	higcit_destroy(it);
	higcit_destroy(itclone1);
	higcit_destroy(itclone2);
	higcit_destroy(itclone3);

	hig_destroy(root);
}

boolean __max_delta(hig_cell *c, void *d) {
	real maxd = *(real *) d;
	Point delta;
	hig_get_delta(c, delta);

	return FLT_LE(delta[0], maxd) && FLT_LE(delta[1], maxd);
}

ATF_TC(filter_max_delta_two_levels);
ATF_TC_HEAD(filter_max_delta_two_levels, tc)
{
    atf_tc_set_md_var(tc, "descr", "Testing Filter: cells of a maximum delta");
}
ATF_TC_BODY(filter_max_delta_two_levels, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid_two_levels(root, nc);
	higcit_celliterator *it = higcit_create_all_leaves(root);
	real maxd = 0.3;

	higcit_celliterator *itf = higcit_create_filter(it, __max_delta, &maxd);

	int pos[] = {
		0,0,-1,
		0,1,-1,
		0,2,-1,
		0,3,-1,
		2,0,-1,
		2,1,-1,
		2,2,-1,
		2,3,-1,
		5,0,-1,
		5,1,-1,
		5,2,-1,
		5,3,-1,
		7,0,-1,
		7,1,-1,
		7,2,-1,
		7,3,-1,
		8,0,-1,
		8,1,-1,
		8,2,-1,
		8,3,-1,
		10,0,-1,
		10,1,-1,
		10,2,-1,
		10,3,-1,
		13,0,-1,
		13,1,-1,
		13,2,-1,
		13,3,-1,
		15,0,-1,
		15,1,-1,
		15,2,-1,
		15,3,-1,
		-2
	};
	__check_cell(pos, root, itf);
	higcit_destroy(itf);

	hig_destroy(root);
}

ATF_TC(filter_max_delta_two_levels_clone);
ATF_TC_HEAD(filter_max_delta_two_levels_clone, tc)
{
    atf_tc_set_md_var(tc, "descr", "Testing Filter: cells of a maximum delta");
}
ATF_TC_BODY(filter_max_delta_two_levels_clone, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid_two_levels(root, nc);
	higcit_celliterator *it = higcit_create_all_leaves(root);
	real maxd = 0.3;

	higcit_celliterator *itf = higcit_create_filter(it, __max_delta, &maxd);
	higcit_celliterator *itclone1 = higcit_clone(itf);
	for(int i = 0; i < 3; i++) {
		higcit_nextcell(itf);
	}
	higcit_celliterator *itclone2 = higcit_clone(itf);
	for(int i = 0; i < 3; i++) {
		higcit_nextcell(itf);
	}
	higcit_celliterator *itclone3 = higcit_clone(itf);

	int pos[] = {
		0,0,-1,
		0,1,-1,
		0,2,-1,
		0,3,-1,
		2,0,-1,
		2,1,-1,
		2,2,-1,
		2,3,-1,
		5,0,-1,
		5,1,-1,
		5,2,-1,
		5,3,-1,
		7,0,-1,
		7,1,-1,
		7,2,-1,
		7,3,-1,
		8,0,-1,
		8,1,-1,
		8,2,-1,
		8,3,-1,
		10,0,-1,
		10,1,-1,
		10,2,-1,
		10,3,-1,
		13,0,-1,
		13,1,-1,
		13,2,-1,
		13,3,-1,
		15,0,-1,
		15,1,-1,
		15,2,-1,
		15,3,-1,
		-2
	};
	__check_cell(&pos[0], root, itclone1);
	__check_cell(&pos[9], root, itclone2);
	__check_cell(&pos[18], root, itclone3);
	higcit_destroy(itf);
	higcit_destroy(itclone1);
	higcit_destroy(itclone2);
	higcit_destroy(itclone3);

	hig_destroy(root);
}

boolean __close_to_origin(hig_cell *c, void *d) {
	real dist = *(real *) d;
	Point center;
	hig_get_center(c, center);

	return FLT_LE(fabs(center[0]), dist) && FLT_LE(fabs(center[1]), dist);
}

ATF_TC(filter_close_to_origin_two_levels);
ATF_TC_HEAD(filter_close_to_origin_two_levels, tc)
{
    atf_tc_set_md_var(tc, "descr", "Testing Filter: close to origin");
}
ATF_TC_BODY(filter_close_to_origin_two_levels, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid_two_levels(root, nc);
	higcit_celliterator *it = higcit_create_all_leaves(root);
	real dist = 0.5;

	higcit_celliterator *itf = higcit_create_filter(it, __close_to_origin, &dist);

	int pos[] = {
		5,0,-1,
		5,1,-1,
		5,2,-1,
		5,3,-1,
		6,-1,
		9,-1,
		10,0,-1,
		10,1,-1,
		10,2,-1,
		10,3,-1,
		-2,
	};
	__check_cell(pos, root, itf);
	higcit_destroy(itf);

	hig_destroy(root);
}

ATF_TC(filter_close_to_origin_and_max_delta_two_levels);
ATF_TC_HEAD(filter_close_to_origin_and_max_delta_two_levels, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(filter_close_to_origin_and_max_delta_two_levels, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid_two_levels(root, nc);
	higcit_celliterator *it = higcit_create_all_leaves(root);
	real dist = 0.5;
	real maxd = 0.3;

	higcit_celliterator *itf1 = higcit_create_filter(it, __close_to_origin, &dist);
	higcit_celliterator *itf2 = higcit_create_filter(itf1, __max_delta, &maxd);

	int pos[] = {
		5,0,-1,
		5,1,-1,
		5,2,-1,
		5,3,-1,
		10,0,-1,
		10,1,-1,
		10,2,-1,
		10,3,-1,
		-2
	};
	__check_cell(pos, root, itf2);

	higcit_destroy(itf2);

	hig_destroy(root);
}

boolean __always_false(hig_cell *c, void *d) {
	return 0;
}

ATF_TC(filter_emtpy_two_levels);
ATF_TC_HEAD(filter_emtpy_two_levels, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(filter_emtpy_two_levels, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid_two_levels(root, nc);
	higcit_celliterator *it = higcit_create_all_leaves(root);

	higcit_celliterator *itf1 = higcit_create_filter(it, __always_false, (void *) NULL);

	int pos[] = {
		-2
	};
	__check_cell(pos, root, itf1);

	higcit_destroy(itf1);

	hig_destroy(root);
}
boolean __visible(hig_cell *c, void *d) {
	real maxd = *(real *)d;
	Point delta;
	hig_get_delta(c, delta);
	return (FLT_GT(delta[0],maxd) && (hig_get_number_of_children(c) == 0 || FLT_LT(delta[0]/2.0,maxd)));
}

ATF_TC(filter_multi_grid_four_levels);
ATF_TC_HEAD(filter_multi_grid_four_levels, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(filter_multi_grid_four_levels, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	assert(root != NULL);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	DEBUG_PASS;
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid_two_levels(root, nc);
	hig_refine_uniform(hig_get_child(hig_get_child(root, 10), 0), nc);
	hig_refine_uniform(hig_get_child(hig_get_child(root, 13), 3), nc);
	hig_refine_uniform(hig_get_child(hig_get_child(hig_get_child(root, 13), 3), 1), nc);
	DEBUG_PASS;
	higcit_celliterator *it = higcit_create_all_higtree(root);
	real maxd = 0.06;
	higcit_celliterator *itf = higcit_create_filter(it, __visible, &maxd);

	int pos1[] = {
		0,0,-1,
		0,1,-1,
		0,2,-1,
		0,3,-1,
		1,-1,
		2,0,-1,
		2,1,-1,
		2,2,-1,
		2,3,-1,
		3,-1,
		4,-1,
		5,0,-1,
		5,1,-1,
		5,2,-1,
		5,3,-1,
		6,-1,
		7,0,-1,
		7,1,-1,
		7,2,-1,
		7,3,-1,
		8,0,-1,
		8,1,-1,
		8,2,-1,
		8,3,-1,
		9,-1,
		10,0,0,-1,
		10,0,1,-1,
		10,0,2,-1,
		10,0,3,-1,
		10,1,-1,
		10,2,-1,
		10,3,-1,
		11,-1,
		12,-1,
		13,0,-1,
		13,1,-1,
		13,2,-1,
		13,3,0,-1,
		13,3,1,0,-1,
		13,3,1,1,-1,
		13,3,1,2,-1,
		13,3,1,3,-1,
		13,3,2,-1,
		13,3,3,-1,
		14,-1,
		15,0,-1,
		15,1,-1,
		15,2,-1,
		15,3,-1,
		-2,
	};
	__check_cell(pos1, root, itf);

	higcit_destroy(itf);

	it = higcit_create_all_higtree(root);
	maxd *= 2.0;
	itf = higcit_create_filter(it, __visible, &maxd);

	int pos2[] = {
		0,0,-1,
		0,1,-1,
		0,2,-1,
		0,3,-1,
		1,-1,
		2,0,-1,
		2,1,-1,
		2,2,-1,
		2,3,-1,
		3,-1,
		4,-1,
		5,0,-1,
		5,1,-1,
		5,2,-1,
		5,3,-1,
		6,-1,
		7,0,-1,
		7,1,-1,
		7,2,-1,
		7,3,-1,
		8,0,-1,
		8,1,-1,
		8,2,-1,
		8,3,-1,
		9,-1,
		10,0,0,-1,
		10,0,1,-1,
		10,0,2,-1,
		10,0,3,-1,
		10,1,-1,
		10,2,-1,
		10,3,-1,
		11,-1,
		12,-1,
		13,0,-1,
		13,1,-1,
		13,2,-1,
		13,3,0,-1,
		13,3,1,-1,
		13,3,2,-1,
		13,3,3,-1,
		14,-1,
		15,0,-1,
		15,1,-1,
		15,2,-1,
		15,3,-1,
		-2,
	};
	__check_cell(pos2, root, itf);

	higcit_destroy(itf);


	it = higcit_create_all_higtree(root);
	maxd *= 2.0;
	itf = higcit_create_filter(it, __visible, &maxd);

	int pos3[] = {
		0,0,-1,
		0,1,-1,
		0,2,-1,
		0,3,-1,
		1,-1,
		2,0,-1,
		2,1,-1,
		2,2,-1,
		2,3,-1,
		3,-1,
		4,-1,
		5,0,-1,
		5,1,-1,
		5,2,-1,
		5,3,-1,
		6,-1,
		7,0,-1,
		7,1,-1,
		7,2,-1,
		7,3,-1,
		8,0,-1,
		8,1,-1,
		8,2,-1,
		8,3,-1,
		9,-1,
		10,0,-1,
		10,1,-1,
		10,2,-1,
		10,3,-1,
		11,-1,
		12,-1,
		13,0,-1,
		13,1,-1,
		13,2,-1,
		13,3,-1,
		14,-1,
		15,0,-1,
		15,1,-1,
		15,2,-1,
		15,3,-1,
		-2,
	};
	__check_cell(pos3, root, itf);

	higcit_destroy(itf);

	it = higcit_create_all_higtree(root);
	maxd *= 2.0;
	itf = higcit_create_filter(it, __visible, &maxd);

	int pos4[] = {
		0,-1,
		1,-1,
		2,-1,
		3,-1,
		4,-1,
		5,-1,
		6,-1,
		7,-1,
		8,-1,
		9,-1,
		10,-1,
		11,-1,
		12,-1,
		13,-1,
		14,-1,
		15,-1,
		-2,
	};
	__check_cell(pos4, root, itf);

	higcit_destroy(itf);

	it = higcit_create_all_higtree(root);
	maxd *= 2.0;
	itf = higcit_create_filter(it, __visible, &maxd);

	int pos5[] = {
		-2,
	};
	__check_cell(pos5, root, itf);

	higcit_destroy(itf);

	hig_destroy(root);
}

ATF_TC(all_higtree_two_levels);
ATF_TC_HEAD(all_higtree_two_levels, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(all_higtree_two_levels, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	assert(root != NULL);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	DEBUG_PASS;
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid_two_levels(root, nc);
	DEBUG_PASS;
	higcit_celliterator *it = higcit_create_all_higtree(root);

	int pos[] = {
		-1,
		0,-1,
		0,0,-1,
		0,1,-1,
		0,2,-1,
		0,3,-1,
		1,-1,
		2,-1,
		2,0,-1,
		2,1,-1,
		2,2,-1,
		2,3,-1,
		3,-1,
		4,-1,
		5,-1,
		5,0,-1,
		5,1,-1,
		5,2,-1,
		5,3,-1,
		6,-1,
		7,-1,
		7,0,-1,
		7,1,-1,
		7,2,-1,
		7,3,-1,
		8,-1,
		8,0,-1,
		8,1,-1,
		8,2,-1,
		8,3,-1,
		9,-1,
		10,-1,
		10,0,-1,
		10,1,-1,
		10,2,-1,
		10,3,-1,
		11,-1,
		12,-1,
		13,-1,
		13,0,-1,
		13,1,-1,
		13,2,-1,
		13,3,-1,
		14,-1,
		15,-1,
		15,0,-1,
		15,1,-1,
		15,2,-1,
		15,3,-1,
		-2
	};
	__check_cell(pos, root, it);

	higcit_destroy(it);

	hig_destroy(root);
}

ATF_TC(all_higtree_two_levels_clone);
ATF_TC_HEAD(all_higtree_two_levels_clone, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(all_higtree_two_levels_clone, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	assert(root != NULL);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	DEBUG_PASS;
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid_two_levels(root, nc);
	DEBUG_PASS;
	higcit_celliterator *it = higcit_create_all_higtree(root);
	higcit_celliterator *itclone1 = higcit_clone(it);
	for(int i = 0; i < 3; i++) {
		higcit_nextcell(it);
	}
	higcit_celliterator *itclone2 = higcit_clone(it);
	for(int i = 0; i < 3; i++) {
		higcit_nextcell(it);
	}
	higcit_celliterator *itclone3 = higcit_clone(it);

	int pos[] = {
		-1,
		0,-1,
		0,0,-1,
		0,1,-1,
		0,2,-1,
		0,3,-1,
		1,-1,
		2,-1,
		2,0,-1,
		2,1,-1,
		2,2,-1,
		2,3,-1,
		3,-1,
		4,-1,
		5,-1,
		5,0,-1,
		5,1,-1,
		5,2,-1,
		5,3,-1,
		6,-1,
		7,-1,
		7,0,-1,
		7,1,-1,
		7,2,-1,
		7,3,-1,
		8,-1,
		8,0,-1,
		8,1,-1,
		8,2,-1,
		8,3,-1,
		9,-1,
		10,-1,
		10,0,-1,
		10,1,-1,
		10,2,-1,
		10,3,-1,
		11,-1,
		12,-1,
		13,-1,
		13,0,-1,
		13,1,-1,
		13,2,-1,
		13,3,-1,
		14,-1,
		15,-1,
		15,0,-1,
		15,1,-1,
		15,2,-1,
		15,3,-1,
		-2
	};
	__check_cell(&pos[0], root, itclone1);
	__check_cell(&pos[6], root, itclone2);
	__check_cell(&pos[15], root, itclone3);

	higcit_destroy(it);
	higcit_destroy(itclone1);
	higcit_destroy(itclone2);
	higcit_destroy(itclone3);

	hig_destroy(root);
}

ATF_TC(all_higtree_four_levels);
ATF_TC_HEAD(all_higtree_four_levels, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(all_higtree_four_levels, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	assert(root != NULL);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	DEBUG_PASS;
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid_two_levels(root, nc);
	hig_refine_uniform(hig_get_child(hig_get_child(root, 10), 0), nc);
	hig_refine_uniform(hig_get_child(hig_get_child(root, 13), 3), nc);
	hig_refine_uniform(hig_get_child(hig_get_child(hig_get_child(root, 13), 3), 1), nc);
	DEBUG_PASS;
	higcit_celliterator *it = higcit_create_all_higtree(root);

	int pos[] = {
		-1,
		0,-1,
		0,0,-1,
		0,1,-1,
		0,2,-1,
		0,3,-1,
		1,-1,
		2,-1,
		2,0,-1,
		2,1,-1,
		2,2,-1,
		2,3,-1,
		3,-1,
		4,-1,
		5,-1,
		5,0,-1,
		5,1,-1,
		5,2,-1,
		5,3,-1,
		6,-1,
		7,-1,
		7,0,-1,
		7,1,-1,
		7,2,-1,
		7,3,-1,
		8,-1,
		8,0,-1,
		8,1,-1,
		8,2,-1,
		8,3,-1,
		9,-1,
		10,-1,
		10,0,-1,
		10,0,0,-1,
		10,0,1,-1,
		10,0,2,-1,
		10,0,3,-1,
		10,1,-1,
		10,2,-1,
		10,3,-1,
		11,-1,
		12,-1,
		13,-1,
		13,0,-1,
		13,1,-1,
		13,2,-1,
		13,3,-1,
		13,3,0,-1,
		13,3,1,-1,
		13,3,1,0,-1,
		13,3,1,1,-1,
		13,3,1,2,-1,
		13,3,1,3,-1,
		13,3,2,-1,
		13,3,3,-1,
		14,-1,
		15,-1,
		15,0,-1,
		15,1,-1,
		15,2,-1,
		15,3,-1,
		-2
	};
	__check_cell(pos, root, it);

	higcit_destroy(it);

	it = higcit_create_all_higtree(hig_get_child(root,10));

	int pos2[] = {
		10,-1,
		10,0,-1,
		10,0,0,-1,
		10,0,1,-1,
		10,0,2,-1,
		10,0,3,-1,
		10,1,-1,
		10,2,-1,
		10,3,-1,
		-2
	};
	__check_cell(pos2, root, it);

	higcit_destroy(it);

	it = higcit_create_all_higtree(hig_get_child(hig_get_child(root,10), 0));

	int pos3[] = {
		10,0,-1,
		10,0,0,-1,
		10,0,1,-1,
		10,0,2,-1,
		10,0,3,-1,
		-2
	};
	__check_cell(pos3, root, it);

	higcit_destroy(it);

	it = higcit_create_all_higtree(hig_get_child(hig_get_child(hig_get_child(root,10), 0), 1));

	int pos4[] = {
		10,0,1,-1,
		-2
	};
	__check_cell(pos4, root, it);

	higcit_destroy(it);

	hig_destroy(root);
}

boolean __closer_to_point(hig_cell *c1, hig_cell *c2, void *d) {
	coordtype *point = (coordtype *)d;
	Point c1center;
	hig_get_center(c1, c1center);
	Point c2center;
	hig_get_center(c2, c2center);
	Point d1;
	POINT_SUB(d1, c1center, point);
	Point d2;
	POINT_SUB(d2, c2center, point);
	POINT_MULT(d1, d1, d1);
	POINT_MULT(d2, d2, d2);
	real squareddist1;
	POINT_FOLD(squareddist1, d1, 0.0, +);
	real squareddist2;
	POINT_FOLD(squareddist2, d2, 0.0, +);
	return FLT_LT(squareddist1, squareddist2);
}

boolean __closer_to_line(hig_cell *c1, hig_cell *c2, void *d) {
	real *x = (real *)d;
	Point c1center;
	hig_get_center(c1, c1center);
	Point c2center;
	hig_get_center(c2, c2center);
	return FLT_LT(fabs(c1center[0] - *x), fabs(c2center[0] - *x));
}

ATF_TC(higcit_sorted);
ATF_TC_HEAD(higcit_sorted, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(higcit_sorted, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	assert(root != NULL);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid_two_levels(root, nc);
	Point origin;
	higcit_celliterator *it;
	higcit_celliterator *its;

	it = higcit_create_all_leaves(root);
	POINT_ASSIGN_REALS(origin, 0.0, 0.0);
	its = higcit_create_sorted(it, 10, __closer_to_point, origin);
	int pos[] = {
		5,3,-1,
		10,0,-1,
		6,-1,
		9,-1,
		5,1,-1,
		5,2,-1,
		10,1,-1,
		10,2,-1,
		5,0,-1,
		10,3,-1,
		-2,
	};
	__check_cell(pos, root, its);
	higcit_destroy(its);

	it = higcit_create_all_leaves(root);
	POINT_ASSIGN_REALS(origin, -3.0/8.0, -5.0/8.0);
	its = higcit_create_sorted(it, 10, __closer_to_point, origin);
	int pos2[] = {
		4,-1,
		0,3,-1,
		5,0,-1,
		0,2,-1,
		5,2,-1,
		0,1,-1,
		5,1,-1,
		8,1,-1,
		1,-1,
		0,0,-1,
		-2,
	};
	__check_cell(pos2, root, its);
	higcit_destroy(its);

	it = higcit_create_all_leaves(root);
	real x = 0.3;
	its = higcit_create_sorted(it, 10, __closer_to_line, &x);
	int pos3[] = {
		9,-1,
		11,-1,
		8,2,-1,
		8,3,-1,
		10,2,-1,
		10,3,-1,
		8,0,-1,
		8,1,-1,
		10,0,-1,
		10,1,-1,
		-2,
	};
	__check_cell(pos3, root, its);
	higcit_destroy(its);

	hig_destroy(root);
}

ATF_TC(higcit_sorted_clone);
ATF_TC_HEAD(higcit_sorted_clone, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(higcit_sorted_clone, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	assert(root != NULL);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid_two_levels(root, nc);
	Point origin;
	higcit_celliterator *it;
	higcit_celliterator *its;

	it = higcit_create_all_leaves(root);
	POINT_ASSIGN_REALS(origin, 0.0, 0.0);
	its = higcit_create_sorted(it, 10, __closer_to_point, origin);
	higcit_celliterator *itclone1 = higcit_clone(its);
	for(int i = 0; i < 3; i++) {
		higcit_nextcell(its);
	}
	higcit_celliterator *itclone2 = higcit_clone(its);
	for(int i = 0; i < 3; i++) {
		higcit_nextcell(its);
	}
	higcit_celliterator *itclone3 = higcit_clone(its);
	int pos[] = {
		5,3,-1,
		10,0,-1,
		6,-1,
		9,-1,
		5,1,-1,
		5,2,-1,
		10,1,-1,
		10,2,-1,
		5,0,-1,
		10,3,-1,
		-2,
	};
	__check_cell(&pos[0], root, itclone1);
	__check_cell(&pos[8], root, itclone2);
	__check_cell(&pos[16], root, itclone3);
	higcit_destroy(its);
	higcit_destroy(itclone1);
	higcit_destroy(itclone2);
	higcit_destroy(itclone3);

	hig_destroy(root);
}

ATF_TC(higcit_concat);
ATF_TC_HEAD(higcit_concat, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(higcit_concat, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	assert(root != NULL);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid_two_levels(root, nc);
	Point origin;
	higcit_celliterator *it;
	higcit_celliterator *its[3];
	POINT_ASSIGN_REALS(lp, -0.8, -0.8);
	POINT_ASSIGN_REALS(hp, -0.1, -0.1);
	its[0] = higcit_create_bounding_box(root, lp, hp);

	POINT_ASSIGN_REALS(lp, -0.8, 0.4);
	POINT_ASSIGN_REALS(hp, 0.6, 0.6);
	its[1] = higcit_create_bounding_box(root, lp, hp);

	POINT_ASSIGN_REALS(lp, 0.9, -0.9);
	POINT_ASSIGN_REALS(hp, 0.95, -0.85);
	its[2] = higcit_create_bounding_box(root, lp, hp);

	it = higcit_create_concat(its, 3);
	int pos[] = {
		0,0,-1,
		0,1,-1,
		0,2,-1,
		0,3,-1,
		1,-1,
		4,-1,
		5,0,-1,
		5,1,-1,
		5,2,-1,
		5,3,-1,
		2,1,-1,
		2,3,-1,
		3,-1,
		6,-1,
		7,0,-1,
		7,2,-1,
		10,1,-1,
		10,3,-1,
		11,-1,
		14,-1,
		15,0,-1,
		12,-1,
		-2,
	};
	__check_cell(pos, root, it);
	higcit_destroy(it);
	hig_destroy(root);
}

ATF_TC(higcit_concat_clone);
ATF_TC_HEAD(higcit_concat_clone, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(higcit_concat_clone, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	assert(root != NULL);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid_two_levels(root, nc);
	Point origin;
	higcit_celliterator *it;
	higcit_celliterator *its[3];
	POINT_ASSIGN_REALS(lp, -0.8, -0.8);
	POINT_ASSIGN_REALS(hp, -0.1, -0.1);
	its[0] = higcit_create_bounding_box(root, lp, hp);

	POINT_ASSIGN_REALS(lp, -0.8, 0.4);
	POINT_ASSIGN_REALS(hp, 0.6, 0.6);
	its[1] = higcit_create_bounding_box(root, lp, hp);

	POINT_ASSIGN_REALS(lp, 0.9, -0.9);
	POINT_ASSIGN_REALS(hp, 0.95, -0.85);
	its[2] = higcit_create_bounding_box(root, lp, hp);

	it = higcit_create_concat(its, 3);
	higcit_celliterator *itclone1 = higcit_clone(it);
	for(int i = 0; i < 3; i++) {
		higcit_nextcell(it);
	}
	higcit_celliterator *itclone2 = higcit_clone(it);
	for(int i = 0; i < 3; i++) {
		higcit_nextcell(it);
	}
	higcit_celliterator *itclone3 = higcit_clone(it);

	int pos[] = {
		0,0,-1,
		0,1,-1,
		0,2,-1,
		0,3,-1,
		1,-1,
		4,-1,
		5,0,-1,
		5,1,-1,
		5,2,-1,
		5,3,-1,
		2,1,-1,
		2,3,-1,
		3,-1,
		6,-1,
		7,0,-1,
		7,2,-1,
		10,1,-1,
		10,3,-1,
		11,-1,
		14,-1,
		15,0,-1,
		12,-1,
		-2,
	};
	__check_cell(&pos[0], root, itclone1);
	__check_cell(&pos[9], root, itclone2);
	__check_cell(&pos[16], root, itclone3);
	higcit_destroy(it);
	higcit_destroy(itclone1);
	higcit_destroy(itclone2);
	higcit_destroy(itclone3);

	hig_destroy(root);
}

ATF_TP_ADD_TCS(tp)
{
    ATF_TP_ADD_TC(tp, all_cells_one_level);
    ATF_TP_ADD_TC(tp, all_cells_two_levels);
    ATF_TP_ADD_TC(tp, all_cells_two_levels_clone);
    ATF_TP_ADD_TC(tp, bounding_box_two_levels1);
    ATF_TP_ADD_TC(tp, bounding_box_two_levels2);
    ATF_TP_ADD_TC(tp, bounding_box_two_levels3);
    ATF_TP_ADD_TC(tp, bounding_box_two_levels_clone);
    ATF_TP_ADD_TC(tp, bounding_box_two_levels_outside);
    ATF_TP_ADD_TC(tp, bounding_box_two_levels_one_cell_with_border);
    ATF_TP_ADD_TC(tp, bounding_box_two_levels_one_point);
    ATF_TP_ADD_TC(tp, neighbours_two_levels);
    ATF_TP_ADD_TC(tp, neighbours_two_levels_clone);
    ATF_TP_ADD_TC(tp, filter_max_delta_two_levels);
    ATF_TP_ADD_TC(tp, filter_max_delta_two_levels_clone);
    ATF_TP_ADD_TC(tp, filter_close_to_origin_two_levels);
    ATF_TP_ADD_TC(tp, filter_close_to_origin_and_max_delta_two_levels);
    ATF_TP_ADD_TC(tp, filter_emtpy_two_levels);
    ATF_TP_ADD_TC(tp, filter_multi_grid_four_levels);
    ATF_TP_ADD_TC(tp, all_higtree_two_levels);
    ATF_TP_ADD_TC(tp, all_higtree_two_levels_clone);
    ATF_TP_ADD_TC(tp, all_higtree_four_levels);
    ATF_TP_ADD_TC(tp, higcit_sorted);
    ATF_TP_ADD_TC(tp, higcit_sorted_clone);
    ATF_TP_ADD_TC(tp, higcit_concat);
    ATF_TP_ADD_TC(tp, higcit_concat_clone);
}

