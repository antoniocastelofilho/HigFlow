#include <stdlib.h>
#include <string.h>
#include <atf-c.h>
#include "mapper.h"
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

ATF_TC(mapper_basic_operations);
ATF_TC_HEAD(mapper_basic_operations, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(mapper_basic_operations, tc)
{
	mp_mapper *m = mp_create();

	mp_assign(m, 0, 100);
	ATF_CHECK(mp_lookup(m, 0) == 100);
	ATF_CHECK(mp_lookup(m, 1) == MP_UNDEF);
	mp_assign(m, 1, 101);
	ATF_CHECK(mp_lookup(m, 1) == 101);
	mp_assign(m, 0, 101);
	ATF_CHECK(mp_lookup(m, 0) == 101);
	ATF_CHECK(mp_get_largest_key(m) == 1);

	mp_destroy(m);
}

ATF_TC(mapper_numbits);
ATF_TC_HEAD(mapper_numbits, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(mapper_numbits, tc)
{
	mp_mapper *m;
	m = mp_create();
	mp_set_numbits(m, MP_1_BIT);

	mp_assign(m, 10, 1);
	ATF_CHECK(mp_lookup(m, 10) == 1);
	ATF_CHECK(mp_lookup(m, 11) == 0);
	mp_assign(m, 11, 1);
	ATF_CHECK(mp_lookup(m, 11) == 1);
	mp_assign(m, 10, 0);
	ATF_CHECK(mp_lookup(m, 10) == 0);
	mp_assign(m, 103, 1);
	ATF_CHECK(mp_lookup(m, 103) == 1);
	ATF_CHECK(mp_lookup(m, 102) == 0);
	ATF_CHECK(mp_lookup(m, 104) == 0);
	mp_assign(m, 103, 0);
	ATF_CHECK(mp_lookup(m, 103) == 0);
	ATF_CHECK(mp_get_largest_key(m) == 103);
	mp_destroy(m);

	m = mp_create();
	mp_set_numbits(m, MP_2_BITS);
	mp_assign(m, 10, 3);
	ATF_CHECK(mp_lookup(m, 10) == 3);
	ATF_CHECK(mp_lookup(m, 11) == 0);
	mp_assign(m, 11, 2);
	ATF_CHECK(mp_lookup(m, 11) == 2);
	mp_assign(m, 10, 1);
	ATF_CHECK(mp_lookup(m, 10) == 1);
	mp_assign(m, 103, 3);
	ATF_CHECK(mp_lookup(m, 103) == 3);
	ATF_CHECK(mp_lookup(m, 102) == 0);
	ATF_CHECK(mp_lookup(m, 104) == 0);
	mp_assign(m, 103, 0);
	ATF_CHECK(mp_lookup(m, 103) == 0);

	mp_destroy(m);
	m = mp_create();
	mp_set_numbits(m, MP_8_BITS);
	mp_assign(m, 10, 3);
	ATF_CHECK(mp_lookup(m, 10) == 3);
	ATF_CHECK(mp_lookup(m, 11) == 0);
	mp_assign(m, 11, 5);
	ATF_CHECK(mp_lookup(m, 11) == 5);
	mp_assign(m, 10, 4);
	ATF_CHECK(mp_lookup(m, 10) == 4);
	mp_assign(m, 103, 7);
	ATF_CHECK(mp_lookup(m, 103) == 7);
	ATF_CHECK(mp_lookup(m, 102) == 0);
	ATF_CHECK(mp_lookup(m, 104) == 0);
	mp_assign(m, 103, 5);
	ATF_CHECK(mp_lookup(m, 103) == 5);

	mp_destroy(m);
}

ATF_TC(mapper_higcit_cell_iterator);
ATF_TC_HEAD(mapper_higcit_cell_iterator, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(mapper_higcit_cell_iterator, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid(root, nc);

	higcit_celliterator *it = higcit_create_all_leaves(root);

	Point p;
	mp_mapper *m = mp_create();
	mp_assign_from_celliterator(m, it, 0);
	higcit_destroy(it);

	POINT_ASSIGN_REALS(p, -0.4, -0.4);
	ATF_CHECK(mp_lookup(m, hig_get_id(hig_get_cell_with_point(root, p), 0)) == 11);

	POINT_ASSIGN_REALS(p, -0.4, 0.4);
	ATF_CHECK(mp_lookup(m, hig_get_id(hig_get_cell_with_point(root, p), 0)) == 15);

	POINT_ASSIGN_REALS(p, 0.4, -0.4);
	ATF_CHECK(mp_lookup(m, hig_get_id(hig_get_cell_with_point(root, p), 0)) == 24);

	POINT_ASSIGN_REALS(p, 0.4, 0.4);
	ATF_CHECK(mp_lookup(m, hig_get_id(hig_get_cell_with_point(root, p), 0)) == 28);

	hig_destroy(root);
	mp_destroy(m);
}

ATF_TC(mapper_many_higcit_celliterators);
ATF_TC_HEAD(mapper_many_higcit_celliterators, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(mapper_many_higcit_celliterators, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -4.0);
	POINT_ASSIGN_SCALAR(hp, 4.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid(root, nc);

	higcit_celliterator *it;
	mp_value_t id = 0;

	Point p;
	mp_mapper *m;
	m = mp_create();
	POINT_ASSIGN_REALS(lp, -3.5, -3.5);
	POINT_ASSIGN_REALS(hp, -1.5, -1.5);
	it = higcit_create_bounding_box(root, lp, hp);
	id = mp_assign_from_celliterator(m, it, id);
	higcit_destroy(it);

	POINT_ASSIGN_REALS(lp, 0.5, 0.5);
	POINT_ASSIGN_REALS(hp, 3.5, 3.5);
	it = higcit_create_bounding_box(root, lp, hp);
	id = mp_assign_from_celliterator(m, it, id);
	higcit_destroy(it);

	POINT_ASSIGN_REALS(lp, -2.5, -2.5);
	POINT_ASSIGN_REALS(hp, 2.5, 2.5);
	it = higcit_create_bounding_box(root, lp, hp);
	id = mp_assign_from_celliterator(m, it, id);
	higcit_destroy(it);

	it = higcit_create_all_leaves(root);
	id = mp_assign_from_celliterator(m, it, id);
	higcit_destroy(it);

	int ids[] = {
		0, 1, 2, 3,
		4,
		32, 33, 17, 18,
		19,

		5, 
		6, 20, 21, 22,
		23,
		24, 34, 25, 35,
		
		36, 26, 37, 27,
		28,
		7, 8, 9, 10,
		11,
		
		29,
		30, 31, 38, 39,
		12,
		13, 14, 15, 16,
	};
	it = higcit_create_all_leaves(root);
	for(int i = 0; !higcit_isfinished(it); higcit_nextcell(it), i++) {
		hig_cell *c = higcit_getcell(it);
		uniqueid cid = hig_get_id(c, 0);
		ATF_CHECK(mp_lookup(m, cid) == ids[i]);
	}
	higcit_destroy(it);
	mp_destroy(m);

	m = mp_create();
	id = 0;
	it = higcit_create_neighbours(hig_get_child(hig_get_child(root, 0), 2));
	id = mp_assign_from_celliterator(m, it, id);
	higcit_destroy(it);
	it = higcit_create_neighbours(hig_get_child(root, 12));
	id = mp_assign_from_celliterator(m, it, id);
	higcit_destroy(it);
	it = higcit_create_neighbours(hig_get_child(root, 11));
	id = mp_assign_from_celliterator(m, it, id);
	higcit_destroy(it);
	it = higcit_create_neighbours(hig_get_child(hig_get_child(root, 5), 1));
	id = mp_assign_from_celliterator(m, it, id);
	higcit_destroy(it);
	it = higcit_create_neighbours(hig_get_child(hig_get_child(root, 7), 1));
	id = mp_assign_from_celliterator(m, it, id);
	higcit_destroy(it);
	it = higcit_create_all_leaves(root);
	id = mp_assign_from_celliterator(m, it, id);
	higcit_destroy(it);
	int ids2[] = {
		0, 1, 24, 2,
		17,
		25, 26, 18, 27,
		22,

		3, 
		19, 28, 20, 21,
		9,
		23, 29, 10, 11,
		
		30, 31, 4, 5,
		6,
		32, 12, 33, 13,
		34,
		
		35,
		7, 36, 8, 37,
		14,
		15, 16, 38, 39,
	};
	it = higcit_create_all_leaves(root);
	for(int i = 0; !higcit_isfinished(it); higcit_nextcell(it), i++) {
		hig_cell *c = higcit_getcell(it);
		uniqueid cid = hig_get_id(c, 0);
		DEBUG_INSPECT(mp_lookup(m, cid), %d);
		ATF_CHECK(mp_lookup(m, cid) == ids2[i]);
	}
	higcit_destroy(it);
	mp_destroy(m);
	hig_destroy(root);
}

ATF_TC(mapper_basic_operations_pointer);
ATF_TC_HEAD(mapper_basic_operations_pointer, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(mapper_basic_operations_pointer, tc)
{
	mpp_mapper *m = mpp_create();

	void *x1 = malloc(1);
	void *x2 = malloc(1);
	mpp_assign(m, 10, x1);
	ATF_CHECK(mpp_lookup(m, 10) == x1);
	ATF_CHECK(mpp_lookup(m, 21) == NULL);
	ATF_CHECK(mpp_lookup(m, 25) == NULL);
	mpp_assign(m, 21, x2);
	ATF_CHECK(mpp_lookup(m, 21) == x2);
	mpp_assign(m, 10, x2);
	ATF_CHECK(mpp_lookup(m, 10) == x2);

	free(x1);
	free(x2);
	mpp_destroy(m);
}


ATF_TP_ADD_TCS(tp)
{
	ATF_TP_ADD_TC(tp, mapper_basic_operations);
	ATF_TP_ADD_TC(tp, mapper_basic_operations_pointer);
	ATF_TP_ADD_TC(tp, mapper_numbits);
	ATF_TP_ADD_TC(tp, mapper_higcit_cell_iterator);
	ATF_TP_ADD_TC(tp, mapper_many_higcit_celliterators);
}

