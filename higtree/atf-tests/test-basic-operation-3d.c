#include <string.h>
#include <atf-c.h>

#include "higtree.h"
#include "higtree-iterator.h"
#include "utils.h"
#define DEBUG
#include "Debug-c.h"

hig_cell *__build_mesh() {
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, 0.0);
	POINT_ASSIGN_SCALAR(hp, 27.0);

	hig_cell *root = hig_create_root(lp, hp);
	
	int numcells[DIM];
	POINT_ASSIGN_SCALAR(numcells, 3);
	DEBUG_INSPECT(numcells[0], %d);
	DEBUG_INSPECT(numcells[1], %d);
	DEBUG_INSPECT(numcells[2], %d);
	hig_refine_uniform(root, numcells);

	hig_cell *c13 = hig_get_child(root, 13);
	hig_refine_uniform(c13, numcells);

	hig_cell *c13_0 = hig_get_child(c13, 0);
	hig_refine_uniform(c13_0, numcells);
	return root;
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

ATF_TC(basic_operations_3d);
ATF_TC_HEAD(basic_operations_3d, tc)
{
	atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(basic_operations_3d, tc)
{
	hig_cell *root = __build_mesh();
	ATF_CHECK(hig_get_number_of_leaves(root) == 79);
	Point p;
	hig_cell *c;

	POINT_ASSIGN_REALS(p, 4.5, 4.5, 4.5);
	c = hig_get_cell_with_point(root, p);
	ATF_CHECK(hig_get_child(root, 0) == c);

	POINT_ASSIGN_REALS(p, 13.5, 4.5, 4.5);
	c = hig_get_cell_with_point(root, p);
	ATF_CHECK(hig_get_child(root, 9) == c);

	POINT_ASSIGN_REALS(p, 22.5, 4.5, 4.5);
	c = hig_get_cell_with_point(root, p);
	ATF_CHECK(hig_get_child(root, 18) == c);

	POINT_ASSIGN_REALS(p, 4.5, 12.5, 4.5);
	c = hig_get_cell_with_point(root, p);
	ATF_CHECK(hig_get_child(root, 3) == c);

	POINT_ASSIGN_REALS(p, 4.5, 22.5, 4.5);
	c = hig_get_cell_with_point(root, p);
	ATF_CHECK(hig_get_child(root, 6) == c);

	POINT_ASSIGN_REALS(p, 4.5, 4.5, 13.5);
	c = hig_get_cell_with_point(root, p);
	ATF_CHECK(hig_get_child(root, 1) == c);

	POINT_ASSIGN_REALS(p, 4.5, 4.5, 22.5);
	c = hig_get_cell_with_point(root, p);
	ATF_CHECK(hig_get_child(root, 2) == c);

	POINT_ASSIGN_REALS(p, 4.5, 13.5, 13.5);
	c = hig_get_cell_with_point(root, p);
	ATF_CHECK(hig_get_child(root, 4) == c);

	POINT_ASSIGN_REALS(p, 4.5, 13.5, 22.5);
	c = hig_get_cell_with_point(root, p);
	ATF_CHECK(hig_get_child(root, 5) == c);

	POINT_ASSIGN_REALS(p, 13.5, 13.5, 13.5);
	c = hig_get_cell_with_point(root, p);
	ATF_CHECK(hig_get_child(hig_get_child(root, 13), 13) == c);

	POINT_ASSIGN_REALS(p, 10.5, 10.5, 10.5);
	c = hig_get_cell_with_point(root, p);
	ATF_CHECK(hig_get_child(hig_get_child(hig_get_child(root, 13), 0), 13) == c);

	POINT_ASSIGN_REALS(p, 10.5, 9.5, 10.5);
	c = hig_get_cell_with_point(root, p);
	ATF_CHECK(hig_get_child(hig_get_child(hig_get_child(root, 13), 0), 10) == c);

	hig_destroy(root);
}

ATF_TC(all_cells_3d);
ATF_TC_HEAD(all_cells_3d, tc)
{
	atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(all_cells_3d, tc)
{
	hig_cell *root = __build_mesh();
	int indices[] = {
		0, -1,
		1, -1,
		2, -1,
		3, -1,
		4, -1,
		5, -1,
		6, -1,
		7, -1,
		8, -1,
		9, -1,
		10, -1,
		11, -1,
		12, -1,
		13, 0, 0, -1,
		13, 0, 1, -1,
		13, 0, 2, -1,
		13, 0, 3, -1,
		13, 0, 4, -1,
		13, 0, 5, -1,
		13, 0, 6, -1,
		13, 0, 7, -1,
		13, 0, 8, -1,
		13, 0, 9, -1,
		13, 0, 10, -1,
		13, 0, 11, -1,
		13, 0, 12, -1,
		13, 0, 13, -1,
		13, 0, 14, -1,
		13, 0, 15, -1,
		13, 0, 16, -1,
		13, 0, 17, -1,
		13, 0, 18, -1,
		13, 0, 19, -1,
		13, 0, 20, -1,
		13, 0, 21, -1,
		13, 0, 22, -1,
		13, 0, 23, -1,
		13, 0, 24, -1,
		13, 0, 25, -1,
		13, 0, 26, -1,
		13, 1, -1,
		13, 2, -1,
		13, 3, -1,
		13, 4, -1,
		13, 5, -1,
		13, 6, -1,
		13, 7, -1,
		13, 8, -1,
		13, 9, -1,
		13, 10, -1,
		13, 11, -1,
		13, 12, -1,
		13, 13, -1,
		13, 14, -1,
		13, 15, -1,
		13, 16, -1,
		13, 17, -1,
		13, 18, -1,
		13, 19, -1,
		13, 20, -1,
		13, 21, -1,
		13, 22, -1,
		13, 23, -1,
		13, 24, -1,
		13, 25, -1,
		13, 26, -1,
		14, -1,
		15, -1,
		16, -1,
		17, -1,
		18, -1,
		19, -1,
		20, -1,
		21, -1,
		22, -1,
		23, -1,
		24, -1,
		25, -1,
		26, -1,
		-2
	};
	higcit_celliterator *it = higcit_create_all_leaves(root);
	__check_cell(indices, root, it);
	higcit_destroy(it);
	hig_destroy(root);
}

ATF_TC(bounding_box_3d);
ATF_TC_HEAD(bounding_box_3d, tc)
{
	atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(bounding_box_3d, tc)
{
	hig_cell *root = __build_mesh();
	hig_cell *c;
	higcit_celliterator *it;

	int indices1[] = {
		13, 26, -1,
		14, -1,
		16, -1,
		17, -1,
		22, -1,
		23, -1,
		25, -1,
		-2
	};
	c = hig_get_child(root, 26);
	it = higcit_create_neighbours(c);
	__check_cell(indices1, root, it);
	higcit_destroy(it);
	
	int indices2[] = {
		0, -1,
		1, -1,
		3, -1,
		4, -1,
		6, -1,
		7, -1,
		9, -1,
		10, -1,
		13, 0, 0, -1,
		13, 0, 3, -1,
		13, 0, 6, -1,
		13, 0, 9, -1,
		13, 0, 12, -1,
		13, 0, 15, -1,
		13, 0, 18, -1,
		13, 0, 21, -1,
		13, 0, 24, -1,
		13, 3, -1,
		13, 6, -1,
		13, 9, -1,
		13, 12, -1,
		13, 15, -1,
		13, 18, -1,
		13, 21, -1,
		13, 24, -1,
		15, -1,
		16, -1,
		18, -1,
		19, -1,
		21, -1,
		22, -1,
		24, -1,
		25, -1,
		-2
	};
	c = hig_get_child(root, 12);
	it = higcit_create_neighbours(c);
	__check_cell(indices2, root, it);
	higcit_destroy(it);

	int indices3[] = {
		3, -1,
		4, -1,
		12, -1,
		13, 0, 3, -1,
		13, 0, 4, -1,
		13, 0, 7, -1,
		13, 0, 12, -1,
		13, 0, 13, -1,
		13, 0, 15, -1,
		13, 0, 16, -1,
		13, 3, -1,
		-2
	};
	c = hig_get_child(hig_get_child(hig_get_child(root, 13), 0), 6);
	it = higcit_create_neighbours(c);
	__check_cell(indices3, root, it);
	higcit_destroy(it);
	hig_destroy(root);
}

ATF_TC(neighbours_3d);
ATF_TC_HEAD(neighbours_3d, tc)
{
	atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(neighbours_3d, tc)
{
	hig_cell *root = __build_mesh();
	hig_cell *c;
	higcit_celliterator *it;

	int indices1[] = {
		0, -1,
		1, -1,
		3, -1,
		4, -1,
		9, -1,
		10, -1,
		12, -1,
		13, 0, 0, -1,
		13, 0, 1, -1,
		13, 0, 2, -1,
		13, 0, 3, -1,
		13, 0, 4, -1,
		13, 0, 5, -1,
		13, 0, 6, -1,
		13, 0, 7, -1,
		13, 0, 8, -1,
		13, 0, 9, -1,
		13, 0, 10, -1,
		13, 0, 11, -1,
		13, 0, 12, -1,
		13, 0, 13, -1,
		13, 0, 14, -1,
		13, 0, 15, -1,
		13, 0, 16, -1,
		13, 0, 17, -1,
		13, 1, -1,
		13, 3, -1,
		13, 4, -1,
		-2
	};
	Point lp, hp;
	POINT_ASSIGN_REALS(lp, 4.5, 4.5, 4.5);
	POINT_ASSIGN_REALS(hp, 10.5, 13.5, 13.5);
	it = higcit_create_bounding_box(root, lp, hp);
	__check_cell(indices1, root, it);
	higcit_destroy(it);
	
	hig_destroy(root);
}

ATF_TP_ADD_TCS(tp)
{
	ATF_TP_ADD_TC(tp, basic_operations_3d);
	ATF_TP_ADD_TC(tp, all_cells_3d);
	ATF_TP_ADD_TC(tp, neighbours_3d);
	ATF_TP_ADD_TC(tp, bounding_box_3d);
}
