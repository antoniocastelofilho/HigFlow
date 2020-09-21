#include <string.h>
#include <atf-c.h>
#include "higtree.h"
#include "utils.h"
#define DEBUG
#include "Debug-c.h"

ATF_TC(basic_operations);
ATF_TC_HEAD(basic_operations, tc)
{
    atf_tc_set_md_var(tc, "descr", "This test case ensures that hig_create_root works");
}
ATF_TC_BODY(basic_operations, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	Point center;
	hig_get_center(root, center);
	int i;
	for(i = 0; i < DIM; i++) {
		ATF_CHECK(FLT_EQ(center[i], 0.0));
	}
	ATF_CHECK(root->children == NULL);
	ATF_CHECK(hig_get_id(root, 0) > 0);

	ATF_CHECK(FLT_EQ(hig_distance_from_center(lp, root), sqrt(2.0)));

	ATF_CHECK(hig_memory_used(root) == sizeof(hig_cell));
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 2);
	hig_refine_uniform(root, nc);
	ATF_CHECK(hig_memory_used(root) == 5*sizeof(hig_cell) + 4*sizeof(hig_cell*));
	
	Point l, h;
	POINT_ASSIGN_SCALAR(l, -0.5);
	POINT_ASSIGN_SCALAR(h, 0.5);
	ATF_CHECK(hig_is_bounding_box_inside_cell(root, l, h));

	hig_destroy(root);
}

ATF_TC(translate_coord);
ATF_TC_HEAD(translate_coord, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(translate_coord, tc)
{
	Point lp, delta, coord;
	int nc[DIM];
	int pp[DIM];

	lp[0] = 1.0;
	lp[1] = 0.5;
	nc[0] = 3;
	nc[1] = 4;
	delta[0] = 1.5;
	delta[1] = 1.0;

	pp[0] = 1;
	pp[1] = 1;
	hig_translate_coord(lp, delta, pp, nc, coord);
	ATF_CHECK(FLT_EQ(coord[0], 1.5));
	ATF_CHECK(FLT_EQ(coord[1], 0.75));

	pp[0] = 2;
	pp[1] = 0;
	hig_translate_coord(lp, delta, pp, nc, coord);
	ATF_CHECK(FLT_EQ(coord[0], 2.0));
	ATF_CHECK(FLT_EQ(coord[1], 0.5));
}

ATF_TC(hig_get_neighbour);
ATF_TC_HEAD(hig_get_neighbour, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(hig_get_neighbour, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);

	hig_cell *c, *n;
	int dir[DIM];

	c = hig_get_child(root, 5);
	dir[0] = 0;
	dir[1] = 0;
	n = hig_get_neighbour(c, dir);
	ATF_CHECK(n == hig_get_child(root, 5));

	dir[0] = 1;
	dir[1] = 0;
	n = hig_get_neighbour(c, dir);
	ATF_CHECK(n == hig_get_child(root, 9));

	dir[0] = -1;
	dir[1] = 0;
	n = hig_get_neighbour(c, dir);
	ATF_CHECK(n == hig_get_child(root, 1));

	dir[0] = -1;
	dir[1] = 1;
	n = hig_get_neighbour(c, dir);
	ATF_CHECK(n == hig_get_child(root, 2));

	dir[0] = 1;
	dir[1] = 1;
	n = hig_get_neighbour(c, dir);
	ATF_CHECK(n == hig_get_child(root, 10));

	dir[0] = -1;
	dir[1] = -1;
	n = hig_get_neighbour(c, dir);
	ATF_CHECK(n == hig_get_child(root, 0));

	c = hig_get_child(root, 1);
	dir[0] = 0;
	dir[1] = -1;
	n = hig_get_neighbour(c, dir);
	ATF_CHECK(n == hig_get_child(root, 0));

	c = hig_get_child(root, 1);
	dir[0] = -1;
	dir[1] = 1;
	n = hig_get_neighbour(c, dir);
	ATF_CHECK(n == NULL);

	hig_destroy(root);
}


ATF_TC(hig_check_conformable);
ATF_TC_HEAD(hig_check_conformable, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(hig_check_conformable, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	int nc[DIM];
	hig_cell *c1, *c2;

	hig_cell *root = hig_create_root(lp, hp);
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	ATF_CHECK(hig_check_conformable(root));
	hig_destroy(root);

	ATF_CHECK(hig_check_conformable(NULL));

	root = hig_create_root(lp, hp);
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	c1 = hig_get_child(root, 0);
	c1->lowpoint[0] += EPSDELTA;
	ATF_CHECK(!hig_check_conformable(root));
	hig_destroy(root);

	root = hig_create_root(lp, hp);
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	c1 = hig_get_child(root, 15);
	c1->highpoint[0] += EPSDELTA;
	ATF_CHECK(!hig_check_conformable(root));
	hig_destroy(root);

	root = hig_create_root(lp, hp);
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	c1 = hig_get_child(root, 15);
	hig_refine_uniform(c1, nc);
	c2 = hig_get_child(c1, 15);
	c2->highpoint[0] += EPSDELTA;
	ATF_CHECK(!hig_check_conformable(root));
	hig_destroy(root);

	root = hig_create_root(lp, hp);
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	c1 = hig_get_child(root, 0);
	c1->highpoint[0] += EPSDELTA;
	ATF_CHECK(!hig_check_conformable(root));
	hig_destroy(root);

}

ATF_TC(baseposition);
ATF_TC_HEAD(baseposition, tc)
{
    atf_tc_set_md_var(tc, "descr", "This test case ensures that hig_toposition and hig_tobase invert each other");
	if (DIM != 2)
		atf_tc_skip("DIM != 2");
}
ATF_TC_BODY(baseposition, tc)
{
	int nc[DIM];
	int pp[DIM];
	int pos;
	int mx = 3, my = 4;
	nc[0] = mx;
	nc[1] = my;
	pp[0] = 1;
	pp[1] = 1;
	pos = hig_toposition(nc, pp);
	ATF_CHECK(pos == 1 + my);
	pp[0] = 2;
	pp[1] = 3;
	pos = hig_toposition(nc, pp);
	ATF_CHECK(pos == 2*my + 3);
	int i, j;
	int count = 0;
	for(i = 0; i < mx; i++) {
		pp[0] = i;
		for(j = 0; j < my; j++) {
			pp[1] = j;
			pos = hig_toposition(nc, pp);
			ATF_CHECK(pos == count);
			count++;
		}
	}
	hig_tobase(9, nc, pp);
	ATF_CHECK(pp[0] == 2);
	ATF_CHECK(pp[1] == 1);
	for(i = 0; i < mx; i++) {
		pp[0] = i;
		for(j = 0; j < my; j++) {
			pp[1] = j;
			pos = hig_toposition(nc, pp);
			int pp2[DIM];
			hig_tobase(pos, nc, pp2);
			for(int i = 0; i < DIM; i++) {
				ATF_CHECK(pp[i] == pp2[i]);
			}
		}
	}
}

ATF_TC(number_of_children);
ATF_TC_HEAD(number_of_children, tc)
{
    atf_tc_set_md_var(tc, "descr", "This test case ensures that hig_get_number_of_children works");
}
ATF_TC_BODY(number_of_children, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	long int numleaves = hig_get_number_of_leaves(root);
	ATF_CHECK(numleaves == 1);
	ATF_CHECK(hig_get_number_of_children(root) == 0);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 5);
	hig_refine_uniform(root, nc);
	numleaves = hig_get_number_of_leaves(root);
	ATF_CHECK(numleaves == 25);
	ATF_CHECK(hig_get_number_of_children(root) == 25);
	Point p;
	POINT_ASSIGN_SCALAR(p, 0.5);
	hig_cell * c = hig_get_cell_with_point(root, p);
	ATF_CHECK(c != NULL);
	hig_refine_uniform(c, nc);
	numleaves = hig_get_number_of_leaves(root);
	ATF_CHECK(numleaves == 49);
	c = hig_get_cell_with_point(root, p);
	ATF_CHECK(c != NULL);
	hig_refine_uniform(c, nc);
	numleaves = hig_get_number_of_leaves(root);
	ATF_CHECK(numleaves == 73);
	hig_destroy(root);
}

ATF_TC(relative_coord);
ATF_TC_HEAD(relative_coord, tc)
{
    atf_tc_set_md_var(tc, "descr", "This test case ensures that relative_coord works");
}
ATF_TC_BODY(relative_coord, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	long int numleaves = hig_get_number_of_leaves(root);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 5);
	hig_refine_uniform(root, nc);
	Point p;
	POINT_ASSIGN_SCALAR(p, 0.5);
	hig_cell * c = hig_get_cell_with_point(root, p);

	Point relcoord;
	relcoord[0] = -0.5;
	relcoord[1] = -0.5;
	Point coords;
	hig_get_relative_coord(c, relcoord, coords);
	ATF_CHECK(FLT_EQ(coords[0], 0.0));
	ATF_CHECK(FLT_EQ(coords[1], 0.0));

	relcoord[0] = 1.25;
	relcoord[1] = -0.5;
	hig_get_relative_coord(c, relcoord, coords);
	ATF_CHECK(FLT_EQ(coords[0], 0.7));
	ATF_CHECK(FLT_EQ(coords[1], 0.0));

	relcoord[0] = 0.5;
	relcoord[1] = 0.5;
	hig_get_relative_coord(c, relcoord, coords);
	ATF_CHECK(FLT_EQ(coords[0], 0.4));
	ATF_CHECK(FLT_EQ(coords[1], 0.4));

	relcoord[0] = 0.0;
	relcoord[1] = 1.0;
	hig_get_relative_coord(c, relcoord, coords);
	ATF_CHECK(FLT_EQ(coords[0], 0.2));
	ATF_CHECK(FLT_EQ(coords[1], 0.6));

	hig_destroy(root);
}


ATF_TC(hig_equal);
ATF_TC_HEAD(hig_equal, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(hig_equal, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 5);
	hig_cell *c1, *c2;
	c1 = hig_create_root(lp, hp);
	hig_refine_uniform(c1, nc);
	c2 = hig_create_root(lp, hp);
	hig_refine_uniform(c2, nc);
	
	ATF_CHECK(hig_equal(c1, c2));
	hig_refine_uniform(hig_get_child(c2, 0), nc);
	ATF_CHECK(!hig_equal(c1, c2));
	hig_destroy(c2);

	POINT_ASSIGN_SCALAR(nc, 4);
	c2 = hig_create_root(lp, hp);
	hig_refine_uniform(c2, nc);
	ATF_CHECK(!hig_equal(c1, c2));
	hig_destroy(c2);

	POINT_ASSIGN_SCALAR(lp, -0.5);
	c2 = hig_create_root(lp, hp);
	hig_refine_uniform(c2, nc);
	ATF_CHECK(!hig_equal(c1, c2));
	hig_destroy(c2);

	hig_destroy(c1);
}

ATF_TC(get_child);
ATF_TC_HEAD(get_child, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(get_child, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 5);
	hig_refine_uniform(root, nc);
	Point p;
	int pp[DIM];
	int cnt = 0;
	for(int i = 0; i < 5; i++) {
		pp[0] = i;
		for (int j = 0; j < 5; j++) {
			pp[1] = j;
			ATF_CHECK(root->children[cnt] == hig_get_child_in_grid(root, pp));
			ATF_CHECK(root->children[cnt] == hig_get_child(root, cnt));
			cnt++;
		}
	}
	hig_destroy(root);
}

ATF_TC(get_cell_coords_of_point);
ATF_TC_HEAD(get_cell_coords_of_point, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(get_cell_coords_of_point, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 5);
	hig_refine_uniform(root, nc);
	Point p;
	int pp[DIM];

	p[0] = 0.0;
	p[1] = 0.0;
	hig_get_cell_coords_of_point(root, p, pp);
	ATF_CHECK(pp[0] == 2);
	ATF_CHECK(pp[1] == 2);

	p[0] = 1.0;
	p[1] = 1.0;
	hig_get_cell_coords_of_point(root, p, pp);
	ATF_CHECK(pp[0] == 4);
	ATF_CHECK(pp[1] == 4);

	p[0] = 1.1;
	p[1] = 1.1;
	hig_get_cell_coords_of_point(root, p, pp);
	ATF_CHECK(pp[0] == 5);
	ATF_CHECK(pp[1] == 5);

	p[0] = -1.1;
	p[1] = -1.1;
	hig_get_cell_coords_of_point(root, p, pp);
	ATF_CHECK(pp[0] == -1);
	ATF_CHECK(pp[1] == -1);

	p[0] = 0.2;
	p[1] = 0.2;
	hig_get_cell_coords_of_point(root, p, pp);
	ATF_CHECK(pp[0] == 2);
	ATF_CHECK(pp[1] == 2);

	p[0] = -0.2;
	p[1] = -0.2;
	hig_get_cell_coords_of_point(root, p, pp);
	ATF_CHECK(pp[0] == 1);
	ATF_CHECK(pp[1] == 1);

	hig_destroy(root);
}

ATF_TC(get_cell_with_point);
ATF_TC_HEAD(get_cell_with_point, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(get_cell_with_point, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 5);
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	for(int i = 0; i < 25; i+=2) {
		hig_refine_uniform(hig_get_child(root, i), nc);
	}
	Point p;
	hig_cell *c;

	p[0] = 0.05;
	p[1] = 0.05;
 	c = hig_get_cell_with_point(root, p);	
	ATF_CHECK(c == hig_get_child(hig_get_child(root, 12), 3));

	p[0] = 0.0;
	p[1] = 0.05;
 	c = hig_get_cell_with_point(root, p);	
	ATF_CHECK(c == hig_get_child(hig_get_child(root, 12), 1));

	p[0] = -0.4;
	p[1] = 0.2;
 	c = hig_get_cell_with_point(root, p);	
	ATF_CHECK(c == hig_get_child(root, 7));

	p[0] = 0.0;
	p[1] = 0.0;
 	c = hig_get_cell_with_point(root, p);	
	ATF_CHECK(c == hig_get_child(hig_get_child(root, 12), 0));

	p[0] = -0.5;
	p[1] = -0.5;
 	c = hig_get_cell_with_point(root, p);	
	ATF_CHECK(c == hig_get_child(hig_get_child(root, 6), 0));

	p[0] = -0.2;
	p[1] = -0.2;
 	c = hig_get_cell_with_point(root, p);	
	ATF_CHECK(c == hig_get_child(hig_get_child(root, 6), 3));

	p[0] = -1.0;
	p[1] = -1.0;
 	c = hig_get_cell_with_point(root, p);	
	ATF_CHECK(c == hig_get_child(hig_get_child(root, 0), 0));

	p[0] = -0.8;
	p[1] = -0.8;
 	c = hig_get_cell_with_point(root, p);	
	ATF_CHECK(c == hig_get_child(hig_get_child(root, 0), 0));

	p[0] = -1.01;
	p[1] = -1.01;
 	c = hig_get_cell_with_point(root, p);	
	ATF_CHECK(c == NULL);

	p[0] = 1.0;
	p[1] = 1.0;
 	c = hig_get_cell_with_point(root, p);	
	ATF_CHECK(c == hig_get_child(hig_get_child(root, 24), 3));

	hig_cell *top = hig_get_child(root, 12);

	p[0] = -0.5;
	p[1] = -0.5;
 	c = hig_get_cell_with_point(top, p);	
	ATF_CHECK(c == NULL);

	p[0] = -0.2;
	p[1] = -0.2;
 	c = hig_get_cell_with_point(top, p);	
	ATF_CHECK(c == hig_get_child(hig_get_child(root, 12), 0));

	hig_destroy(root);
}

ATF_TC(hig_merge_children);
ATF_TC_HEAD(hig_merge_children, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(hig_merge_children, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 5);
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	for(int i = 0; i < 25; i+=2) {
		hig_refine_uniform(hig_get_child(root, i), nc);
	}
	Point p;
	hig_cell *c;

	ATF_CHECK(hig_get_number_of_leaves(root) == 64);
	hig_merge_children(hig_get_child(root, 0));
	ATF_CHECK(hig_get_number_of_leaves(root) == 61);
	hig_merge_children(hig_get_child(root, 4));
	ATF_CHECK(hig_get_number_of_leaves(root) == 58);
	hig_merge_children(hig_get_child(root, 8));
	ATF_CHECK(hig_get_number_of_leaves(root) == 55);
	hig_merge_children(hig_get_child(root, 12));
	ATF_CHECK(hig_get_number_of_leaves(root) == 52);
	hig_merge_children(hig_get_child(root, 16));
	ATF_CHECK(hig_get_number_of_leaves(root) == 49);
	hig_merge_children(hig_get_child(root, 20));
	ATF_CHECK(hig_get_number_of_leaves(root) == 46);
	hig_merge_children(hig_get_child(root, 24));
	ATF_CHECK(hig_get_number_of_leaves(root) == 43);

	for(int i = 0; i < 25; i++) {
		ATF_CHECK(hig_get_number_of_children(hig_get_child(root, i)) == ((i%4==2)?4:0));
	}

	hig_destroy(root);
}

ATF_TP_ADD_TCS(tp)
{
    ATF_TP_ADD_TC(tp, basic_operations);
    ATF_TP_ADD_TC(tp, baseposition);
    ATF_TP_ADD_TC(tp, translate_coord);
    ATF_TP_ADD_TC(tp, number_of_children);
    ATF_TP_ADD_TC(tp, relative_coord);
    ATF_TP_ADD_TC(tp, get_child);
    ATF_TP_ADD_TC(tp, get_cell_coords_of_point);
    ATF_TP_ADD_TC(tp, get_cell_with_point);
    ATF_TP_ADD_TC(tp, hig_get_neighbour);
    ATF_TP_ADD_TC(tp, hig_check_conformable);
    ATF_TP_ADD_TC(tp, hig_equal);
    ATF_TP_ADD_TC(tp, hig_merge_children);
}
