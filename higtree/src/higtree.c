#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "higtree.h"
#include "higtree-iterator.h"
#include "uniqueid.h"
#include "utils.h"

#include "Debug-c.h"

#ifdef HASFACETID
int hig_requires_facet_ids() {
	return 1;
}
#endif

#ifdef HASVERTEXID
int hig_requires_vertex_ids() {
	return 1;
}
#endif

int hig_toposition(const int numcells[DIM], const int pp[DIM]) {
	int res = 0;
	for (int i = 0; i < DIM; i++) {
		res *= numcells[i];
		res += pp[i];
	}
	return res;
}

void hig_tobase(int position, const int numcells[DIM], int pp[DIM]) {
	for (int i = DIM - 1; i >= 0; i--) {
		pp[i] = position % numcells[i];
		position /= numcells[i];
	}
}

void hig_translate_coord(Point lowpoint, Point delta, int pp[DIM], int numcells[DIM], Point chlp) {
	POINT_DIV(chlp, delta, numcells);
	POINT_MULT(chlp, chlp, pp);
	POINT_ADD(chlp, chlp, lowpoint);
}

int __hig_get_number_of_children(int numcells[DIM]) {
	int numChildren;
	POINT_FOLD(numChildren, numcells, 1, *);
	return numChildren;
}

int hig_get_number_of_children(hig_cell *c) {
	assert(c != NULL);

	return __hig_get_number_of_children(c->numcells);
}

void hig_get_cells_per_dim(hig_cell *c, int nc[DIM]) {
	POINT_ASSIGN(nc, c->numcells);
}

void hig_get_center(hig_cell *c, Point center) {
	assert(c != NULL);
	POINT_ADD(center, c->lowpoint, c->highpoint);

	POINT_DIV_SCALAR(center, center, 2.0);
}

void hig_get_delta(hig_cell *c, Point delta) {
	assert(c != NULL);

	POINT_SUB(delta, c->highpoint, c->lowpoint);
}

void hig_get_lowpoint(hig_cell *c, Point p) {
	assert(c != NULL);

	POINT_ASSIGN(p, c->lowpoint);
}

void hig_get_highpoint(hig_cell *c, Point p) {
	assert(c != NULL);

	POINT_ASSIGN(p, c->highpoint);
}

void hig_get_bounding_box(hig_cell *c, Rect *box) {
	hig_get_lowpoint(c, box->lo);
	hig_get_highpoint(c, box->hi);
}

void hig_get_element_coords(hig_cell *c, int elem[DIM], Point lg, Point hg) {
	Point lowpoint;
	Point highpoint;
	hig_get_lowpoint(c, lowpoint);
	hig_get_highpoint(c, highpoint);
	for(int i = 0; i < DIM; i++) {
		if (elem[i] == 0) {
			lg[i] = lowpoint[i];
			hg[i] = highpoint[i];
		} else if (elem[i] == -1) {
			lg[i] = hg[i] = lowpoint[i];
		} else if (elem[i] == 1) {
			lg[i] = hg[i] = highpoint[i];
		}
	}
}

uniqueid hig_get_id(hig_cell *c, int i) {
	assert(i < NUMIDS);
	assert(c != NULL);

	return c->ids[i];
}

uniqueid hig_get_cid(hig_cell *c) {
	assert(c != NULL);

	return hig_get_id(c, FIRSTCELLID);
}

uniqueid hig_get_vid(hig_cell *c, int i) {
	assert(c != NULL);

	return hig_get_id(c, FIRSTVERTEXID + i);
}

uniqueid hig_get_fid_with_dim_dir(hig_cell *c, int d, int i) {
	return c->ids[FIRSTFACETID+(d*2)+i];
}

void hig_set_fid_with_dim_dir(hig_cell *c, int d, int i, uniqueid id) {
	c->ids[FIRSTFACETID+(d*2)+i] = id;
}

hig_cell * hig_get_parent(hig_cell *c) {
	assert(c != NULL);

	return c->parent;
}

hig_cell *hig_get_child_in_grid(hig_cell *c, const int pp[DIM]) {
	assert(c != NULL);
	assert(c->children != NULL);

	int pos = hig_toposition(c->numcells, pp);
	assert(pos < hig_get_number_of_children(c));
	return c->children[pos];
}

hig_cell *hig_get_child(hig_cell *c, int pos) {
	assert(c != NULL);
	assert(c->children != NULL);
	assert(pos < hig_get_number_of_children(c));

	return c->children[pos];
}

hig_cell *hig_get_neighbour(hig_cell *c, int dir[DIM]) {
	assert(c != NULL);

	Point point;
	hig_get_center(c, point);
	for(int i = 0; i < DIM; i++) {
		if (dir[i] > 0) {
			point[i] = c->highpoint[i] + EPSDELTA;
		} else if (dir[i] < 0) {
			point[i] = c->lowpoint[i] - EPSDELTA;
		}
	}
	return hig_get_cell_with_point_from_another_cell(c, point);
}

int hig_get_narrowest_dim(hig_cell *c) {
	real minw = 0.0;
	int dim = -1;
	for(int i = 0; i < DIM; i++) {
		real w = c->highpoint[i] - c->lowpoint[i];
		if (dim == -1 || w < minw) {
			minw = w;
			dim = i;
		}
	}
	return dim;
}

real hig_sq_distance_to_point(hig_cell *c, const Point p)
{
	real accum = 0.0;
	for(unsigned dim = 0; dim < DIM; ++dim) {
		real delta = max(0,
			max(c->lowpoint[dim] - p[dim], p[dim] - c->highpoint[dim]));
		accum += delta * delta;
	}
	return accum;
}

void hig_get_relative_coord(hig_cell *c, Point relcoord, Point coord) {
	assert(c != NULL);

	POINT_ASSIGN(coord, c->lowpoint);
	Point delta;
	hig_get_delta(c, delta);
	POINT_MULT(delta, delta, relcoord);
	POINT_ADD(coord, coord, delta);
}

void hig_fill_empty(hig_cell *cell, Point lowpoint, Point highpoint)
{
	POINT_ASSIGN(cell->lowpoint, lowpoint);
	POINT_ASSIGN(cell->highpoint, highpoint);
	POINT_ASSIGN_SCALAR(cell->numcells, 0);

	cell->children = NULL;
	for(int i = 0; i < NUMIDS; i++) {
		cell->ids[i] = uid_getuniqueid();
	}
}

static hig_cell * __hig_create_cell(Point lowpoint, Point highpoint) {
	DECL_AND_ALLOC(hig_cell, c, 1);
	hig_fill_empty(c, lowpoint, highpoint);
	return c;
}

hig_cell * hig_create_empty_leaf(hig_cell *parent, int pp) {
	assert(parent->children[pp] == NULL);

	ALLOC_INFER(parent->children[pp], 1);
	hig_cell *c = parent->children[pp];
	c->parent = parent;
	c->posinparent = pp;

	return c;
}

hig_cell * hig_create_leaf(Point lowpoint, Point highpoint, hig_cell *cell, int pp) {
	assert(cell->children[pp] == NULL);

	hig_cell *c = __hig_create_cell(lowpoint, highpoint);
	c->posinparent = pp;
	c->parent = cell;
	for(int i = 0; i < NUMIDS; i++) {
		c->ids[i] = uid_getuniqueid();
	}
	cell->children[pp] = c;

	return c;
}

hig_cell * hig_create_root(Point lowpoint, Point highpoint) {
	hig_cell *c = __hig_create_cell(lowpoint, highpoint);
	c->posinparent = -1;
	c->parent = NULL;
	return c;
}

hig_cell * hig_refine_uniform(hig_cell * cell, int numcells[DIM]) {
	assert(cell != NULL);
	assert(cell->children == NULL);

	POINT_ASSIGN(cell->numcells, numcells);
	Point delta;
	hig_get_delta(cell, delta);
	POINT_DIV(delta, delta, cell->numcells);
	POINT_ASSERT_ALL_TRUE(numcells, > 0);

	int numChildren = hig_get_number_of_children(cell);
	ALLOC(hig_cell *, cell->children, numChildren);

	/* Zeroing allocated vector because hig_create_leaf()
	 * asserts if it is NULL... */
	memset(cell->children, 0, numChildren * sizeof *cell->children);

	for (int j = 0; j < numChildren; j++) {
		int pp[DIM];
		hig_tobase(j, numcells, pp);
		Point chlp, chhp;
	  	POINT_MULT(chlp, delta, pp);
	  	POINT_ADD(chlp, chlp, cell->lowpoint);
		POINT_ADD(chhp, chlp, delta);
		hig_create_leaf(chlp, chhp, cell, j);
	}
	return cell;
}

void hig_refine_empty(hig_cell * cell, int numcells[DIM]) {
	assert(cell != NULL);
	assert(cell->children == NULL);

	POINT_ASSIGN(cell->numcells, numcells);
	int numChildren = hig_get_number_of_children(cell);
	if(numChildren) {
		// Allocates a zeroed buffer
		ALLOC_INFER(cell->children, numChildren);
		memset(cell->children, 0, numChildren * sizeof * cell->children);
	}
}

//hig_cell * hig_refine(hig_cell * cell, int numcells[DIM], Point lowpoint[], Point highpoint[]) {
//	assert(cell != NULL);
//	assert(cell->children == NULL);
//
//	POINT_ASSIGN(cell->numcells, numcells);
//	Point delta;
//	hig_get_delta(cell, delta);
//	POINT_DIV(delta, delta, cell->numcells);
//	POINT_ASSERT_ALL_TRUE(numcells, > 0);
//	int numChildren = hig_get_number_of_children(cell);
//	cell->children = calloc(numChildren, sizeof(hig_cell *));
//	assert(cell->children != NULL);
//	for (int j = 0; j < numChildren; j++) {
//		int pp[DIM];
//		hig_tobase(j, numcells, pp);
//		Point chlp, chhp;
//		for (int i = 0; i < DIM; i++) {
//			chlp[i] = lowpoint[pp[i]][i];
//			chhp[i] = highpoint[pp[i]][i];
//		}
//		hig_create_leaf(chlp, chhp, cell, j);
//	}
//	return cell;
//}

void hig_destroy(hig_cell *cell) {
	assert(cell != NULL);

	if (cell->children != NULL) {
		int numChildren = hig_get_number_of_children(cell);
		for (int j = 0; j < numChildren; j++) {
			// There may be incomplete nodes, created with
			// hig_refine_empty(), where not all children were created.
			// So, before destroying each child, we check if is not NULL.
			if(cell->children[j])
				hig_destroy(cell->children[j]);
		}
		free(cell->children);
	}
	free(cell);
}

hig_cell * hig_merge_children(hig_cell *cell) {
	assert(cell != NULL);

	if (cell->children != NULL) {
		int numChildren = hig_get_number_of_children(cell);
		for (int j = 0; j < numChildren; j++) {
			hig_destroy(cell->children[j]);
		}
		free(cell->children);
		cell->children = NULL;
		POINT_ASSIGN_SCALAR(cell->numcells, 0);
	}
	return cell;
}

long hig_get_number_of_leaves(hig_cell *cell) {
	assert(cell != NULL);

	if (cell->children == NULL) {
		return 1;
	} else {
		int numleaves = 0;
		int numChildren = hig_get_number_of_children(cell);
		for (int j = 0; j < numChildren; j++) {
			numleaves += hig_get_number_of_leaves(cell->children[j]);
		}
		return numleaves;
	}
}

void hig_get_cell_coords_of_point(hig_cell *cell, const Point point, int pp[DIM]) {
	assert(cell != NULL);
	assert(cell->children != NULL);

	POINT_ASSIGN_SCALAR(pp, 0);
	for (int dim = 0; dim < DIM; dim++) {
		int first = 0;
		int last = cell->numcells[dim] - 1;
		while (last >= first) {
			pp[dim] = first;
			if (last == first) {
				break;
			} else {
				int pfirst = hig_toposition(cell->numcells, pp);
				pp[dim] = last;
				int plast = hig_toposition(cell->numcells, pp);
				coordtype delta = (cell->children[plast]->highpoint[dim] - cell->children[pfirst]->lowpoint[dim])/(last - first + 1);
				int pos = (int)((point[dim] - cell->children[pfirst]->lowpoint[dim])/delta) + first;
				if (pos > last) {
					pos = last;
				}
				if (pos < first) {
					pos = first;
				}
				pp[dim] = pos;
				int ppos = hig_toposition(cell->numcells, pp);
				if (FLT_LE(point[dim], cell->children[ppos]->lowpoint[dim])) {
					last = pos - 1;
				} else if (FLT_GT(point[dim], cell->children[ppos]->highpoint[dim])) {
					first = pos + 1;
				} else {
					break;
				}
			}
		}
	}
	for (int dim = 0; dim < DIM; dim++) {
		if (FLT_LT(point[dim], cell->lowpoint[dim])) {
			pp[dim] = -1;
		}
		if (FLT_GT(point[dim], cell->highpoint[dim])) {
			pp[dim] = cell->numcells[dim];
		}
	}
}

hig_cell * hig_get_cell_with_point(hig_cell * top, const Point point) {
	assert(top != NULL);

	hig_cell * cell = top;
	for (int i = 0; i < DIM; i++) {
		if (FLT_GT(cell->lowpoint[i], point[i]) || FLT_LT(cell->highpoint[i], point[i])) {
			return NULL;
		}
	}
	while (cell->children != NULL) {
		int pp[DIM];
		hig_get_cell_coords_of_point(cell, point, pp);
		int p = hig_toposition(cell->numcells, pp);
		cell = cell->children[p];
	}
	return cell;
}

hig_cell * hig_get_first_cell_recursive(hig_cell *cell) {
	assert(cell != NULL);

	hig_cell * p = cell;
	while (p->children != NULL) {
		p = p->children[0];
	}
	return p;
}

hig_cell * hig_get_cell_with_point_from_another_cell(hig_cell *cell, Point point) {
	if (cell == NULL) {
		return NULL;
	}
	for (int i = 0; i < DIM; i++) {
		if (FLT_LT(point[i], cell->lowpoint[i]) || FLT_GT(point[i], cell->highpoint[i])) {
			return hig_get_cell_with_point_from_another_cell(cell->parent, point);
		}
	}
	return hig_get_cell_with_point(cell, point);
}

real hig_distance_from_center(Point p, hig_cell *cell) {
	assert(cell != NULL);

	Point center;
	POINT_ADD(center, cell->highpoint, cell->lowpoint);
	POINT_DIV_SCALAR(center, center, 2.0);
	return co_distance(p, center);
}

bool hig_is_bounding_box_inside_cell(hig_cell *cell,
	const Point lowpoint, const Point highpoint)
{
	assert(cell != NULL);

	for (int i = 0; i < DIM; i++) {
		if (FLT_LT(lowpoint[i], cell->lowpoint[i]) || FLT_GT(highpoint[i], cell->highpoint[i])) {
			return 0;
		}
	}
	return 1;
}

bool hig_intersect_bounding_box(hig_cell *cell, const Point lowpoint, const Point highpoint)
{
	assert(cell != NULL);

	for (int i = 0; i < DIM; i++) {
		if ((FLT_GT(lowpoint[i], cell->highpoint[i]) || FLT_LT(highpoint[i], cell->lowpoint[i]))) {
			return 0;
		}
	}
	return 1;
}

bool __hig_equal(hig_cell *c1, hig_cell *c2) {
	for (int i = 0; i < DIM; i++) {
		if (FLT_NE(c1->lowpoint[i], c2->lowpoint[i]) || FLT_NE(c1->highpoint[i], c2->highpoint[i])) {
			return 0;
		}
		if (c1->numcells[i] != c2->numcells[i]) {
			return 0;
		}
	}
	if (c1->children == NULL || c2->children == NULL) {
		return c1->children == c2->children;
	}
	int numChildren = hig_get_number_of_children(c1);
	for (int j = 0; j < numChildren; j++) {
		if (!__hig_equal(c1->children[j], c2->children[j])) {
			return 0;
		}
	}
	return 1;
}

bool hig_equal(hig_cell *c1, hig_cell *c2) {
	return __hig_equal(c1, c2);
}

bool hig_check_conformable(hig_cell *cell) {
	if (cell == NULL) {
		return 1;
	}
	if (cell->children == NULL) {
		return 1;
	}
	int numChildren = hig_get_number_of_children(cell);
	for (int i = 0; i < DIM; i++) {
		if (FLT_NE(cell->lowpoint[i], cell->children[0]->lowpoint[i])) {
			return 0;
		}
		if (FLT_NE(cell->highpoint[i], cell->children[numChildren-1]->highpoint[i])) {
			return 0;
		}
	}
	int pp[DIM];
	for (int j = 0; j < numChildren; j++) {
		hig_tobase(j, cell->numcells, pp);
		for (int i = 0; i < DIM; i++) {
			if (pp[i] < cell->numcells[i] - 1) {
				int pp2[DIM];
				POINT_ASSIGN(pp2, pp);
				pp2[i]++;
				int j2 = hig_toposition(cell->numcells, pp2);
				if (FLT_NE(cell->children[j]->highpoint[i], cell->children[j2]->lowpoint[i])) {
					return 0;
				}
			}
		}
	}
	for (int j = 0; j < numChildren; j++) {
		if (!hig_check_conformable(cell->children[j])) {
			return 0;
		}
	}
	return 1;
}

long int hig_memory_used(hig_cell *parent) {
	long int used = sizeof(hig_cell);
	if (parent->children != NULL) {
		int numChildren = hig_get_number_of_children(parent);
		for (int j = 0; j < numChildren; j++) {
			used += hig_memory_used(parent->children[j]);
		}
		used += numChildren * sizeof(hig_cell *);
	}
	return used;
}

static void __hig_copy_regular_refinement_aux(hig_cell *from, hig_cell *to) {
	int numChildren = hig_get_number_of_children(from);
	if (numChildren > 0) {
		int nc[DIM];
		hig_get_cells_per_dim(from, nc);
		hig_refine_uniform(to, nc);
		for (int j = 0; j < numChildren; j++) {
			__hig_copy_regular_refinement_aux(hig_get_child(from, j), hig_get_child(to, j));
		}
	}
}

void hig_copy_regular_refinement(hig_cell *from, hig_cell *to) {
	__hig_copy_regular_refinement_aux(from, to);
}

hig_cell *hig_clone(hig_cell *from) {
	Point l, h;
	hig_get_lowpoint(from, l);
	hig_get_highpoint(from, h);
	hig_cell *to = hig_create_root(l, h);
	hig_copy_regular_refinement(from, to);
	return to;
}

void hig_adjust_facet_ids(hig_cell *root) {
	higcit_celliterator *it;
	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		Point cl, ch, ccenter;
		hig_get_lowpoint(c, cl);
		hig_get_highpoint(c, ch);
		hig_get_center(c, ccenter);
		for(int dim = 0; dim < DIM; dim++) {
			Point ofcenter;
			POINT_ASSIGN(ofcenter, ccenter);
			for(int dir = 0; dir < 2; dir++) {
				if(dir == 0) {
					ofcenter[dim] = cl[dim]-EPSDELTA;
				} else {
					ofcenter[dim] = ch[dim]+EPSDELTA;
				}
				bool isoppositefacetdifferent = false;
				bool isoppositefacetsmaller = false;
				hig_cell *ofc = hig_get_cell_with_point(root, ofcenter);
				if (ofc != NULL) {
					Point ofl, ofh;
					hig_get_lowpoint(ofc, ofl);
					hig_get_highpoint(ofc, ofh);
					for(int k = 0; k < DIM; k++) {
						if (k != dim) {
							if (FLT_NE(cl[k], ofl[k]) || FLT_NE(ch[k], ofh[k])) {
								isoppositefacetdifferent = true;
								if (FLT_GT((ch[k] - cl[k]), (ofh[k] - ofl[k]))) {
									isoppositefacetsmaller = true;
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
	}
	higcit_destroy(it);
}

int hig_get_facet_dim(hig_facet *facet) {
	return facet->dim;
}

int hig_get_facet_dir(hig_facet *facet) {
	return facet->dir;
}

hig_cell *hig_get_facet_cell(hig_facet *facet) {
	return facet->c;
}

uniqueid hig_get_fid(hig_facet *facet) {
	return hig_get_fid_with_dim_dir(facet->c, facet->dim, facet->dir);
}

void hig_get_facet_lowpoint(hig_facet *facet, Point l) {
	hig_get_lowpoint(facet->c, l);
	if (facet->dir == 1) {
		Point h;
		hig_get_highpoint(facet->c, h);
		l[facet->dim] = h[facet->dim];
	}
}

void hig_get_facet_highpoint(hig_facet *facet, Point h) {
	hig_get_highpoint(facet->c, h);
	if (facet->dir == 0) {
		Point l;
		hig_get_lowpoint(facet->c, l);
		h[facet->dim] = l[facet->dim];
	}
}

void hig_get_facet_center(hig_facet *facet, Point center) {
	hig_get_center(facet->c, center);
	if (facet->dir == 0) {
		Point l;
		hig_get_lowpoint(facet->c, l);
		center[facet->dim] = l[facet->dim];
	} else {
		Point h;
		hig_get_highpoint(facet->c, h);
		center[facet->dim] = h[facet->dim];
	}
}

void hig_get_facet_delta(hig_facet *facet, Point delta) {
	hig_get_delta(facet->c, delta);
}

int hig_intersect_intervals(double l1, double h1, double l2, double h2) {
	return ((l2 < h1) && (l1 < h2));
}

bool hig_facet_intersect_bounding_box(hig_facet *facet, const Point l, const Point h)
{
	Point fl, fh;
	hig_get_facet_lowpoint(facet, fl);
	hig_get_facet_highpoint(facet, fh);
	for (int dim = 0; dim < DIM; dim++) {
		if (!hig_intersect_intervals(l[dim], h[dim], fl[dim], fh[dim])) {
			return false;
		}
	}
	return true;
}

int hig_get_facet_with_point(hig_cell *root, int dim, const Point x, hig_facet *facet) {
	hig_cell *c = hig_get_cell_with_point(root, x);
	if (c == NULL) {
		return 0;
	}
	Point lowpoint, highpoint;
	hig_get_lowpoint(c, lowpoint);
	hig_get_highpoint(c, highpoint);
	int dir = -1;
	if (FLT_EQ(x[dim], lowpoint[dim])) {
		dir = 0;
	} else if (FLT_EQ(x[dim], highpoint[dim])) {
		dir = 1;
	}
	if (dir == -1) {
		return 0;
	}
	int id = hig_get_fid_with_dim_dir(c, dim, dir);
	if (id > 0) {
		facet->c = c;
		facet->dim = dim;
		facet->dir = dir;
		return 1;
	} else {
		Point x2;
		POINT_ASSIGN(x2, x);
		if (dir == 0) {
			x2[dim] -= EPSDELTA;
		} else {
			x2[dim] += EPSDELTA;
		}
		hig_cell *c2 = hig_get_cell_with_point_from_another_cell(c, x2);
		if (c2 == NULL) {
			return 0;
		}
		facet->c = c2;
		facet->dim = dim;
		facet->dir = 1-dir;
		return 1;
	}
	return 0;
}
