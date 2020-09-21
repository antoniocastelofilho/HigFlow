
#include<stdlib.h>
#include<stdio.h>
#include<assert.h>
#include<stdbool.h>
#include "higtree.h"
#include "higtree-iterator.h"
#include "uniqueid.h"
#include "utils.h"

#include "Debug-c.h"

hig_cell * higcit_getcell(higcit_celliterator *it){
	return it->bi.c;
}

bool higcit_isfinished(higcit_celliterator *it){
	return it->bi.c == NULL;
}

hig_cell * higcit_nextcell(higcit_celliterator *it) {
	assert(it->bi.next != NULL);
	if (it->bi.c != NULL) {
		it->bi.c = it->bi.next(it);
	}
	return higcit_getcell(it);
}

void higcit_destroy(higcit_celliterator *it) {
	assert(it->bi.destroy != NULL);
	it->bi.destroy(it);
}

higcit_celliterator * higcit_clone(higcit_celliterator *it) {
	assert(it->bi.clone != NULL);
	return it->bi.clone(it);
}

long higcit_count(higcit_celliterator *it) {
	long cnt = 0;
	for(;!higcit_isfinished(it); higcit_nextcell(it)) {
		cnt++;
	}
	return cnt;
}

long higcit_count_without_advancing(higcit_celliterator *it) {
	higcit_celliterator *itc = higcit_clone(it);
	long cnt = higcit_count(itc);
	higcit_destroy(itc);
	return cnt;
}

static hig_cell * __higcit_all_leaves_inc(higcit_celliterator *it) {
	higcit_all_leaves_extra *itex = (higcit_all_leaves_extra *) it->ex;
	hig_cell *p = it->bi.c;
	it->bi.c = NULL;
	while (p != itex->root && p->parent != NULL) {
		int pos = p->posinparent + 1;
		if (pos < hig_get_number_of_children(p->parent)) {
			it->bi.c = hig_get_first_cell_recursive(p->parent->children[pos]);
			break;
		} else {
			p = p->parent;
		}
	}
	return it->bi.c;
}

static void __higcit_simple_destroy(higcit_celliterator *it) {
	free(it->ex);
	free(it);
}

static struct higcit_celliterator * __higcit_all_leaves_clone(higcit_celliterator *oldit) {
	higcit_all_leaves_extra *olditex = (higcit_all_leaves_extra *) oldit->ex;
	DECL_AND_ALLOC(higcit_celliterator, it, 1);
	DECL_AND_ALLOC(higcit_all_leaves_extra, itex, 1);
	itex->root = olditex->root;
	it->ex = itex;
	it->bi.next = __higcit_all_leaves_inc;
	it->bi.destroy = __higcit_simple_destroy;
	it->bi.clone = __higcit_all_leaves_clone;
	it->bi.c = oldit->bi.c;
	return it;
}

higcit_celliterator * higcit_create_all_leaves(hig_cell *cell) {
	DECL_AND_ALLOC(higcit_celliterator, it, 1);
	DECL_AND_ALLOC(higcit_all_leaves_extra, itex, 1);
	itex->root = cell;
	it->ex = itex;
	it->bi.next = __higcit_all_leaves_inc;
	it->bi.destroy = __higcit_simple_destroy;
	it->bi.clone = __higcit_all_leaves_clone;
	it->bi.c = hig_get_first_cell_recursive(cell);
	return it;
}


static hig_cell *__higcit_neighbours_inc(higcit_celliterator *it) {
	higcit_neighbours_extra *itex = (higcit_neighbours_extra *) it->ex;
	it->bi.c = higcit_nextcell(itex->it);
	if (it->bi.c == itex->cell){
		it->bi.c = higcit_nextcell(itex->it);
	}
	return it->bi.c;
}

static void __higcit_neighbours_destroy(higcit_celliterator *it) {
	higcit_neighbours_extra *itex = (higcit_neighbours_extra *) it->ex;
	higcit_destroy(itex->it);
	free(it->ex);
	free(it);
}

static higcit_celliterator *__higcit_neighbours_clone(higcit_celliterator *oldit) {
	higcit_neighbours_extra *olditex = (higcit_neighbours_extra *) oldit->ex;
	DECL_AND_ALLOC(higcit_celliterator, it, 1);
	DECL_AND_ALLOC(higcit_neighbours_extra, itex, 1);
	it->ex = itex;
	itex->cell = olditex->cell;
	itex->it = higcit_clone(olditex->it);
	it->bi.next = __higcit_neighbours_inc;
	it->bi.destroy = __higcit_neighbours_destroy;
	it->bi.clone = __higcit_neighbours_clone;
	it->bi.c = oldit->bi.c;
	return it;
}

higcit_celliterator * higcit_create_neighbours(hig_cell *cell) {
	DECL_AND_ALLOC(higcit_celliterator, it, 1);
	DECL_AND_ALLOC(higcit_neighbours_extra, itex, 1);
	it->ex = itex;
	it->bi.next = __higcit_neighbours_inc;
	it->bi.destroy = __higcit_neighbours_destroy;
	it->bi.clone = __higcit_neighbours_clone;
	itex->cell = cell;
	Point lowpoint;
	Point highpoint;
	POINT_SUB_SCALAR(lowpoint, cell->lowpoint, EPSDELTA);
	POINT_ADD_SCALAR(highpoint, cell->highpoint, EPSDELTA);
	itex->it = higcit_create_bounding_box(cell, lowpoint, highpoint);
	it->bi.c = higcit_getcell(itex->it);
	if (it->bi.c == itex->cell){
		it->bi.c = higcit_nextcell(itex->it);
	}
	return it;
}

static hig_cell *__higcit_get_first_of_bounding_box(hig_cell *cell, higcit_bounding_box_extra *itex) {
	hig_cell *p = cell;
	while (p->children != NULL) {
		itex->curlevel++;
		int level = itex->curlevel;
		itex->cells[level] = p;
		hig_get_cell_coords_of_point(p, itex->lowpoint, itex->ppl[level]);
		hig_get_cell_coords_of_point(p, itex->highpoint, itex->pph[level]);
		for(int i = 0; i < DIM; i++) {
			if (itex->ppl[level][i] < 0) {
				itex->ppl[level][i] = 0;
			}
			if (itex->pph[level][i] >= p->numcells[i]) {
				itex->pph[level][i] = p->numcells[i] - 1;
			}
		}
		POINT_ASSIGN(itex->pps[level], itex->ppl[level]);
		int j = hig_toposition(p->numcells, itex->pps[level]);
		p = p->children[j];
	}
	return p;
}

static hig_cell *__higcit_bounding_box_inc(higcit_celliterator *it) {
	higcit_bounding_box_extra *itex = (higcit_bounding_box_extra *) it->ex;
	it->bi.c = NULL;
	while (itex->curlevel >= 0) {
		int level = itex->curlevel;
		int i = DIM - 1;
		while (i >= 0 && itex->pps[level][i] >= itex->pph[level][i]) {
	  		itex->pps[level][i] = itex->ppl[level][i];
			i--;
		}
		if (i >= 0) {
			itex->pps[level][i]++;
			hig_cell *p = itex->cells[level];
			int j = hig_toposition(p->numcells, itex->pps[level]);
			p = p->children[j];
			it->bi.c = __higcit_get_first_of_bounding_box(p, itex);
			break;
		} else {
			itex->curlevel--;
		}
	}
	return it->bi.c;
}

static higcit_celliterator *__higcit_bounding_box_clone(higcit_celliterator *oldit) {
	higcit_bounding_box_extra *olditex = (higcit_bounding_box_extra *) oldit->ex;
	DECL_AND_ALLOC(higcit_celliterator, it, 1);
	DECL_AND_ALLOC(higcit_bounding_box_extra, itex, 1);
	it->bi.next = __higcit_bounding_box_inc;
	it->bi.destroy = __higcit_simple_destroy;
	it->bi.clone = __higcit_bounding_box_clone;
	it->ex = itex;
	POINT_ASSIGN(itex->lowpoint, olditex->lowpoint);
	POINT_ASSIGN(itex->highpoint, olditex->highpoint);
	itex->curlevel = olditex->curlevel;
	itex->top = olditex->top;
	for(int i = 0; i < MAXINNERCELLITERATORLEVELS; i++) {
		itex->cells[i] = olditex->cells[i];
		POINT_ASSIGN(itex->ppl[i], olditex->ppl[i]);
		POINT_ASSIGN(itex->pph[i], olditex->pph[i]);
		POINT_ASSIGN(itex->pps[i], olditex->pps[i]);
	}
	it->bi.c = oldit->bi.c;
	return it;
}

higcit_celliterator *higcit_create_bounding_box(hig_cell *cell, const Point lowpoint, const Point highpoint) {
	DECL_AND_ALLOC(higcit_celliterator, it, 1);
	DECL_AND_ALLOC(higcit_bounding_box_extra, itex, 1);

	it->bi.next = __higcit_bounding_box_inc;
	it->bi.destroy = __higcit_simple_destroy;
	it->bi.clone = __higcit_bounding_box_clone;
	it->ex = itex;

	POINT_ASSIGN(itex->lowpoint, lowpoint);
	POINT_ASSIGN(itex->highpoint, highpoint);
	itex->curlevel = -1;
	hig_cell *p = cell;

	if (p != NULL) {
		while (p->parent != NULL && !hig_is_bounding_box_inside_cell(p, lowpoint, highpoint)) {
			p = p->parent;
		}
		itex->top = p;

		if (hig_intersect_bounding_box(p, lowpoint, highpoint)) {
			p = __higcit_get_first_of_bounding_box(p, itex);
		} else {
			p = NULL;
		}
	}

	it->bi.c = p;
	return it;
}

static hig_cell *__higcit_filter_inc(higcit_celliterator *it) {
	higcit_filter_iterator_extra *itex = (higcit_filter_iterator_extra *)it->ex;
	hig_cell *c;
	higcit_filter_func f = itex->f;
	void *data = itex->data;
	do {
		c = higcit_nextcell(itex->it);
	} while (c != NULL && !f(c, data));
	it->bi.c = c;
	return it->bi.c;
}

static void __higcit_filter_destroy(higcit_celliterator *it) {
	higcit_filter_iterator_extra *itex = (higcit_filter_iterator_extra *) it->ex;
	higcit_destroy(itex->it);
	free(it->ex);
	free(it);
}

static higcit_celliterator *__higcit_filter_clone(higcit_celliterator *oldit) {
	higcit_filter_iterator_extra *olditex = (higcit_filter_iterator_extra *)oldit->ex;
	DECL_AND_ALLOC(higcit_filter_iterator_extra, itex, 1);
	DECL_AND_ALLOC(higcit_celliterator, it, 1);
	it->ex = itex;
	it->bi.next = __higcit_filter_inc;
	it->bi.destroy = __higcit_filter_destroy;
	it->bi.clone = __higcit_filter_clone;
	itex->it = higcit_clone(olditex->it);
	itex->f = olditex->f;
	itex->data = olditex->data;
	it->bi.c = oldit->bi.c;
	return it;
}

higcit_celliterator *higcit_create_filter(higcit_celliterator *itin, higcit_filter_func f, void *data) {
	DECL_AND_ALLOC(higcit_filter_iterator_extra, itex, 1);
	DECL_AND_ALLOC(higcit_celliterator, it, 1);
	it->ex = itex;
	it->bi.next = __higcit_filter_inc;
	it->bi.destroy = __higcit_filter_destroy;
	it->bi.clone = __higcit_filter_clone;
	itex->it = itin;
	itex->f = f;
	itex->data = data;
	hig_cell *c = higcit_getcell(itin);
	if (c != NULL && !f(c, data)) {
		c = __higcit_filter_inc(it);
	}
	it->bi.c = c;
	return it;
}

static hig_cell * __higcit_all_higtree_inc(higcit_celliterator *it) {
	higcit_all_higtree_extra *itex = (higcit_all_higtree_extra *) it->ex;
	hig_cell *p = it->bi.c;
	if (p != NULL) {
		it->bi.c = NULL;
		if (hig_get_number_of_children(p) > 0) {
			it->bi.c = hig_get_child(p, 0);
		} else {
			while (p != itex->root && hig_get_parent(p) != NULL) {
				int pos = p->posinparent + 1;
				if (pos < hig_get_number_of_children(hig_get_parent(p))) {
					it->bi.c = hig_get_child(hig_get_parent(p), pos);
					break;
				} else {
					p = p->parent;
				}
			}
		}
	}
	return it->bi.c;
}

static higcit_celliterator * __higcit_all_higtree_clone(higcit_celliterator *oldit) {
	higcit_all_higtree_extra *olditex = (higcit_all_higtree_extra *) oldit->ex;
	DECL_AND_ALLOC(higcit_celliterator, it, 1);
	DECL_AND_ALLOC(higcit_all_higtree_extra, itex, 1);
	it->ex = itex;
	itex->root = olditex->root;
	it->bi.next = __higcit_all_higtree_inc;
	it->bi.destroy = __higcit_simple_destroy;
	it->bi.clone = __higcit_all_higtree_clone;
	it->bi.c = oldit->bi.c;
	return it;
}

higcit_celliterator * higcit_create_all_higtree(hig_cell *cell) {
	DECL_AND_ALLOC(higcit_celliterator, it, 1);
	DECL_AND_ALLOC(higcit_all_higtree_extra, itex, 1);
	itex->root = cell;
	it->ex = itex;
	it->bi.next = __higcit_all_higtree_inc;
	it->bi.destroy = __higcit_simple_destroy;
	it->bi.clone = __higcit_all_higtree_clone;
	it->bi.c = cell;
	return it;
}

static hig_cell * __higcit_breadth_first_inc(higcit_celliterator *it) {
	higcit_breadth_first_extra *itex = (higcit_breadth_first_extra *) it->ex;
	hig_cell *p = it->bi.c;
	if(p == NULL)
		return NULL;

	// If there are children, there is at least one other
	// level that must be visited.
	if(!itex->has_more_level && hig_get_number_of_children(p) > 0) {
		itex->has_more_level = p;
	}

	// Now we need to get the next sibling...
	bool found = false;
	size_t depth = itex->cur_depth;

	// Outer loop: search for next node in correct level
	do {
		// Go right/up in the tree searching for siblings
		if(p == itex->root)
			break; // Has no sibling, it is the root!
		hig_cell *parent = hig_get_parent(p);
		int nextpos = p->posinparent + 1;
		if(nextpos < hig_get_number_of_children(parent)) {
			p = hig_get_child(parent, nextpos);
			// Found a possible sibiling...
			if(depth == itex->cur_depth) {
				found = true;
			} else if(hig_get_number_of_children(p) > 0) {
				// Inner loop: Go down in the tree searching for
				// a node in the correct depth.
				while(depth < itex->cur_depth) {
					if(hig_get_number_of_children(p) > 0) {
						p = hig_get_child(p, 0);
						++depth;
					} else {
						// Dead end... there is no node
						// of the required depth on this branch.
						// break the down search, and continue
						// to the sibling search.
						break;
					}
				}
				if(depth == itex->cur_depth)
					found = true;
			}
			// If not found, the next iteration will continue the
			// search to the next sibling.
		} else {
			// No more siblings. Go up in the tree.
			p = parent;
			--depth;
		}
	} while(!found);

	if(!found) {
		// Within current depth, no new sibilings were found.
		// Depth exausted, continue with next level if needed.
		if(itex->has_more_level) {
			assert(depth == 0);
			assert(p == itex->root);

			depth = ++itex->cur_depth;
			p = hig_get_child(itex->has_more_level, 0);
			itex->has_more_level = NULL;
			found = true;
		}
		else {
			// No more levels to seach for, iteration is done.
			found = false;
		}
	}

	if(found) {
		it->bi.c = p;
	} else {
		// All levels exausted, iteration completed.
		it->bi.c = NULL;
	}

	return it->bi.c;
}

static higcit_celliterator * __higcit_breadth_first_clone(higcit_celliterator *oldit) {
	higcit_breadth_first_extra *olditex = (higcit_breadth_first_extra *) oldit->ex;
	DECL_AND_ALLOC(higcit_celliterator, it, 1);
	DECL_AND_ALLOC(higcit_breadth_first_extra, itex, 1);

	*it = *oldit;
	*itex = *olditex;
	it->ex = itex;
	return it;
}

higcit_celliterator * higcit_create_breadth_first(hig_cell *cell) {
	DECL_AND_ALLOC(higcit_celliterator, it, 1);
	DECL_AND_ALLOC(higcit_breadth_first_extra, itex, 1);

	itex->root = cell;
	itex->cur_depth = 0;
	itex->has_more_level = NULL;
	it->ex = itex;
	it->bi.next = __higcit_breadth_first_inc;
	it->bi.destroy = __higcit_simple_destroy;
	it->bi.clone = __higcit_breadth_first_clone;
	it->bi.c = cell;
	return it;
}

static hig_cell *__higcit_sorted_inc(higcit_celliterator *it) {
	higcit_sorted_iterator_extra *itex = (higcit_sorted_iterator_extra *)it->ex;
	hig_cell *c;
	if (itex->pos < itex->numcells) {
		c = itex->cells[itex->pos];
		itex->pos++;
	} else {
		c = NULL;
	}
	it->bi.c = c;
	return it->bi.c;
}

static void __higcit_sorted_destroy(higcit_celliterator *it) {
	higcit_sorted_iterator_extra *itex = (higcit_sorted_iterator_extra *)it->ex;
	free(itex->cells);
	free(it->ex);
	free(it);
}

static higcit_celliterator *__higcit_sorted_clone(higcit_celliterator *oldit) {
	higcit_sorted_iterator_extra *olditex = (higcit_sorted_iterator_extra *)oldit->ex;
	DECL_AND_ALLOC(higcit_sorted_iterator_extra, itex, 1);
	DECL_AND_ALLOC(higcit_celliterator, it, 1);
	it->ex = itex;
	it->bi.next = __higcit_sorted_inc;
	it->bi.destroy = __higcit_sorted_destroy;
	it->bi.clone = __higcit_sorted_clone;
	itex->maxcells = olditex->maxcells;
	ALLOC(hig_cell *, itex->cells, itex->maxcells + 1);
	for(int i = 0; i < itex->maxcells + 1; i++) {
		itex->cells[i] = olditex->cells[i];
	}
	itex->pos = olditex->pos;
	itex->numcells = olditex->numcells;
	it->bi.c = oldit->bi.c;
	return it;
}

higcit_celliterator *higcit_create_sorted(higcit_celliterator *itin, int maxcells, higcit_lessthan_func f, void *data) {
	DECL_AND_ALLOC(higcit_sorted_iterator_extra, itex, 1);
	DECL_AND_ALLOC(higcit_celliterator, it, 1);
	it->ex = itex;
	it->bi.next = __higcit_sorted_inc;
	it->bi.destroy = __higcit_sorted_destroy;
	it->bi.clone = __higcit_sorted_clone;
	itex->maxcells = maxcells;
	ALLOC(hig_cell *, itex->cells, maxcells + 1);
	itex->pos = 0;
	itex->numcells = 0;
	for(;!higcit_isfinished(itin); higcit_nextcell(itin)) {
		hig_cell *c = higcit_getcell(itin);
		int j = itex->numcells - 1;
		while(j >= 0 && f(c, itex->cells[j], data)) {
			itex->cells[j+1] = itex->cells[j];
			j--;
		}
		itex->cells[j+1] = c;
		if (itex->numcells < itex->maxcells) {
			itex->numcells++;
		}
	}
	higcit_destroy(itin);
	__higcit_sorted_inc(it);
	return it;
}

static hig_cell *__higcit_concat_inc(higcit_celliterator *it) {
	higcit_concat_iterator_extra *itex = (higcit_concat_iterator_extra *)it->ex;
	it->bi.c = NULL;
	if (itex->posit < itex->numits) {
		higcit_nextcell(itex->its[itex->posit]);
		while(itex->posit < itex->numits && higcit_isfinished(itex->its[itex->posit])) {
			itex->posit++;
		}
		if (itex->posit < itex->numits) {
			it->bi.c = higcit_getcell(itex->its[itex->posit]);
		}
	}
	return it->bi.c;
}

static void __higcit_concat_destroy(higcit_celliterator *it) {
	higcit_concat_iterator_extra *itex = (higcit_concat_iterator_extra *)it->ex;
	for(int i = 0; i < itex->numits; i++) {
		higcit_destroy(itex->its[i]);
	}
	free(itex->its);
	free(it->ex);
	free(it);
}

static higcit_celliterator *__higcit_concat_clone(higcit_celliterator *oldit) {
	higcit_concat_iterator_extra *olditex = (higcit_concat_iterator_extra *)oldit->ex;
	DECL_AND_ALLOC(higcit_concat_iterator_extra, itex, 1);
	DECL_AND_ALLOC(higcit_celliterator, it, 1);
	it->ex = itex;
	it->bi.next = __higcit_concat_inc;
	it->bi.destroy = __higcit_concat_destroy;
	it->bi.clone = __higcit_concat_clone;
	itex->numits = olditex->numits;
	ALLOC(higcit_celliterator *, itex->its, itex->numits + 1);
	for(int i = 0; i < itex->numits; i++) {
		itex->its[i] = higcit_clone(olditex->its[i]);
	}
	itex->posit = olditex->posit;
	it->bi.c = oldit->bi.c;
	return it;
}

higcit_celliterator *higcit_create_concat(higcit_celliterator *its[], int numits) {
	DECL_AND_ALLOC(higcit_concat_iterator_extra, itex, 1);
	DECL_AND_ALLOC(higcit_celliterator, it, 1);
	it->ex = itex;
	it->bi.next = __higcit_concat_inc;
	it->bi.destroy = __higcit_concat_destroy;
	it->bi.clone = __higcit_concat_clone;
	itex->numits = numits;
	ALLOC(higcit_celliterator *, itex->its, numits + 1);
	for(int i = 0; i < itex->numits; i++) {
		itex->its[i] = its[i];
	}
	itex->posit = 0;
	it->bi.c = NULL;
	for(; itex->posit < itex->numits; itex->posit++) {
		if (!higcit_isfinished(itex->its[itex->posit])) {
			it->bi.c = higcit_getcell(itex->its[itex->posit]);
			break;
		}
	}
	return it;
}

hig_facet * higfit_getfacet(higfit_facetiterator *fit) {
	return &fit->bi.facet;
}

higfit_facetiterator * higfit_clone(higfit_facetiterator *it) {
	assert(it->bi.clone != NULL);
	return it->bi.clone(it);
}

void higfit_nextfacet(higfit_facetiterator *fit) {
	assert(fit->bi.next != NULL);
	fit->bi.next(fit);
}

void higfit_destroy(higfit_facetiterator *fit) {
	assert(fit->bi.destroy != NULL);
	fit->bi.destroy(fit);
}

int higfit_isfinished(higfit_facetiterator *fit) {
	return fit->bi.facet.c == NULL;
}

int higfit_isvalid_facet(higfit_facetiterator *fit) {
	if (!higcit_isfinished(fit->bi.it)) {
		if (fit->bi.facet.dim < DIM && fit->bi.dimofinterest[fit->bi.facet.dim]) {
			if (fit->bi.facet.dir < 2) {
				if (fit->bi.facet.c != NULL) {
					int fid = hig_get_fid_with_dim_dir(fit->bi.facet.c, fit->bi.facet.dim, fit->bi.facet.dir);
					if (fid > 0) {
						return 1;
					}
				}
			}
		}
	}
	return 0;
}

hig_facet * __higfit_inc_facet(higfit_facetiterator *fit) {
	if (fit->bi.facet.dir < 1) {
		fit->bi.facet.dir++;
	} else {
		fit->bi.facet.dir = 0;
		if (fit->bi.facet.dim < DIM - 1) {
			fit->bi.facet.dim++;
		} else {
			fit->bi.facet.dim = 0;
			if (!higcit_isfinished(fit->bi.it)) {
				higcit_nextcell(fit->bi.it);
				fit->bi.facet.c = higcit_getcell(fit->bi.it);
			}
		}
	}
	return &fit->bi.facet;
}

hig_facet * __higfit_allfacets_next(higfit_facetiterator *fit) {
	do {
		__higfit_inc_facet(fit);
	} while(!higcit_isfinished(fit->bi.it) && !higfit_isvalid_facet(fit));
	return &fit->bi.facet;
}

void __higfit_allfacets_destroy(higfit_facetiterator *fit) {
	higcit_destroy(fit->bi.it);
	free(fit->ex);
	free(fit);
}

higfit_facetiterator *__higfit_allfacets_clone(higfit_facetiterator *oldfit) {
	DECL_AND_ALLOC(higfit_facetiterator, fit, 1);
	DECL_AND_ALLOC(higfit_allfacets_extra, fitex, 1);
	fit->ex = fitex;
	fit->bi.it = higcit_clone(oldfit->bi.it);
	POINT_ASSIGN(fit->bi.dimofinterest, oldfit->bi.dimofinterest);
	fit->bi.facet.dim = oldfit->bi.facet.dim;
	fit->bi.facet.dir = oldfit->bi.facet.dir;
	fit->bi.facet.c = oldfit->bi.facet.c;
	fit->bi.next = oldfit->bi.next;
	fit->bi.destroy = oldfit->bi.destroy;
	fit->bi.clone = oldfit->bi.clone;
	return fit;
}

higfit_facetiterator * higfit_create_allfacets(higcit_celliterator *it, int dimofinterest[DIM]) {
	DECL_AND_ALLOC(higfit_facetiterator, fit, 1);
	ALLOC(higfit_allfacets_extra, fit->ex, 1);
	fit->bi.it = it;
	POINT_ASSIGN(fit->bi.dimofinterest, dimofinterest);
	fit->bi.facet.dim = 0;
	fit->bi.facet.dir = 0;
	fit->bi.facet.c = higcit_getcell(it);
	fit->bi.next = __higfit_allfacets_next;
	fit->bi.destroy = __higfit_allfacets_destroy;
	fit->bi.clone = __higfit_allfacets_clone;
	while(!higcit_isfinished(fit->bi.it) && !higfit_isvalid_facet(fit)) {
		__higfit_inc_facet(fit);
	}
	return fit;
}

hig_facet * __higfit_bounding_box_next(higfit_facetiterator *fit) {
	higfit_bounding_box_extra * fitex = (higfit_bounding_box_extra *) fit->ex;
	do {
		__higfit_inc_facet(fit);
	} while(!higcit_isfinished(fit->bi.it) && !(higfit_isvalid_facet(fit) && hig_facet_intersect_bounding_box(&fit->bi.facet, fitex->lowpoint, fitex->highpoint)));
	return &fit->bi.facet;
}

void __higfit_bounding_box_destroy(higfit_facetiterator *fit) {
	higfit_bounding_box_extra * fitex = (higfit_bounding_box_extra *) fit->ex;
	higcit_destroy(fit->bi.it);
	free(fitex);
	free(fit);
}

higfit_facetiterator *__higfit_bounding_box_clone(higfit_facetiterator *oldfit) {
	higfit_bounding_box_extra * oldfitex = (higfit_bounding_box_extra *) oldfit->ex;
	DECL_AND_ALLOC(higfit_facetiterator, fit, 1);
	DECL_AND_ALLOC(higfit_bounding_box_extra, fitex, 1);
	POINT_ASSIGN(fitex->lowpoint, oldfitex->lowpoint);
	POINT_ASSIGN(fitex->highpoint, oldfitex->highpoint);
	fit->ex = fitex;
	fit->bi.it = higcit_clone(oldfit->bi.it);
	POINT_ASSIGN(fit->bi.dimofinterest, oldfit->bi.dimofinterest);
	fit->bi.facet.dim = oldfit->bi.facet.dim;
	fit->bi.facet.dir = oldfit->bi.facet.dir;
	fit->bi.facet.c = oldfit->bi.facet.c;
	fit->bi.next = oldfit->bi.next;
	fit->bi.destroy = oldfit->bi.destroy;
	fit->bi.clone = oldfit->bi.clone;
	return fit;
}

higfit_facetiterator * higfit_create_bounding_box_facets(hig_cell *root, int dimofinterest[DIM],
	const Point l, const Point h)
{
	DECL_AND_ALLOC(higfit_facetiterator, fit, 1);
	DECL_AND_ALLOC(higfit_bounding_box_extra, fitex, 1);
	POINT_ASSIGN(fitex->lowpoint, l);
	POINT_ASSIGN(fitex->highpoint, h);
	fit->ex = fitex;
	fit->bi.it = higcit_create_bounding_box(root, l, h);
	POINT_ASSIGN(fit->bi.dimofinterest, dimofinterest);
	fit->bi.facet.dim = 0;
	fit->bi.facet.dir = 0;
	fit->bi.facet.c = higcit_getcell(fit->bi.it);
	fit->bi.next = __higfit_bounding_box_next;
	fit->bi.destroy = __higfit_bounding_box_destroy;
	fit->bi.clone = __higfit_bounding_box_clone;
	while(!higcit_isfinished(fit->bi.it) && !(higfit_isvalid_facet(fit) && hig_facet_intersect_bounding_box(&fit->bi.facet, fitex->lowpoint, fitex->highpoint))) {
		__higfit_inc_facet(fit);
	}
	return fit;
}

static hig_facet *__higfit_concat_inc(higfit_facetiterator *it) {
	higfit_concat_iterator_extra *itex = (higfit_concat_iterator_extra *)it->ex;
	it->bi.facet.c = NULL;
	if (itex->posit < itex->numits) {
		higfit_nextfacet(itex->its[itex->posit]);
		while(itex->posit < itex->numits && higfit_isfinished(itex->its[itex->posit])) {
			itex->posit++;
		}
		if (itex->posit < itex->numits) {
			it->bi.facet = *higfit_getfacet(itex->its[itex->posit]);
		}
	}
	return &it->bi.facet;
}

static void __higfit_concat_destroy(higfit_facetiterator *it) {
	higfit_concat_iterator_extra *itex = (higfit_concat_iterator_extra *)it->ex;
	for(int i = 0; i < itex->numits; i++) {
		higfit_destroy(itex->its[i]);
	}
	free(itex->its);
	free(it->ex);
	free(it);
}

static higfit_facetiterator *__higfit_concat_clone(higfit_facetiterator *oldit) {
	higfit_concat_iterator_extra *olditex = (higfit_concat_iterator_extra *)oldit->ex;
	DECL_AND_ALLOC(higfit_concat_iterator_extra, itex, 1);
	DECL_AND_ALLOC(higfit_facetiterator, it, 1);
	it->ex = itex;
	it->bi.next = __higfit_concat_inc;
	it->bi.destroy = __higfit_concat_destroy;
	it->bi.clone = __higfit_concat_clone;
	itex->numits = olditex->numits;
	ALLOC(higfit_facetiterator *, itex->its, itex->numits + 1);
	for(int i = 0; i < itex->numits; i++) {
		itex->its[i] = higfit_clone(olditex->its[i]);
	}
	itex->posit = olditex->posit;
	it->bi.facet = oldit->bi.facet;
	return it;
}

higfit_facetiterator *higfit_create_concat(higfit_facetiterator *its[], int numits) {
	DECL_AND_ALLOC(higfit_concat_iterator_extra, itex, 1);
	DECL_AND_ALLOC(higfit_facetiterator, it, 1);
	it->ex = itex;
	it->bi.next = __higfit_concat_inc;
	it->bi.destroy = __higfit_concat_destroy;
	it->bi.clone = __higfit_concat_clone;
	itex->numits = numits;
	ALLOC(higfit_facetiterator *, itex->its, numits + 1);
	for(int i = 0; i < itex->numits; i++) {
		itex->its[i] = its[i];
	}
	itex->posit = 0;
	it->bi.facet.c = NULL;
	for(; itex->posit < itex->numits; itex->posit++) {
		if (!higfit_isfinished(itex->its[itex->posit])) {
			it->bi.facet = *higfit_getfacet(itex->its[itex->posit]);
			break;
		}
	}
	return it;
}


long higfit_count(higfit_facetiterator *it) {
	long cnt = 0;
	for(;!higfit_isfinished(it); higfit_nextfacet(it)) {
		cnt++;
	}
	return cnt;
}
