#ifndef MTREE_ITERATOR_INTERNAL_H
#define MTREE_ITERATOR_INTERNAL_H

#include <stdbool.h>
#include "higtree.h"

struct higcit_celliterator;

//! Pointer to a function which advances the iterator.
typedef hig_cell * (*higcit_nextfunc)(struct higcit_celliterator *);

//! Pointer to a function which destroys the iterator.
typedef void (*higcit_destroyfunc)(struct higcit_celliterator *);

//! Pointer to a function which clones the iterator.
typedef struct higcit_celliterator * (*higcit_clonefunc)(struct higcit_celliterator *);

//! The data which are shared by all cell iterators.
typedef struct higcit_basiciterator {
	hig_cell *c;
	higcit_nextfunc next;
	higcit_destroyfunc destroy;
	higcit_clonefunc clone;
} higcit_basiciterator;

//! A cell iterator contains a basic iterator plus a set of extra information, which depends on the iterator type.
typedef struct higcit_celliterator {
	higcit_basiciterator bi;
	void *ex;
} higcit_celliterator;

//! \brief Extra data for the iterator which returns all leaves.
typedef struct higcit_all_leaves_extra {
  hig_cell *root;
} higcit_all_leaves_extra;

#define MAXINNERCELLITERATORLEVELS 100
//! \brief Extra data for the iterator which returns all leaves intersecting a bounding box.
typedef struct higcit_bounding_box_extra {
  Point lowpoint;
  Point highpoint;
  hig_cell *top;
  hig_cell *cells[MAXINNERCELLITERATORLEVELS];
  int ppl[MAXINNERCELLITERATORLEVELS][DIM];
  int pph[MAXINNERCELLITERATORLEVELS][DIM];
  int pps[MAXINNERCELLITERATORLEVELS][DIM];
  int curlevel;
} higcit_bounding_box_extra;

//! Extra data for the iterator which returns all neighbours of a cell.
typedef struct higcit_neighbours_extra {
	hig_cell *cell;
	higcit_celliterator *it;
} higcit_neighbours_extra;

//! Type of a function which returns whether a cell should be considered by a filter iterator.
typedef bool (* higcit_filter_func)(hig_cell *, void *);

//! Extra data for the filter iterator.
typedef struct higcit_filter_iterator_extra {
	higcit_celliterator *it;
	higcit_filter_func f;
	void *data;
} higcit_filter_iterator_extra;

//! Extra data for the iterator which returns all cells in the higtree (not only the leaves)
typedef struct higcit_all_higtree_extra {
  hig_cell *root;
} higcit_all_higtree_extra;

//! Extra data for the iterator which returns all cells in the higtree in breadth-first order
typedef struct higcit_breadth_first_extra {
  hig_cell *root;
  hig_cell *has_more_level;
  size_t cur_depth;
} higcit_breadth_first_extra;

//! Type of a function which return which compares to cells.
typedef bool (* higcit_lessthan_func)(hig_cell *, hig_cell *, void *);

//! Extra data for the sorted iterator.
typedef struct higcit_sorted_iterator_extra {
	hig_cell **cells;
	int numcells;
	int maxcells;
	int pos;
} higcit_sorted_iterator_extra;

//! Extra data for the iterator which concatenates the result of other iterators.
typedef struct higcit_concat_iterator_extra {
	int numits;
	int posit;
	higcit_celliterator **its;
} higcit_concat_iterator_extra;

struct higfit_facetiterator;

// Facet iterators
//! Pointer to a function which advances the facet iterator.
typedef hig_facet * (*higfit_nextfacetfunc)(struct higfit_facetiterator *);

//! Pointer to a function which destroys the facet iterator.
typedef void (*higfit_destroyfacetfunc)(struct higfit_facetiterator *);

//! Pointer to a function which clones the facet iterator.
typedef struct higfit_facetiterator * (*higfit_clonefacetfunc)(struct higfit_facetiterator *);

//! Structure which contains the basic information of the facet iterator.
typedef struct higfit_basicfacetiterator {
	higcit_celliterator *it;
	int dimofinterest[DIM];
	hig_facet facet;
	higfit_nextfacetfunc next;
	higfit_destroyfacetfunc destroy;
	higfit_clonefacetfunc clone;
} higfit_basicfacetiterator;

//! A facet iterator contains a basic iterator plus a set of extra information, which depends on the iterator type.
typedef struct higfit_facetiterator {
	higfit_basicfacetiterator bi;
	void *ex;
} higfit_facetiterator;

//! \brief Extra data for the iterator which returns all facets. Currently, it is empty.
typedef struct higfit_allfacets_extra {
} higfit_allfacets_extra;

//! \brief Extra data for the iterator which returns all cells intersecting a bounding box.
typedef struct higfit_bounding_box_extra {
	Point lowpoint, highpoint;
} higfit_bounding_box_extra;


typedef struct higfit_concat_iterator_extra {
	int numits;
	int posit;
	higfit_facetiterator **its;
} higfit_concat_iterator_extra;

#endif
