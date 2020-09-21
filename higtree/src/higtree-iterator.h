
#ifndef MTREE_ITERATOR_H
#define MTREE_ITERATOR_H

#include "types.h"
#include "higtree.h"
#include "higtree-iterator-internal.h"

//! Gets the cell currently associated with the iteration it.
hig_cell * higcit_getcell(higcit_celliterator *it);

//! Advances to the next cell.
hig_cell * higcit_nextcell(higcit_celliterator *it);

//! Clones the iterator. The clone version replicates all relevant information of the cloned iterator.
higcit_celliterator * higcit_clone(higcit_celliterator *it);

//! Checks whether the iterator has already reached the last cell.
bool higcit_isfinished(higcit_celliterator *it);

//! Destroy the cell iterator.
void higcit_destroy(higcit_celliterator *it);

//! Counts how many cells the iterator returns. The iterator is advanced; thus, after calling this function, the iterator is finished.
long higcit_count(higcit_celliterator *it);

//! Counts how many cells a clone of the iterator returns. Thus, the iterator is not advanced.
long higcit_count_without_advancing(higcit_celliterator *it);

//! Creates an iterator which returns all leaves (i.e., cells without children) of the hig-tree rooted at c.
higcit_celliterator *higcit_create_all_leaves(hig_cell *c);

//! Creates an iterator which returns all leaves which intersect with the bounding box defined by lowpoint and highpoint.
higcit_celliterator *higcit_create_bounding_box(hig_cell *cell, const Point lowpoint, const Point highpoint);

//! Creates an iterator which returns all leaves which touch the cell, i.e., its neighbours.
higcit_celliterator *higcit_create_neighbours(hig_cell *cell);

//! Creates an iterator which returns all cells returned by the iterator itin and for that the function f returns true.
higcit_celliterator *higcit_create_filter(higcit_celliterator *itin, higcit_filter_func f, void *data);

//! Creates an iterator which returns all cells of the higtree rooted at c, regardless if it is a leaf or not.
higcit_celliterator *higcit_create_all_higtree(hig_cell *c);

//! Creates an iterator which returns all cells of the higtree rooted at c at breadth-first order.
//
// To save memory, it uses "Iterative deepening depth-first search", that is equivalent to
// a breadth first search.
higcit_celliterator *higcit_create_breadth_first(hig_cell *c);

//! Creates an iterator which returns the cells returned by the iterator itin, sorted according to the function f.
higcit_celliterator *higcit_create_sorted(higcit_celliterator *itin, int maxcells, higcit_lessthan_func f, void *data);

//! Creates an iterator which returns all cells returned by the iterators in its.
higcit_celliterator *higcit_create_concat(higcit_celliterator *its[], int numits);

//! Gets the current facet associated with the facet iterator fit.
hig_facet * higfit_getfacet(higfit_facetiterator *fit);

//! Clones the facet iterator. The clone version replicates all relevant information of the cloned iterator.
higfit_facetiterator * higfit_clone(higfit_facetiterator *it);

//! Advances to the next facet.
void higfit_nextfacet(higfit_facetiterator *fit);

//! Destroys the facet iterator.
void higfit_destroy(higfit_facetiterator *fit);

//! Checks whether the facet iterator has already reached the last facet.
int higfit_isfinished(higfit_facetiterator *fit);

//! Counts how many facets is returned by the facet iterator.
long higfit_count(higfit_facetiterator *fit);

//! \brief Creates a facet iterator which return all facets of the cells returned by the cell iterator it, considering the dimensions of interest.
//! The dimensions of interest contains 1 for the dimension whose facets are to be returned and 0 for the others. More than one dimension
//! can be defined for a given facet iterator.
higfit_facetiterator * higfit_create_allfacets(higcit_celliterator *it, int dimofinterest[DIM]);

//! Creates a facet iterator which return all facets of the cells returned by the cell iterator it and that intersect the bounding box defined by the point l and h, considering the dimensions of interest.
higfit_facetiterator * higfit_create_bounding_box_facets(hig_cell *root, int dimofinterest[DIM], const Point l, const Point h);

//! Creates an iterator which returns all cells returned by the iterators in its.
higfit_facetiterator *higfit_create_concat(higfit_facetiterator *its[], int numits);

#endif
