#ifndef MTREE_H
#define MTREE_H

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#include "coord.h"
#include "rect.h"
#include "types.h"
#include "uniqueid.h"

//! Indicates that each cell has a unique ID
#define HASCELLID
//! Indicates that each facet has a unique ID.
#define HASFACETID
// Indicates that each vertex has a unique ID.
//#define HASVERTEXID


#ifdef HASCELLID
//! there is one id per cell
#define NUMCELLIDS 1
#else
#define NUMCELLIDS 0
#endif

#ifdef HASFACETID
//! there are 2 ids per dimension per cell
#define NUMFACETIDS (2*DIM)
#else
#define NUMFACETIDS 0
#endif

#ifdef HASVERTEXID
//! there are 2^DIM ids per cell
#define NUMVERTEXIDS (1<<DIM)
#else
#define NUMVERTEXIDS 0
#endif

#define FIRSTCELLID 0
#define FIRSTFACETID (FIRSTCELLID + NUMCELLIDS)
#define FIRSTVERTEXID (FIRSTFACETID + NUMFACETIDS)
#define NUMIDS (NUMCELLIDS+NUMFACETIDS+NUMVERTEXIDS)

//! This function can be called in the beginning of a program, to make sure that HASFACETID is defined.
//! It aborts the program if it is not.
int hig_requires_facet_ids();
//! This function can be called in the beginning of a program, to make sure that HASVERTEXID is defined
//! It aborts the program if it is not.
int hig_requires_vertex_ids();

//! \brief hig_cell represents a cell in the mesh. <BR>
//! Each cell is defined by two points (lowpoint and higpoint) and a number of division in each dimension. <BR>
//! It knows which cell is its parent (in the tree) and which is its position among its parent's children. <BR>
//! It also knows its children.<BR>
//! It keeps a set of ids (the number depends on whether facet and vertex ids are required).
typedef struct hig_cell {
	Point lowpoint;
	Point highpoint;
	int numcells[DIM];
	struct hig_cell * parent;
	int posinparent;
	struct hig_cell **children;
	uniqueid ids[NUMIDS];
} hig_cell;

//! \brief hig_facet represents a facet in the mesh. <BR>
//! Each facet is defined by a cell, a dimension and a direction (0 = left/bottom/back, 1 = right/top/front).
typedef struct hig_facet {
	hig_cell *c;
	int dim;
	int dir;
} hig_facet;

//! Given the number of cell division in each dimension and the position of the cell in each dimension, computes the the relative position of the cell among all its parent children.<BR>
//! It can be thought of an enumeration of the cells, allowing to represent a object of DIM dimensions in a unidimensional array.
//! See also hig_tobase.
int hig_toposition(const int numcells[DIM], const int pp[DIM]);

//! Converts the position of a cell in a unidimensional array in the its coordinates in the DIM dimensions.
void hig_tobase(int position, const int numcells[DIM], int pp[DIM]);

//! Computes the lowpoint of a cell in a regular grid.
void hig_translate_coord(Point lowpoint, Point delta, int pp[DIM], int numcells[DIM], Point chlp);

//! Computes the euclidean distance of a point p to the center of a cell.
real hig_distance_from_center(Point p, hig_cell *c);

//! Returns the square of the smallest distance from cell to point p.
// If p is contained in cell, zero is returned.
real hig_sq_distance_to_point(hig_cell *c, const Point p);

//! Computes the center of a cell: center = (lowpoint + highpoint) / 2
void hig_get_center(hig_cell *c, Point center);

//! Computes the delta of a cell: delta = (highpoint - lowpoint)
void hig_get_delta(hig_cell *c, Point delta);

//! Gets the lowpoint of a cell.
void hig_get_lowpoint(hig_cell *c, Point p);

//! Get the highpoint of a cell.
void hig_get_highpoint(hig_cell *c, Point p);

//! Get the bounding box of a cell.
void hig_get_bounding_box(hig_cell *c, Rect *box);

//! Get the i-th id stored in the cell. It can be a cell, facet or vertex id, depending on i.<BR>
//! Use hig_get_cid and hig_get_fid*
uniqueid hig_get_id(hig_cell *c, int i);

//! Get the cell id.
uniqueid hig_get_cid(hig_cell *c);

//! Get the facet id.
uniqueid hig_get_fid(hig_facet *f);

//! Get the facet id, where the facet is represented by a cell, a dimension and a direction.
uniqueid hig_get_fid_with_dim_dir(hig_cell *c, int dim, int dir);

//! Set the facet id, where the facet is represented by a cell, a dimension and a direction.
void hig_set_fid_with_dim_dir(hig_cell *c, int dim, int dir, uniqueid id);

//! Gets the vertex id.
uniqueid hig_get_vid(hig_cell *c, int i);

//! Computes how many children the cell has.
int hig_get_number_of_children(hig_cell *c);

//! Gets the child given its coordinates (among its parent's children).
hig_cell *hig_get_child_in_grid(hig_cell *c, const int pp[DIM]);

//! Gets the pos-th child of a cell.
hig_cell *hig_get_child(hig_cell *c, int pos);

//! Gets the cell's parent.
hig_cell *hig_get_parent(hig_cell *c);

//! Gets the number of division in each direction. How many cells are per dimension.
void hig_get_cells_per_dim(hig_cell *c, int nc[DIM]);

//! TODO
hig_cell *hig_get_neighbour(hig_cell *c, int dir[DIM]);

int hig_get_narrowest_dim(hig_cell *c);

//! Converts a cell relative coordinate (between 0 and 1) to absolute coordinates.
void hig_get_relative_coord(hig_cell *c, Point relcoord, Point coord);

//! Creates a hig_cell, with the given low and high points. The cell has no parent. Thus, it is a root of the tree. It also has no children.
hig_cell * hig_create_root(Point lowpoint, Point highpoint);

//! Creates a hig_cell, with given low point and delta (highpoint = lowpoint + delta), parent and position among its parent's children. The child is assigned to the parent.
hig_cell * hig_create_leaf(Point lowpoint, Point delta, hig_cell * parent, int pp);

//! Creates an empty hig_cell, to be filled with hig_fill_empty().
hig_cell * hig_create_empty_leaf(hig_cell *parent, int pp);

//! Refines a cell, creating a grid where which cell has the same delta (the delta can be different in different dimensions).
hig_cell * hig_refine_uniform(hig_cell * parent, int numcells[DIM]);
//hig_cell * hig_refine(hig_cell * cell, int numcells[DIM], Point lowpoint[], Point highpoint[]);

//! Refines a cell, but does not creates the children.
// The caller is responsible to filling the children with hig_create_leaf() or
// hig_create_empty_leaf() calls, up to the size returned by hig_get_number_of_children().
void hig_refine_empty(hig_cell * parent, int numcells[DIM]);

//! Set the properties of an empty cell, as created by hig_create_empty_leaf().
void hig_fill_empty(hig_cell *cell, Point lowpoint, Point highpoint);

//! Destroys a cell and all its children, recursively.
void hig_destroy(hig_cell *cell);

//! Computes in the coordinates of the children of cell which contains point.
void hig_get_cell_coords_of_point(hig_cell *cell, const Point point, int pp[DIM]);

//! Computes which cell contains a given point. Starts the search from parent and moves down the tree.
hig_cell * hig_get_cell_with_point(hig_cell *parent, const Point point);

//! Computes which cell contains a given point. Starts from a given cell, moving up if necessary, then moving down the tree.
hig_cell * hig_get_cell_with_point_from_another_cell(hig_cell *cell, Point point);

//! Merges all children of a cell, deleting them. After that, p is a leave.
hig_cell * hig_merge_children(hig_cell *p);

//! Computes the number of leaves of the tree, starting from parent.
long hig_get_number_of_leaves(hig_cell *parent);

//! Determines whether the bounding box defined by lowpoint and highpoint is completely inside cell.
bool hig_is_bounding_box_inside_cell(hig_cell *cell, const Point lowpoint, const Point highpoint);

//! Determines whether the bounding box defined by lowpoint and highpoint intersects cell.
bool hig_intersect_bounding_box(hig_cell *cell, const Point lowpoint, const Point highpoint);

//! Obtains the leftmost cell of a tree, recursively.
//! The leftmost cell of a leaf is the leaf itself. Otherwise, it is the leftmost cell of its first child.
hig_cell * hig_get_first_cell_recursive(hig_cell *cell);

//! Verifies whether the tree rooted at root is ok. It checks whether the children defines a grid, there is no overlapping nor gaps between cells.
bool hig_check_conformable(hig_cell *root);

//! Estimates the amount of memory used by a tree rooted at parent.
long int hig_memory_used(hig_cell *parent);

//! Determines whether two trees, rooted at c1 and c2, are equal (only ids are not checked)
bool hig_equal(hig_cell *c1, hig_cell *c2);

//! Copies the refiment of tree rooted at *from* into the tree rooted at *to*. The *to* tree is assumed to be a leaf, i.e., without refiment.
//! Moreover, the refiment should be regular. After executing this function, the trees rooted at *from* and *to* are equal (except for ids).
void hig_copy_regular_refinement(hig_cell *from, hig_cell *to);

//! Creates a copy of the tree rooted at *from*. Ids are not copied; thus, the new tree has unique ids.
hig_cell *hig_clone(hig_cell *from);

//! Adjusts the ids of facets, using the following principles:<BR>
//! 1) if the same facet which is shared by two cells, the cell with the lowest coordinate keeps the id and the other cell sets it to -id. <BR>
//! 2) if a facet f1 of cell c1 includes the facet f2 of cell c2 (i.e., f2 is within the bounding box of f1), then c2 keeps the id of f2 and c1 sets the id of f1 to 0 (an invalid id) <BR>
//! It assumes that a facet area is can partitioned into the areas of a set of facets.
void hig_adjust_facet_ids(hig_cell *root);

//! Gets the dimension the facet. The coordinates of this dimension for low and high points are equal.
int hig_get_facet_dim(hig_facet *facet);

//! Gets the direction of a facet, that is, whether it is the left (low), value 0, or the right (high) facet, value 1.
int hig_get_facet_dir(hig_facet *facet);

//! Gets the cell which contains the facet.
hig_cell *hig_get_facet_cell(hig_facet *facet);

//! Gets the facet id.
uniqueid hig_get_facet_id(hig_facet *facet);

//! Gets the low point of the facet.
void hig_get_facet_lowpoint(hig_facet *facet, Point l);

//! Gets the high point of the facet.
void hig_get_facet_highpoint(hig_facet *facet, Point h);

//! Gets the center of the facet.
void hig_get_facet_center(hig_facet *facet, Point center);

//! Gets the delta of the facet.
void hig_get_facet_delta(hig_facet *facet, Point delta);

//! Auxiliary function. Determines whether two interval intersect.
int hig_intersect_intervals(double l1, double h1, double l2, double h2);

//! Determines whether a facet intersects a bounding box.
bool hig_facet_intersect_bounding_box(hig_facet *facet, const Point l, const Point h);

//! Computes the facet which contains a given point.
int hig_get_facet_with_point(hig_cell *root, int dim, const Point x, hig_facet *facet);


#endif
