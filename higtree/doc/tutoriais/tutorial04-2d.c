#include<stdio.h>
#include<stdlib.h>

#include "higtree.h"
#include "higtree-io.h"
#include "higtree-iterator.h"


int main (int argc, char *argv[]) {
    higtree_initialize(&argc, &argv);
    Point lo, hi;
    POINT_ASSIGN_SCALAR(lo, -1.0);
    POINT_ASSIGN_SCALAR(hi, 1.0);

    //creates the root cell
    hig_cell *root = hig_create_root(lo, hi);

    // refining into 2 cells in both directions
    int nc[2] = {2,2};
    hig_refine_uniform(root, nc);

    //refines the cell with refine_point into more 2 cell in each direction
    Point refine_point;
    POINT_ASSIGN_SCALAR(refine_point, 0.5);
    hig_cell *c = hig_get_cell_with_point(root, refine_point);
    hig_refine_uniform(c, nc);
    FILE *fd = fopen("VTKS/tut04.vtk", "w");
    higio_print_in_vtk2d(fd, root);
    fclose(fd);
    
    higcit_celliterator *cit; /** creates a cell iterator */
    printf("centers of the cells:\n");
    //iterates over all the leaves of the tree, i.e., all the cell centers
    for(cit = higcit_create_all_leaves(root); !higcit_isfinished(cit); higcit_nextcell(cit)) {
        hig_cell *ccell = higcit_getcell(cit); /** cell pointed by the iterator currently */
        Point center;
        hig_get_center(ccell, center); /** gets the center of the cell */
        int id = hig_get_cid(ccell);
        printf("cell: %d \t(x= %.2f, y= %.2f)\n", id, center[0], center[1]);
    }
    /** deallocates iterator */
    higcit_destroy(cit);


    higfit_facetiterator *fit; /** creates a face iterator */
    cit = higcit_create_all_leaves(root); /** the facets an iterator runs through depends on the cells of a cell iterator */

    printf("\ncenters of the facets, \tdirection of the facet\n");

    /** Iterates through all the facets in all directions */
    int interest_dimensions[2] = {1,1};
    for(fit = higfit_create_allfacets(cit, interest_dimensions); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
        hig_facet *f = higfit_getfacet(fit); /** gets the facet of the iterator */
        Point fcenter;
        hig_get_facet_center(f, fcenter); /** gets the coordinates of the facet center */

        int dim = hig_get_facet_dim(f);
        int dir = hig_get_facet_dir(f);
        int id = hig_get_fid(f);

        printf("id: %d \t(x= %.2f, y= %.2f), \tdim: %d \tdir = %d\n", id, fcenter[0], fcenter[1], dim, dir);
    }
    
    /** deallocates iterator (already destroys cell iterator)*/
    higfit_destroy(fit);
}
