#include <stdio.h>
#include<stdlib.h>

#include "coord.h"
#include "higtree-iterator-internal.h"
#include "higtree-iterator.h"
#include "higtree.h"
#include "higtree-io.h"
#include "domain.h"
#include "mapper.h"
#include "utils.h"

/**
 * @file creates a domain with 2 HiG-Trees and dirichlet boundary conditions (Neumann bc are similar)
 * 
 */


int main (int argc, char *argv[]) {
    higtree_initialize(&argc, &argv);

    int nc[2] = {10, 10};

    Point lo, hi;

    /* Creates the first hig-tree*/
    POINT_ASSIGN_SCALAR(lo, -1.0);
    POINT_ASSIGN_SCALAR(hi, 1.0);
    hig_cell *root1 = hig_create_root(lo, hi);
    hig_refine_uniform(root1, nc);

    /* writes to vtk for visualization*/
    FILE *f1 = fopen("VTKS/tut05sd1.vtk", "w");
    higio_print_in_vtk2d(f1, root1);


    /* creates the second hig-tree*/
    POINT_ASSIGN_INTS(nc, 5, 5);
    POINT_ASSIGN_REALS(lo, 1.0,-1.0);
    POINT_ASSIGN_REALS(hi, 2.0, 1.0);
    hig_cell *root2 = hig_create_root(lo, hi);
    hig_refine_uniform(root2, nc);

    /*Writes to vtk for visualization*/
    FILE *f2 = fopen("VTKS/tut05sd2.vtk", "w");
    higio_print_in_vtk2d(f2, root2);
    

    higcit_celliterator *it;

    sim_boundary *sb[4];
    /**
     * Boundary conditions are hig-trees with width defined as the smallest delta of a cell.
     * It's a fine line. Every point inside it is part of the bc. The hig cell of the bc is made in a
     * way that its cell's center coincide with the facet center of desired cells.
    */
    for(int i = 0; i<4; i++){ //iterates in all 4 directions that need a boundary condition
        /* Every instance of sim_boundary (bc) needs a mapper*/
        mp_mapper *bcm = mp_create();
        switch (i) {
            case 0: //left boundary
                POINT_ASSIGN_INTS(nc, 1, 10);
                POINT_ASSIGN_REALS(lo, -1.0 - EPSDELTA,-1.0);
                POINT_ASSIGN_REALS(hi, -1.0, 1.0 + EPSDELTA);
                break;
            case 1: //top boundary
                POINT_ASSIGN_INTS(nc, 15, 1);
                POINT_ASSIGN_REALS(lo, -1.0,1.0 - EPSDELTA);
                POINT_ASSIGN_REALS(hi, 2.0 + EPSDELTA, 1.0);
                break;
            case 2: //right boundary
                POINT_ASSIGN_INTS(nc, 1, 5);
                POINT_ASSIGN_REALS(lo, 2.0 - EPSDELTA,-1.0);
                POINT_ASSIGN_REALS(hi, 2.0, 1.0 + EPSDELTA);
                break;
            case 3:
                POINT_ASSIGN_INTS(nc, 15, 1);
                POINT_ASSIGN_REALS(lo, -1.0,-1.0 - EPSDELTA);
                POINT_ASSIGN_REALS(hi, 2.0 + EPSDELTA, -1.0);
                break;
        }
        hig_cell *sbc = hig_create_root(lo, hi);
        hig_refine_uniform(sbc, nc);
        
        it = higcit_create_all_leaves(sbc);
        mp_assign_from_celliterator(bcm, it, 0); /*assign an id for the cells in the iterator*/
        higcit_destroy(it);

        it = higcit_create_all_leaves(sbc);
        sb[i] = sb_create(sbc, DIRICHLET, bcm); /*creates a bc*/
        for(; !higcit_isfinished(it); higcit_nextcell(it)) {
            hig_cell *bcell = higcit_getcell(it);
			Point bccenter;
			hig_get_center(bcell, bccenter);
			int bcgid = mp_lookup(bcm, hig_get_cid(bcell));
			sb_set_value(sb[i], bcgid, 0.0); /*sets the value of the bc to 0*/
        }
        
        char buffer[100];
        snprintf(buffer, sizeof(buffer), "VTKS/tut05-sb%d.vtk", i);
        FILE *fd = fopen(buffer, "w");
        higio_print_in_vtk2d(fd, sbc); /*more visualization*/
        
        mp_destroy(bcm);
    }
    
    mp_mapper *mp = mp_create();
    /* We create the domain that contains the hig-trees and the boundary conditions*/
    sim_domain *sd = sd_create(mp);
    sd_add_higtree(sd, root1); /*Adds hig-tree*/
    sd_add_higtree(sd, root2); /*Adds hig-tree*/
   
    for(int i =0; i < 4; i++) {
        sd_add_boundary(sd, sb[i]); /*Adds all boundaries*/
    }
    
    /*Assigning an id for every cell in the domain*/
    int last_id = mp_assign_from_celliterator(mp, it, 0);
    higcit_destroy(it);
    it = higcit_create_all_leaves(root2);
    last_id = mp_assign_from_celliterator(mp, it, last_id + 1);
    
}