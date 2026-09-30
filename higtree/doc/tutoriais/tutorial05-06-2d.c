#include <stdio.h>
#include<stdlib.h>


//FOR TUTORIAL 5
#include "higtree-iterator.h"
#include "higtree.h"
#include "higtree-io.h"
#include "domain.h"
#include "mapper.h"

//FOR TUTORIAL 6
#include "solver.h"

#define DEBUG
#include "Debug-c.h"


/**
 * @file Tutorial 5: creates a domain with 2 HiG-Trees and dirichlet boundary conditions (Neumann bc are similar)
 *       Tutorial 6: Uses the domain created to solve a poisson equation
 * All is done in 2D, but the 3D version is very similar
 */



/* Global variable for the selection of the test function */
int func = 0;

/**
 * @brief Analytical function u(x,y) - exact solution
 * @param x x coordinate
 * @param y y coordinate
 * @return function value at point (x,y)
 */
real u(real x, real y) {
	switch(func) {
		case 0:	return sin(x)*cos(y);              // Smooth solution
		case 1: return 4*x*x*y+2*x*x-4*y*x-5*x;      // Polynomial
		case 2: return cos(x*y)+sin(x+y);       // Trigonometric combination
		case 3: return cos(x-y)*sin(x+y);       // Trigonometric product
		case 4: return cos(x+sin(y));           // Trigonometric composition
		case 5: return x*x+y*y-x*y-1.1234;           // Square form
	}
	printf("undefined function! Should be in 0..5\n");
	exit(0);
}

/**
 * @brief source term g(x,y) = ∇²u (RHS of Poisson equation)
 * @param x x coordinate
 * @param y y coordinate
 * @return Laplacian value at point (x,y)
 */
real g(real x, real y) {
	switch(func) {
		case 0: return -2.0*sin(x)*cos(y);                                                        // Laplacian of sin(x)cos(y)
		case 1: return 4 + 8*y;                                                                     // Laplacian of polinomial
		case 2: return (-x*x-y*y)*cos(x*y)-2*sin(x+y);                                         // Complex Laplacian
		case 3: return -4*cos(x-y)*sin(x+y);                                                   // Laplacian of product
		case 4: return -(1+cos(y)*cos(y))*cos(x+sin(y))+sin(y)*(sin(x+sin(y)));    // Laplacian of composition
		case 5: return 4.0;                                                                         // Constant Laplacian
	}
}


int main (int argc, char *argv[]) {
    higtree_initialize(&argc, &argv);

    /**======================
        TUTORIAL 5 PART
    =======================*/
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
            case 0: // left boundary (x = -1.0)
                POINT_ASSIGN_INTS(nc, 1, 10);
                POINT_ASSIGN_REALS(lo, -1.0 - EPSDELTA, -1.0);
                POINT_ASSIGN_REALS(hi, -1.0 + EPSDELTA,  1.0);
                break;
            case 1: // top boundary (y = 1.0)
                POINT_ASSIGN_INTS(nc, 15, 1);
                POINT_ASSIGN_REALS(lo, -1.0, 1.0 - EPSDELTA);
                POINT_ASSIGN_REALS(hi,  2.0, 1.0 + EPSDELTA);
                break;
            case 2: // right boundary (x = 2.0)
                POINT_ASSIGN_INTS(nc, 1, 5);
                POINT_ASSIGN_REALS(lo,  2.0 - EPSDELTA, -1.0);
                POINT_ASSIGN_REALS(hi,  2.0 + EPSDELTA,  1.0);
                break;
            case 3: // bottom boundary (y = -1.0)
                POINT_ASSIGN_INTS(nc, 15, 1);
                POINT_ASSIGN_REALS(lo, -1.0, -1.0 - EPSDELTA);
                POINT_ASSIGN_REALS(hi,  2.0, -1.0 + EPSDELTA);
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
			sb_set_value(sb[i], bcgid, u(bccenter[0], bccenter[1])); /*sets the value of the bc to the actual solution, as we are using Dirichlet BC*/
        }

        char buffer[100];
        snprintf(buffer, sizeof(buffer), "VTKS/tut05-sb%d.vtk", i);
        FILE *fd = fopen(buffer, "w");
        higio_print_in_vtk2d(fd, sbc); /*more visualization*/
    }
    
    mp_mapper *mp = mp_create();
    /* We create the domain that contains the hig-trees and the boundary conditions*/
    sim_domain *sd = sd_create(mp);
    sd_add_higtree(sd, root1); /*Adds hig-tree*/
    sd_add_higtree(sd, root2); /*Adds hig-tree*/
   
    for(int i =0; i < 4; i++) {
        sd_add_boundary(sd, sb[i]); /*Adds all boundaries*/
    }

    
    /*Assigning an id for every cell in the domain.
    Each cell of a hig-tree already has an id, but it is local to that hig-tree. We need to make the ids global, with the mapper*/
    it = higcit_create_all_leaves(root1);
    int size = mp_assign_from_celliterator(mp, it, 0); /*returns the last_id+1 - first_id, that is, the number of cells of a hig-tree*/
    higcit_destroy(it);

    it = higcit_create_all_leaves(root2);
    size = mp_assign_from_celliterator(mp, it, size); /*we input size because size-1 was the last id used*/
    higcit_destroy(it);
    


    /**======================
        TUTORIAL 6 PART
    =======================*/
    /* Preparing to solve linear system */

    /* A stencil is related to a cell and gets the weights of other cells in relation to the main one.*/
	sim_stencil *stn = stn_create();

	printf("generating matrix...\n");
	DEBUG_DIFF_TIME;
	solver *slv = slv_create(SOLVER_ANY, 0, size);
	sd_set_interpolator_order(sd, 3); /* We set the order of the polynomial that will make the interpolations*/
    
    /* the max number of non-zeros numbers in a line. since we're doind a 5-point stencil poisson, there's a max of 4 interpolations.
    Each interpolation (with order 3) uses a max of 60 points, which gives a real max of 241 non zeros per line*/
    slv_set_maxnonzeros(slv, 300); 
    

	/* Making the system matrix using finite difference */

    for(it = sd_get_domain_celliterator(sd); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        Point ccenter;
        Point cdelta;
        hig_get_center(c, ccenter);
        real cx = ccenter[0];
        real cy = ccenter[1];
        hig_get_delta(c, cdelta);
        real dx = cdelta[0];
        real dy = cdelta[1];
        dy *= 1.0;

        /* Defines the points for a 5-point stencil*/
        Point center = {cx     , cy     };
        Point left   = {cx - dx, cy     };
        Point right  = {cx + dx, cy     };
        Point top    = {cx,      cy + dy};
        Point bottom = {cx,      cy - dy};

        /* Makes the discrete Laplacian stencil */
        stn_reset(stn);
        stn_set_rhs(stn, g(cx, cy));  // right side
        
        /* Coefficients of the finite difference scheme:
        ∇²u ≈ (u_{i-1,j} + u_{i+1,j} - 2u_{i,j})/dx² + (u_{i,j-1} + u_{i,j+1} - 2u_{i,j})/dy² */
        sd_get_stencil(sd, center, center, (-2.0/(dx*dx)-2.0/(dy*dy)), stn);  // Center
        sd_get_stencil(sd, center, left,    (1.0/(dx*dx)),              stn);  // Left neighbor
        sd_get_stencil(sd, center, right,   (1.0/(dx*dx)),              stn);  // Right neighbor
        sd_get_stencil(sd, center, top,     (1.0/(dy*dy)),              stn);  // Top neighbor
        sd_get_stencil(sd, center, bottom,  (1.0/(dy*dy)),              stn);  // Bottom neighbor

        /* Adds line in the system matrix */
        int *ids = stn_get_ids(stn);
        real *vals = stn_get_vals(stn);
        int numelems = stn_get_numelems(stn);

        int cgid = mp_lookup(mp, hig_get_cid(c));
        slv_set_Ai(slv, cgid, numelems, ids, vals);
        slv_set_bi(slv, cgid, stn_get_rhs(stn));
    }
    printf("done\n");
    DEBUG_DIFF_TIME;
	higcit_destroy(it);
	
	/* Mount and solve linear system */
	printf("assembling...\n");
	slv_assemble(slv);
	DEBUG_DIFF_TIME;

	printf("solving...\n");
	slv_solve(slv);
	DEBUG_DIFF_TIME;

	/* Calculating errors */
	real norm_inf = -1.0;
	real norm_2 = 0.0;

    for(int i=0; i < sd_get_num_higtrees(sd); i++) {
        hig_cell *root = sd_get_higtree(sd, i);
        for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
            hig_cell *c = higcit_getcell(it);
            int cgid = mp_lookup(mp, hig_get_cid(c));
            Point ccenter;
            hig_get_center(c, ccenter);
            real cx = ccenter[0];
            real cy = ccenter[1];

            /* Absolute error in point */
            real error = fabs(u(cx, cy) - slv_get_xi(slv, cgid));
            norm_2 += error*error;
            norm_inf = ((norm_inf<error) ? error : norm_inf);
        }
        norm_2 /= (real) size;
        norm_2 = sqrt(norm_2);

        /* Show results */
        printf("norm_inf = %f\n", norm_inf);
        printf("norm_2 = %f\n", norm_2);
        higcit_destroy(it);
        return 0;
    }
}