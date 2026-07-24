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
#include "wls.h"
#include "solver.h"

#define DEBUG
#include "Debug-c.h"


/**
 * @file Stokes equation
*/

real u(real x, real y) {
    return x*x - y*y;
}

real v(real x, real y) {
    return -2*x*y;
}

real p(real x, real y) {
    return 0.0;
}


void set_dirichlet_boundary_conditions(sim_domain *sd, Point lo, Point hi, int divide, real (*func)(real,real)) {
    higcit_celliterator *it;
    sim_boundary *sb[4];
    int nc[2];

    for(int i = 0; i<4; i++){ //iterates in all 4 directions that need a boundary condition
        /* Every instance of sim_boundary (bc) needs a mapper*/
        mp_mapper *bcm = mp_create();
        switch (i) {
            case 0: // left boundary (x = -1.0)
                POINT_ASSIGN_INTS(nc, 1, divide);
                POINT_ASSIGN_REALS(lo, -1.0 - EPSDELTA, -1.0);
                POINT_ASSIGN_REALS(hi, -1.0 + EPSDELTA,  1.0);
                break;
            case 1: // top boundary (y = 1.0)
                POINT_ASSIGN_INTS(nc, divide, 1);
                POINT_ASSIGN_REALS(lo, -1.0, 1.0 - EPSDELTA);
                POINT_ASSIGN_REALS(hi,  1.0, 1.0 + EPSDELTA);
                break;
            case 2: // right boundary (x = 2.0)
                POINT_ASSIGN_INTS(nc, 1, divide);
                POINT_ASSIGN_REALS(lo,  1.0 - EPSDELTA, -1.0);
                POINT_ASSIGN_REALS(hi,  1.0 + EPSDELTA,  1.0);
                break;
            case 3: // bottom boundary (y = -1.0)
                POINT_ASSIGN_INTS(nc, divide, 1);
                POINT_ASSIGN_REALS(lo, -1.0, -1.0 - EPSDELTA);
                POINT_ASSIGN_REALS(hi,  1.0, -1.0 + EPSDELTA);
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
			sb_set_value(sb[i], bcgid, func(bccenter[0], bccenter[1])); /*sets the value of the bc to the actual solution, as we are using Dirichlet BC*/
        }
    }

    for(int i =0; i < 4; i++) {
        sd_add_boundary(sd, sb[i]); /*Adds all boundaries*/
    }
}


int main (int argc, char *argv[]) {
    higtree_initialize(&argc, &argv);
    real ni = 1.0;

    int nc[2] = {100, 100};

    Point lo, hi;

    /* Creates the hig-tree*/
    POINT_ASSIGN_SCALAR(lo, -1.0);
    POINT_ASSIGN_SCALAR(hi, 1.0);
    hig_cell *root = hig_create_root(lo, hi);
    hig_refine_uniform(root, nc);


    

    mp_mapper *mp = mp_create();
    /* We create the domain that contains the hig-trees and the boundary conditions*/
    sim_domain *sd = sd_create(mp);
    sd_add_higtree(sd, root); /*Adds hig-tree*/
    //set_dirichlet_boundary_conditions(sd, lo, hi, nc[0], p);
    sd_set_interpolator_order(sd, 3);

    higcit_celliterator *it;

    int sizef[2];
    sim_facet_domain *sfd[2];
    for(int dim = 0; dim < 2; dim++) {
        sfd[dim] = sfd_create(NULL, dim);
        sfd_copy_higtrees_from_center_domain(sfd[dim], sd);
        set_dirichlet_boundary_conditions(sfd[dim]->cdom, lo, hi, nc[0], dim == 0? u : v);
        sfd_set_interpolator_order(sfd[dim], 3);
        sfd_adjust_facet_ids(sfd[dim]);

        int dim_of_interest[2] = {1-dim, dim};
        higcit_celliterator *cit = higcit_create_all_leaves(sfd_get_higtree(sfd[dim],0));
        higfit_facetiterator *fit = higfit_create_allfacets(cit, dim_of_interest);
        sizef[dim] = mp_assign_from_facetiterator(sfd_get_domain_mapper(sfd[dim]), fit, dim == 0? 0 : sizef[0]);
        higfit_destroy(fit);
    }

    it = higcit_create_all_leaves(root);
    int size = mp_assign_from_celliterator(mp, it, sizef[1]); /*returns the last_id+1, that is, the number of cells of a hig-tree*/
    higcit_destroy(it);



    sim_stencil *stn = stn_create();

    solver *slv = slv_create(SOLVER_ANY, 0, size);
    slv_set_maxnonzeros(slv, 400);


    for(int dim =0; dim < 2; dim++) {
        int dim_of_interest[2] = {1-dim, dim};
        higcit_celliterator *cit = higcit_create_all_leaves(sfd_get_higtree(sfd[dim],0));
        higfit_facetiterator *fit = higfit_create_allfacets(cit, dim_of_interest);
        for(; !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            Point fdelta;
            Point fcenter;

            hig_facet *facet = higfit_getfacet(fit);
            hig_get_facet_delta(facet, fdelta);
            real dx = fdelta[0];
            real dy = fdelta[1];
            hig_get_facet_center(facet, fcenter);
            real cx = fcenter[0];
            real cy = fcenter[1];

            Point ftop      = {cx     , cy + dy};
            Point fbottom   = {cx     , cy - dy};
            Point fright    = {cx + dx, cy     };
            Point fleft     = {cx - dx, cy     };

            Point cleft = {cx, cy};
            Point cright = {cx, cy};
            real d;
            if(dim == 0) {
                cleft[0] = fcenter[0] - dx/2.0;
                cright[0] = fcenter[0] + dx/2.0;
                d = dx;
            }
            else {
                cleft[1] = fcenter[1] - dy/2.0;
                cright[1] = fcenter[1] + dy/2.0;
                d = dy;
            }

            /* Get the Global ID of the facet */
            int cgid = mp_lookup(sfd_get_domain_mapper(sfd[dim]), hig_get_fid(facet));
            
            stn_reset(stn);

            stn_set_rhs(stn, 0.0);  // right side
            
            sfd_get_stencil(sfd[dim], fcenter, fcenter, -ni*(-2.0/(dx*dx) -2.0/(dy*dy)), stn);
            sfd_get_stencil(sfd[dim], fcenter, ftop,    -ni/(dy*dy),                     stn);
            sfd_get_stencil(sfd[dim], fcenter, fbottom, -ni/(dy*dy),                     stn);
            sfd_get_stencil(sfd[dim], fcenter, fright,  -ni/(dx*dx),                     stn);
            sfd_get_stencil(sfd[dim], fcenter, fleft,   -ni/(dx*dx),                     stn);
            sd_get_stencil(sd, fcenter, cleft, -1.0/d, stn);
            sd_get_stencil(sd, fcenter, cright, 1.0/d, stn);

            /* Adds line in the system matrix */
            int *ids = stn_get_ids(stn);
            real *vals = stn_get_vals(stn);
            int numelems = stn_get_numelems(stn);

            slv_set_Ai(slv, cgid, numelems, ids, vals);
            slv_set_bi(slv, cgid, stn_get_rhs(stn));
        }
        higfit_destroy(fit);
    }

    int first = 1;
    it = higcit_create_all_leaves(root);
    for(; !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        Point ccenter, cdelta;
        hig_get_center(c, ccenter);
        real cx = ccenter[0];
        real cy = ccenter[1];
        hig_get_delta(c, cdelta);
        real dx = cdelta[0];
        real dy = cdelta[1];
        
        Point ctop      = {cx         , cy + dy/2.0};
        Point cbottom   = {cx         , cy - dy/2.0};
        Point cright    = {cx + dx/2.0, cy         };
        Point cleft     = {cx - dx/2.0, cy         };


        stn_reset(stn);
        
        int cgid = mp_lookup(mp, hig_get_cid(c));
     
        stn_set_rhs(stn, 0.0);
        sd_get_stencil(sd, ccenter, ccenter, 1.0e-8, stn);
        sfd_get_stencil(sfd[0], ccenter, cright,   1.0/dx, stn);
        sfd_get_stencil(sfd[0], ccenter, cleft,   -1.0/dx, stn);
        sfd_get_stencil(sfd[1], ccenter, ctop,     1.0/dy, stn);
        sfd_get_stencil(sfd[1], ccenter, cbottom, -1.0/dy, stn);
        int *ids = stn_get_ids(stn);
        real *vals = stn_get_vals(stn);
        int numelems = stn_get_numelems(stn);

        slv_set_Ai(slv, cgid, numelems, ids, vals);
        slv_set_bi(slv, cgid, stn_get_rhs(stn));
    
    }
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

    // for(int i=0; i < sd_get_num_higtrees(sd); i++) {
    //     for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
    //         hig_cell *c = higcit_getcell(it);
    //         int cgid = mp_lookup(mp, hig_get_cid(c));
    //         Point ccenter;
    //         hig_get_center(c, ccenter);
    //         real cx = ccenter[0];
    //         real cy = ccenter[1];

    //         /* Absolute error in point */
    //         real error = fabs(p(cx, cy) - slv_get_xi(slv, cgid));
    //         norm_2 += error*error;
    //         norm_inf = ((norm_inf<error) ? error : norm_inf);
    //     }

    //     higcit_destroy(it);
    // }

    for(int dim =0; dim < 2; dim++) {
        int dim_of_interest[2] = {1-dim, dim};
        higcit_celliterator *cit = higcit_create_all_leaves(sfd_get_higtree(sfd[dim], 0));
        higfit_facetiterator *fit = higfit_create_allfacets(cit, dim_of_interest);
        for(; !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            hig_facet *f = higfit_getfacet(fit);
            mp_mapper *fmp = sfd_get_domain_mapper(sfd[dim]);
            int cgid = mp_lookup(fmp, hig_get_fid(f));
            Point ccenter;
            hig_get_facet_center(f, ccenter);
            real cx = ccenter[0];
            real cy = ccenter[1];

            real error;
            /* Absolute error in point */
            if(dim == 0) {
                error = fabs(u(cx, cy) - slv_get_xi(slv, cgid));
            }
            else {
                error = fabs(v(cx, cy) - slv_get_xi(slv, cgid));
            }
            
            norm_2 += error*error;
            norm_inf = ((norm_inf<error) ? error : norm_inf);
        }
        higfit_destroy(fit);
    }

    norm_2 /= (real) size;
    norm_2 = sqrt(norm_2);
    /* Show results */
    printf("norm_inf = %f\n", norm_inf);
    printf("norm_2 = %f\n", norm_2);
}