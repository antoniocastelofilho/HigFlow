#include <stdio.h>
#include<stdlib.h>
#include<math.h>

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
 * @file We solve a vector Helmholtz equation in the facets of the domain
 *We can write the decoupled equations: -∇²u(x,y) + k²u(x,y) = f₁(x,y)
 *                                      -∇²v(x,y) + k²v(x,y) = f₂(x,y)
 *u is the solution across x-axis and are stored in the facet of dim 0, and v is the solution in y-axis, stored in the facets of dim 1
*/


/* Analytical solution u*/
real u(real x, real y) {
    return cos((M_PI*x)/2)*cos((M_PI*y)/2);
}

/* Forcing term of u (f₁)*/
real fu(real x, real y, real k) {
    return (((M_PI*M_PI)/2 + k*k)*u(x,y));
}

/* Analytical solution v*/
real v(real x, real y) {
    return sin(M_PI*x)*sin(M_PI*y);
}

/* Forcing term of v (f₂)*/
real fv(real x, real y, real k) {
    return (2*M_PI*M_PI + k*k)*v(x,y);
}

/**
 *@brief Given a function, we set the corresponding dirichlet boundary conditions to ONE higtree with 
 *lo low point and hi high point. Then the boundary conditions are added to the domain sd.
 *divide is an array of size DIM and tells how much to refine the BC in each axis.
*/
void set_dirichlet_boundary_conditions(sim_domain *sd, Point lo, Point hi, int *divide, real (*func)(real,real)) {
    higcit_celliterator *it;
    sim_boundary *sb[4];
    int nc[2];

    for(int i = 0; i<4; i++){ //iterates in all 4 directions that need a boundary condition
        /* Every instance of sim_boundary (bc) needs a mapper*/
        mp_mapper *bcm = mp_create();
        switch (i) {
            case 0: // left boundary (x = -1.0)
                POINT_ASSIGN_INTS(nc, 1, divide[1]);
                POINT_ASSIGN_REALS(lo, -1.0 - EPSDELTA, -1.0);
                POINT_ASSIGN_REALS(hi, -1.0 + EPSDELTA,  1.0);
                break;
            case 1: // top boundary (y = 1.0)
                POINT_ASSIGN_INTS(nc, divide[0], 1);
                POINT_ASSIGN_REALS(lo, -1.0, 1.0 - EPSDELTA);
                POINT_ASSIGN_REALS(hi,  1.0, 1.0 + EPSDELTA);
                break;
            case 2: // right boundary (x = 2.0)
                POINT_ASSIGN_INTS(nc, 1, divide[1]);
                POINT_ASSIGN_REALS(lo,  1.0 - EPSDELTA, -1.0);
                POINT_ASSIGN_REALS(hi,  1.0 + EPSDELTA,  1.0);
                break;
            case 3: // bottom boundary (y = -1.0)
                POINT_ASSIGN_INTS(nc, divide[0], 1);
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
    real k = 1.0;

    int nc[2] = {100, 100}; /*How much we are refining in each direction*/

    Point lo, hi;

    /* Creates the hig-tree*/
    POINT_ASSIGN_SCALAR(lo, -1.0);
    POINT_ASSIGN_SCALAR(hi, 1.0);
    hig_cell *root = hig_create_root(lo, hi);
    hig_refine_uniform(root, nc);


    higcit_celliterator *it;

    int size;
    int start = 0;
    /*We create 2 facet domains, one for each axis*/
    sim_facet_domain *sfd[2];
    for(int dim = 0; dim < 2; dim++) {
        sfd[dim] = sfd_create(NULL, dim);
        sfd_add_higtree(sfd[dim], root);
        set_dirichlet_boundary_conditions(sfd[dim]->cdom, lo, hi, nc, dim == 0? u : v);
        sfd_set_interpolator_order(sfd[dim], 3);
        sfd_adjust_facet_ids(sfd[dim]); /*Facets with same center coordinate are dealt with this function*/

        int dim_of_interest[2] = {1-dim, dim}; /*if the entry of dim_of_interest is 1, then it is a dim of interest!*/
        higcit_celliterator *cit = higcit_create_all_leaves(sfd_get_higtree(sfd[dim],0));
        higfit_facetiterator *fit = higfit_create_allfacets(cit, dim_of_interest); /*Iterates through all the facets in this dim*/
        /* Assign the ids to the mapper. We are doing it so all the u unkowns are the first ids and the v unknowns come after.
         * In this way, the solver creates a decoupled block system, because the equations are decoupled. If the equations were coupled,
         * we'd have the ids in the right way to make the coupled system.
         */
        size = mp_assign_from_facetiterator(sfd_get_domain_mapper(sfd[dim]), fit, start);
        start = size;
        higfit_destroy(fit);
    }



    sim_stencil *stn = stn_create();

    solver *slv = slv_create(SOLVER_ANY, 0, size);
    slv_set_maxnonzeros(slv, 300);

    /* We make the stencil for each direction...*/
    for(int dim =0; dim < 2; dim++) {
        /* ...and each facet in this direction*/
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


            /* Get the Global ID of the facet */
            int cgid = mp_lookup(sfd_get_domain_mapper(sfd[dim]), hig_get_fid(facet));
            
            stn_reset(stn);

            stn_set_rhs(stn, dim == 0 ? fu(cx, cy, k) : fv(cx, cy, k));  // right side
            
            sfd_get_stencil(sfd[dim], fcenter, fcenter, 2.0/(dx*dx) + 2.0/(dy*dy) + k*k, stn);
            sfd_get_stencil(sfd[dim], fcenter, ftop,    -1/(dy*dy),                      stn);
            sfd_get_stencil(sfd[dim], fcenter, fbottom, -1/(dy*dy),                      stn);
            sfd_get_stencil(sfd[dim], fcenter, fright,  -1/(dx*dx),                      stn);
            sfd_get_stencil(sfd[dim], fcenter, fleft,   -1/(dx*dx),                      stn);

            /* Adds line in the system matrix */
            int *ids = stn_get_ids(stn);
            real *vals = stn_get_vals(stn);
            int numelems = stn_get_numelems(stn);

            slv_set_Ai(slv, cgid, numelems, ids, vals);
            slv_set_bi(slv, cgid, stn_get_rhs(stn));
        }
        higfit_destroy(fit);
    }


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