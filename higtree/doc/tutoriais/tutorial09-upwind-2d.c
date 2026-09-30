#include <math.h>
#include<stdio.h>
#include<stdlib.h>

#include "coord.h"
#include "domain.h"
#include "higtree-iterator-internal.h"
#include "higtree-iterator.h"
#include "utils.h"
#include "higtree.h"
#include "higtree-io.h"
#include "mapper.h"
#include "wls.h"
#include "glib.h"
#include "pdomain.h"
#include "solver.h"
#include "lbal.h"

#define DEBUG
#include "Debug-c.h"

/*u_t + (v1,v2)*grad(u)=0
 *u(x,y,0) = f(x,y)
 *v1, v2 > 0
 *Using Upwind to solve it
 */

#define PI 3.1415926535897932385
/*Function of the initial condition*/
real f(real x, real y) {
    if(x <= 1.0 && x >= -1.0 &&
       y <= 1.0 && y >=- 1.0) {
        return cos(PI*x/2)*cos(PI*y/2);
    }
    else return 0.0;
}

/*Analytical solution of the PDE*/
real u(real x, real y, real t, real v1, real v2) {
    return f(x-v1*t, y-v2*t);
}

/*For this problem, we only need 2 BCs, at the bottom and the left*/
void set_dirichlet_boundary_conditions(sim_domain *sd, hig_cell *root) {
    Point lo;
    hig_get_lowpoint(root, lo);
    sim_domain *temp = sd_create(NULL); /*We create a temporary domain*/
    sd_create_boundary(root, temp, DIRICHLET); /*And create the boundary conditions for the tree*/
    for(int i =0; i < 4; i++) {
        sim_boundary *sb = sd_get_bc(temp, DIRICHLET, i); /*We get the BCs*/
        hig_cell *sbc = sb_get_higtree(sb);
        hig_cell *test = hig_get_cell_with_point(sbc, lo);
        if(test != NULL) { /*Testing if the low point is in the BC. If it is, we know it's one of the BCs we wanted*/
            sd_add_boundary(sd, sb); /*Add this BC to the domain we want*/
            higcit_celliterator *it;
            mp_mapper *bm = sb_get_mapper(sb);
            for(it= sb_get_celliterator(sb); !higcit_isfinished(it); higcit_nextcell(it)) {
                hig_cell *bcell = higcit_getcell(it);
                Point bccenter;
                hig_get_center(bcell, bccenter);
                int bcgid = mp_lookup(bm, hig_get_cid(bcell));
                sb_set_value(sb, bcgid, 0.0); /*setting the value for the cells of the BC*/
            }
            higcit_destroy(it);
        }
        else sb_destroy(sb); /*If it is a top/right BC, we just destroy it*/
    }
    /*free the temp domain*/
    mp_destroy(temp->m);
    ptm_destroy(temp->cwls);
    free(temp);
}

/*Get the smallest delta in each direction*/
void get_smallest_delta(hig_cell *root, Point delta) {
    higcit_celliterator *it;
    Point delta_smallest = {1000.0, 1000.0};
    for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        Point delta_tmp;
        hig_get_delta(c, delta_tmp);
        for(int i = 0; i <2; i ++) {
            if(delta_smallest[i] > delta_tmp[i]) delta_smallest[i] = delta_tmp[i];
        }
    }
    higcit_destroy(it);
    delta[0] = delta_smallest[0];
    delta[1] = delta_smallest[1];
}

/*Prints the errors of the numerical solution in relation to the analytical solution at each time*/
void print_errors(sim_domain *sd, distributed_property *dp, real v1, real v2, real time, int myrank) {
    higcit_celliterator *it;
    real norm_1 = 0.0;
    real norm_2 = 0.0;
    real norm_inf = -1.0;
    mp_mapper *m = sd_get_domain_mapper(sd);
    for(it = sd_get_domain_celliterator(sd); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        Point ccenter, delta;
        hig_get_center(c, ccenter);
        hig_get_delta(c, delta);
        int cgid = mp_lookup(m, hig_get_cid(c));
        real error = fabs(u(ccenter[0], ccenter[1], time, v1, v2) - dp_get_value(dp, cgid));
        norm_1 += error * delta[0] * delta[1];
        norm_2 += error * error * delta[0] * delta[1];
        if(error > norm_inf) norm_inf = error;
    }
    higcit_destroy(it);
    norm_2 = sqrt(norm_2);

    printf("rank: %d \t time: %f \tnorm_1: %f \t norm_2: %f \tnorm_inf: %f\n", myrank, time, norm_1, norm_2, norm_inf);
}

int main(int argc, char *argv[]) {
	int size;
	DEBUG_DIFF_TIME;
	int myrank;
	int ntasks;
    int total_iterations= 100;
    real v1 = 1.0; //Parameter of the PDE. Must be greater than 0
    real v2 = 1.0; //Parameter of the PDE. Must be greater than 0

	higtree_initialize(&argc, &argv);

    /*We get the local rank and number of processor used by the MPI*/
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    /*The partition graph, whose function is to store the division of the higtrees of the domain through the different processors*/
	partition_graph *pg = pg_create(MPI_COMM_WORLD);
    pg_set_fringe_size(pg, 4);
	/* our domain*/
    psim_domain *psd;

    /* We load the mesh from an .amr file*/
	FILE *fd = fopen(argv[1], "r");
	hig_cell *root = higio_read_from_amr(fd);
	fclose(fd);
    Point delta;
    get_smallest_delta(root, delta);


    real dt = 0.9*1/(v1/delta[0] + v2/delta[1]); //dt satisfying CFL conditions

    DEBUG_INSPECT(delta[0], %f);
    DEBUG_INSPECT(delta[1], %f);
    DEBUG_INSPECT(dt, %f);
    /*Create the domain to be partioned. Every processor gets the full BC in this tutorial*/
	sim_domain *sd = sd_create(NULL);
    set_dirichlet_boundary_conditions(sd, root);
    sd_set_interpolator_order(sd, 1); /*If it's bigger than 1, the method becomes unstable for some reason*/
    /*With order one, the error keeps growing, but at a smaller rate. I recommend to use this method only with regular meshes*/
    

    /*Load balancer makes the division of the higtrees */
	load_balancer *lb = lb_create(MPI_COMM_WORLD, 1);


	if (myrank == 0) {
		lb_add_input_tree(lb, root, true, 0);
	} else {
		hig_destroy(root);
	}
    /*Function that actually make the partitioned higtrees*/
	lb_calc_partition(lb, pg);

	{
		unsigned numtrees = lb_get_num_local_trees(lb);
		for(unsigned i = 0; i < numtrees; ++i) {
			hig_cell *tree = lb_get_local_tree(lb, i, NULL);
			sd_add_higtree(sd, tree); /*Gets the partitioned higtree*/
		}
	}
	lb_destroy(lb);

	psd = psd_create(sd, pg);

	DEBUG_DIFF_TIME;
	psd_synced_mapper(psd); /*Syncs the mapper. All the work that had to be done by hand, is done here*/
	int localdomainsize = psd_get_local_domain_size(psd);

    /*A distributed property assigns a value for each cell. This makes possible using previous values to calculate new ones at another time*/
    distributed_property *dpu_last = psd_create_property(psd);
    distributed_property *dpu_now = psd_create_property(psd);

    sim_domain *localdomain = psd_get_local_domain(psd);
	higcit_celliterator *it;
	sim_stencil *stn = stn_create();
    real time = dt;
    mp_mapper *m = sd_get_domain_mapper(localdomain);

    /*Setting the initial conditions to the dp*/
    for(it = sd_get_domain_celliterator(localdomain); !higcit_isfinished(it); higcit_nextcell(it)) {
        Point ccenter;
        hig_cell *c = higcit_getcell(it);
        hig_get_center(c, ccenter);
        int cgid = mp_lookup(m, hig_get_cid(c));
        dp_set_value(dpu_last, cgid, f(ccenter[0], ccenter[1]));
    }

    higcit_destroy(it);

    /*We iterate in time*/
    for(int t = 0; t < total_iterations; t++) {
        dp_sync(dpu_last);
        /*for each step of time, we iterate the whole domain*/
        for(it = sd_get_domain_celliterator(localdomain); !higcit_isfinished(it); higcit_nextcell(it)) {
            hig_cell *c = higcit_getcell(it);
            Point ccenter, cdelta;
            hig_get_center(c, ccenter);
            hig_get_delta(c, cdelta);
            real dx = cdelta[0];
            real dy = cdelta[1];
            // Point top = {ccenter[0], ccenter[1] + dy};
            Point bottom = {ccenter[0], ccenter[1] - dy};
            // Point right = {ccenter[0] + dx, ccenter[1]};
            Point left = {ccenter[0] - dx, ccenter[1]};

            int cgid = mp_lookup(m, hig_get_cid(c));
            // real uright = sd_dp_interpolate(localdomain, dpu, ccenter, right, stn);
            real uleft = sd_dp_interpolate(localdomain, dpu_last, ccenter, left, stn);
            // real utop = sd_dp_interpolate(localdomain, dpu, ccenter, top, stn);
            real ubottom = sd_dp_interpolate(localdomain, dpu_last, ccenter, bottom, stn);
            real ucenter = dp_get_value(dpu_last, cgid);

            /*Plug the values in the discretization of the PDE*/
            real next_u = ucenter -((v1*dt)/dx)*(ucenter-uleft) -((v2*dt)/dy)*(ucenter-ubottom);

            
            dp_set_value(dpu_now, cgid, next_u);
        }
        // dp_sync(dpu_now);
        print_errors(localdomain, dpu_now, v1, v2, time, myrank);
        dp_copy_values(dpu_last, dpu_now); /*refresh the dp's to be used in the next iteration*/
        time += dt;
        higcit_destroy(it);
    }

}