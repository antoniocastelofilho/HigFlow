#include<stdio.h>
#include<stdlib.h>

#include "domain.h"
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

int func = 1;

#define PI 3.1415926535897932385
real u(Point p) { /*Analytical solution*/
	real x = p[0];
	real y = p[1];
	switch(func) {
		case 0:	return sin(2.0*PI*x)*sin(2.0*PI*y);
		case 1:	return sin(x)*cos(y);
		case 2:	return 5.0;
	}
	printf("undefined function! Should be in 0..4\n");
	exit(0);
}

real g(Point p) { /*forcing function*/
	real x = p[0];
	real y = p[1];
	switch(func) {
		case 0: return -8.0*PI*PI*sin(2.0*PI*x)*sin(2.0*PI*y);
		case 1: return -2.0*sin(x)*cos(y);
		case 2: return 0.0;
	}
}

/* Same as the other tutorials*/
void set_dirichlet_boundary_conditions(sim_domain *sd) {
	higcit_celliterator *it;
	int numbcs = sd_get_num_bcs(sd, DIRICHLET);
	for (int i = 0; i < numbcs; i++) {
		sim_boundary *bc = sd_get_bc(sd, DIRICHLET, i);
		mp_mapper *bm = sb_get_mapper(bc);
		for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
			hig_cell *bcell = higcit_getcell(it);
			Point bccenter;
			hig_get_center(bcell, bccenter);
			int bcgid = mp_lookup(bm, hig_get_cid(bcell));
			sb_set_value(bc, bcgid, u(bccenter));
		}
		higcit_destroy(it);
	}
}

/*writes the VTK file for each rank and its BC and domain*/
void writevtk(psim_domain *psd, int myrank) {
	higcit_celliterator *it;

	sim_domain *local_domain = psd_get_local_domain(psd);
	char filename[1024];
	sprintf(filename, "./VTKS/tut08-sbd-%d.vtk", myrank);
	FILE *fd = fopen(filename, "w");
	it = sd_get_domain_celliterator(local_domain);
	higio_print_celliterator_in_vtk(fd, it);
	fclose(fd);
	higcit_destroy(it);

	sprintf(filename, "./VTKS/tut08-bc-%d.vtk", myrank);
	fd = fopen(filename, "w");
	it = sd_get_bcs_celliterator(local_domain);
	higio_print_celliterator_in_vtk(fd, it);
	fclose(fd);
	higcit_destroy(it);
}

/*Prints the errors*/
void inspect_norm_from_solver_result(solver *s, sim_domain *domain, mp_mapper *global_mapper, int myrank)
{
	const size_t domainsize = slv_get_local_size(s);
	real fval[domainsize];
    real fval2[domainsize];
	higcit_celliterator *it;
    Point deltas[domainsize];
	// computing error

	for(it = sd_get_domain_celliterator(domain); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *cell = higcit_getcell(it);
		Point center;
		int fgid = mp_lookup(global_mapper, hig_get_cid(cell));
		hig_get_center(cell, center);
		fval[fgid] = u(center);
        hig_get_delta(cell, deltas[fgid]);
	}
	higcit_destroy(it);

	{
		slv_get_x(s, fval2);

		for(size_t i = 0; i < domainsize; ++i) {
			fval[i] -= fval2[i];
		}
	}
	real norm_inf, norm_1, norm_2;

	norm_1 = 0.0;
	norm_2 = 0.0;
	norm_inf = 0.0;
	for(size_t i = 0; i < domainsize; ++i)
	{
		real abs = fabs(fval[i]);
		if(abs > norm_inf)
			norm_inf = abs;
		norm_1 += abs * deltas[i][0] * deltas[i][1];
		norm_2 += abs*abs * deltas[i][0] * deltas[i][1];
	}
	norm_2 = sqrt(norm_2);

    printf("rank: %d \tnorm_1: %g \tnorm_2: %g \tnorm_inf: %g \n", myrank, norm_1, norm_2, norm_inf);
}


int main(int argc, char *argv[]) {
	int size;
	DEBUG_DIFF_TIME;
	int myrank;
	int ntasks;

	higtree_initialize(&argc, &argv);

    /*We get the local rank and number of processor used by the MPI*/
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    /*The partition graph, whose function is to store the division of the higtrees of the domain through the different processors*/
	partition_graph *pg = pg_create(MPI_COMM_WORLD);
	/* our domain*/
    psim_domain *psd;

    /* We load the mesh from an .amr file*/
	FILE *fd = fopen(argv[1], "r");
	hig_cell *root = higio_read_from_amr(fd);
	fclose(fd);

    /*Create the domain to be partioned. Every processor gets the full BC in this tutorial*/
	sim_domain *sd = sd_create(NULL);
	sd_create_boundary(root, sd, DIRICHLET);
	set_dirichlet_boundary_conditions(sd);

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

	writevtk(psd, myrank); /*save the VTKS*/

	higcit_celliterator *it;
	sim_stencil *stn = stn_create();

	DEBUG_DIFF_TIME;
	DEBUG_INSPECT("generating matrix...", %s);
	sim_domain *local_domain = psd_get_local_domain(psd);
	DEBUG_INSPECT(localdomainsize, %d);
	solver *slv = slv_create(SOLVER_ANY, psd_get_first_id(psd), localdomainsize);
	slv_set_maxnonzeros(slv, 300);

	mp_mapper *m = sd_get_domain_mapper(local_domain);
	sd_set_interpolator_order(local_domain, 3);
	for(it = sd_get_domain_celliterator(local_domain); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		Point ccenter;
		Point cdelta;
		hig_get_center(c, ccenter);
		real cx = ccenter[0];
		real cy = ccenter[1];
		hig_get_delta(c, cdelta);
		real dx = cdelta[0];
		real dy = cdelta[1];

		Point center = {cx     , cy     };
		Point left   = {cx - dx, cy     };
		Point right  = {cx + dx, cy     };
		Point top    = {cx,      cy + dy};
		Point bottom = {cx,      cy - dy};


		stn_reset(stn);
		stn_set_rhs(stn, g(ccenter));

		sd_get_stencil(local_domain, center, center, (-2.0/(dx*dx)-2.0/(dy*dy)), stn);

		sd_get_stencil(local_domain, center, left,    (1.0/(dx*dx)),              stn);
		sd_get_stencil(local_domain, center, right,   (1.0/(dx*dx)),              stn);
		sd_get_stencil(local_domain, center, top,     (1.0/(dy*dy)),              stn);
		sd_get_stencil(local_domain, center, bottom,  (1.0/(dy*dy)),              stn);

		int *ids = psd_stn_get_gids(psd, stn);
		real *vals = stn_get_vals(stn);
		int numelems = stn_get_numelems(stn);

		int cgid = psd_get_global_id(psd, c);
		slv_set_Ai(slv, cgid, numelems, ids, vals);
		slv_set_bi(slv, cgid, stn_get_rhs(stn));
	}
	higcit_destroy(it);
	stn_destroy(stn);
	DEBUG_DIFF_TIME;
	DEBUG_INSPECT("assembling...", %s);
	slv_assemble(slv);
	DEBUG_DIFF_TIME;

	DEBUG_INSPECT("solving...", %s);
	slv_solve(slv);
	//MatView(slv->A, PETSC_VIEWER_STDOUT_WORLD);

	inspect_norm_from_solver_result(slv, local_domain, m, myrank);

	pg_destroy(pg);
	sd_destroy(local_domain);
	psd_destroy(psd);
	slv_destroy(slv);
	return 0;
}
