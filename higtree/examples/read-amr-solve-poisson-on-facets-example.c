#include<stdio.h>
#include<stdlib.h>

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
#define PI 3.1415926535897932385

int func = 1;
real f(Point p) {
	real x = p[0];
	real y = p[1];
	switch(func) {
		case 0:	return sin(x)*cos(y);
		case 1:	return 5.0;
	}
	printf("undefined function! Should be in 0..4\n");
	exit(0);
}

real g(Point p) {
	real x = p[0];
	real y = p[1];
	switch(func) {
		case 0: return -2.0*sin(x)*cos(y);
		case 1: return 0.0;
	}
}

real dfdy(real x, real y) {
	switch(func) {
		case 0: return -sin(x)*sin(y);
	}
}

real dfdx(real x, real y) {
	switch(func) {
		case 0: return cos(x)*cos(y);
	}
}

const char * f_str() {
	switch(func) {
		case 0:	return "sin(x)*cos(y)"
#if DIM == 3
			"*sin(z)"
#endif
				;
	}
	return "undefined function!";
}

void set_dirichlet_boundary_conditions(sim_facet_domain *sfd) {
	higcit_celliterator *it;
	int numbcs = sfd_get_num_bcs(sfd, DIRICHLET);
	for (int i = 0; i < numbcs; i++) {
		sim_boundary *bc = sfd_get_bc(sfd, DIRICHLET, i);
		mp_mapper *bm = sb_get_mapper(bc);
		for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
			hig_cell *bcell = higcit_getcell(it);
			Point bccenter;
			hig_get_center(bcell, bccenter);
			int bcgid = mp_lookup(bm, hig_get_cid(bcell));
			sb_set_value(bc, bcgid, f(bccenter));
		}
		higcit_destroy(it);
	}
}

//void writevtk(psim_domain *psfd, int myrank) {
//	higcit_celliterator *it;
//	if (myrank == 0) {
//		char filename[1024];
//		sprintf(filename, "/tmp/d.vtk", myrank);
//		FILE *fd = fopen(filename, "w");
//		sim_domain *global_domain = psd_get_global_domain(psfd);
//		it = sd_get_domain_celliterator(global_domain);
//		higio_print_celliterator_in_vtk2d(fd, it);
//		fclose(fd);
//		higcit_destroy(it);
//	}
//	sim_domain *local_domain = psd_get_local_domain(psfd);
//	char filename[1024];
//	sprintf(filename, "/tmp/sbd-%d.vtk", myrank);
//	FILE *fd = fopen(filename, "w");
//	it = sd_get_domain_celliterator(local_domain);
//	higio_print_celliterator_in_vtk2d(fd, it);
//	fclose(fd);
//	higcit_destroy(it);
//
//	sprintf(filename, "/tmp/sbfd-%d.vtk", myrank);
//	fd = fopen(filename, "w");
//	it = higcit_create_all_leaves(sd_get_higtree(local_domain, 0));
//	higio_print_celliterator_in_vtk2d(fd, it);
//	fclose(fd);
//	higcit_destroy(it);
//
//	sprintf(filename, "/tmp/bc-%d.vtk", myrank);
//	fd = fopen(filename, "w");
//	it = sd_get_bcs_celliterator(local_domain);
//	higio_print_celliterator_in_vtk2d(fd, it);
//	fclose(fd);
//	higcit_destroy(it);
//}

void __print_costs(int numcellsperblock[], int numblocks) {
	printf("\n");
	for(int c = 0; c < numblocks; c++) {
		printf("%d = %d\n", c, numcellsperblock[c]);
	}
}

int __calc_num_points_aux(int order) {
        if (DIM == 2) {
                if (order == 3) {
                        return 40;
                } else if (order == 2) {
                        return 25;
                }
        } else if (DIM == 3) {
                if (order == 3) {
                        return 120;
                } else if (order == 2) {
                        return 45;
                }
        }
        return 2 * DIM * wls_num_min_points(order);
}

wls_interpolator * __sd_create_wls_interpolator_aux(sim_domain *d, int maxpts) {
	if (d->inter.wls == NULL) {
		d->inter.wls = wls_create(d->inter.order, maxpts);
	}
	return d->inter.wls;
}

void inspect_norm_from_solver_result(solver *s, sim_facet_domain *domain, mp_mapper *global_mapper)
{
	const size_t domainsize = slv_get_local_size(s);
	real fval[domainsize];
	higfit_facetiterator *fit;

	// computing error

	slv_get_x(s, fval);
	for(fit = sfd_get_domain_facetiterator(domain); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
		hig_facet *fac = higfit_getfacet(fit);
		Point fcenter;
		int fgid = mp_lookup(global_mapper, hig_get_fid(fac));
		hig_get_facet_center(fac, fcenter);
		fval[fgid] = f(fcenter);
	}
	higfit_destroy(fit);

	{
		real fval2[domainsize];
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
		norm_1 += abs;
		norm_2 += abs*abs;
	}
	norm_2 = sqrt(norm_2);

	DEBUG_INSPECT(norm_2, %g);
	DEBUG_INSPECT(norm_inf, %g);
}

int main(int argc, char *argv[]) {
	assert(DIM == 2 || DIM == 3);

	higtree_initialize(&argc, &argv);

	int myrank;
	int ntasks;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

	partition_graph *pg = pg_create(MPI_COMM_WORLD);
	pg_set_fringe_size(pg, 2);
	psim_facet_domain *psfd;

	DEBUG_DIFF_TIME;
	FILE *fd = fopen(argv[1], "r");
	higio_amr_info *mi = higio_read_amr_info(fd);
	fclose(fd);
	load_balancer *lb = higio_partition_from_amr_info(pg, mi);

	higcit_celliterator *it;
	mp_mapper *m = mp_create();
	int dim = DIM == 2 ? 0 : 1;
	sim_facet_domain *sfd = sfd_create(m, dim);

	{
		unsigned numtrees = lb_get_num_local_trees(lb);
		for(unsigned i = 0; i < numtrees; ++i) {
			hig_cell *root = lb_get_local_tree(lb, i, NULL);
			sfd_add_higtree(sfd, root);
		}
	}
	lb_destroy(lb);

	sfd_adjust_facet_ids(sfd);
	DEBUG_DIFF_TIME;
//
	for(int dim = 0; dim < DIM; dim++) {
		for(int dir = 0; dir < 2; dir++) {
			hig_cell *bcg = higio_read_bc_from_amr(mi, dim, dir);
			mp_mapper *bm = mp_create();
			it = higcit_create_all_leaves(bcg);
			mp_assign_from_celliterator(bm, it, 0);
			higcit_destroy(it);
			sim_boundary *bc = sb_create(bcg, DIRICHLET, bm);
			sfd_add_boundary(sfd, bc);
		}
	}
	DEBUG_DIFF_TIME;
	set_dirichlet_boundary_conditions(sfd);
	DEBUG_DIFF_TIME;
	psfd = psfd_create_from_pg(sfd, pg);
	DEBUG_DIFF_TIME;
        psfd_compute_sfbi(psfd);
	psfd_synced_mapper(psfd);
	DEBUG_DIFF_TIME;
//
//	//writevtk(psfd, myrank);
//	psd_destroy_global_domain(psfd);
//
	sim_stencil *stn = stn_create();

	DEBUG_DIFF_TIME;
	DEBUG_INSPECT("generating matrix...", %s);
	sim_facet_domain *local_domain = psfd_get_local_domain(psfd);
	higfit_facetiterator *fit;
	m = sfd_get_domain_mapper(local_domain);
	psfd_synced_mapper(psfd);
	int localdomainsize = psfd_get_local_domain_size(psfd);
	DEBUG_INSPECT(localdomainsize, %d);
	solver *slv = slv_create(SOLVER_ANY, psfd_get_first_id(psfd), localdomainsize);
	slv_set_maxnonzeros(slv, 100);

	sfd_set_interpolator_order(local_domain, 2);
	for(fit = sfd_get_domain_facetiterator(local_domain); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
		hig_facet *f = higfit_getfacet(fit);
		Point fcenter;
		Point fdelta;
		hig_get_facet_center(f, fcenter);
		hig_get_facet_delta(f, fdelta);

		stn_reset(stn);
		stn_set_rhs(stn, g(fcenter));
		real alpha = 0.0;
		for(int d = 0; d < DIM; d++) {
			Point p;
			real weight = (1.0/(fdelta[d]*fdelta[d]));
			POINT_ASSIGN(p, fcenter);
			p[d] = fcenter[d] + fdelta[d];
			sfd_get_stencil(local_domain, fcenter, p, weight, stn);
			p[d] = fcenter[d] - fdelta[d];
			sfd_get_stencil(local_domain, fcenter, p, weight, stn);
			alpha -= 2.0 * weight;
		}
		sfd_get_stencil(local_domain, fcenter, fcenter, alpha, stn);

		int *ids = psfd_stn_get_gids(psfd, stn);
		real *vals = stn_get_vals(stn);
		int numelems = stn_get_numelems(stn);

		int fgid = psfd_get_global_id(psfd, f);
		slv_set_Ai(slv, fgid, numelems, ids, vals);
		slv_set_bi(slv, fgid, stn_get_rhs(stn));
	}
	higfit_destroy(fit);
	stn_destroy(stn);
	DEBUG_DIFF_TIME;
	DEBUG_INSPECT("assembling...", %s);
	slv_assemble(slv);
	DEBUG_DIFF_TIME;

	DEBUG_INSPECT("solving...", %s);
	int psize[2];
	slv_solve(slv);
	//MatView(slv->A, PETSC_VIEWER_STDOUT_WORLD);

	if (myrank == 0) {
		inspect_norm_from_solver_result(slv, local_domain, m);
	}

	sfd_destroy(local_domain);
	psfd_destroy(psfd);
	slv_destroy(slv);
	pg_destroy(pg);
}

