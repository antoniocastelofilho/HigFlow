#include "higtree.h"
#include "higtree-io.h"
#include "higtree-iterator.h"
#include "pdomain.h"
#include "lbal.h"
#include "solver.h"
#define DEBUG
#include "Debug-c.h"

DECL_CLOCK(total)
DECL_CLOCK(poisson)
DECL_CLOCK(syncprop)
DECL_CLOCK(firstiter)

sim_boundary *make_bc(hig_cell *bcg, int type) {
	mp_mapper *bm = mp_create();
	higcit_celliterator *it;
	it = higcit_create_all_leaves(bcg);
	mp_assign_from_celliterator(bm, it, 0);
	higcit_destroy(it);
	sim_boundary *bc = sb_create(bcg, type, bm);
	return bc;
}

void computeBCP(psim_domain *psd, higio_amr_info *mi) {
	sim_domain *sd = psd_get_local_domain(psd);
	higcit_celliterator *it;
	for(int dim = 0; dim < DIM; dim++) {
		for(int dir = 0; dir < 2; dir++) {
			hig_cell *bcg = higio_read_bc_from_amr(mi, dim, dir);
			sim_boundary *bc = make_bc(bcg, NEUMANN);
			sd_add_boundary(sd, bc);
		}
	}
	int numbcs = sd_get_num_bcs(sd, NEUMANN);
	for (int i = 0; i < numbcs; i++) {
		sim_boundary *bc = sd_get_bc(sd, NEUMANN, i);
		mp_mapper *bm = sb_get_mapper(bc);
		for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
			hig_cell *bcell = higcit_getcell(it);
			Point bccenter;
			hig_get_center(bcell, bccenter);
			int bcgid = mp_lookup(bm, hig_get_cid(bcell));
			sb_set_value(bc, bcgid, 0.0);
		}
		higcit_destroy(it);
	}
	Point l, h, center;
	hig_cell *c = sd_get_higtree(sd, 0);
	hig_get_center(c, center);
	POINT_SUB_SCALAR(l, center, EPSDELTA);
	POINT_ADD_SCALAR(h, center, EPSDELTA);
	hig_cell *bcg = hig_create_root(l, h);
	sim_boundary *bc = make_bc(bcg, DIRICHLET);
	sd_add_boundary(sd, bc);
	mp_mapper *bm = sb_get_mapper(bc);
	int bcgid = mp_lookup(bm, hig_get_cid(bcg));
	sb_set_value(bc, bcgid, 0.0);
}

void computeBCU(psim_facet_domain *psfdu[DIM], higio_amr_info *mi) {
	higcit_celliterator *it;
	sim_facet_domain *sfdu[DIM];
	int hig = 0;
	sim_boundary *bctop = NULL;
	for(int dim = 0; dim < DIM; dim++) {
		sfdu[dim] = psfd_get_local_domain(psfdu[dim]);
		for(int dim2 = 0; dim2 < DIM; dim2++) {
			for(int dir = 0; dir < 2; dir++) {
				hig_cell *bcg = higio_read_bc_from_amr(mi, dim2, dir);
				sim_boundary *bc = make_bc(bcg, DIRICHLET);
				sfd_add_boundary(sfdu[dim], bc);
				if (dim == 0 && dir == 1 && dim2 == 1) {
					bctop = bc;
				}
			}
		}
		sim_domain *sd = sfdu[dim]->cdom;
		int numbcs = sd_get_num_bcs(sd, DIRICHLET);
		DEBUG_INSPECT(numbcs, %d);
		for (int i = 0; i < numbcs; i++) {
			sim_boundary *bc = sd_get_bc(sd, DIRICHLET, i);
			mp_mapper *bm = sb_get_mapper(bc);
			for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
				hig_cell *bcell = higcit_getcell(it);
				Point bccenter;
				hig_get_center(bcell, bccenter);
				int bcgid = mp_lookup(bm, hig_get_cid(bcell));
				sb_set_value(bc, bcgid, 0.0);
			}
			higcit_destroy(it);
		}
	}
	mp_mapper *bm = sb_get_mapper(bctop);
	Point l = {0.0, 0.0, 0.0};
	Point h = {1.0, 1.0, 1.0};
	Point delta;
        POINT_SUB(delta, h, l);

	real dist = delta[0];
	real max = 1.0;
	real a = 16.0 * max;
	real b = -32.0 * max;
	real c = 16.0 * max;
	for(it = sb_get_celliterator(bctop); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *bcell = higcit_getcell(it);
		Point bccenter;
		hig_get_center(bcell, bccenter);
		int bcgid = mp_lookup(bm, hig_get_cid(bcell));
		real pos = (bccenter[0] - l[0]) / dist;
		real v = (a * pos * pos * pos * pos) + (b * pos * pos * pos) + (c * pos * pos);
		v = 1.0;
		sb_set_value(bctop, bcgid, v);
	}
	higcit_destroy(it);
}

void initializeP(psim_domain *psdp, distributed_property *dpp) {
	higcit_celliterator *it;
	sim_domain *sdp = psd_get_local_domain(psdp);
	mp_mapper *m = sd_get_domain_mapper(sdp);
	for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		int cgid = mp_lookup(m, hig_get_cid(c));
		dp_set_value(dpp, cgid, 0.0);
	}
	higcit_destroy(it);
}

real analyticU(real x, real y) {
	return 0.0;
	return 8.0*(x*x*(1.0+x*(-2.0+x)))*(y*(-2.0+4.0*y*y));
}

real analyticV(real x, real y) {
	return 0.0;
	return -8.0*(x*(2.0+x*(-6.0+4.0*x)))*(y*y*(-1.0+y*y));;
}

void initializeU(psim_facet_domain *psfdu[DIM], distributed_property *dpu[DIM]) {
	higfit_facetiterator *fit;
	sim_facet_domain *sfdu[DIM];
	for(int dim = 0; dim < DIM; dim++) {
		sfdu[dim] = psfd_get_local_domain(psfdu[dim]);
		mp_mapper *m = sfd_get_domain_mapper(sfdu[dim]);
		Point l, h, delta;
		hig_cell *dom = sfd_get_higtree(sfdu[dim], 0);
		hig_get_lowpoint(dom, l);
		hig_get_highpoint(dom, h);
		hig_get_delta(dom, delta);
		for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
			hig_facet *f = higfit_getfacet(fit);
			Point fcenter;
			hig_get_facet_center(f, fcenter);
			Point x;
			POINT_SUB(x, fcenter, l);
			POINT_DIV(x, x, delta);
			int fgid = mp_lookup(m, hig_get_fid(f));
			real val;
			if (dim == 0) {
				val = analyticU(x[0], x[1]);
			} else {
				val = analyticV(x[0], x[1]);
			}
			dp_set_value(dpu[dim], fgid, val);
		}
		higfit_destroy(fit);
	}
}

void compute_facet_dudx_and_du2dx2_at_point(sim_facet_domain *sfdu, Point center, Point delta, int dim, real alpha, distributed_property *dpu, sim_stencil *stn, real valueatpoint, real *dudx, real *du2dx2) {
	Point p;
	POINT_ASSIGN(p, center);
	p[dim] = center[dim] - alpha * delta[dim];
	real valuel = sfd_dp_interpolate(sfdu, dpu, center, p, stn);
	p[dim] = center[dim] + alpha * delta[dim];
	real valueh = sfd_dp_interpolate(sfdu, dpu, center, p, stn);
	if (dudx != NULL) {
		*dudx = (valueh - valuel) / (2.0 * alpha * delta[dim]);
	}
	if (du2dx2 != NULL) {
		*du2dx2 = (valuel - 2.0*valueatpoint + valueh)/(delta[dim] * delta[dim] * alpha * alpha);
	}
}

void compute_dpdx_at_point(sim_domain *sdp, Point center, Point delta, int dim, real alpha, distributed_property *dpp, sim_stencil *stn, real *dpdx) {
	Point p;
	POINT_ASSIGN(p, center);
	p[dim] = center[dim] - alpha * delta[dim];
	real valuel = sd_dp_interpolate(sdp, dpp, center, p, stn);
	p[dim] = center[dim] + alpha * delta[dim];
	real valueh = sd_dp_interpolate(sdp, dpp, center, p, stn);
	*dpdx = (valueh - valuel) / (2.0 * alpha * delta[dim]);
}

void compute_ustar(psim_facet_domain *psfdu[DIM], int dim, distributed_property *dpu[DIM], real dt, real Re, distributed_property *dpustar[DIM], sim_stencil *stn) {
	higfit_facetiterator *fit;
	sim_facet_domain *sfdu[DIM];
	for(int dim = 0; dim < DIM; dim++) {
		sfdu[dim] = psfd_get_local_domain(psfdu[dim]);
	}
	mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
	for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
		hig_facet *f = higfit_getfacet(fit);
		int fgid = mp_lookup(mu, hig_get_fid(f));
		Point fcenter;
		Point fdelta;
		hig_get_facet_center(f, fcenter);
		hig_get_facet_delta(f, fdelta);

		real valuec = dp_get_value(dpu[dim], fgid);
		real tot = 0.0;
		for (int dim2 = 0; dim2 < DIM; dim2++) {
			real du2dx2 = 0.0;
			real dudx = 0.0;
			compute_facet_dudx_and_du2dx2_at_point(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], stn, valuec, &dudx, &du2dx2);
			real u;
			// Bug when value at center
			u = sfd_dp_interpolate(sfdu[dim2], dpu[dim2], fcenter, fcenter, stn);
			tot += - u * dudx + du2dx2 / Re;
		}
		real ustar = valuec + dt * tot;
		dp_set_value(dpustar[dim], fgid, ustar);
	}
	higfit_destroy(fit);

	START_CLOCK(syncprop);
	dp_sync(dpustar[dim]);
	STOP_CLOCK(syncprop);
}


solver *setup_p_solver(psim_domain *psdp, sim_stencil *stn) {
        /* This function sets up the solver for use as long as the matrix operaton doesn't change.*/

        /* Initialing the solver */
	solver *slvp = NULL;

        /* Number of elements */
	int localdomainsize = psd_get_local_domain_size(psdp);

        /* Creating the solver local domain */
	slvp = slv_create(SOLVER_ANY, psd_get_first_id(psdp), localdomainsize);

        /* Setting the maximum number of non zeros values */
	slv_set_maxnonzeros(slvp, 300);

	slv_set_maxiteration(slvp, 100);

        /* Declaring the cell iterator */
	higcit_celliterator *it;

        /* Getting the local domain */
	sim_domain *sdp = psd_get_local_domain(psdp);

	/* Mapping the local domain for the cell center */
	mp_mapper *mp = sd_get_domain_mapper(sdp);

	int total_nonzeros = 0;
	int cgid;

        /* Traversing the cells */
	for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {

                /* Getting the cell */
		hig_cell *c = higcit_getcell(it);

                /* Getting the cell center and cell size */
		Point ccenter;
		Point cdelta;
		hig_get_center(c, ccenter);
		hig_get_delta(c, cdelta);

                /* Reset the stencil */
		stn_reset(stn);

                /* Setting the parameters of the discretizated Poisson equation */
		real alpha = 0.0;

                /* Calculanting the parameters and the stencil points*/
		for(int dim = 0; dim < DIM; dim++) {
                        /* Parameter */
			real w = 1.0/(cdelta[dim]*cdelta[dim]);
			alpha -= 2.0 * w;

                        /* Assign the right and left point to the stencil */
			Point p;
			POINT_ASSIGN(p, ccenter);

                        /* Virtual right point */
			p[dim] = ccenter[dim] + cdelta[dim];
			sd_get_stencil(sdp, ccenter, p, w, stn);

                        /* Virtual left point */
			p[dim] = ccenter[dim] - cdelta[dim];
			sd_get_stencil(sdp, ccenter, p, w, stn);
		}
                /* Center point */
		sd_get_stencil(sdp, ccenter, ccenter, alpha, stn);

                /* Getting the number of real points */
		int numelems = stn_get_numelems(stn);

		total_nonzeros += numelems;

                /* Getting the coeficients of real points */
		real *vals = stn_get_vals(stn);

                /* Getting the identifiers of real points */
		int *ids = psd_stn_get_gids(psdp,stn);

                /* Getting the global identifier */
		cgid = psd_get_global_id(psdp, c);

               for(int i = 0;  i< numelems; ++i) {
                       if(ids[i] < 0) {
                               printf("Col ID of %d in cell %d (%lf, %lf, %lf), num of cols: %d n", ids[i], cgid, ccenter[0], ccenter[1], ccenter[2], numelems);
                       }
               }

                /* Setting the matrix line for this point */
		slv_set_Ai(slvp, cgid, numelems, ids, vals);
	}

	printf("Non-zeros: %d\n", total_nonzeros);

        /* Destroying the iterator */
	higcit_destroy(it);

	/* Fix one pressure point. Last point of process 0. */
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(rank == 0) {
		slv_impose_value(slvp, cgid, 0.0);
	}

        /* Assembling the local matrix to the global matrix */
	slv_assemble_matrix(slvp);

	/* Assemble the output vector.
	 * I guess that assembling it only once here, the output
	 * of the previous timestep will be used as the input of
	 * the first one.*/
	slv_assemble_output(slvp);

	return slvp;
}


void compute_p(psim_domain *psdp, psim_facet_domain *psfdu[DIM], distributed_property *dpp, distributed_property *dpustar[DIM], real dt, real Re, sim_stencil *stn) {
	static solver *slvp = NULL;
	if (slvp == NULL) {
		slvp = setup_p_solver(psdp, stn);
	}

	higcit_celliterator *it;
	sim_domain *sdp = psd_get_local_domain(psdp);
	sim_facet_domain *sfdu[DIM];
	for(int dim = 0; dim < DIM; dim++) {
		sfdu[dim] = psfd_get_local_domain(psfdu[dim]);
	}
	mp_mapper *mp = sd_get_domain_mapper(sdp);
	for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		Point ccenter;
		Point cdelta;
		hig_get_center(c, ccenter);
		hig_get_delta(c, cdelta);

		real sumdudx = 0.0;
		for(int dim = 0; dim < DIM; dim++) {
			real dudx;
			compute_facet_dudx_and_du2dx2_at_point(sfdu[dim], ccenter, cdelta, dim, 0.5, dpustar[dim], stn, 0.0, &dudx, NULL);
			sumdudx += dudx;
		}
		real rhs = sumdudx / dt;

		int cgid = psd_get_global_id(psdp, c);
		slv_set_bi(slvp, cgid, rhs);
	}
	higcit_destroy(it);
	slv_assemble_rhs(slvp);
	DEBUG_DIFF_TIME;
	START_CLOCK(poisson);
	slv_solve(slvp);
	STOP_CLOCK(poisson);
	DEBUG_DIFF_TIME;
	dp_slv_load_from_solver(dpp, slvp);

	printf("Solver stats: iteration count = %d; rel residual norm = %e.\n",
			slv_get_iteration_count(slvp),
			slv_get_relative_residual_norm(slvp));
}

void compute_ut(psim_domain *psdp, psim_facet_domain *psfdu[DIM], int dim, distributed_property *dpp, distributed_property *dpustar[DIM], distributed_property *dpu[dim], real dt, real Re, sim_stencil *stn) {
	sim_domain *sdp = psd_get_local_domain(psdp);
	sim_facet_domain *sfdu[DIM];
	higfit_facetiterator *fit;
	for(int dim = 0; dim < DIM; dim++) {
		sfdu[dim] = psfd_get_local_domain(psfdu[dim]);
		mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
		for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
			hig_facet *f = higfit_getfacet(fit);
			int fgid = mp_lookup(mu, hig_get_fid(f));
			Point fcenter;
			Point fdelta;
			hig_get_facet_center(f, fcenter);
			hig_get_facet_delta(f, fdelta);

			real dpdx;
			compute_dpdx_at_point(sdp, fcenter, fdelta, dim, 0.5, dpp, stn, &dpdx);
			real ustar = dp_get_value(dpustar[dim], fgid);
			real utdt = ustar - dt * dpdx;
			dp_set_value(dpu[dim], fgid, utdt);
		}
		higfit_destroy(fit);
		START_CLOCK(syncprop);
		dp_sync(dpu[dim]);
		STOP_CLOCK(syncprop);
	}
}

void updatePandU(psim_domain *psdp, psim_facet_domain *psfdu[DIM], distributed_property *dpp, distributed_property *dpu[DIM], distributed_property *dpustar[DIM], real dt, real Re, sim_stencil *stn) {
	higfit_facetiterator *fit;
	DEBUG_DIFF_TIME;
	for(int dim = 0; dim < DIM; dim++) {
	//	compute_ustar(psfdu, dim, dpu, dt, Re, dpustar, stn);
	}
	DEBUG_DIFF_TIME;
	compute_p(psdp, psfdu, dpp, dpustar, dt, Re, stn);
	DEBUG_DIFF_TIME;

	for(int dim = 0; dim < DIM; dim++) {
		compute_ut(psdp, psfdu, dim, dpp, dpustar, dpu, dt, Re, stn);
	}
	DEBUG_DIFF_TIME;
}

void printResults(psim_domain *psdp, psim_facet_domain *psfdu[DIM], distributed_property *dpp, distributed_property *dpu[DIM]) {
	sim_domain *sdp = psd_get_local_domain(psdp);
	sim_stencil *stn = stn_create();
	higcit_celliterator *it;
	int hig = 0;
	hig_cell *rootp = sd_get_higtree(sdp, hig);
	sim_facet_domain *sfdu[DIM];
	for(int dim = 0; dim < DIM; dim++) {
		sfdu[dim] = psfd_get_local_domain(psfdu[dim]);
	}
        Point l = {0.0, 0.0, 0.0};
        Point h = {1.0, 1.0, 1.0};
	Point center;
	POINT_SUB(center, h, l);
	POINT_DIV_SCALAR(center, center, 2.0);
	POINT_ADD(center, center, l);
	for(int dim = 0; dim < DIM; dim++) {
		if (dim != 1) {
			l[dim] = center[dim] + 3.0 * EPSDELTA;
			h[dim] = center[dim] + 5.0 * EPSDELTA;
		}
	}
	l[1] = l[1] + 2.0 * (h[1] - l[1])/3.0;

	for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		Point ccenter;
		Point cdelta;
		hig_get_center(c, ccenter);
		hig_get_delta(c, cdelta);
		real value[DIM];
		for(int dim = 0; dim < DIM; dim++) {
			value[dim] = sfd_dp_interpolate(sfdu[dim], dpu[dim], ccenter, ccenter, stn);
		}
		printf("v(");
		for(int dim = 0; dim < DIM; dim++) {
			printf("%f ", ccenter[dim]);
		}
		printf(") = ");
		for(int dim = 0; dim < DIM; dim++) {
			printf("%0.15f ", value[dim]);
		}
		printf("\n");
	}
	stn_destroy(stn);
}

void computeextraload(psim_domain *psdp, psim_facet_domain *psfdu[DIM], distributed_property *dpp, distributed_property *dpu[DIM], int extraload) {
	sim_domain *sdp = psd_get_local_domain(psdp);

	mp_mapper *mp = sd_get_domain_mapper(sdp);
	DECL_AND_ALLOC(Point, rndpts, extraload);
	DECL_AND_ALLOC(int, closestpts, extraload);
	higcit_celliterator *it;

	for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		Point lp, hp;
		hig_get_lowpoint(c, lp);
		hig_get_highpoint(c, hp);
		Point x;
		for (int i = 0; i < extraload; i++) {
			randpoint(lp, hp, x);
			POINT_ASSIGN(rndpts[i], x);
		}
		real maxpossibledist = co_distance(lp, hp);
		for (int i = 0; i < extraload; i++) {
			real shortestdist = maxpossibledist;
			int closestpoint = i;
			for (int j = 0; j < extraload; j++) {
				if (i != j) {
					real dist = co_distance(rndpts[i], rndpts[j]);
					if (dist < shortestdist) {
						shortestdist = dist;
						closestpoint = j;
					}
				}
			}
			closestpts[i] = closestpoint;
		}
	}
	free(rndpts);
	free(closestpts);
}


int main (int argc, char *argv[]) {
	if(argc < 6) {
		fputs("Missing parameters: amr_file, time_delta, reynolds, num_steps, extraload\n", stderr);
		return -1;
	}
	higtree_initialize(&argc, &argv);

	int myrank;
	int ntasks;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

	START_CLOCK(total);
	partition_graph *pg = pg_create(MPI_COMM_WORLD);
	pg_set_fringe_size(pg, 4);

	DEBUG_DIFF_TIME;

	char *amrfilename = argv[1];
	real dt = atof(argv[2]);
	real Re = atof(argv[3]);
	int numsteps = atoi(argv[4]);
	int extraload = atoi(argv[5]);
	DEBUG_INSPECT(numsteps, %d);

	/* Start the partitioning process. */
	load_balancer *lb = lb_create(MPI_COMM_WORLD, 1);

	higio_amr_info *mi;
	{
		FILE *fd = fopen(amrfilename, "r");
		mi = higio_read_amr_info(fd);
		fclose(fd);
	}
	if(myrank == 0) {
        	/* Build the initial tree */
		hig_cell *tree_domain = higio_read_from_amr_info(mi);

		/* Set the tree for partitioning. */
		lb_add_input_tree(lb, tree_domain, true, 0);
	}

	lb_calc_partition(lb, pg);

	DEBUG_DIFF_TIME
	sim_domain *sdp = sd_create(NULL);
	sd_set_interpolator_order(sdp, 2);

	unsigned num_trees = lb_get_num_local_trees(lb);
	for(unsigned i = 0; i < num_trees; ++i) {
		hig_cell *root = lb_get_local_tree(lb, i, NULL);
		sd_add_higtree(sdp, root);
	}
	psim_domain *psdp = psd_create(sdp, pg);
	psd_synced_mapper(psdp);
	distributed_property * dpp = psd_create_property(psdp);

	sim_facet_domain *sfdu[DIM];
	psim_facet_domain *psfdu[DIM];
	distributed_property *dpu[DIM];
	distributed_property *dpustar[DIM];
	for(int dim = 0; dim < DIM; dim++) {
		sfdu[dim] = sfd_create(NULL, dim);
		sfd_set_interpolator_order(sfdu[dim], 2);
		sfd_copy_higtrees_from_center_domain(sfdu[dim], sdp);
		for(unsigned i = 0; i < sfdu[dim]->cdom->numhigtrees; ++i)
		{
			char fname[100];
			snprintf(fname, 100, "%d-%d-%u.vtk", myrank, dim, i);
			FILE *fd = fopen(fname, "w");
			higio_print_in_vtk3d(fd, sfdu[dim]->cdom->higtrees[i]);
			fclose(fd);
		}
		sfd_adjust_facet_ids(sfdu[dim]);
		psfdu[dim] = psfd_create(sfdu[dim], psdp);
	}

	DEBUG_DIFF_TIME

	computeBCP(psdp, mi);
	computeBCU(psfdu, mi);

	higio_amr_info_destroy(mi);

	for(int dim = 0; dim < DIM; dim++) {
		psfd_compute_sfbi(psfdu[dim]);
		psfd_synced_mapper(psfdu[dim]);
		dpu[dim] = psfd_create_property(psfdu[dim]);
		dpustar[dim] = psfd_create_property(psfdu[dim]);
	}

	initializeP(psdp, dpp);
	initializeU(psfdu, dpu);
	sim_stencil *stn = stn_create();

	xdmf_output *output_ctx = xdmf_init("cavity_3d", sdp);
	xdmf_register_cell_property(output_ctx, "Pressure", 1, dpp);
	xdmf_register_facet_property(output_ctx, "Velocity", 3,
		sfdu[0], dpu[0], sfdu[1], dpu[1], sfdu[2], dpu[2]);

	for(int step = 0; step < numsteps; step++) {
		//printResults(psdp, psfdu, dpp, dpu);
		if (step == 0) { START_CLOCK(firstiter); }
		updatePandU(psdp, psfdu, dpp, dpu, dpustar, dt, Re, stn);
		if (step == 0) { STOP_CLOCK(firstiter); }
		//computeextraload(psdp, psfdu, dpp, dpu, extraload);
		xdmf_write_timestep(output_ctx, step, step);
		DEBUG_INSPECT(step, %d);
	}
	xdmf_destroy(output_ctx);
	for(int dim = 0; dim < DIM; dim++) {
		dp_destroy(dpustar[dim]);
	}
	stn_destroy(stn);
	STOP_CLOCK(total);
	DEBUG_DIFF_TIME;
	DEBUG_INSPECT(GET_NSEC_CLOCK(total)/1.0e9, %lf);
	DEBUG_INSPECT(GET_NSEC_CLOCK(poisson)/1.0e9, %lf);
	DEBUG_INSPECT(GET_NSEC_CLOCK(syncprop)/1.0e9, %lf);
	DEBUG_INSPECT(GET_NSEC_CLOCK(firstiter)/1.0e9, %lf);
	return 0;
}
