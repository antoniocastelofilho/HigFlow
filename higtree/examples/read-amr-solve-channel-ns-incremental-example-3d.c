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

static void
computeBCP(psim_domain *psd, int numbcs, load_balancer *lb, int pbctypes[], real pbcvalues[], int is_incr, real dt) {
    sim_domain *sd = psd_get_local_domain(psd);

    for(int h = 0; h < numbcs; h++) {
	unsigned num_trees = lb_get_num_local_trees_in_group(lb, h+1);
	for(unsigned i = 0; i < num_trees; ++i) {
		hig_cell *tree =
			lb_get_local_tree_in_group(lb, i, h+1);

		int bc_t = ((pbctypes[h]==0)?DIRICHLET:NEUMANN);
		sim_boundary *bc = make_bc(tree, bc_t);

		/* Adding the boundary condition */
		sd_add_boundary(sd, bc);
		mp_mapper *bm = sb_get_mapper(bc);

		higcit_celliterator *it;

		/* Traversing the cells of the Diriclet boundaries conditions */
		for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {

		   /* Getting the cell */
		    hig_cell *bcell = higcit_getcell(it);

		    /* Set the value 0.0 for the center of the cell */
		    int bcgid = mp_lookup(bm, hig_get_cid(bcell));
		    real bcval;
		    bcval = pbcvalues[h];

		    sb_set_value(bc, bcgid, bcval);
		}
		higcit_destroy(it);
	    }
    }
}

static void
computeBCU(psim_facet_domain *psfdu[DIM], int numbcs, load_balancer *lb, int ubctypes[], real ubcvalues[], int vbctypes[], real vbcvalues[], int wbctypes[], real wbcvalues[]) {
    int *bctypes[3] = {ubctypes, vbctypes, wbctypes};
    real *bcvalues[3] = {ubcvalues, vbcvalues, wbcvalues};

    sim_facet_domain *sfd;
    for(int h = 0; h < numbcs; h++) {
	unsigned num_trees = lb_get_num_local_trees_in_group(lb, h+1);
	for(unsigned i = 0; i < num_trees; ++i) {
		hig_cell *tree =
			lb_get_local_tree_in_group(lb, i, h+1);

		for(int dim = 0; dim < 3; dim++) {
		    sfd = psfd_get_local_domain(psfdu[dim]);

		    int bc_t = ((bctypes[dim][h]==0)?DIRICHLET:NEUMANN);
		    sim_boundary *bc = make_bc(tree, bc_t);

		    /* Adding the boundary condition */
		    sfd_add_boundary(sfd, bc);
		    mp_mapper *bm = sb_get_mapper(bc);

		    higcit_celliterator *it;

		    /* Traversing the cells of the Diriclet boundaries conditions */
		    for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {

		       /* Getting the cell */
			hig_cell *bcell = higcit_getcell(it);

		       /* Set the value 0.0 for the center of the cell */
			int bcgid = mp_lookup(bm, hig_get_cid(bcell));
			real bcval;
			bcval = bcvalues[dim][h];

			sb_set_value(bc, bcgid, bcval);
		    }
		    higcit_destroy(it);
		}
	}
    }
}

void initializeP(psim_domain *psdp, distributed_property *dpp) {

    /* Setting the cell iterator */
    higcit_celliterator *it;
    /* Getting the local domain */
    sim_domain *sdp = psd_get_local_domain(psdp);
    /* Getting the Mapper for the local domain */
    mp_mapper *m = sd_get_domain_mapper(sdp);
    /* Traversing the cells of local domain */
    for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {

        /* Getting the cell */
        hig_cell *c = higcit_getcell(it);

        /* Get the cell identifier */
        int cgid = mp_lookup(m, hig_get_cid(c));

        /* Set the value 0.0 for the pressure in this cell */
        dp_set_value(dpp, cgid, 0.0);
    }
    /* Destroying the iterator */
    higcit_destroy(it);
}


void initializeU(psim_facet_domain *psfdu[DIM], distributed_property *dpu[DIM]) {

    /* Setting the facet-cell iterator */
    higfit_facetiterator *fit;

    /* Setting the volocity values for the domain */
    sim_facet_domain *sfdu[DIM];

    /* Setting the velocity U(dim) */
    for(int dim = 0; dim < DIM; dim++) {
        /* Get the local domain property */
        sfdu[dim] = psfd_get_local_domain(psfdu[dim]);
        /* Get the Mapper for the local domain */
        mp_mapper *m = sfd_get_domain_mapper(sfdu[dim]);

        for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            hig_facet *f = higfit_getfacet(fit);
            int fgid = mp_lookup(m, hig_get_fid(f));
            real val;
            val = 0.0;
            dp_set_value(dpu[dim], fgid, val);
        }

        higfit_destroy(fit);
    }
}

void compute_facet_dudx_and_du2dx2_at_point(sim_facet_domain *sfdu, Point center, Point delta, int dim, real alpha, distributed_property *dpu, sim_stencil *stn, real valueatpoint, real *dudx, real *du2dx2) {

    //DEBUG_PASS;
    Point p;
    POINT_ASSIGN(p, center);
    p[dim] = center[dim] - alpha * delta[dim];
    //DEBUG_PASS;
    real valuel = sfd_dp_interpolate(sfdu, dpu, center, p, stn);
    p[dim] = center[dim] + alpha * delta[dim];
    //DEBUG_PASS;
    real valueh = sfd_dp_interpolate(sfdu, dpu, center, p, stn);
    //DEBUG_PASS;
    if (dudx != NULL) {
        *dudx = (valueh - valuel) / (2.0 * alpha * delta[dim]);
    }
    //DEBUG_PASS;
    if (du2dx2 != NULL) {
        *du2dx2 = (valuel - 2.0*valueatpoint + valueh)/(delta[dim] * delta[dim] * alpha * alpha);
    }
    //DEBUG_PASS;
}

void compute_dpdx_at_point(sim_domain *sdp, Point center, Point delta, int dim, real alpha, distributed_property *dpp, sim_stencil *stn, real *dpdx) {
    Point p;
    POINT_ASSIGN(p, center);
    p[dim] = center[dim] - alpha * delta[dim];
    real valuel = sd_dp_interpolate(sdp, dpp, center, p, stn);
    p[dim] = center[dim] + alpha * delta[dim];
    real valueh = sd_dp_interpolate(sdp, dpp, center, p, stn);
    if (dpdx != NULL) {
        *dpdx = (valueh - valuel) / (2.0 * alpha * delta[dim]);
    }
}

void compute_ustar(psim_facet_domain *psfdu[DIM], int dim, distributed_property *dpu[DIM], real t, real dt, real Re, distributed_property *dpustar[DIM], sim_stencil *stn) {

    higfit_facetiterator *fit;
    sim_facet_domain *sfdu[DIM];
    // Note que o uso de dim esta errado aqui
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
            u = sfd_dp_interpolate(sfdu[dim2], dpu[dim2],fcenter, fcenter, stn);
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

void compute_ustar_incremental(psim_domain *psdp, psim_facet_domain *psfdu[DIM],  int dim, distributed_property *dpp,distributed_property *dpu[DIM], real t, real dt, real Re, distributed_property *dpustar[DIM], sim_stencil *stn) {

    higfit_facetiterator *fit;

    sim_domain *sdp = psd_get_local_domain(psdp);
    sim_facet_domain *sfdu[DIM];

    for(int dim2 = 0; dim2 < DIM; dim2++) {
        sfdu[dim2] = psfd_get_local_domain(psfdu[dim2]);
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
        real dpdx;
        compute_dpdx_at_point(sdp, fcenter, fdelta, dim, 0.5, dpp, stn, &dpdx);
        tot += -dpdx;
        //tot += 0.0;

        for (int dim2 = 0; dim2 < DIM; dim2++) {
            //DEBUG_INSPECT(dim2, %d);
            real du2dx2 = 0.0;
            real dudx = 0.0;
            //DEBUG_PASS;
            compute_facet_dudx_and_du2dx2_at_point(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], stn, valuec, &dudx, &du2dx2);
            real u;
            // Bug when value at center
            //DEBUG_PASS;
            u = sfd_dp_interpolate(sfdu[dim2], dpu[dim2], fcenter, fcenter, stn);
            tot += - u * dudx + du2dx2 / Re;
            //DEBUG_PASS;
        }
        //DEBUG_PASS;
        real ustar = valuec + dt * tot;
        dp_set_value(dpustar[dim], fgid, ustar);
    }
    higfit_destroy(fit);
    START_CLOCK(syncprop);
    dp_sync(dpustar[dim]);
    STOP_CLOCK(syncprop);
}

// Maybe remove_pressure_singularity - it will depend on the boundary situation

void compute_p(psim_domain *psdp, psim_facet_domain *psfdu[DIM], distributed_property *dpp, distributed_property *dpn, distributed_property *dpustar[DIM], real t, real dt, real Re, int is_incr, sim_stencil *stn) {
    static solver *slvp = NULL;
    if (slvp == NULL) {
        int localdomainsize = psd_get_local_domain_size(psdp);
        slvp = slv_create(SOLVER_ANY, psd_get_first_id(psdp), localdomainsize);               //Creates a solver. If islocalsize is 1, then _size is the number of elements in the local domain. If islocalsize is 0, then _size is the number of elements in the global domain.
        slv_set_maxnonzeros(slvp, 300);
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
        stn_reset(stn);
        stn_set_rhs(stn, sumdudx / dt);
        real alpha = 0.0;
        for(int dim = 0; dim < DIM; dim++) {
            real w = 1.0/(cdelta[dim]*cdelta[dim]);
            alpha -= 2.0 * w; // vou ter que colocar a soma dos rhos
            Point p;
            POINT_ASSIGN(p, ccenter);
            p[dim] = ccenter[dim] + cdelta[dim];
            sd_get_stencil(sdp, ccenter, p, w, stn);
            p[dim] = ccenter[dim] - cdelta[dim];
            sd_get_stencil(sdp, ccenter, p, w, stn);
        }
        sd_get_stencil(sdp, ccenter, ccenter, alpha, stn);

        int *ids   = psd_stn_get_gids(psdp, stn);
        real *vals = stn_get_vals(stn);
        int numelems = stn_get_numelems(stn);

        int cgid = psd_get_global_id(psdp, c);

        slv_set_bi(slvp, cgid, stn_get_rhs(stn));
        slv_set_Ai(slvp, cgid, numelems, ids, vals);
    }

    higcit_destroy(it);

    slv_assemble(slvp);

    // aqui é pra remover
    // remove_pressure_singularity(psdp, t, dt, is_incr, slvp);
    //DEBUG_DIFF_TIME;
    START_CLOCK(poisson);
    slv_solve(slvp);
    STOP_CLOCK(poisson);
    //DEBUG_DIFF_TIME;

    dp_slv_load_from_solver(dpp, slvp);
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

void compute_pt(psim_domain *psdp, distributed_property *ddeltap, distributed_property *dpp, real dt, sim_stencil *stn) {

    sim_domain *sdp = psd_get_local_domain(psdp);

    higcit_celliterator *it;
    mp_mapper *mp = sd_get_domain_mapper(sdp);

    for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);

        int cgid = mp_lookup(mp, hig_get_cid(c));
        real p      = dp_get_value(dpp, cgid);
        real deltap = dp_get_value(ddeltap, cgid);
        real newp   = p + deltap;
        dp_set_value(dpp, cgid, p + deltap);

    }
    higcit_destroy(it);
    START_CLOCK(syncprop);
    dp_sync(dpp);
    STOP_CLOCK(syncprop);
}

void updatePandU(psim_domain *psdp, psim_facet_domain *psfdu[DIM], distributed_property *dpp, distributed_property *dpu[DIM], distributed_property *dpustar[DIM], real tn, real dt, real Re,sim_stencil *stn) {

    int is_incr = 0;
    higfit_facetiterator *fit;
    //DEBUG_DIFF_TIME;

    //updateBCU(psfdu, tn+dt);

    for(int dim = 0; dim < DIM; dim++) {
        compute_ustar(psfdu, dim, dpu, tn, dt, Re, dpustar, stn);
    }
    //updateBCP(psdp, t);
    //DEBUG_DIFF_TIME;
   // updateBCP(psdp, tn+dt);
    compute_p(psdp, psfdu, dpp, dpp, dpustar, tn, dt, Re, is_incr, stn);
    //DEBUG_DIFF_TIME;

    for(int dim = 0; dim < DIM; dim++) {
        compute_ut(psdp, psfdu, dim, dpp, dpustar, dpu, dt, Re, stn);
    }
    //DEBUG_DIFF_TIME;
}

void updatePandU_incremental(psim_domain *psdp, psim_facet_domain *psfdu[DIM], distributed_property *dpp, distributed_property *ddeltap, distributed_property *dpu[DIM], distributed_property *dpustar[DIM], real tn, real dt, real Re, sim_stencil *stn) {

    int is_incr = 1;
    higfit_facetiterator *fit;
    //DEBUG_DIFF_TIME;

    // Segundo Guermond tem que ser no tempo n+1
    //updateBCU(psfdu, tn+dt);

    for(int dim = 0; dim < DIM; dim++) {
        compute_ustar_incremental(psdp, psfdu, dim, dpp, dpu, tn, dt, Re, dpustar, stn);
    }
    //updateBCP_incremental(psdp, t, dt);
    //DEBUG_DIFF_TIME;
    //updateBCP_incremental(psdp, tn, dt);
    compute_p(psdp, psfdu, ddeltap, dpp, dpustar, tn, dt, Re, is_incr, stn);
    //DEBUG_DIFF_TIME;
    for(int dim = 0; dim < DIM; dim++) {
        compute_ut(psdp, psfdu, dim, ddeltap, dpustar, dpu, dt, Re, stn);
        // Tem que tomar cuidado com o gradiente de pressao nas faces do contorno e com a condicao de contorno para velocidade
    }
    // updateBCP(psdp, tn+dt);
    compute_pt(psdp, ddeltap, dpp, dt, stn);
    //DEBUG_DIFF_TIME;

}

void write_property_cell(hig_cell *root) {
	higcit_celliterator *it;
	it = higcit_create_all_leaves(root);
	int numcells = higcit_count(it);
	higcit_destroy(it);
	printf("CELL_DATA %d\nSCALARS u float 1\nLOOKUP_TABLE default\n", numcells);

	higcit_destroy(it);
}

void hig_printVTK3dcell(psim_domain *psdp, distributed_property *dpp, distributed_property  *dpu[DIM], psim_facet_domain *psfdu[DIM], int numhigs, const char *basefilename, const char *propertyname, int rank, int frame)
{

    /* Open the VTK file */
    char vtkname[1024];
    sprintf(vtkname, "%s_%s_%d-%d.vtk", basefilename, propertyname, rank, frame);
    FILE *f = fopen(vtkname, "w");
    if (f == NULL) {
        return;
    }

    sim_domain *sdp = psd_get_local_domain(psdp);

	higcit_celliterator *it_t;
	// it_t = higcit_create_all_leaves(root);
	it_t = sd_get_domain_celliterator(sdp);
	long numleafs = higcit_count_without_advancing(it_t);
	fprintf(f, "# vtk DataFile Version 3.0\n");
	fprintf(f, "higtree\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(f, "POINTS %ld float\n", 8*numleafs);
	hig_cell *c;
	for (; !higcit_isfinished(it_t); higcit_nextcell(it_t)) {
		c = higcit_getcell(it_t);
		fprintf(f, "%lf %lf %lf\n", c->lowpoint [0], c->lowpoint [1], c->lowpoint [2]);
		fprintf(f, "%lf %lf %lf\n", c->highpoint[0], c->lowpoint [1], c->lowpoint [2]);
		fprintf(f, "%lf %lf %lf\n", c->highpoint[0], c->highpoint[1], c->lowpoint [2]);
		fprintf(f, "%lf %lf %lf\n", c->lowpoint [0], c->highpoint[1], c->lowpoint [2]);
		fprintf(f, "%lf %lf %lf\n", c->lowpoint [0], c->lowpoint [1], c->highpoint[2]);
		fprintf(f, "%lf %lf %lf\n", c->highpoint[0], c->lowpoint [1], c->highpoint[2]);
		fprintf(f, "%lf %lf %lf\n", c->highpoint[0], c->highpoint[1], c->highpoint[2]);
		fprintf(f, "%lf %lf %lf\n", c->lowpoint [0], c->highpoint[1], c->highpoint[2]);
	}
	fprintf(f, "CELLS %ld %ld\n", numleafs, 9*numleafs);
	for(int i = 0; i < numleafs; i++) {
		fprintf(f, "%d %d %d %d %d %d %d %d %d\n", 8, 8*i+0, 8*i+1, 8*i+2, 8*i+3, 8*i+4, 8*i+5, 8*i+6, 8*i+7); //x == 0
		// fprintf(force, "%d %d %d %d %d\n", 4, 8*i+1, 8*i+2, 8*i+6, 8*i+5); //x == 1
		// fprintf(force, "%d %d %d %d %d\n", 4, 8*i+0, 8*i+1, 8*i+5, 8*i+4); //y == 0
		// fprintf(force, "%d %d %d %d %d\n", 4, 8*i+3, 8*i+2, 8*i+6, 8*i+7); //y == 1
		// fprintf(force, "%d %d %d %d %d\n", 4, 8*i+0, 8*i+1, 8*i+2, 8*i+3); //z == 0
		// fprintf(force, "%d %d %d %d %d\n", 4, 8*i+4, 8*i+5, 8*i+6, 8*i+7); //z == 1
	}
	fprintf(f, "CELL_TYPES %ld\n", numleafs);
	for(int i = 0; i < numleafs; i++) {
		fprintf(f, "12 "); //x == 0
		// fprintf(force, "%d %d %d %d %d\n", 4, 8*i+1, 8*i+2, 8*i+6, 8*i+5); //x == 1
		// fprintf(force, "%d %d %d %d %d\n", 4, 8*i+0, 8*i+1, 8*i+5, 8*i+4); //y == 0
		// fprintf(force, "%d %d %d %d %d\n", 4, 8*i+3, 8*i+2, 8*i+6, 8*i+7); //y == 1
		// fprintf(force, "%d %d %d %d %d\n", 4, 8*i+0, 8*i+1, 8*i+2, 8*i+3); //z == 0
		// fprintf(force, "%d %d %d %d %d\n", 4, 8*i+4, 8*i+5, 8*i+6, 8*i+7); //z == 1
	}
	higcit_destroy(it_t);


	// higcit_celliterator *it = higcit_create_all_leaves(root);
	higcit_celliterator *it = sd_get_domain_celliterator(sdp);
	long numcells = higcit_count(it);

	printf("\n ---- Number of cells: %lu\n", numcells);

	higcit_destroy(it);
	// sim_domain *sdp = psd_get_local_domain(psdp);
	mp_mapper *m = sd_get_domain_mapper(sdp);

	sim_facet_domain *sfdu2[DIM];
	int dimension[DIM]; dimension[0] = 0; dimension[1] = 0; dimension[2] = 0;

	int degree = 2;
	const int maxpts = 2*DIM*wls_num_min_points(degree);
	Point pts[maxpts+1];
	real dist[maxpts+1];
	real w0[maxpts], w1[maxpts], w2[maxpts], w3[maxpts], w4[maxpts], w5[maxpts], w6[maxpts], w7[maxpts];
	uniqueid gids[maxpts];

    fprintf(f, "\nPOINT_DATA %ld\nVECTORS vel FLOAT\n", 8*numleafs);
    for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);

		Point cdelta, ccenter, clowpoint, chightpoint;

		hig_get_delta(c,cdelta);
		hig_get_center(c,ccenter);

		clowpoint[0] = ccenter[0] - 2.0*cdelta[0];
		clowpoint[1] = ccenter[1] - 2.0*cdelta[1];
		clowpoint[2] = ccenter[2] - 2.0*cdelta[2];

		chightpoint[0] = ccenter[0] + 2.0*cdelta[0];
		chightpoint[1] = ccenter[1] + 2.0*cdelta[1];
		chightpoint[2] = ccenter[2] + 2.0*cdelta[2];

		// Pontos onde será interpolada a velocidade
		Point p0, p1, p2, p3, p4, p5, p6, p7;
		p0[0] = c->lowpoint[0];  p0[1] = c->lowpoint[1];  p0[2] = c->lowpoint[2];
		p1[0] = c->highpoint[0]; p1[1] = c->lowpoint[1];  p1[2] = c->lowpoint[2];
		p2[0] = c->highpoint[0]; p2[1] = c->highpoint[1]; p2[2] = c->lowpoint[2];
		p3[0] = c->lowpoint[0];  p3[1] = c->highpoint[1]; p3[2] = c->lowpoint[2];
		p4[0] = c->lowpoint[0];  p4[1] = c->lowpoint[1];  p4[2] = c->highpoint[2];
		p5[0] = c->highpoint[0]; p5[1] = c->lowpoint[1];  p5[2] = c->highpoint[2];
		p6[0] = c->highpoint[0]; p6[1] = c->highpoint[1]; p6[2] = c->highpoint[2];
		p7[0] = c->lowpoint[0];  p7[1] = c->highpoint[1]; p7[2] = c->highpoint[2];

		uniqueid id = hig_get_cid(c);
		int cgid = mp_lookup(m, id);

		double lu0 = 0, lv0 = 0, lw0 = 0;
		double lu1 = 0, lv1 = 0, lw1 = 0;
		double lu2 = 0, lv2 = 0, lw2 = 0;
		double lu3 = 0, lv3 = 0, lw3 = 0;
		double lu4 = 0, lv4 = 0, lw4 = 0;
		double lu5 = 0, lv5 = 0, lw5 = 0;
		double lu6 = 0, lv6 = 0, lw6 = 0;
		double lu7 = 0, lv7 = 0, lw7 = 0;
	wls_interpolator *wls;


            for(int dim = 0; dim < DIM; dim++)
			{

				    int numpts = 0;
				    int cont = 0;

                for (int hig = 0; hig < numhigs; hig++){

                    hig_cell *rootp = sd_get_higtree(sdp, hig);
				    sfdu2[dim] = psfd_get_local_domain(psfdu[dim]);
				    mp_mapper *m = sfd_get_domain_mapper(sfdu2[dim]);

				    higfit_facetiterator *fitn;

    				dimension[dim] = 1;

	       			//verificar passagem root --> rootp
			     	for(fitn = higfit_create_bounding_box_facets(rootp, dimension, clowpoint, chightpoint); !higfit_isfinished(fitn); higfit_nextfacet(fitn))
				    {
					   hig_facet *fn = higfit_getfacet(fitn);
					   Point fcenter;
					   hig_get_facet_center(fn, fcenter);
					   int fgid = mp_lookup(m, hig_get_fid(fn));
					   if (fgid < 0)
					   {
						  continue;
					   }

					   int distance = co_distance(fcenter, ccenter);
					   int j;
					   for(j = cont - 1; j >= 0 && dist[j] > distance; j--) {
					   		gids[j+1] = gids[j];
					   		// Guardando os centros das facetas para calcular os pesos depois.
					   		POINT_ASSIGN(pts[j+1], pts[j]);
					   		dist[j+1] = dist[j];
					   }
					   // Guardando o id da faceta para pegar o valor da velocidade depois.
					   gids[j+1] = fgid;
					   // Guardando os centros das facetas para calcular os pesos depois.
					   POINT_ASSIGN(pts[j+1], fcenter);
					   dist[j+1] = distance;


					   if (cont < maxpts) {
						cont++;
					   }

				    }
				    higfit_destroy(fitn);
				    if (cont > maxpts) {
					cont = maxpts;
				    }
			}

				    wls = wls_create(degree, maxpts);
				    wls_set_points(wls, cont, pts, p0);
				    wls_calc(wls, p0, cont, w0);

				    for(int i = 0; i < cont; i++)
				    {
					   real val = dp_get_value(dpu[dim], gids[i]);
					   switch (dim)
					   {
						  case 0:
							lu0 = lu0 + w0[i]*val;
						  break;
						  case 1:
						      	lv0 = lv0 + w0[i]*val;
						  break;
						  case 2:
							lw0 = lw0 + w0[i]*val;
						  break;
					   }
				    }
    				wls_destroy(wls);

    				wls = wls_create(degree, maxpts);
    				wls_set_points(wls, cont, pts, p1);
    				wls_calc(wls, p1, cont, w1);

    				for(int i = 0; i < cont; i++)
    				{
    					real val = dp_get_value(dpu[dim], gids[i]);
    					switch (dim)
    					{
    						case 0:
    							lu1 = lu1 + w1[i]*val;
    						break;
    						case 1:
    							lv1 = lv1 + w1[i]*val;
    						break;
    						case 2:
    							lw1 = lw1 + w1[i]*val;
    						break;
    					}
    				}
    				wls_destroy(wls);

    				wls = wls_create(degree, maxpts);
    				wls_set_points(wls, cont, pts, p2);
    				wls_calc(wls, p2, cont, w2);

    				for(int i = 0; i < cont; i++)
    				{
    					real val = dp_get_value(dpu[dim], gids[i]);
    					switch (dim)
    					{
    						case 0:
    							lu2 = lu2 + w2[i]*val;
    						break;
    						case 1:
    							lv2 = lv2 + w2[i]*val;
    						break;
    						case 2:
    							lw2 = lw2+ w2[i]*val;
    						break;
    					}
    				}
    				wls_destroy(wls);

    				wls = wls_create(degree, maxpts);
    				wls_set_points(wls, cont, pts, p3);
    				wls_calc(wls, p3, cont, w3);

    				for(int i = 0; i < cont; i++)
    				{
    					real val = dp_get_value(dpu[dim], gids[i]);
    					switch (dim)
    					{
    						case 0:
    							lu3 = lu3 + w3[i]*val;
    						break;
    						case 1:
    							lv3 = lv3 + w3[i]*val;
    						break;
    						case 2:
    							lw3 = lw3 + w3[i]*val;
    						break;
    					}
    				}
    				wls_destroy(wls);

    				wls = wls_create(degree, maxpts);
    				wls_set_points(wls, cont, pts, p4);
    				wls_calc(wls, p4, cont, w4);

    				for(int i = 0; i < cont; i++)
    				{
    					real val = dp_get_value(dpu[dim], gids[i]);
    					switch (dim)
    					{
    						case 0:
    							lu4 = lu4 + w4[i]*val;
    						break;
    						case 1:
    							lv4 = lv4 + w4[i]*val;
    						break;
    						case 2:
    							lw4 = lw4 + w4[i]*val;
    						break;
    					}
    				}
    				wls_destroy(wls);

    				wls = wls_create(degree, maxpts);
    				wls_set_points(wls, cont, pts, p5);
    				wls_calc(wls, p5, cont, w5);

    				for(int i = 0; i < cont; i++)
    				{
    					real val = dp_get_value(dpu[dim], gids[i]);
    					switch (dim)
    					{
    						case 0:
    							lu5 = lu5 + w5[i]*val;
    						break;
    						case 1:
    							lv5 = lv5 + w5[i]*val;
    						break;
    						case 2:
    							lw5 = lw5 + w5[i]*val;
    						break;
    					}
    				}
    				wls_destroy(wls);

    				wls = wls_create(degree, maxpts);
    				wls_set_points(wls, cont, pts, p6);
    				wls_calc(wls, p6, cont, w6);

    				for(int i = 0; i < cont; i++)
    				{
    					real val = dp_get_value(dpu[dim], gids[i]);
    					switch (dim)
    					{
    						case 0:
    							lu6 = lu6 + w6[i]*val;
    						break;
    						case 1:
    							lv6 = lv6 + w6[i]*val;
    						break;
    						case 2:
    							lw6 = lw6 + w6[i]*val;
    						break;
    					}
    				}
    				wls_destroy(wls);

    				wls = wls_create(degree, maxpts);
    				wls_set_points(wls, cont, pts, p7);
    				wls_calc(wls, p7, cont, w7);

    				for(int i = 0; i < cont; i++)
    				{
    					real val = dp_get_value(dpu[dim], gids[i]);
    					switch (dim)
    					{
    						case 0:
    							lu7 = lu7 + w7[i]*val;
    						break;
    						case 1:
    							lv7 = lv7 + w7[i]*val;
    						break;
    						case 2:
    							lw7 = lw7 + w7[i]*val;
    						break;
    					}
    				}
    				wls_destroy(wls);

    				dimension[dim] = 0;

    				//printf("Máximo de pontos --> %d\nMínimo de pontos --> %d\nPontos usados --> %d\n",maxpts,wls_num_min_points(degree),cont);
            }
    		fprintf(f, "%e %e %e\n", lu0, lv0, lw0);
    		fprintf(f, "%e %e %e\n", lu1, lv1, lw1);
    		fprintf(f, "%e %e %e\n", lu2, lv2, lw2);
    		fprintf(f, "%e %e %e\n", lu3, lv3, lw3);
    		fprintf(f, "%e %e %e\n", lu4, lv4, lw4);
    		fprintf(f, "%e %e %e\n", lu5, lv5, lw5);
    		fprintf(f, "%e %e %e\n", lu6, lv6, lw6);
    		fprintf(f, "%e %e %e\n", lu7, lv7, lw7);
        }
        higcit_destroy(it);

        fprintf(f, "\nCELL_DATA %ld\nSCALARS p FLOAT\nLOOKUP_TABLE default\n", numcells);
        for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            hig_cell *c = higcit_getcell(it);

            uniqueid id = hig_get_cid(c);
            int cgid = mp_lookup(m, id);

            real val = dp_get_value(dpp, cgid);

            fprintf(f, "%e\n", val);

            // fprintf(f, "%f\n", val);
        }
        higcit_destroy(it);
}

void writePBC(psim_domain *psdp, const char *basefilename, const char *propertyname, int rank) {
    /* Getting the local domain */
    sim_domain *sdp = psd_get_local_domain(psdp);

    /* Open the VTK format file */
    char vtkname[1024];
    sprintf(vtkname, "%s_%s_%d.vtk", basefilename, propertyname, rank);
    FILE *fdout = fopen(vtkname, "w");
    if (fdout == NULL) {
        return;
    }

    /* Getting the cell for print */
    higcit_celliterator *it;
    it = sd_get_bcs_celliterator(sdp);

    /* Printing the property */
    higio_print_celliterator_in_vtk3d(fdout, it);

    /* Destroying the iterator */
    higcit_destroy(it);

    /* Close the file */
    fclose(fdout);
}


int main (int argc, char *argv[])
{
    /* Initializing Higtree */
    higtree_initialize(&argc, &argv);

    int myrank, ntasks;
    int mpierr = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    /* Partitioning the domain */
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    /* Seeing the success or not of initializing */
    DEBUG_INSPECT(mpierr, %d);
    DEBUG_INSPECT(MPI_SUCCESS, %d);
    DEBUG_INSPECT(myrank, %d);

    START_CLOCK(total);
    /* Taking the sub-domain of the myrank process */
    partition_graph *pg = pg_create(MPI_COMM_WORLD); //Creates a partition table.
    /* Setting the fringe size of the sub-domain */
    pg_set_fringe_size(pg, 5);     //Sets the size of the fringe. The fringe is a buffer around the cells of a given node.

    int numhigs = atoi(argv[1]);
    DEBUG_INSPECT(numhigs, %d);

    argv++;
    char *amrfilename[numhigs];

    for(int h = 0; h < numhigs; h++) {
        amrfilename[h] = argv[1]; argv++;
    }

    int total_numbcs = atoi(argv[1]);
    DEBUG_INSPECT(total_numbcs, %d);

    argv++;
    char *amrBCfilename[total_numbcs];
    int pbctypes[total_numbcs];
    int ubctypes[total_numbcs];
    int vbctypes[total_numbcs];
    int wbctypes[total_numbcs];
    real pbcvalues[total_numbcs];
    real ubcvalues[total_numbcs];
    real vbcvalues[total_numbcs];
    real wbcvalues[total_numbcs];

    for(int h = 0; h < total_numbcs; h++) {
        amrBCfilename[h] = argv[1]; argv++;
        //printf("--- %s\n", amrBCfilename[h] );
        pbctypes[h] = atoi(argv[1]); argv++;
        pbcvalues[h] = atof(argv[1]); argv++;
        ubctypes[h] = atoi(argv[1]); argv++;
        ubcvalues[h] = atof(argv[1]); argv++;
        vbctypes[h] = atoi(argv[1]); argv++;
        vbcvalues[h] = atof(argv[1]); argv++;
        wbctypes[h] = atoi(argv[1]); argv++;
        wbcvalues[h] = atof(argv[1]); argv++;
        //printf("#### %d %.2f %d %.2f %d %.2f %d %.2f \n", pbctypes[h], pbcvalues[h], ubctypes[h], ubcvalues[h], vbctypes[h], vbcvalues[h], wbctypes[h], wbcvalues[h]);
    }

    real dt = atof(argv[1]); argv++;
    real Re = atof(argv[1]); argv++;
    int numsteps = atoi(argv[1]); argv++;
    int incremental = atoi(argv[1]); argv++;
    //int perstep = atoi(argv[1]); argv++;
    char *vtkresults = argv[1]; argv++;
    FILE *fd;

    int perstep;
    perstep = 10;

    assert(numhigs == 1 || ntasks == 1);

    // simulation domain (SD)
    sim_domain *sdp = sd_create(NULL);

    //reuse interpolation, 0 on, 1 off
    sd_use_cache(sdp, 1);


    //Sets the order of the interpolation to bhe used for the SD.
    sd_set_interpolator_order(sdp, 3);

    higio_amr_info *mi[numhigs];

    for(int h = 0; h < numhigs; h++) {
        /* Open the AMR format file */
        FILE *fd = fopen(amrfilename[h], "r");
        /* Reanding the IO information of the file */
        mi[h] = higio_read_amr_info(fd);
        fclose(fd);
    }

    /* Create the partitioning context. */
    unsigned num_groups = 1 + total_numbcs;
    load_balancer *lb = lb_create(MPI_COMM_WORLD, num_groups);
    for(unsigned i = 1; i < num_groups; ++i)
        lb_set_group_weight(lb, i, 0.0001);

    /* Loading domain trees from AMR information */
    for(int h = 0; h < numhigs; h++) {
        /* Creating the HigTree data structure */
        hig_cell *root = higio_read_from_amr_info(mi[h]);

	/* Add tree to partition. */
	lb_add_input_tree(lb, root, true, 0);
    }

    // Load boundary condition trees
    for(int h = myrank; h < total_numbcs; h += ntasks) {
	/* Open the AMR format file */
        fd = fopen(amrBCfilename[h], "r");

        /* Creating the distributed HigTree data structure */
        hig_cell *root = higio_read_from_amr(fd);

        /* Closing the AMR format file */
        fclose(fd);

	lb_add_input_tree(lb, root, true, h+1);
    }

    /*Computes a partition of the hig-tree root in nparts, so that the number of leaves in each partition are as balanced as possible. */
    lb_calc_partition(lb, pg);

    /* Bind the partitioned grid to the domain. */
    {
	    unsigned ntrees = lb_get_num_local_trees_in_group(lb, 0);
	    for(unsigned i = 0; i < ntrees; ++i) {
		hig_cell *tree = lb_get_local_tree_in_group(lb, i, 0);
		sd_add_higtree(sdp, tree);
	    }
    }

    /* Creating the partitioned sub-domain to simulation */
    psim_domain *psdp = psd_create(sdp, pg);

    higcit_celliterator *it;
    it = sd_get_domain_celliterator(sdp);
    fd = fopen("t.vtk", "w");
    higio_print_celliterator_in_vtk3d(fd, it);
    fclose(fd);

    /* Pressure boundary condition */
    computeBCP(psdp, total_numbcs, lb, pbctypes, pbcvalues, incremental, dt);

    ///* Print the pressure */
    char vtkname[1024];
    writePBC(psdp, vtkresults, "bc", myrank);

    /* Creating the distributed properties */
    psd_synced_mapper(psdp);

    /* Distributed property for pressure */
    distributed_property * dpp     = psd_create_property(psdp);
    distributed_property * ddeltap = psd_create_property(psdp);

    /* Creating the vector properties from facets */
    /* Facets list */
    sim_facet_domain  *sfdu[DIM];
    /* Property facet list */
    psim_facet_domain *psfdu[DIM];
    /* Distributed properties from facets in the domain */
    /* Final velocity */
    distributed_property *dpu[DIM];
    /* Intermediate velocity */
    distributed_property *dpustar[DIM];

    /* Linking the property to the facets */
    for(int dim = 0; dim < DIM; dim++) {

        /* Creating the list of cells center in the local domain */
        sfdu[dim] = sfd_create(NULL, dim);
        //sfd_use_cache(sfdu[dim], 0); //---CHECK!!!

        /* Setting the order 2 for the properties interpolation */
        sfd_set_interpolator_order(sfdu[dim], 3);   // ---CHECK ORDER 3

        /* Copying the list of cells center in the domain for pressure */
        sfd_copy_higtrees_from_center_domain(sfdu[dim], sdp);
        sfd_adjust_facet_ids(sfdu[dim]);

        /* Creating property for the facets */
        psfdu[dim] = psfd_create(sfdu[dim], psdp);
    }

    // Calculate the boundary conditions - velocity
    computeBCU(psfdu, total_numbcs, lb, ubctypes, ubcvalues, vbctypes, vbcvalues, wbctypes, wbcvalues);

    /* Parallel load balance context no longer needed. Destroy it. */
    lb_destroy(lb);

    for(int dim = 0; dim < DIM; dim++) {

        psfd_compute_sfbi(psfdu[dim]);
        /* Mapping the properties in the domain */
        psfd_synced_mapper(psfdu[dim]);

        /* Final velocity */
        dpu[dim] = psfd_create_property(psfdu[dim]);
        /* Intermediate velocity */
        dpustar[dim] = psfd_create_property(psfdu[dim]);
    }

    /* Intitializing the properties (velocity and pressure) */
    initializeP(psdp, dpp);
    initializeP(psdp, ddeltap); //verificar, para incremetal
    initializeU(psfdu, dpu);

    /* Creating the stencil for properties interpolation */
    sim_stencil *stn = stn_create();

    /* Initializing the number of frames for print */
    int numframe = 0;
    //  DEBUG_INSPECT(numprobes, %d);

    if (incremental) {
        //* Loop integration of the Navier-Stokes equations */
        for(int step = 0; step < numsteps; step++) {

            if (myrank == 0) printf("****** \n  Step: %d | Results at t = %g \n", step, (step+1)*dt);
            if (step == 0) { START_CLOCK(firstiter); }

            /* Update velocities and pressure using the projection method */
            updatePandU_incremental(psdp, psfdu, dpp, ddeltap, dpu, dpustar, step*dt, dt, Re, stn);

            if (step == 0) { STOP_CLOCK(firstiter); }
            if (myrank == 0) printf("****** \n");
            //DEBUG_INSPECT(step, %d);
            if (step % perstep == 0){
                /* Printing the results to standard output */
                //printResults(psdp, psfdu, numprobes, probe, prop, factor, dpp, dpu);

                /* Printing the velocities to visualize */
                //writeframe(psdp, psfdu[2], dpp, dpu[2], vtkresults, "w", myrank, numframe);
                hig_printVTK3dcell(psdp, dpp, dpu, psfdu, numhigs, vtkresults, "incr", myrank, numframe);

                numframe++;
            }
        }
    }
    else {
        for(int step = 0; step < numsteps; step++) {

            if (myrank == 0) printf("****** \n  Step: %d | Results at t = %g \n", step, (step+1)*dt);
            if (step == 0) { START_CLOCK(firstiter); }
            updatePandU(psdp, psfdu, dpp, dpu, dpustar, step*dt, dt, Re, stn);

            if (step == 0) { STOP_CLOCK(firstiter); }
            if (myrank == 0) printf("****** \n\n\n");
        }
    }

    //writevtk(psdp, myrank);
/*
	FILE *fdout = fopen("saida.vtk", "w");
	if (fdout == NULL)
	{
    	printf("Error opening file!\n");
    	exit(1);
	}
*/
    //Save files - pos processing
    //hig_printVTK3dcell(FILE *f, hig_cell *root, psim_domain *psdp, distributed_property *dpp, distributed_property  *dpu[DIM], distributed_property  *dpf[DIM], psim_facet_domain *psfdu[DIM])
	hig_printVTK3dcell(psdp, dpp, dpu, psfdu, numhigs, vtkresults, "final", myrank, numframe);

    //printResults(psdp, psfdu, dpp, dpu, numhigs);

	//fclose(fdout);

    for(int dim = 0; dim < DIM; dim++) {
        dp_destroy(dpu[dim]);
        dp_destroy(dpustar[dim]);
    }
    dp_destroy(dpp);

    stn_destroy(stn);

    STOP_CLOCK(total);
    /* Getting the time */
    if(myrank == 0)
    {
        DEBUG_INSPECT(GET_NSEC_CLOCK(total)/1.0e9, %lf);
        DEBUG_INSPECT(GET_NSEC_CLOCK(poisson)/1.0e9, %lf);
        DEBUG_INSPECT(GET_NSEC_CLOCK(syncprop)/1.0e9, %lf);
        DEBUG_INSPECT(GET_NSEC_CLOCK(firstiter)/1.0e9, %lf);
    }
}
