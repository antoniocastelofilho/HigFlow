/*********************************************************************************************/
/***    Including files: Begin                                                             ***/
/*********************************************************************************************/

#include "higtree.h"
#include "higtree-io.h"
#include "higtree-iterator.h"
#include "pdomain.h"
#include "lbal.h"
#include "solver.h"
#define DEBUG
#include "Debug-c.h"

/*********************************************************************************************/
/***    Including files: End                                                               ***/
/*********************************************************************************************/






/*********************************************************************************************/
/***    Functions Declaration: Begin                                                       ***/
/*********************************************************************************************/

real analyticU(real x, real y);
real analyticV(real x, real y);

real B(real x, real y, real Re);
real bodyforce(Point p, int dim, real Re);

sim_boundary *make_bc(hig_cell *bcg, int type);
void computeBCP(psim_domain *psd, higio_amr_info *mi, int rank);
void computeBCU(psim_facet_domain *psfdu[DIM], higio_amr_info *mi);

void initializeP(psim_domain *psdp, distributed_property *dpp);
void initializeU(psim_facet_domain *psfdu[DIM], distributed_property *dpu[DIM]);

void compute_facet_dudx_at_point(sim_facet_domain *sfdu, Point center, Point delta, int dim, real alpha, distributed_property *dpu, sim_stencil *stn, real valueatpoint, real *dudx);
void compute_dpdx_at_point(sim_domain *sdp, Point center, Point delta, int dim, real alpha, distributed_property *dpp, sim_stencil *stn, real *dpdx);
void compute_facet_du2dx2_at_point(sim_facet_domain *sfdu, Point center, Point delta, int dim, real alpha, distributed_property *dpu, sim_stencil *stn, real valueatpoint, real *du2dx2);

void compute_p(psim_domain *psdp, psim_facet_domain *psfdu[DIM], distributed_property *dpp, distributed_property *dpustar[DIM], real dt, real Re, sim_stencil *stn);

void compute_ustar(psim_facet_domain *psfdu[DIM], int dim, distributed_property *dpu[DIM], real dt, real Re, distributed_property *dpustar[DIM], sim_stencil *stn);

void compute_ut(psim_domain *psdp, psim_facet_domain *psfdu[DIM], int dim, distributed_property *dpp, distributed_property *dpustar[DIM], distributed_property *dpu[dim], real dt, real Re, sim_stencil *stn);

void updatePandU(psim_domain *psdp, psim_facet_domain *psfdu[DIM], distributed_property *dpp, distributed_property *dpu[DIM], distributed_property *dpustar[DIM], real dt, real Re, sim_stencil *stn);

void writePBC(psim_domain *psdp, const char *basefilename, const char *propertyname, int rank);
void writeframe(psim_domain *psdp, psim_facet_domain *psfdu, distributed_property *dpp, distributed_property *dpu, const char *basefilename, const char *propertyname, int rank, int frame);
void printResults(psim_domain *psdp, psim_facet_domain *psfdu[DIM], distributed_property *dpp, distributed_property *dpu[DIM]);

int main (int argc, char *argv[]);

/*********************************************************************************************/
/***    Functions Declaration: End                                                         ***/
/*********************************************************************************************/









/*********************************************************************************************/
/***    Analytical Solution: Begin                                                         ***/
/*********************************************************************************************/

real analyticU(real x, real y) {
        /* Analitical solution of cavity problem on the report page 5 */

        /* Velocity u */
	return 8.0*(x*x*(1.0+x*(-2.0+x)))*(y*(-2.0+4.0*y*y));
}

real analyticV(real x, real y) {
        /* Analitical solution of cavity problem on the report page 5 */

        /* Velocity v */
	return -8.0*(x*(2.0+x*(-6.0+4.0*x)))*(y*y*(-1.0+y*y));;
}

/*********************************************************************************************/
/***    Analytical Solution: End                                                           ***/
/*********************************************************************************************/











/*********************************************************************************************/
/***    External Forces: Begin                                                             ***/
/*********************************************************************************************/

/* Body Force */
real B(real x, real y, real Re) {
       /* Analitical solution of cavity problem on the report page 5 */

       /* External Force */
       real x2 = x * x;
       real x3 = x2 * x;
       real x4 = x3 * x;
       real x5 = x4 * x;
       real y2 = y * y;
       real y3 = y2 * y;
       real y4 = y3 * y;
       real fx = x4 - 2.0* x3 + x2;
       real fpx = 4.0*x3 - 6.0 *x2 + 2.0*x;
       real fppx = 12.0*x2 - 12.0*x + 2.0;
       real fpppx = 24.0*x - 12.0;
       real gy = y4 - y2;
       real gpy = 4.0*y3 - 2.0*y;
       real gppy = 12.0 *y2 - 2.0;
       real gpppy = 24.0*y;
       real Fx = x5/5.0 - x4/2.0 + x3/3.0;
       real F1x = fx*fppx - fpx*fpx;
       real F2x = fx*fx/2.0;
       real G1y = gy*gpppy - gpy*gppy;
       real res = -(8.0/Re)*(24.0*Fx + 2.0*fpx*gppy + fpppx*gy) - 64.0*(F2x*G1y - gy*gpy*F1x);
       return res;
}

real bodyforce(Point p, int dim, real Re) {
       return 0.0;
}

/*********************************************************************************************/
/***    External Forces: End                                                               ***/
/*********************************************************************************************/








/*********************************************************************************************/
/***    Boundary Conditions: Begin                                                         ***/
/*********************************************************************************************/

sim_boundary *make_bc(hig_cell *bcg, int type) {
        /* Creating boundaries conditions */

	/* Creating the mapper for the boundary */
	mp_mapper *bm = mp_create();

        /* Creating the iterator for the boundary condition */
	higcit_celliterator *it;
	it = higcit_create_all_leaves(bcg);

        /* Assign the mapper for the iterator */
	mp_assign_from_celliterator(bm, it, 0);

        /* Destroying the iterator */
	higcit_destroy(it);

        /* Creating the boundary condition for the mapp */
	sim_boundary *bc = sb_create(bcg, type, bm);

	return bc;
}

void computeBCP(psim_domain *psd, higio_amr_info *mi, int rank) {
	/* Computing the boundary contition for the pressure */

        /* Get the local domain */
	sim_domain *sd = psd_get_local_domain(psd);

        /* Creating the iterator and the grid cells for Dirichet boundary condition */
	higcit_celliterator *it;
	hig_cell *dirichletpoint = NULL;

        /* Reading and setting the boundary conditions */
	for(int dim = 0; dim < DIM; dim++) {
		for(int dir = 0; dir < 2; dir++) {

                        /* Reading the boundary in dim-axe and dir-direction */
			hig_cell *bcg = higio_read_bc_from_amr(mi, dim, dir);

                        /* Creating the Neumann Boundary Condition */
			int bc_t = NEUMANN;
			if (dim == 0 && dir == 1) {
				bc_t = DIRICHLET;
			}
			sim_boundary *bc = make_bc(bcg, bc_t);

                        /* Adding the boundary condition */
			sd_add_boundary(sd, bc);
		}
	}

        /* Getting the number of the Neumann boundaries conditions */
	int numbcs = sd_get_num_bcs(sd, NEUMANN);

        /* Traversing the Neumann boundaries conditions */
	for (int i = 0; i < numbcs; i++) {

                /* Getting the i Neumann boundary condition */
		sim_boundary *bc = sd_get_bc(sd, NEUMANN, i);

                /* Getting the mapper for this boundary condition */
		mp_mapper *bm = sb_get_mapper(bc);

                /* Taversing the cells of the boundary contition */
		for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {

                        /* Get the cell */
			hig_cell *bcell = higcit_getcell(it);

                        /* Get the center of the cell */
			Point bccenter;
			hig_get_center(bcell, bccenter);

                        /* Getting the identifier of this cell in the boundary condition object */
			int bcgid = mp_lookup(bm, hig_get_cid(bcell));

                        /* Setting the value 0.0 for this cell */
			if (bcgid >= 0) {
				sb_set_value(bc, bcgid, 0.0);
			}
		}

                /* Destroying the iterator */
		higcit_destroy(it);
	}

        /* Getting the number of the Dirichlet boundaries conditions */
	numbcs = sd_get_num_bcs(sd, DIRICHLET);

        /* Traversing the Dirichlet boundaries conditions */
	for (int i = 0; i < numbcs; i++) {

                /* Getting the i Dirichlet boundary condition */
		sim_boundary *bc = sd_get_bc(sd, DIRICHLET, i);

                /* Getting the mapper for this boundary condition */
		mp_mapper *bm = sb_get_mapper(bc);

                /* Taversing the cells of the boundary contition */
		for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {

                        /* Get the cell */
			hig_cell *bcell = higcit_getcell(it);

                        /* Get the center of the cell */
			Point bccenter;
			hig_get_center(bcell, bccenter);

                        /* Getting the identifier of this cell in the boundary condition object */
			int bcgid = mp_lookup(bm, hig_get_cid(bcell));

                        /* Setting the value 0.0 for this cell */
			if (bcgid >= 0) {
				sb_set_value(bc, bcgid, 0.0);
			}
		}

                /* Destroying the iterator */
		higcit_destroy(it);
	}
}

void computeBCU(psim_facet_domain *psfdu[DIM], higio_amr_info *mi) {
        /* Creating Boundaries Conditions for Velocity and Setting for the cavity problem */
        /* See the report page xxx */

        /* Creating the iterator and the grid cells for velocity boundary condition */
	higcit_celliterator *it;
	sim_facet_domain *sfdu[DIM];
	sim_boundary *bcinflowu = NULL;
	sim_boundary *bcoutflow[DIM];

        /* Reading and setting the boundary conditions */
	for(int dim = 0; dim < DIM; dim++) {

                /* Getting the local domain */
		sfdu[dim] = psfd_get_local_domain(psfdu[dim]);

                /* Reading the boundary condition for the velocity U(dim) in dim2-axe and dir-direction */
		for(int dim2 = 0; dim2 < DIM; dim2++) {
			for(int dir = 0; dir < 2; dir++) {

                                /* Reading the boundary conditiona from amr structure file */
				hig_cell *bcg = higio_read_bc_from_amr(mi, dim2, dir);

				int bc_t = DIRICHLET;
				if (dim2 == 0 && dir == 1) {
					bc_t = NEUMANN;
				}
				sim_boundary *bc = make_bc(bcg, bc_t);

                                /* Adding this boundary contition to the boundary condition object */
				sfd_add_boundary(sfdu[dim], bc);

				if (dim == 0 && dim2 == 0) {
					if (dir == 0) {
						bcinflowu = bc;
					} else {
						bcoutflow[0] = bc;
					}
				}
				if (dim == 1 && dim2 == 0 && dir == 1) {
					bcoutflow[1] = bc;
				}
			}
		}

                /* Local domain */
		sim_domain *sd = sfdu[dim]->cdom;

                /* Getting the number of Diriclet boundary conditions */
		int numbcs = sd_get_num_bcs(sd, DIRICHLET);

                /* Traversing the Diriclet boundaries conditions */
		for (int i = 0; i < numbcs; i++) {

                        /* Getting the i Driclet boundary condition */
			sim_boundary *bc = sd_get_bc(sd, DIRICHLET, i);

                        /* Mapping this boundary condition */
			mp_mapper *bm = sb_get_mapper(bc);

                        /* Traversing the cells of the Diriclet boundaries conditions */
			for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {

                                /* Getting the cell */
				hig_cell *bcell = higcit_getcell(it);

                                /* Set the value 0.0 for the center of the cell */
				Point bccenter;
				hig_get_center(bcell, bccenter);
				int bcgid = mp_lookup(bm, hig_get_cid(bcell));
				sb_set_value(bc, bcgid, 0.0);
			}

                        /* Destroing the iterator */
			higcit_destroy(it);
		}
	}

        /* Setting the Mapper for the boundary on the inflow */
	mp_mapper *bm = sb_get_mapper(bcinflowu);
	hig_cell *bcincell = sb_get_higtree(bcinflowu);
	Point bcinl, bcinh;
	hig_get_lowpoint(bcincell, bcinl);
	hig_get_highpoint(bcincell, bcinh);

        /* Traversing the cells of boundary on the inflow */
	for(it = sb_get_celliterator(bcinflowu); !higcit_isfinished(it); higcit_nextcell(it)) {

                /* Getting the cell */
		hig_cell *bcell = higcit_getcell(it);

                /* Calculating the velocity at the center of the cell */
		Point bccenter;
		hig_get_center(bcell, bccenter);
		real x = (bccenter[1] - bcinl[1]) / (bcinh[1] - bcinl[1]);

                /* Getting the cell identifier */
		int bcgid = mp_lookup(bm, hig_get_cid(bcell));

                /* Calculating the velocity */
		real v;
		if (x > 0.8) {
			v = 1.0;
		} else {
			v = 0.0;
		}

                /* Setting the velocity at boundary condition object */
		sb_set_value(bcinflowu, bcgid, v);
	}

        /* Destroying the iterator */
	higcit_destroy(it);

	for(int dim = 0; dim < DIM; dim++) {
		mp_mapper *bm = sb_get_mapper(bcinflowu);

		/* Traversing the cells of boundary on the outflow */
		for(it = sb_get_celliterator(bcinflowu); !higcit_isfinished(it); higcit_nextcell(it)) {

			/* Getting the cell */
			hig_cell *bcell = higcit_getcell(it);

			/* Calculating the velocity at the center of the cell */
			Point bccenter;
			hig_get_center(bcell, bccenter);

			/* Getting the cell identifier */
			int bcgid = mp_lookup(bm, hig_get_cid(bcell));

			/* Calculating the velocity */
			real v = 0.0;

			/* Setting the velocity at boundary condition object */
			sb_set_value(bcinflowu, bcgid, v);
		}

		/* Destroying the iterator */
		higcit_destroy(it);
	}

}

/*********************************************************************************************/
/***    Boundary Conditions: End                                                           ***/
/*********************************************************************************************/








/*********************************************************************************************/
/***    Initial Conditions: Begin                                                          ***/
/*********************************************************************************************/

void initializeP(psim_domain *psdp, distributed_property *dpp) {
        /* Setting the initial value for the pressure into the cavity */

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
        /* Setting the initial value for the pressure into the cavity */

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

                /* Calculating the velocity at the facet cell */
		Point l, h, delta;
		hig_cell *dom = sfd_get_higtree(sfdu[dim], 0);
		hig_get_lowpoint(dom, l);
		hig_get_highpoint(dom, h);
		hig_get_delta(dom, delta);
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

/*********************************************************************************************/
/***    Initial Conditions: End                                                            ***/
/*********************************************************************************************/












/*********************************************************************************************/
/***    Computing Properties Derivative at the Point: Begin                                ***/
/*********************************************************************************************/

void compute_facet_dudx_at_point(sim_facet_domain *sfdu, Point center, Point delta, int dim, real alpha, distributed_property *dpu, sim_stencil *stn, real valueatpoint, real *dudx) {
        /* Computing second order aproximation of first derivative of velocity at the point */

        /* Point to calculate the derivative */
	Point p;
	POINT_ASSIGN(p, center);

        /* Left point */
	p[dim] = center[dim] - alpha * delta[dim];
	real valuel = sfd_dp_interpolate(sfdu, dpu, center, p, stn);

        /* Right point */
	p[dim] = center[dim] + alpha * delta[dim];
	real valueh = sfd_dp_interpolate(sfdu, dpu, center, p, stn);

        /* First order derivative of velocity at the point p */
	*dudx = (valueh - valuel) / (2.0 * alpha * delta[dim]);
}

void compute_dpdx_at_point(sim_domain *sdp, Point center, Point delta, int dim, real alpha, distributed_property *dpp, sim_stencil *stn, real *dpdx) {
        /* Computing second order aproximation of first derivative of pressure at the point */

        /* Point to calculate the derivative */
	Point p;
	POINT_ASSIGN(p, center);

        /* Left point */
	p[dim] = center[dim] - alpha * delta[dim];
	real valuel = sd_dp_interpolate(sdp, dpp, center, p, stn);

        /* Right point */
	p[dim] = center[dim] + alpha * delta[dim];
	real valueh = sd_dp_interpolate(sdp, dpp, center, p, stn);

        /* First order derivative of pressure at the point p */
	*dpdx = (valueh - valuel) / (2.0 * alpha * delta[dim]);
}

void compute_facet_du2dx2_at_point(sim_facet_domain *sfdu, Point center, Point delta, int dim, real alpha, distributed_property *dpu, sim_stencil *stn, real valueatpoint, real *du2dx2) {
        /* Computing second order aproximation of second derivative of velocity at the point */

        /* Point to calculate the derivative */
	//DEBUG_INSPECT(valueatpoint, %lf);
	//DEBUG_INSPECT(alpha, %lf);
	//DEBUG_INSPECT(delta[0], %lf);
	//DEBUG_INSPECT(delta[1], %lf);
	Point p;
	POINT_ASSIGN(p, center);

        /* Left point */
	p[dim] = center[dim] - alpha * delta[dim];
	//DEBUG_INSPECT(p[0], %lf);
	//DEBUG_INSPECT(p[1], %lf);
	real valuel = sfd_dp_interpolate(sfdu, dpu, center, p, stn);
	//DEBUG_INSPECT(valuel, %lf);

        /* Right point */
	p[dim] = center[dim] + alpha * delta[dim];
	//DEBUG_INSPECT(p[0], %lf);
	//DEBUG_INSPECT(p[1], %lf);
	real valueh = sfd_dp_interpolate(sfdu, dpu, center, p, stn);

	//DEBUG_INSPECT(valueh, %lf);
        /* Second order derivative of pressure at the point p */
	*du2dx2 = (valuel - 2.0*valueatpoint + valueh)/(delta[dim] * delta[dim] * alpha * alpha);
	//DEBUG_INSPECT(*du2dx2, %lf);
}

/*********************************************************************************************/
/***    Computing Properties Derivative at the Point: End                                  ***/
/*********************************************************************************************/








/*********************************************************************************************/
/***    Computing Pressure Solving Poisson Equation: Begin                                 ***/
/*********************************************************************************************/

void compute_p(psim_domain *psdp, psim_facet_domain *psfdu[DIM], distributed_property *dpp, distributed_property *dpustar[DIM], real dt, real Re, sim_stencil *stn) {
        /* This function solve the Poisson equation to calculate the pressure */

        /* Initialing the solver */
	static solver *slvp = NULL;
	if (slvp == NULL) {

                /* Creating the local domain */

                /* Number of elements */
		int localdomainsize = psd_get_local_domain_size(psdp);

                /* Cranting the solver local domain */
		slvp = slv_create(SOLVER_ANY, psd_get_first_id(psdp), localdomainsize);

                /* Setting the maximum number of non zeros values */
		slv_set_maxnonzeros(slvp, 300);
	}

        /* Declaring the cell iterator */
	higcit_celliterator *it;

        /* Getting the local domain */
	sim_domain *sdp = psd_get_local_domain(psdp);
	sim_facet_domain *sfdu[DIM];

        /* Getting the velocty domain */
	for(int dim = 0; dim < DIM; dim++) {
		sfdu[dim] = psfd_get_local_domain(psfdu[dim]);
	}

        /* Mapping the local domain for the cell center */
	mp_mapper *mp = sd_get_domain_mapper(sdp);

        /* Traversing the cells */
	for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {

                /* Getting the cell */
		hig_cell *c = higcit_getcell(it);

                /* Getting the cell center and cell size */
		Point ccenter;
		Point cdelta;
		hig_get_center(c, ccenter);
		hig_get_delta(c, cdelta);

                /* Calculating the velocity divergent at the center */
		real sumdudx = 0.0;
		for(int dim = 0; dim < DIM; dim++) {
                        /* First order derivative of the velocity */
			real dudx;
			compute_facet_dudx_at_point(sfdu[dim], ccenter, cdelta, dim, 0.5, dpustar[dim], stn, 0.0, &dudx);
			sumdudx += dudx;
		}

                /* Reset the stencil */
		stn_reset(stn);

                /* Setting the right side of the equation */
		stn_set_rhs(stn, sumdudx / dt);

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

                /* Getting the identifiers of real points */
		int *ids = psd_stn_get_gids(psdp, stn);

                /* Getting the coeficients of real points */
		real *vals = stn_get_vals(stn);

                /* Getting the global identifier */
		int cgid = psd_get_global_id(psdp, c);

                /* Setting the matrix line for this point */
		slv_set_Ai(slvp, cgid, numelems, ids, vals);

                /* Setting the independet vector for this point */
		slv_set_bi(slvp, cgid, stn_get_rhs(stn));
	}

        /* Destroying the iterator */
	higcit_destroy(it);

        /* Assembling the local matriz to the global matriz */
	slv_assemble(slvp);

        /* Solving the linear system */
	slv_solve(slvp);

        /* Getting the solution of the linear system */
	dp_slv_load_from_solver(dpp, slvp);
}

/*********************************************************************************************/
/***    Computing Pressure Solving Poisson Equation: End                                   ***/
/*********************************************************************************************/





/*********************************************************************************************/
/***    Computing the Intermediate Velocity by Projection Method: Begin                    ***/
/*********************************************************************************************/

void compute_ustar(psim_facet_domain *psfdu[DIM], int dim, distributed_property *dpu[DIM], real dt, real Re, distributed_property *dpustar[DIM], sim_stencil *stn) {
        /* Computing the intermediate veloicities using the first order projection method  */

        /* Mapping the domain and properties */
	higfit_facetiterator *fit;
	sim_facet_domain *sfdu[DIM];
	for(int dim = 0; dim < DIM; dim++) {
		sfdu[dim] = psfd_get_local_domain(psfdu[dim]);
	}
	mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);

        /* iterating on the facets for compute the intermediate velocity */
	for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
                /* Getting the identifiers and data at the facets */
		hig_facet *f = higfit_getfacet(fit);
		int fgid = mp_lookup(mu, hig_get_fid(f));
		Point fcenter;
		Point fdelta;
		hig_get_facet_center(f, fcenter);
		hig_get_facet_delta(f, fdelta);

                /* Getting the value of the velocity at center of facet */
		real valuec = dp_get_value(dpu[dim], fgid);

                /* Initializing the right side of the Navier-Equation */
		real tot = 0.0;

                /* Computing the derivatives */
		for (int dim2 = 0; dim2 < DIM; dim2++) {

                        /* Computing the first order derivative */
			real dudx = 0.0;
			compute_facet_dudx_at_point(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], stn, valuec, &dudx);

                        /* Computing the first order derivative */
			real du2dx2 = 0.0;
			compute_facet_du2dx2_at_point(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], stn, valuec, &du2dx2);

                        /* Computing the value of the velocity at the center of the facet */
                        /* for the convective term                                        */
			real u;
			Point fakecenter;
			POINT_ASSIGN(fakecenter, fcenter); fakecenter[0] += fdelta[0];
			u = sfd_dp_interpolate(sfdu[dim2], dpu[dim2], fakecenter, fcenter, stn);
                        /* Convective term */
			tot += -u*dudx;

                        /* Diffusive term */
			tot += du2dx2/Re;
		}

                /* External force */
                tot += - bodyforce(fcenter, dim, Re);

                /* Calculating the intermediate velocity */
		real ustar = valuec + dt * tot;

                /* Storing the value of the intermediate velocity */
		dp_set_value(dpustar[dim], fgid, ustar);
	}

        /* Destroying the iterator */
	higfit_destroy(fit);

        /* Sincronizing the intermediate velocity for the domain */
	dp_sync(dpustar[dim]);
}

/*********************************************************************************************/
/***    Computing the Intermediate Velocity by Projection Method: Begin                    ***/
/*********************************************************************************************/








/*********************************************************************************************/
/***    Computing the Final Velocity by the Projection Method: Begin                       ***/
/*********************************************************************************************/

void compute_ut(psim_domain *psdp, psim_facet_domain *psfdu[DIM], int dim, distributed_property *dpp, distributed_property *dpustar[DIM], distributed_property *dpu[dim], real dt, real Re, sim_stencil *stn) {
        /* Computing the final veloicities and pressure using the first order projection method  */

        /* Mapping the domain and properties */
	sim_domain *sdp = psd_get_local_domain(psdp);
	sim_facet_domain *sfdu[DIM];
	higfit_facetiterator *fit;
	sfdu[dim] = psfd_get_local_domain(psfdu[dim]);
	mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);

        /* iterating on the facets for update the final velocity */
	for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {

                /* Getting the identifiers and data at the facets */
		hig_facet *f = higfit_getfacet(fit);
		int fgid = mp_lookup(mu, hig_get_fid(f));
		Point fcenter;
		Point fdelta;
		hig_get_facet_center(f, fcenter);
		hig_get_facet_delta(f, fdelta);

                /* Computing the pressure gradient */
		real dpdx;
		compute_dpdx_at_point(sdp, fcenter, fdelta, dim, 0.5, dpp, stn, &dpdx);

                /* Getting the intermetiate velocity */
		real ustar = dp_get_value(dpustar[dim], fgid);

                /* Computing the final velocity */
		real utdt = ustar - dt * dpdx;

                /* Storing the final velocity */
		dp_set_value(dpu[dim], fgid, utdt);
	}

        /* destroying the iterator */
	higfit_destroy(fit);

        /* Synchronizing distributed proprieties */
	dp_sync(dpu[dim]);
}

/*********************************************************************************************/
/***    Computing the Final Velocity by the Projection Method: End                         ***/
/*********************************************************************************************/








/*********************************************************************************************/
/***    Solve the Incomnpressible Navier-Stokes using the Projetion Method: Begin          ***/
/*********************************************************************************************/

void updatePandU(psim_domain *psdp, psim_facet_domain *psfdu[DIM], distributed_property *dpp, distributed_property *dpu[DIM], distributed_property *dpustar[DIM], real dt, real Re, sim_stencil *stn) {
        /* Computing the Velocities and Pressure using the first order projection method  */

        /* Facet iterator */
	higfit_facetiterator *fit;

        /* Computing the intermediate velocity from projection method */
        /* See report on the projection method - page 1               */
	for(int dim = 0; dim < DIM; dim++) {
		compute_ustar(psfdu, dim, dpu, dt, Re, dpustar, stn);
	}

        /* Computing the pressure by Poisson Equation solver */
        /* See report on the projection method - page 2      */
	compute_p(psdp, psfdu, dpp, dpustar, dt, Re, stn);

        /* Computing the final velocity from the projection method */
        /* See report on the projection method - page 3            */
	for(int dim = 0; dim < DIM; dim++) {
		compute_ut(psdp, psfdu, dim, dpp, dpustar, dpu, dt, Re, stn);
	}
}

/*********************************************************************************************/
/***    Solve the Incomnpressible Navier-Stokes using the Projetion Method: Begin          ***/
/*********************************************************************************************/







/*********************************************************************************************/
/***    Printing proprieties for visualizations: Begin                                     ***/
/*********************************************************************************************/

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
	higio_print_celliterator_in_vtk2d(fdout, it);

        /* Destroying the iterator */
	higcit_destroy(it);

        /* Close the file */
	fclose(fdout);
}


void writeframe(psim_domain *psdp, psim_facet_domain *psfdu, distributed_property *dpp, distributed_property *dpu, const char *basefilename, const char *propertyname, int rank, int frame) {
        /* This function print the properties in VTK format  */

        /* Getting the domain */
	sim_domain *sdp = psd_get_local_domain(psdp);

        /* Creating the stencil */
	sim_stencil *stn = stn_create();

        /* Open the VTK file */
	char vtkname[1024];
	sprintf(vtkname, "%s_%s_%d-%d.vtk", basefilename, propertyname, rank, frame);
	FILE *fdout = fopen(vtkname, "w");
	if (fdout == NULL) {
		return;
	}

        /* Getting the local domain */
	Point l, h;
	hig_cell *rootp = sd_get_higtree(sdp, 0);
	sim_facet_domain *sfdu;
	sfdu = psfd_get_local_domain(psfdu);

        /* Getting the sub-domain iterator */
	higcit_celliterator *it;
	it = sd_get_domain_celliterator(sdp);

        /* Printing the grid data */
	higio_print_celliterator_in_vtk2d(fdout, it);

        /* Destroying the iterator */
	higcit_destroy(it);

        /* Getting the cell center iterator */
	it = sd_get_domain_celliterator(sdp);

        /* Getting the number of cells to print */
	int numcells = higcit_count(it);

        /* Destroying the iterator */
	higcit_destroy(it);

        /* Printing the head of file */
        fprintf(fdout, "CELL_DATA %d\nSCALARS %s float 1\nLOOKUP_TABLE default\n", numcells, propertyname);

        /* Traversing the cell */
	for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {

                /* Getting the cell */
		hig_cell *c = higcit_getcell(it);

                /* Calculating the value in the center of the facet */
		Point ccenter;
		Point cdelta;
		hig_get_center(c, ccenter);
		hig_get_delta(c, cdelta);
		POINT_DIV_SCALAR(cdelta, cdelta, 2.0);
		real value;

                /* Computing the velocity at the facet */
		Point fakecenter;
		POINT_ASSIGN(fakecenter, ccenter); fakecenter[0] += cdelta[0];
		value = sfd_dp_interpolate(sfdu, dpu, fakecenter, ccenter, stn);

                /* Printing the value */
		fprintf(fdout, "%f\n", value);
	}

        /* Destroying the iterator */
	higcit_destroy(it);

        /* Destroying the stencil */
	stn_destroy(stn);

        /* Close the file */
	fclose(fdout);
}

/*********************************************************************************************/
/***    Printing proprieties for visualizations: End                                       ***/
/*********************************************************************************************/







/*********************************************************************************************/
/***    Printing velocytie error: Begin                                                    ***/
/*********************************************************************************************/

void printResults(psim_domain *psdp, psim_facet_domain *psfdu[DIM], distributed_property *dpp, distributed_property *dpu[DIM]) {
        /* This function print the properties in some points at the standad output */

        /* Getting the local domain */
	sim_domain *sdp = psd_get_local_domain(psdp);

        /* Creating the stencil */
	sim_stencil *stn = stn_create();


        /* Getting the local domain */
	int hig = 0;
	Point l, h;
	sim_facet_domain *sfdu[DIM];
	for(int dim = 0; dim < DIM; dim++) {
		sfdu[dim] = psfd_get_local_domain(psfdu[dim]);
	}
        /* Getting the bounding box */

        /* traversing the cells in the bounding box */
	higcit_celliterator *it;
	for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {

                /* Getting the cell */
		hig_cell *c = higcit_getcell(it);

                /* Getting the center and size */
		Point ccenter;
		Point cdelta;
		hig_get_center(c, ccenter);
		hig_get_delta(c, cdelta);

                /* Setting the center for interpolation and compute the value of the property */
		real value[DIM];
		for(int dim = 0; dim < DIM; dim++) {
			Point fakecenter;
			POINT_ASSIGN(fakecenter, ccenter); fakecenter[0] += cdelta[0];
			value[dim] = sfd_dp_interpolate(sfdu[dim], dpu[dim], fakecenter, ccenter, stn);
		}

                /* Printing the value of the velocity */
		printf("v(");

                /* Printing the point */
		for(int dim = 0; dim < DIM; dim++) {
			printf("%f ", ccenter[dim]);
		}
		printf(") = ");

                /* Printing the value */
		for(int dim = 0; dim < DIM; dim++) {
			printf("%20.15f ", value[dim]);
		}

                /* Printing the error */
		for(int i = 0; i < value[0]*20; i++) {
			printf("*");
		}
		printf("\n");
	}

        /* Destroying the stencil */
	stn_destroy(stn);
}

/*********************************************************************************************/
/***    Printing velocytie error: End                                                      ***/
/*********************************************************************************************/








/*********************************************************************************************/
/***    Main Program: Begin                                                                ***/
/*********************************************************************************************/

int main (int argc, char *argv[]) {
        /* Initializing Higtree */
	higtree_initialize(&argc, &argv);

	int myrank;
	int ntasks;
	int mpierr = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

        /* Seeing the success or not of initializing */
	//DEBUG_INSPECT(mpierr, %d);
	//DEBUG_INSPECT(MPI_SUCCESS, %d);
	//DEBUG_INSPECT(myrank, %d);

        /* Partitioning the domain */
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
	//DEBUG_INSPECT(myrank, %d);

        /* Taking the sub-domain of the myrank process */
	partition_graph *pg = pg_create(MPI_COMM_WORLD);

        /* Setting the fringe size of the sub-domain */
	pg_set_fringe_size(pg, 4);

        /* Getting the time */
	DEBUG_DIFF_TIME;

        /* Adapting the data input by the command line */
	char *amrfilename = argv[1];
	real dt = atof(argv[2]);
	real Re = atof(argv[3]);
	int numsteps = atoi(argv[4]);
	int perstep = atoi(argv[5]);

        /* Open the AMR format file */
	FILE *fd = fopen(amrfilename, "r");

        /* Reanding the IO information of the file */
	higio_amr_info *mi = higio_read_amr_info(fd);

        /* Closing the AMR format file */
	fclose(fd);

        /* Partitioning the grid from AMR information */
	load_balancer *lb = higio_partition_from_amr_info(pg, mi);

        /* Getting the time */
	DEBUG_DIFF_TIME

        /* Creating the simulation domain */
	sim_domain *sdp = sd_create(NULL);

        /* Setting the order of the polynomial interpolation */
	sd_set_interpolator_order(sdp, 2);

	/* Getting the numbers of trees assigned to this process. */
	unsigned numtrees = lb_get_num_local_trees(lb);

        /* Linking the HigTree grid to the simulation domain */
	for(unsigned i = 0; i < numtrees; ++i) {
		sd_add_higtree(sdp, lb_get_local_tree(lb, i, NULL));
	}

        /* Creating the partitioned sub-domain to simulation */
	psim_domain *psdp = psd_create(sdp, pg);

        /* Pressure boundary condition */
	computeBCP(psdp, mi, myrank);

        /* Creating the distributed properties */
	psd_synced_mapper(psdp);
	distributed_property * dpp = psd_create_property(psdp);

        /* Creating the vector properties from facets */

        /* Facets list */
	sim_facet_domain *sfdu[DIM];

        /* Property facet list */
	psim_facet_domain *psfdu[DIM];

        /* Distributed properties from facets in the domain */
        /* Final velocity */
	distributed_property * dpu[DIM];
        /* Intermediate velocity */
	distributed_property *dpustar[DIM];

        /* Linking the property to the facets */
	for(int dim = 0; dim < DIM; dim++) {

                /* Creating the list of cells center in the local domain */
		sfdu[dim] = sfd_create(NULL, dim);

                /* Setting the order 2 for the properties interpolation */
		sfd_set_interpolator_order(sfdu[dim], 2);

                /* Copying the list of cells center in the domain for pressure */
		sfd_copy_higtrees_from_center_domain(sfdu[dim], sdp);

		sfd_adjust_facet_ids(sfdu[dim]);

                /* Creating property for the facets */
		psfdu[dim] = psfd_create(sfdu[dim], psdp);


	}

        /* Getting the time */
	DEBUG_DIFF_TIME

        /* Computing the boundary conditions */

        /* Print the pressure */
	writePBC(psdp, "vtks/bc", "p", myrank);

        /* Velocity boundary condition */
	computeBCU(psfdu, mi);

	for(int dim = 0; dim < DIM; dim++) {
		psfd_compute_sfbi(psfdu[dim]);
                /* Mapping the properties in the domain */
		psfd_synced_mapper(psfdu[dim]);

                /* Final velocity */
		dpu[dim] = psfd_create_property(psfdu[dim]);

                /* Intermediate velocity */
		dpustar[dim] = psfd_create_property(psfdu[dim]);
	}

	sim_facet_domain *sfd = psfd_get_local_domain(psfdu[0]);
	mp_mapper *m = sfd_get_domain_mapper(sfd);
	higfit_facetiterator *it;
	for(it = sfd_get_domain_facetiterator(sfd); !higfit_isfinished(it); higfit_nextfacet(it)) {
		hig_facet *f = higfit_getfacet(it);
		Point center;
		hig_get_facet_center(f, center);
		int gid = hig_get_fid(f);
		int cid = mp_lookup(m, gid);
	//	printf("r = %d (%lf, %lf) = %d\n", myrank, center[0], center[1], cid);

	}

        /* Intitializing the properties (velocity and pressure) */
	initializeP(psdp, dpp);
	initializeU(psfdu, dpu);

        /* Setting the number of steps */

        /* Creating the stencil for properties interpolation */
	sim_stencil *stn = stn_create();

        /* Initializing the number of frames for print */
	int numframe = 0;

        /* Loop integration of the Navier-Stokes equations */
	for(int step = 0; step < numsteps; step++) {

                /* Update velocities and pressure using the projection method */
		updatePandU(psdp, psfdu, dpp, dpu, dpustar, dt, Re, stn);

                /* Printing the results */
		if (step % perstep == 0) {

                        /* Getting the time */
			DEBUG_DIFF_TIME;

                        /* Printing the results to standard output */
			printResults(psdp, psfdu, dpp, dpu);

                        /* Printing the velocities to visualize */
			writeframe(psdp, psfdu[0], dpp, dpu[0], "vtks/mesh", "u", myrank, numframe);
			writeframe(psdp, psfdu[1], dpp, dpu[1], "vtks/mesh", "v", myrank, numframe);

                        /* Increasing the number of prints */
			numframe++;
		}
	}

        /* Destroing the properties objects */
	for(int dim = 0; dim < DIM; dim++) {
		dp_destroy(dpustar[dim]);
	}

        /* Destroing the stencil for interpolation */
	stn_destroy(stn);

        /* Getting the time */
	DEBUG_DIFF_TIME;

	return 0;
}

/*********************************************************************************************/
/***    Main Program: End                                                                  ***/
/*********************************************************************************************/

