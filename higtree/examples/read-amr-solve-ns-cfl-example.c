/*********************************************************************************************/
/***    Including files: Begin                                                             ***/
/*********************************************************************************************/

#include <stdio.h>
#include <stdarg.h>
#include <float.h>
#include <string.h>
#include "higtree.h"
#include "higtree-io.h"
#include "higtree-iterator.h"
#include "global-stencil.h"
#include "lbal.h"
#include "solver-petsc.h"
#include "Debug-c.h"

/*********************************************************************************************/
/***    Including files: End                                                               ***/
/*********************************************************************************************/

/*********************************************************************************************/
/***    External Forces: Begin                                                             ***/
/*********************************************************************************************/

real bodyforce(Point p, int dim, real Re)
{
	return 0.0;
}

/*********************************************************************************************/
/***    External Forces: End                                                               ***/
/*********************************************************************************************/

static void
_adjust_outflow(psim_facet_domain * psfdu, distributed_property * dpu,
	real alpha)
{
	sim_facet_domain *sfd = psfd_get_local_domain(psfdu);
	int dim = sfd_get_dim(sfd);
	int numnbcs = sfd_get_num_bcs(sfd, NEUMANN);
	int numhigs = sfd_get_num_higtrees(sfd);
	for (int i = 0; i < numnbcs; i++) {
		sim_boundary *sb = sfd_get_bc(sfd, NEUMANN, i);
		hig_cell *bchig = sb_get_higtree(sb);
		Point bcl, bch;
		hig_get_lowpoint(bchig, bcl);
		hig_get_highpoint(bchig, bch);
		int bcdim = hig_get_narrowest_dim(bchig);
		POINT_SUB_SCALAR(bcl, bcl, EPSDELTA);
		POINT_ADD_SCALAR(bch, bch, EPSDELTA);
		for (int j = 0; j < numhigs; j++) {
			hig_cell *root = sfd_get_higtree(sfd, j);
			mp_mapper *m = sfd_get_domain_mapper(sfd);
			Point delta;
			Point center;
			higcit_celliterator *it =
			    higcit_create_bounding_box(root, bcl, bch);
			if (higcit_isfinished(it)) {
				higcit_destroy(it);
				continue;
			}
			hig_cell *c = higcit_getcell(it);
			hig_get_delta(c, delta);
			hig_get_center(c, center);
			higcit_destroy(it);
			Point backwarddelta;
			POINT_ASSIGN_SCALAR(backwarddelta, 0.0);
			Point flowl, flowh;
			POINT_ASSIGN(flowl, bcl);
			POINT_ASSIGN(flowh, bch);
			if (center[bcdim] < bcl[bcdim]) {
				backwarddelta[bcdim] = -alpha * delta[bcdim];
				flowl[bcdim] -= alpha * delta[bcdim];
			} else {
				backwarddelta[bcdim] = alpha * delta[bcdim];
				flowh[bcdim] += alpha * delta[bcdim];
			}
			int dimofinterest[DIM];
			POINT_ASSIGN_SCALAR(dimofinterest, 0);
			dimofinterest[dim] = 1;
			higfit_facetiterator *fit;
			for (fit =
			     higfit_create_bounding_box_facets(root,
							       dimofinterest,
							       flowl, flowh);
			     !higfit_isfinished(fit); higfit_nextfacet(fit)) {
				hig_facet *f = higfit_getfacet(fit);
				int fgid = mp_lookup(m, hig_get_fid(f));
				Point fcenter;
				hig_get_facet_center(f, fcenter);

				hig_facet fc;
				Point fcopy;
				POINT_ADD(fcopy, fcenter, backwarddelta);

				hig_get_facet_with_point(root, dim, fcopy, &fc);
				int fcgid = mp_lookup(m, hig_get_fid(&fc));
				if (fgid >= 0 && fcgid >= 0) {

					dp_set_value(dpu, fgid,
						     dp_get_value(dpu, fcgid));
				}
			}
		}
	}
}

static real
calc_diffusive_cfl(real max_cfl, sim_domain *sd, real Re)
{
	real min_dx2 = DBL_MAX;

	higcit_celliterator *it;
	for (it = sd_get_domain_celliterator(sd);
	     !higcit_isfinished(it); higcit_nextcell(it))
	{
		/* Getting the cell */
		hig_cell *c = higcit_getcell(it);

		/* Getting the center and size */
		Point cdelta;
		hig_get_delta(c, cdelta);

		/* Compute the value of the property */
		for (int dim = 0; dim < DIM; dim++) {
			real dx2 = cdelta[dim] * cdelta[dim];
			if(dx2 < min_dx2) {
				min_dx2 = dx2;
			}
		}
	}
	higcit_destroy(it);

	real ret;
	MPI_Allreduce(&min_dx2, &ret, 1, MPI_HIGREAL, MPI_MIN, MPI_COMM_WORLD);
	return max_cfl * ret * Re;
}

static void
calc_advective_cfl_for_facet(real *min_dt, hig_facet *f,
	sim_facet_domain *sfdu, distributed_property *dpu)
{
	/* Getting the center and size */
	Point center;
	Point delta;
	hig_get_facet_center(f, center);
	hig_get_facet_delta(f, delta);
	int dim = hig_get_facet_dim(f);

	const int fid = sfd_get_local_id(sfdu, f);
	const real dt = delta[dim] / fabs(
		dp_get_value(dpu, fid)
	);

	if(dt < *min_dt) {
		*min_dt = dt;
	}
}

static real
calc_advective_cfl(real max_cfl, sim_facet_domain **sfdu, distributed_property **dpu)
{
	real min_dt = DBL_MAX;

	for(unsigned dim = 0; dim < DIM; ++dim) {
		higfit_facetiterator *it;
		for (it = sfd_get_domain_facetiterator(sfdu[dim]);
	     		!higfit_isfinished(it); higfit_nextfacet(it))
		{
			/* Getting the cell */
			hig_facet *f = higfit_getfacet(it);

			calc_advective_cfl_for_facet(&min_dt, f, sfdu[dim], dpu[dim]);
		}
		higfit_destroy(it);
	}

	real ret;
	MPI_Allreduce(&min_dt, &ret, 1, MPI_HIGREAL, MPI_MIN, MPI_COMM_WORLD);
	return max_cfl * ret;
}

struct sim_state {
	real Re;
	real cfl;
	real dt;

	psim_domain *psdp;
       	psim_facet_domain *psfdu[DIM];

	distributed_property *dpp;
	distributed_property *dpu[DIM];
	distributed_property *dpustar[DIM];
};
typedef struct sim_state sim_state;

typedef void (*ns_solver_func)(sim_state *ss, sim_stencil * stn);



/*********************************************************************************************/
/***    Framework for Monolytic Solvers: Begin                                             ***/
/*********************************************************************************************/

struct monolytic_system_impl {
	real (*set_mass_eq)(hig_cell *c, int lid, sim_stencil *stn,
		sim_stencil *stn_u[DIM], sim_facet_domain *sfdu[DIM],
		distributed_property *dpu[DIM]);

	real (*set_ns_eq)(hig_facet *f, int flid, int dim, real dt, real Re,
				sim_stencil *stn[2], sim_stencil *stn_p,
				sim_stencil *stn_u[DIM],
				sim_facet_domain *sfdu[DIM],
				distributed_property *dpu[DIM],
				distributed_property *dpustar[DIM],
				sim_domain *sdp, distributed_property *dpp);

	bool (*func_norm_cond)(real mass_norm, real ns_norm);
	bool (*diff_norm_cond)(real diff_norm);
	real (*update_solution)(solver *slv, sim_state *ss, real *dt_out);
};
typedef struct monolytic_system_impl monolytic_system_impl;

static void set_composed_row(solver *slv, int gid,
	psim_domain *psdp, psim_facet_domain *psfdu[DIM],
	sim_stencil *stn_p, sim_stencil *stn_u[DIM])
{
	// Calculate total size
	const size_t p_numelems = stn_get_numelems(stn_p);
	size_t u_numelems[DIM];
	int numelems = p_numelems;

	for (unsigned dim = 0; dim < DIM; ++dim) {
		u_numelems[dim] = stn_get_numelems(stn_u[dim]);
		numelems += u_numelems[dim];
	}

	// Assemble all the properties in a single row
	int ids[numelems];
	real vals[numelems];

	memcpy(ids, psd_stn_get_gids(psdp, stn_p),
		p_numelems * sizeof ids[0]);
	memcpy(vals, stn_get_vals(stn_p),
		p_numelems * sizeof vals[0]);
	size_t start = p_numelems;
	for (unsigned dim = 0; dim < DIM; ++dim) {
		const int *gids = psfd_stn_get_gids(
			psfdu[dim], stn_u[dim]
		);
		memcpy(&ids[start], gids,
			u_numelems[dim] * sizeof ids[0]);
		memcpy(&vals[start],
			stn_get_vals(stn_u[dim]),
			u_numelems[dim] * sizeof vals[0]);

		start += u_numelems[dim];
	}
	assert(start == numelems);

	slv_set_Ai(slv, gid, numelems, ids, vals);
}

static void copy_facet_prop(distributed_property *to[DIM],
		distributed_property *from[DIM])
{
	for(unsigned dim = 0; dim < DIM; ++dim)
	{
		dp_copy_values(to[dim], from[dim]);
	}
}

static real set_mass_eqs(solver *slv, const monolytic_system_impl *impl,
	sim_state *ss, global_stencil *p_eqs[], sim_stencil *stn,
	sim_stencil *stn_u[DIM])
{
	higcit_celliterator *it;

	// Take the local domains.
	sim_facet_domain *sfdu[DIM];
	for (unsigned dim = 0; dim < DIM; ++dim) {
		sfdu[dim] = psfd_get_local_domain(ss->psfdu[dim]);
	}
	sim_domain *sdp = psd_get_local_domain(ss->psdp);

	real inf_norm = 0;

	mp_mapper *mp = sd_get_domain_mapper(sdp);
	for(it = sd_get_domain_celliterator(sdp);
		!higcit_isfinished(it); higcit_nextcell(it))
	{
		hig_cell *c = higcit_getcell(it);
		int lid = mp_lookup(mp, hig_get_cid(c));

		for (unsigned dim2 = 0; dim2 < DIM; ++dim2) {
			stn_reset(stn_u[dim2]);
		}
		const real mass_rhs = impl->set_mass_eq(c, lid, stn, stn_u, sfdu, ss->dpu);

		const int gid = psd_lid_to_gid(ss->psdp, lid);

		if(p_eqs) {
			// Gets the partialy assembled pressure equation.
			global_stencil *row = p_eqs[lid];

			// Add zero to the main diagonal so that the entry exists.
			gs_add_to_element(row, gid, 0.0);

			// Finish pressure equation by adding continuity
			// equation to it.
			const real scale = 1.0;
			for(unsigned dim = 0; dim < DIM; ++dim) {
				gs_add_facet_stencil(row, ss->psfdu[dim], stn_u[dim], scale);
			}
			gs_add_to_rhs(row, scale * mass_rhs);

			// Sets the complete pressure equation in the linear system.
			const real rhs = gs_get_rhs(row);

			// Set linear system righ-hand side
			slv_set_bi(slv, gid, rhs);

			// Set matrix row
			slv_set_Ai(slv, gid,
				gs_get_numelems(row),
				gs_get_ids(row),
				gs_get_vals(row)
			);

			// Reset for use in next iteration
			gs_reset(row);
		} else {
			/* Use empty stn as a placeholder for pressure stencil. */
			stn_reset(stn);

			/* Set linear system matrix row. */
			set_composed_row(slv, gid, ss->psdp, ss->psfdu, stn, stn_u);

			/* Set linear system right hand side. */
			slv_set_bi(slv, gid, mass_rhs);
		}

		const real mag = fabs(mass_rhs);
		if(mag > inf_norm) {
			inf_norm = mag;
		}
	}
	higcit_destroy(it);

	return inf_norm;
}

static void add_to_p_equations(global_stencil *p_eqs[],
	psim_domain *psdp, psim_facet_domain *psfdu[DIM],
	sim_stencil *stn_p, sim_stencil *stn_u[DIM], real rhs)
{
	unsigned size = stn_get_numelems(stn_p);
	const int *lcids = stn_get_ids(stn_p);
	const real *lvals = stn_get_vals(stn_p);
	unsigned lsize = psd_get_local_domain_size(psdp);

	for(unsigned i = 0; i < size; ++i) {
		const int lcid = lcids[i];
		if(lcid >= lsize) {
			continue;
		}

		const real scale = (lvals[i] < 0) ? -1.0 : 1.0;

		global_stencil *row = p_eqs[lcids[i]];

		gs_add_cell_stencil(row, psdp, stn_p, scale);
		for(unsigned dim = 0; dim < DIM; ++dim) {
			gs_add_facet_stencil(row, psfdu[dim], stn_u[dim], scale);
		}
		gs_add_to_rhs(row, scale * rhs);
	}
}

static real set_ns_eqs(solver *slv, const monolytic_system_impl *impl,
	sim_state *ss, global_stencil *p_eqs[], sim_stencil *stn1,
	sim_stencil *stn_p, sim_stencil *stn_u[DIM])
{
	higfit_facetiterator *fit;
	static sim_stencil *stn[2] = {NULL};

	// Anoter stencil is needed here, create it.
	if(!stn[0]) {
		stn[0] = stn_create();
	}
	stn[1] = stn1;

	// Take the local domains.
	sim_facet_domain *sfdu[DIM];
	for (unsigned dim = 0; dim < DIM; ++dim) {
		sfdu[dim] = psfd_get_local_domain(ss->psfdu[dim]);
	}
	sim_domain *sdp = psd_get_local_domain(ss->psdp);

	real inf_norm = 0;

	for(unsigned dim = 0; dim < DIM; ++dim) {
		mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
		for (fit = sfd_get_domain_facetiterator(sfdu[dim]);
		     !higfit_isfinished(fit); higfit_nextfacet(fit))
		{
			/* Getting the identifiers and data at the facets */
			hig_facet *f = higfit_getfacet(fit);
			int flid = mp_lookup(mu, hig_get_fid(f));

			stn_reset(stn_p);
			for (unsigned dim2 = 0; dim2 < DIM; ++dim2) {
				stn_reset(stn_u[dim2]);
			}

			const real rhs = impl->set_ns_eq(f, flid, dim, ss->dt,
				ss->Re, stn, stn_p, stn_u, sfdu, ss->dpu,
				ss->dpustar, sdp, ss->dpp);

			if(p_eqs) {
				// Adds momentum equations to pressure equations in
				// order to make the linear system easier to solve.
				add_to_p_equations(p_eqs, ss->psdp, ss->psfdu, stn_p, stn_u, rhs);
			}

			// Set linear system right-hand side
			int fgid = psfd_lid_to_gid(ss->psfdu[dim], flid);
			slv_set_bi(slv, fgid, rhs);

			/* Set linear system matrix row */
			set_composed_row(slv, fgid, ss->psdp, ss->psfdu, stn_p, stn_u);

			const real mag = fabs(rhs);
			if(mag > inf_norm) {
				inf_norm = mag;
			}
		}
		higfit_destroy(fit);
	}

	return inf_norm;
}

static void monolytic_method(sim_state *ss, sim_stencil * stn,
	const monolytic_system_impl *impl)
{
	/* Computing the Velocities and Pressure using the Newton-Raphson method. */

	/* Facet iterator */
	higfit_facetiterator *fit;

	/* Use pressure equations. */
	const bool use_p_eqs = false;

	/* Initializing the solver */
	static sim_stencil *stn_p = NULL;
	static sim_stencil *stn_u[DIM] = {NULL};
	static global_stencil **p_eqs = NULL;
	static solver *slv = NULL;
	size_t p_size;
	if (slv == NULL) {
		/* Create the stencil for pression elements. */
		stn_p = stn_create();

		/* Calculate total local domain size. */
		size_t u_size[DIM] = {0};
		p_size = psd_get_local_domain_size(ss->psdp);
		size_t localdomainsize = p_size;
		for(unsigned dim = 0; dim < DIM; ++dim) {
			u_size[dim] = psfd_get_local_domain_size(ss->psfdu[dim]);
			localdomainsize += u_size[dim];

			/* Create the stencil for velocity elements. */
			stn_u[dim] = stn_create();
		}

		// Pressure equations will be the continuity equations summed
		// with momentum equations that uses the corresponding pressure
		// point.
		if(use_p_eqs) {
			ALLOC_INFER(p_eqs, p_size);
			for(unsigned i = 0; i < p_size; ++i) {
				p_eqs[i] = gs_create();
			}
		}

		/* We need a global ID common to all physical properties. */
		int firstid = psd_set_joint_gid(1, &ss->psdp, DIM, ss->psfdu);

		/* Create the solver for the local domain */
		slv = slv_create(SOLVER_ANY, firstid, localdomainsize);

		// Set field split up.
		size_t us_size = 0;
		for(unsigned dim = 0; dim < DIM; ++dim) {
			us_size += u_size[dim];
		}
		slv_PETSc_set_field(slv, NULL,
				psfd_get_first_id(ss->psfdu[0]), us_size);

		slv_PETSc_set_field(slv, NULL, psd_get_first_id(ss->psdp), p_size);

		slv_set_maxiteration(slv, 2000);
	}

	/* When recalculating dt, use a smaller CFL, to allow for choosing
	 * for bigger dt. */
	const real rcfl = 0.8 * ss->cfl;

	/* Copy the old velocity field, to be used in temporal discrete term. */
	copy_facet_prop(ss->dpustar, ss->dpu);

	/* Assert we have a way out of the loop. */
	assert(impl->func_norm_cond || impl->diff_norm_cond);

	/* Non-linear solving iterations. */
	unsigned icount = 0;
	for(;;) {
		// Set momentum equations in the linear system.
		const real ns_norm = set_ns_eqs(slv, impl, ss, p_eqs, stn, stn_p, stn_u);

		// Set mass continuity equations in the linear system.
		const real mass_norm = set_mass_eqs(slv, impl, ss, p_eqs, stn, stn_u);

		if(impl->func_norm_cond && impl->func_norm_cond(mass_norm, ns_norm)) {
			break;
		}

		/* Solve the linear system */
		slv_assemble(slv);
		slv_solve(slv);

		print0f("Linear solver residual norm: %g (%d iterations)\n",
			slv_get_relative_residual_norm(slv),
			slv_get_iteration_count(slv));

		// Update current solution with the delta calculated by the
		// linear system.
		real new_dt = DBL_MAX;
		real diff_norm = impl->update_solution(slv, ss, &new_dt);

		{
			const real sendbuf[2] = {new_dt, -diff_norm};
			real recvbuf[2];
			MPI_Allreduce(&sendbuf, recvbuf, 2, MPI_HIGREAL,
				MPI_MIN, MPI_COMM_WORLD);
			new_dt = recvbuf[0];
			diff_norm = -recvbuf[1];
		}

		if(impl->diff_norm_cond && impl->diff_norm_cond(diff_norm)) {
			break;
		}

		if(ss->dt / new_dt > ss->cfl) {
			new_dt *= rcfl;
			print0f("Δt shrunk from velocity increase: %g L/U → %g L/U\n",
				ss->dt, new_dt);
			ss->dt = new_dt;

			// Since the system changed, reset the properties.
			// This is needed to converge difficult time steps
			copy_facet_prop(ss->dpu, ss->dpustar);
			//dp_set_all_values(ss->dpp, 0.0);
		}

		++icount;
	}

	for(unsigned dim = 0; dim < DIM; ++dim) {
		_adjust_outflow(ss->psfdu[dim], ss->dpu[dim], 3.0);
	}
}

/*********************************************************************************************/
/***    Framework for Monolytic Solvers: End                                               ***/
/*********************************************************************************************/



/*********************************************************************************************/
/***    Newton-Raphson Monolytic Solver Implementation: Begin                              ***/
/*********************************************************************************************/

static real newton_set_ns_eq(hig_facet *f, int lid, int dim, real Dt,
		real Re, sim_stencil *stn[2], sim_stencil *df_dp,
		sim_stencil *df_du[DIM], sim_facet_domain *sfdu[],
		distributed_property *dpu[], distributed_property *dpu_prev[],
		sim_domain *sdp, distributed_property *dpp)
{
	Point P;
	Point delta;
	hig_get_facet_center(f, P);
	hig_get_facet_delta(f, delta);

	const real Dx = delta[dim];

	real eval = 0.0;

	// Temporal term ( f = (uP - uP_prev) / dt ):
	const real uP = dp_get_value(dpu[dim], lid);
	{
		real uP_prev = dp_get_value(dpu_prev[dim], lid);

		// - value:
		eval += (uP - uP_prev) / Dt;
		// - d f / d uP:
		stn_add_to_element(df_du[dim], lid,
			1.0 / Dt
		);
	}

	Point p;
	POINT_ASSIGN(p, P);

	// Pressure term (f = (pH - pL) / Dx):
	{
		// High point:
		p[dim] = P[dim] + 0.5 * Dx;
		const real pH = sd_dp_interpolate(sdp, dpp, P, p, stn[0]);
		// High point jacobian elements:
		const real df_dpH = 1.0 / Dx;
		stn_mult_scalar_add(df_dp, stn[0], df_dpH);

		// Low point:
		p[dim] = P[dim] - 0.5 * Dx;
		const real pL = sd_dp_interpolate(sdp, dpp, P, p, stn[0]);
		// Low point jacobian elements:
		const real df_dpL = -1.0 / Dx;
		stn_mult_scalar_add(df_dp, stn[0], df_dpL);

		eval += (pH - pL) / Dx;

		// Restore p to center.
		p[dim] = P[dim];
	}

	// - Velocity terms:
	for(unsigned dim2 = 0; dim2 < DIM; ++dim2) {
		const real Dx = delta[dim2]; // TODO: verify what delta[] to use...
		const real diff_k = 1.0 / (Re * Dx * Dx);

		// Jacobian element for the central point in Diffusive term:
		// Being the diffusive term f = (2uP - uH - uL) / (Re * Dx²),
		// Then:
		// → d f / d uP = 2 / (Re * Dx²)
		stn_add_to_element(df_du[dim], lid,
			2.0 * diff_k
		);

		// Get high points:
		p[dim2] = P[dim2] + Dx;
		const real uH = sfd_dp_interpolate(sfdu[dim], dpu[dim], P, p, stn[0]);

		// Jacobian elements for the high points in the advective term:
		// Being the advective term f = (uH * vH - uL * vL) / (2 * Dx),
		// Then:
		// → d f / d uH = vH / (2 * Dx)
		// → d f / d vH = uH / (2 * Dx)
		real vH;
		if(dim2 != dim) {
			// Handle general case, where u != v:
			vH = sfd_dp_interpolate(sfdu[dim2], dpu[dim2], P, p, stn[1]);

			// Use chain rule of partial derivatives to calculate the
			// derivative relative to the stencil points.
			const real df_duH = 0.5 * vH / Dx;
			stn_mult_scalar_add(df_du[dim], stn[0], df_duH);

			const real df_dvH = 0.5 * uH / Dx;
			stn_mult_scalar_add(df_du[dim2], stn[1], df_dvH);
		} else {
			// Handle special case, where u == v
			// i.e.  f = (uH² - uL²) / (2 * Dx)
			// Then:
			// → d f / d uH = uH / Dx
			const real df_duH = uH / Dx;
			stn_mult_scalar_add(df_du[dim], stn[0], df_duH);
		}

		// Jacobian element for the high point in the diffusive term:
		// Being the diffusive term f = (2uP - uH - uL) / (Re * Dx²),
		// Then:
		// → d f / d uH = -1 / (Re * Dx²)
		stn_mult_scalar_add(df_du[dim], stn[0],
			-diff_k
		);

		// Get low points:
		p[dim2] = P[dim2] - Dx;
		const real uL = sfd_dp_interpolate(sfdu[dim], dpu[dim], P, p, stn[0]);

		// Finish up advective term:
		if(dim2 != dim) {
			// Jacobian elements for the low points in the advective term:
			// → d f / d uL = -vL / (2 * Dx)
			// → d f / d vL = -uL / (2 * Dx)
			const real vL = sfd_dp_interpolate(sfdu[dim2], dpu[dim2], P, p, stn[1]);

			const real df_duL = -0.5 * vL / Dx;
			stn_mult_scalar_add(df_du[dim], stn[0], df_duL);

			const real df_dvL = -0.5 * uL / Dx;
			stn_mult_scalar_add(df_du[dim2], stn[1], df_dvL);

			// Value
			eval += 0.5 * (uH * vH - uL * vL) / Dx;
		} else {
			// Special case, where u == v
			// Jacobian:
			// → d f / d uL = -uL / Dx
			const real df_duL = -uL / Dx;
			stn_mult_scalar_add(df_du[dim], stn[0], df_duL);

			// Value:
			eval += 0.5 * (uH * uH - uL * uL) / Dx;
		}

		// Jacobian element for the low point in the diffusive term:
		// → d f / d uL = -1 / (Re * Dx²)
		stn_mult_scalar_add(df_du[dim], stn[0],
			-diff_k
		);

		// The diffusive term:
		eval += (2*uP - uH - uL) * diff_k;

		// Restore p to center.
		p[dim2] = P[dim2];
	}

	return eval;
}

static real newton_set_mass_eq(hig_cell *c, int lid, sim_stencil *stn,
		sim_stencil *df_du[DIM], sim_facet_domain *sfdu[],
		distributed_property *dpu[])
{
	Point P;
	Point delta;
	hig_get_center(c, P);
	hig_get_delta(c, delta);

	real eval = 0.0;

	Point p;
	POINT_ASSIGN(p, P);

	for(unsigned dim = 0; dim < DIM; ++dim) {
		const real Dx = delta[dim];
		const real iDx = 1.0 / Dx;

		// High point:
		p[dim] = P[dim] + 0.5 * Dx;
		const real uH = sfd_dp_interpolate(sfdu[dim], dpu[dim], P, p, stn);
		// its jacobian:
		stn_mult_scalar_add(df_du[dim], stn, iDx);

		// Low point:
		p[dim] = P[dim] - 0.5 * Dx;
		const real uL = sfd_dp_interpolate(sfdu[dim], dpu[dim], P, p, stn);
		// its jacobian:
		stn_mult_scalar_add(df_du[dim], stn, -iDx);

		// Central difference scheme:
		eval += (uH - uL) * iDx;

		p[dim] = P[dim];
	}

	return eval;
}

static bool newton_func_norm_cond(const real mass_norm, const real ns_norm)
{
	real inf_norm = ns_norm > mass_norm ? ns_norm : mass_norm;
	{
		const real sendbuf = inf_norm;
		MPI_Allreduce(&sendbuf, &inf_norm, 1, MPI_HIGREAL,
			MPI_MAX, MPI_COMM_WORLD);
	}

	print0f("Inf Norm of residual vector: %g\n", inf_norm);

	// TODO: use a runtime parameter as tolerance:
	return inf_norm <= 1e-4;
}

static real newton_update_solution(solver *slv, sim_state *ss, real *dt_out)
{
	// Take the local domains.
	sim_facet_domain *sfdu[DIM];
	for (unsigned dim = 0; dim < DIM; ++dim) {
		sfdu[dim] = psfd_get_local_domain(ss->psfdu[dim]);
	}
	sim_domain *sdp = psd_get_local_domain(ss->psdp);

	// Get velocity.
	for(unsigned dim = 0; dim < DIM; ++dim)
	{
		higfit_facetiterator *fit;
		for(fit = sfd_get_domain_facetiterator(sfdu[dim]);
		     !higfit_isfinished(fit); higfit_nextfacet(fit))
		{
			hig_facet *f = higfit_getfacet(fit);
			int lid = sfd_get_local_id(sfdu[dim], f);

			real val = slv_get_xi(slv,
				psfd_lid_to_gid(ss->psfdu[dim], lid));

			dp_add_value(ss->dpu[dim], lid, -val);

			calc_advective_cfl_for_facet(dt_out, f,
				sfdu[dim], ss->dpu[dim]);
		}
		higfit_destroy(fit);
		dp_sync(ss->dpu[dim]);
	}

	// Get pressure.
	{
		higcit_celliterator *it;
		for(it = sd_get_domain_celliterator(sdp);
			!higcit_isfinished(it); higcit_nextcell(it))
		{
			hig_cell *c = higcit_getcell(it);
			int lid = sd_get_local_id(sdp, c);

			real val = slv_get_xi(slv,
				psd_lid_to_gid(ss->psdp, lid));

			real abs_val = fabs(val);

			dp_add_value(ss->dpp, lid, -val);
		}

		higcit_destroy(it);
	}
	dp_sync(ss->dpp);

	return 0.0;
}

static void newton_method(sim_state *ss, sim_stencil * stn)
{
	static const monolytic_system_impl impl = {
		.set_mass_eq = newton_set_mass_eq,
		.set_ns_eq = newton_set_ns_eq,
		.func_norm_cond = newton_func_norm_cond,
		.diff_norm_cond = NULL,
		.update_solution = newton_update_solution
	};

	monolytic_method(ss, stn, &impl);
}

/*********************************************************************************************/
/***    Newton-Raphson Monolytic Solver Implementation: End                                ***/
/*********************************************************************************************/



/*********************************************************************************************/
/***    Gauß-Siedel Monolytic Solver Implementation: Begin                                      ***/
/*********************************************************************************************/

static real gausss_set_mass_eq(hig_cell *c, int lid, sim_stencil *stn,
	sim_stencil *stn_u[DIM], sim_facet_domain *sfdu[],
	distributed_property *dpu[])
{
	Point P;
	Point delta;
	hig_get_center(c, P);
	hig_get_delta(c, delta);

	real rhs = 0.0;

	Point p;
	POINT_ASSIGN(p, P);

	for(unsigned dim = 0; dim < DIM; ++dim) {
		const real Dx = delta[dim];
		const real iDx = 1.0 / Dx;

		// High point:
		p[dim] = P[dim] + 0.5 * Dx;
		sfd_get_stencil(sfdu[dim], P, p, iDx, stn_u[dim]);

		// Low point:
		p[dim] = P[dim] - 0.5 * Dx;
		sfd_get_stencil(sfdu[dim], P, p, -iDx, stn_u[dim]);

		p[dim] = P[dim];

		rhs += stn_get_rhs(stn_u[dim]);
	}

	return rhs;
}

static real gausss_set_ns_eq(hig_facet *f, int lid, int dim, real Dt,
	real Re, sim_stencil *stn[2], sim_stencil *stn_p,
	sim_stencil *stn_u[DIM], sim_facet_domain *sfdu[],
	distributed_property *dpu[], distributed_property *dpu_prev[],
	sim_domain *sdp, distributed_property *dpp)
{
	Point P;
	Point delta;
	hig_get_facet_center(f, P);
	hig_get_facet_delta(f, delta);

	real rhs = 0;

	Point p;
	POINT_ASSIGN(p, P);

	// Pressure term (f = (pH - pL) / Dx):
	{
		const real Dx = delta[dim];
		const real iDx = 1.0 / Dx;

		// High point:
		p[dim] = P[dim] + 0.5 * Dx;
		sd_get_stencil(sdp, P, p, iDx, stn_p);

		// Low point:
		p[dim] = P[dim] - 0.5 * Dx;
		sd_get_stencil(sdp, P, p, -iDx, stn_p);

		rhs += stn_get_rhs(stn_p);

		// Restore p to center.
		p[dim] = P[dim];
	}

	// - Velocity terms:
	real uP_coef = 0.0;
	for(unsigned dim2 = 0; dim2 < DIM; ++dim2) {
		const real Dy = delta[dim2];

		// Diffusive term:
		// f = (2uP - uH - uL) / (Re * Dy²)

		const real diff_factor = -1.0 / (Dy * Dy * Re);

		// uP term:
		uP_coef += diff_factor;

		// High point in diffusive term:
		p[dim2] = P[dim2] + Dy;
		const real prev_uH = sfd_dp_interpolate(sfdu[dim], dpu[dim], P, p, stn[0]);

		stn_mult_scalar_add(stn_u[dim], stn[0], diff_factor);
		rhs += stn_get_rhs(stn[0]) * diff_factor;

		// Low point in diffusive term:
		p[dim2] = P[dim2] - Dy;
		const real prev_uL = sfd_dp_interpolate(sfdu[dim], dpu[dim], P, p, stn[1]);

		stn_mult_scalar_add(stn_u[dim], stn[1], diff_factor);
		rhs += stn_get_rhs(stn[1]) * diff_factor;

		// Advective term (non-conservative form):
		/*// f = vP * (uH - uL) / (2 * Dy)
		const real adv_factor = (prev_uH - prev_uL) / (2.0 * Dy);
		if(dim != dim2) {
			// vP stencil
			sfd_get_stencil(sfdu[dim2], P, adv_factor, stn_u[dim2]);
		} else {
			// set uP
			stn_add_to_element(stn_u[dim], lid, adv_factor);
		} */

		// Advective term (conservative form):
		// f = (uH * vH - uL * vL) / (2 * Dy)
		const real adv_factor = 0.5 / Dy;

		if(dim != dim2) {
			// uH != vH, so we need to retrieve the stencils again,
			// this time for v.

			// vL stencil
			//p[dim2] = P[dim2] - Dy; // p is already set
			sfd_get_stencil(sfdu[dim2], P, p, -prev_uL * adv_factor, stn_u[dim2]);

			// vH stencil
			p[dim2] = P[dim2] + Dy;
			sfd_get_stencil(sfdu[dim2], P, p, prev_uH * adv_factor, stn_u[dim2]);
		} else {
			// High point in advective term:
			stn_mult_scalar_add(stn_u[dim2], stn[0], prev_uH * adv_factor);

			// Low point in advective term:
			stn_mult_scalar_add(stn_u[dim2], stn[1], -prev_uL * adv_factor);
		}

		// Reset p to center
		p[dim2] = P[dim2];
	}
	// Finish up diffusive factor in uP:
	uP_coef *= -2.0;

	// Retrieve RHS
	for(unsigned dim = 0; dim < DIM; ++dim) {
		rhs += stn_get_rhs(stn_u[dim]);
	}

	// Temporal term:
	// f = (uP - uP_prev) / dt
	{
		const real factor = 1.0 / Dt;

		// stencil coefficient:
		uP_coef += factor;

		// RHS:
		const real uP_prev = dp_get_value(dpu_prev[dim], lid);
		rhs += uP_prev * factor;
	}

	// Store uP coefficient:
	stn_add_to_element(stn_u[dim], lid, uP_coef);

	return rhs;
}

static bool gausss_diff_norm_cond(const real diff_norm)
{
	print0f("Inf Norm of difference vector: %g\n", diff_norm);
	// TODO: use a runtime parameter as tolerance:
	return diff_norm <= 1e-4;
}

static real gausss_update_solution(solver *slv, sim_state *ss, real *dt_out)
{
	real norm = 0.0;

	// Take the local domains.
	sim_facet_domain *sfdu[DIM];
	for (unsigned dim = 0; dim < DIM; ++dim) {
		sfdu[dim] = psfd_get_local_domain(ss->psfdu[dim]);
	}
	sim_domain *sdp = psd_get_local_domain(ss->psdp);

	// Get velocity.
	for(unsigned dim = 0; dim < DIM; ++dim)
	{
		higfit_facetiterator *fit;
		for(fit = sfd_get_domain_facetiterator(sfdu[dim]);
		     !higfit_isfinished(fit); higfit_nextfacet(fit))
		{
			hig_facet *f = higfit_getfacet(fit);
			const int lid = sfd_get_local_id(sfdu[dim], f);

			const real nval = slv_get_xi(slv,
				psfd_lid_to_gid(ss->psfdu[dim], lid));
			const real oval = dp_get_value(ss->dpu[dim], lid);

			dp_set_value(ss->dpu[dim], lid, nval);

			const real diff = fabs(oval - nval);
			if(diff > norm) {
				norm = diff;
			}

			calc_advective_cfl_for_facet(dt_out, f,
				sfdu[dim], ss->dpu[dim]);
		}
		higfit_destroy(fit);
		dp_sync(ss->dpu[dim]);
	}

	// Get pressure.
	{
		int size = psd_get_local_domain_size(ss->psdp);
		for(int lid = 0; lid < size; ++lid)
		{
			const real nval = slv_get_xi(slv,
				psd_lid_to_gid(ss->psdp, lid));
			const real oval = dp_get_value(ss->dpp, lid);

			dp_set_value(ss->dpp, lid, nval);

			const real diff = fabs(oval - nval);
			if(diff > norm) {
				norm = diff;
			}
		}
	}
	dp_sync(ss->dpp);

	return norm;
}

static void gausss_method(sim_state *ss, sim_stencil * stn)
{
	static const monolytic_system_impl impl = {
		.set_mass_eq = gausss_set_mass_eq,
		.set_ns_eq = gausss_set_ns_eq,
		.func_norm_cond = NULL,
		.diff_norm_cond = gausss_diff_norm_cond,
		.update_solution = gausss_update_solution
	};

	monolytic_method(ss, stn, &impl);
}

/*********************************************************************************************/
/***    Jacobi Monolytic Solver Implementation: End                                        ***/
/*********************************************************************************************/



/*********************************************************************************************/
/***    Boundary Conditions: Begin                                                         ***/
/*********************************************************************************************/

sim_boundary *make_bc(hig_cell * bcg, int type)
{
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

void computeBCP(psim_domain * psd, int numbcs, hig_cell * bcg[], int pbctypes[],
		real pbcvalues[])
{
	for (int h = 0; h < numbcs; h++) {
		sim_domain *sd = psd_get_local_domain(psd);

		int bc_t = ((pbctypes[h] == 0) ? DIRICHLET : NEUMANN);
		sim_boundary *bc = make_bc(bcg[h], bc_t);

		/* Adding the boundary condition */
		sd_add_boundary(sd, bc);
		mp_mapper *bm = sb_get_mapper(bc);
		Point bcl, bch;
		hig_get_lowpoint(bcg[h], bcl);
		hig_get_highpoint(bcg[h], bch);
		real bclength = co_distance(bcl, bch);

		higcit_celliterator *it;
		/* Traversing the cells of the Diriclet boundaries conditions */
		for (it = sb_get_celliterator(bc); !higcit_isfinished(it);
		     higcit_nextcell(it)) {

			/* Getting the cell */
			hig_cell *bcell = higcit_getcell(it);

			/* Set the value 0.0 for the center of the cell */
			int bcgid = mp_lookup(bm, hig_get_cid(bcell));
			real bcval;
			if (bc_t == NEUMANN) {
				bcval = pbcvalues[h];
			} else {
				Point bcellcenter;
				hig_get_center(bcell, bcellcenter);
				real dist = co_distance(bcl, bcellcenter);
				real pos = dist / bclength;
				if (pos < 0.5) {
					bcval = pbcvalues[h] * pos * 2.0;
				} else {
					bcval =
					    pbcvalues[h] * (1.0 - pos) * 2.0;
				}
			}
			sb_set_value(bc, bcgid, bcval);
		}
		higcit_destroy(it);
	}
}

/*********************************************************************************************/
/***    Initial Conditions: Begin                                                          ***/
/*********************************************************************************************/

void initializeP(psim_domain * psdp, distributed_property * dpp)
{
	/* Setting the initial value for the pressure into the cavity */

	/* Setting the cell iterator */
	higcit_celliterator *it;

	/* Getting the local domain */
	sim_domain *sdp = psd_get_local_domain(psdp);

	/* Getting the Mapper for the local domain */
	mp_mapper *m = sd_get_domain_mapper(sdp);

	/* Traversing the cells of local domain */
	for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it);
	     higcit_nextcell(it)) {

		/* Getting the cell */
		hig_cell *c = higcit_getcell(it);

		/* Get the cell identifier */
		int cgid = mp_lookup(m, hig_get_cid(c));

		/* Set the value 0.0 for the pressure in this cell */
		dp_set_value(dpp, cgid, 0.0);
	}

	/* Destroying the iterator */
	higcit_destroy(it);

	dp_sync(dpp);
}

void initializeU(psim_facet_domain * psfdu[DIM],
		 distributed_property * dpu[DIM])
{
	/* Setting the initial value for the pressure into the cavity */

	/* Setting the facet-cell iterator */
	higfit_facetiterator *fit;

	/* Setting the volocity values for the domain */
	sim_facet_domain *sfdu[DIM];

	/* Setting the velocity U(dim) */
	for (int dim = 0; dim < DIM; dim++) {

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
		for (fit = sfd_get_domain_facetiterator(sfdu[dim]);
		     !higfit_isfinished(fit); higfit_nextfacet(fit)) {
			hig_facet *f = higfit_getfacet(fit);
			int fgid = mp_lookup(m, hig_get_fid(f));
			real val;
			val = 0.0;
			dp_set_value(dpu[dim], fgid, val);
		}
		higfit_destroy(fit);

		dp_sync(dpu[dim]);
	}
}

/*********************************************************************************************/
/***    Initial Conditions: End                                                            ***/
/*********************************************************************************************/

/*********************************************************************************************/
/***    Computing Properties Derivative at the Point: Begin                                ***/
/*********************************************************************************************/

void compute_facet_dudx_at_point(sim_facet_domain * sfdu, Point center,
				 Point delta, int dim, real alpha,
				 distributed_property * dpu, sim_stencil * stn,
				 real valueatpoint, real * dudx)
{
	/* Computing second order aproximation of first derivative of velocity at the point */

	/* Point to calculate the derivative */
	Point p;
	POINT_ASSIGN(p, center);

	/* Left point */
	p[dim] = center[dim] - alpha * delta[dim];
	real valuel =
	    sfd_dp_interpolate(sfdu, dpu, center, p, stn);

	/* Right point */
	p[dim] = center[dim] + alpha * delta[dim];
	real valueh =
	    sfd_dp_interpolate(sfdu, dpu, center, p, stn);

	/* First order derivative of velocity at the point p */
	*dudx = (valueh - valuel) / (2.0 * alpha * delta[dim]);
}

void compute_dpdx_at_point(sim_domain * sdp, Point center, Point delta, int dim,
			   real alpha, distributed_property * dpp,
			   sim_stencil * stn, real * dpdx)
{
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

void compute_facet_du2dx2_at_point(sim_facet_domain * sfdu, Point center,
				   Point delta, int dim, real alpha,
				   distributed_property * dpu,
				   sim_stencil * stn, real valueatpoint,
				   real * du2dx2)
{
	/* Computing second order aproximation of second derivative of velocity at the point */

	/* Point to calculate the derivative */
	Point p;
	POINT_ASSIGN(p, center);

	/* Left point */
	p[dim] = center[dim] - alpha * delta[dim];
	real valuel =
	    sfd_dp_interpolate(sfdu, dpu, center, p, stn);

	/* Right point */
	p[dim] = center[dim] + alpha * delta[dim];
	real valueh =
	    sfd_dp_interpolate(sfdu, dpu, center, p, stn);

	/* Second order derivative of pressure at the point p */
	*du2dx2 =
	    (valuel - 2.0 * valueatpoint +
	     valueh) / (delta[dim] * delta[dim] * alpha * alpha);
}

/*********************************************************************************************/
/***    Computing Properties Derivative at the Point: End                                  ***/
/*********************************************************************************************/

/*********************************************************************************************/
/***    Computing Pressure Solving Poisson Equation: Begin                                 ***/
/*********************************************************************************************/

static solver *setup_p_solver(sim_state *ss, sim_stencil * stn)
{
	/* This function sets up the solver for use as long as the matrix operaton doesn't change. */

	/* Initialing the solver */
	solver *slvp = NULL;

	/* Number of elements */
	int localdomainsize = psd_get_local_domain_size(ss->psdp);

	/* Creating the solver local domain */
	slvp = slv_create(SOLVER_ANY, psd_get_first_id(ss->psdp), localdomainsize);

	/* Setting the maximum number of non zeros values */
	// TODO: better estimate this number (maybe it is readly available)...
	slv_set_maxnonzeros(slvp, 300);

	/* Declaring the cell iterator */
	higcit_celliterator *it;

	/* Getting the local domain */
	sim_domain *sdp = psd_get_local_domain(ss->psdp);

	int total_nonzeros = 0;
	int cgid;

	/* Traversing the cells */
	for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it);
	     higcit_nextcell(it)) {

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

		/* Calculanting the parameters and the stencil points */
		for (int dim = 0; dim < DIM; dim++) {
			/* Parameter */
			real w = 1.0 / (cdelta[dim] * cdelta[dim]);
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

		/* Getting the local identifiers of real points */
		int *ids = psd_stn_get_gids(ss->psdp, stn);

		/* Getting the coeficients of real points */
		real *vals = stn_get_vals(stn);

		/* Getting the global identifier */
		cgid = psd_get_global_id(ss->psdp, c);

		/* Setting the matrix line for this point */
		slv_set_Ai(slvp, cgid, numelems, ids, vals);
	}

	printf("Non-zeros: %d\n", total_nonzeros);

	/* Destroying the iterator */
	higcit_destroy(it);

	/* Assembling the local matrix to the global matrix */
	slv_assemble_matrix(slvp);

	/* Assemble the output vector.
	 * I guess that assembling it only once here, the output
	 * of the previous timestep will be used as the input of
	 * the first one.*/
	slv_assemble_output(slvp);

	return slvp;
}

static void compute_p(sim_state *ss, sim_stencil * stn)
{
	/* This function solve the Poisson equation to calculate the pressure */

	/* Initialing the solver */
	static solver *slvp = NULL;
	static int myrank;
	if (slvp == NULL) {
		slvp = setup_p_solver(ss, stn);
		MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	}

	/* Declaring the cell iterator */
	higcit_celliterator *it;

	/* Getting the local domain */
	sim_domain *sdp = psd_get_local_domain(ss->psdp);
	sim_facet_domain *sfdu[DIM];

	/* Getting the velocty domain */
	for (int dim = 0; dim < DIM; dim++) {
		sfdu[dim] = psfd_get_local_domain(ss->psfdu[dim]);
	}

	/* Traversing the cells */
	for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it);
	     higcit_nextcell(it)) {

		/* Getting the cell */
		hig_cell *c = higcit_getcell(it);

		/* Getting the cell center and cell size */
		Point ccenter;
		Point cdelta;
		hig_get_center(c, ccenter);
		hig_get_delta(c, cdelta);

		/* Calculating the velocity divergent at the center */
		real sumdudx = 0.0;
		for (int dim = 0; dim < DIM; dim++) {
			/* First order derivative of the velocity */
			real dudx;
			compute_facet_dudx_at_point(sfdu[dim], ccenter, cdelta,
						    dim, 0.5, ss->dpustar[dim], stn,
						    0.0, &dudx);
			sumdudx += dudx;
		}

		/* Setting the right side of the equation */
		real rhs = sumdudx / ss->dt;

		/* Getting the global identifier */
		int cgid = psd_get_global_id(ss->psdp, c);

		/* Setting the independet vector for this point */
		slv_set_bi(slvp, cgid, rhs);
	}

	/* Destroying the iterator */
	higcit_destroy(it);

	/* Assembling the local matrix to the global matrix */
	slv_assemble(slvp);

	/* Solving the linear system */
	slv_solve(slvp);

	/* Getting the solution of the linear system */
	dp_slv_load_from_solver(ss->dpp, slvp);

	/* Display solver statistics */
	if(myrank == 0) {
		printf("Solver stats: iteration count = %d; rel residual norm = %e.\n",
		       slv_get_iteration_count(slvp),
		       slv_get_relative_residual_norm(slvp));
	}
}

/*********************************************************************************************/
/***    Computing Pressure Solving Poisson Equation: End                                   ***/
/*********************************************************************************************/

/*********************************************************************************************/
/***    Computing the Intermediate Velocity by Projection Method: Begin                    ***/
/*********************************************************************************************/

static void compute_ustar(int dim, sim_state *ss, sim_stencil * stn)
{
	/* Computing the intermediate veloicities using the first order projection method  */

	/* Mapping the domain and properties */
	higfit_facetiterator *fit;
	sim_facet_domain *sfdu[DIM];
	for (int dim = 0; dim < DIM; dim++) {
		sfdu[dim] = psfd_get_local_domain(ss->psfdu[dim]);
	}
	mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);

	/* iterating on the facets for compute the intermediate velocity */
	for (fit = sfd_get_domain_facetiterator(sfdu[dim]);
	     !higfit_isfinished(fit); higfit_nextfacet(fit)) {
		/* Getting the identifiers and data at the facets */
		hig_facet *f = higfit_getfacet(fit);
		int fgid = mp_lookup(mu, hig_get_fid(f));
		Point fcenter;
		Point fdelta;
		hig_get_facet_center(f, fcenter);
		hig_get_facet_delta(f, fdelta);

		/* Getting the value of the velocity at center of facet */
		real valuec = dp_get_value(ss->dpu[dim], fgid);

		/* Initializing the right side of the Navier-Equation */
		real tot = 0.0;

		/* Computing the derivatives */
		//              DEBUG_PASS;
		for (int dim2 = 0; dim2 < DIM; dim2++) {

			/* Computing the first order derivative */
			real dudx = 0.0;
			compute_facet_dudx_at_point(sfdu[dim], fcenter, fdelta,
						    dim2, 1.0, ss->dpu[dim], stn,
						    valuec, &dudx);

			/* Computing the first order derivative */
			real du2dx2 = 0.0;
			compute_facet_du2dx2_at_point(sfdu[dim], fcenter,
						      fdelta, dim2, 1.0,
						      ss->dpu[dim], stn, valuec,
						      &du2dx2);

			/* Computing the value of the velocity at the center of the facet */
			/* for the convective term                                        */
			real u;
			u = sfd_dp_interpolate(sfdu[dim2], ss->dpu[dim2], fcenter, fcenter, stn);

			/* Convective term */
			tot += -u * dudx;

			/* Diffusive term */
			tot += du2dx2 / ss->Re;
		}

		/* External force */
		tot += -bodyforce(fcenter, dim, ss->Re);

		/* Calculating the intermediate velocity */
		real ustar = valuec + ss->dt * tot;

		/* Storing the value of the intermediate velocity */
		dp_set_value(ss->dpustar[dim], fgid, ustar);
	}

	/* Destroying the iterator */
	higfit_destroy(fit);

	/* Sincronizing the intermediate velocity for the domain */
	dp_sync(ss->dpustar[dim]);
}

/*********************************************************************************************/
/***    Computing the Intermediate Velocity by Projection Method: End                      ***/
/*********************************************************************************************/

/*********************************************************************************************/
/***    Computing the Final Velocity by the Projection Method: Begin                       ***/
/*********************************************************************************************/

static void compute_ut(int dim, sim_state *ss, sim_stencil * stn)
{
	/* Computing the final veloicities and pressure using the first order projection method  */

	/* Mapping the domain and properties */
	sim_domain *sdp = psd_get_local_domain(ss->psdp);
	sim_facet_domain *sfdu[DIM];
	higfit_facetiterator *fit;
	sfdu[dim] = psfd_get_local_domain(ss->psfdu[dim]);
	mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);

	/* iterating on the facets for update the final velocity */
	for (fit = sfd_get_domain_facetiterator(sfdu[dim]);
	     !higfit_isfinished(fit); higfit_nextfacet(fit)) {

		/* Getting the identifiers and data at the facets */
		hig_facet *f = higfit_getfacet(fit);
		int fgid = mp_lookup(mu, hig_get_fid(f));
		Point fcenter;
		Point fdelta;
		hig_get_facet_center(f, fcenter);
		hig_get_facet_delta(f, fdelta);

		/* Computing the pressure gradient */
		real dpdx;
		compute_dpdx_at_point(sdp, fcenter, fdelta, dim, 0.5, ss->dpp,
				stn, &dpdx);

		/* Getting the intermetiate velocity */
		real ustar = dp_get_value(ss->dpustar[dim], fgid);

		/* Computing the final velocity */
		real utdt = ustar - ss->dt * dpdx;

		/* Storing the final velocity */
		dp_set_value(ss->dpu[dim], fgid, utdt);
	}

	/* destroying the iterator */
	higfit_destroy(fit);

	/* Synchronizing distributed proprieties */
	_adjust_outflow(ss->psfdu[dim], ss->dpu[dim], 3.0);
	dp_sync(ss->dpu[dim]);
}

/*********************************************************************************************/
/***    Computing the Final Velocity by the Projection Method: End                         ***/
/*********************************************************************************************/

/*********************************************************************************************/
/***    Solve the Incomnpressible Navier-Stokes using the Projetion Method: Begin          ***/
/*********************************************************************************************/

static void projection_method(sim_state *ss, sim_stencil *stn)
{
	/* Computing the Velocities and Pressure using the first order projection method  */

	/* Facet iterator */
	higfit_facetiterator *fit;

	/* Computing the intermediate velocity from projection method */
	/* See report on the projection method - page 1               */
	for (int dim = 0; dim < DIM; dim++) {
		compute_ustar(dim, ss, stn);
	}

	/* Computing the pressure by Poisson Equation solver */
	/* See report on the projection method - page 2      */
	compute_p(ss, stn);

	/* Computing the final velocity from the projection method */
	/* See report on the projection method - page 3            */
	for (int dim = 0; dim < DIM; dim++) {
		compute_ut(dim, ss, stn);
	}
}

/*********************************************************************************************/
/***    Solve the Incomnpressible Navier-Stokes using the Projetion Method: End            ***/
/*********************************************************************************************/



/*********************************************************************************************/
/***    Printing velocytie error: Begin                                                    ***/
/*********************************************************************************************/

void printResults(psim_domain * psdp, psim_facet_domain * psfdu[DIM],
		  int numprobes, hig_cell * probe[], distributed_property * dpp,
		  distributed_property * dpu[DIM])
{
	/* This function print the properties in some points at the standad output */

	/* Getting the local domain */
	sim_domain *sdp = psd_get_local_domain(psdp);

	/* Creating the stencil */
	sim_stencil *stn = stn_create();

	/* Getting the local domain */
	int hig = 0;
	Point l, h;
	hig_cell *rootp = sd_get_higtree(sdp, hig);
	sim_facet_domain *sfdu[DIM];
	for (int dim = 0; dim < DIM; dim++) {
		sfdu[dim] = psfd_get_local_domain(psfdu[dim]);
	}

	/* Getting the bounding box */

	/* traversing the cells in the bounding box */
	higcit_celliterator *it;
	for (int p = 0; p < numprobes; p++) {
		for (it = higcit_create_all_leaves(probe[p]);
		     !higcit_isfinished(it); higcit_nextcell(it)) {

			/* Getting the cell */
			hig_cell *c = higcit_getcell(it);

			/* Getting the center and size */
			Point ccenter;
			Point cdelta;
			hig_get_center(c, ccenter);
			hig_get_delta(c, cdelta);

			/* Setting the center for interpolation and compute the value of the property */
			real value[DIM];
			for (int dim = 0; dim < DIM; dim++) {
				value[dim] = sfd_dp_interpolate(sfdu[dim], dpu[dim], ccenter, ccenter, stn);
			}

			/* Printing the value of the velocity */
			printf("v(");

			/* Printing the point */
			for (int dim = 0; dim < DIM; dim++) {
				printf("%14.10f ", ccenter[dim]);
			}
			printf(") = ");

			/* Printing the value */
			for (int dim = 0; dim < DIM; dim++) {
				printf("%20.15f ", value[dim]);
			}

			printf("\n");
		}
		printf("\n");
	}


	/* Destroying the stencil */
	stn_destroy(stn);
}

/*********************************************************************************************/
/***    Printing velocytie error: End                                                      ***/
/*********************************************************************************************/

/*void __print_sfbi(psim_facet_domain *psfd, int hig, sim_facet_block_info *sfbi) {
    partition_graph *pg = psfd_get_partition_graph(psfd);
    for(int j = 0; j < pt->domaingridsize[hig][1]; j++) {
        for(int i = 0; i < pt->domaingridsize[hig][0]; i++) {
            sim_facet_block_info *s = &sfbi[j+(i*pt->domaingridsize[hig][1])];
            printf(" %c%d%c", (s->l?'.':' '), s->type, (s->h?'.':' '));
        }
        printf("\n");
    }
    printf("\n\n");
}*/

/*********************************************************************************************/
/***    Main Program: Begin                                                                ***/
/*********************************************************************************************/

void computeBCU(psim_facet_domain * psfdu[DIM], int numbcs, hig_cell * bcg[numbcs],
	int *bctypes[], real *bcvalues[])
{
	sim_facet_domain *sfd;
	for (int dim = 0; dim < DIM; dim++) {
		for (int h = 0; h < numbcs; h++) {
			sfd = psfd_get_local_domain(psfdu[dim]);

			int bc_t =
			    ((bctypes[dim][h] == 0) ? DIRICHLET : NEUMANN);
			sim_boundary *bc = make_bc(bcg[h], bc_t);
			Point bcl, bch;
			hig_get_lowpoint(bcg[h], bcl);
			hig_get_highpoint(bcg[h], bch);
			real bclength = co_distance(bcl, bch);

			/* Adding the boundary condition */
			sfd_add_boundary(sfd, bc);
			mp_mapper *bm = sb_get_mapper(bc);

			higcit_celliterator *it;
			/* Traversing the cells of the Diriclet boundaries conditions */
			for (it = sb_get_celliterator(bc);
			     !higcit_isfinished(it); higcit_nextcell(it)) {

				/* Getting the cell */
				hig_cell *bcell = higcit_getcell(it);

				/* Set the value 0.0 for the center of the cell */
				int bcgid = mp_lookup(bm, hig_get_cid(bcell));
				real bcval;
				if (bc_t == NEUMANN) {
					bcval = bcvalues[dim][h];
				} else {
					Point bcellcenter;
					hig_get_center(bcell, bcellcenter);
					real dist =
					    co_distance(bcl, bcellcenter);
					real pos = dist / bclength;
					if (pos < 0.5) {
						bcval =
						    bcvalues[dim][h] * pos *
						    2.0;
					} else {
						bcval =
						    bcvalues[dim][h] * (1.0 -
									pos) *
						    2.0;
					}
				}
				sb_set_value(bc, bcgid, bcval);
			}
			higcit_destroy(it);
		}
	}
}

int main(int argc, char *argv[])
{
	/* Initializing Higtree */
	higtree_initialize(&argc, &argv);

	int myrank;
	int ntasks;
	int mpierr = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	/* Partitioning the domain */
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

	if(argc < 9) {
		fprintf(stderr,
			"Error!\n"
			"Missing parameters:\n"
			"  %s \\\n"
			"    num_domain_tress_N {tree_1.amr ... tree_N.amr} \\\n"
			"    num_bc_trees_M {bc_tree_1.amr p_bc_1_type p_bc_1_value \\\n"
			"                       {u[0]_bc_1_type u[0]_bc_1_value ... u[D]_bc_1_type u[D]_bc_1_value} \\\n"
			"                    ... \\\n"
			"                    bc_tree_M.amr p_bc_M_type p_bc_M_value \\\n"
			"                       {u[0]_bc_M_type u[0]_bc_M_value ... u[D]_bc_M_type u[D]_bc_M_value} \\\n"
			"                   } \\\n"
			"    num_probe_trees_P {probe_tree_1.amr ... probe_tree_P.amr} \\\n"
			"    cfl reynolds total_simulation_time output_frequency \\\n"
			"    (projection | monolytic-gs | monolytic-newton) output_prefix\n",
			argv[0]
		);
		return 1;
	}

	/* Define simulation parameters variable. */
	sim_state ss;

	/* Adapting the data input by the command line */
	int numhigs = atoi(argv[1]);
	argv++;
	char *amrfilename[numhigs];
	for (int h = 0; h < numhigs; h++) {
		amrfilename[h] = argv[1];
		argv++;
	}
	int total_numbcs = atoi(argv[1]);
	argv++;
	char *amrBCfilename[total_numbcs];
	int pbctypes[total_numbcs];
	int ubctypes[DIM][total_numbcs];
	real pbcvalues[total_numbcs];
	real ubcvalues[DIM][total_numbcs];
	for (int h = 0; h < total_numbcs; h++) {
		amrBCfilename[h] = argv[1];
		argv++;
		pbctypes[h] = atoi(argv[1]);
		argv++;
		pbcvalues[h] = atof(argv[1]);
		argv++;
		for(unsigned d = 0; d < DIM; ++d) {
			ubctypes[d][h] = atoi(argv[1]);
			argv++;
			ubcvalues[d][h] = atof(argv[1]);
			argv++;
		}
	}

	int total_numprobes = atoi(argv[1]);
	argv++;
	char *probefilename[total_numprobes];
	for (int h = 0; h < total_numprobes; h++) {
		probefilename[h] = argv[1];
		argv++;
	}

	ss.cfl = atof(argv[1]);
	argv++;
	ss.Re = atof(argv[1]);
	argv++;
	real total_time = atof(argv[1]);
	argv++;
	real output_frequency = atof(argv[1]);

	/* Choose the Navier-Stokes solver method to use. */
	argv++;
	const char* method_name = argv[1];
	ns_solver_func updatePandU;
	if(strcmp(method_name, "projection") == 0) {
		updatePandU = projection_method;
	} else if(strcmp(method_name, "monolytic-gs") == 0) {
		updatePandU = gausss_method;
	} else if(strcmp(method_name, "monolytic-newton") == 0) {
		updatePandU = newton_method;
	} else {
		print0f("Choose one of the methods to use:\n"
			"  projection - Projection method\n"
			"  monolytic-gs - Gauß-Siedel non-linear solver\n"
			"  monolytic-newton - Newton-Raphson non-linear solver\n");
		return 1;
	}

	argv++;
	char *results = argv[1];

	/* Print some simulation information. */
	real output_period = 1.0 / output_frequency;

	/* Number of calculated output steps.
	 * Not counting +1 from the initial state.*/
	unsigned num_out_steps = ceil(total_time * output_frequency);

	if(myrank == 0) {
		printf("Simulation parameters, where L * U / ν = Re = %g:\n"
			" - Solution Method: %s\n"
			" - CFL: %g\n"
			" - Output frequency: %g U/L\n"
			" - Output period: %g L/U\n"
			" - Number of output steps: %u\n"
			" - Total time requested: %g L/U\n"
			" - Time of last output (rounded from output frequency): %g L/U\n",
			ss.Re, method_name, ss.cfl, output_frequency,
			output_period, num_out_steps + 1, total_time,
			output_period * num_out_steps
		);
	}

	// Distributed partition context:
	// Each boundary condition will have its own group, so that we can identify
	// its properties later by the group.

	print0f("Creating partitioning context.\n");

	unsigned num_groups = 1 + total_numbcs + total_numprobes;
	load_balancer *part_ctx = lb_create(MPI_COMM_WORLD, num_groups);
	for (unsigned i = 1; i < num_groups; ++i)
		lb_set_group_weight(part_ctx, i, 0.0001);

	// Disable fringe for probes
	for (unsigned i = total_numbcs + 1; i < num_groups; ++i)
		lb_set_group_has_fringe(part_ctx, i, false);

	print0f("Loading domain trees.\n");

	// Load domain trees
	FILE *fd;
	for (int h = myrank; h < numhigs; h += ntasks) {
		/* Open the AMR format file */
		fd = fopen(amrfilename[h], "r");

		/* Creating the distributed HigTree data structure */
		hig_cell *root = higio_read_from_amr(fd);

		/* Closing the AMR format file */
		fclose(fd);

		lb_add_input_tree(part_ctx, root, true, 0);
	}

	print0f("Loading boudary condition trees.\n");

	// Load boundary condition trees
	for (int h = myrank; h < total_numbcs; h += ntasks) {
		/* Open the AMR format file */
		fd = fopen(amrBCfilename[h], "r");

		/* Creating the distributed HigTree data structure */
		hig_cell *root = higio_read_from_amr(fd);

		/* Closing the AMR format file */
		fclose(fd);

		lb_add_input_tree(part_ctx, root, true, h + 1);
	}

	for (int h = myrank; h < total_numprobes; h++) {
		fd = fopen(probefilename[h], "r");
		hig_cell *probe = higio_read_from_amr(fd);
		fclose(fd);

		lb_add_input_tree(part_ctx, probe, true,
					   h + total_numbcs + 1);
	}

	print0f("Partitioning the domain.\n");

	/* Taking the sub-domain of the myrank process */
	partition_graph *pg = pg_create(MPI_COMM_WORLD);

	/* Setting the fringe size of the sub-domain */
	pg_set_fringe_size(pg, 4);

	/* Perform the actual partitioning... */
	lb_calc_partition(part_ctx, pg);

	print0f("Creating cell domains.\n");

	/* Create the local domain from the partitioned trees. */
	sim_domain *sdp = sd_create(NULL);
	//sd_use_cache(sdp, 0);
	size_t tidx = 0;
	numhigs = lb_get_num_local_trees_in_group(part_ctx, 0);
	for (; tidx < numhigs; ++tidx) {
		unsigned group;
		hig_cell *tree =
		    lb_get_local_tree_in_group(part_ctx, tidx, 0);
		sd_add_higtree(sdp, tree);
	}

	/* Creating the partitioned sub-domain to simulation */
	ss.psdp = psd_create(sdp, pg);

	/* Setting the order of the polynomial interpolation */
	sd_set_interpolator_order(sdp, 2);

	/* Creating the vector properties from facets */

	print0f("Creating facet domains.\n");

	/* Facets list */
	sim_facet_domain *sfdu[DIM];

	///* Linking the property to the facets */
	for (int dim = 0; dim < DIM; dim++) {

		/* Creating the list of cells center in the local domain */
		sfdu[dim] = sfd_create(NULL, dim);
		//sfd_use_cache(sfdu[dim], 0);

		/* Setting the order 2 for the properties interpolation */
		sfd_set_interpolator_order(sfdu[dim], 2);

		/* Copying the list of cells center in the domain for pressure */
		sfd_copy_higtrees_from_center_domain(sfdu[dim], sdp);

		sfd_adjust_facet_ids(sfdu[dim]);

		/* Creating property for the facets */
		ss.psfdu[dim] = psfd_create(sfdu[dim], ss.psdp);
	}

	print0f("Assigning boundary conditions trees.\n");

	unsigned max_trees = lb_get_num_local_trees(part_ctx) - 1;

	/* Get the bcs assigned to this process. */
	{
		hig_cell *bctrees[max_trees];
		int l_pbctypes[max_trees];
		real l_pbcvalues[max_trees];

		int *l_ubctypes[DIM];
		real *l_ubcvalues[DIM];
		for(unsigned d = 0; d < DIM; ++d) {
			ALLOC_INFER(l_ubctypes[d], max_trees);
			ALLOC_INFER(l_ubcvalues[d], max_trees);
		}

		unsigned numbcs = 0;
		for (unsigned g = 1; g <= total_numbcs; ++g) {
			unsigned bcid = g - 1;

			unsigned group_count =
			    lb_get_num_local_trees_in_group(part_ctx, g);

			for (unsigned i = 0; i < group_count; ++i) {
				bctrees[numbcs] =
				    lb_get_local_tree_in_group(part_ctx, i, g);

				l_pbctypes[numbcs] = pbctypes[bcid];
				l_pbcvalues[numbcs] = pbcvalues[bcid];
				for(unsigned d = 0; d < DIM; ++d) {
					l_ubctypes[d][numbcs] = ubctypes[d][bcid];
					l_ubcvalues[d][numbcs] = ubcvalues[d][bcid];
				}

				++numbcs;
			}
		}

		/* Pressure boundary condition */
		computeBCP(ss.psdp, numbcs, bctrees, l_pbctypes, l_pbcvalues);

		/* Velocity boundary condition */
		computeBCU(ss.psfdu, numbcs, bctrees, l_ubctypes, l_ubcvalues);

		for(unsigned d = 0; d < DIM; ++d) {
			free(l_ubctypes[d]);
			free(l_ubcvalues[d]);
		}
	}

	print0f("Creating distributed properties.\n");

	/* Creating the distributed properties */
	psd_synced_mapper(ss.psdp);

	/* Distributed property for pressure */
	ss.dpp = psd_create_property(ss.psdp);

	/* Get probes assigned to this proc. */
	DECL_AND_ALLOC(hig_cell *, probe, max_trees);
	unsigned numprobes = 0;
	for (unsigned g = total_numbcs + 1; g < num_groups; ++g) {
		unsigned prid = g - 1;

		unsigned group_count =
		    lb_get_num_local_trees_in_group(part_ctx, g);
		for (unsigned i = 0; i < group_count; ++i) {
			probe[numprobes] =
			    lb_get_local_tree_in_group(part_ctx, i, g);

			++numprobes;
		}
	}
	REALLOC_INFER(probe, numprobes);

	lb_destroy(part_ctx);

	/* Computing the boundary conditions */
	print0f("Computind boundary conditions.\n");

	for (int dim = 0; dim < DIM; dim++) {
		psfd_compute_sfbi(ss.psfdu[dim]);

		/* Mapping the properties in the domain */
		psfd_synced_mapper(ss.psfdu[dim]);

		/* Final velocity */
		ss.dpu[dim] = psfd_create_property(ss.psfdu[dim]);

		/* Intermediate velocity */
		ss.dpustar[dim] = psfd_create_property(ss.psfdu[dim]);
	}

	sim_facet_domain *sfd = psfd_get_local_domain(ss.psfdu[0]);
	mp_mapper *m = sfd_get_domain_mapper(sfd);

	///* Intitializing the properties (velocity and pressure) */
	print0f("Initializing phisical properties.\n");
	initializeP(ss.psdp, ss.dpp);
	initializeU(ss.psfdu, ss.dpu);

	/* Setting the number of steps */

	print0f("Setting up for time integration loop.\n");

	/* Creating the stencil for properties interpolation */
	sim_stencil *stn = stn_create();

	/* Create the xdmf output context if needed. */
	xdmf_output *xdmf_out = xdmf_init(results, sdp);
	xdmf_register_cell_property(xdmf_out, "Pressure", 1, ss.dpp);
	xdmf_register_facet_property_array(xdmf_out, "Velocity", DIM, sfdu, ss.dpu);

	print0f("Calculating diffusive CFL.\n");

	/* Minimum diffusive CFL won't change unless the tree changes, so
	 * calculate it only once. */
	const real diffusive_dt = calc_diffusive_cfl(ss.cfl, sdp, ss.Re);

	/* For the first time step */
	ss.dt = 0.0;

	/* Time accumulated since last output */
	real frame_time = 0.0;

	/* Initializing the number of frames for print */
	int numframe = 0;

	/* Currently output steps. */
	unsigned output_steps = 0;

	/* Used only on process 0, for output. */
	unsigned total_steps = 0;

	/* If next time step is an output frame. */
	bool print_output = true;

	for(;;) {
		/* First half of iteration loop is data output. */
		print0f("\nStep %u: Δt = %g L/U; t = %g L/U", total_steps++, ss.dt, numframe * output_period + frame_time);

		/* Printing the results */
		if (print_output) {
			/* Printing the velocities to visualize */
			xdmf_write_timestep(xdmf_out, numframe, numframe);

			print0f(" (output file %u/%u)", numframe+1, num_out_steps+1);

			/* Printing the results to standard output */
			printResults(ss.psdp, ss.psfdu, numprobes, probe, ss.dpp, ss.dpu);
		}
		print0f("\n");

		/* The loop exit is in the middle, because we want data output of the first
		 * state before the simulation, and of the last step simulated. */

		if(numframe >= num_out_steps) {
			break;
		}

		/* Second half of iteration loop is the time integration step. */

		/* Calculate the next timestep using CFL. */
		ss.dt = calc_advective_cfl(ss.cfl, sfdu, ss.dpu);
		if(diffusive_dt < ss.dt) {
			ss.dt = diffusive_dt;
		}

		/* Clamp to the next outputting frame time, if needed. */
		bool may_output_frame = false;
		real frame_dt;
		if(ss.dt + frame_time >= output_period) {
			ss.dt = output_period - frame_time;
			may_output_frame = true;
			frame_dt = ss.dt;
		}

		/* Update velocities and pressure using the projection method */
		updatePandU(&ss, stn);

		if(may_output_frame && ss.dt >= frame_dt) {
			assert(ss.dt == frame_dt);
			frame_time = 0.0;
			print_output = true;

			/* Increasing the number of prints */
			numframe++;
		} else {
			print_output = false;
			frame_time += ss.dt;
		}

		if (myrank == 0) {
			DEBUG_INSPECT(time, %lf);
		}
	}

	print0f("Cleaning up.\n");

	xdmf_destroy(xdmf_out);

	/* Destroing the properties objects */
	for (int dim = 0; dim < DIM; dim++) {
		dp_destroy(ss.dpu[dim]);
		dp_destroy(ss.dpustar[dim]);
	}
	dp_destroy(ss.dpp);

	/* Destroing the stencil for interpolation */
	stn_destroy(stn);

	///* Getting the time */
	//DEBUG_DIFF_TIME;
}

/*********************************************************************************************/
/***    Main Program: End                                                                  ***/
/*********************************************************************************************/
