#if DIM == 2
// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver geometric split advection  - version 11/2022
// *******************************************************************
// *******************************************************************

// ***********************************************************************
// Volume Fraction Transport Step with PLIC fractional step on direction x
// ***********************************************************************

#include "hig-flow-vof-advection-2D.h"
#include "hig-flow-vof-adap-hf.h"


// ***********************************************************************
// Volume Fraction Transport Step with PLIC fractional step on direction x
// ***********************************************************************
void higflow_plic_advection_volume_fraction_x_direction_2D(higflow_solver *ns, int dim) {
	// Get the local sub-domain for the cells
	sim_domain *sdp = psd_get_local_domain(ns->psdp);
	// Get the local sub-domain for the facets
	sim_facet_domain *sfdu[DIM];
	for(int i = 0; i < DIM; i++) {
		sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
	}
	// Get the map for the domain properties
	mp_mapper *mp = sd_get_domain_mapper(sdp);
	// Loop for each cell
	higcit_celliterator *it;
	real tol_u = 1.0e-8;
	for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
		// Get the cell
		hig_cell *c = higcit_getcell(it);
		// Get the cell identifier
		int clid    = mp_lookup(mp, hig_get_cid(c));
		// Get the center of the cell
		Point ccenter;
		hig_get_center(c, ccenter);
		// Get the delta of the cell
		Point cdelta;
		hig_get_delta(c, cdelta);
		// Get the velocity at facet
		int infacet;
//			// Get the velocity in the left facet center
		real ul = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
//			// Get the velocity in the right facet center
		real ur = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);

		
		//real ur, ul;
		
		////===== Zalesak' test===========================================
		//ur = 0.5 - ccenter[1];
		//ul = ur;
	    ////==============================================================
		
		//real pi = 3.14159265358979323846;
		//real T = 8.0;
		////===== Shearing flow===========================================
		//ul = -pow(sin(pi*ccenter[0]) ,2)*sin(2.0*pi*ccenter[1])*cos(pi*ns->par.t/T);
		//ur = ul;
		////==============================================================
		
		
		//arquivoV(nome_vx,ccenter[0],ccenter[1],ul,ur);
		
		Point Normal, p, Delta_New;
		real d, fracvol,fracr,fracl;
		
		p[0]=ccenter[0];p[1]=ccenter[1];
		p[dim]=p[dim]+0.5*cdelta[dim];
		fracr=compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
		
		p[0]=ccenter[0];p[1]=ccenter[1];
		p[dim]=p[dim]-0.5*cdelta[dim];
		fracl=compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
		
		
		fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
		
		fraction_correction_at_get(&fracr);
		fraction_correction_at_get(&fracl);
		fraction_correction_at_get(&fracvol);
		
		/*if((fracr==1.0&&fracl==1.0&&fracvol==1.0)||(fracr==0.0&&fracl==0.0&&fracvol==0.0)){
			
			dp_set_value(ns->ed.mult.dpfracvolaux, clid, fracvol);
			continue;
			
		}*/

		// Area
		real Ar = 0.0;
		real Al = 0.0;
		fracr = 0.0;
		fracl = 0.0;
		//  Right Facet
		if(fabs(ur)>tol_u){
			if (fabs(ur) * ns->par.dt > 0.5 * cdelta[0]) {
				printf("Time step is large!!!\n");
				exit(1);
			}
			if(ur>0.0){
				p[0] = ccenter[0];
				p[1] = ccenter[1];
				// Normal
				Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
				Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
				// Correction of Normal
				normal_correction_at_get(Normal);
				// Fraction
				fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
				// Fraction correction
				fraction_correction_at_get(&fracvol);
				if(fracvol==1.0 || fracvol==0.0 || (Normal[0]==0.0 && Normal[1]==0.0)){
					Ar = fracvol*cdelta[1]*ur*ns->par.dt;
				} else {
					d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
					d = d - 0.5*(cdelta[0]-ur*ns->par.dt)*Normal[0];
					Delta_New[0] = ur*ns->par.dt;
					Delta_New[1] = cdelta[1];
					Ar = area_left_line_origin_center(Normal, Delta_New, d);
				}
				fracr = Ar/(fabs(ur)*ns->par.dt*cdelta[1]);
			} else {
				p[0] = ccenter[0]+cdelta[0];
				p[1] = ccenter[1];
				// Normal
				Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
				Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
				// Correction of Normal
				normal_correction_at_get(Normal);
				// Fraction
				fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
				// Fraction correction
				fraction_correction_at_get(&fracvol);
				if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
					Ar = fracvol * cdelta[1] * fabs(ur) * ns->par.dt;
				} else {
					d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
					Delta_New[0] = fabs(ur)*ns->par.dt;
					Delta_New[1] = cdelta[1];
					d = d + 0.5 * (cdelta[0] - fabs(ur) * ns->par.dt)*Normal[0];
					Ar = area_left_line_origin_center(Normal, Delta_New, d);
				}
				fracr = Ar/(fabs(ur)*ns->par.dt*cdelta[1]);
			}
		}
		//  Left Facet
		if (fabs(ul) > tol_u) {
			if (fabs(ul) * ns->par.dt > 0.5 * cdelta[0]) {
				printf("Time step is large!!!\n");
				exit(1);
			}
			if (ul > 0.0) {
				p[0] = ccenter[0] - cdelta[0];
				p[1] = ccenter[1];
				// Normal
				Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
				Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
				// Correction of Normal
				normal_correction_at_get(Normal);
				// Fraction
				fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
				// Fraction correction
				fraction_correction_at_get(&fracvol);
				if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
					Al = fracvol * cdelta[1] * ul * ns->par.dt;
				} else {
					d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
					d = d - 0.5 * (cdelta[0] - ul * ns->par.dt) * Normal[0];
					Delta_New[0] = ul*ns->par.dt;
					Delta_New[1] = cdelta[1];
					Al = area_left_line_origin_center(Normal, Delta_New, d);
				}
				fracl = Al/(fabs(ul)*ns->par.dt*cdelta[1]);
			} else {
				p[0] = ccenter[0];
				p[1] = ccenter[1];
				// Normal
				Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
				Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
				// Correction of Normal
				normal_correction_at_get(Normal);
				// Fraction
				fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
				fraction_correction_at_get(&fracvol);
				if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
					Al = fracvol * cdelta[1] * fabs(ul) * ns->par.dt;
				} else {
					d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
					d = d + 0.5 * (cdelta[0] - fabs(ul) * ns->par.dt) * Normal[0];
					Delta_New[0] = fabs(ul)*ns->par.dt;
					Delta_New[1] = cdelta[1];
					Al = area_left_line_origin_center(Normal, Delta_New, d);
				}
				fracl = Al/(fabs(ul)*ns->par.dt*cdelta[1]);
			}
		}
		
		// Fraction
		fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
		real fracvolaux;
		fraction_correction_at_get(&fracr);
		fraction_correction_at_get(&fracl);
		
		fracvolaux = fracvol*(1 + ns->par.dt*(ur - ul)/cdelta[0]) - ns->par.dt*(fracr*ur - fracl*ul)/cdelta[0];
		// Auxiliary fraction correction
		fraction_correction_at_set(&fracvolaux);
		
		dp_set_value(ns->ed.mult.dpfracvolaux, clid, fracvolaux);
	}
	// Destroy the iterator
	higcit_destroy(it);
	// Sync the distributed vol frac aux property
	dp_sync(ns->ed.mult.dpfracvolaux);
}

// ***********************************************************************
// Volume Fraction Transport Step with PLIC fractional step on direction y
// ***********************************************************************
void higflow_plic_advection_volume_fraction_y_direction_2D(higflow_solver *ns, int dim){
	// Get the local sub-domain for the cells
	sim_domain *sdp = psd_get_local_domain(ns->psdp);
	// Get the local sub-domain for the facets
	sim_facet_domain *sfdu[DIM];
	for(int i = 0; i < DIM; i++) {
		sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
	}
	// Get the map for the domain properties
	mp_mapper *mp = sd_get_domain_mapper(sdp);
	// Loop for each cell
	higcit_celliterator *it;
	real tol_u = 1.0e-8;
	for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
		// Get the cell
		hig_cell *c = higcit_getcell(it);
		// Get the cell identifier
		int clid    = mp_lookup(mp, hig_get_cid(c));
		// Get the center of the cell
		Point ccenter;
		hig_get_center(c, ccenter);
		// Get the delta of the cell
		Point cdelta;
		hig_get_delta(c, cdelta);
		// Get the velocity at facet
		int infacet;
//			// Get the velocity in the left facet center
		real ul = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
//			// Get the velocity in the right facet center
		real ur = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);

		
		//real ur, ul;
		
		////===== Zalesak' test===========================================
		//ur = ccenter[0] - 0.5;
		//ul = ur;
		////==============================================================
		
		//real pi = 3.14159265358979323846;
		//real T = 8.0;
		////===== Shearing flow===========================================
		//ul = pow(sin(pi*ccenter[1]) ,2)*sin(2.0*pi*ccenter[0])*cos(pi*ns->par.t/T);
		//ur = ul;
		////==============================================================
		
		
		Point Normal, p, Delta_New;
		real d, fracvol,fracr,fracl;
		
		p[0]=ccenter[0];p[1]=ccenter[1];
		p[dim]=p[dim]+0.5*cdelta[dim];
		fracr=compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
		
		p[0]=ccenter[0];p[1]=ccenter[1];
		p[dim]=p[dim]-0.5*cdelta[dim];
		fracl=compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
		
		
		fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
		
		fraction_correction_at_get(&fracr);
		fraction_correction_at_get(&fracl);
		fraction_correction_at_get(&fracvol);

		// Area
		real Ar = 0.0;
		real Al = 0.0;
		fracr = 0.0;
		fracl = 0.0;
		//  Right Facet (UP)
		if(fabs(ur)>tol_u){
			if (fabs(ur)*ns->par.dt > 0.5*cdelta[0]) {
				printf("Time step is large!!!\n");
				exit(1);
			}
			if(ur>0.0){
				p[0] = ccenter[0];
				p[1] = ccenter[1];
				// Normal
				Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
				Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
				// Correction of Normal
				normal_correction_at_get(Normal);
				// Fraction
				fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
				// Fraction correction
				fraction_correction_at_get(&fracvol);
				if(fracvol==1.0 || fracvol==0.0 || (Normal[0]==0.0 && Normal[1]==0.0)){
					Ar = fracvol*cdelta[0]*ur*ns->par.dt;
				} else {
					d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
					d = d - 0.5*(cdelta[1]-ur*ns->par.dt)*Normal[1];
					Delta_New[0] = cdelta[0];
					Delta_New[1] = ur*ns->par.dt;
					Ar = area_left_line_origin_center(Normal, Delta_New, d);
				}
				fracr = Ar/(fabs(ur)*ns->par.dt*cdelta[0]);
			} else {
				p[0] = ccenter[0];
				p[1] = ccenter[1]+cdelta[1];
				// Normal
				Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
				Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
				// Correction of Normal
				normal_correction_at_get(Normal);
				// Fraction
				fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
				// Fraction correction
				fraction_correction_at_get(&fracvol);
				if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
					Ar = fracvol * cdelta[0] * fabs(ur) * ns->par.dt;
				} else {
					d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
					d = d + 0.5 * (cdelta[1] - fabs(ur) * ns->par.dt)*Normal[1];
					Delta_New[0] = cdelta[0];
					Delta_New[1] = fabs(ur)*ns->par.dt;
					Ar = area_left_line_origin_center(Normal, Delta_New, d);
				}
				fracr = Ar/(fabs(ur)*ns->par.dt*cdelta[0]);
			}
		}
		//  Left Facet (DWON)
		if (fabs(ul) > tol_u) {
			if (fabs(ur)*ns->par.dt > 0.5*cdelta[0]) {
				printf("Time step is large!!!\n");
				exit(1);
			}
			if (ul > 0.0) {
				p[0] = ccenter[0];
				p[1] = ccenter[1]-cdelta[1];
				// Normal
				Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
				Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
				// Correction of Normal
				normal_correction_at_get(Normal);
				// Fraction
				fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
				// Fraction correction
				fraction_correction_at_get(&fracvol);
				if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
					Al = fracvol * cdelta[0] * ul * ns->par.dt;
				} else {
					d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
					d = d - 0.5 * (cdelta[1] - ul * ns->par.dt) * Normal[1];
					Delta_New[0] = cdelta[0];
					Delta_New[1] = ul*ns->par.dt;
					Al = area_left_line_origin_center(Normal, Delta_New, d);
				}
				fracl = Al/(fabs(ul)*ns->par.dt*cdelta[0]);
			} else {
				p[0] = ccenter[0];
				p[1] = ccenter[1];
				// Normal
				Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
				Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
				// Correction of Normal
				normal_correction_at_get(Normal);
				// Fraction
				fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
				fraction_correction_at_get(&fracvol);
				if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
					Al = fracvol * cdelta[0] * fabs(ul) * ns->par.dt;
				} else {
					d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
					d = d + 0.5 * (cdelta[1] - fabs(ul) * ns->par.dt) * Normal[1];
					Delta_New[0] = cdelta[0];
					Delta_New[1] = fabs(ul)*ns->par.dt;
					Al = area_left_line_origin_center(Normal, Delta_New, d);
				}
				fracl = Al/(fabs(ul)*ns->par.dt*cdelta[0]);
			}
		}
		
		// Fraction
		fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
		real fracvolaux;
		fraction_correction_at_get(&fracr);
		fraction_correction_at_get(&fracl);
		
		fracvolaux = fracvol*(1 + ns->par.dt*(ur - ul)/cdelta[1]) - ns->par.dt*(fracr*ur - fracl*ul)/cdelta[1];
		// Auxiliary fraction correction
		fraction_correction_at_set(&fracvolaux);
		
		dp_set_value(ns->ed.mult.dpfracvolaux, clid, fracvolaux);
	}
	// Destroy the iterator
	higcit_destroy(it);
	// Sync the distributed vol frac aux property
	dp_sync(ns->ed.mult.dpfracvolaux);
}

// ***********************************************************************

// ***********************************************************************
// Volume Fraction Transport Step with PLIC fractional step on direction x
// ***********************************************************************
void higflow_plic_advection_volume_fraction_x_direction_imp_2D(higflow_solver *ns, int dim) {
	// Get the local sub-domain for the cells
	sim_domain *sdp = psd_get_local_domain(ns->psdp);
	// Get the local sub-domain for the facets
	sim_facet_domain *sfdu[DIM];
	for(int i = 0; i < DIM; i++) {
		sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
	}
	// Get the map for the domain properties
	mp_mapper *mp = sd_get_domain_mapper(sdp);
	// Loop for each cell
	higcit_celliterator *it;
	real tol_u = 1.0e-8;
	for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
		// Get the cell
		hig_cell *c = higcit_getcell(it);
		// Get the cell identifier
		int clid    = mp_lookup(mp, hig_get_cid(c));
		// Get the center of the cell
		Point ccenter;
		hig_get_center(c, ccenter);
		// Get the delta of the cell
		Point cdelta;
		hig_get_delta(c, cdelta);
		// Get the velocity at facet
		int infacet;
//			// Get the velocity in the left facet center
		real ul = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
//			// Get the velocity in the right facet center
		real ur = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);

		
		//real ur, ul;
		
		////===== Zalesak' test===========================================
		//ur = 0.5 - ccenter[1];
		//ul = ur;
		////===============================================================
		
		//real pi = 3.14159265358979323846;
		//real T = 8.0;
		////===== Shearing flow===========================================
		//ul = -pow(sin(pi*ccenter[0]) ,2)*sin(2.0*pi*ccenter[1])*cos(pi*ns->par.t/T);
		//ur = ul;
		////==============================================================
		

		Point Normal, p, Delta_New;
		real d, fracvol,fracr,fracl;

		p[0]=ccenter[0];p[1]=ccenter[1];
		p[dim]=p[dim]+0.5*cdelta[dim];
		fracr=compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);

		p[0]=ccenter[0];p[1]=ccenter[1];
		p[dim]=p[dim]-0.5*cdelta[dim];
		fracl=compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);


		fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);

		fraction_correction_at_get(&fracr);
		fraction_correction_at_get(&fracl);
		fraction_correction_at_get(&fracvol);


		// Area
		real Ar = 0.0;
		real Al = 0.0;
		fracr = 0.0;
		fracl = 0.0;
		//  Right Facet
		if(fabs(ur)>tol_u){
			if (fabs(ur) * ns->par.dt > 0.5 * cdelta[0]) {
				printf("Time step is large!!!\n");
				exit(1);
			}
			if(ur>0.0){
				p[0] = ccenter[0];
				p[1] = ccenter[1];
				// Normal
				Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
				Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
				// Correction of Normal
				normal_correction_at_get(Normal);
				// Fraction
				fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
				// Fraction correction
				fraction_correction_at_get(&fracvol);
				if(fracvol==1.0 || fracvol==0.0 || (Normal[0]==0.0 && Normal[1]==0.0)){
					Ar = fracvol*cdelta[1]*ur*ns->par.dt;
				} else {
					d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
					d = d - 0.5*(cdelta[0]-ur*ns->par.dt)*Normal[0];
					Delta_New[0] = ur*ns->par.dt;
					Delta_New[1] = cdelta[1];
					Ar = area_left_line_origin_center(Normal, Delta_New, d);
				}
				fracr = Ar/(fabs(ur)*ns->par.dt*cdelta[1]);
			} else {
				p[0] = ccenter[0]+cdelta[0];
				p[1] = ccenter[1];
				// Normal
				Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
				Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
				// Correction of Normal
				normal_correction_at_get(Normal);
				// Fraction
				fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
				// Fraction correction
				fraction_correction_at_get(&fracvol);
				if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
					Ar = fracvol * cdelta[1] * fabs(ur) * ns->par.dt;
				} else {
					d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
					Delta_New[0] = fabs(ur)*ns->par.dt;
					Delta_New[1] = cdelta[1];
					d = d + 0.5 * (cdelta[0] - fabs(ur) * ns->par.dt)*Normal[0];
					Ar = area_left_line_origin_center(Normal, Delta_New, d);
				}
				fracr = Ar/(fabs(ur)*ns->par.dt*cdelta[1]);
			}
		}
		//  Left Facet
		if (fabs(ul) > tol_u) {
			if (fabs(ul) * ns->par.dt > 0.5 * cdelta[0]) {
				printf("Time step is large!!!\n");
				exit(1);
			}
			if (ul > 0.0) {
				p[0] = ccenter[0] - cdelta[0];
				p[1] = ccenter[1];
				// Normal
				Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
				Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
				// Correction of Normal
				normal_correction_at_get(Normal);
				// Fraction
				fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
				// Fraction correction
				fraction_correction_at_get(&fracvol);
				if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
					Al = fracvol * cdelta[1] * ul * ns->par.dt;
				} else {
					d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
					d = d - 0.5 * (cdelta[0] - ul * ns->par.dt) * Normal[0];
					Delta_New[0] = ul*ns->par.dt;
					Delta_New[1] = cdelta[1];
					Al = area_left_line_origin_center(Normal, Delta_New, d);
				}
				fracl = Al/(fabs(ul)*ns->par.dt*cdelta[1]);
			} else {
				p[0] = ccenter[0];
				p[1] = ccenter[1];
				// Normal
				Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
				Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
				// Correction of Normal
				normal_correction_at_get(Normal);
				// Fraction
				fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
				fraction_correction_at_get(&fracvol);
				if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
					Al = fracvol * cdelta[1] * fabs(ul) * ns->par.dt;
				} else {
					d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
					d = d + 0.5 * (cdelta[0] - fabs(ul) * ns->par.dt) * Normal[0];
					Delta_New[0] = fabs(ul)*ns->par.dt;
					Delta_New[1] = cdelta[1];
					Al = area_left_line_origin_center(Normal, Delta_New, d);
				}
				fracl = Al/(fabs(ul)*ns->par.dt*cdelta[1]);
			}
		}

		// Fraction
		fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);

		real fracvolaux;
		fraction_correction_at_get(&fracr);
		fraction_correction_at_get(&fracl);
		
		fracvolaux = (fracvol - ns->par.dt*(fracr*ur - fracl*ul)/cdelta[0])/(1.0 - ns->par.dt*(ur - ul)/cdelta[0]);
		// Auxiliary fraction correction
		fraction_correction_at_set(&fracvolaux);

		dp_set_value(ns->ed.mult.dpfracvolaux, clid, fracvolaux);
	}
	// Destroy the iterator
	higcit_destroy(it);
	// Sync the distributed vol frac aux property
	dp_sync(ns->ed.mult.dpfracvolaux);
}

// ***********************************************************************
// Volume Fraction Transport Step with PLIC fractional step on direction y
// ***********************************************************************
void higflow_plic_advection_volume_fraction_y_direction_imp_2D(higflow_solver *ns, int dim){
	// Get the local sub-domain for the cells
	sim_domain *sdp = psd_get_local_domain(ns->psdp);
	// Get the local sub-domain for the facets
	sim_facet_domain *sfdu[DIM];
	for(int i = 0; i < DIM; i++) {
		sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
	}
	// Get the map for the domain properties
	mp_mapper *mp = sd_get_domain_mapper(sdp);
	// Loop for each cell
	higcit_celliterator *it;
	real tol_u = 1.0e-8;
	for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
		// Get the cell
		hig_cell *c = higcit_getcell(it);
		// Get the cell identifier
		int clid    = mp_lookup(mp, hig_get_cid(c));
		// Get the center of the cell
		Point ccenter;
		hig_get_center(c, ccenter);
		// Get the delta of the cell
		Point cdelta;
		hig_get_delta(c, cdelta);
		// Get the velocity at facet
		int infacet;
//			// Get the velocity in the left facet center
		real ul = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
//			// Get the velocity in the right facet center
		real ur = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);

		
		//real ur, ul;
		
		////===== Zalesak' test===========================================
		//ur = ccenter[0] - 0.5;
		//ul = ur;
		////==============================================================
		
		//real pi = 3.14159265358979323846;
		//real T = 8.0;
		////===== Shearing flow===========================================
		//ul = pow(sin(pi*ccenter[1]) ,2)*sin(2.0*pi*ccenter[0])*cos(pi*ns->par.t/T);
		//ur = ul;
		////==============================================================
		

		//arquivoV(nome_vy,ccenter[0],ccenter[1],ul,ur);

		Point Normal, p, Delta_New;
		real d, fracvol,fracr,fracl;

		p[0]=ccenter[0];p[1]=ccenter[1];
		p[dim]=p[dim]+0.5*cdelta[dim];
		fracr=compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);

		p[0]=ccenter[0];p[1]=ccenter[1];
		p[dim]=p[dim]-0.5*cdelta[dim];
		fracl=compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);


		fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);

		fraction_correction_at_get(&fracr);
		fraction_correction_at_get(&fracl);
		fraction_correction_at_get(&fracvol);

		// Area
		real Ar = 0.0;
		real Al = 0.0;
		fracr = 0.0;
		fracl = 0.0;
		//  Right Facet (UP)
		if(fabs(ur)>tol_u){
			if (fabs(ur)*ns->par.dt > 0.5*cdelta[0]) {
				printf("Time step is large!!!\n");
				exit(1);
			}
			if(ur>0.0){
				p[0] = ccenter[0];
				p[1] = ccenter[1];
				// Normal
				Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
				Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
				// Correction of Normal
				normal_correction_at_get(Normal);
				// Fraction
				fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
				// Fraction correction
				fraction_correction_at_get(&fracvol);
				if(fracvol==1.0 || fracvol==0.0 || (Normal[0]==0.0 && Normal[1]==0.0)){
					Ar = fracvol*cdelta[0]*ur*ns->par.dt;
				} else {
					d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
					d = d - 0.5*(cdelta[1]-ur*ns->par.dt)*Normal[1];
					Delta_New[0] = cdelta[0];
					Delta_New[1] = ur*ns->par.dt;
					Ar = area_left_line_origin_center(Normal, Delta_New, d);
				}
				fracr = Ar/(fabs(ur)*ns->par.dt*cdelta[0]);
			} else {
				p[0] = ccenter[0];
				p[1] = ccenter[1]+cdelta[1];
				// Normal
				Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
				Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
				// Correction of Normal
				normal_correction_at_get(Normal);
				// Fraction
				fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
				// Fraction correction
				fraction_correction_at_get(&fracvol);
				if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
					Ar = fracvol * cdelta[0] * fabs(ur) * ns->par.dt;
				} else {
					d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
					d = d + 0.5 * (cdelta[1] - fabs(ur) * ns->par.dt)*Normal[1];
					Delta_New[0] = cdelta[0];
					Delta_New[1] = fabs(ur)*ns->par.dt;
					Ar = area_left_line_origin_center(Normal, Delta_New, d);
				}
				fracr = Ar/(fabs(ur)*ns->par.dt*cdelta[0]);
			}
		}
		//  Left Facet (DWON)
		if (fabs(ul) > tol_u) {
			if (fabs(ur)*ns->par.dt > 0.5*cdelta[0]) {
				printf("Time step is large!!!\n");
				exit(1);
			}
			if (ul > 0.0) {
				p[0] = ccenter[0];
				p[1] = ccenter[1]-cdelta[1];
				// Normal
				Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
				Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
				// Correction of Normal
				normal_correction_at_get(Normal);
				// Fraction
				fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
				// Fraction correction
				fraction_correction_at_get(&fracvol);
				if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
					Al = fracvol * cdelta[0] * ul * ns->par.dt;
				} else {
					d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
					d = d - 0.5 * (cdelta[1] - ul * ns->par.dt) * Normal[1];
					Delta_New[0] = cdelta[0];
					Delta_New[1] = ul*ns->par.dt;
					Al = area_left_line_origin_center(Normal, Delta_New, d);
				}
				fracl = Al/(fabs(ul)*ns->par.dt*cdelta[0]);
			} else {
				p[0] = ccenter[0];
				p[1] = ccenter[1];
				// Normal
				Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
				Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
				// Correction of Normal
				normal_correction_at_get(Normal);
				// Fraction
				fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
				fraction_correction_at_get(&fracvol);
				if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
					Al = fracvol * cdelta[0] * fabs(ul) * ns->par.dt;
				} else {
					d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
					d = d + 0.5 * (cdelta[1] - fabs(ul) * ns->par.dt) * Normal[1];
					Delta_New[0] = cdelta[0];
					Delta_New[1] = fabs(ul)*ns->par.dt;
					Al = area_left_line_origin_center(Normal, Delta_New, d);
				}
				fracl = Al/(fabs(ul)*ns->par.dt*cdelta[0]);
			}
		}

		// Fraction
		fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
		real fracvolaux;
		fraction_correction_at_get(&fracr);
		fraction_correction_at_get(&fracl);

		fracvolaux = (fracvol - ns->par.dt*(fracr*ur - fracl*ul)/cdelta[1])/(1.0 - ns->par.dt*(ur - ul)/cdelta[1]);
		// Auxiliary fraction correction
		fraction_correction_at_set(&fracvolaux);

		dp_set_value(ns->ed.mult.dpfracvolaux, clid, fracvolaux);
	}
	// Destroy the iterator
	higcit_destroy(it);
	// Sync the distributed vol frac aux property
	dp_sync(ns->ed.mult.dpfracvolaux);
}
//************************************************************************

// ***********************************************************************
//Copy Auxiliary Volume Fraction to Volume Fraction
// ***********************************************************************
void higflow_plic_copy_fractionaux_to_fraction(higflow_solver *ns) {
    if (ns->contr.flowtype == 2) {
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->psdp);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int i = 0; i < DIM; i++) {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdp);
        // Loop for each cell
        higcit_celliterator *it;

        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid    = mp_lookup(mp, hig_get_cid(c));
            // Get the center of the cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            real fracvolaux  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvolaux, ns->ed.stn);
            dp_set_value(ns->ed.mult.dpfracvol, clid, fracvolaux);
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the distributed pressure property
        dp_sync(ns->ed.mult.dpfracvol);
        dp_sync(ns->ed.mult.dpfracvolaux);
  }
}

void fraction_correction_at_set(real *frac){
	
	if(fabs(*frac - 1.0)<1.0e-8){
		*frac = 1.0;
	}else if(fabs(*frac)<1.0e-8){
		*frac=0.0;
	}
	
	//return;
	
	if(*frac>1.0){
		*frac=1.0;
	}else if(*frac<0.0){
		*frac=0.0;
	}
}

// Correction Normal
void normal_correction_at_get(Point Normal){
	real tol_n = 1.0e-6;
	for(int i=0;i<DIM;i++){
		if(fabs(Normal[i])<tol_n) {
			Normal[i] = 0.0;
		}
	}
}
#endif
