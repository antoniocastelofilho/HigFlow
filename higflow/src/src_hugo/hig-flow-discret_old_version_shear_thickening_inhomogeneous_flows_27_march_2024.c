// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Discretization - version 10/11/2016
// *******************************************************************
// *******************************************************************

#include "hig-flow-discret.h"

// *******************************************************************
// Navier-Stokes Discretization
// *******************************************************************

// Compute the first derivative using the central method
real compute_facet_dudxc(Point delta, int dim, real alpha, real uc, real ul, real ur) {
	// Midpoint first derivative
	real dudx = (ur - ul) / (2.0 * alpha * delta[dim]);
	return dudx;
}

// Compute the first derivative using the forward method
real compute_facet_dudxr(Point delta, int dim, real alpha, real uc, real ul, real ur) {
	// Midpoint first derivative
	real dudx = (ur - uc) / (alpha * delta[dim]);
	return dudx;
}

// Compute the first derivative using the backward method
real compute_facet_dudxl(Point delta, int dim, real alpha, real uc, real ul, real ur) {
	// Midpoint first derivative
	real dudx = (uc - ul) / (alpha * delta[dim]);
	return dudx;
}

// Compute the second derivative using the forward method
real compute_facet_du2dx2r(Point delta, int dim, real alpha, real u, real ur, real urr) {
	// Midpoint second derivative
	real du2dx2 = (urr - 2.0*ur + u)/(delta[dim] * delta[dim] * alpha * alpha);
	return du2dx2;
}

// Compute the second derivative using the backward method
real compute_facet_du2dx2l(Point delta, int dim, real alpha, real u, real ul, real ull) {
	// Midpoint second derivative
	real du2dx2 = (u - 2.0*ul + ull)/(delta[dim] * delta[dim] * alpha * alpha);
	return du2dx2;
}

// Compute the second derivative using the central method
real compute_facet_du2dx2(Point delta, int dim, real alpha, real u, real ul, real ur) {
	// Midpoint second derivative
	real du2dx2 = (ul - 2.0*u + ur)/(delta[dim] * delta[dim] * alpha * alpha);
	return du2dx2;
}

// Compute the first derivative using the central method
real compute_dpdx_at_point(Point delta, int dim, real alpha, real valuel, real valueh) {
	// Midpoint first derivative
	real dpdx = (valueh - valuel) / (2.0 * alpha * delta[dim]);
	return dpdx;
}

// Compute the first derivative using the forward method
real compute_dpdxr_at_point(Point delta, int dim, real alpha, real valuec, real valueh) {
	// Midpoint first derivative
	real dpdx = (valueh - valuec) / (alpha * delta[dim]);
	return dpdx;
}

// Compute the first derivative using the backward method
real compute_dpdxl_at_point(Point delta, int dim, real alpha, real valuel, real valuec) {
	// Midpoint first derivative
	real dpdx = (valuec - valuel) / (alpha * delta[dim]);
	return dpdx;
}

// Compute the mid point value
real compute_value_at_mid_point(real valuel, real valuer) {
	// Midpoint value
	real value = (valuer + valuel) / 2.0;
	return value;
}

// *******************************************************************
// Navier-Stokes Computational Cell
// *******************************************************************

// Computing the necessary term for the Navier-Stokes equation
void higflow_computational_cell(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]) {
	// Set the computational cell
	real  uc[DIM], ul[DIM], ur[DIM], Sl[DIM], Sr[DIM], pl, pr, Fl, Fr;
	real  ull[DIM], urr[DIM];
	real  a, b, c, d, e, tol, fi, conv1,conv2;
	a     = 1.7500;
	b     = 0.3750;
	c     = 0.7500;
	d     = 0.1250;
	e     = 0.2500;
	tol   = 1.0e-14;
	conv1 = 0.0;
	conv2 = 0.0;
	int   infacet, infacet_r, infacet_l, infacet_rr, infacet_ll;
	// Spatial discretization
	switch (ns->contr.spatialdiscrtype) {
	// Second order
	case 0:
		// Set the convective method
		ns->cc.convec_type = ns->contr.convecdiscrtype;
		// Get the velocity in the facet center
		ns->cc.ucell  = dp_get_value(dpu[dim], fgid);
		// Get the facet source tern in the facet center
		ns->cc.F      = compute_facet_value_at_point(ns->sfdF[dim], fcenter, fcenter, 1.0, ns->dpFU[dim], ns->stnF);
		// Get the pressure in the left cell
		pl            = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Get the pressure in the right cell
		pr            = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Compute the pressure derivative
		ns->cc.dpdx   = compute_dpdx_at_point(fdelta, dim, 0.5, pl, pr);
		// Get the cell source term in the left cell
		Fl            = compute_center_p_left(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Get the cell source term in the left cell
		Fr            = compute_center_p_right(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Compute the source term at the mid point
		ns->cc.F     += compute_value_at_mid_point(Fl, Fr);
		for (int dim2 = 0; dim2 < DIM; dim2++) {
			// Get the velocity in the center of facet
			if (dim2 == dim) {
				ns->cc.v[dim2] = ns->cc.ucell;
			} else {
				ns->cc.v[dim2] = compute_facet_u_2(sfdu[dim2], fcenter, fdelta, dim, dim2, 1.0, dpu[dim2], ns->stn);
			}
			// Get the velocity in the left facet center
			ul[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet_l);
			// Get the velocity in the right facet center
			ur[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet_r);
			// Compute dudxc at facet center
			ns->cc.dudxc[dim2]  = compute_facet_dudxc(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxr at facet center
			ns->cc.dudxr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxl at facet center
			ns->cc.dudxl[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet_l == 0)){// && (ns->cc.convec_type == 2)) {
				ul[dim2] = compute_facet_u_2_left(sfdu[dim], fcenter, fdelta, dim, dim2, 1.0, dpu[dim2], ns->stn);
				ns->cc.dudxl[dim2]  = compute_facet_dudxl(fdelta, dim2, 0.5, ns->cc.ucell, ul[dim2], ur[dim2]);
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				//ns->cc.convec_type = 1;
			}
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet_r == 0)){// && (ns->cc.convec_type == 2)) {
				ur[dim2] = compute_facet_u_2_right(sfdu[dim], fcenter, fdelta, dim, dim2, 1.0, dpu[dim2], ns->stn);
				//Trocada pela linha de baixo --- ns->cc.dudxr[dim2]  = compute_facet_dudxl(fdelta, dim2, 0.5, ns->cc.ucell, ul[dim2], ur[dim2]);
				ns->cc.dudxr[dim2]  = compute_facet_dudxr(fdelta, dim2, 0.5, ns->cc.ucell, ul[dim2], ur[dim2]);
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				//ns->cc.convec_type = 1;
			}
			// Get the velocity in the right facet center
			urr[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			// Get the velocity in the left facet center
			ull[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			if (ns->cc.convec_type == 2) {
				if (ns->contr.secondconvecdiscrtype == 2) {
					// Compute dudxrr at facet center
					ns->cc.dudxrr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ur[dim2], ns->cc.ucell, urr[dim2]);
					// Compute dudxll at facet center
					ns->cc.dudxll[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ul[dim2], ull[dim2], ns->cc.ucell);
					// Compute terms for second order upwind CUBISTA
				}else if (ns->contr.secondconvecdiscrtype == 1){
					// Get the velocity in the left facet center
					ul[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the left facet center
					ull[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_ll);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity in the right facet center
					ur[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the right facet center
					urr[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_rr);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity v1bar(i+1/2,j+1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_right(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ur[dim]-ul[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv1 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ul[dim]);
								if ((fi >= b) && (fi <= c))
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + b*ur[dim] -d*ul[dim]);
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(e*ns->cc.ucell + c*ur[dim]);
							}
						}
						//v1bar < 0.0
					}else {
						if (fabs(ns->cc.ucell - urr[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ur[dim];
						} else {
							//Trocada pela linha de baixo --- fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							fi = (ur[dim]-urr[dim])/(ns->cc.ucell -urr[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ur[dim];
							}else {
								if (fi < b)
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(a*ur[dim] - c*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if ((fi >= b) && (fi <= c))
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(c*ur[dim] + b*ns->cc.ucell -d*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ur[dim]);
							}
						}
					}
					// Get the velocity  v2bar(i+1/2,j-1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_left(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ns->cc.ucell-ull[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ul[dim];
						}else {
							fi = (ul[dim] - ull[dim])/(ns->cc.ucell-ull[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ul[dim];
							}else {
								if (fi < b)
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(a*ul[dim] - c*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if ((fi >= b) && (fi <= c))
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(b*ns->cc.ucell + c*ul[dim] - d*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ul[dim]);
							}
						}
					}else {
						//v2bar < 0.0
						if (fabs(ul[dim] - ur[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ur[dim])/(ul[dim] - ur[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv2 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ur[dim]);
								if ((fi >= b) && (fi <= c))
									conv2 = ns->cc.vc[dim2]*(b*ul[dim] + c*ns->cc.ucell - d*ur[dim]);
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ul[dim] + e*ns->cc.ucell);
							}
						}
					}
					ns->cc.vc[dim2] = ((conv1-conv2)/fdelta[dim]);
				}
			}
			// Newtonian
			ns->cc.du2dx2[dim2] = compute_facet_du2dx2(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
		}
		break;
	}
}


// Computing the necessary term for the Generalized Newtonian Navier-Stokes equation
void higflow_computational_cell_gen_newt(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]) {
	// Set the computational cell
	real  uc[DIM], ul[DIM], ur[DIM], Sl[DIM], Sr[DIM], pl, pr, Fl, Fr;
	real  ull[DIM], urr[DIM];
	real  a, b, c, d, e, tol, fi, conv1,conv2;
	a     = 1.7500;
	b     = 0.3750;
	c     = 0.7500;
	d     = 0.1250;
	e     = 0.2500;
	tol   = 1.0e-14;
	conv1 = 0.0;
	conv2 = 0.0;
	int   infacet, infacet_rr, infacet_ll;
	// Spatial discretization
	switch (ns->contr.spatialdiscrtype) {
	// Second order
	case 0:
		// Set the convective method
		ns->cc.convec_type = ns->contr.convecdiscrtype;
		// Get the velocity in the facet center
		ns->cc.ucell  = dp_get_value(dpu[dim], fgid);
		// Get the facet source tern in the facet center
		ns->cc.F      = compute_facet_value_at_point(ns->sfdF[dim], fcenter, fcenter, 1.0, ns->dpFU[dim], ns->stnF);
		// Get the pressure in the left cell
		pl            = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Get the pressure in the right cell
		pr            = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Compute the pressure derivative
		ns->cc.dpdx   = compute_dpdx_at_point(fdelta, dim, 0.5, pl, pr);
		// Get the cell source term in the left cell
		Fl            = compute_center_p_left(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Get the cell source term in the left cell
		Fr            = compute_center_p_right(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Compute the source term at the mid point
		ns->cc.F     += compute_value_at_mid_point(Fl, Fr);
		// Get the cell viscolity in the left cell
		ns->cc.viscl = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.gn.dpvisc, ns->ed.stn);
		// Get the cell viscolity in the right cell
		ns->cc.viscr = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.gn.dpvisc, ns->ed.stn);
		for (int dim2 = 0; dim2 < DIM; dim2++) {
			// Get the velocity in the center of facet
			if (dim2 == dim) {
				ns->cc.v[dim2] = ns->cc.ucell;
			} else {
				ns->cc.v[dim2] = compute_facet_u_2(sfdu[dim2], fcenter, fdelta, dim, dim2, 1.0, dpu[dim2], ns->stn);
			}
			// Get the velocity in the left facet center
			ul[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Get the velocity in the right facet center
			ur[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Compute dudxc at facet center
			ns->cc.dudxc[dim2]  = compute_facet_dudxc(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxr at facet center
			ns->cc.dudxr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxl at facet center
			ns->cc.dudxl[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Get the velocity in the right facet center
			urr[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			// Get the velocity in the left facet center
			ull[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			if (ns->cc.convec_type == 2) {
				if (ns->contr.secondconvecdiscrtype == 2) {
					// Compute dudxrr at facet center
					ns->cc.dudxrr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ur[dim2], ns->cc.ucell, urr[dim2]);
					// Compute dudxll at facet center
					ns->cc.dudxll[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ul[dim2], ull[dim2], ns->cc.ucell);
					// Compute terms for second order upwind CUBISTA
				}
				else if (ns->contr.secondconvecdiscrtype == 1){
					// Get the velocity in the left facet center
					ul[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the left facet center
					ull[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_ll);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity in the right facet center
					ur[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the right facet center
					urr[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_rr);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity v1bar(i+1/2,j+1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_right(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ur[dim]-ul[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv1 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ul[dim]);
								if ((fi >= b) && (fi <= c))
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + b*ur[dim] -d*ul[dim]);
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(e*ns->cc.ucell + c*ur[dim]);
							}
						}
						//v1bar < 0.0
					}else {
						if (fabs(ns->cc.ucell - urr[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ur[dim];
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ur[dim];
							}else {
								if (fi < b)
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(a*ur[dim] - c*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if ((fi >= b) && (fi <= c))
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(c*ur[dim] + b*ns->cc.ucell -d*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ur[dim]);
							}
						}
					}
					// Get the velocity  v2bar(i+1/2,j-1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_left(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ns->cc.ucell-ull[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ul[dim];
						}else {
							fi = (ul[dim] - ull[dim])/(ns->cc.ucell-ull[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ul[dim];
							}else {
								if (fi < b)
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(a*ul[dim] - c*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if ((fi >= b) && (fi <= c))
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(b*ns->cc.ucell + c*ul[dim] - d*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ul[dim]);
							}
						}
					}else {
						//v2bar < 0.0
						if (fabs(ul[dim] - ur[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ur[dim])/(ul[dim] - ur[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv2 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ur[dim]);
								if ((fi >= b) && (fi <= c))
									conv2 = ns->cc.vc[dim2]*(b*ul[dim] + c*ns->cc.ucell - d*ur[dim]);
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ul[dim] + e*ns->cc.ucell);
							}
						}
					}
					ns->cc.vc[dim2] = ((conv1-conv2)/fdelta[dim]);
				}
			}
			// Generalized Newtonian
			ns->cc.du2dx2[dim2] = 0.0;
			Point p;
			POINT_ASSIGN(p, fcenter);

			p[dim2]      = fcenter[dim2] + 0.5*fdelta[dim2];

			real ucl     = compute_facet_u_left(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet);
			real ucr     = compute_facet_u_right(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet);
			real duidxjr = compute_facet_dudxc(fdelta, dim2, 0.5, ns->cc.ucell, ucl, ucr);

			real vcl     = compute_facet_u_left(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet);
			real vcr     = compute_facet_u_right(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet);
			real dujdxir = compute_facet_dudxc(fdelta, dim, 0.5, ns->cc.ucell, vcl, vcr);

			p[dim2]      =  fcenter[dim2] - 0.5*fdelta[dim2];

			ucl          = compute_facet_u_left(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet);
			ucr          = compute_facet_u_right(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet);
			real duidxjl = compute_facet_dudxc(fdelta, dim2, 0.5, ns->cc.ucell, ucl, ucr);

			vcl          = compute_facet_u_left(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet);
			vcr          = compute_facet_u_right(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet);
			real dujdxil = compute_facet_dudxc(fdelta, dim, 0.5, ns->cc.ucell, vcl, vcr);

			ns->cc.du2dx2[dim2] += (ns->cc.viscr*duidxjr - ns->cc.viscl*duidxjl)/fdelta[dim];
			ns->cc.du2dx2[dim2] += (ns->cc.viscr*dujdxir - ns->cc.viscl*dujdxil)/fdelta[dim];


			real eps[DIM][DIM];
			//Levi-Civita matrix in 2D
			for (int i = 0; i < DIM; i++) {
				for (int j = i+1; j < DIM; j++) {
					eps[j][1] = 1.0;
					eps[1][j] = -1.0;
				}
				eps[i][i]  = 0.0;
			}

			//Odd Viscosity terms
			//if ((ns->contr.modelflowtype== 5)) {
			//	//We asume by now this value for the Odd Viscosity Oeta
			//	real Oeta = 1.0/ns->par.Re;
				//ns->cc.dSdx[dim2] = 0.0;
			//	Point po;
			//	POINT_ASSIGN(po, fcenter);
			//	po[dim2]      = fcenter[dim2] + 0.5*fdelta[dim2];
				//Calculation of the term (epsik)*dUjdxk
			//	real ucol, ucor, duidxjro;
			//	real epsikdujdxkr = 0.0;
			//	for (int dim3 = 0; dim3 < DIM; dim3++) {
			//		ucol = compute_facet_u_left(sfdu[dim2], po, fdelta, dim3, 0.5, dpu[dim2], ns->stn, &infacet);
			//		ucor = compute_facet_u_right(sfdu[dim2], po, fdelta, dim3, 0.5, dpu[dim2], ns->stn, &infacet);
			//		duidxjro = compute_facet_dudxc(fdelta, dim3, 0.5, ns->cc.ucell, ucol, ucor);
			//		epsikdujdxkr += eps[dim][dim3]*duidxjro;
			//	}

				//Calculation of the term d(epsjk)Uk/dxi
			//	real vcol, vcor, dukdxiro;
			//	real depsjkukdxir = 0.0;
			//	for (int dim3 = 0; dim3 < DIM; dim3++) {
			//		vcol = compute_facet_u_left(sfdu[dim3], po, fdelta, dim, 0.5, dpu[dim3], ns->stn, &infacet);
			//		vcor = compute_facet_u_right(sfdu[dim3], po, fdelta, dim, 0.5, dpu[dim3], ns->stn, &infacet);
			//		dukdxiro = compute_facet_dudxc(fdelta, dim, 0.5, ns->cc.ucell, vcol, vcor);
			//		depsjkukdxir += eps[dim2][dim3]*dukdxiro;
			//	}

			//	po[dim2]      =  fcenter[dim2] - 0.5*fdelta[dim2];
				//Calculation of the term (epsik)*dUjdxk
			//	real duidxjlo; 
			//	real epsikdujdxkl = 0.0;
			//	for (int dim3 = 0; dim3 < DIM; dim3++) {
			//		ucol = compute_facet_u_left(sfdu[dim2], po, fdelta, dim3, 0.5, dpu[dim2], ns->stn, &infacet);
			//		ucor = compute_facet_u_right(sfdu[dim2], po, fdelta, dim3, 0.5, dpu[dim2], ns->stn, &infacet);
			//		duidxjlo = compute_facet_dudxc(fdelta, dim3, 0.5, ns->cc.ucell, ucol, ucor);
			//		epsikdujdxkl += eps[dim][dim3]*duidxjlo;
			//	}

				//Calculation of the term d(epsjk)Uk/dxi
			//	real dukdxilo;
			//	real depsjkukdxil = 0.0;
			//	for (int dim3 = 0; dim3 < DIM; dim3++) {
			//		vcol = compute_facet_u_left(sfdu[dim3], po, fdelta, dim, 0.5, dpu[dim3], ns->stn, &infacet);
			//		vcor = compute_facet_u_right(sfdu[dim3], po, fdelta, dim, 0.5, dpu[dim3], ns->stn, &infacet);
			//		dukdxilo = compute_facet_dudxc(fdelta, dim, 0.5, ns->cc.ucell, vcol, vcor);
			//		depsjkukdxil += eps[dim2][dim3]*dukdxilo;
			//	}

			//	ns->cc.du2dx2[dim2] += ns->par.Re*(epsikdujdxkr - epsikdujdxkl)/(fdelta[dim]*Oeta);
			//	ns->cc.du2dx2[dim2] += ns->par.Re*(depsjkukdxir - depsjkukdxil)/(fdelta[dim]*Oeta);
				//ns->cc.dSdx[dim2] += (epsikdujdxkr - epsikdujdxkl)/(fdelta[dim]*Oeta);
				//ns->cc.dSdx[dim2] += (depsjkukdxir - depsjkukdxil)/(fdelta[dim]*Oeta);

			//}

		}
		break;
	}
}


// Computing the necessary term for the Implicit Generalized Newtonian Navier-Stokes equation
void higflow_computational_cell_imp_gen_newt(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]) {
	// Set the computational cell
	real  uc[DIM], ul[DIM], ur[DIM], Sl[DIM], Sr[DIM], pl, pr, Fl, Fr;
	real  ull[DIM], urr[DIM];
	real  a, b, c, d, e, tol, fi, conv1,conv2;
	a     = 1.7500;
	b     = 0.3750;
	c     = 0.7500;
	d     = 0.1250;
	e     = 0.2500;
	tol   = 1.0e-14;
	conv1 = 0.0;
	conv2 = 0.0;
	int   infacet, infacet_rr, infacet_ll;

	//Odd Viscosity terms
	//***************************************
	//real Oeta = ns->par.Re;
	//***************************************
	// Spatial discretization
	switch (ns->contr.spatialdiscrtype) {
	// Second order
	case 0:
		// Set the convective method
		ns->cc.convec_type = ns->contr.convecdiscrtype;
		// Get the velocity in the facet center
		ns->cc.ucell  = dp_get_value(dpu[dim], fgid);
		// Get the facet source tern in the facet center
		ns->cc.F      = compute_facet_value_at_point(ns->sfdF[dim], fcenter, fcenter, 1.0, ns->dpFU[dim], ns->stnF);
		// Get the pressure in the left cell
		pl            = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Get the pressure in the right cell
		pr            = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Compute the pressure derivative
		ns->cc.dpdx   = compute_dpdx_at_point(fdelta, dim, 0.5, pl, pr);
		// Get the cell source term in the left cell
		Fl            = compute_center_p_left(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Get the cell source term in the left cell
		Fr            = compute_center_p_right(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Compute the source term at the mid point
		ns->cc.F     += compute_value_at_mid_point(Fl, Fr);
		// Get the cell viscolity in the left cell
		ns->cc.viscl        = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.gn.dpvisc, ns->ed.stn);
		// Get the cell viscolity in the right cell
		ns->cc.viscr        = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.gn.dpvisc, ns->ed.stn);
		for (int dim2 = 0; dim2 < DIM; dim2++) {
			// Get the velocity in the center of facet
			if (dim2 == dim) {
				ns->cc.v[dim2] = ns->cc.ucell;
			} else {
				ns->cc.v[dim2] = compute_facet_u_2(sfdu[dim2], fcenter, fdelta, dim, dim2, 1.0, dpu[dim2], ns->stn);
			}
			// Get the velocity in the left facet center
			ul[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Get the velocity in the right facet center
			ur[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Compute dudxc at facet center
			ns->cc.dudxc[dim2]  = compute_facet_dudxc(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxr at facet center
			ns->cc.dudxr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxl at facet center
			ns->cc.dudxl[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Get the velocity in the right facet center
			urr[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			// Get the velocity in the left facet center
			ull[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			if (ns->cc.convec_type == 2) {
				if (ns->contr.secondconvecdiscrtype == 2) {
					// Compute dudxrr at facet center
					ns->cc.dudxrr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ur[dim2], ns->cc.ucell, urr[dim2]);
					// Compute dudxll at facet center
					ns->cc.dudxll[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ul[dim2], ull[dim2], ns->cc.ucell);
					// Compute terms for second order upwind CUBISTA
				}
				else if (ns->contr.secondconvecdiscrtype == 1){
					// Get the velocity in the left facet center
					ul[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the left facet center
					ull[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_ll);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity in the right facet center
					ur[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the right facet center
					urr[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_rr);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity v1bar(i+1/2,j+1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_right(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ur[dim]-ul[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv1 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ul[dim]);
								if ((fi >= b) && (fi <= c))
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + b*ur[dim] -d*ul[dim]);
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(e*ns->cc.ucell + c*ur[dim]);
							}
						}
						//v1bar < 0.0
					}else {
						if (fabs(ns->cc.ucell - urr[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ur[dim];
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ur[dim];
							}else {
								if (fi < b)
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(a*ur[dim] - c*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if ((fi >= b) && (fi <= c))
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(c*ur[dim] + b*ns->cc.ucell -d*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ur[dim]);
							}
						}
					}
					// Get the velocity  v2bar(i+1/2,j-1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_left(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ns->cc.ucell-ull[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ul[dim];
						}else {
							fi = (ul[dim] - ull[dim])/(ns->cc.ucell-ull[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ul[dim];
							}else {
								if (fi < b)
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(a*ul[dim] - c*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if ((fi >= b) && (fi <= c))
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(b*ns->cc.ucell + c*ul[dim] - d*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ul[dim]);
							}
						}
					}else {
						//v2bar < 0.0
						if (fabs(ul[dim] - ur[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ur[dim])/(ul[dim] - ur[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv2 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ur[dim]);
								if ((fi >= b) && (fi <= c))
									conv2 = ns->cc.vc[dim2]*(b*ul[dim] + c*ns->cc.ucell - d*ur[dim]);
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ul[dim] + e*ns->cc.ucell);
							}
						}
					}
					ns->cc.vc[dim2] = ((conv1-conv2)/fdelta[dim]);
				}
			}
			// Implicit Generalized Newtonian
			ns->cc.du2dx2[dim2] = 0.0;
			Point p;
			POINT_ASSIGN(p, fcenter);

			p[dim2]      = fcenter[dim2] + 0.5*fdelta[dim2];

			/*non-mixed*/
			real ucl     = compute_facet_u_left(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet);
			real ucr     = compute_facet_u_right(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet);
			real duidxjr = compute_facet_dudxc(fdelta, dim2, 0.5, ns->cc.ucell, ucl, ucr);

			/*mixed*/
			real vcl     = compute_facet_u_left(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet);
			real vcr     = compute_facet_u_right(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet);
			real dujdxir = compute_facet_dudxc(fdelta, dim, 0.5, ns->cc.ucell, vcl, vcr);

			p[dim2]      =  fcenter[dim2] - 0.5*fdelta[dim2];

			/*non-mixed*/
			ucl          = compute_facet_u_left(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet);
			ucr          = compute_facet_u_right(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet);
			real duidxjl = compute_facet_dudxc(fdelta, dim2, 0.5, ns->cc.ucell, ucl, ucr);

			/*mixed*/
			vcl          = compute_facet_u_left(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet);
			vcr          = compute_facet_u_right(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet);
			real dujdxil = compute_facet_dudxc(fdelta, dim, 0.5, ns->cc.ucell, vcl, vcr);

			//ns->cc.du2dx2[dim2] += (ns->cc.viscr*duidxjr - ns->cc.viscl*duidxjl)/fdelta[dim2];//ns->cc.du2dx2[dim2]=0.0;
			ns->cc.du2dx2[dim2] += (ns->cc.viscr*dujdxir - ns->cc.viscl*dujdxil)/fdelta[dim2];


			real eps[DIM][DIM];
			//Levi-Civita matrix in 2D
			for (int i = 0; i < DIM; i++) {
				for (int j = 0; j < DIM; j++) {
					eps[j][1] = 1.0;
					eps[1][j] = -1.0;
				}
				eps[i][i]  = 0.0;
			}
			//ns->cc.dSdx[dim2] = 0.0;
			//DEBUG_INSPECT(eps[1][1],%lf);


			//Odd Viscosity terms
			//if ((ns->contr.modelflowtype== 5)) {
				

				//DEBUG_INSPECT(TOVL[1][0],%lf);

				//->cc.dSdx[dim2] += 0;
				//ns->cc.du2dx2[dim2] += ns->par.Re*(epsikdujdxkr - epsikdujdxkl)/(fdelta[dim]*Oeta);
				//ns->cc.du2dx2[dim2] += ns->par.Re*(depsjkukdxir - depsjkukdxil)/(fdelta[dim]*Oeta);
				//ns->cc.dSdx[dim2] += (epsikdujdxkr[dim][dim2])/(fdelta[dim]*Oeta);
				//ns->cc.du2dx2[dim2] += ns->par.Re*(epsikdujdxkr[dim][dim2]-epsikdujdxkl[dim][dim2])/(fdelta[dim]*Oeta);
				//ns->cc.du2dx2[dim2] += ns->par.Re*(depsjkukdxir[dim2][dim]-depsjkukdxil[dim2][dim])/(fdelta[dim]*Oeta);
				//ns->cc.dSdx[dim2] += (epsikdujdxkr[dim][dim2]-epsikdujdxkl[dim][dim2])/(fdelta[dim]*Oeta);
				//ns->cc.dSdx[dim2] += (depsjkukdxir[dim2][dim]-depsjkukdxil[dim2][dim])/(fdelta[dim]*Oeta);
				//ns->cc.dSdx[dim2] += (epsikdujdxkr-epsikdujdxkl)/(fdelta[dim]*Oeta);
				//ns->cc.dSdx[dim2] += (depsjkukdxir-depsjkukdxil)/(fdelta[dim]*Oeta);
				//ns->cc.dSdx[dim2] += (epsikdujdxkr - epsikdujdxkl)/fdelta[dim];
				//ns->cc.dSdx[dim2] += (depsjkukdxir - depsjkukdxil)/fdelta[dim];
				//ns->cc.dSdx[dim2] = ns->par.Re*(epsikdujdxkr - epsikdujdxkl + depsjkukdxir - depsjkukdxil)/(fdelta[dim]*Oeta);
				//ns->cc.dSdx[dim2] += 0.0;

			//}



		}
		break;
	}
}


// Computing the necessary term for the Multifase Navier-Stokes equation
void higflow_computational_cell_multiphase(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]) {
	// Set the computational cell
	real  uc[DIM], ul[DIM], ur[DIM], Sl[DIM], Sr[DIM], pl, pr, Fl, Fr;
	real  ull[DIM], urr[DIM];
	real  rhol, rhor;
	real  a, b, c, d, e, tol, fi, conv1, conv2, curvl, curvr;
	a     = 1.7500;
	b     = 0.3750;
	c     = 0.7500;
	d     = 0.1250;
	e     = 0.2500;
	tol   = 1.0e-14;
	conv1 = 0.0;
	conv2 = 0.0;
	int   infacet, infacet_rr, infacet_ll;
	// Spatial discretization
	switch (ns->contr.spatialdiscrtype) {
	// Second order
	case 0:
		// Set the convective method
		ns->cc.convec_type = ns->contr.convecdiscrtype;
		// Get the velocity in the facet center
		ns->cc.ucell  = dp_get_value(dpu[dim], fgid);
		// Get the facet source tern in the facet center
		ns->cc.F      = compute_facet_value_at_point(ns->sfdF[dim], fcenter, fcenter, 1.0, ns->dpFU[dim], ns->stnF);
		// Get the pressure in the left cell
		pl            = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Get the pressure in the right cell
		pr            = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Compute the pressure derivative
		ns->cc.dpdx   = compute_dpdx_at_point(fdelta, dim, 0.5, pl, pr);
		// Get the cell source term in the left cell
		Fl            = compute_center_p_left(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Get the cell source term in the left cell
		Fr            = compute_center_p_right(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Compute the source term at the mid point
		ns->cc.F     += compute_value_at_mid_point(Fl, Fr);
		// Get the cell curvature in the left cell
		curvl = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpcurvature, ns->ed.stn);
		// Get the cell curvature in the right cell
		curvr = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpcurvature, ns->ed.stn);
		ns->cc.curv = compute_value_at_mid_point(curvl, curvr);
		// Get the cell density in the left cell
		rhol          = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpdens, ns->ed.stn);
		// Get the cell density in the right cell
		rhor          = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpdens, ns->ed.stn);
		// Compute the density at the mid point
		ns->cc.dens   = compute_value_at_mid_point(rhol, rhor);
		// Get the cell fraction in the left cell
		real fracl = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpfracvol, ns->ed.stn);
		// Get the cell fraction in the right cell
		real fracr = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpfracvol, ns->ed.stn);
		real curvaux = 0.0;
		//if (((0.0 < fracr)&&(fracr < 1.0)) || ((0.0 < fracl)&&(fracl < 1.0))){
		real wwi   = fracr*(1.0 - fracr);
		real wwim1 = fracl*(1.0 - fracl);
		if (wwi + wwim1 > 0.0){
			curvaux = (wwi*curvr + wwim1*curvl)/(wwi + wwim1);
			//printf("wwi + wwim1 = %.12lf  fracl= %.12lf fracr= %.12lf \n",wwi + wwim1,fracl,fracr);
			//getchar();
		}
		//ns->cc.IF= (fracr - fracl)*ns->cc.curv/(fdelta[dim]);
		ns->cc.IF= (fracr - fracl)*curvaux/(fdelta[dim]);

		//real IFl          = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpIF[dim], ns->ed.stn);
		//// Get the cell interfacial force in the right cell
		//real IFr          = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpIF[dim], ns->ed.stn);
		//// Compute the interfacial force at the mid point
		//ns->cc.IF   = compute_value_at_mid_point(IFl, IFr);

		for (int dim2 = 0; dim2 < DIM; dim2++) {
			// Get the velocity in the center of facet
			if (dim2 == dim) {
				ns->cc.v[dim2] = ns->cc.ucell;
			} else {
				ns->cc.v[dim2] = compute_facet_u_2(sfdu[dim2], fcenter, fdelta, dim, dim2, 1.0, dpu[dim2], ns->stn);
			}
			// Get the velocity in the left facet center
			ul[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Get the velocity in the right facet center
			ur[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Compute dudxc at facet center
			ns->cc.dudxc[dim2]  = compute_facet_dudxc(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxr at facet center
			ns->cc.dudxr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxl at facet center
			ns->cc.dudxl[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Get the velocity in the right facet center
			urr[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			// Get the velocity in the left facet center
			ull[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			if (ns->cc.convec_type == 2) {
				if (ns->contr.secondconvecdiscrtype == 2) {
					// Compute dudxrr at facet center
					ns->cc.dudxrr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ur[dim2], ns->cc.ucell, urr[dim2]);
					// Compute dudxll at facet center
					ns->cc.dudxll[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ul[dim2], ull[dim2], ns->cc.ucell);
					// Compute terms for second order upwind CUBISTA
				}else if (ns->contr.secondconvecdiscrtype == 1){
					// Get the velocity in the left facet center
					ul[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the left facet center
					ull[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_ll);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity in the right facet center
					ur[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the right facet center
					urr[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_rr);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity v1bar(i+1/2,j+1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_right(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ur[dim]-ul[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv1 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ul[dim]);
								if ((fi >= b) && (fi <= c))
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + b*ur[dim] -d*ul[dim]);
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(e*ns->cc.ucell + c*ur[dim]);
							}
						}
						//v1bar < 0.0
					}else {
						if (fabs(ns->cc.ucell - urr[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ur[dim];
						}else {
							fi = (ur[dim]-urr[dim])/(ns->cc.ucell -urr[dim]);
							//Trocada pela linha de cima --- fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ur[dim];
							}else {
								if (fi < b)
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(a*ur[dim] - c*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if ((fi >= b) && (fi <= c))
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(c*ur[dim] + b*ns->cc.ucell -d*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ur[dim]);
							}
						}
					}
					// Get the velocity  v2bar(i+1/2,j-1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_left(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ns->cc.ucell-ull[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ul[dim];
						}else {
							fi = (ul[dim] - ull[dim])/(ns->cc.ucell-ull[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ul[dim];
							}else {
								if (fi < b)
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(a*ul[dim] - c*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if ((fi >= b) && (fi <= c))
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(b*ns->cc.ucell + c*ul[dim] - d*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ul[dim]);
							}
						}
					}else {
						//v2bar < 0.0
						if (fabs(ul[dim] - ur[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ur[dim])/(ul[dim] - ur[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv2 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ur[dim]);
								if ((fi >= b) && (fi <= c))
									conv2 = ns->cc.vc[dim2]*(b*ul[dim] + c*ns->cc.ucell - d*ur[dim]);
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ul[dim] + e*ns->cc.ucell);
							}
						}
					}
					ns->cc.vc[dim2] = ((conv1-conv2)/fdelta[dim]);
				}
			}
			
			
			
			if(dim2==dim){
			
				// Get the cell viscosity in the left cell
				ns->cc.viscl = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.mult.dpvisc, ns->ed.stn);
				// Get the cell viscosity in the right cell
				ns->cc.viscr = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.mult.dpvisc, ns->ed.stn);
			}
			else{
				
				Point p1,p2,p3,p4,p3_,p4_;
				
				POINT_ASSIGN(p1, fcenter);POINT_ASSIGN(p2, fcenter);
				
				p1[dim]=p1[dim]-0.5*fdelta[dim];p2[dim]=p2[dim]+0.5*fdelta[dim];
				
				POINT_ASSIGN(p3, p1);POINT_ASSIGN(p4, p2);
				POINT_ASSIGN(p3_, p1);POINT_ASSIGN(p4_, p2);
				
				p3[dim2]=p3[dim2]+fdelta[dim2];p4[dim2]=p4[dim2]+fdelta[dim2];
				p3_[dim2]=p3_[dim2]-fdelta[dim2];p4_[dim2]=p4_[dim2]-fdelta[dim2];
				
				real v1=compute_value_at_point(ns->ed.sdED,fcenter,p1,1.0,ns->ed.mult.dpvisc,ns->ed.stn);
				real v2=compute_value_at_point(ns->ed.sdED,fcenter,p2,1.0,ns->ed.mult.dpvisc,ns->ed.stn);
				
				real v3=compute_value_at_point(ns->ed.sdED,fcenter,p3,1.0,ns->ed.mult.dpvisc,ns->ed.stn);
				real v4=compute_value_at_point(ns->ed.sdED,fcenter,p4,1.0,ns->ed.mult.dpvisc,ns->ed.stn);
				
				real v3_=compute_value_at_point(ns->ed.sdED,fcenter,p3_,1.0,ns->ed.mult.dpvisc,ns->ed.stn);
				real v4_=compute_value_at_point(ns->ed.sdED,fcenter,p4_,1.0,ns->ed.mult.dpvisc,ns->ed.stn);
				
				ns->cc.viscl=4*v1*v2*v3_*v4_/(v1+v2+v3_+v4_);
				ns->cc.viscr=4*v1*v2*v3*v4/(v1+v2+v3+v4);
				
				
				// Get the cell viscosity in the left cell
				//ns->cc.viscl = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.mult.dpvisc, ns->ed.stn);
				// Get the cell viscosity in the right cell
				//ns->cc.viscr = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.mult.dpvisc, ns->ed.stn);
				
			}
			
			
			// Compute multiphase viscous term
			ns->cc.du2dx2[dim2] = 0.0;
			Point p;
			POINT_ASSIGN(p, fcenter);

			p[dim2]      = fcenter[dim2] + 0.5*fdelta[dim2];

			real ucl     = compute_facet_u_left(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet);
			real ucr     = compute_facet_u_right(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet);
			real duidxjr = compute_facet_dudxc(fdelta, dim2, 0.5, ns->cc.ucell, ucl, ucr);

			real vcl     = compute_facet_u_left(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet);
			real vcr     = compute_facet_u_right(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet);
			real dujdxir = compute_facet_dudxc(fdelta, dim, 0.5, ns->cc.ucell, vcl, vcr);

			p[dim2]      =  fcenter[dim2] - 0.5*fdelta[dim2];

			ucl          = compute_facet_u_left(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet);
			ucr          = compute_facet_u_right(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet);
			real duidxjl = compute_facet_dudxc(fdelta, dim2, 0.5, ns->cc.ucell, ucl, ucr);

			vcl          = compute_facet_u_left(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet);
			vcr          = compute_facet_u_right(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet);
			real dujdxil = compute_facet_dudxc(fdelta, dim, 0.5, ns->cc.ucell, vcl, vcr);

			ns->cc.du2dx2[dim2] += (ns->cc.viscr*duidxjr - ns->cc.viscl*duidxjl)/fdelta[dim2];
			ns->cc.du2dx2[dim2] += (ns->cc.viscr*dujdxir - ns->cc.viscl*dujdxil)/fdelta[dim2];
			
			
			
			// Compute the viscoelastic contribution
			if (dim2 == dim) {
				// Get the tensor in the left cell
				Sl[dim2]          = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.mult.dpS[dim][dim2], ns->ed.stn);
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.mult.dpS[dim][dim2], ns->ed.stn);
				// Compute the tensor derivative
				ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 0.5, Sl[dim2], Sr[dim2]);
			} else {
				// Get the tensor in the left cell
				Sl[dim2]          = compute_center_p_left_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.mult.dpS[dim][dim2], ns->ed.stn);
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.mult.dpS[dim][dim2], ns->ed.stn);
				// Compute the tensor derivative
				ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 1.0, Sl[dim2], Sr[dim2]);
			}
			
			
			
		}
		//break;
	}
}

void higflow_computational_cell_imp_multiphase(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]) {
	// Set the computational cell
	real  uc[DIM], ul[DIM], ur[DIM], Sl[DIM], Sr[DIM], pl, pr, Fl, Fr;
	real  ull[DIM], urr[DIM];
	real  rhol, rhor;
	real  a, b, c, d, e, tol, fi, conv1, conv2, curvl, curvr;
	a     = 1.7500;
	b     = 0.3750;
	c     = 0.7500;
	d     = 0.1250;
	e     = 0.2500;
	tol   = 1.0e-14;
	conv1 = 0.0;
	conv2 = 0.0;
	int   infacet, infacet_rr, infacet_ll;
	// Spatial discretization
	switch (ns->contr.spatialdiscrtype) {
	// Second order
	case 0:
		// Set the convective method
		ns->cc.convec_type = ns->contr.convecdiscrtype;
		// Get the velocity in the facet center
		ns->cc.ucell  = dp_get_value(dpu[dim], fgid);
		// Get the facet source tern in the facet center
		ns->cc.F      = compute_facet_value_at_point(ns->sfdF[dim], fcenter, fcenter, 1.0, ns->dpFU[dim], ns->stnF);
		// Get the pressure in the left cell
		pl            = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Get the pressure in the right cell
		pr            = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Compute the pressure derivative
		ns->cc.dpdx   = compute_dpdx_at_point(fdelta, dim, 0.5, pl, pr);
		// Get the cell source term in the left cell
		Fl            = compute_center_p_left(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Get the cell source term in the left cell
		Fr            = compute_center_p_right(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Compute the source term at the mid point
		ns->cc.F     += compute_value_at_mid_point(Fl, Fr);
		// Get the cell curvature in the left cell
		curvl = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpcurvature, ns->ed.stn);
		// Get the cell curvature in the right cell
		curvr = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpcurvature, ns->ed.stn);
		//ns->cc.curv = compute_value_at_mid_point(curvl, curvr);
		// Get the cell density in the left cell
		rhol          = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpdens, ns->ed.stn);
		// Get the cell density in the right cell
		rhor          = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpdens, ns->ed.stn);
		// Compute the density at the mid point
		ns->cc.dens   = compute_value_at_mid_point(rhol, rhor);
		// Get the cell fraction in the left cell
		real fracl = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpfracvol, ns->ed.stn);
		// Get the cell fraction in the right cell
		real fracr = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpfracvol, ns->ed.stn);
		real curvaux = 0.0;
		//if (((0.0 < fracr)&&(fracr < 1.0)) || ((0.0 < fracl)&&(fracl < 1.0))){
		real wwi   = fracr*(1.0 - fracr);
		real wwim1 = fracl*(1.0 - fracl);
		if (wwi + wwim1 > 0.0){
			curvaux = (wwi*curvr + wwim1*curvl)/(wwi + wwim1);
			//printf("wwi + wwim1 = %.12lf  fracl= %.12lf fracr= %.12lf \n",wwi + wwim1,fracl,fracr);
			//getchar();
		}
		//ns->cc.IF= (fracr - fracl)*ns->cc.curv/(fdelta[dim]);
		ns->cc.IF= (fracr - fracl)*curvaux/(fdelta[dim]);
		//// Get the cell interfacial force in the left cell
		//real IFl          = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpIF[dim], ns->ed.stn);
		//// Get the cell interfacial force in the right cell
		//real IFr          = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpIF[dim], ns->ed.stn);
		//// Compute the interfacial force at the mid point
		//ns->cc.IF   = compute_value_at_mid_point(IFl, IFr);
		for (int dim2 = 0; dim2 < DIM; dim2++) {
			// Get the velocity in the center of facet
			if (dim2 == dim) {
				ns->cc.v[dim2] = ns->cc.ucell;
			} else {
				ns->cc.v[dim2] = compute_facet_u_2(sfdu[dim2], fcenter, fdelta, dim, dim2, 1.0, dpu[dim2], ns->stn);
			}
			// Get the velocity in the left facet center
			ul[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Get the velocity in the right facet center
			ur[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Compute dudxc at facet center
			ns->cc.dudxc[dim2]  = compute_facet_dudxc(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxr at facet center
			ns->cc.dudxr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxl at facet center
			ns->cc.dudxl[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Get the velocity in the right facet center
			urr[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			// Get the velocity in the left facet center
			ull[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			if (ns->cc.convec_type == 2) {
				if (ns->contr.secondconvecdiscrtype == 2) {
					// Compute dudxrr at facet center
					ns->cc.dudxrr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ur[dim2], ns->cc.ucell, urr[dim2]);
					// Compute dudxll at facet center
					ns->cc.dudxll[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ul[dim2], ull[dim2], ns->cc.ucell);
					// Compute terms for second order upwind CUBISTA
				}else if (ns->contr.secondconvecdiscrtype == 1){
					// Get the velocity in the left facet center
					ul[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the left facet center
					ull[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_ll);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity in the right facet center
					ur[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the right facet center
					urr[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_rr);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity v1bar(i+1/2,j+1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_right(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ur[dim]-ul[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv1 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ul[dim]);
								if ((fi >= b) && (fi <= c))
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + b*ur[dim] -d*ul[dim]);
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(e*ns->cc.ucell + c*ur[dim]);
							}
						}
						//v1bar < 0.0
					}else {
						if (fabs(ns->cc.ucell - urr[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ur[dim];
						}else {
							fi = (ur[dim]-urr[dim])/(ns->cc.ucell -urr[dim]);
							//Trocada pela linha de cima --- fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ur[dim];
							}else {
								if (fi < b)
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(a*ur[dim] - c*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if ((fi >= b) && (fi <= c))
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(c*ur[dim] + b*ns->cc.ucell -d*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ur[dim]);
							}
						}
					}
					// Get the velocity  v2bar(i+1/2,j-1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_left(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ns->cc.ucell-ull[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ul[dim];
						}else {
							fi = (ul[dim] - ull[dim])/(ns->cc.ucell-ull[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ul[dim];
							}else {
								if (fi < b)
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(a*ul[dim] - c*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if ((fi >= b) && (fi <= c))
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(b*ns->cc.ucell + c*ul[dim] - d*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ul[dim]);
							}
						}
					}else {
						//v2bar < 0.0
						if (fabs(ul[dim] - ur[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ur[dim])/(ul[dim] - ur[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv2 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ur[dim]);
								if ((fi >= b) && (fi <= c))
									conv2 = ns->cc.vc[dim2]*(b*ul[dim] + c*ns->cc.ucell - d*ur[dim]);
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ul[dim] + e*ns->cc.ucell);
							}
						}
					}
					ns->cc.vc[dim2] = ((conv1-conv2)/fdelta[dim]);
				}
			}

			// Get the cell viscosity in the left cell
			ns->cc.viscl = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.mult.dpvisc, ns->ed.stn);
			// Get the cell viscosity in the right cell
			ns->cc.viscr = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.mult.dpvisc, ns->ed.stn);

			// Compute multiphase viscous term
			ns->cc.du2dx2[dim2] = 0.0;
			Point p;
			POINT_ASSIGN(p, fcenter);

			p[dim2]      = fcenter[dim2] + 0.5*fdelta[dim2];
			/* non mixed */
			real ucl     = compute_facet_u_left(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet);
			real ucr     = compute_facet_u_right(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet);
			real duidxjr = compute_facet_dudxc(fdelta, dim2, 0.5, ns->cc.ucell, ucl, ucr);

			/* mixed */
			real vcl     = compute_facet_u_left(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet);
			real vcr     = compute_facet_u_right(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet);
			real dujdxir = compute_facet_dudxc(fdelta, dim, 0.5, ns->cc.ucell, vcl, vcr);

			p[dim2]      =  fcenter[dim2] - 0.5*fdelta[dim2];

			/* non mixed */
			ucl          = compute_facet_u_left(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet);
			ucr          = compute_facet_u_right(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet);
			real duidxjl = compute_facet_dudxc(fdelta, dim2, 0.5, ns->cc.ucell, ucl, ucr);

			/* mixed */
			vcl          = compute_facet_u_left(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet);
			vcr          = compute_facet_u_right(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet);
			real dujdxil = compute_facet_dudxc(fdelta, dim, 0.5, ns->cc.ucell, vcl, vcr);

			//ns->cc.du2dx2[dim2] += (ns->cc.viscr*duidxjr - ns->cc.viscl*duidxjl)/fdelta[dim];
			ns->cc.du2dx2[dim2] += (ns->cc.viscr*dujdxir - ns->cc.viscl*dujdxil)/fdelta[dim];
			
			
			// Compute the viscoelastic contribution
			if (dim2 == dim) {
				// Get the tensor in the left cell
				Sl[dim2]          = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.mult.dpS[dim][dim2], ns->ed.stn);
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.mult.dpS[dim][dim2], ns->ed.stn);
				// Compute the tensor derivative
				ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 0.5, Sl[dim2], Sr[dim2]);
			} else {
				// Get the tensor in the left cell
				Sl[dim2]          = compute_center_p_left_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.mult.dpS[dim][dim2], ns->ed.stn);
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.mult.dpS[dim][dim2], ns->ed.stn);
				// Compute the tensor derivative
				ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 1.0, Sl[dim2], Sr[dim2]);
			}
			
		}
		break;
	}
}

// Computing the necessary term for the viscoelastic Navier-Stokes equation
void higflow_computational_cell_viscoelastic(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]) {
	// Set the computational cell
	real  uc[DIM], ul[DIM], ur[DIM], Sl[DIM], Sr[DIM], pl, pr, Fl, Fr;
	real  ull[DIM], urr[DIM];
	real  a, b, c, d, e, tol, fi, conv1,conv2;
	a     = 1.7500;
	b     = 0.3750;
	c     = 0.7500;
	d     = 0.1250;
	e     = 0.2500;
	tol   = 1.0e-14;
	conv1 = 0.0;
	conv2 = 0.0;
	int   infacet, infacet_rr, infacet_ll;
	// Spatial discretization
	switch (ns->contr.spatialdiscrtype) {
	// Second order
	case 0:
		// Set the convective method
		ns->cc.convec_type = ns->contr.convecdiscrtype;
		// Get the velocity in the facet center
		ns->cc.ucell  = dp_get_value(dpu[dim], fgid);
		// Get the facet source tern in the facet center
		ns->cc.F      = compute_facet_value_at_point(ns->sfdF[dim], fcenter, fcenter, 1.0, ns->dpFU[dim], ns->stnF);
		// Get the pressure in the left cell
		pl            = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Get the pressure in the right cell
		pr            = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Compute the pressure derivative
		ns->cc.dpdx   = compute_dpdx_at_point(fdelta, dim, 0.5, pl, pr);
		// Get the cell source term in the left cell
		Fl            = compute_center_p_left(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Get the cell source term in the left cell
		Fr            = compute_center_p_right(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Compute the source term at the mid point
		ns->cc.F     += compute_value_at_mid_point(Fl, Fr);
		for (int dim2 = 0; dim2 < DIM; dim2++) {
			// Get the velocity in the center of facet
			if (dim2 == dim) {
				ns->cc.v[dim2] = ns->cc.ucell;
			} else {
				ns->cc.v[dim2] = compute_facet_u_2(sfdu[dim2], fcenter, fdelta, dim, dim2, 1.0, dpu[dim2], ns->stn);
			}
			// Get the velocity in the left facet center
			ul[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Get the velocity in the right facet center
			ur[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Compute dudxc at facet center
			ns->cc.dudxc[dim2]  = compute_facet_dudxc(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxr at facet center
			ns->cc.dudxr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxl at facet center
			ns->cc.dudxl[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Get the velocity in the right facet center
			urr[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			// Get the velocity in the left facet center
			ull[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			if (ns->cc.convec_type == 2) {
				if (ns->contr.secondconvecdiscrtype == 2) {
					// Compute dudxrr at facet center
					ns->cc.dudxrr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ur[dim2], ns->cc.ucell, urr[dim2]);
					// Compute dudxll at facet center
					ns->cc.dudxll[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ul[dim2], ull[dim2], ns->cc.ucell);
					// Compute terms for second order upwind CUBISTA
				}else if (ns->contr.secondconvecdiscrtype == 1){
					// Get the velocity in the left facet center
					ul[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the left facet center
					ull[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_ll);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity in the right facet center
					ur[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the right facet center
					urr[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_rr);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity v1bar(i+1/2,j+1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_right(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ur[dim]-ul[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv1 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ul[dim]);
								if ((fi >= b) && (fi <= c))
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + b*ur[dim] -d*ul[dim]);
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(e*ns->cc.ucell + c*ur[dim]);
							}
						}
						//v1bar < 0.0
					}else {
						if (fabs(ns->cc.ucell - urr[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ur[dim];
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ur[dim];
							}else {
								if (fi < b)
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(a*ur[dim] - c*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if ((fi >= b) && (fi <= c))
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(c*ur[dim] + b*ns->cc.ucell -d*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ur[dim]);
							}
						}
					}
					// Get the velocity  v2bar(i+1/2,j-1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_left(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ns->cc.ucell-ull[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ul[dim];
						}else {
							fi = (ul[dim] - ull[dim])/(ns->cc.ucell-ull[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ul[dim];
							}else {
								if (fi < b)
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(a*ul[dim] - c*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if ((fi >= b) && (fi <= c))
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(b*ns->cc.ucell + c*ul[dim] - d*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ul[dim]);
							}
						}
					}else {
						//v2bar < 0.0
						if (fabs(ul[dim] - ur[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ur[dim])/(ul[dim] - ur[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv2 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ur[dim]);
								if ((fi >= b) && (fi <= c))
									conv2 = ns->cc.vc[dim2]*(b*ul[dim] + c*ns->cc.ucell - d*ur[dim]);
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ul[dim] + e*ns->cc.ucell);
							}
						}
					}
					ns->cc.vc[dim2] = ((conv1-conv2)/fdelta[dim]);
				}
			}
			// Compute viscoelastic viscous term
			ns->cc.du2dx2[dim2] = compute_facet_du2dx2(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute the viscoelastic contribution
			if (dim2 == dim) {
				// Get the tensor in the left cell
				Sl[dim2]          = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.ve.dpS[dim][dim2], ns->ed.stn);
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.ve.dpS[dim][dim2], ns->ed.stn);
				// Compute the tensor derivative
				ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 0.5, Sl[dim2], Sr[dim2]);
			} else {
				// Get the tensor in the left cell
				Sl[dim2]          = compute_center_p_left_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.ve.dpS[dim][dim2], ns->ed.stn);
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.ve.dpS[dim][dim2], ns->ed.stn);
				// Compute the tensor derivative
				ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 1.0, Sl[dim2], Sr[dim2]);
			}
		}
		break;
	}
}

// Computing the necessary term for the viscoelastic Navier-Stokes equation
void higflow_computational_cell_viscoelastic_variable_viscosity(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]) {
	// Set the computational cell
	real  uc[DIM], ul[DIM], ur[DIM], Sl[DIM], Sr[DIM], pl, pr, Fl, Fr;
	real  ull[DIM], urr[DIM];
	real  a, b, c, d, e, tol, fi, conv1,conv2;
	a     = 1.7500;
	b     = 0.3750;
	c     = 0.7500;
	d     = 0.1250;
	e     = 0.2500;
	tol   = 1.0e-14;
	conv1 = 0.0;
	conv2 = 0.0;
	int   infacet, infacet_rr, infacet_ll;
	// Spatial discretization
	switch (ns->contr.spatialdiscrtype) {
	// Second order
	case 0:
		// Set the convective method
		ns->cc.convec_type = ns->contr.convecdiscrtype;
		// Get the velocity in the facet center
		ns->cc.ucell  = dp_get_value(dpu[dim], fgid);
		// Get the facet source tern in the facet center
		ns->cc.F      = compute_facet_value_at_point(ns->sfdF[dim], fcenter, fcenter, 1.0, ns->dpFU[dim], ns->stnF);
		// Get the pressure in the left cell
		pl            = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Get the pressure in the right cell
		pr            = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Compute the pressure derivative
		ns->cc.dpdx   = compute_dpdx_at_point(fdelta, dim, 0.5, pl, pr);
		// Get the cell source term in the left cell
		Fl            = compute_center_p_left(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Get the cell source term in the left cell
		Fr            = compute_center_p_right(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Compute the source term at the mid point
		ns->cc.F     += compute_value_at_mid_point(Fl, Fr);
		for (int dim2 = 0; dim2 < DIM; dim2++) {
			// Get the velocity in the center of facet
			if (dim2 == dim) {
				ns->cc.v[dim2] = ns->cc.ucell;
			} else {
				ns->cc.v[dim2] = compute_facet_u_2(sfdu[dim2], fcenter, fdelta, dim, dim2, 1.0, dpu[dim2], ns->stn);
			}
			// Get the velocity in the left facet center
			ul[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Get the velocity in the right facet center
			ur[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Compute dudxc at facet center
			ns->cc.dudxc[dim2]  = compute_facet_dudxc(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxr at facet center
			ns->cc.dudxr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxl at facet center
			ns->cc.dudxl[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Get the velocity in the right facet center
			urr[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			// Get the velocity in the left facet center
			ull[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			if (ns->cc.convec_type == 2) {
				if (ns->contr.secondconvecdiscrtype == 2) {
					// Compute dudxrr at facet center
					ns->cc.dudxrr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ur[dim2], ns->cc.ucell, urr[dim2]);
					// Compute dudxll at facet center
					ns->cc.dudxll[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ul[dim2], ull[dim2], ns->cc.ucell);
					// Compute terms for second order upwind CUBISTA
				}else if (ns->contr.secondconvecdiscrtype == 1){
					// Get the velocity in the left facet center
					ul[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the left facet center
					ull[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_ll);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity in the right facet center
					ur[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the right facet center
					urr[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_rr);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity v1bar(i+1/2,j+1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_right(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ur[dim]-ul[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv1 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ul[dim]);
								if ((fi >= b) && (fi <= c))
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + b*ur[dim] -d*ul[dim]);
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(e*ns->cc.ucell + c*ur[dim]);
							}
						}
						//v1bar < 0.0
					}else {
						if (fabs(ns->cc.ucell - urr[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ur[dim];
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ur[dim];
							}else {
								if (fi < b)
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(a*ur[dim] - c*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if ((fi >= b) && (fi <= c))
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(c*ur[dim] + b*ns->cc.ucell -d*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ur[dim]);
							}
						}
					}
					// Get the velocity  v2bar(i+1/2,j-1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_left(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ns->cc.ucell-ull[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ul[dim];
						}else {
							fi = (ul[dim] - ull[dim])/(ns->cc.ucell-ull[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ul[dim];
							}else {
								if (fi < b)
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(a*ul[dim] - c*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if ((fi >= b) && (fi <= c))
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(b*ns->cc.ucell + c*ul[dim] - d*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ul[dim]);
							}
						}
					}else {
						//v2bar < 0.0
						if (fabs(ul[dim] - ur[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ur[dim])/(ul[dim] - ur[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv2 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ur[dim]);
								if ((fi >= b) && (fi <= c))
									conv2 = ns->cc.vc[dim2]*(b*ul[dim] + c*ns->cc.ucell - d*ur[dim]);
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ul[dim] + e*ns->cc.ucell);
							}
						}
					}
					ns->cc.vc[dim2] = ((conv1-conv2)/fdelta[dim]);
				}
			}
			// Compute viscoelastic viscous term
			ns->cc.du2dx2[dim2] = compute_facet_du2dx2(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute the viscoelastic contribution
			if (dim2 == dim) {
				// Get the tensor in the left cell
				Sl[dim2]          = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.vevv.dpS[dim][dim2], ns->ed.stn);
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.vevv.dpS[dim][dim2], ns->ed.stn);
				// Compute the tensor derivative
				ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 0.5, Sl[dim2], Sr[dim2]);
			} else {
				// Get the tensor in the left cell
				Sl[dim2]          = compute_center_p_left_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.vevv.dpS[dim][dim2], ns->ed.stn);
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.vevv.dpS[dim][dim2], ns->ed.stn);
				// Compute the tensor derivative
				ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 1.0, Sl[dim2], Sr[dim2]);
			}
		}
		break;
	}
}

// Computing the necessary term for the viscoelastic Navier-Stokes equation
void higflow_computational_cell_viscoelastic_shear_banding(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]) {
	// Set the computational cell
	real  uc[DIM], ul[DIM], ur[DIM], Sl[DIM], Sr[DIM], pl, pr, Fl, Fr;
	real  ull[DIM], urr[DIM];
	real  a, b, c, d, e, tol, fi, conv1,conv2;
	a     = 1.7500;
	b     = 0.3750;
	c     = 0.7500;
	d     = 0.1250;
	e     = 0.2500;
	tol   = 1.0e-14;
	conv1 = 0.0;
	conv2 = 0.0;
	int   infacet, infacet_rr, infacet_ll;
	// Spatial discretization
	switch (ns->contr.spatialdiscrtype) {
	// Second order
	case 0:
		// Set the convective method
		ns->cc.convec_type = ns->contr.convecdiscrtype;
		// Get the velocity in the facet center
		ns->cc.ucell  = dp_get_value(dpu[dim], fgid);
		// Get the facet source tern in the facet center
		ns->cc.F      = compute_facet_value_at_point(ns->sfdF[dim], fcenter, fcenter, 1.0, ns->dpFU[dim], ns->stnF);
		// Get the pressure in the left cell
		pl            = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Get the pressure in the right cell
		pr            = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Compute the pressure derivative
		ns->cc.dpdx   = compute_dpdx_at_point(fdelta, dim, 0.5, pl, pr);
		// Get the cell source term in the left cell
		Fl            = compute_center_p_left(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Get the cell source term in the left cell
		Fr            = compute_center_p_right(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Compute the source term at the mid point
		ns->cc.F     += compute_value_at_mid_point(Fl, Fr);
		for (int dim2 = 0; dim2 < DIM; dim2++) {
			// Get the velocity in the center of facet
			if (dim2 == dim) {
				ns->cc.v[dim2] = ns->cc.ucell;
			} else {
				ns->cc.v[dim2] = compute_facet_u_2(sfdu[dim2], fcenter, fdelta, dim, dim2, 1.0, dpu[dim2], ns->stn);
			}
			// Get the velocity in the left facet center
			ul[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Get the velocity in the right facet center
			ur[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Compute dudxc at facet center
			ns->cc.dudxc[dim2]  = compute_facet_dudxc(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxr at facet center
			ns->cc.dudxr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxl at facet center
			ns->cc.dudxl[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Get the velocity in the right facet center
			urr[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			// Get the velocity in the left facet center
			ull[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			if (ns->cc.convec_type == 2) {
				if (ns->contr.secondconvecdiscrtype == 2) {
					// Compute dudxrr at facet center
					ns->cc.dudxrr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ur[dim2], ns->cc.ucell, urr[dim2]);
					// Compute dudxll at facet center
					ns->cc.dudxll[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ul[dim2], ull[dim2], ns->cc.ucell);
					// Compute terms for second order upwind CUBISTA
				}else if (ns->contr.secondconvecdiscrtype == 1){
					// Get the velocity in the left facet center
					ul[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the left facet center
					ull[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_ll);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity in the right facet center
					ur[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the right facet center
					urr[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_rr);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity v1bar(i+1/2,j+1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_right(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ur[dim]-ul[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv1 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ul[dim]);
								if ((fi >= b) && (fi <= c))
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + b*ur[dim] -d*ul[dim]);
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(e*ns->cc.ucell + c*ur[dim]);
							}
						}
						//v1bar < 0.0
					}else {
						if (fabs(ns->cc.ucell - urr[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ur[dim];
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ur[dim];
							}else {
								if (fi < b)
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(a*ur[dim] - c*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if ((fi >= b) && (fi <= c))
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(c*ur[dim] + b*ns->cc.ucell -d*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ur[dim]);
							}
						}
					}
					// Get the velocity  v2bar(i+1/2,j-1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_left(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ns->cc.ucell-ull[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ul[dim];
						}else {
							fi = (ul[dim] - ull[dim])/(ns->cc.ucell-ull[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ul[dim];
							}else {
								if (fi < b)
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(a*ul[dim] - c*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if ((fi >= b) && (fi <= c))
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(b*ns->cc.ucell + c*ul[dim] - d*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ul[dim]);
							}
						}
					}else {
						//v2bar < 0.0
						if (fabs(ul[dim] - ur[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ur[dim])/(ul[dim] - ur[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv2 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ur[dim]);
								if ((fi >= b) && (fi <= c))
									conv2 = ns->cc.vc[dim2]*(b*ul[dim] + c*ns->cc.ucell - d*ur[dim]);
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ul[dim] + e*ns->cc.ucell);
							}
						}
					}
					ns->cc.vc[dim2] = ((conv1-conv2)/fdelta[dim]);
				}
			}
			// Compute viscoelastic viscous term
			ns->cc.du2dx2[dim2] = compute_facet_du2dx2(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute the viscoelastic contribution
			if (dim2 == dim) {
				// Get the tensor in the left cell
				Sl[dim2]          = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.vesb.dpS[dim][dim2], ns->ed.stn);
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.vesb.dpS[dim][dim2], ns->ed.stn);
				// Compute the tensor derivative
				ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 0.5, Sl[dim2], Sr[dim2]);
			} else {
				// Get the tensor in the left cell
				Sl[dim2]          = compute_center_p_left_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.vesb.dpS[dim][dim2], ns->ed.stn);
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.vesb.dpS[dim][dim2], ns->ed.stn);
				// Compute the tensor derivative
				ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 1.0, Sl[dim2], Sr[dim2]);
			}
		}
		break;
	}
}

//Computing the terms needed for the simulation of viscoelastic flows that exhibit shear-banding behaviour using the VCM model
void higflow_computational_cell_shear_banding_VCM_model(higflow_solver *ns, sim_domain *sdn, int clid, Point ccenter, Point cdelta, int dim, distributed_property *dpn){
	real nc, nr, nl;
	//real  AC, BC;
	int  incell_l, incell_r;
	switch (ns->contr.spatialdiscrtype) {
	// Second order
	case 0:
		// Get the n in the cell center
		nc    = dp_get_value(dpn, clid);
		ns->cc.nABcell = nc;
		//AB  = dp_get_value(ns->ed.vesb.dpA[dim][dim], clid);
		//ns->cc.ABcell = AB;
		// Get the n in the left cell
		nl    = compute_center_p_left_22(sdn, ccenter, cdelta, dim, 1.0, dpn, ns->ed.stn, &incell_l);
		// Get the n in the right cell
		nr    = compute_center_p_right_22(sdn, ccenter, cdelta, dim, 1.0, dpn, ns->ed.stn, &incell_r);
		// Compute the necessary potentials derivative
		ns->cc.dnABdx     = compute_dpdx_at_point(cdelta, dim, 1.0, nl, nr);
		ns->cc.d2nABdx2   = compute_facet_du2dx2(cdelta, dim, 1.0, nc, nl, nr);
		break;
	}
}


//Computing the terms needed for the simulation of viscoelastic flows that exhibit shear-banding behaviour using the VCM model
//void higflow_computational_cell_conformation_tensor_shear_banding_VCM_model(higflow_solver *ns, sim_domain *sdn, int clid, Point ccenter, Point cdelta, int dim, int i, int j, distributed_property *dpA[DIM][DIM]){
	//real nc, nr, nl;
//	real  ABl[DIM][DIM], ABr[DIM][DIM];
//	int  incell_l, incell_r;
//	switch (ns->contr.spatialdiscrtype) {
	// Second order
//	case 0:
		// Get the tensor at the cell 
		//ns->cc.ABcell  = dp_get_value(dpA[i][j], clid);
		// Get the conformation tensor in the left cell
//		ABl = compute_center_p_left_22(sdn, ccenter, cdelta, dim, 1.0, dpA[i][j], ns->ed.stn, &incell_l);
		//compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.vesb.dpS[dim][dim2], ns->ed.stn)
		// Get the conformation tensor in the right cell
//		ABr = compute_center_p_right_22(sdn, ccenter, cdelta, dim, 1.0, dpA[i][j], ns->ed.stn, &incell_r);
		// Compute the necessary derivatives
//		ns->cc.d2ABdx2   = compute_facet_du2dx2(cdelta, dim, 1.0, dpA[i][j], ABl, ABr);
//		break;
//	}
//}

// Computing the necessary term for the Navier-Stokes equation
void higflow_computational_cell_electroosmotic(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]) {
	// Set the computational cell
	real  uc[DIM], ul[DIM], ur[DIM], Sl[DIM], Sr[DIM], pl, pr, Fl, Fr;
	real  ull[DIM], urr[DIM];
	real  a, b, c, d, e, tol, fi, conv1,conv2;
	a     = 1.7500;
	b     = 0.3750;
	c     = 0.7500;
	d     = 0.1250;
	e     = 0.2500;
	tol   = 1.0e-14;
	conv1 = 0.0;
	conv2 = 0.0;
	int   infacet, infacet_r, infacet_l, infacet_rr, infacet_ll;
	// Spatial discretization
	switch (ns->contr.spatialdiscrtype) {
	// Second order
	case 0:
		// Set the convective method
		ns->cc.convec_type = ns->contr.convecdiscrtype;
		// Get the velocity in the facet center
		ns->cc.ucell  = dp_get_value(dpu[dim], fgid);
		// Get the facet source tern in the facet center
		ns->cc.F      = compute_facet_value_at_point(ns->sfdF[dim], fcenter, fcenter, 1.0, ns->dpFU[dim], ns->stnF);
		// Get the pressure in the left cell
		pl            = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Get the pressure in the right cell
		pr            = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Compute the pressure derivative
		ns->cc.dpdx   = compute_dpdx_at_point(fdelta, dim, 0.5, pl, pr);
		// Get the cell source term in the left cell
		Fl            = compute_center_p_left(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Get the cell source term in the left cell
		Fr            = compute_center_p_right(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Compute the source term at the mid point
		ns->cc.F      += compute_value_at_mid_point(Fl, Fr);
		// Compute the EO source term at the mid point
		ns->cc.Feo   = compute_facet_value_at_point(ns->ed.eo.sfdEOFeo[dim], fcenter, fcenter, 1.0, ns->ed.eo.dpFeo[dim], ns->ed.stn);
		for (int dim2 = 0; dim2 < DIM; dim2++) {
			// Get the velocity in the center of facet
			if (dim2 == dim) {
				ns->cc.v[dim2] = ns->cc.ucell;
			} else {
				ns->cc.v[dim2] = compute_facet_u_2(sfdu[dim2], fcenter, fdelta, dim, dim2, 1.0, dpu[dim2], ns->stn);
			}
			// Get the velocity in the left facet center
			ul[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet_l);
			// Get the velocity in the right facet center
			ur[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet_r);
			// Compute dudxc at facet center
			ns->cc.dudxc[dim2]  = compute_facet_dudxc(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxr at facet center
			ns->cc.dudxr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxl at facet center
			ns->cc.dudxl[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet_l == 0)){// && (ns->cc.convec_type == 2)) {
				ul[dim2] = compute_facet_u_2_left(sfdu[dim], fcenter, fdelta, dim, dim2, 1.0, dpu[dim2], ns->stn);
				ns->cc.dudxl[dim2]  = compute_facet_dudxl(fdelta, dim2, 0.5, ns->cc.ucell, ul[dim2], ur[dim2]);
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				//ns->cc.convec_type = 1;
			}
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet_r == 0)){// && (ns->cc.convec_type == 2)) {
				ur[dim2] = compute_facet_u_2_right(sfdu[dim], fcenter, fdelta, dim, dim2, 1.0, dpu[dim2], ns->stn);
				ns->cc.dudxr[dim2]  = compute_facet_dudxl(fdelta, dim2, 0.5, ns->cc.ucell, ul[dim2], ur[dim2]);
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				//ns->cc.convec_type = 1;
			}
			// Get the velocity in the right facet center
			urr[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			// Get the velocity in the left facet center
			ull[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			if (ns->cc.convec_type == 2) {
				if (ns->contr.secondconvecdiscrtype == 2) {
					// Compute dudxrr at facet center
					ns->cc.dudxrr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ur[dim2], ns->cc.ucell, urr[dim2]);
					// Compute dudxll at facet center
					ns->cc.dudxll[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ul[dim2], ull[dim2], ns->cc.ucell);
					// Compute terms for second order upwind CUBISTA
				}else if (ns->contr.secondconvecdiscrtype == 1){
					// Get the velocity in the left facet center
					ul[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the left facet center
					ull[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_ll);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity in the right facet center
					ur[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the right facet center
					urr[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_rr);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity v1bar(i+1/2,j+1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_right(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ur[dim]-ul[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv1 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ul[dim]);
								if ((fi >= b) && (fi <= c))
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + b*ur[dim] -d*ul[dim]);
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(e*ns->cc.ucell + c*ur[dim]);
							}
						}
						//v1bar < 0.0
					}else {
						if (fabs(ns->cc.ucell - urr[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ur[dim];
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ur[dim];
							}else {
								if (fi < b)
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(a*ur[dim] - c*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if ((fi >= b) && (fi <= c))
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(c*ur[dim] + b*ns->cc.ucell -d*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ur[dim]);
							}
						}
					}
					// Get the velocity  v2bar(i+1/2,j-1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_left(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ns->cc.ucell-ull[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ul[dim];
						}else {
							fi = (ul[dim] - ull[dim])/(ns->cc.ucell-ull[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ul[dim];
							}else {
								if (fi < b)
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(a*ul[dim] - c*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if ((fi >= b) && (fi <= c))
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(b*ns->cc.ucell + c*ul[dim] - d*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ul[dim]);
							}
						}
					}else {
						//v2bar < 0.0
						if (fabs(ul[dim] - ur[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ur[dim])/(ul[dim] - ur[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv2 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ur[dim]);
								if ((fi >= b) && (fi <= c))
									conv2 = ns->cc.vc[dim2]*(b*ul[dim] + c*ns->cc.ucell - d*ur[dim]);
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ul[dim] + e*ns->cc.ucell);
							}
						}
					}
					ns->cc.vc[dim2] = ((conv1-conv2)/fdelta[dim]);
				}
			}
			// Newtonian
			ns->cc.du2dx2[dim2] = compute_facet_du2dx2(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute viscoelastic viscous term
			if (ns->contr.flowtype == 3){
				// Compute the viscoelastic contribution
				if (dim2 == dim) {
					// Get the tensor in the left cell
					Sl[dim2]          = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.ve.dpS[dim][dim2], ns->ed.stn);
					// Get the tensor in the right cell
					Sr[dim2]          = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.ve.dpS[dim][dim2], ns->ed.stn);
					// Compute the tensor derivative
					ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 0.5, Sl[dim2], Sr[dim2]);
				} else {
					// Get the tensor in the left cell
					Sl[dim2]          = compute_center_p_left_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.ve.dpS[dim][dim2], ns->ed.stn);
					// Get the tensor in the right cell
					Sr[dim2]          = compute_center_p_right_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.ve.dpS[dim][dim2], ns->ed.stn);
					// Compute the tensor derivative
					ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 1.0, Sl[dim2], Sr[dim2]);
				}
			}
		}
		break;
	}
}

/*void higflow_computational_cell_electroosmotic_ionic(higflow_solver *ns, sim_domain *sdn, sim_domain *sdpsi, sim_domain *sdphi, int clid, Point ccenter, Point cdelta, int dim, distributed_property *dpn, distributed_property *dppsi, distributed_property *dpphi){
              real psic, phic, psir,psil, phir, phil, nc, nr, nl, nplus, nminus; 
              int incell_r, incell_l;
              switch (ns->contr.spatialdiscrtype) {
            // Second order
            case 0:
            // Get the n in the cell center
            nplus           = compute_value_at_point(ns->ed.eo.sdEOnplus, ccenter, ccenter, 1.0, ns->ed.eo.dpnplus, ns->ed.stn);
            // Get the n in the cell center
            nminus          = compute_value_at_point(ns->ed.eo.sdEOnminus, ccenter, ccenter, 1.0, ns->ed.eo.dpnminus, ns->ed.stn);
            // Get the n in the cell center
            nc              = dp_get_value(dpn, clid);
            // Get the n in the left cell
            nl              = compute_center_p_left_22(sdn, ccenter, cdelta, dim, 1.0, dpn, ns->ed.stn, &incell_l); 
            // Get the n in the right cell
            nr              = compute_center_p_right_22(sdn, ccenter, cdelta, dim, 1.0, dpn, ns->ed.stn, &incell_r); 
            // Get the  psi in the cell center
            psic            = dp_get_value(dppsi, clid);
            // Get the  phi in the cell center
            phic            = dp_get_value(dpphi, clid);
            // Get the potential in the left cell
            psil            = compute_center_p_left(sdpsi, ccenter, cdelta, dim, 1.0, dppsi, ns->ed.stn);
            // Get the potential in the right cell
            psir            = compute_center_p_right(sdpsi, ccenter, cdelta, dim, 1.0, dppsi, ns->ed.stn);
            // Get the potential in the left cell
            phil            = compute_center_p_left(sdphi, ccenter, cdelta, dim, 1.0, dpphi, ns->ed.stn);
            // Get the potential in the right cell
            phir            = compute_center_p_right(sdphi, ccenter, cdelta, dim, 1.0, dpphi, ns->ed.stn);
            // Compute the necessary potentials derivative
            ns->cc.dpsidx   = compute_dpdx_at_point(cdelta, dim, 1.0, psil, psir);
            ns->cc.dphidx   = compute_dpdx_at_point(cdelta, dim, 1.0, phil, phir);
            //   ns->cc.d2psidx2 = compute_facet_du2dx2(cdelta, dim, 1.0, psic, psil, psir);
            //   ns->cc.d2phidx2 = compute_facet_du2dx2(cdelta, dim, 1.0, phic, phil, phir);
            ns->cc.d2psidx2 = (nplus - nminus)*ns->ed.eo.par.H*ns->ed.eo.par.H*ns->ed.eo.par.ez*ns->ed.eo.par.n0/ns->ed.eo.par.eps_eo/ns->ed.eo.par.zeta0;
            ns->cc.d2phidx2 = 0.0;
            ns->cc.d2ndx2   = compute_facet_du2dx2(cdelta, dim, 1.0, nc, nl, nr);
            break; 
            }
            }*/

void higflow_computational_cell_electroosmotic_ionic(higflow_solver *ns, sim_domain *sdn, sim_domain *sdpsi, sim_domain *sdphi, int clid, Point ccenter, Point cdelta, int dim, distributed_property *dpn, distributed_property *dppsi, distributed_property *dpphi){
	real psic, phic, psir,psil, phir, phil, nc, nr, nl;
	int  incell_l, incell_r;
	switch (ns->contr.spatialdiscrtype) {
	// Second order
	case 0:
		// Get the n in the cell center
		nc    = dp_get_value(dpn, clid);
		ns->cc.ncell = nc;
		// Get the  psi in the cell center
		psic  = dp_get_value(dppsi, clid);
		// Get the  phi in the cell center
		phic  = dp_get_value(dpphi, clid);
		// Get the n in the left cell
		nl    = compute_center_p_left_22(sdn, ccenter, cdelta, dim, 1.0, dpn, ns->ed.stn, &incell_l);
		// Get the n in the right cell
		nr    = compute_center_p_right_22(sdn, ccenter, cdelta, dim, 1.0, dpn, ns->ed.stn, &incell_r);
		// Get the potential in the left cell
		psil  = compute_center_p_left(sdpsi, ccenter, cdelta, dim, 1.0, dppsi, ns->ed.stn);
		// Get the potential in the right cell
		psir  = compute_center_p_right(sdpsi, ccenter, cdelta, dim, 1.0, dppsi, ns->ed.stn);
		// Get the potential in the left cell
		phil  = compute_center_p_left(sdphi, ccenter, cdelta, dim, 1.0, dpphi, ns->ed.stn);
		// Get the potential in the right cell
		phir  = compute_center_p_right(sdphi, ccenter, cdelta, dim, 1.0, dpphi, ns->ed.stn);
		// Compute the necessary potentials derivative
		ns->cc.dpsidx   = compute_dpdx_at_point(cdelta, dim, 1.0, psil, psir);
		ns->cc.dphidx   = compute_dpdx_at_point(cdelta, dim, 1.0, phil, phir);
		ns->cc.dndx     = compute_dpdx_at_point(cdelta, dim, 1.0, nl, nr);
		ns->cc.d2ndx2   = compute_facet_du2dx2(cdelta, dim, 1.0, nc, nl, nr);
		ns->cc.d2psidx2 = compute_facet_du2dx2(cdelta, dim, 1.0, psic, psil, psir);
		//            ns->cc.d2phidx2 = compute_facet_du2dx2(cdelta, dim, 1.0, phic, phil, phir);
		//            ns->cc.d2psidx2 = ns->ed.eo.par.delta*(dp_get_value(ns->ed.eo.dpnminus, clid) - dp_get_value(ns->ed.eo.dpnplus,clid));
		ns->cc.d2phidx2 = 0.0;
		break;
	}
}


// Computing the necessary term for the viscoelastic Navier-Stokes equation
void higflow_computational_cell_viscoelastic_integral(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]) {
	// Set the computational cell
	real  uc[DIM], ul[DIM], ur[DIM], Sl[DIM], Sr[DIM], pl, pr, Fl, Fr;
	real  ull[DIM], urr[DIM];
	real  a, b, c, d, e, tol, fi, conv1,conv2;
	a     = 1.7500;
	b     = 0.3750;
	c     = 0.7500;
	d     = 0.1250;
	e     = 0.2500;
	tol   = 1.0e-14;
	conv1 = 0.0;
	conv2 = 0.0;
	int   infacet, infacet_rr, infacet_ll;
	// Spatial discretization
	switch (ns->contr.spatialdiscrtype) {
	// Second order
	case 0:
		// Set the convective method
		ns->cc.convec_type = ns->contr.convecdiscrtype;
		// Get the velocity in the facet center
		ns->cc.ucell  = dp_get_value(dpu[dim], fgid);
		// Get the facet source tern in the facet center
		ns->cc.F      = compute_facet_value_at_point(ns->sfdF[dim], fcenter, fcenter, 1.0, ns->dpFU[dim], ns->stnF);
		// Get the pressure in the left cell
		pl            = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Get the pressure in the right cell
		pr            = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Compute the pressure derivative
		ns->cc.dpdx   = compute_dpdx_at_point(fdelta, dim, 0.5, pl, pr);
		// Get the cell source term in the left cell
		//Fl            = compute_center_p_left(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Get the cell source term in the left cell
		//Fr            = compute_center_p_right(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Compute the source term at the mid point
		//ns->cc.F     += compute_value_at_mid_point(Fl, Fr);
		for (int dim2 = 0; dim2 < DIM; dim2++) {
			// Get the velocity in the center of facet
			if (dim2 == dim) {
				ns->cc.v[dim2] = ns->cc.ucell;
			} else {
				ns->cc.v[dim2] = compute_facet_u_2(sfdu[dim2], fcenter, fdelta, dim, dim2, 1.0, dpu[dim2], ns->stn);
			}
			// Get the velocity in the left facet center
			ul[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Get the velocity in the right facet center
			ur[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Compute dudxc at facet center
			ns->cc.dudxc[dim2]  = compute_facet_dudxc(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxr at facet center
			ns->cc.dudxr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxl at facet center
			ns->cc.dudxl[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Get the velocity in the right facet center
			urr[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			// Get the velocity in the left facet center
			ull[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			if (ns->cc.convec_type == 2) {
				if (ns->contr.secondconvecdiscrtype == 2) {
					// Compute dudxrr at facet center
					ns->cc.dudxrr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ur[dim2], ns->cc.ucell, urr[dim2]);
					// Compute dudxll at facet center
					ns->cc.dudxll[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ul[dim2], ull[dim2], ns->cc.ucell);
					// Compute terms for second order upwind CUBISTA
				}else if (ns->contr.secondconvecdiscrtype == 1){
					// Get the velocity in the left facet center
					ul[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the left facet center
					ull[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_ll);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity in the right facet center
					ur[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the right facet center
					urr[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_rr);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity v1bar(i+1/2,j+1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_right(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ur[dim]-ul[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv1 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ul[dim]);
								if ((fi >= b) && (fi <= c))
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + b*ur[dim] -d*ul[dim]);
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(e*ns->cc.ucell + c*ur[dim]);
							}
						}
						//v1bar < 0.0
					}else {
						if (fabs(ns->cc.ucell - urr[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ur[dim];
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ur[dim];
							}else {
								if (fi < b)
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(a*ur[dim] - c*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if ((fi >= b) && (fi <= c))
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(c*ur[dim] + b*ns->cc.ucell -d*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ur[dim]);
							}
						}
					}
					// Get the velocity  v2bar(i+1/2,j-1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_left(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ns->cc.ucell-ull[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ul[dim];
						}else {
							fi = (ul[dim] - ull[dim])/(ns->cc.ucell-ull[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ul[dim];
							}else {
								if (fi < b)
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(a*ul[dim] - c*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if ((fi >= b) && (fi <= c))
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(b*ns->cc.ucell + c*ul[dim] - d*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ul[dim]);
							}
						}
					}else {
						//v2bar < 0.0
						if (fabs(ul[dim] - ur[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ur[dim])/(ul[dim] - ur[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv2 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ur[dim]);
								if ((fi >= b) && (fi <= c))
									conv2 = ns->cc.vc[dim2]*(b*ul[dim] + c*ns->cc.ucell - d*ur[dim]);
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ul[dim] + e*ns->cc.ucell);
							}
						}
					}
					ns->cc.vc[dim2] = ((conv1-conv2)/fdelta[dim]);
				}
			}
			// Compute viscoelastic viscous term
			ns->cc.du2dx2[dim2] = compute_facet_du2dx2(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute the viscoelastic contribution
			if (dim2 == dim) {
				// Get the tensor in the left cell
				Sl[dim2]          = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.im.dpS[dim][dim2], ns->ed.stn);
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.im.dpS[dim][dim2], ns->ed.stn);
				// Compute the tensor derivative
				ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 0.5, Sl[dim2], Sr[dim2]);
			} else {
				// Get the tensor in the left cell
				Sl[dim2]          = compute_center_p_left_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.im.dpS[dim][dim2], ns->ed.stn);
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.im.dpS[dim][dim2], ns->ed.stn);
				// Compute the tensor derivative
				ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 1.0, Sl[dim2], Sr[dim2]);
			}
		}
		break;
	}
}

// Computing the necessary term for the elastoviscoplastic Navier-Stokes equation
void higflow_computational_cell_elastoviscoplastic(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]) {
	// Set the computational cell
	real  uc[DIM], ul[DIM], ur[DIM], Sl[DIM], Sr[DIM], pl, pr, Fl, Fr;
	real  ull[DIM], urr[DIM];
	real  a, b, c, d, e, tol, fi, conv1,conv2;
	a     = 1.7500;
	b     = 0.3750;
	c     = 0.7500;
	d     = 0.1250;
	e     = 0.2500;
	tol   = 1.0e-14;
	conv1 = 0.0;
	conv2 = 0.0;
	int   infacet, infacet_rr, infacet_ll;
	// Spatial discretization
	switch (ns->contr.spatialdiscrtype) {
	// Second order
	case 0:
		// Set the convective method
		ns->cc.convec_type = ns->contr.convecdiscrtype;
		// Get the velocity in the facet center
		ns->cc.ucell  = dp_get_value(dpu[dim], fgid);
		// Get the facet source tern in the facet center
		ns->cc.F      = compute_facet_value_at_point(ns->sfdF[dim], fcenter, fcenter, 1.0, ns->dpFU[dim], ns->stnF);
		// Get the pressure in the left cell
		pl            = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Get the pressure in the right cell
		pr            = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Compute the pressure derivative
		ns->cc.dpdx   = compute_dpdx_at_point(fdelta, dim, 0.5, pl, pr);
		// Get the cell source term in the left cell
		Fl            = compute_center_p_left(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Get the cell source term in the left cell
		Fr            = compute_center_p_right(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Compute the source term at the mid point
		ns->cc.F     += compute_value_at_mid_point(Fl, Fr);
		for (int dim2 = 0; dim2 < DIM; dim2++) {
			// Get the velocity in the center of facet
			if (dim2 == dim) {
				ns->cc.v[dim2] = ns->cc.ucell;
			} else {
				ns->cc.v[dim2] = compute_facet_u_2(sfdu[dim2], fcenter, fdelta, dim, dim2, 1.0, dpu[dim2], ns->stn);
			}
			// Get the velocity in the left facet center
			ul[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Get the velocity in the right facet center
			ur[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Compute dudxc at facet center
			ns->cc.dudxc[dim2]  = compute_facet_dudxc(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxr at facet center
			ns->cc.dudxr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxl at facet center
			ns->cc.dudxl[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Get the velocity in the right facet center
			urr[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			// Get the velocity in the left facet center
			ull[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			if (ns->cc.convec_type == 2) {
				if (ns->contr.secondconvecdiscrtype == 2) {
					// Compute dudxrr at facet center
					ns->cc.dudxrr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ur[dim2], ns->cc.ucell, urr[dim2]);
					// Compute dudxll at facet center
					ns->cc.dudxll[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ul[dim2], ull[dim2], ns->cc.ucell);
					// Compute terms for second order upwind CUBISTA
				}else if (ns->contr.secondconvecdiscrtype == 1){
					// Get the velocity in the left facet center
					ul[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the left facet center
					ull[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_ll);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity in the right facet center
					ur[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the right facet center
					urr[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_rr);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity v1bar(i+1/2,j+1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_right(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ur[dim]-ul[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv1 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ul[dim]);
								if ((fi >= b) && (fi <= c))
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + b*ur[dim] -d*ul[dim]);
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(e*ns->cc.ucell + c*ur[dim]);
							}
						}
						//v1bar < 0.0
					}else {
						if (fabs(ns->cc.ucell - urr[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ur[dim];
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ur[dim];
							}else {
								if (fi < b)
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(a*ur[dim] - c*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if ((fi >= b) && (fi <= c))
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(c*ur[dim] + b*ns->cc.ucell -d*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ur[dim]);
							}
						}
					}
					// Get the velocity  v2bar(i+1/2,j-1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_left(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ns->cc.ucell-ull[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ul[dim];
						}else {
							fi = (ul[dim] - ull[dim])/(ns->cc.ucell-ull[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ul[dim];
							}else {
								if (fi < b)
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(a*ul[dim] - c*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if ((fi >= b) && (fi <= c))
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(b*ns->cc.ucell + c*ul[dim] - d*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ul[dim]);
							}
						}
					}else {
						//v2bar < 0.0
						if (fabs(ul[dim] - ur[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ur[dim])/(ul[dim] - ur[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv2 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ur[dim]);
								if ((fi >= b) && (fi <= c))
									conv2 = ns->cc.vc[dim2]*(b*ul[dim] + c*ns->cc.ucell - d*ur[dim]);
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ul[dim] + e*ns->cc.ucell);
							}
						}
					}
					ns->cc.vc[dim2] = ((conv1-conv2)/fdelta[dim]);
				}
			}
			// Compute viscoelastic viscous term
			ns->cc.du2dx2[dim2] = compute_facet_du2dx2(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute the viscoelastic contribution
			if (dim2 == dim) {
				// Get the tensor in the left cell
				Sl[dim2]          = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.vepl.dpS[dim][dim2], ns->ed.stn);
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.vepl.dpS[dim][dim2], ns->ed.stn);
				// Compute the tensor derivative
				ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 0.5, Sl[dim2], Sr[dim2]);
			} else {
				// Get the tensor in the left cell
				Sl[dim2]          = compute_center_p_left_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.vepl.dpS[dim][dim2], ns->ed.stn);
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.vepl.dpS[dim][dim2], ns->ed.stn);
				// Compute the tensor derivative
				ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 1.0, Sl[dim2], Sr[dim2]);
			}
		}
		break;
	}
}

// Computing the necessary term for the shear-thickening suspensions Navier-Stokes equation
void higflow_computational_cell_shear_thickening_suspensions(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]) {
	// Set the computational cell
	real  uc[DIM], ul[DIM], ur[DIM], Sl[DIM], Sr[DIM], pl, pr, Fl, Fr;
	real  ull[DIM], urr[DIM];
	real  a, b, c, d, e, tol, fi, conv1,conv2;
	a     = 1.7500;
	b     = 0.3750;
	c     = 0.7500;
	d     = 0.1250;
	e     = 0.2500;
	tol   = 1.0e-14;
	conv1 = 0.0;
	conv2 = 0.0;
	int   infacet, infacet_rr, infacet_ll;
	// Spatial discretization
	switch (ns->contr.spatialdiscrtype) {
	// Second order
	case 0:
		// Set the convective method
		ns->cc.convec_type = ns->contr.convecdiscrtype;
		// Get the velocity in the facet center
		ns->cc.ucell  = dp_get_value(dpu[dim], fgid);
		// Get the facet source tern in the facet center
		ns->cc.F      = compute_facet_value_at_point(ns->sfdF[dim], fcenter, fcenter, 1.0, ns->dpFU[dim], ns->stnF);
		// Get the pressure in the left cell
		pl            = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Get the pressure in the right cell
		pr            = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Compute the pressure derivative
		ns->cc.dpdx   = compute_dpdx_at_point(fdelta, dim, 0.5, pl, pr);
		// Get the cell source term in the left cell
		Fl            = compute_center_p_left(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Get the cell source term in the left cell
		Fr            = compute_center_p_right(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Compute the source term at the mid point
		ns->cc.F     += compute_value_at_mid_point(Fl, Fr);
		for (int dim2 = 0; dim2 < DIM; dim2++) {
			// Get the velocity in the center of facet
			if (dim2 == dim) {
				ns->cc.v[dim2] = ns->cc.ucell;
			} else {
				ns->cc.v[dim2] = compute_facet_u_2(sfdu[dim2], fcenter, fdelta, dim, dim2, 1.0, dpu[dim2], ns->stn);
			}
			// Get the velocity in the left facet center
			ul[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Get the velocity in the right facet center
			ur[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			// Compute dudxc at facet center
			ns->cc.dudxc[dim2]  = compute_facet_dudxc(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxr at facet center
			ns->cc.dudxr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxl at facet center
			ns->cc.dudxl[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Get the velocity in the right facet center
			urr[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			// Get the velocity in the left facet center
			ull[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			if (ns->cc.convec_type == 2) {
				if (ns->contr.secondconvecdiscrtype == 2) {
					// Compute dudxrr at facet center
					ns->cc.dudxrr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ur[dim2], ns->cc.ucell, urr[dim2]);
					// Compute dudxll at facet center
					ns->cc.dudxll[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ul[dim2], ull[dim2], ns->cc.ucell);
					// Compute terms for second order upwind CUBISTA
				}else if (ns->contr.secondconvecdiscrtype == 1){
					// Get the velocity in the left facet center
					ul[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the left facet center
					ull[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_ll);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity in the right facet center
					ur[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the right facet center
					urr[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_rr);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity v1bar(i+1/2,j+1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_right(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ur[dim]-ul[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv1 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ul[dim]);
								if ((fi >= b) && (fi <= c))
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + b*ur[dim] -d*ul[dim]);
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(e*ns->cc.ucell + c*ur[dim]);
							}
						}
						//v1bar < 0.0
					}else {
						if (fabs(ns->cc.ucell - urr[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ur[dim];
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ur[dim];
							}else {
								if (fi < b)
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(a*ur[dim] - c*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if ((fi >= b) && (fi <= c))
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(c*ur[dim] + b*ns->cc.ucell -d*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ur[dim]);
							}
						}
					}
					// Get the velocity  v2bar(i+1/2,j-1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_left(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ns->cc.ucell-ull[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ul[dim];
						}else {
							fi = (ul[dim] - ull[dim])/(ns->cc.ucell-ull[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ul[dim];
							}else {
								if (fi < b)
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(a*ul[dim] - c*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if ((fi >= b) && (fi <= c))
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(b*ns->cc.ucell + c*ul[dim] - d*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ul[dim]);
							}
						}
					}else {
						//v2bar < 0.0
						if (fabs(ul[dim] - ur[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ur[dim])/(ul[dim] - ur[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv2 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ur[dim]);
								if ((fi >= b) && (fi <= c))
									conv2 = ns->cc.vc[dim2]*(b*ul[dim] + c*ns->cc.ucell - d*ur[dim]);
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ul[dim] + e*ns->cc.ucell);
							}
						}
					}
					ns->cc.vc[dim2] = ((conv1-conv2)/fdelta[dim]);
				}
			}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			// Compute viscous term
			ns->cc.du2dx2[dim2] = compute_facet_du2dx2(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute the viscoelastic contribution
			if (dim2 == dim) {
				// Get the tensor in the left cell
				Sl[dim2]          = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				// Compute the tensor derivative
				ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 0.5, Sl[dim2], Sr[dim2]);
			} else {
				// Get the tensor in the left cell
				Sl[dim2]          = compute_center_p_left_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				// Compute the tensor derivative
				ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 1.0, Sl[dim2], Sr[dim2]);
			}
		}
		break;
	}
}

//Computational cell used to calculate the terms of the particle migration equation (evolution of the fraction volume)
void higflow_computational_cell_volume_fraction_suspensions(higflow_solver *ns, sim_domain *sdphi, int clid, Point ccenter, Point cdelta, int dim,  distributed_property *dpphi){
	real psic, phic, psir,psil, phir, phil, nc, nr, nl;
	real varphic, varphil, varphir;
	real Sl[DIM], Sr[DIM];
	real Sll[DIM], Srr[DIM];
	real Sc[DIM];
	real DU1l[DIM], DU1r[DIM];
	real DU12l[DIM], DU12r[DIM];
	real DU1ll[DIM], DU1rr[DIM];
	int  incell_l, incell_r;
	real frac_tol   = ns->ed.stsp.par.gdrms;
    real eps = 1.0e-4;
	real phij = ns->ed.stsp.par.phi;
	switch (ns->contr.spatialdiscrtype) {
	// Second order
	case 0:
		// Get the volume fraction value in the cell center
		varphic    = dp_get_value(dpphi, clid);
		if(varphic > 1.0){
			varphic = 1.0-eps;
        }
		//if(varphic > phij- frac_tol){
		//	varphic = phij- frac_tol;
        //}
        if(varphic < eps){
			varphic = 0.0;
        }
		ns->cc.phivfcell = varphic;
		// Get the  psi in the cell center
		//psic  = dp_get_value(dppsi, clid);
		// Get the  phi in the cell center
		//phic  = dp_get_value(dpphi, clid);
		// Get the n in the left cell
		//varphil    = compute_center_p_left_22(sdphi, ccenter, cdelta, dim, 1.0, dpphi, ns->ed.stn, &incell_l);
		// Get the n in the right cell
		//varphir    = compute_center_p_right_22(sdphi, ccenter, cdelta, dim, 1.0, dpphi, ns->ed.stn, &incell_r);
		// Get the n in the left cell

		varphil  = compute_center_p_left_22(sdphi, ccenter, cdelta, dim, 0.5, dpphi, ns->ed.stn, &incell_l);
		if(varphil > 1.0){
			varphil = 1.0-eps;
        }
		//if(varphil > phij- frac_tol){
		//	varphil = phij- frac_tol;
        //}
        if(varphil < eps){
			varphil = 0.0;
        }
		// Get the n in the right cell
		varphir  = compute_center_p_right_22(sdphi, ccenter, cdelta, dim, 0.5, dpphi, ns->ed.stn, &incell_r);
		if(varphir > 1.0){
			varphir = 1.0-eps;
        }
		//if(varphir > phij- frac_tol){
		//	varphir = phij- frac_tol;
        //}
        if(varphir < eps){
			varphir = 0.0;
        }
		// Get the potential in the left cell
		//psil  = compute_center_p_left(sdpsi, ccenter, cdelta, dim, 1.0, dppsi, ns->ed.stn);
		// Get the potential in the right cell
		//psir  = compute_center_p_right(sdpsi, ccenter, cdelta, dim, 1.0, dppsi, ns->ed.stn);
		// Get the potential in the left cell
		//phil  = compute_center_p_left(sdphi, ccenter, cdelta, dim, 1.0, dpphi, ns->ed.stn);
		// Get the potential in the right cell
		//phir  = compute_center_p_right(sdphi, ccenter, cdelta, dim, 1.0, dpphi, ns->ed.stn);
		// Compute the necessary potentials derivative
		//ns->cc.dpsidx   = compute_dpdx_at_point(cdelta, dim, 1.0, psil, psir);
		//ns->cc.dphivfdx   = compute_dpdx_at_point(cdelta, dim, 1.0, varphil, varphir);
		ns->cc.dphivfdx   = compute_dpdx_at_point(cdelta, dim, 0.5, varphil, varphir);
		/*
		if (dim == 1) {
			ns->cc.dphivfdx   = compute_dpdx_at_point(cdelta, dim, 0.5, varphil, varphir);
		} else {
			ns->cc.dphivfdx   = 0.0;
		}
		*/
		//ns->cc.dphivfdx   = compute_dpdx_at_point(cdelta, dim, 0.5, varphil, varphir);
		// Get the velocity derivative tensor Du and the microstructure tensor
		/*
		real Du[DIM][DIM];
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				// Get Du
				Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
			}
		}
        // Rate of strain tensor E
		real E[DIM][DIM];
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				E[i][j] = 0.5 * (Du[i][j] + Du[j][i]);
			}
		}
		*/
		//DEBUG_INSPECT(E[0][1], %lf);
		//real d2Sdx2 = 0.0;
		real d2Tdx2t = 0.0;
		real dTdxpR  = 0.0;
		real dTdxpL = 0.0;
		real d2Tdx2p = 0.0;
		for (int dim2 = 0; dim2 < DIM; dim2++) {
			// Compute the stress contribution
			if (dim2 == dim) {
				// Get the tensor in the left cell
				Sl[dim2]          = compute_center_p_left(ns->ed.sdED, ccenter, cdelta, dim2, 0.5, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				//DU1l[dim2] 		  = compute_center_p_left(ns->ed.sdED, ccenter, cdelta, dim2, 0.5, ns->ed.stsp.dpD[dim][dim2], ns->ed.stn);
				//DU12l[dim2]       =  compute_center_p_left(ns->ed.sdED, ccenter, cdelta, dim2, 0.5, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//DU1ll[dim2] 	  = compute_center_p_left(ns->ed.sdED, ccenter, cdelta, dim2, 0.5, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//Sl[dim2]         += 1.0*(ns->ed.stsp.par.eta0)*(DU1l[dim2]);
				//Sl[dim2]         += 2.0*(ns->ed.stsp.par.eta0)*(DU1ll[dim2]+DU12l[dim2]);
				//DEBUG_INSPECT(Sl[0], %lf);
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right(ns->ed.sdED, ccenter, cdelta, dim2, 0.5, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				//DU1r[dim2] 		  = compute_center_p_right(ns->ed.sdED, ccenter, cdelta, dim2, 0.5, ns->ed.stsp.dpD[dim][dim2], ns->ed.stn);
				//DU12r[dim2]       = compute_center_p_right(ns->ed.sdED, ccenter, cdelta, dim2, 0.5, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//DU1rr[dim2] 	  = compute_center_p_right(ns->ed.sdED, ccenter, cdelta, dim2, 0.5, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//Sr[dim2]         += 1.0*(ns->ed.stsp.par.eta0)*DU1r[dim2];
				//Sr[dim2]         += 2.0*(ns->ed.stsp.par.eta0)*(DU1rr[dim2]+DU12r[dim2]);
				//DEBUG_INSPECT(Sr[0], %lf);
				Sc[dim2] 		  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				//real DUc          = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpD[dim][dim2], ns->ed.stn);
				//real DUc1         = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//DEBUG_INSPECT(Sc, %lf);
				//Sc 		         += 1.0*(ns->ed.stsp.par.eta0)*DUc;
				//Sc 		         += 2.0*(ns->ed.stsp.par.eta0)*(DUc1+DUc);
				// Compute the tensor derivative
				ns->cc.dTdx[dim2] = compute_dpdx_at_point(cdelta, dim2, 0.5, Sl[dim2], Sr[dim2]);
				
				//d2Tdx2t      = compute_facet_du2dx2(cdelta, dim, 1.0, Sc[dim2], Sl[dim2], Sr[dim2]);
				d2Tdx2t = compute_facet_du2dx2(cdelta, dim, 1.0, Sc[dim2], Sl[dim2], Sr[dim2]);
				//DEBUG_INSPECT(ns->cc.dTdx[0], %lf);
				//if (dim2 == 1) {
				//	ns->cc.d2Tdx2 += d2Tdx2t;
				//} else {
				//	ns->cc.d2Tdx2 +=  0.0;
				//}
				//ns->cc.d2Tdx2	  += dim2;
				ns->cc.d2Tdx2	  += d2Tdx2t;
			} else {
				// Get the tensor in the left cell
				//Sl[dim2]          = compute_center_p_left_2(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				//DU1l[dim2]        = compute_center_p_left_2(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpD[dim][dim2], ns->ed.stn);
				Sl[dim2]          = compute_center_p_left_2_r(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				//DU1l[dim2]        = compute_center_p_left_2_r(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpD[dim][dim2], ns->ed.stn);
				//DU12l[dim2]       = compute_center_p_left_2_r(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//DU1ll[dim2]       = compute_center_p_left_2(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//Sl[dim2]         += 2.0*(ns->ed.stsp.par.eta0)*(DU1l[dim2]+DU12l[dim2]);
				//Sl[dim2]         += 2.0*(ns->ed.stsp.par.eta0)*DU1ll[dim2];
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right_2_r(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				//DU1r[dim2]        = compute_center_p_right_2_r(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpD[dim][dim2], ns->ed.stn);
				//DU12r[dim2]        = compute_center_p_right_2_r(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//DU1rr[dim2]       = compute_center_p_right_2(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//Sr[dim2]         += 1.0*(ns->ed.stsp.par.eta0)*DU1r[dim2];
				//Sr[dim2]         += 2.0*(ns->ed.stsp.par.eta0)*(DU1rr[dim2]+DU12r[dim2]);

				//real Sc 		  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				//real DUc          = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpD[dim][dim2], ns->ed.stn);
				//real DUc1         = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//Sc 		         += 1.0*(ns->ed.stsp.par.eta0)*DUc;
				//Sc 		         += 2.0*(ns->ed.stsp.par.eta0)*DUc1;
				// Compute the tensor derivative
				//real dTdxpR        = compute_dpdx_at_point(cdelta, dim2, 1.0, Sl[dim2], Sr[dim2]);
				dTdxpR        = compute_dpdx_at_point(cdelta, dim2, 0.5, Sl[dim2], Sr[dim2]);
				//ns->cc.dTdx[dim2] = 0.0;
				ns->cc.dTdx[dim2] = dTdxpR;

				//Calculation of the double derivatives
				Sll[dim2]          = compute_center_p_left_2_1(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				//DU1ll[dim2]        = compute_center_p_left_2_1(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpD[dim][dim2], ns->ed.stn);
				//DU1ll[dim2]       = compute_center_p_left_2(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//Sll[dim2]         += 1.0*(ns->ed.stsp.par.eta0)*DU1ll[dim2];
				// Get the tensor in the right cell
				Srr[dim2]          = compute_center_p_right_2_1(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				//DU1rr[dim2]        = compute_center_p_right_2_1(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpD[dim][dim2], ns->ed.stn);
				//DU1rr[dim2]       = compute_center_p_right_2(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//Srr[dim2]         += 1.0*(ns->ed.stsp.par.eta0)*DU1rr[dim2];
				// Compute the tensor derivative
				//real dTdxpL       = compute_dpdx_at_point(cdelta, dim2, 1.0, Sll[dim2], Srr[dim2]);
				dTdxpL       = compute_dpdx_at_point(cdelta, dim2, 0.5, Sll[dim2], Srr[dim2]);

				d2Tdx2p      = compute_dpdx_at_point(cdelta, dim, 0.5, dTdxpL, dTdxpR);
				//ns->cc.d2Tdx2 += 0.0;
				ns->cc.d2Tdx2 	  += d2Tdx2p;
			}
			
		}

		//ns->cc.d2Tdx2 = d2Sdx2;
		//DEBUG_INSPECT(ns->cc.dTdx[0], %lf);
		//DEBUG_INSPECT(ns->cc.dTdx[1], %lf);
		//DEBUG_INSPECT(ns->cc.dTdx[2], %lf);
		//DEBUG_INSPECT(ns->cc.d2Tdx2, %lf);

		//ns->cc.dndx     = compute_dpdx_at_point(cdelta, dim, 1.0, nl, nr);
		//ns->cc.d2ndx2   = compute_facet_du2dx2(cdelta, dim, 1.0, nc, nl, nr);
		//ns->cc.d2psidx2 = compute_facet_du2dx2(cdelta, dim, 1.0, psic, psil, psir);
		//            ns->cc.d2phidx2 = compute_facet_du2dx2(cdelta, dim, 1.0, phic, phil, phir);
		//            ns->cc.d2psidx2 = ns->ed.eo.par.delta*(dp_get_value(ns->ed.eo.dpnminus, clid) - dp_get_value(ns->ed.eo.dpnplus,clid));
		//ns->cc.d2phidx2 = 0.0;
		break;
	}
}


//Computing and storing the temperature value to be included in the Bousinessq term of the momentum equation
void higflow_computational_cell_bousinessq_term_momentum_equation(higflow_solver *ns, int clid, Point ccenter, distributed_property *dpn){
	real Tc, Bs;
	//Thermo parameters
	real Ga   = ns->ed.nif.par.Ga;
	real Gr   = ns->ed.nif.par.Gr;
	real T0   = ns->ed.nif.par.T0;
    real T1   = ns->ed.nif.par.T1;
    real alphaT = ns->ed.nif.par.alphaT;
    real A0 = ns->ed.nif.par.A0;
    real B0 = ns->ed.nif.par.B0;
    real c1 = ns->ed.nif.par.c1;
    real c2 = ns->ed.nif.par.c2;
	// Get the Temperature in the cell center
	Tc = dp_get_value(dpn, clid);
	//Bs = (Ga-Gr*(Tc));
	//If the Bousinessq approximation is activated
	if (ns->ed.nif.contr.bousinessq == 1) {
		//ns->cc.Bous = Bs;
		real Bous   =  dp_get_value(ns->ed.nif.dpEF, clid);
		ns->cc.Bous = Bous;
	} else {
		ns->cc.Bous = 0.0;
	}
	//Define the function of temperature f(T)
	//real Temp = Tc*(T1-T0) + T0;
	real Temp = Tc*(T0-T1) + T1;
	real fTemp = 0.0;
	switch (ns->ed.nif.contr.ftmodel) {
		case 0: 
			//Function f(T): Arrhenius
			fTemp = exp(alphaT*( (1.0/Temp)-(1.0/T0) ));
        	break;
        case 1: 
            //Function f(T): Modified-Arrhenius
            fTemp = exp( -1.0*(alphaT*( Temp-T0 )));
            break;
        case 2: 
            //Function f(T): Vogel-Tamman-Fulcher (VTF)
            fTemp = pow(10.0, B0 + (A0/(Temp-T0)));
            break;
        case 3: 
            //Function f(T): Williams-Landel-Ferry
            fTemp = pow(10.0, -1.0*((c1*(Temp-T0))/(c2 + (Temp-T0))));
            break;
        case 4: 
            //Function f(T): Constant
            fTemp = 1.0;
            break;
    }
	//This is the function that will be used to scale the solvent viscosity with temperature
	ns->cc.FT = fTemp;
}

//Computational cell used to calculate the terms needed to solve the energy equation
void higflow_computational_cell_energy_equation(higflow_solver *ns, sim_domain *sdT, int clid, Point ccenter, Point cdelta, int dim, distributed_property *dpT){
	real psic, phic, psir,psil, phir, phil, nc, nr, nl;
	real Tc, Tcd, Tl, Tld, Tr, Trd;
	real KT0, KTC, KTL, KTR;
	int  incell_l, incell_r;
	real k0 = ns->ed.nif.par.k0;
    real k1 = ns->ed.nif.par.k1;
    real k2 = ns->ed.nif.par.k2;
    real k3 = ns->ed.nif.par.k3;
	real T0 = ns->ed.nif.par.T0;
	real T1 = ns->ed.nif.par.T1;
	switch (ns->contr.spatialdiscrtype) {
	// Second order
	case 0:
		// Get the n in the cell center
		Tc    = dp_get_value(dpT, clid);
		//ns->cc.Temp = Tc;
		// Get the T in the left cell
		Tl    = compute_center_p_left_22(sdT, ccenter, cdelta, dim, 1.0, dpT, ns->ed.stn, &incell_l);
		// Get the n in the right cell
		Tr    = compute_center_p_right_22(sdT, ccenter, cdelta, dim, 1.0, dpT, ns->ed.stn, &incell_r);
		//Get the dimensional temperature T in the cell center
		//Tcd = (Tc)*(T1-T0) + T0;
		Tcd = (Tc)*(T0-T1) + T1;
		//real Tcdc = Tcd - 273.15;
		//Tcd = 468.15;
		//Get the dimensional temperature T in the left cell
		//Tld = (Tl)*(T1-T0) + T0;
		Tld = (Tl)*(T0-T1) + T1;
		//Tld = 463.15;
		//real Tldc = Tld - 273.15;
		//DEBUG_INSPECT(Tld,%lf);
		//Tld = 453.15
		//Get the dimensional temperature T in the right cell
		//Trd = (Tr)*(T1-T0) + T0;
		Trd = (Tr)*(T0-T1) + T1;
		//real Trdc = Trd - 273.15;
		//Calculate the thermal conductivity K at the temperature of reference
		//KT0 = k0 + k1*T0 + k2*(pow(T0,2.0)) + k2*(pow(T0,3.0));
		//real LC = 1.0;
		//LC /= KT0;
		//DEBUG_INSPECT(LC,%lf);
		//DEBUG_INSPECT(KT0,%lf);
		//Calculate the thermal conductivity K at the temperaturate in the cell center
		KTC = k0 + k1*Tcd + k2*(pow(Tcd,2.0)) + k2*(pow(Tcd,3.0));
		//KTC /= 5.0;
		//DEBUG_INSPECT(KTC,%lf);
		//KTC /= KT0;
		//Calculate the thermal conductivity K at the temperaturate in the left cell
		KTL = k0 + k1*Tld + k2*(pow(Tld,2.0)) + k2*(pow(Tld,3.0));
		ns->cc.KTL = KTL;
		//DEBUG_INSPECT(ns->cc.KTL,%lf);
		//KTL /= KT0;
		//Calculate the thermal conductivity K at the temperaturate in the right cell
		KTR = k0 + k1*Trd + k2*(pow(Trd,2.0)) + k2*(pow(Trd,3.0));
		ns->cc.KTR = KTR;
		//DEBUG_INSPECT(ns->cc.KTR,%lf);
		//KTR /= KT0;
		// Compute the necessary derivatives
		ns->cc.Temp      = Tc;
		ns->cc.dKTdx     = compute_dpdx_at_point(cdelta, dim, 1.0, KTL, KTR);
		ns->cc.dTempdx   = compute_dpdx_at_point(cdelta, dim, 1.0, Tl, Tr);
		//ns->cc.d2Tempdx2 = compute_facet_du2dx2(cdelta, dim, 1.0, Tc, Tl, Tr);
		ns->cc.KT = KTC;

		//ns->cc.viscl = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.gn.dpvisc, ns->ed.stn);
		// Get the cell viscolity in the right cell
		//ns->cc.viscr = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.gn.dpvisc, ns->ed.stn);
		Point p;
		POINT_ASSIGN(p, ccenter);
		p[dim]      = ccenter[dim] + cdelta[dim];
		//real ucl     = compute_facet_u_left(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet);
		real Tlr = compute_center_p_left_22(sdT, p, cdelta, dim, 1.0, dpT, ns->ed.stn, &incell_l);
		real Trr = compute_center_p_right_22(sdT, p, cdelta, dim, 1.0, dpT, ns->ed.stn, &incell_r);
		real dTdxjr = compute_dpdx_at_point(cdelta, dim, 1.0, Tlr, Trr);

		p[dim]      =  ccenter[dim] - cdelta[dim];
		real Tll = compute_center_p_left_22(sdT, p, cdelta, dim, 1.0, dpT, ns->ed.stn, &incell_l);
		real Trl = compute_center_p_right_22(sdT, p, cdelta, dim, 1.0, dpT, ns->ed.stn, &incell_r);
		real dTdxjl = compute_dpdx_at_point(cdelta, dim, 1.0, Tll, Trl);

		//ns->cc.du2dx2[dim2] += (ns->cc.viscr*duidxjr - ns->cc.viscl*duidxjl)/fdelta[dim];
		ns->cc.d2Tempdx2 += (KTR*dTdxjr - KTL*dTdxjl)/(2.0*cdelta[dim]);
		//real const0 = 5.0;
		//ns->cc.d2Tempdx2 += const0;

		break;
	}
}