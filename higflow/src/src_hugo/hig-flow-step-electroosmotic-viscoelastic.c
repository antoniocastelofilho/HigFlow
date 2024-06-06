// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Step Electro-osmotic - version 02/02/2018
// *******************************************************************
// *******************************************************************

#include "hig-flow-step-electroosmotic-viscoelastic.h"

//********************************************************************
// Viscoelastic functions
//********************************************************************

// Computing the Kernel Tensor
void higflow_compute_kernel_tensor(higflow_solver *ns) {
    if (ns->contr.flowtype == 3) {
        // Get the cosntants
        real Re   = ns->par.Re;
        real De   = ns->ed.ve.par.De;
        real beta = ns->ed.ve.par.beta;
        real tol  = ns->ed.ve.par.kernel_tol;
        //tol       = 1.0e-3;
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
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
            // Get the velocity derivative tensor Du and the S tensor
            real Du[DIM][DIM], S[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get Du
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                    // Get S
                    S[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                }
            }
            // Calculate the tensor A
            real A[DIM][DIM], D[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    D[i][j] = 0.5*(Du[i][j]+Du[j][i]);
                    A[i][j] = Re*De*S[i][j]/(1.0-beta) + 2.0*De*D[i][j];
                }
                A[i][i] += 1.0;
            }
            // Eige-values and eige-vectors of A
            real R[DIM][DIM], lambda[DIM];
            hig_flow_jacobi(A, lambda, R);
            // Calculate the Kernel tansformation matrix
            real Kernel[DIM][DIM];
            hig_flow_calculate_kernel(ns, lambda, R, Kernel, tol);
            // Store the Kernel Tensor
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                   dp_set_value(ns->ed.ve.dpKernel[i][j], clid, Kernel[i][j]);
                }
            }
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed pressure property
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.ve.dpKernel[i][j]);
            }
        }
    }
}

// Computing the Polymeric Tensor
void higflow_compute_polymeric_tensor(higflow_solver *ns) {
    if (ns->contr.flowtype == 3) {
        // Get the constants
        real Re   = ns->par.Re;
        real De   = ns->ed.ve.par.De;
        real beta = ns->ed.ve.par.beta;
        real tol  = ns->ed.ve.par.kernel_tol;
        real Smax[DIM][DIM], Smin[DIM][DIM];
        for (int i = 0; i < DIM; i++) {
           for (int j = 0; j < DIM; j++) {
               Smax[i][j] = -1.0e16;
               Smin[i][j] =  1.0e16;
           }
        }
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdp);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid    = mp_lookup(mp, hig_get_cid(c));
            // Get the inside/outside inflow point cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Get the velocity derivative tensor Du and the Kernel tensor
            real Kernel[DIM][DIM], Du[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get Du
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                    // Get Kernel
                    Kernel[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpKernel[i][j], ns->ed.stn);
                }
            }
            // Eige-values and eige-vectors of Kernel
            real R[DIM][DIM], lambda[DIM];
            hig_flow_jacobi(Kernel, lambda, R);
            // Calculate the Inverse Kernel tansformation matrix
            real B[DIM][DIM], S[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = i+1; j < DIM; j++) {
                    B[i][j] = 0.0;
                    B[j][i] = 0.0;
                }
                B[i][i] = ns->ed.ve.get_kernel_inverse(i, lambda[i], tol);
            }
            // Calculate A matrix >> A = R B R^t
            // D = 0.5*(Du + Du^t)
            real A[DIM][DIM], D[DIM][DIM];
            hig_flow_matrix_transpose_product(B, R, A);
            // Calculate the tensor S
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    D[i][j] = 0.5*(Du[i][j]+Du[j][i]);
                    S[i][j] = (1.0-beta)*(A[i][j]-2.0*De*D[i][j])/(Re*De);
                }
                S[i][i] += -(1.0-beta)/(Re*De);
            }
            // Store the Polymeric Tensor
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                   real T = S[i][j] + 2.0*(1-beta)*D[i][j]/Re;
                   if (T > Smax[i][j]) Smax[i][j] = T;
                   if (T < Smin[i][j]) Smin[i][j] = T;
                   dp_set_value(ns->ed.ve.dpS[i][j], clid, S[i][j]);
                }
            }
        }
        for (int i = 0; i < DIM; i++) {
           for (int j = 0; j < DIM; j++) {
               // Printing the min and max tensor
               printf("===> %d %d: Tmin = %lf <===> Tmax = %lf <===\n",i,j,Smin[i][j],Smax[i][j]);
           }
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed pressure property
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.ve.dpS[i][j]);
            }
        }
    }
}


// *******************************************************************
// Constitutive Equation Step for the Explicit Euler Method
// *******************************************************************
void higflow_explicit_euler_constitutive_equation(higflow_solver *ns) {
    if (ns->contr.flowtype == 3) {
        // Get the cosntants
        real Re    = ns->par.Re;
        real De    = ns->ed.ve.par.De;
        real beta  = ns->ed.ve.par.beta;
        real tol   = ns->ed.ve.par.kernel_tol;
        real small = 1.0e-14;
        //tol        = 1.0e-3;
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
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
            // Get the velocity derivative tensor Du, S and Kernel tensor
            real Du[DIM][DIM], S[DIM][DIM], Kernel[DIM][DIM];
            // Get the S tensor trace
            real tr = 0.0;
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get Du
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                    // Get S
                    S[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                    // Get Kernel
                    Kernel[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpKernel[i][j], ns->ed.stn);
                }
                tr += S[i][i];
            }
            // Calculate the tensor A
            real A[DIM][DIM], D[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    D[i][j] = 0.5*(Du[i][j]+Du[j][i]);
                    A[i][j] = Re*De*S[i][j]/(1.0-beta) + 2.0*De*D[i][j];
                }
                A[i][i] += 1.0;
            }
            // Eige-values and eige-vectors of A
            real R[DIM][DIM], lambda[DIM];
            hig_flow_jacobi(A, lambda, R);
            // Calculate M matrix >> M = R^t Du R
            real M[DIM][DIM];
            hig_flow_matrix_product(Du, R, M);
            // Calculate Omega matrix >> Omega = R Omega_aux R^t
            real Omega[DIM][DIM];
            hig_flow_calculate_omega(lambda, R, M, Omega, small);
            // Calculate the matrix MM and BB for the model
            real BB[DIM][DIM], MM[DIM][DIM];
            switch (ns->ed.ve.contr.model) {
                case 0: 
                    // Oldroyd Model
                    hig_flow_calculate_m_and_b_oldroyd(ns, lambda, R, M, MM, BB, tol);
                    break;
                case 1: 
                    // Giesekus Model
                    hig_flow_calculate_m_and_b_giesekus(ns, lambda, R, M, MM, BB, tol);
                    break;
                case 2: 
                    // LPTT Model
                    hig_flow_calculate_m_and_b_lptt(ns, tr, lambda, R, M, MM, BB, tol);
                    break;
            }
            // Calculate RHS = Omega Kernel - Kernel Omega + 2BB + MM/De
            real RHS[DIM][DIM];
            hig_flow_kernel_rhs(De, Kernel, Omega, BB, MM, RHS);
            // Get the velocity at cell center 
            real u[DIM], dKdx[DIM];
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            // Solving the Constitutive Equation using the Euler Method
            for (int i = 0; i < DIM; i++) {
                for (int j = i; j < DIM; j++) {
                    // Right hand side equation
                    real rhs = 0.0;
                    switch (ns->ed.ve.contr.convecdiscrtype) {
                        case 0: 
                            // Kernel derivative at cell center
                            hig_flow_derivative_kernel_at_center_cell(ns, ccenter, cdelta, i, j, Kernel[i][j], dKdx);
                            for (int dim = 0; dim < DIM; dim++) {
                                //Compute convective tensor term in rhs
                                rhs -= u[dim]*dKdx[dim];
                            }
                            break;
                        case 1: 
                            //Compute convective tensor term CUBISTA in rhs
                            for (int dim = 0; dim < DIM; dim++) {
                                rhs -= hig_flow_convective_tensor_term_cubista(ns, ns->dpu[dim], ns->ed.sdED, ns->ed.stn, Kernel, ccenter, cdelta, dim, i, j);
                            }
                            break;
                    }
                    // Compute the final rhs
                    rhs         += RHS[i][j];
                    // Compute the Kernel at next time
                    real kernel  = Kernel[i][j] + ns->par.dt * rhs;
                    // Store Kernel in S
                    dp_set_value(ns->ed.ve.dpS[i][j], clid, kernel);
                    if (i != j) {
                         dp_set_value(ns->ed.ve.dpS[j][i], clid, kernel);
                    }
                }
            }
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed pressure property
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.ve.dpS[i][j]);
            }
        }
        // Store the Kernel Tensor
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
            // Get the S tensor and store in Kernel
            real S[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get S
                    S[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                    // Store Kernel
                    dp_set_value(ns->ed.ve.dpKernel[i][j], clid, S[i][j]);
                }
            }
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed pressure property
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.ve.dpKernel[i][j]);
            }
        }
    }
}

// *******************************************************************
// Calculate convective tensor term CUBISTA
// *******************************************************************
real hig_flow_convective_tensor_term_cubista(higflow_solver *ns, distributed_property *dpu, sim_domain *sdp, sim_stencil *stn, real K[DIM][DIM], Point ccenter, Point cdelta, int dim, int i, int j) {
    real  vbar[DIM], dKdx[dim], kr, krr, kl, kll, kc, a, b, c, d, e, tol, fi, conv1,conv2;
    a     = 1.7500;
    b     = 0.3750;
    c     = 0.7500;
    d     = 0.1250;
    e     = 0.2500;
    tol   = 1.0e-14;
    conv1 = 0.0;
    conv2 = 0.0;
    int   incell_r, incell_l, incell_ll, incell_rr, infacet;
    // Get the kernel at center cell
    kc  = K[i][j];
    // Get the low, high, lowlow, highhigh component kernel at center cell
    kl  = compute_center_p_left_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.ve.dpKernel[i][j], ns->ed.stn, &incell_l); 
    kr  = compute_center_p_right_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.ve.dpKernel[i][j], ns->ed.stn, &incell_r); 
    kll = compute_center_p_left_22(ns->ed.sdED, ccenter, cdelta, dim, 2.0, ns->ed.ve.dpKernel[i][j], ns->ed.stn, &incell_ll);
    krr = compute_center_p_right_22(ns->ed.sdED, ccenter, cdelta, dim, 2.0, ns->ed.ve.dpKernel[i][j], ns->ed.stn, &incell_rr);
    // Get the velocity  v1bar(i+1/2,j) in the facet center
    vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
    if (vbar[dim] > 0.0){
        if (fabs(kr - kl) <= tol){
            conv1 = vbar[dim]*kc;
        }else {
            fi = (kc - kl)/(kr - kl);
            if ((fi <= 0.0) || (fi >= 1.0)) {
                conv1 = vbar[dim]*kc;
            }else {
                if (fi < b){ 
                    if (incell_l == 1)                    conv1 = vbar[dim]*(a*kc - c*kl);
                    else                                  conv1 = vbar[dim]*kc;
                }
	        if ((fi >= b) && (fi <= c)){
                    if ((incell_l == 1)&&(incell_r == 1)) conv1 = vbar[dim]*(c*kc + b*kr -d*kl);
                    else                                  conv1 = vbar[dim]*kc;
                }
	        if (fi > c){ 
                    if (incell_r == 1)                    conv1 = vbar[dim]*(e*kc + c*kr);
                    else                                  conv1 = vbar[dim]*kc;
                }
                    
            }    
        }
    //v1bar < 0.0
    }else {
        if ((incell_r == 1) && (incell_rr == 1)){
            if (fabs(kc - krr) <= tol){
                conv1 = vbar[dim]*kr;
            }else {
                fi = (kr- krr)/(kc - krr);
                if ((fi <= 0.0) || (fi >= 1.0)) {
                    conv1 = vbar[dim]*kr;
                }else {
		    if (fi < b) 
                        conv1 = vbar[dim]*(a*kr - c*krr);
                    if ((fi >= b) && (fi <= c))
                        conv1 = vbar[dim]*(c*kr + b*kc -d*krr);
	            if (fi > c) 
                        conv1 = vbar[dim]*(c*kc + e*kr);
                }
            }
        //Return upwind value at boundary
        }else if ((incell_r == 1) && (incell_rr == 0)){
            if (fabs(kc - krr) <= tol){
                conv1 = vbar[dim]*kr;
            }else {
                fi = (kr- krr)/(kc - krr);
                if ((fi <= 0.0) || (fi >= 1.0)) {
                    conv1 = vbar[dim]*kr;
                }else {
		    if (fi <= c) 
                        conv1 = vbar[dim]*kr;
	            if (fi > c) 
                        conv1 = vbar[dim]*(c*kc + e*kr);
                }
            }/*
            vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            if (vbar[dim] > 0.0) conv1 = vbar[dim]*kc;
            else                 conv1 = vbar[dim]*kr;
            vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            if (vbar[dim] > 0.0) conv2 = vbar[dim]*kl;
            else                 conv2 = vbar[dim]*kc;
            return ((conv1 - conv2)/cdelta[dim]); */
        }else {
                vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
                if (vbar[dim] > 0.0) conv1 = vbar[dim]*kc;
                else                 conv1 = vbar[dim]*kc;
                vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
                if (vbar[dim] > 0.0) conv2 = vbar[dim]*kl;
                else                 conv2 = vbar[dim]*kc;
                return ((conv1 - conv2)/cdelta[dim]); 
        }
        
    }
    // Get the velocity  v2bar(i-1/2,j) in the facet center
    vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
    if (vbar[dim] > 0.0){
        if ((incell_l == 1) && (incell_ll == 1)){
            if (fabs(kc-kll) <= tol) {
	        conv2 = vbar[dim]*kl;
            }else {
	        fi = (kl - kll)/(kc - kll);
	        if ((fi <= 0.0) || (fi >= 1.0)) {
	            conv2 = vbar[dim]*kl;
	        }else {
	            if (fi < b)
	                conv2 = vbar[dim]*(a*kl - c*kll);
	            if ((fi >= b) && (fi <= c))
	                conv2 = vbar[dim]*(b*kc + c*kl - d*kll);
	            if (fi > c)  
	                conv2 = vbar[dim]*(c*kc + e*kl);
	        }
	    }
        }else if ((incell_l == 1) && (incell_ll == 0)){
            if (fabs(kc-kll) <= tol) {
	        conv2 = vbar[dim]*kl;
            }else {
	        fi = (kl - kll)/(kc - kll);
	        if ((fi <= 0.0) || (fi >= 1.0)) {
	            conv2 = vbar[dim]*kl;
	        }else {
	            if (fi <= c)
	                conv2 = vbar[dim]*kl;
	            if (fi > c)  
	                conv2 = vbar[dim]*(c*kc + e*kl);
	        }
	    }/*
            vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            if (vbar[dim] > 0.0) conv1 = vbar[dim]*kc;
            else                 conv1 = vbar[dim]*kr;
            vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            if (vbar[dim] > 0.0) conv2 = vbar[dim]*kl;
            else                 conv2 = vbar[dim]*kc;
            return ((conv1 - conv2)/cdelta[dim]); */
       }else {
                vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
                if (vbar[dim] > 0.0) conv1 = vbar[dim]*kc;
                else                 conv1 = vbar[dim]*kr;
                vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
                if (vbar[dim] > 0.0) conv2 = vbar[dim]*kc;
                else                 conv2 = vbar[dim]*kc;
                return ((conv1 - conv2)/cdelta[dim]); 
        } 
    }else {
    //v2bar < 0.0 
        if (fabs(kl - kr) <= tol) {
            conv2 = vbar[dim]*kc;
        }else {
            fi = (kc - kr)/(kl - kr);
            if ((fi <= 0.0) || (fi >= 1.0)) {
                conv2 = vbar[dim]*kc;
            }else {
	        if (fi < b){
                    if (incell_r == 1)                    conv2 = vbar[dim]*(a*kc - c*kr);
                    else                                  conv2 = vbar[dim]*kc;
                }
	        if ((fi >= b) && (fi <= c)){
                    if ((incell_l == 1)&&(incell_r == 1)) conv2 = vbar[dim]*(c*kc + b*kl -d*kr);
                    else                                  conv2 = vbar[dim]*kc;
                }
	        if (fi > c){ 
                    if (incell_l == 1)                    conv2 = vbar[dim]*(e*kc + c*kl);
                    else                                  conv2 = vbar[dim]*kc;
                }
	    }
        }
    }
    return ((conv1-conv2)/cdelta[dim]);
}

// *******************************************************************
// Constitutive Equation Step for the Implicit Euler Method
// *******************************************************************
void higflow_implicit_euler_constitutive_equation(higflow_solver *ns) {
    if (ns->contr.flowtype == 3) {
        // Get the cosntants
        real dt    = ns->par.dt;
        real Re    = ns->par.Re;
        real De    = ns->ed.ve.par.De;
        real beta  = ns->ed.ve.par.beta;
        real tol   = ns->ed.ve.par.kernel_tol;
        real small = 1.0e-14;
        //tol        = 1.0e-3;
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
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
            // Get the velocity derivative tensor Du, S and Kernel tensor
            real Du[DIM][DIM], S[DIM][DIM], Kernel[DIM][DIM];
            real tr = 0.0;
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get Du
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                    // Get S
                    S[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                    // Get Kernel
                    Kernel[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpKernel[i][j], ns->ed.stn);
                }
                tr += S[i][i];
            }
            // Calculate the tensor A
            real A[DIM][DIM], D[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    D[i][j] = 0.5*(Du[i][j]+Du[j][i]);
                    A[i][j] = Re*De*S[i][j]/(1.0-beta) + 2.0*De*D[i][j];
                }
                A[i][i] += 1.0;
            }
            // Eige-values and eige-vectors of A
            real R[DIM][DIM], lambda[DIM];
            hig_flow_jacobi(A, lambda, R);
            // Calculate M matrix >> M = R^t Du R
            real M[DIM][DIM];
            hig_flow_matrix_product(Du, R, M);
            // Calculate Omega matrix >> Omega = R Omega_aux R^t
            real Omega[DIM][DIM];
            hig_flow_calculate_omega(lambda, R, M, Omega, small);
            // Calculate the matrix MM and BB for the model
            real BB[DIM][DIM], MM[DIM][DIM];
            switch (ns->ed.ve.contr.model) {
                case 0: 
                    // Oldroyd Model
                    hig_flow_calculate_m_and_b_oldroyd(ns, lambda, R, M, MM, BB, tol);
                    break;
                case 1: 
                    // Giesekus Model
                    hig_flow_calculate_m_and_b_giesekus(ns, lambda, R, M, MM, BB, tol);
                    break;
                case 2: 
                    // LPTT Model
                    hig_flow_calculate_m_and_b_lptt(ns, tr, lambda, R, M, MM, BB, tol);
                    break;
            }
            // Calculate RHS = Omega Kernel - Kernel Omega + 2BB + MM/De
            real RHS[DIM][DIM];
            hig_flow_implicit_kernel_rhs(De, BB, MM, RHS);//retorna RHS =  2B + M/De
            // Get the velocity at cell center 
            real u[DIM], dKdx[DIM];
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            real b[DIM*DIM];
            real w[DIM*DIM][DIM*DIM+1]; 
            //Assign the linear system to solve the constitutive equation K - dtKOmega + dtOmegaK = rhs 
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Kernel derivative at cell center
                    hig_flow_derivative_kernel_at_center_cell(ns, ccenter, cdelta, i, j, Kernel[i][j], dKdx);
                    // Right hand side equation
                    real rhs = 0.0;
                    switch (ns->ed.ve.contr.convecdiscrtype) {
                        case 0: 
                            // Kernel derivative at cell center
                            hig_flow_derivative_kernel_at_center_cell(ns, ccenter, cdelta, i, j, Kernel[i][j], dKdx);
                            for (int dim = 0; dim < DIM; dim++) {
                                //Compute convective tensor term in rhs
                                rhs -= u[dim]*dKdx[dim];
                            }
                            break;
                        case 1: 
                            //Compute convective tensor term CUBISTA in rhs
                            for (int dim = 0; dim < DIM; dim++) {
                                rhs -= hig_flow_convective_tensor_term_cubista(ns, ns->dpu[dim], ns->ed.sdED, ns->ed.stn, Kernel, ccenter, cdelta, dim, i, j);
                            }
                            break;
                    }
                    // Compute the final rhs
                    rhs         += RHS[i][j];
                    rhs         *= dt;
                    rhs         += Kernel[i][j];
                    b[i*DIM + j] = rhs;
                }
            }
            //Calculate de kronecker product Omega*I - I*Omega
            hig_flow_kernel_system_matrix (w, Omega, dt);
            //Solve the linear system
            hig_flow_solve_system_constitutive_equation (DIM*DIM, w, b);
            // Get the solution of linear system
            for (int i = 0; i < DIM; i++) {
                for (int j = i; j < DIM; j++) {
               	    // Get the value of kernel
                    real kernel = b[i*DIM+j];
                    // Set the value of kernel
                    dp_set_value(ns->ed.ve.dpS[i][j], clid, kernel);
                    if (i != j) {
                        dp_set_value(ns->ed.ve.dpS[j][i], clid, kernel);
                    }
                }
            }  
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed pressure property
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.ve.dpS[i][j]);
            }
        }
        // Store the Kernel Tensor
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
            // Get the S tensor and store in Kernel
            real S[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get S
                    S[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                    // Store Kernel
                    dp_set_value(ns->ed.ve.dpKernel[i][j], clid, S[i][j]);
                }
            }
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed pressure property
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.ve.dpKernel[i][j]);
            }
        }
    }
}

// Calculate the eige-value and eige-vectors using the Jacobi method
void hig_flow_jacobi(real A[DIM][DIM], real d[DIM], real V[DIM][DIM]) {
    real b[DIM];
    real z[DIM];
    int  NMAX = 500;
    int  nrot = 0;
    for (int ip = 0; ip < DIM; ip++) {
       for (int iq = 0; iq < DIM; iq++) {
          V[ip][iq] = 0.0;
       }
       V[ip][ip] = 1.0;
       b[ip]     = A[ip][ip];
       d[ip]     = b[ip];
       z[ip]     = 0.0;
    }
    for (int i = 1; i <= NMAX; i++) {      
       real tresh;
       real sm = 0.0;
       for (int ip = 0;  ip < DIM-1; ip++)
          for (int iq = ip+1; iq < DIM; iq++)
             sm += fabs(A[ip][iq]);
       if (sm == 0.0) return;
       if (i < 4) 
          tresh = 0.2*sm/(DIM*DIM);
       else
          tresh = 0.0;
      for (int ip = 0; ip < DIM-1; ip++) 
         for (int iq = ip+1; iq < DIM; iq++) {
            real g = 100.0*fabs(A[ip][iq]);
            if ((i > 4) && (fabs(d[ip])+g == fabs(d[ip])) && 
                        (fabs(d[iq])+g == fabs(d[iq]))) {
               A[ip][iq] = 0.0;
            } else if (fabs(A[ip][iq]) > tresh) {
               real t, theta;
               real h = d[iq] - d[ip];
               if (fabs(h)+g ==  fabs(h)) {
                  t = A[ip][iq]/h;
               } else {
                  theta = 0.5*h/A[ip][iq];
                  t     = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                  if (theta < 0.0) t = -t;
               }
               real c    = 1.0/sqrt(1.0+t*t);
               real s    = t*c;
               real tau  = s/(1.0+c);
               h         = t*A[ip][iq];
               z[ip]     = z[ip] - h;
               z[iq]     = z[iq] + h;
               d[ip]     = d[ip] - h;
               d[iq]     = d[iq] + h;
               A[ip][iq] = 0.0;
               for (int j = 0; j <= ip-1; j++) {
                  g        = A[j][ip];
                  h        = A[j][iq];
                  A[j][ip] = g-s*(h+g*tau);
                  A[j][iq] = h+s*(g-h*tau);
               }
               for (int j = ip+1; j <= iq-1; j++) {
                  g        = A[ip][j];
                  h        = A[j][iq];
                  A[ip][j] = g-s*(h+g*tau);
                  A[j][iq] = h+s*(g-h*tau);
               }
               for (int j = iq+1; j < DIM; j++) {
                  g        = A[ip][j];
                  h        = A[iq][j];
                  A[ip][j] = g-s*(h+g*tau);
                  A[iq][j] = h+s*(g-h*tau);
               } 
               for (int j = 0; j < DIM; j++) {
                  g        = V[j][ip];
                  h        = V[j][iq];
                  V[j][ip] = g-s*(h+g*tau);
                  V[j][iq] = h+s*(g-h*tau);
               }
               nrot++;
           }
        }
        for (int ip = 0; ip < DIM; ip++) {
           b[ip] = b[ip]+z[ip];
           d[ip] = b[ip];
           z[ip] = 0.0;
        }
   }
}

// Calculate the matrix product
void hig_flow_matrix_product (real A[DIM][DIM], real R[DIM][DIM], real B[DIM][DIM]) {
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            B[i][j] = 0.0;
            for (int k = 0; k < DIM; k++) {
                for (int l = 0; l < DIM; l++) {
                    B[i][j] += R[k][i]*A[k][l]*R[l][j];
                }
            }
        }
    }
}

// Calculate the matrix product transpose
void hig_flow_matrix_transpose_product (real A[DIM][DIM], real R[DIM][DIM], real B[DIM][DIM]) {
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            B[i][j] = 0.0;
            for (int k = 0; k < DIM; k++) {
                for (int l = 0; l < DIM; l++) {
                    B[i][j] += R[i][k]*A[k][l]*R[j][l];
                }
            }
        }
    }
}

// Calculate RHS = OK - KO + 2B * M/De
void hig_flow_kernel_rhs (real De, real K[DIM][DIM], real O[DIM][DIM], real B[DIM][DIM], real M[DIM][DIM], real RHS[DIM][DIM]) {
    real OK[DIM][DIM], KO[DIM][DIM];
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            OK[i][j] = 0.0;
            KO[i][j] = 0.0;
            for (int k = 0; k < DIM; k++) {
                OK[i][j] += O[i][k]*K[k][j];
                KO[i][j] += K[i][k]*O[k][j];
            }
        }
    }
    // Calculate RHS = OK - KO + 2B * M/De
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            RHS[i][j] = OK[i][j] - KO[i][j] + 2.0*B[i][j] + M[i][j]/De;
        }
    }
}

// Calculate the Kernel matrix
void hig_flow_calculate_kernel (higflow_solver *ns, real lambda[DIM], real R[DIM][DIM], real Kernel[DIM][DIM], real tol) {
   // Calculate the Kernel tansformation matrix
   real Kernel_aux[DIM][DIM];
   for (int i = 0; i < DIM; i++) {
       for (int j = i+1; j < DIM; j++) {
           Kernel_aux[i][j] = 0.0;
           Kernel_aux[j][i] = 0.0;
       }
       Kernel_aux[i][i] = ns->ed.ve.get_kernel(i, lambda[i], tol);
   }
   // Calculate Kernel matrix >> Kernel = R Kernel_aux R^t
   hig_flow_matrix_transpose_product(Kernel_aux, R, Kernel);
}

// Calculate the Omega matrix
void hig_flow_calculate_omega (real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real Omega[DIM][DIM], real small) {
   // Calculate the Omega and B matrix
   real Omega_aux[DIM][DIM];
   for (int i = 0; i < DIM-1; i++) {
       for (int j = i+1; j < DIM; j++) {
           Omega_aux[i][j] = (M[i][j]*lambda[j]+M[j][i]*lambda[i])/(lambda[j]-lambda[i]+small);
           Omega_aux[j][i] = -Omega_aux[i][j];
       }
       Omega_aux[i][i] = 0.0;
    }
    Omega_aux[DIM-1][DIM-1] = 0.0;
    // Calculate Omega matrix >> Omega = R Omega_aux R^t
    hig_flow_matrix_transpose_product(Omega_aux, R, Omega);
}

// Calculate the matrix MM and BB for Oldroyd model
void hig_flow_calculate_m_and_b_oldroyd (higflow_solver *ns, real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real MM[DIM][DIM], real BB[DIM][DIM], real tol) {
    // Calculate the matrix MM and BB for Oldroyd model
    real M_aux[DIM][DIM], B_aux[DIM][DIM];
    for (int i = 0; i < DIM; i++) {
        for (int j = i+1; j < DIM; j++) {
            B_aux[i][j] = 0.0;
            B_aux[j][i] = 0.0;
            M_aux[i][j] = 0.0;
            M_aux[j][i] = 0.0;
        }
        real jlambda = ns->ed.ve.get_kernel_jacobian(i, lambda[i], tol);
        B_aux[i][i]  = M[i][i]*lambda[i]*jlambda;
        M_aux[i][i]  = (1.0-lambda[i])*jlambda;
    }
    // Calculate Kernel matrix >> BB = R Btilde Lambda JLambda R^t
    hig_flow_matrix_transpose_product(B_aux, R, BB);
    // Calculate Kernel matrix >> MM = R M(Lambda) JLambda R^t
    hig_flow_matrix_transpose_product(M_aux, R, MM);
}

// Calculate the matrix MM and BB for Giesekus model
void hig_flow_calculate_m_and_b_giesekus (higflow_solver *ns, real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real MM[DIM][DIM], real BB[DIM][DIM], real tol) {
    // Calculate the matrix MM and BB for Oldroyd model
    real M_aux[DIM][DIM], B_aux[DIM][DIM];
    for (int i = 0; i < DIM; i++) {
        for (int j = i+1; j < DIM; j++) {
            B_aux[i][j] = 0.0;
            B_aux[j][i] = 0.0;
            M_aux[i][j] = 0.0;
            M_aux[j][i] = 0.0;
        }
        real jlambda = ns->ed.ve.get_kernel_jacobian(i, lambda[i], tol);
        B_aux[i][i]  = M[i][i]*lambda[i]*jlambda;
        real aux     = 1.0-lambda[i];
        M_aux[i][i]  = (aux - ns->ed.ve.par.alpha*aux*aux)*jlambda;
    }
    // Calculate Kernel matrix >> BB = R Btilde Lambda JLambda R^t
    hig_flow_matrix_transpose_product(B_aux, R, BB);
    // Calculate Kernel matrix >> MM = R M(Lambda) JLambda R^t
    hig_flow_matrix_transpose_product(M_aux, R, MM);
}

// Calculate the matrix MM and BB for LPTT model
void hig_flow_calculate_m_and_b_lptt (higflow_solver *ns, real tr, real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real MM[DIM][DIM], real BB[DIM][DIM], real tol) {
    // Calculate the matrix MM and BB for Oldroyd model
    real B[DIM][DIM], jlambda[DIM];
    real M_aux[DIM][DIM], B_aux[DIM][DIM];
    for (int i = 0; i < DIM; i++) {
        for (int j = i+1; j < DIM; j++) {
            B_aux[i][j] = 0.0;
            B_aux[j][i] = 0.0;
        }
        B_aux[i][i]  = M[i][i];
    }
    // Calculate Kernel matrix >> B = R B_aux R^t
    hig_flow_matrix_transpose_product(B_aux, R, B);
    // Calculate the BB and MM matrix
    for (int i = 0; i < DIM; i++) {
        for (int j = i+1; j < DIM; j++) {
            B_aux[i][j] = 0.0;
            B_aux[j][i] = 0.0;
            M_aux[i][j] = 0.0;
            M_aux[j][i] = 0.0;
        }
        jlambda[i]   = ns->ed.ve.get_kernel_jacobian(i, lambda[i], tol);
        B_aux[i][i]  = M[i][i]*lambda[i]*jlambda[i];
        M_aux[i][i]  = (1.0-lambda[i])*(1.0 + ns->ed.ve.par.epsilon*ns->ed.ve.par.De*ns->par.Re*tr/(1.0-ns->ed.ve.par.beta))*jlambda[i];
    }
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            M_aux[i][j] += -2.0*(B[i][j]-B[i][j]*lambda[j])*ns->ed.ve.par.De*ns->ed.ve.par.psi*jlambda[j];
        }
    }
    // Calculate Kernel matrix >> BB = R Btilde Lambda JLambda R^t
    hig_flow_matrix_transpose_product(B_aux, R, BB);
    // Calculate Kernel matrix >> MM = R M(Lambda) JLambda R^t
    hig_flow_matrix_transpose_product(M_aux, R, MM);
}

// Get the velocity at cell center 
void hig_flow_velocity_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, real u[DIM]) {
    for (int dim = 0; dim < DIM; dim++) {
        // Verity if is in facet
        int infacet;
        // Get the velocity in the left facet center
        real ul = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
        // Get the velocity in the right facet center
        real ur = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
        // Setting the velocity at cell center
        u[dim]  = 0.5*(ul + ur);
    }
}

// Get the derivative of Kernel 
void hig_flow_derivative_kernel_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, int i, int j, real Kcenter, real dKdx[DIM]) {
    for (int dim = 0; dim < DIM; dim++) {
        int incell_left, incell_right;
        // Get the Kernel in the left cell
        real Kleft = compute_center_p_left_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.ve.dpKernel[i][j], ns->ed.stn, &incell_left);
        // Get the Kernel in the right cell
        real Kright = compute_center_p_right_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.ve.dpKernel[i][j], ns->ed.stn, &incell_right);
        // Compute the Kernel derivative
        if ((incell_left == 1) && (incell_right == 1)) { 
           dKdx[dim] = compute_dpdx_at_point(cdelta, dim, 1.0, Kleft, Kright);
        } else if (incell_right == 1) { 
           dKdx[dim] = compute_dpdxr_at_point(cdelta, dim, 1.0, Kcenter, Kright);
        } else {
           dKdx[dim] = compute_dpdxl_at_point(cdelta, dim, 1.0, Kleft, Kcenter);
        }
    }
}

//Gauss elimination to solve the constitutive equation
void hig_flow_solve_system_constitutive_equation ( int n, real A[DIM*DIM][DIM*DIM+1], real x[DIM*DIM] ) {
    real s, max, aux;
    int l, k, i, j;
    //Assign the vector b to the expanded matrix of the system Ax = b
    for (j = 0; j < n; j++){
           A[j][n] = x[j];
    }
    // Partial pivot
    for (k = 0; k < n-1; k++){
        max = fabs(A[k][k]);
        l = k;
        for (i = k+1; i < n; i++){
            if (fabs(A[i][k]) > max){
	        max = fabs(A[i][k]);
	        l = i;
             }
        }
        if (l != k){
            for (j = k; j <= n; j++){
	        aux = A[k][j];
	            A[k][j] = A[l][j];
	            A[l][j] = aux;
            }
        }
        // Gauss algorithm
        for (i = k+1; i < n; i++){
                s = A[i][k]/A[k][k];
                A[i][k] = 0;
            for (j = k+1; j <= n; j++){
	        A[i][j] -= s * A[k][j];
            }
        }
    }
    x[n-1] = A[n-1][n] / A[n-1][n-1];
    for (i = n-2; i >= 0; i--){
        s = 0;
        for (j = i+1; j < n; j++){
            s += A[i][j] * x[j];
        }
        x[i] = (A[i][n] - s)/ A[i][i];
    }
}

// Calculate RHS = 2B * M/De
void hig_flow_implicit_kernel_rhs (real De, real B[DIM][DIM], real M[DIM][DIM], real RHS[DIM][DIM]) {
    // Calculate RHS = 2B + M/De
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            RHS[i][j] = 2.0*B[i][j] + M[i][j]/De;
        }
    }
}

// Calculate the matrix product
void hig_flow_kernel_system_matrix (real w[DIM*DIM][DIM*DIM+1], real Omega[DIM][DIM], real dt) {
    real I[DIM][DIM];
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            I[i][j] = 0.0;
	    if (i==j)
                I[i][j] = 1.0;
	}
    }
    for (int i = 0; i < DIM; i++) {
        for (int j = i; j < DIM; j++) {
            for (int k = 0; k < DIM; k++) {
                for (int l = 0; l < DIM; l++) {
                    if (i==j){
                        w[i*DIM + k][j*DIM + l] = I[k][l] - dt*Omega[k][l] - dt*Omega[i][i]*I[k][l];
                    }
                    else{
		        w[i*DIM + k][j*DIM + l] =  dt*Omega[j][i]*I[k][l];
                        w[j*DIM + l][i*DIM + k] = - dt*Omega[j][i]*I[k][l];
                    }
                }
            }
        }
    }
}

// *******************************************************************
// Navier-Stokes step elements
// *******************************************************************

// *******************************************************************
// Ionic transport equation 
// *******************************************************************
void higflow_explicit_euler_ionic_transport_equation_nplus(higflow_solver *ns) {
    if (ns->ed.eo.contr.eo_model == 0) {
        real alphaeo = ns->ed.eo.par.alpha;
        // Get the cosntants
        // Get the local sub-domain for the cells
        sim_domain *sdnplus  = psd_get_local_domain(ns->ed.eo.psdEOnplus);
        sim_domain *sdnminus = psd_get_local_domain(ns->ed.eo.psdEOnminus);
        sim_domain *sdpsi    = psd_get_local_domain(ns->ed.eo.psdEOpsi);
        sim_domain *sdphi    = psd_get_local_domain(ns->ed.eo.psdEOphi);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int i = 0; i < DIM; i++) {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdnplus);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdnplus); !higcit_isfinished(it); higcit_nextcell(it)) {
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
            // Get the velocity at cell center 
            real u[DIM], dnplusdx[DIM], nplus;
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            // compute nplus value at point 
            nplus = compute_value_at_point(ns->ed.eo.sdEOnplus, ccenter, ccenter, 1.0, ns->ed.eo.dpnplus, ns->ed.stn);
            // Solving the Transport Equation using the Euler Method
            // Right hand side equation
            real rhs = 0.0;
            switch (ns->ed.eo.contr.convecdiscrtype) {
                // Central scheme
                case 0: 
                    // Ionic derivative at cell center
                    hig_flow_derivative_nplus_at_center_cell(ns, ccenter, cdelta, nplus, dnplusdx);
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_electroosmotic_ionic(ns, sdnplus, sdpsi, sdphi, clid, ccenter, cdelta, dim, ns->ed.eo.dpnplus, ns->ed.eo.dppsi, ns->ed.eo.dpphi);
                        //Compute convective inonic term in rhs
                        rhs    -= u[dim]*dnplusdx[dim];
                        // Compute the diffusive ionic term rhs
                        rhs    += higflow_diffusive_ionic_term(ns);
                        // Compute the potential ionic term rhs
                        rhs    += higflow_potential_ionic_term(ns);
                    }
                    break;
                // CUBISTA scheme
                case 1: 
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_electroosmotic_ionic(ns, sdnplus, sdpsi, sdphi, clid, ccenter, cdelta, dim, ns->ed.eo.dpnplus, ns->ed.eo.dppsi, ns->ed.eo.dpphi);
                        //Compute convective ionic term CUBISTA in rhs
                        rhs    -= higflow_convective_ionic_term_cubista(ns, ns->dpu[dim], ns->ed.eo.dpnplus, ns->ed.eo.sdEOnplus, ns->ed.stn, nplus, ccenter, cdelta, dim);
                        // Compute the diffusive ionic term rhs
                        rhs    += higflow_diffusive_ionic_term(ns);
                        // Compute the potential ionic term rhs
                        rhs    += higflow_potential_ionic_term(ns);
                    }
                    break;
            }
            // Compute the final value step time
            real newnplus =  nplus + ns->par.dt * rhs;
            // Set property value  
            dp_set_value(ns->ed.eo.dpnplus, clid, newnplus);
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed ionic property
        dp_sync(ns->ed.eo.dpnplus);
    }
}

void higflow_explicit_euler_ionic_transport_equation_nminus(higflow_solver *ns) {
    if (ns->ed.eo.contr.eo_model == 0) {
        real alphaeo = ns->ed.eo.par.alpha;
        // Get the cosntants
        real tol   = ns->ed.ve.par.kernel_tol;
        real small = 1.0e-14;
        // Get the local sub-domain for the cells
        sim_domain *sdnminus = psd_get_local_domain(ns->ed.eo.psdEOnminus);
        sim_domain *sdpsi    = psd_get_local_domain(ns->ed.eo.psdEOpsi);
        sim_domain *sdphi    = psd_get_local_domain(ns->ed.eo.psdEOphi);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int i = 0; i < DIM; i++) {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
        // Get the map for the domain property nminus
        mp_mapper *m = sd_get_domain_mapper(sdnminus);
        // Loop for each cell
        higcit_celliterator *itt;
        for (itt = sd_get_domain_celliterator(sdnminus); !higcit_isfinished(itt); higcit_nextcell(itt)) {
            // Get the cell
            hig_cell *c = higcit_getcell(itt);
            // Get the cell identifier
            int clid    = mp_lookup(m, hig_get_cid(c));
            // Get the center of the cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Get the velocity at cell center 
            real u[DIM], dnminusdx[DIM], nminus;
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            // compute nminus value at point 
            nminus = compute_value_at_point(ns->ed.eo.sdEOnminus, ccenter, ccenter, 1.0, ns->ed.eo.dpnminus, ns->ed.stn);
            // Solving the Transport Equation using the Euler Method
            // Right hand side equation
            real rhs = 0.0;
            switch (ns->ed.eo.contr.convecdiscrtype) {
                // Central scheme
                case 0: 
                    // Ionic derivative at cell center
                    hig_flow_derivative_nminus_at_center_cell(ns, ccenter, cdelta, nminus, dnminusdx);
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_electroosmotic_ionic(ns, sdnminus, sdpsi, sdphi, clid, ccenter, cdelta, dim, ns->ed.eo.dpnminus, ns->ed.eo.dppsi, ns->ed.eo.dpphi);
                        //Compute convective inonic term in rhs
                        rhs    -= u[dim]*dnminusdx[dim];
                        // Compute the diffusive ionic term rhs
                        rhs    += higflow_diffusive_ionic_term(ns);
                        // Compute the potential ionic term rhs
                        rhs    -= higflow_potential_ionic_term(ns);
                    }
                    break;
                // CUBISTA scheme
                case 1: 
                    //Compute convective ionic term CUBISTA in rhs
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_electroosmotic_ionic(ns, sdnminus, sdpsi, sdphi, clid, ccenter, cdelta, dim, ns->ed.eo.dpnminus, ns->ed.eo.dppsi, ns->ed.eo.dpphi);
                        //Compute convective ionic term CUBISTA in rhs
                        rhs    -= higflow_convective_ionic_term_cubista(ns, ns->dpu[dim], ns->ed.eo.dpnminus, ns->ed.eo.sdEOnminus, ns->ed.stn, nminus, ccenter, cdelta, dim);
                        // Compute the diffusive ionic term rhs
                        rhs    += higflow_diffusive_ionic_term(ns);
                        // Compute the potential ionic term rhs
                        rhs    -= higflow_potential_ionic_term(ns);
                    }
                    break;
            }
            // Compute the final value step time
            real newnminus =  nminus + ns->par.dt * rhs;
            // Set ionic property value  
            dp_set_value(ns->ed.eo.dpnminus, clid, newnminus);
        }
        // Destroy the iterator
        higcit_destroy(itt);
        // Sync the ditributed ionic property
        dp_sync(ns->ed.eo.dpnminus);
    }
}

//Semi-implicit ionic transport equation
void higflow_implicit_euler_ionic_transport_equation_nplus(higflow_solver *ns) {
    if (ns->ed.eo.contr.eo_model == 0) {
        // Get the local sub-domain for the cells
        sim_domain *sdnplus  = psd_get_local_domain(ns->ed.eo.psdEOnplus);
        sim_domain *sdpsi    = psd_get_local_domain(ns->ed.eo.psdEOpsi);
        sim_domain *sdphi    = psd_get_local_domain(ns->ed.eo.psdEOphi);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int i = 0; i < DIM; i++) {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdnplus);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdnplus); !higcit_isfinished(it); higcit_nextcell(it)) {
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
            // Get the velocity at cell center 
            real u[DIM], dnplusdx[DIM], nplus;
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            // compute nplus value at point 
            nplus = compute_value_at_point(ns->ed.eo.sdEOnplus, ccenter, ccenter, 1.0, ns->ed.eo.dpnplus, ns->ed.stn);
            // Solving the Transport Equation using the Euler Method
            // Right hand side equation
            real rhs = 0.0;
            switch (ns->ed.eo.contr.convecdiscrtype) {
                // Central scheme
                case 0: 
                    // Ionic derivative at cell center
                    hig_flow_derivative_nplus_at_center_cell(ns, ccenter, cdelta, nplus, dnplusdx);
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_electroosmotic_ionic(ns, sdnplus, sdpsi, sdphi, clid, ccenter, cdelta, dim, ns->ed.eo.dpnplus, ns->ed.eo.dppsi, ns->ed.eo.dpphi);
                        //Compute convective inonic term in rhs
                        rhs    -= u[dim]*dnplusdx[dim];
                        // Compute the potential ionic term rhs
                        rhs    += higflow_potential_ionic_term(ns);
                    }
                    break;
                // CUBISTA scheme
                case 1: 
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_electroosmotic_ionic(ns, sdnplus, sdpsi, sdphi, clid, ccenter, cdelta, dim, ns->ed.eo.dpnplus, ns->ed.eo.dppsi, ns->ed.eo.dpphi);
                        //Compute convective ionic term CUBISTA in rhs
                        rhs    -= higflow_convective_ionic_term_cubista(ns, ns->dpu[dim], ns->ed.eo.dpnplus, ns->ed.eo.sdEOnplus, ns->ed.stn, nplus, ccenter, cdelta, dim);
                        // Compute the potential ionic term rhs
                        rhs    += higflow_potential_ionic_term(ns);
                    }
                    break;
            }
            rhs *= ns->par.dt;
            rhs += ns->cc.ncell;
            // Reset the stencil
            stn_reset(ns->ed.stn);
            // Set the right side of stencil
            stn_set_rhs(ns->ed.stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w  = -ns->par.dt/(ns->ed.eo.par.Pe*cdelta[dim2]*cdelta[dim2]);
                alpha  -= 2.0*w ;
                Point p;
                POINT_ASSIGN(p, ccenter);
                // Stencil point update: right point
                p[dim2] = ccenter[dim2] + cdelta[dim2];
                sd_get_stencil(sdnplus, ccenter, p, w, ns->ed.stn);
                // Stencil point update: left point
                p[dim2] = ccenter[dim2] - cdelta[dim2];
                sd_get_stencil(sdnplus, ccenter, p, w, ns->ed.stn);
            }
            alpha = 1.0 + alpha;
            // Get the stencil
            sd_get_stencil(sdnplus, ccenter, ccenter,alpha, ns->ed.stn);
            // Get the index of the stencil
            int *ids   = psd_stn_get_gids(ns->ed.eo.psdEOnplus, ns->ed.stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->ed.stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->ed.stn);
	    int cgid = psd_get_global_id(ns->ed.eo.psdEOnplus, c);
            // Set the right side of solver linear system
            slv_set_bi(ns->ed.eo.slvnplus, cgid, stn_get_rhs(ns->ed.stn));
            // Set the line of matrix of the solver linear system
            slv_set_Ai(ns->ed.eo.slvnplus, cgid, numelems, ids, vals);
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Assemble the solver
        slv_assemble(ns->ed.eo.slvnplus);
        // Solve the linear system
        slv_solve(ns->ed.eo.slvnplus);
        //Load property from solver
        dp_slv_load_from_solver(ns->ed.eo.dpnplus, ns->ed.eo.slvnplus);
        // Syncing the distributed property
        dp_sync(ns->ed.eo.dpnplus);
    }
}

//Semi-implicit ionic transport equation
void higflow_implicit_euler_ionic_transport_equation_nminus(higflow_solver *ns) {
    if (ns->ed.eo.contr.eo_model == 0) {
        // Get the local sub-domain for the cells
        sim_domain *sdnminus = psd_get_local_domain(ns->ed.eo.psdEOnminus);
        sim_domain *sdpsi    = psd_get_local_domain(ns->ed.eo.psdEOpsi);
        sim_domain *sdphi    = psd_get_local_domain(ns->ed.eo.psdEOphi);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int i = 0; i < DIM; i++) {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
        // Get the map for the domain property nminus
        mp_mapper *m = sd_get_domain_mapper(sdnminus);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdnminus); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid    = mp_lookup(m, hig_get_cid(c));
            // Get the center of the cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Get the velocity at cell center 
            real u[DIM], dnminusdx[DIM], nminus;
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            // compute nminus value at point 
            nminus = compute_value_at_point(ns->ed.eo.sdEOnminus, ccenter, ccenter, 1.0, ns->ed.eo.dpnminus, ns->ed.stn);
            // Right hand side equation
            real rhs = 0.0;
            switch (ns->ed.eo.contr.convecdiscrtype) {
                // Central scheme
                case 0: 
                    // Ionic derivative at cell center
                    hig_flow_derivative_nminus_at_center_cell(ns, ccenter, cdelta, nminus, dnminusdx);
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_electroosmotic_ionic(ns, sdnminus, sdpsi, sdphi, clid, ccenter, cdelta, dim, ns->ed.eo.dpnminus, ns->ed.eo.dppsi, ns->ed.eo.dpphi);
                        //Compute convective inonic term in rhs
                        rhs    -= u[dim]*dnminusdx[dim];
                        // Compute the potential ionic term rhs
                        rhs    -= higflow_potential_ionic_term(ns);
                    }
                    break;
                // CUBISTA scheme
                case 1: 
                    //Compute convective ionic term CUBISTA in rhs
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_electroosmotic_ionic(ns, sdnminus, sdpsi, sdphi, clid, ccenter, cdelta, dim, ns->ed.eo.dpnminus, ns->ed.eo.dppsi, ns->ed.eo.dpphi);
                        //Compute convective ionic term CUBISTA in rhs
                        rhs    -= higflow_convective_ionic_term_cubista(ns, ns->dpu[dim], ns->ed.eo.dpnminus, ns->ed.eo.sdEOnminus, ns->ed.stn, nminus, ccenter, cdelta, dim);
                        // Compute the potential ionic term rhs
                        rhs    -= higflow_potential_ionic_term(ns);
                    }
                    break;
            }
            rhs *= ns->par.dt;
            rhs += ns->cc.ncell;
            // Reset the stencil
            stn_reset(ns->ed.stn);
            // Set the right side of stencil
            stn_set_rhs(ns->ed.stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha2 = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w   = -ns->par.dt/(ns->ed.eo.par.Pe*cdelta[dim2]*cdelta[dim2]);
                alpha2  -= 2.0*w ;
                Point p;
                POINT_ASSIGN(p, ccenter);
                // Stencil point update: right point
                p[dim2] = ccenter[dim2] + cdelta[dim2];
                sd_get_stencil(sdnminus, ccenter, p, w, ns->ed.stn);
                // Stencil point update: left point
                p[dim2] = ccenter[dim2] - cdelta[dim2];
                sd_get_stencil(sdnminus, ccenter, p, w, ns->ed.stn);
            }
            alpha2 = 1.0 + alpha2;
            // Get the stencil
            sd_get_stencil(sdnminus, ccenter, ccenter,alpha2, ns->ed.stn);
            // Get the index of the stencil
            int *ids   = psd_stn_get_gids(ns->ed.eo.psdEOnminus, ns->ed.stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->ed.stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->ed.stn);
	    int cgid = psd_get_global_id(ns->ed.eo.psdEOnminus, c);
            // Set the right side of solver linear system
            slv_set_bi(ns->ed.eo.slvnminus, cgid, stn_get_rhs(ns->ed.stn));
            // Set the line of matrix of the solver linear system
            slv_set_Ai(ns->ed.eo.slvnminus, cgid, numelems, ids, vals);
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Assemble the solver
        slv_assemble(ns->ed.eo.slvnminus);
        // Solve the linear system
        slv_solve(ns->ed.eo.slvnminus);
        //Load property from solver
        dp_slv_load_from_solver(ns->ed.eo.dpnminus, ns->ed.eo.slvnminus);
        // Syncing the distributed property
        dp_sync(ns->ed.eo.dpnminus);
    }
}


// *******************************************************************
// Calculate convective ionic term CUBISTA
// *******************************************************************
real higflow_convective_ionic_term_cubista(higflow_solver *ns, distributed_property *dpu, distributed_property *dpn, sim_domain *sdp, sim_stencil *stn, real n, Point ccenter, Point cdelta, int dim) {
    real  vbar[DIM], dKdx[dim], kr, krr, kl, kll, kc, a, b, c, d, e, tol, fi, conv1,conv2;
    a     = 1.7500;
    b     = 0.3750;
    c     = 0.7500;
    d     = 0.1250;
    e     = 0.2500;
    tol   = 1.0e-14;
    conv1 = 0.0;
    conv2 = 0.0;
    int   incell_r, incell_l, incell_ll, incell_rr, infacet;
    // Get the kernel at center cell
    kc  = n;
    // Get the low, high, lowlow, highhigh component kernel at center cell
    kl  = compute_center_p_left_22(sdp, ccenter, cdelta, dim, 1.0, dpn, stn, &incell_l); 
    kr  = compute_center_p_right_22(sdp, ccenter, cdelta, dim, 1.0, dpn, stn, &incell_r); 
    kll = compute_center_p_left_22(sdp, ccenter, cdelta, dim, 2.0, dpn, stn, &incell_ll);
    krr = compute_center_p_right_22(sdp, ccenter, cdelta, dim, 2.0, dpn, stn, &incell_rr);
    // Get the velocity  v1bar(i+1/2,j) in the facet center
    vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
    if (vbar[dim] > 0.0){
        if (fabs(kr - kl) <= tol){
            conv1 = vbar[dim]*kc;
        }else {
            fi = (kc - kl)/(kr - kl);
            if ((fi <= 0.0) || (fi >= 1.0)) {
                conv1 = vbar[dim]*kc;
            }else {
                if (fi < b){ 
                    conv1 = vbar[dim]*(a*kc - c*kl);
                }
	        if ((fi >= b) && (fi <= c)){
                    conv1 = vbar[dim]*(c*kc + b*kr -d*kl);
                }
	        if (fi > c){ 
                    conv1 = vbar[dim]*(e*kc + c*kr);
                }
                    
            }    
        }
    //v1bar < 0.0
    }else {
        if (incell_r == 1){
            if (fabs(kc - krr) <= tol){
                conv1 = vbar[dim]*kr;
            }else {
                fi = (kr- krr)/(kc - krr);
                if ((fi <= 0.0) || (fi >= 1.0)) {
                    conv1 = vbar[dim]*kr;
                }else {
		    if (fi < b) 
                        conv1 = vbar[dim]*(a*kr - c*krr);
                    if ((fi >= b) && (fi <= c))
                        conv1 = vbar[dim]*(c*kr + b*kc -d*krr);
	            if (fi > c) 
                        conv1 = vbar[dim]*(c*kc + e*kr);
                }
            }
        //Return upwind value at boundary
        }else {
                vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
                if (vbar[dim] > 0.0) conv1 = vbar[dim]*kc;
                else                 conv1 = vbar[dim]*kr;
                vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
                if (vbar[dim] > 0.0) conv2 = vbar[dim]*kl;
                else                 conv2 = vbar[dim]*kc;
                return ((conv1 - conv2)/cdelta[dim]); 
        }
        
    }
    // Get the velocity  v2bar(i-1/2,j) in the facet center
    vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
    if (vbar[dim] > 0.0){
        if (incell_l == 1){
            if (fabs(kc-kll) <= tol) {
	        conv2 = vbar[dim]*kl;
            }else {
	        fi = (kl - kll)/(kc - kll);
	        if ((fi <= 0.0) || (fi >= 1.0)) {
	            conv2 = vbar[dim]*kl;
	        }else {
	            if (fi < b)
	                conv2 = vbar[dim]*(a*kl - c*kll);
	            if ((fi >= b) && (fi <= c))
	                conv2 = vbar[dim]*(b*kc + c*kl - d*kll);
	            if (fi > c)  
	                conv2 = vbar[dim]*(c*kc + e*kl);
	        }
	    }
       }else {
                vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
                if (vbar[dim] > 0.0) conv1 = vbar[dim]*kc;
                else                 conv1 = vbar[dim]*kr;
                vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
                if (vbar[dim] > 0.0) conv2 = vbar[dim]*kl;
                else                 conv2 = vbar[dim]*kc;
                return ((conv1 - conv2)/cdelta[dim]); 
        } 
    }else {
    //v2bar < 0.0 
        if (fabs(kl - kr) <= tol) {
            conv2 = vbar[dim]*kc;
        }else {
            fi = (kc - kr)/(kl - kr);
            if ((fi <= 0.0) || (fi >= 1.0)) {
                conv2 = vbar[dim]*kc;
            }else {
	        if (fi < b){
                    conv2 = vbar[dim]*(a*kc - c*kr);
                }
	        if ((fi >= b) && (fi <= c)){
                    conv2 = vbar[dim]*(c*kc + b*kl -d*kr);
                }
	        if (fi > c){ 
                    conv2 = vbar[dim]*(e*kc + c*kl);
                }
	    }
        }
    }
    return ((conv1-conv2)/cdelta[dim]);
}

// Get the derivative of nminus 
void hig_flow_derivative_nminus_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, real ncenter, real dndx[DIM]) {
    for (int dim = 0; dim < DIM; dim++) {
        int incell_left, incell_right;
        // Get the concentration in the left cell
        real nleft = compute_center_p_left_22(ns->ed.eo.sdEOnminus, ccenter, cdelta, dim, 1.0, ns->ed.eo.dpnminus, ns->ed.stn, &incell_left);
        // Get the concentration in the right cell
        real nright = compute_center_p_right_22(ns->ed.eo.sdEOnminus, ccenter, cdelta, dim, 1.0, ns->ed.eo.dpnminus, ns->ed.stn, &incell_left);
        // Compute the concentration derivative
        if ((incell_left == 1) && (incell_right == 1)) { 
           dndx[dim] = compute_dpdx_at_point(cdelta, dim, 1.0, nleft, nright);
        } else if (incell_right == 1) { 
           dndx[dim] = compute_dpdxr_at_point(cdelta, dim, 1.0, ncenter, nright);
        } else {
           dndx[dim] = compute_dpdxl_at_point(cdelta, dim, 1.0, nleft, ncenter);
        }
    }
}

// Get the derivative of nplus 
void hig_flow_derivative_nplus_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, real ncenter, real dndx[DIM]) {
    for (int dim = 0; dim < DIM; dim++) {
        int incell_left, incell_right;
        // Get the concentration in the left cell
        real nleft = compute_center_p_left_22(ns->ed.eo.sdEOnplus, ccenter, cdelta, dim, 1.0, ns->ed.eo.dpnplus, ns->ed.stn, &incell_left);
        // Get the concentration in the right cell
        real nright = compute_center_p_right_22(ns->ed.eo.sdEOnplus, ccenter, cdelta, dim, 1.0, ns->ed.eo.dpnplus, ns->ed.stn, &incell_left);
        // Compute the concentration derivative
        if ((incell_left == 1) && (incell_right == 1)) { 
           dndx[dim] = compute_dpdx_at_point(cdelta, dim, 1.0, nleft, nright);
        } else if (incell_right == 1) { 
           dndx[dim] = compute_dpdxr_at_point(cdelta, dim, 1.0, ncenter, nright);
        } else {
           dndx[dim] = compute_dpdxl_at_point(cdelta, dim, 1.0, nleft, ncenter);
        }
    }
}

// Electroosmotic induced potential psi 
void higflow_electroosmotic_psi(higflow_solver *ns) {
    real alphaeo = ns->ed.eo.par.alpha;
    real delta   = ns->ed.eo.par.delta;
    real nminus, nplus;
    // Get the local sub-domain for the cells
    sim_domain *sdp = psd_get_local_domain(ns->ed.eo.psdEOpsi);
    // Get the map for the domain property nminus
    mp_mapper *m = sd_get_domain_mapper(sdp);
    // Loop for each cell
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Get the cell
        hig_cell *c = higcit_getcell(it);
        // Get the center of the cell
        Point ccenter;
        hig_get_center(c, ccenter);
        // Get the delta of the cell
        Point cdelta;
        hig_get_delta(c, cdelta);
        int clid    = mp_lookup(m, hig_get_cid(c));
        // Reset the stencil
        stn_reset(ns->ed.stn);
        // Initialize rhs
        real rhs = 0.0;
        // Set the right side of stencil
        if (ns->ed.eo.contr.eo_model == 0) {
            // Poisson-Nernst-Planck model 
            // Get the ionic concentration n- at center cell
            nplus    = dp_get_value(ns->ed.eo.dpnplus, clid);
            nminus   = dp_get_value(ns->ed.eo.dpnminus, clid);
            rhs      = delta*(nminus - nplus);
        } 
        // Calculate the point and weight of the stencil
        stn_set_rhs(ns->ed.stn, rhs);
        real alpha = 0.0;
        for(int dim = 0; dim < DIM; dim++) {
            // Stencil weight update
            real w = 1.0/(cdelta[dim]*cdelta[dim]);
            alpha -= 2.0 * w;
            // Stencil point update
            Point p;
            POINT_ASSIGN(p, ccenter);
            // Stencil point update: right point
            p[dim] = ccenter[dim] + cdelta[dim];
            sd_get_stencil(sdp, ccenter, p, w, ns->ed.stn);
            // Stencil point update: left point
            p[dim] = ccenter[dim] - cdelta[dim];
            sd_get_stencil(sdp, ccenter, p, w, ns->ed.stn);
        }
        switch (ns->ed.eo.contr.eo_model) {
            case 1:
               // Poisson-Boltzmann model 
               printf("xxxxxx Please, you must implement this model xxxxxxx\n");
               exit(1);
               break;
            case 2:
               // Debye-Hückel model (solving poisson equation for psi) 
               alpha -= 2.0*alphaeo*delta;
               break;
            case 3:
               // Debye-Hückel model (using the analytic solution for psi)
               printf("???????This model should not solve poisson equation for potential psi!???????\n");
               exit(1);
               break;
        }
        // Get the stencil
        sd_get_stencil(sdp, ccenter, ccenter, alpha, ns->ed.stn);
        // Get the index of the stencil
        int *ids   = psd_stn_get_gids(ns->ed.eo.psdEOpsi, ns->ed.stn);
        // Get the value of the stencil
        real *vals = stn_get_vals(ns->ed.stn);
        // Get the number of elements of the stencil
        int numelems = stn_get_numelems(ns->ed.stn);
        // Get the cell identifier of the cell
        int cgid = psd_get_global_id(ns->ed.eo.psdEOpsi, c);
        // Set the right side of solver linear system
        slv_set_bi(ns->ed.eo.slvpsi, cgid, stn_get_rhs(ns->ed.stn));
        // Set the line of matrix of the solver linear system
        slv_set_Ai(ns->ed.eo.slvpsi, cgid, numelems, ids, vals);
    }
    // Destroy the iterator
    higcit_destroy(it);
    // Assemble the solver
    slv_assemble(ns->ed.eo.slvpsi);
    // Solve the linear system
    slv_solve(ns->ed.eo.slvpsi);
    // Set the solver solution in distributed property
    dp_slv_load_from_solver(ns->ed.eo.dppsi, ns->ed.eo.slvpsi);
    dp_sync(ns->ed.eo.dppsi);
}

// Electroosmotic applied potential phi 
void higflow_electroosmotic_phi(higflow_solver *ns) {
    // Get the local sub-domain for the cells
    sim_domain *sdp = psd_get_local_domain(ns->ed.eo.psdEOphi);
    // Loop for each cell
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Get the cell
        hig_cell *c = higcit_getcell(it);
        // Get the center of the cell
        Point ccenter;
        hig_get_center(c, ccenter);
        // Get the delta of the cell
        Point cdelta;
        hig_get_delta(c, cdelta);
        // Reset the stencil
        stn_reset(ns->ed.stn);
        // Set the right side of stencil
        stn_set_rhs(ns->ed.stn, 0.0);
        // Calculate the point and weight of the stencil
        real alpha = 0.0;
        for(int dim = 0; dim < DIM; dim++) {
            // Stencil weight update
            real w = 1.0/(cdelta[dim]*cdelta[dim]);
            alpha -= 2.0 * w;
            // Stencil point update
            Point p;
            POINT_ASSIGN(p, ccenter);
            // Stencil point update: right point
            p[dim] = ccenter[dim] + cdelta[dim];
            sd_get_stencil(sdp, ccenter, p, w, ns->ed.stn);
            // Stencil point update: left point
            p[dim] = ccenter[dim] - cdelta[dim];
            sd_get_stencil(sdp, ccenter, p, w, ns->ed.stn);
        }
        // Get the stencil
        sd_get_stencil(sdp, ccenter, ccenter, alpha, ns->ed.stn);
        // Get the index of the stencil
        int *ids   = psd_stn_get_gids(ns->ed.eo.psdEOphi, ns->ed.stn);
        // Get the value of the stencil
        real *vals = stn_get_vals(ns->ed.stn);
        // Get the number of elements of the stencil
        int numelems = stn_get_numelems(ns->ed.stn);
        // Get the cell identifier of the cell
        int cgid = psd_get_global_id(ns->ed.eo.psdEOphi, c);
        // Set the right side of solver linear system
        slv_set_bi(ns->ed.eo.slvphi, cgid, stn_get_rhs(ns->ed.stn));
        // Set the line of matrix of the solver linear system
        slv_set_Ai(ns->ed.eo.slvphi, cgid, numelems, ids, vals);
    }
    // Destroy the iterator
    higcit_destroy(it);
    // Assemble the solver
    slv_assemble(ns->ed.eo.slvphi);
    // Solve the linear system
    slv_solve(ns->ed.eo.slvphi);
    // Set the solver solution in distributed property
    dp_slv_load_from_solver(ns->ed.eo.dpphi, ns->ed.eo.slvphi);
    dp_sync(ns->ed.eo.dpphi);
}

// *******************************************************************
// Electro-osmotic source term
// *******************************************************************
void higflow_calculate_electroosmotic_source_term_pnp( higflow_solver *ns) {
    real psil, psir, nplus, nminus;
    real alphaeo= ns->ed.eo.par.alpha;
    real delta  = ns->ed.eo.par.delta;
    // Get the local sub-domain
    sim_domain *sdpsi = psd_get_local_domain(ns->ed.eo.psdEOpsi);
    sim_domain *sdphi  = psd_get_local_domain(ns->ed.eo.psdEOphi);
    sim_facet_domain *sfdF[DIM];
    // Loop for each dimension
    higfit_facetiterator *fit;
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local partitioned domain for facets
        sfdF[dim] = psfd_get_local_domain(ns->ed.eo.psfdEOFeo[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdF[dim]);
        // Loop for each facet
        for(fit = sfd_get_domain_facetiterator(sfdF[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Get the applied potential in the left cell
            psil        = compute_center_p_left(sdphi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpphi, ns->ed.stn);
            // Get the applied potential in the right cell
            psir        = compute_center_p_right(sdphi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpphi, ns->ed.stn);
            // Set the derivative of applied potential
            real phi    = compute_value_at_mid_point(psil, psir);
            real dphidx = compute_dpdx_at_point(fdelta, dim, 0.5, psil, psir);
            // Get the induced potential in the left cell
            psil        = compute_center_p_left(sdpsi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dppsi, ns->ed.stn);
            // Get the induced potential in the right cell
            psir        = compute_center_p_right(sdpsi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dppsi, ns->ed.stn);
            // Set the induced potential
            real psi    = compute_value_at_mid_point(psil, psir);
            // Set the derivative of induced potential
            real dpsidy = compute_dpdx_at_point(fdelta, dim, 0.5, psil, psir);
            // Get the ionic concentration n+ in the left cell
            psil     = compute_center_p_left(ns->ed.eo.sdEOnplus, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpnplus, ns->ed.stn);
            // Get the ionic concentration n+ in the right cell
            psir     = compute_center_p_right(ns->ed.eo.sdEOnplus, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpnplus, ns->ed.stn);
            // Get the ionic concentration n+ at center cell
            nplus    = compute_value_at_mid_point(psil, psir);
            // Get the ionic concentration n- in the left cell
            psil     = compute_center_p_left(ns->ed.eo.sdEOnminus, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpnminus, ns->ed.stn);
            // Get the ionic concentration n- in the right cell
            psir     = compute_center_p_right(ns->ed.eo.sdEOnminus, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpnminus, ns->ed.stn);
            // Get the ionic concentration n- at center cell
            nminus   = compute_value_at_mid_point(psil, psir);
            // Compute the charge density
            real rho = (nplus - nminus)*delta;   
            // Compute the electro-osmotic source term
            real Feo;
            Feo   = -rho*(dphidx + dpsidy); 
            // Get the electroosmotic source term defined by user
            Feo   += ns->ed.eo.get_electroosmotic_source_term(fcenter, dim, ns->par.t);
            // Set the distributed source term property
            dp_set_value(ns->ed.eo.dpFeo[dim], flid, Feo);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the ditributed velocity property
        dp_sync(ns->ed.eo.dpFeo[dim]);
    }
}

void higflow_calculate_electroosmotic_source_term_pb( higflow_solver *ns) {
    real psil, psir;
    real alphaeo= ns->ed.eo.par.alpha;
    real delta  = ns->ed.eo.par.delta;
    real eps    = 7.0e-10;
    real zeta   = 0.025;
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->ed.eo.psdEOpsi);
    sim_domain *sdphi = psd_get_local_domain(ns->ed.eo.psdEOphi);
    sim_facet_domain *sfdF[DIM];
    // Loop for each dimension
    higfit_facetiterator *fit;
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local partitioned domain for facets
        sfdF[dim] = psfd_get_local_domain(ns->ed.eo.psfdEOFeo[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdF[dim]);
        // Loop for each facet
        for(fit = sfd_get_domain_facetiterator(sfdF[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            real y      = fcenter[1];
            // Get the applied potential in the left cell
            psil        = compute_center_p_left(sdphi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpphi, ns->ed.stn);
            // Get the applied potential in the right cell
            psir        = compute_center_p_right(sdphi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpphi, ns->ed.stn);
            // Set the derivative of applied potential
            real dphidx = compute_dpdx_at_point(fdelta, dim, 0.5, psil, psir);
            // Get the induced potential in the left cell
            psil        = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->ed.eo.dppsi, ns->ed.stn);
            // Get the induced potential in the right cell
            psir        = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->ed.eo.dppsi, ns->ed.stn);
            // Set the induced potential
            real psi    = compute_value_at_mid_point(psil, psir);
            // Set the derivative of induced potential
            real dpsidy = compute_dpdx_at_point(fdelta, dim, 0.5, psil, psir);
            // Set the charge density 
            real rho    = -eps*zeta*2.0*alphaeo*delta*psi;
            // Set the electro-osmotic source term
            real Feo;
            if (dim == 0){
                  Feo   = -rho*dphidx;
            }else Feo   = -rho*dpsidy;
            // Get the electroosmotic source term defined by user
            Feo        += ns->ed.eo.get_electroosmotic_source_term(fcenter, dim, ns->par.t);
            // Set the distributed source term property
            dp_set_value(ns->ed.eo.dpFeo[dim], flid, Feo);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the ditributed velocity property
        dp_sync(ns->ed.eo.dpFeo[dim]);
    }
}

void higflow_calculate_electroosmotic_source_term_pbdh( higflow_solver *ns) {
    real psil, psir;
    real alphaeo= ns->ed.eo.par.alpha;
    real delta  = ns->ed.eo.par.delta;
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->ed.eo.psdEOpsi);
    sim_domain *sdphi = psd_get_local_domain(ns->ed.eo.psdEOphi);
    sim_facet_domain *sfdF[DIM];
    // Loop for each dimension
    higfit_facetiterator *fit;
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local partitioned domain for facets
        sfdF[dim] = psfd_get_local_domain(ns->ed.eo.psfdEOFeo[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdF[dim]);
        // Loop for each facet
        for(fit = sfd_get_domain_facetiterator(sfdF[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            real y      = fcenter[1];
            // Get the applied potential in the left cell
            psil        = compute_center_p_left(sdphi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpphi, ns->ed.stn);
            // Get the applied potential in the right cell
            psir        = compute_center_p_right(sdphi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpphi, ns->ed.stn);
            // Compute the applied potential
            real phi    = compute_value_at_mid_point(psil, psir);
            // Compute the derivative of applied potential
            real dphidx = compute_dpdx_at_point(fdelta, dim, 0.5, psil, psir);
            // Get the induced potential in the left cell
            psil        = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->ed.eo.dppsi, ns->ed.stn);
            // Get the induced potential in the right cell
            psir        = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->ed.eo.dppsi, ns->ed.stn);
            // Compute the induced potential
            real psi    = compute_value_at_mid_point(psil, psir);
            // Set the derivative of induced potential
            real dpsidy = compute_dpdx_at_point(fdelta, dim, 0.5, psil, psir);
            // Set the charge density where  "epsez = epsilon/ez" 
            real rho   = -2.0*alphaeo*delta*psi;
            // Set the electro-osmotic source term
            real Feo;
            Feo   = -rho*(dphidx + dpsidy) ; 
            // Get the electroosmotic source term defined by user
            Feo        += ns->ed.eo.get_electroosmotic_source_term(fcenter, dim, ns->par.t);
            // Set the distributed source term property
            dp_set_value(ns->ed.eo.dpFeo[dim], flid, Feo);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the ditributed velocity property
        dp_sync(ns->ed.eo.dpFeo[dim]);
    }
}

void higflow_calculate_electroosmotic_source_term_analytic_pbdh( higflow_solver *ns) {
    real eps    = 7.0e-10;
    real E      = -0.912;
    real H      = 1.0e-5;
    real zeta   = 0.025;
    real kappa  = 10.0e5;
    // Get the local sub-domain
    //sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdF[DIM];
    // Loop for each dimension
    higfit_facetiterator *fit;
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local partitioned domain for facets
        sfdF[dim] = psfd_get_local_domain(ns->ed.eo.psfdEOFeo[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdF[dim]);
        // Loop for each facet
        for(fit = sfd_get_domain_facetiterator(sfdF[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            real y   = fcenter[1];
            // Set the induced potential
            real psi = cosh(kappa*H*y)/cosh(kappa*H);
            // Set the charge density 
            real rho = -eps*zeta*kappa*kappa*psi;
            // Set the electro-osmotic source term
            real Feo;
            if(dim == 0) Feo = rho*E;
            else         Feo = 0.0;
            // Get the electroosmotic source term defined by user
            Feo     += ns->ed.eo.get_electroosmotic_source_term(fcenter, dim, ns->par.t);
            // Set the distributed source term property
            dp_set_value(ns->ed.eo.dpFeo[dim], flid, Feo);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the ditributed velocity property
        dp_sync(ns->ed.eo.dpFeo[dim]);
    }
}

// Apply the boundary condition for source term 
void higflow_boundary_condition_for_electroosmotic_source_term(higflow_solver *ns) {
    // Facet iterator
    higcit_celliterator *it;
    // Local sub-domain
    sim_facet_domain *sfdFeo[DIM];
    // For each dimension
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local sub-domain
        sfdFeo[dim]      = psfd_get_local_domain(ns->ed.eo.psfdEOFeo[dim]);
        sim_domain *sd = sfdFeo[dim]->cdom;
        // Get the number of Dirichlet boundary
        int numbcs = sd_get_num_bcs(sd, DIRICHLET);
        // For each Dirichlet boundary
        for (int i = 0; i < numbcs; i++) {
            // Get the Dirichlet boundary
            sim_boundary *bc = sd_get_bc(sd, DIRICHLET, i);
            // Get the id defined by the user
            int userid       = sb_get_userid(bc);
            // Get the value type of the boundary condition
            bc_valuetype valuetype = sb_get_valuetype(bc);
            // Apply the boundary condition if time dependent
            if (valuetype == timedependent) {
                // Get the mapper
                mp_mapper *bm    = sb_get_mapper(bc);
                // For each cell of the boundary
                for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
                    // Get the cell
                    hig_cell *bcell = higcit_getcell(it);
                    // Get the cell center
                    Point bccenter;
                    hig_get_center(bcell, bccenter);
                    // Get the id of the cell
                    int bcgid = mp_lookup(bm, hig_get_cid(bcell));
                    // Set the time to apply the boundary condition
                    real t   = ns->par.t + ns->par.dt;
                    // Get the velocity defined by the user
                    real val = ns->ed.eo.get_boundary_electroosmotic_source_term(userid, bccenter, dim, t);
                    // Set the value
                    sb_set_value(bc, bcgid, val);
                }
                // Destroy the iterator
                higcit_destroy(it);
            }
        }
    }
}

// *******************************************************************
// Navier-Stokes Step for the Explicit Euler Method
// *******************************************************************
void higflow_explicit_euler_intermediate_velocity_electroosmotic(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]) {
    // Get the facet iterator
    higfit_facetiterator *fit;
    // Get the local domain for cell
    sim_domain *sdp    = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++) {
        // Get the local domain for facet cell
        for (int dim2 = 0; dim2 < DIM; dim2++) {
            sfdu[dim2] = psfd_get_local_domain(ns->psfdu[dim2]);
        }
        // Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Set the computational cell
            higflow_computational_cell_electroosmotic(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Electo-osmotic source term
            rhs += higflow_electroosmotic_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Difusive term contribution
            rhs += higflow_difusive_term(ns, fdelta);
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Compute the intermediate velocity
            real ustar = ns->cc.ucell + ns->par.dt * rhs;
            // Update the distributed property intermediate velocity
            dp_set_value(dpustar[dim], flid, ustar);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Syncing the intermediate velocity
        dp_sync(dpustar[dim]);
        // Set the velocity at outflow
	//set_outflow(ns->psfdu[dim], ns->dpustar[dim], 20.0);
    }
}


// *******************************************************************
// Navier-Stokes Step for the Second Order Explicit Runge-Kutta Method
// *******************************************************************
void higflow_explicit_runge_kutta_2_intermediate_velocity_electroosmotic(higflow_solver *ns) {
    // Calculate the auxiliar velocity by the explicit Euler method
    higflow_explicit_euler_intermediate_velocity_electroosmotic(ns, ns->dpu, ns->dpuaux);
    // Calculate the star velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_electroosmotic(ns, ns->dpuaux, ns->dpustar);
    // Calculate the order 2 Runge-Kutta method using the euler method
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Loop for each dimension
    higfit_facetiterator *fit;
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local partitioned domain for facets
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Get the intermediate velocity
            real u     = dp_get_value(ns->dpu[dim], flid);
            real ustar = dp_get_value(ns->dpustar[dim], flid);
            // Compute the intermediate velocity
            real urk2  = 0.5*(u + ustar);
            // Set the final velocity in the distributed velocity property
            dp_set_value(ns->dpustar[dim], flid, urk2);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the ditributed velocity property
        dp_sync(ns->dpustar[dim]);
    }
}

// *******************************************************************
// Navier-Stokes Step for third Order Explicit Runge-Kutta Method
// *******************************************************************
void higflow_explicit_runge_kutta_3_intermediate_velocity_electroosmotic(higflow_solver *ns) {
    // Calculate the auxiliar velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_electroosmotic(ns, ns->dpu, ns->dpuaux);
    // Calculate the second stage velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_electroosmotic(ns, ns->dpuaux, ns->dpustar);
    // Calculate the order 2 Runge-Kutta method using the euler method
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Loop for each dimension
    higfit_facetiterator *fit;
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local partitioned domain for facets
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Get the intermediate velocity
            real u     = dp_get_value(ns->dpu[dim], flid);
            real ustar = dp_get_value(ns->dpustar[dim], flid);
            // Compute the intermediate velocity
            real urk3  = 0.75*u + 0.25*ustar;
            // Set the final velocity in the distributed velocity property
            dp_set_value(ns->dpuaux[dim], flid, urk3);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the ditributed velocity property
        dp_sync(ns->dpuaux[dim]);
    }
    // Calculate the order 2 Runge-Kutta method using the euler method
    higflow_explicit_euler_intermediate_velocity_electroosmotic(ns, ns->dpuaux, ns->dpustar);
    // Loop for each dimension
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local partitioned domain for facets
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
    // Calculate the third stage velocity by the explicit euler method
        for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Get the intermediate velocity
            real u     = dp_get_value(ns->dpu[dim], flid);
            real ustar = dp_get_value(ns->dpustar[dim], flid);
            // Compute the intermediate velocity
            real urk3  = u/3.0 + 2.0*ustar/3.0;
            // Set the final velocity in the distributed velocity property
            dp_set_value(ns->dpustar[dim], flid, urk3);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the ditributed velocity property
        dp_sync(ns->dpustar[dim]);
    }
}

// *******************************************************************
// Navier-Stokes Step for the Implicit Euler Method
// *******************************************************************
void higflow_semi_implicit_euler_intermediate_velocity_electroosmotic(higflow_solver *ns) {
    // Get the facet iterator
    higfit_facetiterator *fit;
    // Get the local domain for cell
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Get the local domain for facet cell
    for (int dim2 = 0; dim2 < DIM; dim2++) {
        sfdu[dim2] = psfd_get_local_domain(ns->psfdu[dim2]);
    }
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++) {
	// Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
	    int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Set the computational cell
            higflow_computational_cell_electroosmotic(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Electo-osmotic source term
            rhs += higflow_electroosmotic_source_term(ns);
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Total contribuition terms by delta t
            rhs *= ns->par.dt;
            // Velocity term contribution
            rhs += ns->cc.ucell;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w  = -ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
                alpha  -= 2.0*w ;
                Point p;
                POINT_ASSIGN(p, fcenter);
                // Stencil point update: right point
                p[dim2] = fcenter[dim2] + fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
                // Stencil point update: left point
                p[dim2] = fcenter[dim2] - fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
            }
            alpha = 1.0 + alpha;
            // Get the stencil
            sfd_get_stencil(sfdu[dim], fcenter, fcenter,alpha, ns->stn);
            // Get the index of the stencil
            int *ids   = psfd_stn_get_gids(ns->psfdu[dim], ns->stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->stn);
	    int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Set the right side of solver linear system
            slv_set_bi(ns->slvu[dim], fgid, stn_get_rhs(ns->stn));
            // Set the line of matrix of the solver linear system
            slv_set_Ai(ns->slvu[dim], fgid, numelems, ids, vals);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Assemble the solver
        slv_assemble(ns->slvu[dim]);
        // Solve the linear system
        slv_solve(ns->slvu[dim]);
        // Gets the values of the solution
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
	    int flid = mp_lookup(mu, hig_get_fid(f));
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Get the value of ustar
            real ustar = slv_get_xi(ns->slvu[dim], fgid);
            // Set the value of ustar
            dp_set_value(ns->dpustar[dim], flid, ustar);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Syncing the intermediate velocity
        dp_sync(ns->dpustar[dim]);
    }
}

// *******************************************************************
// Navier-Stokes Step for the Semi-Implicit Crank-Nicolson Method
// *******************************************************************
void higflow_semi_implicit_crank_nicolson_intermediate_velocity_electroosmotic(higflow_solver *ns) {
    // Get the facet iterator
    higfit_facetiterator *fit;
    // Get the local domain for cell
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Get the local domain for facet cell
    for (int dim2 = 0; dim2 < DIM; dim2++) {
        sfdu[dim2] = psfd_get_local_domain(ns->psfdu[dim2]);
    }
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++) {
	// Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
	    int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Set the computational cell
            higflow_computational_cell_electroosmotic(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Electo-osmotic source term
            rhs += higflow_electroosmotic_source_term(ns);
            // Diffusive term term contribution
            rhs += 0.5*higflow_difusive_term(ns, fdelta);
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Total contribuition terms times delta t
            rhs *= ns->par.dt;
            // Velocity term contribution
            rhs += ns->cc.ucell;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w  = -0.5*ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
                alpha  -= 2.0*w ;
                Point p;
                POINT_ASSIGN(p, fcenter);
                // Stencil point update: right point
                p[dim2] = fcenter[dim2] + fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
                // Stencil point update: left point
                p[dim2] = fcenter[dim2] - fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
            }
            alpha = 1.0 + alpha;
            // Get the stencil
            sfd_get_stencil(sfdu[dim], fcenter, fcenter,alpha, ns->stn);
            // Get the index of the stencil
            int *ids   = psfd_stn_get_gids(ns->psfdu[dim], ns->stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->stn);
	    int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Set the right side of solver linear system
            slv_set_bi(ns->slvu[dim], fgid, stn_get_rhs(ns->stn));
            // Set the line of matrix of the solver linear system
            slv_set_Ai(ns->slvu[dim], fgid, numelems, ids, vals);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Assemble the solver
        slv_assemble(ns->slvu[dim]);
        // Solve the linear system
        slv_solve(ns->slvu[dim]);
        // Gets the values of the solution
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
	    int flid = mp_lookup(mu, hig_get_fid(f));
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Get the value of ustar
            real ustar = slv_get_xi(ns->slvu[dim], fgid);
            // Set the value of ustar
            dp_set_value(ns->dpustar[dim], flid, ustar);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Syncing the intermediate velocity
        dp_sync(ns->dpustar[dim]);
    }
}

// *******************************************************************
// Navier-Stokes Step for the Implicit BDF2 Method
// *******************************************************************
void higflow_semi_implicit_bdf2_intermediate_velocity_electroosmotic(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]) {
    // Firt stage of Tr-BDF2 method
    // Get the facet iterator
    higfit_facetiterator *fit;
    // Get the local domain for cell
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Get the local domain for facet cell
    for (int dim2 = 0; dim2 < DIM; dim2++) {
        sfdu[dim2] = psfd_get_local_domain(ns->psfdu[dim2]);
    }
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++) {
	// Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
	    int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Set the computational cell
            higflow_computational_cell_electroosmotic(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Electo-osmotic source term
            rhs += higflow_electroosmotic_source_term(ns);
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Total contribuition terms times delta t
            rhs *= 0.5*ns->par.dt;
            // Difusive term contribution
            rhs += 0.25*ns->par.dt*higflow_difusive_term(ns, fdelta);
            // Velocity term contribution
            rhs += ns->cc.ucell;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w  = -0.25*ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
                alpha  -= 2.0*w; //divide po 4 para usar regra trapezio em t(n+1/2)
                Point p;
                POINT_ASSIGN(p, fcenter);
                // Stencil point update: right point
                p[dim2] = fcenter[dim2] + fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
                // Stencil point update: left point
                p[dim2] = fcenter[dim2] - fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
            }
            alpha = 1.0 + alpha;
            // Get the stencil
            sfd_get_stencil(sfdu[dim], fcenter, fcenter, alpha, ns->stn);
            // Get the index of the stencil
            int *ids   = psfd_stn_get_gids(ns->psfdu[dim], ns->stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->stn);
	    int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Set the right side of solver linear system
            slv_set_bi(ns->slvu[dim], fgid, stn_get_rhs(ns->stn));
            // Set the line of matrix of the solver linear system
            slv_set_Ai(ns->slvu[dim], fgid, numelems, ids, vals);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Assemble the solver
        slv_assemble(ns->slvu[dim]);
        // Solve the linear system
        slv_solve(ns->slvu[dim]);
        // Gets the values of the solution
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
	    int flid = mp_lookup(mu, hig_get_fid(f));
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Get the value of ustar
            real uaux = slv_get_xi(ns->slvu[dim], fgid);
            // Set the value of ustar
            dp_set_value(ns->dpuaux[dim], flid, uaux);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Syncing the intermediate velocity
        dp_sync(ns->dpuaux[dim]);
    }
    // Second Stage of Tr-BDF2
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++) {
	// Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
	    int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Set the computational cell
            higflow_computational_cell_electroosmotic(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            //Get the uaux
            real uaux = dp_get_value(ns->dpuaux[dim], flid);
            // Right hand side equation
            real rhs = 0.0;
            // Electo-osmotic source term
            rhs += higflow_electroosmotic_source_term(ns);
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Total contribuition terms times delta t
            rhs *= ns->par.dt/3.0;
            rhs += (4.0*uaux - ns->cc.ucell)/3.0;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w  = - 1.0/3.0*ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
                alpha  -=  2.0*w ;
                Point p;
                POINT_ASSIGN(p, fcenter);
                // Stencil point update: right point
                p[dim2] = fcenter[dim2] + fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
                // Stencil point update: left point
                p[dim2] = fcenter[dim2] - fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
            }
            alpha = 1.0 + alpha;
            // Get the stencil
            sfd_get_stencil(sfdu[dim], fcenter, fcenter,alpha, ns->stn);
            // Get the index of the stencil
            int *ids   = psfd_stn_get_gids(ns->psfdu[dim], ns->stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->stn);
	    int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Set the right side of solver linear system
            slv_set_bi(ns->slvu[dim], fgid, stn_get_rhs(ns->stn));
            // Set the line of matrix of the solver linear system
            slv_set_Ai(ns->slvu[dim], fgid, numelems, ids, vals);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Assemble the solver
        slv_assemble(ns->slvu[dim]);
        // Solve the linear system
        slv_solve(ns->slvu[dim]);
        // Get the solution of linear system
        //Vec *vecu = slv_get_solution_vec(ns->slvu[dim]);
        // Gets the values of the solution
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
	    int flid = mp_lookup(mu, hig_get_fid(f));
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Get the value of ustar
            real ustar = slv_get_xi(ns->slvu[dim], flid);
            // Set the value of ustar
            dp_set_value(ns->dpustar[dim], flid, ustar);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Syncing the intermediate velocity
        dp_sync(ns->dpustar[dim]);
    }
}

// One step of the Navier-Stokes the projection method for viscoelastic flow
void higflow_solver_step_electroosmotic_viscoelastic(higflow_solver *ns) {
    // Boundary condition for velocity
    higflow_boundary_condition_for_velocity(ns);
    // Set outflow for velocity
    higflow_outflow_u_step(ns);
    // Calculate the electro-osmotic source term
    switch (ns->ed.eo.contr.eo_model) {
        case 0:
           // Poisson-Nernst-Planck model 
           //higflow_electroosmotic_phi(ns); //phi calculado no ns-example-2d.c
           higflow_implicit_euler_ionic_transport_equation_nplus(ns);
           higflow_implicit_euler_ionic_transport_equation_nminus(ns);
           //higflow_explicit_euler_ionic_transport_equation_nplus(ns);
           //higflow_explicit_euler_ionic_transport_equation_nminus(ns);
           higflow_electroosmotic_psi(ns);
           higflow_calculate_electroosmotic_source_term_pnp(ns); 
           break;
        case 1:
           // Poisson-Boltzmann model 
           //higflow_electroosmotic_phi(ns);
           higflow_electroosmotic_psi(ns);
           higflow_calculate_electroosmotic_source_term_pb(ns); 
           break;
        case 2:
           // Debye-Hückel model (solving poisson equation for psi) 
           //higflow_electroosmotic_phi(ns);
           higflow_electroosmotic_psi(ns);
           higflow_calculate_electroosmotic_source_term_pbdh(ns); 
           break;
        case 3:
           // Debye-Hückel model (using the analytic solution for psi)
           higflow_calculate_electroosmotic_source_term_analytic_pbdh(ns); 
           break;
    }
    // Calculate the intermediate velocity
    switch (ns->contr.tempdiscrtype) {
        case 0:
           // Explicit Euler method
           higflow_explicit_euler_intermediate_velocity_electroosmotic(ns, ns->dpu, ns->dpustar);
           break;
        case 1:
           // Explicit RK2 method
           higflow_explicit_runge_kutta_2_intermediate_velocity_electroosmotic(ns);
           break;
        case 2:
           // Explicit RK3 method
           higflow_explicit_runge_kutta_3_intermediate_velocity_electroosmotic(ns);
           break;
        case 3:
           // Semi-Implicit Euler Method
           higflow_semi_implicit_euler_intermediate_velocity_electroosmotic(ns);
           break;
        case 4:
           // Semi-Implicit Crank-Nicolson Method
           higflow_semi_implicit_crank_nicolson_intermediate_velocity_electroosmotic(ns);
           break;
        case 5:
           // Semi-Implicit Crank-Nicolson Method
           higflow_semi_implicit_bdf2_intermediate_velocity_electroosmotic(ns, ns->dpu, ns->dpustar);
           break;
    }
    // Set outflow for ustar velocity 
    higflow_outflow_ustar_step(ns);
    // Boundary condition for pressure
    higflow_boundary_condition_for_pressure(ns);
    // Calculate the pressure
    higflow_pressure(ns);
    // Calculate the final velocity
    higflow_final_velocity(ns);
    // Boundary condition for velocity
    higflow_boundary_condition_for_velocity(ns);
    // Set outflow for u velocity 
    higflow_outflow_u_step(ns);
    // Calculate the final pressure
    higflow_final_pressure(ns);
    // Calculate the velocity derivative tensor
    higflow_compute_velocity_derivative_tensor(ns);
    // Computing the Kernel Tensor
    higflow_compute_kernel_tensor(ns);
    // Constitutive Equation Step for the Explicit Euler Method
    switch (ns->ed.ve.contr.discrtype) {
        case 0:
           // Explicit method
           higflow_explicit_euler_constitutive_equation(ns);
           break;
        case 1: 
           // Implicit method
           higflow_implicit_euler_constitutive_equation(ns);
           break;
    }
    // Computing the Polymeric Tensor
    higflow_compute_polymeric_tensor(ns);
}

