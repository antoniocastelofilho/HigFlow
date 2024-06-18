// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Step Viscoelastic variable viscosity - version 10/11/2016
// *******************************************************************
// *******************************************************************

#include "hig-flow-step-viscoelastic-variable-viscosity.h"
#include "hig-flow-mittag-leffler.h"
// *******************************************************************
// Constitutive Equations
// *******************************************************************


// Computing the Kernel Tensor
void higflow_compute_kernel_tensor_variable_viscosity(higflow_solver *ns) {
    if ((ns->ed.nn_contr.rheotype == PLM) || (ns->ed.nn_contr.rheotype == THIXOTROPIC)){
        // Get the cosntants
        real Re   = ns->par.Re;
        real De   = ns->ed.vevv.par.De;
        real beta = ns->ed.vevv.par.beta;
        real tol  = ns->ed.vevv.par.kernel_tol;
        real Kmax[DIM][DIM], Kmin[DIM][DIM];
        for (int i = 0; i < DIM; i++) {
           for (int j = 0; j < DIM; j++) {
               Kmax[i][j] = -1.0e16;
               Kmin[i][j] =  1.0e16;
           }
        }
        //tol       = 1.0e-3;
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Get the local sub-domain for the cells (viscosity)
        sim_domain *sdvisc = psd_get_local_domain(ns->ed.vevv.psdVisc);
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
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
                    // Get S
                    S[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vevv.dpS[i][j], ns->ed.stn);
                }
            }
            //Get the viscosity value
            //real eta = 1.0;
            //real eta = compute_value_at_point(ns->ed.vevv.sdVisc, ccenter, ccenter, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            // Calculate the tensor A
            real A[DIM][DIM], D[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    D[i][j] = 0.5*(Du[i][j]+Du[j][i]);
					A[i][j] = Re*De*S[i][j]/(1.0-beta) + 2.0*De*D[i][j];
                    //A[i][j] = Re*De*S[i][j]/(1.0-beta) + 2.0*De*D[i][j];
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
                   if (Kernel[i][j] > Kmax[i][j]) Kmax[i][j] = Kernel[i][j];
                   if (Kernel[i][j] < Kmin[i][j]) Kmin[i][j] = Kernel[i][j];
                   dp_set_value(ns->ed.vevv.dpKernel[i][j], clid, Kernel[i][j]);
                }
            }
        }
        //for (int i = 0; i < DIM; i++) {
        //   for (int j = 0; j < DIM; j++) {
               // Printing the min and max tensor
        //       printf("===> %d %d: Kmin = %lf <===> Kmax = %lf <===\n",i,j,Kmin[i][j],Kmax[i][j]);
        //   }
        //}
        //printf("PAUSE 72\n");getchar();
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed pressure property
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.vevv.dpKernel[i][j]);
            }
        }
    }
}

// Computing the Polymeric Tensor
void higflow_compute_polymeric_tensor_variable_viscosity(higflow_solver *ns) {
    if ((ns->ed.nn_contr.rheotype == PLM)||(ns->ed.nn_contr.rheotype == THIXOTROPIC)) {
        // Get the constants
        real Re   = ns->par.Re;
        real De   = ns->ed.vevv.par.De;
        real beta = ns->ed.vevv.par.beta;
        real tol  = ns->ed.vevv.par.kernel_tol;
        real Smax[DIM][DIM], Smin[DIM][DIM];
        for (int i = 0; i < DIM; i++) {
           for (int j = 0; j < DIM; j++) {
               Smax[i][j] = -1.0e16;
               Smin[i][j] =  1.0e16;
           }
        }
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
         // Get the local sub-domain for the cells (viscosity)
        sim_domain *sdvisc = psd_get_local_domain(ns->ed.vevv.psdVisc);       
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
            int inflowcell, outflowcell;
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
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
                    // Get Kernel
                    Kernel[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vevv.dpKernel[i][j], ns->ed.stn);
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
                B[i][i] = ns->ed.vevv.get_kernel_inverse(i, lambda[i], tol);
            }
            //Get the viscosity value
            //real eta = compute_value_at_point(ns->ed.vevv.sdVisc, ccenter, ccenter, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
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
                   dp_set_value(ns->ed.vevv.dpS[i][j], clid, S[i][j]);
                }
            }
        }
        //for (int i = 0; i < DIM; i++) {
        //   for (int j = 0; j < DIM; j++) {
               // Printing the min and max tensor
        //       printf("===> %d %d: Tmin = %lf <===> Tmax = %lf <===\n",i,j,Smin[i][j],Smax[i][j]);
        //   }
        //}
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed pressure property
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.vevv.dpS[i][j]);
            }
        }
    }
}

// *******************************************************************
// Constitutive Equation Step for the Explicit Euler Method
// *******************************************************************
void higflow_explicit_euler_constitutive_equation_variable_viscosity(higflow_solver *ns) {
    if ((ns->ed.nn_contr.rheotype == PLM)||(ns->ed.nn_contr.rheotype == THIXOTROPIC)) {
        // Get the cosntants
        real Re    = ns->par.Re;
        real De    = ns->ed.vevv.par.De;
        real beta  = ns->ed.vevv.par.beta;
        real tol   = ns->ed.vevv.par.kernel_tol;
        real small = 1.0e-14;
        switch (ns->ed.vevv.contr.model) {
                case GPTT:
            ns->ed.vevv.par.gamma_gptt = tgamma(ns->ed.vevv.par.beta_gptt);
        }
        //tol        = 1.0e-3;
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int i = 0; i < DIM; i++) {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
        // Get the local sub-domain for the cells (viscosity)
        sim_domain *sdvisc = psd_get_local_domain(ns->ed.vevv.psdVisc);       
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
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
                    // Get S
                    S[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vevv.dpS[i][j], ns->ed.stn);
                    // Get Kernel
                    Kernel[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vevv.dpKernel[i][j], ns->ed.stn);
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
            // Calculate the matrix BB 
            real BB[DIM][DIM];
            // Calculate the matrix BB 
            hig_flow_calculate_b (ns, lambda, R, M, BB, tol);
            // Calculate the matrix MM for the model
            real MM[DIM][DIM], M_aux[DIM][DIM];
            switch (ns->ed.vevv.contr.model) {
                case USERSET: 
                    // User Model
                    ns->ed.ve.calculate_m_user(ns->par.Re, ns->ed.vevv.par.De, ns->ed.vevv.par.beta, tr, lambda, R, M, M_aux, tol);
                    break;
                case OLDROYD_B: 
                    // Oldroyd Model
                    hig_flow_calculate_m_oldroyd(ns, lambda, M, M_aux, tol);
                    break;
                case GIESEKUS: 
                    // Giesekus Model
                    hig_flow_calculate_m_giesekus(ns, lambda, M, M_aux, tol);
                    break;
                case LPTT: 
                    // LPTT Model
                    hig_flow_calculate_m_lptt(ns, tr, lambda, M, R, M_aux, tol);
                    break;
                case GPTT: 
                    // GPTT Model
                    hig_flow_calculate_m_gptt(ns, tr, lambda, M, R, M_aux, tol);
                    break;
            }
            // Calculate Kernel matrix >> MM = R M(Lambda) JLambda R^t
            hig_flow_matrix_transpose_product(M_aux, R, MM);
            //Get the viscosity value
            //real eta = compute_value_at_point(sdvisc, ccenter, ccenter, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            real eta = compute_value_at_point(ns->ed.vevv.sdVisc, ccenter, ccenter, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            //Calculate kernel matrix with viscosity MMV = MM/eta (Do a loop here)
            real MMV[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    MMV[i][j] = MM[i][j]/eta;
                }
            }
            // Calculate RHS = Omega Kernel - Kernel Omega + 2BB + MMV/De
            real RHS[DIM][DIM];
            hig_flow_kernel_rhs(De, Kernel, Omega, BB, MMV, RHS);
            // Get the velocity at cell center 
            real u[DIM], dKdx[DIM];
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            // Solving the Constitutive Equation using the Euler Method
            for (int i = 0; i < DIM; i++) {
                for (int j = i; j < DIM; j++) {
                    // Right hand side equation
                    real rhs = 0.0;
                    switch (ns->ed.vevv.contr.convecdiscrtype) {
                        case CELL_UPWIND: 
                            // Kernel derivative at cell center
                            hig_flow_derivative_kernel_at_center_cell(ns, ccenter, cdelta, i, j, Kernel[i][j], dKdx);
                            for (int dim = 0; dim < DIM; dim++) {
                                //Compute convective tensor term in rhs
                                rhs -= u[dim]*dKdx[dim];
                            }
                            break;
                        case CELL_CUBISTA: 
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
                    dp_set_value(ns->ed.vevv.dpS[i][j], clid, kernel);
                    if (i != j) {
                         dp_set_value(ns->ed.vevv.dpS[j][i], clid, kernel);
                    }
                }
            }
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed pressure property
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.vevv.dpS[i][j]);
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
                    S[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vevv.dpS[i][j], ns->ed.stn);
                    // Store Kernel
                    dp_set_value(ns->ed.vevv.dpKernel[i][j], clid, S[i][j]);
                }
            }
        }

        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed pressure property
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.vevv.dpKernel[i][j]);
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
    kl  = compute_center_p_left_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.vevv.dpKernel[i][j], ns->ed.stn, &incell_l); 
    kr  = compute_center_p_right_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.vevv.dpKernel[i][j], ns->ed.stn, &incell_r); 
    kll = compute_center_p_left_22(ns->ed.sdED, ccenter, cdelta, dim, 2.0, ns->ed.vevv.dpKernel[i][j], ns->ed.stn, &incell_ll);
    krr = compute_center_p_right_22(ns->ed.sdED, ccenter, cdelta, dim, 2.0, ns->ed.vevv.dpKernel[i][j], ns->ed.stn, &incell_rr);
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
void higflow_implicit_euler_constitutive_equation_variable_viscosity(higflow_solver *ns) {
    if ((ns->ed.nn_contr.rheotype == PLM)||(ns->ed.nn_contr.rheotype == THIXOTROPIC)) {
        // Get the cosntants
        real dt    = ns->par.dt;
        real Re    = ns->par.Re;
        real De    = ns->ed.vevv.par.De;
        real beta  = ns->ed.vevv.par.beta;
        real tol   = ns->ed.vevv.par.kernel_tol;
        real small = 1.0e-14;
        switch (ns->ed.vevv.contr.model) {
                case GPTT:
            ns->ed.vevv.par.gamma_gptt = tgamma(ns->ed.vevv.par.beta_gptt);
        }
        //tol        = 1.0e-3;
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int i = 0; i < DIM; i++) {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
        // Get the local sub-domain for the cells (viscosity)
        sim_domain *sdvisc = psd_get_local_domain(ns->ed.vevv.psdVisc);  
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
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
                    // Get S
                    S[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vevv.dpS[i][j], ns->ed.stn);
                    // Get Kernel
                    Kernel[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vevv.dpKernel[i][j], ns->ed.stn);
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
            // Calculate the matrix BB
            real BB[DIM][DIM];
            hig_flow_calculate_b (ns, lambda, R, M, BB, tol);
            // Calculate the matrix MM for the model
            real MM[DIM][DIM], M_aux[DIM][DIM];
            switch (ns->ed.vevv.contr.model) {
                case USERSET: 
                    // User Model
                    ns->ed.ve.calculate_m_user(ns->par.Re, ns->ed.vevv.par.De, ns->ed.vevv.par.beta, tr, lambda, R, M, M_aux, tol);
                    break;
                case OLDROYD_B: 
                    // Oldroyd Model
                    hig_flow_calculate_m_oldroyd(ns, lambda, M, M_aux, tol);
                    break;
                case GIESEKUS: 
                    // Giesekus Model
                    hig_flow_calculate_m_giesekus(ns, lambda, M, M_aux, tol);
                    break;
                case LPTT: 
                    // LPTT Model
                    hig_flow_calculate_m_lptt(ns, tr, lambda, M, R, M_aux, tol);
                    break;
                case GPTT: 
                    // GPTT Model
                    hig_flow_calculate_m_gptt(ns, tr, lambda, M, R, M_aux, tol);
                    break;
            }
            // Calculate Kernel matrix >> MM = R M(Lambda) JLambda R^t
            hig_flow_matrix_transpose_product(M_aux, R, MM);
            //Get the viscosity value
            //real eta = compute_value_at_point(sdvisc, ccenter, ccenter, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            real eta = compute_value_at_point(ns->ed.vevv.sdVisc, ccenter, ccenter, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            //Calculate kernel matrix with viscosity MMV = MM/eta (Do a loop here)
            real MMV[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    MMV[i][j] = MM[i][j]/eta;
                }
            }
            // Calculate RHS = Omega Kernel - Kernel Omega + 2BB + MM/De
            real RHS[DIM][DIM];
            hig_flow_implicit_kernel_rhs(De, BB, MMV, RHS);//retorna RHS =  2B + M/De
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
                    switch (ns->ed.vevv.contr.convecdiscrtype) {
                        case CELL_UPWIND: 
                            // Kernel derivative at cell center
                            hig_flow_derivative_kernel_at_center_cell(ns, ccenter, cdelta, i, j, Kernel[i][j], dKdx);
                            for (int dim = 0; dim < DIM; dim++) {
                                //Compute convective tensor term in rhs
                                rhs -= u[dim]*dKdx[dim];
                            }
                            break;
                        case CELL_CUBISTA: 
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
                    dp_set_value(ns->ed.vevv.dpS[i][j], clid, kernel);
                    if (i != j) {
                        dp_set_value(ns->ed.vevv.dpS[j][i], clid, kernel);
                    }
                }
            }  
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed pressure property
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.vevv.dpS[i][j]);
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
                    S[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vevv.dpS[i][j], ns->ed.stn);
                    // Store Kernel
                    dp_set_value(ns->ed.vevv.dpKernel[i][j], clid, S[i][j]);
                }
            }
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed pressure property
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.vevv.dpKernel[i][j]);
            }
        }
    }
}


// *******************************************************************
// Navier-Stokes Step for the Explicit Euler Method
// *******************************************************************
void higflow_explicit_euler_intermediate_velocity_viscoelastic_variable_viscosity(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]) {
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
            higflow_computational_cell_viscoelastic_variable_viscosity(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Difusive term contribution
            rhs += higflow_difusive_term(ns, fdelta);
            // Compute the intermediate velocity
            real ustar = ns->cc.ucell + ns->par.dt * rhs;
            // Update the distributed property intermediate velocity
            dp_set_value(dpustar[dim], flid, ustar);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Syncing the intermediate velocity
        dp_sync(dpustar[dim]);
    }
}


// *******************************************************************
// Navier-Stokes Step for the Second Order Explicit Runge-Kutta Method
// *******************************************************************
void higflow_explicit_runge_kutta_2_intermediate_velocity_viscoelastic_variable_viscosity(higflow_solver *ns) {
    // Calculate the auxiliar velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_viscoelastic_variable_viscosity(ns, ns->dpu, ns->dpuaux);
    // Calculate the star velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_viscoelastic_variable_viscosity(ns, ns->dpuaux, ns->dpustar);
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
void higflow_explicit_runge_kutta_3_intermediate_velocity_viscoelastic_variable_viscosity(higflow_solver *ns) {
    // Calculate the auxiliar velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_viscoelastic_variable_viscosity(ns, ns->dpu, ns->dpuaux);
    // Calculate the second stage velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_viscoelastic_variable_viscosity(ns, ns->dpuaux, ns->dpustar);
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
    higflow_explicit_euler_intermediate_velocity_viscoelastic_variable_viscosity(ns, ns->dpuaux, ns->dpustar);
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
void higflow_semi_implicit_euler_intermediate_velocity_viscoelastic_variable_viscosity(higflow_solver *ns) {
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
            higflow_computational_cell_viscoelastic_variable_viscosity(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Total contribuition terms by delta t
            rhs *= ns->par.dt;
            // Velocity term contribution
            rhs += ns->cc.ucell;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0, wr = 0.0, wl = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w = - ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
                alpha -= 2.0 * w ;
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
            int *ids   = psfd_stn_get_gids(ns->psfdu[dim],ns->stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->stn);
            // Get the cell identifier of the cell
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
void higflow_semi_implicit_crank_nicolson_intermediate_velocity_viscoelastic_variable_viscosity(higflow_solver *ns) {
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
            higflow_computational_cell_viscoelastic_variable_viscosity(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Diffusive term term contribution
            rhs += 0.5 * higflow_difusive_term(ns, fdelta);
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
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
                real w = - 0.5 *  ns->par.dt*(1.0-ns->ed.vevv.par.beta)*(ns->cc.viscr + ns->cc.viscl)/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
                alpha -= 2.0 * w ;
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
            int *ids   = psfd_stn_get_gids(ns->psfdu[dim],ns->stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->stn);
            // Get the cell identifier of the cell
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
void higflow_semi_implicit_bdf2_intermediate_velocity_viscoelastic_variable_viscosity(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]) {
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
            higflow_computational_cell_viscoelastic_variable_viscosity(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Total contribuition terms times delta t
            rhs *= 0.25*ns->par.dt;
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
                real w = - 0.25*ns->par.dt*(1.0-ns->ed.vevv.par.beta)*(ns->cc.viscr + ns->cc.viscl)/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
                alpha -= 2.0 * w ; //divide po 4 para usar regra trapezio em t(n+1/2)
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
            int *ids   = psfd_stn_get_gids(ns->psfdu[dim],ns->stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->stn);
            // Get the cell identifier of the cell
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
    //Second Stage of Tr-BDF2
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
            higflow_computational_cell_viscoelastic_variable_viscosity(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            //Get the uaux
            real uaux = dp_get_value(ns->dpuaux[dim], flid);
            // Right hand side equation
            real rhs = 0.0;
            rhs = (4.0*uaux - ns->cc.ucell)/3.0;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w = - 1.0/3.0* ns->par.dt*(1.0-ns->ed.vevv.par.beta)*(ns->cc.viscr + ns->cc.viscl)/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
                alpha -=  2.0 * w ;
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
            int *ids   = psfd_stn_get_gids(ns->psfdu[dim],ns->stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->stn);
            // Get the cell identifier of the cell
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

// One step of the Navier-Stokes the projection method
void higflow_solver_step_viscoelastic_variable_viscosity(higflow_solver *ns) {
    // Boundary condition for velocity
    higflow_boundary_condition_for_velocity(ns);
    // Calculate the source term
    higflow_calculate_source_term(ns);
    // Calculate the facet source term
    higflow_calculate_facet_source_term(ns);
    // Calculate the intermediated velocity
    switch (ns->contr.tempdiscrtype) {
        case 0:
           // Explicit Euler method
           higflow_explicit_euler_intermediate_velocity_viscoelastic_variable_viscosity(ns, ns->dpu, ns->dpustar);
           break;
        case 1: 
           // Explicit RK2 method
           higflow_explicit_runge_kutta_2_intermediate_velocity_viscoelastic_variable_viscosity(ns);
           break;
        case 2: 
           // Explicit RK3 method
           higflow_explicit_runge_kutta_3_intermediate_velocity_viscoelastic_variable_viscosity(ns);
           break;
        case 3: 
           // Semi-Implicit Euler Method
           higflow_semi_implicit_euler_intermediate_velocity_viscoelastic_variable_viscosity(ns);
           break;
        case 4: 
           // Semi-Implicit Crank-Nicolson Method
           higflow_semi_implicit_crank_nicolson_intermediate_velocity_viscoelastic_variable_viscosity(ns);
           break;
        case 5: 
           // Semi-Implicit Crank-Nicolson Method
           higflow_semi_implicit_bdf2_intermediate_velocity_viscoelastic_variable_viscosity(ns, ns->dpu, ns->dpustar);
           break;
    }
    // Set outflow for ustar velocity 
    //higflow_outflow_ustar_step(ns);
    // Boundary condition for pressure
    higflow_boundary_condition_for_pressure(ns);
    // Calculate the pressure
    higflow_pressure(ns);
    // Calculate the final velocity
    higflow_final_velocity(ns);
    // Boundary condition for velocity
    higflow_boundary_condition_for_velocity(ns);
    // Set outflow for velocity
    //higflow_outflow_u_step(ns);
    // Calculate the final pressure
    higflow_final_pressure(ns);
    // Calculate the velocity derivative tensor
    higflow_compute_velocity_derivative_tensor(ns);
    //Calculate the viscosity (model defined by the user)
    if (ns->ed.nn_contr.rheotype == PLM){
        //User defined model
        higflow_compute_viscosity_user_model_vevv(ns);
        // Computing the Kernel Tensor
        higflow_compute_kernel_tensor_variable_viscosity(ns);
        // Constitutive Equation Step for the Explicit Euler Method
        switch (ns->ed.vevv.contr.discrtype) {
            case EXPLICIT:
            // Explicit method
            higflow_explicit_euler_constitutive_equation_variable_viscosity(ns);
            break;
            case IMPLICIT: 
            // Implicit method
            higflow_implicit_euler_constitutive_equation_variable_viscosity(ns);
            break;
        }
        // Computing the Polymeric Tensor
        higflow_compute_polymeric_tensor_variable_viscosity(ns);
    }
    //Calculate the viscosity (BMP model)
    if (ns->ed.nn_contr.rheotype == THIXOTROPIC){
        //BMP Model
        switch (ns->ed.vevv.contr.structpdiscrtype){
            case EXPLICIT:
                //Explicit method
                higflow_explicit_euler_BMP_viscosity_evolution_equation(ns);
                break;
            case IMPLICIT:
                //Implicit method
                higflow_implicit_euler_BMP_viscosity_evolution_equation(ns);
                break;
        }
        // Computing the Kernel Tensor
        higflow_compute_kernel_tensor_variable_viscosity(ns);
        // Constitutive Equation Step for the Explicit Euler Method
        switch (ns->ed.vevv.contr.discrtype) {
            case EXPLICIT:
                // Explicit method
                higflow_explicit_euler_constitutive_equation_variable_viscosity(ns);
                break;
            case IMPLICIT: 
                // Implicit method
                higflow_implicit_euler_constitutive_equation_variable_viscosity(ns);
                break;
        }
        // Computing the Polymeric Tensor
        higflow_compute_polymeric_tensor_variable_viscosity(ns);
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
       Kernel_aux[i][i] = ns->ed.vevv.get_kernel(i, lambda[i], tol);
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

// Calculate the matrix BB
void hig_flow_calculate_b (higflow_solver *ns, real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real BB[DIM][DIM], real tol) {
    // Calculate the matrix BB 
    real B_aux[DIM][DIM];
    for (int i = 0; i < DIM; i++) {
        for (int j = i+1; j < DIM; j++) {
            B_aux[i][j] = 0.0;
            B_aux[j][i] = 0.0;
        }
        real jlambda = ns->ed.vevv.get_kernel_jacobian(i, lambda[i], tol);
        B_aux[i][i]  = M[i][i]*lambda[i]*jlambda;
    }
    // Calculate Kernel matrix >> BB = R Btilde Lambda JLambda R^t
    hig_flow_matrix_transpose_product(B_aux, R, BB);
}

// Calculate the matrix MM for Oldroyd model
void hig_flow_calculate_m_oldroyd (higflow_solver *ns, real lambda[DIM], real M[DIM][DIM], real M_aux[DIM][DIM], real tol) {
    // Calculate the matrix MM for Oldroyd model
    for (int i = 0; i < DIM; i++) {
        for (int j = i+1; j < DIM; j++) {
            M_aux[i][j] = 0.0;
            M_aux[j][i] = 0.0;
        }
        real jlambda = ns->ed.vevv.get_kernel_jacobian(i, lambda[i], tol);
        M_aux[i][i]  = (1.0-lambda[i])*jlambda;
    }
}

// Calculate the matrix MM for Giesekus model
void hig_flow_calculate_m_giesekus (higflow_solver *ns, real lambda[DIM], real M[DIM][DIM], real M_aux[DIM][DIM], real tol) {
    // Calculate the matrix MM for Giesekus model
    for (int i = 0; i < DIM; i++) {
        for (int j = i+1; j < DIM; j++) {
            M_aux[i][j] = 0.0;
            M_aux[j][i] = 0.0;
        }
        real jlambda = ns->ed.vevv.get_kernel_jacobian(i, lambda[i], tol);
        real aux     = 1.0-lambda[i];
        M_aux[i][i]  = (aux - ns->ed.ve.par.alpha*aux*aux)*jlambda;
    }
}

// Calculate the matrix MM for LPTT model
void hig_flow_calculate_m_lptt (higflow_solver *ns, real tr, real lambda[DIM],  real M[DIM][DIM], real R[DIM][DIM], real M_aux[DIM][DIM], real tol) {
    // Calculate the matrix MM for LPTT model
    real B[DIM][DIM], jlambda[DIM];
    real B_aux[DIM][DIM];
    for (int i = 0; i < DIM; i++) {
        for (int j = i+1; j < DIM; j++) {
            B_aux[i][j] = 0.0;
            B_aux[j][i] = 0.0;
        }
        B_aux[i][i]  = M[i][i];
    }
    // Calculate Kernel matrix >> B = R B_aux R^t
    hig_flow_matrix_transpose_product(B_aux, R, B);
    // Calculate the MM matrix
    for (int i = 0; i < DIM; i++) {
        for (int j = i+1; j < DIM; j++) {
            M_aux[i][j] = 0.0;
            M_aux[j][i] = 0.0;
        }
        jlambda[i]   = ns->ed.vevv.get_kernel_jacobian(i, lambda[i], tol);
        M_aux[i][i]  = (1.0-lambda[i])*(1.0+(ns->ed.vevv.par.epsilon*ns->par.Re*ns->ed.vevv.par.De*tr)/(1.0-ns->ed.vevv.par.beta))*jlambda[i];
    }
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            M_aux[i][j] += -2.0*(B[i][j]-B[i][j]*lambda[j])*ns->ed.vevv.par.De*ns->ed.vevv.par.psi*jlambda[j];
        }
    }
}

// Calculate the matrix MM for GPTT model
void hig_flow_calculate_m_gptt (higflow_solver *ns, real tr, real lambda[DIM],  real M[DIM][DIM], real R[DIM][DIM], real M_aux[DIM][DIM], real tol) {
    // Calculate the matrix MM for LPTT model
    real B[DIM][DIM], jlambda[DIM];
    real B_aux[DIM][DIM];
    for (int i = 0; i < DIM; i++) {
        for (int j = i+1; j < DIM; j++) {
            B_aux[i][j] = 0.0;
            B_aux[j][i] = 0.0;
        }
        B_aux[i][i]  = M[i][i];
    }
    // Calculate Kernel matrix >> B = R B_aux R^t
    hig_flow_matrix_transpose_product(B_aux, R, B);
    // Calculate the MM matrix
    for (int i = 0; i < DIM; i++) {
        for (int j = i+1; j < DIM; j++) {
            M_aux[i][j] = 0.0;
            M_aux[j][i] = 0.0;
        }
        jlambda[i]   = ns->ed.vevv.get_kernel_jacobian(i, lambda[i], tol);
	real alfa1 = ns->ed.ve.par.alpha_gptt;
	real beta1 = ns->ed.ve.par.beta_gptt;
	real gama1 = ns->ed.ve.par.gamma_gptt;
	//real gama1, gama1n; 
	//gama1 = 0.572365;
        //gama1n = lgamma(beta1);
        numc z, mitt;
	z.real = ns->ed.ve.par.epsilon*ns->par.Re*ns->ed.ve.par.De*tr/(1.0-ns->ed.ve.par.beta);
	z.imag = 0.0;
	mitt   =  mlfv(alfa1,beta1, z, 6);
        // Exponêncial
        //real mittreal;
        //mittreal = exp(ns->ed.ve.par.epsilon*ns->par.Re*ns->ed.ve.par.De*tr);
        /*if(gama1==0)
         printf("Gamma1=0");
        else
        printf("===> gamma1: %lf=\n",gama1);
        
        if (abs(mittreal-mitt.real)>1.0)
        printf("===> %d: MittExp = %lf <===> MittML = %lf <===\n",i,mittreal,mitt.real);
        */
        // M da exponencial
        //M_aux[i][i]  = (1.0-lambda[i])*mittreal*jlambda[i];
        // M do GPTT
        M_aux[i][i]  = (1.0-lambda[i])*gama1*mitt.real*jlambda[i];
        // M teste
         //M_aux[i][i]  = (1.0-lambda[i])*mitt.real*jlambda[i];
    }
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            M_aux[i][j] += -2.0*(B[i][j]-B[i][j]*lambda[j])*ns->ed.ve.par.De*ns->ed.ve.par.psi*jlambda[j];
        }
    }
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
        real Kleft = compute_center_p_left_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.vevv.dpKernel[i][j], ns->ed.stn, &incell_left);
        // Get the Kernel in the right cell
        real Kright = compute_center_p_right_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.vevv.dpKernel[i][j], ns->ed.stn, &incell_right);
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
// Calculate the viscosity (with User defined model or BMP model)
// *******************************************************************

// Computing the viscosity (User defined model)
void higflow_compute_viscosity_user_model_vevv(higflow_solver *ns) {
    real etamax = -1.0e16;
    real etamin = 1.0e16;
    real qmax = -1.0e16;
    real qmin = 1.0e16;
    if (ns->contr.flowtype == 6) {
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.vevv.psdVisc);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int dim = 0; dim < DIM; dim++) {
            sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
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
            // Get the velocity derivative tensor
            real Du[DIM][DIM];
            for (int dim = 0; dim < DIM; dim++) {
                for (int dim2 = 0; dim2 < DIM; dim2++) {
                    // Get Du
                    Du[dim][dim2] = compute_value_at_point(ns->ed.vevv.sdVisc, ccenter, ccenter, 1.0, ns->ed.vevv.dpD[dim][dim2], ns->ed.stn);
                }
            }
            // Calculate the rate of deformation tensor
            real D[DIM][DIM];
            for (int dim = 0; dim < DIM; dim++) {
                for (int dim2 = 0; dim2 < DIM; dim2++) {
                    // Calculate the tensor
                    D[dim][dim2] = 0.5*(Du[dim][dim2] + Du[dim2][dim]);
                }
            }
            // Calculate the local shear rate
            real D2[DIM][DIM];
            for (int dim = 0; dim < DIM; dim++) {
                for (int dim2 = 0; dim2 < DIM; dim2++) {
                    D2[dim][dim2] = 0.0;
                    for (int dim3 = 0; dim3 < DIM; dim3++) {
                        D2[dim][dim2] += D[dim][dim3] * D[dim3][dim2];
                    }
                }
            }
            // Local shear rate
            real q = 0.0;
            for (int dim = 0; dim < DIM; dim++) {
                q += D2[dim][dim];
            }
            q         = sqrt(2.0*q);
            //Check to see min and max shear rate values
            if (q > qmax) qmax = q;
            if (q < qmin) qmin = q;
            // Calculate the viscosity
            real visc = ns->ed.vevv.get_viscosity(ccenter, q, ns->par.t, ns->ed.vevv.par.beta, 0.0);
            //Check to see min and max viscosity values
            if (visc > etamax) etamax = visc;
            if (visc < etamin) etamin = visc;
            // Set the viscosity in the distributed viscosity property
            dp_set_value(ns->ed.vevv.dpvisc, clid, visc);
        }
        //Printing the min and max shear rate values
        printf("===> qmin = %lf <===> qmax = %lf <===\n", qmin, qmax);
        //Printing the min and max viscosity values
        printf("===> etamin = %lf <===> etamax = %lf <===\n", etamin, etamax);
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed pressure property
        dp_sync(ns->ed.vevv.dpvisc);
    }
}

// Computing the viscosity of the BMP model using an Euler explicit method
void higflow_explicit_euler_BMP_viscosity_evolution_equation(higflow_solver *ns) {
    if (ns->ed.nn_contr.rheotype == THIXOTROPIC) {
        // Get the cosntants
        real Re    = ns->par.Re;
        real De    = ns->ed.vevv.par.De;
        real beta  = ns->ed.vevv.par.beta;
        real Lambda = ns->ed.vevv.par.Lambda;
        real Phi = ns->ed.vevv.par.Phi;
        real Gamma = ns->ed.vevv.par.Gamma;
        real small = 1.0e-14;
        real etamax = -1.0e16;
        real etamin = 1.0e16;
        real spmax = -1.0e16;
        real spmin = 1.0e16;
        real qmax = -1.0e16;
        real qmin = 1.0e16;
        real RHSmax = -1.0e16;
        real RHSmin = 1.0e16;
        real TDmax = -1.0e16;
        real TDmin = 1.0e16;
        real Dmax[DIM][DIM], Dmin[DIM][DIM];
        for (int i = 0; i < DIM; i++) {
           for (int j = 0; j < DIM; j++) {
               Dmax[i][j] = -1.0e16;
               Dmin[i][j] =  1.0e16;
           }
        }
        real TSmax[DIM][DIM], TSmin[DIM][DIM];
        for (int i = 0; i < DIM; i++) {
           for (int j = 0; j < DIM; j++) {
               TSmax[i][j] = -1.0e16;
               TSmin[i][j] =  1.0e16;
           }
        }
        // Get the local sub-domain for the cells (viscosity)
        sim_domain *sdvisc = psd_get_local_domain(ns->ed.vevv.psdVisc);   
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int i = 0; i < DIM; i++) {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }    
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdvisc);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdvisc); !higcit_isfinished(it); higcit_nextcell(it)) {
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
            // Get the velocity derivative tensor Du, and the tensor S
            real Du[DIM][DIM], S[DIM][DIM];
            // Get the tensors Du and S
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get Du
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
                    // Get S
                    S[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vevv.dpS[i][j], ns->ed.stn);
                }
            }

            // Calculate the deformation tensor D
            real D[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    D[i][j] = 0.5*(Du[i][j]+Du[j][i]);
                    if (D[i][j]  > Dmax[i][j]) Dmax[i][j] = D[i][j];
                    if (D[i][j]  < Dmin[i][j]) Dmin[i][j] = D[i][j];
                }
            }

            // Calculate the local shear rate
            real D2[DIM][DIM];
            for (int dim = 0; dim < DIM; dim++) {
                for (int dim2 = 0; dim2 < DIM; dim2++) {
                    D2[dim][dim2] = 0.0;
                    for (int dim3 = 0; dim3 < DIM; dim3++) {
                        D2[dim][dim2] += D[dim][dim3] * D[dim3][dim2];
                    }
                }
            }

            // Local shear rate
            real q = 0.0;
            for (int dim = 0; dim < DIM; dim++) {
                q += D2[dim][dim];
            }
            q         = sqrt(2.0*q);
            //Check to see min and max shear rate values
            if (q > qmax) qmax = q;
            if (q < qmin) qmin = q;

            //Get the viscosity value at the time "n"
            //real eta = ns->ed.vevv.get_viscosity(ccenter, q, ns->par.t, beta, StructParP);
            //real eta = compute_value_at_point(ns->ed.vevv.sdVisc, ccenter, ccenter, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);

            // Calculate the tensor TS
            real TS[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    //For the BMP with solvent (1) and the NM_T models (4) we take into consideration the solvent contribution
                    if ((ns->ed.vevv.contr.structparmodel == BMP_SOLVENT) || (ns->ed.vevv.contr.structparmodel == NM_T)){
                        TS[i][j] = S[i][j] + 2.0*D[i][j]/Re;
                    } else {
                        TS[i][j] = S[i][j] + 2.0*(1.0-beta)*D[i][j]/Re;
                    }
					//TS[i][j] = S[i][j] + 2.0*D[i][j]/Re;
                    if (TS[i][j]  > TSmax[i][j]) TSmax[i][j] = TS[i][j];
                    if (TS[i][j]  < TSmin[i][j]) TSmin[i][j] = TS[i][j];
                }
            }

            // Calculate the first inner product of the tensors TS and D
            real TD1[DIM][DIM];
            for (int dim = 0; dim < DIM; dim++) {
                for (int dim2 = 0; dim2 < DIM; dim2++) {
                    TD1[dim][dim2] = 0.0;
                    for (int dim3 = 0; dim3 < DIM; dim3++) {
                        //TD1[dim][dim2] += TS[dim][dim3] * D[dim3][dim2];
                        TD1[dim][dim2] += fabs(TS[dim][dim3] * D[dim3][dim2]);
                    }
                }
            }
            // Double inner product of the tensors TS and D
            //real TD = TS[0][1]*D[0][1] + TS[1][0]*D[1][0];
            real TD = 0.0;
            for (int dim = 0; dim < DIM; dim++) {
                TD += TD1[dim][dim];
            }
            if (TD > TDmax) TDmax = TD;
            if (TD < TDmin) TDmin = TD;

            //Calculate the value of the structural parameter at time "n"
            real StructParP = compute_value_at_point(ns->ed.vevv.sdVisc, ccenter, ccenter, 1.0, ns->ed.vevv.dpStructPar, ns->ed.stn);
            real StructParT;
            //Calculate varphi_T (solvent + polymer viscosities) at time "n" (only for the BMP model with solvent)
            if (ns->ed.vevv.contr.structparmodel == BMP_SOLVENT) {
                StructParT = StructParP/(1.0 + (beta/(1.0-beta))*StructParP);
            } else {
                StructParT = StructParP;
            }
            //real StructParT = StructParP/(1.0 + (beta/(1.0-beta))*StructParP);

            // Calculate RHS = (Phi-varphi_T)/Lambda + Gamma*(Re/(1-beta))*(1-varphi_T)*TD 
            //real RHS = hig_flow_structural_par_BMP_model_RHS(Lambda, Phi, Gamma, Re, beta, StructParT, TD);
            real RHS = 0.0;
            //real RHS = hig_flow_structural_par_BMP_model2_RHS(Lambda, Phi, Gamma, Re, De, beta, StructParT, TD);
            switch (ns->ed.vevv.contr.structparmodel) {
                case BMP: 
                    //Original BMP model
                    RHS = hig_flow_structural_par_BMP_model_RHS(Lambda, Phi, Gamma, Re, beta, StructParT, TD);
                    break;
                case BMP_SOLVENT: 
                    //BMP model wit solvent viscosity
                    RHS = hig_flow_structural_par_BMP_model_RHS(Lambda, Phi, Gamma, Re, beta, StructParT, TD);
                    break;
                case MBM: 
                    //MBM model
                    RHS = hig_flow_structural_par_BMP_model2_RHS(Lambda, Phi, Gamma, Re, StructParT, TD);
                    break;
                case NM_TAUP: 
                    // NM_taup model
                    RHS = hig_flow_structural_par_BMP_model3_RHS(Lambda, Phi, Gamma, Re, De, StructParT, TD);
                    break;
                case NM_T: 
                    // NM_T model
                    RHS = hig_flow_structural_par_BMP_model3_RHS(Lambda, Phi, Gamma, Re, De, StructParT, TD);
                    break;
            }

            if (RHS > RHSmax) RHSmax = RHS;
            if (RHS < RHSmin) RHSmin = RHS;
            // Get the velocity at cell center 
            real u[DIM], dSpardx[DIM];
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            // Solving the Fluidity Equation using the Euler Method
            // Right hand side equation
            real rhs = 0.0;
            switch (ns->ed.vevv.contr.structpconvecdiscrtype) {
                case CELL_UPWIND: 
                    // Structural parameter derivative at cell center
                    hig_flow_derivative_structural_parameter_center_cell(ns, ccenter, cdelta, StructParP, dSpardx);
                    for (int dim = 0; dim < DIM; dim++) {
                        //Compute convective structural parameter term in rhs
                        rhs    -= u[dim]*dSpardx[dim];
                    }
                    break;
                case CELL_CUBISTA: 
                    //Compute convective tensor term CUBISTA in rhs of the structural parameter equation
                    for (int dim = 0; dim < DIM; dim++) {
                        rhs -= higflow_convective_BMP_structural_parameter_term_cubista(ns, ns->dpu[dim], ns->ed.vevv.dpStructPar, ns->ed.vevv.sdVisc, ns->ed.stn, StructParP, ccenter, cdelta, dim);
                    }
                    break;
            }
            
            // Compute the final rhs (Convective tern + RHS)
            rhs += RHS;
            // Compute the Structural Parameter at next time
            real structparnew = StructParP + ns->par.dt * rhs;
            //Check to see min and max structural parameter values
            if (structparnew > spmax) spmax = structparnew;
            if (structparnew < spmin) spmin = structparnew;
            // Store Structural Parameter
            dp_set_value(ns->ed.vevv.dpStructPar, clid, structparnew);

            //Calculate the new viscosity
            real viscnew = ns->ed.vevv.get_viscosity(ccenter, q, ns->par.t, beta, structparnew);
            //Check to see min and max viscosity values
            if (viscnew > etamax) etamax = viscnew;
            if (viscnew < etamin) etamin = viscnew;
            // Store the viscosity
            dp_set_value(ns->ed.vevv.dpvisc, clid, viscnew);

        }
        //Printing the min and max deformation tensor values
        //for (int i = 0; i < DIM; i++) {
        //    for (int j = 0; j < DIM; j++) {
        //        // Printing the min and max tensor
        //        printf("===> %d %d: Dmin = %lf <===> Dmax = %lf <===\n",i,j,Dmin[i][j],Dmax[i][j]);
        //    }
        //}
        //Printing the min and max elastic + solvent stress tensor values
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                // Printing the min and max tensor
                printf("===> %d %d: TSmin = %lf <===> TSmax = %lf <===\n",i,j,TSmin[i][j],TSmax[i][j]);
            }
        }
        //Printing the min and max TD values
        printf("===> TDmin = %lf <===> TDmax = %lf <===\n", TDmin, TDmax);
        //Printing the min and max shear rate values
        //printf("===> qmin = %lf <===> qmax = %lf <===\n", qmin, qmax);
        //Printing the min and max viscosity values
        printf("===> etamin = %lf <===> etamax = %lf <===\n", etamin, etamax);
        //Printing the min and max structural parameter values
        printf("===> spmin = %lf <===> spmax = %lf <===\n", spmin, spmax);
        //Printing the min and max RHS values
        //printf("===> RHSmin = %lf <===> RHSmax = %lf <===\n", RHSmin, RHSmax);
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed structural parameter property
        dp_sync(ns->ed.vevv.dpStructPar);
        // Sync the ditributed viscosity property
        dp_sync(ns->ed.vevv.dpvisc);

    }
}

// Calculate RHS = (Phi-varphi_T)/Lambda + Gamma*(Re/(1-beta))*(1-varphi_T)*TD 
real hig_flow_structural_par_BMP_model_RHS (real Lambda, real Phi, real Gamma, real Re, real beta, real SPT, real TD) {
    real F1S, F2S, RHS;
    F1S = (1.0-SPT)/Lambda;
    //F1S = (Phi-SPT)/Lambda;
    F2S = Gamma*(Re/(1.0-beta))*(1.0/Phi-SPT)*TD;
    //F2S = Gamma*(Re/(1.0-beta))*(1.0-SPT)*TD;
    RHS = F1S + F2S;
    return RHS;
}

// Calculate RHS = (Phi-varphi_T)/Lambda + Gamma*(Re/(1-beta))*(1-varphi_T)*TD 
real hig_flow_structural_par_BMP_model2_RHS (real Lambda, real Phi, real Gamma, real Re, real SPT, real TD) {
    real F1S, F2S, RHS;
    F1S = (1.0-SPT)/Lambda;
    //F1S = (Phi-SPT)/Lambda;
    //F2S = Gamma*Re*De*TD;
    F2S = Gamma*Re*TD;
    RHS = F1S + F2S;
    return RHS;
}

// Calculate RHS = (Phi-varphi_T)/Lambda + Gamma*(Re/(1-beta))*(1-varphi_T)*TD 
real hig_flow_structural_par_BMP_model3_RHS (real Lambda, real Phi, real Gamma, real Re, real De, real SPT, real TD) {
    real F1S, F2S, RHS;
    F1S = (1.0-SPT)/Lambda;
    //F1S = (Phi-SPT)/Lambda;
    F2S = Gamma*Re*De*TD;
    //F2S = Gamma*Re*TD;
    RHS = F1S + F2S;
    return RHS;
}

//Calculate the derivative of the structural parameter at the cell center
void hig_flow_derivative_structural_parameter_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, real spcenter, real dspdx[DIM]) {
    for (int dim = 0; dim < DIM; dim++) {
        int incell_left, incell_right;
        // Get the structural parameter value in the left cell
        real spleft = compute_center_p_left_22(ns->ed.vevv.sdVisc, ccenter, cdelta, dim, 1.0, ns->ed.vevv.dpStructPar, ns->ed.stn, &incell_left);
        // Get the structural parameter value in the right cell
        real spright = compute_center_p_right_22(ns->ed.vevv.sdVisc, ccenter, cdelta, dim, 1.0, ns->ed.vevv.dpStructPar, ns->ed.stn, &incell_left);
        // Compute the concentration derivative
           dspdx[dim] = compute_dpdx_at_point(cdelta, dim, 1.0, spleft, spright);
    }
}

// *******************************************************************
// Calculate convective structural parameter term CUBISTA of the BMP model
// *******************************************************************
real higflow_convective_BMP_structural_parameter_term_cubista(higflow_solver *ns, distributed_property *dpu, distributed_property *dpSP, sim_domain *sdp, sim_stencil *stn, real sp, Point ccenter, Point cdelta, int dim) {
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
    // Get the structural parameter at the center cell
    kc  = sp;
    // Get the low, high, lowlow, highhigh component kernel at center cell
    kl  = compute_center_p_left_22(sdp, ccenter, cdelta, dim, 1.0, dpSP, stn, &incell_l); 
    kr  = compute_center_p_right_22(sdp, ccenter, cdelta, dim, 1.0, dpSP, stn, &incell_r); 
    kll = compute_center_p_left_22(sdp, ccenter, cdelta, dim, 2.0, dpSP, stn, &incell_ll);
    krr = compute_center_p_right_22(sdp, ccenter, cdelta, dim, 2.0, dpSP, stn, &incell_rr);
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



// Computing the viscosity of the BMP model using an Euler implicit method
void higflow_implicit_euler_BMP_viscosity_evolution_equation(higflow_solver *ns) {
    if (ns->ed.nn_contr.rheotype == THIXOTROPIC) {
        // Get the cosntants
        real Re    = ns->par.Re;
        real De    = ns->ed.vevv.par.De;
        real beta  = ns->ed.vevv.par.beta;
        real Lambda = ns->ed.vevv.par.Lambda;
        real Phi = ns->ed.vevv.par.Phi;
        real Gamma = ns->ed.vevv.par.Gamma;
        real small = 1.0e-14;
        real etamax = -1.0e16;
        real etamin = 1.0e16;
        real spmax = -1.0e16;
        real spmin = 1.0e16;
        real qmax = -1.0e16;
        real qmin = 1.0e16;
        real RHSmax = -1.0e16;
        real RHSmin = 1.0e16;
        real TDmax = -1.0e16;
        real TDmin = 1.0e16;
        real Dmax[DIM][DIM], Dmin[DIM][DIM];
        for (int i = 0; i < DIM; i++) {
           for (int j = 0; j < DIM; j++) {
               Dmax[i][j] = -1.0e16;
               Dmin[i][j] =  1.0e16;
           }
        }
        real TSmax[DIM][DIM], TSmin[DIM][DIM];
        for (int i = 0; i < DIM; i++) {
           for (int j = 0; j < DIM; j++) {
               TSmax[i][j] = -1.0e16;
               TSmin[i][j] =  1.0e16;
           }
        }
        // Get the local sub-domain for the cells (viscosity)
        sim_domain *sdvisc = psd_get_local_domain(ns->ed.vevv.psdVisc);   
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int i = 0; i < DIM; i++) {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }    
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdvisc);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdvisc); !higcit_isfinished(it); higcit_nextcell(it)) {
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
            // Get the velocity derivative tensor Du, and the tensor S
            real Du[DIM][DIM], S[DIM][DIM];
            // Get the tensors Du and S
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get Du
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
                    // Get S
                    S[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vevv.dpS[i][j], ns->ed.stn);
                }
            }

            // Calculate the deformation tensor D
            real D[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    D[i][j] = 0.5*(Du[i][j]+Du[j][i]);
                    if (D[i][j]  > Dmax[i][j]) Dmax[i][j] = D[i][j];
                    if (D[i][j]  < Dmin[i][j]) Dmin[i][j] = D[i][j];
                }
            }

            // Calculate the local shear rate
            real D2[DIM][DIM];
            for (int dim = 0; dim < DIM; dim++) {
                for (int dim2 = 0; dim2 < DIM; dim2++) {
                    D2[dim][dim2] = 0.0;
                    for (int dim3 = 0; dim3 < DIM; dim3++) {
                        D2[dim][dim2] += D[dim][dim3] * D[dim3][dim2];
                    }
                }
            }

            // Local shear rate
            real q = 0.0;
            for (int dim = 0; dim < DIM; dim++) {
                q += D2[dim][dim];
            }
            q         = sqrt(2.0*q);
            //Check to see min and max shear rate values
            if (q > qmax) qmax = q;
            if (q < qmin) qmin = q;

            //Get the viscosity value at the time "n"
            //real eta = ns->ed.vevv.get_viscosity(ccenter, q, ns->par.t, beta, StructParP);
            //real eta = compute_value_at_point(ns->ed.vevv.sdVisc, ccenter, ccenter, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);

            // Calculate the tensor TS (elastic stress + solvent contribution)
            real TS[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    TS[i][j] = S[i][j] + 2.0*(1.0-beta)*D[i][j]/Re;
					//TS[i][j] = S[i][j] + 2.0*D[i][j]/Re;
                    if (TS[i][j]  > TSmax[i][j]) TSmax[i][j] = TS[i][j];
                    if (TS[i][j]  < TSmin[i][j]) TSmin[i][j] = TS[i][j];
                }
            }

            // Calculate the first inner product of the tensors TS and D
            real TD1[DIM][DIM];
            for (int dim = 0; dim < DIM; dim++) {
                for (int dim2 = 0; dim2 < DIM; dim2++) {
                    TD1[dim][dim2] = 0.0;
                    for (int dim3 = 0; dim3 < DIM; dim3++) {
                        //TD1[dim][dim2] += TS[dim][dim3] * D[dim3][dim2];
                        TD1[dim][dim2] += fabs(TS[dim][dim3] * D[dim3][dim2]);
                    }
                }
            }
            // Double inner product of the tensors TS and D
            //real TD = TS[0][1]*D[0][1] + TS[1][0]*D[1][0];
            real TD = 0.0;
            for (int dim = 0; dim < DIM; dim++) {
                TD += TD1[dim][dim];
            }
            if (TD > TDmax) TDmax = TD;
            if (TD < TDmin) TDmin = TD;

            //Calculate the value of the structural parameter at time "n"
            real StructParP = compute_value_at_point(ns->ed.vevv.sdVisc, ccenter, ccenter, 1.0, ns->ed.vevv.dpStructPar, ns->ed.stn);
            //Calculate varphi_T (solvent + polymer viscosities) at time "n"
            real StructParT = StructParP;
            //real StructParT = StructParP/(1.0 + (beta/(1.0-beta))*StructParP);

            // Calculate RHS = Gamma*(Re/(1-beta))*TD 
            //real RHS = hig_flow_structural_par_BMP_model_RHS(Lambda, Phi, Gamma, Re, beta, StructParT, TD);
            real RHS = hig_flow_structural_par_BMP_model_implicit_RHS(Gamma, Re, De, beta, TD);
            if (RHS > RHSmax) RHSmax = RHS;
            if (RHS < RHSmin) RHSmin = RHS;
            // Get the velocity at cell center 
            real u[DIM], dSpardx[DIM];
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            // Solving the Fluidity Equation using the Euler Method
            // Right hand side equation
            real rhs = 0.0;
            switch (ns->ed.vevv.contr.structpconvecdiscrtype) {
                case CELL_UPWIND: 
                    // Structural parameter derivative at cell center
                    hig_flow_derivative_structural_parameter_center_cell(ns, ccenter, cdelta, StructParP, dSpardx);
                    for (int dim = 0; dim < DIM; dim++) {
                        //Compute convective structural parameter term in rhs
                        rhs    -= u[dim]*dSpardx[dim];
                    }
                    break;
                case CELL_CUBISTA: 
                    //Compute convective tensor term CUBISTA in rhs of the structural parameter equation
                    for (int dim = 0; dim < DIM; dim++) {
                        rhs -= higflow_convective_BMP_structural_parameter_term_cubista(ns, ns->dpu[dim], ns->ed.vevv.dpStructPar, ns->ed.vevv.sdVisc, ns->ed.stn, StructParP, ccenter, cdelta, dim);
                    }
                    break;
            }
            
            // Compute the final rhs (Convective tern + RHS)
            rhs += 1.0/Lambda + (1.0/Phi)*RHS;
            // Compute the Structural Parameter at next time
            real structparnew = (StructParP + ns->par.dt * rhs)/(1.0 + ns->par.dt * (1.0/Lambda + RHS));
            //Check to see min and max structural parameter values
            if (structparnew > spmax) spmax = structparnew;
            if (structparnew < spmin) spmin = structparnew;
            // Store Structural Parameter
            dp_set_value(ns->ed.vevv.dpStructPar, clid, structparnew);

            //Calculate the new viscosity
            real viscnew = ns->ed.vevv.get_viscosity(ccenter, q, ns->par.t, beta, structparnew);
            //Check to see min and max viscosity values
            if (viscnew > etamax) etamax = viscnew;
            if (viscnew < etamin) etamin = viscnew;
            // Store the viscosity
            dp_set_value(ns->ed.vevv.dpvisc, clid, viscnew);

        }
        //Printing the min and max deformation tensor values
        //for (int i = 0; i < DIM; i++) {
        //    for (int j = 0; j < DIM; j++) {
        //        // Printing the min and max tensor
        //        printf("===> %d %d: Dmin = %lf <===> Dmax = %lf <===\n",i,j,Dmin[i][j],Dmax[i][j]);
        //    }
        //}
        //Printing the min and max elastic + solvent stress tensor values
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                // Printing the min and max tensor
                printf("===> %d %d: TSmin = %lf <===> TSmax = %lf <===\n",i,j,TSmin[i][j],TSmax[i][j]);
            }
        }
        //Printing the min and max TD values
        printf("===> TDmin = %lf <===> TDmax = %lf <===\n", TDmin, TDmax);
        //Printing the min and max shear rate values
        //printf("===> qmin = %lf <===> qmax = %lf <===\n", qmin, qmax);
        //Printing the min and max viscosity values
        printf("===> etamin = %lf <===> etamax = %lf <===\n", etamin, etamax);
        //Printing the min and max structural parameter values
        printf("===> spmin = %lf <===> spmax = %lf <===\n", spmin, spmax);
        //Printing the min and max RHS values
        //printf("===> RHSmin = %lf <===> RHSmax = %lf <===\n", RHSmin, RHSmax);
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed structural parameter property
        dp_sync(ns->ed.vevv.dpStructPar);
        // Sync the ditributed viscosity property
        dp_sync(ns->ed.vevv.dpvisc);

    }
    
}

// Calculate RHS for the implicit methid = Gamma*(Re/(1-beta))*TD 
real hig_flow_structural_par_BMP_model_implicit_RHS(real Gamma, real Re, real De, real beta, real TD) {
    real RHS;
    RHS = Gamma*(Re/(1.0-beta))*TD;
    return RHS;
}
