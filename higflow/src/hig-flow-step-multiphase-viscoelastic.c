// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Step Multiphase - version 10/11/2016
// *******************************************************************
// *******************************************************************
#include "hig-flow-step-multiphase-viscoelastic.h"

char nome_Frac_Visc[100];
char nome_Txx0[100]; char nome_Txx1[100];
char nome_Txy0[100]; char nome_Txy1[100];
char nome_Tyx0[100]; char nome_Tyx1[100];
char nome_Tyy0[100]; char nome_Tyy1[100];

void arquivoTempo_visc(int step) {
    sprintf(nome_Frac_Visc,"DATA/%d_frac_visc.txt", step);
    sprintf(nome_Txx0,"DATA/%d_Txx.txt", step);
    sprintf(nome_Txx1,"DATA/%d_Txx.txt", step);
    sprintf(nome_Txy0,"DATA/%d_Txy.txt", step);
    sprintf(nome_Txy1,"DATA/%d_Txy.txt", step);
    sprintf(nome_Tyx0,"DATA/%d_Tyx.txt", step);
    sprintf(nome_Tyx1,"DATA/%d_Tyx.txt", step);
    sprintf(nome_Tyy0,"DATA/%d_Tyy.txt", step);
    sprintf(nome_Tyy1,"DATA/%d_Tyy.txt", step);

    
    FILE*fp=fopen(nome_Frac_Visc,"a");
    fclose(fp);
    
    fp=fopen(nome_Txx0,"a");
    fclose(fp);
    
    fp=fopen(nome_Txx1,"a");
    fclose(fp);

    fp=fopen(nome_Txy0,"a");
    fclose(fp);
    
    fp=fopen(nome_Txy1,"a");
    fclose(fp);
    
	 fp=fopen(nome_Tyx0,"a");
    fclose(fp);
    
    fp=fopen(nome_Tyx1,"a");
    fclose(fp);
    
	 fp=fopen(nome_Tyy0,"a");
    fclose(fp);
    
	 fp=fopen(nome_Tyy1,"a");
    fclose(fp);
}

void arquivo_Frac_Visc(char*nome,real x,real y,real frac) {
   FILE *fp=fopen(nome,"a");
   fprintf(fp,"%lf %lf %lf\n",x,y,frac);
   fclose(fp);
}

void arquivo_Frac_Ten(char*nome,real x,real y,real Ten) {
   FILE *fp=fopen(nome,"a");
   fprintf(fp,"%lf %lf %lf\n",x,y,Ten);
   fclose(fp);
}

void save_cell_values_visc(higflow_solver *ns,int aux) {
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
        int infacet;
        if(aux==1) {
            //saving volume fraction
            real fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
            arquivo_Frac_Visc(nome_Frac_Visc,ccenter[0],ccenter[1],fracvol);
            
            //saving Tensor
            real S0[DIM][DIM], S1[DIM][DIM];
            // Get Kernel
            S0[1][1] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.dpS0[1][1], ns->ed.stn);
            S0[1][2] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.dpS0[1][2], ns->ed.stn);
            S0[2][1] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.dpS0[2][1], ns->ed.stn);
            S0[2][2] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.dpS0[2][2], ns->ed.stn);
            
            S1[1][1] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.dpS1[1][1], ns->ed.stn);
            S1[1][2] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.dpS1[1][2], ns->ed.stn);
            S1[2][1] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.dpS1[2][1], ns->ed.stn);
            S1[2][2] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.dpS1[2][2], ns->ed.stn);
            
				arquivo_Frac_Ten(nome_Txx0,ccenter[0],ccenter[1],S0[1][1]); 
				arquivo_Frac_Ten(nome_Txx1,ccenter[0],ccenter[1],S1[1][1]); 
            arquivo_Frac_Ten(nome_Txy0,ccenter[0],ccenter[1],S0[1][2]); 
            arquivo_Frac_Ten(nome_Txy1,ccenter[0],ccenter[1],S1[1][2]); 
            arquivo_Frac_Ten(nome_Tyx0,ccenter[0],ccenter[1],S0[2][1]); 
            arquivo_Frac_Ten(nome_Tyx1,ccenter[0],ccenter[1],S1[2][1]); 
            arquivo_Frac_Ten(nome_Tyy0,ccenter[0],ccenter[1],S0[2][2]); 
            arquivo_Frac_Ten(nome_Tyy1,ccenter[0],ccenter[1],S1[2][2]); 
            
            continue;
            }
        }
    higcit_destroy(it);
}

// Computing beta viscoelastic
void higflow_compute_beta_multiphase_viscoelastic(higflow_solver *ns) {
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
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
            // Calculate the density
            real fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
            // Calculate the density
            real beta0 = ns->ed.mult.par.beta0;
            real beta1 = ns->ed.mult.par.beta1;
            real beta  = (1.0-fracvol)*beta0 + fracvol*beta1;
            // Set the viscosity in the distributed viscosity property
            dp_set_value(ns->ed.mult.dpbeta, clid, beta);
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the distributed beta property
        dp_sync(ns->ed.mult.dpbeta);
}

// Computing S viscoelastic
void higflow_compute_S_multiphase_viscoelastic(higflow_solver *ns) {
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
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
            // Calculate the volume fraction
            real fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
            real S0[DIM][DIM], S1[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get Kernel
                    S0[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.dpS0[i][j], ns->ed.stn);
                    S1[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.dpS1[i][j], ns->ed.stn);
                    S0[i][j]  = (1.0 - fracvol)*S0[i][j];
                    S1[i][j]  = fracvol*S1[i][j];
                    // Set tensor
                    dp_set_value(ns->ed.mult.dpS0[i][j], clid, S0[i][j]);
                    dp_set_value(ns->ed.mult.dpS1[i][j], clid, S1[i][j]);
                }
            }
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed S property
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.mult.dpS0[i][j]);
                dp_sync(ns->ed.mult.dpS1[i][j]);
            }
        }
}

// *******************************************************************
// Constitutive Equations
// Computing the Kernel Tensor
void higflow_compute_kernel_tensor_multiphase_viscoelastic(higflow_solver *ns) {
        // Get the cosntants
        real Re   = ns->par.Re;
        real De0   = ns->ed.mult.par.De0;
        real De1   = ns->ed.mult.par.De1;
        real beta0 = ns->ed.mult.par.beta0;
        real beta1 = ns->ed.mult.par.beta1;
        real tol0  = ns->ed.mult.par.kernel_tol0;
        real tol1  = ns->ed.mult.par.kernel_tol1;
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
            // Multiphase
            real frac  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
            real beta2, De2;
            int casovisc = 1;
            if (casovisc == 0) {
               // Ambas as fases são iguais (viscoelasticas)
               beta2 = (1.0 - frac)*beta + frac*beta;
               De2   = (1.0 - frac)*De    + frac*De;
            } else if (casovisc == 1) {
               beta2 = (1.0 - frac)*1.0 + frac*beta;
               De2   = (1.0 - frac)*0.0 + frac*De;
            } else if (casovisc == 2) {
               // Bolha Newtoniana
               beta2 = (1.0 - frac)*beta + frac*1.0;
               De2   = (1.0 - frac)*De   + frac*0.0;
            }
            // Calculate the tensor A
            real A[DIM][DIM], D[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    D[i][j] = 0.5*(Du[i][j]+Du[j][i]);
                    A[i][j]=0.0;
                    if (beta2 <= 0.999){
                       A[i][j] = Re*De2*S[i][j]/(1.0-beta2) + 2.0*De2*D[i][j];
                    }
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

// Constitutive Equation Step for the Explicit Euler Method
void higflow_explicit_euler_constitutive_equation_multiphase_viscoelastic(higflow_solver *ns) {
        // Get the cosntants
        real Re    = ns->par.Re;
        real De    = ns->ed.ve.par.De;
        real beta  = ns->ed.ve.par.beta;
        real tol   = ns->ed.ve.par.kernel_tol;
        real small = 1.0e-14;
        switch (ns->ed.ve.contr.model) {
            case 3:
               ns->ed.ve.par.gamma_gptt = tgamma(ns->ed.ve.par.beta_gptt);
               break;
        }
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
            // Multiphase
            real frac  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
            real beta2, De2;
            int casovisc = 1;
            if (casovisc == 0) {
               // Ambas as fases são iguais (viscoelasticas)
               beta2 = (1.0 - frac)*beta + frac*beta;
               De2   = (1.0 - frac)*De    + frac*De;
            } else if (casovisc == 1) {
               beta2 = (1.0 - frac)*1.0 + frac*beta;
               De2   = (1.0 - frac)*0.0 + frac*De;
            } else if (casovisc == 2) {
               // Bolha Newtoniana
               beta2 = (1.0 - frac)*beta + frac*1.0;
               De2   = (1.0 - frac)*De   + frac*0.0;
            }
            real A[DIM][DIM], D[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    D[i][j] = 0.5*(Du[i][j]+Du[j][i]);
                    A[i][j]=0.0;
                    if (beta2 <= 0.999){
                       A[i][j] = Re*De2*S[i][j]/(1.0-beta2) + 2.0*De2*D[i][j];
                    }
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
            switch (ns->ed.ve.contr.model) {
                case -1: 
                    // User Model
                    ns->ed.ve.calculate_m_user(ns->par.Re, ns->ed.ve.par.De, ns->ed.ve.par.beta, tr, lambda, R, M, M_aux, tol);
                    break;
                case 0: 
                    // Oldroyd Model
                    hig_flow_calculate_m_oldroyd(ns, lambda, M, M_aux, tol);
                    break;
                case 1: 
                    // Giesekus Model
                    hig_flow_calculate_m_giesekus(ns, lambda, M, M_aux, tol);
                    break;
                case 2: 
                    // LPTT Model
                    hig_flow_calculate_m_lptt(ns, tr, lambda, M, R, M_aux, tol);
                    break;
                case 3: 
                    // GPTT Model
                    hig_flow_calculate_m_gptt(ns, tr, lambda, M, R, M_aux, tol);
                    break;
            }
            // Calculate Kernel matrix >> MM = R M(Lambda) JLambda R^t
            hig_flow_matrix_transpose_product(M_aux, R, MM);
            // Calculate RHS = Omega Kernel - Kernel Omega + 2BB + MM/De
            real RHS[DIM][DIM];
            hig_flow_kernel_rhs_multiphase_viscoelastic(De2, Kernel, Omega, BB, MM, RHS);
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
                                rhs -= hig_flow_conv_tensor_term_cub_mult_visc(ns, ns->dpu[dim], ns->ed.sdED, ns->ed.stn, Kernel, ccenter, cdelta, dim, i, j);
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

// *******************************************************************
// Calculate convective tensor term CUBISTA
// *******************************************************************
real hig_flow_conv_tensor_term_cub_mult_visc(higflow_solver *ns, distributed_property *dpu, sim_domain *sdp, sim_stencil *stn, real K[DIM][DIM], Point ccenter, Point cdelta, int dim, int i, int j) {
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

// Constitutive Equation Step for the Implicit Euler Method
void higflow_implicit_euler_constitutive_equation_multiphase_viscoelastic(higflow_solver *ns) {
        // Get the cosntants
        real dt    = ns->par.dt;
        real Re    = ns->par.Re;
        real De    = ns->ed.ve.par.De;
        real beta  = ns->ed.ve.par.beta;
        real tol   = ns->ed.ve.par.kernel_tol;
        real small = 1.0e-14;
        switch (ns->ed.ve.contr.model) {
            case 3:
               ns->ed.ve.par.gamma_gptt = tgamma(ns->ed.ve.par.beta_gptt);
            break;
        }
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
            // Multiphase
            real frac  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
            real beta2, De2;
            int casovisc = 1;
            if (casovisc == 0) {
               // Ambas as fases são iguais (viscoelasticas)
               beta2 = (1.0 - frac)*beta + frac*beta;
               De2   = (1.0 - frac)*De    + frac*De;
            } else if (casovisc == 1) {
               beta2 = (1.0 - frac)*1.0 + frac*beta;
               De2   = (1.0 - frac)*0.0 + frac*De;
            } else if (casovisc == 2) {
               // Bolha Newtoniana
               beta2 = (1.0 - frac)*beta + frac*1.0;
               De2   = (1.0 - frac)*De   + frac*0.0;
            }
            real A[DIM][DIM], D[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    D[i][j] = 0.5*(Du[i][j]+Du[j][i]);
                    A[i][j]=0.0;
                    if (beta2 <= 0.999){
                       A[i][j] = Re*De2*S[i][j]/(1.0-beta2) + 2.0*De2*D[i][j];
                    }
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
            switch (ns->ed.ve.contr.model) {
                case -1: 
                    // User Model
                    ns->ed.ve.calculate_m_user(ns->par.Re, ns->ed.ve.par.De, ns->ed.ve.par.beta, tr, lambda, R, M, M_aux, tol);
                    break;
                case 0: 
                    // Oldroyd Model
                    hig_flow_calculate_m_oldroyd(ns, lambda, M, M_aux, tol);
                    break;
                case 1: 
                    // Giesekus Model
                    hig_flow_calculate_m_giesekus(ns, lambda, M, M_aux, tol);
                    break;
                case 2: 
                    // LPTT Model
                    hig_flow_calculate_m_lptt(ns, tr, lambda, M, R, M_aux, tol);
                    break;
                case 3: 
                    // GPTT Model
                    hig_flow_calculate_m_gptt(ns, tr, lambda, M, R, M_aux, tol);
                    break;
            }
            // Calculate RHS = Omega Kernel - Kernel Omega + 2BB + MM/De
            real RHS[DIM][DIM];
            hig_flow_implicit_kernel_rhs(De2, BB, MM, RHS);//retorna RHS =  2B + M/De
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
                                rhs -= hig_flow_conv_tensor_term_cub_mult_visc(ns, ns->dpu[dim], ns->ed.sdED, ns->ed.stn, Kernel, ccenter, cdelta, dim, i, j);
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

// Computing the Polymeric Tensor
void higflow_compute_polymeric_tensor_multiphase_viscoelastic(higflow_solver *ns) {
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
            // Multiphase
            real frac  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
            real beta2, De2;
            int casovisc = 1;
            if (casovisc == 0) {
               // Ambas as fases são iguais (viscoelasticas)
               beta2 = (1.0 - frac)*beta + frac*beta;
               De2   = (1.0 - frac)*De    + frac*De;
            } else if (casovisc == 1) {
               beta2 = (1.0 - frac)*1.0 + frac*beta;
               De2   = (1.0 - frac)*0.0 + frac*De;
            } else if (casovisc == 2) {
               // Bolha Newtoniana
               beta2 = (1.0 - frac)*beta + frac*1.0;
               De2   = (1.0 - frac)*De   + frac*0.0;
            }
            // Calculate the tensor S
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    D[i][j] = 0.5*(Du[i][j]+Du[j][i]);
                    S[i][j] = 0.0;
                    if (De2 > 0.0001) {
                       S[i][j] = (1.0-beta2)*(A[i][j]-2.0*De2*D[i][j])/(Re*De2);
                    }
                }
                if (De2 > 0.0001) {
                   S[i][i] += -(1.0-beta2)/(Re*De2);
                }
            }
            // Store the Polymeric Tensor
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                   real T = S[i][j] + 2.0*(1-beta2)*D[i][j]/Re;
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

// One step of the Navier-Stokes the projection method
void higflow_solver_step_multiphase_viscoelastic(higflow_solver *ns) {
    // Boundary condition for velocity
    higflow_boundary_condition_for_velocity(ns);
    // Calculate the source term
    higflow_calculate_source_term(ns);
    // Calculate the facet source term
    higflow_calculate_facet_source_term(ns);
    // Calculate beta
    higflow_compute_beta_multiphase_viscoelastic(ns); //OKKK
    // Calculate S 
    higflow_compute_S_multiphase_viscoelastic(ns);//OKKK
    // Calculate the viscosity
    higflow_compute_viscosity_multiphase(ns);
    // Calculate the viscosity
    higflow_compute_density_multiphase(ns);
    // Calculate the curvature
    higflow_compute_curvature_multiphase(ns);
    // Calculate the intermediated velocity
    switch (ns->contr.tempdiscrtype) {
        case 0:
           // Explicit Euler method
           higflow_explicit_euler_intermediate_velocity_multiphase(ns, ns->dpu, ns->dpustar);
           break;
        case 1: 
           // Explicit RK2 method
           higflow_explicit_runge_kutta_2_intermediate_velocity_multiphase(ns);
           break;
        case 2: 
           // Explicit RK3 method
           higflow_explicit_runge_kutta_3_intermediate_velocity_multiphase(ns);
           break;
        case 3: 
           // Semi-Implicit Euler Method
           higflow_semi_implicit_euler_intermediate_velocity_multiphase(ns);
           break;
        case 4: 
           // Semi-Implicit Crank-Nicolson Method
           higflow_semi_implicit_crank_nicolson_intermediate_velocity_multiphase(ns);
           break;
        case 5: 
           // Semi-Implicit Crank-Nicolson Method
           higflow_semi_implicit_bdf2_intermediate_velocity_multiphase(ns, ns->dpu, ns->dpustar);
           break;
    }
    // Boundary condition for pressure
    higflow_boundary_condition_for_pressure(ns);
    // Calculate the pressure
    higflow_pressure_multiphase(ns);
    // Calculate the final velocity
    higflow_final_velocity_multiphase(ns);
    // Boundary condition for velocity
    higflow_boundary_condition_for_velocity(ns);
    // Calculate the final pressure
    higflow_final_pressure(ns);
    /*if (ns->par.stepaux%1000==0) {
        printf("creating archives at step: %d\n",ns->par.stepaux);
        arquivoTempo_visc(ns->par.stepaux);
        save_cell_values_visc(ns,1);
    }*/
    // Calculate the velocity derivative tensor
    higflow_compute_velocity_derivative_tensor(ns);
    // Computing the Kernel Tensor
    higflow_compute_kernel_tensor_multiphase_viscoelastic(ns); //OKKK
    // Constitutive Equation Step for the Explicit Euler Method
    switch (ns->ed.ve.contr.discrtype) {
        case 0:
           // Explicit method
           higflow_explicit_euler_constitutive_equation_multiphase_viscoelastic(ns);//OKKK
           break;
        case 1: 
           // Implicit method
           higflow_implicit_euler_constitutive_equation_multiphase_viscoelastic(ns);//OKK
           break;
    }
    // Computing the Polymeric Tensor
    higflow_compute_polymeric_tensor_multiphase_viscoelastic(ns); //OKKK
    // Calculate the volume fraction
    higflow_plic_advection_volume_fraction(ns); //OKKK
    //higflow_explicit_euler_volume_fraction(ns);
}

// Calculate RHS = OK - KO + 2B * M/De
void hig_flow_kernel_rhs_multiphase_viscoelastic(real De2, real K[DIM][DIM], real O[DIM][DIM], real B[DIM][DIM], real M[DIM][DIM], real RHS[DIM][DIM]) {
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
            RHS[i][j] = 0.0;
            if (De2 > 0.0001) {
                RHS[i][j] = OK[i][j] - KO[i][j] + 2.0*B[i][j] + M[i][j]/De2;
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
void hig_flow_kernel_rhs (real De2, real K[DIM][DIM], real O[DIM][DIM], real B[DIM][DIM], real M[DIM][DIM], real RHS[DIM][DIM]) {
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
            RHS[i][j] = 0.0;
            if (De2 > 0.0001) {
                RHS[i][j] = OK[i][j] - KO[i][j] + 2.0*B[i][j] + M[i][j]/De2;
            }
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

// Calculate the matrix BB
void hig_flow_calculate_b (higflow_solver *ns, real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real BB[DIM][DIM], real tol) {
    // Calculate the matrix BB 
    real B_aux[DIM][DIM];
    for (int i = 0; i < DIM; i++) {
        for (int j = i+1; j < DIM; j++) {
            B_aux[i][j] = 0.0;
            B_aux[j][i] = 0.0;
        }
        real jlambda = ns->ed.ve.get_kernel_jacobian(i, lambda[i], tol);
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
        real jlambda = ns->ed.ve.get_kernel_jacobian(i, lambda[i], tol);
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
        real jlambda = ns->ed.ve.get_kernel_jacobian(i, lambda[i], tol);
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
        jlambda[i]   = ns->ed.ve.get_kernel_jacobian(i, lambda[i], tol);
        M_aux[i][i]  = (1.0-lambda[i])*(1.0+(ns->ed.ve.par.epsilon*ns->par.Re*ns->ed.ve.par.De*tr)/(1.0-ns->ed.ve.par.beta))*jlambda[i];
    }
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            M_aux[i][j] += -2.0*(B[i][j]-B[i][j]*lambda[j])*ns->ed.ve.par.De*ns->ed.ve.par.psi*jlambda[j];
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
        jlambda[i]   = ns->ed.ve.get_kernel_jacobian(i, lambda[i], tol);
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
void hig_flow_implicit_kernel_rhs (real De2, real B[DIM][DIM], real M[DIM][DIM], real RHS[DIM][DIM]) {
    // Calculate RHS = 2B + M/De
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            RHS[i][j] = 0.0;
            if (De2 > 0.0001) {
                RHS[i][j] = 2.0*B[i][j] + M[i][j]/De2;
            }
        }
    }
}

// Calculate the matrix product
void hig_flow_kernel_system_matrix (real w[DIM*DIM][DIM*DIM+1], real Omega[DIM][DIM], real dt) {
    real I[DIM][DIM];
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            I[i][j] = 0.0;
            if (i==j) I[i][j] = 1.0;
        }
    }
    for (int i = 0; i < DIM; i++) {
        for (int j = i; j < DIM; j++) {
            for (int k = 0; k < DIM; k++) {
                for (int l = 0; l < DIM; l++) {
                    if (i==j){
                        w[i*DIM + k][j*DIM + l] = I[k][l] - dt*Omega[k][l] - dt*Omega[i][i]*I[k][l];
                    } else{
                        w[i*DIM + k][j*DIM + l] =  dt*Omega[j][i]*I[k][l];
                        w[j*DIM + l][i*DIM + k] = - dt*Omega[j][i]*I[k][l];
                    }
                }
            }
        }
    }
}
