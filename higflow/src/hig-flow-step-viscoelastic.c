// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Step Viscoelastic - version 10/11/2016
// *******************************************************************
// *******************************************************************

#include "hig-flow-step-viscoelastic.h"

// *******************************************************************
// Constitutive Equations
// *******************************************************************

// Computing the Kernel Tensor
// void higflow_compute_kernel_tensor(higflow_solver *ns) {
//     if (ns->contr.flowtype == VISCOELASTIC){
//         // Get the cosntants
//         real Re   = ns->par.Re;
//         real De   = ns->ed.ve.par.De;
//         real beta = ns->ed.ve.par.beta;
//         real tol  = ns->ed.ve.par.kernel_tol;
//         //tol       = 1.0e-3;
//         // Get the local sub-domain for the cells
//         sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
//         // Get the map for the domain properties
//         mp_mapper *mp = sd_get_domain_mapper(sdp);
//         // Loop for each cell
//         higcit_celliterator *it;
//         int num_neg_lambda = 0;
//         for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
//             // Get the cell
//             hig_cell *c = higcit_getcell(it);
//             // Get the cell identifier
//             int clid    = mp_lookup(mp, hig_get_cid(c));
//             // Get the center of the cell
//             Point ccenter;
//             hig_get_center(c, ccenter);
//             // Get the velocity derivative tensor Du and the S tensor
//             real Du[DIM][DIM], S[DIM][DIM];
//             for (int i = 0; i < DIM; i++) {
//                 for (int j = 0; j < DIM; j++) {
//                     // Get Du
//                     Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
//                     // Get S
//                     S[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
//                 }
//             }

//             // Calculate the tensor A
//             real A[DIM][DIM], D[DIM][DIM];

//             switch(ns->ed.ve.contr.model) {
//                 case 4: ;////////////////////////////////////////////////// FENE-P
//                     real Rhs[DIM][DIM];
//                     real L2 = ns->ed.ve.par.L2_fene;
//                     //real b_fene = L2;
//                     real a = L2/(L2-3);
//                     real trRhs = 0;
//                     for (int i = 0; i < DIM; i++) {
//                         for (int j = 0; j < DIM; j++) {
//                             D[i][j] = 0.5*(Du[i][j]+Du[j][i]);
//                             Rhs[i][j]=0.0;
//                             Rhs[i][j] = Re*De*S[i][j]/(1.0-beta) + 2.0*De*D[i][j];    //FENE-P
//                         }
//                         Rhs[i][i] += a*1.0;
//                         trRhs += Rhs[i][i];
//                     }
                    
//                     real cRhs = L2/(L2+trRhs);
//                     for (int i = 0; i < DIM; i++) {
//                         for (int j = 0; j < DIM; j++) {
//                             A[i][j] = cRhs*Rhs[i][j];    //FENE-P
//                         }
//                     }
//                     break;
//                 case 5: ;////////////////////////////////////////////////// e-FENE
//                     real Du_prev[DIM][DIM];
//                     // for(int i = 0; i < DIM; i++) {
//                     //         for(int j = 0; j < DIM; j++) {
//                     //             Du_prev[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpD_prev[i][j], ns->ed.stn);
//                     //     }
//                     // }

//                     // if(ns->par.step == 0 && ns->par.t == 0.0){ //initial conformation tensor
//                     //     real Rhs[DIM][DIM], D_prev[DIM][DIM];
//                     //     for (int i = 0; i < DIM; i++) {
//                     //         for (int j = 0; j < DIM; j++) {
//                     //             D_prev[i][j] = 0.5*(Du_prev[i][j]+Du_prev[j][i]);  
//                     //             Rhs[i][j] = Re*De*S[i][j]/(1.0-beta) + 2.0*De*D_prev[i][j]; //previous adimensional t_p
//                     //         }
//                     //         Rhs[i][i] += 1.0;
//                     //     }
//                     //     hig_flow_compute_initial_conformation_e_fene(Rhs,A,ns->ed.ve.par.L2_fene,ns->ed.ve.par.lambda_fene,ns->ed.ve.par.E_fene);
//                     // }
//                     // else{
//                         /// obtain A using kernel, then use the difference of Du to get the new A
//                         real Kernel[DIM][DIM];
//                         for (int i = 0; i < DIM; i++) {
//                             for (int j = 0; j < DIM; j++) {
//                                 Kernel[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpKernel[i][j], ns->ed.stn);
//                             }
//                         }
//                         // Eige-values and eige-vectors of Kernel
//                         real R[DIM][DIM], lambda[DIM];
//                         hig_flow_jacobi(Kernel, lambda, R);
//                         // Calculate the Inverse Kernel tansformation matrix
//                         real B[DIM][DIM];
//                         for (int i = 0; i < DIM; i++) {
//                             for (int j = i+1; j < DIM; j++) {
//                                 B[i][j] = 0.0;
//                                 B[j][i] = 0.0;
//                             }
//                             B[i][i] = ns->ed.ve.get_kernel_inverse(i, lambda[i], tol);
//                         }
//                         // Calculate A matrix >> A = R B R^t
//                         // D = 0.5*(Du + Du^t)
//                         real A[DIM][DIM];
//                         hig_flow_matrix_transpose_product(B, R, A);
//                         // new A
//                         real D_dif[DIM][DIM];
//                         for (int i = 0; i < DIM; i++) {
//                             for (int j = 0; j < DIM; j++) {
//                                 //difference between D and previous D
//                                 D_dif[i][j] = 0.5*(Du[i][j]+Du[j][i])-0.5*(Du_prev[i][j]+Du_prev[j][i]);
//                                 A[i][j] = A[i][j] + 2.0*De*D_dif[i][j]; 
//                             }
//                         }
//                     // }
//                     break;
//                 default: ///////////////////////////////////////////// Outros
//                     for (int i = 0; i < DIM; i++) {
//                         for (int j = 0; j < DIM; j++) {
//                             D[i][j] = 0.5*(Du[i][j]+Du[j][i]);
//                             A[i][j]=0.0;
//                             A[i][j] = Re*De*S[i][j]/(1.0-beta) + 2.0*De*D[i][j]; //Alterar para FENE-P
//                         }
//                         A[i][i] += 1.0;
//                     }
//                     //printf("A = (%lf %lf)\n    (%lf %lf)\n",A[0][0],A[0][1],A[1][0],A[1][1]);
//                     break;
//             }

//             // Eige-values and eige-vectors of A
//             real R[DIM][DIM], lambda[DIM];
//             hig_flow_jacobi(A, lambda, R);
//             for (int dim = 0; dim < DIM; dim++)
//                 if (lambda[dim] < 0.0) num_neg_lambda++;
//             // Calculate the Kernel tansformation matrix
//             real Kernel[DIM][DIM];
//             hig_flow_calculate_kernel(ns, lambda, R, Kernel, tol);
//             // Store the Kernel Tensor
//             for (int i = 0; i < DIM; i++) {
//                 for (int j = 0; j < DIM; j++) {
//                    dp_set_value(ns->ed.ve.dpKernel[i][j], clid, Kernel[i][j]);
//                 }
//             }
//         }
//         print0f("===> negative eigenvalues - CKT = %d | ",num_neg_lambda);
//         // Destroy the iterator
//         higcit_destroy(it);
//         // Sync the distributed pressure property
//         for (int i = 0; i < DIM; i++) {
//             for (int j = 0; j < DIM; j++) {
//                 dp_sync(ns->ed.ve.dpKernel[i][j]);
//             }
//         }
//     }
// }

// void higflow_compute_initial_kernel_tensor(higflow_solver *ns) {
//     if (ns->contr.flowtype == VISCOELASTIC){
//         if(ns->par.step == 0 && ns->par.t == 0.0){ //initial conformation tensor
//             real tol  = ns->ed.ve.par.kernel_tol;
//             // Get the local sub-domain for the cells
//             sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
//             // Get the map for the domain properties
//             mp_mapper *mp = sd_get_domain_mapper(sdp);
//             // Loop for each cell
//             higcit_celliterator *it;
//             real Kernel[DIM][DIM];
//             for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
//                 // Get the cell
//                 hig_cell *c = higcit_getcell(it);
//                 // Get the cell identifier
//                 int clid    = mp_lookup(mp, hig_get_cid(c));
//                 for (int i = 0; i < DIM; i++) {
//                     for (int j = 0; j < DIM; j++) {
//                         Kernel[i][j] = 0.0;
//                     }
//                     Kernel[i][i] = ns->ed.ve.get_kernel(i, 1.0, tol);
//                 }
//                 // Store the Kernel Tensor
//                 for (int i = 0; i < DIM; i++) {
//                     for (int j = 0; j < DIM; j++) {
//                     dp_set_value(ns->ed.ve.dpKernel[i][j], clid, Kernel[i][j]);
//                     }
//                 }
//             }
//             // Destroy the iterator
//             higcit_destroy(it);
//             // Sync the distributed pressure property
//             for (int i = 0; i < DIM; i++) {
//                 for (int j = 0; j < DIM; j++) {
//                     dp_sync(ns->ed.ve.dpKernel[i][j]);
//                 }
//             }
//         }
//     }
// }


// Computing the Polymeric Tensor
void higflow_compute_polymeric_tensor(higflow_solver *ns) {
    if (ns->contr.flowtype == VISCOELASTIC) {
        // Get the constants
        real Re   = ns->par.Re;
        real De   = ns->ed.ve.par.De;
        real beta = ns->ed.ve.par.beta;
        real tol  = ns->ed.ve.par.kernel_tol;
        real Tmax[DIM][DIM], Tmin[DIM][DIM];
        for (int i = 0; i < DIM; i++) {
           for (int j = 0; j < DIM; j++) {
               Tmax[i][j] = -1.0e16;
               Tmin[i][j] =  1.0e16;
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
                    // if(ns->ed.ve.contr.model==5) //e-fene
                    //     dp_set_value(ns->ed.ve.dpD_prev[i][j], clid, Du[i][j]);
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
            
            real a, trA, fA;
            switch (ns->ed.ve.contr.model) {
                case FENE_P: ;///////////////////////////////////// FENE-P
                    real L2 = ns->ed.ve.par.L2_fene;
                    trA = 0;
                    for (int i = 0; i < DIM; i++)
                        trA += A[i][i];
                    fA = L2/(L2 - trA);
                    a = L2/(L2-3);
                    break;
                case E_FENE: ;///////////////////////////////////// e-FENE
                    real b_fene = ns->ed.ve.par.L2_fene;
                    real lambda_fene = ns->ed.ve.par.lambda_fene;
                    real E = ns->ed.ve.par.E_fene;
                    trA = 0;
                    for (int i = 0; i < DIM; i++)
                        trA += A[i][i];
                    fA = b_fene/(b_fene-trA) - E*sqrt(b_fene)*exp(-sqrt(trA)/lambda_fene)*(1.0/(trA*lambda_fene)+1/(trA*sqrt(trA)));
                    a  = 1.0;
                default: ///////////////////////////////////// Outros
                    fA = 1.0;
                    a  = 1.0;
                    break;
            } 

            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    D[i][j] = 0.5*(Du[i][j]+Du[j][i]);
                    S[i][j] = (1.0-beta)*(fA*A[i][j]-2.0*De*D[i][j])/(Re*De);
                }
                S[i][i] += -a*(1.0-beta)/(Re*De);
            }

            // Store the Polymeric Tensor
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                   real T = S[i][j] + 2.0*(1-beta)*D[i][j]/Re;
                   if (T > Tmax[i][j]) Tmax[i][j] = T;
                   if (T < Tmin[i][j]) Tmin[i][j] = T;
                   dp_set_value(ns->ed.ve.dpS[i][j], clid, S[i][j]);
                }
            }
        }
        for (int i = 0; i < DIM; i++) {
           for (int j = 0; j <= i; j++) {
               // Printing the min and max tensor
               real tmin_global, tmax_global;
               MPI_Allreduce(&Tmin[i][j], &tmin_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
               MPI_Allreduce(&Tmax[i][j], &tmax_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
               print0f("===> %d %d: Tmin = %15.10lf <===> Tmax = %15.10lf <===\n",i,j,tmin_global,tmax_global);
           }
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the distributed pressure property
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
    if (ns->contr.flowtype == VISCOELASTIC) {
        // Get the cosntants
        real Re          = ns->par.Re;
        real De          = ns->ed.ve.par.De;
        real tol         = ns->ed.ve.par.kernel_tol;
        real small       = EPSMACH;
        switch (ns->ed.ve.contr.model) {
            case GPTT:
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
            real Du[DIM][DIM], S[DIM][DIM], Kernel[DIM][DIM], KernelCopy[DIM][DIM];
            // Get the S tensor trace
            real trS = 0.0;
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get Du
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                    // Get S
                    S[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                    // Get Kernel
                    Kernel[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpKernel[i][j], ns->ed.stn);
                    KernelCopy[i][j] = Kernel[i][j];
                }
                trS += S[i][i];
            }

            // Eige-values and eige-vectors of A
            real R[DIM][DIM], Klambda[DIM], lambda[DIM];
            hig_flow_jacobi(KernelCopy, Klambda, R);
            for (int dim = 0; dim < DIM; dim++)
                lambda[dim] = ns->ed.ve.get_kernel_inverse(dim, Klambda[dim], tol);

            // Calculate M matrix >> M = R^t Du R
            real M[DIM][DIM];
            hig_flow_matrix_product(Du, R, M);
            // Calculate Omega matrix >> Omega = R Omega_aux R^t
            real Omega[DIM][DIM];
            hig_flow_calculate_omega(lambda, R, M, Omega, small);
            // Calculate the Kernel jacobian
            real jlambda[DIM];
            for(int i = 0; i < DIM; i++)
                jlambda[i] = ns->ed.ve.get_kernel_jacobian(i, lambda[i], tol);
            // Calculate the matrix BB and the matrix B
            real BB[DIM][DIM]; real B[DIM][DIM];
            hig_flow_calculate_bs (lambda, jlambda, R, M, BB, B);
            // Calculate the matrix MM for the model
            real MM[DIM][DIM], M_aux[DIM][DIM];
            switch (ns->ed.ve.contr.model) {
                case USERSET: 
                    // User Model
                    ns->ed.ve.calculate_m_user(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.ve.par));
                    break;
                case OLDROYD_B: 
                    // Oldroyd-B Model
                    hig_flow_calculate_m_oldroyd(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.ve.par));
                    break;
                case GIESEKUS: 
                    // Giesekus Model
                    hig_flow_calculate_m_giesekus(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.ve.par));
                    break;
                case LPTT: 
                    // LPTT Model
                    hig_flow_calculate_m_lptt(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.ve.par));
                    break;
                case GPTT: 
                    // GPTT Model
                    hig_flow_calculate_m_gptt(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.ve.par));
                    break;
                case FENE_P: 
                    // FENE-P Model
                    hig_flow_calculate_m_fene_p(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.ve.par));
                    break;
                case E_FENE:
                    // e-FENE Model
                    hig_flow_calculate_m_e_fene(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.ve.par));
                    break;
            }
            // Calculate Kernel matrix >> MM = R M(Lambda) JLambda R^t
            hig_flow_matrix_transpose_product(M_aux, R, MM);
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
                                rhs -= hig_flow_convective_cell_term_cubista(ns->dpu[dim], ns->sfdu[dim], ns->stn, ns->ed.ve.dpKernel[i][j], ns->ed.sdED, ns->ed.stn, Kernel[i][j], ccenter, cdelta, dim);
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
        // Sync the distributed pressure property
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.ve.dpS[i][j]);
            }
        }
        // Store the Kernel Tensor
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                    // Get the cell
                    hig_cell *c = higcit_getcell(it);
                    // Get the cell identifier
                    int clid    = mp_lookup(mp, hig_get_cid(c));
                    // Get the center of the cell
                    Point ccenter;
                    hig_get_center(c, ccenter);
                    // Get the S tensor and store in Kernel
                    real S = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);

                    if(j>=i) UPDATE_RESIDUAL_BUFFER_CELL(ns, dp_get_value(ns->ed.ve.dpKernel[i][j], clid), S, c, ccenter)
                
                    // Store Kernel
                    dp_set_value(ns->ed.ve.dpKernel[i][j], clid, S);   
                }
                // Destroy the iterator
                higcit_destroy(it);

                if(j>=i) UPDATE_RESIDUALS(ns, ns->residuals->Kernel[i][j])
            }
        }
        
        // Sync the distributed pressure property
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.ve.dpKernel[i][j]);
            }
        }
    }
}

// // *******************************************************************
// // Calculate convective tensor term CUBISTA
// // *******************************************************************
// real hig_flow_convective_cell_term_cubista(distributed_property *dpu, sim_facet_domain *sfdu, sim_stencil *stn, distributed_property *dpK, sim_domain *sdED, sim_stencil *stnED, real kc, Point ccenter, Point cdelta, int dim) {
//     real  vbar, kr, krr, kl, kll, kc, a, b, c, d, e, fi, conv1,conv2;
//     a     = 1.7500;
//     b     = 0.3750;
//     c     = 0.7500;
//     d     = 0.1250;
//     e     = 0.2500;
//     conv1 = 0.0;
//     conv2 = 0.0;
//     int   incell_r, incell_l, incell_ll, incell_rr, infacet;

//     // Get the low, high, lowlow, highhigh component kernel at center cell
//     kl  = compute_center_p_left_22(sdED, ccenter, cdelta, dim, 1.0, dpK, stnED, &incell_l); 
//     kr  = compute_center_p_right_22(sdED, ccenter, cdelta, dim, 1.0, dpK, stnED, &incell_r); 
//     kll = compute_center_p_left_22(sdED, ccenter, cdelta, dim, 2.0, dpK, stnED, &incell_ll);
//     krr = compute_center_p_right_22(sdED, ccenter, cdelta, dim, 2.0, dpK, stnED, &incell_rr);
//     // Get the velocity  v1bar(i+1/2,j) in the facet center
//     vbar = compute_facet_u_right(sfdu, ccenter, cdelta, dim, 0.5, dpu, stn, &infacet);
//     if (vbar > 0.0){
//         if (FLT_EQ(kr, kl)){
//             conv1 = vbar*kc;
//         }else {
//             fi = (kc - kl)/(kr - kl);
//             if ((fi <= 0.0) || (fi >= 1.0)) {
//                 conv1 = vbar*kc;
//             }else {
//                 if (fi < b){ 
//                     if (incell_l == 1)                    conv1 = vbar*(a*kc - c*kl);
//                     else                                  conv1 = vbar*kc;
//                 }
//            if ((fi >= b) && (fi <= c)){
//                     if ((incell_l == 1)&&(incell_r == 1)) conv1 = vbar*(c*kc + b*kr -d*kl);
//                     else                                  conv1 = vbar*kc;
//                 }
//            if (fi > c){ 
//                     if (incell_r == 1)                    conv1 = vbar*(e*kc + c*kr);
//                     else                                  conv1 = vbar*kc;
//                 }
                    
//             }    
//         }
//     //v1bar < 0.0
//     }else {
//         if ((incell_r == 1) && (incell_rr == 1)){
//             if (FLT_EQ(kc, krr)){
//                 conv1 = vbar*kr;
//             }else {
//                 fi = (kr- krr)/(kc - krr);
//                 if ((fi <= 0.0) || (fi >= 1.0)) {
//                     conv1 = vbar*kr;
//                 }else {
//           if (fi < b) 
//                         conv1 = vbar*(a*kr - c*krr);
//                     if ((fi >= b) && (fi <= c))
//                         conv1 = vbar*(c*kr + b*kc -d*krr);
//                if (fi > c) 
//                         conv1 = vbar*(c*kc + e*kr);
//                 }
//             }
//         //Return upwind value at boundary
//         }else if ((incell_r == 1) && (incell_rr == 0)){
//             if (FLT_EQ(kc, krr)){
//                 conv1 = vbar*kr;
//             }else {
//                 fi = (kr- krr)/(kc - krr);
//                 if ((fi <= 0.0) || (fi >= 1.0)) {
//                     conv1 = vbar*kr;
//                 }else {
//           if (fi <= c) 
//                         conv1 = vbar*kr;
//                if (fi > c) 
//                         conv1 = vbar*(c*kc + e*kr);
//                 }
//             }/*
//             vbar = compute_facet_u_right(sfdu, ccenter, cdelta, dim, 0.5, dpu, stn, &infacet);
//             if (vbar > 0.0) conv1 = vbar*kc;
//             else                 conv1 = vbar*kr;
//             vbar = compute_facet_u_left(sfdu, ccenter, cdelta, dim, 0.5, dpu, stn, &infacet);
//             if (vbar > 0.0) conv2 = vbar*kl;
//             else                 conv2 = vbar*kc;
//             return ((conv1 - conv2)/cdelta[dim]); */
//         }else {
//                 vbar = compute_facet_u_right(sfdu, ccenter, cdelta, dim, 0.5, dpu, stn, &infacet);
//                 if (vbar > 0.0) conv1 = vbar*kc;
//                 else                 conv1 = vbar*kc;
//                 vbar = compute_facet_u_left(sfdu, ccenter, cdelta, dim, 0.5, dpu, stn, &infacet);
//                 if (vbar > 0.0) conv2 = vbar*kl;
//                 else                 conv2 = vbar*kc;
//                 return ((conv1 - conv2)/cdelta[dim]); 
//         }
        
//     }
//     // Get the velocity  v2bar(i-1/2,j) in the facet center
//     vbar = compute_facet_u_left(sfdu, ccenter, cdelta, dim, 0.5, dpu, stn, &infacet);
//     if (vbar > 0.0){
//         if ((incell_l == 1) && (incell_ll == 1)){
//             if (FLT_EQ(kc, kll)) {
//            conv2 = vbar*kl;
//             }else {
//            fi = (kl - kll)/(kc - kll);
//            if ((fi <= 0.0) || (fi >= 1.0)) {
//                conv2 = vbar*kl;
//            }else {
//                if (fi < b)
//                    conv2 = vbar*(a*kl - c*kll);
//                if ((fi >= b) && (fi <= c))
//                    conv2 = vbar*(b*kc + c*kl - d*kll);
//                if (fi > c)  
//                    conv2 = vbar*(c*kc + e*kl);
//            }
//        }
//         }else if ((incell_l == 1) && (incell_ll == 0)){
//             if (FLT_EQ(kc, kll)) {
//            conv2 = vbar*kl;
//             }else {
//            fi = (kl - kll)/(kc - kll);
//            if ((fi <= 0.0) || (fi >= 1.0)) {
//                conv2 = vbar*kl;
//            }else {
//                if (fi <= c)
//                    conv2 = vbar*kl;
//                if (fi > c)  
//                    conv2 = vbar*(c*kc + e*kl);
//            }
//        }/*
//             vbar = compute_facet_u_right(sfdu, ccenter, cdelta, dim, 0.5, dpu, stn, &infacet);
//             if (vbar > 0.0) conv1 = vbar*kc;
//             else                 conv1 = vbar*kr;
//             vbar = compute_facet_u_left(sfdu, ccenter, cdelta, dim, 0.5, dpu, stn, &infacet);
//             if (vbar > 0.0) conv2 = vbar*kl;
//             else                 conv2 = vbar*kc;
//             return ((conv1 - conv2)/cdelta[dim]); */
//        }else {
//                 vbar = compute_facet_u_right(sfdu, ccenter, cdelta, dim, 0.5, dpu, stn, &infacet);
//                 if (vbar > 0.0) conv1 = vbar*kc;
//                 else                 conv1 = vbar*kr;
//                 vbar = compute_facet_u_left(sfdu, ccenter, cdelta, dim, 0.5, dpu, stn, &infacet);
//                 if (vbar > 0.0) conv2 = vbar*kc;
//                 else                 conv2 = vbar*kc;
//                 return ((conv1 - conv2)/cdelta[dim]); 
//         } 
//     }else {
//     //v2bar < 0.0 
//         if (FLT_EQ(kl, kr)) {
//             conv2 = vbar*kc;
//         }else {
//             fi = (kc - kr)/(kl - kr);
//             if ((fi <= 0.0) || (fi >= 1.0)) {
//                 conv2 = vbar*kc;
//             }else {
//            if (fi < b){
//                     if (incell_r == 1)                    conv2 = vbar*(a*kc - c*kr);
//                     else                                  conv2 = vbar*kc;
//                 }
//            if ((fi >= b) && (fi <= c)){
//                     if ((incell_l == 1)&&(incell_r == 1)) conv2 = vbar*(c*kc + b*kl -d*kr);
//                     else                                  conv2 = vbar*kc;
//                 }
//            if (fi > c){ 
//                     if (incell_l == 1)                    conv2 = vbar*(e*kc + c*kl);
//                     else                                  conv2 = vbar*kc;
//                 }
//        }
//         }
//     }
//     return ((conv1-conv2)/cdelta[dim]);
// }

// *******************************************************************
// Constitutive Equation Step for the Implicit Euler Method
// *******************************************************************
void higflow_implicit_euler_constitutive_equation(higflow_solver *ns) {
    if (ns->contr.flowtype == VISCOELASTIC) {
        // Get the cosntants
        real dt    = ns->par.dt;
        real Re          = ns->par.Re;
        real De          = ns->ed.ve.par.De;
        real tol         = ns->ed.ve.par.kernel_tol;
        real small       = EPSMACH;
        switch (ns->ed.ve.contr.model) {
            case GPTT:
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
            real Du[DIM][DIM], S[DIM][DIM], Kernel[DIM][DIM], KernelCopy[DIM][DIM];
            real trS = 0.0;
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get Du
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                    // Get S
                    S[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                    // Get Kernel
                    Kernel[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpKernel[i][j], ns->ed.stn);
                    KernelCopy[i][j] = Kernel[i][j];
                }
                trS += S[i][i];
            }
            
            // Eige-values and eige-vectors of A
            real R[DIM][DIM], Klambda[DIM], lambda[DIM];
            hig_flow_jacobi(KernelCopy, Klambda, R);
            for (int dim = 0; dim < DIM; dim++)
                lambda[dim] = ns->ed.ve.get_kernel_inverse(dim, Klambda[dim], tol);

            // Calculate M matrix >> M = R^t Du R
            real M[DIM][DIM];
            hig_flow_matrix_product(Du, R, M);
            // Calculate Omega matrix >> Omega = R Omega_aux R^t
            real Omega[DIM][DIM];
            hig_flow_calculate_omega(lambda, R, M, Omega, small);
            // Calculate the Kernel jacobian
            real jlambda[DIM];
            for(int i = 0; i < DIM; i++)
                jlambda[i] = ns->ed.ve.get_kernel_jacobian(i, lambda[i], tol);
            // Calculate the matrix BB and the matrix B
            real BB[DIM][DIM]; real B[DIM][DIM];
            hig_flow_calculate_bs (lambda, jlambda, R, M, BB, B);
            // Calculate the matrix MM for the model
            real MM[DIM][DIM], M_aux[DIM][DIM];
            switch (ns->ed.ve.contr.model) {
                case USERSET: 
                    // User Model
                    ns->ed.ve.calculate_m_user(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.ve.par));
                    break;
                case OLDROYD_B: 
                    // Oldroyd-B Model
                    hig_flow_calculate_m_oldroyd(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.ve.par));
                    break;
                case GIESEKUS: 
                    // Giesekus Model
                    hig_flow_calculate_m_giesekus(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.ve.par));
                    break;
                case LPTT: 
                    // LPTT Model
                    hig_flow_calculate_m_lptt(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.ve.par));
                    break;
                case GPTT: 
                    // GPTT Model
                    hig_flow_calculate_m_gptt(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.ve.par));
                    break;
                case FENE_P: 
                    // FENE-P Model
                    hig_flow_calculate_m_fene_p(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.ve.par));
                    break;
                case E_FENE:
                    // e-FENE Model
                    hig_flow_calculate_m_e_fene(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.ve.par));
                    break;
            }
            // Calculate Kernel matrix >> MM = R M(Lambda) JLambda R^t
            hig_flow_matrix_transpose_product(M_aux, R, MM);
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
                    // Right hand side equation
                    real rhs = 0.0;
                    switch (ns->ed.ve.contr.convecdiscrtype) {
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
                                rhs -= hig_flow_convective_cell_term_cubista(ns->dpu[dim], ns->sfdu[dim], ns->stn, ns->ed.ve.dpKernel[i][j], ns->ed.sdED, ns->ed.stn, Kernel[i][j], ccenter, cdelta, dim);
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
        // Sync the distributed pressure property
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.ve.dpS[i][j]);
            }
        }
        // Store the Kernel Tensor
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                    // Get the cell
                    hig_cell *c = higcit_getcell(it);
                    // Get the cell identifier
                    int clid    = mp_lookup(mp, hig_get_cid(c));
                    // Get the center of the cell
                    Point ccenter;
                    hig_get_center(c, ccenter);
                    // Get the S tensor and store in Kernel
                    real S = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);

                    if(j>=i) UPDATE_RESIDUAL_BUFFER_CELL(ns, dp_get_value(ns->ed.ve.dpKernel[i][j], clid), S, c, ccenter)

                    // Store Kernel
                    dp_set_value(ns->ed.ve.dpKernel[i][j], clid, S);   
                }
                // Destroy the iterator
                higcit_destroy(it);

                if(j>=i) UPDATE_RESIDUALS(ns, ns->residuals->Kernel[i][j])
            }
        }
        
        // Sync the distributed pressure property
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.ve.dpKernel[i][j]);
            }
        }
    }
}

// *******************************************************************
// Navier-Stokes Step for the Explicit Euler Method
// *******************************************************************
void higflow_explicit_euler_intermediate_velocity_viscoelastic(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]) {
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
            higflow_computational_cell_viscoelastic(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
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
            rhs += higflow_diffusive_term(ns, fdelta);
            // Compute the intermediate velocity
            real ustar = ns->cc.ufacet + ns->par.dt * rhs;
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
void higflow_explicit_runge_kutta_2_intermediate_velocity_viscoelastic(higflow_solver *ns) {
    // Calculate the auxiliar velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_viscoelastic(ns, ns->dpu, ns->dpuaux);
    // Calculate the star velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_viscoelastic(ns, ns->dpuaux, ns->dpustar);
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
        // Sync the distributed velocity property
        dp_sync(ns->dpustar[dim]);
    }
}

// *******************************************************************
// Navier-Stokes Step for third Order Explicit Runge-Kutta Method
// *******************************************************************
void higflow_explicit_runge_kutta_3_intermediate_velocity_viscoelastic(higflow_solver *ns) {
    // Calculate the auxiliar velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_viscoelastic(ns, ns->dpu, ns->dpuaux);
    // Calculate the second stage velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_viscoelastic(ns, ns->dpuaux, ns->dpustar);
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
        // Sync the distributed velocity property
        dp_sync(ns->dpuaux[dim]);
    }
    // Calculate the order 2 Runge-Kutta method using the euler method
    higflow_explicit_euler_intermediate_velocity_viscoelastic(ns, ns->dpuaux, ns->dpustar);
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
        // Sync the distributed velocity property
        dp_sync(ns->dpustar[dim]);
    }
}

// *******************************************************************
// Navier-Stokes Step for the Implicit Euler Method
// *******************************************************************
void higflow_semi_implicit_euler_intermediate_velocity_viscoelastic(higflow_solver *ns) {
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
            higflow_computational_cell_viscoelastic(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
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
            rhs += ns->cc.ufacet;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
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
void higflow_semi_implicit_crank_nicolson_intermediate_velocity_viscoelastic(higflow_solver *ns) {
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
            higflow_computational_cell_viscoelastic(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Diffusive term term contribution
            rhs += 0.5 * higflow_diffusive_term(ns, fdelta);
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
            rhs += ns->cc.ufacet;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w = - 0.5 * ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
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
void higflow_semi_implicit_bdf2_intermediate_velocity_viscoelastic(higflow_solver *ns) {
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
            higflow_computational_cell_viscoelastic(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
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
            rhs *= 0.5*ns->par.dt;
            // Difusive term contribution
            rhs += 0.25*ns->par.dt*higflow_diffusive_term(ns, fdelta);
            // Velocity term contribution
            rhs += ns->cc.ufacet;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w = - 0.25*ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
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
            higflow_computational_cell_viscoelastic(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            //Get the uaux
            real uaux = dp_get_value(ns->dpuaux[dim], flid);
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
            rhs *= 1.0/3.0*ns->par.dt;
            rhs += (4.0*uaux - ns->cc.ufacet)/3.0;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w = - 1.0/3.0*ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
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
// void hig_flow_calculate_kernel (higflow_solver *ns, real lambda[DIM], real R[DIM][DIM], real Kernel[DIM][DIM], real tol) {
//    // Calculate the Kernel tansformation matrix
//    real Kernel_aux[DIM][DIM];
//    for (int i = 0; i < DIM; i++) {
//        for (int j = i+1; j < DIM; j++) {
//            Kernel_aux[i][j] = 0.0;
//            Kernel_aux[j][i] = 0.0;
//        }
//        Kernel_aux[i][i] = ns->ed.ve.get_kernel(i, lambda[i], tol);
//    }
//    // Calculate Kernel matrix >> Kernel = R Kernel_aux R^t
//    hig_flow_matrix_transpose_product(Kernel_aux, R, Kernel);
// }

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

// Calculate the matrix BB and the matrix B
void hig_flow_calculate_bs (real lambda[DIM], real jlambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real BB[DIM][DIM], real B[DIM][DIM]) {
    // Calculate the matrix BB 
    real B_aux[DIM][DIM];
    for (int i = 0; i < DIM; i++) {
        for (int j = i+1; j < DIM; j++) {
            B_aux[i][j] = 0.0;
            B_aux[j][i] = 0.0;
        }
        B_aux[i][i]  = M[i][i]*lambda[i]*jlambda[i];
    }
    // Calculate Kernel matrix >> BB = R Btilde Lambda JLambda R^t
    hig_flow_matrix_transpose_product(B_aux, R, BB);

    for (int i = 0; i < DIM; i++) {
        for (int j = i+1; j < DIM; j++) {
            B_aux[i][j] = 0.0;
            B_aux[j][i] = 0.0;
        }
        B_aux[i][i]  = M[i][i];
    }
    // Calculate Kernel matrix >> B = R B_aux R^t
    hig_flow_matrix_transpose_product(B_aux, R, B);
}



// Calculate the matrix MM for Oldroyd-B model
void hig_flow_calculate_m_oldroyd (real lambda[DIM], real jlambda[DIM],  real B[DIM][DIM], real M_aux[DIM][DIM], real Re, real trS, ve_parameters *par) {
    // Calculate the matrix MM for Oldroyd-B model
    for (int i = 0; i < DIM; i++) {
        for (int j = i+1; j < DIM; j++) {
            M_aux[i][j] = 0.0;
            M_aux[j][i] = 0.0;
        }
        M_aux[i][i]  = (1.0-lambda[i])*jlambda[i];
    }
}

// Calculate the matrix MM for Giesekus model
void hig_flow_calculate_m_giesekus (real lambda[DIM], real jlambda[DIM],  real B[DIM][DIM], real M_aux[DIM][DIM], real Re, real trS, ve_parameters *par) {
    // Calculate the matrix MM for Giesekus model
    real alpha = par->alpha;

    for (int i = 0; i < DIM; i++) {
        for (int j = i+1; j < DIM; j++) {
            M_aux[i][j] = 0.0;
            M_aux[j][i] = 0.0;
        }
        real aux     = 1.0-lambda[i];
        M_aux[i][i]  = (aux - alpha*aux*aux)*jlambda[i];
    }
}

// Calculate the matrix MM for LPTT model
void hig_flow_calculate_m_lptt (real lambda[DIM], real jlambda[DIM],  real B[DIM][DIM], real M_aux[DIM][DIM], real Re, real trS, ve_parameters *par) {
    // Calculate the matrix MM for LPTT model
    // real B[DIM][DIM], jlambda[DIM];
    // real B_aux[DIM][DIM];
    // for (int i = 0; i < DIM; i++) {
    //     for (int j = i+1; j < DIM; j++) {
    //         B_aux[i][j] = 0.0;
    //         B_aux[j][i] = 0.0;
    //     }
    //     B_aux[i][i]  = M[i][i];
    // }
    // // Calculate Kernel matrix >> B = R B_aux R^t
    // hig_flow_matrix_transpose_product(B_aux, R, B);

    real De      = par->De;
    real beta    = par->beta;
    real epsilon = par->epsilon;
    real xi     = par->xi;
    // Calculate the MM matrix
    for (int i = 0; i < DIM; i++) {
        for (int j = i+1; j < DIM; j++) {
            M_aux[i][j] = 0.0;
            M_aux[j][i] = 0.0;
        }
        M_aux[i][i]  = (1.0-lambda[i])*(1.0+(epsilon*Re*De*trS)/(1.0-beta))*jlambda[i];
    }
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            M_aux[i][j] += -2.0*(B[i][j]-B[i][j]*lambda[j])*De*xi*jlambda[j];
        }
    }
}

// Calculate the matrix MM for GPTT model
void hig_flow_calculate_m_gptt (real lambda[DIM], real jlambda[DIM],  real B[DIM][DIM], real M_aux[DIM][DIM], real Re, real trS, ve_parameters *par) {
    // Calculate the matrix MM for LPTT model
    // real B[DIM][DIM], jlambda[DIM];
    // real B_aux[DIM][DIM];
    // for (int i = 0; i < DIM; i++) {
    //     for (int j = i+1; j < DIM; j++) {
    //         B_aux[i][j] = 0.0;
    //         B_aux[j][i] = 0.0;
    //     }
    //     B_aux[i][i]  = M[i][i];
    // }
    // // Calculate Kernel matrix >> B = R B_aux R^t
    // hig_flow_matrix_transpose_product(B_aux, R, B);

    real De      = par->De;
    real beta    = par->beta;
    real epsilon = par->epsilon;
    real xi     = par->xi;
    real alpha_gptt = par->alpha_gptt;
    real beta_gptt = par->beta_gptt;
    real gamma_gptt = par->gamma_gptt;

    // Calculate the MM matrix
    for (int i = 0; i < DIM; i++) {
        for (int j = i+1; j < DIM; j++) {
            M_aux[i][j] = 0.0;
            M_aux[j][i] = 0.0;
        }
        //real gama1, gama1n; 
        //gama1 = 0.572365;
                //gama1n = lgamma(beta1);
                numc z, mitt;
        z.real = epsilon*Re*De*trS/(1.0-beta);
        z.imag = 0.0;
        mitt   =  mlfv(alpha_gptt, beta_gptt, z, 6);
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
        M_aux[i][i]  = (1.0-lambda[i])*gamma_gptt*mitt.real*jlambda[i];
        // M teste
         //M_aux[i][i]  = (1.0-lambda[i])*mitt.real*jlambda[i];
    }
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            M_aux[i][j] += -2.0*(B[i][j]-B[i][j]*lambda[j])*De*xi*jlambda[j];
        }
    }
}

// Calculate the matrix MM for FENE-P model
void hig_flow_calculate_m_fene_p (real lambda[DIM], real jlambda[DIM],  real B[DIM][DIM], real M_aux[DIM][DIM], real Re, real trS, ve_parameters *par) {

    real L2 = par->L2_fene;
    //real b_fene = L2;
    real a = L2/(L2-3);
    real trA = 0;
    for (int i = 0; i < DIM; i++) {
        trA += lambda[i];
    }
    real fA = L2/(L2 - trA);

    // Calculate the MM matrix
    for (int i = 0; i < DIM; i++) {
        for (int j = i+1; j < DIM; j++) {
            M_aux[i][j] = 0.0;
            M_aux[j][i] = 0.0;
        }
        M_aux[i][i]  = (a - fA*lambda[i])*jlambda[i];
        //M_aux[i][i]  = (1.0-lambda[i])*(1.0+(ns->ed.ve.par.epsilon*ns->par.Re*ns->ed.ve.par.De*tr)/(1.0-ns->ed.ve.par.beta))*jlambda[i];
    }
}

// Calculate the matrix MM for e-FENE model
void hig_flow_calculate_m_e_fene (real lambda[DIM], real jlambda[DIM],  real B[DIM][DIM], real M_aux[DIM][DIM], real Re, real trS, ve_parameters *par) {

    real L2 = par->L2_fene;
    real lambda_fene = par->lambda_fene;
    real E = par->E_fene;

    real b = L2;
    real trA = 0;
    for (int i = 0; i < DIM; i++) {
        trA += lambda[i];
    }
    real fA = b/(b-trA) - E*sqrt(b)*exp(-sqrt(trA)/lambda_fene)*(1.0/(trA*lambda_fene)+1/(trA*sqrt(trA)));

    // Calculate the MM matrix
    for (int i = 0; i < DIM; i++) {
        for (int j = i+1; j < DIM; j++) {
            M_aux[i][j] = 0.0;
            M_aux[j][i] = 0.0;
        }
        M_aux[i][i]  = (1.0 - fA*lambda[i])*jlambda[i];
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

void hig_flow_compute_initial_conformation_e_fene(real Rhs[DIM][DIM], real A[DIM][DIM], real b, real l, real E){
    int maxiter = 20;
    real eps = 1.0e-6;
    real num,den,step;

    //take the trace of T
    real r = 0.0;
    for (int i = 0; i < DIM; i++) {
        r += Rhs[i][i];
    }
    //printf("r = %15.10lf\n", r);

    for (int i = 0; i < DIM; i++) {
        for(int j = 0; j < DIM; j++) {
            if(fabs(Rhs[i][j])<eps) {
                A[i][j] = 0.0;
                continue;
            }
            // newton's method
            real tr = r;
            real tr_next;
            for (int i = 0; i < maxiter; i++) {
                num = 2*l*tr*(tr-b)*(sqrt(b)*E*(sqrt(tr)+l)*(-tr+b) - exp(sqrt(tr)/l)*l*sqrt(tr)*((b+r)*tr-b*r));
                den = 2*b*b*exp(sqrt(tr)/l)*l*l*sqrt(tr)*tr + sqrt(b)*E*(-tr+b)*(-tr+b)*(tr+l*sqrt(tr)+l*l);
                step = num/den;
                //printf("step = %lf\n", step);
                tr_next = tr - step;
                tr = tr_next;
                if(fabs(step) < eps) break;
            }

            real test = (-sqrt(b)*E*exp(-sqrt(tr)/l)*(1.0/(sqrt(tr)*tr)+1.0/(l*tr)) + b/(b-tr))*tr - r;
            //printf("test = %lf\n", test);

            real mul = -sqrt(b)*E*exp(-sqrt(tr)/l)*(1.0/(sqrt(tr)*tr)+1.0/(l*tr)) + b/(b-tr);
            real val = Rhs[i][j]/mul;
            //printf("val = %lf", val);
            A[i][j] = val;
        }
    }
}

// One step of the Navier-Stokes the projection method
void higflow_solver_step_viscoelastic(higflow_solver *ns) {

    // Calculate the velocity derivative tensor
    higflow_compute_velocity_derivative_tensor(ns);
    // Computing the Kernel Tensor
    //higflow_compute_initial_kernel_tensor(ns);
    //higflow_compute_kernel_tensor(ns);

    // Constitutive Equation Step for the Explicit Euler Method
    switch (ns->ed.ve.contr.discrtype) {
        case EXPLICIT:
        // Explicit method
        higflow_explicit_euler_constitutive_equation(ns);
        break;
        case IMPLICIT: 
        // Implicit method
        higflow_implicit_euler_constitutive_equation(ns);
        break;
    }
    // Computing the Polymeric Tensor
    higflow_compute_polymeric_tensor(ns);

    // Boundary condition for velocity
    higflow_boundary_condition_for_velocity(ns);
    // Boundary conditions for source term
    higflow_boundary_condition_for_cell_source_term(ns);
    higflow_boundary_condition_for_facet_source_term(ns);
    // Boundary condition for pressure
    higflow_boundary_condition_for_pressure(ns);
    // // Initial velocity derivative
    // if(ns->par.step == 0 && ns->par.t == 0.0) hig_flow_compute_initial_velocity_derivative_tensor(ns);
    // Calculate the source term
    higflow_calculate_source_term(ns);
    // Calculate the facet source term
    higflow_calculate_facet_source_term(ns);
    // Calculate the intermediated velocity
    switch (ns->contr.tempdiscrtype) {
        case EXPLICIT_EULER:
           // Explicit Euler method
           higflow_explicit_euler_intermediate_velocity_viscoelastic(ns, ns->dpu, ns->dpustar);
           break;
        case EXPLICIT_RK2: 
           // Explicit RK2 method
           higflow_explicit_runge_kutta_2_intermediate_velocity_viscoelastic(ns);
           break;
        case EXPLICIT_RK3: 
           // Explicit RK3 method
           higflow_explicit_runge_kutta_3_intermediate_velocity_viscoelastic(ns);
           break;
        case SEMI_IMPLICIT_EULER: 
           // Semi-Implicit Euler Method
           higflow_semi_implicit_euler_intermediate_velocity_viscoelastic(ns);
           break;
        case SEMI_IMPLICIT_CN: 
           // Semi-Implicit Crank-Nicolson Method
           higflow_semi_implicit_crank_nicolson_intermediate_velocity_viscoelastic(ns);
           break;
        case SEMI_IMPLICIT_BDF2: 
           // Semi-Implicit Crank-Nicolson Method
           higflow_semi_implicit_bdf2_intermediate_velocity_viscoelastic(ns);
           break;
    }
    // Set outflow for ustar velocity 
    //higflow_outflow_ustar_step(ns);
    // Calculate the pressure
    higflow_pressure(ns);
    // Calculate the final velocity
    higflow_final_velocity(ns);
    // Set outflow for velocity
    //higflow_outflow_u_step(ns);
    // Calculate the final pressure
    higflow_final_pressure(ns);
}


// void hig_flow_compute_initial_velocity_derivative_tensor(higflow_solver *ns){
//     int infacet;
//     // Get the local sub-domain for the cells
//     sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
//     // Get the local sub-domain for the facets
//     sim_facet_domain *sfdu[DIM];
//     for(int dim = 0; dim < DIM; dim++) {
//         sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
//     }
//     // Get the map for the domain properties
//     mp_mapper *mp = sd_get_domain_mapper(sdp);
//     // Loop for each cell
//     higcit_celliterator *it;

//     for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
//         // Get the cell
//         hig_cell *c = higcit_getcell(it);
//         // Get the cell identifier
//         int clid    = mp_lookup(mp, hig_get_cid(c));
//         // Get the center of the cell
//         Point ccenter;
//         hig_get_center(c, ccenter);
//         // Get the delta of the cell
//         Point cdelta;
//         hig_get_delta(c, cdelta);
//         // Calculate the velocity derivative tensor
//         for (int dim = 0; dim < DIM; dim++) {
//             for (int dim2 = 0; dim2 < DIM; dim2++) {
//                 real dudx, ul, ur;
//                 if (dim == dim2) {
//                     // Get the velocity in the left facet center
//                     ul   = compute_facet_u_left(sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
//                     // Get the velocity in the right facet center
//                     ur   = compute_facet_u_right(sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
//                 } else {
//                     // Get the velocity in the left facet center
//                     ul = compute_facet_u_4_left(sfdu[dim], ccenter, cdelta, dim, dim2, 1.0, ns->dpu[dim], ns->stn);
//                     // Get the velocity in the right facet center
//                     ur = compute_facet_u_4_right(sfdu[dim], ccenter, cdelta, dim, dim2, 1.0, ns->dpu[dim], ns->stn);
//                 }
//                 dudx = compute_facet_dudxc(cdelta, dim2, 0.5, ul, ul, ur);
//                 dp_set_value(ns->ed.ve.dpD_prev[dim][dim2], clid, dudx);
//             }
//         }
//     }
//     // Destroy the iterator
//     higcit_destroy(it);
//     // Sync the distributed pressure property
//     for (int dim = 0; dim < DIM; dim++) {
//         for (int dim2 = 0; dim2 < DIM; dim2++) {
//             dp_sync(ns->ed.ve.dpD_prev[dim][dim2]);
//         }
//     }
// }
