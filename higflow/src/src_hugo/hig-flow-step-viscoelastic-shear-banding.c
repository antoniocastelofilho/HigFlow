// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Step Electro-osmotic - version 02/02/2018
// *******************************************************************
// *******************************************************************

#include "hig-flow-step-viscoelastic-shear-banding.h"

//********************************************************************
// Viscoelastic functions
//********************************************************************

// Computing the Polymeric Tensor
void higflow_compute_polymeric_tensor_shear_banding(higflow_solver *ns) {
    if (ns->contr.flowtype == 7) {
        // Get the constants
        real Re   = ns->par.Re;
        real De   = ns->ed.vesb.par.De;
        real beta = ns->ed.vesb.par.beta;
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
            real A[DIM][DIM], B[DIM][DIM], Du[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get Du
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpD[i][j], ns->ed.stn);
                    // Get conformation tensor of specie A
                    A[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpA[i][j], ns->ed.stn);
                    // Get conformation tensor of specie B
                    B[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpB[i][j], ns->ed.stn);
                }
            }

            //Get the density number nA
            real nA = compute_value_at_point(ns->ed.vesb.sdSBnA, ccenter, ccenter, 1.0, ns->ed.vesb.dpnA, ns->ed.stn);
            //Get the density number nB
            real nB = compute_value_at_point(ns->ed.vesb.sdSBnB, ccenter, ccenter, 1.0, ns->ed.vesb.dpnB, ns->ed.stn);

            // Calculate the tensor S
            real S[DIM][DIM], D[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    D[i][j] = 0.5*(Du[i][j]+Du[j][i]);
                    S[i][j] = De*(A[i][j]+2.0*B[i][j]-2.0*D[i][j])/Re;
                }
                S[i][i] += -De*(nA+nB)/Re;
            }
            // Store the Polymeric Tensor
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                   real T = S[i][j] + 2.0*De*D[i][j]/Re;
                   if (T > Smax[i][j]) Smax[i][j] = T;
                   if (T < Smin[i][j]) Smin[i][j] = T;
                   dp_set_value(ns->ed.vesb.dpS[i][j], clid, S[i][j]);
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
                dp_sync(ns->ed.vesb.dpS[i][j]);
            }
        }
    }
}

// Computing the breakage rate of specie A (cA)
void higflow_compute_viscoelastic_shear_banding_cA_VCM(higflow_solver *ns) {
    real CAeq   = ns->ed.vesb.par.CAeq;
    real chi   = ns->ed.vesb.par.chi;
    real CAmax = -1.0e16;
    real CAmin = 1.0e16;
    real Dmax[DIM][DIM], Dmin[DIM][DIM];
        for (int i = 0; i < DIM; i++) {
           for (int j = 0; j < DIM; j++) {
               Dmax[i][j] = -1.0e16;
               Dmin[i][j] =  1.0e16;
           }
        }
    if (ns->contr.rheotype == VCM) {
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
            // Get the velocity derivative tensor Du, S and B tensor
            real Du[DIM][DIM], A[DIM][DIM];
            // Get the tensors
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get Du
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpD[i][j], ns->ed.stn);
                    // Get conformation tensor of specie A
                    A[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpA[i][j], ns->ed.stn);
                }
            }

            //Get the density number nA
            real nA = compute_value_at_point(ns->ed.vesb.sdSBnA, ccenter, ccenter, 1.0, ns->ed.vesb.dpnA, ns->ed.stn);

            // Calculate the deformation tensor D
            real D[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    //D[i][j] = 0.5*(Du[i][j]+Du[j][i]);
                    D[i][j] = Du[i][j]+Du[j][i];
                    if (D[i][j]  > Dmax[i][j]) Dmax[i][j] = D[i][j];
                    if (D[i][j]  < Dmin[i][j]) Dmin[i][j] = D[i][j];
                }
            }

            //Calculate tensor ANA = A/nA (Do a loop here)
            real ANA[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    ANA[i][j] = A[i][j]/nA;
                }
            }

            // Calculate the first inner product of the tensors A/nA and D
            real ANAD1[DIM][DIM];
            for (int dim = 0; dim < DIM; dim++) {
                for (int dim2 = 0; dim2 < DIM; dim2++) {
                    ANAD1[dim][dim2] = 0.0;
                    for (int dim3 = 0; dim3 < DIM; dim3++) {
                        //TD1[dim][dim2] += TS[dim][dim3] * D[dim3][dim2];
                        //ANAD1[dim][dim2] += fabs(D[dim][dim3] * ANA[dim3][dim2]);
                        ANAD1[dim][dim2] += D[dim][dim3] * ANA[dim3][dim2];
                    }
                }
            }

            // Double inner product of the tensors A/NA and D
            //real TD = TS[0][1]*D[0][1] + TS[1][0]*D[1][0];
            real ANAD = 0.0;
            for (int dim = 0; dim < DIM; dim++) {
                ANAD += ANAD1[dim][dim];
            }

            // Calculate the breakage rate of specie A
            real cA = ns->ed.vesb.get_cA(ccenter, ns->par.t, CAeq, chi, ANAD);
            //Check to see min and max viscosity values
            if (cA > CAmax) CAmax = cA;
            if (cA < CAmin) CAmin = cA;
            // Set the viscosity in the distributed viscosity property
            dp_set_value(ns->ed.vesb.dpcA, clid, cA);
        }
        //Printing the min and max deformation tensor values
        //for (int i = 0; i < DIM; i++) {
            //for (int j = 0; j < DIM; j++) {
        //        // Printing the min and max tensor
             //   printf("===> %d %d: Dmin = %lf <===> Dmax = %lf <===\n",i,j,Dmin[i][j],Dmax[i][j]);
            //}
        //}
        //Printing the min and max breakage rate values
        printf("===> CAmin = %lf <===> CAmax = %lf <===\n", CAmin, CAmax);
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed pressure property
        dp_sync(ns->ed.vesb.dpcA);
    }
}

// Computing the breakage rate of specie B (cB)
void higflow_compute_viscoelastic_shear_banding_cB_VCM(higflow_solver *ns) {
    real CBeq   = ns->ed.vesb.par.CBeq;
    real chi   = ns->ed.vesb.par.chi;
    real CBMaxMin;
    if (ns->contr.rheotype == VCM) {
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
            // Double inner product of the tensors A/NA and D
            //real TD = TS[0][1]*D[0][1] + TS[1][0]*D[1][0];
            real ANAD = 0.0;
            // Calculate the breakage rate of specie A
            real cB = ns->ed.vesb.get_cB(ccenter, ns->par.t, CBeq, chi, ANAD);
            CBMaxMin = cB;
            // Set the viscosity in the distributed viscosity property
            dp_set_value(ns->ed.vesb.dpcB, clid, cB);
        }
        //Printing the min and max breakage rate values
        printf("===> CBeq = %lf <===> CBeq = %lf <===\n", CBMaxMin, CBMaxMin);
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed pressure property
        dp_sync(ns->ed.vesb.dpcB);
    }
}


// *******************************************************************
// Constitutive Equation Step for the Explicit Euler Method
// *******************************************************************

// Solve the Constitutive Equation of specie A using the Explicit Euler Method
void higflow_explicit_euler_conformation_tensor_A(higflow_solver *ns) {
    if (ns->contr.rheotype == VCM) {
        // Get the cosntants
        real Re    = ns->par.Re;
        real DeA    = ns->ed.vesb.par.DeA;
        real PeA   = ns->ed.vesb.par.PeA;
        real CAeq   = ns->ed.vesb.par.CAeq;
        real CBeq   = ns->ed.vesb.par.CBeq;
        real chi   = ns->ed.vesb.par.chi;
        real small = 1.0e-14;
        real Amax[DIM][DIM], Amin[DIM][DIM];
        for (int i = 0; i < DIM; i++) {
           for (int j = 0; j < DIM; j++) {
               Amax[i][j] = -1.0e16;
               Amin[i][j] =  1.0e16;
           }
        }
        real RHSmax[DIM][DIM], RHSmin[DIM][DIM];
        for (int i = 0; i < DIM; i++) {
           for (int j = 0; j < DIM; j++) {
               RHSmax[i][j] = -1.0e16;
               RHSmin[i][j] =  1.0e16;
           }
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
            // Get the velocity derivative tensor Du, S and B tensor
            real Du[DIM][DIM], S[DIM][DIM], A[DIM][DIM], B[DIM][DIM];
            // Get the tensors
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get Du
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpD[i][j], ns->ed.stn);
                    // Get S
                    S[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpS[i][j], ns->ed.stn);
                    // Get conformation tensor of specie A
                    A[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpA[i][j], ns->ed.stn);
                    // Get conformation tensor of specie B
                    B[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpB[i][j], ns->ed.stn);
                }
            }

            //Get the density number nA
            real nA = compute_value_at_point(ns->ed.vesb.sdSBnA, ccenter, ccenter, 1.0, ns->ed.vesb.dpnA, ns->ed.stn);
            //Get the density number nB
            real nB = compute_value_at_point(ns->ed.vesb.sdSBnB, ccenter, ccenter, 1.0, ns->ed.vesb.dpnB, ns->ed.stn);

            //Get the breakage/reformation rate of specie A
            real CA = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpcA, ns->ed.stn);
            //Get the breakage/reformation rate of specie B
            real CB = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpcB, ns->ed.stn);

            // Calculate the deformation tensor D
            real D[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    D[i][j] = 0.5*(Du[i][j]+Du[j][i]);
                }
            }

            // Eige-values and eige-vectors of A
            //real R[DIM][DIM], lambda[DIM];
            //hig_flow_jacobi(A, lambda, R);
            // Calculate M matrix >> M = R^t Du R
            //real M[DIM][DIM];
            //hig_flow_matrix_product(Du, R, M);
            // Calculate Omega matrix >> Omega = R Omega_aux R^t
            //real Omega[DIM][DIM];
            //hig_flow_calculate_omega(lambda, R, M, Omega, small);
            // Calculate the matrix BB 
            //real BB[DIM][DIM];
            // Calculate the matrix BB 
            //hig_flow_calculate_b (lambda, R, M, BB);
            // Calculate RHS = Omega Kernel - Kernel Omega + 2BB + MM/De
            real RHS[DIM][DIM];
            // Calculate RHS = A*Du + (Du)*T A + (nA*I-A + CB*nB*B -CA*A)/De
            hig_flow_conformation_tensor_A_VCM_rhs(DeA, A, Du, B, CA, CB, nA, nB, RHS);
            // Get the velocity at cell center 
            real u[DIM], dAdx[DIM], d2Adx2[DIM];
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            // Solving the Constitutive Equation using the Euler Method
            for (int i = 0; i < DIM; i++) {
                for (int j = i; j < DIM; j++) {
                    // Right hand side equation
                    real rhs = 0.0;
                    switch (ns->ed.vesb.contr.convecdiscrtype) {
                        case 0: 
                            // Tensor derivative at cell center
                            hig_flow_derivative_tensor_A_at_center_cell (ns, ccenter, cdelta, i, j, A[i][j], dAdx);
                            for (int dim = 0; dim < DIM; dim++) {
                                //Compute convective tensor term in rhs
                                rhs -= u[dim]*dAdx[dim];
                                rhs += higflow_computational_cell_conformation_tensor_A_shear_banding_VCM_model(ns, ns->ed.sdED, ns->ed.stn, A, ccenter, cdelta, dim, i, j);
                            }
                            break;
                        case 1: 
                            //Compute convective tensor term CUBISTA in rhs
                            for (int dim = 0; dim < DIM; dim++) {
                                rhs -= hig_flow_convective_tensor_A_term_cubista(ns, ns->dpu[dim], ns->ed.sdED, ns->ed.stn, A, ccenter, cdelta, dim, i, j);
                                rhs += higflow_computational_cell_conformation_tensor_A_shear_banding_VCM_model(ns, ns->ed.sdED, ns->ed.stn, A, ccenter, cdelta, dim, i, j);
                            }
                            break;
                    }
                    // Compute the diffusive term rhs
                    //for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        //higflow_computational_cell_conformation_tensor_shear_banding_VCM_model(ns, sdp, clid, ccenter, cdelta, dim, i, j, ns->ed.vesb.dpA[i][j]);
                        // Compute the diffusive term rhs
                        //rhs    += hig_flow_second_derivative_tensor_at_center_cell (ns, ccenter, cdelta, dim, i, j, A[i][j], d2Adx2[i][j])/PeA;
                    //}
                    //rhs    += higflow_diffusive_shear_banding_conformation_tensor_A_term(ns);
                    //if (rhs > RHSmax[i][j]) RHSmax[i][j] = rhs;
                    //if (rhs < RHSmin[i][j]) RHSmin[i][j] = rhs;
                    // Compute the final rhs
                    rhs         += RHS[i][j];
                    //Calculate RHS max and min
                    // Compute A at next time
                    real Anew  = A[i][j] + ns->par.dt * rhs;
                    // Store A in S
                    dp_set_value(ns->ed.vesb.dpA[i][j], clid, Anew);
                    if (i != j) {
                         dp_set_value(ns->ed.vesb.dpA[j][i], clid, Anew);
                    }
                    //Calculate max and min A
                    if (Anew > Amax[i][j]) Amax[i][j] = Anew;
                    if (Anew < Amin[i][j]) Amin[i][j] = Anew;
                }
            }
        }

        // Destroy the iterator
        //higcit_destroy(it);
        // Sync the ditributed pressure property
        //for (int i = 0; i < DIM; i++) {
        //    for (int j = 0; j < DIM; j++) {
        //        dp_sync(ns->ed.vesb.dpS[i][j]);
        //    }
        //}

        // Store the Kernel Tensor
        //for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
        //    hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
        //    int clid    = mp_lookup(mp, hig_get_cid(c));
            // Get the center of the cell
        //    Point ccenter;
        //    hig_get_center(c, ccenter);
            // Get the delta of the cell
        //    Point cdelta;
        //    hig_get_delta(c, cdelta);
            // Get the S tensor and store in Kernel
        //    real S[DIM][DIM];
        //    for (int i = 0; i < DIM; i++) {
        //        for (int j = 0; j < DIM; j++) {
                    // Get S
        //            S[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpS[i][j], ns->ed.stn);
                    // Store Kernel
        //            dp_set_value(ns->ed.vesb.dpA[i][j], clid, S[i][j]);
        //        }
        //    }
        //}

        //for (int i = 0; i < DIM; i++) {
        //   for (int j = 0; j < DIM; j++) {
                //Printing the min and max tensor
        //       printf("===> %d %d: Amin = %lf <===> Amax = %lf <===\n",i,j, Amin[i][j], Amax[i][j]);
        //   }
        //}

        //Printing the min and max deformation tensor values
        //for (int i = 0; i < DIM; i++) {
            //for (int j = 0; j < DIM; j++) {
                // Printing the min and max tensor
                //printf("===> %d %d: Dmin = %lf <===> Dmax = %lf <===\n",i,j,Dmin[i][j],Dmax[i][j]);
        //    }
        //}

        //Printing the min and max deformation tensor values
        //for (int i = 0; i < DIM; i++) {
            //for (int j = 0; j < DIM; j++) {
                // Printing the min and max tensor
                //printf("===> %d %d: ANAmin = %lf <===> ANAmax = %lf <===\n",i,j,ANAmin[i][j],ANAmax[i][j]);
            //}
        //}

        //for (int i = 0; i < DIM; i++) {
           //for (int j = 0; j < DIM; j++) {
               // Printing the min and max tensor
               //printf("===> %d %d: RHSAmin = %lf <===> RHSAmax = %lf <===\n",i,j, RHSmin[i][j], RHSmax[i][j]);
           //}
        //}

        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed pressure property
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.vesb.dpA[i][j]);
            }
        }
    }
}

// Solve the Constitutive Equation of specie B using the Explicit Euler Method
void higflow_explicit_euler_conformation_tensor_B(higflow_solver *ns) {
    if (ns->contr.rheotype == VCM) {
        // Get the cosntants
        real Re    = ns->par.Re;
        real DeA    = ns->ed.vesb.par.DeA;
        real epsilon   = ns->ed.vesb.par.epsilon;
        real PeB   = ns->ed.vesb.par.PeB;
        real CAeq   = ns->ed.vesb.par.CAeq;
        real CBeq   = ns->ed.vesb.par.CBeq;
        real chi   = ns->ed.vesb.par.chi;
        real small = 1.0e-14;
        real Bmax[DIM][DIM], Bmin[DIM][DIM];
        for (int i = 0; i < DIM; i++) {
           for (int j = 0; j < DIM; j++) {
               Bmax[i][j] = -1.0e16;
               Bmin[i][j] =  1.0e16;
           }
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
            // Get the velocity derivative tensor Du, S and B tensor
            real Du[DIM][DIM], S[DIM][DIM], A[DIM][DIM], B[DIM][DIM];
            // Get the tensors
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get Du
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpD[i][j], ns->ed.stn);
                    // Get S
                    S[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpS[i][j], ns->ed.stn);
                    // Get conformation tensor of specie A
                    A[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpA[i][j], ns->ed.stn);
                    // Get conformation tensor of specie B
                    B[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpB[i][j], ns->ed.stn);
                }
            }

            //Get the density number nA
            real nA = compute_value_at_point(ns->ed.vesb.sdSBnA, ccenter, ccenter, 1.0, ns->ed.vesb.dpnA, ns->ed.stn);
            //Get the density number nB
            real nB = compute_value_at_point(ns->ed.vesb.sdSBnB, ccenter, ccenter, 1.0, ns->ed.vesb.dpnB, ns->ed.stn);

            //Get the breakage/reformation rate of specie A
            real CA = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpcA, ns->ed.stn);
            //Get the breakage/reformation rate of specie B
            real CB = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpcB, ns->ed.stn);

            // Calculate the rate of deformation tensor
            real D[DIM][DIM];
            for (int dim = 0; dim < DIM; dim++) {
                for (int dim2 = 0; dim2 < DIM; dim2++) {
                    // Calculate the tensor
                    D[dim][dim2] = 0.5*(Du[dim][dim2] + Du[dim2][dim]);
                }
            }

            //Calculate DeB = epsilon*DeA
            real DeB = epsilon*DeA;
            // Eige-values and eige-vectors of B
            //real R[DIM][DIM], lambda[DIM];
            //hig_flow_jacobi(B, lambda, R);
            // Calculate M matrix >> M = R^t Du R
            //real M[DIM][DIM];
            //hig_flow_matrix_product(Du, R, M);
            // Calculate Omega matrix >> Omega = R Omega_aux R^t
            //real Omega[DIM][DIM];
            //hig_flow_calculate_omega(lambda, R, M, Omega, small);
            // Calculate the matrix BB 
            //real BB[DIM][DIM];
            // Calculate the matrix BB 
            //hig_flow_calculate_b (lambda, R, M, BB);
            // Calculate RHS = A*Du + (Du)*T A + (nA*I-A + CB*nB*B -CA*A)/De
            //hig_flow_conformation_tensor_A_VCM_rhs(DeA, A, Omega, BB, B, CA, CB, nA, nB, RHS);
            // Calculate RHS = Omega Kernel - Kernel Omega + 2BB + MM/De
            real RHS[DIM][DIM];
            // Calculate RHS = A*Du + (Du)*T A + (nA*I-A + CB*nB*B -CA*A)/De
            hig_flow_conformation_tensor_B_VCM_rhs (DeA, DeB, B, Du, A, CA, CB, nA, nB, RHS);
            // Get the velocity at cell center 
            real u[DIM], dBdx[DIM], d2Bdx2[DIM][DIM];
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            // Solving the Constitutive Equation using the Euler Method
            for (int i = 0; i < DIM; i++) {
                for (int j = i; j < DIM; j++) {
                    // Right hand side equation
                    real rhs = 0.0;
                    switch (ns->ed.vesb.contr.convecdiscrtype) {
                        case 0: 
                            // Tensor derivative at cell center
                            hig_flow_derivative_tensor_B_at_center_cell (ns, ccenter, cdelta, i, j, B[i][j], dBdx);
                            for (int dim = 0; dim < DIM; dim++) {
                                //Compute convective tensor term in rhs
                                rhs -= u[dim]*dBdx[dim];
                                rhs += higflow_computational_cell_conformation_tensor_B_shear_banding_VCM_model(ns, ns->ed.sdED, ns->ed.stn, B, ccenter, cdelta, dim, i, j);
                            }
                            break;
                        case 1: 
                            //Compute convective tensor term CUBISTA in rhs
                            for (int dim = 0; dim < DIM; dim++) {
                                rhs -= hig_flow_convective_tensor_B_term_cubista(ns, ns->dpu[dim], ns->ed.sdED, ns->ed.stn, B, ccenter, cdelta, dim, i, j);
                                rhs += higflow_computational_cell_conformation_tensor_B_shear_banding_VCM_model(ns, ns->ed.sdED, ns->ed.stn, B, ccenter, cdelta, dim, i, j);
                            }
                            break;
                    }
                    //Compute the diffusive term rhs
                    //for (int dim = 0; dim < DIM; dim++) {
                        // Compute the diffusive term rhs
                        //rhs    += hig_flow_second_derivative_tensor_at_center_cell (ns, ccenter, cdelta, dim, i, j, B[i][j], d2Bdx2[i][j])/PeB;
                    //}
                    // Compute the final rhs
                    rhs         += RHS[i][j];
                    // Compute B at next time
                    real Bnew  = B[i][j] + ns->par.dt * rhs;
                    // Store A in S
                    dp_set_value(ns->ed.vesb.dpB[i][j], clid, Bnew);
                    if (i != j) {
                         dp_set_value(ns->ed.vesb.dpB[j][i], clid, Bnew);
                    }
                    //Calculate max and min B
                    if (Bnew > Bmax[i][j]) Bmax[i][j] = Bnew;
                    if (Bnew < Bmin[i][j]) Bmin[i][j] = Bnew;
                }
            }
        }

        // Destroy the iterator
        //higcit_destroy(it);
        // Sync the ditributed pressure property
        //for (int i = 0; i < DIM; i++) {
        //    for (int j = 0; j < DIM; j++) {
        //        dp_sync(ns->ed.vesb.dpS[i][j]);
        //    }
        //}

        // Store the Kernel Tensor
        //for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
        //    hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
        //    int clid    = mp_lookup(mp, hig_get_cid(c));
            // Get the center of the cell
        //    Point ccenter;
        //    hig_get_center(c, ccenter);
            // Get the delta of the cell
        //    Point cdelta;
        //    hig_get_delta(c, cdelta);
            // Get the S tensor and store in Kernel
        //    real S[DIM][DIM];
        //    for (int i = 0; i < DIM; i++) {
        //        for (int j = 0; j < DIM; j++) {
                    // Get S
        //            S[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpS[i][j], ns->ed.stn);
                    // Store Kernel
        //            dp_set_value(ns->ed.vesb.dpB[i][j], clid, S[i][j]);
        //        }
        //    }
        //}
        
        //for (int i = 0; i < DIM; i++) {
        //   for (int j = 0; j < DIM; j++) {
               // Printing the min and max tensor
        //       printf("===> %d %d: Bmin = %lf <===> Bmax = %lf <===\n",i,j, Bmin[i][j], Bmax[i][j]);
        //   }
        //}

        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed pressure property
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.vesb.dpB[i][j]);
            }
        }

    }
}


// Solve the Constitutive Equation of specie A using the Implicit Euler Method
void higflow_implicit_euler_conformation_tensor_A(higflow_solver *ns) {
    if (ns->contr.rheotype == VCM) {
        // Get the cosntants
        real dt    = ns->par.dt;
        real Re    = ns->par.Re;
        real DeA    = ns->ed.vesb.par.DeA;
        real PeA   = ns->ed.vesb.par.PeA;
        real CAeq   = ns->ed.vesb.par.CAeq;
        real CBeq   = ns->ed.vesb.par.CBeq;
        real chi   = ns->ed.vesb.par.chi;
        real small = 1.0e-14;
        real Amax[DIM][DIM], Amin[DIM][DIM];
        for (int i = 0; i < DIM; i++) {
           for (int j = 0; j < DIM; j++) {
               Amax[i][j] = -1.0e16;
               Amin[i][j] =  1.0e16;
           }
        }
        real RHSmax[DIM][DIM], RHSmin[DIM][DIM];
        for (int i = 0; i < DIM; i++) {
           for (int j = 0; j < DIM; j++) {
               RHSmax[i][j] = -1.0e16;
               RHSmin[i][j] =  1.0e16;
           }
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
            // Get the velocity derivative tensor Du, S and B tensor
            real Du[DIM][DIM], S[DIM][DIM], A[DIM][DIM], B[DIM][DIM];
            // Get the tensors
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get Du
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpD[i][j], ns->ed.stn);
                    // Get S
                    S[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpS[i][j], ns->ed.stn);
                    // Get conformation tensor of specie A
                    A[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpA[i][j], ns->ed.stn);
                    // Get conformation tensor of specie B
                    B[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpB[i][j], ns->ed.stn);
                }
            }

            //Get the density number nA
            real nA = compute_value_at_point(ns->ed.vesb.sdSBnA, ccenter, ccenter, 1.0, ns->ed.vesb.dpnA, ns->ed.stn);
            //Get the density number nB
            real nB = compute_value_at_point(ns->ed.vesb.sdSBnB, ccenter, ccenter, 1.0, ns->ed.vesb.dpnB, ns->ed.stn);

            //Get the breakage/reformation rate of specie A
            real CA = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpcA, ns->ed.stn);
            //Get the breakage/reformation rate of specie B
            real CB = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpcB, ns->ed.stn);

            // Calculate the deformation tensor D and Omega (for the Euler implicit method, omega=DU)
            real D[DIM][DIM];
            real Omega[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    D[i][j] = 0.5*(Du[i][j]+Du[j][i]);
                    //OmegaT[i][j] = Du[j][i];
                    Omega[i][j] = Du[i][j];
                    //Omega[i][j] = Du[j][i];
                }
            }

            // Eige-values and eige-vectors of A
            //real R[DIM][DIM], lambda[DIM];
            //hig_flow_jacobi(A, lambda, R);
            // Calculate M matrix >> M = R^t Du R
            //real M[DIM][DIM];
            //hig_flow_matrix_product(Du, R, M);
            // Calculate Omega matrix >> Omega = R Omega_aux R^t
            //real Omega[DIM][DIM];
            //hig_flow_calculate_omega(lambda, R, M, Omega, small);
            // Calculate the matrix BB 
            //real BB[DIM][DIM];
            // Calculate the matrix BB 
            //hig_flow_calculate_b (lambda, R, M, BB);
            // Calculate RHS = Omega Kernel - Kernel Omega + 2BB + MM/De
            real RHS[DIM][DIM];
            // Calculate RHS =(nA*I-A + CB*nB*B -CA*A)/De
            hig_flow_implicit_conformation_tensor_A_VCM_rhs(DeA, A, B, CA, CB, nA, nB, RHS);
            // Get the velocity at cell center 
            real u[DIM], dAdx[DIM], d2Adx2[DIM];
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            real bS[DIM*DIM];
            real w[DIM*DIM][DIM*DIM+1];
            // Solving the Constitutive Equation using the Euler Method
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Right hand side equation
                    real rhs = 0.0;
                    switch (ns->ed.vesb.contr.convecdiscrtype) {
                        case 0: 
                            // Tensor derivative at cell center
                            hig_flow_derivative_tensor_A_at_center_cell (ns, ccenter, cdelta, i, j, A[i][j], dAdx);
                            for (int dim = 0; dim < DIM; dim++) {
                                //Compute convective tensor term in rhs
                                rhs -= u[dim]*dAdx[dim];
                                rhs += higflow_computational_cell_conformation_tensor_A_shear_banding_VCM_model(ns, ns->ed.sdED, ns->ed.stn, A, ccenter, cdelta, dim, i, j);
                            }
                            break;
                        case 1: 
                            //Compute convective tensor term CUBISTA in rhs
                            for (int dim = 0; dim < DIM; dim++) {
                                rhs -= hig_flow_convective_tensor_A_term_cubista(ns, ns->dpu[dim], ns->ed.sdED, ns->ed.stn, A, ccenter, cdelta, dim, i, j);
                                rhs += higflow_computational_cell_conformation_tensor_A_shear_banding_VCM_model(ns, ns->ed.sdED, ns->ed.stn, A, ccenter, cdelta, dim, i, j);
                            }
                            break;
                    }
                    // Compute the diffusive term rhs
                    //for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        //higflow_computational_cell_conformation_tensor_shear_banding_VCM_model(ns, sdp, clid, ccenter, cdelta, dim, i, j, ns->ed.vesb.dpA[i][j]);
                        // Compute the diffusive term rhs
                        //rhs    += hig_flow_second_derivative_tensor_at_center_cell (ns, ccenter, cdelta, dim, i, j, A[i][j], d2Adx2[i][j])/PeA;
                    //}
                    //rhs    += higflow_diffusive_shear_banding_conformation_tensor_A_term(ns);
                    //if (rhs > RHSmax[i][j]) RHSmax[i][j] = rhs;
                    //if (rhs < RHSmin[i][j]) RHSmin[i][j] = rhs;
                    // Compute the final rhs
                    rhs         += RHS[i][j];
                    rhs         *= dt;
                    rhs         += A[i][j];
                    bS[i*DIM + j] = rhs;
                    //Calculate RHS max and min
                    // Compute A at next time
                    //real Anew  = A[i][j] + ns->par.dt * rhs;
                    // Store A in S
                    //dp_set_value(ns->ed.vesb.dpA[i][j], clid, Anew);
                    //if (i != j) {
                    //     dp_set_value(ns->ed.vesb.dpA[j][i], clid, Anew);
                    //}
                    //Calculate max and min A
                    //if (Anew > Amax[i][j]) Amax[i][j] = Anew;
                    //if (Anew < Amin[i][j]) Amin[i][j] = Anew;
                }
            }
            //Calculate de kronecker product Omega*I - I*Omega
            hig_flow_kernel_system_matrix_VCM (w, Omega,dt);
            //Solve the linear system
            hig_flow_solve_system_constitutive_equation (DIM*DIM, w, bS);
            // Get the solution of linear system
            for (int i = 0; i < DIM; i++) {
                for (int j = i; j < DIM; j++) {
               	    // Get the value of kernel
                    real Anew = bS[i*DIM+j];
                    // Set the value of kernel
                    dp_set_value(ns->ed.vesb.dpA[i][j], clid, Anew);
                    if (i != j) {
                        dp_set_value(ns->ed.vesb.dpA[j][i], clid, Anew);
                    }
                    //Calculate max and min A
                    if (Anew > Amax[i][j]) Amax[i][j] = Anew;
                    if (Anew < Amin[i][j]) Amin[i][j] = Anew;
                }
            } 
        }

        // Destroy the iterator
        //higcit_destroy(it);
        // Sync the ditributed pressure property
        //for (int i = 0; i < DIM; i++) {
        //    for (int j = 0; j < DIM; j++) {
        //        dp_sync(ns->ed.vesb.dpS[i][j]);
        //    }
        //}

        // Store the Kernel Tensor
        //for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
        //    hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
        //    int clid    = mp_lookup(mp, hig_get_cid(c));
            // Get the center of the cell
        //    Point ccenter;
        //    hig_get_center(c, ccenter);
            // Get the delta of the cell
        //    Point cdelta;
        //    hig_get_delta(c, cdelta);
            // Get the S tensor and store in Kernel
        //    real S[DIM][DIM];
        //    for (int i = 0; i < DIM; i++) {
        //        for (int j = 0; j < DIM; j++) {
                    // Get S
        //            S[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpS[i][j], ns->ed.stn);
                    // Store Kernel
        //            dp_set_value(ns->ed.vesb.dpA[i][j], clid, S[i][j]);
        //        }
        //    }
        //}

        for (int i = 0; i < DIM; i++) {
           for (int j = 0; j < DIM; j++) {
                //Printing the min and max tensor
               printf("===> %d %d: Amin = %lf <===> Amax = %lf <===\n",i,j, Amin[i][j], Amax[i][j]);
           }
        }

        //Printing the min and max deformation tensor values
        //for (int i = 0; i < DIM; i++) {
            //for (int j = 0; j < DIM; j++) {
                // Printing the min and max tensor
                //printf("===> %d %d: Dmin = %lf <===> Dmax = %lf <===\n",i,j,Dmin[i][j],Dmax[i][j]);
        //    }
        //}

        //Printing the min and max deformation tensor values
        //for (int i = 0; i < DIM; i++) {
            //for (int j = 0; j < DIM; j++) {
                // Printing the min and max tensor
                //printf("===> %d %d: ANAmin = %lf <===> ANAmax = %lf <===\n",i,j,ANAmin[i][j],ANAmax[i][j]);
            //}
        //}

        //for (int i = 0; i < DIM; i++) {
           //for (int j = 0; j < DIM; j++) {
               // Printing the min and max tensor
               //printf("===> %d %d: RHSAmin = %lf <===> RHSAmax = %lf <===\n",i,j, RHSmin[i][j], RHSmax[i][j]);
           //}
        //}

        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed pressure property
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.vesb.dpA[i][j]);
            }
        }
    }
}

// Solve the Constitutive Equation of specie B using the Implicit Euler Method
void higflow_implicit_euler_conformation_tensor_B(higflow_solver *ns) {
    if (ns->contr.rheotype == VCM) {
        // Get the cosntants
        real dt    = ns->par.dt;
        real Re    = ns->par.Re;
        real DeA    = ns->ed.vesb.par.DeA;
        real epsilon   = ns->ed.vesb.par.epsilon;
        real PeB   = ns->ed.vesb.par.PeB;
        real CAeq   = ns->ed.vesb.par.CAeq;
        real CBeq   = ns->ed.vesb.par.CBeq;
        real chi   = ns->ed.vesb.par.chi;
        real small = 1.0e-14;
        real Bmax[DIM][DIM], Bmin[DIM][DIM];
        for (int i = 0; i < DIM; i++) {
           for (int j = 0; j < DIM; j++) {
               Bmax[i][j] = -1.0e16;
               Bmin[i][j] =  1.0e16;
           }
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
            // Get the velocity derivative tensor Du, S and B tensor
            real Du[DIM][DIM], S[DIM][DIM], A[DIM][DIM], B[DIM][DIM];
            // Get the tensors
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get Du
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpD[i][j], ns->ed.stn);
                    // Get S
                    S[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpS[i][j], ns->ed.stn);
                    // Get conformation tensor of specie A
                    A[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpA[i][j], ns->ed.stn);
                    // Get conformation tensor of specie B
                    B[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpB[i][j], ns->ed.stn);
                }
            }

            //Get the density number nA
            real nA = compute_value_at_point(ns->ed.vesb.sdSBnA, ccenter, ccenter, 1.0, ns->ed.vesb.dpnA, ns->ed.stn);
            //Get the density number nB
            real nB = compute_value_at_point(ns->ed.vesb.sdSBnB, ccenter, ccenter, 1.0, ns->ed.vesb.dpnB, ns->ed.stn);

            //Get the breakage/reformation rate of specie A
            real CA = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpcA, ns->ed.stn);
            //Get the breakage/reformation rate of specie B
            real CB = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpcB, ns->ed.stn);


            // Calculate the deformation tensor D and Omega (for the Euler implicit method, omega=DU T)
            real D[DIM][DIM];
            real Omega[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    D[i][j] = 0.5*(Du[i][j]+Du[j][i]);
                    //OmegaT[i][j] = Du[j][i];
                    Omega[i][j] = Du[i][j];
                    //Omega[i][j] = Du[j][i];
                }
            }


            //Calculate DeB = epsilon*DeA
            real DeB = epsilon*DeA;
            // Eige-values and eige-vectors of B
            //real R[DIM][DIM], lambda[DIM];
            //hig_flow_jacobi(B, lambda, R);
            // Calculate M matrix >> M = R^t Du R
            //real M[DIM][DIM];
            //hig_flow_matrix_product(Du, R, M);
            // Calculate Omega matrix >> Omega = R Omega_aux R^t
            //real Omega[DIM][DIM];
            //hig_flow_calculate_omega(lambda, R, M, Omega, small);
            // Calculate the matrix BB 
            //real BB[DIM][DIM];
            // Calculate the matrix BB 
            //hig_flow_calculate_b (lambda, R, M, BB);
            // Calculate RHS = A*Du + (Du)*T A + (nA*I-A + CB*nB*B -CA*A)/De
            //hig_flow_conformation_tensor_A_VCM_rhs(DeA, A, Omega, BB, B, CA, CB, nA, nB, RHS);
            // Calculate RHS = Omega Kernel - Kernel Omega + 2BB + MM/De
            real RHS[DIM][DIM];
            // Calculate RHS = (nA*I-A + CB*nB*B -CA*A)/De
            hig_flow_implicit_conformation_tensor_B_VCM_rhs(DeA, DeB, B, A, CA, CB, nA, nB, RHS);
            // Get the velocity at cell center 
            real u[DIM], dBdx[DIM], d2Bdx2[DIM][DIM];
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            real bS[DIM*DIM];
            real w[DIM*DIM][DIM*DIM+1];
            // Solving the Constitutive Equation using the Euler Method
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Right hand side equation
                    real rhs = 0.0;
                    switch (ns->ed.vesb.contr.convecdiscrtype) {
                        case 0: 
                            // Tensor derivative at cell center
                            hig_flow_derivative_tensor_B_at_center_cell (ns, ccenter, cdelta, i, j, B[i][j], dBdx);
                            for (int dim = 0; dim < DIM; dim++) {
                                //Compute convective tensor term in rhs
                                rhs -= u[dim]*dBdx[dim];
                                rhs += higflow_computational_cell_conformation_tensor_B_shear_banding_VCM_model(ns, ns->ed.sdED, ns->ed.stn, B, ccenter, cdelta, dim, i, j);
                            }
                            break;
                        case 1: 
                            //Compute convective tensor term CUBISTA in rhs
                            for (int dim = 0; dim < DIM; dim++) {
                                rhs -= hig_flow_convective_tensor_B_term_cubista(ns, ns->dpu[dim], ns->ed.sdED, ns->ed.stn, B, ccenter, cdelta, dim, i, j);
                                rhs += higflow_computational_cell_conformation_tensor_B_shear_banding_VCM_model(ns, ns->ed.sdED, ns->ed.stn, B, ccenter, cdelta, dim, i, j);
                            }
                            break;
                    }
                    //Compute the diffusive term rhs
                    //for (int dim = 0; dim < DIM; dim++) {
                        // Compute the diffusive term rhs
                        //rhs    += hig_flow_second_derivative_tensor_at_center_cell (ns, ccenter, cdelta, dim, i, j, B[i][j], d2Bdx2[i][j])/PeB;
                    //}
                    // Compute the final rhs
                    rhs         += RHS[i][j];
                    rhs         *= dt;
                    rhs         += B[i][j];
                    bS[i*DIM + j] = rhs;
                    // Compute B at next time
                    //real Bnew  = B[i][j] + ns->par.dt * rhs;
                    // Store A in S
                    //dp_set_value(ns->ed.vesb.dpB[i][j], clid, Bnew);
                    //if (i != j) {
                    //     dp_set_value(ns->ed.vesb.dpB[j][i], clid, Bnew);
                    //}
                    //Calculate max and min B
                    //if (Bnew > Bmax[i][j]) Bmax[i][j] = Bnew;
                    //if (Bnew < Bmin[i][j]) Bmin[i][j] = Bnew;
                }
            }
            //Calculate de kronecker product Omega*I - I*Omega
            hig_flow_kernel_system_matrix_VCM (w, Omega,dt);
            //Solve the linear system
            hig_flow_solve_system_constitutive_equation (DIM*DIM, w, bS);
            for (int i = 0; i < DIM; i++) {
                for (int j = i; j < DIM; j++) {
               	    // Get the value of kernel
                    real Bnew = bS[i*DIM+j];
                    // Set the value of B
                    dp_set_value(ns->ed.vesb.dpB[i][j], clid, Bnew);
                    if (i != j) {
                        dp_set_value(ns->ed.vesb.dpB[j][i], clid, Bnew);
                    }
                    //Calculate max and min B
                    if (Bnew > Bmax[i][j]) Bmax[i][j] = Bnew;
                    if (Bnew < Bmin[i][j]) Bmin[i][j] = Bnew;
                }
            }
        }

        // Destroy the iterator
        //higcit_destroy(it);
        // Sync the ditributed pressure property
        //for (int i = 0; i < DIM; i++) {
        //    for (int j = 0; j < DIM; j++) {
        //        dp_sync(ns->ed.vesb.dpS[i][j]);
        //    }
        //}

        // Store the Kernel Tensor
        //for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
        //    hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
        //    int clid    = mp_lookup(mp, hig_get_cid(c));
            // Get the center of the cell
        //    Point ccenter;
        //    hig_get_center(c, ccenter);
            // Get the delta of the cell
        //    Point cdelta;
        //    hig_get_delta(c, cdelta);
            // Get the S tensor and store in Kernel
        //    real S[DIM][DIM];
        //    for (int i = 0; i < DIM; i++) {
        //        for (int j = 0; j < DIM; j++) {
                    // Get S
        //            S[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpS[i][j], ns->ed.stn);
                    // Store Kernel
        //            dp_set_value(ns->ed.vesb.dpB[i][j], clid, S[i][j]);
        //        }
        //    }
        //}
        
        for (int i = 0; i < DIM; i++) {
           for (int j = 0; j < DIM; j++) {
               // Printing the min and max tensor
               printf("===> %d %d: Bmin = %lf <===> Bmax = %lf <===\n",i,j, Bmin[i][j], Bmax[i][j]);
           }
        }

        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed pressure property
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.vesb.dpB[i][j]);
            }
        }

    }
}


// *******************************************************************
// Calculate convective tensor term CUBISTA
// *******************************************************************
real hig_flow_convective_tensor_A_term_cubista(higflow_solver *ns, distributed_property *dpu, sim_domain *sdp, sim_stencil *stn, real K[DIM][DIM], Point ccenter, Point cdelta, int dim, int i, int j) {
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
    // Get the tensor at center cell
    kc  = K[i][j];
    // Get the low, high, lowlow, highhigh component kernel at center cell
    kl  = compute_center_p_left_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.vesb.dpA[i][j], ns->ed.stn, &incell_l); 
    kr  = compute_center_p_right_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.vesb.dpA[i][j], ns->ed.stn, &incell_r); 
    kll = compute_center_p_left_22(ns->ed.sdED, ccenter, cdelta, dim, 2.0, ns->ed.vesb.dpA[i][j], ns->ed.stn, &incell_ll);
    krr = compute_center_p_right_22(ns->ed.sdED, ccenter, cdelta, dim, 2.0, ns->ed.vesb.dpA[i][j], ns->ed.stn, &incell_rr);
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

real hig_flow_convective_tensor_B_term_cubista(higflow_solver *ns, distributed_property *dpu, sim_domain *sdp, sim_stencil *stn, real K[DIM][DIM], Point ccenter, Point cdelta, int dim, int i, int j) {
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
    // Get the tensor at center cell
    kc  = K[i][j];
    // Get the low, high, lowlow, highhigh component kernel at center cell
    kl  = compute_center_p_left_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.vesb.dpB[i][j], ns->ed.stn, &incell_l); 
    kr  = compute_center_p_right_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.vesb.dpB[i][j], ns->ed.stn, &incell_r); 
    kll = compute_center_p_left_22(ns->ed.sdED, ccenter, cdelta, dim, 2.0, ns->ed.vesb.dpB[i][j], ns->ed.stn, &incell_ll);
    krr = compute_center_p_right_22(ns->ed.sdED, ccenter, cdelta, dim, 2.0, ns->ed.vesb.dpB[i][j], ns->ed.stn, &incell_rr);
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



//Computing the terms needed for the simulation of viscoelastic flows that exhibit shear-banding behaviour using the VCM model
real higflow_computational_cell_conformation_tensor_A_shear_banding_VCM_model(higflow_solver *ns, sim_domain *sdp, sim_stencil *stn, real A[DIM][DIM], Point ccenter, Point cdelta, int dim, int i, int j){
    real Ac, Ar, Al;
    //real  ABl[DIM][DIM], ABr[DIM][DIM];
    int  incell_l, incell_r;
    // Get the tensor at the cell 
    Ac  = A[i][j];
	// Get the conformation tensor in the left cell
    //Al = compute_center_p_left_22(sdp, ccenter, cdelta, dim, 1.0, dpA[i][j], ns->ed.stn, &incell_l);
    Al  = compute_center_p_left_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.vesb.dpA[i][j], ns->ed.stn, &incell_l); 
    Ar  = compute_center_p_right_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.vesb.dpA[i][j], ns->ed.stn, &incell_r); 
    // Compute the necessary derivatives
	real d2Adx2   = compute_facet_du2dx2(cdelta, dim, 1.0, Ac, Al, Ar);
    d2Adx2 /= ns->ed.vesb.par.PeA;
    return d2Adx2;
//		break;
//	}
}

//Computing the terms needed for the simulation of viscoelastic flows that exhibit shear-banding behaviour using the VCM model
real higflow_computational_cell_conformation_tensor_B_shear_banding_VCM_model(higflow_solver *ns, sim_domain *sdp, sim_stencil *stn, real B[DIM][DIM], Point ccenter, Point cdelta, int dim, int i, int j){
    real Bc, Br, Bl;
    //real  ABl[DIM][DIM], ABr[DIM][DIM];
    int  incell_l, incell_r;
    // Get the tensor at the cell 
    Bc  = B[i][j];
	// Get the conformation tensor in the left cell
    //Al = compute_center_p_left_22(sdp, ccenter, cdelta, dim, 1.0, dpA[i][j], ns->ed.stn, &incell_l);
    Bl  = compute_center_p_left_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.vesb.dpB[i][j], ns->ed.stn, &incell_l); 
    Br  = compute_center_p_right_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.vesb.dpB[i][j], ns->ed.stn, &incell_r); 
    // Compute the necessary derivatives
	real d2Bdx2   = compute_facet_du2dx2(cdelta, dim, 1.0, Bc, Bl, Br);
    d2Bdx2 /= ns->ed.vesb.par.PeB;
    return d2Bdx2;
//		break;
//	}
}


// *******************************************************************
// Constitutive Equation Step for the Implicit Euler Method
// *******************************************************************

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
void hig_flow_calculate_b (real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real BB[DIM][DIM]) {
    // Calculate the matrix BB 
    real B_aux[DIM][DIM];
    for (int i = 0; i < DIM; i++) {
        for (int j = i+1; j < DIM; j++) {
            B_aux[i][j] = 0.0;
            B_aux[j][i] = 0.0;
        }
        B_aux[i][i]  = M[i][i]*lambda[i];
    }
    // Calculate Kernel matrix >> BB = R Btilde Lambda JLambda R^t
    hig_flow_matrix_transpose_product(B_aux, R, BB);
}

// Calculate RHS = A*Du + (Du)*T A + (nA*I-A + CB*nB*B -CA*A)/De
void hig_flow_conformation_tensor_A_VCM_rhs (real De, real A[DIM][DIM], real DU[DIM][DIM], real B[DIM][DIM], real CA, real CB, real nA, real nB, real RHS[DIM][DIM]) {
    real DUA[DIM][DIM], ADU[DIM][DIM];
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            DUA[i][j] = 0.0;
            ADU[i][j] = 0.0;
            for (int k = 0; k < DIM; k++) {
                ADU[i][j] += DU[i][k]*A[k][j];
                DUA[i][j] += DU[j][k]*A[i][k];
            }
        }
    }
    // Calculate RHS = A*Du + (Du)*T A + (nA*I-A + CB*nB*B -CA*A)/De
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            RHS[i][j] = ADU[i][j] + DUA[i][j] + (CB*nB*B[i][j] -CA*A[i][j] -A[i][j])/De;
        }
        RHS[i][i] += nA/De;
    }
}

// Calculate RHS = B*Du + (Du)*T B + (nB*I/2-B)/De +2*(CA*A-CB*nB*B)/DeB
void hig_flow_conformation_tensor_B_VCM_rhs (real DeA, real DeB, real B[DIM][DIM], real DU[DIM][DIM], real A[DIM][DIM], real CA, real CB, real nA, real nB, real RHS[DIM][DIM]) {
    real DUB[DIM][DIM], BDU[DIM][DIM];
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            DUB[i][j] = 0.0;
            BDU[i][j] = 0.0;
            for (int k = 0; k < DIM; k++) {
                BDU[i][j] += DU[i][k]*B[k][j];
                DUB[i][j] += DU[j][k]*B[i][k];
            }
        }
    }
    // Calculate RHS = B*Du + (Du)*T B + (nB*I/2-B + CA*A - CB*nB*B)/De
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            RHS[i][j] = BDU[i][j] + DUB[i][j] + 2.0*(CA*A[i][j] - CB*nB*B[i][j])/DeA -B[i][j]/DeB;
        }
        RHS[i][i] += (0.5*nB)/DeB;
    }
}

// Calculate RHS = (nA*I-A + CB*nB*B -CA*A)/De
void hig_flow_implicit_conformation_tensor_A_VCM_rhs (real De, real A[DIM][DIM], real B[DIM][DIM], real CA, real CB, real nA, real nB, real RHS[DIM][DIM]) {
    // Calculate RHS = (nA*I-A + CB*nB*B -CA*A)/De
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            RHS[i][j] = (CB*nB*B[i][j] -CA*A[i][j] -A[i][j])/De;
        }
        RHS[i][i] += nA/De;
    }
}

// Calculate RHS = (nB*I/2-B)/De +2*(CA*A-CB*nB*B)/DeB
void hig_flow_implicit_conformation_tensor_B_VCM_rhs (real DeA, real DeB, real B[DIM][DIM], real A[DIM][DIM], real CA, real CB, real nA, real nB, real RHS[DIM][DIM]) {
    // Calculate RHS = B*Du + (Du)*T B + (nB*I/2-B + CA*A - CB*nB*B)/De
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            RHS[i][j] = 2.0*(CA*A[i][j] - CB*nB*B[i][j])/DeA -B[i][j]/DeB;
        }
        RHS[i][i] += (0.5*nB)/DeB;
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

// Get the derivative of conformation tensor A
void hig_flow_derivative_tensor_A_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, int i, int j, real Acenter, real dAdx[DIM]) {
    for (int dim = 0; dim < DIM; dim++) {
        int incell_left, incell_right;
        // Get the Kernel in the left cell
        real Aleft = compute_center_p_left_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.vesb.dpA[i][j], ns->ed.stn, &incell_left);
        // Get the Kernel in the right cell
        real Aright = compute_center_p_right_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.vesb.dpA[i][j], ns->ed.stn, &incell_right);
        // Compute the Kernel derivative
        if ((incell_left == 1) && (incell_right == 1)) { 
           dAdx[dim] = compute_dpdx_at_point(cdelta, dim, 1.0, Aleft, Aright);
        } else if (incell_right == 1) { 
           dAdx[dim] = compute_dpdxr_at_point(cdelta, dim, 1.0, Acenter, Aright);
        } else {
           dAdx[dim] = compute_dpdxl_at_point(cdelta, dim, 1.0, Aleft, Acenter);
        }
    }
}

// Get the derivative of conformation tensor B
void hig_flow_derivative_tensor_B_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, int i, int j, real Bcenter, real dBdx[DIM]) {
    for (int dim = 0; dim < DIM; dim++) {
        int incell_left, incell_right;
        // Get the Kernel in the left cell
        real Bleft = compute_center_p_left_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.vesb.dpB[i][j], ns->ed.stn, &incell_left);
        // Get the Kernel in the right cell
        real Bright = compute_center_p_right_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.vesb.dpB[i][j], ns->ed.stn, &incell_right);
        // Compute the Kernel derivative
        if ((incell_left == 1) && (incell_right == 1)) { 
           dBdx[dim] = compute_dpdx_at_point(cdelta, dim, 1.0, Bleft, Bright);
        } else if (incell_right == 1) { 
           dBdx[dim] = compute_dpdxr_at_point(cdelta, dim, 1.0, Bcenter, Bright);
        } else {
           dBdx[dim] = compute_dpdxl_at_point(cdelta, dim, 1.0, Bleft, Bcenter);
        }
    }
}

// Get the second derivative of conformation tensor A
//void hig_flow_second_derivative_tensor_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, int dim, int i, int j, real Acenter, real d2Adx2[DIM][DIM]) {
    //real ABl[DIM][DIM], ABr[DIM][DIM];
    //int incell_l, incell_r;
    //real ABl = compute_center_p_left_22(sdn, ccenter, cdelta, dim, 1.0, Acenter, ns->ed.stn, &incell_l);
	//compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.vesb.dpS[dim][dim2], ns->ed.stn)
	// Get the conformation tensor in the right cell
	//real ABr = compute_center_p_right_22(sdn, ccenter, cdelta, dim, 1.0, Acenter, ns->ed.stn, &incell_r);
    //d2Adx2[i][j] = compute_facet_du2dx2(cdelta, dim, 1.0, Acenter, ABl, ABr) ;
//}

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


// Calculate the matrix product
void hig_flow_kernel_system_matrix_VCM (real w[DIM*DIM][DIM*DIM+1], real Omega[DIM][DIM], real dt) {
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
		        w[i*DIM + k][j*DIM + l] =  - dt*Omega[i][j]*I[k][l];
                        w[j*DIM + l][i*DIM + k] = - dt*Omega[j][i]*I[k][l];
                    }
                }
            }
        }
    }
}

// Calculate the matrix product
//void hig_flow_kernel_system_matrix_VCM (real w[DIM*DIM][DIM*DIM+1], real Omega[DIM][DIM], real dt) {
//    real I[DIM][DIM];
//    for (int i = 0; i < DIM; i++) {
//        for (int j = 0; j < DIM; j++) {
//            I[i][j] = 0.0;
//	    if (i==j)
//                I[i][j] = 1.0;
//	}
//    }
//    for (int i = 0; i < DIM; i++) {
//        for (int j = i; j < DIM; j++) {
//            for (int k = 0; k < DIM; k++) {
//                for (int l = 0; l < DIM; l++) {
//                    if (i==j){
//                        w[i*DIM + k][j*DIM + l] = I[k][l] - dt*Omega[k][l] - dt*Omega[l][k];
                        //w[i*DIM + k][j*DIM + l] = I[k][l] - dt*OmegaT[i][i]*I[k][l] - dt*Omega[i][i]*I[k][l];
                        //w[i*DIM + k][j*DIM + l] = I[k][l] - 2.0*dt*Omega[i][i]*I[k][l];
//                    }
//                    else{
		        //w[i*DIM + k][j*DIM + l] =  - dt*Omega[i][j]*I[k][l] - dt*Omega[j][i]*I[k][l];
//                w[i*DIM + k][j*DIM + l] =  0.0;
                        //w[j*DIM + l][i*DIM + k] = - dt*Omega[j][i]*I[k][l] - dt*Omega[i][j]*I[k][l];
//                        w[j*DIM + l][i*DIM + k] = 0.0;
//                    }
//                }
//            }
//        }
//    }
//}

// *******************************************************************
// Navier-Stokes step elements
// *******************************************************************

// *******************************************************************
// Transport equations of the density numbers of species A and B
// *******************************************************************

// Transport equation of the density numbers of specie A using Euler explicit method
void higflow_explicit_euler_shear_banding_transport_equation_nA(higflow_solver *ns) {
    if (ns->contr.rheotype == VCM) {
        // Get the cosntants
        real DeA    = ns->ed.vesb.par.DeA;
        real epsilon   = ns->ed.vesb.par.epsilon;
        real PeA   = ns->ed.vesb.par.PeA;
        real CAeq   = ns->ed.vesb.par.CAeq;
        real CBeq   = ns->ed.vesb.par.CBeq;
        real chi   = ns->ed.vesb.par.chi;
        real nAmax = -1.0e16;
        real nAmin = 1.0e16;
        // Get the local sub-domain for the cells
        sim_domain *sdnA  = psd_get_local_domain(ns->ed.vesb.psdSBnA);
        sim_domain *sdnB = psd_get_local_domain(ns->ed.vesb.psdSBnB);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int i = 0; i < DIM; i++) {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdnA);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdnA); !higcit_isfinished(it); higcit_nextcell(it)) {
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
            // Get the velocity derivative tensor Du, S and B tensor
            //real Du[DIM][DIM], A[DIM][DIM];
            // Get Du and A tensors
            //for (int i = 0; i < DIM; i++) {
            //    for (int j = 0; j < DIM; j++) {
                    // Get Du
            //        Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpD[i][j], ns->ed.stn);
                    // Get conformation tensor of specie A
            //        A[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpA[i][j], ns->ed.stn);
            //    }
            //}

            //Get the density number nA
            real nA = compute_value_at_point(ns->ed.vesb.sdSBnA, ccenter, ccenter, 1.0, ns->ed.vesb.dpnA, ns->ed.stn);
            //Get the density number nB
            real nB = compute_value_at_point(ns->ed.vesb.sdSBnB, ccenter, ccenter, 1.0, ns->ed.vesb.dpnB, ns->ed.stn);

            //Get the breakage/reformation rate of specie A
            real CA = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpcA, ns->ed.stn);
            //Get the breakage/reformation rate of specie B
            real CB = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpcB, ns->ed.stn);

            // Calculate the rate of deformation tensor
            //real D[DIM][DIM];
            //for (int dim = 0; dim < DIM; dim++) {
            //    for (int dim2 = 0; dim2 < DIM; dim2++) {
                    // Calculate the tensor
            //        D[dim][dim2] = 0.5*(Du[dim][dim2] + Du[dim2][dim]);
            //    }
            //}

            // Calculate RHS
            real RHS;
            // Calculate RHS = cB*nB**2/(2*DeA) -cA*nA/DeA
            RHS = hig_flow_nA_VCM_RHS (CA, CB, nA, nB, DeA);
            // Get the velocity at cell center 
            real u[DIM], dnAdx[DIM];
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            //Get the density number nA
            //real nA = compute_value_at_point(ns->ed.vesb.sdSBnA, ccenter, ccenter, 1.0, ns->ed.vesb.dpnA, ns->ed.stn);
            //Get the density number nB
            //real nB = compute_value_at_point(ns->ed.vesb.sdSBnB, ccenter, ccenter, 1.0, ns->ed.vesb.dpnB, ns->ed.stn);
            // Solving the Transport Equation using the Euler Method
            // Right hand side equation
            real rhs = 0.0;
            switch (ns->ed.vesb.contr.nAnBconvecdiscrtype) {
                // Central scheme
                case 0: 
                    // Derivative at cell center
                    hig_flow_derivative_nA_at_center_cell(ns, ccenter, cdelta, nA, dnAdx);
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_shear_banding_VCM_model(ns, sdnA, clid, ccenter, cdelta, dim, ns->ed.vesb.dpnA);
                        //Compute convective term in rhs
                        rhs    -= u[dim]*dnAdx[dim];
                        // Compute the diffusive term rhs
                        rhs    += higflow_diffusive_shear_banding_nA_term(ns);
                    }
                    break;
                // CUBISTA scheme
                case 1: 
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_shear_banding_VCM_model(ns, sdnA, clid, ccenter, cdelta, dim, ns->ed.vesb.dpnA);
                        //Compute convective ionic term CUBISTA in rhs
                        rhs    -= higflow_convective_term_shear_banding_VCM_cubista(ns, ns->dpu[dim], ns->ed.vesb.dpnA, ns->ed.vesb.sdSBnA, ns->ed.stn, nA, ccenter, cdelta, dim);
                        // Compute the diffusive term rhs
                        rhs    += higflow_diffusive_shear_banding_nA_term(ns);
                    }
                    break;
            }
            // Compute the final rhs (Convective tern + RHS + Diffusive term)
            rhs += RHS;
            // Compute the final value step time
            real newnA =  nA + ns->par.dt * rhs;
            if (newnA > nAmax) nAmax = newnA;
            if (newnA < nAmin) nAmin = newnA;
            // Set property value  
            dp_set_value(ns->ed.vesb.dpnA, clid, newnA);
        }

        //Printing the min and max structural parameter values
        printf("===> nAmin = %lf <===> nAmax = %lf <===\n", nAmin, nAmax);
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed ionic property
        dp_sync(ns->ed.vesb.dpnA);
    }
}


// Transport equation of the density numbers of specie B using Euler explicit method
void higflow_explicit_euler_shear_banding_transport_equation_nB(higflow_solver *ns) {
    if (ns->contr.rheotype == VCM) {
        // Get the cosntants
        real DeA    = ns->ed.vesb.par.DeA;
        real epsilon   = ns->ed.vesb.par.epsilon;
        real PeA   = ns->ed.vesb.par.PeB;
        real CAeq   = ns->ed.vesb.par.CAeq;
        real CBeq   = ns->ed.vesb.par.CBeq;
        real chi   = ns->ed.vesb.par.chi;
        real nBmax = -1.0e16;
        real nBmin = 1.0e16;
        // Get the local sub-domain for the cells
        sim_domain *sdnA  = psd_get_local_domain(ns->ed.vesb.psdSBnA);
        sim_domain *sdnB = psd_get_local_domain(ns->ed.vesb.psdSBnB);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int i = 0; i < DIM; i++) {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdnB);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdnB); !higcit_isfinished(it); higcit_nextcell(it)) {
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
            // Get the velocity derivative tensor Du, S and B tensor
            //real Du[DIM][DIM], A[DIM][DIM];
            // Get Du and A tensors
            //for (int i = 0; i < DIM; i++) {
            //    for (int j = 0; j < DIM; j++) {
                    // Get Du
            //        Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpD[i][j], ns->ed.stn);
                    // Get conformation tensor of specie A
            //        A[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpA[i][j], ns->ed.stn);
            //    }
            //}

            //Get the density number nA
            real nA = compute_value_at_point(ns->ed.vesb.sdSBnA, ccenter, ccenter, 1.0, ns->ed.vesb.dpnA, ns->ed.stn);
            //Get the density number nB
            real nB = compute_value_at_point(ns->ed.vesb.sdSBnB, ccenter, ccenter, 1.0, ns->ed.vesb.dpnB, ns->ed.stn);

            //Get the breakage/reformation rate of specie A
            real CA = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpcA, ns->ed.stn);
            //Get the breakage/reformation rate of specie B
            real CB = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpcB, ns->ed.stn);

            // Calculate the rate of deformation tensor
            //real D[DIM][DIM];
            //for (int dim = 0; dim < DIM; dim++) {
            //    for (int dim2 = 0; dim2 < DIM; dim2++) {
                    // Calculate the tensor
            //        D[dim][dim2] = 0.5*(Du[dim][dim2] + Du[dim2][dim]);
            //    }
            //}

            // Calculate RHS
            real RHS;
            // Calculate RHS = cB*nB**2/(2*DeA) -cA*nA/DeA
            RHS = hig_flow_nB_VCM_RHS (CA, CB, nA, nB, DeA);
            // Get the velocity at cell center 
            real u[DIM], dnBdx[DIM];
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            // Solving the Transport Equation using the Euler Method
            // Right hand side equation
            real rhs = 0.0;
            switch (ns->ed.vesb.contr.nAnBconvecdiscrtype) {
                // Central scheme
                case 0: 
                    // Ionic derivative at cell center
                    hig_flow_derivative_nB_at_center_cell(ns, ccenter, cdelta, nB, dnBdx);
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_shear_banding_VCM_model(ns, sdnB, clid, ccenter, cdelta, dim, ns->ed.vesb.dpnB);
                        //Compute convective inonic term in rhs
                        rhs    -= u[dim]*dnBdx[dim];
                        // Compute the diffusive term rhs
                        rhs    += higflow_diffusive_shear_banding_nB_term(ns);
                    }
                    break;
                // CUBISTA scheme
                case 1: 
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_shear_banding_VCM_model(ns, sdnB, clid, ccenter, cdelta, dim, ns->ed.vesb.dpnB);
                        //Compute convective ionic term CUBISTA in rhs
                        rhs    -= higflow_convective_term_shear_banding_VCM_cubista(ns, ns->dpu[dim], ns->ed.vesb.dpnB, ns->ed.vesb.sdSBnB, ns->ed.stn, nB, ccenter, cdelta, dim);
                        // Compute the diffusive term rhs
                        rhs    += higflow_diffusive_shear_banding_nB_term(ns);
                    }
                    break;
            }
            // Compute the final rhs (Convective tern + RHS + Diffusive term)
            rhs += RHS;
            // Compute the final value step time
            real newnB =  nB + ns->par.dt * rhs;
            if (newnB > nBmax) nBmax = newnB;
            if (newnB < nBmin) nBmin = newnB;
            // Set property value  
            dp_set_value(ns->ed.vesb.dpnB, clid, newnB);
        }

        //Printing the min and max structural parameter values
        printf("===> nBmin = %lf <===> nBmax = %lf <===\n", nBmin, nBmax);
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed ionic property
        dp_sync(ns->ed.vesb.dpnB);
    }
}


// Transport equation of the density numbers of specie A using Euler implicit method
void higflow_implicit_euler_shear_banding_transport_equation_nA(higflow_solver *ns) {
    if (ns->contr.rheotype == VCM) {
        // Get the cosntants
        real DeA    = ns->ed.vesb.par.DeA;
        real epsilon   = ns->ed.vesb.par.epsilon;
        real PeA   = ns->ed.vesb.par.PeB;
        real CAeq   = ns->ed.vesb.par.CAeq;
        real CBeq   = ns->ed.vesb.par.CBeq;
        real chi   = ns->ed.vesb.par.chi;
        real nAmax = -1.0e16;
        real nAmin = 1.0e16;
        // Get the local sub-domain for the cells
        sim_domain *sdnA  = psd_get_local_domain(ns->ed.vesb.psdSBnA);
        sim_domain *sdnB = psd_get_local_domain(ns->ed.vesb.psdSBnB);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int i = 0; i < DIM; i++) {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdnA);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdnA); !higcit_isfinished(it); higcit_nextcell(it)) {
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
            // Get the velocity derivative tensor Du, S and B tensor
            //real Du[DIM][DIM], A[DIM][DIM];
            // Get Du and A tensors
            //for (int i = 0; i < DIM; i++) {
            //    for (int j = 0; j < DIM; j++) {
                    // Get Du
            //        Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpD[i][j], ns->ed.stn);
                    // Get conformation tensor of specie A
            //        A[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpA[i][j], ns->ed.stn);
            //    }
            //}

            //Get the density number nA
            real nA = compute_value_at_point(ns->ed.vesb.sdSBnA, ccenter, ccenter, 1.0, ns->ed.vesb.dpnA, ns->ed.stn);
            //Get the density number nB
            real nB = compute_value_at_point(ns->ed.vesb.sdSBnB, ccenter, ccenter, 1.0, ns->ed.vesb.dpnB, ns->ed.stn);

            //Get the breakage/reformation rate of specie A
            real CA = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpcA, ns->ed.stn);
            //Get the breakage/reformation rate of specie B
            real CB = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpcB, ns->ed.stn);

            // Calculate the rate of deformation tensor
            //real D[DIM][DIM];
            //for (int dim = 0; dim < DIM; dim++) {
            //    for (int dim2 = 0; dim2 < DIM; dim2++) {
                    // Calculate the tensor
            //        D[dim][dim2] = 0.5*(Du[dim][dim2] + Du[dim2][dim]);
            //    }
            //}

            // Calculate RHS
            real RHS = 0.0;
            // Calculate RHS = cB*nB**2/(2*DeA) -cA*nA/DeA
            RHS = hig_flow_nA_VCM_RHS (CA, CB, nA, nB, DeA);
            // Get the velocity at cell center 
            real u[DIM], dnAdx[DIM];
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            //Get the density number nA
            //real nA = compute_value_at_point(ns->ed.vesb.sdSBnA, ccenter, ccenter, 1.0, ns->ed.vesb.dpnA, ns->ed.stn);
            //Get the density number nB
            //real nB = compute_value_at_point(ns->ed.vesb.sdSBnB, ccenter, ccenter, 1.0, ns->ed.vesb.dpnB, ns->ed.stn);
            // Solving the Transport Equation using the Euler Method
            // Right hand side equation
            real rhs = 0.0;
            switch (ns->ed.vesb.contr.nAnBconvecdiscrtype) {
                // Central scheme
                case 0: 
                    // Derivative at cell center
                    hig_flow_derivative_nA_at_center_cell(ns, ccenter, cdelta, nA, dnAdx);
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_shear_banding_VCM_model(ns, sdnA, clid, ccenter, cdelta, dim, ns->ed.vesb.dpnA);
                        //Compute convective term in rhs
                        rhs    -= u[dim]*dnAdx[dim];
                    }
                    break;
                // CUBISTA scheme
                case 1: 
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_shear_banding_VCM_model(ns, sdnA, clid, ccenter, cdelta, dim, ns->ed.vesb.dpnA);
                        //Compute convective ionic term CUBISTA in rhs
                        rhs    -= higflow_convective_term_shear_banding_VCM_cubista(ns, ns->dpu[dim], ns->ed.vesb.dpnA, ns->ed.vesb.sdSBnA, ns->ed.stn, nA, ccenter, cdelta, dim);
                    }
                    break;
            }
            // Compute the final rhs (Convective tern + RHS + Diffusive term)
            rhs += RHS;
            rhs *= ns->par.dt;
            rhs += ns->cc.nABcell;
            // Reset the stencil
            stn_reset(ns->ed.stn);
            // Set the right side of stencil
            stn_set_rhs(ns->ed.stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w  = -ns->par.dt/(ns->ed.vesb.par.PeA*cdelta[dim2]*cdelta[dim2]);
                alpha  -= 4.0*w ;
                Point p;
                POINT_ASSIGN(p, ccenter);
                // Stencil point update: right point
                p[dim2] = ccenter[dim2] + cdelta[dim2];
                sd_get_stencil(sdnA, ccenter, p, w, ns->ed.stn);
                // Stencil point update: left point
                p[dim2] = ccenter[dim2] - cdelta[dim2];
                sd_get_stencil(sdnA, ccenter, p, w, ns->ed.stn);
            }
            alpha = 1.0 + alpha;
            // Get the stencil
            sd_get_stencil(sdnA, ccenter, ccenter,alpha, ns->ed.stn);
            // Get the index of the stencil
            int *ids   = psd_stn_get_gids(ns->ed.vesb.psdSBnA, ns->ed.stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->ed.stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->ed.stn);
	        int cgid = psd_get_global_id(ns->ed.vesb.psdSBnA, c);
            // Set the right side of solver linear system
            slv_set_bi(ns->ed.vesb.slvnA, cgid, stn_get_rhs(ns->ed.stn));
            // Set the line of matrix of the solver linear system
            slv_set_Ai(ns->ed.vesb.slvnA, cgid, numelems, ids, vals);
        }

        // Destroy the iterator
        higcit_destroy(it);
        // Assemble the solver
        slv_assemble(ns->ed.vesb.slvnA);
        // Solve the linear system
        slv_solve(ns->ed.vesb.slvnA);
        //Load property from solver
        dp_slv_load_from_solver(ns->ed.vesb.dpnA, ns->ed.vesb.slvnA);
        // Syncing the distributed property
        dp_sync(ns->ed.vesb.dpnA);
    }
}


// Transport equation of the density numbers of specie B using Euler implicit method
void higflow_implicit_euler_shear_banding_transport_equation_nB(higflow_solver *ns) {
    if (ns->contr.rheotype == VCM) {
        // Get the cosntants
        real DeA    = ns->ed.vesb.par.DeA;
        real epsilon   = ns->ed.vesb.par.epsilon;
        real PeA   = ns->ed.vesb.par.PeB;
        real CAeq   = ns->ed.vesb.par.CAeq;
        real CBeq   = ns->ed.vesb.par.CBeq;
        real chi   = ns->ed.vesb.par.chi;
        real nBmax = -1.0e16;
        real nBmin = 1.0e16;
        // Get the local sub-domain for the cells
        sim_domain *sdnA  = psd_get_local_domain(ns->ed.vesb.psdSBnA);
        sim_domain *sdnB = psd_get_local_domain(ns->ed.vesb.psdSBnB);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int i = 0; i < DIM; i++) {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdnB);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdnB); !higcit_isfinished(it); higcit_nextcell(it)) {
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
            // Get the velocity derivative tensor Du, S and B tensor
            //real Du[DIM][DIM], A[DIM][DIM];
            // Get Du and A tensors
            //for (int i = 0; i < DIM; i++) {
            //    for (int j = 0; j < DIM; j++) {
                    // Get Du
            //        Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpD[i][j], ns->ed.stn);
                    // Get conformation tensor of specie A
            //        A[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpA[i][j], ns->ed.stn);
            //    }
            //}

            //Get the density number nA
            real nA = compute_value_at_point(ns->ed.vesb.sdSBnA, ccenter, ccenter, 1.0, ns->ed.vesb.dpnA, ns->ed.stn);
            //Get the density number nB
            real nB = compute_value_at_point(ns->ed.vesb.sdSBnB, ccenter, ccenter, 1.0, ns->ed.vesb.dpnB, ns->ed.stn);

            //Get the breakage/reformation rate of specie A
            real CA = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpcA, ns->ed.stn);
            //Get the breakage/reformation rate of specie B
            real CB = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.vesb.dpcB, ns->ed.stn);

            // Calculate the rate of deformation tensor
            //real D[DIM][DIM];
            //for (int dim = 0; dim < DIM; dim++) {
            //    for (int dim2 = 0; dim2 < DIM; dim2++) {
                    // Calculate the tensor
            //        D[dim][dim2] = 0.5*(Du[dim][dim2] + Du[dim2][dim]);
            //    }
            //}

            // Calculate RHS
            real RHS = 0.0;
            // Calculate RHS = cB*nB**2/(2*DeA) -cA*nA/DeA
            RHS = hig_flow_nB_VCM_RHS (CA, CB, nA, nB, DeA);
            // Get the velocity at cell center 
            real u[DIM], dnBdx[DIM];
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            // Solving the Transport Equation using the Euler Method
            // Right hand side equation
            real rhs = 0.0;
            switch (ns->ed.vesb.contr.nAnBconvecdiscrtype) {
                // Central scheme
                case 0: 
                    // Ionic derivative at cell center
                    hig_flow_derivative_nB_at_center_cell(ns, ccenter, cdelta, nB, dnBdx);
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_shear_banding_VCM_model(ns, sdnB, clid, ccenter, cdelta, dim, ns->ed.vesb.dpnB);
                        //Compute convective inonic term in rhs
                        rhs    -= u[dim]*dnBdx[dim];
                    }
                    break;
                // CUBISTA scheme
                case 1: 
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_shear_banding_VCM_model(ns, sdnB, clid, ccenter, cdelta, dim, ns->ed.vesb.dpnB);
                        //Compute convective ionic term CUBISTA in rhs
                        rhs    -= higflow_convective_term_shear_banding_VCM_cubista(ns, ns->dpu[dim], ns->ed.vesb.dpnB, ns->ed.vesb.sdSBnB, ns->ed.stn, nB, ccenter, cdelta, dim);
                    }
                    break;
            }
            // Compute the final rhs (Convective tern + RHS)
            rhs += RHS;
            rhs *= ns->par.dt;
            rhs += ns->cc.nABcell;
            // Reset the stencil
            stn_reset(ns->ed.stn);
            // Set the right side of stencil
            stn_set_rhs(ns->ed.stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w  = -ns->par.dt/(ns->ed.vesb.par.PeB*cdelta[dim2]*cdelta[dim2]);
                alpha  -= 4.0*w ;
                Point p;
                POINT_ASSIGN(p, ccenter);
                // Stencil point update: right point
                p[dim2] = ccenter[dim2] + cdelta[dim2];
                sd_get_stencil(sdnB, ccenter, p, w, ns->ed.stn);
                // Stencil point update: left point
                p[dim2] = ccenter[dim2] - cdelta[dim2];
                sd_get_stencil(sdnB, ccenter, p, w, ns->ed.stn);
            }
            alpha = 1.0 + alpha;
            // Get the stencil
            sd_get_stencil(sdnB, ccenter, ccenter,alpha, ns->ed.stn);
            // Get the index of the stencil
            int *ids   = psd_stn_get_gids(ns->ed.vesb.psdSBnB, ns->ed.stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->ed.stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->ed.stn);
	        int cgid = psd_get_global_id(ns->ed.vesb.psdSBnB, c);
            // Set the right side of solver linear system
            slv_set_bi(ns->ed.vesb.slvnB, cgid, stn_get_rhs(ns->ed.stn));
            // Set the line of matrix of the solver linear system
            slv_set_Ai(ns->ed.vesb.slvnB, cgid, numelems, ids, vals);

        }

        // Destroy the iterator
        higcit_destroy(it);
        // Assemble the solver
        slv_assemble(ns->ed.vesb.slvnB);
        // Solve the linear system
        slv_solve(ns->ed.vesb.slvnB);
        //Load property from solver
        dp_slv_load_from_solver(ns->ed.vesb.dpnB, ns->ed.vesb.slvnB);
        // Syncing the distributed property
        dp_sync(ns->ed.vesb.dpnB);
    }
}


// *******************************************************************
// Calculate convective term of the transport equation of species A and B using the CUBISTA method
// *******************************************************************
real higflow_convective_term_shear_banding_VCM_cubista(higflow_solver *ns, distributed_property *dpu, distributed_property *dpn, sim_domain *sdp, sim_stencil *stn, real n, Point ccenter, Point cdelta, int dim) {
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

// Get the derivative of nA 
void hig_flow_derivative_nA_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, real ncenter, real dndx[DIM]) {
    for (int dim = 0; dim < DIM; dim++) {
        int incell_left, incell_right;
        // Get the concentration in the left cell
        real nleft = compute_center_p_left_22(ns->ed.vesb.sdSBnA, ccenter, cdelta, dim, 1.0, ns->ed.vesb.dpnA, ns->ed.stn, &incell_left);
        // Get the concentration in the right cell
        real nright = compute_center_p_right_22(ns->ed.vesb.sdSBnA, ccenter, cdelta, dim, 1.0, ns->ed.vesb.dpnA, ns->ed.stn, &incell_left);
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

// Get the derivative of nB 
void hig_flow_derivative_nB_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, real ncenter, real dndx[DIM]) {
    for (int dim = 0; dim < DIM; dim++) {
        int incell_left, incell_right;
        // Get the concentration in the left cell
        real nleft = compute_center_p_left_22(ns->ed.vesb.sdSBnB, ccenter, cdelta, dim, 1.0, ns->ed.vesb.dpnB, ns->ed.stn, &incell_left);
        // Get the concentration in the right cell
        real nright = compute_center_p_right_22(ns->ed.vesb.sdSBnB, ccenter, cdelta, dim, 1.0, ns->ed.vesb.dpnB, ns->ed.stn, &incell_left);
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




// Calculate RHS = cB*nB**2/(2*DeA) -cA*nA/DeA
real hig_flow_nA_VCM_RHS (real CA, real CB, real nA, real nB, real De) {
    real F1S, F2S, RHS;
    F1S = 0.5*CB*(pow(nB,2))/De;
    //F1S = 0.5*CB*nB*nB/De;
    F2S = CA*nA/De;
    RHS = F1S - F2S;
    return RHS;
}

// Calculate RHS = 2*cA*nA/DeA - cB*nB**2/DeA
real hig_flow_nB_VCM_RHS (real CA, real CB, real nA, real nB, real De) {
    real F1S, F2S, RHS;
    F1S = 2.0*CA*nA/De;
    F2S = CB*(pow(nB,2))/De;
    //F2S = CB*nB*nB/De;
    RHS = F1S - F2S;
    return RHS;
}

// *******************************************************************
// Navier-Stokes Step for the Explicit Euler Method
// *******************************************************************
void higflow_explicit_euler_intermediate_velocity_viscoelastic_shear_banding(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]) {
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
            higflow_computational_cell_viscoelastic_shear_banding(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
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
void higflow_explicit_runge_kutta_2_intermediate_velocity_shear_banding(higflow_solver *ns) {
    // Calculate the auxiliar velocity by the explicit Euler method
    higflow_explicit_euler_intermediate_velocity_viscoelastic_shear_banding(ns, ns->dpu, ns->dpuaux);
    // Calculate the star velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_viscoelastic_shear_banding(ns, ns->dpuaux, ns->dpustar);
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
void higflow_explicit_runge_kutta_3_intermediate_velocity_shear_banding(higflow_solver *ns) {
    // Calculate the auxiliar velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_viscoelastic_shear_banding(ns, ns->dpu, ns->dpuaux);
    // Calculate the second stage velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_viscoelastic_shear_banding(ns, ns->dpuaux, ns->dpustar);
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
    higflow_explicit_euler_intermediate_velocity_viscoelastic_shear_banding(ns, ns->dpuaux, ns->dpustar);
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
void higflow_semi_implicit_euler_intermediate_velocity_shear_banding(higflow_solver *ns) {
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
            higflow_computational_cell_viscoelastic_shear_banding(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
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
                real w  = -(ns->par.dt)*(1.0 + ns->ed.vesb.par.beta)*(ns->ed.vesb.par.De)/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
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
void higflow_semi_implicit_crank_nicolson_intermediate_velocity_shear_banding(higflow_solver *ns) {
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
            higflow_computational_cell_viscoelastic_shear_banding(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
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
void higflow_semi_implicit_bdf2_intermediate_velocity_shear_banding(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]) {
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
            higflow_computational_cell_viscoelastic_shear_banding(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
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
            higflow_computational_cell_viscoelastic_shear_banding(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            //Get the uaux
            real uaux = dp_get_value(ns->dpuaux[dim], flid);
            // Right hand side equation
            real rhs = 0.0;
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
void higflow_solver_step_viscoelastic_shear_banding(higflow_solver *ns) {
    // Boundary condition for velocity
    higflow_boundary_condition_for_velocity(ns);
        // Calculate the source term
    higflow_calculate_source_term(ns);
    // Calculate the facet source term
    higflow_calculate_facet_source_term(ns);
    // Set outflow for velocity
    //higflow_outflow_u_step(ns);
    // Calculate the intermediate velocity
    switch (ns->contr.tempdiscrtype) {
        case 0:
           // Explicit Euler method
           higflow_explicit_euler_intermediate_velocity_viscoelastic_shear_banding(ns, ns->dpu, ns->dpustar);
           break;
        case 1:
           // Explicit RK2 method
           higflow_explicit_runge_kutta_2_intermediate_velocity_shear_banding(ns);
           break;
        case 2:
           // Explicit RK3 method
           higflow_explicit_runge_kutta_3_intermediate_velocity_shear_banding(ns);
           break;
        case 3:
           // Semi-Implicit Euler Method
           higflow_semi_implicit_euler_intermediate_velocity_shear_banding(ns);
           break;
        case 4:
           // Semi-Implicit Crank-Nicolson Method
           higflow_semi_implicit_crank_nicolson_intermediate_velocity_shear_banding(ns);
           break;
        case 5:
           // Semi-Implicit Crank-Nicolson Method
           higflow_semi_implicit_bdf2_intermediate_velocity_shear_banding(ns, ns->dpu, ns->dpustar);
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
    // Set outflow for u velocity 
    //higflow_outflow_u_step(ns);
    // Calculate the final pressure
    higflow_final_pressure(ns);
    // Calculate the velocity derivative tensor
    higflow_compute_velocity_derivative_tensor(ns);
    // Solve the transport equations of the density number of species A and B
    if (ns->contr.rheotype == VCM) {
        //Calculate the breakage and reformation rates of species A and B
        higflow_compute_viscoelastic_shear_banding_cA_VCM(ns);
        higflow_compute_viscoelastic_shear_banding_cB_VCM(ns);
        //Type of the Discretization of the transport equations to be selected
        switch (ns->ed.vesb.contr.nAnBdiscrtype) {
            case 0: //Explicit Euler Method
                //VCM model
                switch (ns->ed.contr.rheotype) {
                    case 0:
                    //Standard VCM model
                    higflow_explicit_euler_shear_banding_transport_equation_nA(ns);
                    higflow_explicit_euler_shear_banding_transport_equation_nB(ns);
                    break;
                    case 1:
                    //New VCM model that could be implemented here
                    higflow_explicit_euler_shear_banding_transport_equation_nA(ns);
                    higflow_explicit_euler_shear_banding_transport_equation_nB(ns);           
                    break;
                }
            break;
            case 1: //Implicit Euler Method
                //VCM model
                switch (ns->ed.contr.rheotype) {
                    case 0:
                    //Standard VCM model
                    higflow_implicit_euler_shear_banding_transport_equation_nA(ns);
                    higflow_implicit_euler_shear_banding_transport_equation_nB(ns);
                    break;
                    case 1:
                    //New VCM model that could be implemented here
                    higflow_explicit_euler_shear_banding_transport_equation_nA(ns);
                    higflow_explicit_euler_shear_banding_transport_equation_nB(ns);           
                    break;
                }
            break;    
        }
        //Calculate the new breakage rate of specie A
        //higflow_compute_viscoelastic_shear_banding_cA_VCM(ns);

    }
    // Constitutive Equation Step for the Explicit Euler Method
    switch (ns->ed.vesb.contr.discrtype) {
        case 0:
           // Explicit method
           //Solving the conformation tensor equation of specie A
           higflow_explicit_euler_conformation_tensor_A(ns);
           //Solving the conformation tensor equation of specie B
           higflow_explicit_euler_conformation_tensor_B(ns);
           break;
        case 1: 
           // Implicit method
           //Solving the conformation tensor equation of specie A
           higflow_implicit_euler_conformation_tensor_A(ns);
           //Solving the conformation tensor equation of specie B
           higflow_implicit_euler_conformation_tensor_B(ns);
           break;
    }
    // Computing the Polymeric Tensor
    higflow_compute_polymeric_tensor_shear_banding(ns);
}

