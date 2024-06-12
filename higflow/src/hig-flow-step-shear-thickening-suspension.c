// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Step Viscoelastic - version 10/11/2016
// *******************************************************************
// *******************************************************************

#include "hig-flow-step-shear-thickening-suspension.h"
#include "hig-flow-mittag-leffler.h"
// *******************************************************************
// Constitutive Equations
// *******************************************************************

// Computing the polymeric Tensor S = Tau - 2*eta0*D
void higflow_compute_polymeric_tensor_shear_thickening_suspension(higflow_solver *ns)
{
    if ((ns->contr.flowtype == 9))
    {
        // Get the constants
        // real Re   = ns->par.Re;
        real alpha = ns->ed.stsp.par.alpha;
        real eta0 = ns->ed.stsp.par.eta0;
        real chij1 = ns->ed.stsp.par.chij1;
        real chij2 = ns->ed.stsp.par.chij2;
        real X0 = ns->ed.stsp.par.X0;
        real Pic = ns->ed.stsp.par.Pic;
        real phircp = ns->ed.stsp.par.phi;
        real fmax = -1.0e16;
        real fmin = 1.0e16;
        real chimax = -1.0e16;
        real chimin = 1.0e16;
        real chiJmax = -1.0e16;
        real chiJmin = 1.0e16;
        real alphamax = -1.0e16;
        real alphamin = 1.0e16;
        real small = 1.0e-10;
        real small2 = 1.0e-6;
        real Smax[DIM][DIM], Smin[DIM][DIM];
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                Smax[i][j] = -1.0e16;
                Smin[i][j] = 1.0e16;
            }
        }
        real EEmax[DIM][DIM], EEmin[DIM][DIM];
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                EEmax[i][j] = -1.0e16;
                EEmin[i][j] = 1.0e16;
            }
        }
        real ECmax[DIM][DIM], ECmin[DIM][DIM];
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                ECmax[i][j] = -1.0e16;
                ECmin[i][j] = 1.0e16;
            }
        }
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdp);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid = mp_lookup(mp, hig_get_cid(c));
            // Get the inside/outside inflow point cell
            int inflowcell, outflowcell;
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Get the velocity derivative tensor Du and the microstructure tensor
            real A[DIM][DIM], Du[DIM][DIM], S0[DIM][DIM];
            // AM is the Trace of the microstructure tensor
            real AM = 0.0;
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    // Get Du
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
                    // Get the microstructure tensor
                    A[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpA[i][j], ns->ed.stn);
                    // Get the previous value of polymeric tensor, S0
                    S0[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpS[i][j], ns->ed.stn);
                }
                AM += A[i][i];
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
            // DEBUG_INSPECT(Du[0][1],%lf);
            // DEBUG_INSPECT(E[1][0],%lf);
            // Calculate T of previous step, T0
            real T0[DIM][DIM];
            real Pis = 0.0;
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    if (ns->ed.stsp.contr.model == GW_WC_IF){
                        T0[i][j] = S0[i][j];
                    }
                    else {
                        T0[i][j] = S0[i][j] + 2.0 * eta0 * E[i][j];
                    }
                }
                Pis += T0[i][i];
            }
            Pis = -1.0 * Pis / 3.0;
            // Pis = -1.0/3.0;
            // DEBUG_INSPECT(Pis,%lf);
            //  Eigen-values and eigen-vectors of the deformation tensor E
            real R[DIM][DIM], lambda[DIM];
            // DEBUG_INSPECT(E[0][1],%lf);
            // DEBUG_INSPECT(E[1][0],%lf);
            hig_flow_jacobi(E, lambda, R);
            // DEBUG_INSPECT(lambda[0],%lf);
            //DEBUG_INSPECT(R[0][0],%lf);
            //DEBUG_INSPECT(R[0][1],%lf);
            //DEBUG_INSPECT(R[0][2],%lf);

            // DEBUG_INSPECT(lambda[1],%lf);
            // Calculation of the compressive rate of strain Ec and the extensional rate of strain Ee
            real E_C[DIM][DIM], E_E[DIM][DIM];
            hig_flow_compressive_extensional_tensors(lambda, R, E_C, E_E);
            // DEBUG_INSPECT(E_E[0][0],%lf);
            //DEBUG_INSPECT(E_E[0][0],%lf);
            //DEBUG_INSPECT(E_E[0][1],%lf);
            //DEBUG_INSPECT(E_E[0][2],%lf);
            //  Calculate the max and min values of T
            // for (int i = 0; i < DIM; i++) {
            //     for (int j = 0; j < DIM; j++) {
            //         real ETEST = E_E[i][j];
            //         if (ETEST > EEmax[i][j]) EEmax[i][j] = ETEST;
            //         if (ETEST < EEmin[i][j]) EEmin[i][j] = ETEST;
            //         if (E_C[i][j] > ECmax[i][j]) ECmax[i][j] = E_C[i][j];
            //         if (E_C[i][j] < ECmin[i][j]) ECmin[i][j] = E_C[i][j];
            //     }
            // }
            // for (int i = 0; i < DIM; i++) {
            //     for (int j = 0; j < DIM; j++) {
            //         real ETEST = E_C[i][j];
            //         if (ETEST > ECmax[i][j]) ECmax[i][j] = ETEST;
            //         if (ETEST < ECmin[i][j]) ECmin[i][j] = ETEST;
            //         if (E_C[i][j] > ECmax[i][j]) ECmax[i][j] = E_C[i][j];
            //         if (E_C[i][j] < ECmin[i][j]) ECmin[i][j] = E_C[i][j];
            //     }
            // }
            // DEBUG_INSPECT(E_C[1][0],%lf);
            // Kronecker's delta function
            real KD[DIM][DIM];
            // real KM;
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    if (i == j)
                    {
                        KD[i][i] = 1.0;
                    }
                    else
                    {
                        KD[i][j] = 0.0;
                    }
                    // printf("%lf     ", KD[i][j]);
                }
                // KM += KD[i][i];
                // printf("\n")
            }
            // DEBUG_INSPECT(KD[1][0],%lf);

            //******************************************************************************************+
            // Tests and debugging
            // real lambda1[DIM];
            // real gd1 = 1.0;
            // lambda1[0] = -1.0*gd1;
            // lambda1[1] = 1.0*gd1;
            // lambda1[2] = 0.0;

            // real R0[DIM][DIM];
            // for (int i =0; i<DIM;i++){
            //     for (int j=0;j<DIM;j++){
            //         R0[i][j] = 0.0;
            //     }
            // }

            // R0[0][0]= 1.0/sqrt(2);
            // R0[0][1]= 1.0/sqrt(2);
            // R0[1][0] = -R0[0][1];
            // R0[1][1] = 1.0/sqrt(2);

            // Calculation of the compressive rate of strain Ec and the extensional rate of strain Ee
            // real E_C1[DIM][DIM], E_E1[DIM][DIM];
            // hig_flow_compressive_extensional_tensors (lambda1, R0, E_C1, E_E1);
            // DEBUG_INSPECT(E_E1[1][1],%lf);
            // DEBUG_INSPECT(E_E1[2][2],%lf);
            // real p[DIM];
            // p[0] = -1.0/sqrt(2.0);
            // p[1] = -p[0];
            // p[2] = 0.0;
            // real PP[DIM][DIM];
            // real E0[DIM][DIM];
            // real PP0 =0.0;
            // for (int i=0;i<DIM;i++){
            //     for (int j=0;j<DIM;j++){
            //         PP[i][j] = p[i]*p[j];
            //     }
            //     PP0 += PP[i][i];
            // }
            // DEBUG_INSPECT(PP0,%lf);

            // for (int i=0;i<DIM;i++){
            //     for (int j=0;j<DIM;j++){
            //         E0[i][j] = 0.0;
            //     }
            // }
            // E0[0][1] = 1.0;
            // E0[1][0] = 1.0;
            // Calculation of the fourth order orientation moment <pppp>
            // real FOM0[DIM][DIM][DIM][DIM];
            // hig_flow_calculate_4th_order_orientation_moment (PP0, KD, PP, FOM0);
            // for (int i=0;i<DIM;i++){
            //     for (int j=0;j<DIM;j++){
            //         for (int k=0;k<DIM;k++){
            //             for (int l=0;l<DIM;l++){
            // FOM0[i][j][k][l] = PP[i][j]*PP[k][l];
            //                FOM0[i][j][k][l] = 0.0;
            // FOM0[i][j][k][l] = -(1.0/35.0)*AM*(KD[i][j]*KD[k][l] + KD[i][k]*KD[j][l] + KD[i][l]*KD[j][k]) + (1.0/7.0)*(KD[i][j]*a[k][l] + KD[i][k]*a[j][l] + KD[i][l]*a[j][k] + a[i][j]*KD[k][l] + a[i][k]*KD[j][l] + a[i][l]*KD[j][k]);
            //                FOM0[i][j][k][l] += (-1.0)*(1.0/35.0)*PP0*(KD[i][j]*KD[k][l] + KD[i][k]*KD[j][l] + KD[i][l]*KD[j][k]);
            //                FOM0[i][j][k][l] += (1.0/7.0)*(KD[i][j]*PP[k][l] + KD[i][k]*PP[j][l] + KD[i][l]*PP[j][k]);
            //                FOM0[i][j][k][l] += (1.0/7.0)*(PP[i][j]*KD[k][l] + PP[i][k]*KD[j][l] + PP[i][l]*KD[j][k]);
            //            }
            //        }
            //    }
            //}

            // Calculate the suspension stress tensor
            // real TC[DIM][DIM];
            // hig_flow_calculate_suspension_stress_tensor (E0, FOM0, TC);
            // for (int m=0;m<DIM;m++){
            //    for (int n=0;n<DIM;n++){
            //        TC[m][n] = 0.0;
            //        for (int i=0;i<DIM;i++){
            //            for (int j=0;j<DIM;j++){
            // FOM2[m][n] += cons*E[i][j]*FOM[i][j][m][n];
            //                TC[m][n] += E0[i][j]*FOM0[m][n][i][j];
            //            }
            //        }
            //    }
            //}

            // DEBUG_INSPECT(TC[0][2],%lf);
            // DEBUG_INSPECT(TC[2][2],%lf);

            //******************************************************************************************+

            // Calculation of the fourth order orientation moment <pppp>
            real FOM[DIM][DIM][DIM][DIM];
            hig_flow_calculate_4th_order_orientation_moment(AM, KD, A, FOM);
            real tau_M[DIM][DIM];
            if (ns->ed.stsp.contr.model == GW){
                // For GW model
                real E2[DIM][DIM];
                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        E2[i][j] = 0.5 * (Du[i][j] + Du[j][i]);
                    }
                }
                // DEBUG_INSPECT(E2[0][1],%lf);
                // DEBUG_INSPECT(E2[1][0],%lf);
                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        //tau_M[i][j] = 1.0*(ns->ed.stsp.par.alpha) * (ns->ed.stsp.par.eta0) * E2[i][j];
                        tau_M[i][j] = 2.0*(ns->ed.stsp.par.alpha) * (ns->ed.stsp.par.eta0) * E2[i][j];
                    }
                }
                // DEBUG_INSPECT(tau_M[0][1],%lf);
                // DEBUG_INSPECT(tau_M[1][0],%lf);
            } else if (ns->ed.stsp.contr.model == GW_WC_IF) {
                // GW-WC model for inhomogeneous flows
                real E2[DIM][DIM];
                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        E2[i][j] = 0.5 * (Du[i][j] + Du[j][i]);
                    }
                }
                //DEBUG_INSPECT(E2[0][1],%lf);
                //DEBUG_INSPECT(E2[1][0],%lf);
                // Calculation of chi value
                real chi = hig_flow_calculate_jamming_coordinate_chi(A, E_C);
                if (chi > chimax)
                    chimax = chi;
                if (chi < chimin)
                    chimin = chi;
                // DEBUG_INSPECT(chi,%lf);
                // chi = fabs(chi);
                Pis = fabs(Pis);
                // real fPi;
                real fPi = exp(-Pic / (Pis + small2));
                // if (Pis <= small2) {
                //     fPi = exp(-Pic/(Pis+small2));
                // } else {
                //     fPi = exp(-Pic/Pis);
                // }
                // Calculation of function f(Pi)
                // real fPi = exp(-Pic/(Pis+small2));
                // Check to see min and max values of f
                if (fPi > fmax)
                    fmax = fPi;
                if (fPi < fmin)
                    fmin = fPi;
                // DEBUG_INSPECT(fPi,%lf);
                // Calculation of the jamming point Chi_J
                real chi_J = (1.0 - fPi) * chij1 + fPi * chij2;
                /*
                real tol_chi = 0.1;
                real chi2chij = chi/chi_J;
                real tol2chi = 1.0-chi2chij;
                //real tol_chi = 0.977*chi_J;
                //if (chi >= phircp - tol_chi){
                //    chi = phircp - tol_chi;
                //}
                //if (chi >= tol_chi){
                //    chi = tol_chi;
                //}
                if (tol2chi <= tol_chi){
                    chi = (1.0-tol_chi)*chi_J;
                }*/
                //Tolerance value that works
                real tol_chi = 0.85*chi_J;
                //real tol_chi = 0.86*chi_J;
                if (chi >= tol_chi){
                    chi = tol_chi;
                }

                real Xchi = ns->ed.stsp.get_X(ccenter, ns->par.t, X0, chi, chi_J);
                // Check to see min and max values of Xchi
                if (Xchi > chiJmax)
                    chiJmax = Xchi;
                if (Xchi < chiJmin)
                    chiJmin = Xchi;
                // DEBUG_INSPECT(Xchi,%lf);
                Xchi = fabs(Xchi);
                real Xchi2;
                if (Xchi < small2)
                {
                    Xchi2 = small2;
                }
                else
                {
                    Xchi2 = Xchi;
                }

                // DEBUG_INSPECT(E_C[1][0],%lf);
                // DEBUG_INSPECT(Xchi,%lf);
                // DEBUG_INSPECT(E2CT[0][1],%lf);
                // for (int i = 0; i < DIM; i++) {ñ
                //     for (int j = 0; j < DIM; j++) {
                //         tau_M[i][j] = 0.0;
                //     }
                // }

                // Get the volume fraction value
                real varphi = compute_value_at_point(ns->ed.stsp.sdphi, ccenter, ccenter, 1.0, ns->ed.stsp.dpphi, ns->ed.stn);
                real frac_tol = ns->ed.stsp.par.gdrms;
                real eps = 1.0e-4;
                if(varphi > 1.0){
                    varphi = 1.0-eps;
                }
                //if(varphi > phircp - frac_tol){
                //    varphi = phircp - frac_tol;
                //}
                if(varphi < eps){
                    varphi = 0.0;
                }
                //DEBUG_INSPECT(varphi,%lf);
                real alphaVF = ns->ed.stsp.get_alpha(ccenter, ns->par.t, alpha, varphi, phircp);
                //DEBUG_INSPECT(alphaVF,%lf);
                if (alphaVF > alphamax)
                    alphamax = alphaVF;
                if (alphaVF < alphamin)
                    alphamin = alphaVF;
                alphaVF = fabs(alphaVF);
                real alphaVF2;
                if (alphaVF < small2)
                {
                    alphaVF2 = small2;
                }
                else
                {
                    alphaVF2 = alphaVF;
                }

                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        tau_M[i][j] = 0.0;
                        // tau_M[i][j] = eta0*(alpha*E2[i][j] + Xchi*E_C[i][j]);
                        // tau_M[i][j] += eta0*alpha*E2[i][j];
                        // tau_M[i][j] += eta0*Xchi2*E_C[i][j];
                        tau_M[i][j] += eta0 * alphaVF2 * E2[i][j];
                        tau_M[i][j] += eta0 * Xchi2 * E_C[i][j];
                        //tau_M[i][j] *= 2.0;
                        // tau_M[i][j] = (ns->ed.stsp.par.eta0)*(ns->ed.stsp.par.alpha)*E2[i][j];
                        // tau_M[i][j] += (ns->ed.stsp.par.eta0)*(small + Xchi)*E_C[i][j];
                        // tau_M[i][j] = (ns->ed.stsp.par.eta0)*((ns->ed.stsp.par.alpha)*E2[i][j] + E2CT[i][j]);
                    }
                }
            } else {
                // for GW-WC model
                real E2[DIM][DIM];
                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        E2[i][j] = 0.5 * (Du[i][j] + Du[j][i]);
                    }
                }
                // DEBUG_INSPECT(E2[0][1],%lf);
                // DEBUG_INSPECT(E2[1][0],%lf);
                // Calculation of chi value
                real chi = hig_flow_calculate_jamming_coordinate_chi(A, E_C);
                if (chi > chimax)
                    chimax = chi;
                if (chi < chimin)
                    chimin = chi;
                // DEBUG_INSPECT(chi,%lf);
                // chi = fabs(chi);
                Pis = fabs(Pis);
                // real fPi;
                real fPi = exp(-Pic / (Pis + small2));
                // if (Pis <= small2) {
                //     fPi = exp(-Pic/(Pis+small2));
                // } else {
                //     fPi = exp(-Pic/Pis);
                // }
                // Calculation of function f(Pi)
                // real fPi = exp(-Pic/(Pis+small2));
                // Check to see min and max values of f
                if (fPi > fmax)
                    fmax = fPi;
                if (fPi < fmin)
                    fmin = fPi;
                // DEBUG_INSPECT(fPi,%lf);
                // Calculation of the jamming point Chi_J
                real chi_J = (1.0 - fPi) * chij1 + fPi * chij2;
                real tol_chi = 0.15;
                //real tol_chi = 0.01876;
                real chi2chij = chi/chi_J;
                real tol2chi = 1.0-chi2chij;
                //real tol_chi = 0.977*chi_J;
                //if (chi >= phircp - tol_chi){
                //    chi = phircp - tol_chi;
                //}
                //if (chi >= tol_chi){
                //    chi = tol_chi;
                //}
                if (tol2chi <= tol_chi){
                    chi = (1.0-tol_chi)*chi_J;
                }
                real Xchi = ns->ed.stsp.get_X(ccenter, ns->par.t, X0, chi, chi_J);
                // Check to see min and max values of Xchi
                if (Xchi > chiJmax)
                    chiJmax = Xchi;
                if (Xchi < chiJmin)
                    chiJmin = Xchi;
                // DEBUG_INSPECT(Xchi,%lf);
                Xchi = fabs(Xchi);
                real Xchi2;
                if (Xchi < small2)
                {
                    Xchi2 = small2;
                }
                else
                {
                    Xchi2 = Xchi;
                }

                // DEBUG_INSPECT(E_C[1][0],%lf);
                // DEBUG_INSPECT(Xchi,%lf);
                // DEBUG_INSPECT(E2CT[0][1],%lf);
                // for (int i = 0; i < DIM; i++) {ñ
                //     for (int j = 0; j < DIM; j++) {
                //         tau_M[i][j] = 0.0;
                //     }
                // }

                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        tau_M[i][j] = 0.0;
                        // tau_M[i][j] = eta0*(alpha*E2[i][j] + Xchi*E_C[i][j]);
                        // tau_M[i][j] += eta0*alpha*E2[i][j];
                        // tau_M[i][j] += eta0*Xchi2*E_C[i][j];
                        tau_M[i][j] += eta0 * alpha * E2[i][j];
                        tau_M[i][j] += eta0 * Xchi2 * E_C[i][j];
                        // tau_M[i][j] = (ns->ed.stsp.par.eta0)*(ns->ed.stsp.par.alpha)*E2[i][j];
                        // tau_M[i][j] += (ns->ed.stsp.par.eta0)*(small + Xchi)*E_C[i][j];
                        // tau_M[i][j] = (ns->ed.stsp.par.eta0)*((ns->ed.stsp.par.alpha)*E2[i][j] + E2CT[i][j]);
                    }
                }
            }

            /*here
            if ((ns->ed.stsp.contr.model == GW)){
                //For GW model
                real E2[DIM][DIM];
                for (int i = 0; i < DIM; i++) {
                    for (int j = 0; j < DIM; j++) {
                        E2[i][j] = 0.5*(Du[i][j]+Du[j][i]);
                    }
                }
                //DEBUG_INSPECT(E2[0][1],%lf);
                //DEBUG_INSPECT(E2[1][0],%lf);
                for (int i = 0; i < DIM; i++) {
                    for (int j = 0; j < DIM; j++) {
                       tau_M[i][j] = (ns->ed.stsp.par.alpha)*(ns->ed.stsp.par.eta0)*E2[i][j];
                    }
                }
                //DEBUG_INSPECT(tau_M[0][1],%lf);
                //DEBUG_INSPECT(tau_M[1][0],%lf);
            } else {
                //for GW-WC model
                real E2[DIM][DIM];
                for (int i = 0; i < DIM; i++) {
                    for (int j = 0; j < DIM; j++) {
                        E2[i][j] = 0.5*(Du[i][j]+Du[j][i]);
                    }
                }
                //DEBUG_INSPECT(E2[0][1],%lf);
                //DEBUG_INSPECT(E2[1][0],%lf);
                //Calculation of chi value
                real chi = hig_flow_calculate_jamming_coordinate_chi (A, E_C);
                if (chi > chimax) chimax = chi;
                if (chi < chimin) chimin = chi;
                //DEBUG_INSPECT(chi,%lf);
                //chi = fabs(chi);
                Pis = fabs(Pis);
                //real fPi;
                real fPi = exp(-Pic/(Pis+small2));
                //if (Pis <= small2) {
                //    fPi = exp(-Pic/(Pis+small2));
                //} else {
                //    fPi = exp(-Pic/Pis);
                //}
                //Calculation of function f(Pi)
                //real fPi = exp(-Pic/(Pis+small2));
                //Check to see min and max values of f
                if (fPi > fmax) fmax = fPi;
                if (fPi < fmin) fmin = fPi;
                //DEBUG_INSPECT(fPi,%lf);
                //Calculation of the jamming point Chi_J
                real chi_J = (1.0-fPi)*chij1 + fPi*chij2;
                real Xchi = ns->ed.stsp.get_X(ccenter, ns->par.t, X0, chi, chi_J);
                //Check to see min and max values of Xchi
                if (Xchi > chiJmax) chiJmax = Xchi;
                if (Xchi < chiJmin) chiJmin = Xchi;
                //DEBUG_INSPECT(Xchi,%lf);
                Xchi = fabs(Xchi);
                real Xchi2;
                if (Xchi < small2) {
                    Xchi2 = small2;
                } else {
                    Xchi2 = Xchi;
                }

                //DEBUG_INSPECT(E_C[1][0],%lf);
                //DEBUG_INSPECT(Xchi,%lf);
                //DEBUG_INSPECT(E2CT[0][1],%lf);
                //for (int i = 0; i < DIM; i++) {ñ
                //    for (int j = 0; j < DIM; j++) {
                //        tau_M[i][j] = 0.0;
                //    }
                //}

                for (int i = 0; i < DIM; i++) {
                    for (int j = 0; j < DIM; j++) {
                        tau_M[i][j] = 0.0;
                        //tau_M[i][j] = eta0*(alpha*E2[i][j] + Xchi*E_C[i][j]);
                        //tau_M[i][j] += eta0*alpha*E2[i][j];
                        //tau_M[i][j] += eta0*Xchi2*E_C[i][j];
                        tau_M[i][j] += eta0*alpha*E2[i][j];
                        tau_M[i][j] += eta0*Xchi2*E_C[i][j];
                       //tau_M[i][j] = (ns->ed.stsp.par.eta0)*(ns->ed.stsp.par.alpha)*E2[i][j];
                       //tau_M[i][j] += (ns->ed.stsp.par.eta0)*(small + Xchi)*E_C[i][j];
                       //tau_M[i][j] = (ns->ed.stsp.par.eta0)*((ns->ed.stsp.par.alpha)*E2[i][j] + E2CT[i][j]);
                    }
                }
            }
            */
            //DEBUG_INSPECT(tau_M[0][0],%lf);
            //DEBUG_INSPECT(tau_M[1][1],%lf);

            // DEBUG_INSPECT(tau_M[0][0],%lf);
            //  Calculate the suspension stress tensor
            real T[DIM][DIM], S[DIM][DIM];
            hig_flow_calculate_suspension_stress_tensor(tau_M, FOM, T);
            // DEBUG_INSPECT(T[0][1],%lf);
            // DEBUG_INSPECT(T[1][0],%lf);
            //  Calculate the max and min values of T
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    if (T[i][j] > Smax[i][j])
                        Smax[i][j] = T[i][j];
                    if (T[i][j] < Smin[i][j])
                        Smin[i][j] = T[i][j];
                    S[i][j] = T[i][j] - 2.0 * (ns->ed.stsp.par.eta0) * 0.5 * (Du[i][j] + Du[j][i]);
                    // S[i][j] = T[i][j]-2.0*(ns->ed.stsp.par.eta0)*E[i][j];
                    //  Store the Polymeric Tensor
                    dp_set_value(ns->ed.stsp.dpS[i][j], clid, S[i][j]);
                }
            }
            // real S[DIM][DIM];
            //  Calculate the tensor S
            /*
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    //if (T[i][j] > Smax[i][j]) Smax[i][j] = T[i][j];
                    //if (T[i][j] < Smin[i][j]) Smin[i][j] = T[i][j];
                    S[i][j] = T[i][j]-2.0*(ns->ed.stsp.par.eta0)*E[i][j];
                    // Store the Polymeric Tensor
                    dp_set_value(ns->ed.stsp.dpS[i][j], clid, S[i][j]);
                }
            } */
            // DEBUG_INSPECT(Du[0][1],%lf);
            // DEBUG_INSPECT(Du[1][0],%lf);
            //  Calculate the tensor S
            // for (int i = 0; i < DIM; i++) {
            //     for (int j = 0; j < DIM; j++) {
            // if (T[i][j] > Smax[i][j]) Smax[i][j] = T[i][j];
            // if (T[i][j] < Smin[i][j]) Smin[i][j] = T[i][j];
            // S[i][j] = T[i][j]-2.0*eta0*E[i][j];
            //        S[i][j] = T[i][j]-2.0*eta0*0.5*(Du[i][j]+Du[j][i]);
            // Store the Polymeric Tensor
            //        dp_set_value(ns->ed.stsp.dpS[i][j], clid, S[i][j]);
            // if (i != j) {
            //     dp_set_value(ns->ed.stsp.dpS[j][i], clid, S[i][j]);
            //}
            //    }
            //}
        }
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                // Printing the min and max tensor
                printf("===> %d %d: Tmin = %lf <===> Tmax = %lf <===\n", i, j, Smin[i][j], Smax[i][j]);
            }
        }
        // for (int i = 0; i < DIM; i++) {
        //    for (int j = 0; j < DIM; j++) {
        //  Printing the min and max tensor
        //       printf("===> %d %d: E_Emin = %lf <===> E_Emax = %lf <===\n",i,j,EEmin[i][j],EEmax[i][j]);
        // printf("===> %d %d: E_Cmin = %lf <===> E_Cmax = %lf <===\n",i,j,ECmin[i][j],ECmax[i][j]);
        //   }
        //}
        // for (int i = 0; i < DIM; i++) {
        //   for (int j = 0; j < DIM; j++) {
        // Printing the min and max tensor
        // printf("===> %d %d: E_Emin = %lf <===> E_Emax = %lf <===\n",i,j,EEmin[i][j],EEmax[i][j]);
        //       printf("===> %d %d: E_Cmin = %lf <===> E_Cmax = %lf <===\n",i,j,ECmin[i][j],ECmax[i][j]);
        //   }
        //}

        if ((ns->ed.stsp.contr.model == GW_WC))
        {
            // Printing the min and max values of f
            printf("===> fmin = %lf <===> fmax = %lf <===\n", fmin, fmax);
            // Printing the min and max values of chi
            printf("===> chimin = %lf <===> chimax = %lf <===\n", chimin, chimax);
            // Printing the min and max values of Xchi
            printf("===> Xchimin = %lf <===> Xchimax = %lf <===\n", chiJmin, chiJmax);
        }
        if ((ns->ed.stsp.contr.model == GW_WC_IF))
        {
            // Printing the min and max values of f
            printf("===> fmin = %lf <===> fmax = %lf <===\n", fmin, fmax);
            // Printing the min and max values of chi
            printf("===> chimin = %lf <===> chimax = %lf <===\n", chimin, chimax);
            // Printing the min and max values of Xchi
            printf("===> Xchimin = %lf <===> Xchimax = %lf <===\n", chiJmin, chiJmax);
            // Printing the min and max values of alpha
            printf("===> alphamin = %lf <===> alphamax = %lf <===\n", alphamin, alphamax);
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed pressure property
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                dp_sync(ns->ed.stsp.dpS[i][j]);
            }
        }
    }
}


// Computing the particle stress TensorTau
void higflow_compute_particle_stress_tensor_shear_thickening_suspension(higflow_solver *ns)
{
    if ((ns->contr.flowtype == 9))
    {
        // Get the constants
        // real Re   = ns->par.Re;
        real alpha = ns->ed.stsp.par.alpha;
        real eta0 = ns->ed.stsp.par.eta0;
        real chij1 = ns->ed.stsp.par.chij1;
        real chij2 = ns->ed.stsp.par.chij2;
        real X0 = ns->ed.stsp.par.X0;
        real Pic = ns->ed.stsp.par.Pic;
        real phircp = ns->ed.stsp.par.phi;
        real fmax = -1.0e16;
        real fmin = 1.0e16;
        real chimax = -1.0e16;
        real chimin = 1.0e16;
        real chiJmax = -1.0e16;
        real chiJmin = 1.0e16;
        real alphamax = -1.0e16;
        real alphamin = 1.0e16;
        real small = 1.0e-10;
        real small2 = 1.0e-6;
        real Smax[DIM][DIM], Smin[DIM][DIM];
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                Smax[i][j] = -1.0e16;
                Smin[i][j] = 1.0e16;
            }
        }
        real EEmax[DIM][DIM], EEmin[DIM][DIM];
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                EEmax[i][j] = -1.0e16;
                EEmin[i][j] = 1.0e16;
            }
        }
        real ECmax[DIM][DIM], ECmin[DIM][DIM];
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                ECmax[i][j] = -1.0e16;
                ECmin[i][j] = 1.0e16;
            }
        }
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdp);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid = mp_lookup(mp, hig_get_cid(c));
            // Get the inside/outside inflow point cell
            int inflowcell, outflowcell;
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Get the velocity derivative tensor Du and the microstructure tensor
            real A[DIM][DIM], Du[DIM][DIM], S0[DIM][DIM];
            // AM is the Trace of the microstructure tensor
            real AM = 0.0;
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    // Get Du
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
                    // Get the microstructure tensor
                    A[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpA[i][j], ns->ed.stn);
                    // Get the previous value of polymeric tensor, S0
                    S0[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpS[i][j], ns->ed.stn);
                }
                AM += A[i][i];
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
            // DEBUG_INSPECT(Du[0][1],%lf);
            // DEBUG_INSPECT(E[1][0],%lf);
            // Calculate T of previous step, T0
            real T0[DIM][DIM];
            real Pis = 0.0;
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    T0[i][j] = S0[i][j] + 2.0 * eta0 * E[i][j];
                }
                Pis += T0[i][i];
            }
            Pis = -1.0 * Pis / 3.0;
            // Pis = -1.0/3.0;
            // DEBUG_INSPECT(Pis,%lf);
            //  Eigen-values and eigen-vectors of the deformation tensor E
            real R[DIM][DIM], lambda[DIM];
            // DEBUG_INSPECT(E[0][1],%lf);
            // DEBUG_INSPECT(E[1][0],%lf);
            hig_flow_jacobi(E, lambda, R);
            // DEBUG_INSPECT(lambda[0],%lf);

            // DEBUG_INSPECT(lambda[1],%lf);
            // Calculation of the compressive rate of strain Ec and the extensional rate of strain Ee
            real E_C[DIM][DIM], E_E[DIM][DIM];
            hig_flow_compressive_extensional_tensors(lambda, R, E_C, E_E);
            // Kronecker's delta function
            real KD[DIM][DIM];
            // real KM;
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    if (i == j)
                    {
                        KD[i][i] = 1.0;
                    }
                    else
                    {
                        KD[i][j] = 0.0;
                    }
                    // printf("%lf     ", KD[i][j]);
                }
                // KM += KD[i][i];
                // printf("\n")
            }

            //******************************************************************************************+

            // Calculation of the fourth order orientation moment <pppp>
            real FOM[DIM][DIM][DIM][DIM];
            hig_flow_calculate_4th_order_orientation_moment(AM, KD, A, FOM);
            real tau_M[DIM][DIM];
            if (ns->ed.stsp.contr.model == GW){
                // For GW model
                real E2[DIM][DIM];
                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        E2[i][j] = 0.5 * (Du[i][j] + Du[j][i]);
                    }
                }
                // DEBUG_INSPECT(E2[0][1],%lf);
                // DEBUG_INSPECT(E2[1][0],%lf);
                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        tau_M[i][j] = 2.0*(ns->ed.stsp.par.alpha) * (ns->ed.stsp.par.eta0) * E2[i][j];
                    }
                }
                // DEBUG_INSPECT(tau_M[0][1],%lf);
                // DEBUG_INSPECT(tau_M[1][0],%lf);
            } else if (ns->ed.stsp.contr.model == GW_WC_IF) {
                // GW-WC model for inhomogeneous flows
                real E2[DIM][DIM];
                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        E2[i][j] = 0.5 * (Du[i][j] + Du[j][i]);
                    }
                }
                //DEBUG_INSPECT(E2[0][1],%lf);
                //DEBUG_INSPECT(E2[1][0],%lf);
                // Calculation of chi value
                real chi = hig_flow_calculate_jamming_coordinate_chi(A, E_C);
                if (chi > chimax)
                    chimax = chi;
                if (chi < chimin)
                    chimin = chi;
                // DEBUG_INSPECT(chi,%lf);
                // chi = fabs(chi);
                Pis = fabs(Pis);
                // real fPi;
                real fPi = exp(-Pic / (Pis + small2));
                // if (Pis <= small2) {
                //     fPi = exp(-Pic/(Pis+small2));
                // } else {
                //     fPi = exp(-Pic/Pis);
                // }
                // Calculation of function f(Pi)
                // real fPi = exp(-Pic/(Pis+small2));
                // Check to see min and max values of f
                if (fPi > fmax)
                    fmax = fPi;
                if (fPi < fmin)
                    fmin = fPi;
                // DEBUG_INSPECT(fPi,%lf);
                // Calculation of the jamming point Chi_J
                real chi_J = (1.0 - fPi) * chij1 + fPi * chij2;
                /*
                real tol_chi = 0.1;
                real chi2chij = chi/chi_J;
                real tol2chi = 1.0-chi2chij;
                //real tol_chi = 0.977*chi_J;
                //if (chi >= phircp - tol_chi){
                //    chi = phircp - tol_chi;
                //}
                //if (chi >= tol_chi){
                //    chi = tol_chi;
                //}
                if (tol2chi <= tol_chi){
                    chi = (1.0-tol_chi)*chi_J;
                }*/
                //Tolerance value that works
                
                real tol_chi = 0.85*chi_J;
                //real tol_chi = 0.86*chi_J;
                if (chi >= tol_chi){
                    chi = tol_chi;
                }
                real Xchi = ns->ed.stsp.get_X(ccenter, ns->par.t, X0, chi, chi_J);
                // Check to see min and max values of Xchi
                if (Xchi > chiJmax)
                    chiJmax = Xchi;
                if (Xchi < chiJmin)
                    chiJmin = Xchi;
                // DEBUG_INSPECT(Xchi,%lf);
                Xchi = fabs(Xchi);
                real Xchi2;
                if (Xchi < small2)
                {
                    Xchi2 = small2;
                }
                else
                {
                    Xchi2 = Xchi;
                }

                // DEBUG_INSPECT(E_C[1][0],%lf);
                // DEBUG_INSPECT(Xchi,%lf);
                // DEBUG_INSPECT(E2CT[0][1],%lf);
                // for (int i = 0; i < DIM; i++) {ñ
                //     for (int j = 0; j < DIM; j++) {
                //         tau_M[i][j] = 0.0;
                //     }
                // }

                // Get the volume fraction value
                real varphi = compute_value_at_point(ns->ed.stsp.sdphi, ccenter, ccenter, 1.0, ns->ed.stsp.dpphi, ns->ed.stn);
                real frac_tol = ns->ed.stsp.par.gdrms;
                real eps = 1.0e-4;
                if(varphi > 1.0){
                    varphi = 1.0-eps;
                }
                //if(varphi > phircp - frac_tol){
                //    varphi = phircp - frac_tol;
                //}
                if(varphi < eps){
                    varphi = 0.0;
                }
                //DEBUG_INSPECT(varphi,%lf);
                real alphaVF = ns->ed.stsp.get_alpha(ccenter, ns->par.t, alpha, varphi, phircp);
                //DEBUG_INSPECT(alphaVF,%lf);
                if (alphaVF > alphamax)
                    alphamax = alphaVF;
                if (alphaVF < alphamin)
                    alphamin = alphaVF;
                alphaVF = fabs(alphaVF);
                real alphaVF2;
                if (alphaVF < small2)
                {
                    alphaVF2 = small2;
                }
                else
                {
                    alphaVF2 = alphaVF;
                }


                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        tau_M[i][j] = 0.0;
                        // tau_M[i][j] = eta0*(alpha*E2[i][j] + Xchi*E_C[i][j]);
                        // tau_M[i][j] += eta0*alpha*E2[i][j];
                        // tau_M[i][j] += eta0*Xchi2*E_C[i][j];
                        tau_M[i][j] += eta0 * alphaVF2 * E2[i][j];
                        tau_M[i][j] += eta0 * Xchi2 * E_C[i][j];
                        //tau_M[i][j] *= 2.0;
                        // tau_M[i][j] = (ns->ed.stsp.par.eta0)*(ns->ed.stsp.par.alpha)*E2[i][j];
                        // tau_M[i][j] += (ns->ed.stsp.par.eta0)*(small + Xchi)*E_C[i][j];
                        // tau_M[i][j] = (ns->ed.stsp.par.eta0)*((ns->ed.stsp.par.alpha)*E2[i][j] + E2CT[i][j]);
                    }
                }
            } else {
                // for GW-WC model
                real E2[DIM][DIM];
                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        E2[i][j] = 0.5 * (Du[i][j] + Du[j][i]);
                    }
                }
                // DEBUG_INSPECT(E2[0][1],%lf);
                // DEBUG_INSPECT(E2[1][0],%lf);
                // Calculation of chi value
                real chi = hig_flow_calculate_jamming_coordinate_chi(A, E_C);
                if (chi > chimax)
                    chimax = chi;
                if (chi < chimin)
                    chimin = chi;
                // DEBUG_INSPECT(chi,%lf);
                // chi = fabs(chi);
                Pis = fabs(Pis);
                // real fPi;
                real fPi = exp(-Pic / (Pis + small2));
                // if (Pis <= small2) {
                //     fPi = exp(-Pic/(Pis+small2));
                // } else {
                //     fPi = exp(-Pic/Pis);
                // }
                // Calculation of function f(Pi)
                // real fPi = exp(-Pic/(Pis+small2));
                // Check to see min and max values of f
                if (fPi > fmax)
                    fmax = fPi;
                if (fPi < fmin)
                    fmin = fPi;
                // DEBUG_INSPECT(fPi,%lf);
                // Calculation of the jamming point Chi_J
                real chi_J = (1.0 - fPi) * chij1 + fPi * chij2;
                real tol_chi = chi_J - 0.01;
                //real tol_chi = 0.88*chi_J;
                if (chi >= tol_chi){
                    chi = tol_chi;
                }
                real Xchi = ns->ed.stsp.get_X(ccenter, ns->par.t, X0, chi, chi_J);
                // Check to see min and max values of Xchi
                if (Xchi > chiJmax)
                    chiJmax = Xchi;
                if (Xchi < chiJmin)
                    chiJmin = Xchi;
                // DEBUG_INSPECT(Xchi,%lf);
                Xchi = fabs(Xchi);
                real Xchi2;
                if (Xchi < small2)
                {
                    Xchi2 = small2;
                }
                else
                {
                    Xchi2 = Xchi;
                }

                // DEBUG_INSPECT(E_C[1][0],%lf);
                // DEBUG_INSPECT(Xchi,%lf);
                // DEBUG_INSPECT(E2CT[0][1],%lf);
                // for (int i = 0; i < DIM; i++) {ñ
                //     for (int j = 0; j < DIM; j++) {
                //         tau_M[i][j] = 0.0;
                //     }
                // }

                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        tau_M[i][j] = 0.0;
                        // tau_M[i][j] = eta0*(alpha*E2[i][j] + Xchi*E_C[i][j]);
                        // tau_M[i][j] += eta0*alpha*E2[i][j];
                        // tau_M[i][j] += eta0*Xchi2*E_C[i][j];
                        tau_M[i][j] += eta0 * alpha * E2[i][j];
                        tau_M[i][j] += eta0 * Xchi2 * E_C[i][j];
                        // tau_M[i][j] = (ns->ed.stsp.par.eta0)*(ns->ed.stsp.par.alpha)*E2[i][j];
                        // tau_M[i][j] += (ns->ed.stsp.par.eta0)*(small + Xchi)*E_C[i][j];
                        // tau_M[i][j] = (ns->ed.stsp.par.eta0)*((ns->ed.stsp.par.alpha)*E2[i][j] + E2CT[i][j]);
                    }
                }
            }


            // DEBUG_INSPECT(tau_M[0][0],%lf);
            //  Calculate the suspension stress tensor
            real T[DIM][DIM], S[DIM][DIM];
            hig_flow_calculate_suspension_stress_tensor(tau_M, FOM, T);
            // DEBUG_INSPECT(T[0][1],%lf);
            // DEBUG_INSPECT(T[1][0],%lf);
            //  Calculate the max and min values of T
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    if (T[i][j] > Smax[i][j])
                        Smax[i][j] = T[i][j];
                    if (T[i][j] < Smin[i][j])
                        Smin[i][j] = T[i][j];
                    S[i][j] = T[i][j];
                    // S[i][j] = T[i][j]-2.0*(ns->ed.stsp.par.eta0)*E[i][j];
                    //  Store the Polymeric Tensor
                    dp_set_value(ns->ed.stsp.dpS[i][j], clid, S[i][j]);
                }
            }
            // real S[DIM][DIM];
            //  Calculate the tensor S
            /*
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    //if (T[i][j] > Smax[i][j]) Smax[i][j] = T[i][j];
                    //if (T[i][j] < Smin[i][j]) Smin[i][j] = T[i][j];
                    S[i][j] = T[i][j]-2.0*(ns->ed.stsp.par.eta0)*E[i][j];
                    // Store the Polymeric Tensor
                    dp_set_value(ns->ed.stsp.dpS[i][j], clid, S[i][j]);
                }
            } */
            // DEBUG_INSPECT(Du[0][1],%lf);
            // DEBUG_INSPECT(Du[1][0],%lf);
            //  Calculate the tensor S
            // for (int i = 0; i < DIM; i++) {
            //     for (int j = 0; j < DIM; j++) {
            // if (T[i][j] > Smax[i][j]) Smax[i][j] = T[i][j];
            // if (T[i][j] < Smin[i][j]) Smin[i][j] = T[i][j];
            // S[i][j] = T[i][j]-2.0*eta0*E[i][j];
            //        S[i][j] = T[i][j]-2.0*eta0*0.5*(Du[i][j]+Du[j][i]);
            // Store the Polymeric Tensor
            //        dp_set_value(ns->ed.stsp.dpS[i][j], clid, S[i][j]);
            // if (i != j) {
            //     dp_set_value(ns->ed.stsp.dpS[j][i], clid, S[i][j]);
            //}
            //    }
            //}
        }
        /*
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                // Printing the min and max tensor
                printf("===> %d %d: Tmin = %lf <===> Tmax = %lf <===\n", i, j, Smin[i][j], Smax[i][j]);
            }
        }
        // for (int i = 0; i < DIM; i++) {
        //    for (int j = 0; j < DIM; j++) {
        //  Printing the min and max tensor
        //       printf("===> %d %d: E_Emin = %lf <===> E_Emax = %lf <===\n",i,j,EEmin[i][j],EEmax[i][j]);
        // printf("===> %d %d: E_Cmin = %lf <===> E_Cmax = %lf <===\n",i,j,ECmin[i][j],ECmax[i][j]);
        //   }
        //}
        // for (int i = 0; i < DIM; i++) {
        //   for (int j = 0; j < DIM; j++) {
        // Printing the min and max tensor
        // printf("===> %d %d: E_Emin = %lf <===> E_Emax = %lf <===\n",i,j,EEmin[i][j],EEmax[i][j]);
        //       printf("===> %d %d: E_Cmin = %lf <===> E_Cmax = %lf <===\n",i,j,ECmin[i][j],ECmax[i][j]);
        //   }
        //}

        if ((ns->ed.stsp.contr.model == GW_WC))
        {
            // Printing the min and max values of f
            printf("===> fmin = %lf <===> fmax = %lf <===\n", fmin, fmax);
            // Printing the min and max values of chi
            printf("===> chimin = %lf <===> chimax = %lf <===\n", chimin, chimax);
            // Printing the min and max values of Xchi
            printf("===> Xchimin = %lf <===> Xchimax = %lf <===\n", chiJmin, chiJmax);
        }
        if ((ns->ed.stsp.contr.model == GW_WC_IF))
        {
            // Printing the min and max values of f
            printf("===> fmin = %lf <===> fmax = %lf <===\n", fmin, fmax);
            // Printing the min and max values of chi
            printf("===> chimin = %lf <===> chimax = %lf <===\n", chimin, chimax);
            // Printing the min and max values of Xchi
            printf("===> Xchimin = %lf <===> Xchimax = %lf <===\n", chiJmin, chiJmax);
            // Printing the min and max values of alpha
            printf("===> alphamin = %lf <===> alphamax = %lf <===\n", alphamin, alphamax);
        }
        */
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed pressure property
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                dp_sync(ns->ed.stsp.dpS[i][j]);
            }
        }
    }
}



// *******************************************************************
// Constitutive Equation Step for the Explicit Euler Method
// *******************************************************************
// Solve the evolution equation of the microstructure tensor <pp>
void higflow_explicit_euler_evolution_equation_microstructure_tensor(higflow_solver *ns)
{
    if ((ns->contr.flowtype == 9))
    {
        // Get the constants
        real beta = ns->ed.stsp.par.beta;
        real small = 1.0e-14;
        real phi;
        if (ns->ed.stsp.contr.model == GW)
        {
            // Value of phi for the GW model
            phi = 1.0;
        }
        else
        {
            // Value of phi for the GW-WC model
            phi = ns->ed.stsp.par.phi;
        }
        // standard deviation of the shear rate fluctuations
        real gdrms = ns->ed.stsp.par.gdrms;

        real Amax[DIM][DIM], Amin[DIM][DIM];
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                Amax[i][j] = -1.0e16;
                Amin[i][j] = 1.0e16;
            }
        }
        real RHSmax[DIM][DIM], RHSmin[DIM][DIM];
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                RHSmax[i][j] = -1.0e16;
                RHSmin[i][j] = 1.0e16;
            }
        }
        // tol        = 1.0e-3;
        //  Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for (int i = 0; i < DIM; i++)
        {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdp);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid = mp_lookup(mp, hig_get_cid(c));
            // Get the center of the cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Get the velocity derivative tensor Du and the microstructure tensor A
            real Du[DIM][DIM], A[DIM][DIM];
            // AM is the trace of the microstructure tensor
            real AM = 0.0;
            // Get the tensors
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    // Get Du
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
                    // Get microstructure tensor A =  <pp>
                    A[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpA[i][j], ns->ed.stn);
                }
                AM += A[i][i];
            }

            // Calculate the rate of strain tensor E
            real E[DIM][DIM], L[DIM][DIM];
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    L[i][j] = Du[j][i];
                    E[i][j] = 0.5 * (Du[i][j] + Du[j][i]);
                }
            }

            // Eigen-values and eigen-vectors of the deformation tensor E
            real R[DIM][DIM], lambda[DIM];
            hig_flow_jacobi(E, lambda, R);
            // Calculation of the compressive rate of strain Ec and the extensional rate of strain Ee
            real E_C[DIM][DIM], E_E[DIM][DIM];
            hig_flow_compressive_extensional_tensors(lambda, R, E_C, E_E);
            // Kronecker's delta function
            real KD[DIM][DIM];
            // real KM;
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    if (i == j)
                    {
                        KD[i][i] = 1.0;
                    }
                    else
                    {
                        KD[i][j] = 0.0;
                    }
                    // printf("%lf     ", KD[i][j]);
                }
                // KM += KD[i][i];
                // printf("\n")
            }
            // Calculation of the fourth order orientation moment <pppp>
            real FOM[DIM][DIM][DIM][DIM];
            hig_flow_calculate_4th_order_orientation_moment(AM, KD, A, FOM);
            // Calculate the trace of tensor E_C
            real TrEC = 0.0;
            for (int i = 0; i < DIM; i++)
            {
                TrEC += E_C[i][i];
            }
            real LP[DIM][DIM], EP[DIM][DIM];
            // Calculation of the term (Du)^T : <pppp>
            hig_flow_calculate_double_dot_product_tensors(L, FOM, LP);
            // Calculation of the E_E tensor with fluctuations
            real E_ERMS[DIM][DIM], E_CRMS[DIM][DIM];
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    E_ERMS[i][j] = E_E[i][j] + 0.25 * gdrms * KD[i][j];
                    //E_CRMS[i][j] = E_C[i][j] + 0.25 * gdrms * KD[i][j];
                    E_CRMS[i][j] = E_C[i][j] - 0.25 * gdrms * KD[i][j];
                }
            }
            //DEBUG_INSPECT(E_CRMS[1][1],%lf);
            //DEBUG_INSPECT(E_C[1][1],%lf);
            // Calculate the trace of tensor E_C
            real TrECRMS = 0.0;
            for (int i = 0; i < DIM; i++)
            {
                TrECRMS += E_CRMS[i][i];
            }
            //DEBUG_INSPECT(TrECRMS,%lf);
            //DEBUG_INSPECT(TrEC,%lf);
            // Calculation of the term E_E : <pppp>
            // hig_flow_calculate_double_dot_product_tensors (E_E, FOM, EP);
            //  Calculate RHS = A*Du + (Du)^T * A- 2*(Du)^T : <pppp> - beta*( E_E:<pppp> + (phi/15)*(2*E_C+Tr(E_C)*KD) )
            real RHS[DIM][DIM];
            if (ns->ed.stsp.contr.model == GW_WC_IF)
            {
                // Model with fluctuations
                hig_flow_calculate_double_dot_product_tensors(E_ERMS, FOM, EP);
                real varphi = compute_value_at_point(ns->ed.stsp.sdphi, ccenter, ccenter, 1.0, ns->ed.stsp.dpphi, ns->ed.stn);
                real frac_tol = ns->ed.stsp.par.gdrms;
                real eps = 1.0e-4;
                if(varphi > 1.0){
                    varphi = 1.0-eps;
                }
                //if(varphi > phi- frac_tol){
                //    varphi = phi- frac_tol;
                //}
                if(varphi < eps){
                    varphi = 0.0;
                }
                //DEBUG_INSPECT(varphi,%lf);
                hig_flow_evolution_equation_microstructure_tensor_rhs(beta, varphi, A, L, E_CRMS, TrECRMS, KD, FOM, LP, EP, RHS);
            }
            else
            {
                hig_flow_calculate_double_dot_product_tensors(E_E, FOM, EP);
                // hig_flow_evolution_equation_microstructure_tensor_rhs (beta, phi, A, Du, E_C, TrEC, KD, FOM, LP, EP, RHS);
                hig_flow_evolution_equation_microstructure_tensor_rhs(beta, phi, A, L, E_C, TrEC, KD, FOM, LP, EP, RHS);
            }
            // Get the velocity at cell center
            real u[DIM], dAdx[DIM], d2Adx2[DIM];
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            // Solving the Constitutive Equation using the Euler Method
            for (int i = 0; i < DIM; i++)
            {
                for (int j = i; j < DIM; j++)
                {
                    // Right hand side equation
                    real rhs = 0.0;
                    switch (ns->ed.stsp.contr.convecdiscrtype)
                    {
                    case CELL_CENTRAL:
                        // Tensor derivative at cell center
                        hig_flow_derivative_tensor_A_at_center_cell(ns, ccenter, cdelta, i, j, A[i][j], dAdx);
                        for (int dim = 0; dim < DIM; dim++)
                        {
                            // Compute convective tensor term in rhs
                            rhs -= u[dim] * dAdx[dim];
                        }
                        break;
                    case CELL_CUBISTA:
                        // Compute convective tensor term CUBISTA in rhs
                        for (int dim = 0; dim < DIM; dim++)
                        {
                            rhs -= hig_flow_convective_tensor_A_term_cubista(ns, ns->dpu[dim], ns->ed.sdED, ns->ed.stn, A, ccenter, cdelta, dim, i, j);
                            //rhs -= 0.0;
                        }
                        break;
                    }
                    // if (rhs > RHSmax[i][j]) RHSmax[i][j] = rhs;
                    // if (rhs < RHSmin[i][j]) RHSmin[i][j] = rhs;
                    //  Compute the final rhs
                    rhs += RHS[i][j];
                    // Calculate RHS max and min
                    //  Compute A at next time
                    real Anew = A[i][j] + ns->par.dt * rhs;
                    // Store A in S
                    dp_set_value(ns->ed.stsp.dpA[i][j], clid, Anew);
                    if (i != j)
                    {
                        dp_set_value(ns->ed.stsp.dpA[j][i], clid, Anew);
                    }
                    // Calculate max and min A
                    if (Anew > Amax[i][j])
                        Amax[i][j] = Anew;
                    if (Anew < Amin[i][j])
                        Amin[i][j] = Anew;
                }
            }
        }

        // Destroy the iterator
        // higcit_destroy(it);
        // Sync the ditributed pressure property
        // for (int i = 0; i < DIM; i++) {
        //    for (int j = 0; j < DIM; j++) {
        //        dp_sync(ns->ed.vesb.dpS[i][j]);
        //    }
        //}

        // Store the Kernel Tensor
        // for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
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

        // for (int i = 0; i < DIM; i++) {
        //    for (int j = 0; j < DIM; j++) {
        // Printing the min and max tensor
        //       printf("===> %d %d: Amin = %lf <===> Amax = %lf <===\n",i,j, Amin[i][j], Amax[i][j]);
        //   }
        //}

        // Printing the min and max deformation tensor values
        // for (int i = 0; i < DIM; i++) {
        // for (int j = 0; j < DIM; j++) {
        //  Printing the min and max tensor
        // printf("===> %d %d: Dmin = %lf <===> Dmax = %lf <===\n",i,j,Dmin[i][j],Dmax[i][j]);
        //    }
        //}

        // Printing the min and max deformation tensor values
        // for (int i = 0; i < DIM; i++) {
        // for (int j = 0; j < DIM; j++) {
        //  Printing the min and max tensor
        // printf("===> %d %d: ANAmin = %lf <===> ANAmax = %lf <===\n",i,j,ANAmin[i][j],ANAmax[i][j]);
        //}
        //}

        // for (int i = 0; i < DIM; i++) {
        // for (int j = 0; j < DIM; j++) {
        //  Printing the min and max tensor
        // printf("===> %d %d: RHSAmin = %lf <===> RHSAmax = %lf <===\n",i,j, RHSmin[i][j], RHSmax[i][j]);
        //}
        //}

        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed pressure property
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                dp_sync(ns->ed.stsp.dpA[i][j]);
            }
        }
    }
}

// *******************************************************************
// Calculate convective tensor term CUBISTA
// *******************************************************************
real hig_flow_convective_tensor_A_term_cubista(higflow_solver *ns, distributed_property *dpu, sim_domain *sdp, sim_stencil *stn, real K[DIM][DIM], Point ccenter, Point cdelta, int dim, int i, int j)
{
    real vbar[DIM], dKdx[dim], kr, krr, kl, kll, kc, a, b, c, d, e, tol, fi, conv1, conv2;
    a = 1.7500;
    b = 0.3750;
    c = 0.7500;
    d = 0.1250;
    e = 0.2500;
    tol = 1.0e-14;
    conv1 = 0.0;
    conv2 = 0.0;
    int incell_r, incell_l, incell_ll, incell_rr, infacet;
    // Get the tensor at center cell
    kc = K[i][j];
    // Get the low, high, lowlow, highhigh component kernel at center cell
    kl = compute_center_p_left_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.stsp.dpA[i][j], ns->ed.stn, &incell_l);
    kr = compute_center_p_right_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.stsp.dpA[i][j], ns->ed.stn, &incell_r);
    kll = compute_center_p_left_22(ns->ed.sdED, ccenter, cdelta, dim, 2.0, ns->ed.stsp.dpA[i][j], ns->ed.stn, &incell_ll);
    krr = compute_center_p_right_22(ns->ed.sdED, ccenter, cdelta, dim, 2.0, ns->ed.stsp.dpA[i][j], ns->ed.stn, &incell_rr);
    // Get the velocity  v1bar(i+1/2,j) in the facet center
    vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
    if (vbar[dim] > 0.0)
    {
        if (fabs(kr - kl) <= tol)
        {
            conv1 = vbar[dim] * kc;
        }
        else
        {
            fi = (kc - kl) / (kr - kl);
            if ((fi <= 0.0) || (fi >= 1.0))
            {
                conv1 = vbar[dim] * kc;
            }
            else
            {
                if (fi < b)
                {
                    if (incell_l == 1)
                        conv1 = vbar[dim] * (a * kc - c * kl);
                    else
                        conv1 = vbar[dim] * kc;
                }
                if ((fi >= b) && (fi <= c))
                {
                    if ((incell_l == 1) && (incell_r == 1))
                        conv1 = vbar[dim] * (c * kc + b * kr - d * kl);
                    else
                        conv1 = vbar[dim] * kc;
                }
                if (fi > c)
                {
                    if (incell_r == 1)
                        conv1 = vbar[dim] * (e * kc + c * kr);
                    else
                        conv1 = vbar[dim] * kc;
                }
            }
        }
        // v1bar < 0.0
    }
    else
    {
        if ((incell_r == 1) && (incell_rr == 1))
        {
            if (fabs(kc - krr) <= tol)
            {
                conv1 = vbar[dim] * kr;
            }
            else
            {
                fi = (kr - krr) / (kc - krr);
                if ((fi <= 0.0) || (fi >= 1.0))
                {
                    conv1 = vbar[dim] * kr;
                }
                else
                {
                    if (fi < b)
                        conv1 = vbar[dim] * (a * kr - c * krr);
                    if ((fi >= b) && (fi <= c))
                        conv1 = vbar[dim] * (c * kr + b * kc - d * krr);
                    if (fi > c)
                        conv1 = vbar[dim] * (c * kc + e * kr);
                }
            }
            // Return upwind value at boundary
        }
        else if ((incell_r == 1) && (incell_rr == 0))
        {
            if (fabs(kc - krr) <= tol)
            {
                conv1 = vbar[dim] * kr;
            }
            else
            {
                fi = (kr - krr) / (kc - krr);
                if ((fi <= 0.0) || (fi >= 1.0))
                {
                    conv1 = vbar[dim] * kr;
                }
                else
                {
                    if (fi <= c)
                        conv1 = vbar[dim] * kr;
                    if (fi > c)
                        conv1 = vbar[dim] * (c * kc + e * kr);
                }
            } /*
             vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
             if (vbar[dim] > 0.0) conv1 = vbar[dim]*kc;
             else                 conv1 = vbar[dim]*kr;
             vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
             if (vbar[dim] > 0.0) conv2 = vbar[dim]*kl;
             else                 conv2 = vbar[dim]*kc;
             return ((conv1 - conv2)/cdelta[dim]); */
        }
        else
        {
            vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            if (vbar[dim] > 0.0)
                conv1 = vbar[dim] * kc;
            else
                conv1 = vbar[dim] * kc;
            vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            if (vbar[dim] > 0.0)
                conv2 = vbar[dim] * kl;
            else
                conv2 = vbar[dim] * kc;
            return ((conv1 - conv2) / cdelta[dim]);
        }
    }
    // Get the velocity  v2bar(i-1/2,j) in the facet center
    vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
    if (vbar[dim] > 0.0)
    {
        if ((incell_l == 1) && (incell_ll == 1))
        {
            if (fabs(kc - kll) <= tol)
            {
                conv2 = vbar[dim] * kl;
            }
            else
            {
                fi = (kl - kll) / (kc - kll);
                if ((fi <= 0.0) || (fi >= 1.0))
                {
                    conv2 = vbar[dim] * kl;
                }
                else
                {
                    if (fi < b)
                        conv2 = vbar[dim] * (a * kl - c * kll);
                    if ((fi >= b) && (fi <= c))
                        conv2 = vbar[dim] * (b * kc + c * kl - d * kll);
                    if (fi > c)
                        conv2 = vbar[dim] * (c * kc + e * kl);
                }
            }
        }
        else if ((incell_l == 1) && (incell_ll == 0))
        {
            if (fabs(kc - kll) <= tol)
            {
                conv2 = vbar[dim] * kl;
            }
            else
            {
                fi = (kl - kll) / (kc - kll);
                if ((fi <= 0.0) || (fi >= 1.0))
                {
                    conv2 = vbar[dim] * kl;
                }
                else
                {
                    if (fi <= c)
                        conv2 = vbar[dim] * kl;
                    if (fi > c)
                        conv2 = vbar[dim] * (c * kc + e * kl);
                }
            } /*
                 vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
                 if (vbar[dim] > 0.0) conv1 = vbar[dim]*kc;
                 else                 conv1 = vbar[dim]*kr;
                 vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
                 if (vbar[dim] > 0.0) conv2 = vbar[dim]*kl;
                 else                 conv2 = vbar[dim]*kc;
                 return ((conv1 - conv2)/cdelta[dim]); */
        }
        else
        {
            vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            if (vbar[dim] > 0.0)
                conv1 = vbar[dim] * kc;
            else
                conv1 = vbar[dim] * kr;
            vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            if (vbar[dim] > 0.0)
                conv2 = vbar[dim] * kc;
            else
                conv2 = vbar[dim] * kc;
            return ((conv1 - conv2) / cdelta[dim]);
        }
    }
    else
    {
        // v2bar < 0.0
        if (fabs(kl - kr) <= tol)
        {
            conv2 = vbar[dim] * kc;
        }
        else
        {
            fi = (kc - kr) / (kl - kr);
            if ((fi <= 0.0) || (fi >= 1.0))
            {
                conv2 = vbar[dim] * kc;
            }
            else
            {
                if (fi < b)
                {
                    if (incell_r == 1)
                        conv2 = vbar[dim] * (a * kc - c * kr);
                    else
                        conv2 = vbar[dim] * kc;
                }
                if ((fi >= b) && (fi <= c))
                {
                    if ((incell_l == 1) && (incell_r == 1))
                        conv2 = vbar[dim] * (c * kc + b * kl - d * kr);
                    else
                        conv2 = vbar[dim] * kc;
                }
                if (fi > c)
                {
                    if (incell_l == 1)
                        conv2 = vbar[dim] * (e * kc + c * kl);
                    else
                        conv2 = vbar[dim] * kc;
                }
            }
        }
    }
    return ((conv1 - conv2) / cdelta[dim]);
}

// *******************************************************************
// Constitutive Equation Step for the Implicit Euler Method
// *******************************************************************
// Solve the evolution equation of the microstructure tensor <pp>
void higflow_implicit_euler_evolution_equation_microstructure_tensor(higflow_solver *ns)
{
    if ((ns->contr.flowtype == 9))
    {
        // Get the constants
        real beta = ns->ed.stsp.par.beta;
        real small = 1.0e-14;
        real dt = ns->par.dt;
        real phi;
        if (ns->ed.stsp.contr.model == GW)
        {
            // Value of phi for the GW model
            phi = 1.0;
        }
        else
        {
            // Value of phi for the GW-WC model
            phi = ns->ed.stsp.par.phi;
        }

        real Amax[DIM][DIM], Amin[DIM][DIM];
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                Amax[i][j] = -1.0e16;
                Amin[i][j] = 1.0e16;
            }
        }
        real RHSmax[DIM][DIM], RHSmin[DIM][DIM];
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                RHSmax[i][j] = -1.0e16;
                RHSmin[i][j] = 1.0e16;
            }
        }
        // tol        = 1.0e-3;
        //  Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for (int i = 0; i < DIM; i++)
        {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdp);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid = mp_lookup(mp, hig_get_cid(c));
            // Get the center of the cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Get the velocity derivative tensor Du and the microstructure tensor A
            real Du[DIM][DIM], A[DIM][DIM];
            real AM = 0.0;
            // Get the tensors
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    // Get Du
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
                    // Get microstructure tensor A =  <pp>
                    A[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpA[i][j], ns->ed.stn);
                }
                AM += A[i][i];
            }

            // Calculate the rate of strain tensor E
            real E[DIM][DIM], L[DIM][DIM];
            real Omega[DIM][DIM];
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    L[i][j] = Du[i][j];
                    // L[i][j] = Du[j][i];
                    E[i][j] = 0.5 * (Du[i][j] + Du[j][i]);
                    Omega[i][j] = Du[i][j];
                }
            }

            // Eigen-values and eigen-vectors of the deformation tensor E
            real R[DIM][DIM], lambda[DIM];
            hig_flow_jacobi(E, lambda, R);
            // Calculation of the compressive rate of strain Ec and the extensional rate of strain Ee
            real E_C[DIM][DIM], E_E[DIM][DIM];
            hig_flow_compressive_extensional_tensors(lambda, R, E_C, E_E);
            // Kronecker's delta function
            real KD[DIM][DIM];
            // real KM;
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    if (i == j)
                    {
                        KD[i][i] = 1.0;
                    }
                    else
                    {
                        KD[i][j] = 0.0;
                    }
                    // printf("%lf     ", KD[i][j]);
                }
                // KM += KD[i][i];
                // printf("\n")
            }
            // Calculation of the fourth order orientation moment <pppp>
            real FOM[DIM][DIM][DIM][DIM];
            hig_flow_calculate_4th_order_orientation_moment(AM, KD, A, FOM);
            // Calculate the trace of tensor E_C
            real LP[DIM][DIM], EP[DIM][DIM];
            // Calculation of the term (Du)^T : <pppp>
            hig_flow_calculate_double_dot_product_tensors(L, FOM, LP);
            // Calculation of the term E_E : <pppp>
            hig_flow_calculate_double_dot_product_tensors(E_E, FOM, EP);
            // Calculate the trace of tensor E_C
            real TrEC = 0.0;
            for (int i = 0; i < DIM; i++)
            {
                TrEC += E_C[i][i];
            }
            // Calculate RHS = - 2*(Du)^T : <pppp> - beta*( E_E:<pppp> + (phi/15)*(2*E_C+Tr(E_C)*KD) )
            real RHS[DIM][DIM];
            hig_flow_implicit_evolution_equation_microstructure_tensor_rhs(beta, phi, A, E_C, TrEC, KD, FOM, LP, EP, RHS);
            // Get the velocity at cell center
            real u[DIM], dAdx[DIM], d2Adx2[DIM];
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            real bS[DIM * DIM];
            real w[DIM * DIM][DIM * DIM + 1];
            // Solving the Constitutive Equation using the Euler Method
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    // Right hand side equation
                    real rhs = 0.0;
                    switch (ns->ed.stsp.contr.convecdiscrtype)
                    {
                    case CELL_CENTRAL:
                        // Tensor derivative at cell center
                        hig_flow_derivative_tensor_A_at_center_cell(ns, ccenter, cdelta, i, j, A[i][j], dAdx);
                        for (int dim = 0; dim < DIM; dim++)
                        {
                            // Compute convective tensor term in rhs
                            rhs -= u[dim] * dAdx[dim];
                        }
                        break;
                    case CELL_CUBISTA:
                        // Compute convective tensor term CUBISTA in rhs
                        for (int dim = 0; dim < DIM; dim++)
                        {
                            // rhs -= hig_flow_convective_tensor_A_term_cubista(ns, ns->dpu[dim], ns->ed.sdED, ns->ed.stn, A, ccenter, cdelta, dim, i, j);
                            rhs -= 0.0;
                        }
                        break;
                    }
                    // rhs    += higflow_diffusive_shear_banding_conformation_tensor_A_term(ns);
                    // if (rhs > RHSmax[i][j]) RHSmax[i][j] = rhs;
                    // if (rhs < RHSmin[i][j]) RHSmin[i][j] = rhs;
                    //  Compute the final rhs
                    rhs += RHS[i][j];
                    rhs *= dt;
                    rhs += A[i][j];
                    bS[i * DIM + j] = rhs;
                    // Calculate RHS max and min
                    //  Compute A at next time
                    // real Anew  = A[i][j] + ns->par.dt * rhs;
                    //  Store A in S
                    // dp_set_value(ns->ed.vesb.dpA[i][j], clid, Anew);
                    // if (i != j) {
                    //      dp_set_value(ns->ed.vesb.dpA[j][i], clid, Anew);
                    // }
                    // Calculate max and min A
                    // if (Anew > Amax[i][j]) Amax[i][j] = Anew;
                    // if (Anew < Amin[i][j]) Amin[i][j] = Anew;
                }
            }
            // Calculate de kronecker product Omega*I - I*Omega
            // Omega = Du
            hig_flow_kernel_system_matrix_shear_thickening_suspensions(w, Omega, dt);
            // Solve the linear system
            hig_flow_solve_system_constitutive_equation(DIM * DIM, w, bS);
            // Get the solution of linear system
            for (int i = 0; i < DIM; i++)
            {
                for (int j = i; j < DIM; j++)
                {
                    // Get the value of kernel
                    real Anew = bS[i * DIM + j];
                    // Set the value of kernel
                    dp_set_value(ns->ed.stsp.dpA[i][j], clid, Anew);
                    if (i != j)
                    {
                        dp_set_value(ns->ed.stsp.dpA[j][i], clid, Anew);
                    }
                    // Calculate max and min A
                    if (Anew > Amax[i][j])
                        Amax[i][j] = Anew;
                    if (Anew < Amin[i][j])
                        Amin[i][j] = Anew;
                }
            }
        }

        // Destroy the iterator
        // higcit_destroy(it);
        // Sync the ditributed pressure property
        // for (int i = 0; i < DIM; i++) {
        //    for (int j = 0; j < DIM; j++) {
        //        dp_sync(ns->ed.vesb.dpS[i][j]);
        //    }
        //}

        // Store the Kernel Tensor
        // for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
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

        // for (int i = 0; i < DIM; i++) {
        //    for (int j = 0; j < DIM; j++) {
        // Printing the min and max tensor
        //       printf("===> %d %d: Amin = %lf <===> Amax = %lf <===\n",i,j, Amin[i][j], Amax[i][j]);
        //   }
        //}

        // Printing the min and max deformation tensor values
        // for (int i = 0; i < DIM; i++) {
        // for (int j = 0; j < DIM; j++) {
        //  Printing the min and max tensor
        // printf("===> %d %d: Dmin = %lf <===> Dmax = %lf <===\n",i,j,Dmin[i][j],Dmax[i][j]);
        //    }
        //}

        // Printing the min and max deformation tensor values
        // for (int i = 0; i < DIM; i++) {
        // for (int j = 0; j < DIM; j++) {
        //  Printing the min and max tensor
        // printf("===> %d %d: ANAmin = %lf <===> ANAmax = %lf <===\n",i,j,ANAmin[i][j],ANAmax[i][j]);
        //}
        //}

        // for (int i = 0; i < DIM; i++) {
        // for (int j = 0; j < DIM; j++) {
        //  Printing the min and max tensor
        // printf("===> %d %d: RHSAmin = %lf <===> RHSAmax = %lf <===\n",i,j, RHSmin[i][j], RHSmax[i][j]);
        //}
        //}

        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed pressure property
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                dp_sync(ns->ed.stsp.dpA[i][j]);
            }
        }
    }
}

// *******************************************************************
// Navier-Stokes Step for the Explicit Euler Method
// *******************************************************************
void higflow_explicit_euler_intermediate_velocity_shear_thickening_suspensions(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM])
{
    // Get the facet iterator
    higfit_facetiterator *fit;
    // Get the local domain for cell
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++)
    {
        // Get the local domain for facet cell
        for (int dim2 = 0; dim2 < DIM; dim2++)
        {
            sfdu[dim2] = psfd_get_local_domain(ns->psfdu[dim2]);
        }
        // Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit))
        {
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
            higflow_computational_cell_shear_thickening_suspensions(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Convective term contribution: we neglect the convective term of the Navier-Stokes equation in the simulation of shear-thickening suspensions
            // rhs -= higflow_convective_term(ns, fdelta, dim);
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
void higflow_explicit_runge_kutta_2_intermediate_velocity_shear_thickening_suspensions(higflow_solver *ns)
{
    // Calculate the auxiliar velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_shear_thickening_suspensions(ns, ns->dpu, ns->dpuaux);
    // Calculate the star velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_shear_thickening_suspensions(ns, ns->dpuaux, ns->dpustar);
    // Calculate the order 2 Runge-Kutta method using the euler method
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Loop for each dimension
    higfit_facetiterator *fit;
    for (int dim = 0; dim < DIM; dim++)
    {
        // Get the local partitioned domain for facets
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit))
        {
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
            real u = dp_get_value(ns->dpu[dim], flid);
            real ustar = dp_get_value(ns->dpustar[dim], flid);
            // Compute the intermediate velocity
            real urk2 = 0.5 * (u + ustar);
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
void higflow_explicit_runge_kutta_3_intermediate_velocity_shear_thickening_suspensions(higflow_solver *ns)
{
    // Calculate the auxiliar velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_shear_thickening_suspensions(ns, ns->dpu, ns->dpuaux);
    // Calculate the second stage velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_shear_thickening_suspensions(ns, ns->dpuaux, ns->dpustar);
    // Calculate the order 2 Runge-Kutta method using the euler method
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Loop for each dimension
    higfit_facetiterator *fit;
    for (int dim = 0; dim < DIM; dim++)
    {
        // Get the local partitioned domain for facets
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit))
        {
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
            real u = dp_get_value(ns->dpu[dim], flid);
            real ustar = dp_get_value(ns->dpustar[dim], flid);
            // Compute the intermediate velocity
            real urk3 = 0.75 * u + 0.25 * ustar;
            // Set the final velocity in the distributed velocity property
            dp_set_value(ns->dpuaux[dim], flid, urk3);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the ditributed velocity property
        dp_sync(ns->dpuaux[dim]);
    }
    // Calculate the order 2 Runge-Kutta method using the euler method
    higflow_explicit_euler_intermediate_velocity_shear_thickening_suspensions(ns, ns->dpuaux, ns->dpustar);
    // Loop for each dimension
    for (int dim = 0; dim < DIM; dim++)
    {
        // Get the local partitioned domain for facets
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        // Calculate the third stage velocity by the explicit euler method
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit))
        {
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
            real u = dp_get_value(ns->dpu[dim], flid);
            real ustar = dp_get_value(ns->dpustar[dim], flid);
            // Compute the intermediate velocity
            real urk3 = u / 3.0 + 2.0 * ustar / 3.0;
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
void higflow_semi_implicit_euler_intermediate_velocity_shear_thickening_suspensions(higflow_solver *ns)
{
    // Get the facet iterator
    higfit_facetiterator *fit;
    // Get the local domain for cell
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Get the local domain for facet cell
    for (int dim2 = 0; dim2 < DIM; dim2++)
    {
        sfdu[dim2] = psfd_get_local_domain(ns->psfdu[dim2]);
    }
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++)
    {
        // Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit))
        {
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
            higflow_computational_cell_shear_thickening_suspensions(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // We divide the past terms over the density, which has the value of Re
            rhs /= ns->par.Re;
            // Convective term contribution: we neglect the convective term of the Navier-Stokes in the simulation of shear-thickening suspensions
            // This term is only considered when solving the lid-driven cavity problem and is no divided over the density
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Total contribuition terms by delta t
            rhs *= ns->par.dt;
            // Velocity term contribution
            rhs += ns->cc.ucell;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn, rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for (int dim2 = 0; dim2 < DIM; dim2++)
            {
                // Stencil weight update
                // real w = - ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
                // We set Re=1 for the simulation of shear-thickening suspensions, and we always assume rho = 1.0
                // value of eta_s used in the simulation of suspensions: eta_S = 0.00001
                // real w = - 1.0*(ns->par.dt)*(0.01 + ns->ed.stsp.par.eta0)/(fdelta[dim2]*fdelta[dim2]);
                // real w = - 1.0*(ns->par.dt)*(0.1 + ns->ed.stsp.par.eta0)/(fdelta[dim2]*fdelta[dim2]);
                // Here, Re represents the density value, which is a constant
                //real w = -1.0 * (ns->par.dt) * (0.01 + ns->ed.stsp.par.eta0) / (ns->par.Re * fdelta[dim2] * fdelta[dim2]);
                //real w = -1.0 * (ns->par.dt) * (ns->ed.stsp.par.eta0) / (ns->par.Re * fdelta[dim2] * fdelta[dim2]);
                real w = -1.0 * (ns->par.dt) * (ns->ed.stsp.par.eta0 + ns->ed.stsp.par.eta0) / (ns->par.Re * fdelta[dim2] * fdelta[dim2]);
                // real w = - 2.0*(ns->par.dt)*(ns->ed.stsp.par.eta0)/(fdelta[dim2]*fdelta[dim2]);
                alpha -= 2.0 * w;
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
            int *ids = psfd_stn_get_gids(ns->psfdu[dim], ns->stn);
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
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit))
        {
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
void higflow_semi_implicit_crank_nicolson_intermediate_velocity_shear_thickening_suspensions(higflow_solver *ns)
{
    // Get the facet iterator
    higfit_facetiterator *fit;
    // Get the local domain for cell
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Get the local domain for facet cell
    for (int dim2 = 0; dim2 < DIM; dim2++)
    {
        sfdu[dim2] = psfd_get_local_domain(ns->psfdu[dim2]);
    }
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++)
    {
        // Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit))
        {
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
            higflow_computational_cell_shear_thickening_suspensions(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
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
            // Convective term contribution: we neglect the convective term of the N-S equations in the simulation of shear-thickening suspensions
            // rhs -= higflow_convective_term(ns, fdelta, dim);
            // Total contribuition terms times delta t
            rhs *= ns->par.dt;
            // Velocity term contribution
            rhs += ns->cc.ucell;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn, rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for (int dim2 = 0; dim2 < DIM; dim2++)
            {
                // Stencil weight update
                // real w = - 0.5 * ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
                // We set Re = 1.0 and we always consider rho = 1.0 in the simulation of shear-thickening suspensions
                real w = -0.5 * 2.0 * (ns->par.dt) * (ns->ed.stsp.par.eta0) / (fdelta[dim2] * fdelta[dim2]);
                alpha -= 2.0 * w;
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
            int *ids = psfd_stn_get_gids(ns->psfdu[dim], ns->stn);
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
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit))
        {
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
void higflow_semi_implicit_bdf2_intermediate_velocity_shear_thickening_suspensions(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM])
{
    // Firt stage of Tr-BDF2 method
    // Get the facet iterator
    higfit_facetiterator *fit;
    // Get the local domain for cell
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Get the local domain for facet cell
    for (int dim2 = 0; dim2 < DIM; dim2++)
    {
        sfdu[dim2] = psfd_get_local_domain(ns->psfdu[dim2]);
    }
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++)
    {
        // Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit))
        {
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
            higflow_computational_cell_shear_thickening_suspensions(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Convective term contribution: we neglect the convective term of the N-S equations in the simulation of shear-thickening suspensions
            // rhs -= higflow_convective_term(ns, fdelta, dim);
            // Total contribuition terms times delta t
            rhs *= 0.25 * ns->par.dt;
            // Velocity term contribution
            rhs += ns->cc.ucell;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn, rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for (int dim2 = 0; dim2 < DIM; dim2++)
            {
                // Stencil weight update
                // real w = - 0.25*ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
                real w = -0.25 * 2.0 * (ns->par.dt) * (ns->ed.stsp.par.eta0) / (fdelta[dim2] * fdelta[dim2]);
                alpha -= 2.0 * w; // divide po 4 para usar regra trapezio em t(n+1/2)
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
            int *ids = psfd_stn_get_gids(ns->psfdu[dim], ns->stn);
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
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit))
        {
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
    //  Looping for the velocity
    for (int dim = 0; dim < DIM; dim++)
    {
        // Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit))
        {
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
            higflow_computational_cell_shear_thickening_suspensions(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Get the uaux
            real uaux = dp_get_value(ns->dpuaux[dim], flid);
            // Right hand side equation
            real rhs = 0.0;
            rhs = (4.0 * uaux - ns->cc.ucell) / 3.0;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn, rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for (int dim2 = 0; dim2 < DIM; dim2++)
            {
                // Stencil weight update
                // real w = - 1.0/3.0*ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
                real w = -1.0 / 3.0 * 2.0 * (ns->par.dt) * (ns->ed.stsp.par.eta0) / (fdelta[dim2] * fdelta[dim2]);
                alpha -= 2.0 * w;
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
            int *ids = psfd_stn_get_gids(ns->psfdu[dim], ns->stn);
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
        // Vec *vecu = slv_get_solution_vec(ns->slvu[dim]);
        // Gets the values of the solution
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit))
        {
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
void higflow_solver_step_shear_thickening_suspensions(higflow_solver *ns)
{
    // Boundary condition for velocity
    higflow_boundary_condition_for_velocity(ns);
    // Calculate the source term
    higflow_calculate_source_term(ns);
    // Calculate the facet source term
    higflow_calculate_facet_source_term(ns);
    // Calculate the intermediated velocity
    switch (ns->contr.tempdiscrtype)
    {
    case 0:
        // Explicit Euler method
        higflow_explicit_euler_intermediate_velocity_shear_thickening_suspensions(ns, ns->dpu, ns->dpustar);
        break;
    case 1:
        // Explicit RK2 method
        higflow_explicit_runge_kutta_2_intermediate_velocity_shear_thickening_suspensions(ns);
        break;
    case 2:
        // Explicit RK3 method
        higflow_explicit_runge_kutta_3_intermediate_velocity_shear_thickening_suspensions(ns);
        break;
    case 3:
        // Semi-Implicit Euler Method
        higflow_semi_implicit_euler_intermediate_velocity_shear_thickening_suspensions(ns);
        break;
    case 4:
        // Semi-Implicit Crank-Nicolson Method
        higflow_semi_implicit_crank_nicolson_intermediate_velocity_shear_thickening_suspensions(ns);
        break;
    case 5:
        // Semi-Implicit Crank-Nicolson Method
        higflow_semi_implicit_bdf2_intermediate_velocity_shear_thickening_suspensions(ns, ns->dpu, ns->dpustar);
        break;
    }
    // Set outflow for ustar velocity
    // higflow_outflow_ustar_step(ns);
    // Boundary condition for pressure
    higflow_boundary_condition_for_pressure(ns);
    // Calculate the pressure
    higflow_pressure(ns);
    // Calculate the final velocity
    higflow_final_velocity(ns);
    // Boundary condition for velocity
    higflow_boundary_condition_for_velocity(ns);
    // Set outflow for velocity
    // higflow_outflow_u_step(ns);
    // Calculate the final pressure
    higflow_final_pressure(ns);
    // Calculate the velocity derivative tensor
    higflow_compute_velocity_derivative_tensor(ns);
    //solving the particle migration equation
    switch (ns->ed.stsp.contr.discrtype)
    {
    case EXPLICIT_EULER:
        // Explicit method
        higflow_explicit_euler_evolution_equation_microstructure_tensor(ns);
        break;
    case SEMI_IMPLICIT_EULER:
        // Implicit method
        higflow_implicit_euler_evolution_equation_microstructure_tensor(ns);
        break;
    }
    // Constitutive Equation Step for the Explicit Euler Method
    if (ns->ed.stsp.contr.model == GW_WC_IF){
        higflow_compute_particle_stress_tensor_shear_thickening_suspension(ns);
        higflow_explicit_euler_volume_fraction_equation(ns);
    }
    // Computing the Polymeric Tensor
    higflow_compute_polymeric_tensor_shear_thickening_suspension(ns);
}

// Calculate the eige-value and eige-vectors using the Jacobi method
void hig_flow_jacobi(real A[DIM][DIM], real d[DIM], real V[DIM][DIM])
{
    real b[DIM];
    real z[DIM];
    int NMAX = 500;
    int nrot = 0;
    for (int ip = 0; ip < DIM; ip++)
    {
        for (int iq = 0; iq < DIM; iq++)
        {
            V[ip][iq] = 0.0;
        }
        V[ip][ip] = 1.0;
        b[ip] = A[ip][ip];
        d[ip] = b[ip];
        z[ip] = 0.0;
    }
    for (int i = 1; i <= NMAX; i++)
    {
        real tresh;
        real sm = 0.0;
        for (int ip = 0; ip < DIM - 1; ip++)
            for (int iq = ip + 1; iq < DIM; iq++)
                sm += fabs(A[ip][iq]);
        if (sm == 0.0)
            return;
        if (i < 4)
            tresh = 0.2 * sm / (DIM * DIM);
        else
            tresh = 0.0;
        for (int ip = 0; ip < DIM - 1; ip++)
            for (int iq = ip + 1; iq < DIM; iq++)
            {
                real g = 100.0 * fabs(A[ip][iq]);
                if ((i > 4) && (fabs(d[ip]) + g == fabs(d[ip])) &&
                    (fabs(d[iq]) + g == fabs(d[iq])))
                {
                    A[ip][iq] = 0.0;
                }
                else if (fabs(A[ip][iq]) > tresh)
                {
                    real t, theta;
                    real h = d[iq] - d[ip];
                    if (fabs(h) + g == fabs(h))
                    {
                        t = A[ip][iq] / h;
                    }
                    else
                    {
                        theta = 0.5 * h / A[ip][iq];
                        t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                        if (theta < 0.0)
                            t = -t;
                    }
                    real c = 1.0 / sqrt(1.0 + t * t);
                    real s = t * c;
                    real tau = s / (1.0 + c);
                    h = t * A[ip][iq];
                    z[ip] = z[ip] - h;
                    z[iq] = z[iq] + h;
                    d[ip] = d[ip] - h;
                    d[iq] = d[iq] + h;
                    A[ip][iq] = 0.0;
                    for (int j = 0; j <= ip - 1; j++)
                    {
                        g = A[j][ip];
                        h = A[j][iq];
                        A[j][ip] = g - s * (h + g * tau);
                        A[j][iq] = h + s * (g - h * tau);
                    }
                    for (int j = ip + 1; j <= iq - 1; j++)
                    {
                        g = A[ip][j];
                        h = A[j][iq];
                        A[ip][j] = g - s * (h + g * tau);
                        A[j][iq] = h + s * (g - h * tau);
                    }
                    for (int j = iq + 1; j < DIM; j++)
                    {
                        g = A[ip][j];
                        h = A[iq][j];
                        A[ip][j] = g - s * (h + g * tau);
                        A[iq][j] = h + s * (g - h * tau);
                    }
                    for (int j = 0; j < DIM; j++)
                    {
                        g = V[j][ip];
                        h = V[j][iq];
                        V[j][ip] = g - s * (h + g * tau);
                        V[j][iq] = h + s * (g - h * tau);
                    }
                    nrot++;
                }
            }
        for (int ip = 0; ip < DIM; ip++)
        {
            b[ip] = b[ip] + z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }
}

// Calculate the matrix product
void hig_flow_matrix_product(real A[DIM][DIM], real R[DIM][DIM], real B[DIM][DIM])
{
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            B[i][j] = 0.0;
            for (int k = 0; k < DIM; k++)
            {
                for (int l = 0; l < DIM; l++)
                {
                    B[i][j] += R[k][i] * A[k][l] * R[l][j];
                }
            }
        }
    }
}

// Calculate the matrix product transpose
void hig_flow_matrix_transpose_product(real A[DIM][DIM], real R[DIM][DIM], real B[DIM][DIM])
{
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            B[i][j] = 0.0;
            for (int k = 0; k < DIM; k++)
            {
                for (int l = 0; l < DIM; l++)
                {
                    B[i][j] += R[i][k] * A[k][l] * R[j][l];
                }
            }
        }
    }
}

// Step function 1
real step_function_1(real l)
{
    real value;
    if (l > 0.0)
    {
        value = 1.0;
    }
    else
    {
        value = 0.0;
    }
    return value;
}

// Step function 2
real step_function_2(real l)
{
    real value;
    if (l > 0.0)
    {
        value = 0.0;
    }
    else
    {
        value = 1.0;
    }
    return value;
}

// Calculates the tensors E_E and E_C
void hig_flow_compressive_extensional_tensors(real lambda[DIM], real R[DIM][DIM], real E_C[DIM][DIM], real E_E[DIM][DIM])
{
    // Calculation of E_E: the extensional (positive eigenvalues) part of the deformation tensor D
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            E_E[i][j] = 0.0;
            for (int k = 0; k < DIM; k++)
            {
                E_E[i][j] += lambda[k] * step_function_1(lambda[k]) * R[i][k] * R[j][k];
            }
            // printf("%lf     ", E_E[i][j]);
        }
        // printf("\n");
    }

    // Calculation of E_C: the compressive (positive eigenvalues) part of the deformation tensor D
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            E_C[i][j] = 0.0;
            for (int k = 0; k < DIM; k++)
            {
                E_C[i][j] += lambda[k] * step_function_2(lambda[k]) * R[i][k] * R[j][k];
            }
            // printf("%lf     ", E_C[i][j]);
        }
        // printf("\n");
    }
}

// Calculate the fourth-order orientation moment <pppp>
void hig_flow_calculate_4th_order_orientation_moment(real AM, real KD[DIM][DIM], real A[DIM][DIM], real FOM[DIM][DIM][DIM][DIM])
{
    // Calculation of GW using the linear closure model
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            for (int k = 0; k < DIM; k++)
            {
                for (int l = 0; l < DIM; l++)
                {
                    // FOM[i][j][k][l] = a[i][j]*a[k][l];
                    FOM[i][j][k][l] = 0.0;
                    // FOM[i][j][k][l] = -(1.0/35.0)*AM*(KD[i][j]*KD[k][l] + KD[i][k]*KD[j][l] + KD[i][l]*KD[j][k]) + (1.0/7.0)*(KD[i][j]*a[k][l] + KD[i][k]*a[j][l] + KD[i][l]*a[j][k] + a[i][j]*KD[k][l] + a[i][k]*KD[j][l] + a[i][l]*KD[j][k]);
                    FOM[i][j][k][l] += (-1.0) * (1.0 / 35.0) * AM * (KD[i][j] * KD[k][l] + KD[i][k] * KD[j][l] + KD[i][l] * KD[j][k]);
                    // FOM[i][j][k][l] = (KD[i][j]*KD[k][l] + KD[i][k]*KD[j][l] + KD[i][l]*KD[j][k]);
                    FOM[i][j][k][l] += (1.0 / 7.0) * (KD[i][j] * A[k][l] + KD[i][k] * A[j][l] + KD[i][l] * A[j][k]);
                    FOM[i][j][k][l] += (1.0 / 7.0) * (A[i][j] * KD[k][l] + A[i][k] * KD[j][l] + A[i][l] * KD[j][k]);
                    // printf("%lf     ", FOM[i][j][k][l]);
                }
                // printf("\n");
            }
        }
    }
}

// Calculate the suspension stress tensor
void hig_flow_calculate_suspension_stress_tensor(real tau_M[DIM][DIM], real FOM[DIM][DIM][DIM][DIM], real T[DIM][DIM])
{
    // Calculation of T (suspension stress tensor)
    for (int m = 0; m < DIM; m++)
    {
        for (int n = 0; n < DIM; n++)
        {
            T[m][n] = 0.0;
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    // FOM2[m][n] += cons*E[i][j]*FOM[i][j][m][n];
                    T[m][n] += tau_M[i][j] * FOM[m][n][i][j];
                }
            }
        }
    }
}

// Calculate the double dot product between two tensors (one is a second order tensor and the other one is a 4th order tensor)
void hig_flow_calculate_double_dot_product_tensors(real M[DIM][DIM], real FOM[DIM][DIM][DIM][DIM], real P[DIM][DIM])
{
    // Calculation of T (suspension stress tensor)
    for (int m = 0; m < DIM; m++)
    {
        for (int n = 0; n < DIM; n++)
        {
            P[m][n] = 0.0;
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    // FOM2[m][n] += cons*E[i][j]*FOM[i][j][m][n];
                    P[m][n] += M[i][j] * FOM[m][n][i][j];
                }
            }
        }
    }
}

// Calculate the jamming coordinate value chi
real hig_flow_calculate_jamming_coordinate_chi(real A[DIM][DIM], real E_C[DIM][DIM])
{
    real chi_v;
    // Calculate the first inner product of the tensors A and E_C
    real AE1[DIM][DIM];
    // Calculate the first inner product of the tensors E_C
    real EC2X1[DIM][DIM];
    for (int dim = 0; dim < DIM; dim++)
    {
        for (int dim2 = 0; dim2 < DIM; dim2++)
        {
            AE1[dim][dim2] = 0.0;
            EC2X1[dim][dim2] = 0.0;
            for (int dim3 = 0; dim3 < DIM; dim3++)
            {
                AE1[dim][dim2] += A[dim][dim3] * E_C[dim3][dim2];
                EC2X1[dim][dim2] += E_C[dim][dim3] * E_C[dim3][dim2];
            }
        }
    }

    // Double inner product of the tensors A and E_C
    real AE = 0.0;
    // Double inner product of the tensors E_C
    real EC2X = 0.0;
    for (int dim = 0; dim < DIM; dim++)
    {
        AE += AE1[dim][dim];
        EC2X += EC2X1[dim][dim];
    }
    EC2X = sqrt(EC2X);
    chi_v = -1.0 * AE / EC2X;
    return chi_v;
}

// Calculate RHS = A*Du + (Du)^T * A- 2*(Du)^T : <pppp> - beta*( E_E:<pppp> + (phi/15)*(2*E_C+Tr(E_C)*KD) )
void hig_flow_evolution_equation_microstructure_tensor_rhs(real beta, real phi, real A[DIM][DIM], real DU[DIM][DIM], real E_C[DIM][DIM], real TrEC, real KD[DIM][DIM], real FOM[DIM][DIM][DIM][DIM], real LP[DIM][DIM], real EP[DIM][DIM], real RHS[DIM][DIM])
{
    real DUA[DIM][DIM], ADU[DIM][DIM];
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            DUA[i][j] = 0.0;
            ADU[i][j] = 0.0;
            for (int k = 0; k < DIM; k++)
            {
                ADU[i][j] += DU[i][k] * A[k][j];
                DUA[i][j] += DU[j][k] * A[i][k];
            }
        }
    }

    // Calculate RHS = A*Du + (Du)^T * A- 2*(Du)^T : <pppp> - beta*( E_E:<pppp> + (phi/15)*(2*E_C+Tr(E_C)*KD) )
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            RHS[i][j] = ADU[i][j] + DUA[i][j] - 2.0 * LP[i][j] - beta * (EP[i][j] + phi * (2.0 * E_C[i][j] + TrEC * KD[i][j]) / 15.0);
        }
    }
}

// Calculate RHS for the implicit method RHS = - 2*(Du)^T : <pppp> - beta*( E_E:<pppp> + (phi/15)*(2*E_C+Tr(E_C)*KD) )
void hig_flow_implicit_evolution_equation_microstructure_tensor_rhs(real beta, real phi, real A[DIM][DIM], real E_C[DIM][DIM], real TrEC, real KD[DIM][DIM], real FOM[DIM][DIM][DIM][DIM], real LP[DIM][DIM], real EP[DIM][DIM], real RHS[DIM][DIM])
{
    // Calculate RHS = A*Du + (Du)^T * A- 2*(Du)^T : <pppp> - beta*( E_E:<pppp> + (phi/15)*(2*E_C+Tr(E_C)*KD) )
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            RHS[i][j] = -2.0 * LP[i][j] - beta * (EP[i][j] + phi * (2.0 * E_C[i][j] + TrEC * KD[i][j]) / 15.0);
        }
    }
}

// Get the velocity at cell center
void hig_flow_velocity_at_center_cell(higflow_solver *ns, Point ccenter, Point cdelta, real u[DIM])
{
    for (int dim = 0; dim < DIM; dim++)
    {
        // Verity if is in facet
        int infacet;
        // Get the velocity in the left facet center
        real ul = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
        // Get the velocity in the right facet center
        real ur = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
        // Setting the velocity at cell center
        u[dim] = 0.5 * (ul + ur);
    }
}

// Get the derivative of microstructure tensor A
void hig_flow_derivative_tensor_A_at_center_cell(higflow_solver *ns, Point ccenter, Point cdelta, int i, int j, real Acenter, real dAdx[DIM])
{
    for (int dim = 0; dim < DIM; dim++)
    {
        int incell_left, incell_right;
        // Get the Kernel in the left cell
        real Aleft = compute_center_p_left_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.stsp.dpA[i][j], ns->ed.stn, &incell_left);
        // Get the Kernel in the right cell
        real Aright = compute_center_p_right_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.stsp.dpA[i][j], ns->ed.stn, &incell_right);
        // Compute the Kernel derivative
        if ((incell_left == 1) && (incell_right == 1))
        {
            dAdx[dim] = compute_dpdx_at_point(cdelta, dim, 1.0, Aleft, Aright);
        }
        else if (incell_right == 1)
        {
            dAdx[dim] = compute_dpdxr_at_point(cdelta, dim, 1.0, Acenter, Aright);
        }
        else
        {
            dAdx[dim] = compute_dpdxl_at_point(cdelta, dim, 1.0, Aleft, Acenter);
        }
    }
}

// Gauss elimination to solve the constitutive equation
void hig_flow_solve_system_constitutive_equation(int n, real A[DIM * DIM][DIM * DIM + 1], real x[DIM * DIM])
{
    real s, max, aux;
    int l, k, i, j;
    // Assign the vector b to the expanded matrix of the system Ax = b
    for (j = 0; j < n; j++)
    {
        A[j][n] = x[j];
    }
    // Partial pivot
    for (k = 0; k < n - 1; k++)
    {
        max = fabs(A[k][k]);
        l = k;
        for (i = k + 1; i < n; i++)
        {
            if (fabs(A[i][k]) > max)
            {
                max = fabs(A[i][k]);
                l = i;
            }
        }
        if (l != k)
        {
            for (j = k; j <= n; j++)
            {
                aux = A[k][j];
                A[k][j] = A[l][j];
                A[l][j] = aux;
            }
        }
        // Gauss algorithm
        for (i = k + 1; i < n; i++)
        {
            s = A[i][k] / A[k][k];
            A[i][k] = 0;
            for (j = k + 1; j <= n; j++)
            {
                A[i][j] -= s * A[k][j];
            }
        }
    }
    x[n - 1] = A[n - 1][n] / A[n - 1][n - 1];
    for (i = n - 2; i >= 0; i--)
    {
        s = 0;
        for (j = i + 1; j < n; j++)
        {
            s += A[i][j] * x[j];
        }
        x[i] = (A[i][n] - s) / A[i][i];
    }
}

// Calculate the matrix product
void hig_flow_kernel_system_matrix_shear_thickening_suspensions(real w[DIM * DIM][DIM * DIM + 1], real Omega[DIM][DIM], real dt)
{
    real I[DIM][DIM];
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            I[i][j] = 0.0;
            if (i == j)
                I[i][j] = 1.0;
        }
    }
    for (int i = 0; i < DIM; i++)
    {
        for (int j = i; j < DIM; j++)
        {
            for (int k = 0; k < DIM; k++)
            {
                for (int l = 0; l < DIM; l++)
                {
                    if (i == j)
                    {
                        w[i * DIM + k][j * DIM + l] = I[k][l] - dt * Omega[k][l] - dt * Omega[i][i] * I[k][l];
                    }
                    else
                    {
                        w[i * DIM + k][j * DIM + l] = dt * Omega[j][i] * I[k][l];
                        w[j * DIM + l][i * DIM + k] = -dt * Omega[j][i] * I[k][l];
                    }
                }
            }
        }
    }
}

// Get the derivative of the volume fraction
void hig_flow_derivative_volfrac_at_center_cell(higflow_solver *ns, Point ccenter, Point cdelta, real ncenter, real dndx[DIM])
{
    for (int dim = 0; dim < DIM; dim++)
    {
        int incell_left, incell_right;
        // Get the fraction value in the left cell
        real nleft = compute_center_p_left_22(ns->ed.stsp.sdphi, ccenter, cdelta, dim, 1.0, ns->ed.stsp.dpphi, ns->ed.stn, &incell_left);
        // Get the fraction value in the right cell
        real nright = compute_center_p_right_22(ns->ed.stsp.sdphi, ccenter, cdelta, dim, 1.0, ns->ed.stsp.dpphi, ns->ed.stn, &incell_left);
        // Compute the volume fraction derivative
        dndx[dim] = compute_dpdx_at_point(cdelta, dim, 1.0, nleft, nright);
    }
}

// Calculate convective term CUBISTA for volume fraction
// *******************************************************************
real hig_flow_fraction_volume_suspensions_term_cubista(higflow_solver *ns, distributed_property *dpu, sim_domain *sdp, sim_stencil *stn, real varphi, Point ccenter, Point cdelta, int dim)
{
    real vbar[DIM], dKdx[dim], kr, krr, kl, kll, kc, a, b, c, d, e, tol, frac_tol, fi, conv1, conv2;
    a = 1.7500;
    b = 0.3750;
    c = 0.7500;
    d = 0.1250;
    e = 0.2500;
    tol = 1.0e-14;
    conv1 = 0.0;
    conv2 = 0.0;
    int incell_r, incell_l, incell_ll, incell_rr, infacet;
    // Get the kernel at center cell
    kc = varphi;
    // Get the low, high, lowlow, highhigh component kernel at center cell
    kl = compute_center_p_left_22(sdp, ccenter, cdelta, dim, 1.0, ns->ed.stsp.dpphi, ns->ed.stn, &incell_l);
    kr = compute_center_p_right_22(sdp, ccenter, cdelta, dim, 1.0, ns->ed.stsp.dpphi, ns->ed.stn, &incell_r);
    kll = compute_center_p_left_22(sdp, ccenter, cdelta, dim, 2.0, ns->ed.stsp.dpphi, ns->ed.stn, &incell_ll);
    krr = compute_center_p_right_22(sdp, ccenter, cdelta, dim, 2.0, ns->ed.stsp.dpphi, ns->ed.stn, &incell_rr);
    // Get the velocity  v1bar(i+1/2,j) in the facet center
    vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
    if (vbar[dim] > 0.0)
    {
        if (fabs(kr - kl) <= tol)
        {
            conv1 = vbar[dim] * kc;
        }
        else
        {
            fi = (kc - kl) / (kr - kl);
            if ((fi <= 0.0) || (fi >= 1.0))
            {
                conv1 = vbar[dim] * kc;
            }
            else
            {
                if (fi < b)
                {
                    if (incell_l == 1)
                        conv1 = vbar[dim] * (a * kc - c * kl);
                    else
                        conv1 = vbar[dim] * kc;
                }
                if ((fi >= b) && (fi <= c))
                {
                    if ((incell_l == 1) && (incell_r == 1))
                        conv1 = vbar[dim] * (c * kc + b * kr - d * kl);
                    else
                        conv1 = vbar[dim] * kc;
                }
                if (fi > c)
                {
                    if (incell_r == 1)
                        conv1 = vbar[dim] * (e * kc + c * kr);
                    else
                        conv1 = vbar[dim] * kc;
                }
            }
        }
        // v1bar < 0.0
    }
    else
    {
        if ((incell_r == 1) && (incell_rr == 1))
        {
            if (fabs(kc - krr) <= tol)
            {
                conv1 = vbar[dim] * kr;
            }
            else
            {
                fi = (kr - krr) / (kc - krr);
                if ((fi <= 0.0) || (fi >= 1.0))
                {
                    conv1 = vbar[dim] * kr;
                }
                else
                {
                    if (fi < b)
                        conv1 = vbar[dim] * (a * kr - c * krr);
                    if ((fi >= b) && (fi <= c))
                        conv1 = vbar[dim] * (c * kr + b * kc - d * krr);
                    if (fi > c)
                        conv1 = vbar[dim] * (c * kc + e * kr);
                }
            }
            // Return upwind value at boundary
        }
        else if ((incell_r == 1) && (incell_rr == 0))
        {
            if (fabs(kc - krr) <= tol)
            {
                conv1 = vbar[dim] * kr;
            }
            else
            {
                fi = (kr - krr) / (kc - krr);
                if ((fi <= 0.0) || (fi >= 1.0))
                {
                    conv1 = vbar[dim] * kr;
                }
                else
                {
                    if (fi <= c)
                        conv1 = vbar[dim] * kr;
                    if (fi > c)
                        conv1 = vbar[dim] * (c * kc + e * kr);
                }
            } /*
             vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
             if (vbar[dim] > 0.0) conv1 = vbar[dim]*kc;
             else                 conv1 = vbar[dim]*kr;
             vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
             if (vbar[dim] > 0.0) conv2 = vbar[dim]*kl;
             else                 conv2 = vbar[dim]*kc;
             return ((conv1 - conv2)/cdelta[dim]); */
        }
        else
        {
            vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            if (vbar[dim] > 0.0)
                conv1 = vbar[dim] * kc;
            else
                conv1 = vbar[dim] * kc;
            vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            if (vbar[dim] > 0.0)
                conv2 = vbar[dim] * kl;
            else
                conv2 = vbar[dim] * kc;

            real value_tol = ((conv1 - conv2) / cdelta[dim]);
            return value_tol;
        }
    }
    // Get the velocity  v2bar(i-1/2,j) in the facet center
    vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
    if (vbar[dim] > 0.0)
    {
        if ((incell_l == 1) && (incell_ll == 1))
        {
            if (fabs(kc - kll) <= tol)
            {
                conv2 = vbar[dim] * kl;
            }
            else
            {
                fi = (kl - kll) / (kc - kll);
                if ((fi <= 0.0) || (fi >= 1.0))
                {
                    conv2 = vbar[dim] * kl;
                }
                else
                {
                    if (fi < b)
                        conv2 = vbar[dim] * (a * kl - c * kll);
                    if ((fi >= b) && (fi <= c))
                        conv2 = vbar[dim] * (b * kc + c * kl - d * kll);
                    if (fi > c)
                        conv2 = vbar[dim] * (c * kc + e * kl);
                }
            }
        }
        else if ((incell_l == 1) && (incell_ll == 0))
        {
            if (fabs(kc - kll) <= tol)
            {
                conv2 = vbar[dim] * kl;
            }
            else
            {
                fi = (kl - kll) / (kc - kll);
                if ((fi <= 0.0) || (fi >= 1.0))
                {
                    conv2 = vbar[dim] * kl;
                }
                else
                {
                    if (fi <= c)
                        conv2 = vbar[dim] * kl;
                    if (fi > c)
                        conv2 = vbar[dim] * (c * kc + e * kl);
                }
            } /*
                  vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
                  if (vbar[dim] > 0.0) conv1 = vbar[dim]*kc;
                  else                 conv1 = vbar[dim]*kr;
                  vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
                  if (vbar[dim] > 0.0) conv2 = vbar[dim]*kl;
                  else                 conv2 = vbar[dim]*kc;
                  return ((conv1 - conv2)/cdelta[dim]); */
        }
        else
        {
            vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            if (vbar[dim] > 0.0)
                conv1 = vbar[dim] * kc;
            else
                conv1 = vbar[dim] * kr;
            vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            if (vbar[dim] > 0.0)
                conv2 = vbar[dim] * kc;
            else
                conv2 = vbar[dim] * kc;

            real value_tol = ((conv1 - conv2) / cdelta[dim]);
            return value_tol;
        }
    }
    else
    {
        // v2bar < 0.0
        if (fabs(kl - kr) <= tol)
        {
            conv2 = vbar[dim] * kc;
        }
        else
        {
            fi = (kc - kr) / (kl - kr);
            if ((fi <= 0.0) || (fi >= 1.0))
            {
                conv2 = vbar[dim] * kc;
            }
            else
            {
                if (fi < b)
                {
                    if (incell_r == 1)
                        conv2 = vbar[dim] * (a * kc - c * kr);
                    else
                        conv2 = vbar[dim] * kc;
                }
                if ((fi >= b) && (fi <= c))
                {
                    if ((incell_l == 1) && (incell_r == 1))
                        conv2 = vbar[dim] * (c * kc + b * kl - d * kr);
                    else
                        conv2 = vbar[dim] * kc;
                }
                if (fi > c)
                {
                    if (incell_l == 1)
                        conv2 = vbar[dim] * (e * kc + c * kl);
                    else
                        conv2 = vbar[dim] * kc;
                }
            }
        }
    }

    real value_tol = ((conv1 - conv2) / cdelta[dim]);
    return value_tol;
}

// *******************************************************************
// Volume fraction evolution equation
// *******************************************************************
void higflow_explicit_euler_volume_fraction_equation(higflow_solver *ns)
{
    if (ns->ed.stsp.contr.model == GW_WC_IF)
    {
        real apsize = ns->ed.stsp.par.apsize;
        real rho = ns->par.Re;
        real eta0 = ns->ed.stsp.par.eta0;
        real volfracmax = -1.0e16;
        real volfracmin = 1.0e16;
        real phij = ns->ed.stsp.par.phi;
        real frac_tol  = ns->ed.stsp.par.gdrms;
        real eps = 1.0e-4;
        //DEBUG_INSPECT(apsize, %lf);
        // Get the constants
        // Get the local sub-domain for the cells
        sim_domain *sdphi = psd_get_local_domain(ns->ed.stsp.psdphi);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for (int i = 0; i < DIM; i++)
        {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdphi);
        // Loop for each cell
        higcit_celliterator *it;
        //printf("=+=+=+= WE ARE HERE before loop =+=+=+=\n");
        for (it = sd_get_domain_celliterator(sdphi); !higcit_isfinished(it); higcit_nextcell(it)) {
            //printf("=+=+=+= WE ARE HERE after loop =+=+=+=\n");
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid = mp_lookup(mp, hig_get_cid(c));
            // Get the center of the cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Get the velocity at cell center
            real u[DIM], dphidx[DIM];
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            // compute nplus value at point
            real varphi = compute_value_at_point(ns->ed.stsp.sdphi, ccenter, ccenter, 1.0, ns->ed.stsp.dpphi, ns->ed.stn);
            if(varphi > phij- frac_tol){
                varphi = phij- frac_tol;
            }
            if(varphi < eps){
                varphi = 0.0;
            }
            //DEBUG_INSPECT(varphi, %lf);
            // Solving the Transport Equation using the Euler Method
            // Right hand side equation
            real rhs = 0.0;
            ns->cc.d2Tdx2 = 0.0;
            switch (ns->ed.stsp.contr.volfracconvecdiscrtype)
            {
            // Central scheme
            case CELL_CENTRAL:
                // volume fraction derivative at cell center
                hig_flow_derivative_volfrac_at_center_cell(ns, ccenter, cdelta, varphi, dphidx);
                for (int dim = 0; dim < DIM; dim++)
                {
                    // Set the computational cell
                    //higflow_computational_cell_electroosmotic_ionic(ns, sdnplus, sdpsi, sdphi, clid, ccenter, cdelta, dim, ns->ed.eo.dpnplus, ns->ed.eo.dppsi, ns->ed.eo.dpphi);
                    higflow_computational_cell_volume_fraction_suspensions(ns, sdphi, clid, ccenter, cdelta, dim, ns->ed.stsp.dpphi);
                    // Compute convective inonic term in rhs
                    //rhs = 5.0;
                    rhs -= u[dim] * dphidx[dim];
                    // Compute the diffusive ionic term rhs
                    //rhs += higflow_diffusive_ionic_term(ns);
                    // Compute the potential ionic term rhs
                    //rhs += higflow_potential_ionic_term(ns);
                    //compute the term -4*(1-phi)3*(dphidx)*(dTdx)
                    //rhs += higflow_vol_frac_term1(ns);
                    //compute the second diffusive term: laplacian of the stress tensor
                    //rhs += higflow_vol_frac_term2(ns);
                }
                //rhs += higflow_vol_frac_term2(ns, varphi);
                break;
            // CUBISTA scheme
            case CELL_CUBISTA:
                for (int dim = 0; dim < DIM; dim++)
                {
                    // Set the computational cell
                    // higflow_computational_cell_electroosmotic_ionic(ns, sdnplus, sdpsi, sdphi, clid, ccenter, cdelta, dim, ns->ed.eo.dpnplus, ns->ed.eo.dppsi, ns->ed.eo.dpphi);
                    //higflow_computational_cell_volume_fraction_suspensions(ns, ns->ed.stsp.sdphi, clid, ccenter, cdelta, dim, ns->ed.stsp.dpphi);
                    //higflow_computational_cell_volume_fraction_suspensions(ns, ns->ed.stsp.sdphi, clid, ccenter, cdelta, dim, ns->ed.stsp.dpphi);
                    higflow_computational_cell_volume_fraction_suspensions(ns, sdphi, clid, ccenter, cdelta, dim, ns->ed.stsp.dpphi);
                    // Compute convective ionic term CUBISTA in rhs
                    rhs -= hig_flow_fraction_volume_suspensions_term_cubista(ns, ns->dpu[dim], ns->ed.stsp.sdphi, ns->ed.stn, varphi, ccenter, cdelta, dim);
                    //rhs    -= higflow_convective_ionic_term_cubista(ns, ns->dpu[dim], ns->ed.eo.dpnplus, ns->ed.eo.sdEOnplus, ns->ed.stn, nplus, ccenter, cdelta, dim);
                    //  Compute the diffusive ionic term rhs
                    //rhs += higflow_diffusive_ionic_term(ns);
                    // Compute the potential ionic term rhs
                    //rhs += higflow_potential_ionic_term(ns);
                    //compute the term -4*(1-phi)3*(dphidx)*(dTdx)
                    //rhs += higflow_vol_frac_term1(ns);
                    //rhs += 1.0;
                    //compute the second diffusive term: laplacian of the stress tensor
                }
                //rhs += higflow_vol_frac_term2(ns, varphi);
                break;
            }
            rhs += higflow_vol_frac_term2(ns, varphi);
            //real rhst = fabs(rhs);
            /*
            if(ns->par.t <= ns->par.dtp){
                //rhs = 0.0;
                if ((rhs < -1.0e16) ||  (rhs > 1.0e4)){
                    rhs = 0.0;
                } 
            }
            */
            //rhs = fabs(rhs);
            //if (rhs < 0.00001) rhs = 0.00001; 
            //DEBUG_INSPECT(rhs, %lf);
            // Compute the final value step time
            real newvarphi = varphi + ns->par.dt * rhs;
            if(newvarphi > 1.0){
                newvarphi = 1.0-eps;
            }
            //if(newvarphi > phij- frac_tol){
            //    newvarphi = phij- frac_tol;
            //}
            if(newvarphi < eps){
                newvarphi = eps;
            }
            //Check to see min and max volume fraction values
            if (newvarphi > volfracmax) volfracmax = newvarphi;
            if (newvarphi < volfracmin) volfracmin = newvarphi;
            // Set property value
            dp_set_value(ns->ed.stsp.dpphi, clid, newvarphi);
        }
        //Printing the min and max viscosity values
        printf("===> volfracmin = %lf <===> volfracmax = %lf <===\n", volfracmin, volfracmax);
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed ionic property
        dp_sync(ns->ed.stsp.dpphi);
    }
}