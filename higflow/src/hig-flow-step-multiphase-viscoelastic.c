// *******************************************************************
//  HiG-Flow Solver Step Multiphase - version 20/01/2022
// *******************************************************************

#include "hig-flow-step-multiphase-viscoelastic.h"


char nome_Frac_Visc[100];
char nome_Txx[100];
char nome_Txy[100];
char nome_Tyx[100];
char nome_Tyy[100];
char nome_d[100];     char nome_N[100];

void arquivoTempo_visc(int step) {
    sprintf(nome_Frac_Visc, "DATA/%d_frac_visc.txt", step);
    sprintf(nome_Txx, "DATA/%d_Txx.txt", step);
    sprintf(nome_Txy, "DATA/%d_Txy.txt", step);
    sprintf(nome_Tyx, "DATA/%d_Tyx.txt", step);
    sprintf(nome_Tyy, "DATA/%d_Tyy.txt", step);
    sprintf(nome_d, "DATA/%d_d.txt", step);
    sprintf(nome_N, "DATA/%d_norm.txt", step);

    FILE* fp = fopen(nome_Frac_Visc, "a"); fclose(fp);

    fp = fopen(nome_Txx, "a"); fclose(fp);

    fp = fopen(nome_Txy, "a"); fclose(fp);

    fp = fopen(nome_Tyx, "a"); fclose(fp);

    fp = fopen(nome_Tyy, "a"); fclose(fp);

    fp = fopen(nome_d, "a"); fclose(fp);

    fp = fopen(nome_N, "a"); fclose(fp);
}

void arquivo_Frac_Visc(char* nome, real x, real y, real frac) {
    FILE* fp = fopen(nome, "a");
    fprintf(fp, "%lf %lf %lf\n", x, y, frac);
    fclose(fp);
}

void arquivo_Frac_Ten(char* nome, real x, real y, real Ten) {
    FILE* fp = fopen(nome, "a");
    fprintf(fp, "%lf %lf %lf\n", x, y, Ten);
    fclose(fp);
}

void arquivo_d(char* nome, real d) {
    FILE* fp = fopen(nome, "a");
    fprintf(fp, "%lf\n", d);
    fclose(fp);
}

void arquivo_N(char* nome, real nx, real ny) {
    FILE* fp = fopen(nome, "a");
    fprintf(fp, "%lf %lf\n", nx, ny);
    fclose(fp);
}

void save_cell_values_visc(higflow_solver* ns, int aux) {
    // Get the local sub-domain for the cells
    sim_domain* sdm = psd_get_local_domain(ns->ed.mult.psdmult);
    // Get the local sub-domain for the facets
    sim_facet_domain* sfdu[DIM];
    for (int i = 0; i < DIM; i++) {
        sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
    }
    // Get the map for the domain properties
    mp_mapper* mp = sd_get_domain_mapper(sdm);
    // Loop for each cell
    higcit_celliterator* it;
    for (it = sd_get_domain_celliterator(sdm); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Get the cell
        hig_cell* c = higcit_getcell(it);
        // Get the cell identifier
        int clid = mp_lookup(mp, hig_get_cid(c));
        // Get the center of the cell
        Point ccenter;
        hig_get_center(c, ccenter);
        // Get the delta of the cell
        Point cdelta;
        hig_get_delta(c, cdelta);
        int infacet;
        if (aux == 1) {
            //saving volume fraction
            real fracvol = compute_value_at_point(sdm, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
            arquivo_Frac_Visc(nome_Frac_Visc, ccenter[0], ccenter[1], fracvol);

            real d = compute_value_at_point(sdm, ccenter, ccenter, 1.0, ns->ed.mult.dpdistance, ns->ed.mult.stn);
            arquivo_d(nome_d, d);

            Point Normal;
            Normal[0] = compute_value_at_point(sdm, ccenter, ccenter, 1.0, ns->ed.mult.dpnormal[0], ns->ed.mult.stn);
            Normal[1] = compute_value_at_point(sdm, ccenter, ccenter, 1.0, ns->ed.mult.dpnormal[1], ns->ed.mult.stn);
            arquivo_N(nome_N, Normal[0], Normal[1]);

            //saving Tensor
            real S[DIM][DIM];
            // Get Kernel
            S[0][0] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.ve.dpS[1][1], ns->ed.stn);
            S[0][1] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.ve.dpS[1][2], ns->ed.stn);
            S[1][0] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.ve.dpS[2][1], ns->ed.stn);
            S[1][1] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.ve.dpS[2][2], ns->ed.stn);

            arquivo_Frac_Ten(nome_Txx, ccenter[0], ccenter[1], S[0][0]);
            arquivo_Frac_Ten(nome_Txy, ccenter[0], ccenter[1], S[0][1]);
            arquivo_Frac_Ten(nome_Tyx, ccenter[0], ccenter[1], S[1][0]);
            arquivo_Frac_Ten(nome_Tyy, ccenter[0], ccenter[1], S[1][1]);
            //continue;
        }
    }
    higcit_destroy(it);
}

real higflow_interp_beta_multiphase_viscoelastic(real beta0, real beta1, real visc0, real visc1, real fracvol) {
    real visc = (1.0 - fracvol) * visc0 + fracvol * visc1;
    return ((1 - fracvol) * visc0 * beta0 + fracvol * visc1 * beta1) / visc;
}

real higflow_interp_De_multiphase_viscoelastic(real De0, real De1, real fracvol) {
    return (1.0 - fracvol) * De0 + fracvol * De1;
}

void higlow_interp_MM_multiphase_viscoelastic(real MM0[DIM][DIM], real MM1[DIM][DIM], real MM[DIM][DIM], real fracvol) {
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            MM[i][j] = (1.0 - fracvol) * MM0[i][j] + fracvol * MM1[i][j];
}

typedef enum phase_location {
    PHASE0 = 0,
    PHASE1 = 1,
    INTERFACE = 2
} phase_location;

// Constitutive Equation Step for the Explicit Euler Method
void higflow_explicit_euler_constitutive_equation_multiphase_viscoelastic(higflow_solver* ns) {
    // Get the cosntants
    real Re = ns->par.Re;
    real De;
    real De0 = ns->ed.mult.ve.par0.De;
    real De1 = ns->ed.mult.ve.par1.De;
    flow_type flowtype0 = ns->ed.mult.contr.flowtype0;
    flow_type flowtype1 = ns->ed.mult.contr.flowtype1;
    visc_model_type model0 = ns->ed.mult.ve.contr.model0;
    visc_model_type model1 = ns->ed.mult.ve.contr.model1;
   
    switch (model0) {
        case GPTT:
            ns->ed.mult.ve.par0.gamma_gptt = tgamma(ns->ed.mult.ve.par0.beta_gptt);
            break;
    }
    switch (model1) {
        case GPTT:
            ns->ed.mult.ve.par1.gamma_gptt = tgamma(ns->ed.mult.ve.par1.beta_gptt);
            break;
    }
    real tol         = ns->ed.mult.ve.par0.kernel_tol;
    real small       = 1.0e-14;
    // Get the local sub-domain for the cells
    sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
    // Get the local sub-domain for the facets
    sim_facet_domain *sfdu[DIM];
    for (int i = 0; i < DIM; i++) {
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
        int clid = mp_lookup(mp, hig_get_cid(c));
        // Get the center of the cell
        Point ccenter;
        hig_get_center(c, ccenter);
        // Get the delta of the cell
        Point cdelta;
        hig_get_delta(c, cdelta);

        // Get Volume fraction
        real fracvol = compute_value_at_point(ns->ed.mult.sdmult, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
        phase_location phase_loc;
        if (fracvol < small) phase_loc = PHASE0;
        else if (fracvol > 1.0 - small) phase_loc = PHASE1;
        else phase_loc = INTERFACE;

        // pure newtonian in this phase
        if ((phase_loc == PHASE0 && flowtype0 != VISCOELASTIC) ||
            (phase_loc == PHASE1 && flowtype1 != VISCOELASTIC)) {
            for (int i = 0; i < DIM; i++) {
                dp_set_value(ns->ed.mult.ve.dpS[i][i], clid, 1.0);
                for (int j = i + 1; j < DIM; j++) {
                    dp_set_value(ns->ed.mult.ve.dpS[i][j], clid, 0.0);
                    dp_set_value(ns->ed.mult.ve.dpS[j][i], clid, 0.0);
                }
            }
            continue;
        }

        // Get the velocity derivative tensor Du, S and Kernel tensor
        real Du[DIM][DIM], S[DIM][DIM], Kernel[DIM][DIM], KernelCopy[DIM][DIM];
        // Get the S tensor trace
        real trS = 0.0;
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                // Get Du
                Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.ve.dpD[i][j], ns->ed.stn);
                // Get S
                S[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.ve.dpS[i][j], ns->ed.stn);
                // Get Kernel
                Kernel[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.ve.dpKernel[i][j], ns->ed.stn);
                KernelCopy[i][j] = Kernel[i][j];
            }
            trS += S[i][i];
        }

        // Eige-values and eige-vectors of A
        real R[DIM][DIM], Klambda[DIM], lambda[DIM];
        hig_flow_jacobi(KernelCopy, Klambda, R);
        for (int dim = 0; dim < DIM; dim++)
            lambda[dim] = ns->ed.mult.ve.get_kernel_inverse(dim, Klambda[dim], tol);

        // Calculate M matrix >> M = R^t Du R
        real M[DIM][DIM];
        hig_flow_matrix_product(Du, R, M);
        // Calculate Omega matrix >> Omega = R Omega_aux R^t
        real Omega[DIM][DIM];
        hig_flow_calculate_omega(lambda, R, M, Omega, small);
        // Calculate the Kernel jacobian
        real jlambda[DIM];
        for(int i = 0; i < DIM; i++)
            jlambda[i] = ns->ed.mult.ve.get_kernel_jacobian(i, lambda[i], tol);
        // Calculate the matrix BB and the matrix B
        real BB[DIM][DIM]; real B[DIM][DIM];
        hig_flow_calculate_bs (lambda, jlambda, R, M, BB, B);
        // Calculate the matrix MM for the model
        real MM0[DIM][DIM], MM1[DIM][DIM];
        real MM[DIM][DIM], M_aux[DIM][DIM];
        if (phase_loc == PHASE0 || (phase_loc == INTERFACE && flowtype0 == VISCOELASTIC)) {
            // Phase 0
            switch (model0) {
                case USERSET: // User Model
                    ns->ed.mult.ve.calculate_m_user_multiphase(fracvol, lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par0), &(ns->ed.mult.ve.par1));
                    break;
                case OLDROYD_B: // Oldroyd-B Model
                    hig_flow_calculate_m_oldroyd(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par0));
                    break;
                case GIESEKUS: // Giesekus Model
                    hig_flow_calculate_m_giesekus(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par0));
                    break;
                case LPTT: // LPTT Model
                    hig_flow_calculate_m_lptt(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par0));
                    break;
                case GPTT: // GPTT Model
                    hig_flow_calculate_m_gptt(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par0));
                    break;
                case FENE_P: // FENE-P Model
                    hig_flow_calculate_m_fene_p(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par0));
                    break;
                case E_FENE: // e-FENE Model
                    hig_flow_calculate_m_e_fene(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par0));
                    break;
                }
            // Calculate Kernel matrix >> MM = R M(Lambda) JLambda R^t
            hig_flow_matrix_transpose_product(M_aux, R, MM0);
        } else {
            for (int i = 0; i < DIM; i++)
                for (int j = 0; j < DIM; j++)
                    MM0[i][j] = 0.0;
        }
        if (phase_loc == PHASE1 || (phase_loc == INTERFACE && flowtype1 == VISCOELASTIC)) {
            // Phase 1
            switch (model1) {
                case USERSET: // User Model
                    ns->ed.mult.ve.calculate_m_user_multiphase(fracvol, lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par0), &(ns->ed.mult.ve.par1));
                    break;
                case OLDROYD_B: // Oldroyd-B Model
                    hig_flow_calculate_m_oldroyd(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par1));
                    break;
                case GIESEKUS: // Giesekus Model
                    hig_flow_calculate_m_giesekus(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par1));
                    break;
                case LPTT: // LPTT Model
                    hig_flow_calculate_m_lptt(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par1));
                    break;
                case GPTT: // GPTT Model
                    hig_flow_calculate_m_gptt(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par1));
                    break;
                case FENE_P: // FENE-P Model
                    hig_flow_calculate_m_fene_p(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par1));
                    break;
                case E_FENE: // e-FENE Model
                    hig_flow_calculate_m_e_fene(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par1));
                    break;
            }
            // Calculate Kernel matrix >> MM = R M(Lambda) JLambda R^t
            hig_flow_matrix_transpose_product(M_aux, R, MM1);
        } else {
            for (int i = 0; i < DIM; i++)
                for (int j = 0; j < DIM; j++)
                    MM1[i][j] = 0.0;
        }

        // Calculate RHS = Omega Kernel - Kernel Omega + 2BB + MM/De
        real RHS[DIM][DIM];

        De = higflow_interp_De_multiphase_viscoelastic(De0, De1, fracvol);
        higlow_interp_MM_multiphase_viscoelastic(MM0, MM1, MM, fracvol);
        hig_flow_kernel_rhs(De, Kernel, Omega, BB, MM, RHS);

        // Get the velocity at cell center 
        real u[DIM], dKdx[DIM];
        hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
        // Solving the Constitutive Equation using the Euler Method
        for (int i = 0; i < DIM; i++) {
            for (int j = i; j < DIM; j++) {
                // Right hand side equation
                real rhs = 0.0;
                switch (ns->ed.mult.ve.contr.convecdiscrtype) {
                case CELL_CENTRAL:
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
                        rhs -= hig_flow_convective_cell_term_cubista(ns->dpu[dim], ns->sfdu[dim], ns->stn, ns->ed.mult.ve.dpKernel[i][j], ns->ed.sdED, ns->ed.stn, Kernel[i][j], ccenter, cdelta, dim);
                    }
                    break;
                }
                // Compute the final rhs
                rhs         += RHS[i][j];
                // Compute the Kernel at next time
                real kernel = Kernel[i][j] + ns->par.dt * rhs;
                // if (phase_loc == INTERFACE) { // interpolate with identity
                //     if (flowtype0 != VISCOELASTIC) kernel = (1.0 - fracvol) * (i == j) + fracvol * kernel;
                //     if (flowtype1 != VISCOELASTIC) kernel = (1.0 - fracvol) * kernel + fracvol * (i == j);
                // }
                // Store the Kernel
                dp_set_value(ns->ed.mult.ve.dpS[i][j], clid, kernel);
                if (i != j) {
                    dp_set_value(ns->ed.mult.ve.dpS[j][i], clid, kernel);
                }
            }
        }
    }
    // Destroy the iterator
    higcit_destroy(it);
    // Sync the distributed pressure property
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            dp_sync(ns->ed.mult.ve.dpS[i][j]);
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
                real S = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.ve.dpS[i][j], ns->ed.stn);

                if (j>=i) UPDATE_RESIDUAL_BUFFER_CELL(ns, dp_get_value(ns->ed.mult.ve.dpKernel[i][j], clid), S, c, ccenter)

                    // Store Kernel
                    dp_set_value(ns->ed.mult.ve.dpKernel[i][j], clid, S);
            }
            // Destroy the iterator
            higcit_destroy(it);

            if (j>=i) UPDATE_RESIDUALS(ns, ns->residuals->Kernel[i][j])
        }
    }

    // Sync the distributed pressure property
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            dp_sync(ns->ed.mult.ve.dpKernel[i][j]);
        }
    }
}


// Constitutive Equation Step for the Implicit Euler Method
void higflow_implicit_euler_constitutive_equation_multiphase_viscoelastic(higflow_solver* ns) {
    // Get the cosntants
    real Re = ns->par.Re;
    real De;
    real De0 = ns->ed.mult.ve.par0.De;
    real De1 = ns->ed.mult.ve.par1.De;
    flow_type flowtype0 = ns->ed.mult.contr.flowtype0;
    flow_type flowtype1 = ns->ed.mult.contr.flowtype1;
    visc_model_type model0 = ns->ed.mult.ve.contr.model0;
    visc_model_type model1 = ns->ed.mult.ve.contr.model1;
   
    switch (model0) {
        case GPTT:
            ns->ed.mult.ve.par0.gamma_gptt = tgamma(ns->ed.mult.ve.par0.beta_gptt);
            break;
    }
    switch (model1) {
        case GPTT:
            ns->ed.mult.ve.par1.gamma_gptt = tgamma(ns->ed.mult.ve.par1.beta_gptt);
            break;
    }
    real tol         = ns->ed.mult.ve.par0.kernel_tol;
    real small       = 1.0e-14;
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

        // Get Volume fraction
        real fracvol = compute_value_at_point(ns->ed.mult.sdmult, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
        phase_location phase_loc;
        if (fracvol < small) phase_loc = PHASE0;
        else if (fracvol > 1.0 - small) phase_loc = PHASE1;
        else phase_loc = INTERFACE;

        // pure newtonian in this phase
        if ((phase_loc == PHASE0 && flowtype0 != VISCOELASTIC) ||
            (phase_loc == PHASE1 && flowtype1 != VISCOELASTIC)) {
            for (int i = 0; i < DIM; i++) {
                dp_set_value(ns->ed.mult.ve.dpS[i][i], clid, 1.0);
                for (int j = i + 1; j < DIM; j++) {
                    dp_set_value(ns->ed.mult.ve.dpS[i][j], clid, 0.0);
                    dp_set_value(ns->ed.mult.ve.dpS[j][i], clid, 0.0);
                }
            }
            continue;
        }

        // Get the velocity derivative tensor Du, S and Kernel tensor
        real Du[DIM][DIM], S[DIM][DIM], Kernel[DIM][DIM], KernelCopy[DIM][DIM];
        // Get the S tensor trace
        real trS = 0.0;
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                // Get Du
                Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.ve.dpD[i][j], ns->ed.stn);
                // Get S
                S[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.ve.dpS[i][j], ns->ed.stn);
                // Get Kernel
                Kernel[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.ve.dpKernel[i][j], ns->ed.stn);
                KernelCopy[i][j] = Kernel[i][j];
            }
            trS += S[i][i];
        }

        // Eige-values and eige-vectors of A
        real R[DIM][DIM], Klambda[DIM], lambda[DIM];
        hig_flow_jacobi(KernelCopy, Klambda, R);
        for (int dim = 0; dim < DIM; dim++)
            lambda[dim] = ns->ed.mult.ve.get_kernel_inverse(dim, Klambda[dim], tol);

        // Calculate M matrix >> M = R^t Du R
        real M[DIM][DIM];
        hig_flow_matrix_product(Du, R, M);
        // Calculate Omega matrix >> Omega = R Omega_aux R^t
        real Omega[DIM][DIM];
        hig_flow_calculate_omega(lambda, R, M, Omega, small);
        // Calculate the Kernel jacobian
        real jlambda[DIM];
        for(int i = 0; i < DIM; i++)
            jlambda[i] = ns->ed.mult.ve.get_kernel_jacobian(i, lambda[i], tol);
        // Calculate the matrix BB and the matrix B
        real BB[DIM][DIM]; real B[DIM][DIM];
        hig_flow_calculate_bs (lambda, jlambda, R, M, BB, B);
        // Calculate the matrix MM for the model
        real MM0[DIM][DIM], MM1[DIM][DIM];
        real MM[DIM][DIM], M_aux[DIM][DIM];
        if (phase_loc == PHASE0 || (phase_loc == INTERFACE && flowtype0 == VISCOELASTIC)) {
            // Phase 0
            switch (model0) {
                case USERSET: // User Model
                    ns->ed.mult.ve.calculate_m_user_multiphase(fracvol, lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par0), &(ns->ed.mult.ve.par1));
                    break;
                case OLDROYD_B: // Oldroyd-B Model
                    hig_flow_calculate_m_oldroyd(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par0));
                    break;
                case GIESEKUS: // Giesekus Model
                    hig_flow_calculate_m_giesekus(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par0));
                    break;
                case LPTT: // LPTT Model
                    hig_flow_calculate_m_lptt(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par0));
                    break;
                case GPTT: // GPTT Model
                    hig_flow_calculate_m_gptt(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par0));
                    break;
                case FENE_P: // FENE-P Model
                    hig_flow_calculate_m_fene_p(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par0));
                    break;
                case E_FENE: // e-FENE Model
                    hig_flow_calculate_m_e_fene(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par0));
                    break;
                }
            // Calculate Kernel matrix >> MM = R M(Lambda) JLambda R^t
            hig_flow_matrix_transpose_product(M_aux, R, MM0);
        } else {
            for (int i = 0; i < DIM; i++)
                for (int j = 0; j < DIM; j++)
                    MM0[i][j] = 0.0;
        }
        if (phase_loc == PHASE1 || (phase_loc == INTERFACE && flowtype1 == VISCOELASTIC)) {
            // Phase 1
            switch (model1) {
                case USERSET: // User Model
                    ns->ed.mult.ve.calculate_m_user_multiphase(fracvol, lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par0), &(ns->ed.mult.ve.par1));
                    break;
                case OLDROYD_B: // Oldroyd-B Model
                    hig_flow_calculate_m_oldroyd(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par1));
                    break;
                case GIESEKUS: // Giesekus Model
                    hig_flow_calculate_m_giesekus(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par1));
                    break;
                case LPTT: // LPTT Model
                    hig_flow_calculate_m_lptt(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par1));
                    break;
                case GPTT: // GPTT Model
                    hig_flow_calculate_m_gptt(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par1));
                    break;
                case FENE_P: // FENE-P Model
                    hig_flow_calculate_m_fene_p(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par1));
                    break;
                case E_FENE: // e-FENE Model
                    hig_flow_calculate_m_e_fene(lambda, jlambda, B, M_aux, Re, trS, &(ns->ed.mult.ve.par1));
                    break;
            }
            // Calculate Kernel matrix >> MM = R M(Lambda) JLambda R^t
            hig_flow_matrix_transpose_product(M_aux, R, MM1);
        } else {
            for (int i = 0; i < DIM; i++)
                for (int j = 0; j < DIM; j++)
                    MM1[i][j] = 0.0;
        }

        // Calculate RHS = Omega Kernel - Kernel Omega + 2BB + MM/De
        real RHS[DIM][DIM];

        De = higflow_interp_De_multiphase_viscoelastic(De0, De1, fracvol);
        higlow_interp_MM_multiphase_viscoelastic(MM0, MM1, MM, fracvol);
        hig_flow_implicit_kernel_rhs(De, BB, MM, RHS);

        // Get the velocity at cell center 
        real u[DIM], dKdx[DIM];
        hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
        real b[DIM*DIM];
        real w[DIM*DIM][DIM*DIM+1]; 
        // Solving the Constitutive Equation using the Euler Method
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                // Right hand side equation
                real rhs = 0.0;
                switch (ns->ed.mult.ve.contr.convecdiscrtype) {
                case CELL_CENTRAL:
                    // Kernel derivative at cell center
                    hig_flow_derivative_kernel_at_center_cell(ns, ccenter, cdelta, i, j, Kernel[i][j], dKdx);
                    for (int dim = 0; dim < DIM; dim++) {
                        //Compute convective tensor term in rhs
                        rhs -= u[dim] * dKdx[dim];
                    }
                    break;
                case CELL_CUBISTA:
                    //Compute convective tensor term CUBISTA in rhs
                    for (int dim = 0; dim < DIM; dim++) {
                        rhs -= hig_flow_convective_cell_term_cubista(ns->dpu[dim], ns->sfdu[dim], ns->stn, ns->ed.mult.ve.dpKernel[i][j], ns->ed.sdED, ns->ed.stn, Kernel[i][j], ccenter, cdelta, dim);
                    }
                    break;
                }
                // Compute the final rhs
                rhs += RHS[i][j];
                rhs         *= ns->par.dt;
                rhs         += Kernel[i][j];
                b[i*DIM + j] = rhs;
            }
        }
        //Calculate de kronecker product Omega*I - I*Omega
        hig_flow_kernel_system_matrix (w, Omega, ns->par.dt);
        //Solve the linear system
        hig_flow_solve_system_constitutive_equation (DIM*DIM, w, b);
        // Get the solution of linear system
        for (int i = 0; i < DIM; i++) {
            for (int j = i; j < DIM; j++) {
                // Get the value of kernel
                real kernel = b[i*DIM+j];
                // if (phase_loc == INTERFACE) { // interpolate with identity
                //     if (flowtype0 != VISCOELASTIC) kernel = (1.0 - fracvol) * (i == j) + fracvol * kernel;
                //     if (flowtype1 != VISCOELASTIC) kernel = (1.0 - fracvol) * kernel + fracvol * (i == j);
                // }
                // Set the value of kernel
                dp_set_value(ns->ed.mult.ve.dpS[i][j], clid, kernel);
                if (i != j) {
                    dp_set_value(ns->ed.mult.ve.dpS[j][i], clid, kernel);
                }
            }
        }  
    }
    // Destroy the iterator
    higcit_destroy(it);
    // Sync the distributed pressure property
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            dp_sync(ns->ed.mult.ve.dpS[i][j]);
        }
    }
    // Store the Kernel Tensor
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                // Get the cell
                hig_cell* c = higcit_getcell(it);
                // Get the cell identifier
                int clid = mp_lookup(mp, hig_get_cid(c));
                // Get the center of the cell
                Point ccenter;
                hig_get_center(c, ccenter);
                // Get the S tensor and store in Kernel
                real S = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.ve.dpS[i][j], ns->ed.stn);

                if (j >= i) UPDATE_RESIDUAL_BUFFER_CELL(ns, dp_get_value(ns->ed.mult.ve.dpKernel[i][j], clid), S, c, ccenter)

                    // Store Kernel
                    dp_set_value(ns->ed.mult.ve.dpKernel[i][j], clid, S);
            }
            // Destroy the iterator
            higcit_destroy(it);

            if (j >= i) UPDATE_RESIDUALS(ns, ns->residuals->Kernel[i][j])
        }
    }

    // Sync the distributed pressure property
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            dp_sync(ns->ed.mult.ve.dpKernel[i][j]);
        }
    }
}

real higflow_interp_fA_multiphase_viscoelastic(real fA0, real fA1, real fracvol) {
    return (1.0 - fracvol) * fA0 + fracvol * fA1;
}

real higflow_interp_a_multiphase_viscoelastic(real a0, real a1, real fracvol) {
    return (1.0 - fracvol) * a0 + fracvol * a1;
}


// Computing the Polymeric Tensor
void higflow_compute_polymeric_tensor_multiphase_viscoelastic(higflow_solver* ns) {
    if (ns->contr.flowtype == VISCOELASTIC) {
        // Get the constants
        real Re   = ns->par.Re;
        real De, De0, De1;
        real beta, beta0, beta1, visc0, visc1;
        real fA, fA0, fA1;
        real a, a0, a1;
        real trA;
        int flowtype0 = ns->ed.mult.contr.flowtype0;
        int flowtype1 = ns->ed.mult.contr.flowtype1;
        real tol = ns->ed.mult.ve.par0.kernel_tol;
        real Tmax[DIM][DIM], Tmin[DIM][DIM];
        for (int i = 0; i < DIM; i++) {
           for (int j = 0; j < DIM; j++) {
               Tmax[i][j] = -1.0e16;
               Tmin[i][j] =  1.0e16;
           }
        }
        real small = 1.0e-14;
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
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);

            // Get Volume fraction
            real fracvol = compute_value_at_point(ns->ed.mult.sdmult, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
            phase_location phase_loc;
            if (fracvol < small) phase_loc = PHASE0;
            else if (fracvol > 1.0 - small) phase_loc = PHASE1;
            else phase_loc = INTERFACE;

            // Get the velocity derivative tensor Du and the Kernel tensor
            real Kernel[DIM][DIM], Du[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get Du
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.ve.dpD[i][j], ns->ed.stn);
                    // Get Kernel
                    Kernel[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.ve.dpKernel[i][j], ns->ed.stn);
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
                B[i][i] = ns->ed.mult.ve.get_kernel_inverse(i, lambda[i], tol);
            }
            // Calculate A matrix >> A = R B R^t
            // D = 0.5*(Du + Du^t)
            real A[DIM][DIM], D[DIM][DIM];
            hig_flow_matrix_transpose_product(B, R, A);
            // Calculate the tensor S
            
            if(flowtype0 == VISCOELASTIC) {
                switch (ns->ed.mult.ve.contr.model0) {
                    case 4: ;///////////////////////////////////// FENE-P
                        real L2 = ns->ed.mult.ve.par0.L2_fene;
                        trA = 0.0;
                        for (int i = 0; i < DIM; i++)
                            trA += A[i][i];
                        fA0 = L2/(L2 - trA);
                        a0 = L2/(L2-3);
                        break;
                    case 5: ;///////////////////////////////////// e-FENE
                        real b_fene = ns->ed.mult.ve.par0.L2_fene;
                        real lambda_fene = ns->ed.mult.ve.par0.lambda_fene;
                        real E = ns->ed.mult.ve.par0.E_fene;
                        trA = 0.0;
                        for (int i = 0; i < DIM; i++)
                            trA += A[i][i];
                        fA0 = b_fene/(b_fene-trA) - E*sqrt(b_fene)*exp(-sqrt(trA)/lambda_fene)*(1.0/(trA*lambda_fene)+1/(trA*sqrt(trA)));
                        a0  = 1.0;
                    default: ///////////////////////////////////// Outros
                        fA0 = 1.0;
                        a0  = 1.0;
                        break;
                } 
            } else{
                fA0 = 1.0;
                a0  = 1.0;
            }
            if(flowtype1 == VISCOELASTIC) {
                switch (ns->ed.mult.ve.contr.model1) {
                    case 4: ;///////////////////////////////////// FENE-P
                        real L2 = ns->ed.mult.ve.par1.L2_fene;
                        trA = 0.0;
                        for (int i = 0; i < DIM; i++)
                            trA += A[i][i];
                        fA1 = L2/(L2 - trA);
                        a1 = L2/(L2-3);
                        break;
                    case 5: ;///////////////////////////////////// e-FENE
                        real b_fene = ns->ed.mult.ve.par1.L2_fene;
                        real lambda_fene = ns->ed.mult.ve.par1.lambda_fene;
                        real E = ns->ed.mult.ve.par1.E_fene;
                        trA = 0.0;
                        for (int i = 0; i < DIM; i++)
                            trA += A[i][i];
                        fA1 = b_fene/(b_fene-trA) - E*sqrt(b_fene)*exp(-sqrt(trA)/lambda_fene)*(1.0/(trA*lambda_fene)+1/(trA*sqrt(trA)));
                        a1  = 1.0;
                    default: ///////////////////////////////////// Outros
                        fA1 = 1.0;
                        a1  = 1.0;
                        break;
                } 
            } else{
                fA1 = 1.0;
                a1  = 1.0;
            }

            // interpolate parameters
            fA = higflow_interp_fA_multiphase_viscoelastic(fA0, fA1, fracvol);
            a = higflow_interp_a_multiphase_viscoelastic(a0, a1, fracvol);
            De = higflow_interp_De_multiphase_viscoelastic(De0, De1, fracvol);
            visc0 = ns->ed.mult.get_viscosity0(ccenter, ns->par.t);
            visc1 = ns->ed.mult.get_viscosity1(ccenter, ns->par.t);
            beta = higflow_interp_beta_multiphase_viscoelastic(beta0, beta1, visc0, visc1, fracvol);

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
                   dp_set_value(ns->ed.mult.ve.dpS[i][j], clid, S[i][j]);
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
                dp_sync(ns->ed.mult.ve.dpS[i][j]);
            }
        }
    }
}


// // Set velocity for the test: Zalesak and single vortex
// void higflow_set_velocity_multiphase(higflow_solver *ns) {
//     // Get the local sub-domain
//     sim_domain *sdp = psd_get_local_domain(ns->psdp);
//     sim_facet_domain *sfdu[DIM];
//     // Loop for each dimension
//     higfit_facetiterator *fit;
//     //real absvelmin =  1.0e16;
//     for (int dim = 0; dim < DIM; dim++) {
//         // Get the local partitioned domain for facets
//         sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
//         // Get the map of the distributed properties in the facets
//         mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
//         // Loop for each facet
//         for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
//             // Get the facet
//             hig_facet *f = higfit_getfacet(fit);
//             int flid = mp_lookup(mu, hig_get_fid(f));
//             // Get the center of the facet
//             Point fcenter;
//             hig_get_facet_center(f, fcenter);
//             // Get the delta of the facet
//             Point fdelta;
//             hig_get_facet_delta(f, fdelta);
//             // Get the velocity
//             real uset, pii;
//             pii = 3.14159265359;
//             if (dim == 0) {
//                // Zalesak
//                uset = 0.5 - fcenter[1];
//                // Single vortex
//                /*if (ns->par.t < 1.0) {
//                   uset = -2.0*pow(sin(pii*fcenter[0]),2)*sin(pii*fcenter[1])*cos(pii*fcenter[1]);
//                } else {
//                   uset = 2.0*pow(sin(pii*fcenter[0]),2)*sin(pii*fcenter[1])*cos(pii*fcenter[1]);
//                }*/
//             } else if (dim == 1) {
//                // Zalesak
//                uset =  fcenter[0] - 0.5;
//                // Single vortex
//                /*if (ns->par.t < 1.0) {
//                   uset = 2.0*pow(sin(pii*fcenter[1]),2)*sin(pii*fcenter[0])*cos(pii*fcenter[0]);
//                } else {
//                   uset = -2.0*pow(sin(pii*fcenter[1]),2)*sin(pii*fcenter[0])*cos(pii*fcenter[0]);
//                }*/           
//             }
//             // Set the final velocity in the distributed velocity property
//             dp_set_value(ns->dpu[dim], flid, uset);
//         }
//         // Destroy the iterator
//         higfit_destroy(fit);
//         // Sync the distributed velocity property
//         dp_sync(ns->dpu[dim]);
//     }
// }


// One step of the Navier-Stokes the projection method
void higflow_solver_step_multiphase_viscoelastic(higflow_solver* ns) {
    // Calculate the velocity derivative tensor
    higflow_compute_velocity_derivative_tensor(ns);
    // Constitutive Equation Step for the Explicit Euler Method
    switch (ns->ed.mult.ve.contr.discrtype) {
        case EXPLICIT:
            // Explicit method
            higflow_explicit_euler_constitutive_equation_multiphase_viscoelastic(ns);
            break;
        case IMPLICIT:
            // Implicit method
            higflow_implicit_euler_constitutive_equation_multiphase_viscoelastic(ns);
            break;
        }
    // Computing the Polymeric Tensor
    higflow_compute_polymeric_tensor_multiphase_viscoelastic(ns);

    // Boundary condition for velocity
    higflow_boundary_condition_for_velocity(ns);
    // Boundary conditions for source term
    higflow_boundary_condition_for_cell_source_term(ns);
    higflow_boundary_condition_for_facet_source_term(ns);
    // Calculate the source term
    higflow_calculate_source_term(ns);
    // Calculate the facet source term
    higflow_calculate_facet_source_term(ns);
    // Boundary condition for pressure
    higflow_boundary_condition_for_pressure(ns);
    // Calculate the viscosity
    higflow_compute_viscosity_multiphase(ns);
    // Calculate the viscosity
    higflow_compute_density_multiphase(ns);
    // Calculate the curvature
    higflow_compute_curvature_multiphase(ns);
    // Calculate the intermediated velocity
    switch (ns->contr.tempdiscrtype) {
        case EXPLICIT_EULER:
            // Explicit Euler method
            higflow_explicit_euler_intermediate_velocity_multiphase(ns, ns->dpu, ns->dpustar);
            break;
        case EXPLICIT_RK2:
            // Explicit RK2 method
            higflow_explicit_runge_kutta_2_intermediate_velocity_multiphase(ns);
            break;
        case EXPLICIT_RK3:
            // Explicit RK3 method
            higflow_explicit_runge_kutta_3_intermediate_velocity_multiphase(ns);
            break;
        case SEMI_IMPLICIT_EULER:
            // Semi-Implicit Euler Method
            higflow_semi_implicit_euler_intermediate_velocity_multiphase(ns);
            break;
        case SEMI_IMPLICIT_CN:
            // Semi-Implicit Crank-Nicolson Method
            higflow_semi_implicit_crank_nicolson_intermediate_velocity_multiphase(ns);
            break;
        case SEMI_IMPLICIT_BDF2:
            // Semi-Implicit Crank-Nicolson Method
            higflow_semi_implicit_bdf2_intermediate_velocity_multiphase(ns, ns->dpu, ns->dpustar);
            break;
    }
    // Calculate the pressure
    higflow_pressure_multiphase(ns);
    // Calculate the final velocity
    higflow_final_velocity_multiphase(ns);
    // Calculate the final pressure
    higflow_final_pressure(ns);
    // if (ns->par.stepaux % 20000 == 0) {
    //     printf("creating archives at step: %d\n", ns->par.stepaux);
    //     arquivoTempo_visc(ns->par.stepaux);
    //     save_cell_values_visc(ns, 1);
    // }
    // Calculate the volume fraction
    higflow_plic_advection_volume_fraction(ns);
    //higflow_explicit_euler_volume_fraction(ns);
}
