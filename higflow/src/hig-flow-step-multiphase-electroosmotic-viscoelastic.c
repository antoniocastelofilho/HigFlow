// *******************************************************************
//  HiG-Flow Solver Step Electro-osmotic - version 12/05/2024
// *******************************************************************

#include "hig-flow-step-multiphase-electroosmotic-viscoelastic.h"

void higflow_solver_step_multiphase_electroosmotic_viscoelastic(higflow_solver *ns) {
    // Calculate the deformation rate tensor
    higflow_compute_velocity_derivative_tensor(ns);
    // Constitutive Equation Step for viscoelastic flow
    switch (ns->ed.mult.ve.contr.discrtype) {
        case EXPLICIT:
            higflow_explicit_euler_constitutive_equation_multiphase_viscoelastic(ns);
            break;
        case IMPLICIT:
            higflow_implicit_euler_constitutive_equation_multiphase_viscoelastic(ns);
            break;
    }
    // Calculate the elastic tensor to be used in the momentum equation
    higflow_compute_polymeric_tensor_multiphase_viscoelastic(ns);

    // Boundary conditions and source terms
    higflow_boundary_condition_for_velocity(ns);
    higflow_boundary_condition_for_pressure(ns);
    higflow_boundary_condition_for_cell_source_term(ns);
    higflow_boundary_condition_for_facet_source_term(ns);
    higflow_calculate_source_term(ns);
    higflow_calculate_facet_source_term(ns);

    // Calculate the electro-osmotic terms
    higflow_boundary_condition_for_phi(ns);
    higflow_boundary_condition_for_psi(ns);
    switch (ns->ed.mult.eo.contr.eo_model) {
    case PNP: // Poisson-Nernst-Planck model
        for(int k=0; k<3; k++) {
            higflow_boundary_condition_for_electroosmotic_nplus(ns);
            higflow_boundary_condition_for_electroosmotic_nminus(ns);
            if( (ns->par.step == 0) || ns->ed.mult.eo.contr.is_phibc_timedependent == true)
            higflow_electroosmotic_phi(ns);
            switch (ns->ed.mult.eo.contr.tempdiscrtype) {
                case EXPLICIT_EULER:
                    higflow_explicit_euler_ionic_transport_equation_nplus_multiphase(ns);
                    higflow_explicit_euler_ionic_transport_equation_nminus_multiphase(ns);
                    break;
                case SEMI_IMPLICIT_EULER:
                    higflow_semi_implicit_euler_ionic_transport_equation_nplus_multiphase(ns);
                    higflow_semi_implicit_euler_ionic_transport_equation_nminus_multiphase(ns);
                    break;
            }
            real max_psi_res = higflow_multiphase_electroosmotic_psi(ns);
            if(k>0) print0f("=+=+=+ psi residual in inner iteration %d: %15.10lf =+=+=+\n", k+1, max_psi_res);
        }
        dp_copy_values(ns->ed.eo.dpnplus, ns->ed.eo.dpnplus_temp);
        dp_copy_values(ns->ed.eo.dpnminus, ns->ed.eo.dpnminus_temp);
        higflow_calculate_multiphase_electroosmotic_source_term(ns);
        break;
    case PB: // Poisson-Boltzmann model 
        if( (ns->par.step == 0) || ns->ed.mult.eo.contr.is_phibc_timedependent == true)
        higflow_electroosmotic_phi(ns);
        if( (ns->par.step == 0) || ns->ed.eo.contr.is_psibc_timedependent == true)
        higflow_multiphase_electroosmotic_solve_pb(ns);
        if( (ns->par.step == 0) || ns->ed.mult.eo.contr.is_phibc_timedependent == true
                              || ns->ed.mult.eo.contr.is_psibc_timedependent == true)
        higflow_calculate_multiphase_electroosmotic_source_term(ns);
        break;
    case PBDH: // Debye-Hückel model (solving poisson equation for psi) 
        if( (ns->par.step == 0) || ns->ed.mult.eo.contr.is_phibc_timedependent == true)
        higflow_electroosmotic_phi(ns);
        if( (ns->par.step == 0) || ns->ed.mult.eo.contr.is_psibc_timedependent == true)
            higflow_multiphase_electroosmotic_psi(ns);
        if( (ns->par.step == 0) || ns->ed.mult.eo.contr.is_phibc_timedependent == true
                              || ns->ed.mult.eo.contr.is_psibc_timedependent == true)
            higflow_calculate_multiphase_electroosmotic_source_term(ns);
        break;
    case PBDH_ANALYTIC:
        printf("PBDH_ANALYTIC not implemented\n");
        break;
    }

    // Interpolate the viscosity and density
    higflow_compute_viscosity_multiphase(ns);
    higflow_compute_density_multiphase(ns);
    // Calculate the curvature, interfacial force and normal
    higflow_compute_curvature_interfacial_force_normal_multiphase_2D_hf_shirani(ns);
    higflow_compute_distance_multiphase_2D(ns);
    higflow_compute_plic_lines_2d(ns);

    // Calculate the intermediate velocity
    switch (ns->contr.tempdiscrtype) {
    case EXPLICIT_EULER:
        higflow_explicit_euler_intermediate_velocity_multiphase_electroosmotic(ns, ns->dpu, ns->dpustar);
        break;
    case EXPLICIT_RK2:
        higflow_explicit_runge_kutta_2_intermediate_velocity_multiphase_electroosmotic(ns);
        break;
    case EXPLICIT_RK3:
        higflow_explicit_runge_kutta_3_intermediate_velocity_multiphase_electroosmotic(ns);
        break;
    case SEMI_IMPLICIT_EULER:
        higflow_semi_implicit_euler_intermediate_velocity_multiphase_electroosmotic(ns);
        break;
    case SEMI_IMPLICIT_CN:
        higflow_semi_implicit_crank_nicolson_intermediate_velocity_multiphase_electroosmotic(ns);
        break;
    case SEMI_IMPLICIT_BDF2:
        higflow_semi_implicit_bdf2_intermediate_velocity_multiphase_electroosmotic(ns);
        break;
    }

    // Projection
    higflow_pressure_multiphase(ns);
    higflow_final_velocity_multiphase(ns);
    higflow_final_pressure(ns);

    // Calculate the volume fraction
    higflow_plic_advection_volume_fraction(ns);

}
