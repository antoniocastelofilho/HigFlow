// *******************************************************************
//  HiG-Flow Solver Step Electro-osmotic - version 12/05/2024
// *******************************************************************

#include "hig-flow-step-multiphase-electroosmotic-viscoelastic.h"

void higflow_solver_step_multiphase_electroosmotic_viscoelastic(higflow_solver *ns) {
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
    // Boundary condition for n+ 
    higflow_boundary_condition_for_electroosmotic_nplus(ns);
    // Boundary condition for n- 
    higflow_boundary_condition_for_electroosmotic_nminus(ns);
    // Boundary condition for pressure
    higflow_boundary_condition_for_pressure(ns);
    // Calculate the source term
    higflow_calculate_source_term(ns);
    // Calculate the facet source term
    higflow_calculate_facet_source_term(ns);
    // Calculate the electro-osmotic source term
    switch (ns->ed.mult.eo.contr.eo_model) {
    case PNP:
        // Poisson-Nernst-Planck model
        if( (ns->par.step == 0) || ns->ed.mult.eo.contr.is_phibc_timedependent == true) {
            higflow_boundary_condition_for_phi(ns);
            higflow_electroosmotic_phi(ns);
        }
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
        higflow_boundary_condition_for_psi(ns);
        higflow_multiphase_electroosmotic_psi(ns);
        higflow_calculate_multiphase_electroosmotic_source_term(ns);
        break;
    case PB:
        // Poisson-Boltzmann model 
        if( (ns->par.step == 0) || ns->ed.mult.eo.contr.is_phibc_timedependent == true) {
            higflow_boundary_condition_for_phi(ns);
            higflow_electroosmotic_phi(ns);
        }
        higflow_boundary_condition_for_psi(ns);
        higflow_multiphase_electroosmotic_psi(ns);
        higflow_calculate_multiphase_electroosmotic_source_term(ns);
        break;
    case PBDH:
        // Debye-Hückel model (solving poisson equation for psi) 
        // Poisson-Boltzmann model 
        if( (ns->par.step == 0) || ns->ed.mult.eo.contr.is_phibc_timedependent == true) {
            higflow_boundary_condition_for_phi(ns);
            higflow_electroosmotic_phi(ns);
        }
        if( (ns->par.step == 0) || ns->ed.mult.eo.contr.is_psibc_timedependent == true) {
            higflow_boundary_condition_for_psi(ns);
            higflow_multiphase_electroosmotic_psi(ns);
        }
        if( (ns->par.step == 0) || ns->ed.mult.eo.contr.is_phibc_timedependent == true
                              || ns->ed.mult.eo.contr.is_psibc_timedependent == true) {
            higflow_calculate_multiphase_electroosmotic_source_term(ns);
        }
        break;
    case PBDH_ANALYTIC:
        printf("PBDH_ANALYTIC not implemented\n");
        break;
    }

    // Calculate the viscosity
    higflow_compute_viscosity_multiphase(ns);
    // Calculate the density
    higflow_compute_density_multiphase(ns);
    // Calculate the curvature, interfacial force and normal
    higflow_compute_curvature_interfacial_force_normal_multiphase_2D_hf_shirani(ns);
    higflow_compute_distance_multiphase_2D(ns);

    // Calculate the intermediate velocity
    switch (ns->contr.tempdiscrtype) {
    case EXPLICIT_EULER:
        // Explicit Euler method
        higflow_explicit_euler_intermediate_velocity_multiphase_electroosmotic(ns, ns->dpu, ns->dpustar);
        break;
    case EXPLICIT_RK2:
        // Explicit RK2 method
        higflow_explicit_runge_kutta_2_intermediate_velocity_multiphase_electroosmotic(ns);
        break;
    case EXPLICIT_RK3:
        // Explicit RK3 method
        higflow_explicit_runge_kutta_3_intermediate_velocity_multiphase_electroosmotic(ns);
        break;
    case SEMI_IMPLICIT_EULER:
        // Semi-Implicit Euler Method
        higflow_semi_implicit_euler_intermediate_velocity_multiphase_electroosmotic(ns);
        break;
    case SEMI_IMPLICIT_CN:
        // Semi-Implicit Crank-Nicolson Method
        higflow_semi_implicit_crank_nicolson_intermediate_velocity_multiphase_electroosmotic(ns);
        break;
    case SEMI_IMPLICIT_BDF2:
        // Semi-Implicit Crank-Nicolson Method
        higflow_semi_implicit_bdf2_intermediate_velocity_multiphase_electroosmotic(ns);
        break;
    }
    // Set outflow for ustar velocity 
    // higflow_outflow_ustar_step(ns);
    // Calculate the pressure
    higflow_pressure_multiphase(ns);
    // Calculate the final velocity
    higflow_final_velocity_multiphase(ns);
    // Set outflow for u velocity 
    // higflow_outflow_u_step(ns);
    // Calculate the final pressure
    higflow_final_pressure(ns);
    // Calculate the volume fraction
    higflow_plic_advection_volume_fraction(ns);

}
