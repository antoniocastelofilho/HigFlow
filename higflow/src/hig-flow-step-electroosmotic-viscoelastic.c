// *******************************************************************
//  HiG-Flow Solver Step Electro-osmotic - version 03/2024
// *******************************************************************

#include "hig-flow-step-electroosmotic-viscoelastic.h"

// One step of the Navier-Stokes the projection method for viscoelastic flow
void higflow_solver_step_electroosmotic_viscoelastic(higflow_solver *ns) {

    // Compute deformation rate tensor
    higflow_compute_velocity_derivative_tensor(ns);
    // Constitutive Equation Step for viscoelastic flow
    switch (ns->ed.ve.contr.discrtype) {
    case EXPLICIT:
        higflow_explicit_euler_constitutive_equation(ns);
        break;
    case IMPLICIT:
        higflow_implicit_euler_constitutive_equation(ns);
        break;
    }
    // Calculate the elastic tensor to be used in the momentum equation
    higflow_compute_polymeric_tensor(ns);

    // Boundary conditions and source terms
    higflow_boundary_condition_for_velocity(ns);
    higflow_boundary_condition_for_pressure(ns);
    higflow_calculate_source_term(ns);
    higflow_calculate_facet_source_term(ns);

    // Calculate the electro-osmotic terms
    higflow_boundary_condition_for_phi(ns);
    higflow_boundary_condition_for_psi(ns);
    switch (ns->ed.eo.contr.eo_model) {
    case PNP: // Poisson-Nernst-Planck model
        higflow_boundary_condition_for_electroosmotic_nplus(ns);
        higflow_boundary_condition_for_electroosmotic_nminus(ns);
        if( (ns->par.step == 0) || ns->ed.eo.contr.is_phibc_timedependent == true)
            higflow_electroosmotic_phi(ns);
        switch (ns->ed.eo.contr.tempdiscrtype) {
            case EXPLICIT_EULER:
                higflow_explicit_euler_ionic_transport_equation_nplus(ns);
                higflow_explicit_euler_ionic_transport_equation_nminus(ns);
                break;
            case SEMI_IMPLICIT_EULER:
                higflow_semi_implicit_euler_ionic_transport_equation_nplus(ns);
                higflow_semi_implicit_euler_ionic_transport_equation_nminus(ns);
                break;
        }
        higflow_electroosmotic_psi(ns);
        higflow_calculate_electroosmotic_source_term(ns);
        break;
    case PB: // Poisson-Boltzmann model 
        if( (ns->par.step == 0) || ns->ed.eo.contr.is_phibc_timedependent == true)
            higflow_electroosmotic_phi(ns);
        if( (ns->par.step == 0) || ns->ed.eo.contr.is_psibc_timedependent == true)
            higflow_electroosmotic_solve_pb(ns);
        if( (ns->par.step == 0) || ns->ed.eo.contr.is_phibc_timedependent == true
                                || ns->ed.eo.contr.is_psibc_timedependent == true)
            higflow_calculate_electroosmotic_source_term(ns);
        break;
    case PBDH: // Debye-Hückel model (solving poisson equation for psi) 
        if( (ns->par.step == 0) || ns->ed.eo.contr.is_phibc_timedependent == true)
            higflow_electroosmotic_phi(ns);
        if( (ns->par.step == 0) || ns->ed.eo.contr.is_psibc_timedependent == true)
            higflow_electroosmotic_psi(ns);
        if( (ns->par.step == 0) || ns->ed.eo.contr.is_phibc_timedependent == true
                                || ns->ed.eo.contr.is_psibc_timedependent == true)
            higflow_calculate_electroosmotic_source_term(ns);
        break;
    case PBDH_ANALYTIC: // Debye-Hückel model (using the analytic solution for psi)
        higflow_calculate_electroosmotic_source_term_analytic_pbdh(ns);
        break;
    } 

    // Calculate the intermediate velocity
    switch (ns->contr.tempdiscrtype) {
    case EXPLICIT_EULER:
        higflow_explicit_euler_intermediate_velocity_electroosmotic(ns, ns->dpu, ns->dpustar);
        break;
    case EXPLICIT_RK2:
        higflow_explicit_runge_kutta_2_intermediate_velocity_electroosmotic(ns);
        break;
    case EXPLICIT_RK3:
        higflow_explicit_runge_kutta_3_intermediate_velocity_electroosmotic(ns);
        break;
    case SEMI_IMPLICIT_EULER:
        higflow_semi_implicit_euler_intermediate_velocity_electroosmotic(ns);
        break;
    case SEMI_IMPLICIT_CN:
        higflow_semi_implicit_crank_nicolson_intermediate_velocity_electroosmotic(ns);
        break;
    case SEMI_IMPLICIT_BDF2:
        higflow_semi_implicit_bdf2_intermediate_velocity_electroosmotic(ns);
        break;
    }

   // Projection
   higflow_pressure(ns);
   higflow_final_velocity(ns);
   higflow_final_pressure(ns);
}

