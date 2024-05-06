// *******************************************************************
//  HiG-Flow Solver Step Electro-osmotic - version 25/01/2022
// *******************************************************************

#include "hig-flow-step-electroosmotic-viscoelastic.h"

// One step of the Navier-Stokes the projection method for viscoelastic flow
void higflow_solver_step_electroosmotic_viscoelastic(higflow_solver *ns) {
    // Calculate the velocity derivative tensor
    higflow_compute_velocity_derivative_tensor(ns);

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
   switch (ns->ed.eo.contr.eo_model) {
    case PNP:
        // Poisson-Nernst-Planck model
        if( (ns->par.step==0) || ns->ed.eo.contr.is_phibc_timedependent == true) {
            higflow_boundary_condition_for_phi(ns);
            higflow_electroosmotic_phi(ns);
        }
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
        higflow_boundary_condition_for_psi(ns);
        higflow_electroosmotic_psi(ns);
        higflow_calculate_electroosmotic_source_term_pnp(ns);
        break;
    case PB:
        // Poisson-Boltzmann model 
        if( (ns->par.step==0) || ns->ed.eo.contr.is_phibc_timedependent == true) {
            higflow_boundary_condition_for_phi(ns);
            higflow_electroosmotic_phi(ns);
        }
        if( (ns->par.step==0) || ns->ed.eo.contr.is_psibc_timedependent == true) {
            higflow_boundary_condition_for_psi(ns);
            higflow_electroosmotic_psi(ns);
        }
        if( (ns->par.step==0) || ns->ed.eo.contr.is_phibc_timedependent == true
                              || ns->ed.eo.contr.is_psibc_timedependent == true) {
            higflow_calculate_electroosmotic_source_term_pb(ns);
        }
        break;
    case PBDH:
        // Debye-Hückel model (solving poisson equation for psi) 
        // Poisson-Boltzmann model 
        if( (ns->par.step==0) || ns->ed.eo.contr.is_phibc_timedependent == true) {
            higflow_boundary_condition_for_phi(ns);
            higflow_electroosmotic_phi(ns);
        }
        if( (ns->par.step==0) || ns->ed.eo.contr.is_psibc_timedependent == true) {
            higflow_boundary_condition_for_psi(ns);
            higflow_electroosmotic_psi(ns);
        }
        if( (ns->par.step==0) || ns->ed.eo.contr.is_phibc_timedependent == true
                              || ns->ed.eo.contr.is_psibc_timedependent == true) {
            higflow_calculate_electroosmotic_source_term_pbdh(ns);
        }
        break;
    case PBDH_ANALYTIC:
        // Debye-Hückel model (using the analytic solution for psi)
        higflow_calculate_electroosmotic_source_term_analytic_pbdh(ns);
        break;
    }
    // Calculate the intermediated velocity
   // Calculate the intermediate velocity
   switch (ns->contr.tempdiscrtype) {
   case EXPLICIT_EULER:
      // Explicit Euler method
      higflow_explicit_euler_intermediate_velocity_electroosmotic(ns, ns->dpu, ns->dpustar);
      break;
   case EXPLICIT_RK2:
      // Explicit RK2 method
      higflow_explicit_runge_kutta_2_intermediate_velocity_electroosmotic(ns);
      break;
   case EXPLICIT_RK3:
      // Explicit RK3 method
      higflow_explicit_runge_kutta_3_intermediate_velocity_electroosmotic(ns);
      break;
   case SEMI_IMPLICIT_EULER:
      // Semi-Implicit Euler Method
      higflow_semi_implicit_euler_intermediate_velocity_electroosmotic(ns);
      break;
   case SEMI_IMPLICIT_CN:
      // Semi-Implicit Crank-Nicolson Method
      higflow_semi_implicit_crank_nicolson_intermediate_velocity_electroosmotic(ns);
      break;
   case SEMI_IMPLICIT_BDF2:
      // Semi-Implicit Crank-Nicolson Method
      higflow_semi_implicit_bdf2_intermediate_velocity_electroosmotic(ns, ns->dpu, ns->dpustar);
      break;
   }
   // Set outflow for ustar velocity 
   // higflow_outflow_ustar_step(ns);
   // Calculate the pressure
   higflow_pressure(ns);
   // Calculate the final velocity
   higflow_final_velocity(ns);
   // Set outflow for u velocity 
   // higflow_outflow_u_step(ns);
   // Calculate the final pressure
   higflow_final_pressure(ns);
}

