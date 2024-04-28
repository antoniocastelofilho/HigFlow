// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Step Electro-osmotic - version 02/02/2018
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_STEP_ELECTROOSMOTIC
#define HIG_FLOW_STEP_ELECTROOSMOTIC

#include "hig-flow-step-viscoelastic.h"

typedef enum charge_type {
    POSITIVE = 1,
    NEGATIVE = -1,
} charge_type;

// *******************************************************************
// Electro-osmotic source term
// *******************************************************************

//Semi-implicit ionic transport equation for nplus
void higflow_semi_implicit_euler_ionic_transport_equation_nplus(higflow_solver *ns); 

//Semi-implicit ionic transport equation for nminus
void higflow_semi_implicit_euler_ionic_transport_equation_nminus(higflow_solver *ns);

// Time step of ionic transport equation for n+
void higflow_explicit_euler_ionic_transport_equation_nplus(higflow_solver *ns); 

// Time step of ionic transport equation for n-
void higflow_explicit_euler_ionic_transport_equation_nminus(higflow_solver *ns); 

// Compute the electric convective ionic term
real hig_flow_convective_ionic_cell_term_cubista(higflow_solver *ns, real nc, Point ccenter, Point cdelta, int dim, charge_type charge_sign);

// Souce term by analytical PBDH model
void higflow_calculate_electroosmotic_source_term_analytic_pbdh( higflow_solver *ns);

// Souce term by PBDH model
void higflow_calculate_electroosmotic_source_term_pbdh( higflow_solver *ns);

// Souce term by PB model
void higflow_calculate_electroosmotic_source_term_pb( higflow_solver *ns);

// Souce term by PNP model
void higflow_calculate_electroosmotic_source_term_pnp( higflow_solver *ns);

// Boundary conditiono for electro-osmotic source term
void higflow_boundary_condition_for_electroosmotic_source_term(higflow_solver *ns); 

// Electroosmotic induced potential psi 
void higflow_electroosmotic_psi(higflow_solver *ns);

// Electroosmotic applied potential psi 
void higflow_electroosmotic_phi(higflow_solver *ns);

// Navier-Stokes Step for the Explicit Euler Method
void higflow_explicit_euler_intermediate_velocity_electroosmotic(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]); 

// Navier-Stokes Step for the Second Order Explicit Runge-Kutta Method
void higflow_explicit_runge_kutta_2_intermediate_velocity_electroosmotic(higflow_solver *ns); 

// Navier-Stokes Step for third Order Explicit Runge-Kutta Method
void higflow_explicit_runge_kutta_3_intermediate_velocity_electroosmotic(higflow_solver *ns); 

// Navier-Stokes Step for the Implicit Euler Method
void higflow_semi_implicit_euler_intermediate_velocity_electroosmotic(higflow_solver *ns); 

// Navier-Stokes Step for the Semi-Implicit Crank-Nicolson Method
void higflow_semi_implicit_crank_nicolson_intermediate_velocity_electroosmotic(higflow_solver *ns); 

// Navier-Stokes Step for the Implicit BDF2 Method
void higflow_semi_implicit_bdf2_intermediate_velocity_electroosmotic(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]); 

// print minimum and maximum values of distributed properties
void print_minmax_properties(higflow_solver *ns);

// One step of the Navier-Stokes the projection method
void higflow_solver_step_electroosmotic(higflow_solver *ns); 

// Apply the boundary condition for source term 
void higflow_boundary_condition_for_electroosmotic_nplus(higflow_solver *ns);
void higflow_boundary_condition_for_electroosmotic_nminus(higflow_solver *ns);
void higflow_boundary_condition_for_phi(higflow_solver *ns);
void higflow_boundary_condition_for_psi(higflow_solver *ns);

#endif
