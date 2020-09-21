// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Step - version 10/11/2016
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_STEP
#define HIG_FLOW_STEP

#include "hig-flow-kernel.h"
#include "hig-flow-eval.h"
#include "hig-flow-discret.h"
#include "hig-flow-terms.h"

// *******************************************************************
// Navier-Stokes Step
// *******************************************************************

// Maybe remove_pressure_singularity - it will depend on the boundary situation
void remove_pressure_singularity(higflow_solver *ns, solver *slvp); 

// Set the velocity at outflow
void set_outflow(psim_facet_domain * psfdu, distributed_property * dpu, real alpha); 

// Navier-Stokes pressure using the projection method
void higflow_pressure(higflow_solver *ns); 

// Navier-Stokes final velocity using the projection method 
void higflow_final_velocity(higflow_solver *ns); 

// Navier-Stokes final pressure using the projection method
void higflow_final_pressure(higflow_solver *ns); 

// Computing the necessary term for the Navier-Stokes equation
void higflow_compute_rate_of_deformation(higflow_solver *ns); 
    
// Navier-Stokes calculate the source term
void higflow_calculate_source_term(higflow_solver *ns); 

// Navier-Stokes calculate the facet source term
void higflow_calculate_facet_source_term(higflow_solver *ns); 

// Apply the boundary condition for the velocity
void higflow_boundary_condition_for_velocity(higflow_solver *ns); 

// Apply the boundary condition for the pressure
void higflow_boundary_condition_for_pressure(higflow_solver *ns);

// Navier-Stokes outflow for u velocity 
void higflow_outflow_u_step(higflow_solver *ns);

// Navier-Stokes outflow for ustar velocity 
void higflow_outflow_ustar_step(higflow_solver *ns);

// One step of the explicit Euler method - Intermediate Velocity
void higflow_explicit_euler_intermediate_velocity(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *pdustar[DIM]); 

// One step of the explicit order 2 Runge-Kutta method - Intermediate Velocity
void higflow_explicit_runge_kutta_2_intermediate_velocity(higflow_solver *ns);

// One step of the explicit order 3 Runge-Kutta method - Intermediate Velocity
void higflow_explicit_runge_kutta_3_intermediate_velocity(higflow_solver *ns);

// One step of the Semi-Implicit Euler method - Intermediate Velocity
void higflow_semi_implicit_euler_intermediate_velocity(higflow_solver *ns); 

// One step of the semi-implicit Crank-Nicolson method - Intermediate Velocity
void higflow_semi_implicit_crank_nicolson_intermediate_velocity(higflow_solver *ns);

// One step of the semi-implicit BDF2 method - Intermediate Velocity
void higflow_semi_implicit_bdf2_intermediate_velocity(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]); 

// One step of the Navier-Stokes the projection method
void higflow_solver_step(higflow_solver *ns); 

#endif
