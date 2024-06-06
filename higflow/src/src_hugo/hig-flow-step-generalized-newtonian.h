// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Step Generalized Newtonian - version 10/11/2016
// *******************************************************************
// *******************************************************************

#ifndef HIGFLOW_STEP_GENERALIZED_NEWTONIAN
#define HIGFLOW_STEP_GENERALIZED_NEWTONIAN

#include "hig-flow-step.h"

// *******************************************************************
// Navier-Stokes Step
// *******************************************************************

// Calculate the viscosity
void higflow_compute_viscosity_gn(higflow_solver *ns); 

// Computing the necessary term for the Navier-Stokes equation
void higflow_compute_velocity_derivative_tensor(higflow_solver *ns);
    
// One step of the explicit Euler method - Intermediate Velocity
void higflow_explicit_euler_intermediate_velocity_gen_newt(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *pdustar[DIM]); 

// One step of the explicit order 2 Runge-Kutta method - Intermediate Velocity
void higflow_explicit_runge_kutta_2_intermediate_velocity_gen_newt(higflow_solver *ns);

// One step of the explicit order 3 Runge-Kutta method - Intermediate Velocity
void higflow_explicit_runge_kutta_3_intermediate_velocity_gen_newt(higflow_solver *ns);

// One step of the Semi-Implicit Euler method - Intermediate Velocity
void higflow_semi_implicit_euler_intermediate_velocity_gen_newt(higflow_solver *ns); 

// One step of the semi-implicit Crank-Nicolson method - Intermediate Velocity
void higflow_semi_implicit_crank_nicolson_intermediate_velocity_gen_newt(higflow_solver *ns);

// One step of the semi-implicit BDF2 method - Intermediate Velocity
void higflow_semi_implicit_bdf2_intermediate_velocity_gen_newt(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]); 

// One step of the Navier-Stokes the projection method
void higflow_solver_step_gen_newt(higflow_solver *ns); 

#endif
