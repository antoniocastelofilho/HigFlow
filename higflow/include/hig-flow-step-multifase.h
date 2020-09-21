// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Step Multifase - version 10/11/2016
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_STEP_MULTIFASE
#define HIG_FLOW_STEP_MULTIFASE

#include "hig-flow-step.h"

// *******************************************************************
// Navier-Stokes Step
// *******************************************************************

// Calculate the viscosity
void higflow_compute_viscosity_multifase(higflow_solver *ns); 

// Calculate the density
void higflow_compute_density_multifase(higflow_solver *ns); 

// Navier-Stokes final velocity using the projection method for multifase flow
void higflow_final_velocity_multifase(higflow_solver *ns); 

// Navier-Stokes pressure with variable density using the projection method
void higflow_pressure_multifase(higflow_solver *ns); 

// One step of the explicit Euler method - Intermediate Velocity
void higflow_explicit_euler_intermediate_velocity_multifase(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *pdustar[DIM]); 

// One step of the explicit order 2 Runge-Kutta method - Intermediate Velocity
void higflow_explicit_runge_kutta_2_intermediate_velocity_multifase(higflow_solver *ns);

// One step of the explicit order 3 Runge-Kutta method - Intermediate Velocity
void higflow_explicit_runge_kutta_3_intermediate_velocity_multifase(higflow_solver *ns);

// One step of the Semi-Implicit Euler method - Intermediate Velocity
void higflow_semi_implicit_euler_intermediate_velocity_multifase(higflow_solver *ns); 

// One step of the semi-implicit Crank-Nicolson method - Intermediate Velocity
void higflow_semi_implicit_crank_nicolson_intermediate_velocity_multifase(higflow_solver *ns);

// One step of the semi-implicit BDF2 method - Intermediate Velocity
void higflow_semi_implicit_bdf2_intermediate_velocity_multifase(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]); 

// One step of the Navier-Stokes the projection method
void higflow_solver_step_multifase(higflow_solver *ns); 

#endif
