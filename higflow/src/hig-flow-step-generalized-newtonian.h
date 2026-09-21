// Generalized Newtonian: viscosity is a function of the shear rate, with no elastic
// history -- no tensor to advance, just a viscosity recomputed per cell.
//
// NO EXAMPLE SELECTS THIS PATH.  No case in the tree sets a generalized-Newtonian
// flow type, so nothing here is covered by the regression suite; the nearest
// exercised relative is the variable-viscosity path reached by example2d_BMP.
// Changes here are unverified until a case is written or the code is injected into
// a case that does run.
//
// WHAT A STEP FILE CONTAINS: the intermediate-velocity stage of the projection for
// one rheology, plus whatever constitutive equation that rheology has to advance.
// The shared stages (pressure solve, final velocity) stay in hig-flow-step.h.
//
// Each intermediate-velocity variant is offered in several time discretizations
// (explicit Euler, RK2, RK3, semi-implicit, implicit).  WHICH ONE RUNS COMES FROM
// THE YAML, so a variant present here is not necessarily a variant any case selects.

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
void higflow_semi_implicit_bdf2_intermediate_velocity_gen_newt(higflow_solver *ns); 

// One step of the Navier-Stokes the projection method
void higflow_solver_step_gen_newt(higflow_solver *ns); 

#endif
