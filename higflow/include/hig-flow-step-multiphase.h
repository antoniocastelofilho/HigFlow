// *******************************************************************
//  HiG-Flow Solver Step multiphase - version 20/01/2022
// *******************************************************************

#ifndef HIG_FLOW_STEP_MULTIPHASE
#define HIG_FLOW_STEP_MULTIPHASE

#include "hig-flow-step.h"
#include "hig-flow-vof-plic.h"
#include "hig-flow-vof-9-cells.h"
#include "hig-flow-vof-elvira.h"
#include "hig-flow-vof-adap-hf.h"

// *******************************************************************
// Navier-Stokes Step
// *******************************************************************
// Computing the filter for the curvature
void higflow_compute_filter_fracvol_multiphase(higflow_solver *ns);

// Computing the curvature
void higflow_compute_curvature_multiphase(higflow_solver *ns);

// Computing the curvature for DIM = 3
void higflow_compute_curvature_multiphase_3D(higflow_solver *ns);

// Calculate the viscosity
void higflow_compute_viscosity_multiphase(higflow_solver *ns); 

// Calculate the density
void higflow_compute_density_multiphase(higflow_solver *ns);

// ******************************************************************
// Plic advection
// ******************************************************************
void higflow_plic_advection_volume_fraction(higflow_solver *ns);

void higflow_plic_advection_volume_fraction_x_direction(higflow_solver *ns, int dim);

void higflow_plic_advection_volume_fraction_y_direction(higflow_solver *ns, int dim);

void higflow_plic_advection_volume_fraction_x_direction_imp(higflow_solver *ns, int dim);

void higflow_plic_advection_volume_fraction_y_direction_imp(higflow_solver *ns, int dim);

void higflow_plic_copy_fractionaux_to_fraction(higflow_solver *ns);

void normal_correction_at_get(Point Normal);

// *******************************************************************
// Volume Fraction Transport Step for the Explicit Euler Method
// *******************************************************************
void higflow_explicit_euler_volume_fraction(higflow_solver *ns); 

// Get the derivative of Kernel 
void hig_flow_derivative_fracvol_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, real FVcenter, real dfracvoldx[DIM]); 

// Get the velocity at cell center 
void hig_flow_velocity_at_center_cell_multiphase (higflow_solver *ns, Point ccenter, Point cdelta, real u[DIM]); 

// *******************************************************************
// Calculate convective term CUBISTA for volume fraction
// *******************************************************************
real hig_flow_fracvol_term_cubista(higflow_solver *ns, distributed_property *dpu, sim_domain *sdp, sim_stencil *stn, real fracvol, Point ccenter, Point cdelta, int dim); 

// Navier-Stokes final velocity using the projection method for multiphase flow
void higflow_final_velocity_multiphase(higflow_solver *ns); 

// Navier-Stokes pressure with variable density using the projection method
void higflow_pressure_multiphase(higflow_solver *ns); 

// One step of the explicit Euler method - Intermediate Velocity
void higflow_explicit_euler_intermediate_velocity_multiphase(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *pdustar[DIM]); 

// One step of the explicit order 2 Runge-Kutta method - Intermediate Velocity
void higflow_explicit_runge_kutta_2_intermediate_velocity_multiphase(higflow_solver *ns);

// One step of the explicit order 3 Runge-Kutta method - Intermediate Velocity
void higflow_explicit_runge_kutta_3_intermediate_velocity_multiphase(higflow_solver *ns);

// One step of the Semi-Implicit Euler method - Intermediate Velocity
void higflow_semi_implicit_euler_intermediate_velocity_multiphase(higflow_solver *ns); 

// One step of the semi-implicit Crank-Nicolson method - Intermediate Velocity
void higflow_semi_implicit_crank_nicolson_intermediate_velocity_multiphase(higflow_solver *ns);

// One step of the semi-implicit BDF2 method - Intermediate Velocity
void higflow_semi_implicit_bdf2_intermediate_velocity_multiphase(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]); 

// One step of the Navier-Stokes the projection method
void higflow_solver_step_multiphase(higflow_solver *ns); 

#endif
