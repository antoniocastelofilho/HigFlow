// Two-phase flow: a volume fraction tracks the interface, and density, viscosity and
// curvature are derived from it per cell.  The interface reconstruction itself is in
// the hig-flow-vof-*.h family.
//
// EXERCISED by example2d_VOF, example2d_VOF_Oldroyd and example2d_VOF_Gptt.
//
// UNDER A MULTIPHASE RUN THE SINGLE-PHASE flowtype IS INERT -- flowtype0 and
// flowtype1, the per-phase rheologies, are what select behaviour.  The converse also
// holds and is the commoner confusion: every single-phase case in this tree still
// carries flowtype0/flowtype1 values in its YAML, and they do nothing there.  Do not
// read a case's rheology off the wrong key.
//
// WHAT A STEP FILE CONTAINS: the intermediate-velocity stage of the projection for
// one rheology, plus whatever constitutive equation that rheology has to advance.
// The shared stages (pressure solve, final velocity) stay in hig-flow-step.h.
//
// Each intermediate-velocity variant is offered in several time discretizations
// (explicit Euler, RK2, RK3, semi-implicit, implicit).  WHICH ONE RUNS COMES FROM
// THE YAML, so a variant present here is not necessarily a variant any case selects.

// *******************************************************************
//  HiG-Flow Solver Step multiphase - version 20/01/2022
// *******************************************************************

#ifndef HIG_FLOW_STEP_MULTIPHASE
#define HIG_FLOW_STEP_MULTIPHASE

#include "hig-flow-step.h"
#include "hig-flow-vof-plic.h"
#include "hig-flow-vof-plic-3D.h"
#include "hig-flow-vof-9-cells.h"
#include "hig-flow-vof-elvira.h"
#include "hig-flow-vof-adap-hf.h"
#include "hig-flow-vof-HF-3D.h"
#include "hig-flow-vof-advection-3D.h"
#include "hig-flow-vof-mehta.h"
#include "hig-flow-vof-finite-difference-normal-curvature_3D.h"

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


// *******************************************************************
// Volume Fraction Transport Step for the Explicit Euler Method
// *******************************************************************
void higflow_explicit_euler_volume_fraction(higflow_solver *ns); 

// Get the derivative of Kernel 
void hig_flow_derivative_fracvol_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, real FVcenter, real dfracvoldx[DIM]); 

// Get the velocity at cell center 

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
void higflow_semi_implicit_bdf2_intermediate_velocity_multiphase(higflow_solver *ns); 

// One step of the Navier-Stokes the projection method
void higflow_solver_step_multiphase(higflow_solver *ns); 

#endif
