// Elastoviscoplastic flow: elastic below a yield stress, viscoplastic above it.  The
// yield criterion is evaluated per cell, so the same run has both regimes at once.
//
// NO EXAMPLE SELECTS THIS PATH.  Uncovered by the suite; it shares the conformation-
// tensor machinery in hig-flow-viscoelastic-kernel.h with the paths that are
// covered, so a change THERE is witnessed, a change HERE is not.
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
//  HiG-Flow Solver Step Viscoelastic - version 10/11/2016
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_STEP_ELASTOVISCOPLASTIC
#define HIG_FLOW_STEP_ELASTOVISCOPLASTIC
#include "hig-flow-viscoelastic-kernel.h"

#include "hig-flow-linear-algebra.h"

#include "hig-flow-step-generalized-newtonian.h"

// *******************************************************************
// Navier-Stokes Step
// *******************************************************************

// One step of the explicit Euler method - Intermediate Velocity
void higflow_explicit_euler_intermediate_velocity_elastoviscoplastic(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *pdustar[DIM]); 

// One step of the explicit order 2 Runge-Kutta method - Intermediate Velocity
void higflow_explicit_runge_kutta_2_intermediate_velocity_elastoviscoplastic(higflow_solver *ns);

// One step of the explicit order 3 Runge-Kutta method - Intermediate Velocity
void higflow_explicit_runge_kutta_3_intermediate_velocity_elastoviscoplastic(higflow_solver *ns);

// One step of the Semi-Implicit Euler method - Intermediate Velocity
void higflow_semi_implicit_euler_intermediate_velocity_elastoviscoplastic(higflow_solver *ns); 

// One step of the semi-implicit Crank-Nicolson method - Intermediate Velocity
void higflow_semi_implicit_crank_nicolson_intermediate_velocity_elastoviscoplastic(higflow_solver *ns);

// One step of the semi-implicit BDF2 method - Intermediate Velocity
void higflow_semi_implicit_bdf2_intermediate_velocity_elastoviscoplastic(higflow_solver *ns); 

// One step of the Navier-Stokes the projection method
void higflow_solver_step_elastoviscoplastic(higflow_solver *ns); 

// Calculate the eige-value and eige-vectors using the Jacobi method

// Calculate the matrix product

// Calculate the matrix product transpose

// Calculate RHS = OK - KO + 2B * M/De

// Calculate the Kernel matrix
void hig_flow_calculate_kernel (higflow_solver *ns, real lambda[DIM], real R[DIM][DIM], real Kernel[DIM][DIM], real tol); 

// Calculate the Omega matrix

// Calculate the matrix BB
void hig_flow_calculate_b (higflow_solver *ns, real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real BB[DIM][DIM], real tol); 

// Calculate the matrix MM for Oldroyd-B-Bingham model
void hig_flow_calculate_m_oldroyd_bingham (higflow_solver *ns, real lambda[DIM], real M[DIM][DIM], real R[DIM][DIM], real M_aux[DIM][DIM], real tol, real smallTD, real SD); 

// Calculate the matrix MM for Oldroyd-B-Herschel-Bulkley
void hig_flow_calculate_m_oldroyd_HB (higflow_solver *ns, real lambda[DIM], real M[DIM][DIM], real R[DIM][DIM], real M_aux[DIM][DIM], real tol, real smallTD, real SD); 

// Calculate the matrix MM for LPTT-Bingham model
void hig_flow_calculate_m_lptt_bingham (higflow_solver *ns, real tr, real lambda[DIM], real M[DIM][DIM], real R[DIM][DIM], real M_aux[DIM][DIM], real tol, real smallTD, real SD); 

// Calculate the matrix MM for EPTT-Bingham model
void hig_flow_calculate_m_eptt_bingham(higflow_solver *ns, real tr, real lambda[DIM],  real M[DIM][DIM], real R[DIM][DIM], real M_aux[DIM][DIM], real tol, real smallTD, real SD); 

// Get the velocity at cell center 

// Get the derivative of Kernel 
void hig_flow_derivative_kernel_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, int i, int j, real Kcenter, real dKdx[DIM]); 

// Calculate the kronecker product 

// Calculate RHS = 2B * M/De

// Solve linear system for constitutive equation

//Calculate convective tensor term CUBISTA
real hig_flow_convective_tensor_term_cubista(higflow_solver *ns, distributed_property *dpu, sim_domain *sdp, sim_stencil *stn, real K[DIM][DIM], Point ccenter, Point cdelta, int dim, int i, int j);

// Computing the Polymeric Tensor
void higflow_compute_polymeric_tensor_elastoviscoplastic(higflow_solver *ns);

// Constitutive Equation Step for the Explicit Euler Method
void higflow_explicit_euler_constitutive_equation_elastoviscoplastic(higflow_solver *ns);

// Constitutive Equation Step for the Implicit Euler Method
void higflow_implicit_euler_constitutive_equation_elastoviscoplastic(higflow_solver *ns);

// Computing the Kernel Tensor
void higflow_compute_kernel_tensor_elastoviscoplastic(higflow_solver *ns);

#endif
