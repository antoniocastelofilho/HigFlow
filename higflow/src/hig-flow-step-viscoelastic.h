// Viscoelastic flow: the conformation-tensor path, and the reference implementation
// the other viscoelastic variants were derived from.
//
// EXERCISED by example2d_Oldroyd and example2d_Gptt.  Oldroyd is the consistent case
// used as an injection target when a suspicious pattern has to be tested in a solver
// that no example reaches.
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

#ifndef HIG_FLOW_STEP_VISCOELASTIC
#define HIG_FLOW_STEP_VISCOELASTIC
#include "hig-flow-viscoelastic-kernel.h"

#include "hig-flow-linear-algebra.h"

#include "hig-flow-step-generalized-newtonian.h"
#include "hig-flow-mittag-leffler.h"

// *******************************************************************
// Navier-Stokes Step
// *******************************************************************

// One step of the explicit Euler method - Intermediate Velocity
void higflow_explicit_euler_intermediate_velocity_viscoelastic(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *pdustar[DIM]); 

// One step of the explicit order 2 Runge-Kutta method - Intermediate Velocity
void higflow_explicit_runge_kutta_2_intermediate_velocity_viscoelastic(higflow_solver *ns);

// One step of the explicit order 3 Runge-Kutta method - Intermediate Velocity
void higflow_explicit_runge_kutta_3_intermediate_velocity_viscoelastic(higflow_solver *ns);

// One step of the Semi-Implicit Euler method - Intermediate Velocity
void higflow_semi_implicit_euler_intermediate_velocity_viscoelastic(higflow_solver *ns); 

// One step of the semi-implicit Crank-Nicolson method - Intermediate Velocity
void higflow_semi_implicit_crank_nicolson_intermediate_velocity_viscoelastic(higflow_solver *ns);

// One step of the semi-implicit BDF2 method - Intermediate Velocity
void higflow_semi_implicit_bdf2_intermediate_velocity_viscoelastic(higflow_solver *ns); 

// One step of the Navier-Stokes the projection method
void higflow_solver_step_viscoelastic(higflow_solver *ns); 

// Calculate the eige-value and eige-vectors using the Jacobi method

// Calculate the matrix product

// Calculate the matrix product transpose

// Calculate RHS = OK - KO + 2B * M/De

// Calculate the Kernel matrix
void hig_flow_calculate_kernel (higflow_solver *ns, real lambda[DIM], real R[DIM][DIM], real Kernel[DIM][DIM], real tol); 

// Calculate the Omega matrix

// Calculate the matrix BB and the matrix B
void hig_flow_calculate_b (real lambda[DIM], real jlambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real BB[DIM][DIM]);

// Calculate the matrix MM for Oldroyd-B model
void hig_flow_calculate_m_oldroyd (real lambda[DIM], real jlambda[DIM],real M_aux[DIM][DIM], real Re, real trS, ve_parameters *par);

// Calculate the matrix MM for Giesekus model
void hig_flow_calculate_m_giesekus (real lambda[DIM], real jlambda[DIM],real M_aux[DIM][DIM], real Re, real trS, ve_parameters *par);

// Calculate the matrix MM for LPTT model
void hig_flow_calculate_m_lptt (real lambda[DIM], real jlambda[DIM],real M_aux[DIM][DIM], real Re, real trS, ve_parameters *par);

// Calculate the matrix MM for GPTT model
void hig_flow_calculate_m_gptt (real lambda[DIM], real jlambda[DIM],real M_aux[DIM][DIM], real Re, real trS, ve_parameters *par); 

// Calculate the matrix MM for FENE-P model
void hig_flow_calculate_m_fene_p (real lambda[DIM], real jlambda[DIM],real M_aux[DIM][DIM], real Re, real trS, ve_parameters *par);

// Calculate the matrix MM for e-FENE model
void hig_flow_calculate_m_e_fene (real lambda[DIM], real jlambda[DIM],real M_aux[DIM][DIM], real Re, real trS, ve_parameters *par);

// Get the velocity at cell center 

// Get the derivative of Kernel 
void hig_flow_derivative_kernel_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, int i, int j, real Kcenter, real dKdx[DIM]); 

// Calculate the kronecker product 

// Calculate RHS = 2B * M/De

// Solve linear system for constitutive equation

// //Calculate convective tensor term CUBISTA
// real hig_flow_convective_tensor_term_cubista(distributed_property *dpu, sim_facet_domain *sfdu, sim_stencil *stn, distributed_property *dpK, sim_domain *sdED, sim_stencil *stnED, real kc, Point ccenter, Point cdelta, int dim);

// // Computing the Kernel Tensor
// void higflow_compute_kernel_tensor(higflow_solver *ns);

// // Computing the Initial Kernel Tensor
// void higflow_compute_initial_kernel_tensor(higflow_solver *ns);

// Constitutive Equation Step for the Explicit Euler Method
void higflow_explicit_euler_constitutive_equation(higflow_solver *ns);

// Constitutive Equation Step for the Implicit Euler Method
void higflow_implicit_euler_constitutive_equation(higflow_solver *ns);

// Computing the Polymeric Tensor
void higflow_compute_polymeric_tensor(higflow_solver *ns);

// Compute initial conformation tensor for e-FENE model
void hig_flow_compute_initial_conformation_e_fene(real Rhs[DIM][DIM], real A[DIM][DIM], real b, real l, real E);

// Compute initial velocity derivative tensor for e-FENE model
void hig_flow_compute_initial_velocity_derivative_tensor(higflow_solver *ns);

#endif
