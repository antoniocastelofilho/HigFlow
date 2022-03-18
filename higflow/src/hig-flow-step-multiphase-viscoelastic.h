// *******************************************************************
//  HiG-Flow Solver Step multiphase - version 10/11/2016
// *******************************************************************

#ifndef HIG_FLOW_STEP_MULTIPHASE_VISC
#define HIG_FLOW_STEP_MULTIPHASE_VISC
#include "hig-flow-step.h"
#include "hig-flow-vof-plic.h"
#include "hig-flow-vof-elvira.h"
#include "hig-flow-vof-adap-hf.h"
#include "hig-flow-mittag-leffler.h"
#include "hig-flow-step-multiphase.h"
#include "hig-flow-step-generalized-newtonian.h"

// *******************************************************************
// Navier-Stokes Step
// *******************************************************************
// Computing beta viscoelastic
void higflow_compute_beta_multiphase_viscoelastic(higflow_solver *ns);

// Computing tensor viscoelastic
void higflow_compute_S_multiphase_viscoelastic(higflow_solver *ns) ;

// Computing the Kernel Tensor
void higflow_compute_kernel_tensor_multiphase_viscoelastic(higflow_solver *ns);

// Computing the Polymeric Tensor
void higflow_compute_polymeric_tensor_multiphase_viscoelastic(higflow_solver *ns) ;

// Constitutive Equation Step for the Explicit Euler Method
void higflow_explicit_euler_constitutive_equation_multiphase_viscoelastic(higflow_solver *ns);

// Calculate convective tensor term CUBISTA - phase 0 
real hig_flow_conv_tensor_term_cub_mult_visc_phase0(higflow_solver *ns, distributed_property *dpu, sim_domain *sdp, sim_stencil *stn, real K[DIM][DIM], Point ccenter, Point cdelta, int dim, int i, int j);

// Calculate convective tensor term CUBISTA - phase 1
real hig_flow_conv_tensor_term_cub_mult_visc_phase1(higflow_solver *ns, distributed_property *dpu, sim_domain *sdp, sim_stencil *stn, real K[DIM][DIM], Point ccenter, Point cdelta, int dim, int i, int j);

// Constitutive Equation Step for the Implicit Euler Method
void higflow_implicit_euler_constitutive_equation_multiphase_viscoelastic(higflow_solver *ns);

// Calculate RHS = OK - KO + 2B * M/De
void hig_flow_kernel_rhs_multiphase_viscoelastic(real De, real K[DIM][DIM], real O[DIM][DIM], real B[DIM][DIM], real M[DIM][DIM], real RHS[DIM][DIM]) ;

// Calculate RHS = OK - KO + 2B * M/De
void hig_flow_kernel_rhs_multiphase_viscoelastic(real De, real K[DIM][DIM], real O[DIM][DIM], real B[DIM][DIM], real M[DIM][DIM], real RHS[DIM][DIM]) ;

// Calculate the eige-value and eige-vectors using the Jacobi method
void hig_flow_jacobi(real A[DIM][DIM], real d[DIM], real V[DIM][DIM]) ;

// Calculate the matrix product
void hig_flow_matrix_product (real A[DIM][DIM], real R[DIM][DIM], real B[DIM][DIM]) ;

// Calculate the matrix product transpose
void hig_flow_matrix_transpose_product (real A[DIM][DIM], real R[DIM][DIM], real B[DIM][DIM]) ;

// Calculate RHS = OK - KO + 2B * M/De
void hig_flow_kernel_rhs (real De, real K[DIM][DIM], real O[DIM][DIM], real B[DIM][DIM], real M[DIM][DIM], real RHS[DIM][DIM]) ;

// Calculate the Kernel matrix - phase 0
void hig_flow_calculate_kernel_phase0(higflow_solver *ns, real lambda[DIM], real R[DIM][DIM], real Kernel[DIM][DIM], real tol) ;

// Calculate the Kernel matrix - phase 1
void hig_flow_calculate_kernel_phase1(higflow_solver *ns, real lambda[DIM], real R[DIM][DIM], real Kernel[DIM][DIM], real tol) ;

// Calculate the Omega matrix
void hig_flow_calculate_omega (real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real Omega[DIM][DIM], real small) ;

// Calculate the matrix BB
void hig_flow_calculate_b (higflow_solver *ns, real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real BB[DIM][DIM], real tol) ;

// Calculate the matrix MM for Oldroyd model
void hig_flow_calculate_m_oldroyd (higflow_solver *ns, real lambda[DIM], real M[DIM][DIM], real M_aux[DIM][DIM], real tol) ;

// Calculate the matrix MM for Giesekus model - phase 0
void hig_flow_calculate_m_giesekus_phase0 (higflow_solver *ns, real lambda[DIM], real M[DIM][DIM], real M_aux[DIM][DIM], real tol) ;

// Calculate the matrix MM for Giesekus model - phase 1
void hig_flow_calculate_m_giesekus_phase1 (higflow_solver *ns, real lambda[DIM], real M[DIM][DIM], real M_aux[DIM][DIM], real tol) ;

// Calculate the matrix MM for LPTT model - phase 0
void hig_flow_calculate_m_lptt_phase0 (higflow_solver *ns, real tr, real lambda[DIM],  real M[DIM][DIM], real R[DIM][DIM], real M_aux[DIM][DIM], real tol) ;

// Calculate the matrix MM for LPTT model - phase 1
void hig_flow_calculate_m_lptt_phase1 (higflow_solver *ns, real tr, real lambda[DIM],  real M[DIM][DIM], real R[DIM][DIM], real M_aux[DIM][DIM], real tol) ;

// Calculate the matrix MM for GPTT model - phase 0
void hig_flow_calculate_m_gptt_phase0 (higflow_solver *ns, real tr, real lambda[DIM],  real M[DIM][DIM], real R[DIM][DIM], real M_aux[DIM][DIM], real tol) ;

// Calculate the matrix MM for GPTT model - phase 1
void hig_flow_calculate_m_gptt_phase1 (higflow_solver *ns, real tr, real lambda[DIM],  real M[DIM][DIM], real R[DIM][DIM], real M_aux[DIM][DIM], real tol) ;

// Get the velocity at cell center 
void hig_flow_velocity_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, real u[DIM]) ;

// Get the derivative of Kernel - phase 0
void hig_flow_derivative_kernel_at_center_cell_phase0 (higflow_solver *ns, Point ccenter, Point cdelta, int i, int j, real Kcenter, real dKdx[DIM]) ;

// Get the derivative of Kernel - phase 1
void hig_flow_derivative_kernel_at_center_cell_phase1 (higflow_solver *ns, Point ccenter, Point cdelta, int i, int j, real Kcenter, real dKdx[DIM]) ;

//Gauss elimination to solve the constitutive equation
void hig_flow_solve_system_constitutive_equation ( int n, real A[DIM*DIM][DIM*DIM+1], real x[DIM*DIM] ) ;

// Calculate RHS = 2B * M/De
void hig_flow_implicit_kernel_rhs (real De, real B[DIM][DIM], real M[DIM][DIM], real RHS[DIM][DIM]) ;

// Calculate the matrix product
void hig_flow_kernel_system_matrix (real w[DIM*DIM][DIM*DIM+1], real Omega[DIM][DIM], real dt) ;

// Set velocity for the test: Zalesak and single vortex
void higflow_set_velocity_multiphase(higflow_solver *ns) ;

// One step of the Navier-Stokes the projection method
void higflow_solver_step_multiphase_viscoelastic(higflow_solver *ns); 

#endif
