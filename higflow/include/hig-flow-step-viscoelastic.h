// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Step Viscoelastic - version 10/11/2016
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_STEP_VISCOELASTIC
#define HIG_FLOW_STEP_VISCOELASTIC

#include "hig-flow-step-generalized-newtonian.h"

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
void higflow_semi_implicit_bdf2_intermediate_velocity_viscoelastic(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]); 

// One step of the Navier-Stokes the projection method
void higflow_solver_step_viscoelastic(higflow_solver *ns); 

// Calculate the eige-value and eige-vectors using the Jacobi method
void hig_flow_jacobi(real A[DIM][DIM], real d[DIM], real V[DIM][DIM]); 

// Calculate the matrix product
void hig_flow_matrix_product( real A[DIM][DIM], real R[DIM][DIM], real B[DIM][DIM]); 

// Calculate the matrix product transpose
void hig_flow_matrix_transpose_product( real A[DIM][DIM], real R[DIM][DIM], real B[DIM][DIM]); 

// Calculate RHS = OK - KO + 2B * M/De
void hig_flow_kernel_rhs (real De, real K[DIM][DIM], real O[DIM][DIM], real B[DIM][DIM], real M[DIM][DIM], real RHS[DIM][DIM]); 

// Calculate the Kernel matrix
void hig_flow_calculate_kernel (higflow_solver *ns, real lambda[DIM], real R[DIM][DIM], real Kernel[DIM][DIM], real tol); 

// Calculate the Omega matrix
void hig_flow_calculate_omega (real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real Omega[DIM][DIM], real small); 

// Calculate the matrix BB
void hig_flow_calculate_b (higflow_solver *ns, real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real BB[DIM][DIM], real tol); 

// Calculate the matrix MM for Oldroyd model
void hig_flow_calculate_m_oldroyd (higflow_solver *ns, real lambda[DIM], real M[DIM][DIM], real M_aux[DIM][DIM], real tol); 

// Calculate the matrix MM for Giesekus model
void hig_flow_calculate_m_giesekus (higflow_solver *ns, real lambda[DIM], real M[DIM][DIM], real M_aux[DIM][DIM], real tol); 

// Calculate the matrix MM for LPTT model
void hig_flow_calculate_m_lptt (higflow_solver *ns, real tr, real lambda[DIM], real M[DIM][DIM], real R[DIM][DIM], real M_aux[DIM][DIM], real tol); 

// Get the velocity at cell center 
void hig_flow_velocity_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, real u[DIM]); 

// Get the derivative of Kernel 
void hig_flow_derivative_kernel_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, int i, int j, real Kcenter, real dKdx[DIM]); 

// Calculate the kronecker product 
void hig_flow_kernel_system_matrix (real w[DIM*DIM+1][DIM*DIM+1], real Omega[DIM][DIM], real dt);

// Calculate RHS = 2B * M/De
void hig_flow_implicit_kernel_rhs (real De, real B[DIM][DIM], real M[DIM][DIM], real RHS[DIM][DIM]);

// Solve linear system for constitutive equation
void hig_flow_solve_system_constitutive_equation ( int n, real A[DIM*DIM][DIM*DIM+1], real x[DIM*DIM] ); 

//Calculate convective tensor term CUBISTA
real hig_flow_convective_tensor_term_cubista(higflow_solver *ns, distributed_property *dpu, sim_domain *sdp, sim_stencil *stn, real K[DIM][DIM], Point ccenter, Point cdelta, int dim, int i, int j);

// Computing the Kernel Tensor
void higflow_compute_kernel_tensor(higflow_solver *ns);

// Constitutive Equation Step for the Explicit Euler Method
void higflow_explicit_euler_constitutive_equation(higflow_solver *ns);

// Constitutive Equation Step for the Implicit Euler Method
void higflow_implicit_euler_constitutive_equation(higflow_solver *ns);

// Computing the Polymeric Tensor
void higflow_compute_polymeric_tensor(higflow_solver *ns);

#endif
