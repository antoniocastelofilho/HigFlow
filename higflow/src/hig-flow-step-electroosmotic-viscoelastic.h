// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Step Electro-osmotic - version 02/02/2018
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_STEP_ELECTROOSMOTIC_VISCOELASTIC
#define HIG_FLOW_STEP_ELECTROOSMOTIC_VISCOELASTIC

#include "hig-flow-step-generalized-newtonian.h"

// *******************************************************************
// Electro-osmotic source term

//Semi-implicit ionic transport equation for nplus
void higflow_implicit_euler_ionic_transport_equation_nplus(higflow_solver *ns); 

//Semi-implicit ionic transport equation for nminus
void higflow_implicit_euler_ionic_transport_equation_nminus(higflow_solver *ns);

// Get velocity at center cell
void hig_flow_velocity_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, real u[DIM]); 

// Time step of ionic transport equation for n+
void higflow_explicit_euler_ionic_transport_equation_nplus(higflow_solver *ns); 

// Time step of ionic transport equation for n-
void higflow_explicit_euler_ionic_transport_equation_nminus(higflow_solver *ns); 

// Compute the convective ionic term
real higflow_convective_ionic_term_cubista(higflow_solver *ns, distributed_property *dpu, distributed_property *dpn, sim_domain *sdp, sim_stencil *stn, real n, Point ccenter, Point cdelta, int dim); 

// Compute the ionic concentration derivative
void hig_flow_derivative_nminus_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, real ncenter, real dndx[DIM]); 
void hig_flow_derivative_nplus_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, real ncenter, real dndx[DIM]) ;

// Souce term by analytical PBDH model
void higflow_calculate_electroosmotic_source_term_analytic_pbdh( higflow_solver *ns);

// Souce term by PBDH model
void higflow_calculate_electroosmotic_source_term_pbdh( higflow_solver *ns);

// Souce term by PB model
void higflow_calculate_electroosmotic_source_term_pb( higflow_solver *ns);

// Souce term by PNP model
void higflow_calculate_electroosmotic_source_term_pnp( higflow_solver *ns);

// Boundary conditiono for electro-osmotic source term
void higflow_boundary_condition_for_electroosmotic_source_term(higflow_solver *ns); 

// Electroosmotic induced potential psi 
void higflow_electroosmotic_psi(higflow_solver *ns);

// Electroosmotic applied potential psi 
void higflow_electroosmotic_phi(higflow_solver *ns);

// Navier-Stokes Step for the Explicit Euler Method
void higflow_explicit_euler_intermediate_velocity_electroosmotic(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]); 

// Navier-Stokes Step for the Second Order Explicit Runge-Kutta Method
void higflow_explicit_runge_kutta_2_intermediate_velocity_electroosmotic(higflow_solver *ns); 

// Navier-Stokes Step for third Order Explicit Runge-Kutta Method
void higflow_explicit_runge_kutta_3_intermediate_velocity_electroosmotic(higflow_solver *ns); 

// Navier-Stokes Step for the Implicit Euler Method
void higflow_semi_implicit_euler_intermediate_velocity_electroosmotic(higflow_solver *ns); 

// Navier-Stokes Step for the Semi-Implicit Crank-Nicolson Method
void higflow_semi_implicit_crank_nicolson_intermediate_velocity_electroosmotic(higflow_solver *ns); 

// Navier-Stokes Step for the Implicit BDF2 Method
void higflow_semi_implicit_bdf2_intermediate_velocity_electroosmotic(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]); 

// One step of the Navier-Stokes the projection method
void higflow_solver_step_electroosmotic_viscoelastic(higflow_solver *ns); 

// Apply the boundary condition for source term 
void higflow_boundary_condition_for_electroosmotic_nplus(higflow_solver *ns);
void higflow_boundary_condition_for_electroosmotic_nminus(higflow_solver *ns);

//********************************************************************
// Viscoelastic functions
//********************************************************************
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

// Calculate the matrix MM and BB for Oldroyd model
void hig_flow_calculate_m_and_b_oldroyd (higflow_solver *ns, real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real MM[DIM][DIM], real BB[DIM][DIM], real tol);

// Calculate the matrix MM and BB for Giesekus model
void hig_flow_calculate_m_and_b_giesekus (higflow_solver *ns, real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real MM[DIM][DIM], real BB[DIM][DIM], real tol);

// Calculate the matrix MM and BB for LPTT model
void hig_flow_calculate_m_and_b_lptt (higflow_solver *ns, real tr, real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real MM[DIM][DIM], real BB[DIM][DIM], real tol);

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


#endif
