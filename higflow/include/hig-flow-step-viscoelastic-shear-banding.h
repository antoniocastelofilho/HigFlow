// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Step Electro-osmotic - version 02/02/2018
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_STEP_VISCOELASTIC_SHEAR_BANDING
#define HIG_FLOW_STEP_VISCOELASTIC_SHEAR_BANDING

#include "hig-flow-step-generalized-newtonian.h"

// *******************************************************************

// Get velocity at center cell
void hig_flow_velocity_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, real u[DIM]); 

// Navier-Stokes Step for the Explicit Euler Method
void higflow_explicit_euler_intermediate_velocity_viscoelastic_shear_banding(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]); 

// Navier-Stokes Step for the Second Order Explicit Runge-Kutta Method
void higflow_explicit_runge_kutta_2_intermediate_velocity_shear_banding(higflow_solver *ns); 

// Navier-Stokes Step for third Order Explicit Runge-Kutta Method
void higflow_explicit_runge_kutta_3_intermediate_velocity_shear_banding(higflow_solver *ns); 

// Navier-Stokes Step for the Implicit Euler Method
void higflow_semi_implicit_euler_intermediate_velocity_shear_banding(higflow_solver *ns); 

// Navier-Stokes Step for the Semi-Implicit Crank-Nicolson Method
void higflow_semi_implicit_crank_nicolson_intermediate_velocity_shear_banding(higflow_solver *ns); 

// Navier-Stokes Step for the Implicit BDF2 Method
void higflow_semi_implicit_bdf2_intermediate_velocity_shear_banding(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]); 

// One step of the Navier-Stokes the projection method for viscoelastic flow
void higflow_solver_step_viscoelastic_shear_banding(higflow_solver *ns) ; 

//********************************************************************
// Viscoelastic functions
//********************************************************************


// Calculate the eige-value and eige-vectors using the Jacobi method
void hig_flow_jacobi(real A[DIM][DIM], real d[DIM], real V[DIM][DIM]); 

// Calculate the matrix product
void hig_flow_matrix_product( real A[DIM][DIM], real R[DIM][DIM], real B[DIM][DIM]); 

// Calculate the matrix product transpose
void hig_flow_matrix_transpose_product( real A[DIM][DIM], real R[DIM][DIM], real B[DIM][DIM]); 

// Calculate the Omega matrix
void hig_flow_calculate_omega (real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real Omega[DIM][DIM], real small);

// Calculate the matrix BB
void hig_flow_calculate_b (real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real BB[DIM][DIM]);

// Calculate the kronecker product 
void hig_flow_kernel_system_matrix_VCM (real w[DIM*DIM+1][DIM*DIM+1], real Omega[DIM][DIM], real dt);

// Calculate the matrix product
//void hig_flow_kernel_system_matrix_VCM (real w[DIM*DIM][DIM*DIM+1], real Omega[DIM][DIM], real OmegaT[DIM][DIM], real dt);

// Solve linear system for constitutive equation
void hig_flow_solve_system_constitutive_equation ( int n, real A[DIM*DIM][DIM*DIM+1], real x[DIM*DIM] ); 


//Viscoelastic fluids with shear-banding functions

// Computing the Polymeric Tensor
void higflow_compute_polymeric_tensor_shear_banding(higflow_solver *ns);

// Computing the breakage rate of specie A (cA)
void higflow_compute_viscoelastic_shear_banding_cA_VCM(higflow_solver *ns);

// Computing the breakage rate of specie B (cB)
void higflow_compute_viscoelastic_shear_banding_cB_VCM(higflow_solver *ns);

// Solve the Constitutive Equation of specie A using the Explicit Euler Method
void higflow_explicit_euler_conformation_tensor_A(higflow_solver *ns);

// Solve the Constitutive Equation of specie B using the Explicit Euler Method
void higflow_explicit_euler_conformation_tensor_B(higflow_solver *ns);

// Solve the Constitutive Equation of specie A using the Implicit Euler Method
void higflow_implicit_euler_conformation_tensor_A(higflow_solver *ns);

// Solve the Constitutive Equation of specie B using the Implicit Euler Method
void higflow_implicit_euler_conformation_tensor_B(higflow_solver *ns);

// Calculate convective tensor term CUBISTA (specie A)
real hig_flow_convective_tensor_A_term_cubista(higflow_solver *ns, distributed_property *dpu, sim_domain *sdp, sim_stencil *stn, real K[DIM][DIM], Point ccenter, Point cdelta, int dim, int i, int j);

// Calculate convective tensor term CUBISTA (specie B)
real hig_flow_convective_tensor_B_term_cubista(higflow_solver *ns, distributed_property *dpu, sim_domain *sdp, sim_stencil *stn, real K[DIM][DIM], Point ccenter, Point cdelta, int dim, int i, int j);

//Computing the terms needed for the simulation of viscoelastic flows that exhibit shear-banding behaviour using the VCM model
real higflow_computational_cell_conformation_tensor_A_shear_banding_VCM_model(higflow_solver *ns, sim_domain *sdp, sim_stencil *stn, real A[DIM][DIM], Point ccenter, Point cdelta, int dim, int i, int j);

//Computing the terms needed for the simulation of viscoelastic flows that exhibit shear-banding behaviour using the VCM model
real higflow_computational_cell_conformation_tensor_B_shear_banding_VCM_model(higflow_solver *ns, sim_domain *sdp, sim_stencil *stn, real B[DIM][DIM], Point ccenter, Point cdelta, int dim, int i, int j);

// Calculate RHS = A*Du + (Du)*T A + (nA*I-A + CB*nB*B -CA*A)/De
void hig_flow_conformation_tensor_A_VCM_rhs (real De, real A[DIM][DIM], real DU[DIM][DIM], real B[DIM][DIM], real CA, real CB, real nA, real nB, real RHS[DIM][DIM]);

// Calculate RHS = B*Du + (Du)*T B + (nB*I/2-B)/De +2*(CA*A-CB*nB*B)/DeB
void hig_flow_conformation_tensor_B_VCM_rhs (real DeA, real DeB, real B[DIM][DIM], real DU[DIM][DIM], real A[DIM][DIM], real CA, real CB, real nA, real nB, real RHS[DIM][DIM]);

// Calculate RHS = (nA*I-A + CB*nB*B -CA*A)/De
void hig_flow_implicit_conformation_tensor_A_VCM_rhs (real De, real A[DIM][DIM], real B[DIM][DIM], real CA, real CB, real nA, real nB, real RHS[DIM][DIM]);

// Calculate RHS = (nB*I/2-B)/De +2*(CA*A-CB*nB*B)/DeB
void hig_flow_implicit_conformation_tensor_B_VCM_rhs (real DeA, real DeB, real B[DIM][DIM], real A[DIM][DIM], real CA, real CB, real nA, real nB, real RHS[DIM][DIM]);

// Get the derivative of conformation tensor A
void hig_flow_derivative_tensor_A_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, int i, int j, real Acenter, real dAdx[DIM]);

// Get the derivative of conformation tensor B
void hig_flow_derivative_tensor_B_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, int i, int j, real Bcenter, real dBdx[DIM]);

// Get the second derivative of conformation tensor
//void hig_flow_second_derivative_tensor_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, int dim, int i, int j, real Acenter, real d2Adx2[DIM][DIM]);

// Transport equation of the density numbers of specie A using Euler explicit method
void higflow_explicit_euler_shear_banding_transport_equation_nA(higflow_solver *ns);

// Transport equation of the density numbers of specie B using Euler explicit method
void higflow_explicit_euler_shear_banding_transport_equation_nB(higflow_solver *ns);

// Transport equation of the density numbers of specie A using Euler implicit method
void higflow_implicit_euler_shear_banding_transport_equation_nA(higflow_solver *ns);

// Transport equation of the density numbers of specie B using Euler implicit method
void higflow_implicit_euler_shear_banding_transport_equation_nB(higflow_solver *ns); 

// Calculate RHS = cB*nB**2/(2*DeA) -cA*nA/DeA
real hig_flow_nA_VCM_RHS (real CA, real CB, real nA, real nB, real De);

// Calculate RHS = 2*cA*nA/DeA - cB*nB**2/DeA
real hig_flow_nB_VCM_RHS (real CA, real CB, real nA, real nB, real De);

// Get the derivative of nA 
void hig_flow_derivative_nA_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, real ncenter, real dndx[DIM]);

// Get the derivative of nB 
void hig_flow_derivative_nB_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, real ncenter, real dndx[DIM]);

// Calculate convective term of the transport equation of species A and B using the CUBISTA method
real higflow_convective_term_shear_banding_VCM_cubista(higflow_solver *ns, distributed_property *dpu, distributed_property *dpn, sim_domain *sdp, sim_stencil *stn, real n, Point ccenter, Point cdelta, int dim);



#endif
