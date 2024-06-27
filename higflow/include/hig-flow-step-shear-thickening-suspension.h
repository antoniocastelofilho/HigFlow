// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Step Viscoelastic - version 10/11/2016
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_STEP_SHEAR_THICKENING_SUSPENSION
#define HIG_FLOW_STEP_SHEAR_THICKENING_SUSPENSION

#include "hig-flow-step-generalized-newtonian.h"

// *******************************************************************
// Navier-Stokes Step
// *******************************************************************

// One step of the explicit Euler method - Intermediate Velocity
void higflow_explicit_euler_intermediate_velocity_shear_thickening_suspensions(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *pdustar[DIM]); 

// One step of the explicit order 2 Runge-Kutta method - Intermediate Velocity
void higflow_explicit_runge_kutta_2_intermediate_velocity_shear_thickening_suspensions(higflow_solver *ns);

// One step of the explicit order 3 Runge-Kutta method - Intermediate Velocity
void higflow_explicit_runge_kutta_3_intermediate_velocity_shear_thickening_suspensions(higflow_solver *ns);

// One step of the Semi-Implicit Euler method - Intermediate Velocity
void higflow_semi_implicit_euler_intermediate_velocity_shear_thickening_suspensions(higflow_solver *ns); 

// One step of the semi-implicit Crank-Nicolson method - Intermediate Velocity
void higflow_semi_implicit_crank_nicolson_intermediate_velocity_shear_thickening_suspensions(higflow_solver *ns);

// One step of the semi-implicit BDF2 method - Intermediate Velocity
void higflow_semi_implicit_bdf2_intermediate_velocity_shear_thickening_suspensions(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]); 

// One step of the Navier-Stokes the projection method
void higflow_solver_step_shear_thickening_suspensions(higflow_solver *ns); 

// Calculate the eige-value and eige-vectors using the Jacobi method
void hig_flow_jacobi(real A[DIM][DIM], real d[DIM], real V[DIM][DIM]); 

// Calculate the matrix product
void hig_flow_matrix_product( real A[DIM][DIM], real R[DIM][DIM], real B[DIM][DIM]); 

// Calculate the matrix product transpose
void hig_flow_matrix_transpose_product( real A[DIM][DIM], real R[DIM][DIM], real B[DIM][DIM]); 

//Step function 1
real step_function_1(real l);

//Step function 2
real step_function_2(real l);

// Calculates the tensors E_E and E_C 
void hig_flow_compressive_extensional_tensors (real lambda[DIM], real R[DIM][DIM], real E_C[DIM][DIM], real E_E[DIM][DIM]);

// Calculate the fourth-order orientation moment <pppp>
void hig_flow_calculate_4th_order_orientation_moment (real AM, real KD[DIM][DIM], real A[DIM][DIM], real FOM[DIM][DIM][DIM][DIM]);

// Calculate the suspension stress tensor
void hig_flow_calculate_suspension_stress_tensor (real tau_M[DIM][DIM], real FOM[DIM][DIM][DIM][DIM], real T[DIM][DIM]);

// Calculate the double dot product between two tensors (one is a second order tensor and the other one is a 4th order tensor)
void hig_flow_calculate_double_dot_product_tensors (real M[DIM][DIM], real FOM[DIM][DIM][DIM][DIM], real P[DIM][DIM]);

// Calculate the jamming coordinate value chi
real hig_flow_calculate_jamming_coordinate_chi (real A[DIM][DIM], real E_C[DIM][DIM]);

// Get the velocity at cell center 
void hig_flow_velocity_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, real u[DIM]); 

// Get the derivative of microstructure tensor A
void hig_flow_derivative_tensor_A_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, int i, int j, real Acenter, real dAdx[DIM]);

// Calculate the matrix product
void hig_flow_kernel_system_matrix_shear_thickening_suspensions (real w[DIM*DIM][DIM*DIM+1], real Omega[DIM][DIM], real dt);

// Calculate RHS = A*Du + (Du)^T * A- 2*(Du)^T : <pppp> - beta*( E_E:<pppp> + (phi/15)*(2*E_C+Tr(E_C)*KD) )
void hig_flow_evolution_equation_microstructure_tensor_rhs (real beta, real phi, real A[DIM][DIM], real DU[DIM][DIM], real E_C[DIM][DIM], real TrEC, real KD[DIM][DIM], real FOM[DIM][DIM][DIM][DIM], real LP[DIM][DIM], real EP[DIM][DIM], real RHS[DIM][DIM]);

// Calculate RHS for the implicit method RHS = - 2*(Du)^T : <pppp> - beta*( E_E:<pppp> + (phi/15)*(2*E_C+Tr(E_C)*KD) )
void hig_flow_implicit_evolution_equation_microstructure_tensor_rhs (real beta, real phi, real A[DIM][DIM], real E_C[DIM][DIM], real TrEC, real KD[DIM][DIM], real FOM[DIM][DIM][DIM][DIM], real LP[DIM][DIM], real EP[DIM][DIM], real RHS[DIM][DIM]);

// Solve linear system for constitutive equation
void hig_flow_solve_system_constitutive_equation ( int n, real A[DIM*DIM][DIM*DIM+1], real x[DIM*DIM] ); 

//Calculate convective tensor term CUBISTA
real hig_flow_convective_tensor_A_term_cubista(higflow_solver *ns, distributed_property *dpu, sim_domain *sdp, sim_stencil *stn, real K[DIM][DIM], Point ccenter, Point cdelta, int dim, int i, int j);
//real hig_flow_convective_tensor_term_cubista(higflow_solver *ns, distributed_property *dpu, sim_domain *sdp, sim_stencil *stn, real K[DIM][DIM], Point ccenter, Point cdelta, int dim, int i, int j);

// Solve the evolution equation of the microstructure tensor <pp>
void higflow_explicit_euler_evolution_equation_microstructure_tensor(higflow_solver *ns);

// Solve the evolution equation of the microstructure tensor <pp>
void higflow_implicit_euler_evolution_equation_microstructure_tensor(higflow_solver *ns);

// Get the derivative of the volume fraction 
void hig_flow_derivative_volfrac_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, real ncenter, real dndx[DIM]);

// Calculate convective term CUBISTA for volume fraction
real hig_flow_fraction_volume_suspensions_term_cubista(higflow_solver *ns, distributed_property *dpu, sim_domain *sdp, sim_stencil *stn, real varphi, Point ccenter, Point cdelta, int dim);

// Volume fraction evolution equation
void higflow_explicit_euler_volume_fraction_equation(higflow_solver *ns);

// Computing the particle stress TensorTau
void higflow_compute_particle_stress_tensor_shear_thickening_suspension(higflow_solver *ns);

#endif
