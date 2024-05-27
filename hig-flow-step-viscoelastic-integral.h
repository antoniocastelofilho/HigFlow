// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Step Viscoelastic - version 10/11/2016
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_STEP_VISCOELASTIC_INTEGRAL
#define HIG_FLOW_STEP_VISCOELASTIC_INTEGRAL

#include "hig-flow-step-generalized-newtonian.h"
#include "hig-flow-mittag-leffler.h"

// *******************************************************************
// Constitutive Equation Step for the Explicit Euler Method
// *******************************************************************
void higflow_explicit_euler_constitutive_equation_integral(higflow_solver *ns); 

void hig_flow_integral_equation (higflow_solver *ns);

void hig_flow_integral_equation_KBKZ (higflow_solver *ns);

void hig_flow_integral_equation_KBKZ_Fractional (higflow_solver *ns);
 
// Calculate RHS = Du^t B + B Du
void hig_flow_b_rhs (real B[DIM][DIM], real Du[DIM][DIM], real RHS[DIM][DIM]); 


// Get the velocity at cell center 
void hig_flow_velocity_at_center_cell_integral (higflow_solver *ns, Point ccenter, Point cdelta, real u[DIM]); 

// Get the derivative of B tensor 
void hig_flow_derivative_b_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, int k, int i, int j, real Bcenter, real dBdx[DIM]); 


// *******************************************************************
// Calculate convective tensor term CUBISTA
// *******************************************************************
real hig_flow_convective_tensor_term_b_cubista(higflow_solver *ns, distributed_property *dpu, sim_domain *sdp, sim_stencil *stn, real B[DIM][DIM], Point ccenter, Point cdelta, int dim, int k, int i, int j); 


// *******************************************************************
// Navier-Stokes Step for the Explicit Euler Method
// *******************************************************************
void higflow_explicit_euler_intermediate_velocity_viscoelastic_integral(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]); 


// *******************************************************************
// Navier-Stokes Step for the Second Order Explicit Runge-Kutta Method
// *******************************************************************
void higflow_explicit_runge_kutta_2_intermediate_velocity_viscoelastic_integral(higflow_solver *ns); 

// *******************************************************************
// Navier-Stokes Step for third Order Explicit Runge-Kutta Method
// *******************************************************************
void higflow_explicit_runge_kutta_3_intermediate_velocity_viscoelastic_integral(higflow_solver *ns); 

// *******************************************************************
// Navier-Stokes Step for the Implicit Euler Method
// *******************************************************************
void higflow_semi_implicit_euler_intermediate_velocity_viscoelastic_integral(higflow_solver *ns); 


// *******************************************************************
// Navier-Stokes Step for the Semi-Implicit Crank-Nicolson Method
// *******************************************************************
void higflow_semi_implicit_crank_nicolson_intermediate_velocity_viscoelastic_integral(higflow_solver *ns); 

// *******************************************************************
// Navier-Stokes Step for the Implicit BDF2 Method
// *******************************************************************
void higflow_semi_implicit_bdf2_intermediate_velocity_viscoelastic_integral(higflow_solver *ns); 

// One step of the Navier-Stokes the projection method
void higflow_solver_step_viscoelastic_integral(higflow_solver *ns); 

// Print the Polymeric Tensor
void higflow_print_polymeric_tensor_integral(higflow_solver *ns);


#endif
