// Integral viscoelastic models (KBKZ): stress from an integral over deformation
// history rather than from a differential equation, which is why this path carries
// the Finger tensor and a stored history the differential ones do not have.
//
// REACHED ONLY BY example2d_KBKZ, and that case is marked known_broken in the suite.
// This is therefore the one solver family with an example that does not currently
// produce a usable result -- in practice it is uncovered.
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
