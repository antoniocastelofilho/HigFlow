// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Step multiphase - version 10/11/2016
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_STEP_MULTIPHASE_VISC
#define HIG_FLOW_STEP_MULTIPHASE_VISC

#include "hig-flow-step.h"
#include "hig-flow-vof-plic.h"
#include "hig-flow-vof-adap-hf.h"
#include "hig-flow-mittag-leffler.h"
#include "hig-flow-step-multiphase.h"
#include "hig-flow-step-viscoelastic.h"

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

// Calculate convective tensor term CUBISTA
real hig_flow_conv_tensor_term_cub_mult_visc(higflow_solver *ns, distributed_property *dpu, sim_domain *sdp, sim_stencil *stn, real K[DIM][DIM], Point ccenter, Point cdelta, int dim, int i, int j);

// Constitutive Equation Step for the Implicit Euler Method
void higflow_implicit_euler_constitutive_equation_multiphase_viscoelastic(higflow_solver *ns);

// Calculate RHS = OK - KO + 2B * M/De
void hig_flow_kernel_rhs_multiphase_viscoelastic(real De, real K[DIM][DIM], real O[DIM][DIM], real B[DIM][DIM], real M[DIM][DIM], real RHS[DIM][DIM]) ;

// One step of the Navier-Stokes the projection method
void higflow_solver_step_multiphase_viscoelastic(higflow_solver *ns); 

#endif
