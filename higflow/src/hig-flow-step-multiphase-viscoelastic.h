// *******************************************************************
//  HiG-Flow Solver Step multiphase - version 28/03/2024
// *******************************************************************

#ifndef HIG_FLOW_STEP_MULTIPHASE_VISC
#define HIG_FLOW_STEP_MULTIPHASE_VISC

#include "hig-flow-step-multiphase.h"
#include "hig-flow-step-viscoelastic.h"

// *******************************************************************
// Navier-Stokes Step
// *******************************************************************
// Computing beta viscoelastic
real higflow_interp_beta_multiphase_viscoelastic(real beta0, real beta1, real visc0, real visc1, real fracvol);

// Computing the Deborah Number
real higflow_interp_De_multiphase_viscoelastic(real De0, real De1, real fracvol);

// Computing the rhs model specific matrix
void higlow_interp_MM_multiphase_viscoelastic(real MM0[DIM][DIM], real MM1[DIM][DIM], real MM[DIM][DIM], real fracvol);

// Computing the FENE function
real higflow_interp_fA_multiphase_viscoelastic(real fA0, real fA1, real fracvol);

// Computing the FENE-P a parameter
real higflow_interp_a_multiphase_viscoelastic(real a0, real a1, real fracvol);

// Computing the Polymeric Tensor
void higflow_compute_polymeric_tensor_multiphase_viscoelastic(higflow_solver *ns) ;

// Constitutive Equation Step for the Explicit Euler Method
void higflow_explicit_euler_constitutive_equation_multiphase_viscoelastic(higflow_solver *ns);

// Constitutive Equation Step for the Implicit Euler Method
void higflow_implicit_euler_constitutive_equation_multiphase_viscoelastic(higflow_solver *ns);

// One step of the Navier-Stokes the projection method
void higflow_solver_step_multiphase_viscoelastic(higflow_solver *ns); 

#endif
