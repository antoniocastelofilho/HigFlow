// *******************************************************************
//  HiG-Flow Solver Step Multiphase Electroosmotic - version 07/05/2024
// *******************************************************************

#ifndef HIG_FLOW_STEP_MULTIPHASE_EO
#define HIG_FLOW_STEP_MULTIPHASE_EO

#include "hig-flow-step-multiphase.h"
#include "hig-flow-step-electroosmotic.h"


// *******************************************************************
// PB and PNP Equations
// *******************************************************************

void higflow_explicit_euler_ionic_transport_equation_nplus_multiphase(higflow_solver *ns);

void higflow_explicit_euler_ionic_transport_equation_nminus_multiphase(higflow_solver *ns);

void higflow_semi_implicit_euler_ionic_transport_equation_nplus_multiphase(higflow_solver *ns);

void higflow_semi_implicit_euler_ionic_transport_equation_nminus_multiphase(higflow_solver *ns);

void higflow_multiphase_electroosmotic_psi(higflow_solver *ns);

void higflow_calculate_multiphase_electroosmotic_source_term( higflow_solver *ns);


// *******************************************************************
// Navier-Stokes Step
// *******************************************************************

void higflow_explicit_euler_intermediate_velocity_multiphase_electroosmotic(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]);

void higflow_explicit_runge_kutta_2_intermediate_velocity_multiphase_electroosmotic(higflow_solver *ns);

void higflow_explicit_runge_kutta_3_intermediate_velocity_multiphase_electroosmotic(higflow_solver *ns);

void higflow_semi_implicit_euler_intermediate_velocity_multiphase_electroosmotic(higflow_solver *ns);

void higflow_semi_implicit_crank_nicolson_intermediate_velocity_multiphase_electroosmotic(higflow_solver *ns);

void higflow_semi_implicit_bdf2_intermediate_velocity_multiphase_electroosmotic(higflow_solver *ns);


// *******************************************************************
// One step of multiphase electroosmotic solver
// *******************************************************************

void higflow_solver_step_multiphase_electroosmotic(higflow_solver *ns);



#endif