// Two-phase electroosmotic flow: the ionic transport of the electroosmotic path
// carried across an interface that moves.
//
// NO EXAMPLE SELECTS THIS PATH.  Uncovered by the suite -- both of its parents are
// covered separately, their combination is not.
//
// WHAT A STEP FILE CONTAINS: the intermediate-velocity stage of the projection for
// one rheology, plus whatever constitutive equation that rheology has to advance.
// The shared stages (pressure solve, final velocity) stay in hig-flow-step.h.
//
// Each intermediate-velocity variant is offered in several time discretizations
// (explicit Euler, RK2, RK3, semi-implicit, implicit).  WHICH ONE RUNS COMES FROM
// THE YAML, so a variant present here is not necessarily a variant any case selects.

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

void higflow_semi_implicit_crank_nicolson_ionic_transport_equation_nplus_multiphase(higflow_solver *ns);

void higflow_semi_implicit_crank_nicolson_ionic_transport_equation_nminus_multiphase(higflow_solver *ns);

void higflow_semi_implicit_bdf2_ionic_transport_equation_nplus_multiphase(higflow_solver *ns);

void higflow_semi_implicit_bdf2_ionic_transport_equation_nminus_multiphase(higflow_solver *ns);

void higflow_multiphase_electroosmotic_phi(higflow_solver *ns);

real higflow_multiphase_electroosmotic_psi(higflow_solver *ns);

void higflow_multiphase_electroosmotic_solve_pb(higflow_solver *ns);

void higflow_calculate_multiphase_electroosmotic_source_term( higflow_solver *ns);


// *******************************************************************
// Navier-Stokes Step
// *******************************************************************

void higflow_explicit_euler_intermediate_velocity_multiphase_electroosmotic(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]);



void higflow_semi_implicit_euler_intermediate_velocity_multiphase_electroosmotic(higflow_solver *ns);

void higflow_semi_implicit_crank_nicolson_intermediate_velocity_multiphase_electroosmotic(higflow_solver *ns);

void higflow_semi_implicit_bdf2_intermediate_velocity_multiphase_electroosmotic(higflow_solver *ns);


// *******************************************************************
// One step of multiphase electroosmotic solver
// *******************************************************************

void higflow_solver_step_multiphase_electroosmotic(higflow_solver *ns);



#endif