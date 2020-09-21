// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Terms - version 10/11/2016
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_TERMS
#define HIG_FLOW_TERMS

#include "hig-flow-kernel.h"
#include "hig-flow-eval.h"
#include "hig-flow-discret.h"

// *******************************************************************
// Navier-Stokes Equation Components
// *******************************************************************

// Cell term contribution for the interfacial tension
real higflow_interfacial_tension_term(higflow_solver *ns);

// Cell term contribution for the gravity
real higflow_gravity_term(higflow_solver *ns);

// Pressure term contribution for the Navier-Stokes equation
real higflow_pressure_term(higflow_solver *ns); 

// Source term contribution for the Navier-Stokes equation
real higflow_source_term(higflow_solver *ns); 

// Tensor term contribution for the Navier-Stokes equation
real higflow_tensor_term(higflow_solver *ns); 

// Convective term: Quick scheme
real higflow_convective_term_quick(higflow_solver *ns, Point delta, int dim); 

// Convective term: Cubista scheme
real higflow_convective_term_cubista(higflow_solver *ns, Point delta, int dim); 

// Convective term: Modified Coefficient Upwind scheme
real higflow_convective_term_mcupwind(higflow_solver *ns, Point delta, int dim);

// Convective term: second order scheme
real higflow_convective_term_second_order(higflow_solver *ns, Point delta, int dim); 

// Convective term: first order scheme
real higflow_convective_term_first_order(higflow_solver *ns, Point delta, int dim); 

// Convective term: central scheme
real higflow_convective_term_central(higflow_solver *ns, Point delta, int dim); 

// Convective term contribution for the Navier-Stokes equation
real higflow_convective_term(higflow_solver *ns, Point delta, int dim); 

// Difusive term contribution for the Navier-Stokes equation
real higflow_difusive_term(higflow_solver *ns, Point delta); 

// Cell electroosmotic source term contribution for the Navier-Stokes equation
real higflow_electroosmotic_source_term(higflow_solver *ns);

// Cell electroosmotic diffusive ionic term contribution for the Navier-Stokes equation
real higflow_diffusive_ionic_term(higflow_solver *ns); 

// Cell electroosmotic potential ionic term contribution for the Navier-Stokes equation
real higflow_potential_ionic_term(higflow_solver *ns); 
#endif
