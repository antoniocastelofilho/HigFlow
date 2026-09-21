// Composition layer only: it declares the single step function that runs the
// electroosmotic path with a viscoelastic rheology, and delegates everything to the
// two parent step files.  There is no formulation of its own here.
//
// NO EXAMPLE SELECTS THIS COMBINATION.  Uncovered by the suite.

// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Step Electro-osmotic Viscoelastic - version 03/2024
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_STEP_ELECTROOSMOTIC_VISCOELASTIC
#define HIG_FLOW_STEP_ELECTROOSMOTIC_VISCOELASTIC

#include "hig-flow-step-viscoelastic.h"
#include "hig-flow-step-electroosmotic.h"

// *******************************************************************

// One step of the Navier-Stokes the projection method
void higflow_solver_step_electroosmotic_viscoelastic(higflow_solver *ns); 

#endif
