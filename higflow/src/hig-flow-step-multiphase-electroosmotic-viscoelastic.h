// Composition layer only: the three-way combination (two-phase + electroosmotic +
// viscoelastic), declaring its step function and delegating to the parents.
//
// NO EXAMPLE SELECTS THIS COMBINATION.  Uncovered by the suite, and the deepest
// composition in the tree -- the path most likely to read a cc field that none of
// its parents filled.  See rule 2 in hig-flow-kernel.h.

// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Step Multiphase Electro-osmotic Viscoelastic - version 12/05/2024
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_STEP_MULTIPHASE_ELECTROOSMOTIC_VISCOELASTIC
#define HIG_FLOW_STEP_MULTIPHASE_ELECTROOSMOTIC_VISCOELASTIC

#include "hig-flow-step-multiphase-electroosmotic.h"
#include "hig-flow-step-multiphase-viscoelastic.h"

// *******************************************************************

// One step of the Navier-Stokes the projection method
void higflow_solver_step_multiphase_electroosmotic_viscoelastic(higflow_solver *ns); 

#endif