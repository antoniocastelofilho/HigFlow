// *******************************************************************
// *******************************************************************
//  HiG-Flow Time Step Controller - version 01/08/24
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_TIMESTEP
#define HIG_FLOW_TIMESTEP

#include "hig-flow-kernel.h"
#include "hig-flow-eval.h"

#define MAX_CFL 0.1
#define MAX_DT_FACTOR 1.1

real higflow_compute_CFL(higflow_solver *ns);

real higflow_compute_CFL_PNP(higflow_solver *ns);

real higflow_compute_CFL_PNP_multiphase(higflow_solver *ns);

void higflow_adjust_timestep(higflow_solver *ns);

#endif // HIG_FLOW_TIMESTEP