// Adaptive time step: measure the CFL, then grow or shrink dt.
//
// MAX_CFL (0.25) and MAX_DT_FACTOR (1.05) are COMPILE-TIME, not YAML settings: the
// stability target and the per-step growth cap cannot be changed per case without
// rebuilding.  The asymmetry is deliberate -- dt grows by at most 5% per step and
// may be cut without limit.
//
// There is one CFL variant per physics that adds its own transport speed (PNP adds
// ionic migration, the multiphase variant the interface velocity).  A path that uses
// the plain variant while carrying a faster mechanism gets a dt that satisfies the
// wrong stability condition.

// *******************************************************************
// *******************************************************************
//  HiG-Flow Time Step Controller - version 01/08/24
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_TIMESTEP
#define HIG_FLOW_TIMESTEP

#include "hig-flow-kernel.h"
#include "hig-flow-eval.h"

#define MAX_CFL 0.25
#define MAX_DT_FACTOR 1.05

real higflow_compute_CFL(higflow_solver *ns);

real higflow_compute_CFL_PNP(higflow_solver *ns);

real higflow_compute_CFL_PNP_multiphase(higflow_solver *ns);

void higflow_adjust_timestep(higflow_solver *ns);

#endif // HIG_FLOW_TIMESTEP