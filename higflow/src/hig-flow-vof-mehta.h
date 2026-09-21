// Mehta's interface treatment, called from the multiphase step.
//
// Unlike hig-flow-vof-shirani.h, which is its copy-paste sibling, this one IS live:
// hig-flow-step-multiphase.c calls mehta().  The commented-out shirani declaration
// below is a leftover of that copy and means nothing here.
//
// PART OF THE VOF FAMILY: a volume fraction per cell says how much of each phase is
// there, and the interface is RECONSTRUCTED from it rather than tracked.  The steps
// are: estimate a normal, place an interface from the normal and the fraction, then
// advect the fraction.  These files are the alternative methods for those steps, not
// a pipeline -- several are alternatives to each other.

#ifndef HIG_FLOW_VOF_MEHTA
#define HIG_FLOW_VOF_MEHTA

#include "hig-flow-discret.h"

//void shirani_125_cells(sim_domain *sdp, higflow_solver *ns, int clid,Point center, Point p, Point delta);
void mehta(higflow_solver *ns);

#endif
