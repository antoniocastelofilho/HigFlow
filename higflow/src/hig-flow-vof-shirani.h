// Shirani normal estimation over a 5x5x5 neighbourhood.
//
// DEAD IN THIS TREE: nothing includes this header, and the only calls to
// shirani_125_cells() outside its own .c are commented out in hig-flow-vof-HF-3D.c.
// The function is compiled and linked, and reachable only if someone wires it up.
// The commented-out declaration below is copy-paste shared with
// hig-flow-vof-mehta.h; it is not a second entry point.
//
// PART OF THE VOF FAMILY: a volume fraction per cell says how much of each phase is
// there, and the interface is RECONSTRUCTED from it rather than tracked.  The steps
// are: estimate a normal, place an interface from the normal and the fraction, then
// advect the fraction.  These files are the alternative methods for those steps, not
// a pipeline -- several are alternatives to each other.

#ifndef HIG_FLOW_VOF_SHIRANI
#define HIG_FLOW_VOF_SHIRANI

#include "hig-flow-discret.h"

//void shirani_125_cells(sim_domain *sdp, higflow_solver *ns, int clid,Point center, Point p, Point delta);
void shirani_125_cells(higflow_solver *ns);

#endif
