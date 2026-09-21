// Normal estimation from the 3x3 neighbourhood (Shirani's method at the smallest
// stencil).  Used by the height-function files as a fallback where their columns do
// not resolve.
//
// THIS HEADER IS NOT SELF-CONTAINED: it includes nothing, yet its declaration names
// higflow_solver, sim_domain and Point.  It compiles only when included AFTER the
// headers that define them, so it cannot be the first include in a file.  That is a
// property of this header, not a convention -- include it late.
//
// PART OF THE VOF FAMILY: a volume fraction per cell says how much of each phase is
// there, and the interface is RECONSTRUCTED from it rather than tracked.  The steps
// are: estimate a normal, place an interface from the normal and the fraction, then
// advect the fraction.  These files are the alternative methods for those steps, not
// a pipeline -- several are alternatives to each other.

#ifndef HIG_FLOW_VOF_9_CELLS
#define HIG_FLOW_VOF_9_CELLS

void higflow_compute_normal_multiphase_2D_shirani_9_cells(higflow_solver *ns, sim_domain *sdp, int clid, Point center, Point p, Point delta);

#endif
