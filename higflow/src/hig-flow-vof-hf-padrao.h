// The standard height-function method: sum the volume fraction along a column of
// cells to get the interface height, then differentiate the heights to get curvature.
//
// Integrating before differentiating is what makes this the most accurate curvature
// estimator here for a WELL-RESOLVED interface -- and what makes it fail when the
// column does not span the interface cleanly, which is what the adaptive variant in
// hig-flow-vof-adap-hf.h addresses.
//
// PART OF THE VOF FAMILY: a volume fraction per cell says how much of each phase is
// there, and the interface is RECONSTRUCTED from it rather than tracked.  The steps
// are: estimate a normal, place an interface from the normal and the fraction, then
// advect the fraction.  These files are the alternative methods for those steps, not
// a pipeline -- several are alternatives to each other.

#ifndef HIG_FLOW_VOF_HF_PADRAO
#define HIG_FLOW_VOF_HF_PADRAO

#include "hig-flow-discret.h"

real direction_hf_horizontal(sim_domain *sdp, higflow_solver *ns, Point center, Point delta, int *aux_dir);

real direction_hf_vertical(sim_domain *sdp, higflow_solver *ns, Point center, Point delta, int *aux_dir);

void vertical_collumn_padrao(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *vertical, int *aux);

void horizontal_row_padrao(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *horizontal, int *aux);

void calculate_interfacial_force(sim_domain *sdp, higflow_solver *ns, int clid, Point center, Point IF);

void higflow_compute_curvature_and_normal_with_HF_2D(higflow_solver *ns);

#endif
