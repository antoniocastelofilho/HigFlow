// Height function in 3D: columns along x, y or z, with direction_hf_* choosing which.
//
// The choice is not cosmetic -- a column nearly parallel to the interface gives a
// height that is not single-valued, and the curvature computed from it is wrong
// rather than merely imprecise.
//
// PART OF THE VOF FAMILY: a volume fraction per cell says how much of each phase is
// there, and the interface is RECONSTRUCTED from it rather than tracked.  The steps
// are: estimate a normal, place an interface from the normal and the fraction, then
// advect the fraction.  These files are the alternative methods for those steps, not
// a pipeline -- several are alternatives to each other.

#ifndef HIG_FLOW_VOF_HF_3D
#define HIG_FLOW_VOF_HF_3D

#include "hig-flow-discret.h"

void fraction_correction_at_get_3D(real *fracvol);

int get_frac_vol_3D(sim_domain *sd, higflow_solver *ns, int dim, Point Center, Point P,Point Delta,real *fracvol);

real direction_hf_x(sim_domain *sdp, higflow_solver *ns, Point center, Point delta);

real direction_hf_y(sim_domain *sdp, higflow_solver *ns, Point center, Point delta);

real direction_hf_z(sim_domain *sdp, higflow_solver *ns, Point center, Point delta);

void x_block(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *x_row, int *aux);

void y_block(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *x_row, int *aux);

void z_block(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *x_row, int *aux);

void higflow_compute_curvature_interfacial_force_normal_multiphase_3D_HF_padrao(higflow_solver *ns);

#endif
