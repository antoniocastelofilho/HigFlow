// Advecting the volume fraction in 3D, one direction at a time.  Same split scheme
// and same need for fraction correction as the 2D case.
//
// The _imp variants are the implicit form of a directional sweep, for the steps
// where an explicit sweep would violate the CFL on that axis alone.
//
// PART OF THE VOF FAMILY: a volume fraction per cell says how much of each phase is
// there, and the interface is RECONSTRUCTED from it rather than tracked.  The steps
// are: estimate a normal, place an interface from the normal and the fraction, then
// advect the fraction.  These files are the alternative methods for those steps, not
// a pipeline -- several are alternatives to each other.

#ifndef HIG_FLOW_VOF_ADVECTION_3D
#define HIG_FLOW_VOF_ADVECTION_3D

#include "hig-flow-discret.h"
#include "hig-flow-vof-HF-3D.h"
#include "hig-flow-vof-plic-3D.h"


void fraction_correction_at_set(real *frac);

void higflow_plic_advection_volume_fraction_x_direction(higflow_solver *ns, int dim);

void higflow_plic_advection_volume_fraction_x_direction_imp(higflow_solver *ns, int dim);

void higflow_plic_advection_volume_fraction_y_direction(higflow_solver *ns, int dim);

void higflow_plic_advection_volume_fraction_y_direction_imp(higflow_solver *ns, int dim);

void higflow_plic_advection_volume_fraction_z_direction(higflow_solver *ns, int dim);

void higflow_plic_advection_volume_fraction_z_direction_imp(higflow_solver *ns, int dim);

#endif
