// Advecting the volume fraction in 2D, one direction at a time.
//
// DIRECTIONALLY SPLIT ADVECTION DOES NOT CONSERVE THE FRACTION BY CONSTRUCTION --
// a sweep can leave a cell slightly outside [0,1].  That is why the
// fraction_correction_at_* functions exist and why they are part of the method
// rather than a safety net: removing them does not make the scheme stricter, it
// makes the fraction field invalid.
//
// PART OF THE VOF FAMILY: a volume fraction per cell says how much of each phase is
// there, and the interface is RECONSTRUCTED from it rather than tracked.  The steps
// are: estimate a normal, place an interface from the normal and the fraction, then
// advect the fraction.  These files are the alternative methods for those steps, not
// a pipeline -- several are alternatives to each other.

#ifndef HIG_FLOW_VOF_ADVECTION_2D
#define HIG_FLOW_VOF_ADVECTION_2D

#include "hig-flow-discret.h"
#include "hig-flow-vof-adap-hf.h"
#include "hig-flow-vof-plic.h"

void higflow_plic_copy_fractionaux_to_fraction(higflow_solver *ns);

void fraction_correction_at_set(real *frac);


void higflow_plic_advection_volume_fraction_x_direction_2D(higflow_solver *ns, int dim);

void higflow_plic_advection_volume_fraction_x_direction_imp_2D(higflow_solver *ns, int dim);

void higflow_plic_advection_volume_fraction_y_direction_2D(higflow_solver *ns, int dim);

void higflow_plic_advection_volume_fraction_y_direction_imp_2D(higflow_solver *ns, int dim);

#endif
