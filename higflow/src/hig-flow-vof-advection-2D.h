#ifndef HIG_FLOW_VOF_ADVECTION_2D
#define HIG_FLOW_VOF_ADVECTION_2D

#include "hig-flow-discret.h"
#include "hig-flow-vof-adap-hf.h"
#include "hig-flow-vof-plic.h"

void higflow_plic_copy_fractionaux_to_fraction(higflow_solver *ns);

void fraction_correction_at_set(real *frac);

void normal_correction_at_get(Point Normal);

void higflow_plic_advection_volume_fraction_x_direction_2D(higflow_solver *ns, int dim);

void higflow_plic_advection_volume_fraction_x_direction_imp_2D(higflow_solver *ns, int dim);

void higflow_plic_advection_volume_fraction_y_direction_2D(higflow_solver *ns, int dim);

void higflow_plic_advection_volume_fraction_y_direction_imp_2D(higflow_solver *ns, int dim);

#endif
