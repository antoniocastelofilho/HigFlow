#ifndef HIG_FLOW_VOF_ADVECTION_3D
#define HIG_FLOW_VOF_ADVECTION_3D

#include "hig-flow-discret.h"
#include "hig-flow-vof-HF-3D.h"
#include "hig-flow-vof-plic-3D.h"

void normal_correction_at_get(Point Normal);

void fraction_correction_at_set(real *frac);

void higflow_plic_advection_volume_fraction_x_direction(higflow_solver *ns, int dim);

void higflow_plic_advection_volume_fraction_x_direction_imp(higflow_solver *ns, int dim);

void higflow_plic_advection_volume_fraction_y_direction(higflow_solver *ns, int dim);

void higflow_plic_advection_volume_fraction_y_direction_imp(higflow_solver *ns, int dim);

void higflow_plic_advection_volume_fraction_z_direction(higflow_solver *ns, int dim);

void higflow_plic_advection_volume_fraction_z_direction_imp(higflow_solver *ns, int dim);

#endif
