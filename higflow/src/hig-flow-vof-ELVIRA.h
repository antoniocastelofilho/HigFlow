#ifndef HIG_FLOW_VOF_ELVIRA
#define HIG_FLOW_VOF_ELVIRA

#include "hig-flow-discret.h"

void fraction_correction_at_get(real *fracvol);

void ELVIRA_vertical_collumn(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *vertical, int *aux);

void ELVIRA_horizontal_row(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *horizontal, int *aux);

void higflow_compute_curvature_interfacial_force_normal_multiphase_2D_EL(higflow_solver *ns);


#endif
