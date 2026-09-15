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
