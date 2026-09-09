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
