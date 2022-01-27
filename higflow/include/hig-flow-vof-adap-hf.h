#ifndef HIG_FLOW_VOF_ADAP_HF
#define HIG_FLOW_VOF_ADAP_HF

#include "hig-flow-discret.h"

void fraction_correction_at_get(real *fracvol);

int get_frac_vol(sim_domain *sd, higflow_solver *ns, int dim, Point Center, Point P,Point Delta,real *fracvol);

void vertical_collumn(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *vertical, int *aux, real *orig);

void horizontal_row(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *horizontal, int *aux, real *orig);

real real_max_vec(real *vec,int n);

real real_min_vec(real *vec,int n);

void real_times_vec(real c,real *vec,int n);

void set_vec(real orig,real *vec,int n);

void vec_minus_vec(real *vec1, real *vec2, real *vec_new, int n);

void vec_plus_vec(real *vec1, real *vec2, real *vec_new, int n);

void set_common_orig_vertical(real *V,int auxv,real *origv,int n,real dy);

void set_common_orig_horizontal(real *H,int auxh,real *origh,int n,real dx);

void calculate_interfacial_force(sim_domain *sdp, higflow_solver *ns, int clid, Point center, Point IF);

void higflow_compute_curvature_interfacial_force_normal_multiphase_2D(higflow_solver *ns);

void higflow_compute_curvature_interfacial_force_normal_multiphase_2D_hf_elvira(higflow_solver *ns);

void higflow_compute_curvature_interfacial_force_normal_multiphase_2D_shirami(higflow_solver *ns);

#endif
