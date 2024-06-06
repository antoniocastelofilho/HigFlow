#ifndef HIG_FLOW_VOF_PLIC
#define HIG_FLOW_VOF_PLIC

#include "hig-flow-discret.h"

void solve_equation(int Rx, int Ry, Point Normal, Point Delta, real area, Point Value);

real parallel_left_line_origin_center_distance(Point Normal,Point Delta,real area,real tol_n);

real distance_from_center(Point Normal,Point Delta,real AREA);

real parallel_left_line_origin_center_area(Point Normal,Point Delta,real d_from_center,real tol_n);

real area_left_line_origin_center(Point Normal,Point Delta,real d_from_center);

int left_right(Point P,Point Normal,real d);

real trans_center_2_bl(Point Delta,Point Normal,real d_from_center);

real trans_bl_2_center(Point Delta,Point Normal,real d_from_bl);

real trans_center_2_br(Point Delta,Point Normal,real d_from_center);

real trans_br_2_center(Point Delta,Point Normal,real d_from_br);

real R(real x);

int sign(real value);

void higflow_compute_distance_multiphase_2D(higflow_solver *ns);

void higflow_compute_area_fraction_multiphase_2D(higflow_solver *ns);

#endif
