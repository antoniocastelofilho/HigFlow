// PLIC in 3D: the same inversion as the 2D case, but cutting a box with a PLANE, so
// the volume-as-a-function-of-offset relation is piecewise cubic instead of
// piecewise quadratic.  That is why this file is larger than its 2D counterpart and
// why it carries its own root finder rather than sharing one.
//
// PART OF THE VOF FAMILY: a volume fraction per cell says how much of each phase is
// there, and the interface is RECONSTRUCTED from it rather than tracked.  The steps
// are: estimate a normal, place an interface from the normal and the fraction, then
// advect the fraction.  These files are the alternative methods for those steps, not
// a pipeline -- several are alternatives to each other.

#ifndef HIG_FLOW_VOF_PLIC_3D
#define HIG_FLOW_VOF_PLIC_3D

#include "hig-flow-discret.h"


real solver_equation_n_nonzero(real volume, real n1, real dy, real dz);

real H(real x);

int left_right_of_the_plane(Point P,Point Normal,real d);

real trans_center_to_p0(Point Delta,Point Normal,real d_from_center);

real trans_p0_to_center(Point Delta, Point Normal, real d_from_bl);

real trans_center_to_p1(Point Delta,Point Normal,real d_from_center);

real trans_p1_to_center(Point Delta,Point Normal,real d_from_br);

real trans_center_to_p2(Point Delta,Point Normal,real d_from_center);

real trans_p2_to_center(Point Delta,Point Normal,real d_from_br);

real trans_center_to_p3(Point Delta,Point Normal,real d_from_center);

real trans_p3_to_center(Point Delta,Point Normal,real d_from_br);

real distance_from_center_3D(Point Normal,Point Delta,real VOLUME);

real parallel_case_volume(Point Normal,Point Delta,real d_from_center,real tol_n);

real volume_one_n_zero(real dx, real dy, real dz, real n_y, real n_z, real ny, Point Normal, real d);

real volume_3D(real n_x, real n_y, real n_z, real dx, real dy, real dz, real d, real aux_d);

real volume_left_line_origin_center(Point Normal,Point Delta,real d_from_center);

void higflow_compute_distance_multiphase_3D(higflow_solver *ns);

void higflow_compute_volume_fraction_multiphase_3D(higflow_solver *ns);

real minimo(real vol1, real vol2);

real maximo(real nx, real ny, real nz);

#endif
