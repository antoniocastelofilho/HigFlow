// Normal and curvature in 3D by differencing the volume fraction.
//
// The _dominant variants pick which axis to treat as the interface's principal
// direction before differencing -- differencing along an axis the interface is
// nearly parallel to is what makes this estimator fail, and choosing the dominant
// axis avoids that case.
//
// PART OF THE VOF FAMILY: a volume fraction per cell says how much of each phase is
// there, and the interface is RECONSTRUCTED from it rather than tracked.  The steps
// are: estimate a normal, place an interface from the normal and the fraction, then
// advect the fraction.  These files are the alternative methods for those steps, not
// a pipeline -- several are alternatives to each other.

#ifndef HIG_FLOW_VOF_FINITE_DIFFERENCE_NORMAL_CURVATURE_3D
#define HIG_FLOW_VOF_FINITE_DIFFERENCE_NORMAL_CURVATURE_3D

#include "hig-flow-discret.h"

void calculate_HF_curvature_x_dominant(higflow_solver *ns, int clid, real Xl, real Xr, real Yb, real Yu, real XY, real Xul, real Xur, real Xbl, real Xbr, real dx, real dy, real dz, int aux);

void calculate_HF_curvature_y_dominant(higflow_solver *ns, int clid, real Xl, real Xr, real Yb, real Yu, real XY, real Xul, real Xur, real Xbl, real Xbr, real dx, real dy, real dz, int aux);

void calculate_HF_curvature_z_dominant(higflow_solver *ns, int clid, real Xl, real Xr, real Yb, real Yu, real XY, real Xul, real Xur, real Xbl, real Xbr, real dx, real dy, real dz, int aux);

void calculate_HF_normal_x_dominant(higflow_solver *ns, int clid, real Xl, real Xr, real Yb, real Yu,real dx,real dy, real dz,int auxh);

void calculate_HF_normal_y_dominant(higflow_solver *ns, int clid, real Xl, real Xr, real Yb, real Yu,real dx,real dy, real dz,int auxh);

void calculate_HF_normal_z_dominant(higflow_solver *ns, int clid, real Xl, real Xr, real Yb, real Yu,real dx,real dy, real dz,int auxh);

void calculate_exact_normal_x_dominant(higflow_solver *ns, int clid, Point p, int auxh);

void calculate_exact_normal_y_dominant(higflow_solver *ns, int clid, Point p, int auxh);

void calculate_exact_normal_z_dominant(higflow_solver *ns, int clid, Point p, int auxh);

void calculate_progressive_HF_curvature_x_dominant(higflow_solver *ns, int clid, real Xl, real Xr, real Yuu, real Yu, real XY, real Xul, real Xur, real Xuul, real Xuur, real dx, real dy, real dz, int aux);

void calculate_progressive_HF_normal_x_dominant(higflow_solver *ns, int clid, real Xl, real Xr, real XY, real Yu, real Yuu, real dx,real dy, real dz,int auxh);

void calculate_regressive_HF_curvature_x_dominant(higflow_solver *ns, int clid, real Xl, real Xr, real XY, real Xbl, real Xbr, real Xb, real Xbbl, real Xbbr, real Xbb, real dx, real dy, real dz, int aux);

void calculate_regressive_HF_normal_x_dominant(higflow_solver *ns, int clid, real Xl, real Xr, real XY, real Yu, real Yuu, real dx,real dy, real dz,int auxh);

#endif
