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
