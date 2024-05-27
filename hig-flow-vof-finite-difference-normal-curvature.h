#ifndef HIG_FLOW_VOF_FINITE_DIFFERENCE_NORMAL_CURVATURE
#define HIG_FLOW_VOF_FINITE_DIFFERENCE_NORMAL_CURVATURE

#include "hig-flow-discret.h"

real norm_vec(real *vec, int n);

// Calculate Norm-2
real norm_vec(real *vec, int n);

// Calculate Normal Elvira: Regressive_Horizontal
void elvira_calculate_normal_cell_regressive_1st_order_finite_difference_Horizontal(higflow_solver *ns, int clid, Point Normal, real Hm, real Hb, real dx, real dy, int auxh);

// Calculate Normal Elvira: Central_Horizontal
void elvira_calculate_normal_cell_central_2nd_order_finite_difference_Horizontal(higflow_solver *ns, int clid, Point Normal, real Ht, real Hb, real dx, real dy, int auxh);

// Calculate Normal Elvira: Progressive_Horizontal
void elvira_calculate_normal_cell_progressive_1st_order_finite_difference_Horizontal(higflow_solver *ns, int clid, Point Normal, real Ht, real Hm, real dx, real dy, int auxh);

// Calculate Normal Elvira: Regressive_Vertical
void elvira_calculate_normal_cell_regressive_1st_order_finite_difference_Vertical(higflow_solver *ns, int clid, Point Normal, real Vm, real Vl, real dx, real dy, int auxv);

// Calculate Normal Elvira: Central_Vertical
void elvira_calculate_normal_cell_central_2nd_order_finite_difference_Vertical(higflow_solver *ns, int clid, Point Normal, real Vr, real Vl, real dx, real dy, int auxv);

// Calculate Normal Elvira: Progressive_Vertical
void elvira_calculate_normal_cell_progressive_1st_order_finite_difference_Vertical(higflow_solver *ns, int clid, Point Normal, real Vr, real Vm, real dx, real dy, int auxv);

/////////////////////////////////////////////
// Calculate Normal HF: Regressive_Horizontal
void calculate_normal_cell_regressive_2nd_order_finite_difference_Horizontal(higflow_solver *ns, int clid, real Hk_k_, real Hk_, real H, real dx, real dy, int auxh);

// Calculate Normal HF: Central_Horizontal
void calculate_normal_cell_central_2nd_order_finite_difference_Horizontal(higflow_solver *ns, int clid, real Hk_, real H, real Hk, real dx, real dy, int auxh);

// Calculate Normal HF: Progressive_Horizontal
void calculate_normal_cell_progressive_2nd_order_finite_difference_Horizontal(higflow_solver *ns, int clid, real H, real Hk, real Hkk, real dx, real dy, int auxh);

// Calculate Normal HF: Regressive_Vertical
void calculate_normal_cell_regressive_2nd_order_finite_difference_Vertical(higflow_solver *ns, int clid, real Vk_k_, real Vk_, real V, real dx, real dy, int auxv);

// Calculate Normal HF: Central_Vertical
void calculate_normal_cell_central_2nd_order_finite_difference_Vertical(higflow_solver *ns, int clid, real Vk_, real V, real Vk, real dx, real dy, int auxv);

// Calculate Normal HF: Progressive_Vertical
void calculate_normal_cell_progressive_2nd_order_finite_difference_Vertical(higflow_solver *ns, int clid, real V,real Vk,real Vkk,real dx,real dy,int auxv);

// Calculate Curvature HF: Regressive_Horizontal
void calculate_curvature_cell_regressive_1st_order_finite_difference_Horizontal(higflow_solver *ns, int clid, real Hk_k_,real Hk_,real H,real dx,real dy,int auxh);

// Calculate Curvature HF: Central_Horizontal
void calculate_curvature_cell_central_2nd_order_finite_difference_Horizontal(higflow_solver *ns, int clid, real Hk_,real H,real Hk,real dx,real dy,int auxh);

// Calculate Curvature HF: Progressive_Horizontal
void calculate_curvature_cell_progressive_1st_order_finite_difference_Horizontal(higflow_solver *ns, int clid, real H,real Hk,real Hkk,real dx,real dy,int auxh);

// Calculate Curvature HF: Regressive_Vertical
void calculate_curvature_cell_regressive_1st_order_finite_difference_Vertical(higflow_solver *ns, int clid, real Vk_k_,real Vk_,real V,real dx,real dy,int auxv);

// Calculate Curvature HF: Central_Vertical
void calculate_curvature_cell_central_2nd_order_finite_difference_Vertical(higflow_solver *ns, int clid, real Vk_,real V,real Vk,real dx,real dy,int auxv);

// Calculate Curvature HF: Progressive_Vertical
void calculate_curvature_cell_progressive_1st_order_finite_difference_Vertical(higflow_solver *ns, int clid, real V,real Vk,real Vkk,real dx,real dy,int auxv);

//// Normal Exact
void calculate_exact_normal_x_dominant_2D(higflow_solver *ns, int clid, Point p);

void calculate_exact_normal_y_dominant_2D(higflow_solver *ns, int clid, Point p);

#endif
