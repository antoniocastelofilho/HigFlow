#include "hig-flow-vof-elvira.h"

real area_correction_at_get(Point delta, real area){
   if (FLT_EQ(area, delta[0]*delta[1])) {
      area = delta[0]*delta[1];
   } else if (FLT_EQ(area, 0.0)) {
      area = 0.0;
   }
   return area;
}

void _elvira_vertical_collumn(sim_domain *sdm, higflow_solver *ns, Point center, Point p, Point delta, real *vertical, int *aux){
   int status;
   real fracvol, fracvol_aux;
   Point pp, ppt, ppb;
   *aux  = 0;
   // Vertical Up
   status = get_frac_vol(sdm, ns, 1, center, p, delta, &fracvol);
   if (status != 1) {
      if (status == -1) {
         //printf("Up - Vertical cells with different sizes \n");
      } else {
         //printf("Up - Vertical cell out of domain \n");
      }
      return;
   }
   pp[0] = p[0];
   *vertical = fracvol;
   int i = 0;
   do { i++;
      pp[1] = center[1] + i * delta[1];
      status = get_frac_vol(sdm, ns, 1, center, pp, delta, &fracvol);
      if (status == 1) {
         *vertical += fracvol;
      } else if (status == 0) {
         //printf("Up - Vertical cell out of domain \n");
         return;
      } else {
         //printf("Up - Vertical cells with different sizes \n");
         return;
      }
   } while (i < 1 && status == 1);
   fracvol_aux = fracvol;
   ppt[0] = pp[0]; 
   ppt[1] = pp[1] + 0.5*delta[1]; 
   real fracvol_top = fracvol;             
   // Vertical Down
   i = 0;
   do { i++;
      pp[1] = p[1] - i * delta[1];
      status = get_frac_vol(sdm, ns, 1, center, pp, delta, &fracvol);
      if (status == 1) {
         *vertical += fracvol;
      } else if (status == 0) {
         //printf("Down - Vertical cell out of domain \n");
         return;
      } else {
         //printf("Down - Vertical cells with different sizes \n");
         return;
      }
   } while (i < 1 && status == 1);
   if (fracvol == fracvol_aux) {
      //printf("The phases are not different \n");
      //return;
   }
   if (fracvol == 1.0 && fracvol_top == 0.0) {
      *aux = -1;
   } else {
      *aux = 1;
   }
}

void elvira_vertical_collumn(sim_domain *sdm, higflow_solver *ns, Point center, Point p, Point delta, real *vertical, int *aux){
   int status;
   real fracvol, fracvol_aux;
   Point pp, ppt, ppb;
   *aux  = 0;
   // Vertical Up
   status = get_frac_vol(sdm, ns, 1, center, p, delta, &fracvol);
   if (status != 1) {
      if (status == -1) {
         //printf("Up - Vertical cells with different sizes \n");
      } else {
         //printf("Up - Vertical cell out of domain \n");
      }
      return;
   }
   pp[0] = p[0];
   *vertical = fracvol;
   int i = 0;
   do { i++;
      pp[1] = center[1] + i * delta[1];
      status = get_frac_vol(sdm, ns, 1, center, pp, delta, &fracvol);
      if (status == 1) {
         *vertical += fracvol;
      } else if (status == 0) {
         //printf("Up - Vertical cell out of domain \n");
         return;
      } else {
         //printf("Up - Vertical cells with different sizes \n");
         return;
      }
   } while (i < 1 && status == 1);
   fracvol_aux = fracvol;
   ppt[0] = pp[0]; 
   ppt[1] = pp[1] + 0.5*delta[1]; 
   real fracvol_top = fracvol;             
   // Vertical Down
   i = 0;
   do { i++;
      pp[1] = p[1] - i * delta[1];
      status = get_frac_vol(sdm, ns, 1, center, pp, delta, &fracvol);
      if (status == 1) {
         *vertical += fracvol;
      } else if (status == 0) {
         //printf("Down - Vertical cell out of domain \n");
         return;
      } else {
         //printf("Down - Vertical cells with different sizes \n");
         return;
      }
   } while (i < 1 && status == 1);
   if (fracvol == fracvol_aux) {
      //printf("The phases are not different \n");
      //return;
   }
   if (fracvol == 1.0 && fracvol_top == 0.0) {
      *aux = -1;
   } else {
      *aux = 1;
   }
}

void _elvira_horizontal_row(sim_domain *sdm, higflow_solver *ns, Point center, Point p, Point delta, real *horizontal, int *aux){
   int status;
   real fracvol, fracvol_aux;
   Point pp, ppr, ppl;
   *aux=0;
   // Horizontal Right
   status = get_frac_vol(sdm, ns, 0, center, p, delta, &fracvol);
   if (status != 1) {
      if (status == -1) {
         //printf("Right - Horizontal cells with different sizes \n");
      } else {
         //printf("Right - Horizontal cell out of domain \n");
      }
      return;
   }
   pp[1] = p[1];
   *horizontal = fracvol;
   int i = 0;
   do { i++;
      pp[0] = center[0] + i*delta[0];
      status = get_frac_vol(sdm, ns, 0, center, pp, delta, &fracvol);
      if (status == 1) {
         *horizontal += fracvol;
      } else if (status == 0) {
         //printf("Right - Horizontal cell out of domain \n");
         return;
      } else {
         //printf("Right - Horizontal cells with different sizes \n");
         return;
      }
   } while (i < 1 && status == 1);
   fracvol_aux = fracvol;
   ppr[0]=pp[0] + 0.5*delta[0];
   ppr[1]=pp[1];
   real fracvol_right = fracvol;
   // Horizontal Left
   i = 0;
   // fracvol_aux = 0;
   do { i++;
      pp[0] = p[0] - i*delta[0];
      status = get_frac_vol(sdm, ns, 0, center, pp, delta, &fracvol);
      if (status == 1) {
         *horizontal += fracvol;
      } else if (status == 0) {
         //printf("Left - Horizontal cell out of domain \n");
         return;
      } else {
         //printf("Left - Horizontal cells with different sizes \n");
         return;
      }
   } while (i < 1 && status == 1);
   if (fracvol == fracvol_aux){
      //printf("The phases are not different:\n");
      //return;
   }
   if (fracvol == 1.0 && fracvol_right == 0.0) {
      *aux = -1;
   } else {
      *aux = 1;
   }
}

void elvira_horizontal_row(sim_domain *sdm, higflow_solver *ns, Point center, Point p, Point delta, real *horizontal, int *aux){
   int status;
   real fracvol, fracvol_aux;
   Point pp, ppr, ppl;
   *aux=0;
   // Horizontal Right
   status = get_frac_vol(sdm, ns, 0, center, p, delta, &fracvol);
   if (status != 1) {
      if (status == -1) {
         //printf("Right - Horizontal cells with different sizes \n");
      } else {
         //printf("Right - Horizontal cell out of domain \n");
      }
      return;
   }
   pp[1] = p[1];
   *horizontal = fracvol;
   int i = 0;
   do { i++;
      pp[0] = center[0] + i*delta[0];
      status = get_frac_vol(sdm, ns, 0, center, pp, delta, &fracvol);
      if (status == 1) {
         *horizontal += fracvol;
      } else if (status == 0) {
         //printf("Right - Horizontal cell out of domain \n");
         return;
      } else {
         //printf("Right - Horizontal cells with different sizes \n");
         return;
      }
   } while (i < 1 && status == 1);
   fracvol_aux = fracvol;
   ppr[0]=pp[0] + 0.5*delta[0];
   ppr[1]=pp[1];
   real fracvol_right = fracvol;
   // Horizontal Left
   i = 0;
   // fracvol_aux = 0;
   do { i++;
      pp[0] = p[0] - i*delta[0];
      status = get_frac_vol(sdm, ns, 0, center, pp, delta, &fracvol);
      if (status == 1) {
         *horizontal += fracvol;
      } else if (status == 0) {
         //printf("Left - Horizontal cell out of domain \n");
         return;
      } else {
         //printf("Left - Horizontal cells with different sizes \n");
         return;
      }
   } while (i < 1 && status == 1);
   if (fracvol == fracvol_aux){
      //printf("The phases are not different:\n");
      //return;
   }
   if (fracvol == 1.0 && fracvol_right == 0.0) {
      *aux = -1;
   } else {
      *aux = 1;
   }
}

void elvira_vertical_collumn_adap(sim_domain *sdm, higflow_solver *ns, Point center, Point p, Point delta, real *vertical, int *aux, real *orig){
   int status;
   real fracvol, fracvol_aux;
   Point pp, ppt, ppb;
   *aux  = 0;
   *orig = 0.0;
   // Vertical Up
   status = get_frac_vol(sdm, ns, 1, center, p, delta, &fracvol);
   if (status != 1) {
      if (status == -1) {
         //printf("Up - Vertical cells with different sizes \n");
      } else {
         //printf("Up - Vertical cell out of domain \n");
      }
      return;
   }
   pp[0] = p[0];
   *vertical = fracvol;
   int i = 0;
   do { i++;
      pp[1] = center[1] + i * delta[1];
      status = get_frac_vol(sdm, ns, 1, center, pp, delta, &fracvol);
      if (status == 1) {
         *vertical += fracvol;
      } else if (status == 0) {
         //printf("Up - Vertical cell out of domain \n");
         return;
      } else {
         //printf("Up - Vertical cells with different sizes \n");
         return;
      }
   } while (fracvol > 0.0 && fracvol < 1.0 && status == 1);
   fracvol_aux = fracvol;
   ppt[0]=pp[0];ppt[1]=pp[1]+0.5*delta[1]; 
   real fracvol_top = fracvol;             
   // Vertical Down
   i = 0;
   do { i++;
      pp[1] = p[1] - i * delta[1];
      status = get_frac_vol(sdm, ns, 1, center, pp, delta, &fracvol);
      if (status == 1) {
         *vertical += fracvol;
      } else if (status == 0) {
         //printf("Down - Vertical cell out of domain \n");
         return;
      } else {
         //printf("Down - Vertical cells with different sizes \n");
         return;
      }
   } while (fracvol > 0.0 && fracvol < 1.0 && FLT_NE(fracvol, fracvol_aux) && status == 1);

   if (fracvol == fracvol_aux) {
      //printf("The phases are not different \n");
      //return;
   }
   ppb[0]=pp[0];ppb[1]=pp[1]-0.5*delta[1];

   if (fracvol == 1.0 && fracvol_top == 0.0) {
      *aux = -1;
      *orig=ppb[1];
   } else {
      *aux = 1;
      *orig=ppt[1];
   }
}

void elvira_horizontal_row_adap(sim_domain *sdm, higflow_solver *ns, Point center, Point p, Point delta, real *horizontal, int *aux, real *orig){
   int status;
   real fracvol, fracvol_aux;
   Point pp, ppr, ppl;
   *aux=0;
   *orig=0.0;
   // Horizontal Right
   status = get_frac_vol(sdm, ns, 0, center, p, delta, &fracvol);
   if (status != 1) {
      if (status == -1) {
         //printf("Right - Horizontal cells with different sizes \n");
      } else {
         //printf("Right - Horizontal cell out of domain \n");
      }
      return;
   }
   pp[1] = p[1];
   *horizontal=fracvol;
   int i = 0;
   do { i++;
      pp[0] = center[0] + i*delta[0];
      status=get_frac_vol(sdm, ns, 0, center, pp, delta, &fracvol);
      if (status == 1) {
         *horizontal += fracvol;
      } else if (status == 0) {
         //printf("Right - Horizontal cell out of domain \n");
         return;
      } else {
         //printf("Right - Horizontal cells with different sizes \n");
         return;
      }
   } while (fracvol > 0.0 && fracvol < 1.0 && status == 1);
   fracvol_aux = fracvol;
   ppr[0]=pp[0]+0.5*delta[0];
   ppr[1]=pp[1];
   real fracvol_right = fracvol;
   // Horizontal Left
   i = 0;
   //fracvol_aux=0;
   do { i++;
      pp[0] = p[0] - i*delta[0];
      status = get_frac_vol(sdm, ns, 0, center, pp, delta, &fracvol);
      if (status == 1) {
         *horizontal += fracvol;
      } else if (status == 0) {
         //printf("Left - Horizontal cell out of domain \n");
         return;
      } else {
         //printf("Left - Horizontal cells with different sizes \n");
         return;
      }
   } while (fracvol > 0.0 && fracvol < 1.0 && FLT_NE(fracvol, fracvol_aux) && status == 1);

   if (fracvol == fracvol_aux){
      //printf("The phases are not different:\n");
      //return;
   }
   ppl[0]=pp[0]-0.5*delta[0];
   ppl[1]=pp[1];

   if (fracvol == 1.0 && fracvol_right == 0.0) {
      *aux = -1;
      *orig=ppl[0];
   } else {
      *aux = 1;
      *orig=ppr[0];
   }
}

real real_max_vec_elvira(real *vec,int n){
   real max = vec[1];
   for(int i=1;i<n;i++){
      if(vec[i]>max){
         max=vec[i];
      }
   }
   return max;
}

real real_min_vec_elvira(real *vec,int n){
   real min = vec[1];
   for(int i=1;i<n;i++){
      if(vec[i]<min){
         min=vec[i];
      }
   }
   return min;
}

void real_times_vec_elvira(real c,real *vec,int n){
   for(int i=0;i<n;i++){
      vec[i]=c*vec[i];
   }
}

void set_vec_elvira(real orig,real *vec,int n){
   for(int i=0;i<n;i++){
      vec[i] = orig;
   }
}

void vec_minus_vec_elvira(real *vec1, real *vec2, real *vec_new, int n){
   for(int i=0;i<n;i++){
      vec_new[i] = vec1[i] - vec2[i];
   }
}
void vec_plus_vec_elvira(real *vec1, real *vec2, real *vec_new, int n){
   for(int i=0;i<n;i++){
      vec_new[i] = vec1[i] + vec2[i];
   }
}

void set_common_orig_vertical_elvira(real *V,int auxv,real *origv,int n,real dy){
   real orig;
   real cons[n],delta[n];
   if(auxv>0){
      orig=real_max_vec_elvira(origv,n);
      set_vec_elvira(orig,cons,n);
      vec_minus_vec_elvira(cons,origv,delta,n);
      real_times_vec_elvira(1/dy,delta,n);
      vec_plus_vec_elvira(V,delta,V,n);
   } else {
      orig=real_min_vec_elvira(origv,n);
      set_vec_elvira(-orig,cons,n);
      vec_plus_vec_elvira(cons,origv,delta,n);
      real_times_vec_elvira(1/dy,delta,n);
      vec_plus_vec_elvira(V,delta,V,n);
   }
}

void set_common_orig_horizontal_elvira(real *H,int auxh,real *origh,int n,real dx){
   real orig;
   real cons[n],delta[n];
   if(auxh>0){
      orig=real_max_vec_elvira(origh,n);
      set_vec_elvira(orig,cons,n);
      vec_minus_vec_elvira(cons,origh,delta,n);
      real_times_vec_elvira(1/dx,delta,n);
      vec_plus_vec_elvira(H,delta,H,n);
   } else {
      orig=real_min_vec_elvira(origh,n);
      set_vec_elvira(-orig,cons,n);
      vec_plus_vec_elvira(cons,origh,delta,n);
      real_times_vec_elvira(1/dx,delta,n);
      vec_plus_vec_elvira(H,delta,H,n);
   }
}

real elvira_interface_error(real area[6], real areaE[6][9], int s) {
  real error = 0.0;
  for (int j=0; j<9; j++) {
     error += pow((area[j] - areaE[s][j]),2);
  }
  return sqrt(error);
}

void higflow_compute_normal_multiphase_2D_elvira(higflow_solver *ns, sim_domain *sdm, mp_mapper *mp, higcit_celliterator *it, hig_cell *c,int clid, Point center, Point delta, Point p) {
   // Case bi-dimensional
   real IF[DIM];
   Point pp;
   real area[9]; real fracvol[9];
   pp[0]      = p[0];
   pp[1]      = p[1];
   fracvol[4] = compute_value_at_point(sdm, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
   area[4]    = fracvol[4]*delta[0]*delta[1];
   // y
   // ^
   // |-dx[0]-|-dx[1]-|-dx[2]-|
   // |       |       |       |
   // +-------+-------+-------+--+
   // |       |       |       |  |
   // | a[6]  | a[7]  | a[8]  |  dY[2]
   // |       |       |       |  |
   // +-------+-------+-------+--+
   // |       |       |       |  |
   // | a[3]  | a[4]  | a[5]  |  dY[1]
   // |       |       |       |  |
   // +-------+-------+-------+--+
   // |       |       |       |  |
   // | a[0]  | a[1]  | a[2]  |  dY[0]
   // |       |       |       |  |
   // +-------+-------+-------+--+--> 
   pp[0]      = p[0] - delta[0]; 
   pp[1]      = p[1] - delta[1];
   fracvol[0] = compute_value_at_point(sdm, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
   area[0]    = fracvol[0]*delta[0]*delta[1];
   
   pp[0]      = p[0]; 
   pp[1]      = p[1] - delta[1];
   fracvol[1] = compute_value_at_point(sdm, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
   area[1]    = fracvol[1]*delta[0]*delta[1];
   
   pp[0]      = p[0] + delta[0]; 
   pp[1]      = p[1] - delta[1];
   fracvol[2] = compute_value_at_point(sdm, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
   area[2]    = fracvol[2]*delta[0]*delta[1];
   
   pp[0]      = p[0] - delta[0]; 
   pp[1]      = p[1];
   fracvol[3] = compute_value_at_point(sdm, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
   area[3]    = fracvol[3]*delta[0]*delta[1];
   
   pp[0]      = p[0] + delta[0]; 
   pp[1]      = p[1];
   fracvol[5] = compute_value_at_point(sdm, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
   area[5]    = fracvol[5]*delta[0]*delta[1];
   
   pp[0]      = p[0] - delta[0];
   pp[1]      = p[1] + delta[1];
   fracvol[6] = compute_value_at_point(sdm, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
   area[6]    = fracvol[6]*delta[0]*delta[1];
   
   pp[0]      = p[0];
   pp[1]      = p[1] + delta[1];
   fracvol[7] = compute_value_at_point(sdm, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
   area[7]    = fracvol[7]*delta[0]*delta[1];
   
   pp[0]      = p[0] + delta[0];
   pp[1]      = p[1] + delta[1];
   fracvol[8] = compute_value_at_point(sdm, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
   area[8]    = fracvol[8]*delta[0]*delta[1];
   
   //Height functions
   real hm, hb, ht, vm, vl, vr;
   int aux_mh = 0, aux_b = 0, aux_t = 0;
   int aux_mv = 0, aux_l = 0, aux_r = 0;
   //// Horizontal ////
   pp[0]      = p[0];
   pp[1]      = p[1];
   // Middle
   elvira_horizontal_row(sdm, ns, center, pp, delta, &hm, &aux_mh);
   ////===========================================================
   // Top
   pp[1]      = p[1]+delta[1];
   elvira_horizontal_row(sdm, ns, center, pp, delta, &ht, &aux_t);
   ////===========================================================
   //Bottom
   pp[1]      = p[1]-delta[1];
   elvira_horizontal_row(sdm, ns, center, pp, delta, &hb, &aux_b);
   ////===========================================================
   real H[3];
   H[0]       = hb;
   H[1]       = hm;
   H[2]       = ht;
   
   //// Vertical ////
   pp[0]      = p[0]; 
   pp[1]      = p[1];
   // Middle
   elvira_vertical_collumn(sdm, ns, center, pp, delta, &vm, &aux_mv);
   ////===========================================================
   //Right
   pp[0]      = p[0]+delta[0];
   elvira_vertical_collumn(sdm, ns, center, pp, delta, &vr, &aux_r);
   ////===========================================================
   //Left
   pp[0]      = p[0]-delta[0];
   elvira_vertical_collumn(sdm, ns, center, pp, delta, &vl, &aux_l);
   ////===========================================================
   real V[3];
   V[0]       = vl; 
   V[1]       = vm; 
   V[2]       = vr;
   
   Point Normal; 
   real areaE[6][9], distanceE[6][9], nx[6], ny[6], d[9], aux_x[9], aux_y[9];

   aux_x[0] = -1; aux_x[1] = 0; aux_x[2] = 1;
   aux_x[3] = -1; aux_x[4] = 0; aux_x[5] = 1;
   aux_x[6] = -1; aux_x[7] = 0; aux_x[8] = 1;

   aux_y[0] = -1; aux_y[1] = -1; aux_y[2] = -1;
   aux_y[3] = 0;  aux_y[4] = 0;  aux_y[5] = 0;
   aux_y[6] = 1;  aux_y[7] = 1;  aux_y[8] = 1;

   //=============================================================
   //// Normals, distances and area PLIC
   //=============================================================
   ////horizontal-regressive finite difference 
   elvira_calculate_normal_cell_regressive_1st_order_finite_difference_Horizontal(ns, clid, Normal, H[1], H[0], delta[0], delta[1], aux_mh);
   nx[0] = Normal[0]; ny[0] = Normal[1];
   distanceE[0][4] = distance_from_center(Normal,delta,area[4]);
   for (int j=0; j<9; j++) {
      //distanceE[0][j] = distance_from_center(Normal,delta,area[j]);
      d[j]        = distanceE[0][4] - (aux_x[j]*Normal[0]*delta[0] + aux_y[j]*Normal[1]*delta[1]);
      areaE[0][j] = area_left_line_origin_center(Normal, delta, d[j]);
      areaE[0][j] = area_correction_at_get(delta, areaE[0][j]);
   }
   //=============================================================
   ////horizontal-central finite difference 
   elvira_calculate_normal_cell_central_2nd_order_finite_difference_Horizontal(ns, clid, Normal, H[2], H[0], delta[0], delta[1], aux_mh);
   nx[1] = Normal[0]; ny[1] = Normal[1];
   distanceE[1][4] = distance_from_center(Normal,delta,area[4]);
   for (int j=0; j<9; j++) {
      //distanceE[1][j] = distance_from_center(Normal,delta,area[j]);
      d[j]        = distanceE[1][4] - (aux_x[j]*Normal[0]*delta[0] + aux_y[j]*Normal[1]*delta[1]);
      areaE[1][j] = area_left_line_origin_center(Normal, delta, d[j]);
      areaE[1][j] = area_correction_at_get(delta, areaE[1][j]);
   }
   //=============================================================
   ////horizontal-progressive finite difference 
   elvira_calculate_normal_cell_progressive_1st_order_finite_difference_Horizontal(ns, clid, Normal, H[2], H[1], delta[0], delta[1], aux_mh);
   nx[2] = Normal[0]; ny[2] = Normal[1];
   distanceE[2][4] = distance_from_center(Normal,delta,area[4]);
   for (int j=0; j<9; j++) {
      //distanceE[2][j] = distance_from_center(Normal,delta,area[j]);
      d[j]        = distanceE[2][4] - (aux_x[j]*Normal[0]*delta[0] + aux_y[j]*Normal[1]*delta[1]);
      areaE[2][j] = area_left_line_origin_center(Normal, delta, d[j]);
      areaE[2][j] = area_correction_at_get(delta, areaE[2][j]);
   }
   //=============================================================
   ////vertical-progressive finite difference 
   elvira_calculate_normal_cell_progressive_1st_order_finite_difference_Vertical(ns, clid, Normal, V[2], V[1], delta[0], delta[1], aux_mv);
   nx[3] = Normal[0]; ny[3] = Normal[1];
   distanceE[3][4] = distance_from_center(Normal,delta,area[4]);
   for (int j=0; j<9; j++) {
      //distanceE[3][j] = distance_from_center(Normal,delta,area[j]);
      d[j]        = distanceE[3][4] - (aux_x[j]*Normal[0]*delta[0] + aux_y[j]*Normal[1]*delta[1]);
      areaE[3][j] = area_left_line_origin_center(Normal, delta, d[j]);
      areaE[3][j] = area_correction_at_get(delta, areaE[3][j]);
   }
   ////===========================================================
   ////vertical-central finite difference 
   elvira_calculate_normal_cell_central_2nd_order_finite_difference_Vertical(ns, clid, Normal, V[2], V[0], delta[0], delta[1], aux_mv);
   nx[4] = Normal[0]; ny[4] = Normal[1];
   distanceE[4][4] = distance_from_center(Normal,delta,area[4]);
   for (int j=0; j<9; j++) {
      //distanceE[4][j] = distance_from_center(Normal,delta,area[j]);
      d[j]        = distanceE[4][4] - (aux_x[j]*Normal[0]*delta[0] + aux_y[j]*Normal[1]*delta[1]);
      areaE[4][j] = area_left_line_origin_center(Normal, delta, d[j]);
      areaE[4][j] = area_correction_at_get(delta, areaE[4][j]);
   }
   ////===========================================================
   ////vertical-regressive finite difference 
   elvira_calculate_normal_cell_regressive_1st_order_finite_difference_Vertical(ns, clid, Normal, V[1], V[0], delta[0], delta[1], aux_mv);
   nx[5] = Normal[0]; ny[5] = Normal[1];
   distanceE[5][4] = distance_from_center(Normal,delta,area[4]);
   for (int j=0; j<9; j++) {
      //distanceE[5][j] = distance_from_center(Normal,delta,area[j]);
      d[j]        = distanceE[5][4] - (aux_x[j]*Normal[0]*delta[0] + aux_y[j]*Normal[1]*delta[1]);
      areaE[5][j] = area_left_line_origin_center(Normal, delta, d[j]);
      areaE[5][j] = area_correction_at_get(delta, areaE[5][j]);
   }
   //============================================================
   real err_min = 1.0e32;
   real errors;
   int posi = 0;
   for (int s=0; s<6; s++){
      errors = elvira_interface_error(area, areaE, s);
      if (errors <= err_min) {
         err_min = errors;
         Normal[0] = nx[s];
         Normal[1] = ny[s];
         posi = s;
      }
      //printf("s = %d \t nx = %.16lf \t ny = %.16lf \t errors = %.16lf\n", s, nx[s], ny[s],errors);
    }
    //printf("#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^#^\n");
    //printf("x = %.16lf \t y = %.16lf\n", p[0], p[1]);
    //printf("position = %d \t nx = %.16lf \t ny = %.16lf\n", posi, nx[posi], ny[posi]);
    //getchar();
    for (int i=0; i<DIM; i++){
      dp_set_value(ns->ed.mult.dpnormal[i], clid, Normal[i]);
    }
}


void higflow_compute_normal_multiphase_2D_elvira_adap(higflow_solver *ns, sim_domain *sdm, mp_mapper *mp, higcit_celliterator *it, hig_cell *c, int clid, Point center, Point delta, Point p) {
    // Case bi-dimensional
    real IF[DIM];
    Point pp;
    real area[9]; real fracvol[9];
    pp[0]      = p[0];
    pp[1]      = p[1];
    fracvol[4] = compute_value_at_point(sdm, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
    area[4]    = fracvol[4]*delta[0]*delta[1];
    // y
    // ^
    // |-dx[0]-|-dx[1]-|-dx[2]-|
    // |       |       |       |
    // +-------+-------+-------+--+
    // |       |       |       |  |
    // | a[6]  | a[7]  | a[8]  |  dY[2]
    // |       |       |       |  |
    // +-------+-------+-------+--+
    // |       |       |       |  |
    // | a[3]  | a[4]  | a[5]  |  dY[1]
    // |       |       |       |  |
    // +-------+-------+-------+--+
    // |       |       |       |  |
    // | a[0]  | a[1]  | a[2]  |  dY[0]
    // |       |       |       |  |
    // +-------+-------+-------+--+--> 
    pp[0]      = p[0] - delta[0]; 
    pp[1]      = p[1] - delta[1];
    fracvol[0] = compute_value_at_point(sdm, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
    area[0]    = fracvol[0]*delta[0]*delta[1];
    
    pp[0]      = p[0]; 
    pp[1]      = p[1] - delta[1];
    fracvol[1] = compute_value_at_point(sdm, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
    area[1]    = fracvol[1]*delta[0]*delta[1];
    
    pp[0]      = p[0] + delta[0]; 
    pp[1]      = p[1] - delta[1];
    fracvol[2] = compute_value_at_point(sdm, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
    area[2]    = fracvol[2]*delta[0]*delta[1];
    
    pp[0]      = p[0] - delta[0]; 
    pp[1]      = p[1];
    fracvol[3] = compute_value_at_point(sdm, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
    area[3]    = fracvol[3]*delta[0]*delta[1];
    
    pp[0]      = p[0] + delta[0]; 
    pp[1]      = p[1];
    fracvol[5] = compute_value_at_point(sdm, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
    area[5]    = fracvol[5]*delta[0]*delta[1];
    
    pp[0]      = p[0] - delta[0];
    pp[1]      = p[1] + delta[1];
    fracvol[6] = compute_value_at_point(sdm, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
    area[6]    = fracvol[6]*delta[0]*delta[1];
    
    pp[0]      = p[0];
    pp[1]      = p[1] + delta[1];
    fracvol[7] = compute_value_at_point(sdm, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
    area[7]    = fracvol[7]*delta[0]*delta[1];
    
    pp[0]      = p[0] + delta[0];
    pp[1]      = p[1] + delta[1];
    fracvol[8] = compute_value_at_point(sdm, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
    area[8]    = fracvol[8]*delta[0]*delta[1];
    
    //Height functions
    real hm, hb, ht, vm, vl, vr;
    int aux_mh   = 0, aux_b  = 0, aux_t  = 0;
    int aux_mv   = 0, aux_l  = 0, aux_r  = 0;
    real orig_mh = 0, orig_b = 0, orig_t = 0;
    real orig_mv = 0, orig_l = 0, orig_r = 0;

    //// Horizontal ////
    pp[0] = p[0];
    pp[1] = p[1];
    // Middle
    //printf("Antes  HM: aux_mh = %d\n",aux_mh);
    elvira_horizontal_row_adap(sdm, ns, center, pp, delta, &hm, &aux_mh, &orig_mh);
    //printf("Depois HM: aux_mh = %d\n",aux_mh);
    ////===========================================================
    // Top
    pp[1] = p[1]+delta[1];
    //printf("Antes  HT: aux_t = %d\n",aux_t);
    elvira_horizontal_row_adap(sdm, ns, center, pp, delta, &ht, &aux_t, &orig_t);
    //printf("Depois HT: aux_t = %d\n",aux_t);
    ////===========================================================
    //Bottom
    pp[1] = p[1]-delta[1];
    //printf("Antes  HB: aux_b = %d\n",aux_b);
    elvira_horizontal_row_adap(sdm, ns, center, pp, delta, &hb, &aux_b, &orig_b);
    //printf("Depois HB: aux_b = %d\n",aux_b);
    ////===========================================================
    real H[3], origh[3];
    H[0] = hb;
    H[1] = hm;
    H[2] = ht;
    
    origh[0] = orig_b;
    origh[1] = orig_mh;
    origh[2] = orig_t;

    set_common_orig_horizontal_elvira(H, aux_mh, origh, 3, delta[0]);

    //// Vertical ////
    pp[0] = p[0]; 
    pp[1] = p[1];
    // Middle
    elvira_vertical_collumn_adap(sdm, ns, center, pp, delta, &vm, &aux_mv, &orig_mv);
    ////===========================================================
    //Right
    pp[0]    = p[0]+delta[0];
    elvira_vertical_collumn_adap(sdm, ns, center, pp, delta, &vr, &aux_r, &orig_r);
    ////===========================================================
    //Left
    pp[0]    = p[0]-delta[0];
    elvira_vertical_collumn_adap(sdm, ns, center, pp, delta, &vl, &aux_l, &orig_l);
    ////===========================================================
    real V[3], origv[3];
    V[0] = vl; 
    V[1] = vm; 
    V[2] = vr;
    
    origv[0] = orig_l;
    origv[1] = orig_mv;
    origv[2] = orig_r;
    
    set_common_orig_vertical_elvira(V, aux_mv, origv, 3, delta[1]);

    //printf("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n");

    Point Normal; 
    real areaE[6][9], distanceE[6][9], nx[6], ny[6], d[9], aux_x[9], aux_y[9];

    aux_x[0] = -1; aux_x[1] = 0; aux_x[2] = 1;
    aux_x[3] = -1; aux_x[4] = 0; aux_x[5] = 1;
    aux_x[6] = -1; aux_x[7] = 0; aux_x[8] = 1;

    aux_y[0] = -1; aux_y[1] = -1; aux_y[2] = -1;
    aux_y[3] = 0;  aux_y[4] = 0;  aux_y[5] = 0;
    aux_y[6] = 1;  aux_y[7] = 1;  aux_y[8] = 1;

    //=============================================================
    //// Normals, distances and area PLIC
    //=============================================================
    ////horizontal-regressive finite difference 
    elvira_calculate_normal_cell_regressive_1st_order_finite_difference_Horizontal(ns, clid, Normal, H[1], H[0], delta[0], delta[1], aux_mh);
    nx[0] = Normal[0]; ny[0] = Normal[1];
    distanceE[0][4] = distance_from_center(Normal,delta,area[4]);
    for (int j=0; j<9; j++) {
       d[j]        = distanceE[0][4] - (aux_x[j]*Normal[0]*delta[0] + aux_y[j]*Normal[1]*delta[1]);
       areaE[0][j] = area_left_line_origin_center(Normal, delta, d[j]);
       areaE[0][j] = area_correction_at_get(delta, areaE[0][j]);
    }
    //=============================================================
    ////horizontal-central finite difference 
    elvira_calculate_normal_cell_central_2nd_order_finite_difference_Horizontal(ns, clid, Normal, H[2], H[0], delta[0], delta[1], aux_mh);
    nx[1] = Normal[0]; ny[1] = Normal[1];
    distanceE[1][4] = distance_from_center(Normal,delta,area[4]);
    for (int j=0; j<9; j++) {
       d[j]        = distanceE[1][4] - (aux_x[j]*Normal[0]*delta[0] + aux_y[j]*Normal[1]*delta[1]);
       areaE[1][j] = area_left_line_origin_center(Normal, delta, d[j]);
       areaE[1][j] = area_correction_at_get(delta, areaE[1][j]);
    }
    //=============================================================
    ////horizontal-progressive finite difference 
    elvira_calculate_normal_cell_progressive_1st_order_finite_difference_Horizontal(ns, clid, Normal, H[2], H[1], delta[0], delta[1], aux_mh);
    nx[2] = Normal[0]; ny[2] = Normal[1];
    distanceE[2][4] = distance_from_center(Normal,delta,area[4]);
    for (int j=0; j<9; j++) {
       d[j]        = distanceE[2][4] - (aux_x[j]*Normal[0]*delta[0] + aux_y[j]*Normal[1]*delta[1]);
       areaE[2][j] = area_left_line_origin_center(Normal, delta, d[j]);
       areaE[2][j] = area_correction_at_get(delta, areaE[2][j]);
    }
    //=============================================================
    ////vertical-progressive finite difference 
    elvira_calculate_normal_cell_progressive_1st_order_finite_difference_Vertical(ns, clid, Normal, V[2], V[1], delta[0], delta[1], aux_mv);
    nx[3] = Normal[0]; ny[3] = Normal[1];
    distanceE[3][4] = distance_from_center(Normal,delta,area[4]);
    for (int j=0; j<9; j++) {
       d[j]        = distanceE[3][4] - (aux_x[j]*Normal[0]*delta[0] + aux_y[j]*Normal[1]*delta[1]);
       areaE[3][j] = area_left_line_origin_center(Normal, delta, d[j]);
       areaE[3][j] = area_correction_at_get(delta, areaE[3][j]);
    }
    ////===========================================================
    ////vertical-central finite difference 
    elvira_calculate_normal_cell_central_2nd_order_finite_difference_Vertical(ns, clid, Normal, V[2], V[0], delta[0], delta[1], aux_mv);
    nx[4] = Normal[0]; ny[4] = Normal[1];
    distanceE[4][4] = distance_from_center(Normal,delta,area[4]);
    for (int j=0; j<9; j++) {
       d[j]        = distanceE[4][4] - (aux_x[j]*Normal[0]*delta[0] + aux_y[j]*Normal[1]*delta[1]);
       areaE[4][j] = area_left_line_origin_center(Normal, delta, d[j]);
       areaE[4][j] = area_correction_at_get(delta, areaE[4][j]);
    }
    ////===========================================================
    ////vertical-regressive finite difference 
    elvira_calculate_normal_cell_regressive_1st_order_finite_difference_Vertical(ns, clid, Normal, V[1], V[0], delta[0], delta[1], aux_mv);
    nx[5] = Normal[0]; ny[5] = Normal[1];
    distanceE[5][4] = distance_from_center(Normal,delta,area[4]);
    for (int j=0; j<9; j++) {
       d[j]        = distanceE[5][4] - (aux_x[j]*Normal[0]*delta[0] + aux_y[j]*Normal[1]*delta[1]);
       areaE[5][j] = area_left_line_origin_center(Normal, delta, d[j]);
       areaE[5][j] = area_correction_at_get(delta, areaE[5][j]);
    }
    //============================================================
    real err_min = 1.0e32;
    real errors;
    int posi=0;
    for (int s=0; s<6; s++){
       errors = elvira_interface_error(area, areaE, s);
       if (errors < err_min){
          err_min = errors;
          Normal[0] = nx[s];
          Normal[1] = ny[s];
          posi = s;
       }
       //printf("s = %d \t nx = %.16lf \t ny = %.16lf\n", s, nx[s], ny[s]);
    }
    //printf("######################################################################################\n");
    //printf("x = %.16lf \t y = %.16lf\n", p[0], p[1]);
    //printf("position = %d \t nx = %.16lf \t ny = %.16lf\n", posi, nx[posi], ny[posi]);
    //getchar();
    for (int i=0; i<DIM; i++){
       dp_set_value(ns->ed.mult.dpnormal[i], clid, Normal[i]);
    }
}
