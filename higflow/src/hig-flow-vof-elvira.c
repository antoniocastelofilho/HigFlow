#include "hig-flow-vof-elvira.h"

real area_correction_at_get(Point delta, real area){
   if(fabs(area - delta[0]*delta[1])<1.0e-8){
      area = delta[0]*delta[1];
   }else if(fabs(area)<1.0e-8){
      area = 0.0;
   }
   return area;
}

real direction_hf_horizontal(sim_domain *sdp, higflow_solver *ns, Point center, Point delta){
   real frac, fracvol;
   Point p;
   real hf_horizontal=0.0;
   
   p[0] = center[0]+delta[0];
   //p[1] = center[1];
   real hf_hor_dir = 0.0;
   int i = 0;
   do {
      p[1] = center[1] + (i*delta[1]-delta[1]);
      real frac = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
      hf_hor_dir += frac;
      i++;
   } while (i < 3);
   p[0] = center[0]-delta[0];
   real hf_hor_esq = 0.0;
   i = 0;
   do {
      p[1] = center[1] + (i*delta[1]-delta[1]);
      frac = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
      hf_hor_esq += frac;
      i++;
   } while (i < 3);
   
   hf_horizontal = fabs(hf_hor_dir - hf_hor_esq);
   //printf("%lf %lf %lf\n",hf_hor_dir, hf_hor_esq, hf_horizontal);
   return hf_horizontal;
}

real direction_hf_vertical(sim_domain *sdp, higflow_solver *ns, Point center, Point delta){
   real frac, fracvol;
   Point p;
   real hf_vertical=0.0;
   
   p[1] = center[1]+delta[1];
   //p[1] = center[1];
   real hf_ver_sup = 0.0;
   int i = 0;
   do {
      p[0] = center[0] + (i*delta[0]-delta[0]);
      real frac = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
      hf_ver_sup += frac;
      i++;
   } while (i < 3);
   p[1] = center[1]-delta[1];
   real hf_ver_inf = 0.0;
   i = 0;
   do {
      p[0] = center[0] + (i*delta[0]-delta[0]);
      frac = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
      hf_ver_inf += frac;
      i++;
   } while (i < 3);
   
   hf_vertical = fabs(hf_ver_sup - hf_ver_inf);
   //printf("%lf %lf %lf\n",hf_ver_dir, hf_ver_esq, hf_vertical);
   return hf_vertical;
}

void ELVIRA_vertical_collumn(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *vertical, int *aux){
   int status;
   real fracvol, fracvol_aux;
   Point pp, ppt, ppb;
   *aux  = 0;
   // Vertical Up
   status = get_frac_vol(sdp, ns, 1, center, p, delta, &fracvol);
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
   do { 
      i++;
      pp[1] = center[1] + i * delta[1];
      status = get_frac_vol(sdp, ns, 1, center, pp, delta, &fracvol);
      if (status == 1) {
         *vertical += fracvol;
      } else if (status == 0) {
         //printf("Up - Vertical cell out of domain \n");
         return;
      } else {
         //printf("Up - Vertical cells with different sizes \n");
         return;
      }
   } while (i < 3 && status == 1);
   fracvol_aux = fracvol;
   ppt[0] = pp[0]; 
   ppt[1] = pp[1] + 0.5*delta[1]; 
   real fracvol_top = fracvol;             
   // Vertical Down
   i = 0;
   do {
      i++;
      pp[1] = p[1] - i * delta[1];
      status = get_frac_vol(sdp, ns, 1, center, pp, delta, &fracvol);
      if (status == 1) {
         *vertical += fracvol;
      } else if (status == 0) {
         //printf("Down - Vertical cell out of domain \n");
         return;
      } else {
         //printf("Down - Vertical cells with different sizes \n");
         return;
      }
   } while (i < 3 && status == 1);
   if (fracvol == fracvol_aux) {
      //printf("The phases are not different \n");
      return;
   }
   if (fracvol == 1.0 && fracvol_top == 0.0) {
      *aux = -1;
   } else {
      *aux = 1;
   }
}

void ELVIRA_vertical_collumn_adap(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *vertical, int *aux, real *orig){
   int status;
   real fracvol, fracvol_aux;
   Point pp, ppt, ppb;
   *aux  = 0;
   *orig = 0.0;
   // Vertical Up
   status = get_frac_vol(sdp, ns, 1, center, p, delta, &fracvol);
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
   do { 
      i++;
      pp[1] = center[1] + i * delta[1];
      status = get_frac_vol(sdp, ns, 1, center, pp, delta, &fracvol);
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
   do {
      i++;
      pp[1] = p[1] - i * delta[1];
      status = get_frac_vol(sdp, ns, 1, center, pp, delta, &fracvol);
      if (status == 1) {
         *vertical += fracvol;
      } else if (status == 0) {
         printf("Down - Vertical cell out of domain \n");
         return;
      } else {
         printf("Down - Vertical cells with different sizes \n");
         return;
      }
   } while (fracvol > 0.0 && fracvol < 1.0 && fabs(fracvol - fracvol_aux)>1.0e-12 && status == 1);

   if (fracvol == fracvol_aux) {
      //printf("The phases are not different \n");
      return;
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

void ELVIRA_horizontal_row(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *horizontal, int *aux){
   int status;
   real fracvol, fracvol_aux;
   Point pp, ppr, ppl;
   *aux=0;
   // Horizontal Right
   status = get_frac_vol(sdp, ns, 0, center, p, delta, &fracvol);
   if (status != 1) {
      if (status == -1) {
         printf("Right - Horizontal cells with different sizes \n");
      } else {
         //printf("Right - Horizontal cell out of domain \n");
      }
      return;
   }
   pp[1] = p[1];
   *horizontal = fracvol;
   int i = 0;
   do {
      i++;
      pp[0] = center[0] + i*delta[0];
      status = get_frac_vol(sdp, ns, 0, center, pp, delta, &fracvol);
      if (status == 1) {
         *horizontal += fracvol;
      } else if (status == 0) {
         //printf("Right - Horizontal cell out of domain \n");
         return;
      } else {
         printf("Right - Horizontal cells with different sizes \n");
         return;
      }
   } while (i < 3 && status == 1);
   fracvol_aux = fracvol;
   ppr[0]=pp[0] + 0.5*delta[0];
   ppr[1]=pp[1];
   real fracvol_right = fracvol;
   // Horizontal Left
   i = 0;
   // fracvol_aux = 0;
   do {
      i++;
      pp[0] = p[0] - i*delta[0];
      status = get_frac_vol(sdp, ns, 0, center, pp, delta, &fracvol);
      if (status == 1) {
         *horizontal += fracvol;
      } else if (status == 0) {
         //printf("Left - Horizontal cell out of domain \n");
         return;
      } else {
         printf("Left - Horizontal cells with different sizes \n");
         return;
      }
   } while (i < 3 && status == 1);
   if (fracvol == fracvol_aux){
      //printf("The phases are not different:\n");
      return;
   }
   if (fracvol == 1.0 && fracvol_right == 0.0) {
      *aux = -1;
   } else {
      *aux = 1;
   }
}

void ELVIRA_horizontal_row_adap(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *horizontal, int *aux, real *orig){
   int status;
   real fracvol, fracvol_aux;
   Point pp, ppr, ppl;
   *aux=0;
   *orig=0.0;
   // Horizontal Right
   status = get_frac_vol(sdp, ns, 0, center, p, delta, &fracvol);
   if (status != 1) {
      if (status == -1) {
         printf("Right - Horizontal cells with different sizes \n");
      } else {
         printf("Right - Horizontal cell out of domain \n");
      }
      return;
   }
   pp[1] = p[1];
   *horizontal=fracvol;
   int i = 0;
   do {
      i++;
      pp[0] = center[0] + i*delta[0];
      status=get_frac_vol(sdp, ns, 0, center, pp, delta, &fracvol);
      if (status == 1) {
         *horizontal += fracvol;
      } else if (status == 0) {
         printf("Right - Horizontal cell out of domain \n");
         return;
      } else {
         printf("Right - Horizontal cells with different sizes \n");
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
   do {
      i++;
      pp[0] = p[0] - i*delta[0];
      status = get_frac_vol(sdp, ns, 0, center, pp, delta, &fracvol);
      if (status == 1) {
         *horizontal += fracvol;
      } else if (status == 0) {
         printf("Left - Horizontal cell out of domain \n");
         return;
      } else {
         printf("Left - Horizontal cells with different sizes \n");
         return;
      }
   }while (fracvol > 0.0 && fracvol < 1.0 && fabs(fracvol - fracvol_aux)>1.0e-12 && status == 1);

   if (fracvol == fracvol_aux){
      //printf("The phases are not different:\n");
      return;
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
  //printf("error  = %.16e\n",error);
  return sqrt(error);
}

void higflow_compute_curvature_interfacial_force_normal_multiphase_2D_EL(higflow_solver *ns) {
      real IF[DIM];
      // Get the local sub-domain for the cells
      sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
      // Get the map for the domain properties
      mp_mapper *mp = sd_get_domain_mapper(sdp);
      // Loop for each cell
      higcit_celliterator *it;
      for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
         // Get the cell
         hig_cell *c = higcit_getcell(it);
         // Get the cell identifier
         int clid = mp_lookup(mp, hig_get_cid(c));
         // Get the center of the cell
         Point center;
         hig_get_center(c, center);
         // Get the delta of the cell
         Point delta;
         hig_get_delta(c, delta);
         // Case bi-dimensional
         Point p;
         real area[9];
         p[0] = center[0];   p[1] = center[1];
         real fracvol = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
         
         //dp_set_value(ns->ed.mult.dpcurvature, clid, 0.0);/* Set the curvature in the distributed curvature property*/
         for (int i = 0; i < DIM; i++) {
            dp_set_value(ns->ed.mult.dpIF[i], clid, 0.0);
            dp_set_value(ns->ed.mult.dpnormal[i], clid, 0.0);
         }
         if (fracvol < 1.0e-14 || fabs(fracvol - 1.0) < 1.0e-14){
            continue;
         }
         
         real hf_horizontal = direction_hf_horizontal(sdp, ns, center, delta); 
         //height deriction verification
         real hf_vertical = direction_hf_vertical(sdp, ns, center, delta); //height deriction verification
         
         /*if(hf_horizontal > hf_vertical){
            calculate_exact_normal_x_dominant_2D(ns, clid, p);
         } else{
            calculate_exact_normal_y_dominant_2D(ns, clid, p);
         }*/
         
         area[4] = fracvol*delta[0]*delta[1];
         
         //normal calculation ======================================
         
         real hm, hb, ht, vm, vl, vr;
         int aux_mh = 0, aux_b = 0, aux_t = 0;
         int aux_mv = 0, aux_l = 0, aux_r = 0;
         
         // Middle
         // Horizontal: going middle
         ELVIRA_horizontal_row(sdp, ns, center, p, delta, &hm, &aux_mh);
         
         // Top
         p[1] = center[1] + delta[1];
         // Volume Faction and Area (required for error calculation)
         fracvol = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
         area[7] = fracvol*delta[0]*delta[1];
         // Horizontal: going top
         ELVIRA_horizontal_row(sdp, ns, center, p, delta, &ht, &aux_t);
         
         // Bottom
         p[1] = center[1] - delta[1];
         // Volume Faction and Area (required for error calculation)
         fracvol = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
         area[1] = fracvol*delta[0]*delta[1];
         //Horizontal: going BOTTOM
         ELVIRA_horizontal_row(sdp, ns, center, p, delta, &hb, &aux_b);
         
         real H[3];
         H[0] = hb; H[1] = hm; H[2] = ht;
         p[0] = center[0]; p[1] = center[1];
         // Middle
         // Vertical: going middle
         ELVIRA_vertical_collumn(sdp, ns, center, p, delta, &vm, &aux_mv);
         // Right
         p[0] = center[0] + delta[0];
         // Volume Faction and Area (required for error calculation)
         fracvol = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
         area[5] = fracvol*delta[0]*delta[1];
         // Vertical: going Right
         ELVIRA_vertical_collumn(sdp, ns, center, p, delta, &vr, &aux_r);
         // Left
         p[0] = center[0] - delta[0];
         // Volume Faction and Area (required for error calculation)
         fracvol = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
         area[3] = fracvol*delta[0]*delta[1];
         // Vertical: going Left
         ELVIRA_vertical_collumn(sdp, ns, center, p, delta, &vl, &aux_l);
         
         real V[3];
         V[0] = vl; V[1] = vm; V[2] = vr;
         
         // Volume Faction and Area in the Corns (required for error calculation)
         // Left-Botton
         p[0] = center[0] - delta[0];
         p[1] = center[1] - delta[1];
         fracvol = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
         area[0] = fracvol*delta[0]*delta[1];
         // Left-Top
         p[1] = center[1] + delta[1];
         fracvol = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
         area[6] = fracvol*delta[0]*delta[1];
         // Right-Botton
         p[0] = center[0] + delta[0];
         p[1] = center[1] - delta[1];
         fracvol = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
         area[2] = fracvol*delta[0]*delta[1];
         // Right-Top
         p[1] = center[1] + delta[1];
         fracvol = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
         area[8] = fracvol*delta[0]*delta[1];
         
         Point Normal;
         real nx[6], ny[6];
         real areaE[6][9], distanceE[6][9], d[9];
         real dl;
         real sinalXl[9], sinalYl[9];
         sinalXl[0] = -1;sinalXl[1] = 0;sinalXl[2] = 1;
         sinalXl[3] = -1;sinalXl[4] = 0;sinalXl[5] = 1;
         sinalXl[6] = -1;sinalXl[7] = 0;sinalXl[8] = 1;

         sinalYl[0] = -1; sinalYl[1] = -1; sinalYl[2] = -1;
         sinalYl[3] = 0;  sinalYl[4] = 0;  sinalYl[5] = 0;
         sinalYl[6] = 1;  sinalYl[7] = 1;  sinalYl[8] = 1;
         //==========================================================
         ////horizontal-regressive finite difference 
         ELVIRA_calculate_normal_cell_regressive_1st_order_finite_difference_Horizontal(ns, clid, Normal, hm, hb, delta[0], delta[1], aux_mh);
         normal_correction_at_get(Normal);
         nx[0] = Normal[0]; ny[0] = Normal[1];
         distanceE[0][4] = distance_from_center(Normal,delta,area[4]);
         for (int j=0; j<9; j++) {
            d[j]        = distanceE[0][4] - (sinalXl[j]*Normal[0]*delta[0] + sinalYl[j]*Normal[1]*delta[1]);
            areaE[0][j] = area_left_line_origin_center(Normal, delta, d[j]);
            areaE[0][j] = area_correction_at_get(delta , areaE[0][j]);
         }
         //==========================================================
         ////horizontal-central finite difference 
         ELVIRA_calculate_normal_cell_central_2nd_order_finite_difference_Horizontal(ns, clid, Normal, ht, hb, delta[0], delta[1], aux_mh);
         normal_correction_at_get(Normal);
         nx[1] = Normal[0]; ny[1] = Normal[1];
         distanceE[1][4] = distance_from_center(Normal,delta,area[4]);
         for (int j=0; j<9; j++) {
            d[j]        = distanceE[1][4] - (sinalXl[j]*Normal[0]*delta[0] + sinalYl[j]*Normal[1]*delta[1]);
            //distanceE[0][j] = distance_from_center(Normal,delta,area[j]);
            areaE[1][j] = area_left_line_origin_center(Normal, delta, d[j]);
            areaE[1][j] = area_correction_at_get(delta , areaE[1][j]);
         }
         //==========================================================
         ////horizontal-progressive finite difference 
         ELVIRA_calculate_normal_cell_progressive_1st_order_finite_difference_Horizontal(ns, clid, Normal, ht, hm, delta[0], delta[1], aux_mh);
         normal_correction_at_get(Normal);
         nx[2] = Normal[0]; ny[2] = Normal[1];
         distanceE[2][4] = distance_from_center(Normal,delta,area[4]);
         for (int j=0; j<9; j++) {
            d[j]        = distanceE[2][4] - (sinalXl[j]*Normal[0]*delta[0] + sinalYl[j]*Normal[1]*delta[1]);
            areaE[2][j] = area_left_line_origin_center(Normal, delta, d[j]);
            areaE[2][j] = area_correction_at_get(delta , areaE[2][j]);
         }
         //==========================================================
         ////vertical-progressive finite difference 
         ELVIRA_calculate_normal_cell_progressive_1st_order_finite_difference_Vertical(ns, clid, Normal, vr, vm, delta[0], delta[1], aux_mv);
         normal_correction_at_get(Normal);
         nx[3] = Normal[0]; ny[3] = Normal[1];
         distanceE[3][4] = distance_from_center(Normal,delta,area[4]);
         for (int j=0; j<9; j++) {
            d[j]        = distanceE[3][4] - (sinalXl[j]*Normal[0]*delta[0] + sinalYl[j]*Normal[1]*delta[1]);
            areaE[3][j] = area_left_line_origin_center(Normal, delta, d[j]);
            areaE[3][j] = area_correction_at_get(delta , areaE[3][j]);
         }
         //==========================================================
         //vertical-central finite difference 
         ELVIRA_calculate_normal_cell_central_2nd_order_finite_difference_Vertical(ns, clid, Normal, vr, vl, delta[0], delta[1], aux_mv);
         normal_correction_at_get(Normal);
         nx[4] = Normal[0]; ny[4] = Normal[1];
         distanceE[4][4] = distance_from_center(Normal,delta,area[4]);
         for (int j=0; j<9; j++) {
            d[j]        = distanceE[4][4] - (sinalXl[j]*Normal[0]*delta[0] + sinalYl[j]*Normal[1]*delta[1]);
            areaE[4][j] = area_left_line_origin_center(Normal, delta, d[j]);
            areaE[4][j] = area_correction_at_get(delta , areaE[4][j]);
         }
         //==========================================================
         ////vertical-regressive finite difference 
         ELVIRA_calculate_normal_cell_regressive_1st_order_finite_difference_Vertical(ns, clid, Normal, vm, vl, delta[0], delta[1], aux_mv);
         normal_correction_at_get(Normal);
         nx[5] = Normal[0]; ny[5] = Normal[1];
         distanceE[5][4] = distance_from_center(Normal,delta,area[4]);
         for (int j=0; j<9; j++) {
            d[j]        = distanceE[5][4] - (sinalXl[j]*Normal[0]*delta[0] + sinalYl[j]*Normal[1]*delta[1]);
            areaE[5][j] = area_left_line_origin_center(Normal, delta, d[j]);
            areaE[5][j] = area_correction_at_get(delta , areaE[5][j]);
         }
         //==========================================================
         real err_min = 1.0e32;
         real errors;
         for (int s=0; s<6; s++){
            errors = elvira_interface_error(area, areaE, s);
            if (errors < err_min){
               err_min = errors;
               Normal[0] = nx[s];
               Normal[1] = ny[s];
            }
         }
         for (int i=0; i<DIM; i++){
            dp_set_value(ns->ed.mult.dpnormal[i], clid, Normal[i]);
         }
      }
      // Destroy the iterator
      higcit_destroy(it);
      // Sync the distributed pressure property
      //dp_sync(ns->ed.mult.dpcurvature);
      for (int i = 0; i < DIM; i++) {
         //dp_sync(ns->ed.mult.dpIF[i]);
         dp_sync(ns->ed.mult.dpnormal[i]);
      }
}

void higflow_compute_normal_multiphase_2D_elvira(higflow_solver *ns, sim_domain *sdp, mp_mapper *mp, higcit_celliterator *it, hig_cell *c,int clid, Point center, Point delta) {
      real IF[DIM];
      // Get the local sub-domain for the cells
      //sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
      // Get the map for the domain properties
      //mp_mapper *mp = sd_get_domain_mapper(sdp);
      // Loop for each cell
      //higcit_celliterator *it;
//      for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
         // Get the cell
         //hig_cell *c = higcit_getcell(it);
         // Get the cell identifier
         //int clid = mp_lookup(mp, hig_get_cid(c));
         // Get the center of the cell
         //Point center;
         //hig_get_center(c, center);
         // Get the delta of the cell
         //Point delta;
         //hig_get_delta(c, delta);
         // Case bi-dimensional
         Point p;
         real area[9];
         p[0] = center[0];   p[1] = center[1];
         real fracvol = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
         
         real hf_horizontal = direction_hf_horizontal(sdp, ns, center, delta); 
         //height deriction verification
         real hf_vertical = direction_hf_vertical(sdp, ns, center, delta); //height deriction verification
         
         area[4] = fracvol*delta[0]*delta[1];
         
         //normal calculation ======================================
         real hm, hb, ht, vm, vl, vr;
         int aux_mh = 0, aux_b = 0, aux_t = 0;
         int aux_mv = 0, aux_l = 0, aux_r = 0;
         
         // Middle
         // Horizontal: going middle
         ELVIRA_horizontal_row(sdp, ns, center, p, delta, &hm, &aux_mh);
         
         // Top
         p[1] = center[1] + delta[1];
         // Volume Faction and Area (required for error calculation)
         fracvol = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
         area[7] = fracvol*delta[0]*delta[1];
         // Horizontal: going top
         ELVIRA_horizontal_row(sdp, ns, center, p, delta, &ht, &aux_t);
         
         // Bottom
         p[1] = center[1] - delta[1];
         // Volume Faction and Area (required for error calculation)
         fracvol = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
         area[1] = fracvol*delta[0]*delta[1];
         //Horizontal: going BOTTOM
         ELVIRA_horizontal_row(sdp, ns, center, p, delta, &hb, &aux_b);
         
         real H[3];
         H[0] = hb; H[1] = hm; H[2] = ht;
         p[0] = center[0]; p[1] = center[1];
         // Middle
         // Vertical: going middle
         ELVIRA_vertical_collumn(sdp, ns, center, p, delta, &vm, &aux_mv);
         // Right
         p[0] = center[0] + delta[0];
         // Volume Faction and Area (required for error calculation)
         fracvol = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
         area[5] = fracvol*delta[0]*delta[1];
         // Vertical: going Right
         ELVIRA_vertical_collumn(sdp, ns, center, p, delta, &vr, &aux_r);
         // Left
         p[0] = center[0] - delta[0];
         // Volume Faction and Area (required for error calculation)
         fracvol = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
         area[3] = fracvol*delta[0]*delta[1];
         // Vertical: going Left
         ELVIRA_vertical_collumn(sdp, ns, center, p, delta, &vl, &aux_l);
         
         real V[3];
         V[0] = vl; V[1] = vm; V[2] = vr;
         
         // Volume Faction and Area in the Corns (required for error calculation)
         // Left-Botton
         p[0] = center[0] - delta[0];
         p[1] = center[1] - delta[1];
         fracvol = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
         area[0] = fracvol*delta[0]*delta[1];
         // Left-Top
         p[1] = center[1] + delta[1];
         fracvol = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
         area[6] = fracvol*delta[0]*delta[1];
         // Right-Botton
         p[0] = center[0] + delta[0];
         p[1] = center[1] - delta[1];
         fracvol = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
         area[2] = fracvol*delta[0]*delta[1];
         // Right-Top
         p[1] = center[1] + delta[1];
         fracvol = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
         area[8] = fracvol*delta[0]*delta[1];
         
         Point Normal;
         real nx[6], ny[6];
         real areaE[6][9], distanceE[6][9], d[9];
         real dl;
         real sinalXl[9], sinalYl[9];
         sinalXl[0] = -1;sinalXl[1] = 0;sinalXl[2] = 1;
         sinalXl[3] = -1;sinalXl[4] = 0;sinalXl[5] = 1;
         sinalXl[6] = -1;sinalXl[7] = 0;sinalXl[8] = 1;

         sinalYl[0] = -1; sinalYl[1] = -1; sinalYl[2] = -1;
         sinalYl[3] = 0;  sinalYl[4] = 0;  sinalYl[5] = 0;
         sinalYl[6] = 1;  sinalYl[7] = 1;  sinalYl[8] = 1;
         //==========================================================
         ////horizontal-regressive finite difference 
         ELVIRA_calculate_normal_cell_regressive_1st_order_finite_difference_Horizontal(ns, clid, Normal, hm, hb, delta[0], delta[1], aux_mh);
         normal_correction_at_get(Normal);
         nx[0] = Normal[0]; ny[0] = Normal[1];
         distanceE[0][4] = distance_from_center(Normal,delta,area[4]);
         for (int j=0; j<9; j++) {
            d[j]        = distanceE[0][4] - (sinalXl[j]*Normal[0]*delta[0] + sinalYl[j]*Normal[1]*delta[1]);
            areaE[0][j] = area_left_line_origin_center(Normal, delta, d[j]);
            areaE[0][j] = area_correction_at_get(delta , areaE[0][j]);
         }
         //==========================================================
         ////horizontal-central finite difference 
         ELVIRA_calculate_normal_cell_central_2nd_order_finite_difference_Horizontal(ns, clid, Normal, ht, hb, delta[0], delta[1], aux_mh);
         normal_correction_at_get(Normal);
         nx[1] = Normal[0]; ny[1] = Normal[1];
         distanceE[1][4] = distance_from_center(Normal,delta,area[4]);
         for (int j=0; j<9; j++) {
            d[j]        = distanceE[1][4] - (sinalXl[j]*Normal[0]*delta[0] + sinalYl[j]*Normal[1]*delta[1]);
            //distanceE[0][j] = distance_from_center(Normal,delta,area[j]);
            areaE[1][j] = area_left_line_origin_center(Normal, delta, d[j]);
            areaE[1][j] = area_correction_at_get(delta , areaE[1][j]);
         }
         //==========================================================
         ////horizontal-progressive finite difference 
         ELVIRA_calculate_normal_cell_progressive_1st_order_finite_difference_Horizontal(ns, clid, Normal, ht, hm, delta[0], delta[1], aux_mh);
         normal_correction_at_get(Normal);
         nx[2] = Normal[0]; ny[2] = Normal[1];
         distanceE[2][4] = distance_from_center(Normal,delta,area[4]);
         for (int j=0; j<9; j++) {
            d[j]        = distanceE[2][4] - (sinalXl[j]*Normal[0]*delta[0] + sinalYl[j]*Normal[1]*delta[1]);
            areaE[2][j] = area_left_line_origin_center(Normal, delta, d[j]);
            areaE[2][j] = area_correction_at_get(delta , areaE[2][j]);
         }
         //==========================================================
         ////vertical-progressive finite difference 
         ELVIRA_calculate_normal_cell_progressive_1st_order_finite_difference_Vertical(ns, clid, Normal, vr, vm, delta[0], delta[1], aux_mv);
         normal_correction_at_get(Normal);
         nx[3] = Normal[0]; ny[3] = Normal[1];
         distanceE[3][4] = distance_from_center(Normal,delta,area[4]);
         for (int j=0; j<9; j++) {
            d[j]        = distanceE[3][4] - (sinalXl[j]*Normal[0]*delta[0] + sinalYl[j]*Normal[1]*delta[1]);
            areaE[3][j] = area_left_line_origin_center(Normal, delta, d[j]);
            areaE[3][j] = area_correction_at_get(delta , areaE[3][j]);
         }
         //==========================================================
         //vertical-central finite difference 
         ELVIRA_calculate_normal_cell_central_2nd_order_finite_difference_Vertical(ns, clid, Normal, vr, vl, delta[0], delta[1], aux_mv);
         normal_correction_at_get(Normal);
         nx[4] = Normal[0]; ny[4] = Normal[1];
         distanceE[4][4] = distance_from_center(Normal,delta,area[4]);
         for (int j=0; j<9; j++) {
            d[j]        = distanceE[4][4] - (sinalXl[j]*Normal[0]*delta[0] + sinalYl[j]*Normal[1]*delta[1]);
            areaE[4][j] = area_left_line_origin_center(Normal, delta, d[j]);
            areaE[4][j] = area_correction_at_get(delta , areaE[4][j]);
         }
         //==========================================================
         ////vertical-regressive finite difference 
         ELVIRA_calculate_normal_cell_regressive_1st_order_finite_difference_Vertical(ns, clid, Normal, vm, vl, delta[0], delta[1], aux_mv);
         normal_correction_at_get(Normal);
         nx[5] = Normal[0]; ny[5] = Normal[1];
         distanceE[5][4] = distance_from_center(Normal,delta,area[4]);
         for (int j=0; j<9; j++) {
            d[j]        = distanceE[5][4] - (sinalXl[j]*Normal[0]*delta[0] + sinalYl[j]*Normal[1]*delta[1]);
            areaE[5][j] = area_left_line_origin_center(Normal, delta, d[j]);
            areaE[5][j] = area_correction_at_get(delta , areaE[5][j]);
         }
         //==========================================================
         real err_min = 1.0e32;
         real errors;
         for (int s=0; s<6; s++){
            errors = elvira_interface_error(area, areaE, s);
            if (errors < err_min){
               err_min = errors;
               Normal[0] = nx[s];
               Normal[1] = ny[s];
            }
         }
         for (int i=0; i<DIM; i++){
            dp_set_value(ns->ed.mult.dpnormal[i], clid, Normal[i]);
         }
      //}
      // Destroy the iterator
      //higcit_destroy(it);
      // Sync the distributed pressure property
      for (int i = 0; i < DIM; i++) {
         dp_sync(ns->ed.mult.dpnormal[i]);
      }
}

void higflow_compute_curvature_interfacial_force_normal_multiphase_2D_EL_adap(higflow_solver *ns) {
   real IF[DIM];
   // Get the local sub-domain for the cells
   sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
   // Get the map for the domain properties
   mp_mapper *mp = sd_get_domain_mapper(sdp);
   // Loop for each cell
   higcit_celliterator *it;
   for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
      // Get the cell
      hig_cell *c = higcit_getcell(it);
      // Get the cell identifier
      int clid = mp_lookup(mp, hig_get_cid(c));
      // Get the center of the cell
      Point center;
      hig_get_center(c, center);
      // Get the delta of the cell
      Point delta;
      hig_get_delta(c, delta);
      
      // Case bi-dimensional
      Point p;
      real area[9]; real fracvol[9];
      p[0] = center[0];   p[1] = center[1];
      fracvol[4] = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
      
      dp_set_value(ns->ed.mult.dpcurvature, clid, 0.0);/* Set the curvature in the distributed curvature property*/
      for (int i = 0; i < DIM; i++) {
         dp_set_value(ns->ed.mult.dpIF[i], clid, 0.0);
         dp_set_value(ns->ed.mult.dpnormal[i], clid, 0.0);
      }
      if (fracvol[4] < 1.0e-14 || fabs(fracvol[4] - 1.0) < 1.0e-14){
         continue;
      }
      
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
      // +-------+-------+-------+--+--> x
   
      area[4] = fracvol[4]*delta[0]*delta[1];
      
      p[0] = center[0] - delta[0]; p[1] = center[1] - delta[1];
      fracvol[0] = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
      area[0] = fracvol[0]*delta[0]*delta[1];
      
      p[0] = center[0]; p[1] = center[1] - delta[1];
      fracvol[1] = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
      area[1] = fracvol[1]*delta[0]*delta[1];
      
      p[0] = center[0] + delta[0]; p[1] = center[1] - delta[1];
      fracvol[2] = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
      area[2] = fracvol[2]*delta[0]*delta[1];
      
      p[0] = center[0] - delta[0]; p[1] = center[1];
      fracvol[3] = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
      area[3] = fracvol[3]*delta[0]*delta[1];
      
      p[0] = center[0] + delta[0]; p[1] = center[1];
      fracvol[5] = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
      area[5] = fracvol[5]*delta[0]*delta[1];
      
      p[0] = center[0] - delta[0]; p[1] = center[1] + delta[1];
      fracvol[6] = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
      area[6] = fracvol[6]*delta[0]*delta[1];
      
      p[0] = center[0]; p[1] = center[1] + delta[1];
      fracvol[7] = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
      area[7] = fracvol[7]*delta[0]*delta[1];
      
      p[0] = center[0] + delta[0]; p[1] = center[1] + delta[1];
      fracvol[8] = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
      area[8] = fracvol[8]*delta[0]*delta[1];
      
      //Height functions
      real hm, hb, ht, vm, vl, vr;
      int aux_mh = 0, aux_b = 0, aux_t = 0;
      int aux_mv = 0, aux_l = 0, aux_r = 0;
      real orig_mh, orig_b, orig_t;
      real orig_mv, orig_l, orig_r;
      
      //// Horizontal
      p[0] = center[0]; p[1] = center[1];
      // Middle
      ELVIRA_horizontal_row_adap(sdp, ns, center, p, delta, &hm, &aux_mh, &orig_mh);
      ////===========================================================
      // Top
      p[1] = center[1]+delta[1];
      ELVIRA_horizontal_row_adap(sdp, ns, center, p, delta, &ht, &aux_t, &orig_t);
      ////===========================================================
      //Bottom
      p[1] = center[1]-delta[1];
      ELVIRA_horizontal_row_adap(sdp, ns, center, p, delta, &hb, &aux_t, &orig_b);
      ////===========================================================
      real H[3];
      H[0] = hb; H[1] = hm; H[2] = ht;
      
      real origh[3];
      origh[0] = orig_b;
      origh[1] = orig_mh;
      origh[2] = orig_t;
      set_common_orig_horizontal_elvira(H, aux_mh, origh, 3, delta[0]);
      
      //// Vertical
      p[0] = center[0]; p[1] = center[1];
      // Middle
      ELVIRA_vertical_collumn_adap(sdp, ns, center, p, delta, &vm, &aux_mv, &orig_mv);
      ////===========================================================
      //Right
      p[0] = center[0]+delta[0];
      ELVIRA_vertical_collumn_adap(sdp, ns, center, p, delta, &vr, &aux_mv, &orig_r);
      ////===========================================================
      //Left
      p[0] = center[0]-delta[0];
      ELVIRA_vertical_collumn_adap(sdp, ns, center, p, delta, &vl, &aux_mv, &orig_l);
      ////===========================================================
      real V[3];
      V[0] = vl; V[1] = vm; V[2] = vr;
      
      real origv[3];
      origv[0] = orig_l;
      origv[1] = orig_mv;
      origv[2] = orig_r;
      set_common_orig_vertical_elvira(V, aux_mv, origv, 3, delta[1]);
      
      
      Point Normal; real areaE[6][9], distanceE[6][9];
      real nx[6], ny[6];
      //=============================================================
      //// Normals, distances and area PLIC
      //=============================================================
      ////horizontal-regressive finite difference 
      ELVIRA_calculate_normal_cell_regressive_1st_order_finite_difference_Horizontal(ns, 
      clid, Normal, H[1], H[0], delta[0], delta[1], aux_mh);
      nx[0] = Normal[0]; ny[0] = Normal[1];
      
      for (int j=0; j<9; j++) {
         distanceE[0][j] = distance_from_center(Normal,delta,area[j]);
         areaE[0][j]     = area_left_line_origin_center(Normal, delta, distanceE[0][j]);
         //areaE[0][j]     = area_left_line_origin_center(Normal, delta, distanceE[0][j]+Normal[0]*delta[0]-Normal[1]*delta[1]);
         areaE[0][j]     = area_correction_at_get(delta , areaE[0][j]);
      }
      //=============================================================
      ////horizontal-central finite difference 
      ELVIRA_calculate_normal_cell_central_2nd_order_finite_difference_Horizontal(ns, 
      clid, Normal, H[2], H[0], delta[0], delta[1], aux_mh);
      nx[1] = Normal[0]; ny[1] = Normal[1];
      
      for (int j=0; j<9; j++) {
         distanceE[1][j] = distance_from_center(Normal,delta,area[j]);
         areaE[1][j]     = area_left_line_origin_center(Normal, delta, distanceE[1][j]);
         areaE[1][j]     = area_correction_at_get(delta , areaE[1][j]);
      }
      //=============================================================
      ////horizontal-progressive finite difference 
      ELVIRA_calculate_normal_cell_progressive_1st_order_finite_difference_Horizontal(ns, 
      clid, Normal, H[2], H[1], delta[0], delta[1], aux_mh);
      nx[2] = Normal[0]; ny[2] = Normal[1];
      
      for (int j=0; j<9; j++) {
         distanceE[2][j] = distance_from_center(Normal,delta,area[j]);
         areaE[2][j]     = area_left_line_origin_center(Normal, delta, distanceE[2][j]);
         areaE[2][j]     = area_correction_at_get(delta , areaE[2][j]);
      }
      //=============================================================
      ////vertical-progressive finite difference 
      ELVIRA_calculate_normal_cell_progressive_1st_order_finite_difference_Vertical(ns, 
      clid, Normal, V[2], V[1], delta[0], delta[1], aux_mv);
      nx[3] = Normal[0]; ny[3] = Normal[1];
      
      for (int j=0; j<9; j++) {
         distanceE[3][j] = distance_from_center(Normal,delta,area[j]);
         areaE[3][j]     = area_left_line_origin_center(Normal, delta, distanceE[3][j]);
         areaE[3][j]     = area_correction_at_get(delta , areaE[3][j]);
      }
      ////===========================================================
      ////vertical-central finite difference 
      ELVIRA_calculate_normal_cell_central_2nd_order_finite_difference_Vertical(ns, 
      clid, Normal, V[2], V[0], delta[0], delta[1], aux_mv);
      nx[4] = Normal[0]; ny[4] = Normal[1];
      
      for (int j=0; j<9; j++) {
         distanceE[4][j] = distance_from_center(Normal,delta,area[j]);
         areaE[4][j]     = area_left_line_origin_center(Normal, delta, distanceE[4][j]);
         areaE[4][j]     = area_correction_at_get(delta , areaE[4][j]);
      }
      //==========================================================
      ////vertical-regressive finite difference 
      ELVIRA_calculate_normal_cell_regressive_1st_order_finite_difference_Vertical(ns, clid, Normal, V[1], V[0], delta[0], delta[1], aux_mv);
      nx[5] = Normal[0]; ny[5] = Normal[1];
      for (int j=0; j<9; j++) {
         distanceE[5][j] = distance_from_center(Normal,delta,area[j]);
         areaE[5][j]     = area_left_line_origin_center(Normal, delta, distanceE[5][j]);
         areaE[5][j]     = area_correction_at_get(delta , areaE[5][j]);
      }
      //============================================================
      real err_min = 1.0e32;
      real errors;
      for (int s=0; s<6; s++){
        errors = elvira_interface_error(area, areaE, s);
         if (errors < err_min){
            err_min = errors;
            Normal[0] = nx[s];
            Normal[1] = ny[s];
         }
      }
      for (int i=0; i<DIM; i++){
         dp_set_value(ns->ed.mult.dpnormal[i], clid, Normal[i]);
      }
   }
   // Destroy the iterator
   higcit_destroy(it);
   // Sync the distributed pressure property
   //dp_sync(ns->ed.mult.dpcurvature);
   for (int i = 0; i < DIM; i++) {
      //dp_sync(ns->ed.mult.dpIF[i]);
      dp_sync(ns->ed.mult.dpnormal[i]);
   }
}
