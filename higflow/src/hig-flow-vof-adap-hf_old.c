#include "hig-flow-vof-finite-difference-normal-curvature.h"
#include "hig-flow-vof-adap-hf.h"
#include "hig-flow-vof-9-cells.h"
#include "hig-flow-vof-plic.h"
#include "hig-flow-vof-elvira.h"
#include "hig-flow-discret.h"

void fraction_correction_at_get(real *fracvol){
   real fracvol_tol = 1.0e-16;
   if (fabs(*fracvol - 1.0) < fracvol_tol) {
      *fracvol = 1.0;
   } else if (fabs(*fracvol) < fracvol_tol) {
      *fracvol = 0.0;
   }
}

int get_frac_vol(sim_domain *sd, higflow_solver *ns, int dim, Point Center, Point P,Point Delta,real *fracvol){
   hig_cell *c = sd_get_cell_with_point(sd,P);
   if (c == NULL) {
      return 0;
   } else {
      Point Delta2;
      hig_get_delta(c, Delta2);
      if (fabs(Delta[dim] - Delta2[dim])>1.0e-12){
         printf("Different sizes dc=[%.18lf %.18lf] dp= [%.18lf %.18lf] \t ",Delta[0],Delta[1],Delta2[0],Delta2[1]);
         return -1;
      } else {
         *fracvol = compute_value_at_point(sd, Center, P, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
         fraction_correction_at_get(fracvol);
         return 1;
      }
   }
}

void vertical_collumn(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *vertical, int *aux, real *orig){
   int status;
   real fracvol, fracvol_aux;
   Point pp, ppt, ppb;
   *aux=0;
   *orig=0.0;
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
   do { i++;
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
   do { i++;
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
   } while (fracvol > 0.0 && fracvol < 1.0 && fabs(fracvol - fracvol_aux) > 1.0e-14 && status == 1);

   if (fracvol == fracvol_aux) {
      //printf("The phases are not different \n");
      return;
   }
   ppb[0]=pp[0];ppb[1]=pp[1]-0.5*delta[1];

   if (fracvol == 1.0 && fracvol_top == 0.0) {
      *aux  =  -1;
      *orig = ppb[1];
   } else {
      *aux  =  1;
      *orig = ppt[1];
   }
}

void horizontal_row(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *horizontal, int *aux, real *orig){
   int status;
   real fracvol, fracvol_aux;
   Point pp, ppr, ppl;
   *aux=0;
   *orig=0.0;
   // Horizontal Right
   status = get_frac_vol(sdp, ns, 0, center, p, delta, &fracvol);
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
      status=get_frac_vol(sdp, ns, 0, center, pp, delta, &fracvol);
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
      status = get_frac_vol(sdp, ns, 0, center, pp, delta, &fracvol);
      if (status == 1) {
         *horizontal += fracvol;
      } else if (status == 0) {
         //printf("Left - Horizontal cell out of domain \n");
         return;
      } else {
         //printf("Left - Horizontal cells with different sizes \n");
         return;
      }
   } while (fracvol > 0.0 && fracvol < 1.0 && fabs(fracvol - fracvol_aux)>1.0e-14 && status == 1);

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

real real_max_vec(real *vec,int n){
   real max = vec[1];
   for(int i=1;i<n;i++){
      if(vec[i]>max){
         max=vec[i];
      }
   }
   return max;
}

real real_min_vec(real *vec,int n) {
   real min = vec[1];
   for (int i=1;i<n;i++) {
      if (vec[i]<min) {
         min=vec[i];
      }
   }
   return min;
}

void real_times_vec(real c,real *vec,int n) {
   for (int i=0;i<n;i++) {
      vec[i]=c*vec[i];
   }
}

void set_vec(real orig,real *vec,int n){
   for(int i=0;i<n;i++){
      vec[i] = orig;
   }
}

void vec_minus_vec(real *vec1, real *vec2, real *vec_new, int n){
   for(int i=0;i<n;i++){
      vec_new[i] = vec1[i] - vec2[i];
   }
}

void vec_plus_vec(real *vec1, real *vec2, real *vec_new, int n){
   for(int i=0;i<n;i++){
      vec_new[i] = vec1[i] + vec2[i];
   }
}

void set_common_orig_vertical(real *V,int auxv,real *origv,int n,real dy){
   real orig;
   real cons[n],delta[n];
   if(auxv>0){
      orig=real_max_vec(origv,n);
      set_vec(orig,cons,n);
      vec_minus_vec(cons,origv,delta,n);
      real_times_vec(1/dy,delta,n);
      vec_plus_vec(V,delta,V,n);
   } else {
      orig=real_min_vec(origv,n);
      set_vec(-orig,cons,n);
      vec_plus_vec(cons,origv,delta,n);
      real_times_vec(1/dy,delta,n);
      vec_plus_vec(V,delta,V,n);
   }
}

void set_common_orig_horizontal(real *H,int auxh,real *origh,int n,real dx){
   real orig;
   real cons[n],delta[n];
   if(auxh>0){
      orig=real_max_vec(origh,n);
      set_vec(orig,cons,n);
      vec_minus_vec(cons,origh,delta,n);
      real_times_vec(1/dx,delta,n);
      vec_plus_vec(H,delta,H,n);
   } else {
      orig=real_min_vec(origh,n);
      set_vec(-orig,cons,n);
      vec_plus_vec(cons,origh,delta,n);
      real_times_vec(1/dx,delta,n);
      vec_plus_vec(H,delta,H,n);
   }
}

void calculate_interfacial_force(sim_domain *sdp, higflow_solver *ns, int clid, Point center, Point IF){
   // Set the curvature in the distributed curvature property
   real curvature = compute_value_at_point(sdp, center, center, 1.0, ns->ed.mult.dpcurvature, ns->ed.stn);
   real Normalx = compute_value_at_point(sdp, center, center, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
   real Normaly = compute_value_at_point(sdp, center, center, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
   IF[0] = curvature * Normalx;
   IF[1] = curvature * Normaly;
   //   arquivoCurv(NULL,center[0],center[1],curvature);
   for (int i = 0; i < DIM; i++) {
      dp_set_value(ns->ed.mult.dpIF[i], clid, IF[i]);
   }

}

void higflow_compute_curvature_interfacial_force_normal_multiphase_2D(higflow_solver *ns) {
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
         p[0] = center[0];
         p[1] = center[1];
         real fracvol = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
         fraction_correction_at_get(&fracvol);
         // As funçoes foram alteradas de lugar
         if (fracvol == 0.0 || fracvol == 1.0){
            /* Set the curvature in the distributed curvature property*/
            dp_set_value(ns->ed.mult.dpcurvature, clid, 0.0);
            for (int i = 0; i < DIM; i++) {
               dp_set_value(ns->ed.mult.dpIF[i], clid, 0.0);
               dp_set_value(ns->ed.mult.dpnormal[i], clid, 0.0);
            }
            continue;
         }
         real hm, hb, ht;
         int aux_mh, aux_b, aux_t;
         real orig_mh, orig_b, orig_t;
         // Middle
         horizontal_row(sdp, ns, center, p, delta, &hm, &aux_mh, &orig_mh);
         // Top
         p[1] = center[1] + delta[1];
         horizontal_row(sdp, ns, center, p, delta, &ht, &aux_t, &orig_t);
         // Bottom
         p[1] = center[1] - delta[1];
         horizontal_row(sdp, ns, center, p, delta, &hb, &aux_b, &orig_b);
         if (abs(aux_mh + aux_t + aux_b) == 3) {
            real H[3];
            H[0] = hb;
            H[1] = hm;
            H[2] = ht;
            real origh[3];
            origh[0] = orig_b;
            origh[1] = orig_mh;
            origh[2] = orig_t;
            set_common_orig_horizontal(H, aux_mh, origh, 3, delta[0]);
            calculate_normal_cell_central_2nd_order_finite_difference_Horizontal(ns,clid,H[0],H[1],H[2],delta[0],delta[1],aux_mh);
            calculate_curvature_cell_central_2nd_order_finite_difference_Horizontal(ns,clid,H[0],H[1],H[2],delta[0],delta[1],aux_mh);
            continue;
         }
         // Set p cell
         p[0] = center[0];
         p[1] = center[1];
         real vm, vl, vr;
         int aux_mv, aux_l, aux_r;
         real orig_mv, orig_l, orig_r;
         // Middle
         vertical_collumn(sdp, ns, center, p, delta, &vm, &aux_mv, &orig_mv);
         // Right
         p[0] = center[0] + delta[0];
         vertical_collumn(sdp, ns, center, p, delta, &vr, &aux_r, &orig_r);
         // Left
         p[0] = center[0] - delta[0];
         vertical_collumn(sdp, ns, center, p, delta, &vl, &aux_l, &orig_l);
         if (abs(aux_mv + aux_r + aux_l) == 3){
            real V[3], origv[3];
            V[0] = vl;
            V[1] = vm;
            V[2] = vr;
            origv[0] = orig_l;
            origv[1] = orig_mv;
            origv[2] = orig_r;
            set_common_orig_vertical(V, aux_mv, origv, 3, delta[1]);
            calculate_normal_cell_central_2nd_order_finite_difference_Vertical(ns, clid, V[0], V[1], V[2], delta[0], delta[1], aux_mv);
            calculate_curvature_cell_central_2nd_order_finite_difference_Vertical(ns, clid, V[0], V[1], V[2], delta[0], delta[1], aux_mv);
            continue;
         }
         if (abs(aux_mh+aux_t)==2) {
            // Set p cell
            p[0] = center[0];
            p[1] = center[1];
            real htt;
            int aux_tt;
            real orig_tt;
            //top-top
            p[1] = center[1] + 2*delta[1];
            horizontal_row(sdp, ns, center, p, delta, &htt, &aux_tt, &orig_tt);
            if(abs(aux_mh+aux_t+aux_tt)==3){
               real H[3];
               H[0] = hm;
               H[1] = ht;
               H[2] = htt;
               real origh[3];
               origh[0] = orig_mh;
               origh[1] = orig_t; 
               origh[2] = orig_tt;
               set_common_orig_horizontal(H, aux_mh, origh, 3, delta[0]);
               calculate_normal_cell_progressive_2nd_order_finite_difference_Horizontal(ns, clid, H[0], H[1], H[2], delta[0], delta[1], aux_mh);
               calculate_curvature_cell_progressive_1st_order_finite_difference_Horizontal(ns, clid, H[0], H[1], H[2], delta[0], delta[1], aux_mh);
               continue;
            }
         }
         if(abs(aux_mh+aux_b)==2){
            // Set p cell
            p[0] = center[0];
            p[1] = center[1];
            real hbb;
            int aux_bb;
            real orig_bb;
            //bottom-bottom
            p[1] = center[1] - 2*delta[1];
            horizontal_row(sdp, ns, center, p, delta, &hbb, &aux_bb, &orig_bb);
            if(abs(aux_mh+aux_b+aux_bb)==3){
               real H[3];
               H[0] = hbb;
               H[1] = hb;
               H[2] = hm;
               real origh[3];
               origh[0] = orig_bb;
               origh[1] = orig_b;
               origh[2] = orig_mh;
               set_common_orig_horizontal(H, aux_mh, origh, 3, delta[0]);
               calculate_normal_cell_regressive_2nd_order_finite_difference_Horizontal(ns, clid, H[0], H[1], H[2], delta[0], delta[1], aux_mh);
               calculate_curvature_cell_regressive_1st_order_finite_difference_Horizontal(ns, clid, H[0], H[1], H[2], delta[0], delta[1], aux_mh);
               continue;
            }
         }
         if(abs(aux_mv+aux_r)==2){
            // Set p cell
            p[0] = center[0];
            p[1] = center[1];
            real vrr;
            int aux_rr;
            real orig_rr;
            //right-right
            p[0] = center[0] + 2*delta[0];
            vertical_collumn(sdp, ns, center, p, delta, &vrr, &aux_rr, &orig_rr);
            if(abs(aux_mv+aux_r+aux_rr)==3){
               real V[3], origv[3];
               V[0] = vm;
               V[1] = vr;
               V[2] = vrr;
               origv[0] = orig_mv;
               origv[1] = orig_r; 
               origv[2] = orig_rr;
               set_common_orig_vertical(V, aux_mv, origv, 3, delta[1]);
               calculate_normal_cell_progressive_2nd_order_finite_difference_Vertical(ns, clid, V[0], V[1], V[2], delta[0], delta[1], aux_mv);
               calculate_curvature_cell_progressive_1st_order_finite_difference_Vertical(ns, clid, V[0], V[1], V[2], delta[0], delta[1], aux_mv);
               continue;
            }
         }
         if (abs(aux_mv + aux_l) == 2) {
            // Set p cell
            p[0] = center[0];
            p[1] = center[1];
            real vll;
            int aux_ll;
            real orig_ll;

            //left-left
            p[0] = center[0] - 2 * delta[0];
            vertical_collumn(sdp, ns, center, p, delta, &vll, &aux_ll,&orig_ll);

            if (abs(aux_mv + aux_l + aux_ll) == 3) {
               real V[3];
               V[0] = vll;
               V[1] = vl;
               V[2] = vm;

               real origv[3];
               origv[0] = orig_ll;
               origv[1] = orig_l;
               origv[2] = orig_mv;

               set_common_orig_vertical(V, aux_mv, origv, 3, delta[1]);
               calculate_normal_cell_regressive_2nd_order_finite_difference_Vertical(ns, clid, V[0], V[1], V[2], delta[0], delta[1], aux_mv);
               calculate_curvature_cell_regressive_1st_order_finite_difference_Vertical(ns, clid, V[0], V[1], V[2], delta[0], delta[1], aux_mv);
               continue;
            }
         }
      }
      // Destroy the iterator
      higcit_destroy(it);
      // Sync the distributed pressure property
      dp_sync(ns->ed.mult.dpcurvature);
      for (int i = 0; i < DIM; i++) {
         dp_sync(ns->ed.mult.dpnormal[i]);
      }
}

void higflow_compute_curvature_interfacial_force_normal_multiphase_2D_hf_elvira(higflow_solver *ns) {
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
         p[0] = center[0];
         p[1] = center[1];
         real fracvol = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
         fraction_correction_at_get(&fracvol);
         // As funçoes foram alteradas de lugar
         if (fracvol == 0.0 || fracvol == 1.0){
            dp_set_value(ns->ed.mult.dpcurvature, clid, 0.0);/* Set the curvature in the distributed curvature property*/
            for (int i = 0; i < DIM; i++) {
               dp_set_value(ns->ed.mult.dpIF[i], clid, 0.0);
               dp_set_value(ns->ed.mult.dpnormal[i], clid, 0.0);
            }
            continue;
         }
         // Calculate Normal
         higflow_compute_normal_multiphase_2D_elvira(ns, sdp, mp, it, c, clid, center, delta, p);
         real hm, hb, ht;
         int aux_mh, aux_b, aux_t;
         real orig_mh, orig_b, orig_t;
         // Middle
         horizontal_row(sdp, ns, center, p, delta, &hm, &aux_mh, &orig_mh);
         // Top
         p[1] = center[1] + delta[1];
         horizontal_row(sdp, ns, center, p, delta, &ht, &aux_t, &orig_t);
         // Bottom
         p[1] = center[1] - delta[1];
         horizontal_row(sdp, ns, center, p, delta, &hb, &aux_b, &orig_b);
         if (abs(aux_mh + aux_t + aux_b) == 3) {
            real H[3], origh[3];
            H[0] = hb;
            H[1] = hm;
            H[2] = ht;
            origh[0] = orig_b;
            origh[1] = orig_mh;
            origh[2] = orig_t;
            set_common_orig_horizontal(H, aux_mh, origh, 3, delta[0]);
            calculate_curvature_cell_central_2nd_order_finite_difference_Horizontal(ns,clid,H[0],H[1],H[2],delta[0],delta[1],aux_mh);
            //calculate_interfacial_force(sdp, ns, clid, center, IF);
            continue;
         }
         // Set p cell
         p[0] = center[0];
         p[1] = center[1];
         real vm, vl, vr;
         int aux_mv, aux_l, aux_r;
         real orig_mv, orig_l, orig_r;
         // Middle
         vertical_collumn(sdp, ns, center, p, delta, &vm, &aux_mv, &orig_mv);
         // Right
         p[0] = center[0] + delta[0];
         vertical_collumn(sdp, ns, center, p, delta, &vr, &aux_r, &orig_r);
         // Left
         p[0] = center[0] - delta[0];
         vertical_collumn(sdp, ns, center, p, delta, &vl, &aux_l, &orig_l);
         if (abs(aux_mv + aux_r + aux_l) == 3){
            real V[3];
            V[0] = vl;
            V[1] = vm;
            V[2] = vr;
            real origv[3];
            origv[0] = orig_l;
            origv[1] = orig_mv;
            origv[2] = orig_r;
            set_common_orig_vertical(V, aux_mv, origv, 3, delta[1]);
            calculate_curvature_cell_central_2nd_order_finite_difference_Vertical(ns, clid, V[0], V[1], V[2], delta[0], delta[1], aux_mv);
            //calculate_interfacial_force(sdp, ns, clid, center, IF);
            continue;
         }
         if(abs(aux_mh+aux_t)==2){
            // Set p cell
            p[0] = center[0];
            p[1] = center[1];
            real htt;
            int aux_tt;
            real orig_tt;
            //top-top
            p[1] = center[1] + 2*delta[1];
            horizontal_row(sdp, ns, center, p, delta, &htt, &aux_tt, &orig_tt);
            if(abs(aux_mh+aux_t+aux_tt)==3){
               real H[3], origh[3];
               H[0] = hm;
               H[1] = ht;
               H[2] = htt;
               origh[0] = orig_mh;
               origh[1] = orig_t; 
               origh[2] = orig_tt;
               set_common_orig_horizontal(H, aux_mh, origh, 3, delta[0]);
               calculate_curvature_cell_progressive_1st_order_finite_difference_Horizontal(ns, clid, H[0], H[1], H[2], delta[0], delta[1], aux_mh);
               //calculate_interfacial_force(sdp, ns, clid, center, IF);
               continue;
            }
         }
         if(abs(aux_mh+aux_b)==2){
            // Set p cell
            p[0] = center[0];
            p[1] = center[1];
            real hbb;
            int aux_bb;
            real orig_bb;
            //bottom-bottom
            p[1] = center[1] - 2*delta[1];
            horizontal_row(sdp, ns, center, p, delta, &hbb, &aux_bb, &orig_bb);
            if(abs(aux_mh+aux_b+aux_bb)==3){
               real H[3], origh[3];
               H[0] = hbb;
               H[1] = hb;
               H[2] = hm;
               origh[0] = orig_bb;
               origh[1] = orig_b;
               origh[2] = orig_mh;
               set_common_orig_horizontal(H, aux_mh, origh, 3, delta[0]);
               calculate_curvature_cell_regressive_1st_order_finite_difference_Horizontal(ns, clid, H[0], H[1], H[2], delta[0], delta[1], aux_mh);
               //calculate_interfacial_force(sdp, ns, clid, center, IF);
               continue;
            }
         }
         if(abs(aux_mv+aux_r)==2){
            // Set p cell
            p[0] = center[0];
            p[1] = center[1];
            real vrr, orig_rr;
            int aux_rr;
            //right-right
            p[0] = center[0] + 2*delta[0];
            vertical_collumn(sdp, ns, center, p, delta, &vrr, &aux_rr, &orig_rr);
            if(abs(aux_mv+aux_r+aux_rr)==3){
               real V[3];
               V[0] = vm;
               V[1] = vr;
               V[2] = vrr;
               real origv[3];
               origv[0] = orig_mv;
               origv[1] = orig_r; 
               origv[2] = orig_rr;
               set_common_orig_vertical(V, aux_mv, origv, 3, delta[1]);
               calculate_curvature_cell_progressive_1st_order_finite_difference_Vertical(ns, clid, V[0], V[1], V[2], delta[0], delta[1], aux_mv);
               //calculate_interfacial_force(sdp, ns, clid, center, IF);
               continue;
            }
         }
         if (abs(aux_mv + aux_l) == 2) {
            // Set p cell
            p[0] = center[0];
            p[1] = center[1];
            real vll;
            int aux_ll;
            real orig_ll;
            //left-left
            p[0] = center[0] - 2 * delta[0];
            vertical_collumn(sdp, ns, center, p, delta, &vll, &aux_ll,&orig_ll);
            if (abs(aux_mv + aux_l + aux_ll) == 3) {
               real V[3];
               V[0] = vll;
               V[1] = vl;
               V[2] = vm;
               real origv[3];
               origv[0] = orig_ll;
               origv[1] = orig_l;
               origv[2] = orig_mv;
               set_common_orig_vertical(V, aux_mv, origv, 3, delta[1]);
               calculate_curvature_cell_regressive_1st_order_finite_difference_Vertical(ns, clid, V[0], V[1], V[2], delta[0], delta[1], aux_mv);
               //calculate_interfacial_force(sdp, ns, clid, center, IF);
               continue;
            }
         }
      }
      // Destroy the iterator
      higcit_destroy(it);
      // Sync the distributed pressure property
      dp_sync(ns->ed.mult.dpcurvature);
      for (int i = 0; i < DIM; i++) {
         //dp_sync(ns->ed.mult.dpIF[i]);
         dp_sync(ns->ed.mult.dpnormal[i]);
      }
}

void higflow_compute_curvature_interfacial_force_normal_multiphase_2D_hf_elvira_adap(higflow_solver *ns) {
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
         p[0] = center[0];
         p[1] = center[1];
         real frac = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
         // As funçoes foram alteradas de lugar
         if (frac == 0.0 || frac == 1.0){
            dp_set_value(ns->ed.mult.dpcurvature, clid, 0.0);/* Set the curvature in the distributed curvature property*/
            for (int i = 0; i < DIM; i++) {
               dp_set_value(ns->ed.mult.dpIF[i], clid, 0.0);
               dp_set_value(ns->ed.mult.dpnormal[i], clid, 0.0);
            }
            continue;
         }
          
         higflow_compute_normal_multiphase_2D_elvira_adap(ns, sdp, mp, it, c, clid, center, delta, p);
          
         real hm, hb, ht;
         int aux_mh, aux_b, aux_t;
         real orig_mh, orig_b, orig_t;
         // Middle
         horizontal_row(sdp, ns, center, p, delta, &hm, &aux_mh, &orig_mh);
         // Top
         p[1] = center[1] + delta[1];
         horizontal_row(sdp, ns, center, p, delta, &ht, &aux_t, &orig_t);
         // Bottom
         p[1] = center[1] - delta[1];
         horizontal_row(sdp, ns, center, p, delta, &hb, &aux_b, &orig_b);
         if (abs(aux_mh + aux_t + aux_b) == 3) {
            real H[3], origh[3];
            H[0] = hb;
            H[1] = hm;
            H[2] = ht;
            origh[0] = orig_b;
            origh[1] = orig_mh;
            origh[2] = orig_t;
            set_common_orig_horizontal(H, aux_mh, origh, 3, delta[0]);
            calculate_curvature_cell_central_2nd_order_finite_difference_Horizontal(ns,clid,H[0],H[1],H[2],delta[0],delta[1],aux_mh);
            //calculate_interfacial_force(sdp, ns, clid, center, IF);
            continue;
         }
         // Set p cell
         p[0] = center[0];
         p[1] = center[1];
         real vm, vl, vr;
         int aux_mv, aux_l, aux_r;
         real orig_mv, orig_l, orig_r;
         // Middle
         vertical_collumn(sdp, ns, center, p, delta, &vm, &aux_mv, &orig_mv);
         // Right
         p[0] = center[0] + delta[0];
         vertical_collumn(sdp, ns, center, p, delta, &vr, &aux_r, &orig_r);
         // Left
         p[0] = center[0] - delta[0];
         vertical_collumn(sdp, ns, center, p, delta, &vl, &aux_l, &orig_l);
         if (abs(aux_mv + aux_r + aux_l) == 3){
            real V[3];
            V[0] = vl;
            V[1] = vm;
            V[2] = vr;
            real origv[3];
            origv[0] = orig_l;
            origv[1] = orig_mv;
            origv[2] = orig_r;
            set_common_orig_vertical(V, aux_mv, origv, 3, delta[1]);
            calculate_curvature_cell_central_2nd_order_finite_difference_Vertical(ns, clid, V[0], V[1], V[2], delta[0], delta[1], aux_mv);
            //calculate_interfacial_force(sdp, ns, clid, center, IF);
            continue;
         }
         if(abs(aux_mh+aux_t)==2){
            // Set p cell
            p[0] = center[0];
            p[1] = center[1];
            real htt;
            int aux_tt;
            real orig_tt;
            //top-top
            p[1] = center[1] + 2*delta[1];
            horizontal_row(sdp, ns, center, p, delta, &htt, &aux_tt, &orig_tt);
            if(abs(aux_mh+aux_t+aux_tt)==3){
               real H[3], origh[3];
               H[0] = hm;
               H[1] = ht;
               H[2] = htt;
               origh[0] = orig_mh;
               origh[1] = orig_t; 
               origh[2] = orig_tt;
               set_common_orig_horizontal(H, aux_mh, origh, 3, delta[0]);
               calculate_curvature_cell_progressive_1st_order_finite_difference_Horizontal(ns, clid, H[0], H[1], H[2], delta[0], delta[1], aux_mh);
               //calculate_interfacial_force(sdp, ns, clid, center, IF);
               continue;
            }
         }
         if(abs(aux_mh+aux_b)==2){
            // Set p cell
            p[0] = center[0];
            p[1] = center[1];
            real hbb;
            int aux_bb;
            real orig_bb;
            //bottom-bottom
            p[1] = center[1] - 2*delta[1];
            horizontal_row(sdp, ns, center, p, delta, &hbb, &aux_bb, &orig_bb);
            if(abs(aux_mh+aux_b+aux_bb)==3){
               real H[3], origh[3];
               H[0] = hbb;
               H[1] = hb;
               H[2] = hm;
               origh[0] = orig_bb;
               origh[1] = orig_b;
               origh[2] = orig_mh;
               set_common_orig_horizontal(H, aux_mh, origh, 3, delta[0]);
               calculate_curvature_cell_regressive_1st_order_finite_difference_Horizontal(ns, clid, H[0], H[1], H[2], delta[0], delta[1], aux_mh);
               //calculate_interfacial_force(sdp, ns, clid, center, IF);
               continue;
            }
         }
         if(abs(aux_mv+aux_r)==2){
            // Set p cell
            p[0] = center[0];
            p[1] = center[1];
            real vrr, orig_rr;
            int aux_rr;
            //right-right
            p[0] = center[0] + 2*delta[0];
            vertical_collumn(sdp, ns, center, p, delta, &vrr, &aux_rr, &orig_rr);
            if(abs(aux_mv+aux_r+aux_rr)==3){
               real V[3];
               V[0] = vm;
               V[1] = vr;
               V[2] = vrr;
               real origv[3];
               origv[0] = orig_mv;
               origv[1] = orig_r; 
               origv[2] = orig_rr;
               set_common_orig_vertical(V, aux_mv, origv, 3, delta[1]);
               calculate_curvature_cell_progressive_1st_order_finite_difference_Vertical(ns, clid, V[0], V[1], V[2], delta[0], delta[1], aux_mv);
               //calculate_interfacial_force(sdp, ns, clid, center, IF);
               continue;
            }
         }
         if (abs(aux_mv + aux_l) == 2) {
            // Set p cell
            p[0] = center[0];
            p[1] = center[1];
            real vll;
            int aux_ll;
            real orig_ll;
            //left-left
            p[0] = center[0] - 2 * delta[0];
            vertical_collumn(sdp, ns, center, p, delta, &vll, &aux_ll,&orig_ll);
            if (abs(aux_mv + aux_l + aux_ll) == 3) {
               real V[3];
               V[0] = vll;
               V[1] = vl;
               V[2] = vm;
               real origv[3];
               origv[0] = orig_ll;
               origv[1] = orig_l;
               origv[2] = orig_mv;
               set_common_orig_vertical(V, aux_mv, origv, 3, delta[1]);
               calculate_curvature_cell_regressive_1st_order_finite_difference_Vertical(ns, clid, V[0], V[1], V[2], delta[0], delta[1], aux_mv);
               //calculate_interfacial_force(sdp, ns, clid, center, IF);
               continue;
            }
         }
      }
      // Destroy the iterator
      higcit_destroy(it);
      // Sync the distributed pressure property
      dp_sync(ns->ed.mult.dpcurvature);
      for (int i = 0; i < DIM; i++) {
         //dp_sync(ns->ed.mult.dpIF[i]);
         dp_sync(ns->ed.mult.dpnormal[i]);
      }
}

void higflow_compute_curvature_interfacial_force_normal_multiphase_2D_hf_shirani(higflow_solver *ns) {
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
         p[0] = center[0];
         p[1] = center[1];
         real frac = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
         // As funçoes foram alteradas de lugar
         if (frac == 0.0 || frac == 1.0){
            dp_set_value(ns->ed.mult.dpcurvature, clid, 0.0);/* Set the curvature in the distributed curvature property*/
            for (int i = 0; i < DIM; i++) {
               dp_set_value(ns->ed.mult.dpIF[i], clid, 0.0);
               dp_set_value(ns->ed.mult.dpnormal[i], clid, 0.0);
            }
            continue;
         }
          
         higflow_compute_normal_multiphase_2D_shirani_9_cells(ns, sdp, clid, center, p, delta);
         
         real hm, hb, ht;
         int aux_mh, aux_b, aux_t;
         real orig_mh, orig_b, orig_t;
         // Middle
         horizontal_row(sdp, ns, center, p, delta, &hm, &aux_mh, &orig_mh);
         // Top
         p[1] = center[1] + delta[1];
         horizontal_row(sdp, ns, center, p, delta, &ht, &aux_t, &orig_t);
         // Bottom
         p[1] = center[1] - delta[1];
         horizontal_row(sdp, ns, center, p, delta, &hb, &aux_b, &orig_b);

         if (abs(aux_mh + aux_t + aux_b) == 3){
            real H[3];
            H[0] = hb;
            H[1] = hm;
            H[2] = ht;
            real origh[3];
            origh[0] = orig_b;
            origh[1] = orig_mh;
            origh[2] = orig_t;
            set_common_orig_horizontal(H, aux_mh, origh, 3, delta[0]);
            //calculate_normal_cell_central_2nd_order_finite_difference_Horizontal(ns,clid,H[0],H[1],H[2],delta[0],delta[1],aux_mh);
            calculate_curvature_cell_central_2nd_order_finite_difference_Horizontal(ns,clid,H[0],H[1],H[2],delta[0],delta[1],aux_mh);
            //calculate_interfacial_force(sdp, ns, clid, center, IF);
//            real curvaturech = compute_value_at_point(sdp, center, center,1.0, ns->ed.mult.dpcurvature, ns->ed.stn);
            continue;
         }
//         getchar();
         // Set p cell
         p[0] = center[0];
         p[1] = center[1];

         real vm, vl, vr;
         int aux_mv, aux_l, aux_r;
         real orig_mv, orig_l, orig_r;
         // Middle
         vertical_collumn(sdp, ns, center, p, delta, &vm, &aux_mv, &orig_mv);
         // Right
         p[0] = center[0] + delta[0];
         vertical_collumn(sdp, ns, center, p, delta, &vr, &aux_r, &orig_r);
         // Left
         p[0] = center[0] - delta[0];
         vertical_collumn(sdp, ns, center, p, delta, &vl, &aux_l, &orig_l);
         if (abs(aux_mv + aux_r + aux_l) == 3){
            real V[3];
            V[0] = vl;
            V[1] = vm;
            V[2] = vr;
            real origv[3];
            origv[0] = orig_l;
            origv[1] = orig_mv;
            origv[2] = orig_r;
            set_common_orig_vertical(V, aux_mv, origv, 3, delta[1]);
            //calculate_normal_cell_central_2nd_order_finite_difference_Vertical(ns, clid, V[0], V[1], V[2], delta[0], delta[1], aux_mv);
            calculate_curvature_cell_central_2nd_order_finite_difference_Vertical(ns, clid, V[0], V[1], V[2], delta[0], delta[1], aux_mv);
            //calculate_interfacial_force(sdp, ns, clid, center, IF);
//            real curvaturecv = compute_value_at_point(sdp, center, center,1.0, ns->ed.mult.dpcurvature, ns->ed.stn);
            continue;
         }

         if(abs(aux_mh+aux_t)==2){
            // Set p cell
            p[0] = center[0];
            p[1] = center[1];
            
            real htt;
            int aux_tt;
            real orig_tt;

            //top-top
            p[1] = center[1] + 2*delta[1];
//            printf("Horizontal: going top top\n");
            horizontal_row(sdp, ns, center, p, delta, &htt, &aux_tt, &orig_tt);
//            printf("Horizontal: TOP TOP: htt=%lf auxtt=%d origtt=%lf\n",htt,aux_tt,orig_tt);
            
            if(abs(aux_mh+aux_t+aux_tt)==3){
               real H[3];
               H[0] = hm;
               H[1] = ht;
               H[2] = htt;
               
               real origh[3];
               origh[0] = orig_mh;
               origh[1] = orig_t; 
               origh[2] = orig_tt;
               
               set_common_orig_horizontal(H, aux_mh, origh, 3, delta[0]);
               //calculate_normal_cell_progressive_2nd_order_finite_difference_Horizontal(ns, clid, H[0], H[1], H[2], delta[0], delta[1], aux_mh);
               calculate_curvature_cell_progressive_1st_order_finite_difference_Horizontal(ns, clid, H[0], H[1], H[2], delta[0], delta[1], aux_mh);
               //calculate_interfacial_force(sdp, ns, clid, center, IF);
//               real curvaturecp = compute_value_at_point(sdp, center, center,1.0, ns->ed.mult.dpcurvature, ns->ed.stn);
               continue;
            }
         }

         if(abs(aux_mh+aux_b)==2){
            // Set p cell
            p[0] = center[0];
            p[1] = center[1];
            
            real hbb;
            int aux_bb;
            real orig_bb;
            
            //bottom-bottom
            p[1] = center[1] - 2*delta[1];
//            printf("Horizontal: going botton botton\n");
            horizontal_row(sdp, ns, center, p, delta, &hbb, &aux_bb, &orig_bb);
//            printf("Horizontal: BOTTON BOTTON: hbb=%lf auxbb=%d origbb=%lf\n",hbb,aux_bb,orig_bb);
            
            if(abs(aux_mh+aux_b+aux_bb)==3){
               real H[3];
               H[0] = hbb;
               H[1] = hb;
               H[2] = hm;
               
               real origh[3];
               origh[0] = orig_bb;
               origh[1] = orig_b;
               origh[2] = orig_mh;
               
               set_common_orig_horizontal(H, aux_mh, origh, 3, delta[0]);
               //calculate_normal_cell_regressive_2nd_order_finite_difference_Horizontal(ns, clid, H[0], H[1], H[2], delta[0], delta[1], aux_mh);
               calculate_curvature_cell_regressive_1st_order_finite_difference_Horizontal(ns, clid, H[0], H[1], H[2], delta[0], delta[1], aux_mh);
               //calculate_interfacial_force(sdp, ns, clid, center, IF);
//               real curvaturerh = compute_value_at_point(sdp, center, center, 1.0, ns->ed.mult.dpcurvature, ns->ed.stn);
               continue;
            }
         }

         if(abs(aux_mv+aux_r)==2){
            // Set p cell
            p[0] = center[0];
            p[1] = center[1];
            
            real vrr;
            int aux_rr;
            real orig_rr;
            
            //right-right
            p[0] = center[0] + 2*delta[0];
//            printf("Vertical: going right right\n");
            vertical_collumn(sdp, ns, center, p, delta, &vrr, &aux_rr, &orig_rr);
//            printf("Vertical: RIGHT RIGHT: vrr=%lf auxrr=%d origrr=%lf\n",vrr,aux_rr,orig_rr);
            
            if(abs(aux_mv+aux_r+aux_rr)==3){
               real V[3];
               V[0] = vm;
               V[1] = vr;
               V[2] = vrr;
               
               real origv[3];
               origv[0] = orig_mv;
               origv[1] = orig_r; 
               origv[2] = orig_rr;
               
               set_common_orig_vertical(V, aux_mv, origv, 3, delta[1]);
               //calculate_normal_cell_progressive_2nd_order_finite_difference_Vertical(ns, clid, V[0], V[1], V[2], delta[0], delta[1], aux_mv);
               calculate_curvature_cell_progressive_1st_order_finite_difference_Vertical(ns, clid, V[0], V[1], V[2], delta[0], delta[1], aux_mv);
               //calculate_interfacial_force(sdp, ns, clid, center, IF);
               continue;
            }
         }
         
         if (abs(aux_mv + aux_l) == 2) {
            // Set p cell
            p[0] = center[0];
            p[1] = center[1];

            real vll;
            int aux_ll;
            real orig_ll;

            //left-left
            p[0] = center[0] - 2 * delta[0];
//            printf("Vertical: going left left\n");
            vertical_collumn(sdp, ns, center, p, delta, &vll, &aux_ll,&orig_ll);
//            printf("Vertical: LEFT LEFT: vll=%lf auxll=%d origll=%lf\n",vll,aux_ll,orig_ll);

            if (abs(aux_mv + aux_l + aux_ll) == 3) {
               real V[3];
               V[0] = vll;
               V[1] = vl;
               V[2] = vm;

               real origv[3];
               origv[0] = orig_ll;
               origv[1] = orig_l;
               origv[2] = orig_mv;

               set_common_orig_vertical(V, aux_mv, origv, 3, delta[1]);
               //calculate_normal_cell_regressive_2nd_order_finite_difference_Vertical(ns, clid, V[0], V[1], V[2], delta[0], delta[1], aux_mv);
               calculate_curvature_cell_regressive_1st_order_finite_difference_Vertical(ns, clid, V[0], V[1], V[2], delta[0], delta[1], aux_mv);
               //calculate_interfacial_force(sdp, ns, clid, center, IF);
               continue;
            }
         }
      }
      // Destroy the iterator
      higcit_destroy(it);
      // Sync the distributed pressure property
      dp_sync(ns->ed.mult.dpcurvature);
      for (int i = 0; i < DIM; i++) {
         //dp_sync(ns->ed.mult.dpIF[i]);
         dp_sync(ns->ed.mult.dpnormal[i]);
      }
}
