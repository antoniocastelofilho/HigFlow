#include "hig-flow-vof-plic.h"

void arquivoDist(char*nome_dist,real x,real y,real distance){
   FILE *fp=fopen("dist.txt","a");
   fprintf(fp,"%lf %lf %lf\n",x,y,distance);
   fclose(fp);
}

void arquivoFracAux(char*nome_dist,real x,real y,real frac){
   FILE *fp=fopen("fracaux.txt","a");
   fprintf(fp,"%lf %lf %lf\n",x,y,frac);
   fclose(fp);
}

void solve_equation(int Rx, int Ry, Point Normal, Point Delta, real area, Point Value){real A, B, C;
  
  A = (1.0-Rx-Ry)/fabs(Normal[0]*Normal[1]);
  B = 2.0*(Rx*Delta[0]/fabs(Normal[1])+Ry*Delta[1]/fabs(Normal[0]));
  C = -2.0*area-Rx*Delta[0]*Delta[0]*fabs(Normal[0]/Normal[1])-Ry*Delta[1]*Delta[1]*fabs(Normal[1]/Normal[0]);
  
  if(A==0){
     Value[0] = -C/B;
     Value[1] = -C/B;
     return;
  }
  Value[0] = (-B+sqrt(B*B-4.0*A*C))/(2.0*A);
  Value[1] = (-B-sqrt(B*B-4.0*A*C))/(2.0*A);
  return;
}

real parallel_left_line_origin_center_distance(Point Normal,Point Delta,real area,real tol_n){
   real n_x = fabs(Normal[0]);
   real n_y = fabs(Normal[1]);
   real nx = Normal[0];
   real ny = Normal[1];
   real dx = Delta[0];
   real dy = Delta[1];
   real d;
   
   if (n_x < tol_n) {
      n_x = 0;
      if(ny > 0){
         n_y = 1;
      } else {
         n_y = -1;
      }
   } else if (n_y < tol_n) {
      n_y = 0;
      if(nx > 0){
         n_x = 1;
      } else {
         n_x = -1;
      }
   }
   
   if (n_x == -1) {
      d = -area/dy;
   } else if (n_x == 1) {
      d = -(area/dy - dx);
   } else if (n_y == -1) {
      d = -area/dx;
   } else if (n_y == 1) {
      d = -(area/dx - dy);
   }
   return d;
}

real distance_from_center(Point Normal,Point Delta,real AREA){
   real nx       = Normal[0];
   real ny       = Normal[1];
   real dx       = Delta[0];
   real dy       = Delta[1];
   Point Value;
   real area;
   real d;

   if (FLT_EQ(nx, 0.0) || FLT_EQ(ny, 0.0)) {
      d = parallel_left_line_origin_center_distance(Normal, Delta, AREA, EPSMACH);
      d = trans_bl_2_center(Delta, Normal, d);
      return d;
   } else if (nx * ny > 0) {
      real dLT = nx * (-0.5 * dx) + ny * (0.5 * dy);
      real aLT = area_left_line_origin_center(Normal, Delta, dLT);

      real dRB = nx * (0.5 * dx) + ny * (-0.5 * dy);
      real aRB = area_left_line_origin_center(Normal, Delta, dRB);

      area = AREA;
      if (FLT_EQ(area, aLT)) {
         d = dLT;
         return d;
      } else if (FLT_EQ(area, aRB)) {
         d = dRB;
         return d;
      }

      if (nx > 0) {
         area = dx * dy - AREA;
         aLT = dx * dy - aLT;
         aRB = dx * dy - aRB;
      }

      real amax = fmax(aLT, aRB);
      real amin = fmin(aLT, aRB);
      Point d_vec;
      int Rx, Ry;
      if (area > amax) {
         Rx = 1;
         Ry = 1;
         solve_equation(Rx, Ry, Normal, Delta, area, d_vec);
         d = fmin(d_vec[0], d_vec[1]);
      } else if (area < amax && area > amin) {
         if (aRB > aLT) {
            Rx = 0;
            Ry = 1;
         } else {
            Rx = 1;
            Ry = 0;
         }
         solve_equation(Rx, Ry, Normal, Delta, area, d_vec);
         d = d_vec[0];
      } else if (area < amin) {
         Rx = 0;
         Ry = 0;
         solve_equation(Rx, Ry, Normal, Delta, area, d_vec);
         d = fmax(d_vec[0], d_vec[1]);
      }
      d = sign(nx) * d;
      d = trans_bl_2_center(Delta, Normal, d);
      return d;
   } else {
      real dLB = nx * (0.5 * dx) + ny * (0.5 * dy);
      real aLB = area_left_line_origin_center(Normal, Delta, dLB);

      real dRT = nx * (-0.5 * dx) + ny * (-0.5 * dy);
      real aRT = area_left_line_origin_center(Normal, Delta, dRT);

      area = AREA;
      if (FLT_EQ(area, aRT)) {
         d = dRT;
         return d;
      } else if (FLT_EQ(area, aLB)) {
         d = dLB;
         return d;
      }

      if (nx < 0) {
         area = dx * dy - AREA;
         aRT = dx * dy - aRT;
         aLB = dx * dy - aLB;
      }

      real amax = fmax(aRT, aLB);
      real amin = fmin(aRT, aLB);
      int Rx, Ry;
      Point d_vec;
      if (area > amax) {
         Rx = 1;
         Ry = 1;
         solve_equation(Rx, Ry, Normal, Delta, area, d_vec);
         d = fmin(d_vec[0], d_vec[1]);
      } else if (area > amin && area < amax) {
         if (aRT > aLB) {
            Rx = 0;
            Ry = 1;
         } else {
            Rx = 1;
            Ry = 0;
         }
         solve_equation(Rx, Ry, Normal, Delta, area, d_vec);
         d = d_vec[0];
      } else if (area < amin) {
         Rx = 0;
         Ry = 0;
         solve_equation(Rx, Ry, Normal, Delta, area, d_vec);
         d = fmax(d_vec[0], d_vec[1]);
      }
      d = -sign(nx) * d;
      d = trans_br_2_center(Delta, Normal, d);
      return d;
   }
}

real parallel_left_line_origin_center_area(Point Normal,Point Delta,real d_from_center,real tol_n){
   real nx = Normal[0];
   real ny = Normal[1];
   real dx = Delta[0];
   real dy = Delta[1];
   int n_x, n_y;
   real area, d;
   
   if (fabs(nx) < tol_n) {
      n_x = 0;
      if(ny > 0) {
         n_y = 1;
      } else {
         n_y = -1;
      }
   } else if (fabs(ny) < tol_n) {
      n_y = 0;
      if (nx > 0) {
         n_x = 1;
      } else {
         n_x = -1;
      }
   }
   
   d = trans_center_2_bl(Delta,Normal,d_from_center);
   
   if (n_x==-1) {
      area=-d*dy;
   } else if (n_x==1) {
      area=dy*(dx-d);
   } else if (n_y==-1) {
      area=-d*dx;
   } else if (n_y==1) {
      area=dx*(dy-d);
   }
   return area;
}

real area_left_line_origin_center(Point Normal,Point Delta,real d_from_center){

   real d, area;
   real dx  = Delta[0];
   real dy  = Delta[1];
   real nx  = Normal[0];
   real ny  = Normal[1];

   real n_x   = fabs(Normal[0]);
   real n_y   = fabs(Normal[1]);
   
   Point Prt, Prb, Plt, Plb;
   Prt[0] = 0.5*dx;Prt[1]  = 0.5*dy;
   Prb[0] = 0.5*dx;Prb[1]  = -0.5*dy;
   Plt[0] = -0.5*dx;Plt[1] = 0.5*dy;
   Plb[0] = -0.5*dx;Plb[1] = -0.5*dy;
   
   int srt = left_right(Prt,Normal,d_from_center);
   int srb = left_right(Prb,Normal,d_from_center);
   int slt = left_right(Plt,Normal,d_from_center);
   int slb = left_right(Plb,Normal,d_from_center);

   if (srt <= 0 && srb <= 0 && slt <= 0 && slb <= 0) {
      return area = dx*dy;
   } else if (srt >= 0 && srb >= 0 && slt >= 0 && slb >= 0) {
      return area = 0.0;
   }
   
   if (FLT_EQ(n_x, 0.0) || FLT_EQ(n_y, 0.0)) {
      return area = parallel_left_line_origin_center_area(Normal,Delta,d_from_center,EPSMACH);
   }
   
   if (nx*ny>0) {
      d = trans_center_2_bl(Delta,Normal,d_from_center);
   } else {
      d = trans_center_2_br(Delta,Normal,d_from_center);
   }
   
   d = fabs(d);
   
   real aux1   = 1/(2.0*n_x*n_y);
   real auxx   = d - n_x*dx;
   real auxx_1 = auxx*auxx;
   real auxy   = d - n_y*dy;
   real auxy_1 = auxy*auxy;
   
   area   = aux1*(d*d - R(auxx)*auxx_1 - R(auxy)*auxy_1);
   
   if (ny>0) {
      area=dx*dy-area;
   }
   return area;
}

int left_right(Point P,Point Normal,real d) {
   int side;
   return side = sign(-Normal[0]*P[0] - Normal[1]*P[1] + d);
}

real trans_center_2_bl(Point Delta,Point Normal,real d_from_center){
   real d_from_bl = d_from_center + Normal[0]*0.5*Delta[0] + Normal[1]*0.5*Delta[1];
   return d_from_bl;
}

real trans_bl_2_center(Point Delta, Point Normal, real d_from_bl) {
   real d_from_center = d_from_bl - Normal[0] * 0.5 * Delta[0] - Normal[1] * 0.5 * Delta[1];
   return d_from_center;
}

real trans_center_2_br(Point Delta,Point Normal,real d_from_center){
   real d_from_br = d_from_center - Normal[0]*0.5*Delta[0] + Normal[1]*0.5*Delta[1];
   return d_from_br;
}

real trans_br_2_center(Point Delta,Point Normal,real d_from_br){
   real d_from_center = d_from_br + Normal[0]*0.5*Delta[0] - Normal[1]*0.5*Delta[1];
   return d_from_center;
}

real R(real x){
   real value;
   if (x <= 0.0) {
      value = 0.0;
   } else {
      value = 1.0;
   }
   return value;
}

int sign(real value) {
   if (value > 0.0) {
      return 1;
   } else if (value < 0.0) {
      return -1;
   } else {
      return 0;
   }
}

void higflow_compute_distance_multiphase_2D(higflow_solver *ns) {
    real IF[DIM];
    // Get the local sub-domain for the cells
    sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
    // Get the map for the domain properties
    mp_mapper *mp = sd_get_domain_mapper(sdm);
    // Loop for each cell
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(sdm); !higcit_isfinished(it); higcit_nextcell(it)) {
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
       Point Normal;
       Normal[0] = compute_value_at_point(sdm, center, center, 1.0, ns->ed.mult.dpnormal[0], ns->ed.mult.stn);
       Normal[1] = compute_value_at_point(sdm, center, center, 1.0, ns->ed.mult.dpnormal[1], ns->ed.mult.stn);
       // If the cell is not interfacial, the distance is not calculated 
       if (FLT_EQ(Normal[0], 0.0) && FLT_EQ(Normal[1], 0.0)) {
          continue;
       }
       real fracvol = compute_value_at_point(sdm, center, center, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
       real area = fracvol*delta[0]*delta[1];
       // Calculate distance from center
       real distance = distance_from_center(Normal,delta,area);
       // Set value distance
       dp_set_value(ns->ed.mult.dpdistance, clid, distance);
    }
    // Destroy the iterator
    higcit_destroy(it);
    // Sync the distributed pressure property
    dp_sync(ns->ed.mult.dpdistance);
}

// compute the plic lines of the interface for the 2D case
void higflow_compute_plic_lines_2d(higflow_solver *ns) {
   sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
   mp_mapper *mp = sd_get_domain_mapper(sdm);
   higcit_celliterator *it; hig_cell *c;
   int clid;
   real fracvol;
   int interface_counter = 0;
   for (it = sd_get_domain_celliterator(sdm); !higcit_isfinished(it); higcit_nextcell(it)) {
      c = higcit_getcell(it);
      clid = mp_lookup(mp, hig_get_cid(c));
      fracvol = dp_get_value(ns->ed.mult.dpfracvol, clid);
      if(FLT_GT(fracvol, 0.0) && FLT_LT(fracvol, 1.0)) interface_counter++;
   } higcit_destroy(it);

   if(ns->ed.mult.num_plic_lines !=0) free(ns->ed.mult.plic_lines);
   ns->ed.mult.num_plic_lines = interface_counter;

   if(interface_counter == 0) return;
   ns->ed.mult.plic_lines = (Line*)malloc(sizeof(Line)*interface_counter);

   int line_index = 0;
   Point ccenter, cdelta;
   real nx, ny, dist;
   Line l;
   for (it = sd_get_domain_celliterator(sdm); !higcit_isfinished(it); higcit_nextcell(it)) {
      c = higcit_getcell(it);
      clid = mp_lookup(mp, hig_get_cid(c));
      fracvol = dp_get_value(ns->ed.mult.dpfracvol, clid);
      if(FLT_GT(fracvol, 0.0) && FLT_LT(fracvol, 1.0)) {
         hig_get_center(c, ccenter);
         hig_get_delta(c, cdelta);
         nx= dp_get_value(ns->ed.mult.dpnormal[0], clid);
         ny = dp_get_value(ns->ed.mult.dpnormal[1], clid);
         dist = dp_get_value(ns->ed.mult.dpdistance, clid);

         if(FLT_EQ(nx, 0.0)) {
            l[0][0] = ccenter[0] - 0.5*cdelta[0];
            l[0][1] = ccenter[1] + dist;
            l[1][0] = ccenter[0] + 0.5*cdelta[0];
            l[1][1] = ccenter[1] + dist;
         } else if(FLT_EQ(ny, 0.0)) {
            l[0][0] = ccenter[0] + dist;
            l[0][1] = ccenter[1] - 0.5*cdelta[1];
            l[1][0] = ccenter[0] + dist;
            l[1][1] = ccenter[1] + 0.5*cdelta[1];
         } else {
            int pt_index = 0;
            real x,y;
            x = -0.5*cdelta[0]; y = (dist - x*nx)/ny;
            if(FLT_LE(fabs(y),0.5*cdelta[1])) { // if it intersects the left wall
               l[pt_index][0] = x + ccenter[0]; l[pt_index][1] = y + ccenter[1]; 
               pt_index++;
            }
            x = 0.5*cdelta[0]; y = (dist - x*nx)/ny;
            if(FLT_LE(fabs(y),0.5*cdelta[1])) { // if it intersects the right wall
               l[pt_index][0] = x + ccenter[0]; l[pt_index][1] = y + ccenter[1]; 
               pt_index++;
            }
            y = -0.5*cdelta[1]; x = (dist - y*ny)/nx;
            if(FLT_LE(fabs(x),0.5*cdelta[0])) { // if it intersects the bottom wall
               l[pt_index][0] = x + ccenter[0]; l[pt_index][1] = y + ccenter[1]; 
               pt_index++;
            }
            y = 0.5*cdelta[1]; x = (dist - y*ny)/nx;
            if(FLT_LE(fabs(x),0.5*cdelta[0])) { // if it intersects the top wall
               l[pt_index][0] = x + ccenter[0]; l[pt_index][1] = y + ccenter[1]; 
               pt_index++;
            }
         }

         // check orientation
         real vx = l[1][0] - l[0][0]; real vy = l[1][1] - l[0][1];
         if(FLT_GE(nx*vy - ny*vx, 0.0)) {
            ns->ed.mult.plic_lines[line_index][0][0] = l[0][0];
            ns->ed.mult.plic_lines[line_index][0][1] = l[0][1];
            ns->ed.mult.plic_lines[line_index][1][0] = l[1][0];
            ns->ed.mult.plic_lines[line_index][1][1] = l[1][1];
         } else {
            ns->ed.mult.plic_lines[line_index][0][0] = l[1][0];
            ns->ed.mult.plic_lines[line_index][0][1] = l[1][1];
            ns->ed.mult.plic_lines[line_index][1][0] = l[0][0];
            ns->ed.mult.plic_lines[line_index][1][1] = l[0][1];
         }
         line_index++;
      }
      
   } higcit_destroy(it);
}

void higflow_compute_area_fraction_multiphase_2D(higflow_solver *ns) {
      real IF[DIM];
      // Get the local sub-domain for the cells
      sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);

      // Get the map for the domain properties
      mp_mapper *mp = sd_get_domain_mapper(sdm);

      // Loop for each cell
      higcit_celliterator *it;

      for (it = sd_get_domain_celliterator(sdm); !higcit_isfinished(it); higcit_nextcell(it)) {
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
         Point Normal;
         Normal[0] = compute_value_at_point(sdm, center, center, 1.0, ns->ed.mult.dpnormal[0], ns->ed.mult.stn);
         Normal[1] = compute_value_at_point(sdm, center, center, 1.0, ns->ed.mult.dpnormal[1], ns->ed.mult.stn);

         if (FLT_EQ(Normal[0], 0.0) && FLT_EQ(Normal[1], 0.0)){
            continue;
         }

         real d_from_center = compute_value_at_point(sdm, center, center, 1.0, ns->ed.mult.dpdistance, ns->ed.mult.stn);
         real area = area_left_line_origin_center(Normal, delta, d_from_center);
         real frac = area/(delta[0]*delta[1]);

         // arquivoFracAux(NULL,center[0],center[1],frac);
      }
      // Destroy the iterator
      higcit_destroy(it);
}

