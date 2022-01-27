#include "hig-flow-vof-9-cells.h"
#include "hig-flow-vof-finite-difference-normal-curvature.h"

void higflow_compute_normal_multiphase_2D_shirani_9_cells(higflow_solver *ns, sim_domain *sdp, mp_mapper *mp, higcit_celliterator *it, hig_cell *c,int clid, Point center, Point delta) {
   Point p, pp, Normal;
   real fracvol, gradFij1, gradFij2, gradFij3, gradFij, gradFik1, gradFik2;
   real gradFik3, gradFik, gradFji1, gradFji2, gradFji3, gradFji, gradFjk1;
   real gradFjk2, gradFjk3, gradFjk, gradFki1, gradFki2, gradFki3, gradFki;
   real gradFkj1, gradFkj2, gradFkj3, gradFkj, gradFi, gradFj, gradFk;

   p[0] = center[0]; pp[0]=p[0];
   p[1] = center[1]; pp[1]=p[1];
   
   real f_c = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
   
   pp[0]=p[0]+delta[0];
   pp[1]=p[1];
   real f_r = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
   
   pp[0]=p[0]-delta[0];
   pp[1]=p[1];
   real f_l = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
   
   pp[0]=p[0];
   pp[1]=p[1]+delta[1];
   real f_u = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
   
   pp[0]=p[0]+delta[0];
   pp[1]=p[1]+delta[1];
   real f_ru = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
   
   pp[0]=p[0]-delta[0];
   pp[1]=p[1]+delta[1];
   real f_lu = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
   
   pp[0]=p[0];
   pp[1]=p[1]-delta[1];
   real f_d = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
   
   pp[0]=p[0]+delta[0];
   pp[1]=p[1]-delta[1];
   real f_rd = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
   
   pp[0]=p[0]-delta[0];
   pp[1]=p[1]-delta[1];
   real f_ld = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
   
   //at x direction
   //==================================================================
   //plane i-j
   gradFi = (f_rd-f_ld + 2.0*(f_r-f_l) + f_ru-f_lu)/(delta[0]);
       
    //at y direction
   //==================================================================
   //plane j-i
   gradFj = (f_lu-f_ld + 2.0*(f_u-f_d) + f_ru-f_rd)/delta[1];
   
   //Normal
   //=================================================================
   real norm_grad = sqrt(pow(gradFi,2)+pow(gradFj,2));
    Normal[0] = gradFi/(norm_grad+1.0e-16);
    Normal[1] = gradFj/(norm_grad+1.0e-16);
    
    //printf("%lf %lf\n", Normal[0], Normal[1]);
   
   for (int i=0; i<DIM; i++){
      dp_set_value(ns->ed.mult.dpnormal[i], clid, Normal[i]);
   }
   
}
