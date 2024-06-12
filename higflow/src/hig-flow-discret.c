// *******************************************************************
//  HiG-Flow Solver Discretization - version 10/11/2016
// *******************************************************************

#include "hig-flow-discret.h"

// *******************************************************************
// Navier-Stokes Discretization
// *******************************************************************

// Compute the first derivative using the central method
real compute_facet_dudxc(Point delta, int dim, real alpha, real uc, real ul, real ur) {
   // Midpoint first derivative
   real dudx = (ur - ul) / (2.0 * alpha * delta[dim]);
   return dudx;
}

// Compute the first derivative using the forward method
real compute_facet_dudxr(Point delta, int dim, real alpha, real uc, real ul, real ur) {
   // Midpoint first derivative
   real dudx = (ur - uc) / (alpha * delta[dim]);
   return dudx;
}

// Compute the first derivative using the backward method
real compute_facet_dudxl(Point delta, int dim, real alpha, real uc, real ul, real ur) {
   // Midpoint first derivative
   real dudx = (uc - ul) / (alpha * delta[dim]);
   return dudx;
}

// Compute the second derivative using the forward method
real compute_facet_du2dx2r(Point delta, int dim, real alpha, real u, real ur, real urr) {
   // Midpoint second derivative
   real du2dx2 = (urr - 2.0*ur + u)/(delta[dim] * delta[dim] * alpha * alpha);
   return du2dx2;
}

// Compute the second derivative using the backward method
real compute_facet_du2dx2l(Point delta, int dim, real alpha, real u, real ul, real ull) {
   // Midpoint second derivative
   real du2dx2 = (u - 2.0*ul + ull)/(delta[dim] * delta[dim] * alpha * alpha);
   return du2dx2;
}

// Compute the second derivative using the central method
real compute_facet_du2dx2(Point delta, int dim, real alpha, real u, real ul, real ur) {
   // Midpoint second derivative
   real du2dx2 = (ul - 2.0*u + ur)/(delta[dim] * delta[dim] * alpha * alpha);
   return du2dx2;
}

// Compute the first derivative using the central method
real compute_dpdx_at_point(Point delta, int dim, real alpha, real valuel, real valueh) {
   // Midpoint first derivative
   real dpdx = (valueh - valuel) / (2.0 * alpha * delta[dim]);
   return dpdx;
}

// Compute the first derivative using the forward method
real compute_dpdxr_at_point(Point delta, int dim, real alpha, real valuec, real valuer) {
   // Midpoint first derivative
   real dpdx = (valuer - valuec) / (alpha * delta[dim]);
   return dpdx;
}

// Compute the first derivative using the backward method
real compute_dpdxl_at_point(Point delta, int dim, real alpha, real valuel, real valuec) {
   // Midpoint first derivative
   real dpdx = (valuec - valuel) / (alpha * delta[dim]);
   return dpdx;
}

// Compute the mid point value
real compute_value_at_mid_point(real valuel, real valuer) {
   // Midpoint value
   real value = (valuer + valuel) / 2.0;
   return value;
}

// *******************************************************************
// Navier-Stokes Computational Cell
// *******************************************************************

// Computing the necessary term for the Navier-Stokes equation
void higflow_computational_cell(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int flid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]) {
   // Set the computational cell
   real  uc[DIM], ul[DIM], ur[DIM], Sl[DIM], Sr[DIM], pl, pr, Fl, Fr;
   real  ull[DIM], urr[DIM];
   real  a, b, c, d, e, tol, fi, conv1,conv2;
   a     = 1.7500;
   b     = 0.3750;
   c     = 0.7500;
   d     = 0.1250;
   e     = 0.2500;
   conv1 = 0.0;
   conv2 = 0.0;
   int   infacet, infacet_r, infacet_l, infacet_rr, infacet_ll;
   // Spatial discretization
   switch (ns->contr.spatialdiscrtype) {
   // Second order
   case ORDER2:
      // Set the convective method
      ns->cc.convec_type = ns->contr.convecdiscrtype;
      // Get the velocity in the facet center
      ns->cc.ufacet  = dp_get_value(dpu[dim], flid);
      // Get the facet source tern in the facet center
      ns->cc.F      = compute_facet_value_at_point(ns->sfdF[dim], fcenter, fcenter, 1.0, ns->dpFU[dim], ns->stnF);
      // Get the pressure in the left cell
      pl            = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
      // Get the pressure in the right cell
      pr            = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
      // Compute the pressure derivative
      ns->cc.dpdx   = compute_dpdx_at_point(fdelta, dim, 0.5, pl, pr);
      // Get the cell source term in the left cell
      Fl            = compute_center_p_left(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
      // Get the cell source term in the left cell
      Fr            = compute_center_p_right(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
      // Compute the source term at the mid point
      ns->cc.F     += compute_value_at_mid_point(Fl, Fr);
      for (int dim2 = 0; dim2 < DIM; dim2++) {
         // Get the velocity in the center of facet
         if (dim2 == dim) {
            ns->cc.v[dim2] = ns->cc.ufacet;
         } else {
            ns->cc.v[dim2] = compute_facet_u_2(sfdu[dim2], fcenter, fdelta, dim, dim2, 1.0, dpu[dim2], ns->stn);
         }
         // Get the velocity in the left facet center
         ul[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet_l);
         //Changes to first ordem upwind if there is a boundary cell
         if ((infacet_l == 0) && (ns->cc.convec_type == SECOND_ORDER)) {
            ns->cc.convec_type = FIRST_ORDER;
         }
          // Get the velocity in the right facet center
         ur[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet_r);
         //Changes to first ordem upwind if there is a boundary cell
         if ((infacet_r == 0) && (ns->cc.convec_type == SECOND_ORDER)) {
            ns->cc.convec_type = FIRST_ORDER;
         }
         // Compute dudxc at facet center
         ns->cc.dudxc[dim2]  = compute_facet_dudxc(fdelta, dim2, 1.0, ns->cc.ufacet, ul[dim2], ur[dim2]);
         // Compute dudxr at facet center
         ns->cc.dudxr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ns->cc.ufacet, ul[dim2], ur[dim2]);
         // Compute dudxl at facet center
         ns->cc.dudxl[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ns->cc.ufacet, ul[dim2], ur[dim2]);
         
         // Get the velocity in the right facet center
         urr[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_rr);
         // Get the velocity in the left facet center
         ull[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_ll);
         if (ns->cc.convec_type == SECOND_ORDER) {
            if (ns->contr.secondconvecdiscrtype == QUICK) {
               // Compute dudxrr at facet center
               ns->cc.dudxrr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ur[dim2], ns->cc.ufacet, urr[dim2]);
               // Compute dudxll at facet center
               ns->cc.dudxll[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ul[dim2], ull[dim2], ns->cc.ufacet);
               // Compute terms for second order upwind CUBISTA
            }else if (ns->contr.secondconvecdiscrtype == CUBISTA){
               // Get the velocity in the left facet center
               ul[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet_l);
               // Get the velocity in the left facet center
               ull[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_ll);
               //Changes to central scheme if there is a outside boundary cell
               if (infacet_l == 0)  {
                  ns->cc.convec_type = FIRST_ORDER;
               }
               // Get the velocity in the right facet center
               ur[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet_r);
               // Get the velocity in the right facet center
               urr[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_rr);
               //Changes to central scheme if there is a outside boundary cell
               if (infacet_r == 0)  {
                  ns->cc.convec_type = FIRST_ORDER;
               }
               // Get the velocity v1bar(i+1/2,j+1/2) in the facet center
               ns->cc.vc[dim2]  = compute_facet_u_right(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
               if (ns->cc.vc[dim2] > 0.0){
                  if (FLT_EQ(ur[dim], ul[dim])){
                     conv1 = ns->cc.vc[dim2]*ns->cc.ufacet;
                  }else {
                     fi = (ns->cc.ufacet - ul[dim])/(ur[dim] - ul[dim]);
                     if ((fi <= 0.0) || (fi >= 1.0)) {
                        conv1 = ns->cc.vc[dim2]*ns->cc.ufacet;
                     }else {
                        if (fi < b)
                           conv1 = ns->cc.vc[dim2]*(a*ns->cc.ufacet - c*ul[dim]);
                        if ((fi >= b) && (fi <= c))
                           conv1 = ns->cc.vc[dim2]*(c*ns->cc.ufacet + b*ur[dim] -d*ul[dim]);
                        if (fi > c)
                           conv1 = ns->cc.vc[dim2]*(e*ns->cc.ufacet + c*ur[dim]);
                     }
                  }
                  //v1bar < 0.0
               }else {
                  if (FLT_EQ(ns->cc.ufacet, urr[dim])){
                     conv1 = ns->cc.vc[dim2]*ur[dim];
                  } else {
                     fi = (ur[dim]-urr[dim])/(ns->cc.ufacet -urr[dim]);
                     if ((fi <= 0.0) || (fi >= 1.0)) {
                        conv1 = ns->cc.vc[dim2]*ur[dim];
                     }else {
                        if (fi < b)
                           if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(a*ur[dim] - c*urr[dim]);
                           else                 conv1 = ns->cc.vc[dim2]*ur[dim];
                        if ((fi >= b) && (fi <= c))
                           if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(c*ur[dim] + b*ns->cc.ufacet -d*urr[dim]);
                           else                 conv1 = ns->cc.vc[dim2]*ur[dim];
                        if (fi > c)
                           conv1 = ns->cc.vc[dim2]*(c*ns->cc.ufacet + e*ur[dim]);
                     }
                  }
               }
               // Get the velocity  v2bar(i+1/2,j-1/2) in the facet center
               ns->cc.vc[dim2]  = compute_facet_u_left(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
               if (ns->cc.vc[dim2] > 0.0){
                  if (FLT_EQ(ns->cc.ufacet, ull[dim])) {
                     conv2 = ns->cc.vc[dim2]*ul[dim];
                  }else {
                     fi = (ul[dim] - ull[dim])/(ns->cc.ufacet-ull[dim]);
                     if ((fi <= 0.0) || (fi >= 1.0)) {
                        conv2 = ns->cc.vc[dim2]*ul[dim];
                     }else {
                        if (fi < b)
                           if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(a*ul[dim] - c*ull[dim]);
                           else                  conv2 = ns->cc.vc[dim2]*ul[dim];
                        if ((fi >= b) && (fi <= c))
                           if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(b*ns->cc.ufacet + c*ul[dim] - d*ull[dim]);
                           else                  conv2 = ns->cc.vc[dim2]*ul[dim];
                        if (fi > c)
                           conv2 = ns->cc.vc[dim2]*(c*ns->cc.ufacet + e*ul[dim]);
                     }
                  }
               }else {
                  //v2bar < 0.0
                  if (FLT_EQ(ul[dim], ur[dim])) {
                     conv2 = ns->cc.vc[dim2]*ns->cc.ufacet;
                  }else {
                     fi = (ns->cc.ufacet - ur[dim])/(ul[dim] - ur[dim]);
                     if ((fi <= 0.0) || (fi >= 1.0)) {
                        conv2 = ns->cc.vc[dim2]*ns->cc.ufacet;
                     }else {
                        if (fi < b)
                           conv2 = ns->cc.vc[dim2]*(a*ns->cc.ufacet - c*ur[dim]);
                        if ((fi >= b) && (fi <= c))
                           conv2 = ns->cc.vc[dim2]*(b*ul[dim] + c*ns->cc.ufacet - d*ur[dim]);
                        if (fi > c)
                           conv2 = ns->cc.vc[dim2]*(c*ul[dim] + e*ns->cc.ufacet);
                     }
                  }
               }
               ns->cc.vc[dim2] = ((conv1-conv2)/fdelta[dim]);
            }
         }
         // Newtonian
         ns->cc.du2dx2[dim2] = compute_facet_du2dx2(fdelta, dim2, 1.0, ns->cc.ufacet, ul[dim2], ur[dim2]);
      }
      break;
   }
}


// Computing the necessary term for the Generalized Newtonian Navier-Stokes equation
void higflow_computational_cell_gen_newt(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int flid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]) {
   // Set the computational cell
   higflow_computational_cell(ns, sdp, sfdu, flid, fcenter, fdelta, dim, dpu);
   int infacet_l, infacet_r;
   // Spatial discretization
   switch (ns->contr.spatialdiscrtype) {
   // Second order
   case ORDER2:
      // Get the cell viscolity in the left cell
      ns->cc.viscl = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.gn.dpvisc, ns->ed.stn);
      // Get the cell viscolity in the right cell
      ns->cc.viscr = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.gn.dpvisc, ns->ed.stn);
      for (int dim2 = 0; dim2 < DIM; dim2++) {
         // Generalized Newtonian
         ns->cc.du2dx2[dim2] = 0.0;
         Point p;
         POINT_ASSIGN(p, fcenter);

         p[dim2]      = fcenter[dim2] + 0.5*fdelta[dim2];

         real ucl     = compute_facet_u_left(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet_l);
         real ucr     = compute_facet_u_right(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet_r);
         real duidxjr = compute_facet_dudxc(fdelta, dim2, 0.5, ns->cc.ufacet, ucl, ucr);

         real vcl     = compute_facet_u_left(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet_l);
         real vcr     = compute_facet_u_right(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet_r);
         real dujdxir = compute_facet_dudxc(fdelta, dim, 0.5, ns->cc.ufacet, vcl, vcr);

         p[dim2]      =  fcenter[dim2] - 0.5*fdelta[dim2];

         ucl          = compute_facet_u_left(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet_l);
         ucr          = compute_facet_u_right(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet_r);
         real duidxjl = compute_facet_dudxc(fdelta, dim2, 0.5, ns->cc.ufacet, ucl, ucr);

         vcl          = compute_facet_u_left(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet_l);
         vcr          = compute_facet_u_right(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet_r);
         real dujdxil = compute_facet_dudxc(fdelta, dim, 0.5, ns->cc.ufacet, vcl, vcr);

         ns->cc.du2dx2[dim2] += (ns->cc.viscr*duidxjr - ns->cc.viscl*duidxjl)/fdelta[dim2];
         ns->cc.du2dx2[dim2] += (ns->cc.viscr*dujdxir - ns->cc.viscl*dujdxil)/fdelta[dim2];
      }
      break;
   }
}

// Computing the necessary term for the Multifase Navier-Stokes equation
void higflow_computational_cell_multiphase(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int flid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]) {
   // Set the computational cell
   higflow_computational_cell(ns, sdp, sfdu, flid, fcenter, fdelta, dim, dpu);
   real  Sl[DIM], Sr[DIM];
   int infacet_l, infacet_r;
   real  rhol, rhor, curvl, curvr;
   // Spatial discretization
   switch (ns->contr.spatialdiscrtype) {
   // Second order
   case ORDER2:
      // Get the cell curvature in the left cell
      curvl = compute_center_p_left(ns->ed.mult.sdmult, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpcurvature, ns->ed.mult.stn);
      // Get the cell curvature in the right cell
      curvr = compute_center_p_right(ns->ed.mult.sdmult, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpcurvature, ns->ed.mult.stn);
      //ns->cc.curv = compute_value_at_mid_point(curvl, curvr);
      // Get the cell density in the left cell
      rhol          = compute_center_p_left(ns->ed.mult.sdmult, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpdens, ns->ed.mult.stn);
      // Get the cell density in the right cell
      rhor          = compute_center_p_right(ns->ed.mult.sdmult, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpdens, ns->ed.mult.stn);
      // Compute the density at the mid point
      ns->cc.dens   = compute_value_at_mid_point(rhol, rhor);
      // Get the cell fraction in the left cell
      real fracl = compute_center_p_left(ns->ed.mult.sdmult, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
      // Get the cell fraction in the right cell
      real fracr = compute_center_p_right(ns->ed.mult.sdmult, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
      real curvaux = 0.0;
      //if (((0.0 < fracr)&&(fracr < 1.0)) || ((0.0 < fracl)&&(fracl < 1.0))){
      real wwi   = fracr*(1.0 - fracr);
      real wwim1 = fracl*(1.0 - fracl);
      if (wwi + wwim1 > 0.0){
         curvaux = (wwi*curvr + wwim1*curvl)/(wwi + wwim1);
         //printf("wwi + wwim1 = %.12lf  fracl= %.12lf fracr= %.12lf \n",wwi + wwim1,fracl,fracr);
         //getchar();
      }

      // Get the cell viscosity in the left cell
      ns->cc.viscl = compute_center_p_left(ns->ed.mult.sdmult, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpvisc, ns->ed.mult.stn);
      // Get the cell viscosity in the right cell
      ns->cc.viscr = compute_center_p_right(ns->ed.mult.sdmult, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpvisc, ns->ed.mult.stn);
      
      //ns->cc.IF= (fracr - fracl)*ns->cc.curv/(fdelta[dim]);
      ns->cc.IF = (fracr - fracl)*curvaux/(fdelta[dim]);

      //real IFl          = compute_center_p_left(ns->ed.mult.sdmult, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpIF[dim], ns->ed.mult.stn);
      //// Get the cell interfacial force in the right cell
      //real IFr          = compute_center_p_right(ns->ed.mult.sdmult, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpIF[dim], ns->ed.mult.stn);
      //// Compute the interfacial force at the mid point
      //ns->cc.IF   = compute_value_at_mid_point(IFl, IFr);

      if(ns->ed.mult.contr.eoflow_either == true) {
         // Compute the EO source term at the mid point
         ns->cc.Feo   = compute_facet_value_at_point(ns->ed.eo.sfdEOFeo[dim], fcenter, fcenter, 1.0, ns->ed.eo.dpFeo[dim], ns->ed.eo.stnpsi);
      }

      for (int dim2 = 0; dim2 < DIM; dim2++) {
         if(dim2==dim){
            // Get the cell viscosity in the left cell
            ns->cc.viscl = compute_center_p_left(ns->ed.mult.sdmult, fcenter, fdelta, dim2, 0.5, ns->ed.mult.dpvisc, ns->ed.mult.stn);
            // Get the cell viscosity in the right cell
            ns->cc.viscr = compute_center_p_right(ns->ed.mult.sdmult, fcenter, fdelta, dim2, 0.5, ns->ed.mult.dpvisc, ns->ed.mult.stn);
         } else {
            Point p1,p2,p3,p4,p3_,p4_;
            
            POINT_ASSIGN(p1, fcenter);POINT_ASSIGN(p2, fcenter);
            
            p1[dim]=p1[dim]-0.5*fdelta[dim];p2[dim]=p2[dim]+0.5*fdelta[dim];
            
            POINT_ASSIGN(p3, p1);POINT_ASSIGN(p4, p2);
            POINT_ASSIGN(p3_, p1);POINT_ASSIGN(p4_, p2);
            
            p3[dim2]=p3[dim2]+fdelta[dim2];p4[dim2]=p4[dim2]+fdelta[dim2];
            p3_[dim2]=p3_[dim2]-fdelta[dim2];p4_[dim2]=p4_[dim2]-fdelta[dim2];
            
            real v1=compute_value_at_point(ns->ed.mult.sdmult,fcenter,p1,1.0,ns->ed.mult.dpvisc,ns->ed.mult.stn);
            real v2=compute_value_at_point(ns->ed.mult.sdmult,fcenter,p2,1.0,ns->ed.mult.dpvisc,ns->ed.mult.stn);
            
            real v3=compute_value_at_point(ns->ed.mult.sdmult,fcenter,p3,1.0,ns->ed.mult.dpvisc,ns->ed.mult.stn);
            real v4=compute_value_at_point(ns->ed.mult.sdmult,fcenter,p4,1.0,ns->ed.mult.dpvisc,ns->ed.mult.stn);
            
            real v3_=compute_value_at_point(ns->ed.mult.sdmult,fcenter,p3_,1.0,ns->ed.mult.dpvisc,ns->ed.mult.stn);
            real v4_=compute_value_at_point(ns->ed.mult.sdmult,fcenter,p4_,1.0,ns->ed.mult.dpvisc,ns->ed.mult.stn);
            
            ns->cc.viscl=4.0 / (1.0/v1 + 1.0/v2 + 1.0/v3_ + 1.0/v4_);
            ns->cc.viscr=4.0 / (1.0/v1 + 1.0/v2 + 1.0/v3  + 1.0/v4);
         }

         // Compute multiphase viscous term
         ns->cc.du2dx2[dim2] = 0.0;
         Point p;
         POINT_ASSIGN(p, fcenter);

         p[dim2]      = fcenter[dim2] + 0.5*fdelta[dim2];

         real ucl     = compute_facet_u_left(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet_l);
         real ucr     = compute_facet_u_right(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet_r);
         real duidxjr = compute_facet_dudxc(fdelta, dim2, 0.5, ns->cc.ucell, ucl, ucr);

         real vcl     = compute_facet_u_left(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet_l);
         real vcr     = compute_facet_u_right(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet_r);
         real dujdxir = compute_facet_dudxc(fdelta, dim, 0.5, ns->cc.ucell, vcl, vcr);

         p[dim2]      =  fcenter[dim2] - 0.5*fdelta[dim2];

         ucl          = compute_facet_u_left(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet_l);
         ucr          = compute_facet_u_right(sfdu[dim], p, fdelta, dim2, 0.5, dpu[dim], ns->stn, &infacet_r);
         real duidxjl = compute_facet_dudxc(fdelta, dim2, 0.5, ns->cc.ucell, ucl, ucr);

         vcl          = compute_facet_u_left(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet_l);
         vcr          = compute_facet_u_right(sfdu[dim2], p, fdelta, dim, 0.5, dpu[dim2], ns->stn, &infacet_r);
         real dujdxil = compute_facet_dudxc(fdelta, dim, 0.5, ns->cc.ucell, vcl, vcr);

         ns->cc.du2dx2[dim2] += (ns->cc.viscr*duidxjr - ns->cc.viscl*duidxjl)/fdelta[dim2];
         ns->cc.du2dx2[dim2] += (ns->cc.viscr*dujdxir - ns->cc.viscl*dujdxil)/fdelta[dim2];
         
         if(ns->ed.mult.contr.viscoelastic_either == true) {
            // Compute the viscoelastic contribution
            if (dim2 == dim) {
               // Get the tensor in the left cell
               Sl[dim2]          = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.ve.dpS[dim][dim2], ns->ed.stn);
               // Get the tensor in the right cell
               Sr[dim2]          = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.ve.dpS[dim][dim2], ns->ed.stn);
               // Compute the tensor derivative
               ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 0.5, Sl[dim2], Sr[dim2]);
            } else {
               // Get the tensor in the left cell
               Sl[dim2]          = compute_center_p_left_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.ve.dpS[dim][dim2], ns->ed.stn);
               // Get the tensor in the right cell
               Sr[dim2]          = compute_center_p_right_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.ve.dpS[dim][dim2], ns->ed.stn);
               // Compute the tensor derivative
               ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 1.0, Sl[dim2], Sr[dim2]);
            }
         }
         
      }
      break;
   }
 }
}

// Computing the necessary term for the viscoelastic Navier-Stokes equation
void higflow_computational_cell_viscoelastic(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int flid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]) {
   // Set the computational cell
   higflow_computational_cell(ns, sdp, sfdu, flid, fcenter, fdelta, dim, dpu);
   real  Sl[DIM], Sr[DIM];
   // Spatial discretization
   switch (ns->contr.spatialdiscrtype) {
   // Second order
   case ORDER2:
      for (int dim2 = 0; dim2 < DIM; dim2++) {
         // Compute the viscoelastic contribution
         if (dim2 == dim) {
            // Get the tensor in the left cell
            Sl[dim2]          = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.ve.dpS[dim][dim2], ns->ed.stn);
            // Get the tensor in the right cell
            Sr[dim2]          = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.ve.dpS[dim][dim2], ns->ed.stn);
            // Compute the tensor derivative
            ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 0.5, Sl[dim2], Sr[dim2]);
         } else {
            // Get the tensor in the left cell
            Sl[dim2]          = compute_center_p_left_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.ve.dpS[dim][dim2], ns->ed.stn);
            // Get the tensor in the right cell
            Sr[dim2]          = compute_center_p_right_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.ve.dpS[dim][dim2], ns->ed.stn);
            // Compute the tensor derivative
            ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 1.0, Sl[dim2], Sr[dim2]);
         }
      }
      break;
   }
}

//Computing the terms needed for the simulation of viscoelastic flows that exhibit shear-banding behaviour using the VCM model
void higflow_computational_cell_shear_banding_VCM_model(higflow_solver *ns, sim_domain *sdn, int clid, Point ccenter, Point cdelta, int dim, distributed_property *dpn){
	real nc, nr, nl;
	//real  AC, BC;
	int  incell_l, incell_r;
	switch (ns->contr.spatialdiscrtype) {
	// Second order
	case 0:
		// Get the n in the cell center
		nc    = dp_get_value(dpn, clid);
		ns->cc.nABcell = nc;
		//AB  = dp_get_value(ns->ed.vesb.dpA[dim][dim], clid);
		//ns->cc.ABcell = AB;
		// Get the n in the left cell
		nl    = compute_center_p_left_22(sdn, ccenter, cdelta, dim, 1.0, dpn, ns->ed.stn, &incell_l);
		// Get the n in the right cell
		nr    = compute_center_p_right_22(sdn, ccenter, cdelta, dim, 1.0, dpn, ns->ed.stn, &incell_r);
		// Compute the necessary potentials derivative
		ns->cc.dnABdx     = compute_dpdx_at_point(cdelta, dim, 1.0, nl, nr);
		ns->cc.d2nABdx2   = compute_facet_du2dx2(cdelta, dim, 1.0, nc, nl, nr);
		break;
	}
}

// Computing the necessary term for the Navier-Stokes equation
void higflow_computational_cell_electroosmotic(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int flid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]) {
  // Set the computational cell
   higflow_computational_cell(ns, sdp, sfdu, flid, fcenter, fdelta, dim, dpu);
   real  Sl[DIM], Sr[DIM];
   // Spatial discretization
   switch (ns->contr.spatialdiscrtype) {
   // Second order
   case ORDER2:
      // Compute the EO source term at the mid point
      ns->cc.Feo   = compute_facet_value_at_point(ns->ed.eo.sfdEOFeo[dim], fcenter, fcenter, 1.0, ns->ed.eo.dpFeo[dim], ns->ed.eo.stnpsi);

      if (ns->contr.flowtype == VISCOELASTIC){
         // Compute the viscoelastic contribution
         for (int dim2 = 0; dim2 < DIM; dim2++) {
            if (dim2 == dim) {
               // Get the tensor in the left cell
               Sl[dim2]          = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.ve.dpS[dim][dim2], ns->ed.stn);
               // Get the tensor in the right cell
               Sr[dim2]          = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.ve.dpS[dim][dim2], ns->ed.stn);
               // Compute the tensor derivative
               ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 0.5, Sl[dim2], Sr[dim2]);
            } else {
               // Get the tensor in the left cell
               Sl[dim2]          = compute_center_p_left_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.ve.dpS[dim][dim2], ns->ed.stn);
               // Get the tensor in the right cell
               Sr[dim2]          = compute_center_p_right_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.ve.dpS[dim][dim2], ns->ed.stn);
               // Compute the tensor derivative
               ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 1.0, Sl[dim2], Sr[dim2]);
            }
         }
      }
   break;
   }
}

void higflow_computational_cell_electroosmotic_ionic(higflow_solver *ns, sim_domain *sdn, int clid, Point ccenter, Point cdelta, int dim, distributed_property *dpn, sim_stencil *stnn) {
   
   real psic, phic, psir,psil, phir, phil, nc, nr, nl, ul, ur;
   int  incell_l, incell_r;
   int  infacet_l, infacet_r;
   switch (ns->contr.spatialdiscrtype) {
   // Second order
   case ORDER2:
      // Get the n in the cell center
      nc    = dp_get_value(dpn, clid);
      ns->cc.ncell = nc;
      // Get the  psi in the cell center
      psic  = dp_get_value(ns->ed.eo.dppsi, clid);
      // Get the  phi in the cell center
      phic  = dp_get_value(ns->ed.eo.dpphi, clid);

      // Get the n in the left cell
      nl    = compute_center_p_left_22(sdn, ccenter, cdelta, dim, 1.0, dpn, stnn, &incell_l);
      // Get the n in the right cell
      nr    = compute_center_p_right_22(sdn, ccenter, cdelta, dim, 1.0, dpn, stnn, &incell_r);
      // Get the potential in the left cell
      psil  = compute_center_p_left(ns->ed.eo.sdEOpsi, ccenter, cdelta, dim, 1.0, ns->ed.eo.dppsi, ns->ed.eo.stnpsi);
      // Get the potential in the right cell
      psir  = compute_center_p_right(ns->ed.eo.sdEOpsi, ccenter, cdelta, dim, 1.0, ns->ed.eo.dppsi, ns->ed.eo.stnpsi);
      // Get the potential in the left cell
      phil  = compute_center_p_left(ns->ed.eo.sdEOphi, ccenter, cdelta, dim, 1.0, ns->ed.eo.dpphi, ns->ed.eo.stnphi);
      // Get the potential in the right cell
      phir  = compute_center_p_right(ns->ed.eo.sdEOphi, ccenter, cdelta, dim, 1.0, ns->ed.eo.dpphi, ns->ed.eo.stnphi);
      // Compute the necessary potentials derivative
      ns->cc.dpsidx   = compute_dpdx_at_point(cdelta, dim, 1.0, psil, psir);
      ns->cc.dphidx   = compute_dpdx_at_point(cdelta, dim, 1.0, phil, phir);
      ns->cc.dndx     = compute_dpdx_at_point(cdelta, dim, 1.0, nl, nr);
      ns->cc.d2ndx2   = compute_facet_du2dx2(cdelta, dim, 1.0, nc, nl, nr);
      ns->cc.d2psidx2 = compute_facet_du2dx2(cdelta, dim, 1.0, psic, psil, psir);
      //            ns->cc.d2phidx2 = compute_facet_du2dx2(cdelta, dim, 1.0, phic, phil, phir);
      ns->cc.d2phidx2 = 0.0;
      // // Get psi in the cell center
      // ns->cc.psicell = psic;

      switch (ns->ed.eo.contr.convecdiscrtype) {
         case CELL_CENTRAL:
            ul = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet_l);
            ur = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet_r);
            ns->cc.ucell = 0.5*(ul + ur);
         break;
      }

      break;
   }
}

// Computing the necessary term for the viscoelastic Navier-Stokes equation
void higflow_computational_cell_viscoelastic_integral(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int flid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]) {
   // Set the computational cell
   higflow_computational_cell(ns, sdp, sfdu, flid, fcenter, fdelta, dim, dpu);
   real  Sl[DIM], Sr[DIM];
   // Spatial discretization
   switch (ns->contr.spatialdiscrtype) {
   // Second order
   case ORDER2:
      for (int dim2 = 0; dim2 < DIM; dim2++) {
         // Compute the viscoelastic contribution
         if (dim2 == dim) {
            // Get the tensor in the left cell
            Sl[dim2]          = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.im.dpS[dim][dim2], ns->ed.stn);
            // Get the tensor in the right cell
            Sr[dim2]          = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.im.dpS[dim][dim2], ns->ed.stn);
            // Compute the tensor derivative
            ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 0.5, Sl[dim2], Sr[dim2]);
         } else {
            // Get the tensor in the left cell
            Sl[dim2]          = compute_center_p_left_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.im.dpS[dim][dim2], ns->ed.stn);
            // Get the tensor in the right cell
            Sr[dim2]          = compute_center_p_right_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.im.dpS[dim][dim2], ns->ed.stn);
            // Compute the tensor derivative
            ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 1.0, Sl[dim2], Sr[dim2]);
         }
      }
      break;
   }
}

// Computing the necessary term for the viscoelastic Navier-Stokes equation
void higflow_computational_cell_viscoelastic_variable_viscosity(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]) {
	// Set the computational cell
	real  uc[DIM], ul[DIM], ur[DIM], Sl[DIM], Sr[DIM], pl, pr, Fl, Fr;
	real  ull[DIM], urr[DIM];
	real  a, b, c, d, e, tol, fi, conv1,conv2;
	a     = 1.7500;
	b     = 0.3750;
	c     = 0.7500;
	d     = 0.1250;
	e     = 0.2500;
	tol   = 1.0e-14;
	conv1 = 0.0;
	conv2 = 0.0;
	int   infacet, infacet_rr, infacet_ll;
	// Spatial discretization
	switch (ns->contr.spatialdiscrtype) {
	// Second order
	case 0:
		// Set the convective method
		ns->cc.convec_type = ns->contr.convecdiscrtype;
		// Get the velocity in the facet center
		ns->cc.ucell  = dp_get_value(dpu[dim], fgid);
		// Get the facet source tern in the facet center
		ns->cc.F      = compute_facet_value_at_point(ns->sfdF[dim], fcenter, fcenter, 1.0, ns->dpFU[dim], ns->stnF);
		// Get the pressure in the left cell
		pl            = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Get the pressure in the right cell
		pr            = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Compute the pressure derivative
		ns->cc.dpdx   = compute_dpdx_at_point(fdelta, dim, 0.5, pl, pr);
		// Get the cell source term in the left cell
		Fl            = compute_center_p_left(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Get the cell source term in the left cell
		Fr            = compute_center_p_right(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Compute the source term at the mid point
		ns->cc.F     += compute_value_at_mid_point(Fl, Fr);
		for (int dim2 = 0; dim2 < DIM; dim2++) {
			// Get the velocity in the center of facet
			if (dim2 == dim) {
				ns->cc.v[dim2] = ns->cc.ucell;
			} else {
				ns->cc.v[dim2] = compute_facet_u_2(sfdu[dim2], fcenter, fdelta, dim, dim2, 1.0, dpu[dim2], ns->stn);
			}
			// Get the velocity in the left facet center
			ul[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Get the velocity in the right facet center
			ur[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Compute dudxc at facet center
			ns->cc.dudxc[dim2]  = compute_facet_dudxc(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxr at facet center
			ns->cc.dudxr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxl at facet center
			ns->cc.dudxl[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Get the velocity in the right facet center
			urr[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			// Get the velocity in the left facet center
			ull[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			if (ns->cc.convec_type == 2) {
				if (ns->contr.secondconvecdiscrtype == 2) {
					// Compute dudxrr at facet center
					ns->cc.dudxrr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ur[dim2], ns->cc.ucell, urr[dim2]);
					// Compute dudxll at facet center
					ns->cc.dudxll[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ul[dim2], ull[dim2], ns->cc.ucell);
					// Compute terms for second order upwind CUBISTA
				}else if (ns->contr.secondconvecdiscrtype == 1){
					// Get the velocity in the left facet center
					ul[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the left facet center
					ull[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_ll);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity in the right facet center
					ur[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the right facet center
					urr[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_rr);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity v1bar(i+1/2,j+1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_right(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ur[dim]-ul[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv1 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ul[dim]);
								if ((fi >= b) && (fi <= c))
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + b*ur[dim] -d*ul[dim]);
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(e*ns->cc.ucell + c*ur[dim]);
							}
						}
						//v1bar < 0.0
					}else {
						if (fabs(ns->cc.ucell - urr[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ur[dim];
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ur[dim];
							}else {
								if (fi < b)
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(a*ur[dim] - c*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if ((fi >= b) && (fi <= c))
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(c*ur[dim] + b*ns->cc.ucell -d*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ur[dim]);
							}
						}
					}
					// Get the velocity  v2bar(i+1/2,j-1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_left(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ns->cc.ucell-ull[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ul[dim];
						}else {
							fi = (ul[dim] - ull[dim])/(ns->cc.ucell-ull[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ul[dim];
							}else {
								if (fi < b)
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(a*ul[dim] - c*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if ((fi >= b) && (fi <= c))
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(b*ns->cc.ucell + c*ul[dim] - d*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ul[dim]);
							}
						}
					}else {
						//v2bar < 0.0
						if (fabs(ul[dim] - ur[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ur[dim])/(ul[dim] - ur[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv2 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ur[dim]);
								if ((fi >= b) && (fi <= c))
									conv2 = ns->cc.vc[dim2]*(b*ul[dim] + c*ns->cc.ucell - d*ur[dim]);
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ul[dim] + e*ns->cc.ucell);
							}
						}
					}
					ns->cc.vc[dim2] = ((conv1-conv2)/fdelta[dim]);
				}
			}
			// Compute viscoelastic viscous term
			ns->cc.du2dx2[dim2] = compute_facet_du2dx2(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute the viscoelastic contribution
			if (dim2 == dim) {
				// Get the tensor in the left cell
				Sl[dim2]          = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.vevv.dpS[dim][dim2], ns->ed.stn);
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.vevv.dpS[dim][dim2], ns->ed.stn);
				// Compute the tensor derivative
				ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 0.5, Sl[dim2], Sr[dim2]);
			} else {
				// Get the tensor in the left cell
				Sl[dim2]          = compute_center_p_left_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.vevv.dpS[dim][dim2], ns->ed.stn);
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.vevv.dpS[dim][dim2], ns->ed.stn);
				// Compute the tensor derivative
				ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 1.0, Sl[dim2], Sr[dim2]);
			}
		}
		break;
	}
}

// Computing the necessary term for the elastoviscoplastic Navier-Stokes equation
void higflow_computational_cell_elastoviscoplastic(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]) {
	// Set the computational cell
	real  uc[DIM], ul[DIM], ur[DIM], Sl[DIM], Sr[DIM], pl, pr, Fl, Fr;
	real  ull[DIM], urr[DIM];
	real  a, b, c, d, e, tol, fi, conv1,conv2;
	a     = 1.7500;
	b     = 0.3750;
	c     = 0.7500;
	d     = 0.1250;
	e     = 0.2500;
	tol   = 1.0e-14;
	conv1 = 0.0;
	conv2 = 0.0;
	int   infacet, infacet_rr, infacet_ll;
	// Spatial discretization
	switch (ns->contr.spatialdiscrtype) {
	// Second order
	case 0:
		// Set the convective method
		ns->cc.convec_type = ns->contr.convecdiscrtype;
		// Get the velocity in the facet center
		ns->cc.ucell  = dp_get_value(dpu[dim], fgid);
		// Get the facet source tern in the facet center
		ns->cc.F      = compute_facet_value_at_point(ns->sfdF[dim], fcenter, fcenter, 1.0, ns->dpFU[dim], ns->stnF);
		// Get the pressure in the left cell
		pl            = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Get the pressure in the right cell
		pr            = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Compute the pressure derivative
		ns->cc.dpdx   = compute_dpdx_at_point(fdelta, dim, 0.5, pl, pr);
		// Get the cell source term in the left cell
		Fl            = compute_center_p_left(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Get the cell source term in the left cell
		Fr            = compute_center_p_right(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Compute the source term at the mid point
		ns->cc.F     += compute_value_at_mid_point(Fl, Fr);
		for (int dim2 = 0; dim2 < DIM; dim2++) {
			// Get the velocity in the center of facet
			if (dim2 == dim) {
				ns->cc.v[dim2] = ns->cc.ucell;
			} else {
				ns->cc.v[dim2] = compute_facet_u_2(sfdu[dim2], fcenter, fdelta, dim, dim2, 1.0, dpu[dim2], ns->stn);
			}
			// Get the velocity in the left facet center
			ul[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Get the velocity in the right facet center
			ur[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Compute dudxc at facet center
			ns->cc.dudxc[dim2]  = compute_facet_dudxc(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxr at facet center
			ns->cc.dudxr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxl at facet center
			ns->cc.dudxl[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Get the velocity in the right facet center
			urr[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			// Get the velocity in the left facet center
			ull[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			if (ns->cc.convec_type == 2) {
				if (ns->contr.secondconvecdiscrtype == 2) {
					// Compute dudxrr at facet center
					ns->cc.dudxrr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ur[dim2], ns->cc.ucell, urr[dim2]);
					// Compute dudxll at facet center
					ns->cc.dudxll[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ul[dim2], ull[dim2], ns->cc.ucell);
					// Compute terms for second order upwind CUBISTA
				}else if (ns->contr.secondconvecdiscrtype == 1){
					// Get the velocity in the left facet center
					ul[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the left facet center
					ull[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_ll);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity in the right facet center
					ur[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the right facet center
					urr[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_rr);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity v1bar(i+1/2,j+1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_right(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ur[dim]-ul[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv1 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ul[dim]);
								if ((fi >= b) && (fi <= c))
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + b*ur[dim] -d*ul[dim]);
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(e*ns->cc.ucell + c*ur[dim]);
							}
						}
						//v1bar < 0.0
					}else {
						if (fabs(ns->cc.ucell - urr[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ur[dim];
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ur[dim];
							}else {
								if (fi < b)
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(a*ur[dim] - c*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if ((fi >= b) && (fi <= c))
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(c*ur[dim] + b*ns->cc.ucell -d*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ur[dim]);
							}
						}
					}
					// Get the velocity  v2bar(i+1/2,j-1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_left(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ns->cc.ucell-ull[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ul[dim];
						}else {
							fi = (ul[dim] - ull[dim])/(ns->cc.ucell-ull[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ul[dim];
							}else {
								if (fi < b)
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(a*ul[dim] - c*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if ((fi >= b) && (fi <= c))
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(b*ns->cc.ucell + c*ul[dim] - d*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ul[dim]);
							}
						}
					}else {
						//v2bar < 0.0
						if (fabs(ul[dim] - ur[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ur[dim])/(ul[dim] - ur[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv2 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ur[dim]);
								if ((fi >= b) && (fi <= c))
									conv2 = ns->cc.vc[dim2]*(b*ul[dim] + c*ns->cc.ucell - d*ur[dim]);
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ul[dim] + e*ns->cc.ucell);
							}
						}
					}
					ns->cc.vc[dim2] = ((conv1-conv2)/fdelta[dim]);
				}
			}
			// Compute viscoelastic viscous term
			ns->cc.du2dx2[dim2] = compute_facet_du2dx2(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute the viscoelastic contribution
			if (dim2 == dim) {
				// Get the tensor in the left cell
				Sl[dim2]          = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.vepl.dpS[dim][dim2], ns->ed.stn);
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.vepl.dpS[dim][dim2], ns->ed.stn);
				// Compute the tensor derivative
				ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 0.5, Sl[dim2], Sr[dim2]);
			} else {
				// Get the tensor in the left cell
				Sl[dim2]          = compute_center_p_left_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.vepl.dpS[dim][dim2], ns->ed.stn);
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.vepl.dpS[dim][dim2], ns->ed.stn);
				// Compute the tensor derivative
				ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 1.0, Sl[dim2], Sr[dim2]);
			}
		}
		break;
	}
}

// Computing the necessary term for the shear-thickening suspensions Navier-Stokes equation
void higflow_computational_cell_shear_thickening_suspensions(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]) {
	// Set the computational cell
	real  uc[DIM], ul[DIM], ur[DIM], Sl[DIM], Sr[DIM], pl, pr, Fl, Fr;
	real  ull[DIM], urr[DIM];
	real  a, b, c, d, e, tol, fi, conv1,conv2;
	a     = 1.7500;
	b     = 0.3750;
	c     = 0.7500;
	d     = 0.1250;
	e     = 0.2500;
	tol   = 1.0e-14;
	conv1 = 0.0;
	conv2 = 0.0;
	int   infacet, infacet_rr, infacet_ll;
	// Spatial discretization
	switch (ns->contr.spatialdiscrtype) {
	// Second order
	case 0:
		// Set the convective method
		ns->cc.convec_type = ns->contr.convecdiscrtype;
		// Get the velocity in the facet center
		ns->cc.ucell  = dp_get_value(dpu[dim], fgid);
		// Get the facet source tern in the facet center
		ns->cc.F      = compute_facet_value_at_point(ns->sfdF[dim], fcenter, fcenter, 1.0, ns->dpFU[dim], ns->stnF);
		// Get the pressure in the left cell
		pl            = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Get the pressure in the right cell
		pr            = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
		// Compute the pressure derivative
		ns->cc.dpdx   = compute_dpdx_at_point(fdelta, dim, 0.5, pl, pr);
		// Get the cell source term in the left cell
		Fl            = compute_center_p_left(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Get the cell source term in the left cell
		Fr            = compute_center_p_right(ns->sdF, fcenter, fdelta, dim, 0.5, ns->dpF, ns->stnF);
		// Compute the source term at the mid point
		ns->cc.F     += compute_value_at_mid_point(Fl, Fr);
		for (int dim2 = 0; dim2 < DIM; dim2++) {
			// Get the velocity in the center of facet
			if (dim2 == dim) {
				ns->cc.v[dim2] = ns->cc.ucell;
			} else {
				ns->cc.v[dim2] = compute_facet_u_2(sfdu[dim2], fcenter, fdelta, dim, dim2, 1.0, dpu[dim2], ns->stn);
			}
			// Get the velocity in the left facet center
			ul[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}
			// Get the velocity in the right facet center
			ur[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
			//Changes to first ordem upwind if there is a boundary cell
			if ((infacet == 0) && (ns->cc.convec_type == 2)) {
				//DEBUG_INSPECT(fcenter[0],%lf);
				//DEBUG_INSPECT(fcenter[1],%lf);
				//DEBUG_INSPECT(dim,%d);
				//DEBUG_INSPECT(dim2,%d);
				//printf("\n\n");
				ns->cc.convec_type = 1;
			}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			// Compute dudxc at facet center
			ns->cc.dudxc[dim2]  = compute_facet_dudxc(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxr at facet center
			ns->cc.dudxr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute dudxl at facet center
			ns->cc.dudxl[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Get the velocity in the right facet center
			urr[dim2]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			// Get the velocity in the left facet center
			ull[dim2]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet);
			if (ns->cc.convec_type == 2) {
				if (ns->contr.secondconvecdiscrtype == 2) {
					// Compute dudxrr at facet center
					ns->cc.dudxrr[dim2]  = compute_facet_dudxr(fdelta, dim2, 1.0, ur[dim2], ns->cc.ucell, urr[dim2]);
					// Compute dudxll at facet center
					ns->cc.dudxll[dim2]  = compute_facet_dudxl(fdelta, dim2, 1.0, ul[dim2], ull[dim2], ns->cc.ucell);
					// Compute terms for second order upwind CUBISTA
				}else if (ns->contr.secondconvecdiscrtype == 1){
					// Get the velocity in the left facet center
					ul[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the left facet center
					ull[dim]  = compute_facet_u_left(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_ll);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity in the right facet center
					ur[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 1.0, dpu[dim], ns->stn, &infacet);
					// Get the velocity in the right facet center
					urr[dim]  = compute_facet_u_right(sfdu[dim], fcenter, fdelta, dim2, 2.0, dpu[dim], ns->stn, &infacet_rr);
					//Changes to central scheme if there is a outside boundary cell
					if (infacet == 0)  {
						ns->cc.convec_type = 1;
					}
					// Get the velocity v1bar(i+1/2,j+1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_right(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ur[dim]-ul[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv1 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ul[dim]);
								if ((fi >= b) && (fi <= c))
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + b*ur[dim] -d*ul[dim]);
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(e*ns->cc.ucell + c*ur[dim]);
							}
						}
						//v1bar < 0.0
					}else {
						if (fabs(ns->cc.ucell - urr[dim]) <= tol){
							conv1 = ns->cc.vc[dim2]*ur[dim];
						}else {
							fi = (ns->cc.ucell - ul[dim])/(ur[dim] - ul[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv1 = ns->cc.vc[dim2]*ur[dim];
							}else {
								if (fi < b)
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(a*ur[dim] - c*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if ((fi >= b) && (fi <= c))
									if(infacet_rr == 1)  conv1 = ns->cc.vc[dim2]*(c*ur[dim] + b*ns->cc.ucell -d*urr[dim]);
									else                 conv1 = ns->cc.vc[dim2]*ur[dim];
								if (fi > c)
									conv1 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ur[dim]);
							}
						}
					}
					// Get the velocity  v2bar(i+1/2,j-1/2) in the facet center
					ns->cc.vc[dim2]  = compute_facet_u_left(sfdu[dim2], fcenter, fdelta, dim2, 0.5, dpu[dim2], ns->stn, &infacet);
					if (ns->cc.vc[dim2] > 0.0){
						if (fabs(ns->cc.ucell-ull[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ul[dim];
						}else {
							fi = (ul[dim] - ull[dim])/(ns->cc.ucell-ull[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ul[dim];
							}else {
								if (fi < b)
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(a*ul[dim] - c*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if ((fi >= b) && (fi <= c))
									if (infacet_ll == 1)  conv2 = ns->cc.vc[dim2]*(b*ns->cc.ucell + c*ul[dim] - d*ull[dim]);
									else                  conv2 = ns->cc.vc[dim2]*ul[dim];
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ns->cc.ucell + e*ul[dim]);
							}
						}
					}else {
						//v2bar < 0.0
						if (fabs(ul[dim] - ur[dim]) <= tol) {
							conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
						}else {
							fi = (ns->cc.ucell - ur[dim])/(ul[dim] - ur[dim]);
							if ((fi <= 0.0) || (fi >= 1.0)) {
								conv2 = ns->cc.vc[dim2]*ns->cc.ucell;
							}else {
								if (fi < b)
									conv2 = ns->cc.vc[dim2]*(a*ns->cc.ucell - c*ur[dim]);
								if ((fi >= b) && (fi <= c))
									conv2 = ns->cc.vc[dim2]*(b*ul[dim] + c*ns->cc.ucell - d*ur[dim]);
								if (fi > c)
									conv2 = ns->cc.vc[dim2]*(c*ul[dim] + e*ns->cc.ucell);
							}
						}
					}
					ns->cc.vc[dim2] = ((conv1-conv2)/fdelta[dim]);
				}
			}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			// Compute viscous term
			ns->cc.du2dx2[dim2] = compute_facet_du2dx2(fdelta, dim2, 1.0, ns->cc.ucell, ul[dim2], ur[dim2]);
			// Compute the viscoelastic contribution
			if (dim2 == dim) {
				// Get the tensor in the left cell
				Sl[dim2]          = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				// Compute the tensor derivative
				ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 0.5, Sl[dim2], Sr[dim2]);
			} else {
				// Get the tensor in the left cell
				Sl[dim2]          = compute_center_p_left_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				// Get the tensor in the right cell
				Sr[dim2]          = compute_center_p_right_2(ns->ed.sdED, fcenter, fdelta, dim, dim2, 1.0, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				// Compute the tensor derivative
				ns->cc.dSdx[dim2] = compute_dpdx_at_point(fdelta, dim2, 1.0, Sl[dim2], Sr[dim2]);
			}
		}
		break;
	}
}

//Computational cell used to calculate the terms of the particle migration equation (evolution of the fraction volume)
void higflow_computational_cell_volume_fraction_suspensions(higflow_solver *ns, sim_domain *sdphi, int clid, Point ccenter, Point cdelta, int dim,  distributed_property *dpphi){
	real psic, phic, psir,psil, phir, phil, nc, nr, nl;
	real varphic, varphil, varphir;
	real Sl[DIM], Sr[DIM];
	real Sll[DIM], Srr[DIM];
	real Sc[DIM];

	real Sl0[DIM], Sr0[DIM], Sll0[DIM], Srl0[DIM];
	real dTdxr0[DIM], dTdxl0[DIM];

	real Sl1[DIM], Sr1[DIM], Sll1[DIM], Srl1[DIM];
	real dTdxr1[DIM], dTdxl1[DIM];

	int  incell_l, incell_r;
	real frac_tol   = ns->ed.stsp.par.gdrms;
    real eps = 1.0e-4;
	real phij = ns->ed.stsp.par.phi;
	switch (ns->contr.spatialdiscrtype) {
	// Second order
	case 0:
		// Get the volume fraction value in the cell center
		varphic    = dp_get_value(dpphi, clid);
		if(varphic > 1.0){
			varphic = 1.0-eps;
        }
		//if(varphic > phij- frac_tol){
		//	varphic = phij- frac_tol;
        //}
        if(varphic < eps){
			varphic = 0.0;
        }
		ns->cc.phivfcell = varphic;
		// Get the  psi in the cell center
		//psic  = dp_get_value(dppsi, clid);
		// Get the  phi in the cell center
		//phic  = dp_get_value(dpphi, clid);
		// Get the n in the left cell
		//varphil    = compute_center_p_left_22(sdphi, ccenter, cdelta, dim, 1.0, dpphi, ns->ed.stn, &incell_l);
		// Get the n in the right cell
		//varphir    = compute_center_p_right_22(sdphi, ccenter, cdelta, dim, 1.0, dpphi, ns->ed.stn, &incell_r);
		/*
		if (dim == 1) {
			ns->cc.dphivfdx   = compute_dpdx_at_point(cdelta, dim, 0.5, varphil, varphir);
		} else {
			ns->cc.dphivfdx   = 0.0;
		}
		*/
		//ns->cc.dphivfdx   = compute_dpdx_at_point(cdelta, dim, 0.5, varphil, varphir);
		// Get the velocity derivative tensor Du and the microstructure tensor
		/*
		real Du[DIM][DIM];
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				// Get Du
				Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
			}
		}
        // Rate of strain tensor E
		real E[DIM][DIM];
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				E[i][j] = 0.5 * (Du[i][j] + Du[j][i]);
			}
		}
		
		*/
		real varphil0  = compute_center_p_left_22(sdphi, ccenter, cdelta, dim, 1.0, dpphi, ns->ed.stn, &incell_l);
		if(varphil0 > 1.0){
			varphil0 = 1.0-eps;
		}
		if(varphil0 < eps){
			varphil0 = 0.0;
		}
		// Get the n in the right cell
		real varphir0  = compute_center_p_right_22(sdphi, ccenter, cdelta, dim, 1.0, dpphi, ns->ed.stn, &incell_r);
		if(varphir0 > 1.0){
			varphir0 = 1.0-eps;
		}
		if(varphir0 < eps){
			varphir0 = 0.0;
		}
		real KR = pow((1.0-varphir0),4.0);
		real KL = pow((1.0-varphil0),4.0);

		//DEBUG_INSPECT(E[0][1], %lf);
		//real d2Sdx2 = 0.0;
		real d2Tdx2t = 0.0;
		real dTdxpR  = 0.0;
		real dTdxpL = 0.0;
		real d2Tdx2p = 0.0;

		for (int dim2 = 0; dim2 < DIM; dim2++) {
			// Compute the stress contribution
			if (dim2 == dim) {
				Point p;
				POINT_ASSIGN(p, ccenter);
				p[dim]      = ccenter[dim] + cdelta[dim];

				// Get the tensor in the left cell
				Sl0[dim2]          = compute_center_p_left(ns->ed.sdED, p, cdelta, dim2, 0.5, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				//DU1l[dim2] 		  = compute_center_p_left(ns->ed.sdED, ccenter, cdelta, dim2, 0.5, ns->ed.stsp.dpD[dim][dim2], ns->ed.stn);
				//DU12l[dim2]       =  compute_center_p_left(ns->ed.sdED, ccenter, cdelta, dim2, 0.5, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//DU1ll[dim2] 	  = compute_center_p_left(ns->ed.sdED, ccenter, cdelta, dim2, 0.5, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//Sl[dim2]         += 1.0*(ns->ed.stsp.par.eta0)*(DU1l[dim2]);
				//Sl[dim2]         += 2.0*(ns->ed.stsp.par.eta0)*(DU1ll[dim2]+DU12l[dim2]);
				//DEBUG_INSPECT(Sl[0], %lf);
				// Get the tensor in the right cell
				Sr0[dim2]          = compute_center_p_right(ns->ed.sdED, p, cdelta, dim2, 0.5, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				//DU1r[dim2] 		  = compute_center_p_right(ns->ed.sdED, ccenter, cdelta, dim2, 0.5, ns->ed.stsp.dpD[dim][dim2], ns->ed.stn);
				//DU12r[dim2]       = compute_center_p_right(ns->ed.sdED, ccenter, cdelta, dim2, 0.5, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//DU1rr[dim2] 	  = compute_center_p_right(ns->ed.sdED, ccenter, cdelta, dim2, 0.5, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//Sr[dim2]         += 1.0*(ns->ed.stsp.par.eta0)*DU1r[dim2];
				//Sr[dim2]         += 2.0*(ns->ed.stsp.par.eta0)*(DU1rr[dim2]+DU12r[dim2]);
				//DEBUG_INSPECT(Sr[0], %lf);
	
				//Sc[dim2] 		  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				//real DUc          = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpD[dim][dim2], ns->ed.stn);
				//real DUc1         = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//DEBUG_INSPECT(Sc, %lf);
				//Sc 		         += 1.0*(ns->ed.stsp.par.eta0)*DUc;
				//Sc 		         += 2.0*(ns->ed.stsp.par.eta0)*(DUc1+DUc);
				// Compute the tensor derivative
				dTdxr0[dim2] = compute_dpdx_at_point(cdelta, dim2, 0.5, Sl0[dim2], Sr0[dim2]);

				p[dim]      = ccenter[dim] - cdelta[dim];
				//ns->cc.dTdx[dim2] = compute_dpdx_at_point(cdelta, dim2, 0.5, Sl[dim2], Sr[dim2]);
				Sll0[dim2]          = compute_center_p_left(ns->ed.sdED, p, cdelta, dim2, 0.5, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				Srl0[dim2]          = compute_center_p_right(ns->ed.sdED, p, cdelta, dim2, 0.5, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				// Compute the tensor derivative
				dTdxl0[dim2] = compute_dpdx_at_point(cdelta, dim2, 0.5, Sll0[dim2], Srl0[dim2]);
				//d2Tdx2t      = compute_facet_du2dx2(cdelta, dim, 1.0, Sc[dim2], Sl[dim2], Sr[dim2]);
				//d2Tdx2t = compute_facet_du2dx2(cdelta, dim, 1.0, Sc[dim2], Sl[dim2], Sr[dim2]);
				//DEBUG_INSPECT(ns->cc.dTdx[0], %lf);
				d2Tdx2t = (KR*dTdxr0[dim2]-KL*dTdxl0[dim2])/(2.0*cdelta[dim2]);
				//if (dim2 == 1) {
				//	ns->cc.d2Tdx2 += d2Tdx2t;
				//} else {
				//	ns->cc.d2Tdx2 +=  0.0;
				//}
				//ns->cc.d2Tdx2	  += dim2;
				ns->cc.d2Tdx2	  += d2Tdx2t;
			} else {
				Point p1;
				POINT_ASSIGN(p1, ccenter);
				p1[dim]      = ccenter[dim] + cdelta[dim];
				// Get the tensor in the left cell
				//Sl[dim2]          = compute_center_p_left_2(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				//DU1l[dim2]        = compute_center_p_left_2(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpD[dim][dim2], ns->ed.stn);
				Sl1[dim2]          = compute_center_p_left_2(ns->ed.sdED, p1, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				//DU1l[dim2]        = compute_center_p_left_2_r(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpD[dim][dim2], ns->ed.stn);
				//DU12l[dim2]       = compute_center_p_left_2_r(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//DU1ll[dim2]       = compute_center_p_left_2(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//Sl[dim2]         += 2.0*(ns->ed.stsp.par.eta0)*(DU1l[dim2]+DU12l[dim2]);
				//Sl[dim2]         += 2.0*(ns->ed.stsp.par.eta0)*DU1ll[dim2];
				// Get the tensor in the right cell
				Sr1[dim2]          = compute_center_p_right_2(ns->ed.sdED, p1, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				//DU1r[dim2]        = compute_center_p_right_2_r(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpD[dim][dim2], ns->ed.stn);
				//DU12r[dim2]        = compute_center_p_right_2_r(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//DU1rr[dim2]       = compute_center_p_right_2(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//Sr[dim2]         += 1.0*(ns->ed.stsp.par.eta0)*DU1r[dim2];
				//Sr[dim2]         += 2.0*(ns->ed.stsp.par.eta0)*(DU1rr[dim2]+DU12r[dim2]);
				// Compute the tensor derivative
				dTdxr1[dim2] = compute_dpdx_at_point(cdelta, dim2, 1.0, Sl1[dim2], Sr1[dim2]);
				//real Sc 		  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				//real DUc          = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpD[dim][dim2], ns->ed.stn);
				//real DUc1         = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//Sc 		         += 1.0*(ns->ed.stsp.par.eta0)*DUc;
				
				//Sc 		         += 2.0*(ns->ed.stsp.par.eta0)*DUc1;
				p1[dim]      = ccenter[dim] - cdelta[dim];
				Sll1[dim2]          = compute_center_p_left_2(ns->ed.sdED, p1, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				//DU1l[dim2]        = compute_center_p_left_2_r(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpD[dim][dim2], ns->ed.stn);
				//DU12l[dim2]       = compute_center_p_left_2_r(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//DU1ll[dim2]       = compute_center_p_left_2(ns->ed.sdED, ccenter, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpD[dim2][dim], ns->ed.stn);
				//Sl[dim2]         += 2.0*(ns->ed.stsp.par.eta0)*(DU1l[dim2]+DU12l[dim2]);
				//Sl[dim2]         += 2.0*(ns->ed.stsp.par.eta0)*DU1ll[dim2];
				// Get the tensor in the right cell
				Srl1[dim2]          = compute_center_p_right_2(ns->ed.sdED, p1, cdelta, dim, dim2, 1.0, ns->ed.stsp.dpS[dim][dim2], ns->ed.stn);
				// Compute t// Compute the tensor derivative
				dTdxl1[dim2] = compute_dpdx_at_point(cdelta, dim2, 1.0, Sll1[dim2], Srl1[dim2]);
				//real dTdxpL       = compute_dpdx_at_point(cdelta, dim2, 1.0, Sll[dim2], Srr[dim2]);

				d2Tdx2p = (KR*dTdxr1[dim2]-KL*dTdxl1[dim2])/(2.0*cdelta[dim]);
				//ns->cc.d2Tdx2 += 0.0;
				ns->cc.d2Tdx2 	  += d2Tdx2p;
			}
			
		}

		//ns->cc.d2Tdx2 = d2Sdx2;
		//DEBUG_INSPECT(ns->cc.dTdx[0], %lf);
		//DEBUG_INSPECT(ns->cc.dTdx[1], %lf);
		//DEBUG_INSPECT(ns->cc.dTdx[2], %lf);
		//DEBUG_INSPECT(ns->cc.d2Tdx2, %lf);

		//ns->cc.dndx     = compute_dpdx_at_point(cdelta, dim, 1.0, nl, nr);
		//ns->cc.d2ndx2   = compute_facet_du2dx2(cdelta, dim, 1.0, nc, nl, nr);
		//ns->cc.d2psidx2 = compute_facet_du2dx2(cdelta, dim, 1.0, psic, psil, psir);
		//            ns->cc.d2phidx2 = compute_facet_du2dx2(cdelta, dim, 1.0, phic, phil, phir);
		//            ns->cc.d2psidx2 = ns->ed.eo.par.delta*(dp_get_value(ns->ed.eo.dpnminus, clid) - dp_get_value(ns->ed.eo.dpnplus,clid));
		//ns->cc.d2phidx2 = 0.0;
		break;
	}
}
