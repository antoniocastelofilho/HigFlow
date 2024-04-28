// *******************************************************************
//  HiG-Flow Solver Step Multiphase - version 20/01/2022
// *******************************************************************
//#define DEBUG

#include "Debug-c.h"
#include "hig-flow-step-multiphase.h"

char nome_frac[100]; char nome_fracaux[100];
char nome_vx[100]; char nome_vy[100];
char nome_dmult[100]; char nome_Nmult[100];
char nome_Kappa[100]; char nome_press[100];
char nome_IF[100]; char nome_visc[100];
char nome_beta[100]; int count;

void arquivoTempo(int step) {
   
   //sprintf(nome_curv,"DATA/%f_curv.txt", t);
   //sprintf(nome_norm,"DATA/%f_norm.txt", t);
   sprintf(nome_frac,"DATA/%d_frac.txt", step);
   sprintf(nome_fracaux,"DATA/%d_fracaux.txt", step);
   sprintf(nome_vx,"DATA/%d_vx.txt", step);
   sprintf(nome_vy,"DATA/%d_vy.txt", step);
   sprintf(nome_dmult,"DATA/%d_d.txt", step);
   sprintf(nome_Nmult,"DATA/%d_norm.txt", step);
   sprintf(nome_Kappa,"DATA/%d_curv.txt", step);
   sprintf(nome_press,"DATA/%d_press.txt", step);
   sprintf(nome_IF,"DATA/%d_if.txt", step);
   sprintf(nome_visc,"DATA/%d_visc.txt", step);
   sprintf(nome_beta,"DATA/%d_beta.txt", step);
   
   FILE*fp=fopen(nome_frac,"a");
   fclose(fp);
   fp=fopen(nome_fracaux,"a");
   fclose(fp);
   
   fp=fopen(nome_vx,"a");
   fclose(fp);
   
   fp=fopen(nome_vy,"a");
   fclose(fp);
   
   fp=fopen(nome_dmult,"a");
   fclose(fp);
   
   fp=fopen(nome_Nmult,"a");
   fclose(fp);
   
   fp=fopen(nome_Kappa,"a");
   fclose(fp);
   
   fp=fopen(nome_press,"a");
   fclose(fp);
   
   fp=fopen(nome_IF,"a");
   fclose(fp);
   
   fp=fopen(nome_visc,"a");
   fclose(fp);
   
   fp=fopen(nome_beta,"a");
   fclose(fp);
}

void arquivoFrac(char*nome,real x,real y,real frac) {
   FILE *fp=fopen(nome,"a");
   fprintf(fp,"%lf %lf %lf\n",x,y,frac);
   fclose(fp);
}

void arquivoV(char*nome,real vl,real vr) {
   FILE *fp=fopen(nome,"a");
   fprintf(fp,"%lf %lf\n",vl,vr);
   fclose(fp);
}

void arquivod(char*nome,real d) {
   FILE *fp=fopen(nome,"a");
   fprintf(fp,"%lf\n",d);
   fclose(fp);
}

void arquivoN(char*nome,real nx,real ny) {
   FILE *fp=fopen(nome,"a");
   fprintf(fp,"%lf %lf\n",nx,ny);
   fclose(fp);
}

void arquivoIF(char*nome,real ifx,real ify) {
   FILE *fp=fopen(nome,"a");
   fprintf(fp,"%lf %lf\n",ifx,ify);
   fclose(fp);
}

void arquivoK(char*nome,real k) {
   FILE *fp=fopen(nome,"a");
   fprintf(fp,"%lf\n",k);
   fclose(fp);
}

void arquivobeta(char*nome,real k) {
   FILE *fp=fopen(nome,"a");
   fprintf(fp,"%lf\n",k);
   fclose(fp);
}

void arquivovisc(char*nome,real k) {
   FILE *fp=fopen(nome,"a");
   fprintf(fp,"%lf\n",k);
   fclose(fp);
}

void arquivopress(char*nome,real press) {
   FILE *fp=fopen(nome,"a");
   fprintf(fp,"%lf\n",press);
   fclose(fp);
}

void save_cell_values(higflow_solver *ns,int aux) {
    //count=0;
    // Get the local sub-domain for the cells
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    // Get the local sub-domain for the facets
    sim_facet_domain *sfdu[DIM];
    for(int i = 0; i < DIM; i++) {
        sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
    }
    // Get the map for the domain properties
    mp_mapper *mp = sd_get_domain_mapper(sdp);
    // Loop for each cell
    higcit_celliterator *it;
    //printf("*******saving cell properties at %s*************\n",nome_frac);
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Get the cell
        hig_cell *c = higcit_getcell(it);
        // Get the cell identifier
        int clid    = mp_lookup(mp, hig_get_cid(c));
        // Get the center of the cell
        Point ccenter;
        hig_get_center(c, ccenter);
        // Get the delta of the cell
        Point cdelta;
        hig_get_delta(c, cdelta);
        int infacet;
        if(aux==1) {
            real fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
            arquivoFrac(nome_frac,ccenter[0],ccenter[1],fracvol);
            real d = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
            arquivod(nome_dmult,d);
            Point Normal;
            Normal[0] = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
            Normal[1] = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
            
            arquivoN(nome_Nmult,Normal[0],Normal[1]);
            real k = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpcurvature, ns->ed.stn);
            //dpp ou ddeltap?
            real press = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->dpp, ns->ed.stn);
            arquivopress(nome_press,press);
            
            // real beta = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpbeta, ns->ed.stn);
            // arquivobeta(nome_beta,beta);
            
            real visc = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpvisc, ns->ed.stn);
            arquivopress(nome_visc,visc);
            //saving interfacial force
            real dens=compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpdens, ns->ed.stn);
            
            Point Pr;Pr[0]=ccenter[0]+cdelta[0];Pr[1]=ccenter[1];
            Point Pl;Pl[0]=ccenter[0]-cdelta[0];Pl[1]=ccenter[1];
            
            real fracr  = compute_value_at_point(sdp, ccenter, Pr, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
            real fracl  = compute_value_at_point(sdp, ccenter, Pl, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
            
            real fx=(fracr-fracl)/(2*cdelta[0]);
            Point Pt;Pt[1]=ccenter[1]+cdelta[1];Pt[0]=ccenter[0];
            Point Pb;Pb[1]=ccenter[1]-cdelta[1];Pb[0]=ccenter[0];
            
            real fract  = compute_value_at_point(sdp, ccenter, Pt, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
            real fracb  = compute_value_at_point(sdp, ccenter, Pb, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
            real fy=(fract-fracb)/(2*cdelta[1]);
            real ifx=k*fx/(10.0*dens);
            real ify=k*fy/(10.0*dens);
            
            arquivoIF(nome_IF,ifx,ify);
            int dim=0;
            real ul = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            real ur = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            dim=1;
            
            ul = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            ur = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            continue;
        }
        real fracvolaux  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvolaux, ns->ed.stn);
        //arquivoFrac(nome_fracaux,ccenter[0],ccenter[1],fracvolaux);
    }
    // Destroy the iterator
    higcit_destroy(it);//printf("count=%d\n",count);
    // Sync the distributed pressure property
    //dp_sync(ns->ed.mult.dpfracvol);
    //dp_sync(ns->ed.mult.dpfracvolaux);
}

// *******************************************************************
// Navier-Stokes step elements
// *******************************************************************
// Computing the curvature
void higflow_compute_curvature_multiphase(higflow_solver *ns) {
    if (DIM == 2) {
        //higflow_compute_curvature_interfacial_force_normal_multiphase_2D(ns);
        //higflow_compute_curvature_interfacial_force_normal_multiphase_2D_hf_elvira(ns);
        //higflow_compute_curvature_interfacial_force_normal_multiphase_2D_hf_elvira_adap(ns);
        higflow_compute_curvature_interfacial_force_normal_multiphase_2D_hf_shirani(ns);
        higflow_compute_distance_multiphase_2D(ns);
    } else if (DIM == 3) {
        //higflow_compute_curvature_multiphase_3D(ns);
    } else {
        printf("Dimension out of range %d\n",DIM);
        exit(1);
    }
}

// *******************************************************************
// Volume Fraction Transport Step with PLIC fractional step
// *******************************************************************
void higflow_plic_advection_volume_fraction(higflow_solver *ns) {
   if (ns->par.step % 2) {
      higflow_plic_advection_volume_fraction_x_direction_imp(ns, 0);
      higflow_plic_copy_fractionaux_to_fraction(ns);
      //higflow_compute_curvature_interfacial_force_normal_multiphase_2D(ns);
      //higflow_compute_curvature_interfacial_force_normal_multiphase_2D_hf_elvira(ns);
      //higflow_compute_curvature_interfacial_force_normal_multiphase_2D_hf_elvira_adap(ns);
      higflow_compute_curvature_interfacial_force_normal_multiphase_2D_hf_shirani(ns);
      higflow_compute_distance_multiphase_2D(ns);
      higflow_plic_advection_volume_fraction_y_direction(ns, 1);
      higflow_plic_copy_fractionaux_to_fraction(ns);
   } else {
      higflow_plic_advection_volume_fraction_y_direction_imp(ns, 1);
      higflow_plic_copy_fractionaux_to_fraction(ns);
      //higflow_compute_curvature_interfacial_force_normal_multiphase_2D(ns);
      //higflow_compute_curvature_interfacial_force_normal_multiphase_2D_hf_elvira(ns);
      //higflow_compute_curvature_interfacial_force_normal_multiphase_2D_hf_elvira_adap(ns);
      higflow_compute_curvature_interfacial_force_normal_multiphase_2D_hf_shirani(ns);
      higflow_compute_distance_multiphase_2D(ns);
      higflow_plic_advection_volume_fraction_x_direction(ns, 0);
      higflow_plic_copy_fractionaux_to_fraction(ns);
   }
   // DEBUG_INSPECT(ns->par.step,%d);
   // Sync the distributed pressure property
   dp_sync(ns->ed.mult.dpfracvol);
}

// Computing viscosity
void higflow_compute_viscosity_multiphase(higflow_solver *ns) {
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
        int clid    = mp_lookup(mp, hig_get_cid(c));
        // Get the center of the cell
        Point ccenter;
        hig_get_center(c, ccenter);
        // Get the delta of the cell
        Point cdelta;
        hig_get_delta(c, cdelta);
        real fracvol = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
        // Calculate the viscosity
        real visc0 = ns->ed.mult.get_viscosity0(ccenter, ns->par.t);
        real visc1 = ns->ed.mult.get_viscosity1(ccenter, ns->par.t);
        //real visc  = (1.0-fracvol) + fracvol*visc1;
        real visc  = (1.0 - fracvol)*visc0 + fracvol*visc1;
        // Set the viscosity in the distributed viscosity property
        dp_set_value(ns->ed.mult.dpvisc, clid, visc);
    }
    // Destroy the iterator
    higcit_destroy(it);
    // Sync the distributed pressure property
    dp_sync(ns->ed.mult.dpvisc);
}

// Computing density
void higflow_compute_density_multiphase(higflow_solver *ns) {
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int dim = 0; dim < DIM; dim++) {
            sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        }
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdp);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid    = mp_lookup(mp, hig_get_cid(c));
            // Get the center of the cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Calculate the density
            real fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
            // Calculate the density
            real dens0 = ns->ed.mult.get_density0(ccenter, ns->par.t);
            real dens1 = ns->ed.mult.get_density1(ccenter, ns->par.t);
            real dens  = (1.0 - fracvol)*dens0 + fracvol*dens1;
            // Set the viscosity in the distributed viscosity property
            dp_set_value(ns->ed.mult.dpdens, clid, dens);
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the distributed pressure property
        dp_sync(ns->ed.mult.dpdens);
}

// // Computing beta viscoelastic
// void higflow_compute_beta_visc_multiphase(higflow_solver *ns) {
//         // Get the local sub-domain for the cells
//         sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
//         // Get the local sub-domain for the facets
//         sim_facet_domain *sfdu[DIM];
//         for(int dim = 0; dim < DIM; dim++) {
//             sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
//         }
//         // Get the map for the domain properties
//         mp_mapper *mp = sd_get_domain_mapper(sdp);
//         // Loop for each cell
//         higcit_celliterator *it;
//         for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
//             // Get the cell
//             hig_cell *c = higcit_getcell(it);
//             // Get the cell identifier
//             int clid    = mp_lookup(mp, hig_get_cid(c));
//             // Get the center of the cell
//             Point ccenter;
//             hig_get_center(c, ccenter);
//             // Get the delta of the cell
//             Point cdelta;
//             hig_get_delta(c, cdelta);
//             // Calculate the density
//             real fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
//             // Calculate the density
//             real beta0 = 1.0;
//             real beta1 = 0.5;
//             real beta  = (1.0-fracvol)*beta0 + fracvol*beta1;
//             // Set the viscosity in the distributed viscosity property
//             dp_set_value(ns->ed.mult.dpbeta, clid, beta);
//         }
//         // Destroy the iterator
//         higcit_destroy(it);
//         // Sync the distributed beta property
//         dp_sync(ns->ed.mult.dpbeta);
// }

// // Computing S viscoelastic
// void higflow_compute_S_visc_multiphase(higflow_solver *ns) {
//        // Get the local sub-domain for the cells
//         sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
//         // Get the local sub-domain for the facets
//         sim_facet_domain *sfdu[DIM];
//         for(int dim = 0; dim < DIM; dim++) {
//             sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
//         }
//         // Get the map for the domain properties
//         mp_mapper *mp = sd_get_domain_mapper(sdp);
//         // Loop for each cell
//         higcit_celliterator *it;
//         for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
//             // Get the cell
//             hig_cell *c = higcit_getcell(it);
//             // Get the cell identifier
//             int clid    = mp_lookup(mp, hig_get_cid(c));
//             // Get the center of the cell
//             Point ccenter;
//             hig_get_center(c, ccenter);
//             // Get the delta of the cell
//             Point cdelta;
//             hig_get_delta(c, cdelta);
//             // Calculate the density
//             real fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
//             // Calculate the density
            
//             //real dens  = (1.0-fracvol)*dens0 + fracvol*dens1;
//             // Set the viscosity in the distributed viscosity property
            
//             real S[DIM][DIM];
//             for (int i = 0; i < DIM; i++) {
//                 for (int j = 0; j < DIM; j++) {
//                     // Get Kernel
//                     S[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
//                     //S[i][j]  = (1.0 - fracvol)*S[i][j];
//                     S[i][j]  = fracvol*S[i][j];
//                     //printf("S[%d][%d]=%lf   ",i,j,S[i][j]);
//                     dp_set_value(ns->ed.mult.dpS[i][j], clid, S[i][j]);
//                 }
                
//             }
//             //printf("\n");
//         }
//         // Destroy the iterator
//         higcit_destroy(it);
//         // Sync the distributed beta property
//         // Sync the distributed pressure property
//         for (int i = 0; i < DIM; i++) {
//             for (int j = 0; j < DIM; j++) {
//                 dp_sync(ns->ed.mult.dpS[i][j]);
//             }
//         }
// }

// Fraction Correction  
void fraction_correction_at_set(real *fracvol){
   if (fabs(*fracvol - 1.0) < 1.0e-10) {
      *fracvol = 1.0;
   } else if (fabs(*fracvol) < 1.0e-10) {
      *fracvol = 0.0;
   }
   // Para os casos onde as fracoes forem muito grande
   if (*fracvol > 1.0) {
      *fracvol = 1.0;
   } else if (*fracvol < 0.0) {
      *fracvol = 0.0;
   }
}

// ***********************************************************************
// Volume Fraction Transport Step with PLIC fractional step on direction x
// ***********************************************************************
void higflow_plic_advection_volume_fraction_x_direction(higflow_solver *ns, int dim) {
   // Get the local sub-domain for the cells
   sim_domain *sdp = psd_get_local_domain(ns->psdp);
   // Get the local sub-domain for the facets
   sim_facet_domain *sfdu[DIM];
   for(int i = 0; i < DIM; i++) {
      sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
   }
   // Get the map for the domain properties
   mp_mapper *mp = sd_get_domain_mapper(sdp);
   // Loop for each cell
   higcit_celliterator *it;
   real tol_u = 1.0e-14;
   for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
      // Get the cell
      hig_cell *c = higcit_getcell(it);
      // Get the cell identifier
      int clid    = mp_lookup(mp, hig_get_cid(c));
      // Get the center of the cell
      Point ccenter;
      hig_get_center(c, ccenter);
      // Get the delta of the cell
      Point cdelta;
      hig_get_delta(c, cdelta);
      // Get the velocity at facet
      int infacet;
      // Get the velocity in the left facet center
      real ul = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
      // Get the velocity in the right facet center
      real ur = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
      //real ul = 1.0;ul=sin(3.1416*ns->par.t);ul=-4.0*(ccenter[1] - 1.0)*ccenter[1];
      //real ur = 1.0;ur=sin(3.1416*ns->par.t);ul=-4.0*(ccenter[1] - 1.0)*ccenter[1];
      //real ul = -4.0*(ccenter[1] - 1.0)*ccenter[1];
      //real ur = -4.0*(ccenter[1] - 1.0)*ccenter[1];
      /*real piii = 3.1415926535897932384626433832795029;
      real ul, ur;
      if (ns->par.t < 1.0) {
         ul=-2.0*pow(sin(piii*ccenter[0]),2.0)*sin(piii*ccenter[1])*cos(piii*ccenter[1]);
         ur = ul;
         //ur=2.0*pow(sin(piii*ccenter[1]),2.0)*sin(piii*ccenter[0])*cos(piii*ccenter[0]);
         //arquivoV(nome_vx,ccenter[0],ccenter[1],ul,ur);
      } else {
         ul=2.0*pow(sin(piii*ccenter[0]),2.0)*sin(piii*ccenter[1])*cos(piii*ccenter[1]);
         ur = ul;
         //ur=-2.0*pow(sin(piii*ccenter[1]),2.0)*sin(piii*ccenter[0])*cos(piii*ccenter[0]);
      }*/
      // x
      //real ur = 0.5 - ccenter[1];
      //real ul = 0.5 - ccenter[1];
      // y
      //real ur = ccenter[0]-0.5;
      //real ul = ccenter[0]-0.5;
      
      Point Normal, p, Delta_New;
      real d, fracvol,fracr,fracl;
      
      p[0]=ccenter[0];p[1]=ccenter[1];
      p[dim]=p[dim]+0.5*cdelta[dim];
      fracr=compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
      
      p[0]=ccenter[0];p[1]=ccenter[1];
      p[dim]=p[dim]-0.5*cdelta[dim];
      fracl=compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
      
      fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
      
      fraction_correction_at_get(&fracr);
      fraction_correction_at_get(&fracl);
      fraction_correction_at_get(&fracvol);
      
      // Area
      real Ar = 0.0;
      real Al = 0.0;
      fracr = 0.0;
      fracl = 0.0;
      //  Right Facet
      if(fabs(ur)>tol_u) {
         if (fabs(ur) * ns->par.dt > 0.5 * cdelta[0]) {
            printf("Time step is large!!!\n");
            exit(1);
         }
         if(ur>0.0){
            p[0] = ccenter[0];
            p[1] = ccenter[1];
            // Normal
            Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
            Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
            // Correction of Normal
            normal_correction_at_get(Normal);
            // Fraction
            fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
            // Fraction correction
            fraction_correction_at_get(&fracvol);
            if(fracvol==1.0 || fracvol==0.0 || (Normal[0]==0.0 && Normal[1]==0.0)){
               Ar = fracvol*cdelta[1]*ur*ns->par.dt;
            } else {
               d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
               d = d - 0.5*(cdelta[0]-ur*ns->par.dt)*Normal[0];
               Delta_New[0] = ur*ns->par.dt;
               Delta_New[1] = cdelta[1];
               Ar = area_left_line_origin_center(Normal, Delta_New, d);
            }
            fracr = Ar/(fabs(ur)*ns->par.dt*cdelta[1]);
         } else {
            p[0] = ccenter[0]+cdelta[0];
            p[1] = ccenter[1];
            // Normal
            Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
            Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
            // Correction of Normal
            normal_correction_at_get(Normal);
            // Fraction
            fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
            // Fraction correction
            fraction_correction_at_get(&fracvol);
            if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
               Ar = fracvol * cdelta[1] * fabs(ur) * ns->par.dt;
            } else {
               d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
               Delta_New[0] = fabs(ur)*ns->par.dt;
               Delta_New[1] = cdelta[1];
               d = d + 0.5 * (cdelta[0] - fabs(ur) * ns->par.dt)*Normal[0];
               Ar = area_left_line_origin_center(Normal, Delta_New, d);
            }
            fracr = Ar/(fabs(ur)*ns->par.dt*cdelta[1]);
         }
      }
      //  Left Facet
      if (fabs(ul) > tol_u) {
         if (fabs(ul) * ns->par.dt > 0.5 * cdelta[0]) {
            printf("Time step is large!!!\n");
            exit(1);
         }
         if (ul > 0.0) {
            p[0] = ccenter[0] - cdelta[0];
            p[1] = ccenter[1];
            // Normal
            Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
            Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
            // Correction of Normal
            normal_correction_at_get(Normal);
            // Fraction
            fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
            // Fraction correction
            fraction_correction_at_get(&fracvol);
            if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
               Al = fracvol * cdelta[1] * ul * ns->par.dt;
            } else {
               d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
               d = d - 0.5 * (cdelta[0] - ul * ns->par.dt) * Normal[0];
               Delta_New[0] = ul*ns->par.dt;
               Delta_New[1] = cdelta[1];
               Al = area_left_line_origin_center(Normal, Delta_New, d);
            }
            fracl = Al/(fabs(ul)*ns->par.dt*cdelta[1]);
         } else {
            p[0] = ccenter[0];
            p[1] = ccenter[1];
            // Normal
            Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
            Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
            // Correction of Normal
            normal_correction_at_get(Normal);
            // Fraction
            fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
            fraction_correction_at_get(&fracvol);
            if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
               Al = fracvol * cdelta[1] * fabs(ul) * ns->par.dt;
            } else {
               d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
               d = d + 0.5 * (cdelta[0] - fabs(ul) * ns->par.dt) * Normal[0];
               Delta_New[0] = fabs(ul)*ns->par.dt;
               Delta_New[1] = cdelta[1];
               Al = area_left_line_origin_center(Normal, Delta_New, d);
            }
            fracl = Al/(fabs(ul)*ns->par.dt*cdelta[1]);
         }
      }
         
      // Fraction
      fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
      
      real fracvolaux;
      fraction_correction_at_get(&fracr);
      fraction_correction_at_get(&fracl);
      
      /*if(fracr==1.0 && fracl==1.0){
         fracvolaux=fracvol;
      }else{*/
      fracvolaux = fracvol*(1 + ns->par.dt*(ur - ul)/cdelta[0]) - ns->par.dt*(fracr*ur - fracl*ul)/cdelta[0];
      //}
      
      // Auxiliary fraction correction
      fraction_correction_at_set(&fracvolaux);

      UPDATE_RESIDUAL_BUFFER_CELL(ns, dp_get_value(ns->ed.mult.dpfracvolaux, clid), fracvolaux, c, ccenter)
      
      dp_set_value(ns->ed.mult.dpfracvolaux, clid, fracvolaux);
   }
   // Destroy the iterator
   higcit_destroy(it);

   UPDATE_RESIDUALS(ns, ns->residuals->fracvol_adv[0])

   // Sync the distributed vol frac aux property
   dp_sync(ns->ed.mult.dpfracvolaux);
}

// ***********************************************************************
// Volume Fraction Transport Step with PLIC fractional step on direction y
// ***********************************************************************
void higflow_plic_advection_volume_fraction_y_direction(higflow_solver *ns, int dim){
    if (ns->contr.flowtype == MULTIPHASE) {
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->psdp);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int i = 0; i < DIM; i++) {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdp);
        // Loop for each cell
        higcit_celliterator *it;
        real tol_u = 1.0e-14;
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid    = mp_lookup(mp, hig_get_cid(c));
            // Get the center of the cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Get the velocity at facet
            int infacet;
            // Get the velocity in the left facet center
            real ul = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            // Get the velocity in the right facet center
            real ur = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            //real ul = 0.0;//ul=0.4*cos(3.1416*ns->par.t);
            //real ur = 0.0;//ur=0.4*cos(3.1416*ns->par.t);
            /*real piii = 3.1415926535897932384626433832795029;
            real ul, ur;
            if (ns->par.t < 1.0) {
               //ul=-2.0*pow(sin(piii*ccenter[0]),2.0)*sin(piii*ccenter[1])*cos(piii*ccenter[1]);
               ur=2.0*pow(sin(piii*ccenter[1]),2.0)*sin(piii*ccenter[0])*cos(piii*ccenter[0]);
               ul = ur;
               //arquivoV(nome_vx,ccenter[0],ccenter[1],ul,ur);
            } else {
               //ul=2.0*pow(sin(piii*ccenter[0]),2.0)*sin(piii*ccenter[1])*cos(piii*ccenter[1]);
               ur=-2.0*pow(sin(piii*ccenter[1]),2.0)*sin(piii*ccenter[0])*cos(piii*ccenter[0]);
               ul = ur;
            }*/
            // x
            //real ur = 0.5 - ccenter[1];
            //real ul = 0.5 - ccenter[1];
            // y
            //real ur = ccenter[0]-0.5;
            //real ul = ccenter[0]-0.5;
            //arquivoV(nome_vy,ccenter[0],ccenter[1],ul,ur);
            Point Normal, p, Delta_New;
            real d, fracvol,fracr,fracl;
         
            p[0]=ccenter[0];p[1]=ccenter[1];
            p[dim]=p[dim]+0.5*cdelta[dim];
            fracr=compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
            
            p[0]=ccenter[0];p[1]=ccenter[1];
            p[dim]=p[dim]-0.5*cdelta[dim];
            fracl=compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
            
            
            fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
            
            fraction_correction_at_get(&fracr);
            fraction_correction_at_get(&fracl);
            fraction_correction_at_get(&fracvol);
            
            /*if((fracr==1.0&&fracl==1.0&&fracvol==1.0)||(fracr==0.0&&fracl==0.0&&fracvol==0.0)){
               
               dp_set_value(ns->ed.mult.dpfracvolaux, clid, fracvol);
               continue;
               
            }*/

            // Area
            real Ar = 0.0;
            real Al = 0.0;
            fracr = 0.0;
            fracl = 0.0;
            //  Right Facet (UP)
            if(fabs(ur)>tol_u){
               if (fabs(ur)*ns->par.dt > 0.5*cdelta[0]) {
                  printf("Time step is large!!!\n");
                  exit(1);
               }
               if(ur>0.0){
                  p[0] = ccenter[0];
                  p[1] = ccenter[1];
                  // Normal
                  Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
                  Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
                  // Correction of Normal
                  normal_correction_at_get(Normal);
                  // Fraction
                  fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
                  // Fraction correction
                  fraction_correction_at_get(&fracvol);
                  if(fracvol==1.0 || fracvol==0.0 || (Normal[0]==0.0 && Normal[1]==0.0)){
                     Ar = fracvol*cdelta[0]*ur*ns->par.dt;
                  } else {
                     d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
                     d = d - 0.5*(cdelta[1]-ur*ns->par.dt)*Normal[1];
                     Delta_New[0] = cdelta[0];
                     Delta_New[1] = ur*ns->par.dt;
                     Ar = area_left_line_origin_center(Normal, Delta_New, d);
                  }
                  fracr = Ar/(fabs(ur)*ns->par.dt*cdelta[0]);
               } else {
                  p[0] = ccenter[0];
                  p[1] = ccenter[1]+cdelta[1];
                  // Normal
                  Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
                  Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
                  // Correction of Normal
                  normal_correction_at_get(Normal);
                  // Fraction
                  fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
                  // Fraction correction
                  fraction_correction_at_get(&fracvol);
                  if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
                     Ar = fracvol * cdelta[0] * fabs(ur) * ns->par.dt;
                  } else {
                     d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
                     d = d + 0.5 * (cdelta[1] - fabs(ur) * ns->par.dt)*Normal[1];
                     Delta_New[0] = cdelta[0];
                     Delta_New[1] = fabs(ur)*ns->par.dt;
                     Ar = area_left_line_origin_center(Normal, Delta_New, d);
                  }
                  fracr = Ar/(fabs(ur)*ns->par.dt*cdelta[0]);
               }
            }
            //  Left Facet (DWON)
            if (fabs(ul) > tol_u) {
               if (fabs(ur)*ns->par.dt > 0.5*cdelta[0]) {
                  printf("Time step is large!!!\n");
                  exit(1);
               }
               if (ul > 0.0) {
                  p[0] = ccenter[0];
                  p[1] = ccenter[1]-cdelta[1];
                  // Normal
                  Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
                  Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
                  // Correction of Normal
                  normal_correction_at_get(Normal);
                  // Fraction
                  fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
                  // Fraction correction
                  fraction_correction_at_get(&fracvol);
                  if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
                     Al = fracvol * cdelta[0] * ul * ns->par.dt;
                  } else {
                     d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
                     d = d - 0.5 * (cdelta[1] - ul * ns->par.dt) * Normal[1];
                     Delta_New[0] = cdelta[0];
                     Delta_New[1] = ul*ns->par.dt;
                     Al = area_left_line_origin_center(Normal, Delta_New, d);
                  }
                  fracl = Al/(fabs(ul)*ns->par.dt*cdelta[0]);
               } else {
                  p[0] = ccenter[0];
                  p[1] = ccenter[1];
                  // Normal
                  Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
                  Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
                  // Correction of Normal
                  normal_correction_at_get(Normal);
                  // Fraction
                  fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
                  fraction_correction_at_get(&fracvol);
                  if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
                     Al = fracvol * cdelta[0] * fabs(ul) * ns->par.dt;
                  } else {
                     d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
                     d = d + 0.5 * (cdelta[1] - fabs(ul) * ns->par.dt) * Normal[1];
                     Delta_New[0] = cdelta[0];
                     Delta_New[1] = fabs(ul)*ns->par.dt;
                     Al = area_left_line_origin_center(Normal, Delta_New, d);
                  }
                  fracl = Al/(fabs(ul)*ns->par.dt*cdelta[0]);
               }
            }
            
            // Fraction
            fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
            real fracvolaux;
            fraction_correction_at_get(&fracr);
            fraction_correction_at_get(&fracl);
            
            /*if(fracr==1.0 && fracl==1.0){
               fracvolaux=fracvol;
            }else{*/
            fracvolaux = fracvol*(1 + ns->par.dt*(ur - ul)/cdelta[1]) - ns->par.dt*(fracr*ur - fracl*ul)/cdelta[1];
            //}

            // Auxiliary fraction correction
            fraction_correction_at_set(&fracvolaux);

            UPDATE_RESIDUAL_BUFFER_CELL(ns, dp_get_value(ns->ed.mult.dpfracvolaux, clid), fracvolaux, c, ccenter)
            
            dp_set_value(ns->ed.mult.dpfracvolaux, clid, fracvolaux);
        }
        // Destroy the iterator
        higcit_destroy(it);

        UPDATE_RESIDUALS(ns, ns->residuals->fracvol_adv[1])

        // Sync the distributed vol frac aux property
        dp_sync(ns->ed.mult.dpfracvolaux);
    }
}

// ***********************************************************************
// ***********************************************************************
// Volume Fraction Transport Step with PLIC fractional step on direction x
// ***********************************************************************
void higflow_plic_advection_volume_fraction_x_direction_imp(higflow_solver *ns, int dim) {
    if (ns->contr.flowtype == MULTIPHASE) {
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->psdp);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int i = 0; i < DIM; i++) {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdp);
        // Loop for each cell
        higcit_celliterator *it;
        real tol_u = 1.0e-14;
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid    = mp_lookup(mp, hig_get_cid(c));
            // Get the center of the cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Get the velocity at facet
            int infacet;
            // Get the velocity in the left facet center
            real ul = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            //Get the velocity in the right facet center
            real ur = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            //real ul = 1.0;ul=sin(3.1416*ns->par.t);ul=-4.0*(ccenter[1] - 1.0)*ccenter[1];
            //real ur = 1.0;ur=sin(3.1416*ns->par.t);ul=-4.0*(ccenter[1] - 1.0)*ccenter[1];
            //real ul = -4.0*(ccenter[1] - 1.0)*ccenter[1];
            //real ur = -4.0*(ccenter[1] - 1.0)*ccenter[1];
            /*real piii = 3.1415926535897932384626433832795029;
            real ul, ur;
            if (ns->par.t < 1.0) {
               ul=-2.0*pow(sin(piii*ccenter[0]),2.0)*sin(piii*ccenter[1])*cos(piii*ccenter[1]);
               ur = ul;
               //ur=2.0*pow(sin(piii*ccenter[1]),2.0)*sin(piii*ccenter[0])*cos(piii*ccenter[0]);
               //arquivoV(nome_vx,ccenter[0],ccenter[1],ul,ur);
            } else {
               ul=2.0*pow(sin(piii*ccenter[0]),2.0)*sin(piii*ccenter[1])*cos(piii*ccenter[1]);
               ur = ul;
               //ur=-2.0*pow(sin(piii*ccenter[1]),2.0)*sin(piii*ccenter[0])*cos(piii*ccenter[0]);
            }*/
            // x
            //real ur = 0.5 - ccenter[1];
            //real ul = 0.5 - ccenter[1];
            // y
            //real ur = ccenter[0]-0.5;
            //real ul = ccenter[0]-0.5;
            //arquivoV(nome_vx,ccenter[0],ccenter[1],ul,ur);
            Point Normal, p, Delta_New;
            real d, fracvol,fracr,fracl;

            p[0]=ccenter[0];p[1]=ccenter[1];
            p[dim]=p[dim]+0.5*cdelta[dim];
            fracr=compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);

            p[0]=ccenter[0];p[1]=ccenter[1];
            p[dim]=p[dim]-0.5*cdelta[dim];
            fracl=compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);


            fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);

            fraction_correction_at_get(&fracr);
            fraction_correction_at_get(&fracl);
            fraction_correction_at_get(&fracvol);

            /*if((fracr==1.0&&fracl==1.0&&fracvol==1.0)||(fracr==0.0&&fracl==0.0&&fracvol==0.0)){

               dp_set_value(ns->ed.mult.dpfracvolaux, clid, fracvol);
               continue;

            }*/

            // Area
            real Ar = 0.0;
            real Al = 0.0;
            fracr = 0.0;
            fracl = 0.0;
            //  Right Facet
            if(fabs(ur)>tol_u){
               if (fabs(ur) * ns->par.dt > 0.5 * cdelta[0]) {
                  printf("Time step is large!!!\n");
                  exit(1);
               }
               if(ur>0.0){
                  p[0] = ccenter[0];
                  p[1] = ccenter[1];
                  // Normal
                  Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
                  Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
                  // Correction of Normal
                  normal_correction_at_get(Normal);
                  // Fraction
                  fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
                  // Fraction correction
                  fraction_correction_at_get(&fracvol);
                  if(fracvol==1.0 || fracvol==0.0 || (Normal[0]==0.0 && Normal[1]==0.0)){
                     Ar = fracvol*cdelta[1]*ur*ns->par.dt;
                  } else {
//                     printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
                     d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
//                     printf("d = %lf\n",d);
//                     printf("ur*dt = %lf\n",ur*ns->par.dt);
//                     printf("cdelta[0] = %lf, cdelta[1] = %lf\n",cdelta[0],cdelta[1]);
                     d = d - 0.5*(cdelta[0]-ur*ns->par.dt)*Normal[0];
//                     printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
                     Delta_New[0] = ur*ns->par.dt;
                     Delta_New[1] = cdelta[1];
                     Ar = area_left_line_origin_center(Normal, Delta_New, d);
//                     printf("d_c_r = %lf, Ar = %lf, n_x_r = %lf, n_y_r = %lf\n",d,Ar,Normal[0],Normal[1]);
//                     printf("Frac = %lf, x = %lf, y = %lf, ur = %lf\n",fracvol,ccenter[0],ccenter[1],ur);
//                     getchar();
                  }
                  fracr = Ar/(fabs(ur)*ns->par.dt*cdelta[1]);
               } else {
                  p[0] = ccenter[0]+cdelta[0];
                  p[1] = ccenter[1];
                  // Normal
                  Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
                  Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
                  // Correction of Normal
                  normal_correction_at_get(Normal);
                  // Fraction
                  fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
                  // Fraction correction
                  fraction_correction_at_get(&fracvol);
                  if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
                     Ar = fracvol * cdelta[1] * fabs(ur) * ns->par.dt;
                  } else {
                     d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
                     Delta_New[0] = fabs(ur)*ns->par.dt;
                     Delta_New[1] = cdelta[1];
                     d = d + 0.5 * (cdelta[0] - fabs(ur) * ns->par.dt)*Normal[0];
                     Ar = area_left_line_origin_center(Normal, Delta_New, d);
                  }
                  fracr = Ar/(fabs(ur)*ns->par.dt*cdelta[1]);
               }
            }
            //  Left Facet
            if (fabs(ul) > tol_u) {
               if (fabs(ul) * ns->par.dt > 0.5 * cdelta[0]) {
                  printf("Time step is large!!!\n");
                  exit(1);
               }
               if (ul > 0.0) {
                  p[0] = ccenter[0] - cdelta[0];
                  p[1] = ccenter[1];
                  // Normal
                  Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
                  Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
                  // Correction of Normal
                  normal_correction_at_get(Normal);
                  // Fraction
                  fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
                  // Fraction correction
                  fraction_correction_at_get(&fracvol);
                  if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
                     Al = fracvol * cdelta[1] * ul * ns->par.dt;
                  } else {
                     d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
                     d = d - 0.5 * (cdelta[0] - ul * ns->par.dt) * Normal[0];
                     Delta_New[0] = ul*ns->par.dt;
                     Delta_New[1] = cdelta[1];
                     Al = area_left_line_origin_center(Normal, Delta_New, d);
                  }
                  fracl = Al/(fabs(ul)*ns->par.dt*cdelta[1]);
               } else {
                  p[0] = ccenter[0];
                  p[1] = ccenter[1];
                  // Normal
                  Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
                  Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
                  // Correction of Normal
                  normal_correction_at_get(Normal);
                  // Fraction
                  fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
                  fraction_correction_at_get(&fracvol);
                  if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
                     Al = fracvol * cdelta[1] * fabs(ul) * ns->par.dt;
                  } else {
                     d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
                     d = d + 0.5 * (cdelta[0] - fabs(ul) * ns->par.dt) * Normal[0];
                     Delta_New[0] = fabs(ul)*ns->par.dt;
                     Delta_New[1] = cdelta[1];
                     Al = area_left_line_origin_center(Normal, Delta_New, d);
                  }
                  fracl = Al/(fabs(ul)*ns->par.dt*cdelta[1]);
               }
            }

            // Fraction
            fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);

            real fracvolaux;
            fraction_correction_at_get(&fracr);
            fraction_correction_at_get(&fracl);

            /*if(fracr==1.0 && fracl==1.0){
               fracvolaux=fracvol;
            }else{*/
            fracvolaux = (fracvol - ns->par.dt*(fracr*ur - fracl*ul)/cdelta[0])/(1.0 - ns->par.dt*(ur - ul)/cdelta[0]);
            //}

            // Auxiliary fraction correction
            fraction_correction_at_set(&fracvolaux);

            UPDATE_RESIDUAL_BUFFER_CELL(ns, dp_get_value(ns->ed.mult.dpfracvolaux, clid), fracvolaux, c, ccenter)

            dp_set_value(ns->ed.mult.dpfracvolaux, clid, fracvolaux);
        }
        // Destroy the iterator
        higcit_destroy(it);

        UPDATE_RESIDUALS(ns, ns->residuals->fracvol_adv[0])

        // Sync the distributed vol frac aux property
        dp_sync(ns->ed.mult.dpfracvolaux);
    }
}

// ***********************************************************************
// Volume Fraction Transport Step with PLIC fractional step on direction y
// ***********************************************************************
void higflow_plic_advection_volume_fraction_y_direction_imp(higflow_solver *ns, int dim){
    if (ns->contr.flowtype == MULTIPHASE) {
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->psdp);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int i = 0; i < DIM; i++) {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdp);
        // Loop for each cell
        higcit_celliterator *it;
        real tol_u = 1.0e-14;
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid    = mp_lookup(mp, hig_get_cid(c));
            // Get the center of the cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Get the velocity at facet
            int infacet;
            // Get the velocity in the left facet center
            real ul = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            // Get the velocity in the right facet center
            real ur = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            //real ul = 0.0;//ul=0.4*cos(3.1416*ns->par.t);
            //real ur = 0.0;//ur=0.4*cos(3.1416*ns->par.t);
            /*real piii = 3.1415926535897932384626433832795029;
            real ul, ur;
            if (ns->par.t < 1.0) {
               //ul=-2.0*pow(sin(piii*ccenter[0]),2.0)*sin(piii*ccenter[1])*cos(piii*ccenter[1]);
               ur=2.0*pow(sin(piii*ccenter[1]),2.0)*sin(piii*ccenter[0])*cos(piii*ccenter[0]);
               ul = ur;
               //arquivoV(nome_vx,ccenter[0],ccenter[1],ul,ur);
            } else {
               //ul=2.0*pow(sin(piii*ccenter[0]),2.0)*sin(piii*ccenter[1])*cos(piii*ccenter[1]);
               ur=-2.0*pow(sin(piii*ccenter[1]),2.0)*sin(piii*ccenter[0])*cos(piii*ccenter[0]);
               ul = ur;
            }*/
            // x
            //real ur = 0.5 - ccenter[1];
            //real ul = 0.5 - ccenter[1];
            // y
            //real ur = ccenter[0]-0.5;
            //real ul = ccenter[0]-0.5;
            //arquivoV(nome_vy,ccenter[0],ccenter[1],ul,ur);
            Point Normal, p, Delta_New;
            real d, fracvol,fracr,fracl;

            p[0]=ccenter[0];p[1]=ccenter[1];
            p[dim]=p[dim]+0.5*cdelta[dim];
            fracr=compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);

            p[0]=ccenter[0];p[1]=ccenter[1];
            p[dim]=p[dim]-0.5*cdelta[dim];
            fracl=compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);


            fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);

            fraction_correction_at_get(&fracr);
            fraction_correction_at_get(&fracl);
            fraction_correction_at_get(&fracvol);

            /*if((fracr==1.0&&fracl==1.0&&fracvol==1.0)||(fracr==0.0&&fracl==0.0&&fracvol==0.0)){

               dp_set_value(ns->ed.mult.dpfracvolaux, clid, fracvol);
               continue;

            }*/

            // Area
            real Ar = 0.0;
            real Al = 0.0;
            fracr = 0.0;
            fracl = 0.0;
            //  Right Facet (UP)
            if(fabs(ur)>tol_u){
               if (fabs(ur)*ns->par.dt > 0.5*cdelta[0]) {
                  printf("Time step is large!!!\n");
                  exit(1);
               }
               if(ur>0.0){
                  p[0] = ccenter[0];
                  p[1] = ccenter[1];
                  // Normal
                  Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
                  Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
                  // Correction of Normal
                  normal_correction_at_get(Normal);
                  // Fraction
                  fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
                  // Fraction correction
                  fraction_correction_at_get(&fracvol);
                  if(fracvol==1.0 || fracvol==0.0 || (Normal[0]==0.0 && Normal[1]==0.0)){
                     Ar = fracvol*cdelta[0]*ur*ns->par.dt;
                  } else {
                     d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
                     d = d - 0.5*(cdelta[1]-ur*ns->par.dt)*Normal[1];
                     Delta_New[0] = cdelta[0];
                     Delta_New[1] = ur*ns->par.dt;
                     Ar = area_left_line_origin_center(Normal, Delta_New, d);
                  }
                  fracr = Ar/(fabs(ur)*ns->par.dt*cdelta[0]);
               } else {
                  p[0] = ccenter[0];
                  p[1] = ccenter[1]+cdelta[1];
                  // Normal
                  Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
                  Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
                  // Correction of Normal
                  normal_correction_at_get(Normal);
                  // Fraction
                  fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
                  // Fraction correction
                  fraction_correction_at_get(&fracvol);
                  if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
                     Ar = fracvol * cdelta[0] * fabs(ur) * ns->par.dt;
                  } else {
                     d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
                     d = d + 0.5 * (cdelta[1] - fabs(ur) * ns->par.dt)*Normal[1];
                     Delta_New[0] = cdelta[0];
                     Delta_New[1] = fabs(ur)*ns->par.dt;
                     Ar = area_left_line_origin_center(Normal, Delta_New, d);
                  }
                  fracr = Ar/(fabs(ur)*ns->par.dt*cdelta[0]);
               }
            }
            //  Left Facet (DWON)
            if (fabs(ul) > tol_u) {
               if (fabs(ur)*ns->par.dt > 0.5*cdelta[0]) {
                  printf("Time step is large!!!\n");
                  exit(1);
               }
               if (ul > 0.0) {
                  p[0] = ccenter[0];
                  p[1] = ccenter[1]-cdelta[1];
                  // Normal
                  Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
                  Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
                  // Correction of Normal
                  normal_correction_at_get(Normal);
                  // Fraction
                  fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
                  // Fraction correction
                  fraction_correction_at_get(&fracvol);
                  if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
                     Al = fracvol * cdelta[0] * ul * ns->par.dt;
                  } else {
                     d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
                     d = d - 0.5 * (cdelta[1] - ul * ns->par.dt) * Normal[1];
                     Delta_New[0] = cdelta[0];
                     Delta_New[1] = ul*ns->par.dt;
                     Al = area_left_line_origin_center(Normal, Delta_New, d);
                  }
                  fracl = Al/(fabs(ul)*ns->par.dt*cdelta[0]);
               } else {
                  p[0] = ccenter[0];
                  p[1] = ccenter[1];
                  // Normal
                  Normal[0] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
                  Normal[1] = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
                  // Correction of Normal
                  normal_correction_at_get(Normal);
                  // Fraction
                  fracvol  = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
                  fraction_correction_at_get(&fracvol);
                  if (fracvol == 1.0 || fracvol == 0.0 || (Normal[0] == 0.0 && Normal[1] == 0.0)) {
                     Al = fracvol * cdelta[0] * fabs(ul) * ns->par.dt;
                  } else {
                     d = compute_value_at_point(sdp, ccenter, p, 1.0, ns->ed.mult.dpdistance, ns->ed.stn);
                     d = d + 0.5 * (cdelta[1] - fabs(ul) * ns->par.dt) * Normal[1];
                     Delta_New[0] = cdelta[0];
                     Delta_New[1] = fabs(ul)*ns->par.dt;
                     Al = area_left_line_origin_center(Normal, Delta_New, d);
                  }
                  fracl = Al/(fabs(ul)*ns->par.dt*cdelta[0]);
               }
            }

            // Fraction
            fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
            real fracvolaux;
            fraction_correction_at_get(&fracr);
            fraction_correction_at_get(&fracl);

            /*if(fracr==1.0 && fracl==1.0){
               fracvolaux=fracvol;
            }else{*/
            fracvolaux = (fracvol - ns->par.dt*(fracr*ur - fracl*ul)/cdelta[1])/(1.0 - ns->par.dt*(ur - ul)/cdelta[1]);
            //}

            // Auxiliary fraction correction
            fraction_correction_at_set(&fracvolaux);

            UPDATE_RESIDUAL_BUFFER_CELL(ns, dp_get_value(ns->ed.mult.dpfracvolaux, clid), fracvolaux, c, ccenter)

            dp_set_value(ns->ed.mult.dpfracvolaux, clid, fracvolaux);
        }
        // Destroy the iterator
        higcit_destroy(it);

        UPDATE_RESIDUALS(ns, ns->residuals->fracvol_adv[1])

        // Sync the distributed vol frac aux property
        dp_sync(ns->ed.mult.dpfracvolaux);
    }
}
//************************************************************************

// ***********************************************************************
//Copy Auxiliary Volume Fraction to Volume Fraction
// ***********************************************************************
void higflow_plic_copy_fractionaux_to_fraction(higflow_solver *ns) {
    if (ns->contr.flowtype == MULTIPHASE) {
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->psdp);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int i = 0; i < DIM; i++) {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdp);
        // Loop for each cell
        higcit_celliterator *it;

        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid    = mp_lookup(mp, hig_get_cid(c));
            // Get the center of the cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
//            real fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
//            arquivoFrac(nome_frac,ccenter[0],ccenter[1],fracvol);
            real fracvolaux  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvolaux, ns->ed.stn);
//            arquivoFrac(nome_fracaux,ccenter[0],ccenter[1],fracvolaux);
            dp_set_value(ns->ed.mult.dpfracvol, clid, fracvolaux);
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the distributed pressure property
        dp_sync(ns->ed.mult.dpfracvol);
        dp_sync(ns->ed.mult.dpfracvolaux);
  }
}

// Correction Normal
void normal_correction_at_get(Point Normal) {
   real tol_n = 1.0e-14;
   for(int i=0;i<DIM;i++){
      if (fabs(Normal[i]) < tol_n) {
         Normal[i] = 0.0;
      }
   }
}

// *******************************************************************
// Volume Fraction Transport Step for the Explicit Euler Method
// *******************************************************************
void higflow_explicit_euler_volume_fraction(higflow_solver *ns) {
    if (ns->contr.flowtype == MULTIPHASE) {
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->psdp);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int i = 0; i < DIM; i++) {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdp);
        // Loop for each cell
        higcit_celliterator *it;
       real volume_old = 0.0;
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid    = mp_lookup(mp, hig_get_cid(c));
            // Get the center of the cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
           real  volcell = 1.0;
       for (int i = 0; i < DIM; i++) volcell *= cdelta[i];
       // Calculate the volume fraction
       real fracvol  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
       volume_old += fracvol*volcell;
       // Get the velocity at cell center 
       real u[DIM], dfracvoldx[DIM];
       hig_flow_velocity_at_center_cell_multiphase(ns, ccenter, cdelta, u);
       // Solving the Transport Equation using the Euler Method
       // Right hand side equation
       real rhs = 0.0;
       int  convecdiscrtype = 1;
           switch (convecdiscrtype) {
              case CELL_CENTRAL: 
                 // Kernel derivative at cell center
                 hig_flow_derivative_fracvol_at_center_cell(ns, ccenter, cdelta, fracvol, dfracvoldx);
                 for (int dim = 0; dim < DIM; dim++) {
                    //Compute convective tensor term in rhs
                    rhs -= u[dim]*dfracvoldx[dim];
                 }
                 break;
              case CELL_CUBISTA: 
                 //Compute convective fracvol term CUBISTA in rhs
                 for (int dim = 0; dim < DIM; dim++) {
                    rhs -= hig_flow_fracvol_term_cubista(ns, ns->dpu[dim], ns->ed.sdED, ns->ed.stn, fracvol, ccenter, cdelta, dim);
                 }
                 break;
           }
           // Compute the Kernel at next time
           fracvol  = fracvol + ns->par.dt * rhs;
           real frac_tol   = 1.0e-14;
           if(fracvol < frac_tol){
              fracvol = 0.0;
           }
           if(fracvol > 1.0 - frac_tol){
              fracvol = 1.0;
           }
           // Store Kernel in S
           dp_set_value(ns->ed.mult.dpcurvature, clid, fracvol);
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the distributed pressure property
        dp_sync(ns->ed.mult.dpcurvature);

        real volume_new = 0.0;
        // Store the Volume Fraction 
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid    = mp_lookup(mp, hig_get_cid(c));
            // Get the center of the cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            real  volcell = 1.0;
            for (int i = 0; i < DIM; i++) volcell *= cdelta[i];
               // Get the volume fraction stored in dpcurvature
               real fracvol  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.dpcurvature, ns->ed.stn);
               // Store volume fraction in dpfracvol
               dp_set_value(ns->ed.mult.dpfracvol, clid, fracvol);
               volume_new += fracvol*volcell;
            }
            // Destroy the iterator
            higcit_destroy(it);
            // Sync the distributed pressure property
            dp_sync(ns->ed.mult.dpfracvol);
            // Print the volume
            printf("===> Volume = %16.10lf <====> Volume Error = %16.13lf <===\n",volume_new,volume_new-volume_old);
        }
}

// Get the velocity at cell center 
void hig_flow_velocity_at_center_cell_multiphase (higflow_solver *ns, Point ccenter, Point cdelta, real u[DIM]) {
    for (int dim = 0; dim < DIM; dim++) {
        // Verity if is in facet
        int infacet;
        // Get the velocity in the left facet center
        real ul = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
        // Get the velocity in the right facet center
        real ur = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
        // Setting the velocity at cell center
        u[dim]  = 0.5*(ul + ur);
    }
}

// Get the derivative of Kernel 
void hig_flow_derivative_fracvol_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, real FVcenter, real dfracvoldx[DIM]) {
    for (int dim = 0; dim < DIM; dim++) {
        int incell_left, incell_right;
        // Get the Kernel in the left cell
        real FVleft = compute_center_p_left_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn, &incell_left);
        // Get the Kernel in the right cell
        real FVright = compute_center_p_right_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn, &incell_right);
        // Compute the Kernel derivative
        if ((incell_left == 1) && (incell_right == 1)) { 
           dfracvoldx[dim] = compute_dpdx_at_point(cdelta, dim, 1.0, FVleft, FVright);
        } else if (incell_right == 1) { 
           dfracvoldx[dim] = compute_dpdxr_at_point(cdelta, dim, 1.0, FVcenter, FVright);
        } else {
           dfracvoldx[dim] = compute_dpdxl_at_point(cdelta, dim, 1.0, FVleft, FVcenter);
        }
    }
}

// *******************************************************************
// Calculate convective term CUBISTA for volume fraction
// *******************************************************************
real hig_flow_fracvol_term_cubista(higflow_solver *ns, distributed_property *dpu, sim_domain *sdp, sim_stencil *stn, real fracvol, Point ccenter, Point cdelta, int dim) {
    real  vbar[DIM], dKdx[dim], kr, krr, kl, kll, kc, a, b, c, d, e, tol, frac_tol, fi, conv1,conv2;
    a     = 1.7500;
    b     = 0.3750;
    c     = 0.7500;
    d     = 0.1250;
    e     = 0.2500;
    tol   = 1.0e-14;
    conv1 = 0.0;
    conv2 = 0.0;
    int   incell_r, incell_l, incell_ll, incell_rr, infacet;
    // Get the kernel at center cell
    kc  = fracvol;
    // Get the low, high, lowlow, highhigh component kernel at center cell
    kl  = compute_center_p_left_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn, &incell_l); 
    kr  = compute_center_p_right_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn, &incell_r); 
    kll = compute_center_p_left_22(ns->ed.sdED, ccenter, cdelta, dim, 2.0, ns->ed.mult.dpfracvol, ns->ed.stn, &incell_ll);
    krr = compute_center_p_right_22(ns->ed.sdED, ccenter, cdelta, dim, 2.0, ns->ed.mult.dpfracvol, ns->ed.stn, &incell_rr);
    // Get the velocity  v1bar(i+1/2,j) in the facet center
    vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
    if (vbar[dim] > 0.0){
        if (fabs(kr - kl) <= tol){
            conv1 = vbar[dim]*kc;
        }else {
            fi = (kc - kl)/(kr - kl);
            if ((fi <= 0.0) || (fi >= 1.0)) {
                conv1 = vbar[dim]*kc;
            }else {
                if (fi < b){ 
                    if (incell_l == 1)                    conv1 = vbar[dim]*(a*kc - c*kl);
                    else                                  conv1 = vbar[dim]*kc;
                }
           if ((fi >= b) && (fi <= c)){
                    if ((incell_l == 1)&&(incell_r == 1)) conv1 = vbar[dim]*(c*kc + b*kr -d*kl);
                    else                                  conv1 = vbar[dim]*kc;
                }
           if (fi > c){ 
                    if (incell_r == 1)                    conv1 = vbar[dim]*(e*kc + c*kr);
                    else                                  conv1 = vbar[dim]*kc;
                }
            }
        }
    //v1bar < 0.0
    }else {
        if ((incell_r == 1) && (incell_rr == 1)){
            if (fabs(kc - krr) <= tol){
                conv1 = vbar[dim]*kr;
            }else {
                fi = (kr - krr)/(kc - krr);
                if ((fi <= 0.0) || (fi >= 1.0)) {
                    conv1 = vbar[dim]*kr;
                }else {
          if (fi < b) 
                        conv1 = vbar[dim]*(a*kr - c*krr);
                    if ((fi >= b) && (fi <= c))
                        conv1 = vbar[dim]*(c*kr + b*kc -d*krr);
               if (fi > c) 
                        conv1 = vbar[dim]*(c*kc + e*kr);
                }
            }
        //Return upwind value at boundary
        }else if ((incell_r == 1) && (incell_rr == 0)){
            if (fabs(kc - krr) <= tol){
                conv1 = vbar[dim]*kr;
            }else {
                fi = (kr- krr)/(kc - krr);
                if ((fi <= 0.0) || (fi >= 1.0)) {
                    conv1 = vbar[dim]*kr;
                }else {
          if (fi <= c) 
                        conv1 = vbar[dim]*kr;
               if (fi > c) 
                        conv1 = vbar[dim]*(c*kc + e*kr);
                }
            }/*
            vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            if (vbar[dim] > 0.0) conv1 = vbar[dim]*kc;
            else                 conv1 = vbar[dim]*kr;
            vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            if (vbar[dim] > 0.0) conv2 = vbar[dim]*kl;
            else                 conv2 = vbar[dim]*kc;
            return ((conv1 - conv2)/cdelta[dim]); */
        }else {
           vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
           if (vbar[dim] > 0.0) conv1 = vbar[dim]*kc;
           else                 conv1 = vbar[dim]*kc;
           vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
           if (vbar[dim] > 0.0) conv2 = vbar[dim]*kl;
           else                 conv2 = vbar[dim]*kc;

           real value_tol = ((conv1 - conv2)/cdelta[dim]);
           return value_tol; 
        }
        
    }
    // Get the velocity  v2bar(i-1/2,j) in the facet center
    vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
    if (vbar[dim] > 0.0){
        if ((incell_l == 1) && (incell_ll == 1)){
            if (fabs(kc-kll) <= tol) {
           conv2 = vbar[dim]*kl;
            }else {
           fi = (kl - kll)/(kc - kll);
           if ((fi <= 0.0) || (fi >= 1.0)) {
               conv2 = vbar[dim]*kl;
           }else {
               if (fi < b)
                   conv2 = vbar[dim]*(a*kl - c*kll);
               if ((fi >= b) && (fi <= c))
                   conv2 = vbar[dim]*(b*kc + c*kl - d*kll);
               if (fi > c)  
                   conv2 = vbar[dim]*(c*kc + e*kl);
           }
       }
        }else if ((incell_l == 1) && (incell_ll == 0)){
            if (fabs(kc-kll) <= tol) {
           conv2 = vbar[dim]*kl;
            }else {
           fi = (kl - kll)/(kc - kll);
           if ((fi <= 0.0) || (fi >= 1.0)) {
               conv2 = vbar[dim]*kl;
           }else {
               if (fi <= c)
                   conv2 = vbar[dim]*kl;
               if (fi > c)  
                   conv2 = vbar[dim]*(c*kc + e*kl);
           }
       }/*
            vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            if (vbar[dim] > 0.0) conv1 = vbar[dim]*kc;
            else                 conv1 = vbar[dim]*kr;
            vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            if (vbar[dim] > 0.0) conv2 = vbar[dim]*kl;
            else                 conv2 = vbar[dim]*kc;
            return ((conv1 - conv2)/cdelta[dim]); */
       }else {
                vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
                if (vbar[dim] > 0.0) conv1 = vbar[dim]*kc;
                else                 conv1 = vbar[dim]*kr;
                vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
                if (vbar[dim] > 0.0) conv2 = vbar[dim]*kc;
                else                 conv2 = vbar[dim]*kc;

                real value_tol = ((conv1 - conv2)/cdelta[dim]);
                return value_tol; 
        } 
    }else {
    //v2bar < 0.0 
        if (fabs(kl - kr) <= tol) {
            conv2 = vbar[dim]*kc;
        }else {
            fi = (kc - kr)/(kl - kr);
            if ((fi <= 0.0) || (fi >= 1.0)) {
                conv2 = vbar[dim]*kc;
            }else {
           if (fi < b){
                    if (incell_r == 1)                    conv2 = vbar[dim]*(a*kc - c*kr);
                    else                                  conv2 = vbar[dim]*kc;
                }
           if ((fi >= b) && (fi <= c)){
                    if ((incell_l == 1)&&(incell_r == 1)) conv2 = vbar[dim]*(c*kc + b*kl -d*kr);
                    else                                  conv2 = vbar[dim]*kc;
                }
           if (fi > c){ 
                    if (incell_l == 1)                    conv2 = vbar[dim]*(e*kc + c*kl);
                    else                                  conv2 = vbar[dim]*kc;
                }
       }
        }
    }

    real value_tol = ((conv1 - conv2)/cdelta[dim]);
    return value_tol;
}

// Navier-Stokes final velocity using the projection method for multiphase flow
void higflow_final_velocity_multiphase(higflow_solver *ns) {
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Loop for each dimension
    higfit_facetiterator *fit;
    for (int dim = 0; dim < DIM; dim++) {
        // Initialize the min and max velocity
        real velmax    = -1.0e16;
        real velmin    =  1.0e16;
        // Get the local partitioned domain for facets
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the map of the distributed properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Get the density in the left cell
            real rhol  = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpdens, ns->ed.stn);
            //real rhol  = 1.0;
            // Get the density in the right cell
            real rhor  = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim, 0.5, ns->ed.mult.dpdens, ns->ed.stn);
            //real rhor  = 1.0;
            // Compute the density in the facet
            real rho = 0.5*(rhol + rhor);
            real pl, pr;
            if (ns->contr.projtype == INCREMENTAL) {
                // Get the pressure in the left cell
                pl    = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->ddeltap, ns->stn);
                // Get the pressure in the right cell
                pr    = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->ddeltap, ns->stn);
            } else {
                // Get the pressure in the left cell
                pl    = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
                // Get the pressure in the right cell
                pr    = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
            }
            // Compute the pressure derivative
            real dpdx  = compute_dpdx_at_point(fdelta, dim, 0.5, pl, pr);
            // Get the intermediate velocity
            real ustar = dp_get_value(ns->dpustar[dim], flid);
            // Compute the final velocity
            real u  = ustar - ns->par.dt*dpdx/rho;

            UPDATE_RESIDUAL_BUFFER_FACET(ns, dp_get_value(ns->dpu[dim], flid), u, f, fcenter)

            // Set the final velocity in the distributed velocity property
            dp_set_value(ns->dpu[dim], flid, u);
            // Update the min and max velocity
            if (u > velmax) velmax = u;
            if (u < velmin) velmin = u;
        }
        // Destroy the iterator
        higfit_destroy(fit);

        UPDATE_RESIDUALS(ns, ns->residuals->u[dim])

        // Sync the distributed velocity property
        dp_sync(ns->dpu[dim]);
        // Printing the min and max velocity
        real velmin_global, velmax_global;
        MPI_Allreduce(&velmin, &velmin_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&velmax, &velmax_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        print0f("===> %d: Vmin = %15.10lf <===> Vmax = %15.10lf <===\n",dim,velmin_global,velmax_global);
    }
}

// Navier-Stokes pressure with variable density using the projection method
void higflow_pressure_multiphase(higflow_solver *ns) {
    // Get the local sub-domain for the cells
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    // Get the local sub-domain for the facets
    sim_facet_domain *sfdu[DIM];
    for(int dim = 0; dim < DIM; dim++) {
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
    }

    // Get the map for the domain properties
    //mp_mapper *mp = sd_get_domain_mapper(sdp);
    remove_pressure_singularity(ns, ns->slvp);
    // Loop for each cell
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Get the cell
        hig_cell *c = higcit_getcell(it);
        // Get the center of the cell
        Point ccenter;
        hig_get_center(c, ccenter);
        // Get the delta of the cell
        Point cdelta;
        hig_get_delta(c, cdelta);
        // Calculate the divergence of the intermediate velocity
        real sumdudx = 0.0;
        for(int dim = 0; dim < DIM; dim++) {
            int infacet;
            // Get the left facet intermediate velocity
            real ustarl = compute_facet_u_left(sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpustar[dim], ns->stn, &infacet);
            // Get the right facet intermediate velocity
            real ustarr = compute_facet_u_right(sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpustar[dim], ns->stn, &infacet);
            // Compute the derivative of intermediate velocity
            real dudx   = compute_facet_dudxc(cdelta, dim, 0.5, ustarl, ustarl, ustarr);
            // Compute the divergence of intermediate velocity
            sumdudx += dudx;
        }
        // Reset the stencil
        stn_reset(ns->stn);
        // Set the right side of stencil
        stn_set_rhs(ns->stn, sumdudx / ns->par.dt);
        // Calculate the point and weight of the stencil
        real alpha = 0.0;
        for(int dim = 0; dim < DIM; dim++) {
            // Stencil point update
            Point p;
            POINT_ASSIGN(p, ccenter);
            // Get the value at the point
            real rho    = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.dpdens, ns->ed.stn);
            // Get the density in the left cell
            real rhol   = compute_center_p_left(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.mult.dpdens, ns->ed.stn);
            // Get the density in the right cell
            real rhor   = compute_center_p_right(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.mult.dpdens, ns->ed.stn);
            // Stencil weight update
            real rho_ll = 2.0/(rho + rhol);
            real rho_rr = 2.0/(rho + rhor);
            // Stencil point update: right point
            real w      = rho_rr/(cdelta[dim]*cdelta[dim]);
            p[dim]      = ccenter[dim] + cdelta[dim];
            sd_get_stencil(sdp, ccenter, p, w, ns->stn);
            // Stencil point update: left point
            w           = rho_ll/(cdelta[dim]*cdelta[dim]);
            p[dim]      = ccenter[dim] - cdelta[dim];
            sd_get_stencil(sdp, ccenter, p, w, ns->stn);
            // Calculate the central coefficient
            alpha      -= (rho_rr + rho_ll)/(cdelta[dim]*cdelta[dim]);
        }
        // Get the stencil
        sd_get_stencil(sdp, ccenter, ccenter, alpha, ns->stn);
        // Get the index of the stencil
        int *ids   = psd_stn_get_gids(ns->psdp, ns->stn);
        // Get the value of the stencil
        real *vals = stn_get_vals(ns->stn);
        // Get the number of elements of the stencil
        int numelems = stn_get_numelems(ns->stn);
        // Get the cell identifier of the cell
        int cgid = psd_get_global_id(ns->psdp, c);
        if (ns->contr.desingpressure == true) { // cell with pressure desingularized
            if(cgid==ns->slvp->imposed_line) continue;
        }
        // Set the right side of solver linear system
        slv_set_bi(ns->slvp, cgid, stn_get_rhs(ns->stn));
        // Set the line of matrix of the solver linear system
        slv_set_Ai(ns->slvp, cgid, numelems, ids, vals);
    }
    // Destroy the iterator
    higcit_destroy(it);
    // Assemble the solver
    slv_assemble(ns->slvp);
    // Solve the linear system
    slv_solve(ns->slvp);
    // Set the solver solution in the pressure or pressure difference distributed property
    distributed_property *dp = (ns->contr.projtype == INCREMENTAL) ? ns->ddeltap : ns->dpp;
    dp_slv_load_from_solver(dp, ns->slvp);
    dp_sync(dp);
}

// *******************************************************************
// Navier-Stokes Step for the Explicit Euler Method
// *******************************************************************
void higflow_explicit_euler_intermediate_velocity_multiphase(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]) {
    // Get the facet iterator
    higfit_facetiterator *fit;
    // Get the local domain for cell
    sim_domain *sdp    = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Get the local domain for facet cell
    for (int dim2 = 0; dim2 < DIM; dim2++) {
        sfdu[dim2] = psfd_get_local_domain(ns->psfdu[dim2]);
    }
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++) {
        // Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Set the computational cell
            higflow_computational_cell_multiphase(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Interfacial force contribution
            rhs += higflow_interfacial_tension_term(ns);
            // Cell term contribution for the gravity
            if (dim == 1) rhs -= higflow_gravity_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Difusive term contribution
            rhs += higflow_difusive_term(ns, fdelta);
            // Compute the intermediate velocity
            real ustar = ns->cc.ufacet + ns->par.dt * rhs;
            // Update the distributed property intermediate velocity
            dp_set_value(dpustar[dim], flid, ustar);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Syncing the intermediate velocity
        dp_sync(dpustar[dim]);
    }
}

// *******************************************************************
// Navier-Stokes Step for the Second Order Explicit Runge-Kutta Method
// *******************************************************************
void higflow_explicit_runge_kutta_2_intermediate_velocity_multiphase(higflow_solver *ns) {
    // Calculate the auxiliar velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_multiphase(ns, ns->dpu, ns->dpuaux);
    // Calculate the star velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_multiphase(ns, ns->dpuaux, ns->dpustar);
    // Calculate the order 2 Runge-Kutta method using the euler method
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Loop for each dimension
    higfit_facetiterator *fit;
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local partitioned domain for facets
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the map of the distributed properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Get the intermediate velocity
            real u     = dp_get_value(ns->dpu[dim], flid);
            real ustar = dp_get_value(ns->dpustar[dim], flid);
            // Compute the intermediate velocity
            real urk2  = 0.5*(u + ustar);
            // Set the final velocity in the distributed velocity property
            dp_set_value(ns->dpustar[dim], flid, urk2);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the distributed velocity property
        dp_sync(ns->dpustar[dim]);
    }
}

// *******************************************************************
// Navier-Stokes Step for third Order Explicit Runge-Kutta Method
// *******************************************************************
void higflow_explicit_runge_kutta_3_intermediate_velocity_multiphase(higflow_solver *ns) {
    // Calculate the auxiliar velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_multiphase(ns, ns->dpu, ns->dpuaux);
    // Calculate the second stage velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_multiphase(ns, ns->dpuaux, ns->dpustar);
    // Calculate the order 2 Runge-Kutta method using the euler method
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Loop for each dimension
    higfit_facetiterator *fit;
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local partitioned domain for facets
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Get the intermediate velocity
            real u     = dp_get_value(ns->dpu[dim], flid);
            real ustar = dp_get_value(ns->dpustar[dim], flid);
            // Compute the intermediate velocity
            real urk3  = 0.75*u + 0.25*ustar;
            // Set the final velocity in the distributed velocity property
            dp_set_value(ns->dpuaux[dim], flid, urk3);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the distributed velocity property
        dp_sync(ns->dpuaux[dim]);
    }
    // Calculate the order 2 Runge-Kutta method using the euler method
    higflow_explicit_euler_intermediate_velocity_multiphase(ns, ns->dpuaux, ns->dpustar);
    // Loop for each dimension
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local partitioned domain for facets
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
    // Calculate the third stage velocity by the explicit euler method
        for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Get the intermediate velocity
            real u     = dp_get_value(ns->dpu[dim], flid);
            real ustar = dp_get_value(ns->dpustar[dim], flid);
            // Compute the intermediate velocity
            real urk3  = u/3.0 + 2.0*ustar/3.0;
            // Set the final velocity in the distributed velocity property
            dp_set_value(ns->dpustar[dim], flid, urk3);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the distributed velocity property
        dp_sync(ns->dpustar[dim]);
    }
}

// *******************************************************************
// Navier-Stokes Step for the Implicit Euler Method
// *******************************************************************
void higflow_semi_implicit_euler_intermediate_velocity_multiphase(higflow_solver *ns) {
    // Get the facet iterator
    higfit_facetiterator *fit;
    // Get the local domain for cell
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Get the local domain for facet cell
    for (int dim2 = 0; dim2 < DIM; dim2++) {
        sfdu[dim2] = psfd_get_local_domain(ns->psfdu[dim2]);
    }
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++) {
        // Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Set the computational cell
            higflow_computational_cell_multiphase(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Interfacial tension contribution
            rhs += higflow_interfacial_tension_term(ns);
            // Cell term contribution for the gravity
            if (dim == 1) rhs -= higflow_gravity_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Total contribuition terms by delta t
            rhs *= ns->par.dt;
            // Velocity term contribution
            rhs += ns->cc.ufacet;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // if(dim2==dim) {
                //     // Get the cell viscosity in the left cell
                //     ns->cc.viscl = compute_center_p_left(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.mult.dpvisc, ns->ed.stn);
                //     // Get the cell viscosity in the right cell
                //     ns->cc.viscr = compute_center_p_right(ns->ed.sdED, fcenter, fdelta, dim2, 0.5, ns->ed.mult.dpvisc, ns->ed.stn);
                // } else {
                //     Point p1,p2,p3,p4,p3_,p4_;
                //     // p1 and p2 in facet center
                //     POINT_ASSIGN(p1, fcenter);POINT_ASSIGN(p2, fcenter);
                //     // p1 and p2 in cell center
                //     p1[dim]=p1[dim]-0.5*fdelta[dim];p2[dim]=p2[dim]+0.5*fdelta[dim];
                //     // copy p1 and p2 in p3 and p4
                //     POINT_ASSIGN(p3, p1);POINT_ASSIGN(p4, p2);
                //     POINT_ASSIGN(p3_, p1);POINT_ASSIGN(p4_, p2);
                //     // p3 and p4 in cell center
                //     p3[dim2]=p3[dim2]+fdelta[dim2];p4[dim2]=p4[dim2]+fdelta[dim2];
                //     p3_[dim2]=p3_[dim2]-fdelta[dim2];p4_[dim2]=p4_[dim2]-fdelta[dim2];
                //     // viscosity 
                //     real v1=compute_value_at_point(ns->ed.sdED,fcenter,p1,1.0,ns->ed.mult.dpvisc,ns->ed.stn);
                //     real v2=compute_value_at_point(ns->ed.sdED,fcenter,p2,1.0,ns->ed.mult.dpvisc,ns->ed.stn);
                //     real v3=compute_value_at_point(ns->ed.sdED,fcenter,p3,1.0,ns->ed.mult.dpvisc,ns->ed.stn);
                //     real v4=compute_value_at_point(ns->ed.sdED,fcenter,p4,1.0,ns->ed.mult.dpvisc,ns->ed.stn);
                //     real v3_=compute_value_at_point(ns->ed.sdED,fcenter,p3_,1.0,ns->ed.mult.dpvisc,ns->ed.stn);
                //     real v4_=compute_value_at_point(ns->ed.sdED,fcenter,p4_,1.0,ns->ed.mult.dpvisc,ns->ed.stn);
                //     ns->cc.viscl=4.0/(1.0/v1+1.0/v2+1.0/v3_+1.0/v4_);
                //     ns->cc.viscr=4.0/(1.0/v1+1.0/v2+1.0/v3+1.0/v4);
                // }
                // Stencil weight update
                real wr = - ns->par.dt*ns->cc.viscr/(ns->par.Re*fdelta[dim2]*fdelta[dim2])/ns->cc.dens;
                real wl = - ns->par.dt*ns->cc.viscl/(ns->par.Re*fdelta[dim2]*fdelta[dim2])/ns->cc.dens;
                alpha -= (wr + wl);
                Point p;
                POINT_ASSIGN(p, fcenter);
                // Stencil point update: right point
                p[dim2] = fcenter[dim2] + fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, wr, ns->stn);
                // Stencil point update: left point
                p[dim2] = fcenter[dim2] - fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, wl, ns->stn);
            }
            alpha = 1.0 + alpha;
            // Get the stencil
            sfd_get_stencil(sfdu[dim], fcenter, fcenter, alpha, ns->stn);
            // Get the index of the stencil
            int *ids   = psfd_stn_get_gids(ns->psfdu[dim], ns->stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->stn);
            // Get the cell identifier of the cell
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Set the right side of solver linear system
            slv_set_bi(ns->slvu[dim], fgid, stn_get_rhs(ns->stn));
            // Set the line of matrix of the solver linear system
            slv_set_Ai(ns->slvu[dim], fgid, numelems, ids, vals);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Assemble the solver
        slv_assemble(ns->slvu[dim]);
        // Solve the linear system
        slv_solve(ns->slvu[dim]);
        // Get the solution of linear system
        // Gets the values of the solution
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Get the value of ustar
            real ustar = slv_get_xi(ns->slvu[dim], fgid);
            // Set the value of ustar
            dp_set_value(ns->dpustar[dim], flid, ustar);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Syncing the intermediate velocity
        dp_sync(ns->dpustar[dim]);
    }
}

// *******************************************************************
// Navier-Stokes Step for the Semi-Implicit Crank-Nicolson Method
// *******************************************************************
void higflow_semi_implicit_crank_nicolson_intermediate_velocity_multiphase(higflow_solver *ns) {
    // Get the facet iterator
    higfit_facetiterator *fit;
    // Get the local domain for cell
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Get the local domain for facet cell
    for (int dim2 = 0; dim2 < DIM; dim2++) {
        sfdu[dim2] = psfd_get_local_domain(ns->psfdu[dim2]);
    }
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++) {
        // Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Set the computational cell
            higflow_computational_cell_multiphase(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Diffusive term term contribution
            rhs += 0.5 * higflow_difusive_term(ns, fdelta);
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Total contribuition terms times delta t
            rhs *= ns->par.dt;
            // Velocity term contribution
            rhs += ns->cc.ufacet;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real wr = - 0.5*ns->par.dt*ns->cc.viscr/(ns->par.Re*fdelta[dim2]*fdelta[dim2])/ns->cc.dens;
                real wl = - 0.5*ns->par.dt*ns->cc.viscl/(ns->par.Re*fdelta[dim2]*fdelta[dim2])/ns->cc.dens;
                alpha -= (wr + wl);
                Point p;
                POINT_ASSIGN(p, fcenter);
                // Stencil point update: right point
                p[dim2] = fcenter[dim2] + fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, wr, ns->stn);
                // Stencil point update: left point
                p[dim2] = fcenter[dim2] - fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, wl, ns->stn);
            }
            alpha = 1.0 + alpha;
            // Get the stencil
            sfd_get_stencil(sfdu[dim], fcenter, fcenter, alpha, ns->stn);
            // Get the index of the stencil
            int *ids   = psfd_stn_get_gids(ns->psfdu[dim], ns->stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->stn);
            // Get the cell identifier of the cell
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Set the right side of solver linear system
            slv_set_bi(ns->slvu[dim], fgid, stn_get_rhs(ns->stn));
            // Set the line of matrix of the solver linear system
            slv_set_Ai(ns->slvu[dim], fgid, numelems, ids, vals);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Assemble the solver
        slv_assemble(ns->slvu[dim]);
        // Solve the linear system
        slv_solve(ns->slvu[dim]);
        // Gets the values of the solution
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Get the value of ustar
            real ustar = slv_get_xi(ns->slvu[dim], fgid);
            // Set the value of ustar
            dp_set_value(ns->dpustar[dim], flid, ustar);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Syncing the intermediate velocity
        dp_sync(ns->dpustar[dim]);
    }
}

// *******************************************************************
// Navier-Stokes Step for the Implicit BDF2 Method
// *******************************************************************
void higflow_semi_implicit_bdf2_intermediate_velocity_multiphase(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]) {
    // Firt stage of Tr-BDF2 method
    // Get the facet iterator
    higfit_facetiterator *fit;
    // Get the local domain for cell
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Get the local domain for facet cell
    for (int dim2 = 0; dim2 < DIM; dim2++) {
        sfdu[dim2] = psfd_get_local_domain(ns->psfdu[dim2]);
    }
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++) {
        // Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Set the computational cell
            higflow_computational_cell_multiphase(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Total contribuition terms times delta t
            rhs *= 0.5*ns->par.dt;
            // Difusive term contribution
            rhs += 0.25*ns->par.dt*higflow_difusive_term(ns, fdelta);
            // Velocity term contribution
            rhs += ns->cc.ufacet;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real wr = - 0.25*ns->par.dt*ns->cc.viscr/(ns->par.Re*fdelta[dim2]*fdelta[dim2])/ns->cc.dens;
                real wl = - 0.25*ns->par.dt*ns->cc.viscl/(ns->par.Re*fdelta[dim2]*fdelta[dim2])/ns->cc.dens;
                alpha -= (wr + wl);
                Point p;
                POINT_ASSIGN(p, fcenter);
                // Stencil point update: right point
                p[dim2] = fcenter[dim2] + fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, wr, ns->stn);
                // Stencil point update: left point
                p[dim2] = fcenter[dim2] - fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, wl, ns->stn);
            }
            alpha = 1.0 + alpha;
            // Get the stencil
            sfd_get_stencil(sfdu[dim], fcenter, fcenter, alpha, ns->stn);
            // Get the index of the stencil
            int *ids   = psfd_stn_get_gids(ns->psfdu[dim], ns->stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->stn);
            // Get the cell identifier of the cell
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Set the right side of solver linear system
            slv_set_bi(ns->slvu[dim], fgid, stn_get_rhs(ns->stn));
            // Set the line of matrix of the solver linear system
            slv_set_Ai(ns->slvu[dim], fgid, numelems, ids, vals);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Assemble the solver
        slv_assemble(ns->slvu[dim]);
        // Solve the linear system
        slv_solve(ns->slvu[dim]);
        // Gets the values of the solution
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Get the value of ustar
            real uaux = slv_get_xi(ns->slvu[dim], fgid);
            // Set the value of ustar
            dp_set_value(ns->dpuaux[dim], flid, uaux);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Syncing the intermediate velocity
        dp_sync(ns->dpuaux[dim]);
    }
    //Second Stage of Tr-BDF2
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++) {
        // Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Set the computational cell
            higflow_computational_cell_multiphase(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            //Get the uaux
            real uaux = dp_get_value(ns->dpuaux[dim], flid);
            // Right hand side equation
            real rhs = 0.0;
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Total contribuition terms times delta t
            rhs *= 1.0/3.0*ns->par.dt;
            rhs += (4.0*uaux - ns->cc.ufacet)/3.0;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real wr = - 1.0/3.0*ns->par.dt*ns->cc.viscr/(ns->par.Re*fdelta[dim2]*fdelta[dim2])/ns->cc.dens;
                real wl = - 1.0/3.0*ns->par.dt*ns->cc.viscl/(ns->par.Re*fdelta[dim2]*fdelta[dim2])/ns->cc.dens;
                alpha -= (wr + wl);
                Point p;
                POINT_ASSIGN(p, fcenter);
                // Stencil point update: right point
                p[dim2] = fcenter[dim2] + fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, wr, ns->stn);
                // Stencil point update: left point
                p[dim2] = fcenter[dim2] - fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, wl, ns->stn);
            }
            alpha = 1.0 + alpha;
            // Get the stencil
            sfd_get_stencil(sfdu[dim], fcenter, fcenter, alpha, ns->stn);
            // Get the index of the stencil
            int *ids   = psfd_stn_get_gids(ns->psfdu[dim], ns->stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->stn);
            // Get the cell identifier of the cell
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Set the right side of solver linear system
            slv_set_bi(ns->slvu[dim], fgid, stn_get_rhs(ns->stn));
            // Set the line of matrix of the solver linear system
            slv_set_Ai(ns->slvu[dim], fgid, numelems, ids, vals);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Assemble the solver
        slv_assemble(ns->slvu[dim]);
        // Solve the linear system
        slv_solve(ns->slvu[dim]);
        // Get the solution of linear system
        //Vec *vecu = slv_get_solution_vec(ns->slvu[dim]);
        // Gets the values of the solution
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Get the value of ustar
            real ustar = slv_get_xi(ns->slvu[dim], fgid);
            // Set the value of ustar
            dp_set_value(ns->dpustar[dim], flid, ustar);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Syncing the intermediate velocity
        dp_sync(ns->dpustar[dim]);
    }
}

// One step of the Navier-Stokes the projection method
void higflow_solver_step_multiphase(higflow_solver *ns) {
    // Boundary condition for velocity
    higflow_boundary_condition_for_velocity(ns);
    // Boundary conditions for source term
    higflow_boundary_condition_for_cell_source_term(ns);
    higflow_boundary_condition_for_facet_source_term(ns);
    // Calculate the source term
    higflow_calculate_source_term(ns);
    // Calculate the facet source term
    higflow_calculate_facet_source_term(ns);
    // Boundary condition for pressure
    higflow_boundary_condition_for_pressure(ns);
    // Calculate beta
    //higflow_compute_beta_visc_multiphase(ns); //used to be uncommented
    // Calculate S
    //higflow_compute_S_visc_multiphase(ns); //used to be uncommented
    // Calculate the viscosity
    higflow_compute_viscosity_multiphase(ns);
    // Calculate the viscosity
    higflow_compute_density_multiphase(ns);
    // Calculate the curvature
    higflow_compute_curvature_multiphase(ns);

    // Calculate the intermediated velocity
    switch (ns->contr.tempdiscrtype) {
        case EXPLICIT_EULER:
           // Explicit Euler method
           higflow_explicit_euler_intermediate_velocity_multiphase(ns, ns->dpu, ns->dpustar);
           break;
        case EXPLICIT_RK2: 
           // Explicit RK2 method
           higflow_explicit_runge_kutta_2_intermediate_velocity_multiphase(ns);
           break;
        case EXPLICIT_RK3: 
           // Explicit RK3 method
           higflow_explicit_runge_kutta_3_intermediate_velocity_multiphase(ns);
           break;
        case SEMI_IMPLICIT_EULER: 
           // Semi-Implicit Euler Method
           higflow_semi_implicit_euler_intermediate_velocity_multiphase(ns);
           break;
        case SEMI_IMPLICIT_CN: 
           // Semi-Implicit Crank-Nicolson Method
           higflow_semi_implicit_crank_nicolson_intermediate_velocity_multiphase(ns);
           break;
        case SEMI_IMPLICIT_BDF2: 
           // Semi-Implicit Crank-Nicolson Method
           higflow_semi_implicit_bdf2_intermediate_velocity_multiphase(ns, ns->dpu, ns->dpustar);
           break;
    }

    // Boundary condition for pressure
    higflow_boundary_condition_for_pressure(ns);
    // Calculate the pressure
    higflow_pressure_multiphase(ns);
    // Calculate the final velocity
    higflow_final_velocity_multiphase(ns);
    // Boundary condition for velocity
    higflow_boundary_condition_for_velocity(ns);
    // Calculate the final pressure
    higflow_final_pressure(ns);
    // if (ns->par.stepaux%1000==0) {
    //    printf("creating archives at step: %d\n",ns->par.stepaux);
    //    arquivoTempo(ns->par.stepaux);
    //    save_cell_values(ns,1);
    // }

    // Calculate the velocity derivative tensor
    /*higflow_compute_velocity_derivative_tensor(ns);
    // Computing the Kernel Tensor
    higflow_compute_kernel_tensor(ns);
    // Constitutive Equation Step for the Explicit Euler Method
    switch (ns->ed.ve.contr.discrtype) {
        case EXPLICIT:
        // Explicit method
         higflow_explicit_euler_constitutive_equation(ns);
         break;
        case IMPLICIT: 
           // Implicit method
           higflow_implicit_euler_constitutive_equation(ns);
           break;
    }
    // Computing the Polymeric Tensor
    higflow_compute_polymeric_tensor(ns);*/
    // Calculate the volume fraction
    higflow_plic_advection_volume_fraction(ns);
    //higflow_explicit_euler_volume_fraction(ns);
}
