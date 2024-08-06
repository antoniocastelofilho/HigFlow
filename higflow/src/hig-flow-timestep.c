// *******************************************************************
//  HiG-Flow Solver - version 01/08/24
// *******************************************************************

#include "hig-flow-timestep.h"

real higflow_compute_CFL(higflow_solver *ns) {
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    higcit_celliterator *cit; 
    hig_cell *c; Point ccenter; Point cdelta;
    real ul, ur, u, u_dx_norm, maxu_dx_norm = 0.0;
    real facet_area;
    int infacet;
    for(cit = sd_get_domain_celliterator(sdp); !higcit_isfinished(cit); higcit_nextcell(cit)) {
        c = higcit_getcell(cit);
        hig_get_center(c, ccenter);
        hig_get_delta(c, cdelta);
        u_dx_norm = 0.0;
        for(int dim = 0; dim < DIM; dim++) {
            ul = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            ur = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            u  = (ul + ur) / 2.0;
            facet_area = 1.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                if(dim2 != dim) facet_area *= cdelta[dim2];
            }
            u_dx_norm += fabs(u) / facet_area;
        }
        if(u_dx_norm > maxu_dx_norm) maxu_dx_norm = u_dx_norm;
    }
    higcit_destroy(cit);

    real dt = ns->par.dt;
    real cfl = maxu_dx_norm * dt;

    real cfl_global;
    MPI_Allreduce(&cfl, &cfl_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return cfl_global;
}


real higflow_compute_CFL_PNP(higflow_solver *ns) {
    sim_domain *sdpsi = psd_get_local_domain(ns->ed.eo.psdEOpsi);
    higcit_celliterator *cit; 
    hig_cell *c; Point ccenter; Point cdelta;
    real psil, psir, phil, phir, dpsidx, dphidx, E_dx_norm, maxE_dx_norm = 0.0;
    real facet_area;
    real alphaeo = ns->ed.eo.par.alpha; real Pe = ns->ed.eo.par.Pe;
    for(cit = sd_get_domain_celliterator(sdpsi); !higcit_isfinished(cit); higcit_nextcell(cit)) {
        c = higcit_getcell(cit);
        hig_get_center(c, ccenter);
        hig_get_delta(c, cdelta);
        maxE_dx_norm = 0.0;
        for(int dim = 0; dim < DIM; dim++) {
            psil = compute_center_p_left(ns->ed.eo.sdEOpsi, ccenter, cdelta, dim, 1.0, ns->ed.eo.dppsi, ns->ed.eo.stnpsi);
            psir = compute_center_p_right(ns->ed.eo.sdEOpsi, ccenter, cdelta, dim, 1.0, ns->ed.eo.dppsi, ns->ed.eo.stnpsi);
            phil = compute_center_p_left(ns->ed.eo.sdEOphi, ccenter, cdelta, dim, 1.0, ns->ed.eo.dpphi, ns->ed.eo.stnphi);
            phir = compute_center_p_right(ns->ed.eo.sdEOphi, ccenter, cdelta, dim, 1.0, ns->ed.eo.dpphi, ns->ed.eo.stnphi);
            dphidx = (phir - phil) / cdelta[dim];
            dpsidx = (psir - psil) / cdelta[dim];
            facet_area = 1.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                if(dim2 != dim) facet_area *= cdelta[dim2];
            }
            E_dx_norm += fabs(dphidx + dpsidx) / facet_area;
        }
        E_dx_norm *= alphaeo / Pe;
        if(E_dx_norm > maxE_dx_norm) maxE_dx_norm = E_dx_norm;
    }
    higcit_destroy(cit);

    real dt = ns->par.dt;
    real cfl = maxE_dx_norm * dt;

    real cfl_global;
    MPI_Allreduce(&cfl, &cfl_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return cfl_global;
}

real higflow_compute_CFL_PNP_multiphase(higflow_solver *ns) {
    sim_domain *sdpsi = psd_get_local_domain(ns->ed.eo.psdEOpsi);
    higcit_celliterator *cit; 
    hig_cell *c; Point ccenter; Point cdelta;
    real psil, psir, phil, phir, dpsidx, dphidx, E_dx_norm, maxE_dx_norm = 0.0;
    real facet_area, fracvol;
    real alphaeo0 = ns->ed.mult.eo.par0.alpha; real Pe0 = ns->ed.mult.eo.par0.Pe;
    real alphaeo1 = ns->ed.mult.eo.par1.alpha; real Pe1 = ns->ed.mult.eo.par1.Pe;
    real alphaeo, Pe;
    for(cit = sd_get_domain_celliterator(sdpsi); !higcit_isfinished(cit); higcit_nextcell(cit)) {
        c = higcit_getcell(cit);
        hig_get_center(c, ccenter);
        hig_get_delta(c, cdelta);
        fracvol = compute_value_at_point(ns->ed.mult.sdmult, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
        alphaeo = (1.0 - fracvol) * alphaeo0 + fracvol * alphaeo1;
        Pe = (1.0 - fracvol) * Pe0 + fracvol * Pe1;
        maxE_dx_norm = 0.0;
        for(int dim = 0; dim < DIM; dim++) {
            psil = compute_center_p_left(ns->ed.eo.sdEOpsi, ccenter, cdelta, dim, 1.0, ns->ed.eo.dppsi, ns->ed.eo.stnpsi);
            psir = compute_center_p_right(ns->ed.eo.sdEOpsi, ccenter, cdelta, dim, 1.0, ns->ed.eo.dppsi, ns->ed.eo.stnpsi);
            phil = compute_center_p_left(ns->ed.eo.sdEOphi, ccenter, cdelta, dim, 1.0, ns->ed.eo.dpphi, ns->ed.eo.stnphi);
            phir = compute_center_p_right(ns->ed.eo.sdEOphi, ccenter, cdelta, dim, 1.0, ns->ed.eo.dpphi, ns->ed.eo.stnphi);
            dphidx = (phir - phil) / cdelta[dim];
            dpsidx = (psir - psil) / cdelta[dim];
            facet_area = 1.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                if(dim2 != dim) facet_area *= cdelta[dim2];
            }
            E_dx_norm += fabs(dphidx + dpsidx) / facet_area;
        }
        E_dx_norm *= alphaeo / Pe;
        if(E_dx_norm > maxE_dx_norm) maxE_dx_norm = E_dx_norm;
    }
    higcit_destroy(cit);

    real dt = ns->par.dt;
    real cfl = maxE_dx_norm * dt;

    real cfl_global;
    MPI_Allreduce(&cfl, &cfl_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return cfl_global;
}

void higflow_adjust_timestep(higflow_solver *ns) {
    real dt0 = ns->par.dt0;
    // real dtmin = dt0/32.0;
    real dtmax = dt0*32.0;

    real cfl = higflow_compute_CFL(ns);
    if(ns->contr.eoflow == true){
        if(ns->ed.eo.contr.eo_model == PNP) cfl += higflow_compute_CFL_PNP(ns);
    } else if(ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true) {
        if(ns->ed.mult.eo.contr.eo_model == PNP) cfl += higflow_compute_CFL_PNP_multiphase(ns);
    }
    print0f("CFL = %f\n", cfl);

    real dt_factor = min(MAX_CFL / (cfl+EPSMACH), MAX_DT_FACTOR);
    real newdt = dt_factor*ns->par.dt;
    ns->par.dt = min(newdt, dtmax);
    print0f("dt = %f\n", ns->par.dt);
}