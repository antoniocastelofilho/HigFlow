// *******************************************************************
//  HiG-Flow Solver IO - version 04/05/2021
// *******************************************************************

#include "hig-flow-io.h"
#include <string.h>
#include <libfyaml.h>

// *******************************************************************
// Navier-Stokes Print for Visualize
// *******************************************************************

// Print the electro-osmotic proporties
void higflow_print_electroosmotic_properties(higflow_solver *ns, FILE *data, int dimprint, real pprint) {
    if (ns->contr.eoflow == true) {
        // Get the cosntants
        real epsprint = 1.0e-8;
        real deltaeo = ns->ed.eo.par.delta;
        // Necessary properties
        real phi, psi, np, nm, rhoe;
        // Get the local sub-domain for the cells
        sim_domain *sdphi   = psd_get_local_domain(ns->ed.eo.psdEOphi);
        sim_domain *sdpsi   = psd_get_local_domain(ns->ed.eo.psdEOpsi);
        sim_domain *sdnplus = psd_get_local_domain(ns->ed.eo.psdEOnplus);
        sim_domain *sdnminus = psd_get_local_domain(ns->ed.eo.psdEOnminus);
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdpsi);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdpsi); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid    = mp_lookup(mp, hig_get_cid(c));
            // Get the inside/outside inflow point cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Get the electro-osmotic properties 
            phi  = compute_value_at_point(ns->ed.eo.sdEOphi, ccenter, ccenter, 1.0, ns->ed.eo.dpphi, ns->ed.stn);
            psi  = compute_value_at_point(ns->ed.eo.sdEOpsi, ccenter, ccenter, 1.0, ns->ed.eo.dppsi, ns->ed.stn);
            np   = compute_value_at_point(ns->ed.eo.sdEOnplus, ccenter, ccenter, 1.0, ns->ed.eo.dpnplus, ns->ed.stn);
            nm   = compute_value_at_point(ns->ed.eo.sdEOnminus, ccenter, ccenter, 1.0, ns->ed.eo.dpnminus, ns->ed.stn);
            rhoe = deltaeo*(np - nm);
            //Print data file
            if ((ccenter[dimprint] < pprint + epsprint) && (ccenter[dimprint] > pprint - epsprint)){
                   fprintf(data, "%15.10lf  %15.10lf  %15.12lf  %15.12lf  %15.12lf  %15.12lf  %15.12lf\n", ccenter[0], ccenter[1], phi, psi, np, nm, rhoe);
            }
        }
        // Destroy the iterator
        higcit_destroy(it);
    }
}

// Print the Polymeric Tensor
void higflow_print_polymeric_tensor(higflow_solver *ns, FILE *data, int dimprint, real pprint) {
    if (ns->contr.flowtype == VISCOELASTIC) {
        // Get the cosntants
        real Re   = ns->par.Re;
        real De   = ns->ed.ve.par.De;
        real beta = ns->ed.ve.par.beta;
        real tol  = ns->ed.ve.par.kernel_tol;
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
            // Get the inside/outside inflow point cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Get the velocity derivative tensor Du and the Kernel tensor
            real Kernel[DIM][DIM], Du[DIM][DIM], S[DIM][DIM], D[DIM][DIM];
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get Du
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                    D[i][j]  = 0.5*(Du[i][j]+Du[j][i]);
                    // Get S tensor
                    S[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                }
            }
            //Print polymeric stress data file
            if (ccenter[dimprint] == pprint){
                   fprintf(data, "%lf  %lf  %15.12lf  %15.12lf  %15.12lf  %15.12lf\n", ccenter[0], ccenter[1], S[0][0] + 2.0 * (1 - beta) * D[0][0] / Re, S[0][1] + 2.0 * (1 - beta) * D[0][1] / Re, S[1][0] + 2.0 * (1 - beta) * D[1][0] / Re, S[1][1] + 2.0 * (1 - beta) * D[1][1] / Re);
            }
        }
        // Destroy the iterator
        higcit_destroy(it);
    }
}

// Print the velocity
void higflow_print_velocity(higflow_solver *ns, FILE *data, int dimprint, real pprint) {
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Loop for each dimension
    higfit_facetiterator *fit;
    for (int dim = 0; dim < DIM; dim++) {
        // Initialize the min and max velocity
        // Get the local partitioned domain for facets
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the map of the distributed properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Get the velocity
            real u = dp_get_value(ns->dpu[dim], flid);
            if(fcenter[dimprint] == pprint){
                fprintf(data, "%15.06lf  %15.06lf  %15.10lf\n", fcenter[0], fcenter[1], u);
            }
        }
        // Destroy the iterator
        higfit_destroy(fit);
    }
}

// Print the VTK file for visualize
void higflow_print_vtk(higflow_solver *ns, int rank) {
    switch (DIM) {
        case 2:
            // 2D case
            higflow_print_vtk2D(ns, rank);
            break;
        case 3:
            // 3D case
             higflow_print_vtk3D(ns, rank);
        break;
    }
}

// Print the kinetic energy
void higflow_kinetic_energy(higflow_solver *ns, FILE *data) {
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
    real volume  = 0.0;
    real kinetic = 0.0;
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
        // 
        real norm2vel = 0.0;
        real volcel   = 1.0;
        for(int dim = 0; dim < DIM; dim++) {
            //Get the velocity in the left facet center
            real ul = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            //Get the velocity in the right facet center
            real ur = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            real u = 0.5*(ur + ul);
            norm2vel += u*u;
            volcel   *= cdelta[dim];
        }
        volume  += volcel;
        kinetic += norm2vel*volcel;
    }
    kinetic = 0.5*kinetic/volume;
    // Destroy the iterator
    higcit_destroy(it);
    // Printing the min and max velocity
    printf("<===> t = %.12lf <===> kinetic = %.12lf <===\n",ns->par.t,kinetic);
    fprintf(data,"%.12lf  %.12lf \n",ns->par.t,kinetic);
}

// Print the VTK file for visualize 2D
void higflow_print_vtk2D(higflow_solver *ns, int rank) {
    //  Open the VTK file
    real Re   = ns->par.Re;
    real beta = ns->ed.ve.par.beta;
    char vtkname[1024];
    snprintf(vtkname, sizeof vtkname, "%s_%d-%d.vtk", ns->par.nameprint, rank, ns->par.frame);
    FILE *f = fopen(vtkname, "w");
    if (f == NULL) {
        return;
    }
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    int numhigs = sd_get_num_higtrees(sdp);
    higcit_celliterator *it;
    // it = higcit_create_all_leaves(root);
    it = sd_get_domain_celliterator(sdp);
    long numleafs = higcit_count_without_advancing(it);
    fprintf(f, "# vtk DataFile Version 3.0\n");
    fprintf(f, "higtree\n");
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(f, "\nPOINTS %ld float\n", 4*numleafs);
    hig_cell *c;
    for (; !higcit_isfinished(it); higcit_nextcell(it)) {
        c = higcit_getcell(it);
        fprintf(f, "%e %e 0\n%e %e 0\n%e %e 0\n%e %e 0\n", 
        c->lowpoint[0], c->lowpoint[1], c->highpoint[0], c->lowpoint[1], 
        c->highpoint[0], c->highpoint[1], c->lowpoint[0], c->highpoint[1]);
    }
    higcit_destroy(it);
    fprintf(f, "\nCELLS %ld %ld\n", numleafs, 5*numleafs);
    for (int i = 0; i < numleafs; i++) {
        fprintf(f, "%d %d %d %d %d\n", 4, 4*i, 4*i+1, 4*i+2, 4*i+3);
    }
    fprintf(f, "\nCELL_TYPES %ld\n", numleafs);
    for (int i = 0; i < numleafs; i++) {
        fprintf(f, "7 "); //vtk cell types 7 = VTK_POLYGON // 9
    }
    
    // higcit_celliterator *it = higcit_create_all_leaves(root);
    // higcit_celliterator *it = sd_get_domain_celliterator(sdp);
    // long numcells = higcit_count(it);
    //printf("\n ---- Number of cells: %lu\n", numcells);
    // higcit_destroy(it);
    // sim_domain *sdp = psd_get_local_domain(psdp);

    // int degree = 2;
    // const int maxpts = 2*DIM*wls_num_min_points(DIM,degree);
    // //const int maxpts = 2*DIM*wls_num_min_points(degree);
    // Point pts[maxpts+1];
    // real dist[maxpts+1];
    // real w0[maxpts], w1[maxpts], w2[maxpts], w3[maxpts], w4[maxpts], w5[maxpts], w6[maxpts], w7[maxpts];
    // uniqueid gids[maxpts];
    // Saving vector properties of cell faces
    sim_facet_domain *sfdu2[DIM];
    sfdu2[0] = psfd_get_local_domain(ns->psfdu[0]); sfdu2[1] = psfd_get_local_domain(ns->psfdu[1]);
    distributed_property *dpu[DIM];
    dpu[0] = ns->dpu[0]; dpu[1] = ns->dpu[1];

    fprintf(f, "\n\nPOINT_DATA %ld\nVECTORS vel FLOAT\n", 4*numleafs);

    // Saving vector properties of cell faces
    for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        Point ccenter;
        hig_get_center(c,ccenter);
        // Pontos onde será interpolada a velocidade
        Point p0, p1, p2, p3;
        p0[0] = c->lowpoint[0];  p0[1] = c->lowpoint[1];
        p1[0] = c->highpoint[0]; p1[1] = c->lowpoint[1];
        p2[0] = c->highpoint[0]; p2[1] = c->highpoint[1];
        p3[0] = c->lowpoint[0];  p3[1] = c->highpoint[1];
        Point vel0, vel1, vel2, vel3;

        for(int dim = 0; dim < DIM; dim++) {
            vel0[dim]= compute_facet_value_at_point(sfdu2[dim], ccenter, p0, 1.0, dpu[dim], ns->stn);
            vel1[dim] = compute_facet_value_at_point(sfdu2[dim], ccenter, p1, 1.0, dpu[dim], ns->stn);
            vel2[dim] = compute_facet_value_at_point(sfdu2[dim], ccenter, p2, 1.0, dpu[dim], ns->stn);
            vel3[dim] = compute_facet_value_at_point(sfdu2[dim], ccenter, p3, 1.0, dpu[dim], ns->stn);
        }
        
        fprintf(f, "%e %e 0\n%e %e 0\n%e %e 0\n%e %e 0\n", 
        vel0[0], vel0[1], vel1[0], vel1[1], vel2[0], vel2[1], vel3[0], vel3[1]);
    }
    higcit_destroy(it);

    /*
    switch (ns->contr.eoflow) {
        case false:
        break;

        case true:
            sim_facet_domain *sfdEOFeo2[DIM];
            sfdEOFeo2[0] = psfd_get_local_domain(ns->ed.eo.psfdEOFeo[0]); sfdEOFeo2[1] = psfd_get_local_domain(ns->ed.eo.psfdEOFeo[1]);
            distributed_property *dpEOFeo2[DIM];
            dpEOFeo2[0] = ns->ed.eo.dpFeo[0]; dpEOFeo2[1] = ns->ed.eo.dpFeo[1];
            Point p0, p1, p2, p3;
            fprintf(f, "VECTORS F\u2091\u2092 FLOAT\n");
            for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                hig_cell *c = higcit_getcell(it);
                Point ccenter;
                hig_get_center(c,ccenter);
                // Pontos onde será interpolada a força eletrica
                p0[0] = c->lowpoint[0];  p0[1] = c->lowpoint[1];
                p1[0] = c->highpoint[0]; p1[1] = c->lowpoint[1];
                p2[0] = c->highpoint[0]; p2[1] = c->highpoint[1];
                p3[0] = c->lowpoint[0];  p3[1] = c->highpoint[1];
                Point Feo0, Feo1, Feo2, Feo3;

                for(int dim = 0; dim < DIM; dim++) {
                    Feo0[dim] = compute_facet_value_at_point(sfdEOFeo2[dim], ccenter, p0, 1.0, ns->ed.eo.dpFeo[dim], ns->ed.stn);
                    Feo1[dim] = compute_facet_value_at_point(sfdEOFeo2[dim], ccenter, p1, 1.0, ns->ed.eo.dpFeo[dim], ns->ed.stn);
                    Feo2[dim] = compute_facet_value_at_point(sfdEOFeo2[dim], ccenter, p2, 1.0, ns->ed.eo.dpFeo[dim], ns->ed.stn);
                    Feo3[dim] = compute_facet_value_at_point(sfdEOFeo2[dim], ccenter, p3, 1.0, ns->ed.eo.dpFeo[dim], ns->ed.stn);
                }

                fprintf(f, "%e %e 0\n%e %e 0\n%e %e 0\n%e %e 0\n", 
                Feo0[0], Feo0[1], Feo1[0], Feo1[1], Feo2[0], Feo2[1], Feo3[0], Feo3[1]);
            }
            higcit_destroy(it);
        break;
    }
    */

    // Saving tensorial properties 
    switch (ns->contr.flowtype) {
        case NEWTONIAN:
        break;
        case GENERALIZED_NEWTONIAN:
        break;
        case MULTIPHASE:
            if(ns->ed.mult.contr.flowtype_either == VISCOELASTIC) {
                real beta0 = ns->ed.mult.ve.par0.beta;
                real beta1 = ns->ed.mult.ve.par1.beta;
                real t = ns->par.t;
                real beta_interp;
                real fracvol;
                real visc, visc0, visc1;

                fprintf(f, "\nTENSORS \u03C4\u209A FLOAT\n");
                for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                    hig_cell *c = higcit_getcell(it);
                    Point ccenter;
                    hig_get_center(c,ccenter);

                    // Pontos onde será interpolada a velocidade
                    Point p0, p1, p2, p3;
                    p0[0] = c->lowpoint[0];  p0[1] = c->lowpoint[1];
                    p1[0] = c->highpoint[0]; p1[1] = c->lowpoint[1];
                    p2[0] = c->highpoint[0]; p2[1] = c->highpoint[1];
                    p3[0] = c->lowpoint[0];  p3[1] = c->highpoint[1];

                    real taup0[DIM+1][DIM+1], taup1[DIM+1][DIM+1], taup2[DIM+1][DIM+1], taup3[DIM+1][DIM+1];
                    real Dp0[DIM][DIM], Dp1[DIM][DIM], Dp2[DIM][DIM], Dp3[DIM][DIM];
                    for (int i = 0; i < DIM; i++) {
                        for (int j = 0; j < DIM; j++) {
                            // Get Du
                            Dp0[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p0, 1.0, ns->ed.mult.ve.dpD[i][j], ns->ed.stn);
                            Dp1[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p1, 1.0, ns->ed.mult.ve.dpD[i][j], ns->ed.stn);
                            Dp2[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p2, 1.0, ns->ed.mult.ve.dpD[i][j], ns->ed.stn);
                            Dp3[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p3, 1.0, ns->ed.mult.ve.dpD[i][j], ns->ed.stn);
                        }
                    }
                        
                    for (int i = 0; i < DIM; i++) {
                        for (int j = 0; j < DIM; j++) {
                            
                            fracvol = compute_value_at_point(sdp, ccenter, p0, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
                            visc0 = ns->ed.mult.get_density0(p0, t);
                            visc1 = ns->ed.mult.get_density1(p0, t);
                            visc = (1.0 - fracvol) * visc0 + fracvol * visc1;
                            beta_interp =  ((1 - fracvol) * visc0 * beta0 + fracvol * visc1 * beta1) / visc;
                            taup0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.mult.ve.dpS[i][j], ns->ed.stn);
                            taup0[i][j]+= 2.0*(1-beta_interp)*0.5*(Dp0[i][j]+Dp0[j][i])/Re;

                            fracvol = compute_value_at_point(sdp, ccenter, p1, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
                            visc0 = ns->ed.mult.get_density0(p1, t);
                            visc1 = ns->ed.mult.get_density1(p1, t);
                            visc = (1.0 - fracvol) * visc0 + fracvol * visc1;
                            beta_interp =  ((1 - fracvol) * visc0 * beta0 + fracvol * visc1 * beta1) / visc;
                            taup1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.mult.ve.dpS[i][j], ns->ed.stn);
                            taup1[i][j]+= 2.0*(1-beta_interp)*0.5*(Dp1[i][j]+Dp1[j][i])/Re;

                            fracvol = compute_value_at_point(sdp, ccenter, p2, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
                            visc0 = ns->ed.mult.get_density0(p2, t);
                            visc1 = ns->ed.mult.get_density1(p2, t);
                            visc = (1.0 - fracvol) * visc0 + fracvol * visc1;
                            beta_interp =  ((1 - fracvol) * visc0 * beta0 + fracvol * visc1 * beta1) / visc;
                            taup2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.mult.ve.dpS[i][j], ns->ed.stn);
                            taup2[i][j]+= 2.0*(1-beta_interp)*0.5*(Dp2[i][j]+Dp2[j][i])/Re;
                            
                            fracvol = compute_value_at_point(sdp, ccenter, p3, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
                            visc0 = ns->ed.mult.get_density0(p3, t);
                            visc1 = ns->ed.mult.get_density1(p3, t);
                            visc = (1.0 - fracvol) * visc0 + fracvol * visc1;
                            beta_interp =  ((1 - fracvol) * visc0 * beta0 + fracvol * visc1 * beta1) / visc;
                            taup3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.mult.ve.dpS[i][j], ns->ed.stn);
                            taup3[i][j]+= 2.0*(1-beta_interp)*0.5*(Dp3[i][j]+Dp3[j][i])/Re;
                        } 
                    }  

                    fprintf(f, "%e %e 0\n%e %e 0\n0 0 0\n%e %e 0\n%e %e 0\n0 0 0\n%e %e 0\n%e %e 0\n0 0 0\n%e %e 0\n%e %e 0\n0 0 0\n",
                    taup0[0][0], taup0[0][1], taup0[1][0], taup0[1][1],
                    taup1[0][0], taup1[0][1], taup1[1][0], taup1[1][1],
                    taup2[0][0], taup2[0][1], taup2[1][0], taup2[1][1],
                    taup3[0][0], taup3[0][1], taup3[1][0], taup3[1][1]);
                }
                    
            }
        break;
        case VISCOELASTIC:
           fprintf(f, "\nTENSORS \u03C4\u209A FLOAT\n");
           for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
              hig_cell *c = higcit_getcell(it);
              Point ccenter;
              hig_get_center(c,ccenter);
              // Pontos onde será interpolada a velocidade
              Point p0, p1, p2, p3;
              p0[0] = c->lowpoint[0];  p0[1] = c->lowpoint[1];
              p1[0] = c->highpoint[0]; p1[1] = c->lowpoint[1];
              p2[0] = c->highpoint[0]; p2[1] = c->highpoint[1];
              p3[0] = c->lowpoint[0];  p3[1] = c->highpoint[1];
              
                real taup0[DIM+1][DIM+1], taup1[DIM+1][DIM+1], taup2[DIM+1][DIM+1], taup3[DIM+1][DIM+1];
                real Dp0[DIM][DIM], Dp1[DIM][DIM], Dp2[DIM][DIM], Dp3[DIM][DIM];
                for (int i = 0; i < DIM; i++) {
                    for (int j = 0; j < DIM; j++) {
                        // Get Du
                        Dp0[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p0, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                        Dp1[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p1, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                        Dp2[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p2, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                        Dp3[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p3, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                    }
                }
                    
                for (int i = 0; i < DIM; i++) {
                    for (int j = 0; j < DIM; j++) {
                        taup0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                        taup0[i][j]+= 2.0*(1-beta)*0.5*(Dp0[i][j]+Dp0[j][i])/Re;
                        taup1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                        taup1[i][j]+= 2.0*(1-beta)*0.5*(Dp1[i][j]+Dp1[j][i])/Re;
                        taup2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                        taup2[i][j]+= 2.0*(1-beta)*0.5*(Dp2[i][j]+Dp2[j][i])/Re;
                        taup3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                        taup3[i][j]+= 2.0*(1-beta)*0.5*(Dp3[i][j]+Dp3[j][i])/Re;
                    } 
                }  

                fprintf(f, "%e %e 0\n%e %e 0\n0 0 0\n%e %e 0\n%e %e 0\n0 0 0\n%e %e 0\n%e %e 0\n0 0 0\n%e %e 0\n%e %e 0\n0 0 0\n",
                taup0[0][0], taup0[0][1], taup0[1][0], taup0[1][1],
                taup1[0][0], taup1[0][1], taup1[1][0], taup1[1][1],
                taup2[0][0], taup2[0][1], taup2[1][0], taup2[1][1],
                taup3[0][0], taup3[0][1], taup3[1][0], taup3[1][1]);
            }
            higcit_destroy(it);
        break;
        
        case VISCOELASTIC_INTEGRAL:
        // Integral
            fprintf(f, "\nTENSORS \u03C4\u209A FLOAT\n");
            //real Re   = ns->par.Re;
            for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                hig_cell *c = higcit_getcell(it);
                Point ccenter;
                hig_get_center(c,ccenter);
                // Pontos onde será interpolada a velocidade
                Point p0, p1, p2, p3;
                p0[0] = c->lowpoint[0];  p0[1] = c->lowpoint[1];
                p1[0] = c->highpoint[0]; p1[1] = c->lowpoint[1];
                p2[0] = c->highpoint[0]; p2[1] = c->highpoint[1];
                p3[0] = c->lowpoint[0];  p3[1] = c->highpoint[1];
                
                real taup0[DIM+1][DIM+1], taup1[DIM+1][DIM+1], taup2[DIM+1][DIM+1], taup3[DIM+1][DIM+1];
                real Dp0[DIM][DIM], Dp1[DIM][DIM], Dp2[DIM][DIM], Dp3[DIM][DIM];
                for (int i = 0; i < DIM; i++) {
                    for (int j = 0; j < DIM; j++) {
                        // Get Du
                        Dp0[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p0, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
                        Dp1[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p1, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
                        Dp2[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p2, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
                        Dp3[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p3, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
                    }
                }
                    
                for (int i = 0; i < DIM; i++) {
                    for (int j = 0; j < DIM; j++) {
                        taup0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
                        taup0[i][j]+= 2.0*(1-beta)*0.5*(Dp0[i][j]+Dp0[j][i])/Re;
                        taup1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
                        taup1[i][j]+= 2.0*(1-beta)*0.5*(Dp1[i][j]+Dp1[j][i])/Re;
                        taup2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
                        taup2[i][j]+= 2.0*(1-beta)*0.5*(Dp2[i][j]+Dp2[j][i])/Re;
                        taup3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
                        taup3[i][j]+= 2.0*(1-beta)*0.5*(Dp3[i][j]+Dp3[j][i])/Re;
                    } 
                }
                
                fprintf(f, "%e %e 0\n%e %e 0\n0 0 0\n%e %e 0\n%e %e 0\n0 0 0\n%e %e 0\n%e %e 0\n0 0 0\n%e %e 0\n%e %e 0\n0 0 0\n",
                taup0[0][0], taup0[0][1], taup0[1][0], taup0[1][1],
                taup1[0][0], taup1[0][1], taup1[1][0], taup1[1][1],
                taup2[0][0], taup2[0][1], taup2[1][0], taup2[1][1],
                taup3[0][0], taup3[0][1], taup3[1][0], taup3[1][1]);
            }
            higcit_destroy(it);
        break;
    }

    Point ccenter;
    // Saving pressure from the center of the cell 
    fprintf(f, "\nCELL_DATA %ld\nSCALARS p FLOAT\nLOOKUP_TABLE default\n", numleafs);
    for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        hig_get_center(c,ccenter);
        real val  = compute_value_at_point(ns->sdp, ccenter, ccenter, 1.0, ns->dpp, ns->stn);
        fprintf(f, "%e\n", val);
    }
    higcit_destroy(it);
    
    // Saving scalar properties from the center of the cell 
    switch (ns->contr.flowtype) {
        case NEWTONIAN:
        break;
        
        case GENERALIZED_NEWTONIAN:
        break;
        
        case MULTIPHASE:
            fprintf(f,"\nSCALARS FracVol FLOAT\nLOOKUP_TABLE default\n");
            sfdu2[0] = psfd_get_local_domain(ns->psfdu[0]);
            for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                hig_cell *c = higcit_getcell(it);
                hig_get_center(c, ccenter);
                real value  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
               fprintf(f, "%e\n", value);
            }
            higcit_destroy(it);
         break;
            
        case VISCOELASTIC:
        break;
        
        case VISCOELASTIC_INTEGRAL:
        break;
    }

    switch (ns->contr.eoflow) {
        case false:
        break;

        case true: ;
            // Get the constant
            real deltaeo = ns->ed.eo.par.delta;
            // Necessary properties
            real phi, psi, np, nm, rhoe;
            // loop variables
            Point ccenter;
            hig_cell *c;
            fprintf(f, "\nSCALARS \u03D5 FLOAT\nLOOKUP_TABLE default\n");
            for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                c = higcit_getcell(it);
                hig_get_center(c, ccenter);
                phi  = compute_value_at_point(ns->ed.eo.sdEOphi, ccenter, ccenter, 1.0, ns->ed.eo.dpphi, ns->ed.stn);
                fprintf(f, "%e\n", phi);
            }
            higcit_destroy(it);
            fprintf(f, "\nSCALARS \u03C8 FLOAT\nLOOKUP_TABLE default\n");
            for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                c = higcit_getcell(it);
                hig_get_center(c, ccenter);
                psi  = compute_value_at_point(ns->ed.eo.sdEOpsi, ccenter, ccenter, 1.0, ns->ed.eo.dppsi, ns->ed.stn);
                fprintf(f, "%e\n", psi);
            }
            higcit_destroy(it);
            fprintf(f, "\nSCALARS n\u207A FLOAT\nLOOKUP_TABLE default\n");
            for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                c = higcit_getcell(it);
                hig_get_center(c, ccenter);
                np   = compute_value_at_point(ns->ed.eo.sdEOnplus, ccenter, ccenter, 1.0, ns->ed.eo.dpnplus, ns->ed.stn);
                fprintf(f, "%e\n", np);
            }
            higcit_destroy(it);
            fprintf(f, "\nSCALARS n\u207B FLOAT\nLOOKUP_TABLE default\n");
            for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                c = higcit_getcell(it);
                hig_get_center(c, ccenter);
                nm   = compute_value_at_point(ns->ed.eo.sdEOnminus, ccenter, ccenter, 1.0, ns->ed.eo.dpnminus, ns->ed.stn);
                fprintf(f, "%e\n", nm);
            }
            higcit_destroy(it);
            fprintf(f, "\nSCALARS \u03C1\u2091 FLOAT\nLOOKUP_TABLE default\n");
            for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                c = higcit_getcell(it);
                hig_get_center(c, ccenter);
                np   = compute_value_at_point(ns->ed.eo.sdEOnplus, ccenter, ccenter, 1.0, ns->ed.eo.dpnplus, ns->ed.stn);
                nm   = compute_value_at_point(ns->ed.eo.sdEOnminus, ccenter, ccenter, 1.0, ns->ed.eo.dpnminus, ns->ed.stn);
                rhoe = deltaeo*(np - nm);
                fprintf(f, "%e\n", rhoe);
            }
            higcit_destroy(it);
        break;
    }
    fclose(f);
}


static long get_offset_cummulative(long proc_block_size){
    long offset_cum;
    MPI_Scan(&proc_block_size, &offset_cum, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    offset_cum -= proc_block_size;
    return offset_cum;
}

static long get_offset_sum(long proc_block_size){
    long offset_sum;
    MPI_Allreduce(&proc_block_size, &offset_sum, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    return offset_sum;
}

static void paddn_before_last(char *str, int N) {
    int len = strlen(str);
    assert(len <= N);
    int numSpaces = N - len;
    if(numSpaces == 0) return;

    char curr_last = str[len - 1];
    assert(curr_last == '\n'); // check if last character is newline
    str[N - 1] = curr_last;
    // Insert spaces before the last character
    memset(str + len - 1, ' ', numSpaces); // len - 1 to exclude the last
    str[N] = '\0'; // Ensure null termination at the end
}

static void update_buffer_write(MPI_File *f, long offset, char *write_buff, int buff_count, char *local_str, int local_str_size, int it_count) {
    int irem = it_count % buff_count;
    strncpy(write_buff + irem * local_str_size, local_str, local_str_size);
    if(irem == buff_count - 1) {
        int local_offset = (it_count - irem)*local_str_size;
        MPI_File_write_at(*f, offset + local_offset, write_buff, buff_count * local_str_size, MPI_CHAR, MPI_STATUS_IGNORE);
    }
}

static void write_remainder(MPI_File *f, long offset, char *write_buff, int buff_size, long proc_block_size) {
    int offset_remainder = proc_block_size % buff_size;
    if(offset_remainder != 0) { // write the remaining values
        int local_offset = proc_block_size - offset_remainder;
        MPI_File_write_at(*f, offset + local_offset, write_buff, offset_remainder, MPI_CHAR, MPI_STATUS_IGNORE);
    }
}


void higflow_print_vtk2D_parallel_single(higflow_solver *ns, int rank, int nprocs) {
    //  Open the VTK file
    real Re   = ns->par.Re;
    real beta = ns->ed.ve.par.beta;
    char vtkname[1024];
    snprintf(vtkname, sizeof vtkname, "%s_%d.vtk", ns->par.nameprint, ns->par.frame);
    //FILE *f;
    sim_domain *sdp = psd_get_local_domain(ns->psdp);

    //////////////////// Get number of leafs /////////////////////////////////////////////
    higcit_celliterator *it;
    // it_t = higcit_create_all_leaves(root);
    it = sd_get_domain_celliterator(sdp);
    long numleafs = higcit_count_without_advancing(it);
    long numleafs_global;
    MPI_Allreduce(&numleafs, &numleafs_global, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    /////////////////// Buffers and variables ////////////////////////////////////////////

    long curr_file_ptr_pos = 0, proc_block_size, proc_offset;
    char local_str[1024];
    int local_str_size, lf_max_size = 20, d_max_size, e_max_size = 13;
    long it_count;
    long buff_count = min(numleafs, 1048576), buff_size;
    char *write_buff;

    //////////////////////////////////////////////////////////////////////////////////////

    MPI_File f;
    int openerr = MPI_File_open(MPI_COMM_WORLD, vtkname, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f);
    if(openerr != MPI_SUCCESS) {
        printf("Error opening file %s\n", vtkname);
        return;
    }
    
    /////////////////////////////// Header /////////////////////////////////////////////////
    
    // header
    sprintf(local_str, "# vtk DataFile Version 3.0\nhigtree\nASCII\nDATASET UNSTRUCTURED_GRID\n\nPOINTS %ld float\n",
    4*numleafs_global);
    if(rank == 0) MPI_File_write_at(f, curr_file_ptr_pos, local_str, strlen(local_str), MPI_CHAR, MPI_STATUS_IGNORE);
    curr_file_ptr_pos += strlen(local_str);

    ////////////////////////////// Points //////////////////////////////////////////////////

    // hardcoded buffer size
    local_str_size = 4 * (2 * e_max_size + 4) * sizeof(char);
    write_buff = (char *)malloc(buff_count * local_str_size);
    proc_block_size = local_str_size * numleafs;
    proc_offset = get_offset_cummulative(proc_block_size);
    it_count = 0;
    
    for (; !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        // update buffer
        sprintf(local_str, "%e %e 0\n%e %e 0\n%e %e 0\n%e %e 0\n", 
        c->lowpoint[0], c->lowpoint[1], c->highpoint[0], c->lowpoint[1], 
        c->highpoint[0], c->highpoint[1], c->lowpoint[0], c->highpoint[1]);
        paddn_before_last(local_str, local_str_size);
        update_buffer_write(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count, local_str, local_str_size, it_count);
        it_count++;
    }
    write_remainder(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count*local_str_size, proc_block_size);
    higcit_destroy(it);

    curr_file_ptr_pos += get_offset_sum(proc_block_size);
    free(write_buff);
    
    ///////////////////////////////// Cells /////////////////////////////////////////////////

    // header
    sprintf(local_str, "\nCELLS %ld %ld\n", numleafs_global, 5*numleafs_global);
    if(rank == 0) MPI_File_write_at(f, curr_file_ptr_pos, local_str, strlen(local_str), MPI_CHAR, MPI_STATUS_IGNORE);
    curr_file_ptr_pos += strlen(local_str);

    // hardcoded buffer size
    int num_leafs_cummulative;
    MPI_Scan(&numleafs, &num_leafs_cummulative, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    num_leafs_cummulative -= numleafs;
    char dummy_str[20];
    sprintf(dummy_str, "%ld", 4*(num_leafs_cummulative + numleafs));
    d_max_size = strlen(dummy_str);

    local_str_size = (d_max_size * 4 + 6) * sizeof(char);
    write_buff = (char *)malloc(buff_count * local_str_size);
    proc_block_size = local_str_size * numleafs;
    proc_offset = get_offset_cummulative(proc_block_size);
    it_count = 0;
    
    for (long i = num_leafs_cummulative; i < num_leafs_cummulative + numleafs; i++) {
        sprintf(local_str, "%d %ld %ld %ld %ld\n", 4, 4*i, 4*i+1, 4*i+2, 4*i+3);
        paddn_before_last(local_str, local_str_size);
        update_buffer_write(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count, local_str, local_str_size, it_count);
        it_count++;
    }
    write_remainder(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count*local_str_size, proc_block_size);
    curr_file_ptr_pos += get_offset_sum(proc_block_size);
    free(write_buff);

    ///////////////////////////////// Cell types /////////////////////////////////////////////////

    // header
    sprintf(local_str, "\nCELL_TYPES %ld\n", numleafs_global);
    if(rank == 0) MPI_File_write_at(f, curr_file_ptr_pos, local_str, strlen(local_str), MPI_CHAR, MPI_STATUS_IGNORE);
    curr_file_ptr_pos += strlen(local_str);

    // hardcoded buffer size
    local_str_size = 2 * sizeof(char);
    write_buff = (char *)malloc(buff_count * local_str_size);
    proc_block_size = local_str_size * numleafs;
    proc_offset = get_offset_cummulative(proc_block_size);
    it_count = 0;

    for (int i = num_leafs_cummulative; i < num_leafs_cummulative + numleafs; i++) {
        sprintf(local_str, "7 "); //vtk cell types 7 = VTK_POLYGON // 9
        update_buffer_write(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count, local_str, local_str_size, it_count);
        it_count++;
    }
    write_remainder(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count*local_str_size, proc_block_size);
    curr_file_ptr_pos += get_offset_sum(proc_block_size);
    free(write_buff);
    
    ///////////////////////////////// Vectors /////////////////////////////////////////////////

    // header
    sprintf(local_str, "\n\nPOINT_DATA %ld\nVECTORS vel FLOAT\n", 4*numleafs_global);
    if(rank == 0) MPI_File_write_at(f, curr_file_ptr_pos, local_str, strlen(local_str), MPI_CHAR, MPI_STATUS_IGNORE);
    curr_file_ptr_pos += strlen(local_str);

    // hardcoded buffer size
    local_str_size = 4 * (2 * e_max_size + 4) * sizeof(char);
    write_buff = (char *)malloc(buff_count * local_str_size);
    proc_block_size = local_str_size * numleafs;
    proc_offset = get_offset_cummulative(proc_block_size);
    it_count = 0;

    sim_facet_domain *sfdu2[DIM];
    sfdu2[0] = psfd_get_local_domain(ns->psfdu[0]); sfdu2[1] = psfd_get_local_domain(ns->psfdu[1]);
    distributed_property *dpu[DIM];
    dpu[0] = ns->dpu[0]; dpu[1] = ns->dpu[1];

    // Saving vector properties of cell faces
    for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        Point ccenter;
        hig_get_center(c,ccenter);
        // Pontos onde será interpolada a velocidade
        Point p0, p1, p2, p3;
        p0[0] = c->lowpoint[0];  p0[1] = c->lowpoint[1];
        p1[0] = c->highpoint[0]; p1[1] = c->lowpoint[1];
        p2[0] = c->highpoint[0]; p2[1] = c->highpoint[1];
        p3[0] = c->lowpoint[0];  p3[1] = c->highpoint[1];
        Point vel0, vel1, vel2, vel3;

        for(int dim = 0; dim < DIM; dim++) {
            vel0[dim]= compute_facet_value_at_point(sfdu2[dim], ccenter, p0, 1.0, dpu[dim], ns->stn);
            vel1[dim] = compute_facet_value_at_point(sfdu2[dim], ccenter, p1, 1.0, dpu[dim], ns->stn);
            vel2[dim] = compute_facet_value_at_point(sfdu2[dim], ccenter, p2, 1.0, dpu[dim], ns->stn);
            vel3[dim] = compute_facet_value_at_point(sfdu2[dim], ccenter, p3, 1.0, dpu[dim], ns->stn);
        }
        
        sprintf(local_str, "%e %e 0\n%e %e 0\n%e %e 0\n%e %e 0\n", 
        vel0[0], vel0[1], vel1[0], vel1[1], vel2[0], vel2[1], vel3[0], vel3[1]);
        paddn_before_last(local_str, local_str_size);
        update_buffer_write(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count, local_str, local_str_size, it_count);
        it_count++;
    }
    write_remainder(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count*local_str_size, proc_block_size);
    higcit_destroy(it);

    curr_file_ptr_pos += get_offset_sum(proc_block_size);

    
    switch (ns->contr.eoflow) {
        case false:
        break;

        case true:

            sprintf(local_str, "\nVECTORS F\u2091\u2092 FLOAT\n");
            if(rank == 0) MPI_File_write_at(f, curr_file_ptr_pos, local_str, strlen(local_str), MPI_CHAR, MPI_STATUS_IGNORE);
            curr_file_ptr_pos += strlen(local_str);

            // hardcoded buffer size
            local_str_size = 4 * (2 * e_max_size + 4) * sizeof(char);
            proc_block_size = local_str_size * numleafs;
            proc_offset = get_offset_cummulative(proc_block_size);
            it_count = 0;

            sim_facet_domain *sfdEOFeo2[DIM];
            sfdEOFeo2[0] = psfd_get_local_domain(ns->ed.eo.psfdEOFeo[0]); sfdEOFeo2[1] = psfd_get_local_domain(ns->ed.eo.psfdEOFeo[1]);
            distributed_property *dpEOFeo2[DIM];
            dpEOFeo2[0] = ns->ed.eo.dpFeo[0]; dpEOFeo2[1] = ns->ed.eo.dpFeo[1];
            Point p0, p1, p2, p3;
            
            for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                hig_cell *c = higcit_getcell(it);
                Point ccenter;
                hig_get_center(c,ccenter);
                // Pontos onde será interpolada a força eletrica
                p0[0] = c->lowpoint[0];  p0[1] = c->lowpoint[1];
                p1[0] = c->highpoint[0]; p1[1] = c->lowpoint[1];
                p2[0] = c->highpoint[0]; p2[1] = c->highpoint[1];
                p3[0] = c->lowpoint[0];  p3[1] = c->highpoint[1];
                Point Feo0, Feo1, Feo2, Feo3;

                for(int dim = 0; dim < DIM; dim++) {
                    Feo0[dim] = compute_facet_value_at_point(sfdEOFeo2[dim], ccenter, p0, 1.0, ns->ed.eo.dpFeo[dim], ns->ed.stn);
                    Feo1[dim] = compute_facet_value_at_point(sfdEOFeo2[dim], ccenter, p1, 1.0, ns->ed.eo.dpFeo[dim], ns->ed.stn);
                    Feo2[dim] = compute_facet_value_at_point(sfdEOFeo2[dim], ccenter, p2, 1.0, ns->ed.eo.dpFeo[dim], ns->ed.stn);
                    Feo3[dim] = compute_facet_value_at_point(sfdEOFeo2[dim], ccenter, p3, 1.0, ns->ed.eo.dpFeo[dim], ns->ed.stn);
                }

                sprintf(local_str, "%e %e 0\n%e %e 0\n%e %e 0\n%e %e 0\n",
                Feo0[0], Feo0[1], Feo1[0], Feo1[1], Feo2[0], Feo2[1], Feo3[0], Feo3[1]);
                paddn_before_last(local_str, local_str_size);
                update_buffer_write(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count, local_str, local_str_size, it_count);
                it_count++;
            }
            write_remainder(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count*local_str_size, proc_block_size);
            higcit_destroy(it);

            curr_file_ptr_pos += get_offset_sum(proc_block_size);
        break;
    }

    free(write_buff);
    
    ///////////////////////////////// Tensors /////////////////////////////////////////////////

    // Saving tensorial properties 
    switch (ns->contr.flowtype) {
        case NEWTONIAN:
        break;
        case GENERALIZED_NEWTONIAN:
        break;
        case MULTIPHASE:
            if(ns->ed.mult.contr.flowtype_either == VISCOELASTIC) {
                sprintf(local_str, "\nTENSORS \u03C4\u209A FLOAT\n");
                if(rank == 0) MPI_File_write_at(f, curr_file_ptr_pos, local_str, strlen(local_str), MPI_CHAR, MPI_STATUS_IGNORE);
                curr_file_ptr_pos += strlen(local_str);

                // hardcoded buffer size
                local_str_size = 4 * (4 * e_max_size + 14) * sizeof(char);
                write_buff = (char *)malloc(buff_count * local_str_size);
                proc_block_size = local_str_size * numleafs;
                proc_offset = get_offset_cummulative(proc_block_size);
                it_count = 0;

                
                real beta0 = ns->ed.mult.ve.par0.beta;
                real beta1 = ns->ed.mult.ve.par1.beta;
                real t = ns->par.t;
                real beta_interp;
                real fracvol;
                real visc, visc0, visc1;

                for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                    hig_cell *c = higcit_getcell(it);
                    Point ccenter;
                    hig_get_center(c,ccenter);

                    // Pontos onde será interpolada a velocidade
                    Point p0, p1, p2, p3;
                    p0[0] = c->lowpoint[0];  p0[1] = c->lowpoint[1];
                    p1[0] = c->highpoint[0]; p1[1] = c->lowpoint[1];
                    p2[0] = c->highpoint[0]; p2[1] = c->highpoint[1];
                    p3[0] = c->lowpoint[0];  p3[1] = c->highpoint[1];
                
                    real taup0[DIM][DIM], taup1[DIM][DIM], taup2[DIM][DIM], taup3[DIM][DIM];
                    real Dp0[DIM][DIM], Dp1[DIM][DIM], Dp2[DIM][DIM], Dp3[DIM][DIM];
                    for (int i = 0; i < DIM; i++) {
                        for (int j = 0; j < DIM; j++) {
                            // Get Du
                            Dp0[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p0, 1.0, ns->ed.mult.ve.dpD[i][j], ns->ed.stn);
                            Dp1[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p1, 1.0, ns->ed.mult.ve.dpD[i][j], ns->ed.stn);
                            Dp2[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p2, 1.0, ns->ed.mult.ve.dpD[i][j], ns->ed.stn);
                            Dp3[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p3, 1.0, ns->ed.mult.ve.dpD[i][j], ns->ed.stn);
                        }
                    }
                        
                    for (int i = 0; i < DIM; i++) {
                        for (int j = 0; j < DIM; j++) {
                            
                            fracvol = compute_value_at_point(sdp, ccenter, p0, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
                            visc0 = ns->ed.mult.get_density0(p0, t);
                            visc1 = ns->ed.mult.get_density1(p0, t);
                            visc = (1.0 - fracvol) * visc0 + fracvol * visc1;
                            beta_interp =  ((1 - fracvol) * visc0 * beta0 + fracvol * visc1 * beta1) / visc;
                            taup0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.mult.ve.dpS[i][j], ns->ed.stn);
                            taup0[i][j] += 2.0*(1-beta_interp)*0.5*(Dp0[i][j]+Dp0[j][i])/Re;

                            fracvol = compute_value_at_point(sdp, ccenter, p1, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
                            visc0 = ns->ed.mult.get_density0(p1, t);
                            visc1 = ns->ed.mult.get_density1(p1, t);
                            visc = (1.0 - fracvol) * visc0 + fracvol * visc1;
                            beta_interp =  ((1 - fracvol) * visc0 * beta0 + fracvol * visc1 * beta1) / visc;
                            taup1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.mult.ve.dpS[i][j], ns->ed.stn);
                            taup1[i][j] += 2.0*(1-beta_interp)*0.5*(Dp1[i][j]+Dp1[j][i])/Re;

                            fracvol = compute_value_at_point(sdp, ccenter, p2, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
                            visc0 = ns->ed.mult.get_density0(p2, t);
                            visc1 = ns->ed.mult.get_density1(p2, t);
                            visc = (1.0 - fracvol) * visc0 + fracvol * visc1;
                            beta_interp =  ((1 - fracvol) * visc0 * beta0 + fracvol * visc1 * beta1) / visc;
                            taup2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.mult.ve.dpS[i][j], ns->ed.stn);
                            taup2[i][j] += 2.0*(1-beta_interp)*0.5*(Dp2[i][j]+Dp2[j][i])/Re;

                            fracvol = compute_value_at_point(sdp, ccenter, p3, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
                            visc0 = ns->ed.mult.get_density0(p3, t);
                            visc1 = ns->ed.mult.get_density1(p3, t);
                            visc = (1.0 - fracvol) * visc0 + fracvol * visc1;
                            beta_interp =  ((1 - fracvol) * visc0 * beta0 + fracvol * visc1 * beta1) / visc;
                            taup3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.mult.ve.dpS[i][j], ns->ed.stn);
                            taup3[i][j] += 2.0*(1-beta_interp)*0.5*(Dp3[i][j]+Dp3[j][i])/Re;
                        } 
                    }

                    sprintf(local_str, "%e %e 0\n%e %e 0\n0 0 0\n%e %e 0\n%e %e 0\n0 0 0\n%e %e 0\n%e %e 0\n0 0 0\n%e %e 0\n%e %e 0\n0 0 0\n",
                    taup0[0][0], taup0[0][1], taup0[1][0], taup0[1][1],
                    taup1[0][0], taup1[0][1], taup1[1][0], taup1[1][1],
                    taup2[0][0], taup2[0][1], taup2[1][0], taup2[1][1],
                    taup3[0][0], taup3[0][1], taup3[1][0], taup3[1][1]);
                    paddn_before_last(local_str, local_str_size);
                    update_buffer_write(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count, local_str, local_str_size, it_count);
                    it_count++;
                }
                write_remainder(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count*local_str_size, proc_block_size);
                higcit_destroy(it);

                curr_file_ptr_pos += get_offset_sum(proc_block_size);
                free(write_buff);
            }

        break;
        case VISCOELASTIC:
            sprintf(local_str, "\nTENSORS \u03C4\u209A FLOAT\n");
            if(rank == 0) MPI_File_write_at(f, curr_file_ptr_pos, local_str, strlen(local_str), MPI_CHAR, MPI_STATUS_IGNORE);
            curr_file_ptr_pos += strlen(local_str);

            // hardcoded buffer size
            local_str_size = 4 * (4 * e_max_size + 14) * sizeof(char);
            write_buff = (char *)malloc(buff_count * local_str_size);
            proc_block_size = local_str_size * numleafs;
            proc_offset = get_offset_cummulative(proc_block_size);
            it_count = 0;
      
            for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                hig_cell *c = higcit_getcell(it);
                Point ccenter;
                hig_get_center(c,ccenter);
                // Pontos onde será interpolada a velocidade
                Point p0, p1, p2, p3;
                p0[0] = c->lowpoint[0];  p0[1] = c->lowpoint[1];
                p1[0] = c->highpoint[0]; p1[1] = c->lowpoint[1];
                p2[0] = c->highpoint[0]; p2[1] = c->highpoint[1];
                p3[0] = c->lowpoint[0];  p3[1] = c->highpoint[1];
            
                real taup0[DIM][DIM], taup1[DIM][DIM], taup2[DIM][DIM], taup3[DIM][DIM];
                real Dp0[DIM][DIM], Dp1[DIM][DIM], Dp2[DIM][DIM], Dp3[DIM][DIM];
                for (int i = 0; i < DIM; i++) {
                    for (int j = 0; j < DIM; j++) {
                        // Get Du
                        Dp0[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p0, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                        Dp1[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p1, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                        Dp2[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p2, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                        Dp3[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p3, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                    }
                }
                    
                for (int i = 0; i < DIM; i++) {
                    for (int j = 0; j < DIM; j++) {
                        taup0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                        taup0[i][j]+= 2.0*(1-beta)*0.5*(Dp0[i][j]+Dp0[j][i])/Re;
                        taup1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                        taup1[i][j]+= 2.0*(1-beta)*0.5*(Dp1[i][j]+Dp1[j][i])/Re;
                        taup2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                        taup2[i][j]+= 2.0*(1-beta)*0.5*(Dp2[i][j]+Dp2[j][i])/Re;
                        taup3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                        taup3[i][j]+= 2.0*(1-beta)*0.5*(Dp3[i][j]+Dp3[j][i])/Re;
                        // taup0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.ve.dpKernel[i][j], ns->ed.stn);
                        // taup1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.ve.dpKernel[i][j], ns->ed.stn);
                        // taup2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.ve.dpKernel[i][j], ns->ed.stn);
                        // taup3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.ve.dpKernel[i][j], ns->ed.stn);
                    } 
                }  

                sprintf(local_str, "%e %e 0\n%e %e 0\n0 0 0\n%e %e 0\n%e %e 0\n0 0 0\n%e %e 0\n%e %e 0\n0 0 0\n%e %e 0\n%e %e 0\n0 0 0\n",
                taup0[0][0], taup0[0][1], taup0[1][0], taup0[1][1],
                taup1[0][0], taup1[0][1], taup1[1][0], taup1[1][1],
                taup2[0][0], taup2[0][1], taup2[1][0], taup2[1][1],
                taup3[0][0], taup3[0][1], taup3[1][0], taup3[1][1]);
                paddn_before_last(local_str, local_str_size);
                update_buffer_write(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count, local_str, local_str_size, it_count);
                it_count++;
            }
            write_remainder(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count*local_str_size, proc_block_size);
            higcit_destroy(it);

            curr_file_ptr_pos += get_offset_sum(proc_block_size);
            free(write_buff);
                 
            break;
        
        case VISCOELASTIC_INTEGRAL:
            // Integral
            sprintf(local_str, "\nTENSORS \u03C4\u209A FLOAT\n");
            if(rank == 0) MPI_File_write_at(f, curr_file_ptr_pos, local_str, strlen(local_str), MPI_CHAR, MPI_STATUS_IGNORE);
            curr_file_ptr_pos += strlen(local_str);
            
            // hardcoded buffer size
            local_str_size = 4 * (4 * e_max_size + 14) * sizeof(char);
            write_buff = (char *)malloc(buff_count * local_str_size);
            proc_block_size = local_str_size * numleafs;
            proc_offset = get_offset_cummulative(proc_block_size);
            it_count = 0;

            //real Re   = ns->par.Re;
            for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                hig_cell *c = higcit_getcell(it);
                Point ccenter;
                hig_get_center(c,ccenter);
                // Pontos onde será interpolada a velocidade
                Point p0, p1, p2, p3;
                p0[0] = c->lowpoint[0];  p0[1] = c->lowpoint[1];
                p1[0] = c->highpoint[0]; p1[1] = c->lowpoint[1];
                p2[0] = c->highpoint[0]; p2[1] = c->highpoint[1];
                p3[0] = c->lowpoint[0];  p3[1] = c->highpoint[1];
                
                real taup0[DIM][DIM], taup1[DIM][DIM], taup2[DIM][DIM], taup3[DIM][DIM];
                real Dp0[DIM][DIM], Dp1[DIM][DIM], Dp2[DIM][DIM], Dp3[DIM][DIM];
                for (int i = 0; i < DIM; i++) {
                    for (int j = 0; j < DIM; j++) {
                        // Get Du
                        Dp0[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p0, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
                        Dp1[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p1, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
                        Dp2[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p2, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
                        Dp3[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p3, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
                    }
                }
                    
                for (int i = 0; i < DIM; i++) {
                    for (int j = 0; j < DIM; j++) {
                        taup0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
                        taup0[i][j]+= 2.0*(1-beta)*0.5*(Dp0[i][j]+Dp0[j][i])/Re;
                        taup1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
                        taup1[i][j]+= 2.0*(1-beta)*0.5*(Dp1[i][j]+Dp1[j][i])/Re;
                        taup2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
                        taup2[i][j]+= 2.0*(1-beta)*0.5*(Dp2[i][j]+Dp2[j][i])/Re;
                        taup3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
                        taup3[i][j]+= 2.0*(1-beta)*0.5*(Dp3[i][j]+Dp3[j][i])/Re;
                    } 
                }
                
                sprintf(local_str, "%e %e 0\n%e %e 0\n0 0 0\n%e %e 0\n%e %e 0\n0 0 0\n%e %e 0\n%e %e 0\n0 0 0\n%e %e 0\n%e %e 0\n0 0 0\n",
                taup0[0][0], taup0[0][1], taup0[1][0], taup0[1][1],
                taup1[0][0], taup1[0][1], taup1[1][0], taup1[1][1],
                taup2[0][0], taup2[0][1], taup2[1][0], taup2[1][1],
                taup3[0][0], taup3[0][1], taup3[1][0], taup3[1][1]);
                paddn_before_last(local_str, local_str_size);
                update_buffer_write(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count, local_str, local_str_size, it_count);
                it_count++;
            }
            write_remainder(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count*local_str_size, proc_block_size);
            higcit_destroy(it);

            curr_file_ptr_pos += get_offset_sum(proc_block_size);
            free(write_buff);

            break;
    }


    /////////////////////////// Scalars ///////////////////////////

    sprintf(local_str, "\nCELL_DATA %ld\nSCALARS p FLOAT\nLOOKUP_TABLE default\n", numleafs_global);
    if(rank == 0) MPI_File_write_at(f, curr_file_ptr_pos, local_str, strlen(local_str), MPI_CHAR, MPI_STATUS_IGNORE);
    curr_file_ptr_pos += strlen(local_str);

    // hardcoded buffer size
    local_str_size = (e_max_size + 1) * sizeof(char);
    write_buff = (char *)malloc(buff_count * local_str_size);
    proc_block_size = local_str_size * numleafs;
    proc_offset = get_offset_cummulative(proc_block_size);
    it_count = 0;

    Point ccenter;

    mp_mapper *mp = sd_get_domain_mapper(sdp);

    // Saving pressure from the center of the cell 
    for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        hig_get_center(c,ccenter);
        real val  = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->dpp, ns->stn);
        // int clid = mp_lookup(mp, hig_get_cid(c));
        // real val = dp_get_value(ns->dpp, clid);
        
        sprintf(local_str, "%e\n", val);
        paddn_before_last(local_str, local_str_size);
        update_buffer_write(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count, local_str, local_str_size, it_count);
        it_count++;
    }
    write_remainder(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count*local_str_size, proc_block_size);
    higcit_destroy(it);

    curr_file_ptr_pos += get_offset_sum(proc_block_size);
    
    // Saving scalar properties from the center of the cell 
    switch (ns->contr.flowtype) {
        case NEWTONIAN:
        break;
        
        case GENERALIZED_NEWTONIAN:
        break;
        
        case MULTIPHASE:

            sprintf(local_str, "\nSCALARS FracVol FLOAT\nLOOKUP_TABLE default\n");
            if(rank == 0) MPI_File_write_at(f, curr_file_ptr_pos, local_str, strlen(local_str), MPI_CHAR, MPI_STATUS_IGNORE);
            curr_file_ptr_pos += strlen(local_str);

            // hardcoded buffer size
            local_str_size = (e_max_size + 1) * sizeof(char);
            proc_block_size = local_str_size * numleafs;
            proc_offset = get_offset_cummulative(proc_block_size);
            it_count = 0;

            for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                hig_cell *c = higcit_getcell(it);
                hig_get_center(c, ccenter);
                real value  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
                
                sprintf(local_str, "%e\n", value);
                paddn_before_last(local_str, local_str_size);
                update_buffer_write(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count, local_str, local_str_size, it_count);
                it_count++;
            }
            write_remainder(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count*local_str_size, proc_block_size);
            higcit_destroy(it);

            curr_file_ptr_pos += get_offset_sum(proc_block_size);
                    
        break;
            
        case VISCOELASTIC:
        break;
        
        case VISCOELASTIC_INTEGRAL:
        break;
    }

    switch (ns->contr.eoflow) {
        case false:
        break;

        case true: ;
            // Get the constant
            real deltaeo = ns->ed.eo.par.delta;
            // Necessary properties
            real phi, psi, np, nm, rhoe;
            // loop variables
            Point ccenter;
            hig_cell *c;

            ////////////// phi //////////////

            sprintf(local_str, "\nSCALARS \u03D5 FLOAT\nLOOKUP_TABLE default\n");
            if(rank == 0) MPI_File_write_at(f, curr_file_ptr_pos, local_str, strlen(local_str), MPI_CHAR, MPI_STATUS_IGNORE);
            curr_file_ptr_pos += strlen(local_str);
            
            // hardcoded buffer size
            local_str_size = (e_max_size + 1) * sizeof(char);
            proc_block_size = local_str_size * numleafs;
            proc_offset = get_offset_cummulative(proc_block_size);
            it_count = 0;
            
            for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                c = higcit_getcell(it);
                hig_get_center(c, ccenter);
                phi  = compute_value_at_point(ns->ed.eo.sdEOphi, ccenter, ccenter, 1.0, ns->ed.eo.dpphi, ns->ed.stn);
                
                sprintf(local_str, "%e\n", phi);
                paddn_before_last(local_str, local_str_size);
                update_buffer_write(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count, local_str, local_str_size, it_count);
                it_count++;
            }
            write_remainder(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count*local_str_size, proc_block_size);
            higcit_destroy(it);

            curr_file_ptr_pos += get_offset_sum(proc_block_size);

            ////////////// psi //////////////

            sprintf(local_str, "\nSCALARS \u03C8 FLOAT\nLOOKUP_TABLE default\n");
            if(rank == 0) MPI_File_write_at(f, curr_file_ptr_pos, local_str, strlen(local_str), MPI_CHAR, MPI_STATUS_IGNORE);
            curr_file_ptr_pos += strlen(local_str);
            
            // hardcoded buffer size
            local_str_size = (e_max_size + 1) * sizeof(char);
            proc_block_size = local_str_size * numleafs;
            proc_offset = get_offset_cummulative(proc_block_size);
            it_count = 0;
                    
            for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                c = higcit_getcell(it);
                hig_get_center(c, ccenter);
                psi  = compute_value_at_point(ns->ed.eo.sdEOpsi, ccenter, ccenter, 1.0, ns->ed.eo.dppsi, ns->ed.stn);
                
                sprintf(local_str, "%e\n", psi);
                paddn_before_last(local_str, local_str_size);
                update_buffer_write(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count, local_str, local_str_size, it_count);
                it_count++;
            }
            write_remainder(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count*local_str_size, proc_block_size);
            higcit_destroy(it);

            curr_file_ptr_pos += get_offset_sum(proc_block_size);

            /////////////// nplus ///////////////

            sprintf(local_str, "\nSCALARS n\u207A FLOAT\nLOOKUP_TABLE default\n");
            if(rank == 0) MPI_File_write_at(f, curr_file_ptr_pos, local_str, strlen(local_str), MPI_CHAR, MPI_STATUS_IGNORE);
            curr_file_ptr_pos += strlen(local_str);
            
            // hardcoded buffer size
            local_str_size = (e_max_size + 1) * sizeof(char);
            proc_block_size = local_str_size * numleafs;
            proc_offset = get_offset_cummulative(proc_block_size);
            it_count = 0;

            for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                c = higcit_getcell(it);
                hig_get_center(c, ccenter);
                np   = compute_value_at_point(ns->ed.eo.sdEOnplus, ccenter, ccenter, 1.0, ns->ed.eo.dpnplus, ns->ed.stn);
                
                sprintf(local_str, "%e\n", np);
                paddn_before_last(local_str, local_str_size);
                update_buffer_write(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count, local_str, local_str_size, it_count);
                it_count++;
            }
            write_remainder(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count*local_str_size, proc_block_size);
            higcit_destroy(it);

            curr_file_ptr_pos += get_offset_sum(proc_block_size);

            /////////////// nminus ///////////////

            sprintf(local_str, "\nSCALARS n\u207B FLOAT\nLOOKUP_TABLE default\n");
            if(rank == 0) MPI_File_write_at(f, curr_file_ptr_pos, local_str, strlen(local_str), MPI_CHAR, MPI_STATUS_IGNORE);
            curr_file_ptr_pos += strlen(local_str);
            
            // hardcoded buffer size
            local_str_size = (e_max_size + 1) * sizeof(char);
            proc_block_size = local_str_size * numleafs;
            proc_offset = get_offset_cummulative(proc_block_size);
            it_count = 0;
                    
            for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                c = higcit_getcell(it);
                hig_get_center(c, ccenter);
                nm   = compute_value_at_point(ns->ed.eo.sdEOnminus, ccenter, ccenter, 1.0, ns->ed.eo.dpnminus, ns->ed.stn);
                
                sprintf(local_str, "%e\n", nm);
                paddn_before_last(local_str, local_str_size);
                update_buffer_write(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count, local_str, local_str_size, it_count);
                it_count++;
            }
            write_remainder(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count*local_str_size, proc_block_size);
            higcit_destroy(it);

            curr_file_ptr_pos += get_offset_sum(proc_block_size);

            /////////////// rhoe ///////////////

            sprintf(local_str, "\nSCALARS \u03C1\u2091 FLOAT\nLOOKUP_TABLE default\n");
            if(rank == 0) MPI_File_write_at(f, curr_file_ptr_pos, local_str, strlen(local_str), MPI_CHAR, MPI_STATUS_IGNORE);
            curr_file_ptr_pos += strlen(local_str);
            
            // hardcoded buffer size
            local_str_size = (e_max_size + 1) * sizeof(char);
            proc_block_size = local_str_size * numleafs;
            proc_offset = get_offset_cummulative(proc_block_size);
            it_count = 0;
                    
            for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                c = higcit_getcell(it);
                hig_get_center(c, ccenter);
                np   = compute_value_at_point(ns->ed.eo.sdEOnplus, ccenter, ccenter, 1.0, ns->ed.eo.dpnplus, ns->ed.stn);
                nm   = compute_value_at_point(ns->ed.eo.sdEOnminus, ccenter, ccenter, 1.0, ns->ed.eo.dpnminus, ns->ed.stn);
                rhoe = deltaeo*(np - nm);
                
                sprintf(local_str, "%e\n", rhoe);
                paddn_before_last(local_str, local_str_size);
                update_buffer_write(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count, local_str, local_str_size, it_count);
                it_count++;
            }
            write_remainder(&f, curr_file_ptr_pos + proc_offset, write_buff, buff_count*local_str_size, proc_block_size);
            higcit_destroy(it);

            curr_file_ptr_pos += get_offset_sum(proc_block_size);
      
        break;
    }

    free(write_buff);
    MPI_File_close(&f);
}


// Print the VTK file for visualize 3D
void higflow_print_vtk3D(higflow_solver *ns, int rank) {
    // Open the VTK file
    real Re   = ns->par.Re;
    real beta = ns->ed.ve.par.beta;
    char vtkname[1024];
    snprintf(vtkname, sizeof vtkname, "%s_%d-%d.vtk", ns->par.nameprint, rank, ns->par.frame);
    FILE *f = fopen(vtkname, "w");
    if (f == NULL) {
        return;
    }
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    int numhigs = sd_get_num_higtrees(sdp);
    higcit_celliterator *it_t;
    // it_t = higcit_create_all_leaves(root);
    it_t = sd_get_domain_celliterator(sdp);
    long numleafs = higcit_count_without_advancing(it_t);
    fprintf(f, "# vtk DataFile Version 3.0\n");
    fprintf(f, "higtree\n");
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(f, "POINTS %ld float\n", 8 * numleafs);
    hig_cell *c;
    for (; !higcit_isfinished(it_t); higcit_nextcell(it_t)) {
        c = higcit_getcell(it_t);
        fprintf(f, "%lf %lf %lf\n", c->lowpoint[0], c->lowpoint[1], c->lowpoint[2]);
        fprintf(f, "%lf %lf %lf\n", c->highpoint[0], c->lowpoint[1], c->lowpoint[2]);
        fprintf(f, "%lf %lf %lf\n", c->highpoint[0], c->highpoint[1], c->lowpoint[2]);
        fprintf(f, "%lf %lf %lf\n", c->lowpoint[0], c->highpoint[1], c->lowpoint[2]);
        fprintf(f, "%lf %lf %lf\n", c->lowpoint[0], c->lowpoint[1], c->highpoint[2]);
        fprintf(f, "%lf %lf %lf\n", c->highpoint[0], c->lowpoint[1], c->highpoint[2]);
        fprintf(f, "%lf %lf %lf\n", c->highpoint[0], c->highpoint[1], c->highpoint[2]);
        fprintf(f, "%lf %lf %lf\n", c->lowpoint[0], c->highpoint[1], c->highpoint[2]);
    }
    fprintf(f, "CELLS %ld %ld\n", numleafs, 9 * numleafs);
    for (int i = 0; i < numleafs; i++) {
        fprintf(f, "%d %d %d %d %d %d %d %d %d\n", 8, 8 * i + 0, 8 * i + 1, 8 * i + 2, 8 * i + 3, 8 * i + 4, 8 * i + 5, 8 * i + 6, 8 * i + 7); //x == 0
    }
    fprintf(f, "CELL_TYPES %ld\n", numleafs);
    for (int i = 0; i < numleafs; i++) {
        fprintf(f, "12 "); //x == 0
    }
    higcit_destroy(it_t);
    // higcit_celliterator *it = higcit_create_all_leaves(root);
    higcit_celliterator *it = sd_get_domain_celliterator(sdp);
    long numcells = higcit_count(it);
    //printf("\n ---- Number of cells: %lu\n", numcells);
    higcit_destroy(it);
    // sim_domain *sdp = psd_get_local_domain(psdp);
    mp_mapper *m = sd_get_domain_mapper(sdp);
    sim_facet_domain *sfdu2[DIM];
    int dimension[DIM];
    dimension[0] = 0; dimension[1] = 0; dimension[2] = 0;
    int degree = 2;
    const int maxpts = 2 * DIM * wls_num_min_points(DIM, degree);
    //const int maxpts = 2*DIM*wls_num_min_points(degree);
    Point pts[maxpts + 1];
    real dist[maxpts + 1];
    real w0[maxpts], w1[maxpts], w2[maxpts], w3[maxpts], w4[maxpts], w5[maxpts], w6[maxpts], w7[maxpts];
    uniqueid gids[maxpts];
    fprintf(f, "\nPOINT_DATA %ld\nVECTORS vel FLOAT\n", 8 * numleafs);
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        Point cdelta, ccenter, clowpoint, chightpoint;
        hig_get_delta(c, cdelta);
        hig_get_center(c, ccenter);
        // Pontos onde será interpolada a velocidade
        Point p0, p1, p2, p3, p4, p5, p6, p7;
        p0[0] = c->lowpoint[0]; p0[1] = c->lowpoint[1]; p0[2] = c->lowpoint[2];
        p1[0] = c->highpoint[0]; p1[1] = c->lowpoint[1]; p1[2] = c->lowpoint[2];
        p2[0] = c->highpoint[0]; p2[1] = c->highpoint[1]; p2[2] = c->lowpoint[2];
        p3[0] = c->lowpoint[0]; p3[1] = c->highpoint[1]; p3[2] = c->lowpoint[2];
        p4[0] = c->lowpoint[0]; p4[1] = c->lowpoint[1]; p4[2] = c->highpoint[2];
        p5[0] = c->highpoint[0]; p5[1] = c->lowpoint[1]; p5[2] = c->highpoint[2];
        p6[0] = c->highpoint[0]; p6[1] = c->highpoint[1]; p6[2] = c->highpoint[2];
        p7[0] = c->lowpoint[0]; p7[1] = c->highpoint[1]; p7[2] = c->highpoint[2];
        uniqueid id = hig_get_cid(c);
        int clid = mp_lookup(m, id);
        double lu0 = 0, lv0 = 0, lw0 = 0;
        double lu1 = 0, lv1 = 0, lw1 = 0;
        double lu2 = 0, lv2 = 0, lw2 = 0;
        double lu3 = 0, lv3 = 0, lw3 = 0;
        double lu4 = 0, lv4 = 0, lw4 = 0;
        double lu5 = 0, lv5 = 0, lw5 = 0;
        double lu6 = 0, lv6 = 0, lw6 = 0;
        double lu7 = 0, lv7 = 0, lw7 = 0;
        for (int dim = 0; dim < DIM; dim++) {
            // int numpts = 0;
            // int cont = 0;
            sfdu2[dim] = psfd_get_local_domain(ns->psfdu[dim]);
            mp_mapper *m = sfd_get_domain_mapper(sfdu2[dim]);
            real value0 = compute_facet_value_at_point(sfdu2[dim], ccenter, p0, 1.0, ns->dpu[dim], ns->stn);
            switch (dim) {
                case 0:
                    lu0 = value0;
                break;
                case 1:
                    lv0 = value0;
                break;
                case 2:
                    lw0 = value0;
                break;
            }
            real value1 = compute_facet_value_at_point(sfdu2[dim], ccenter, p1, 1.0, ns->dpu[dim], ns->stn);
            switch (dim) {
                case 0:
                    lu1 = value1;
                break;
                case 1:
                    lv1 = value1;
                break;
                case 2:
                    lw1 = value1;
                break;
            }
            real value2 = compute_facet_value_at_point(sfdu2[dim], ccenter, p2, 1.0, ns->dpu[dim], ns->stn);
            switch (dim) {
                case 0:
                    lu2 = value2;
                break;
                case 1:
                    lv2 = value2;
                break;
                case 2:
                    lw2 = value2;
                break;
            }
            real value3 = compute_facet_value_at_point(sfdu2[dim], ccenter, p3, 1.0, ns->dpu[dim], ns->stn);
            switch (dim) {
                case 0:
                    lu3 = value3;
                break;
                case 1:
                    lv3 = value3;
                break;
                case 2:
                    lw3 = value3;
                break;
            }
            real value4 = compute_facet_value_at_point(sfdu2[dim], ccenter, p4, 1.0, ns->dpu[dim], ns->stn);
            switch (dim) {
                case 0:
                    lu4 = value4;
                break;
                case 1:
                    lv4 = value4;
                break;
                case 2:
                    lw4 = value4;
                break;
            }
            real value5 = compute_facet_value_at_point(sfdu2[dim], ccenter, p5, 1.0, ns->dpu[dim], ns->stn);
            switch (dim) {
                case 0:
                    lu5 = value5;
                break;
                case 1:
                    lv5 = value5;
                break;
                case 2:
                    lw5 = value5;
                break;
            }
            real value6 = compute_facet_value_at_point(sfdu2[dim], ccenter, p6, 1.0, ns->dpu[dim], ns->stn);
            switch (dim) {
                case 0:
                    lu6 = value6;
                break;
                case 1:
                    lv6 = value6;
                break;
                case 2:
                    lw6 = value6;
                break;
            }
            real value7 = compute_facet_value_at_point(sfdu2[dim], ccenter, p7, 1.0, ns->dpu[dim], ns->stn);
            switch (dim) {
                case 0:
                    lu7 = value7;
                break;
                case 1:
                    lv7 = value7;
                break;
                case 2:
                    lw7 = value7;
                break;
            }
            // printf("Máximo de pontos --> %d\nMínimo de pontos --> %d\nPontos usados --> %d\n",maxpts,wls_num_min_points(DIM,degree),cont);
        }
        fprintf(f, "%e %e %e\n", lu0, lv0, lw0);
        fprintf(f, "%e %e %e\n", lu1, lv1, lw1);
        fprintf(f, "%e %e %e\n", lu2, lv2, lw2);
        fprintf(f, "%e %e %e\n", lu3, lv3, lw3);
        fprintf(f, "%e %e %e\n", lu4, lv4, lw4);
        fprintf(f, "%e %e %e\n", lu5, lv5, lw5);
        fprintf(f, "%e %e %e\n", lu6, lv6, lw6);
        fprintf(f, "%e %e %e\n", lu7, lv7, lw7);
    }
    higcit_destroy(it);
    // Saving vectorial properties 
    switch (ns->contr.flowtype) {
        case NEWTONIAN:
        break;
        
        case GENERALIZED_NEWTONIAN:
        break;
        
        case MULTIPHASE:
         break;
            
        case VISCOELASTIC:
            fprintf(f, "\nTENSORS Tensor FLOAT\n");
            for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                hig_cell *c = higcit_getcell(it);
                Point cdelta, ccenter, clowpoint, chightpoint;
                hig_get_delta(c,cdelta);
                hig_get_center(c,ccenter);
                // Pontos onde será interpolada a velocidade
                Point p0, p1, p2, p3, p4, p5, p6, p7;
                p0[0] = c->lowpoint[0];  p0[1] = c->lowpoint[1];  p0[2] = c->lowpoint[2];
                p1[0] = c->highpoint[0]; p1[1] = c->lowpoint[1];  p1[2] = c->lowpoint[2];
                p2[0] = c->highpoint[0]; p2[1] = c->highpoint[1]; p2[2] = c->lowpoint[2];
                p3[0] = c->lowpoint[0];  p3[1] = c->highpoint[1]; p3[2] = c->lowpoint[2];
                p4[0] = c->lowpoint[0];  p4[1] = c->lowpoint[1];  p4[2] = c->highpoint[2];
                p5[0] = c->highpoint[0]; p5[1] = c->lowpoint[1];  p5[2] = c->highpoint[2];
                p6[0] = c->highpoint[0]; p6[1] = c->highpoint[1]; p6[2] = c->highpoint[2];
                p7[0] = c->lowpoint[0];  p7[1] = c->highpoint[1]; p7[2] = c->highpoint[2];
                uniqueid id = hig_get_cid(c);
                int clid = mp_lookup(m, id);
                real taup0[DIM+1][DIM+1], taup1[DIM+1][DIM+1], taup2[DIM+1][DIM+1], taup3[DIM+1][DIM+1];
                real taup4[DIM+1][DIM+1], taup5[DIM+1][DIM+1], taup6[DIM+1][DIM+1], taup7[DIM+1][DIM+1];
                for (int i = 0; i <= DIM; i++) {
                    for (int j = 0; j <= DIM; j++) {
                        taup0[i][j]=0.0; taup1[i][j]=0.0; taup2[i][j]=0.0; taup3[i][j]=0.0;
                        taup4[i][j]=0.0; taup5[i][j]=0.0; taup6[i][j]=0.0; taup7[i][j]=0.0;
                    }
                }
                sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
                real Dp0[DIM][DIM], Dp1[DIM][DIM], Dp2[DIM][DIM], Dp3[DIM][DIM];
                real Dp4[DIM][DIM], Dp5[DIM][DIM], Dp6[DIM][DIM], Dp7[DIM][DIM];
                for (int i = 0; i < DIM; i++) {
                    for (int j = 0; j < DIM; j++) {
                        // Get Du
                        Dp0[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p0, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                        Dp1[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p1, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                        Dp2[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p2, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                        Dp3[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p3, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                        Dp4[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p4, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                        Dp5[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p5, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                        Dp6[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p6, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                        Dp7[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p7, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                    }
                }
                 for (int i = 0; i < DIM; i++) {
                   for (int j = 0; j < DIM; j++) {
                     taup0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                     taup0[i][j]+= 2.0*(1-beta)*0.5*(Dp0[i][j]+Dp0[j][i])/Re;
                     taup1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                     taup1[i][j]+= 2.0*(1-beta)*0.5*(Dp1[i][j]+Dp1[j][i])/Re;
                     taup2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                     taup2[i][j]+= 2.0*(1-beta)*0.5*(Dp2[i][j]+Dp2[j][i])/Re;
                     taup3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                     taup3[i][j]+= 2.0*(1-beta)*0.5*(Dp3[i][j]+Dp3[j][i])/Re;
                     taup4[i][j] = compute_value_at_point(ns->ed.sdED, p4, p4, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                     taup4[i][j]+= 2.0*(1-beta)*0.5*(Dp4[i][j]+Dp4[j][i])/Re;
                     taup5[i][j] = compute_value_at_point(ns->ed.sdED, p5, p5, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                     taup5[i][j]+= 2.0*(1-beta)*0.5*(Dp5[i][j]+Dp5[j][i])/Re;
                     taup6[i][j] = compute_value_at_point(ns->ed.sdED, p6, p6, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                     taup6[i][j]+= 2.0*(1-beta)*0.5*(Dp6[i][j]+Dp6[j][i])/Re;
                     taup7[i][j] = compute_value_at_point(ns->ed.sdED, p7, p7, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                     taup7[i][j]+= 2.0*(1-beta)*0.5*(Dp7[i][j]+Dp7[j][i])/Re;
                  }
                }
                for (int i = 0; i <= DIM; i++) {
                   for (int j = 0; j <= DIM; j++) {
                      fprintf(f, "%e ", taup0[i][j]);
                   }
                   fprintf(f, "\n");
                }
                for (int i = 0; i <= DIM; i++) {
                   for (int j = 0; j <= DIM; j++) {
                      fprintf(f, "%e ", taup1[i][j]);
                   }
                   fprintf(f, "\n");
                }
                for (int i = 0; i <= DIM; i++) {
                   for (int j = 0; j <= DIM; j++) {
                      fprintf(f, "%e ", taup2[i][j]);
                   }
                   fprintf(f, "\n");
                }
                for (int i = 0; i <= DIM; i++) {
                   for (int j = 0; j <= DIM; j++) {
                      fprintf(f, "%e ", taup3[i][j]);
                   }
                   fprintf(f, "\n");
                }
                for (int i = 0; i <= DIM; i++) {
                   for (int j = 0; j <= DIM; j++) {
                      fprintf(f, "%e ", taup4[i][j]);
                   }
                   fprintf(f, "\n");
                }
                for (int i = 0; i <= DIM; i++) {
                   for (int j = 0; j <= DIM; j++) {
                      fprintf(f, "%e ", taup5[i][j]);
                   }
                   fprintf(f, "\n");
                }
                for (int i = 0; i <= DIM; i++) {
                   for (int j = 0; j <= DIM; j++) {
                      fprintf(f, "%e ", taup6[i][j]);
                   }
                   fprintf(f, "\n");
                }
                for (int i = 0; i <= DIM; i++) {
                   for (int j = 0; j <= DIM; j++) {
                      fprintf(f, "%e ", taup7[i][j]);
                   }
                   fprintf(f, "\n");
                }
            }
            higcit_destroy(it);
        break;
        
        case VISCOELASTIC_INTEGRAL:
        // Integral
                   fprintf(f, "\nTENSORS Tensor FLOAT\n");
                   for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                        hig_cell *c = higcit_getcell(it);
                        Point cdelta, ccenter, clowpoint, chightpoint;
                        hig_get_delta(c,cdelta);
                        hig_get_center(c,ccenter);
                        // Pontos onde será interpolada a velocidade
                        Point p0, p1, p2, p3, p4, p5, p6, p7;
                        p0[0] = c->lowpoint[0];  p0[1] = c->lowpoint[1];  p0[2] = c->lowpoint[2];
                        p1[0] = c->highpoint[0]; p1[1] = c->lowpoint[1];  p1[2] = c->lowpoint[2];
                        p2[0] = c->highpoint[0]; p2[1] = c->highpoint[1]; p2[2] = c->lowpoint[2];
                        p3[0] = c->lowpoint[0];  p3[1] = c->highpoint[1]; p3[2] = c->lowpoint[2];
                        p4[0] = c->lowpoint[0];  p4[1] = c->lowpoint[1];  p4[2] = c->highpoint[2];
                        p5[0] = c->highpoint[0]; p5[1] = c->lowpoint[1];  p5[2] = c->highpoint[2];
                        p6[0] = c->highpoint[0]; p6[1] = c->highpoint[1]; p6[2] = c->highpoint[2];
                        p7[0] = c->lowpoint[0];  p7[1] = c->highpoint[1]; p7[2] = c->highpoint[2];

                        uniqueid id = hig_get_cid(c);
                        int clid = mp_lookup(m, id);
                        
                        real taup0[DIM+1][DIM+1], taup1[DIM+1][DIM+1], taup2[DIM+1][DIM+1], taup3[DIM+1][DIM+1];
                        real taup4[DIM+1][DIM+1], taup5[DIM+1][DIM+1], taup6[DIM+1][DIM+1], taup7[DIM+1][DIM+1];
                        
                        for (int i = 0; i <= DIM; i++) {
                          for (int j = 0; j <= DIM; j++) {
                            taup0[i][j]=0.0; taup1[i][j]=0.0; taup2[i][j]=0.0; taup3[i][j]=0.0;
                            taup4[i][j]=0.0; taup5[i][j]=0.0; taup6[i][j]=0.0; taup7[i][j]=0.0;
                          }
                        }
                        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
                        
                        real Dp0[DIM][DIM], Dp1[DIM][DIM], Dp2[DIM][DIM], Dp3[DIM][DIM];
                        real Dp4[DIM][DIM], Dp5[DIM][DIM], Dp6[DIM][DIM], Dp7[DIM][DIM];
                        for (int i = 0; i < DIM; i++) {
                          for (int j = 0; j < DIM; j++) {
                             // Get Du
                             Dp0[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p0, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
                             Dp1[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p1, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
                             Dp2[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p2, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
                             Dp3[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p3, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
                             Dp4[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p4, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
                             Dp5[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p5, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
                             Dp6[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p6, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
                             Dp7[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p7, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
                          }
                         }
                         
                         for (int i = 0; i < DIM; i++) {
                           for (int j = 0; j < DIM; j++) {
                             taup0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
                             taup0[i][j]+= (Dp0[i][j]+Dp0[j][i])/Re;
                             
                             taup1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
                             taup1[i][j]+= (Dp1[i][j]+Dp1[j][i])/Re;
                             
                             taup2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
                             taup2[i][j]+= (Dp2[i][j]+Dp2[j][i])/Re;
                             
                             taup3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
                             taup3[i][j]+= (Dp3[i][j]+Dp3[j][i])/Re;
                             
                             taup4[i][j] = compute_value_at_point(ns->ed.sdED, p4, p4, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
                             taup4[i][j]+= (Dp4[i][j]+Dp4[j][i])/Re;
                             
                             taup5[i][j] = compute_value_at_point(ns->ed.sdED, p5, p5, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
                             taup5[i][j]+= (Dp5[i][j]+Dp5[j][i])/Re;
                             
                             taup6[i][j] = compute_value_at_point(ns->ed.sdED, p6, p6, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
                             taup6[i][j]+= (Dp6[i][j]+Dp6[j][i])/Re;
                             
                             taup7[i][j] = compute_value_at_point(ns->ed.sdED, p7, p7, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
                             taup7[i][j]+= (Dp7[i][j]+Dp7[j][i])/Re;
                          }
                        }
                        
                        for (int i = 0; i <= DIM; i++) {
                           for (int j = 0; j <= DIM; j++) {
                              fprintf(f, "%e ", taup0[i][j]);
                           }
                           fprintf(f, "\n");
                        }
                        
                        for (int i = 0; i <= DIM; i++) {
                           for (int j = 0; j <= DIM; j++) {
                              fprintf(f, "%e ", taup1[i][j]);
                           }
                           fprintf(f, "\n");
                        }
                        
                        for (int i = 0; i <= DIM; i++) {
                           for (int j = 0; j <= DIM; j++) {
                              fprintf(f, "%e ", taup2[i][j]);
                           }
                           fprintf(f, "\n");
                        }
                        
                        for (int i = 0; i <= DIM; i++) {
                           for (int j = 0; j <= DIM; j++) {
                              fprintf(f, "%e ", taup3[i][j]);
                           }
                           fprintf(f, "\n");
                        }
                        
                        for (int i = 0; i <= DIM; i++) {
                           for (int j = 0; j <= DIM; j++) {
                              fprintf(f, "%e ", taup4[i][j]);
                           }
                           fprintf(f, "\n");
                        }
                        
                        for (int i = 0; i <= DIM; i++) {
                           for (int j = 0; j <= DIM; j++) {
                              fprintf(f, "%e ", taup5[i][j]);
                           }
                           fprintf(f, "\n");
                        }
                        
                        for (int i = 0; i <= DIM; i++) {
                           for (int j = 0; j <= DIM; j++) {
                              fprintf(f, "%e ", taup6[i][j]);
                           }
                           fprintf(f, "\n");
                        }
                        
                        for (int i = 0; i <= DIM; i++) {
                           for (int j = 0; j <= DIM; j++) {
                              fprintf(f, "%e ", taup7[i][j]);
                           }
                           fprintf(f, "\n");
                        }
                    }
                    higcit_destroy(it);
        break;
        
    }
    fprintf(f, "\nCELL_DATA %ld\nSCALARS p FLOAT\nLOOKUP_TABLE default\n", numcells);
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        uniqueid id = hig_get_cid(c);
        int clid = mp_lookup(m, id);
        real val = dp_get_value(ns->dpp, clid);
        fprintf(f, "%e\n", val);
        // fprintf(f, "%e\n", val);
    }
    higcit_destroy(it);
    
    // Saving scalar properties from the center of the cell 
    switch (ns->contr.flowtype) {
        case NEWTONIAN:
        break;
        
        case GENERALIZED_NEWTONIAN:
        break;
        
        case MULTIPHASE:
            fprintf(f,"\nSCALARS FracVol FLOAT\nLOOKUP_TABLE default\n");
            sfdu2[0] = psfd_get_local_domain(ns->psfdu[0]);
            for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
               hig_cell *c = higcit_getcell(it);
               uniqueid id = hig_get_cid(c);
               int clid = mp_lookup(m, id);
               Point cdelta, ccenter, p0, p1;
               hig_get_delta(c, cdelta);
               hig_get_center(c, ccenter);
               real value = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
               fprintf(f, "%e\n", value);
            }
            higcit_destroy(it);
         break;
            
        case VISCOELASTIC:
        break;
        
        case VISCOELASTIC_INTEGRAL:
        break;

    }
    
    fclose(f);
}

// Print the VTK file for vilusalize 3D
void higflow_print_vtk3D_viscoelastic(higflow_solver *ns, int rank) {
    // Open the VTK file
    char vtkname[1024];
    snprintf(vtkname, sizeof vtkname, "%s_%d-%d.vtk", ns->par.nameprint, rank, ns->par.frame);
    FILE *f = fopen(vtkname, "w");
    if (f == NULL) {
        return;
    }
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    int numhigs = sd_get_num_higtrees(sdp);
    higcit_celliterator *it_t;
    // it_t = higcit_create_all_leaves(root);
    it_t = sd_get_domain_celliterator(sdp);
    long numleafs = higcit_count_without_advancing(it_t);
    fprintf(f, "# vtk DataFile Version 3.0\n");
    fprintf(f, "higtree\n");
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(f, "POINTS %ld float\n", 8*numleafs);
    hig_cell *c;
    for (; !higcit_isfinished(it_t); higcit_nextcell(it_t)) {
        c = higcit_getcell(it_t);
        fprintf(f, "%lf %lf %lf\n", c->lowpoint[0], c->lowpoint[1], c->lowpoint[2]);
        fprintf(f, "%lf %lf %lf\n", c->highpoint[0], c->lowpoint[1], c->lowpoint[2]);
        fprintf(f, "%lf %lf %lf\n", c->highpoint[0], c->highpoint[1], c->lowpoint[2]);
        fprintf(f, "%lf %lf %lf\n", c->lowpoint[0], c->highpoint[1], c->lowpoint[2]);
        fprintf(f, "%lf %lf %lf\n", c->lowpoint[0], c->lowpoint[1], c->highpoint[2]);
        fprintf(f, "%lf %lf %lf\n", c->highpoint[0], c->lowpoint[1], c->highpoint[2]);
        fprintf(f, "%lf %lf %lf\n", c->highpoint[0], c->highpoint[1], c->highpoint[2]);
        fprintf(f, "%lf %lf %lf\n", c->lowpoint[0], c->highpoint[1], c->highpoint[2]);
    }
    fprintf(f, "CELLS %ld %ld\n", numleafs, 9*numleafs);
    for (int i = 0; i < numleafs; i++) {
        fprintf(f, "%d %d %d %d %d %d %d %d %d\n", 8, 8*i+0, 8*i+1, 8*i+2, 8*i+3, 8*i+4, 8*i+5, 8*i+6, 8*i+7); //x == 0
    }
    fprintf(f, "CELL_TYPES %ld\n", numleafs);
    for (int i = 0; i < numleafs; i++) {
        fprintf(f, "12 "); //x == 0
    }
    higcit_destroy(it_t);
    // higcit_celliterator *it = higcit_create_all_leaves(root);
    higcit_celliterator *it = sd_get_domain_celliterator(sdp);
    long numcells = higcit_count(it);
    //printf("\n ---- Number of cells: %lu\n", numcells);
    higcit_destroy(it);
    // sim_domain *sdp = psd_get_local_domain(psdp);
    mp_mapper *m = sd_get_domain_mapper(sdp);
    sim_facet_domain *sfdu2[DIM];
    int dimension[DIM];
    dimension[0] = 0;
    dimension[1] = 0;
    dimension[2] = 0;
    int degree = 2;
    const int maxpts = 2*DIM*wls_num_min_points(DIM,degree);
    //const int maxpts = 2*DIM*wls_num_min_points(degree);
    Point pts[maxpts+1];
    real dist[maxpts+1];
    real w0[maxpts], w1[maxpts], w2[maxpts], w3[maxpts], w4[maxpts], w5[maxpts], w6[maxpts], w7[maxpts];
    uniqueid gids[maxpts];
    fprintf(f, "\nPOINT_DATA %ld\nVECTORS vel FLOAT\n", 8*numleafs);
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        Point cdelta, ccenter, clowpoint, chightpoint;
        hig_get_delta(c, cdelta);
        hig_get_center(c, ccenter);
        // Pontos onde será interpolada a velocidade
        Point p0, p1, p2, p3, p4, p5, p6, p7;
        p0[0] = c->lowpoint[0];
        p0[1] = c->lowpoint[1];
        p0[2] = c->lowpoint[2];
        p1[0] = c->highpoint[0];
        p1[1] = c->lowpoint[1];
        p1[2] = c->lowpoint[2];
        p2[0] = c->highpoint[0];
        p2[1] = c->highpoint[1];
        p2[2] = c->lowpoint[2];
        p3[0] = c->lowpoint[0];
        p3[1] = c->highpoint[1];
        p3[2] = c->lowpoint[2];
        p4[0] = c->lowpoint[0];
        p4[1] = c->lowpoint[1];
        p4[2] = c->highpoint[2];
        p5[0] = c->highpoint[0];
        p5[1] = c->lowpoint[1];
        p5[2] = c->highpoint[2];
        p6[0] = c->highpoint[0];
        p6[1] = c->highpoint[1];
        p6[2] = c->highpoint[2];
        p7[0] = c->lowpoint[0];
        p7[1] = c->highpoint[1];
        p7[2] = c->highpoint[2];
        uniqueid id = hig_get_cid(c);
        int clid = mp_lookup(m, id);
        double lu0 = 0, lv0 = 0, lw0 = 0;
        double lu1 = 0, lv1 = 0, lw1 = 0;
        double lu2 = 0, lv2 = 0, lw2 = 0;
        double lu3 = 0, lv3 = 0, lw3 = 0;
        double lu4 = 0, lv4 = 0, lw4 = 0;
        double lu5 = 0, lv5 = 0, lw5 = 0;
        double lu6 = 0, lv6 = 0, lw6 = 0;
        double lu7 = 0, lv7 = 0, lw7 = 0;
        for (int dim = 0; dim < DIM; dim++) {
            // int numpts = 0;
            // int cont = 0;
            sfdu2[dim] = psfd_get_local_domain(ns->psfdu[dim]);
            mp_mapper *m = sfd_get_domain_mapper(sfdu2[dim]);
            real value0 = compute_facet_value_at_point(sfdu2[dim], ccenter, p0, 1.0, ns->dpu[dim], ns->stn);
            switch (dim) {
                case 0:
                    lu0 = value0;
                break;
                case 1:
                    lv0 = value0;
                break;
                case 2:
                    lw0 = value0;
                break;
            }
            real value1 = compute_facet_value_at_point(sfdu2[dim], ccenter, p1, 1.0, ns->dpu[dim], ns->stn);
            switch (dim) {
                case 0:
                    lu1 = value1;
                break;
                case 1:
                    lv1 = value1;
                break;
                case 2:
                    lw1 = value1;
                break;
            }
            real value2 = compute_facet_value_at_point(sfdu2[dim], ccenter, p2, 1.0, ns->dpu[dim], ns->stn);
            switch (dim) {
                case 0:
                    lu2 = value2;
                break;
                case 1:
                    lv2 = value2;
                break;
                case 2:
                    lw2 = value2;
                break;
            }
            real value3 = compute_facet_value_at_point(sfdu2[dim], ccenter, p3, 1.0, ns->dpu[dim], ns->stn);
            switch (dim) {
                case 0:
                    lu3 = value3;
                break;
                case 1:
                    lv3 = value3;
                break;
                case 2:
                    lw3 = value3;
                break;
            }
            real value4 = compute_facet_value_at_point(sfdu2[dim], ccenter, p4, 1.0, ns->dpu[dim], ns->stn);
            switch (dim) {
                case 0:
                    lu4 = value4;
                break;
                case 1:
                    lv4 = value4;
                break;
                case 2:
                    lw4 = value4;
                break;
            }
            real value5 = compute_facet_value_at_point(sfdu2[dim], ccenter, p5, 1.0, ns->dpu[dim], ns->stn);
            switch (dim) {
                case 0:
                    lu5 = value5;
                break;
                case 1:
                    lv5 = value5;
                break;
                case 2:
                    lw5 = value5;
                break;
            }
            real value6 = compute_facet_value_at_point(sfdu2[dim], ccenter, p6, 1.0, ns->dpu[dim], ns->stn);
            switch (dim) {
                case 0:
                    lu6 = value6;
                break;
                case 1:
                    lv6 = value6;
                break;
                case 2:
                    lw6 = value6;
                break;
            }
            real value7 = compute_facet_value_at_point(sfdu2[dim], ccenter, p7, 1.0, ns->dpu[dim], ns->stn);
            switch (dim) {
                case 0:
                    lu7 = value7;
                break;
                case 1:
                    lv7 = value7;
                break;
                case 2:
                    lw7 = value7;
                break;
            }
            // printf("Máximo de pontos --> %d\nMínimo de pontos --> %d\nPontos usados --> %d\n",maxpts,wls_num_min_points(DIM,degree),cont);
        }
        fprintf(f, "%e %e %e\n", lu0, lv0, lw0);
        fprintf(f, "%e %e %e\n", lu1, lv1, lw1);
        fprintf(f, "%e %e %e\n", lu2, lv2, lw2);
        fprintf(f, "%e %e %e\n", lu3, lv3, lw3);
        fprintf(f, "%e %e %e\n", lu4, lv4, lw4);
        fprintf(f, "%e %e %e\n", lu5, lv5, lw5);
        fprintf(f, "%e %e %e\n", lu6, lv6, lw6);
        fprintf(f, "%e %e %e\n", lu7, lv7, lw7);
    }
    higcit_destroy(it);

    // Pressure
/*
    fprintf(f, "\nCELL_DATA %ld\nSCALARS p FLOAT\nLOOKUP_TABLE default\n", numcells);
    for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        uniqueid id = hig_get_cid(c);
        int clid = mp_lookup(m, id);
        real val = dp_get_value(ns->dpp, clid);
        fprintf(f, "%e\n", val);
        // fprintf(f, "%e\n", val);
    }
    higcit_destroy(it);
*/

    fprintf(f, "\nCELL_DATA %ld\nSCALARS S00 FLOAT\nLOOKUP_TABLE default\n", numcells);
    for (it = sd_get_domain_celliterator(ns->ed.sdED); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        uniqueid id = hig_get_cid(c);
        int clid = mp_lookup(m, id);
        real val = dp_get_value(ns->ed.ve.dpS[0][0], clid);
        fprintf(f, "%e\n", val);
    }
    higcit_destroy(it);
    fclose(f);
}




// *******************************************************************
// Navier-Stokes Load and Save
// *******************************************************************

void get_nameres(char **source, char *destination) {
    // Find the second to last '/' character
    char *second_last_slash = NULL;
    char *last_slash = NULL;
    char *current = *source;
    
    while (*current != '\0') {
        if (*current == '/') {
            second_last_slash = last_slash;
            last_slash = current;
        }
        current++;
    }

    // Copy the portion of the string up to the second to last '/' character
    if (second_last_slash != NULL) {
        strncpy(destination, *source, second_last_slash - *source);
        destination[second_last_slash - *source] = '\0';
    } else {
        strcpy(destination, *source); // No '/' found, copy the entire string
    }

    // Append "/res/" to the destination string
    strcat(destination, "/res/");
}

// Loading the data files
void higflow_load_data_files(int argc, char *argv[], higflow_solver *ns) {
    // Name of load file
    ns->par.nameload = argv[1];
    argv++;
    // Name of save file
    ns->par.namesave = argv[1];
    argv++;
    // Name of print file
    ns->par.nameprint = argv[1];
    argv++;
    // Name of residuals folder
    get_nameres(&ns->par.namesave, ns->par.nameres);
}

// Loading the controllers and parameters
void higflow_load_controllers_and_parameters_yaml(higflow_solver* ns, int myrank) {
    // Parameters file name
    char ParContr[1024], InitRestart[1024];

    snprintf(InitRestart, sizeof InitRestart, "%s.init.restart.yaml", ns->par.nameload);
    struct fy_document* fydini = NULL;

    snprintf(ParContr, sizeof ParContr, "%s.par.contr.yaml", ns->par.nameload);
    struct fy_document* fyd = NULL;

    fyd = fy_document_build_from_file(NULL, ParContr);
    fydini = fy_document_build_from_file(NULL, InitRestart);

    if (fyd != NULL && fydini != NULL) {
        // Loading the controllers
        int ifd = 0;
        // Load Step, Initial Frame, Initial Time of the Simulation
        ifd += fy_document_scanf(fydini, "/init_par/step %d", &(ns->par.step));
        ns->par.initstep = ns->par.step;
        ifd += fy_document_scanf(fydini, "/init_par/t %lf", &(ns->par.t));
        ifd += fy_document_scanf(fydini, "/init_par/frame %d", &(ns->par.frame));
        ifd += fy_document_scanf(fydini, "/init_par/ts %lf", &(ns->par.ts));
        ifd += fy_document_scanf(fydini, "/init_par/tp %lf", &(ns->par.tp));
        if (ifd != 5) {
            print0f("=+=+=+= Initial Parameters are Missing!!! =+=+=+=\n");
            exit(1);
        }
        char auxchar[1024];
        // Projtype
        ifd += fy_document_scanf(fyd, "/simulation_contr/projtype %s", auxchar);
        if (strcmp(auxchar, "non_incremental") == 0) {
            ns->contr.projtype = NON_INCREMENTAL;
            print0f("=+=+=+= Projection Method: Non Incremental =+=+=+=\n");
        }
        else if (strcmp(auxchar, "incremental") == 0) {
            ns->contr.projtype = INCREMENTAL;
            print0f("=+=+=+= Projection Method: Incremental =+=+=+=\n");
        }
        else {
            print0f("=+=+=+= Projection Method: Invalid =+=+=+=\n");
            exit(1);
        }
        // Flowtype
        ifd += fy_document_scanf(fyd, "/simulation_contr/phasetype %s", auxchar);
        if (strcmp(auxchar, "singlephase") == 0) {
            // Flowtype for Phase 0 (Case singlephase)
            ifd += fy_document_scanf(fyd, "/simulation_contr/flowtype0 %s", auxchar);
            if (strcmp(auxchar, "newtonian") == 0) {
                ns->contr.flowtype = NEWTONIAN;
                print0f("=+=+=+= Flow Type: Newtonian =+=+=+=\n");
            }
            else if (strcmp(auxchar, "generalized_newtonian") == 0) {
                ns->contr.flowtype = GENERALIZED_NEWTONIAN;
                print0f("=+=+=+= Flow Type: Generalized Newtonian =+=+=+=\n");
            }
            else if (strcmp(auxchar, "viscoelastic") == 0) {
                ns->contr.flowtype = VISCOELASTIC;
                print0f("=+=+=+= Flow Type: Viscoelastic =+=+=+=\n");
                
                ifd += fy_document_scanf(fyd, "/simulation_contr/phase0/visc_contr/model %s", auxchar);
                if (strcmp(auxchar, "oldroyd_b") == 0) {
                    ns->ed.ve.contr.model = 0;
                    print0f("=+=+=+= Constitutive Equation Model: Oldroyd-B =+=+=+=\n");
                }
                else if (strcmp(auxchar, "giesekus") == 0) {
                    ns->ed.ve.contr.model = 1;
                    print0f("=+=+=+= Constitutive Equation Model: Giesekus =+=+=+=\n");
                    ifd += fy_document_scanf(fyd, "/simulation_contr/phase0/model_giesekus/alpha %lf", &(ns->ed.ve.par.alpha));
                    print0f("=+=+=+= Alpha: %f =+=+=+=\n", ns->ed.ve.par.alpha);
                }
                else if (strcmp(auxchar, "lptt") == 0) {
                    ns->ed.ve.contr.model = 2;
                    print0f("=+=+=+= Constitutive Equation Model: LPTT =+=+=+=\n");
                    ifd += fy_document_scanf(fyd, "/simulation_contr/phase0/model_lptt/epsilon %lf", &(ns->ed.ve.par.epsilon));
                    ifd += fy_document_scanf(fyd, "/simulation_contr/phase0/model_lptt/psi %lf", &(ns->ed.ve.par.psi));
                    print0f("=+=+=+= Epsilon: %f =+=+=+=\n", ns->ed.ve.par.epsilon);
                    print0f("=+=+=+= Psi: %f =+=+=+=\n", ns->ed.ve.par.psi);
                }
                else if (strcmp(auxchar, "gptt") == 0) {
                    ns->ed.ve.contr.model = 3;
                    print0f("=+=+=+= Constitutive Equation Model: GPTT =+=+=+=\n");
                    ifd += fy_document_scanf(fyd, "/simulation_contr/phase0/model_lptt/alpha_gptt %lf", &(ns->ed.ve.par.alpha_gptt));
                    ifd += fy_document_scanf(fyd, "/simulation_contr/phase0/model_lptt/beta_gptt %lf", &(ns->ed.ve.par.beta_gptt));
                    print0f("=+=+=+= Alpha GPTT: %f =+=+=+=\n", ns->ed.ve.par.alpha_gptt);
                    print0f("=+=+=+= Beta GPTT: %f =+=+=+=\n", ns->ed.ve.par.beta_gptt);
                    //    } else if (strcmp(auxchar,"integral") == 0) {
                    //       ns->ed.ve.contr.model = 4;
                    //       print0f("=+=+=+= Constitutive Equation Model: Integral =+=+=+=\n");
                }
                else if (strcmp(auxchar, "user_set") == 0) {
                    ns->ed.ve.contr.model = -1;
                    print0f("=+=+=+= Constitutive Equation Model: User Set Model =+=+=+=\n");
                }
                else {
                    print0f("=+=+=+= Constitutive Equation Model Invalid =+=+=+=\n");
                    exit(1);
                }

                // Discret Model Flow Type Viscoelastic
                ifd += fy_document_scanf(fyd, "/simulation_contr/phase0/visc_contr/discrtype %s", auxchar);
                if (strcmp(auxchar, "explicit") == 0) {
                    ns->ed.ve.contr.discrtype = EXPLICIT;
                    print0f("=+=+=+= Constitutive Equation Discretization: Explicit =+=+=+=\n");
                }
                else if (strcmp(auxchar, "implicit") == 0) {
                    ns->ed.ve.contr.discrtype = IMPLICIT;
                    print0f("=+=+=+= Constitutive Equation Discretization: Implicit =+=+=+=\n");
                }
                else {
                    print0f("=+=+=+= Constitutive Equation Discretization Invalid =+=+=+=\n");
                    exit(1);
                }

                // Convective Discret Model Flow Type Viscoelastic
                ifd += fy_document_scanf(fyd, "/simulation_contr/phase0/visc_contr/convecdiscrtype %s", auxchar);
                if (strcmp(auxchar, "central") == 0) {
                    ns->ed.ve.contr.convecdiscrtype = CELL_CENTRAL;
                    print0f("=+=+=+= Constitutive Equation Convective Term: Central  =+=+=+=\n");
                }
                else if (strcmp(auxchar, "cubista") == 0) {
                    ns->ed.ve.contr.convecdiscrtype = CELL_CUBISTA;
                    print0f("=+=+=+= Constitutive Equation Convective Term: CUBISTA =+=+=+=\n");
                }
                else {
                    print0f("=+=+=+= Constitutive Equation Convective Term Invalid =+=+=+=\n");
                    exit(1);
                }

                // Common Viscoelastic Parameters
                ifd += fy_document_scanf(fyd, "/simulation_contr/phase0/adimensional/De %lf", &(ns->ed.ve.par.De));
                ifd += fy_document_scanf(fyd, "/simulation_contr/phase0/adimensional/beta %lf", &(ns->ed.ve.par.beta));
                ifd += fy_document_scanf(fyd, "/simulation_contr/phase0/adimensional/kernel_tol %lf", &(ns->ed.ve.par.kernel_tol));

                print0f("=+=+=+= Deborah Number: %f =+=+=+=\n", ns->ed.ve.par.De);
                print0f("=+=+=+= Beta: %f =+=+=+=\n", ns->ed.ve.par.beta);
            }
            //    else {
            //        print0f("=+=+=+= Flow Model Invalid =+=+=+=\n");
            //        exit(1);
            //    }
        }
        else if (strcmp(auxchar, "twophase") == 0) {
            // Flowtype for Phase 0
            ns->contr.flowtype = MULTIPHASE;
            ifd += fy_document_scanf(fyd, "/multiphase_par/Ca %lf", &(ns->ed.mult.par.Ca));
            printf("=+=+=+= Capillary Number: %f =+=+=+=\n", ns->ed.mult.par.Ca);
            ifd += fy_document_scanf(fyd, "/simulation_contr/flowtype0 %s", auxchar);
            printf("------------------------------------------------------ \n");
            printf("=*+=*+=*+=*+=*+=*+=*+= Phase 0 =*+=*+=*+=*+=*+=*+=*+= \n");
            printf("------------------------------------------------------ \n");
            if (strcmp(auxchar, "newtonian") == 0) {
                //ns->contr.flowtype = NEWTONIAN;
                print0f("=+=+=+= Flow Type in Phase 0: Newtonian =+=+=+=\n");
            }
            else if (strcmp(auxchar, "generalized_newtonian") == 0) {
                //ns->contr.flowtype = GENERALIZED_NEWTONIAN;
                print0f("=+=+=+= Flow Type in Phase 0: Generalized Newtonian =+=+=+=\n");
            }
            else if (strcmp(auxchar, "viscoelastic") == 0) {
                //ns->contr.flowtype = VISCOELASTIC;
                print0f("=+=+=+= Flow Type in Phase 0: Viscoelastic =+=+=+=\n");
                // Model Flow Type Viscoelastic
                ifd += fy_document_scanf(fyd, "/simulation_contr/phase0/visc_contr/model %s", auxchar);
                if (strcmp(auxchar, "oldroyd_b") == 0) {
                    ns->ed.ve.contr.model = 0;
                    print0f("=+=+=+= Constitutive Equation Model in Phase 0: Oldroyd-B =+=+=+=\n");
                }
                else if (strcmp(auxchar, "giesekus") == 0) {
                    ns->ed.ve.contr.model = 1;
                    print0f("=+=+=+= Constitutive Equation Model in Phase 0: Giesekus =+=+=+=\n");
                    ifd += fy_document_scanf(fyd, "/simulation_contr/phase0/model_giesekus/alpha %lf", &(ns->ed.ve.par.alpha));
                    print0f("=+=+=+= Alpha: %f =+=+=+=\n", ns->ed.ve.par.alpha);
                }
                else if (strcmp(auxchar, "lptt") == 0) {
                    ns->ed.ve.contr.model = 2;
                    print0f("=+=+=+= Constitutive Equation Model in Phase 0: LPTT =+=+=+=\n");
                    ifd += fy_document_scanf(fyd, "/simulation_contr/phase0/model_lptt/epsilon %lf", &(ns->ed.ve.par.epsilon));
                    ifd += fy_document_scanf(fyd, "/simulation_contr/phase0/model_lptt/psi %lf", &(ns->ed.ve.par.psi));
                    print0f("=+=+=+= Epsilon in Phase 0: %f =+=+=+=\n", ns->ed.ve.par.epsilon);
                    print0f("=+=+=+= Psi in Phase 0: %f =+=+=+=\n", ns->ed.ve.par.psi);
                }
                else if (strcmp(auxchar, "gptt") == 0) {
                    ns->ed.ve.contr.model = 3;
                    print0f("=+=+=+= Constitutive Equation Model in Phase 0: GPTT =+=+=+=\n");
                    // } else if (strcmp(auxchar,"integral") == 0) {
                    //   ns->ed.ve.contr.model = 4;
                    //   print0f("=+=+=+= Constitutive Equation Model in Phase 0: Integral =+=+=+=\n");
                }
                else if (strcmp(auxchar, "user_set") == 0) {
                    ns->ed.ve.contr.model = -1;
                    print0f("=+=+=+= Constitutive Equation Model in Phase 0: User Set Model =+=+=+=\n");
                }
                // else {
                //   print0f("=+=+=+= Constitutive Equation Model Invalid in Phase 0 =+=+=+=\n");
                //   exit(1);
                // }

               // Discret Model Flow Type Viscoelastic
                ifd += fy_document_scanf(fyd, "/simulation_contr/phase0/visc_contr/discrtype %s", auxchar);
                if (strcmp(auxchar, "explicit") == 0) {
                    ns->ed.ve.contr.discrtype = EXPLICIT;
                    print0f("=+=+=+= Constitutive Equation Discretization in Phase 0: Explicit =+=+=+=\n");
                }
                else if (strcmp(auxchar, "giesekus") == 0) {
                    ns->ed.ve.contr.discrtype = IMPLICIT;
                    print0f("=+=+=+= Constitutive Equation Discretization in Phase 0: Implicit =+=+=+=\n");
                }
                else {
                    print0f("=+=+=+= Constitutive Equation Discretization Invalid in Phase 0 =+=+=+=\n");
                    exit(1);
                }

                // Convective Discret Model Flow Type Viscoelastic
                ifd += fy_document_scanf(fyd, "/simulation_contr/phase0/visc_contr/convecdiscrtype %s", auxchar);
                if (strcmp(auxchar, "central") == 0) {
                    ns->ed.ve.contr.convecdiscrtype = CELL_CENTRAL;
                    print0f("=+=+=+= Constitutive Equation Convective Term in Phase 0: Central  =+=+=+=\n");
                }
                else if (strcmp(auxchar, "cubista") == 0) {
                    ns->ed.ve.contr.convecdiscrtype = CELL_CUBISTA;
                    print0f("=+=+=+= Constitutive Equation Convective Term in Phase 0: CUBISTA =+=+=+=\n");
                }
                else {
                    print0f("=+=+=+= Constitutive Equation Convective Term Invalid in Phase 0 =+=+=+=\n");
                    exit(1);
                }

                // Common Viscoelastic Parameters
                ifd += fy_document_scanf(fyd, "/simulation_contr/phase0/adimensional/De %lf", &(ns->ed.ve.par.De));
                ifd += fy_document_scanf(fyd, "/simulation_contr/phase0/adimensional/beta %lf", &(ns->ed.ve.par.beta));
                ifd += fy_document_scanf(fyd, "/simulation_contr/phase0/adimensional/kernel_tol %lf", &(ns->ed.ve.par.kernel_tol));

                print0f("=+=+=+= Deborah Number in Phase 0: %f =+=+=+=\n", ns->ed.ve.par.De);
                print0f("=+=+=+= Beta in Phase 0: %f =+=+=+=\n", ns->ed.ve.par.beta);
                print0f("=+=+=+= Kernel_tol in Phase 0: %f =+=+=+=\n", ns->ed.ve.par.kernel_tol);

            }
            else {
                if (myrank == 0) {
                    printf("=+=+=+= Flow Type Invalid in Phase 0 =+=+=+=\n");
                }
                exit(1);
            }

            // Flowtype for Phase 1
            printf("------------------------------------------------------ \n");
            printf("=*+=*+=*+=*+=*+=*+=*+= Phase 1 =*+=*+=*+=*+=*+=*+=*+= \n");
            printf("------------------------------------------------------ \n");
            ifd += fy_document_scanf(fyd, "/simulation_contr/flowtype1 %s", auxchar);
            if (strcmp(auxchar, "newtonian") == 0) {
                //ns->contr.flowtype = NEWTONIAN;
                print0f("=+=+=+= Flow Type in Phase 1: Newtonian =+=+=+=\n");
            }
            else if (strcmp(auxchar, "generalized_newtonian") == 0) {
                //ns->contr.flowtype = GENERALIZED_NEWTONIAN;
                print0f("=+=+=+= Flow Type in Phase 1: Generalized Newtonian =+=+=+=\n");
            }
            else if (strcmp(auxchar, "viscoelastic") == 0) {
                //ns->contr.flowtype = VISCOELASTIC;
                print0f("=+=+=+= Flow Type in Phase 1: Viscoelastic =+=+=+=\n");
                // Model Flow Type Viscoelastic
                ifd += fy_document_scanf(fyd, "/simulation_contr/phase1/visc_contr/model %s", auxchar);
                if (strcmp(auxchar, "oldroyd_b") == 0) {
                    ns->ed.ve.contr.model = 0;
                    print0f("=+=+=+= Constitutive Equation Model in Phase 1: Oldroyd-B =+=+=+=\n");
                }
                else if (strcmp(auxchar, "giesekus") == 0) {
                    ns->ed.ve.contr.model = 1;
                    print0f("=+=+=+= Constitutive Equation Model in Phase 1: Giesekus =+=+=+=\n");
                    ifd += fy_document_scanf(fyd, "/simulation_contr/phase1/model_giesekus/alpha %lf", &(ns->ed.ve.par.alpha));
                    print0f("=+=+=+= Alpha: %f =+=+=+=\n", ns->ed.ve.par.alpha);
                }
                else if (strcmp(auxchar, "lptt") == 0) {
                    ns->ed.ve.contr.model = 2;
                    print0f("=+=+=+= Constitutive Equation Model in Phase 1: LPTT =+=+=+=\n");
                    ifd += fy_document_scanf(fyd, "/simulation_contr/phase1/model_lptt/epsilon %lf", &(ns->ed.ve.par.epsilon));
                    ifd += fy_document_scanf(fyd, "/simulation_contr/phase1/model_lptt/psi %lf", &(ns->ed.ve.par.psi));
                    print0f("=+=+=+= Epsilon in Phase 1: %f =+=+=+=\n", ns->ed.ve.par.epsilon);
                    print0f("=+=+=+= Psi in Phase 1: %f =+=+=+=\n", ns->ed.ve.par.psi);
                }
                else if (strcmp(auxchar, "gptt") == 0) {
                    ns->ed.ve.contr.model = 3;
                    print0f("=+=+=+= Constitutive Equation Model in Phase 1: GPTT =+=+=+=\n");
                    //    } else if (strcmp(auxchar,"integral") == 0) {
                    //       ns->ed.ve.contr.model = 4;
                    //       print0f("=+=+=+= Constitutive Equation Model in Phase 1: Integral =+=+=+=\n");
                }
                else if (strcmp(auxchar, "user_set") == 0) {
                    ns->ed.ve.contr.model = -1;
                    print0f("=+=+=+= Constitutive Equation Model in Phase 1: User Set Model =+=+=+=\n");
                }
                //    else {
                //       print0f("=+=+=+= Constitutive Equation Model Invalid in Phase 1 =+=+=+=\n");
                //       exit(1);
                //    }

                   // Discret Model Flow Type Viscoelastic
                ifd += fy_document_scanf(fyd, "/simulation_contr/phase1/visc_contr/discrtype %s", auxchar);
                if (strcmp(auxchar, "explicit") == 0) {
                    ns->ed.ve.contr.discrtype = EXPLICIT;
                    print0f("=+=+=+= Constitutive Equation Discretization in Phase 1: Explicit =+=+=+=\n");
                }
                else if (strcmp(auxchar, "giesekus") == 0) {
                    ns->ed.ve.contr.discrtype = IMPLICIT;
                    print0f("=+=+=+= Constitutive Equation Discretization in Phase 1: Implicit =+=+=+=\n");
                }
                else {
                    print0f("=+=+=+= Constitutive Equation Discretization Invalid in Phase 1 =+=+=+=\n");
                    exit(1);
                }

                // Convective Discret Model Flow Type Viscoelastic
                ifd += fy_document_scanf(fyd, "/simulation_contr/phase1/visc_contr/convecdiscrtype %s", auxchar);
                if (strcmp(auxchar, "central") == 0) {
                    ns->ed.ve.contr.convecdiscrtype = CELL_CENTRAL;
                    print0f("=+=+=+= Constitutive Equation Convective Term in Phase 1: Central  =+=+=+=\n");
                }
                else if (strcmp(auxchar, "cubista") == 0) {
                    ns->ed.ve.contr.convecdiscrtype = CELL_CUBISTA;
                    print0f("=+=+=+= Constitutive Equation Convective Term in Phase 1: CUBISTA =+=+=+=\n");
                }
                else {
                    print0f("=+=+=+= Constitutive Equation Convective Term Invalid in Phase 1 =+=+=+=\n");
                    exit(1);
                }

                // Common Viscoelastic Parameters
                ifd += fy_document_scanf(fyd, "/simulation_contr/phase1/adimensional/De %lf", &(ns->ed.ve.par.De));
                ifd += fy_document_scanf(fyd, "/simulation_contr/phase1/adimensional/beta %lf", &(ns->ed.ve.par.beta));
                ifd += fy_document_scanf(fyd, "/simulation_contr/phase1/adimensional/kernel_tol %lf", &(ns->ed.ve.par.kernel_tol));

                print0f("=+=+=+= Deborah Number in Phase 1: %f =+=+=+=\n", ns->ed.ve.par.De);
                print0f("=+=+=+= Beta in Phase 1: %f =+=+=+=\n", ns->ed.ve.par.beta);
                print0f("=+=+=+= Kernel_tol in Phase 1: %f =+=+=+=\n", ns->ed.ve.par.kernel_tol);
                printf("------------------------------------------------------ \n");
            }
            else {
                print0f("=+=+=+= Flow Type Invalid in Phase 1 =+=+=+=\n");
                exit(1);
            }
        }
        else {
            print0f("=+=+=+= Flow Type Invalid =+=+=+=\n");
            exit(1);
        }

        // Model Flow Type
        ifd += fy_document_scanf(fyd, "/simulation_contr/eoflow %s", auxchar);
        if (strcmp(auxchar, "none") == 0) {
            ns->contr.eoflow = false;
            print0f("    =+= Flow Model Type: ---------------- =+=\n");
        }
        else if (strcmp(auxchar, "electroosmotic") == 0) {
            ns->contr.eoflow = true;
            print0f("=+= Flow Model Type: Electroosmotic =+=\n");
        }
        // Time Discret Type
        ifd += fy_document_scanf(fyd, "/simulation_contr/tempdiscrtype %s", auxchar);

        if (strcmp(auxchar, "explicit_euler") == 0) {
            ns->contr.tempdiscrtype = EXPLICIT_EULER;
            print0f("=+=+=+= Temporal Discretization: Explicit Euler =+=+=+=\n");
        }
        else if (strcmp(auxchar, "explicit_rk2") == 0) {
            ns->contr.tempdiscrtype = EXPLICIT_RK2;
            print0f("=+=+=+= Temporal Discretization: Runge-Kutta 2 =+=+=+=\n");
        }
        else if (strcmp(auxchar, "explicit_rk3") == 0) {
            ns->contr.tempdiscrtype = EXPLICIT_RK3;
            print0f("=+=+=+= Temporal Discretization: Runge-Kutta 3 =+=+=+=\n");
        }
        else if (strcmp(auxchar, "semi_implicit_euler") == 0) {
            ns->contr.tempdiscrtype = SEMI_IMPLICIT_EULER;
            print0f("=+=+=+= Temporal Discretization: Semi-Implicit Euler =+=+=+=\n");
        }
        else if (strcmp(auxchar, "semi_implicit_crank_nicolson") == 0) {
            ns->contr.tempdiscrtype = SEMI_IMPLICIT_CN;
            print0f("=+=+=+= Temporal Discretization: Semi-Implicit Crank-Nicolson =+=+=+=\n");
        }
        else if (strcmp(auxchar, "semi_implicit_bdf2") == 0) {
            ns->contr.tempdiscrtype = SEMI_IMPLICIT_BDF2;
            print0f("=+=+=+= Temporal Discretization: Semi-Implicit BDF2 =+=+=+=\n");
        }
        // Spatial Discret Type
        ifd += fy_document_scanf(fyd, "/simulation_contr/spatialdiscrtype %s", auxchar);
        if (strcmp(auxchar, "second_order") == 0) {
            ns->contr.spatialdiscrtype = ORDER2;
            print0f("=+=+=+= Spatial Discretization: Second Order =+=+=+=\n");
        }
        else if (strcmp(auxchar, "forth_order") == 0) {
            ns->contr.spatialdiscrtype = ORDER4;
            print0f("=+=+=+= Spatial Discretization: Forth Order =+=+=+=\n");
        }
        // Convective Discret Type
        ifd += fy_document_scanf(fyd, "/simulation_contr/convecdiscrtype %s", auxchar);
        if (strcmp(auxchar, "central") == 0) {
            ns->contr.convecdiscrtype = CENTRAL;
            print0f("=+=+=+= Convective Scheme: Central =+=+=+=\n");
        }
        else if (strcmp(auxchar, "first_order") == 0) {
            ns->contr.convecdiscrtype = FIRST_ORDER;
            print0f("=+=+=+= Convective Scheme: First Order Upwind =+=+=+=\n");
        }
        else if (strcmp(auxchar, "second_order") == 0) {
            ns->contr.convecdiscrtype = ORDER2;
            print0f("=+=+=+= Convective Scheme: Second Order Upwind =+=+=+=\n");
            ifd += fy_document_scanf(fyd, "/simulation_contr/secondconvecdiscrtype %s", auxchar);
            if (strcmp(auxchar, "modified_coefficient_upwind") == 0) {
                ns->contr.secondconvecdiscrtype = 0;
                printf("    =+= Convective Scheme Type : Modified Coefficient Upwind =+=\n");
            }
            else if (strcmp(auxchar, "cubista") == 0) {
                ns->contr.secondconvecdiscrtype = 1;
                printf("    =+= Convective Scheme Type : Cubista =+=\n");
            }
            else if (strcmp(auxchar, "quick") == 0) {
                ns->contr.secondconvecdiscrtype = 2;
                printf("    =+= Convective Scheme Type : Quick =+=\n");
            }
        }
        // Reynold's Number
        ifd += fy_document_scanf(fyd, "/adimensional_par/Re %lf", &(ns->par.Re));
        print0f("=+=+=+= Reynolds Number: %f =+=+=+=\n", ns->par.Re);
        // The Parameters of the Simulation
        ifd += fy_document_scanf(fyd, "/simulation_par/numsteps %d", &(ns->par.finalstep));
        ifd += fy_document_scanf(fyd, "/simulation_par/dt %lf", &(ns->par.dt));
        ifd += fy_document_scanf(fyd, "/simulation_par/dts %lf", &(ns->par.dts));
        ifd += fy_document_scanf(fyd, "/simulation_par/dtp %lf", &(ns->par.dtp));
    }
    else {
        // Error in open the file
        print0f("=+=+=+= Error loading file %s or %s =+=+=+=\n", ParContr, InitRestart);
        exit(1);
    }
    fy_document_destroy(fyd);
    fy_document_destroy(fydini);
}

// Loading the controllers
void higflow_load_controllers(higflow_solver *ns, int myrank) {
    // Controllers file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.contr", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL) {
        // Loading the controllers
        int ifd;
        ifd = fscanf(fd, "%d", (int *)&(ns->contr.projtype));
        ifd = fscanf(fd, "%d", (int *)&(ns->contr.flowtype));
        ifd = fscanf(fd, "%d", &(ns->contr.eoflow));
        ifd = fscanf(fd, "%d", (int *)&(ns->contr.tempdiscrtype));
        ifd = fscanf(fd, "%d", (int *)&(ns->contr.spatialdiscrtype));
        ifd = fscanf(fd, "%d", (int *)&(ns->contr.convecdiscrtype));
        ifd = fscanf(fd, "%d", (int *)&(ns->contr.secondconvecdiscrtype));
        fclose(fd);

        // Printing the controllers
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (myrank == 0) {
           printf("=+=+=+= Controllers =+=+=+=\n");
           switch (ns->contr.projtype) {
               case NON_INCREMENTAL:
                    printf("=+=+=+= Projection Method: Non Incremental =+=+=+=\n");
                    break;
                case INCREMENTAL:
                    printf("=+=+=+= Projection Method: Incremental =+=+=+=\n");
                    break;
            }
            switch (ns->contr.flowtype) {
                case NEWTONIAN:
                    printf("=+=+=+= Flow Model: Newtonian =+=+=+=\n");
                    break;
                case GENERALIZED_NEWTONIAN:
                    printf("=+= Flow Model: Generalized Newtonian =+=\n");
                    break;
                case MULTIPHASE:
                    printf("=+=+=+= Flow Model: Multifase =+=+=+=\n");
                    break;
                case VISCOELASTIC:
                    printf("=+=+=+= Flow Model: Viscoelastic =+=+=+=\n");
                    break;
            }
            switch (ns->contr.eoflow) {
                case true:
                    printf("=+= Flow Model Type: Eletroosmotic =+=\n");
                    break;
                case false:
                    break;
                default:
                    printf("=+=+=+= Unsupported eoflow!! =+=+=+=\n");
                    exit(1);
                    break;
            }
            switch (ns->contr.tempdiscrtype) {
                case EXPLICIT_EULER:
                    printf("=+=+=+= Temporal Discretization: Explicit Euler =+=+=+=\n");
                    break;
                case EXPLICIT_RK2:
                    printf("=+=+=+= Temporal Discretization: Runge-Kutta 2 =+=+=+=\n");
                    break;
                case EXPLICIT_RK3:
                    printf("=+=+=+= Temporal Discretization: Runge-Kutta 3 =+=+=+=\n");
                    break;
                case SEMI_IMPLICIT_EULER:
                    printf("=+=+=+= Temporal Discretization: Semi-Implicit Euler =+=+=+=\n");
                    break;
                case SEMI_IMPLICIT_CN:
                    printf("=+=+=+= Temporal Discretization: Semi-Implicit Crank-Nicolson =+=+=+=\n");
                    break;
                case SEMI_IMPLICIT_BDF2:
                    printf("=+=+=+= Temporal Discretization: Semi-Implicit BDF2 =+=+=+=\n");
                    break;
            }
            switch (ns->contr.spatialdiscrtype) {
                case ORDER2:
                    printf("=+=+=+= Spatial Discretization: Second Order =+=+=+=\n");
                    break;
                case ORDER4:
                    printf("=+=+=+= Spatial Discretization: Forth Order =+=+=+=\n");
                    break;
            }
            switch (ns->contr.convecdiscrtype) {
                case CENTRAL:
                    printf("=+=+=+= Convective Scheme: Central =+=+=+=\n");
                    break;
                case FIRST_ORDER:
                    printf("=+=+=+= Convective Scheme: First Order Upwind =+=+=+=\n");
                    break;
                case SECOND_ORDER:
                    printf("=+=+=+= Convective Scheme: Second Order Upwind =+=+=+=\n");
                    switch (ns->contr.secondconvecdiscrtype) {
                        case MCU:
                            printf("=+= Convective Scheme Type : Modified Coefficient Upwind =+=\n");
                            break;
                        case CUBISTA:
                            printf("=+= Convective Scheme Type : Cubista =+=\n");
                            break;
                        case QUICK:
                            printf("=+= Convective Scheme Type : Quick =+=\n");
                            break;
                    }
                    break;
            }
        }
    } else {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the controllers
void higflow_save_controllers(higflow_solver *ns, int myrank) {
    // Controllers file name
    if (myrank == 0) {
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.contr", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL) {
            // Saving the controllers
            fprintf(fd, "%d\n", (int)(ns->contr.projtype));
            fprintf(fd, "%d\n", (int)(ns->contr.flowtype));
            fprintf(fd, "%d\n", (ns->contr.eoflow));
            fprintf(fd, "%d\n", (int)(ns->contr.tempdiscrtype));
            fprintf(fd, "%d\n", (int)(ns->contr.spatialdiscrtype));
            fprintf(fd, "%d\n", (int)(ns->contr.convecdiscrtype));
            fprintf(fd, "%d\n", (int)(ns->contr.secondconvecdiscrtype));
            fclose(fd);
        } else {
            // Error in open the file
            printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}

// Loading the parameters
void higflow_load_parameters(higflow_solver *ns, int myrank) {
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.par", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL) {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%d", &(ns->par.step));
        ns->par.initstep = ns->par.step;
        ifd = fscanf(fd, "%d", &(ns->par.finalstep));
        ifd = fscanf(fd, "%lf", &(ns->par.t));
        ifd = fscanf(fd, "%lf", &(ns->par.dt));
        ifd = fscanf(fd, "%lf", &(ns->par.Re));
        ifd = fscanf(fd, "%lf", &(ns->par.dts));
        ifd = fscanf(fd, "%lf", &(ns->par.dtp));
        ifd = fscanf(fd, "%d", &(ns->par.frame));
        ifd = fscanf(fd, "%lf", &(ns->par.ts));
        ifd = fscanf(fd, "%lf", &(ns->par.tp));
        fclose(fd);
        if (myrank == 0) {
            printf("=+=+=+= Reynolds Number: %f =+=+=+=\n", ns->par.Re);
        }
    } else {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the parameters
void higflow_save_parameters(higflow_solver *ns, int myrank) {
    if (myrank == 0) {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.par", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL) {
            // Saving the parameters
            fprintf(fd, "%d\n", (ns->par.step));
            fprintf(fd, "%d\n", (ns->par.finalstep));
            fprintf(fd, "%lf\n", (ns->par.t));
            fprintf(fd, "%lf\n", (ns->par.dt));
            fprintf(fd, "%lf\n", (ns->par.Re));
            fprintf(fd, "%lf\n", (ns->par.dts));
            fprintf(fd, "%lf\n", (ns->par.dtp));
            fprintf(fd, "%d\n", (ns->par.frame));
            fprintf(fd, "%lf\n", (ns->par.ts));
            fprintf(fd, "%lf\n", (ns->par.tp));
            fclose(fd);
        } else {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}

// Loading the properties
void higflow_load_properties(higflow_solver *ns, int myrank, int ntasks) {

    char namefile[1024];
    
    // Velocity filename base
    snprintf(namefile, sizeof namefile, "%s._dpu", ns->par.namesave);
    load_fdp_vec(ns->psfdu, ns->dpu, namefile, myrank, ntasks);

    // Pressure file name
    snprintf(namefile, sizeof namefile, "%s._dpp", ns->par.namesave);
    load_dp_scalar(ns->psdp, ns->dpp, namefile, myrank, ntasks);

    // Cell source term file name
    snprintf(namefile, sizeof namefile, "%s._dpF", ns->par.namesave);
    load_dp_scalar(ns->psdF, ns->dpF, namefile, myrank, ntasks);

    // Facet source term file name base
    snprintf(namefile, sizeof namefile, "%s._dpFU", ns->par.namesave);
    load_fdp_vec(ns->psfdF, ns->dpFU, namefile, myrank, ntasks);

    switch (ns->contr.flowtype) {
        case GENERALIZED_NEWTONIAN:
            // Generalized Newtonian Viscosity file name
            snprintf(namefile, sizeof namefile, "%s._dpvisc_gn", ns->par.namesave);
            load_dp_scalar(ns->ed.psdED, ns->ed.gn.dpvisc, namefile, myrank, ntasks);
        break;
        case MULTIPHASE:
            // Volume Fraction file name
            snprintf(namefile, sizeof namefile, "%s._dpfracvol", ns->par.namesave);
            load_dp_scalar(ns->ed.psdED, ns->ed.mult.dpfracvol, namefile, myrank, ntasks);

            // Multiphase Viscosity file name
            snprintf(namefile, sizeof namefile, "%s._dpvisc_mult", ns->par.namesave);
            load_dp_scalar(ns->ed.psdED, ns->ed.mult.dpvisc, namefile, myrank, ntasks);

            // Multiphase Density file name
            snprintf(namefile, sizeof namefile, "%s._dpdens", ns->par.namesave);
            load_dp_scalar(ns->ed.psdED, ns->ed.mult.dpdens, namefile, myrank, ntasks);

            if(ns->ed.mult.contr.flowtype_either == VISCOELASTIC) {
                // Multiphase Kernel Tensor file name base
                snprintf(namefile, sizeof namefile, "%s._dpKernel_mult", ns->par.namesave);
                load_dp_tensor(ns->ed.psdED, ns->ed.mult.ve.dpKernel, namefile, myrank, ntasks);
            }
        break;
        case VISCOELASTIC:
            // Viscoelastic Kernel Tensor file name
            snprintf(namefile, sizeof namefile, "%s._dpKernel", ns->par.namesave);
            load_dp_tensor(ns->ed.psdED, ns->ed.ve.dpKernel, namefile, myrank, ntasks);
        break;
        case VISCOELASTIC_INTEGRAL:
            // Viscoelastic Integral Elastic (S) Tensor file name
            snprintf(namefile, sizeof namefile, "%s._dpS_im", ns->par.namesave);
            load_dp_tensor(ns->ed.psdED, ns->ed.im.dpS, namefile, myrank, ntasks);

            for(int k=0; k <= NDT; k++){
                // Viscoelastic Integral Finger Tensor (B) file name base at a certain time frame k
                snprintf(namefile, sizeof namefile, "%s._dpB%d", ns->par.namesave, k);
                load_dp_tensor(ns->ed.psdED, ns->ed.im.dpB[k], namefile, myrank, ntasks);
            }
        break;
    }

    if (ns->contr.eoflow == true) {

        // Electroosmotic Applied Potential Phi file name
        snprintf(namefile, sizeof namefile, "%s._dpphi", ns->par.namesave);
        load_dp_scalar(ns->ed.eo.psdEOphi, ns->ed.eo.dpphi, namefile, myrank, ntasks);

        // Electroosmotic Induced Potential (zeta potential) Psi file name
        snprintf(namefile, sizeof namefile, "%s._dppsi", ns->par.namesave);
        load_dp_scalar(ns->ed.eo.psdEOpsi, ns->ed.eo.dppsi, namefile, myrank, ntasks);

        // Electroosmotic Positive charge concentration Nplus file name
        snprintf(namefile, sizeof namefile, "%s._dpnplus", ns->par.namesave);
        load_dp_scalar(ns->ed.eo.psdEOnplus, ns->ed.eo.dpnplus, namefile, myrank, ntasks);

        // Electroosmotic Negative charge concentration Nminus file name
        snprintf(namefile, sizeof namefile, "%s._dpnminus", ns->par.namesave);
        load_dp_scalar(ns->ed.eo.psdEOnminus, ns->ed.eo.dpnminus, namefile, myrank, ntasks);

        // Electroosmotic source term file name
        snprintf(namefile, sizeof namefile, "%s._dpFeo", ns->par.namesave);
        load_fdp_vec(ns->ed.eo.psfdEOFeo, ns->ed.eo.dpFeo, namefile, myrank, ntasks);

    }

}

// Saving the properties
void higflow_save_properties(higflow_solver *ns, int myrank, int ntasks) {
    
    char namefile[1024];

    // Velocity filename base
    snprintf(namefile, sizeof namefile, "%s._dpu", ns->par.namesave);
    save_fdp_vec(ns->psfdu, ns->dpu, namefile, myrank, ntasks);

    // Pressure file name
    snprintf(namefile, sizeof namefile, "%s._dpp", ns->par.namesave);
    save_dp_scalar(ns->psdp, ns->dpp, namefile, myrank, ntasks);

    // Cell source term file name
    snprintf(namefile, sizeof namefile, "%s._dpF", ns->par.namesave);
    save_dp_scalar(ns->psdF, ns->dpF, namefile, myrank, ntasks);

    // Facet source term file name base
    snprintf(namefile, sizeof namefile, "%s._dpFU", ns->par.namesave);
    save_fdp_vec(ns->psfdF, ns->dpFU, namefile, myrank, ntasks);

    switch (ns->contr.flowtype) {
        case GENERALIZED_NEWTONIAN:
            // Generalized Newtonian Viscosity file name
            snprintf(namefile, sizeof namefile, "%s._dpvisc_gn", ns->par.namesave);
            save_dp_scalar(ns->ed.psdED, ns->ed.gn.dpvisc, namefile, myrank, ntasks);
        break;
        case MULTIPHASE:
            // Volume Fraction file name
            snprintf(namefile, sizeof namefile, "%s._dpfracvol", ns->par.namesave);
            save_dp_scalar(ns->ed.psdED, ns->ed.mult.dpfracvol, namefile, myrank, ntasks);

            // Multiphase Viscosity file name
            snprintf(namefile, sizeof namefile, "%s._dpvisc_mult", ns->par.namesave);
            save_dp_scalar(ns->ed.psdED, ns->ed.mult.dpvisc, namefile, myrank, ntasks);

            // Multiphase Density file name
            snprintf(namefile, sizeof namefile, "%s._dpdens", ns->par.namesave);
            save_dp_scalar(ns->ed.psdED, ns->ed.mult.dpdens, namefile, myrank, ntasks);

            if(ns->ed.mult.contr.flowtype_either == VISCOELASTIC) {
                // Multiphase Kernel Tensor file name base
                snprintf(namefile, sizeof namefile, "%s._dpKernel_mult", ns->par.namesave);
                save_dp_tensor(ns->ed.psdED, ns->ed.mult.ve.dpKernel, namefile, myrank, ntasks);
            }
        break;
        case VISCOELASTIC:
            // Viscoelastic Kernel Tensor file name
            snprintf(namefile, sizeof namefile, "%s._dpKernel", ns->par.namesave);
            save_dp_tensor(ns->ed.psdED, ns->ed.ve.dpKernel, namefile, myrank, ntasks);
        break;
        case VISCOELASTIC_INTEGRAL:
            // Viscoelastic Integral Elastic (S) Tensor file name
            snprintf(namefile, sizeof namefile, "%s._dpS_im", ns->par.namesave);
            save_dp_tensor(ns->ed.psdED, ns->ed.im.dpS, namefile, myrank, ntasks);

            for(int k=0; k <= NDT; k++){
                // Viscoelastic Integral Finger Tensor (B) file name base at a certain time frame k
                snprintf(namefile, sizeof namefile, "%s._dpB%d", ns->par.namesave, k);
                save_dp_tensor(ns->ed.psdED, ns->ed.im.dpB[k], namefile, myrank, ntasks);
            }
        break;
    }

    if (ns->contr.eoflow == true) {

        // Electroosmotic Applied Potential Phi file name
        snprintf(namefile, sizeof namefile, "%s._dpphi", ns->par.namesave);
        save_dp_scalar(ns->ed.eo.psdEOphi, ns->ed.eo.dpphi, namefile, myrank, ntasks);

        // Electroosmotic Induced Potential (zeta potential) Psi file name
        snprintf(namefile, sizeof namefile, "%s._dppsi", ns->par.namesave);
        save_dp_scalar(ns->ed.eo.psdEOpsi, ns->ed.eo.dppsi, namefile, myrank, ntasks);

        // Electroosmotic Positive charge concentration Nplus file name
        snprintf(namefile, sizeof namefile, "%s._dpnplus", ns->par.namesave);
        save_dp_scalar(ns->ed.eo.psdEOnplus, ns->ed.eo.dpnplus, namefile, myrank, ntasks);

        // Electroosmotic Negative charge concentration Nminus file name
        snprintf(namefile, sizeof namefile, "%s._dpnminus", ns->par.namesave);
        save_dp_scalar(ns->ed.eo.psdEOnminus, ns->ed.eo.dpnminus, namefile, myrank, ntasks);

        // Electroosmotic source term file name
        snprintf(namefile, sizeof namefile, "%s._dpFeo", ns->par.namesave);
        save_fdp_vec(ns->ed.eo.psfdEOFeo, ns->ed.eo.dpFeo, namefile, myrank, ntasks);

    }

}


// create the destination name by concatenating the save name up to the last slash 
// with the load name after the last slash
void get_amrfilename_save(char (*load)[1024], char (*save)[1024], char *destination) {
    char *lastslash_save = strrchr(*save, '/');

    char *lastslash_load = strrchr(*load, '/');

    // // get second to last slash position in load
    // char *secondlastslash_load = NULL;
    // char *lastslash_load = NULL;
    // char *current = *load;
    
    // while (*current != '\0') {
    //     if (*current == '/') {
    //         secondlastslash_load = lastslash_load;
    //         lastslash_load = current;
    //     }
    //     current++;
    // }
    
    if (lastslash_save != NULL && lastslash_load != NULL) {
        int len = lastslash_save - *save;
        strncpy(destination, *save, len);
        destination[len] = '\0';
        strcat(destination, lastslash_load);
    } else {
        printf("=+=+=+= Error in get_amrfilename_save =+=+=+=\n");
        exit(1);
    }
}


void higflow_save_domain(higflow_solver *ns, int myrank, int ntasks) {
    if(myrank == 0) { 
        // Loading the domain data
        char namefile_load[1024];
        sprintf(namefile_load,"%s.domain",ns->par.nameload);
        char namefile_save[1024];
        sprintf(namefile_save, "%s.domain", ns->par.namesave);

        int samefile = 0;
        if(strcmp(namefile_load, namefile_save) == 0) {
            printf("=+=+=+= Load and Save names are the same - not changing the domain file =+=+=+=\n");
            samefile = 1;
        }


        FILE *fdomain_load = fopen(namefile_load, "r");
        if (fdomain_load == NULL) {
            // Error in open the file
            printf("=+=+=+= Error loading file %s =+=+=+=\n",namefile_load);
            exit(1);
        }

        FILE *fdomain_save;
        if(samefile == 0) {
            fdomain_save = fopen(namefile_save, "w");
            if (fdomain_save == NULL) {
                // Error in open the file
                printf("=+=+=+= Error saving file %s =+=+=+=\n",namefile_save);
                exit(1);
            }
        }
        
        // Number of HigTrees
        int numhigs;
        int ifd = fscanf(fdomain_load,"%d\n",&numhigs);

        // write number of higtrees
        if(samefile == 0) fprintf(fdomain_save,"%d\n",numhigs);

        for(int h = 0; h < numhigs; h++) {
            // Name of the HigTree file
            char amrfilename_load[1024];
            __higflow_readstring(amrfilename_load,1024,fdomain_load);

            // Open the AMR format file
            FILE *fd_load = fopen(amrfilename_load, "r");
            if (fd_load == NULL) {
                // Error in open the file
                printf("=+=+=+= Error loading file %s =+=+=+=\n",amrfilename_load);
                exit(1);
            }
            
            char ch;
            char amrfilename_save[1024];
            get_amrfilename_save(&amrfilename_load, &namefile_save, amrfilename_save);

            int samefile_amr = 0;
            if(strcmp(amrfilename_load, amrfilename_save) == 0) {
                printf("=+=+=+= Load and Save names are the same - not changing the domain %d amr file =+=+=+=\n", h);
                samefile_amr = 1;
            }

            // copy the corresponding filename of the hig
            if(samefile == 0) fprintf(fdomain_save,"%s\n",amrfilename_save);

            // Copying the higtree information to the file 
            if(samefile_amr == 0) {
                printf("=+=+=+= Saving domain %d to %s =+=+=+=\n", h, amrfilename_save);

                FILE *fd_save = fopen(amrfilename_save, "w");
                if (fd_save == NULL) {
                    // Error in open the file
                    printf("=+=+=+= Error saving file %s =+=+=+=\n",amrfilename_save);
                    exit(1);
                }

                while ((ch = fgetc(fd_load)) != EOF) {
                    fputc(ch, fd_save);
                }

                fclose(fd_save);
            }
            
            // Close the AMR format files
            fclose(fd_load);
            
        }
        fclose(fdomain_load);
        if(samefile == 0) fclose(fdomain_save);

    }
}

void higflow_save_boundaries(higflow_solver *ns, int myrank, int ntasks) {
    if(myrank == 0) { 
        // Loading the boundary data
        char namefile_load[1024];
        sprintf(namefile_load,"%s.bc",ns->par.nameload);
        char namefile_save[1024];
        sprintf(namefile_save, "%s.bc", ns->par.namesave);

        int samefile = 0;
        if(strcmp(namefile_load, namefile_save) == 0) {
            printf("=+=+=+= Load and Save names are the same - not changing the bc file =+=+=+=\n");
            samefile = 1;
        }

        FILE *fboundary_load = fopen(namefile_load, "r");
        if (fboundary_load == NULL) {
            // Error in open the file
            printf("=+=+=+= Error loading file %s =+=+=+=\n",namefile_load);
            exit(1);
        }

        FILE *fboundary_save;
        if(samefile == 0) {
            fboundary_save = fopen(namefile_save, "w");
            if (fboundary_save == NULL) {
                // Error in open the file
                printf("=+=+=+= Error saving file %s =+=+=+=\n",namefile_save);
                exit(1);
            }
        }
        
        // Number of HigTrees
        int numbcs;
        int ifd = fscanf(fboundary_load,"%d\n",&numbcs);

        // write number of higtrees
        if(samefile == 0) fprintf(fboundary_save,"%d\n",numbcs);

        for(int h = 0; h < numbcs; h++) {
            // write boundary id
            int bcid;
            ifd = fscanf(fboundary_load,"%d\n",&bcid);
            if(samefile == 0) fprintf(fboundary_save,"%d\n",bcid);

            // Name of the HigTree file
            char amrfilename_load[1024];
            __higflow_readstring(amrfilename_load,1024,fboundary_load);

            // Open the AMR format file
            FILE *fd_load = fopen(amrfilename_load, "r");
            if (fd_load == NULL) {
                // Error in open the file
                printf("=+=+=+= Error loading file %s =+=+=+=\n",amrfilename_load);
                exit(1);
            }
            
            char ch;
            char amrfilename_save[1024];
            get_amrfilename_save(&amrfilename_load, &namefile_save, amrfilename_save);

            int samefile_amr = 0;
            if(strcmp(amrfilename_load, amrfilename_save) == 0) {
                printf("=+=+=+= Load and Save names are the same - not changing the boundary %d amr file =+=+=+=\n", h);
                samefile_amr = 1;
            }

            // copy the corresponding filename of the hig
            if(samefile == 0) fprintf(fboundary_save,"%s\n",amrfilename_save);

            if(samefile_amr == 0) {
                printf("=+=+=+= Saving boundary %d to %s =+=+=+=\n", h, amrfilename_save);

                // Copy
                FILE *fd_save = fopen(amrfilename_save, "w");
                if (fd_save == NULL) {
                    // Error in open the file
                    printf("=+=+=+= Error saving file %s =+=+=+=\n",amrfilename_save);
                    exit(1);
                }

                while ((ch = fgetc(fd_load)) != EOF) {
                    fputc(ch, fd_save);
                }

                fclose(fd_save);
            }

            // Close the AMR format files
            fclose(fd_load);

            // write boundary conditions
            int pbc_type; int pbc_valuetype;
            ifd = fscanf(fboundary_load,"%d %d ",&pbc_type, &pbc_valuetype);
            if(samefile == 0) fprintf(fboundary_save,"%d %d ",pbc_type, pbc_valuetype);
            for(int dim=0; dim < DIM; dim++){
                int ubc_type; int ubc_valuetype;
                ifd = fscanf(fboundary_load,"%d %d ",&ubc_type, &ubc_valuetype);
                if(samefile == 0) fprintf(fboundary_save,"%d %d ",ubc_type, ubc_valuetype);
            }
            if(samefile == 0) fprintf(fboundary_save,"\n");
        }
        fclose(fboundary_load);
        if(samefile == 0) fclose(fboundary_save);

    }
}

void higflow_save_boundaries_electroosmotic(higflow_solver *ns, int myrank, int ntasks) {
    if(myrank == 0) { 
        // Loading the boundary data
        char namefile_load[1024];
        sprintf(namefile_load,"%s.bcelectroosmotic",ns->par.nameload);
        char namefile_save[1024];
        sprintf(namefile_save, "%s.bcelectroosmotic", ns->par.namesave);

        int samefile = 0;
        if(strcmp(namefile_load, namefile_save) == 0) {
            printf("=+=+=+= Load and Save names are the same - not changing the bcelectroosmotic file =+=+=+=\n");
            samefile = 1;
        }

        FILE *fboundary_load = fopen(namefile_load, "r");
        if (fboundary_load == NULL) {
            // Error in open the file
            printf("=+=+=+= Error loading file %s =+=+=+=\n",namefile_load);
            exit(1);
        }

        FILE *fboundary_save;
        if(samefile == 0) {
            fboundary_save = fopen(namefile_save, "w");
            if (fboundary_save == NULL) {
                // Error in open the file
                printf("=+=+=+= Error saving file %s =+=+=+=\n",namefile_save);
                exit(1);
            }
        }
        
        // Number of HigTrees
        int numbcs;
        int ifd = fscanf(fboundary_load,"%d\n",&numbcs);

        // write number of higtrees
        if(samefile == 0) fprintf(fboundary_save,"%d\n",numbcs);

        for(int h = 0; h < numbcs; h++) {
            // write boundary id
            int bcid;
            ifd = fscanf(fboundary_load,"%d\n",&bcid);
            if(samefile == 0) fprintf(fboundary_save,"%d\n",bcid);

            // Name of the HigTree file
            char amrfilename_load[1024];
            __higflow_readstring(amrfilename_load,1024,fboundary_load);

            // Open the AMR format file
            FILE *fd_load = fopen(amrfilename_load, "r");
            if (fd_load == NULL) {
                // Error in open the file
                printf("=+=+=+= Error loading file %s =+=+=+=\n",amrfilename_load);
                exit(1);
            }
            
            char ch;
            char amrfilename_save[1024];
            get_amrfilename_save(&amrfilename_load, &namefile_save, amrfilename_save);

            int samefile_amr = 0;
            if(strcmp(amrfilename_load, amrfilename_save) == 0) {
                printf("=+=+=+= Load and Save names are the same - not changing the boundary %d amr file =+=+=+=\n", h);
                samefile_amr = 1;
            }

            // copy the corresponding filename of the hig
            if(samefile == 0) fprintf(fboundary_save,"%s\n",amrfilename_save);

            if(samefile_amr == 0) {
                printf("=+=+=+= Saving boundary %d to %s =+=+=+=\n", h, amrfilename_save);

                // Copy
                FILE *fd_save = fopen(amrfilename_save, "w");
                if (fd_save == NULL) {
                    // Error in open the file
                    printf("=+=+=+= Error saving file %s =+=+=+=\n",amrfilename_save);
                    exit(1);
                }

                while ((ch = fgetc(fd_load)) != EOF) {
                    fputc(ch, fd_save);
                }
                fclose(fd_save);
            }

            // Close the AMR format files
            fclose(fd_load);

            // write boundary conditions
            int phibc_type; int phibc_valuetype; int psi_type; int psi_valuetype;
            int nplusbc_type; int nplusbc_valuetype; int nminusbc_type; int nminusbc_valuetype;
            ifd = fscanf(fboundary_load,"%d %d %d %d %d %d %d %d\n",&phibc_type, &phibc_valuetype,
                                                                    &psi_type, &psi_valuetype, 
                                                                    &nplusbc_type, &nplusbc_valuetype, 
                                                                    &nminusbc_type, &nminusbc_valuetype);
            if(samefile == 0) fprintf(fboundary_save, "%d %d %d %d %d %d %d %d\n", phibc_type, phibc_valuetype,
                                                                                   psi_type, psi_valuetype, 
                                                                                   nplusbc_type, nplusbc_valuetype, 
                                                                                   nminusbc_type, nminusbc_valuetype);
        }
        fclose(fboundary_load);
        if(samefile == 0) fclose(fboundary_save);

    }
}


// Loading the multiphase parameters
void higflow_load_multiphase_parameters(higflow_solver *ns, int myrank) {
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.mult.par", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL) {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%lf", &(ns->ed.mult.par.Ca));
        fclose(fd);
        if (myrank == 0) {
            printf("=+=+=+= Capillary Number: %f =+=+=+=\n", ns->ed.mult.par.Ca);
        }
    } else {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the multiphase parameters
void higflow_save_multiphase_parameters(higflow_solver *ns, int myrank) {
    if (myrank == 0) {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.mult.par", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL) {
            // Saving the parameters
            fprintf(fd, "%lf\n", (ns->ed.mult.par.Ca));
            fclose(fd);
        } else {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", ns->par.namesave);
            exit(1);
        }
    }
}

// Loading the multiphase controllers
void higflow_load_multiphase_controllers(higflow_solver *ns, int myrank) {
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.mult.contr", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL) {
        // Loading the parameters
        int ifd;
        // Phase 0
        ifd = fscanf(fd, "%d", (int *)&(ns->ed.mult.contr.flowtype0));
        ifd = fscanf(fd, "%d", &(ns->ed.mult.contr.eoflow0));
        // Phase 1
        ifd = fscanf(fd, "%d", (int *)&(ns->ed.mult.contr.flowtype1));
        ifd = fscanf(fd, "%d", &(ns->ed.mult.contr.eoflow1));
        fclose(fd);
        if (myrank == 0) {
            // Flow Type - Phase 0
            switch (ns->ed.mult.contr.flowtype0) {
                case NEWTONIAN:
                    printf("=+=+=+= Flow Type - Phase 0: Newtonian =+=+=+=\n");
                    break;
                case VISCOELASTIC:
                    printf("=+=+=+= Flow Type - Phase 0: Viscoelastic =+=+=+=\n");
                    break;
                default:
                    printf("=+=+=+= Unsupported Flow Type of Phase 0!! =+=+=+=\n");
                    exit(1);
                    break;
            }
            // Electroosmotic Flow - Phase 0
            switch (ns->ed.mult.contr.eoflow0) {
                case true:
                    printf("=+=+=+= Phase 0: Electroosmotic Flow =+=+=+=\n");
                    break;
                case false:
                    break;
                default:
                    printf("=+=+=+= Unsupported eoflow of Phase 0!! =+=+=+=\n");
                    exit(1);
                    break;
            }
            // Flow Type - Phase 1
            switch (ns->ed.mult.contr.flowtype1) {
                case NEWTONIAN:
                    printf("=+=+=+= Flow Type - Phase 1: Newtonian =+=+=+=\n");
                    break;
                case VISCOELASTIC:
                    printf("=+=+=+= Flow Type - Phase 1: Viscoelastic =+=+=+=\n");
                    break;
                default:
                    printf("=+=+=+= Unsupported Flow Type of Phase 1!! =+=+=+=\n");
                    exit(1);
                    break;
            }
            // Electroosmotic Flow - Phase 1
            switch (ns->ed.mult.contr.eoflow1) {
                case true:
                    printf("=+=+=+= Phase 1: Electroosmotic Flow =+=+=+=\n");
                    break;
                case false:
                    break;
                default:
                    printf("=+=+=+= Unsupported eoflow of Phase 1!! =+=+=+=\n");
                    exit(1);
                    break;
            }
        }
    } else {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }

    if (ns->ed.mult.contr.flowtype0 == VISCOELASTIC || ns->ed.mult.contr.flowtype1 == VISCOELASTIC) {
        ns->ed.mult.contr.flowtype_either = VISCOELASTIC;
    }
    else {
        ns->ed.mult.contr.flowtype_either = NEWTONIAN;
    }

    if (ns->ed.mult.contr.eoflow0 == true || ns->ed.mult.contr.eoflow1 == true) {
        ns->ed.mult.contr.eoflow_either = true;
    }
    else {
        ns->ed.mult.contr.eoflow_either = false;
    }
}

// Saving the multiphase controllers
void higflow_save_multiphase_controllers(higflow_solver *ns, int myrank) {
    if (myrank == 0) {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.mult.contr", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL) {
            // Saving the parameters
            fprintf(fd, "%d\n", (int)(ns->ed.mult.contr.flowtype0));
            fprintf(fd, "%d\n", (ns->ed.mult.contr.eoflow0));
            fprintf(fd, "%d\n", (int)(ns->ed.mult.contr.flowtype1));
            fprintf(fd, "%d\n", (ns->ed.mult.contr.eoflow1));
            fclose(fd);
        } else {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", ns->par.namesave);
            exit(1);
        }
    }
}

// Loading the multiphase viscoelastic parameters
void higflow_load_multiphase_viscoelastic_parameters(higflow_solver *ns, int myrank) {
    if (ns->ed.mult.contr.flowtype0 == VISCOELASTIC) {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.mult.viscpar0", ns->par.nameload);
        FILE *fd = fopen(namefile, "r");
        if (fd != NULL) {
            // Loading the parameters
            int ifd;
            // Phase 0
            ifd = fscanf(fd, "%lf", &(ns->ed.mult.ve.par0.De));
            ifd = fscanf(fd, "%lf", &(ns->ed.mult.ve.par0.beta));
            ifd = fscanf(fd, "%lf", &(ns->ed.mult.ve.par0.epsilon));
            ifd = fscanf(fd, "%lf", &(ns->ed.mult.ve.par0.psi));
            ifd = fscanf(fd, "%lf", &(ns->ed.mult.ve.par0.alpha));
            ifd = fscanf(fd, "%lf", &(ns->ed.mult.ve.par0.kernel_tol));
            ifd = fscanf(fd, "%lf", &(ns->ed.mult.ve.par0.alpha_gptt));
            ifd = fscanf(fd, "%lf", &(ns->ed.mult.ve.par0.beta_gptt));
            ifd = fscanf(fd, "%lf", &(ns->ed.mult.ve.par0.L2_fene));
            ifd = fscanf(fd, "%lf", &(ns->ed.mult.ve.par0.lambda_fene));
            ifd = fscanf(fd, "%lf", &(ns->ed.mult.ve.par0.E_fene));
            fclose(fd);
            if (myrank == 0) {
                printf("=+=+=+= Phase 0: =+=+=+=\n");
                printf("=+=+=+= Deborah Number: %f =+=+=+=\n", ns->ed.mult.ve.par0.De);
                printf("=+=+=+= Beta: %f =+=+=+=\n", ns->ed.mult.ve.par0.beta);
                printf("=+=+=+= Epsilon: %f =+=+=+=\n", ns->ed.mult.ve.par0.epsilon);
                printf("=+=+=+= Psi: %f =+=+=+=\n", ns->ed.mult.ve.par0.psi);
                printf("=+=+=+= Alpha Giesekus: %f =+=+=+=\n", ns->ed.mult.ve.par0.alpha);
                printf("=+=+=+= Kernel Tolerance: %f =+=+=+=\n", ns->ed.mult.ve.par0.kernel_tol);
                printf("=+=+=+= Alpha GPTT: %f =+=+=+=\n", ns->ed.mult.ve.par0.alpha_gptt);
                printf("=+=+=+= Beta GPTT: %f =+=+=+=\n", ns->ed.mult.ve.par0.beta_gptt);
                printf("=+=+=+= L2 FENE: %f =+=+=+=\n", ns->ed.mult.ve.par0.L2_fene);
                printf("=+=+=+= Lambda FENE: %f =+=+=+=\n", ns->ed.mult.ve.par0.lambda_fene);
                printf("=+=+=+= E FENE: %f =+=+=+=\n", ns->ed.mult.ve.par0.E_fene);
            }
        } else {
            // Error in open the file
            printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
    else{
        ns->ed.mult.ve.par0.De = 0.0;
        ns->ed.mult.ve.par0.beta = 1.0;
        ns->ed.mult.ve.par0.epsilon = 0.0;
        ns->ed.mult.ve.par0.psi = 0.0;
        ns->ed.mult.ve.par0.alpha = 0.0;
        ns->ed.mult.ve.par0.kernel_tol = 1.0;
        ns->ed.mult.ve.par0.alpha_gptt = 0.0;
        ns->ed.mult.ve.par0.beta_gptt = 0.0;
        ns->ed.mult.ve.par0.L2_fene = 0.0;
        ns->ed.mult.ve.par0.lambda_fene = 0.0;
        ns->ed.mult.ve.par0.E_fene = 0.0;
        if (myrank == 0) {
            printf("=+=+=+= Phase 0 is NEWTONIAN - considering parameters below: =+=+=+=\n");
            printf("=+=+=+= Deborah Number: %f =+=+=+=\n", ns->ed.mult.ve.par0.De);
            printf("=+=+=+= Beta: %f =+=+=+=\n", ns->ed.mult.ve.par0.beta);
            printf("=+=+=+= Epsilon: %f =+=+=+=\n", ns->ed.mult.ve.par0.epsilon);
            printf("=+=+=+= Psi: %f =+=+=+=\n", ns->ed.mult.ve.par0.psi);
            printf("=+=+=+= Alpha Giesekus: %f =+=+=+=\n", ns->ed.mult.ve.par0.alpha);
            printf("=+=+=+= Kernel Tolerance: %f =+=+=+=\n", ns->ed.mult.ve.par0.kernel_tol);
            printf("=+=+=+= Alpha GPTT: %f =+=+=+=\n", ns->ed.mult.ve.par0.alpha_gptt);
            printf("=+=+=+= Beta GPTT: %f =+=+=+=\n", ns->ed.mult.ve.par0.beta_gptt);
            printf("=+=+=+= L2 FENE: %f =+=+=+=\n", ns->ed.mult.ve.par0.L2_fene);
            printf("=+=+=+= Lambda FENE: %f =+=+=+=\n", ns->ed.mult.ve.par0.lambda_fene);
            printf("=+=+=+= E FENE: %f =+=+=+=\n", ns->ed.mult.ve.par0.E_fene);
        }
    }

    if (ns->ed.mult.contr.flowtype1 == VISCOELASTIC) {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.mult.viscpar1", ns->par.nameload);
        FILE *fd = fopen(namefile, "r");
        if (fd != NULL) {
            // Loading the parameters
            int ifd;
            // Phase 1
            ifd = fscanf(fd, "%lf", &(ns->ed.mult.ve.par1.De));
            ifd = fscanf(fd, "%lf", &(ns->ed.mult.ve.par1.beta));
            ifd = fscanf(fd, "%lf", &(ns->ed.mult.ve.par1.epsilon));
            ifd = fscanf(fd, "%lf", &(ns->ed.mult.ve.par1.psi));
            ifd = fscanf(fd, "%lf", &(ns->ed.mult.ve.par1.alpha));
            ifd = fscanf(fd, "%lf", &(ns->ed.mult.ve.par1.kernel_tol));
            ifd = fscanf(fd, "%lf", &(ns->ed.mult.ve.par1.alpha_gptt));
            ifd = fscanf(fd, "%lf", &(ns->ed.mult.ve.par1.beta_gptt));
            ifd = fscanf(fd, "%lf", &(ns->ed.mult.ve.par1.L2_fene));
            ifd = fscanf(fd, "%lf", &(ns->ed.mult.ve.par1.lambda_fene));
            ifd = fscanf(fd, "%lf", &(ns->ed.mult.ve.par1.E_fene));
            fclose(fd);
            if (myrank == 0) {
                printf("=+=+=+= Phase 1: =+=+=+=\n");
                printf("=+=+=+= Deborah Number: %f =+=+=+=\n", ns->ed.mult.ve.par1.De);
                printf("=+=+=+= Beta: %f =+=+=+=\n", ns->ed.mult.ve.par1.beta);
                printf("=+=+=+= Epsilon: %f =+=+=+=\n", ns->ed.mult.ve.par1.epsilon);
                printf("=+=+=+= Psi: %f =+=+=+=\n", ns->ed.mult.ve.par1.psi);
                printf("=+=+=+= Alpha Giesekus: %f =+=+=+=\n", ns->ed.mult.ve.par1.alpha);
                printf("=+=+=+= Kernel Tolerance: %f =+=+=+=\n", ns->ed.mult.ve.par1.kernel_tol);
                printf("=+=+=+= Alpha GPTT: %f =+=+=+=\n", ns->ed.mult.ve.par1.alpha_gptt);
                printf("=+=+=+= Beta GPTT: %f =+=+=+=\n", ns->ed.mult.ve.par1.beta_gptt);
                printf("=+=+=+= L2 FENE: %f =+=+=+=\n", ns->ed.mult.ve.par1.L2_fene);
                printf("=+=+=+= Lambda FENE: %f =+=+=+=\n", ns->ed.mult.ve.par1.lambda_fene);
                printf("=+=+=+= E FENE: %f =+=+=+=\n", ns->ed.mult.ve.par1.E_fene);
            }
        } else {
            // Error in open the file
            printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
    else{
        ns->ed.mult.ve.par1.De = 0.0;
        ns->ed.mult.ve.par1.beta = 1.0;
        ns->ed.mult.ve.par1.epsilon = 0.0;
        ns->ed.mult.ve.par1.psi = 0.0;
        ns->ed.mult.ve.par1.alpha = 0.0;
        ns->ed.mult.ve.par1.kernel_tol = 1.0;
        ns->ed.mult.ve.par1.alpha_gptt = 0.0;
        ns->ed.mult.ve.par1.beta_gptt = 0.0;
        ns->ed.mult.ve.par1.L2_fene = 0.0;
        ns->ed.mult.ve.par1.lambda_fene = 0.0;
        ns->ed.mult.ve.par1.E_fene = 0.0;
        if (myrank == 0) {
            printf("=+=+=+= Phase 1 is NEWTONIAN - considering parameters below: =+=+=+=\n");
            printf("=+=+=+= Deborah Number: %f =+=+=+=\n", ns->ed.mult.ve.par1.De);
            printf("=+=+=+= Beta: %f =+=+=+=\n", ns->ed.mult.ve.par1.beta);
            printf("=+=+=+= Epsilon: %f =+=+=+=\n", ns->ed.mult.ve.par1.epsilon);
            printf("=+=+=+= Psi: %f =+=+=+=\n", ns->ed.mult.ve.par1.psi);
            printf("=+=+=+= Alpha Giesekus: %f =+=+=+=\n", ns->ed.mult.ve.par1.alpha);
            printf("=+=+=+= Kernel Tolerance: %f =+=+=+=\n", ns->ed.mult.ve.par1.kernel_tol);
            printf("=+=+=+= Alpha GPTT: %f =+=+=+=\n", ns->ed.mult.ve.par1.alpha_gptt);
            printf("=+=+=+= Beta GPTT: %f =+=+=+=\n", ns->ed.mult.ve.par1.beta_gptt);
            printf("=+=+=+= L2 FENE: %f =+=+=+=\n", ns->ed.mult.ve.par1.L2_fene);
            printf("=+=+=+= Lambda FENE: %f =+=+=+=\n", ns->ed.mult.ve.par1.lambda_fene);
            printf("=+=+=+= E FENE: %f =+=+=+=\n", ns->ed.mult.ve.par1.E_fene);
        }
    }

    // Set general kernel tolerance
    ns->ed.mult.ve.par0.kernel_tol = min(ns->ed.mult.ve.par0.kernel_tol,ns->ed.mult.ve.par1.kernel_tol);
    ns->ed.mult.ve.par1.kernel_tol = min(ns->ed.mult.ve.par0.kernel_tol,ns->ed.mult.ve.par1.kernel_tol);
}

// Saving the multiphase viscoelastic parameters
void higflow_save_multiphase_viscoelastic_parameters(higflow_solver *ns, int myrank) {
    if (myrank == 0) {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.mult.viscpar0", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL) {
            // Saving the parameters
            fprintf(fd, "%lf\n", (ns->ed.mult.ve.par0.De));
            fprintf(fd, "%lf\n", (ns->ed.mult.ve.par0.beta));
            fprintf(fd, "%lf\n", (ns->ed.mult.ve.par0.epsilon));
            fprintf(fd, "%lf\n", (ns->ed.mult.ve.par0.psi));
            fprintf(fd, "%lf\n", (ns->ed.mult.ve.par0.alpha));
            fprintf(fd, "%lf\n", (ns->ed.mult.ve.par0.kernel_tol));
            fprintf(fd, "%lf\n", (ns->ed.mult.ve.par0.alpha_gptt));
            fprintf(fd, "%lf\n", (ns->ed.mult.ve.par0.beta_gptt));
            fprintf(fd, "%lf\n", (ns->ed.mult.ve.par0.L2_fene));
            fprintf(fd, "%lf\n", (ns->ed.mult.ve.par0.lambda_fene));
            fprintf(fd, "%lf\n", (ns->ed.mult.ve.par0.E_fene));
            fclose(fd);
        } else {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", ns->par.namesave);
            exit(1);
        }
        // Parameters file name
        snprintf(namefile, sizeof namefile, "%s.mult.viscpar1", ns->par.namesave);
        fd = fopen(namefile, "w");
        if (fd != NULL) {
            // Saving the parameters
            fprintf(fd, "%lf\n", (ns->ed.mult.ve.par1.De));
            fprintf(fd, "%lf\n", (ns->ed.mult.ve.par1.beta));
            fprintf(fd, "%lf\n", (ns->ed.mult.ve.par1.epsilon));
            fprintf(fd, "%lf\n", (ns->ed.mult.ve.par1.psi));
            fprintf(fd, "%lf\n", (ns->ed.mult.ve.par1.alpha));
            fprintf(fd, "%lf\n", (ns->ed.mult.ve.par1.kernel_tol));
            fprintf(fd, "%lf\n", (ns->ed.mult.ve.par1.alpha_gptt));
            fprintf(fd, "%lf\n", (ns->ed.mult.ve.par1.beta_gptt));
            fprintf(fd, "%lf\n", (ns->ed.mult.ve.par1.L2_fene));
            fprintf(fd, "%lf\n", (ns->ed.mult.ve.par1.lambda_fene));
            fprintf(fd, "%lf\n", (ns->ed.mult.ve.par1.E_fene));
            fclose(fd);
        } else {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", ns->par.namesave);
            exit(1);
        }
    }
}

// Loading the multiphase viscoelastic controllers
void higflow_load_multiphase_viscoelastic_controllers(higflow_solver *ns, int myrank) {
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.mult.visccontr", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL) {
        // Loading the parameters
        int ifd;
        // Phase 0
        ifd = fscanf(fd, "%d", (int *)&(ns->ed.mult.ve.contr.model0));
        // Phase 1
        ifd = fscanf(fd, "%d", (int *)&(ns->ed.mult.ve.contr.model1));
        // Both Phases
        ifd = fscanf(fd, "%d", (int *)&(ns->ed.mult.ve.contr.discrtype));
        ifd = fscanf(fd, "%d", (int *)&(ns->ed.mult.ve.contr.convecdiscrtype));
        fclose(fd);
        if (myrank == 0) {
            //Constitutive Equation Model - Phase 0
            switch (ns->ed.mult.ve.contr.model0) {
                case USERSET:
                    printf("=+=+=+= Constitutive Equation Model - Phase 0: User Defined =+=+=+=\n");
                    break;
                case OLDROYD_B:
                    printf("=+=+=+= Constitutive Equation Model - Phase 0: Oldroyd-B =+=+=+=\n");
                    break;
                case GIESEKUS:
                    printf("=+=+=+= Constitutive Equation Model - Phase 0: Giesekus =+=+=+=\n");
                    break;
                case LPTT:
                    printf("=+=+=+= Constitutive Equation Model - Phase 0: LPTT =+=+=+=\n");
                    break;
                case GPTT:
                    printf("=+=+=+= Constitutive Equation Model - Phase 0: GPTT =+=+=+=\n");
                    break;
                case FENE_P:
                    printf("=+=+=+= Constitutive Equation Model - Phase 0: FENE-P =+=+=+=\n");
                    break;
                case E_FENE:
                    printf("=+=+=+= Constitutive Equation Model - Phase 0: e-FENE =+=+=+=\n");
                    break;
            }
            //Constitutive Equation Model - Phase 1
            switch (ns->ed.mult.ve.contr.model1) {
                case USERSET:
                    printf("=+=+=+= Constitutive Equation Model - Phase 1: User Defined =+=+=+=\n");
                    break;
                case OLDROYD_B:
                    printf("=+=+=+= Constitutive Equation Model - Phase 1: Oldroyd-B =+=+=+=\n");
                    break;
                case GIESEKUS:
                    printf("=+=+=+= Constitutive Equation Model - Phase 1: Giesekus =+=+=+=\n");
                    break;
                case LPTT:
                    printf("=+=+=+= Constitutive Equation Model - Phase 1: LPTT =+=+=+=\n");
                    break;
                case GPTT:
                    printf("=+=+=+= Constitutive Equation Model - Phase 1 GPTT =+=+=+=\n");
                    break;
                case FENE_P:
                    printf("=+=+=+= Constitutive Equation Model - Phase 1 FENE-P =+=+=+=\n");
                    break;
                case E_FENE:
                    printf("=+=+=+= Constitutive Equation Model - Phase 1 e-FENE =+=+=+=\n");
                    break;
            }
            // Constitutive Equation Discretization
            switch (ns->ed.mult.ve.contr.discrtype) {
                case EXPLICIT:
                    printf("=+=+=+= Constitutive Equation Discretization : Explicit =+=+=+=\n");
                    break;
                case IMPLICIT:
                    printf("=+=+=+= Constitutive Equation Discretization : Implicit =+=+=+=\n");
                    break;
            }
            // Constitutive Equation Convective Term
            switch (ns->ed.mult.ve.contr.convecdiscrtype) {
                case CELL_CENTRAL:
                    printf("=+=+=+= Constitutive Equation Convective Term : Central  =+=+=+=\n");
                    break;
                case CELL_CUBISTA:
                    printf("=+=+=+= Constitutive Equation Convective Term : CUBISTA =+=+=+=\n");
                    break;
            }
          
        }
    } else {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the multiphase viscoelastic controllers
void higflow_save_multiphase_viscoelastic_controllers(higflow_solver *ns, int myrank) {
    if (myrank == 0) {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.mult.visccontr", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL) {
            // Saving the parameters
            fprintf(fd, "%d\n", (int)(ns->ed.mult.ve.contr.model0));
            fprintf(fd, "%d\n", (int)(ns->ed.mult.ve.contr.model1));
            fprintf(fd, "%d\n", (int)(ns->ed.mult.ve.contr.discrtype));
            fprintf(fd, "%d\n", (int)(ns->ed.mult.ve.contr.convecdiscrtype));
            fclose(fd);
        } else {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", ns->par.namesave);
            exit(1);
        }
    }
}


// Loading the viscoelastic parameters
void higflow_load_viscoelastic_parameters(higflow_solver *ns, int myrank) {
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.viscpar", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL) {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%lf", &(ns->ed.ve.par.De));
        ifd = fscanf(fd, "%lf", &(ns->ed.ve.par.beta));
        ifd = fscanf(fd, "%lf", &(ns->ed.ve.par.epsilon));
        ifd = fscanf(fd, "%lf", &(ns->ed.ve.par.psi));
        ifd = fscanf(fd, "%lf", &(ns->ed.ve.par.alpha));
        ifd = fscanf(fd, "%lf", &(ns->ed.ve.par.kernel_tol));
        ifd = fscanf(fd,"%lf",&(ns->ed.ve.par.alpha_gptt));
        ifd = fscanf(fd,"%lf",&(ns->ed.ve.par.beta_gptt));
        ifd = fscanf(fd,"%lf",&(ns->ed.ve.par.L2_fene));
        ifd = fscanf(fd,"%lf",&(ns->ed.ve.par.lambda_fene));
        ifd = fscanf(fd,"%lf",&(ns->ed.ve.par.E_fene));
        fclose(fd);
        if (myrank == 0) {
            printf("=+=+=+= Deborah Number: %f =+=+=+=\n", ns->ed.ve.par.De);
            printf("=+=+=+= Beta: %f =+=+=+=\n", ns->ed.ve.par.beta);
        }
    } else {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}



// Saving the viscoelastic parameters
void higflow_save_viscoelastic_parameters(higflow_solver *ns, int myrank) {
    if (myrank == 0) {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.viscpar", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL) {
            // Saving the parameters
            fprintf(fd, "%lf\n", (ns->ed.ve.par.De));
            fprintf(fd, "%lf\n", (ns->ed.ve.par.beta));
            fprintf(fd, "%lf\n", (ns->ed.ve.par.epsilon));
            fprintf(fd, "%lf\n", (ns->ed.ve.par.psi));
            fprintf(fd, "%lf\n", (ns->ed.ve.par.alpha));
            fprintf(fd, "%lf\n", (ns->ed.ve.par.kernel_tol));
            fprintf(fd,"%lf\n",(ns->ed.ve.par.alpha_gptt));
            fprintf(fd,"%lf\n",(ns->ed.ve.par.beta_gptt));
            fprintf(fd,"%lf\n",(ns->ed.ve.par.L2_fene));
            fprintf(fd,"%lf\n",(ns->ed.ve.par.lambda_fene));
            fprintf(fd,"%lf\n",(ns->ed.ve.par.E_fene));
            fclose(fd);
        } else {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}

// Loading the viscoelastic controllers
void higflow_load_viscoelastic_controllers(higflow_solver *ns, int myrank) {
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.visccontr", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL) {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%d", (int *)&(ns->ed.ve.contr.model));
        ifd = fscanf(fd, "%d", (int *)&(ns->ed.ve.contr.discrtype));
        ifd = fscanf(fd, "%d", (int *)&(ns->ed.ve.contr.convecdiscrtype));
        fclose(fd);
        if (myrank == 0) {
            switch (ns->ed.ve.contr.model) {
                case USERSET:
                    printf("=+=+=+= Constitutive Equation Model: User Defined =+=+=+=\n");
                    break;
                case OLDROYD_B:
                    printf("=+=+=+= Constitutive Equation Model: Oldroyd-B =+=+=+=\n");
                    break;
                case GIESEKUS:
                    printf("=+=+=+= Constitutive Equation Model: Giesekus =+=+=+=\n");
                    break;
                case LPTT:
                    printf("=+=+=+= Constitutive Equation Model: LPTT =+=+=+=\n");
                    break;
                case GPTT:
                    printf("=+=+=+= Constitutive Equation Model: GPTT =+=+=+=\n");
                    break;
                case FENE_P:
                    printf("=+=+=+= Constitutive Equation Model: FENE-P =+=+=+=\n");
                    break;
                case E_FENE:
                    printf("=+=+=+= Constitutive Equation Model: e-FENE =+=+=+=\n");
                    break;
            }
            switch (ns->ed.ve.contr.discrtype) {
                case EXPLICIT:
                    printf("=+=+=+= Constitutive Equation Discretization: Explicit =+=+=+=\n");
                    break;
                case IMPLICIT:
                    printf("=+=+=+= Constitutive Equation Discretization: Implicit =+=+=+=\n");
                    break;
            }
            switch (ns->ed.ve.contr.convecdiscrtype) {
                case CELL_CENTRAL:
                    printf("=+=+=+= Constitutive Equation Convective Term: Central  =+=+=+=\n");
                    break;
                case CELL_CUBISTA:
                    printf("=+=+=+= Constitutive Equation Convective Term: CUBISTA =+=+=+=\n");
                    break;
            }
        }
    } else {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the viscoelastic controllers
void higflow_save_viscoelastic_controllers(higflow_solver *ns, int myrank) {
    if (myrank == 0) {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.visccontr", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL) {
            // Saving the parameters
            fprintf(fd, "%d\n", (int)(ns->ed.ve.contr.model));
            fprintf(fd, "%d\n", (int)(ns->ed.ve.contr.discrtype));
            fprintf(fd, "%d\n", (int)(ns->ed.ve.contr.convecdiscrtype));
            fclose(fd);
        } else {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}

// Loading the electroosmotic parameters
void higflow_load_electroosmotic_parameters(higflow_solver *ns, int myrank) {
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.eopar", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL) {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%lf", &(ns->ed.eo.par.alpha));
        ifd = fscanf(fd, "%lf", &(ns->ed.eo.par.delta));
        ifd = fscanf(fd, "%lf", &(ns->ed.eo.par.Pe));
        ifd = fscanf(fd, "%lf", &(ns->ed.eo.par.Ex));
        fclose(fd);
        if (myrank == 0) {
            printf("=+=+=+= alpha: %e =+=+=+=\n", ns->ed.eo.par.alpha);
            printf("=+=+=+= delta: %e =+=+=+=\n", ns->ed.eo.par.delta);
            printf("=+=+=+= Peclet number: %f =+=+=+=\n", ns->ed.eo.par.Pe);
            printf("=+=+=+= Externel field Ex : %e =+=+=+=\n", ns->ed.eo.par.Ex);
        }
    } else {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the electroosmotic parameters
void higflow_save_electroosmotic_parameters(higflow_solver *ns, int myrank) {
    if (myrank == 0) {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.eopar", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL) {
            // Saving the parameters
            fprintf(fd, "%lf\n", (ns->ed.eo.par.alpha));
            fprintf(fd, "%lf\n", (ns->ed.eo.par.delta));
            fprintf(fd, "%lf\n", (ns->ed.eo.par.Pe));
            fprintf(fd, "%lf\n", (ns->ed.eo.par.Ex));
            fclose(fd);
        } else {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}

// Loading the electroosmotic controllers
void higflow_load_electroosmotic_controllers(higflow_solver *ns, int myrank) {
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.eocontr", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL) {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%d", (int *)&(ns->ed.eo.contr.eo_model));
        ifd = fscanf(fd, "%d", (int *)&(ns->ed.eo.contr.tempdiscrtype));
        ifd = fscanf(fd, "%d", (int *)&(ns->ed.eo.contr.convecdiscrtype));
        // ifd = fscanf(fd, "%d", (int *)&(ns->ed.eo.contr.electric_convecdiscrtype));
        // ifd = fscanf(fd, "%d", (int *)&(ns->ed.eo.contr.electricdiv_discrtype));
        fclose(fd);
        if (myrank == 0) {
            switch (ns->ed.eo.contr.eo_model) {
                case PNP:
                    printf("=+=+=+= Electroosmotic Model: PNP =+=+=+=\n");
                    break;
                case PB:
                    printf("=+=+=+= Electroosmotic Model: PB =+=+=+=\n");
                    break;
                case PBDH:
                    printf("=+=+=+= Electroosmotic Model: PBDH =+=+=+=\n");
                    break;
                case PBDH_ANALYTIC:
                    printf("=+=+=+= Electroosmotic Model: analytic PBDH =+=+=+=\n");
                    break;
            }
            switch (ns->ed.eo.contr.tempdiscrtype) {
                case EXPLICIT_EULER:
                    printf("=+=+=+= Ionic concentration temporal discretization: Euler Explicit =+=+=+=\n");
                    break;
                case SEMI_IMPLICIT_EULER:
                    printf("=+=+=+= Ionic concentration temporal discretization: Euler Semi-Implicit =+=+=+=\n");
                    break;
                default:
                    printf("=+=+=+= Unsupported Temporal Discretization!! =+=+=+=\n");
                    exit(1);
                    break;
            }
            switch (ns->ed.eo.contr.convecdiscrtype) {
                case CELL_CENTRAL:
                    printf("=+=+=+= Ionic concentration convective term: Central =+=+=+=\n");
                    printf("=+=+=+= Electric Field terms will be source          =+=+=+=\n");
                    break;
                case CELL_CUBISTA:
                    printf("=+=+=+= Ionic concentration convective term: Cubista =+=+=+=\n");
                    break;
            }
        }
    } else {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the electroosmotic controllers
void higflow_save_electroosmotic_controllers(higflow_solver *ns, int myrank) {
    if (myrank == 0) {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.eocontr", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL) {
            // Saving the parameters
            fprintf(fd, "%d\n", (int)(ns->ed.eo.contr.eo_model));
            fprintf(fd, "%d\n", (int)(ns->ed.eo.contr.tempdiscrtype));
            fprintf(fd, "%d\n", (int)(ns->ed.eo.contr.convecdiscrtype));
            // fprintf(fd, "%d\n", (ns->ed.eo.contr.electric_convecdiscrtype));
            // fprintf(fd, "%d\n", (ns->ed.eo.contr.electricdiv_discrtype));
            fclose(fd);
        } else {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}

// Loading the viscoelastic parameters
void higflow_load_viscoelastic_integral_parameters(higflow_solver *ns, int myrank) {
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.viscintpar",ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL) {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd,"%lf",&(ns->ed.im.par.De));
        ifd = fscanf(fd,"%lf",&(ns->ed.im.par.alpha)); //for damping function
        ifd = fscanf(fd,"%lf",&(ns->ed.im.par.beta));  //for damping function
        ifd = fscanf(fd,"%lf",&(ns->ed.im.par.scorte));
        ifd = fscanf(fd,"%lf",&(ns->ed.im.par.rho));
        ifd = fscanf(fd,"%lf",&(ns->ed.im.par.v_ref));
        ifd = fscanf(fd,"%lf",&(ns->ed.im.par.l_ref));
       if (ns->ed.im.contr.model == KBKZ) {
            ifd = fscanf(fd,"%d",&(ns->ed.im.par.M));
            for (int i = 0; i < 8; i++) {
                ifd = fscanf(fd,"%lf",&(ns->ed.im.par.a[i]));
            }
            for (int i = 0; i < 8; i++) {
                ifd = fscanf(fd,"%lf",&(ns->ed.im.par.lambda[i]));
            }
   } else if (ns->ed.im.contr.model == KBKZ_FRACTIONAL) {
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.alpha_frac));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.beta_frac));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.Phi1));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.Phi2));
        }
        fclose(fd);
        if (myrank == 0) {
       printf("=+=+=+= Deborah Number: %f =+=+=+=\n",ns->ed.im.par.De);
       printf("=+=+=+= Alpha: %f =+=+=+=\n",ns->ed.im.par.alpha);
       printf("=+=+=+= Beta: %f =+=+=+=\n",ns->ed.im.par.beta);
       if (ns->ed.im.contr.model == KBKZ) {
           printf("=+=+=+= M: %d =+=+=+=\n",ns->ed.im.par.M);
                for (int i = 0; i < 8; i++) {
               printf("=+=+=+= a(%d): %f =+=+=+=\n",i,ns->ed.im.par.a[i]);
                }
                for (int i = 0; i < 8; i++) {
               printf("=+=+=+= Lambda(%d): %f =+=+=+=\n",i,ns->ed.im.par.lambda[i]);
                }
       } else if (ns->ed.im.contr.model == KBKZ_FRACTIONAL) {
           printf("=+=+=+= Alpha Fractional: %f =+=+=+=\n",ns->ed.im.par.alpha_frac);
           printf("=+=+=+= Beta Fractional: %f =+=+=+=\n",ns->ed.im.par.beta_frac);
           printf("=+=+=+= Phi1: %f =+=+=+=\n",ns->ed.im.par.Phi1);
           printf("=+=+=+= Phi2: %f =+=+=+=\n",ns->ed.im.par.Phi2);
            }
        }
    } else {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n",namefile);
        exit(1);
    }
}

// Saving the viscoelastic integral parameters
void higflow_save_viscoelastic_integral_parameters(higflow_solver *ns, int myrank) {
    if (myrank == 0) {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.viscintpar", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL) {
            // Saving the parameters
            fprintf(fd, "%lf\n", (ns->ed.im.par.De));
            fprintf(fd, "%lf\n", (ns->ed.im.par.alpha)); // for damping function
            fprintf(fd, "%lf\n", (ns->ed.im.par.beta));  // for damping function
            fprintf(fd,"%lf\n",(ns->ed.im.par.scorte));
            fprintf(fd,"%lf\n",(ns->ed.im.par.rho));
            fprintf(fd,"%lf\n",(ns->ed.im.par.v_ref));
            fprintf(fd,"%lf\n",(ns->ed.im.par.l_ref));
            if (ns->ed.im.contr.model == KBKZ) {
                fprintf(fd, "%d\n", (ns->ed.im.par.M));
                for (int k = 0; k < 8; k++) {
                    fprintf(fd, "%lf\n", (ns->ed.im.par.a[k]));
                }
                for (int k = 0; k < 8; k++) {
                    fprintf(fd, "%lf\n", (ns->ed.im.par.lambda[k]));
                }
            } else if (ns->ed.im.contr.model == KBKZ_FRACTIONAL) {
                fprintf(fd, "%lf\n", (ns->ed.im.par.alpha_frac));
                fprintf(fd, "%lf\n", (ns->ed.im.par.beta_frac));
                fprintf(fd, "%lf\n", (ns->ed.im.par.Phi1));
                fprintf(fd, "%lf\n", (ns->ed.im.par.Phi2));
            }
            fclose(fd);
        } else {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}

// Loading the viscoelastic integral controllers
void higflow_load_viscoelastic_integral_controllers(higflow_solver *ns, int myrank) {
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.viscintcontr",ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL) {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd,"%d",(int *)&(ns->ed.im.contr.model));
        ifd = fscanf(fd,"%d",(int *)&(ns->ed.im.contr.model_H));
        ifd = fscanf(fd,"%d",(int *)&(ns->ed.im.contr.discrtype));
        ifd = fscanf(fd,"%d",(int *)&(ns->ed.im.contr.convecdiscrtype));
        fclose(fd);
        if (myrank == 0) {
            switch (ns->ed.im.contr.model) {
                case KBKZ:
                    printf("=+=+=+= Constitutive Equation Model: KBKZ-PSM =+=+=+=\n");
                    break;
                case KBKZ_FRACTIONAL:
                    printf("=+=+=+= Constitutive Equation Model: Fractional =+=+=+=\n");
                    break;
            }
            switch (ns->ed.im.contr.model_H) {
                case PSM:
                  printf("=+=+=+= Constitutive Equation Model: PSM =+=+=+=\n");
                break;
                case UCM:
                  printf("=+=+=+= Constitutive Equation Model: UCM =+=+=+=\n");
                break;
            }
            switch (ns->ed.im.contr.discrtype) {
                case EXPLICIT:
                    printf("=+=+=+= Constitutive Equation Discretization: Explicit =+=+=+=\n");
                    break;
                case IMPLICIT:
                    printf("=+=+=+= Constitutive Equation Discretization: Implicit - not implemented =+=+=+=\n");
                    exit(1);
                    break;
            }
            switch (ns->ed.im.contr.convecdiscrtype) {
                case CELL_CENTRAL:
                    printf("=+=+=+= Constitutive Equation Convective Term: Central  =+=+=+=\n");
                    break;
                case CELL_CUBISTA:
                    printf("=+=+=+= Constitutive Equation Convective Term: CUBISTA =+=+=+=\n");
                    break;
            }
        }
    } else {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the viscoelastic integral controllers
void higflow_save_viscoelastic_integral_controllers(higflow_solver *ns, int myrank) {
    if (myrank == 0) {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.viscintcontr",ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL) {
            // Saving the parameters
            fprintf(fd,"%d\n",(int)(ns->ed.im.contr.model));
            fprintf(fd,"%d\n",(int)(ns->ed.im.contr.model_H));
            fprintf(fd,"%d\n",(int)(ns->ed.im.contr.discrtype));
            fprintf(fd,"%d\n",(int)(ns->ed.im.contr.convecdiscrtype));
            fclose(fd);
        } else {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n",namefile);
            exit(1);
        }
    }
}
