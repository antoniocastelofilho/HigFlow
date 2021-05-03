// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver IO - version 10/11/2016
// *******************************************************************
// *******************************************************************

#include "hig-flow-io.h"

// *******************************************************************
// Navier-Stokes Print for Visualize
// *******************************************************************

// Print the electro-osmotic proporties
void higflow_print_electroosmotic_properties(higflow_solver *ns, FILE *data, int dimprint, real pprint) {
    if (ns->contr.modelflowtype == 1) {
        // Get the cosntants
        real epsprint = 1.0e-8;
        real deltaeo = ns->ed.eo.par.delta;
        // Necessary properties
        real phi, psi, np, nm, rhoe;
        // Get the local sub-domain for the cells
        sim_domain *sdphi = psd_get_local_domain(ns->ed.eo.psdEOphi);
        sim_domain *sdpsi = psd_get_local_domain(ns->ed.eo.psdEOpsi);
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
            int clid = mp_lookup(mp, hig_get_cid(c));
            // Get the inside/outside inflow point cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Get the electro-osmotic properties 
            phi = compute_value_at_point(ns->ed.eo.sdEOphi, ccenter, ccenter, 1.0, ns->ed.eo.dpphi, ns->ed.stn);
            psi = compute_value_at_point(ns->ed.eo.sdEOpsi, ccenter, ccenter, 1.0, ns->ed.eo.dppsi, ns->ed.stn);
            np = compute_value_at_point(ns->ed.eo.sdEOnplus, ccenter, ccenter, 1.0, ns->ed.eo.dpnplus, ns->ed.stn);
            nm = compute_value_at_point(ns->ed.eo.sdEOnminus, ccenter, ccenter, 1.0, ns->ed.eo.dpnminus, ns->ed.stn);
            rhoe = deltaeo * (np - nm);
            //Print data file
            if ((ccenter[dimprint] < pprint + epsprint) && (ccenter[dimprint] > pprint - epsprint)) {
                fprintf(data, "%15.10lf  %15.10lf  %15.12lf  %15.12lf  %15.12lf  %15.12lf  %15.12lf\n", ccenter[0],
                        ccenter[1], phi, psi, np, nm, rhoe);
            }
        }
        // Destroy the iterator
        higcit_destroy(it);
    }
}

// Print the Polymeric Tensor
void higflow_print_polymeric_tensor(higflow_solver *ns, FILE *data, int dimprint, real pprint) {
    if (ns->contr.flowtype == 3) {
        // Get the cosntants
        real Re = ns->par.Re;
        real De = ns->ed.ve.par.De;
        real beta = ns->ed.ve.par.beta;
        real tol = ns->ed.ve.par.kernel_tol;
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
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpD[i][j],
                                                      ns->ed.stn);
                    D[i][j] = 0.5 * (Du[i][j] + Du[j][i]);
                    // Get S tensor
                    S[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpS[i][j],
                                                     ns->ed.stn);
                }
            }
            //Print polymeric stress data file
            if (ccenter[dimprint] == pprint) {
                fprintf(data, "%lf  %lf  %15.12lf  %15.12lf  %15.12lf  %15.12lf\n", ccenter[0], ccenter[1],
                        S[0][0] + 2.0 * (1 - beta) * D[0][0] / Re, S[0][1] + 2.0 * (1 - beta) * D[0][1] / Re,
                        S[1][0] + 2.0 * (1 - beta) * D[1][0] / Re, S[1][1] + 2.0 * (1 - beta) * D[1][1] / Re);
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
            if (fcenter[dimprint] == pprint) {
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
            switch (ns->contr.flowtype) {
                case 0:
                    // Newtonian
                    higflow_print_vtk3D(ns, rank);
                    break;
                case 1:
                    // Generalized Newtonian
                    higflow_print_vtk3D(ns, rank);
                    break;
                case 2:
                    // Multiphase
                    higflow_print_vtk3D(ns, rank);
                    break;
                case 3:
                    // Viscoelastic
                    higflow_print_vtk3D_viscoelastic(ns, rank);
                    break;
            }
            break;
    }
}

// Print the VTK file for visualize 2D
void higflow_print_vtk2D(higflow_solver *ns, int rank) {
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
    fprintf(f, "POINTS %ld float\n", 4 * numleafs);
    hig_cell *c;
    for (; !higcit_isfinished(it_t); higcit_nextcell(it_t)) {
        c = higcit_getcell(it_t);
        fprintf(f, "%lf %lf 0\n", c->lowpoint[0], c->lowpoint[1]);
        fprintf(f, "%lf %lf 0\n", c->highpoint[0], c->lowpoint[1]);
        fprintf(f, "%lf %lf 0\n", c->highpoint[0], c->highpoint[1]);
        fprintf(f, "%lf %lf 0\n", c->lowpoint[0], c->highpoint[1]);
    }
    fprintf(f, "CELLS %ld %ld\n", numleafs, 5 * numleafs);
    for (int i = 0; i < numleafs; i++) {
        fprintf(f, "%d %d %d %d %d\n", 4, 4 * i, 4 * i + 1, 4 * i + 2, 4 * i + 3);
    }
    fprintf(f, "CELL_TYPES %ld\n", numleafs);
    for (int i = 0; i < numleafs; i++) {
        fprintf(f, "7 "); //vtk cell types 7 = VTK_POLYGON // 9
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

    int degree = 2;
    const int maxpts = 2 * DIM * wls_num_min_points(DIM, degree);
    //const int maxpts = 2*DIM*wls_num_min_points(degree);
    Point pts[maxpts + 1];
    real dist[maxpts + 1];
    real w0[maxpts], w1[maxpts], w2[maxpts], w3[maxpts], w4[maxpts], w5[maxpts], w6[maxpts], w7[maxpts];
    uniqueid gids[maxpts];
	 // Saving vector properties of cell faces
    fprintf(f, "\nPOINT_DATA %ld\nVECTORS vel FLOAT\n", 4 * numleafs);
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        Point cdelta, ccenter, clowpoint, chightpoint;
        hig_get_delta(c, cdelta);
        hig_get_center(c, ccenter);
        Point p0, p1, p2, p3;
        p0[0] = c->lowpoint[0];
        p0[1] = c->lowpoint[1];
        p1[0] = c->highpoint[0];
        p1[1] = c->lowpoint[1];
        p2[0] = c->highpoint[0];
        p2[1] = c->highpoint[1];
        p3[0] = c->lowpoint[0];
        p3[1] = c->highpoint[1];

        uniqueid id = hig_get_cid(c);
        int cgid = mp_lookup(m, id);
        double lu0 = 0, lv0 = 0, lw0 = 0;
        double lu1 = 0, lv1 = 0, lw1 = 0;
        double lu2 = 0, lv2 = 0, lw2 = 0;
        double lu3 = 0, lv3 = 0, lw3 = 0;

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
            }
            real value1 = compute_facet_value_at_point(sfdu2[dim], ccenter, p1, 1.0, ns->dpu[dim], ns->stn);
            switch (dim) {
                case 0:
                    lu1 = value1;
                    break;
                case 1:
                    lv1 = value1;
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
            }
            real value3 = compute_facet_value_at_point(sfdu2[dim], ccenter, p3, 1.0, ns->dpu[dim], ns->stn);
            switch (dim) {
                case 0:
                    lu3 = value3;
                    break;
                case 1:
                    lv3 = value3;
                    break;
            }
        }
        fprintf(f, "%e %e %e\n", lu0, lv0, lw0);
        fprintf(f, "%e %e %e\n", lu1, lv1, lw1);
        fprintf(f, "%e %e %e\n", lu2, lv2, lw2);
        fprintf(f, "%e %e %e\n", lu3, lv3, lw3);
    }
    higcit_destroy(it);
    
    // Saving pressure from the center of the cell 
    fprintf(f, "\nCELL_DATA %ld\nSCALARS p FLOAT\nLOOKUP_TABLE default\n", numcells);
    for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        uniqueid id = hig_get_cid(c);
        int cgid = mp_lookup(m, id);
        
        Point cdelta, ccenter, clowpoint, chightpoint;
        hig_get_delta(c,cdelta);
        hig_get_center(c,ccenter);
        
        real val  = dp_get_value(ns->dpp, cgid);
        fprintf(f, "%f\n", val);
    }
    higcit_destroy(it);
    
    // Saving scalar properties from the center of the cell 
    switch (ns->contr.flowtype) {
        case 0:
        break;
        
        case 1:
        break;
        
        case 2:
            fprintf(f,"\nSCALARS FracVol FLOAT\nLOOKUP_TABLE default\n");
            sfdu2[0] = psfd_get_local_domain(ns->psfdu[0]);
            for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
               hig_cell *c = higcit_getcell(it);
               uniqueid id = hig_get_cid(c);
               int cgid = mp_lookup(m, id);
               Point cdelta, ccenter, p0, p1;
               hig_get_delta(c, cdelta);
               hig_get_center(c, ccenter);
               real value = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
               fprintf(f, "%f\n", value);
            }
            higcit_destroy(it);
            
            break;
            
            case 3:
            break;
    }
    fclose(f);
}

// Print the VTK file for visualize 3D
void higflow_print_vtk3D(higflow_solver *ns, int rank) {
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
        fprintf(f, "%d %d %d %d %d %d %d %d %d\n", 8, 8 * i + 0, 8 * i + 1, 8 * i + 2, 8 * i + 3, 8 * i + 4, 8 * i + 5,
                8 * i + 6, 8 * i + 7); //x == 0
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
        int cgid = mp_lookup(m, id);
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
    fprintf(f, "\nCELL_DATA %ld\nSCALARS p FLOAT\nLOOKUP_TABLE default\n", numcells);
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        uniqueid id = hig_get_cid(c);
        int cgid = mp_lookup(m, id);
        real val = dp_get_value(ns->dpp, cgid);
        fprintf(f, "%e\n", val);
        // fprintf(f, "%f\n", val);
    }
    higcit_destroy(it);
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
        fprintf(f, "%d %d %d %d %d %d %d %d %d\n", 8, 8 * i + 0, 8 * i + 1, 8 * i + 2, 8 * i + 3, 8 * i + 4, 8 * i + 5,
                8 * i + 6, 8 * i + 7); //x == 0
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
        int cgid = mp_lookup(m, id);
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
        int cgid = mp_lookup(m, id);
        real val = dp_get_value(ns->dpp, cgid);
        fprintf(f, "%e\n", val);
        // fprintf(f, "%f\n", val);
    }
    higcit_destroy(it);
*/

    fprintf(f, "\nCELL_DATA %ld\nSCALARS S00 FLOAT\nLOOKUP_TABLE default\n", numcells);
    for (it = sd_get_domain_celliterator(ns->ed.sdED); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        uniqueid id = hig_get_cid(c);
        int cgid = mp_lookup(m, id);
        real val = dp_get_value(ns->ed.ve.dpS[0][0], cgid);
        fprintf(f, "%e\n", val);
    }
    higcit_destroy(it);
    fclose(f);
}

// *******************************************************************
// Navier-Stokes Load and Save
// *******************************************************************

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
        ifd = fscanf(fd, "%d", &(ns->contr.projtype));
        ifd = fscanf(fd, "%d", &(ns->contr.flowtype));
        ifd = fscanf(fd, "%d", &(ns->contr.modelflowtype));
        ifd = fscanf(fd, "%d", &(ns->contr.tempdiscrtype));
        ifd = fscanf(fd, "%d", &(ns->contr.spatialdiscrtype));
        ifd = fscanf(fd, "%d", &(ns->contr.convecdiscrtype));
        ifd = fscanf(fd, "%d", &(ns->contr.secondconvecdiscrtype));
        fclose(fd);

        // Printing the controllers
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (myrank == 0) {
            printf("=+=+=+= Controllers =+=+=+=\n");
            switch (ns->contr.projtype) {
                case 0:
                    printf("=+=+=+= Projection Method: Non Incremental =+=+=+=\n");
                    break;
                case 1:
                    printf("=+=+=+= Projection Method: Incremental =+=+=+=\n");
                    break;
            }
            switch (ns->contr.flowtype) {
                case 0:
                    printf("=+=+=+= Flow Model: Newtonian =+=+=+=\n");
                    break;
                case 1:
                    printf("=+= Flow Model: Generalized Newtonian =+=\n");
                    break;
                case 2:
                    printf("=+=+=+= Flow Model: Multifase =+=+=+=\n");
                    break;
                case 3:
                    printf("=+=+=+= Flow Model: Viscoelastic =+=+=+=\n");
                    break;
            }
            switch (ns->contr.modelflowtype) {
                case 1:
                    printf("=+= Flow Model Type: Eletroosmotic =+=\n");
                    break;
            }
            switch (ns->contr.tempdiscrtype) {
                case 0:
                    printf("=+=+=+= Temporal Discretization: Explicit Euler =+=+=+=\n");
                    break;
                case 1:
                    printf("=+=+=+= Temporal Discretization: Runge-Kutta 2 =+=+=+=\n");
                    break;
                case 2:
                    printf("=+=+=+= Temporal Discretization: Runge-Kutta 3 =+=+=+=\n");
                    break;
                case 3:
                    printf("=+=+=+= Temporal Discretization: Semi-Implicit Euler =+=+=+=\n");
                    break;
                case 4:
                    printf("=+=+=+= Temporal Discretization: Semi-Implicit Crank-Nicolson =+=+=+=\n");
                    break;
                case 5:
                    printf("=+=+=+= Temporal Discretization: Semi-Implicit BDF2 =+=+=+=\n");
                    break;
            }
            switch (ns->contr.spatialdiscrtype) {
                case 0:
                    printf("=+=+=+= Spatial Discretization: Second Order =+=+=+=\n");
                    break;
                case 1:
                    printf("=+=+=+= Spatial Discretization: Forth Order =+=+=+=\n");
                    break;
            }
            switch (ns->contr.convecdiscrtype) {
                case 0:
                    printf("=+=+=+= Convective Scheme: Central =+=+=+=\n");
                    break;
                case 1:
                    printf("=+=+=+= Convective Scheme: First Order Upwind =+=+=+=\n");
                    break;
                case 2:
                    printf("=+=+=+= Convective Scheme: Second Order Upwind =+=+=+=\n");
                    switch (ns->contr.secondconvecdiscrtype) {
                        case 0:
                            printf("=+= Convective Scheme Type : Modified Coefficient Upwind =+=\n");
                            break;
                        case 1:
                            printf("=+= Convective Scheme Type : Cubista =+=\n");
                            break;
                        case 2:
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
            fprintf(fd, "%d\n", (ns->contr.projtype));
            fprintf(fd, "%d\n", (ns->contr.flowtype));
            fprintf(fd, "%d\n", (ns->contr.modelflowtype));
            fprintf(fd, "%d\n", (ns->contr.tempdiscrtype));
            fprintf(fd, "%d\n", (ns->contr.spatialdiscrtype));
            fprintf(fd, "%d\n", (ns->contr.convecdiscrtype));
            fprintf(fd, "%d\n", (ns->contr.secondconvecdiscrtype));
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
    // Velocity file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.vel", ns->par.namesave);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL) {
        // Loading the velocity
        loadUV(ns->psfdu, ns->dpu, fd, myrank, ntasks);
        fclose(fd);
    } else {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
    // Pressure file name
    snprintf(namefile, sizeof namefile, "%s.pres", ns->par.namesave);
    fd = fopen(namefile, "r");
    if (fd != NULL) {
        // Loading the pressure
        loadP(ns->psdp, ns->dpp, fd, myrank, ntasks);
        fclose(fd);
    } else {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
    // FracVol file name
    if (ns->contr.flowtype == 2) {
        snprintf(namefile, sizeof namefile, "%s.fracvol", ns->par.namesave);
        fd = fopen(namefile, "r");
        if (fd != NULL) {
            // Loading the pressure
            loadVF(ns->psdp, ns->ed.mult.dpfracvol, fd, myrank, ntasks);
            fclose(fd);
        } else {
            // Error in open the file
            printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}

// Saving the properties
void higflow_save_properties(higflow_solver *ns, int myrank, int ntasks) {
    // Velocity file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.vel", ns->par.namesave);
    FILE *fd = fopen(namefile, "w");
    if (fd != NULL) {
        // Saving the velocity
        saveUV(ns->psfdu, ns->dpu, fd, myrank, ntasks);
        fclose(fd);
    } else {
        // Error in open the file
        printf("=+=+=+= Error saving file %s =+=+=+=\n", ns->par.namesave);
        exit(1);
    }
    // Pressure file name
    snprintf(namefile, sizeof namefile, "%s.pres", ns->par.namesave);
    fd = fopen(namefile, "w");
    if (fd != NULL) {
        // Saving the pressure
        saveP(ns->psdp, ns->dpp, fd, myrank, ntasks);
        fclose(fd);
    } else {
        // Error in open the file
        printf("=+=+=+= Error saving file %s =+=+=+=\n", ns->par.namesave);
        exit(1);
    }
    // FracVol file name
    if (ns->contr.flowtype == 2) {
        snprintf(namefile, sizeof namefile, "%s.fracvol", ns->par.namesave);
        fd = fopen(namefile, "w");
        if (fd != NULL) {
            // Saving the Volume Fraction
            saveVF(ns->psdp, ns->ed.mult.dpfracvol, fd, myrank, ntasks);
            fclose(fd);
        } else {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", ns->par.namesave);
            exit(1);
        }
        // Saving Tensor
        snprintf(namefile, sizeof namefile, "%s.tens", ns->par.namesave);
        fd = fopen(namefile, "w");
        if (fd != NULL) {
            // Saving the Volume Fraction
            saveS(ns->psdp, ns->ed.mult.dpS, fd, myrank, ntasks);
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
        ifd = fscanf(fd, "%d", &(ns->ed.ve.contr.model));
        ifd = fscanf(fd, "%d", &(ns->ed.ve.contr.discrtype));
        ifd = fscanf(fd, "%d", &(ns->ed.ve.contr.convecdiscrtype));
        fclose(fd);
        if (myrank == 0) {
            switch (ns->ed.ve.contr.model) {
                case -1:
                    printf("=+=+=+= Constitutive Equation Model: User Defined =+=+=+=\n");
                    break;
                case 0:
                    printf("=+=+=+= Constitutive Equation Model: Oldroyd =+=+=+=\n");
                    break;
                case 1:
                    printf("=+=+=+= Constitutive Equation Model: Giesekus =+=+=+=\n");
                    break;
                case 2:
                    printf("=+=+=+= Constitutive Equation Model: LPTT =+=+=+=\n");
                    break;
            }
            switch (ns->ed.ve.contr.discrtype) {
                case 0:
                    printf("=+=+=+= Constitutive Equation Discretization: Explicit =+=+=+=\n");
                    break;
                case 1:
                    printf("=+=+=+= Constitutive Equation Discretization: Implicit =+=+=+=\n");
                    break;
            }
            switch (ns->ed.ve.contr.convecdiscrtype) {
                case 0:
                    printf("=+=+=+= Constitutive Equation Convective Term: Upwind  =+=+=+=\n");
                    break;
                case 1:
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
            fprintf(fd, "%d\n", (ns->ed.ve.contr.model));
            fprintf(fd, "%d\n", (ns->ed.ve.contr.discrtype));
            fprintf(fd, "%d\n", (ns->ed.ve.contr.convecdiscrtype));
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
        ifd = fscanf(fd, "%lf", &(ns->ed.eo.par.dphidx));
        fclose(fd);
        if (myrank == 0) {
            printf("=+=+=+= alpha: %e =+=+=+=\n", ns->ed.eo.par.alpha);
            printf("=+=+=+= delta: %e =+=+=+=\n", ns->ed.eo.par.delta);
            printf("=+=+=+= Peclet number: %f =+=+=+=\n", ns->ed.eo.par.Pe);
            printf("=+=+=+= Externel field dphidx : %e =+=+=+=\n", ns->ed.eo.par.dphidx);
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
        snprintf(namefile, sizeof namefile, "%s.viscpar", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL) {
            // Saving the parameters
            fprintf(fd, "%lf\n", (ns->ed.eo.par.alpha));
            fprintf(fd, "%lf\n", (ns->ed.eo.par.delta));
            fprintf(fd, "%lf\n", (ns->ed.eo.par.Pe));
            fprintf(fd, "%lf\n", (ns->ed.eo.par.dphidx));
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
        ifd = fscanf(fd, "%d", &(ns->ed.eo.contr.eo_model));
        ifd = fscanf(fd, "%d", &(ns->ed.eo.contr.convecdiscrtype));
        fclose(fd);
        if (myrank == 0) {
            switch (ns->ed.eo.contr.eo_model) {
                case 0:
                    printf("=+=+=+= Electroosmotic Model: PNP =+=+=+=\n");
                    break;
                case 1:
                    printf("=+=+=+= Electroosmotic Model: PB =+=+=+=\n");
                    break;
                case 2:
                    printf("=+=+=+= Electroosmotic Model: PBDH =+=+=+=\n");
                    break;
                case 3:
                    printf("=+=+=+= Electroosmotic Model: analytic PBDH =+=+=+=\n");
                    break;
            }
            switch (ns->ed.eo.contr.convecdiscrtype) {
                case 0:
                    printf("=+=+=+= Ionic concentration convective term: Central =+=+=+=\n");
                    break;
                case 1:
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
            fprintf(fd, "%d\n", (ns->ed.eo.contr.eo_model));
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
    snprintf(namefile, sizeof namefile, "%s.viscintpar", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL) {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%lf", &(ns->ed.im.par.De));
        ifd = fscanf(fd, "%lf", &(ns->ed.im.par.alpha)); //for damping function
        ifd = fscanf(fd, "%lf", &(ns->ed.im.par.beta));  //for damping function
        if (ns->ed.im.contr.model == 0) {
            ifd = fscanf(fd, "%d", &(ns->ed.im.par.M));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.a[0]));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.a[1]));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.a[2]));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.a[3]));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.a[4]));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.a[5]));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.a[6]));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.a[7]));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.lambda[0]));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.lambda[1]));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.lambda[2]));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.lambda[3]));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.lambda[4]));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.lambda[5]));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.lambda[6]));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.lambda[7]));
        } else if (ns->ed.im.contr.model == 1) {
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.alpha_frac));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.beta_frac));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.Phi1));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.Phi2));
        }
        fclose(fd);
        if (myrank == 0) {
            printf("=+=+=+= Deborah Number: %f =+=+=+=\n", ns->ed.im.par.De);
            printf("=+=+=+= Alpha: %f =+=+=+=\n", ns->ed.im.par.alpha);
            printf("=+=+=+= Beta: %f =+=+=+=\n", ns->ed.im.par.beta);
            if (ns->ed.im.contr.model == 0) {
                printf("=+=+=+= M: %d =+=+=+=\n", ns->ed.im.par.M);
                for (int i = 0; i < 8; i++) {
                    printf("=+=+=+= a(%d): %f =+=+=+=\n", i, ns->ed.im.par.a[i]);
                }
                for (int i = 0; i < 8; i++) {
                    printf("=+=+=+= Lambda(%d): %f =+=+=+=\n", i, ns->ed.im.par.lambda[i]);
                }
            } else if (ns->ed.im.contr.model == 1) {
                printf("=+=+=+= Alpha Fractional: %f =+=+=+=\n", ns->ed.im.par.alpha_frac);
                printf("=+=+=+= Beta Fractional: %f =+=+=+=\n", ns->ed.im.par.beta_frac);
                printf("=+=+=+= Phi1: %f =+=+=+=\n", ns->ed.im.par.Phi1);
                printf("=+=+=+= Phi2: %f =+=+=+=\n", ns->ed.im.par.Phi2);
            }
        }
    } else {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
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
            if (ns->ed.im.contr.model == 0) {
                fprintf(fd, "%d\n", (ns->ed.im.par.M));
                for (int k = 0; k < 8; k++) {
                    fprintf(fd, "%lf\n", (ns->ed.im.par.a[k]));
                }
                for (int k = 0; k < 8; k++) {
                    fprintf(fd, "%lf\n", (ns->ed.im.par.lambda[k]));
                }
            } else if (ns->ed.im.contr.model == 1) {
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
    snprintf(namefile, sizeof namefile, "%s.viscintcontr", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL) {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%d", &(ns->ed.im.contr.model));
        ifd = fscanf(fd, "%d", &(ns->ed.im.contr.discrtype));
        ifd = fscanf(fd, "%d", &(ns->ed.im.contr.convecdiscrtype));
        fclose(fd);
        if (myrank == 0) {
            switch (ns->ed.im.contr.model) {
                case 0:
                    printf("=+=+=+= Constitutive Equation Model: KBKZ-PSM =+=+=+=\n");
                    break;
                case 1:
                    printf("=+=+=+= Constitutive Equation Model: UCM =+=+=+=\n");
                    break;
            }
            switch (ns->ed.im.contr.discrtype) {
                case 0:
                    printf("=+=+=+= Constitutive Equation Discretization: Explicit =+=+=+=\n");
                    break;
                case 1:
                    printf("=+=+=+= Constitutive Equation Discretization: Implicit =+=+=+=\n");
                    break;
            }
            switch (ns->ed.im.contr.convecdiscrtype) {
                case 0:
                    printf("=+=+=+= Constitutive Equation Convective Term: Central  =+=+=+=\n");
                    break;
                case 1:
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
        snprintf(namefile, sizeof namefile, "%s.viscintcontr", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL) {
            // Saving the parameters
            fprintf(fd, "%d\n", (ns->ed.im.contr.model));
            fprintf(fd, "%d\n", (ns->ed.im.contr.discrtype));
            fprintf(fd, "%d\n", (ns->ed.im.contr.convecdiscrtype));
            fclose(fd);
        } else {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}
