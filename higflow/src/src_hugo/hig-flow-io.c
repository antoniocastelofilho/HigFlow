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
void higflow_print_electroosmotic_properties(higflow_solver *ns, FILE *data, int dimprint, real pprint)
{
    if (ns->contr.modelflowtype == 1)
    {
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
        for (it = sd_get_domain_celliterator(sdpsi); !higcit_isfinished(it); higcit_nextcell(it))
        {
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
            if ((ccenter[dimprint] < pprint + epsprint) && (ccenter[dimprint] > pprint - epsprint))
            {
                fprintf(data, "%15.10lf  %15.10lf  %15.12lf  %15.12lf  %15.12lf  %15.12lf  %15.12lf\n", ccenter[0], ccenter[1], phi, psi, np, nm, rhoe);
            }
        }
        // Destroy the iterator
        higcit_destroy(it);
    }
}

// Print the Polymeric Tensor
void higflow_print_polymeric_tensor(higflow_solver *ns, FILE *data, int dimprint, real pprint)
{
    if (ns->contr.flowtype == 3)
    {
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
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
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
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    // Get Du
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                    D[i][j] = 0.5 * (Du[i][j] + Du[j][i]);
                    // Get S tensor
                    S[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                }
            }
            //Print polymeric stress data file
            if (ccenter[dimprint] == pprint)
            {
                fprintf(data, "%lf  %lf  %15.12lf  %15.12lf  %15.12lf  %15.12lf\n", ccenter[0], ccenter[1], S[0][0] + 2.0 * (1 - beta) * D[0][0] / Re, S[0][1] + 2.0 * (1 - beta) * D[0][1] / Re, S[1][0] + 2.0 * (1 - beta) * D[1][0] / Re, S[1][1] + 2.0 * (1 - beta) * D[1][1] / Re);
            }
        }
        // Destroy the iterator
        higcit_destroy(it);
    }
}

// Print the velocity
void higflow_print_velocity(higflow_solver *ns, FILE *data, int dimprint, real pprint)
{
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Loop for each dimension
    higfit_facetiterator *fit;
    for (int dim = 0; dim < DIM; dim++)
    {
        // Initialize the min and max velocity
        // Get the local partitioned domain for facets
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the map of the distributed properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit))
        {
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
            if (fcenter[dimprint] == pprint)
            {
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
            case 6:
                //Viscoelastic with variable viscosity
                higflow_print_vtk3D_viscoelastic_variable_viscosity(ns, rank);
                break;
            case 7:
                //Viscoelastic with shear-banding
                higflow_print_vtk3D_viscoelastic_shear_banding(ns, rank);
                break;
            case 8:
                //Elastoviscoplastic
                higflow_print_vtk3D_elastoviscoplastic(ns, rank);
                break;
            case 9:
                //Shear-thickening suspensions
                higflow_print_vtk3D_shear_thickening_suspensions(ns, rank);
                break;
            }
            break;
    }
}

// Print the VTK file for visualize 2D
void higflow_print_vtk2D(higflow_solver *ns, int rank) {
    // Open the VTK file
    real Re = ns->par.Re;
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
    fprintf(f, "POINTS %ld float\n", 4 * numleafs);
    hig_cell *c;
    for (; !higcit_isfinished(it_t); higcit_nextcell(it_t))
    {
        c = higcit_getcell(it_t);
        fprintf(f, "%lf %lf 0\n", c->lowpoint[0], c->lowpoint[1]);
        fprintf(f, "%lf %lf 0\n", c->highpoint[0], c->lowpoint[1]);
        fprintf(f, "%lf %lf 0\n", c->highpoint[0], c->highpoint[1]);
        fprintf(f, "%lf %lf 0\n", c->lowpoint[0], c->highpoint[1]);
    }
    fprintf(f, "CELLS %ld %ld\n", numleafs, 5 * numleafs);
    for (int i = 0; i < numleafs; i++)
    {
        fprintf(f, "%d %d %d %d %d\n", 4, 4 * i, 4 * i + 1, 4 * i + 2, 4 * i + 3);
    }
    fprintf(f, "CELL_TYPES %ld\n", numleafs);
    for (int i = 0; i < numleafs; i++)
    {
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
    fprintf(f, "\nPOINT_DATA %ld\nVECTORS vel FLOAT\n", 4 * numleafs);
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
    {
        hig_cell *c = higcit_getcell(it);
        Point cdelta, ccenter, clowpoint, chightpoint;
        hig_get_delta(c, cdelta);
        hig_get_center(c, ccenter);
        // Pontos onde será interpolada a velocidade
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

        for (int dim = 0; dim < DIM; dim++)
        {
            // int numpts = 0;
            // int cont = 0;
            sfdu2[dim] = psfd_get_local_domain(ns->psfdu[dim]);
            mp_mapper *m = sfd_get_domain_mapper(sfdu2[dim]);
            real value0 = compute_facet_value_at_point(sfdu2[dim], ccenter, p0, 1.0, ns->dpu[dim], ns->stn);
            switch (dim)
            {
            case 0:
                lu0 = value0;
                break;
            case 1:
                lv0 = value0;
                break;
            }
            real value1 = compute_facet_value_at_point(sfdu2[dim], ccenter, p1, 1.0, ns->dpu[dim], ns->stn);
            switch (dim)
            {
            case 0:
                lu1 = value1;
                break;
            case 1:
                lv1 = value1;
                break;
            }
            real value2 = compute_facet_value_at_point(sfdu2[dim], ccenter, p2, 1.0, ns->dpu[dim], ns->stn);
            switch (dim)
            {
            case 0:
                lu2 = value2;
                break;
            case 1:
                lv2 = value2;
                break;
            }
            real value3 = compute_facet_value_at_point(sfdu2[dim], ccenter, p3, 1.0, ns->dpu[dim], ns->stn);
            switch (dim)
            {
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

    switch (ns->contr.flowtype)
    {
    case 0:
        // Newtonian
        fprintf(f, "\nCELL_DATA %ld\nSCALARS p FLOAT\nLOOKUP_TABLE default\n", numcells);
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
            hig_cell *c = higcit_getcell(it);
            uniqueid id = hig_get_cid(c);
            int cgid = mp_lookup(m, id);
            real val = dp_get_value(ns->dpp, cgid);
            fprintf(f, "%e\n", val);
        }
        higcit_destroy(it);

        //Non-isothermal flows
        if (ns->contr.modelflowtype == 5){
            real T0 = ns->ed.nif.par.T0;
            real T1 = ns->ed.nif.par.T1;
            //Printing the temperature
            sim_domain *sdna =  psd_get_local_domain(ns->ed.nif.psdT);
            fprintf(f, "\nSCALARS nA FLOAT\nLOOKUP_TABLE default\n");
            for (it = sd_get_domain_celliterator(sdna); !higcit_isfinished(it); higcit_nextcell(it))
            {
                hig_cell *c = higcit_getcell(it);
                uniqueid id = hig_get_cid(c);
                int cgid = mp_lookup(m, id);
                Point cdelta, ccenter, clowpoint, chightpoint;
                hig_get_delta(c, cdelta);
                //Point ccenter;
                hig_get_center(c, ccenter);
                //Calculation of the dimensionless temperature val0 = (T-T0)/(T1-T0)
                real val0 = compute_value_at_point(ns->ed.nif.sdT, ccenter, ccenter, 1.0, ns->ed.nif.dpT, ns->ed.stn);
                //Printing the dimensional temperature T = (T1-T0)*val0 + T0
                real val = (T1-T0)*val0 + T0;
                fprintf(f, "%e\n", val);
            }
            higcit_destroy(it);

        }

        break;

    case 1:
        // Generalized Newtonian
        fprintf(f, "\nCELL_DATA %ld\nSCALARS p FLOAT\nLOOKUP_TABLE default\n", numcells);
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
            hig_cell *c = higcit_getcell(it);
            uniqueid id = hig_get_cid(c);
            int cgid = mp_lookup(m, id);
            real val = dp_get_value(ns->dpp, cgid);
            fprintf(f, "%e\n", val);
        }
        higcit_destroy(it);

        sim_domain *sdviscgn =  psd_get_local_domain(ns->ed.psdED);
        fprintf(f, "\nSCALARS viscosity FLOAT\nLOOKUP_TABLE default\n");
        for (it = sd_get_domain_celliterator(sdviscgn); !higcit_isfinished(it); higcit_nextcell(it))
        {
            hig_cell *c = higcit_getcell(it);
            uniqueid id = hig_get_cid(c);
            int cgid = mp_lookup(m, id);
            Point cdelta, ccenter, clowpoint, chightpoint;
            hig_get_delta(c, cdelta);
            //Point ccenter;
            hig_get_center(c, ccenter);
            real val = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.gn.dpvisc, ns->ed.stn);
            //real val = dp_get_value(ns->ed.vevv.dpvisc, cgid);
            fprintf(f, "%e\n", val);
        }
        higcit_destroy(it);

        //Non-isothermal flows
        if (ns->contr.modelflowtype == 5){
            real T0 = ns->ed.nif.par.T0;
            real T1 = ns->ed.nif.par.T1;
            //Printing the temperature
            sim_domain *sdna =  psd_get_local_domain(ns->ed.nif.psdT);
            fprintf(f, "\nSCALARS Temperature FLOAT\nLOOKUP_TABLE default\n");
            for (it = sd_get_domain_celliterator(sdna); !higcit_isfinished(it); higcit_nextcell(it))
            {
                hig_cell *c = higcit_getcell(it);
                uniqueid id = hig_get_cid(c);
                int cgid = mp_lookup(m, id);
                Point cdelta, ccenter, clowpoint, chightpoint;
                hig_get_delta(c, cdelta);
                //Point ccenter;
                hig_get_center(c, ccenter);
                //Calculation of the dimensionless temperature val0 = (T-T0)/(T1-T0)
                real val0 = compute_value_at_point(ns->ed.nif.sdT, ccenter, ccenter, 1.0, ns->ed.nif.dpT, ns->ed.stn);
                //Printing the dimensional temperature T = (T1-T0)*val0 + T0
                real val = (T1-T0)*val0 + T0;
                fprintf(f, "%e\n", val);
            }
            higcit_destroy(it);

        }

        break;

    case 2:
        // Multiphase
        fprintf(f, "\nCELL_DATA %ld\nSCALARS p FLOAT\nLOOKUP_TABLE default\n", numcells);
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
            hig_cell *c = higcit_getcell(it);
            uniqueid id = hig_get_cid(c);
            int cgid = mp_lookup(m, id);
            real val = dp_get_value(ns->dpp, cgid);
            fprintf(f, "%e\n", val);
        }
        higcit_destroy(it);

        fprintf(f, "\nSCALARS FracVol FLOAT\nLOOKUP_TABLE default\n");
        sfdu2[0] = psfd_get_local_domain(ns->psfdu[0]);
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
            hig_cell *c = higcit_getcell(it);
            uniqueid id = hig_get_cid(c);
            int cgid = mp_lookup(m, id);
            Point cdelta, ccenter, p0, p1;
            hig_get_delta(c, cdelta);
            hig_get_center(c, ccenter);
            real value = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
            fprintf(f, "%e\n", value);
        }
        higcit_destroy(it);
        break;

    case 3:
        // Viscoelastic
        /* fprintf(f, "\nCELL_DATA %ld\nTENSORS Stensor FLOAT\n", numcells);
                   for(it = sd_get_domain_celliterator(ns->ed.sdED); !higcit_isfinished(it); higcit_nextcell(it)) {
                       hig_cell *c = higcit_getcell(it);
                       uniqueid id = hig_get_cid(c);
                       int cgid = mp_lookup(m, id);
                       real val00 = dp_get_value(ns->ed.im.dpS[0][0], cgid);
                       real val01 = dp_get_value(ns->ed.im.dpS[0][1], cgid);
                       real val02 = 0.0;
                       real val10 = dp_get_value(ns->ed.im.dpS[1][0], cgid);
                       real val11 = dp_get_value(ns->ed.im.dpS[1][1], cgid);
                       real val12 = 0.0;
                       real val20 = 0.0;
                       real val21 = 0.0;
                       real val22 = 0.0;
                       fprintf(f, "%e %e %e\n", val00, val01, val02);
                       fprintf(f, "%e %e %e\n", val10, val11, val12);
                       fprintf(f, "%e %e %e\n", val20, val21, val22);
                   }
                   higcit_destroy(it); */

        fprintf(f, "\nTENSORS stress FLOAT\n");

        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
            hig_cell *c = higcit_getcell(it);
            Point cdelta, ccenter, clowpoint, chightpoint;
            hig_get_delta(c, cdelta);
            hig_get_center(c, ccenter);
            // Pontos onde será interpolada a velocidade
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

            real taup0[DIM + 1][DIM + 1], taup1[DIM + 1][DIM + 1], taup2[DIM + 1][DIM + 1], taup3[DIM + 1][DIM + 1];

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    taup0[i][j] = 0.0;
                    taup1[i][j] = 0.0;
                    taup2[i][j] = 0.0;
                    taup3[i][j] = 0.0;
                }
            }

            //  double Sxx0 = 0, Sxy0 = 0, Sxz0 = 0, Syx0 = 0, Syy0 = 0, Syz0 = 0, Szx0 = 0, Szy0 = 0, Szz0 = 0;
            //  double Sxx1 = 0, Sxy1 = 0, Sxz1 = 0, Syx1 = 0, Syy1 = 0, Syz1 = 0, Szx1 = 0, Szy1 = 0, Szz1 = 0;
            //  double Sxx2 = 0, Sxy2 = 0, Sxz2 = 0, Syx2 = 0, Syy2 = 0, Syz2 = 0, Szx2 = 0, Szy2 = 0, Szz2 = 0;
            //  double Sxx3 = 0, Sxy3 = 0, Sxz3 = 0, Syx3 = 0, Syy3 = 0, Syz3 = 0, Szx3 = 0, Szy3 = 0, Szz3 = 0;

            sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);

            real Dp0[DIM][DIM], Dp1[DIM][DIM], Dp2[DIM][DIM], Dp3[DIM][DIM];
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    // Get Du
                    Dp0[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p0, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                    Dp1[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p1, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                    Dp2[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p2, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                    Dp3[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p3, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
                }
            }

            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    taup0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                    taup0[i][j] += 2.0 * (1 - beta) * 0.5 * (Dp0[i][j] + Dp0[j][i]) / Re;
                    taup1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                    taup1[i][j] += 2.0 * (1 - beta) * 0.5 * (Dp1[i][j] + Dp1[j][i]) / Re;
                    taup2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                    taup2[i][j] += 2.0 * (1 - beta) * 0.5 * (Dp2[i][j] + Dp2[j][i]) / Re;
                    taup3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
                    taup3[i][j] += 2.0 * (1 - beta) * 0.5 * (Dp3[i][j] + Dp3[j][i]) / Re;
                }
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup0[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup1[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup2[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup3[i][j]);
                }
                fprintf(f, "\n");
            }

            /*  
                        Sxx0 = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.ve.dpS[0][0], ns->ed.stn);
                        Sxx0 += 2.0*(1-beta)*Dp0[0][0]/Re; 
                        Sxy0 = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.ve.dpS[0][1], ns->ed.stn);
                        Sxy0 += 2.0*(1-beta)*Dp0[1][0]/Re;
                        Syx0 = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.ve.dpS[1][0], ns->ed.stn);
                        Syx0 += 2.0*(1-beta)*Dp0[0][1]/Re;
                        Syy0 = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.ve.dpS[1][1], ns->ed.stn);
                        Syy0 += 2.0*(1-beta)*Dp0[1][1]/Re;
                        
                        Sxx1 = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.ve.dpS[0][0], ns->ed.stn);
                        Sxx1 += 2.0*(1-beta)*Dp1[0][0]/Re; 
                        Sxy1 = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.ve.dpS[0][1], ns->ed.stn);
                        Sxy1 += 2.0*(1-beta)*Dp1[1][0]/Re;
                        Syx1 = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.ve.dpS[1][0], ns->ed.stn);
                        Syx1 += 2.0*(1-beta)*Dp1[0][1]/Re;
                        Syy1 = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.ve.dpS[1][1], ns->ed.stn);
                        Syy1 += 2.0*(1-beta)*Dp1[1][1]/Re;

                        Sxx2 = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.ve.dpS[0][0], ns->ed.stn);
                        Sxx2 += 2.0*(1-beta)*Dp2[0][0]/Re;
                        Sxy2 = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.ve.dpS[0][1], ns->ed.stn);
                        Sxy2 += 2.0*(1-beta)*Dp2[1][0]/Re;
                        Syx2 = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.ve.dpS[1][0], ns->ed.stn);
                        Syx2 += 2.0*(1-beta)*Dp2[0][1]/Re;
                        Syy2 = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.ve.dpS[1][1], ns->ed.stn);
                        Syy2 += 2.0*(1-beta)*Dp2[1][1]/Re;
                        
                        Sxx3 = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.ve.dpS[0][0], ns->ed.stn);
                        Sxx3 += 2.0*(1-beta)*Dp3[0][0]/Re; 
                        Sxy3 = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.ve.dpS[0][1], ns->ed.stn);
                        Sxy3 += 2.0*(1-beta)*Dp3[1][0]/Re;
                        Syx3 = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.ve.dpS[1][0], ns->ed.stn);
                        Syx3 += 2.0*(1-beta)*Dp3[0][1]/Re;
                        Syy3 = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.ve.dpS[1][1], ns->ed.stn);
                        Syy3 += 2.0*(1-beta)*Dp3[1][1]/Re;
    
                        fprintf(f, "%e %e %e\n", Sxx0, Sxy0, Sxz0);
                        fprintf(f, "%e %e %e\n", Syx0, Syy0, Syz0);
                        fprintf(f, "%e %e %e\n", Szx0, Szy0, Szz0);
                        fprintf(f, "%e %e %e\n", Sxx1, Sxy1, Sxz1);
                        fprintf(f, "%e %e %e\n", Syx1, Syy1, Syz1);
                        fprintf(f, "%e %e %e\n", Szx1, Szy1, Szz1);
                        fprintf(f, "%e %e %e\n", Sxx2, Sxy2, Sxz2);
                        fprintf(f, "%e %e %e\n", Syx2, Syy2, Syz2);
                        fprintf(f, "%e %e %e\n", Szx2, Szy2, Szz2);
                        fprintf(f, "%e %e %e\n", Sxx3, Sxy3, Sxz3);
                        fprintf(f, "%e %e %e\n", Syx3, Syy3, Syz3);
                        fprintf(f, "%e %e %e\n", Szx3, Szy3, Szz3);   */
        }
        higcit_destroy(it);

        //fprintf(f, "\nSCALARS p FLOAT\nLOOKUP_TABLE default\n");
        fprintf(f, "\nCELL_DATA %ld\nSCALARS p FLOAT\nLOOKUP_TABLE default\n", numcells);
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
            hig_cell *c = higcit_getcell(it);
            uniqueid id = hig_get_cid(c);
            int cgid = mp_lookup(m, id);
            real val = dp_get_value(ns->dpp, cgid);
            fprintf(f, "%e\n", val);
        }
        higcit_destroy(it);

        //Non-isothermal flows
        if (ns->contr.modelflowtype == 5){
            real T0 = ns->ed.nif.par.T0;
            real T1 = ns->ed.nif.par.T1;
            //Printing the temperature
            sim_domain *sdna =  psd_get_local_domain(ns->ed.nif.psdT);
            fprintf(f, "\nSCALARS Temperature FLOAT\nLOOKUP_TABLE default\n");
            for (it = sd_get_domain_celliterator(sdna); !higcit_isfinished(it); higcit_nextcell(it))
            {
                hig_cell *c = higcit_getcell(it);
                uniqueid id = hig_get_cid(c);
                int cgid = mp_lookup(m, id);
                Point cdelta, ccenter, clowpoint, chightpoint;
                hig_get_delta(c, cdelta);
                //Point ccenter;
                hig_get_center(c, ccenter);
                //Calculation of the dimensionless temperature val0 = (T-T0)/(T1-T0)
                real val0 = compute_value_at_point(ns->ed.nif.sdT, ccenter, ccenter, 1.0, ns->ed.nif.dpT, ns->ed.stn);
                //Printing the dimensional temperature T = (T1-T0)*val0 + T0
                //real val = (T1-T0)*val0 + T0;
                real val = (T0-T1)*val0 + T1;
                fprintf(f, "%e\n", val);
                //fprintf(f, "%e\n", val0);
            }
            higcit_destroy(it);

        }

        break;

    case 4:
        // Integral
        /* fprintf(f, "\nCELL_DATA %ld\nTENSORS Stensor FLOAT\n", numcells);
                   for(it = sd_get_domain_celliterator(ns->ed.sdED); !higcit_isfinished(it); higcit_nextcell(it)) {
                       hig_cell *c = higcit_getcell(it);
                       uniqueid id = hig_get_cid(c);
                       int cgid = mp_lookup(m, id);
                       real val00 = dp_get_value(ns->ed.im.dpS[0][0], cgid);
                       real val01 = dp_get_value(ns->ed.im.dpS[0][1], cgid);
                       real val02 = 0.0;
                       real val10 = dp_get_value(ns->ed.im.dpS[1][0], cgid);
                       real val11 = dp_get_value(ns->ed.im.dpS[1][1], cgid);
                       real val12 = 0.0;
                       real val20 = 0.0;
                       real val21 = 0.0;
                       real val22 = 0.0;
                       fprintf(f, "%e %e %e\n", val00, val01, val02);
                       fprintf(f, "%e %e %e\n", val10, val11, val12);
                       fprintf(f, "%e %e %e\n", val20, val21, val22);
                   }
                   higcit_destroy(it); */

        fprintf(f, "\nTENSORS Stensor FLOAT\n");
        //real Re   = ns->par.Re;
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
            hig_cell *c = higcit_getcell(it);
            Point cdelta, ccenter, clowpoint, chightpoint;
            hig_get_delta(c, cdelta);
            hig_get_center(c, ccenter);
            // Pontos onde será interpolada a velocidade
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

            real taup0[DIM + 1][DIM + 1], taup1[DIM + 1][DIM + 1], taup2[DIM + 1][DIM + 1], taup3[DIM + 1][DIM + 1];

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    taup0[i][j] = 0.0;
                    taup1[i][j] = 0.0;
                    taup2[i][j] = 0.0;
                    taup3[i][j] = 0.0;
                }
            }

            //double Sxx0 = 0, Sxy0 = 0, Sxz0 = 0, Syx0 = 0, Syy0 = 0, Syz0 = 0, Szx0 = 0, Szy0 = 0, Szz0 = 0;
            //double Sxx1 = 0, Sxy1 = 0, Sxz1 = 0, Syx1 = 0, Syy1 = 0, Syz1 = 0, Szx1 = 0, Szy1 = 0, Szz1 = 0;
            //double Sxx2 = 0, Sxy2 = 0, Sxz2 = 0, Syx2 = 0, Syy2 = 0, Syz2 = 0, Szx2 = 0, Szy2 = 0, Szz2 = 0;
            //double Sxx3 = 0, Sxy3 = 0, Sxz3 = 0, Syx3 = 0, Syy3 = 0, Syz3 = 0, Szx3 = 0, Szy3 = 0, Szz3 = 0;

            sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);

            real Dp0[DIM][DIM], Dp1[DIM][DIM], Dp2[DIM][DIM], Dp3[DIM][DIM];
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    // Get Du
                    Dp0[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p0, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
                    Dp1[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p1, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
                    Dp2[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p2, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
                    Dp3[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p3, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
                }
            }

            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    taup0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
                    taup0[i][j] -= (Dp0[i][j] + Dp0[j][i]) / Re;
                    taup1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
                    taup1[i][j] -= (Dp1[i][j] + Dp1[j][i]) / Re;
                    taup2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
                    taup2[i][j] -= (Dp2[i][j] + Dp2[j][i]) / Re;
                    taup3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
                    taup3[i][j] -= (Dp3[i][j] + Dp3[j][i]) / Re;
                }
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup0[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup1[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup2[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup3[i][j]);
                }
                fprintf(f, "\n");
            }
            /* Sxx0 = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.im.dpS[0][0], ns->ed.stn);
                        Sxy0 = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.im.dpS[0][1], ns->ed.stn);    
                        Syx0 = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.im.dpS[1][0], ns->ed.stn);
                        Syy0 = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.im.dpS[1][1], ns->ed.stn);
                        
                        Sxx1 = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.im.dpS[0][0], ns->ed.stn);
                        Sxy1 = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.im.dpS[0][1], ns->ed.stn);    
                        Syx1 = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.im.dpS[1][0], ns->ed.stn);
                        Syy1 = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.im.dpS[1][1], ns->ed.stn);
                        
                        Sxx2 = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.im.dpS[0][0], ns->ed.stn);
                        Sxy2 = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.im.dpS[0][1], ns->ed.stn);    
                        Syx2 = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.im.dpS[1][0], ns->ed.stn);
                        Syy2 = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.im.dpS[1][1], ns->ed.stn);
                        
                        Sxx3 = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.im.dpS[0][0], ns->ed.stn);
                        Sxy3 = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.im.dpS[0][1], ns->ed.stn);    
                        Syx3 = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.im.dpS[1][0], ns->ed.stn);
                        Syy3 = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.im.dpS[1][1], ns->ed.stn);
    
                        fprintf(f, "%e %e %e\n", Sxx0, Sxy0, Sxz0);
                        fprintf(f, "%e %e %e\n", Syx0, Syy0, Syz0);
                        fprintf(f, "%e %e %e\n", Szx0, Szy0, Szz0);
                        fprintf(f, "%e %e %e\n", Sxx1, Sxy1, Sxz1);
                        fprintf(f, "%e %e %e\n", Syx1, Syy1, Syz1);
                        fprintf(f, "%e %e %e\n", Szx1, Szy1, Szz1);
                        fprintf(f, "%e %e %e\n", Sxx2, Sxy2, Sxz2);
                        fprintf(f, "%e %e %e\n", Syx2, Syy2, Syz2);
                        fprintf(f, "%e %e %e\n", Szx2, Szy2, Szz2);
                        fprintf(f, "%e %e %e\n", Sxx3, Sxy3, Sxz3);
                        fprintf(f, "%e %e %e\n", Syx3, Syy3, Syz3);
                        fprintf(f, "%e %e %e\n", Szx3, Szy3, Szz3);  */
        }
        higcit_destroy(it);

        //fprintf(f, "\nSCALARS p FLOAT\nLOOKUP_TABLE default\n");
        fprintf(f, "\nCELL_DATA %ld\nSCALARS p FLOAT\nLOOKUP_TABLE default\n", numcells);
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
            hig_cell *c = higcit_getcell(it);
            uniqueid id = hig_get_cid(c);
            int cgid = mp_lookup(m, id);
            real val = dp_get_value(ns->dpp, cgid);
            fprintf(f, "%e\n", val);
        }
        higcit_destroy(it);
        break;

    case 6:
        //sim_domain *sdte =  psd_get_local_domain(ns->ed.psdED);
        fprintf(f, "\nTENSORS stress FLOAT\n");

        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
            hig_cell *c = higcit_getcell(it);
            Point cdelta, ccenter, clowpoint, chightpoint;
            hig_get_delta(c, cdelta);
            hig_get_center(c, ccenter);
            // Pontos onde será interpolada a velocidade
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

            real taup0[DIM + 1][DIM + 1], taup1[DIM + 1][DIM + 1], taup2[DIM + 1][DIM + 1], taup3[DIM + 1][DIM + 1];

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    taup0[i][j] = 0.0;
                    taup1[i][j] = 0.0;
                    taup2[i][j] = 0.0;
                    taup3[i][j] = 0.0;
                }
            }

            //  double Sxx0 = 0, Sxy0 = 0, Sxz0 = 0, Syx0 = 0, Syy0 = 0, Syz0 = 0, Szx0 = 0, Szy0 = 0, Szz0 = 0;
            //  double Sxx1 = 0, Sxy1 = 0, Sxz1 = 0, Syx1 = 0, Syy1 = 0, Syz1 = 0, Szx1 = 0, Szy1 = 0, Szz1 = 0;
            //  double Sxx2 = 0, Sxy2 = 0, Sxz2 = 0, Syx2 = 0, Syy2 = 0, Syz2 = 0, Szx2 = 0, Szy2 = 0, Szz2 = 0;
            //  double Sxx3 = 0, Sxy3 = 0, Sxz3 = 0, Syx3 = 0, Syy3 = 0, Syz3 = 0, Szx3 = 0, Szy3 = 0, Szz3 = 0;

            sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);

            real Dp0[DIM][DIM], Dp1[DIM][DIM], Dp2[DIM][DIM], Dp3[DIM][DIM];
            //for (int i = 0; i < DIM; i++)
            //{
            //    for (int j = 0; j < DIM; j++)
            //    {
            //        // Get Du
            //        Dp0[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p0, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
            //        Dp1[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p1, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
            //        Dp2[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p2, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
            //        Dp3[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p3, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
            //    }
            //}
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    // Get Du
                    Dp0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
                    Dp1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
                    Dp2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
                    Dp3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
                }
            }

            //Get the viscosity value
            //real eta0 = compute_value_at_point(ns->ed.vevv.sdVisc, ccenter, p0, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            //real eta1 = compute_value_at_point(ns->ed.vevv.sdVisc, ccenter, p1, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            //real eta2 = compute_value_at_point(ns->ed.vevv.sdVisc, ccenter, p2, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            //real eta3 = compute_value_at_point(ns->ed.vevv.sdVisc, ccenter, p3, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);

            //real eta0 = compute_value_at_point(ns->ed.vevv.sdVisc, p0, p0, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            //real eta1 = compute_value_at_point(ns->ed.vevv.sdVisc, p1, p1, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            //real eta2 = compute_value_at_point(ns->ed.vevv.sdVisc, p2, p2, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            //real eta3 = compute_value_at_point(ns->ed.vevv.sdVisc, p3, p3, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            

            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    taup0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.vevv.dpS[i][j], ns->ed.stn);
                    taup0[i][j] += 2.0 * (1.0 - ns->ed.vevv.par.beta) * 0.5 * (Dp0[i][j] + Dp0[j][i]) / Re;
                    //taup0[i][j] += 2.0* (1.0- ns->ed.vevv.par.beta)* eta0 * 0.5 * (Dp0[i][j] + Dp0[j][i])/Re -2.0*(ns->ed.vevv.par.beta)/Re* 0.5 * (Dp0[i][j] + Dp0[j][i]);
                    taup1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.vevv.dpS[i][j], ns->ed.stn);
                    taup1[i][j] += 2.0 * (1.0 - ns->ed.vevv.par.beta) * 0.5 * (Dp1[i][j] + Dp1[j][i]) / Re;
                    //taup1[i][j] += 2.0* (1.0- ns->ed.vevv.par.beta)* eta1 * 0.5 * (Dp1[i][j] + Dp1[j][i])/Re -2.0*(ns->ed.vevv.par.beta)/Re* 0.5 * (Dp1[i][j] + Dp1[j][i]);
                    taup2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.vevv.dpS[i][j], ns->ed.stn);
                    taup2[i][j] += 2.0 * (1.0 - ns->ed.vevv.par.beta) * 0.5 * (Dp2[i][j] + Dp2[j][i]) / Re;
                    //taup2[i][j] += 2.0* (1.0- ns->ed.vevv.par.beta)* eta2 * 0.5 * (Dp2[i][j] + Dp2[j][i])/Re -2.0*(ns->ed.vevv.par.beta)/Re* 0.5 * (Dp2[i][j] + Dp2[j][i]);
                    taup3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.vevv.dpS[i][j], ns->ed.stn);
                    taup3[i][j] += 2.0 * (1.0 - ns->ed.vevv.par.beta) * 0.5 * (Dp3[i][j] + Dp3[j][i]) / Re;
                    //taup3[i][j] += 2.0* (1.0- ns->ed.vevv.par.beta)* eta3 * 0.5 * (Dp3[i][j] + Dp3[j][i])/Re -2.0*(ns->ed.vevv.par.beta)/Re* 0.5 * (Dp3[i][j] + Dp3[j][i]);
                }
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup0[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup1[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup2[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup3[i][j]);
                }
                fprintf(f, "\n");
            }

            /*  
                        Sxx0 = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.ve.dpS[0][0], ns->ed.stn);
                        Sxx0 += 2.0*(1-beta)*Dp0[0][0]/Re; 
                        Sxy0 = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.ve.dpS[0][1], ns->ed.stn);
                        Sxy0 += 2.0*(1-beta)*Dp0[1][0]/Re;
                        Syx0 = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.ve.dpS[1][0], ns->ed.stn);
                        Syx0 += 2.0*(1-beta)*Dp0[0][1]/Re;
                        Syy0 = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.ve.dpS[1][1], ns->ed.stn);
                        Syy0 += 2.0*(1-beta)*Dp0[1][1]/Re;
                        
                        Sxx1 = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.ve.dpS[0][0], ns->ed.stn);
                        Sxx1 += 2.0*(1-beta)*Dp1[0][0]/Re; 
                        Sxy1 = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.ve.dpS[0][1], ns->ed.stn);
                        Sxy1 += 2.0*(1-beta)*Dp1[1][0]/Re;
                        Syx1 = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.ve.dpS[1][0], ns->ed.stn);
                        Syx1 += 2.0*(1-beta)*Dp1[0][1]/Re;
                        Syy1 = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.ve.dpS[1][1], ns->ed.stn);
                        Syy1 += 2.0*(1-beta)*Dp1[1][1]/Re;

                        Sxx2 = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.ve.dpS[0][0], ns->ed.stn);
                        Sxx2 += 2.0*(1-beta)*Dp2[0][0]/Re;
                        Sxy2 = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.ve.dpS[0][1], ns->ed.stn);
                        Sxy2 += 2.0*(1-beta)*Dp2[1][0]/Re;
                        Syx2 = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.ve.dpS[1][0], ns->ed.stn);
                        Syx2 += 2.0*(1-beta)*Dp2[0][1]/Re;
                        Syy2 = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.ve.dpS[1][1], ns->ed.stn);
                        Syy2 += 2.0*(1-beta)*Dp2[1][1]/Re;
                        
                        Sxx3 = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.ve.dpS[0][0], ns->ed.stn);
                        Sxx3 += 2.0*(1-beta)*Dp3[0][0]/Re; 
                        Sxy3 = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.ve.dpS[0][1], ns->ed.stn);
                        Sxy3 += 2.0*(1-beta)*Dp3[1][0]/Re;
                        Syx3 = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.ve.dpS[1][0], ns->ed.stn);
                        Syx3 += 2.0*(1-beta)*Dp3[0][1]/Re;
                        Syy3 = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.ve.dpS[1][1], ns->ed.stn);
                        Syy3 += 2.0*(1-beta)*Dp3[1][1]/Re;
    
                        fprintf(f, "%e %e %e\n", Sxx0, Sxy0, Sxz0);
                        fprintf(f, "%e %e %e\n", Syx0, Syy0, Syz0);
                        fprintf(f, "%e %e %e\n", Szx0, Szy0, Szz0);
                        fprintf(f, "%e %e %e\n", Sxx1, Sxy1, Sxz1);
                        fprintf(f, "%e %e %e\n", Syx1, Syy1, Syz1);
                        fprintf(f, "%e %e %e\n", Szx1, Szy1, Szz1);
                        fprintf(f, "%e %e %e\n", Sxx2, Sxy2, Sxz2);
                        fprintf(f, "%e %e %e\n", Syx2, Syy2, Syz2);
                        fprintf(f, "%e %e %e\n", Szx2, Szy2, Szz2);
                        fprintf(f, "%e %e %e\n", Sxx3, Sxy3, Sxz3);
                        fprintf(f, "%e %e %e\n", Syx3, Syy3, Syz3);
                        fprintf(f, "%e %e %e\n", Szx3, Szy3, Szz3);   */
        }
        higcit_destroy(it);

        //fprintf(f, "\nSCALARS p FLOAT\nLOOKUP_TABLE default\n");
        fprintf(f, "\nCELL_DATA %ld\nSCALARS p FLOAT\nLOOKUP_TABLE default\n", numcells);
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
            hig_cell *c = higcit_getcell(it);
            uniqueid id = hig_get_cid(c);
            int cgid = mp_lookup(m, id);
            real val = dp_get_value(ns->dpp, cgid);
            fprintf(f, "%e\n", val);
        }
        higcit_destroy(it);

        //Printing the viscosity
        sim_domain *sdvisc =  psd_get_local_domain(ns->ed.vevv.psdVisc);
        //fprintf(f, "\nCELL_DATA %ld\nSCALARS viscosity FLOAT\nLOOKUP_TABLE default\n", numcells);
        //fprintf(f, "\nPOINT_DATA %ld\nSCALARS viscosity FLOAT\n", 4 * numleafs);
        fprintf(f, "\nSCALARS viscosity FLOAT\nLOOKUP_TABLE default\n");
        for (it = sd_get_domain_celliterator(sdvisc); !higcit_isfinished(it); higcit_nextcell(it))
        {
            hig_cell *c = higcit_getcell(it);
            uniqueid id = hig_get_cid(c);
            int cgid = mp_lookup(m, id);
            Point cdelta, ccenter, clowpoint, chightpoint;
            hig_get_delta(c, cdelta);
            //Point ccenter;
            hig_get_center(c, ccenter);
            //Point p0, p1, p2, p3;
            //p0[0] = c->lowpoint[0];
            //p0[1] = c->lowpoint[1];
            //p1[0] = c->highpoint[0];
            //p1[1] = c->lowpoint[1];
            //p2[0] = c->highpoint[0];
            //p2[1] = c->highpoint[1];
            //p3[0] = c->lowpoint[0];
            //p3[1] = c->highpoint[1];
            //real etav0 = 0.0; 
            //real etav1 = 0.0;
            //real etav2 = 0.0;
            //real etav3 = 0.0;
            //etav0 = compute_value_at_point(ns->ed.vevv.sdVisc, p0, p0, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            //etav1 = compute_value_at_point(ns->ed.vevv.sdVisc, p1, p1, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            //etav2 = compute_value_at_point(ns->ed.vevv.sdVisc, p2, p2, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            //etav3 = compute_value_at_point(ns->ed.vevv.sdVisc, p3, p3, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            real val = compute_value_at_point(ns->ed.vevv.sdVisc, ccenter, ccenter, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            //real val = dp_get_value(ns->ed.vevv.dpvisc, cgid);
            fprintf(f, "%e\n", val);
            //fprintf(f, "%e\n", etav0);
            //fprintf(f, "%e\n", etav1);
            //fprintf(f, "%e\n", etav2);
            //fprintf(f, "%e\n", etav3);
        }
        higcit_destroy(it);

        //Non-isothermal flows
        if (ns->contr.modelflowtype == 5){
            real T0 = ns->ed.nif.par.T0;
            real T1 = ns->ed.nif.par.T1;
            //Printing the temperature
            sim_domain *sdna =  psd_get_local_domain(ns->ed.nif.psdT);
            fprintf(f, "\nSCALARS Temperature FLOAT\nLOOKUP_TABLE default\n");
            for (it = sd_get_domain_celliterator(sdna); !higcit_isfinished(it); higcit_nextcell(it))
            {
                hig_cell *c = higcit_getcell(it);
                uniqueid id = hig_get_cid(c);
                int cgid = mp_lookup(m, id);
                Point cdelta, ccenter, clowpoint, chightpoint;
                hig_get_delta(c, cdelta);
                //Point ccenter;
                hig_get_center(c, ccenter);
                //Calculation of the dimensionless temperature val0 = (T-T0)/(T1-T0)
                real val0 = compute_value_at_point(ns->ed.nif.sdT, ccenter, ccenter, 1.0, ns->ed.nif.dpT, ns->ed.stn);
                //Printing the dimensional temperature T = (T1-T0)*val0 + T0
                real val = (T1-T0)*val0 + T0;
                fprintf(f, "%e\n", val);
            }
            higcit_destroy(it);

        }

        break;

    case 7:
        //sim_domain *sdte =  psd_get_local_domain(ns->ed.psdED);
        fprintf(f, "\nTENSORS stress FLOAT\n");

        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
            hig_cell *c = higcit_getcell(it);
            Point cdelta, ccenter, clowpoint, chightpoint;
            hig_get_delta(c, cdelta);
            hig_get_center(c, ccenter);
            // Pontos onde será interpolada a velocidade
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

            real taup0[DIM + 1][DIM + 1], taup1[DIM + 1][DIM + 1], taup2[DIM + 1][DIM + 1], taup3[DIM + 1][DIM + 1];

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    taup0[i][j] = 0.0;
                    taup1[i][j] = 0.0;
                    taup2[i][j] = 0.0;
                    taup3[i][j] = 0.0;
                }
            }

            //  double Sxx0 = 0, Sxy0 = 0, Sxz0 = 0, Syx0 = 0, Syy0 = 0, Syz0 = 0, Szx0 = 0, Szy0 = 0, Szz0 = 0;
            //  double Sxx1 = 0, Sxy1 = 0, Sxz1 = 0, Syx1 = 0, Syy1 = 0, Syz1 = 0, Szx1 = 0, Szy1 = 0, Szz1 = 0;
            //  double Sxx2 = 0, Sxy2 = 0, Sxz2 = 0, Syx2 = 0, Syy2 = 0, Syz2 = 0, Szx2 = 0, Szy2 = 0, Szz2 = 0;
            //  double Sxx3 = 0, Sxy3 = 0, Sxz3 = 0, Syx3 = 0, Syy3 = 0, Syz3 = 0, Szx3 = 0, Szy3 = 0, Szz3 = 0;

            sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);

            real Dp0[DIM][DIM], Dp1[DIM][DIM], Dp2[DIM][DIM], Dp3[DIM][DIM];
            //for (int i = 0; i < DIM; i++)
            //{
            //    for (int j = 0; j < DIM; j++)
            //    {
            //        // Get Du
            //        Dp0[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p0, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
            //        Dp1[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p1, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
            //        Dp2[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p2, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
            //        Dp3[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p3, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
            //    }
            //}
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    // Get Du
                    Dp0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.vesb.dpD[i][j], ns->ed.stn);
                    Dp1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.vesb.dpD[i][j], ns->ed.stn);
                    Dp2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.vesb.dpD[i][j], ns->ed.stn);
                    Dp3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.vesb.dpD[i][j], ns->ed.stn);
                }
            }

            //Get the viscosity value
            //real eta0 = compute_value_at_point(ns->ed.vevv.sdVisc, ccenter, p0, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            //real eta1 = compute_value_at_point(ns->ed.vevv.sdVisc, ccenter, p1, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            //real eta2 = compute_value_at_point(ns->ed.vevv.sdVisc, ccenter, p2, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            //real eta3 = compute_value_at_point(ns->ed.vevv.sdVisc, ccenter, p3, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);

            //real eta0 = compute_value_at_point(ns->ed.vevv.sdVisc, p0, p0, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            //real eta1 = compute_value_at_point(ns->ed.vevv.sdVisc, p1, p1, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            //real eta2 = compute_value_at_point(ns->ed.vevv.sdVisc, p2, p2, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            //real eta3 = compute_value_at_point(ns->ed.vevv.sdVisc, p3, p3, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            

            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    taup0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.vesb.dpS[i][j], ns->ed.stn);
                    taup0[i][j] += 2.0 * (ns->ed.vesb.par.De) * 0.5 * (Dp0[i][j] + Dp0[j][i]) / Re;
                    //taup0[i][j] += 2.0* (1.0- ns->ed.vevv.par.beta)* eta0 * 0.5 * (Dp0[i][j] + Dp0[j][i])/Re -2.0*(ns->ed.vevv.par.beta)/Re* 0.5 * (Dp0[i][j] + Dp0[j][i]);
                    taup1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.vesb.dpS[i][j], ns->ed.stn);
                    taup1[i][j] += 2.0 * (ns->ed.vesb.par.De) * 0.5 * (Dp1[i][j] + Dp1[j][i]) / Re;
                    //taup1[i][j] += 2.0* (1.0- ns->ed.vevv.par.beta)* eta1 * 0.5 * (Dp1[i][j] + Dp1[j][i])/Re -2.0*(ns->ed.vevv.par.beta)/Re* 0.5 * (Dp1[i][j] + Dp1[j][i]);
                    taup2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.vesb.dpS[i][j], ns->ed.stn);
                    taup2[i][j] += 2.0 * (ns->ed.vesb.par.De) * 0.5 * (Dp2[i][j] + Dp2[j][i]) / Re;
                    //taup2[i][j] += 2.0* (1.0- ns->ed.vevv.par.beta)* eta2 * 0.5 * (Dp2[i][j] + Dp2[j][i])/Re -2.0*(ns->ed.vevv.par.beta)/Re* 0.5 * (Dp2[i][j] + Dp2[j][i]);
                    taup3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.vesb.dpS[i][j], ns->ed.stn);
                    taup3[i][j] += 2.0 * (ns->ed.vesb.par.De) * 0.5 * (Dp3[i][j] + Dp3[j][i]) / Re;
                    //taup3[i][j] += 2.0* (1.0- ns->ed.vevv.par.beta)* eta3 * 0.5 * (Dp3[i][j] + Dp3[j][i])/Re -2.0*(ns->ed.vevv.par.beta)/Re* 0.5 * (Dp3[i][j] + Dp3[j][i]);
                }
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup0[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup1[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup2[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup3[i][j]);
                }
                fprintf(f, "\n");
            }

            /*  
                        Sxx0 = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.ve.dpS[0][0], ns->ed.stn);
                        Sxx0 += 2.0*(1-beta)*Dp0[0][0]/Re; 
                        Sxy0 = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.ve.dpS[0][1], ns->ed.stn);
                        Sxy0 += 2.0*(1-beta)*Dp0[1][0]/Re;
                        Syx0 = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.ve.dpS[1][0], ns->ed.stn);
                        Syx0 += 2.0*(1-beta)*Dp0[0][1]/Re;
                        Syy0 = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.ve.dpS[1][1], ns->ed.stn);
                        Syy0 += 2.0*(1-beta)*Dp0[1][1]/Re;
                        
                        Sxx1 = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.ve.dpS[0][0], ns->ed.stn);
                        Sxx1 += 2.0*(1-beta)*Dp1[0][0]/Re; 
                        Sxy1 = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.ve.dpS[0][1], ns->ed.stn);
                        Sxy1 += 2.0*(1-beta)*Dp1[1][0]/Re;
                        Syx1 = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.ve.dpS[1][0], ns->ed.stn);
                        Syx1 += 2.0*(1-beta)*Dp1[0][1]/Re;
                        Syy1 = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.ve.dpS[1][1], ns->ed.stn);
                        Syy1 += 2.0*(1-beta)*Dp1[1][1]/Re;

                        Sxx2 = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.ve.dpS[0][0], ns->ed.stn);
                        Sxx2 += 2.0*(1-beta)*Dp2[0][0]/Re;
                        Sxy2 = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.ve.dpS[0][1], ns->ed.stn);
                        Sxy2 += 2.0*(1-beta)*Dp2[1][0]/Re;
                        Syx2 = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.ve.dpS[1][0], ns->ed.stn);
                        Syx2 += 2.0*(1-beta)*Dp2[0][1]/Re;
                        Syy2 = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.ve.dpS[1][1], ns->ed.stn);
                        Syy2 += 2.0*(1-beta)*Dp2[1][1]/Re;
                        
                        Sxx3 = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.ve.dpS[0][0], ns->ed.stn);
                        Sxx3 += 2.0*(1-beta)*Dp3[0][0]/Re; 
                        Sxy3 = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.ve.dpS[0][1], ns->ed.stn);
                        Sxy3 += 2.0*(1-beta)*Dp3[1][0]/Re;
                        Syx3 = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.ve.dpS[1][0], ns->ed.stn);
                        Syx3 += 2.0*(1-beta)*Dp3[0][1]/Re;
                        Syy3 = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.ve.dpS[1][1], ns->ed.stn);
                        Syy3 += 2.0*(1-beta)*Dp3[1][1]/Re;
    
                        fprintf(f, "%e %e %e\n", Sxx0, Sxy0, Sxz0);
                        fprintf(f, "%e %e %e\n", Syx0, Syy0, Syz0);
                        fprintf(f, "%e %e %e\n", Szx0, Szy0, Szz0);
                        fprintf(f, "%e %e %e\n", Sxx1, Sxy1, Sxz1);
                        fprintf(f, "%e %e %e\n", Syx1, Syy1, Syz1);
                        fprintf(f, "%e %e %e\n", Szx1, Szy1, Szz1);
                        fprintf(f, "%e %e %e\n", Sxx2, Sxy2, Sxz2);
                        fprintf(f, "%e %e %e\n", Syx2, Syy2, Syz2);
                        fprintf(f, "%e %e %e\n", Szx2, Szy2, Szz2);
                        fprintf(f, "%e %e %e\n", Sxx3, Sxy3, Sxz3);
                        fprintf(f, "%e %e %e\n", Syx3, Syy3, Syz3);
                        fprintf(f, "%e %e %e\n", Szx3, Szy3, Szz3);   */
        }
        higcit_destroy(it);

        fprintf(f, "\nTENSORS A FLOAT\n");

        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
            hig_cell *c = higcit_getcell(it);
            Point cdelta, ccenter, clowpoint, chightpoint;
            hig_get_delta(c, cdelta);
            hig_get_center(c, ccenter);
            // Pontos onde será interpolada a velocidade
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

            real taup0[DIM + 1][DIM + 1], taup1[DIM + 1][DIM + 1], taup2[DIM + 1][DIM + 1], taup3[DIM + 1][DIM + 1];

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    taup0[i][j] = 0.0;
                    taup1[i][j] = 0.0;
                    taup2[i][j] = 0.0;
                    taup3[i][j] = 0.0;
                }
            }

            sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);

            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    taup0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.vesb.dpA[i][j], ns->ed.stn);
                    taup1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.vesb.dpA[i][j], ns->ed.stn);
                    taup2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.vesb.dpA[i][j], ns->ed.stn);
                    taup3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.vesb.dpA[i][j], ns->ed.stn);
                }
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup0[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup1[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup2[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup3[i][j]);
                }
                fprintf(f, "\n");
            }
        }
        higcit_destroy(it);


        fprintf(f, "\nTENSORS B FLOAT\n");

        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
            hig_cell *c = higcit_getcell(it);
            Point cdelta, ccenter, clowpoint, chightpoint;
            hig_get_delta(c, cdelta);
            hig_get_center(c, ccenter);
            // Pontos onde será interpolada a velocidade
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

            real taup0[DIM + 1][DIM + 1], taup1[DIM + 1][DIM + 1], taup2[DIM + 1][DIM + 1], taup3[DIM + 1][DIM + 1];

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    taup0[i][j] = 0.0;
                    taup1[i][j] = 0.0;
                    taup2[i][j] = 0.0;
                    taup3[i][j] = 0.0;
                }
            }

            sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);

            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    taup0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.vesb.dpB[i][j], ns->ed.stn);
                    taup1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.vesb.dpB[i][j], ns->ed.stn);
                    taup2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.vesb.dpB[i][j], ns->ed.stn);
                    taup3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.vesb.dpB[i][j], ns->ed.stn);
                }
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup0[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup1[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup2[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup3[i][j]);
                }
                fprintf(f, "\n");
            }
        }
        higcit_destroy(it);

        //fprintf(f, "\nSCALARS p FLOAT\nLOOKUP_TABLE default\n");
        fprintf(f, "\nCELL_DATA %ld\nSCALARS p FLOAT\nLOOKUP_TABLE default\n", numcells);
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
            hig_cell *c = higcit_getcell(it);
            uniqueid id = hig_get_cid(c);
            int cgid = mp_lookup(m, id);
            real val = dp_get_value(ns->dpp, cgid);
            fprintf(f, "%e\n", val);
        }
        higcit_destroy(it);

        if (ns->contr.rheotype == VCM){
            //Printing the density number nA
            sim_domain *sdna =  psd_get_local_domain(ns->ed.vesb.psdSBnA);
            fprintf(f, "\nSCALARS nA FLOAT\nLOOKUP_TABLE default\n");
            for (it = sd_get_domain_celliterator(sdna); !higcit_isfinished(it); higcit_nextcell(it))
            {
                hig_cell *c = higcit_getcell(it);
                uniqueid id = hig_get_cid(c);
                int cgid = mp_lookup(m, id);
                Point cdelta, ccenter, clowpoint, chightpoint;
                hig_get_delta(c, cdelta);
                //Point ccenter;
                hig_get_center(c, ccenter);
                real val = compute_value_at_point(ns->ed.vesb.sdSBnA, ccenter, ccenter, 1.0, ns->ed.vesb.dpnA, ns->ed.stn);
                //real val = dp_get_value(ns->ed.vevv.dpvisc, cgid);
                fprintf(f, "%e\n", val);
            }
            higcit_destroy(it);

            //Printing the density number nB
            sim_domain *sdnb =  psd_get_local_domain(ns->ed.vesb.psdSBnB);
            fprintf(f, "\nSCALARS nB FLOAT\nLOOKUP_TABLE default\n");
            for (it = sd_get_domain_celliterator(sdnb); !higcit_isfinished(it); higcit_nextcell(it))
            {
                hig_cell *c = higcit_getcell(it);
                uniqueid id = hig_get_cid(c);
                int cgid = mp_lookup(m, id);
                Point cdelta, ccenter, clowpoint, chightpoint;
                hig_get_delta(c, cdelta);
                //Point ccenter;
                hig_get_center(c, ccenter);
                real val = compute_value_at_point(ns->ed.vesb.sdSBnB, ccenter, ccenter, 1.0, ns->ed.vesb.dpnB, ns->ed.stn);
                //real val = dp_get_value(ns->ed.vevv.dpvisc, cgid);
                fprintf(f, "%e\n", val);
            }
            higcit_destroy(it);
        }
        break;

    case 8:
        // Elastoviscoplastic

        fprintf(f, "\nTENSORS stress FLOAT\n");

        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
            hig_cell *c = higcit_getcell(it);
            Point cdelta, ccenter, clowpoint, chightpoint;
            hig_get_delta(c, cdelta);
            hig_get_center(c, ccenter);
            // Pontos onde será interpolada a velocidade
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

            real taup0[DIM + 1][DIM + 1], taup1[DIM + 1][DIM + 1], taup2[DIM + 1][DIM + 1], taup3[DIM + 1][DIM + 1];

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    taup0[i][j] = 0.0;
                    taup1[i][j] = 0.0;
                    taup2[i][j] = 0.0;
                    taup3[i][j] = 0.0;
                }
            }

            sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);

            real Dp0[DIM][DIM], Dp1[DIM][DIM], Dp2[DIM][DIM], Dp3[DIM][DIM];
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    // Get Du
                    Dp0[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p0, 1.0, ns->ed.vepl.dpD[i][j], ns->ed.stn);
                    Dp1[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p1, 1.0, ns->ed.vepl.dpD[i][j], ns->ed.stn);
                    Dp2[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p2, 1.0, ns->ed.vepl.dpD[i][j], ns->ed.stn);
                    Dp3[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p3, 1.0, ns->ed.vepl.dpD[i][j], ns->ed.stn);
                }
            }

            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    taup0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.vepl.dpS[i][j], ns->ed.stn);
                    taup0[i][j] += 2.0 * (1.0 - ns->ed.vepl.par.beta) * 0.5 * (Dp0[i][j] + Dp0[j][i]) / Re;
                    taup1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.vepl.dpS[i][j], ns->ed.stn);
                    taup1[i][j] += 2.0 * (1.0 - ns->ed.vepl.par.beta) * 0.5 * (Dp1[i][j] + Dp1[j][i]) / Re;
                    taup2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.vepl.dpS[i][j], ns->ed.stn);
                    taup2[i][j] += 2.0 * (1.0- ns->ed.vepl.par.beta) * 0.5 * (Dp2[i][j] + Dp2[j][i]) / Re;
                    taup3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.vepl.dpS[i][j], ns->ed.stn);
                    taup3[i][j] += 2.0 * (1.0 - ns->ed.vepl.par.beta) * 0.5 * (Dp3[i][j] + Dp3[j][i]) / Re;
                }
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup0[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup1[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup2[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup3[i][j]);
                }
                fprintf(f, "\n");
            }

        }
        higcit_destroy(it);

        //fprintf(f, "\nSCALARS p FLOAT\nLOOKUP_TABLE default\n");
        fprintf(f, "\nCELL_DATA %ld\nSCALARS p FLOAT\nLOOKUP_TABLE default\n", numcells);
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
            hig_cell *c = higcit_getcell(it);
            uniqueid id = hig_get_cid(c);
            int cgid = mp_lookup(m, id);
            real val = dp_get_value(ns->dpp, cgid);
            fprintf(f, "%e\n", val);
        }
        higcit_destroy(it);
        break;

    case 9:
        // Shear-thickening suspensions
        fprintf(f, "\nTENSORS stress FLOAT\n");

        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
            hig_cell *c = higcit_getcell(it);
            Point cdelta, ccenter, clowpoint, chightpoint;
            hig_get_delta(c, cdelta);
            hig_get_center(c, ccenter);
            // Pontos onde será interpolada a velocidade
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

            real taup0[DIM + 1][DIM + 1], taup1[DIM + 1][DIM + 1], taup2[DIM + 1][DIM + 1], taup3[DIM + 1][DIM + 1];

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    taup0[i][j] = 0.0;
                    taup1[i][j] = 0.0;
                    taup2[i][j] = 0.0;
                    taup3[i][j] = 0.0;
                }
            }

            sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);

            real Dp0[DIM][DIM], Dp1[DIM][DIM], Dp2[DIM][DIM], Dp3[DIM][DIM];
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    // Get Du
                    Dp0[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p0, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
                    Dp1[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p1, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
                    Dp2[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p2, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
                    Dp3[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, p3, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
                }
            }

            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    //Tau = S + 2*(eta0)*D
                    taup0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.stsp.dpS[i][j], ns->ed.stn);
                    taup0[i][j] += 2.0 * (ns->ed.stsp.par.eta0) * 0.5 * (Dp0[i][j] + Dp0[j][i]) / Re;
                    taup1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.stsp.dpS[i][j], ns->ed.stn);
                    taup1[i][j] += 2.0 * (ns->ed.stsp.par.eta0) * 0.5 * (Dp1[i][j] + Dp1[j][i]) / Re;
                    taup2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.stsp.dpS[i][j], ns->ed.stn);
                    taup2[i][j] += 2.0 * (ns->ed.stsp.par.eta0) * 0.5 * (Dp2[i][j] + Dp2[j][i]) / Re;
                    taup3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.stsp.dpS[i][j], ns->ed.stn);
                    taup3[i][j] += 2.0 * (ns->ed.stsp.par.eta0) * 0.5 * (Dp3[i][j] + Dp3[j][i]) / Re;
                }
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup0[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup1[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup2[i][j]);
                }
                fprintf(f, "\n");
            }

            for (int i = 0; i <= DIM; i++)
            {
                for (int j = 0; j <= DIM; j++)
                {
                    fprintf(f, "%e ", taup3[i][j]);
                }
                fprintf(f, "\n");
            }

        }
        higcit_destroy(it);

        //Pressure
        fprintf(f, "\nCELL_DATA %ld\nSCALARS p FLOAT\nLOOKUP_TABLE default\n", numcells);
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
            hig_cell *c = higcit_getcell(it);
            uniqueid id = hig_get_cid(c);
            int cgid = mp_lookup(m, id);
            real val = dp_get_value(ns->dpp, cgid);
            fprintf(f, "%e\n", val);
        }
        higcit_destroy(it);

        break;


    }
    fclose(f);
}

// Print the VTK file for visualize 3D
void higflow_print_vtk3D(higflow_solver *ns, int rank)
{
    // Open the VTK file
    char vtkname[1024];
    snprintf(vtkname, sizeof vtkname, "%s_%d-%d.vtk", ns->par.nameprint, rank, ns->par.frame);
    FILE *f = fopen(vtkname, "w");
    if (f == NULL)
    {
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
    for (; !higcit_isfinished(it_t); higcit_nextcell(it_t))
    {
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
    for (int i = 0; i < numleafs; i++)
    {
        fprintf(f, "%d %d %d %d %d %d %d %d %d\n", 8, 8 * i + 0, 8 * i + 1, 8 * i + 2, 8 * i + 3, 8 * i + 4, 8 * i + 5, 8 * i + 6, 8 * i + 7); //x == 0
    }
    fprintf(f, "CELL_TYPES %ld\n", numleafs);
    for (int i = 0; i < numleafs; i++)
    {
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
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
    {
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
        for (int dim = 0; dim < DIM; dim++)
        {
            // int numpts = 0;
            // int cont = 0;
            sfdu2[dim] = psfd_get_local_domain(ns->psfdu[dim]);
            mp_mapper *m = sfd_get_domain_mapper(sfdu2[dim]);
            real value0 = compute_facet_value_at_point(sfdu2[dim], ccenter, p0, 1.0, ns->dpu[dim], ns->stn);
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
    {
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
void higflow_print_vtk3D_viscoelastic(higflow_solver *ns, int rank)
{
    // Open the VTK file
    char vtkname[1024];
    snprintf(vtkname, sizeof vtkname, "%s_%d-%d.vtk", ns->par.nameprint, rank, ns->par.frame);
    FILE *f = fopen(vtkname, "w");
    if (f == NULL)
    {
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
    for (; !higcit_isfinished(it_t); higcit_nextcell(it_t))
    {
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
    for (int i = 0; i < numleafs; i++)
    {
        fprintf(f, "%d %d %d %d %d %d %d %d %d\n", 8, 8 * i + 0, 8 * i + 1, 8 * i + 2, 8 * i + 3, 8 * i + 4, 8 * i + 5, 8 * i + 6, 8 * i + 7); //x == 0
    }
    fprintf(f, "CELL_TYPES %ld\n", numleafs);
    for (int i = 0; i < numleafs; i++)
    {
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
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
    {
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
        for (int dim = 0; dim < DIM; dim++)
        {
            // int numpts = 0;
            // int cont = 0;
            sfdu2[dim] = psfd_get_local_domain(ns->psfdu[dim]);
            mp_mapper *m = sfd_get_domain_mapper(sfdu2[dim]);
            real value0 = compute_facet_value_at_point(sfdu2[dim], ccenter, p0, 1.0, ns->dpu[dim], ns->stn);
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
    for (it = sd_get_domain_celliterator(ns->ed.sdED); !higcit_isfinished(it); higcit_nextcell(it))
    {
        hig_cell *c = higcit_getcell(it);
        uniqueid id = hig_get_cid(c);
        int cgid = mp_lookup(m, id);
        real val = dp_get_value(ns->ed.ve.dpS[0][0], cgid);
        fprintf(f, "%e\n", val);
    }
    higcit_destroy(it);
    fclose(f);
}

// Print the VTK file for vilusalize 3D viscoelastic flows with variable viscosity
void higflow_print_vtk3D_viscoelastic_variable_viscosity(higflow_solver *ns, int rank)
{
    // Open the VTK file
    real Re = ns->par.Re;
    char vtkname[1024];
    snprintf(vtkname, sizeof vtkname, "%s_%d-%d.vtk", ns->par.nameprint, rank, ns->par.frame);
    FILE *f = fopen(vtkname, "w");
    if (f == NULL)
    {
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
    for (; !higcit_isfinished(it_t); higcit_nextcell(it_t))
    {
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
    for (int i = 0; i < numleafs; i++)
    {
        fprintf(f, "%d %d %d %d %d %d %d %d %d\n", 8, 8 * i + 0, 8 * i + 1, 8 * i + 2, 8 * i + 3, 8 * i + 4, 8 * i + 5, 8 * i + 6, 8 * i + 7); //x == 0
    }
    fprintf(f, "CELL_TYPES %ld\n", numleafs);
    for (int i = 0; i < numleafs; i++)
    {
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
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
    {
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
        for (int dim = 0; dim < DIM; dim++)
        {
            // int numpts = 0;
            // int cont = 0;
            sfdu2[dim] = psfd_get_local_domain(ns->psfdu[dim]);
            mp_mapper *m = sfd_get_domain_mapper(sfdu2[dim]);
            real value0 = compute_facet_value_at_point(sfdu2[dim], ccenter, p0, 1.0, ns->dpu[dim], ns->stn);
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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

    //Printing the tensors
    fprintf(f, "\nTENSORS stress FLOAT\n");
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
    {
        hig_cell *c = higcit_getcell(it);
        Point cdelta, ccenter, clowpoint, chightpoint;
        hig_get_delta(c, cdelta);
        hig_get_center(c, ccenter);
        // Points where the tensors will be interpolated
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
        real taup0[DIM][DIM], taup1[DIM][DIM], taup2[DIM][DIM], taup3[DIM][DIM];
        real taup4[DIM][DIM], taup5[DIM][DIM], taup6[DIM][DIM], taup7[DIM][DIM];

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                taup0[i][j] = 0.0;
                taup1[i][j] = 0.0;
                taup2[i][j] = 0.0;
                taup3[i][j] = 0.0;
                taup4[i][j] = 0.0;
                taup5[i][j] = 0.0;
                taup6[i][j] = 0.0;
                taup7[i][j] = 0.0;
            }
        }
        
        real Dp0[DIM][DIM], Dp1[DIM][DIM], Dp2[DIM][DIM], Dp3[DIM][DIM];
        real Dp4[DIM][DIM], Dp5[DIM][DIM], Dp6[DIM][DIM], Dp7[DIM][DIM];
        
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                // Get Du
                Dp0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
                Dp1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
                Dp2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
                Dp3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
                Dp4[i][j] = compute_value_at_point(ns->ed.sdED, p4, p4, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
                Dp5[i][j] = compute_value_at_point(ns->ed.sdED, p5, p5, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
                Dp6[i][j] = compute_value_at_point(ns->ed.sdED, p6, p6, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
                Dp7[i][j] = compute_value_at_point(ns->ed.sdED, p7, p7, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
            }
        }
        
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                taup0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.vevv.dpS[i][j], ns->ed.stn);
                taup0[i][j] += 2.0 * (1.0 - ns->ed.vevv.par.beta) * 0.5 * (Dp0[i][j] + Dp0[j][i]) / Re;
                taup1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.vevv.dpS[i][j], ns->ed.stn);
                taup1[i][j] += 2.0 * (1.0 - ns->ed.vevv.par.beta) * 0.5 * (Dp1[i][j] + Dp1[j][i]) / Re;
                taup2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.vevv.dpS[i][j], ns->ed.stn);
                taup2[i][j] += 2.0 * (1.0 - ns->ed.vevv.par.beta) * 0.5 * (Dp2[i][j] + Dp2[j][i]) / Re;
                taup3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.vevv.dpS[i][j], ns->ed.stn);
                taup3[i][j] += 2.0 * (1.0 - ns->ed.vevv.par.beta) * 0.5 * (Dp3[i][j] + Dp3[j][i]) / Re;
                taup4[i][j] = compute_value_at_point(ns->ed.sdED, p4, p4, 1.0, ns->ed.vevv.dpS[i][j], ns->ed.stn);
                taup4[i][j] += 2.0 * (1.0 - ns->ed.vevv.par.beta) * 0.5 * (Dp4[i][j] + Dp4[j][i]) / Re;
                taup5[i][j] = compute_value_at_point(ns->ed.sdED, p5, p5, 1.0, ns->ed.vevv.dpS[i][j], ns->ed.stn);
                taup5[i][j] += 2.0 * (1.0 - ns->ed.vevv.par.beta) * 0.5 * (Dp5[i][j] + Dp5[j][i]) / Re;
                taup6[i][j] = compute_value_at_point(ns->ed.sdED, p6, p6, 1.0, ns->ed.vevv.dpS[i][j], ns->ed.stn);
                taup6[i][j] += 2.0 * (1.0 - ns->ed.vevv.par.beta) * 0.5 * (Dp6[i][j] + Dp6[j][i]) / Re;
                taup7[i][j] = compute_value_at_point(ns->ed.sdED, p7, p7, 1.0, ns->ed.vevv.dpS[i][j], ns->ed.stn);
                taup7[i][j] += 2.0 * (1.0 - ns->ed.vevv.par.beta) * 0.5 * (Dp7[i][j] + Dp7[j][i]) / Re;
            }
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup0[i][j]);
            }
            fprintf(f, "\n");
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup1[i][j]);
            }
            fprintf(f, "\n");
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup2[i][j]);
            }
            fprintf(f, "\n");
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup3[i][j]);
            }
            fprintf(f, "\n");
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup4[i][j]);
            }
            fprintf(f, "\n");
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup5[i][j]);
            }
            fprintf(f, "\n");
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup6[i][j]);
            }
            fprintf(f, "\n");
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup7[i][j]);
            }
            fprintf(f, "\n");
        }

    }
    higcit_destroy(it);

    // Pressure
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
    
    //Printing the viscosity
    //sim_domain *sdvisc =  psd_get_local_domain(ns->ed.vevv.psdVisc);
    fprintf(f, "\nSCALARS viscosity FLOAT\nLOOKUP_TABLE default\n");
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
    {
        hig_cell *c = higcit_getcell(it);
        uniqueid id = hig_get_cid(c);
        int cgid = mp_lookup(m, id);
        Point cdelta, ccenter, clowpoint, chightpoint;
        hig_get_delta(c, cdelta);
        hig_get_center(c, ccenter);
        real val = compute_value_at_point(ns->ed.vevv.sdVisc, ccenter, ccenter, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
        //real val = dp_get_value(ns->ed.vevv.dpvisc, cgid);
        fprintf(f, "%e\n", val);
    }
    higcit_destroy(it);

    fclose(f);
}

// Print the VTK file for vilusalize 3D viscoelastic flows with shear-banding
void higflow_print_vtk3D_viscoelastic_shear_banding(higflow_solver *ns, int rank)
{
    // Open the VTK file
    char vtkname[1024];
    snprintf(vtkname, sizeof vtkname, "%s_%d-%d.vtk", ns->par.nameprint, rank, ns->par.frame);
    FILE *f = fopen(vtkname, "w");
    if (f == NULL)
    {
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
    for (; !higcit_isfinished(it_t); higcit_nextcell(it_t))
    {
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
    for (int i = 0; i < numleafs; i++)
    {
        fprintf(f, "%d %d %d %d %d %d %d %d %d\n", 8, 8 * i + 0, 8 * i + 1, 8 * i + 2, 8 * i + 3, 8 * i + 4, 8 * i + 5, 8 * i + 6, 8 * i + 7); //x == 0
    }
    fprintf(f, "CELL_TYPES %ld\n", numleafs);
    for (int i = 0; i < numleafs; i++)
    {
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
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
    {
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
        for (int dim = 0; dim < DIM; dim++)
        {
            // int numpts = 0;
            // int cont = 0;
            sfdu2[dim] = psfd_get_local_domain(ns->psfdu[dim]);
            mp_mapper *m = sfd_get_domain_mapper(sfdu2[dim]);
            real value0 = compute_facet_value_at_point(sfdu2[dim], ccenter, p0, 1.0, ns->dpu[dim], ns->stn);
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
    for (it = sd_get_domain_celliterator(ns->ed.sdED); !higcit_isfinished(it); higcit_nextcell(it))
    {
        hig_cell *c = higcit_getcell(it);
        uniqueid id = hig_get_cid(c);
        int cgid = mp_lookup(m, id);
        real val = dp_get_value(ns->ed.vesb.dpS[0][0], cgid);
        fprintf(f, "%e\n", val);
    }
    higcit_destroy(it);
    fclose(f);
}


// Print the VTK file for vilusalize 3D
void higflow_print_vtk3D_elastoviscoplastic(higflow_solver *ns, int rank)
{
    // Open the VTK file
    char vtkname[1024];
    snprintf(vtkname, sizeof vtkname, "%s_%d-%d.vtk", ns->par.nameprint, rank, ns->par.frame);
    FILE *f = fopen(vtkname, "w");
    if (f == NULL)
    {
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
    for (; !higcit_isfinished(it_t); higcit_nextcell(it_t))
    {
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
    for (int i = 0; i < numleafs; i++)
    {
        fprintf(f, "%d %d %d %d %d %d %d %d %d\n", 8, 8 * i + 0, 8 * i + 1, 8 * i + 2, 8 * i + 3, 8 * i + 4, 8 * i + 5, 8 * i + 6, 8 * i + 7); //x == 0
    }
    fprintf(f, "CELL_TYPES %ld\n", numleafs);
    for (int i = 0; i < numleafs; i++)
    {
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
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
    {
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
        for (int dim = 0; dim < DIM; dim++)
        {
            // int numpts = 0;
            // int cont = 0;
            sfdu2[dim] = psfd_get_local_domain(ns->psfdu[dim]);
            mp_mapper *m = sfd_get_domain_mapper(sfdu2[dim]);
            real value0 = compute_facet_value_at_point(sfdu2[dim], ccenter, p0, 1.0, ns->dpu[dim], ns->stn);
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
    for (it = sd_get_domain_celliterator(ns->ed.sdED); !higcit_isfinished(it); higcit_nextcell(it))
    {
        hig_cell *c = higcit_getcell(it);
        uniqueid id = hig_get_cid(c);
        int cgid = mp_lookup(m, id);
        real val = dp_get_value(ns->ed.vepl.dpS[0][0], cgid);
        fprintf(f, "%e\n", val);
    }
    higcit_destroy(it);
    fclose(f);
}


// Print the VTK file for vilusalize 3D
void higflow_print_vtk3D_shear_thickening_suspensions(higflow_solver *ns, int rank)
{
    // Open the VTK file
    char vtkname[1024];
    snprintf(vtkname, sizeof vtkname, "%s_%d-%d.vtk", ns->par.nameprint, rank, ns->par.frame);
    FILE *f = fopen(vtkname, "w");
    if (f == NULL)
    {
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
    for (; !higcit_isfinished(it_t); higcit_nextcell(it_t))
    {
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
    for (int i = 0; i < numleafs; i++)
    {
        fprintf(f, "%d %d %d %d %d %d %d %d %d\n", 8, 8 * i + 0, 8 * i + 1, 8 * i + 2, 8 * i + 3, 8 * i + 4, 8 * i + 5, 8 * i + 6, 8 * i + 7); //x == 0
    }
    fprintf(f, "CELL_TYPES %ld\n", numleafs);
    for (int i = 0; i < numleafs; i++)
    {
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
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
    {
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
        for (int dim = 0; dim < DIM; dim++)
        {
            // int numpts = 0;
            // int cont = 0;
            sfdu2[dim] = psfd_get_local_domain(ns->psfdu[dim]);
            mp_mapper *m = sfd_get_domain_mapper(sfdu2[dim]);
            real value0 = compute_facet_value_at_point(sfdu2[dim], ccenter, p0, 1.0, ns->dpu[dim], ns->stn);
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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
            switch (dim)
            {
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

    //Printing the tensors
    fprintf(f, "\nTENSORS stress FLOAT\n");
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
    {
        hig_cell *c = higcit_getcell(it);
        Point cdelta, ccenter, clowpoint, chightpoint;
        hig_get_delta(c, cdelta);
        hig_get_center(c, ccenter);
        // Points where the tensors will be interpolated
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
        real taup0[DIM][DIM], taup1[DIM][DIM], taup2[DIM][DIM], taup3[DIM][DIM];
        real taup4[DIM][DIM], taup5[DIM][DIM], taup6[DIM][DIM], taup7[DIM][DIM];

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                taup0[i][j] = 0.0;
                taup1[i][j] = 0.0;
                taup2[i][j] = 0.0;
                taup3[i][j] = 0.0;
                taup4[i][j] = 0.0;
                taup5[i][j] = 0.0;
                taup6[i][j] = 0.0;
                taup7[i][j] = 0.0;
            }
        }
        
        real Dp0[DIM][DIM], Dp1[DIM][DIM], Dp2[DIM][DIM], Dp3[DIM][DIM];
        real Dp4[DIM][DIM], Dp5[DIM][DIM], Dp6[DIM][DIM], Dp7[DIM][DIM];
        
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                // Get Du
                Dp0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
                Dp1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
                Dp2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
                Dp3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
                Dp4[i][j] = compute_value_at_point(ns->ed.sdED, p4, p4, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
                Dp5[i][j] = compute_value_at_point(ns->ed.sdED, p5, p5, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
                Dp6[i][j] = compute_value_at_point(ns->ed.sdED, p6, p6, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
                Dp7[i][j] = compute_value_at_point(ns->ed.sdED, p7, p7, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
            }
        }
        
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {

                taup0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.stsp.dpS[i][j], ns->ed.stn);
                taup0[i][j] += 2.0 * (ns->ed.stsp.par.eta0) * 0.5 * (Dp0[i][j] + Dp0[j][i]);
                taup1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.stsp.dpS[i][j], ns->ed.stn);
                taup1[i][j] += 2.0 * (ns->ed.stsp.par.eta0) * 0.5 * (Dp1[i][j] + Dp1[j][i]);
                taup2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.stsp.dpS[i][j], ns->ed.stn);
                taup2[i][j] += 2.0 * (ns->ed.stsp.par.eta0) * 0.5 * (Dp2[i][j] + Dp2[j][i]);
                taup3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.stsp.dpS[i][j], ns->ed.stn);
                taup3[i][j] += 2.0 * (ns->ed.stsp.par.eta0) * 0.5 * (Dp3[i][j] + Dp3[j][i]);
                taup4[i][j] = compute_value_at_point(ns->ed.sdED, p4, p4, 1.0, ns->ed.stsp.dpS[i][j], ns->ed.stn);
                taup4[i][j] += 2.0 * (ns->ed.stsp.par.eta0) * 0.5 * (Dp4[i][j] + Dp4[j][i]);
                taup5[i][j] = compute_value_at_point(ns->ed.sdED, p5, p5, 1.0, ns->ed.stsp.dpS[i][j], ns->ed.stn);
                taup5[i][j] += 2.0 * (ns->ed.stsp.par.eta0) * 0.5 * (Dp5[i][j] + Dp5[j][i]);
                taup6[i][j] = compute_value_at_point(ns->ed.sdED, p6, p6, 1.0, ns->ed.stsp.dpS[i][j], ns->ed.stn);
                taup6[i][j] += 2.0 * (ns->ed.stsp.par.eta0) * 0.5 * (Dp6[i][j] + Dp6[j][i]);
                taup7[i][j] = compute_value_at_point(ns->ed.sdED, p7, p7, 1.0, ns->ed.stsp.dpS[i][j], ns->ed.stn);
                taup7[i][j] += 2.0 * (ns->ed.stsp.par.eta0) * 0.5 * (Dp7[i][j] + Dp7[j][i]);
            }
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup0[i][j]);
            }
            fprintf(f, "\n");
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup1[i][j]);
            }
            fprintf(f, "\n");
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup2[i][j]);
            }
            fprintf(f, "\n");
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup3[i][j]);
            }
            fprintf(f, "\n");
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup4[i][j]);
            }
            fprintf(f, "\n");
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup5[i][j]);
            }
            fprintf(f, "\n");
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup6[i][j]);
            }
            fprintf(f, "\n");
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup7[i][j]);
            }
            fprintf(f, "\n");
        }

    }
    higcit_destroy(it);

    //Printing the tensors
    fprintf(f, "\nTENSORS A FLOAT\n");
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
    {
        hig_cell *c = higcit_getcell(it);
        Point cdelta, ccenter, clowpoint, chightpoint;
        hig_get_delta(c, cdelta);
        hig_get_center(c, ccenter);
        // Points where the tensors will be interpolated
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
        real taup0[DIM][DIM], taup1[DIM][DIM], taup2[DIM][DIM], taup3[DIM][DIM];
        real taup4[DIM][DIM], taup5[DIM][DIM], taup6[DIM][DIM], taup7[DIM][DIM];

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                taup0[i][j] = 0.0;
                taup1[i][j] = 0.0;
                taup2[i][j] = 0.0;
                taup3[i][j] = 0.0;
                taup4[i][j] = 0.0;
                taup5[i][j] = 0.0;
                taup6[i][j] = 0.0;
                taup7[i][j] = 0.0;
            }
        }
        
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {

                taup0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.stsp.dpA[i][j], ns->ed.stn);
                taup1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.stsp.dpA[i][j], ns->ed.stn);
                taup2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.stsp.dpA[i][j], ns->ed.stn);
                taup3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.stsp.dpA[i][j], ns->ed.stn);
                taup4[i][j] = compute_value_at_point(ns->ed.sdED, p4, p4, 1.0, ns->ed.stsp.dpA[i][j], ns->ed.stn);
                taup5[i][j] = compute_value_at_point(ns->ed.sdED, p5, p5, 1.0, ns->ed.stsp.dpA[i][j], ns->ed.stn);
                taup6[i][j] = compute_value_at_point(ns->ed.sdED, p6, p6, 1.0, ns->ed.stsp.dpA[i][j], ns->ed.stn);
                taup7[i][j] = compute_value_at_point(ns->ed.sdED, p7, p7, 1.0, ns->ed.stsp.dpA[i][j], ns->ed.stn);
            }
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup0[i][j]);
            }
            fprintf(f, "\n");
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup1[i][j]);
            }
            fprintf(f, "\n");
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup2[i][j]);
            }
            fprintf(f, "\n");
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup3[i][j]);
            }
            fprintf(f, "\n");
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup4[i][j]);
            }
            fprintf(f, "\n");
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup5[i][j]);
            }
            fprintf(f, "\n");
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup6[i][j]);
            }
            fprintf(f, "\n");
        }

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fprintf(f, "%e ", taup7[i][j]);
            }
            fprintf(f, "\n");
        }

    }
    higcit_destroy(it);


    // Pressure
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

    //only for the model that considers particle migration
    if (ns->ed.stsp.contr.model == GW_WC_IF) {
        //Printing the viscosity
        //sim_domain *sdvisc =  psd_get_local_domain(ns->ed.vevv.psdVisc);
        fprintf(f, "\nSCALARS VolFrac FLOAT\nLOOKUP_TABLE default\n");
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
        {
            hig_cell *c = higcit_getcell(it);
            uniqueid id = hig_get_cid(c);
            int cgid = mp_lookup(m, id);
            Point cdelta, ccenter, clowpoint, chightpoint;
            hig_get_delta(c, cdelta);
            hig_get_center(c, ccenter);
            real val = compute_value_at_point(ns->ed.stsp.sdphi, ccenter, ccenter, 1.0, ns->ed.stsp.dpphi, ns->ed.stn);
            real frac_tol = ns->ed.stsp.par.gdrms;
            real eps = 1.0e-4;
            real phij = ns->ed.stsp.par.phi;
            if(val > phij- frac_tol){
                val = phij- frac_tol;
            }
            if(val < eps){
                val = 0.0;
            }
            //real val = dp_get_value(ns->ed.vevv.dpvisc, cgid);
            fprintf(f, "%e\n", val);
        }
        higcit_destroy(it);
    }
    /*
    //Printing the viscosity
    //sim_domain *sdvisc =  psd_get_local_domain(ns->ed.vevv.psdVisc);
    //fprintf(f, "\nCELL_DATA %ld\nSCALARS viscosity FLOAT\nLOOKUP_TABLE default\n", numcells);
    fprintf(f, "\nSCALARS viscosity FLOAT\nLOOKUP_TABLE default\n");
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it))
    {
        hig_cell *c = higcit_getcell(it);
        Point cdelta, ccenter, clowpoint, chightpoint;
        hig_get_delta(c, cdelta);
        hig_get_center(c, ccenter);
        // Points where the tensors will be interpolated
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
        real taup0 = 0.0, taup1 = 0.0, taup2 = 0.0, taup3 = 0.0;
        real taup4 = 0.0, taup5 = 0.0, taup6 = 0.0, taup7 = 0.0;
        
        
        real Dp0[DIM][DIM], Dp1[DIM][DIM], Dp2[DIM][DIM], Dp3[DIM][DIM];
        real Dp4[DIM][DIM], Dp5[DIM][DIM], Dp6[DIM][DIM], Dp7[DIM][DIM];
        

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                // Get Du
                Dp0[i][j] = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
                Dp1[i][j] = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
                Dp2[i][j] = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
                Dp3[i][j] = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
                Dp4[i][j] = compute_value_at_point(ns->ed.sdED, p4, p4, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
                Dp5[i][j] = compute_value_at_point(ns->ed.sdED, p5, p5, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
                Dp6[i][j] = compute_value_at_point(ns->ed.sdED, p6, p6, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
                Dp7[i][j] = compute_value_at_point(ns->ed.sdED, p7, p7, 1.0, ns->ed.stsp.dpD[i][j], ns->ed.stn);
            }
        }
       
        // Calculate the shear rate q
        real D2p0[DIM][DIM], D2p1[DIM][DIM], D2p2[DIM][DIM], D2p3[DIM][DIM];
        real D2p4[DIM][DIM], D2p5[DIM][DIM], D2p6[DIM][DIM], D2p7[DIM][DIM];
        for (int dim = 0; dim < DIM; dim++) {
            for (int dim2 = 0; dim2 < DIM; dim2++) {
                D2p0[dim][dim2] = 0.0;
                for (int dim3 = 0; dim3 < DIM; dim3++) {
                    D2p0[dim][dim2] += Dp0[dim][dim3] * Dp0[dim3][dim2];
                }
            }
        }

        for (int dim = 0; dim < DIM; dim++) {
            for (int dim2 = 0; dim2 < DIM; dim2++) {
                D2p1[dim][dim2] = 0.0;
                for (int dim3 = 0; dim3 < DIM; dim3++) {
                    D2p1[dim][dim2] += Dp1[dim][dim3] * Dp1[dim3][dim2];
                }
            }
        }

        for (int dim = 0; dim < DIM; dim++) {
            for (int dim2 = 0; dim2 < DIM; dim2++) {
                D2p2[dim][dim2] = 0.0;
                for (int dim3 = 0; dim3 < DIM; dim3++) {
                    D2p2[dim][dim2] += Dp2[dim][dim3] * Dp2[dim3][dim2];
                }
            }
        }

        for (int dim = 0; dim < DIM; dim++) {
            for (int dim2 = 0; dim2 < DIM; dim2++) {
                D2p3[dim][dim2] = 0.0;
                for (int dim3 = 0; dim3 < DIM; dim3++) {
                    D2p3[dim][dim2] += Dp3[dim][dim3] * Dp3[dim3][dim2];
                }
            }
        }

        for (int dim = 0; dim < DIM; dim++) {
            for (int dim2 = 0; dim2 < DIM; dim2++) {
                D2p4[dim][dim2] = 0.0;
                for (int dim3 = 0; dim3 < DIM; dim3++) {
                    D2p4[dim][dim2] += Dp4[dim][dim3] * Dp4[dim3][dim2];
                }
            }
        }

        for (int dim = 0; dim < DIM; dim++) {
            for (int dim2 = 0; dim2 < DIM; dim2++) {
                D2p5[dim][dim2] = 0.0;
                for (int dim3 = 0; dim3 < DIM; dim3++) {
                    D2p5[dim][dim2] += Dp5[dim][dim3] * Dp5[dim3][dim2];
                }
            }
        }

        for (int dim = 0; dim < DIM; dim++) {
            for (int dim2 = 0; dim2 < DIM; dim2++) {
                D2p6[dim][dim2] = 0.0;
                for (int dim3 = 0; dim3 < DIM; dim3++) {
                    D2p6[dim][dim2] += Dp6[dim][dim3] * Dp6[dim3][dim2];
                }
            }
        } 

         for (int dim = 0; dim < DIM; dim++) {
            for (int dim2 = 0; dim2 < DIM; dim2++) {
                D2p7[dim][dim2] = 0.0;
                for (int dim3 = 0; dim3 < DIM; dim3++) {
                    D2p7[dim][dim2] += Dp7[dim][dim3] * Dp7[dim3][dim2];
                }
            }
        }                
        //shear rate calculation
        real qp0 = 0.0, qp1 = 0.0, qp2 = 0.0, qp3 = 0.0;
        real qp4 = 0.0, qp5 = 0.0, qp6 = 0.0, qp7 = 0.0;
        /*
        for (int dim = 0; dim < DIM; dim++) {
                qp0 += D2p0[dim][dim];
        }
        qp0         = sqrt(2.0*qp0);
        for (int dim = 0; dim < DIM; dim++) {
                qp1 += D2p1[dim][dim];
        }
        qp1         = sqrt(2.0*qp1);
        for (int dim = 0; dim < DIM; dim++) {
                qp2 += D2p2[dim][dim];
        }
        qp2         = sqrt(2.0*qp3);
        for (int dim = 0; dim < DIM; dim++) {
                qp3 += D2p3[dim][dim];
        }
        qp3         = sqrt(2.0*qp3);
        for (int dim = 0; dim < DIM; dim++) {
                qp4 += D2p4[dim][dim];
        }
        qp4         = sqrt(2.0*qp4);
        for (int dim = 0; dim < DIM; dim++) {
                qp5 += D2p5[dim][dim];
        }
        qp5         = sqrt(2.0*qp5);
        for (int dim = 0; dim < DIM; dim++) {
                qp6 += D2p6[dim][dim];
        }
        qp6         = sqrt(2.0*qp6);
        for (int dim = 0; dim < DIM; dim++) {
                qp7 += D2p7[dim][dim];
        }
        qp7         = sqrt(2.0*qp7);
       qp0 = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.stsp.dpD[0][1], ns->ed.stn);
       qp1 = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.stsp.dpD[0][1], ns->ed.stn);
       qp2 = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.stsp.dpD[0][1], ns->ed.stn);
       qp3 = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.stsp.dpD[0][1], ns->ed.stn);
       qp4 = compute_value_at_point(ns->ed.sdED, p4, p4, 1.0, ns->ed.stsp.dpD[0][1], ns->ed.stn);
       qp5 = compute_value_at_point(ns->ed.sdED, p5, p5, 1.0, ns->ed.stsp.dpD[0][1], ns->ed.stn);
       qp6 = compute_value_at_point(ns->ed.sdED, p6, p6, 1.0, ns->ed.stsp.dpD[0][1], ns->ed.stn);
       qp7 = compute_value_at_point(ns->ed.sdED, p7, p7, 1.0, ns->ed.stsp.dpD[0][1], ns->ed.stn);

        //calculating the shear stress
        taup0 = compute_value_at_point(ns->ed.sdED, p0, p0, 1.0, ns->ed.stsp.dpS[0][1], ns->ed.stn);
        taup0 += 2.0 * (ns->ed.stsp.par.eta0) * 0.5 * (Dp0[0][1] + Dp0[1][0]);
        taup1 = compute_value_at_point(ns->ed.sdED, p1, p1, 1.0, ns->ed.stsp.dpS[0][1], ns->ed.stn);
        taup1 += 2.0 * (ns->ed.stsp.par.eta0) * 0.5 * (Dp1[0][1] + Dp1[1][0]);
        taup2 = compute_value_at_point(ns->ed.sdED, p2, p2, 1.0, ns->ed.stsp.dpS[0][1], ns->ed.stn);
        taup2 += 2.0 * (ns->ed.stsp.par.eta0) * 0.5 * (Dp2[0][1] + Dp2[1][0]);
        taup3 = compute_value_at_point(ns->ed.sdED, p3, p3, 1.0, ns->ed.stsp.dpS[0][1], ns->ed.stn);
        taup3 += 2.0 * (ns->ed.stsp.par.eta0) * 0.5 * (Dp3[0][1] + Dp3[1][0]);
        taup4 = compute_value_at_point(ns->ed.sdED, p4, p4, 1.0, ns->ed.stsp.dpS[0][1], ns->ed.stn);
        taup4 += 2.0 * (ns->ed.stsp.par.eta0) * 0.5 * (Dp4[0][1] + Dp4[1][0]);
        taup5 = compute_value_at_point(ns->ed.sdED, p5, p5, 1.0, ns->ed.stsp.dpS[0][1], ns->ed.stn);
        taup5 += 2.0 * (ns->ed.stsp.par.eta0) * 0.5 * (Dp5[0][1] + Dp5[1][0]);
        taup6 = compute_value_at_point(ns->ed.sdED, p6, p6, 1.0, ns->ed.stsp.dpS[0][1], ns->ed.stn);
        taup6 += 2.0 * (ns->ed.stsp.par.eta0) * 0.5 * (Dp6[0][1] + Dp6[1][0]);
        taup7 = compute_value_at_point(ns->ed.sdED, p7, p7, 1.0, ns->ed.stsp.dpS[0][1], ns->ed.stn);
        taup7 += 2.0 * (ns->ed.stsp.par.eta0) * 0.5 * (Dp7[0][1] + Dp7[1][0]);

        //Calculate the viscosity eta = Txy/(q*eta0)
        real etap0 = 0.0, etap1 = 0.0, etap2 = 0.0, etap3 = 0.0;
        real etap4 = 0.0, etap5 = 0.0, etap6 = 0.0, etap7 = 0.0;

        etap0 = fabs((taup0)/((ns->ed.stsp.par.eta0)*qp0));
        fprintf(f, "%e\n", etap0);
        etap1 = fabs((taup1)/((ns->ed.stsp.par.eta0)*qp1));
        fprintf(f, "%e\n", etap1);
        etap2 = fabs((taup2)/((ns->ed.stsp.par.eta0)*qp2));
        fprintf(f, "%e\n", etap2);
        etap3 = fabs((taup3)/((ns->ed.stsp.par.eta0)*qp3));
        fprintf(f, "%e\n", etap3);
        etap4 = fabs((taup4)/((ns->ed.stsp.par.eta0)*qp4));
        fprintf(f, "%e\n", etap4);
        etap5 = fabs((taup5)/((ns->ed.stsp.par.eta0)*qp5));
        fprintf(f, "%e\n", etap5);
        etap6 = fabs((taup6)/((ns->ed.stsp.par.eta0)*qp6));
        fprintf(f, "%e\n", etap6);
        etap7 = fabs((taup7)/((ns->ed.stsp.par.eta0)*qp7));
        fprintf(f, "%e\n", etap7);

     //   real val = compute_value_at_point(ns->ed.vevv.sdVisc, ccenter, ccenter, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
        //real val = dp_get_value(ns->ed.vevv.dpvisc, cgid);
        //fprintf(f, "%e\n", val);
    }
    higcit_destroy(it);
    */

    fclose(f);
}


// *******************************************************************
// Navier-Stokes Load and Save
// *******************************************************************

// Loading the data files
void higflow_load_data_files(int argc, char *argv[], higflow_solver *ns)
{
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
void higflow_load_controllers(higflow_solver *ns, int myrank)
{
    // Controllers file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.contr", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL)
    {
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
        if (myrank == 0)
        {
            printf("=+=+=+= Controllers =+=+=+=\n");
            switch (ns->contr.projtype)
            {
            case 0:
                printf("=+=+=+= Projection Method: Non Incremental =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Projection Method: Incremental =+=+=+=\n");
                break;
            }
            switch (ns->contr.flowtype)
            {
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
            case 6:
                printf("=+=+=+= Flow Model: Viscoelastic with variable viscosity=+=+=+=\n");
                break;
            case 7:
                printf("=+=+=+= Flow Model: Viscoelastic with shear-banding=+=+=+=\n");
                break;
            case 8:
                printf("=+=+=+= Flow Model: Elastoviscoplastic =+=+=+=\n");
                break;
            case 9:
                printf("=+=+=+= Flow Model: Suspensions of spherical particles =+=+=+=\n");
                break;
            }
            switch (ns->contr.modelflowtype)
            {
            case 1:
                printf("=+= Flow Model Type: Eletroosmotic =+=\n");
                break;
            case 2:
                printf("=+=+=+= Flow Model Type: User defined model for variable viscosity =+=+=+=\n");
                break;
            case 3:
                printf("=+=+=+= Flow Model: Thixotropic-viscoelastic models (Oldroyd-B model + structural parameter equation) =+=+=+=\n");
                break;
            case 4:
                printf("=+=+=+= Flow Model: Two species model for shear-banding fluids =+=+=+=\n");
                break;
            case 5:
                printf("=+=+=+= Flow Model: Non-isothermal flows =+=+=+=\n");
                break;
            }
            switch (ns->contr.tempdiscrtype)
            {
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
            switch (ns->contr.spatialdiscrtype)
            {
            case 0:
                printf("=+=+=+= Spatial Discretization: Second Order =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Spatial Discretization: Forth Order =+=+=+=\n");
                break;
            }
            switch (ns->contr.convecdiscrtype)
            {
            case 0:
                printf("=+=+=+= Convective Scheme: Central =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Convective Scheme: First Order Upwind =+=+=+=\n");
                break;
            case 2:
                printf("=+=+=+= Convective Scheme: Second Order Upwind =+=+=+=\n");
                switch (ns->contr.secondconvecdiscrtype)
                {
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
    }
    else
    {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the controllers
void higflow_save_controllers(higflow_solver *ns, int myrank)
{
    // Controllers file name
    if (myrank == 0)
    {
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.contr", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL)
        {
            // Saving the controllers
            fprintf(fd, "%d\n", (ns->contr.projtype));
            fprintf(fd, "%d\n", (ns->contr.flowtype));
            fprintf(fd, "%d\n", (ns->contr.modelflowtype));
            fprintf(fd, "%d\n", (ns->contr.tempdiscrtype));
            fprintf(fd, "%d\n", (ns->contr.spatialdiscrtype));
            fprintf(fd, "%d\n", (ns->contr.convecdiscrtype));
            fprintf(fd, "%d\n", (ns->contr.secondconvecdiscrtype));
            fclose(fd);
        }
        else
        {
            // Error in open the file
            printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}

// Loading the parameters
void higflow_load_parameters(higflow_solver *ns, int myrank)
{
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.par", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL)
    {
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
        if (myrank == 0)
        {
            printf("=+=+=+= Reynolds Number: %f =+=+=+=\n", ns->par.Re);
        }
    }
    else
    {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the parameters
void higflow_save_parameters(higflow_solver *ns, int myrank)
{
    if (myrank == 0)
    {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.par", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL)
        {
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
        }
        else
        {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}

// Loading the properties
void higflow_load_properties(higflow_solver *ns)
{
    // Velocity file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.vel", ns->par.nameload);
    /*FILE *fd = fopen(namefile, "r");
    if (fd != NULL) {
        // Loading the velocity
        loadUV(ns->psfdu, ns->dpu, fd);
        fclose(fd);
    } else */
    {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
    // Pressure file name
    snprintf(namefile, sizeof namefile, "%s.pres", ns->par.nameload);
    /*fd = fopen(namefile, "r");
    if (fd != NULL) {
        // Loading the pressure
        load_property_cell(ns->psdp, ns->dpp, fd);
        fclose(fd);
    } else */
    {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the properties
void higflow_save_properties(higflow_solver *ns, int myrank, int ntasks)
{
    // Velocity file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.vel", ns->par.namesave);
    /*FILE *fd = fopen(namefile, "w");
    if (fd != NULL) {
        // Saving the velocity
        saveUV(ns->psfdu, ns->dpu, fd, myrank, ntasks);
        fclose(fd);
    } else */
    {
        // Error in open the file
        printf("=+=+=+= Error saving file %s =+=+=+=\n", ns->par.namesave);
        exit(1);
    }
    // Pressure file name
    snprintf(namefile, sizeof namefile, "%s.pres", ns->par.namesave);
    /*fd = fopen(namefile, "w");
    if (fd != NULL) {
        // Saving the pressure
        save_property_cell(ns->psdp, ns->dpp, fd, myrank, ntasks);
        fclose(fd);
    } else */
    {
        // Error in open the file
        printf("=+=+=+= Error saving file %s =+=+=+=\n", ns->par.namesave);
        exit(1);
    }
}

// Loading the viscoelastic parameters
void higflow_load_viscoelastic_parameters(higflow_solver *ns, int myrank)
{
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.viscpar", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL)
    {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%lf", &(ns->ed.ve.par.De));
        ifd = fscanf(fd, "%lf", &(ns->ed.ve.par.beta));
        ifd = fscanf(fd, "%lf", &(ns->ed.ve.par.epsilon));
        ifd = fscanf(fd, "%lf", &(ns->ed.ve.par.psi));
        ifd = fscanf(fd, "%lf", &(ns->ed.ve.par.alpha));
        ifd = fscanf(fd, "%lf", &(ns->ed.ve.par.kernel_tol));
        ifd = fscanf(fd, "%lf", &(ns->ed.ve.par.alpha_gptt));
        ifd = fscanf(fd, "%lf", &(ns->ed.ve.par.beta_gptt));
        fclose(fd);
        if (myrank == 0)
        {
            printf("=+=+=+= Deborah Number: %f =+=+=+=\n", ns->ed.ve.par.De);
            printf("=+=+=+= Beta: %f =+=+=+=\n", ns->ed.ve.par.beta);
        }
    }
    else
    {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the viscoelastic parameters
void higflow_save_viscoelastic_parameters(higflow_solver *ns, int myrank)
{
    if (myrank == 0)
    {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.viscpar", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL)
        {
            // Saving the parameters
            fprintf(fd, "%lf\n", (ns->ed.ve.par.De));
            fprintf(fd, "%lf\n", (ns->ed.ve.par.beta));
            fprintf(fd, "%lf\n", (ns->ed.ve.par.epsilon));
            fprintf(fd, "%lf\n", (ns->ed.ve.par.psi));
            fprintf(fd, "%lf\n", (ns->ed.ve.par.alpha));
            fprintf(fd, "%lf\n", (ns->ed.ve.par.kernel_tol));
            fprintf(fd, "%lf\n", (ns->ed.ve.par.alpha_gptt));
            fprintf(fd, "%lf\n", (ns->ed.ve.par.beta_gptt));
            fclose(fd);
        }
        else
        {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}

// Loading the viscoelastic controllers
void higflow_load_viscoelastic_controllers(higflow_solver *ns, int myrank)
{
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.visccontr", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL)
    {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%d", &(ns->ed.ve.contr.model));
        ifd = fscanf(fd, "%d", &(ns->ed.ve.contr.discrtype));
        ifd = fscanf(fd, "%d", &(ns->ed.ve.contr.convecdiscrtype));
        fclose(fd);
        if (myrank == 0)
        {
            switch (ns->ed.ve.contr.model)
            {
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
            case 3:
                printf("=+=+=+= Constitutive Equation Model: GPTT =+=+=+=\n");
                break;
            }
            switch (ns->ed.ve.contr.discrtype)
            {
            case 0:
                printf("=+=+=+= Constitutive Equation Discretization: Explicit =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Constitutive Equation Discretization: Implicit =+=+=+=\n");
                break;
            }
            switch (ns->ed.ve.contr.convecdiscrtype)
            {
            case 0:
                printf("=+=+=+= Constitutive Equation Convective Term: Upwind  =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Constitutive Equation Convective Term: CUBISTA =+=+=+=\n");
                break;
            }
        }
    }
    else
    {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the viscoelastic controllers
void higflow_save_viscoelastic_controllers(higflow_solver *ns, int myrank)
{
    if (myrank == 0)
    {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.visccontr", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL)
        {
            // Saving the parameters
            fprintf(fd, "%d\n", (ns->ed.ve.contr.model));
            fprintf(fd, "%d\n", (ns->ed.ve.contr.discrtype));
            fprintf(fd, "%d\n", (ns->ed.ve.contr.convecdiscrtype));
            fclose(fd);
        }
        else
        {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}

// Loading the electroosmotic parameters
void higflow_load_electroosmotic_parameters(higflow_solver *ns, int myrank)
{
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.eopar", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL)
    {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%lf", &(ns->ed.eo.par.alpha));
        ifd = fscanf(fd, "%lf", &(ns->ed.eo.par.delta));
        ifd = fscanf(fd, "%lf", &(ns->ed.eo.par.Pe));
        ifd = fscanf(fd, "%lf", &(ns->ed.eo.par.dphidx));
        fclose(fd);
        if (myrank == 0)
        {
            printf("=+=+=+= alpha: %e =+=+=+=\n", ns->ed.eo.par.alpha);
            printf("=+=+=+= delta: %e =+=+=+=\n", ns->ed.eo.par.delta);
            printf("=+=+=+= Peclet number: %f =+=+=+=\n", ns->ed.eo.par.Pe);
            printf("=+=+=+= Externel field dphidx : %e =+=+=+=\n", ns->ed.eo.par.dphidx);
        }
    }
    else
    {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the electroosmotic parameters
void higflow_save_electroosmotic_parameters(higflow_solver *ns, int myrank)
{
    if (myrank == 0)
    {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.viscpar", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL)
        {
            // Saving the parameters
            fprintf(fd, "%lf\n", (ns->ed.eo.par.alpha));
            fprintf(fd, "%lf\n", (ns->ed.eo.par.delta));
            fprintf(fd, "%lf\n", (ns->ed.eo.par.Pe));
            fprintf(fd, "%lf\n", (ns->ed.eo.par.dphidx));
            fclose(fd);
        }
        else
        {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}

// Loading the electroosmotic controllers
void higflow_load_electroosmotic_controllers(higflow_solver *ns, int myrank)
{
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.eocontr", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL)
    {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%d", &(ns->ed.eo.contr.eo_model));
        ifd = fscanf(fd, "%d", &(ns->ed.eo.contr.convecdiscrtype));
        fclose(fd);
        if (myrank == 0)
        {
            switch (ns->ed.eo.contr.eo_model)
            {
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
                printf("=+=+=+= Function f(T): Modified-Arrhenius =+=+=+=\n");
                break;
            }
            switch (ns->ed.eo.contr.convecdiscrtype)
            {
            case 0:
                printf("=+=+=+= Ionic concentration convective term: Central =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Ionic concentration convective term: Cubista =+=+=+=\n");
                break;
            }
        }
    }
    else
    {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the electroosmotic controllers
void higflow_save_electroosmotic_controllers(higflow_solver *ns, int myrank)
{
    if (myrank == 0)
    {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.eocontr", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL)
        {
            // Saving the parameters
            fprintf(fd, "%d\n", (ns->ed.eo.contr.eo_model));
            fclose(fd);
        }
        else
        {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}

// Loading the viscoelastic parameters
void higflow_load_viscoelastic_integral_parameters(higflow_solver *ns, int myrank)
{
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.viscintpar", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL)
    {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%lf", &(ns->ed.im.par.De));
        ifd = fscanf(fd, "%lf", &(ns->ed.im.par.alpha)); //for damping function
        ifd = fscanf(fd, "%lf", &(ns->ed.im.par.beta));  //for damping function
        ifd = fscanf(fd, "%lf", &(ns->ed.im.par.scorte));
        ifd = fscanf(fd, "%lf", &(ns->ed.im.par.rho));
        ifd = fscanf(fd, "%lf", &(ns->ed.im.par.v_ref));
        ifd = fscanf(fd, "%lf", &(ns->ed.im.par.l_ref));
        if (ns->ed.im.contr.model == 0)
        {
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
        }
        else if (ns->ed.im.contr.model == 1)
        {
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.alpha_frac));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.beta_frac));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.Phi1));
            ifd = fscanf(fd, "%lf", &(ns->ed.im.par.Phi2));
        }
        fclose(fd);
        if (myrank == 0)
        {
            printf("=+=+=+= Deborah Number: %f =+=+=+=\n", ns->ed.im.par.De);
            printf("=+=+=+= Alpha: %f =+=+=+=\n", ns->ed.im.par.alpha);
            printf("=+=+=+= Beta: %f =+=+=+=\n", ns->ed.im.par.beta);
            printf("=+=+=+= Scorte: %f =+=+=+=\n", ns->ed.im.par.scorte);
            printf("=+=+=+= Rho: %f =+=+=+=\n", ns->ed.im.par.rho);
            printf("=+=+=+= v_ref: %f =+=+=+=\n", ns->ed.im.par.v_ref);
            printf("=+=+=+= l_ref: %f =+=+=+=\n", ns->ed.im.par.l_ref);
            if (ns->ed.im.contr.model == 0)
            {
                printf("=+=+=+= M: %d =+=+=+=\n", ns->ed.im.par.M);
                for (int i = 0; i < 8; i++)
                {
                    printf("=+=+=+= a(%d): %f =+=+=+=\n", i, ns->ed.im.par.a[i]);
                }
                for (int i = 0; i < 8; i++)
                {
                    printf("=+=+=+= Lambda(%d): %f =+=+=+=\n", i, ns->ed.im.par.lambda[i]);
                }
            }
            else if (ns->ed.im.contr.model == 1)
            {
                printf("=+=+=+= Alpha Fractional: %f =+=+=+=\n", ns->ed.im.par.alpha_frac);
                printf("=+=+=+= Beta Fractional: %f =+=+=+=\n", ns->ed.im.par.beta_frac);
                printf("=+=+=+= Phi1: %f =+=+=+=\n", ns->ed.im.par.Phi1);
                printf("=+=+=+= Phi2: %f =+=+=+=\n", ns->ed.im.par.Phi2);
            }
        }
    }
    else
    {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the viscoelastic integral parameters
void higflow_save_viscoelastic_integral_parameters(higflow_solver *ns, int myrank)
{
    if (myrank == 0)
    {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.viscintpar", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL)
        {
            // Saving the parameters
            fprintf(fd, "%lf\n", (ns->ed.im.par.De));
            fprintf(fd, "%lf\n", (ns->ed.im.par.alpha)); // for damping function
            fprintf(fd, "%lf\n", (ns->ed.im.par.beta));  // for damping function
            fprintf(fd, "%lf\n", (ns->ed.im.par.scorte));
            fprintf(fd, "%lf\n", (ns->ed.im.par.rho));
            fprintf(fd, "%lf\n", (ns->ed.im.par.v_ref));
            fprintf(fd, "%lf\n", (ns->ed.im.par.l_ref));
            if (ns->ed.im.contr.model == 0)
            {
                fprintf(fd, "%d\n", (ns->ed.im.par.M));
                for (int k = 0; k < 8; k++)
                {
                    fprintf(fd, "%lf\n", (ns->ed.im.par.a[k]));
                }
                for (int k = 0; k < 8; k++)
                {
                    fprintf(fd, "%lf\n", (ns->ed.im.par.lambda[k]));
                }
            }
            else if (ns->ed.im.contr.model == 1)
            {
                fprintf(fd, "%lf\n", (ns->ed.im.par.alpha_frac));
                fprintf(fd, "%lf\n", (ns->ed.im.par.beta_frac));
                fprintf(fd, "%lf\n", (ns->ed.im.par.Phi1));
                fprintf(fd, "%lf\n", (ns->ed.im.par.Phi2));
            }

            fclose(fd);
        }
        else
        {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}

// Loading the viscoelastic integral controllers
void higflow_load_viscoelastic_integral_controllers(higflow_solver *ns, int myrank)
{
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.viscintcontr", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL)
    {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%d", &(ns->ed.im.contr.model));
        ifd = fscanf(fd, "%d", &(ns->ed.im.contr.model_H));
        ifd = fscanf(fd, "%d", &(ns->ed.im.contr.discrtype));
        ifd = fscanf(fd, "%d", &(ns->ed.im.contr.convecdiscrtype));
        fclose(fd);
        if (myrank == 0)
        {
            switch (ns->ed.im.contr.model)
            {
            case 0:
                printf("=+=+=+= Constitutive Equation Model: KBKZ =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Constitutive Equation Model: Fractional =+=+=+=\n");
                break;
            }
            switch (ns->ed.im.contr.model_H)
            {
            case 0:
                printf("=+=+=+= Constitutive Equation Model: PSM =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Constitutive Equation Model: UCM =+=+=+=\n");
                break;
            }
            switch (ns->ed.im.contr.discrtype)
            {
            case 0:
                printf("=+=+=+= Constitutive Equation Discretization: Explicit =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Constitutive Equation Discretization: Implicit =+=+=+=\n");
                break;
            }
            switch (ns->ed.im.contr.convecdiscrtype)
            {
            case 0:
                printf("=+=+=+= Constitutive Equation Convective Term: Central  =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Constitutive Equation Convective Term: CUBISTA =+=+=+=\n");
                break;
            }
        }
    }
    else
    {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the viscoelastic integral controllers
void higflow_save_viscoelastic_integral_controllers(higflow_solver *ns, int myrank)
{
    if (myrank == 0)
    {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.viscintcontr", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL)
        {
            // Saving the parameters
            fprintf(fd, "%d\n", (ns->ed.im.contr.model));
            fprintf(fd, "%d\n", (ns->ed.im.contr.model_H));
            fprintf(fd, "%d\n", (ns->ed.im.contr.discrtype));
            fprintf(fd, "%d\n", (ns->ed.im.contr.convecdiscrtype));
            fclose(fd);
        }
        else
        {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}


// Loading the parameters of viscoelastic flows with variable viscosity
void higflow_load_viscoelastic_variable_viscosity_parameters(higflow_solver *ns, int myrank)
{
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.viscvvpar", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL)
    {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%lf", &(ns->ed.vevv.par.De));
        ifd = fscanf(fd, "%lf", &(ns->ed.vevv.par.beta));
        ifd = fscanf(fd, "%lf", &(ns->ed.vevv.par.epsilon));
        ifd = fscanf(fd, "%lf", &(ns->ed.vevv.par.psi));
        ifd = fscanf(fd, "%lf", &(ns->ed.vevv.par.alpha));
        ifd = fscanf(fd, "%lf", &(ns->ed.vevv.par.kernel_tol));
        ifd = fscanf(fd, "%lf", &(ns->ed.vevv.par.alpha_gptt));
        ifd = fscanf(fd, "%lf", &(ns->ed.vevv.par.beta_gptt));
        if (ns->contr.rheotype == THIXOTROPIC)
        {
            ifd = fscanf(fd, "%lf", &(ns->ed.vevv.par.Lambda));
            ifd = fscanf(fd, "%lf", &(ns->ed.vevv.par.Phi));
            ifd = fscanf(fd, "%lf", &(ns->ed.vevv.par.Gamma));            
        }
        fclose(fd);
        if (myrank == 0)
        {
            printf("=+=+=+= Deborah Number: %f =+=+=+=\n", ns->ed.vevv.par.De);
            printf("=+=+=+= Beta: %f =+=+=+=\n", ns->ed.vevv.par.beta);
            if (ns->contr.rheotype == THIXOTROPIC)
            {
                printf("=+=+=+= Lambda: %f =+=+=+=\n", ns->ed.vevv.par.Lambda);
                printf("=+=+=+= Phi: %f =+=+=+=\n", ns->ed.vevv.par.Phi);
                printf("=+=+=+= Gamma: %f =+=+=+=\n", ns->ed.vevv.par.Gamma);
            }
        }
    }
    else
    {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the parameters of the viscoelastic flows with variable viscosity 
void higflow_save_viscoelastic_variable_viscosity_parameters(higflow_solver *ns, int myrank)
{
    if (myrank == 0)
    {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.viscvvpar", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL)
        {
            // Saving the parameters
            fprintf(fd, "%lf\n", (ns->ed.vevv.par.De));
            fprintf(fd, "%lf\n", (ns->ed.vevv.par.beta));
            fprintf(fd, "%lf\n", (ns->ed.vevv.par.epsilon));
            fprintf(fd, "%lf\n", (ns->ed.vevv.par.psi));
            fprintf(fd, "%lf\n", (ns->ed.vevv.par.alpha));
            fprintf(fd, "%lf\n", (ns->ed.vevv.par.kernel_tol));
            fprintf(fd, "%lf\n", (ns->ed.vevv.par.alpha_gptt));
            fprintf(fd, "%lf\n", (ns->ed.vevv.par.beta_gptt));
            if (ns->contr.rheotype == THIXOTROPIC)
            {
                fprintf(fd, "%lf\n", (ns->ed.vevv.par.Lambda));
                fprintf(fd, "%lf\n", (ns->ed.vevv.par.Phi));
                fprintf(fd, "%lf\n", (ns->ed.vevv.par.Gamma));
            }
            fclose(fd);
        }
        else
        {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}

// Loading the controllers for viscoelastic flows with variable viscosity
void higflow_load_viscoelastic_variable_viscosity_controllers(higflow_solver *ns, int myrank)
{
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.viscvvcontr", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL)
    {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%d", &(ns->ed.vevv.contr.model));
        ifd = fscanf(fd, "%d", &(ns->ed.vevv.contr.discrtype));
        ifd = fscanf(fd, "%d", &(ns->ed.vevv.contr.convecdiscrtype));
        ifd = fscanf(fd, "%d", &(ns->ed.vevv.contr.structpdiscrtype));
        ifd = fscanf(fd, "%d", &(ns->ed.vevv.contr.structpconvecdiscrtype));
        ifd = fscanf(fd, "%d", &(ns->ed.vevv.contr.structparmodel));
        fclose(fd);
        if (myrank == 0)
        {
            switch (ns->ed.vevv.contr.model)
            {
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
            case 3:
                printf("=+=+=+= Constitutive Equation Model: GPTT =+=+=+=\n");
                break;
            }
            switch (ns->ed.vevv.contr.discrtype)
            {
            case 0:
                printf("=+=+=+= Constitutive Equation Discretization: Explicit =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Constitutive Equation Discretization: Implicit =+=+=+=\n");
                break;
            }
            switch (ns->ed.vevv.contr.convecdiscrtype)
            {
            case 0:
                printf("=+=+=+= Constitutive Equation Convective Term: Upwind  =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Constitutive Equation Convective Term: CUBISTA =+=+=+=\n");
                break;
            }
            if (ns->contr.rheotype == THIXOTROPIC)
            {
                if (ns->ed.vevv.contr.structparmodel == 0)
                {
                    switch (ns->ed.vevv.contr.structpdiscrtype)
                    {
                    case 0:
                        printf("=+=+=+= BMP Model Viscosity Evolution Equation Discretization: Explicit  =+=+=+=\n"); 
                        break;
                    case 1:
                        printf("=+=+=+= BMP Model Viscosity Evolution Equation Discretization: Implicit  =+=+=+=\n");
                        break;
                    }
                    switch (ns->ed.vevv.contr.structpconvecdiscrtype)
                    {
                    case 0:
                        printf("=+=+=+= BMP Model Viscosity Evolution Equation Convective Term: Upwind  =+=+=+=\n");
                        break;
                    case 1:
                        printf("=+=+=+= BMP Model Viscosity Evolution Equation Convective Term: CUBISTA  =+=+=+=\n");
                        break;            
                    }
                }
                
                if (ns->ed.vevv.contr.structparmodel == 1)
                {
                    switch (ns->ed.vevv.contr.structpdiscrtype)
                    {
                    case 0:
                        printf("=+=+=+= BMP Model with solvent Viscosity Evolution Equation Discretization: Explicit  =+=+=+=\n"); 
                        break;
                    case 1:
                        printf("=+=+=+= BMP Model with solvent Viscosity Evolution Equation Discretization: Implicit  =+=+=+=\n");
                        break;
                    }
                    switch (ns->ed.vevv.contr.structpconvecdiscrtype)
                    {
                    case 0:
                        printf("=+=+=+= BMP Model with solvent Viscosity Evolution Equation Convective Term: Upwind  =+=+=+=\n");
                        break;
                    case 1:
                        printf("=+=+=+= BMP Model with solvent Viscosity Evolution Equation Convective Term: CUBISTA  =+=+=+=\n");
                        break;            
                    }
                }

                if (ns->ed.vevv.contr.structparmodel == 2)
                {
                    switch (ns->ed.vevv.contr.structpdiscrtype)
                    {
                    case 0:
                        printf("=+=+=+= MBM Model Viscosity Evolution Equation Discretization: Explicit  =+=+=+=\n"); 
                        break;
                    case 1:
                        printf("=+=+=+= MBM Model Viscosity Evolution Equation Discretization: Implicit  =+=+=+=\n");
                        break;
                    }
                    switch (ns->ed.vevv.contr.structpconvecdiscrtype)
                    {
                    case 0:
                        printf("=+=+=+= MBM Model Viscosity Evolution Equation Convective Term: Upwind  =+=+=+=\n");
                        break;
                    case 1:
                        printf("=+=+=+= MBM Model Viscosity Evolution Equation Convective Term: CUBISTA  =+=+=+=\n");
                        break;            
                    }
                }

                if (ns->ed.vevv.contr.structparmodel == 3)
                {
                    switch (ns->ed.vevv.contr.structpdiscrtype)
                    {
                    case 0:
                        printf("=+=+=+= NM_taup Model Viscosity Evolution Equation Discretization: Explicit  =+=+=+=\n"); 
                        break;
                    case 1:
                        printf("=+=+=+= NM_taup Model Viscosity Evolution Equation Discretization: Implicit  =+=+=+=\n");
                        break;
                    }
                    switch (ns->ed.vevv.contr.structpconvecdiscrtype)
                    {
                    case 0:
                        printf("=+=+=+= NM_taup Model Viscosity Evolution Equation Convective Term: Upwind  =+=+=+=\n");
                        break;
                    case 1:
                        printf("=+=+=+= NM_taup Model Viscosity Evolution Equation Convective Term: CUBISTA  =+=+=+=\n");
                        break;            
                    }
                }

                if (ns->ed.vevv.contr.structparmodel == 4)
                {
                    switch (ns->ed.vevv.contr.structpdiscrtype)
                    {
                    case 0:
                        printf("=+=+=+= NM_T Model Viscosity Evolution Equation Discretization: Explicit  =+=+=+=\n"); 
                        break;
                    case 1:
                        printf("=+=+=+= NM_T Model Viscosity Evolution Equation Discretization: Implicit  =+=+=+=\n");
                        break;
                    }
                    switch (ns->ed.vevv.contr.structpconvecdiscrtype)
                    {
                    case 0:
                        printf("=+=+=+= NM_T Model Viscosity Evolution Equation Convective Term: Upwind  =+=+=+=\n");
                        break;
                    case 1:
                        printf("=+=+=+= NM_T Model Viscosity Evolution Equation Convective Term: CUBISTA  =+=+=+=\n");
                        break;            
                    }
                }
            }
        }
    }
    else
    {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the controllers for viscoelastic flows with variable viscosity
void higflow_save_viscoelastic_variable_viscosity_controllers(higflow_solver *ns, int myrank)
{
    if (myrank == 0)
    {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.viscvvcontr", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL)
        {
            // Saving the parameters
            fprintf(fd, "%d\n", (ns->ed.vevv.contr.model));
            fprintf(fd, "%d\n", (ns->ed.vevv.contr.discrtype));
            fprintf(fd, "%d\n", (ns->ed.vevv.contr.convecdiscrtype));
            if (ns->contr.rheotype == THIXOTROPIC)
            {
                fprintf(fd, "%d\n", (ns->ed.vevv.contr.structpdiscrtype));
                fprintf(fd, "%d\n", (ns->ed.vevv.contr.structpconvecdiscrtype));
                fprintf(fd, "%d\n", (ns->ed.vevv.contr.structparmodel));
            }
            fclose(fd);
        }
        else
        {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}


// Loading the parameters of viscoelastic flows with shear-banding
void higflow_load_viscoelastic_variable_shear_banding_parameters(higflow_solver *ns, int myrank)
{
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.viscsbpar", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL)
    {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%lf", &(ns->ed.vesb.par.De));
        ifd = fscanf(fd, "%lf", &(ns->ed.vesb.par.beta));
        if (ns->contr.rheotype == VCM)
        {
            ifd = fscanf(fd, "%lf", &(ns->ed.vesb.par.DeA));
            ifd = fscanf(fd, "%lf", &(ns->ed.vesb.par.epsilon));
            ifd = fscanf(fd, "%lf", &(ns->ed.vesb.par.PeA));
            ifd = fscanf(fd, "%lf", &(ns->ed.vesb.par.PeB));
            ifd = fscanf(fd, "%lf", &(ns->ed.vesb.par.chi));
            ifd = fscanf(fd, "%lf", &(ns->ed.vesb.par.CAeq));
            ifd = fscanf(fd, "%lf", &(ns->ed.vesb.par.CBeq));           
        }

        fclose(fd);
        if (myrank == 0)
        {
            printf("=+=+=+= Deborah Number: %f =+=+=+=\n", ns->ed.vesb.par.De);
            printf("=+=+=+= Beta: %f =+=+=+=\n", ns->ed.vesb.par.beta);
            if (ns->contr.rheotype == VCM)
            {
                printf("=+=+=+= Deborah Number of specie A: %f =+=+=+=\n", ns->ed.vesb.par.DeA);
                printf("=+=+=+= Epsilon = lambdaB/lambdaA: %f =+=+=+=\n", ns->ed.vesb.par.epsilon);
                printf("=+=+=+= Peclet number of specie A: %f =+=+=+=\n", ns->ed.vesb.par.PeA);
                printf("=+=+=+= Peclet number of specie B: %f =+=+=+=\n", ns->ed.vesb.par.PeB);
                printf("=+=+=+= Chi parameter: %f =+=+=+=\n", ns->ed.vesb.par.chi);
                printf("=+=+=+= Equilibrium concentration of A: %f =+=+=+=\n", ns->ed.vesb.par.CAeq);
                printf("=+=+=+= Equilibrium concentration of B: %f =+=+=+=\n", ns->ed.vesb.par.CBeq);
            }
        }
    }
    else
    {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the parameters of the viscoelastic flows with shear-banding
void higflow_save_viscoelastic_shear_banding_parameters(higflow_solver *ns, int myrank)
{
    if (myrank == 0)
    {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.viscsbpar", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL)
        {
            // Saving the parameters
            fprintf(fd, "%lf\n", (ns->ed.vesb.par.De));
            fprintf(fd, "%lf\n", (ns->ed.vesb.par.beta));
            if (ns->contr.rheotype == VCM)
            {
                fprintf(fd, "%lf\n", (ns->ed.vesb.par.DeA));
                fprintf(fd, "%lf\n", (ns->ed.vesb.par.epsilon));
                fprintf(fd, "%lf\n", (ns->ed.vesb.par.PeA));
                fprintf(fd, "%lf\n", (ns->ed.vesb.par.PeB));
                fprintf(fd, "%lf\n", (ns->ed.vesb.par.chi));
                fprintf(fd, "%lf\n", (ns->ed.vesb.par.CAeq));
                fprintf(fd, "%lf\n", (ns->ed.vesb.par.CBeq));  
            }
            fclose(fd);
        }
        else
        {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}

// Loading the controllers for viscoelastic flows with shear-banding
void higflow_load_viscoelastic_shear_banding_controllers(higflow_solver *ns, int myrank)
{
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.viscsbcontr", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL)
    {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%d", &(ns->ed.vesb.contr.model));
        ifd = fscanf(fd, "%d", &(ns->ed.vesb.contr.discrtype));
        ifd = fscanf(fd, "%d", &(ns->ed.vesb.contr.convecdiscrtype));
        if (ns->contr.rheotype == VCM)
        {
            ifd = fscanf(fd, "%d", &(ns->ed.vesb.contr.nAnBdiscrtype));
            ifd = fscanf(fd, "%d", &(ns->ed.vesb.contr.nAnBconvecdiscrtype));     
        }
        fclose(fd);
        if (myrank == 0)
        {
            switch (ns->ed.vesb.contr.model)
            {
            case 0:
                printf("=+=+=+= The Vazquez-McKinley-Cook (VCM) model =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= The modified Vazquez-McKinley-Cook (VCM) model =+=+=+=\n");
                break;
            }
            switch (ns->ed.vesb.contr.discrtype)
            {
            case 0:
                printf("=+=+=+= Viscoelastic Equation Discretization: Explicit =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Viscoelastic Equation Discretization: Implicit =+=+=+=\n");
                break;
            }
            switch (ns->ed.vesb.contr.convecdiscrtype)
            {
            case 0:
                printf("=+=+=+= Viscoelastic Equation Convective Term: Upwind  =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Viscoelastic Equation Convective Term: CUBISTA =+=+=+=\n");
                break;
            }
            if (ns->contr.rheotype == VCM)
            {
                switch (ns->ed.vesb.contr.nAnBdiscrtype)
                {
                    case 0:
                        printf("=+=+=+= Discretization of the Equations of the Density Numbers (nA and nB): Explicit  =+=+=+=\n"); 
                        break;
                    case 1:
                        printf("=+=+=+= Discretization of the Equations of the Density Numbers (nA and nB): Implicit  =+=+=+=\n");
                        break;
                }
                switch (ns->ed.vesb.contr.nAnBconvecdiscrtype)
                {
                    case 0:
                        printf("=+=+=+= Discretization of the Convective Term of of the Equations of the Density Numbers (nA and nB): Upwind  =+=+=+=\n");
                        break;
                    case 1:
                        printf("=+=+=+= Discretization of the Convective Term of of the Equations of the Density Numbers (nA and nB): CUBISTA  =+=+=+=\n");
                        break;            
                }                
            }
        }
    }
    else
    {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the controllers for viscoelastic flows with shear_banding
void higflow_save_viscoelastic_shear_banding_controllers(higflow_solver *ns, int myrank)
{
    if (myrank == 0)
    {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.viscsbcontr", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL)
        {
            // Saving the parameters
            fprintf(fd, "%d\n", (ns->ed.vesb.contr.model));
            fprintf(fd, "%d\n", (ns->ed.vesb.contr.discrtype));
            fprintf(fd, "%d\n", (ns->ed.vesb.contr.convecdiscrtype));
            if (ns->contr.rheotype == VCM)
            {
                fprintf(fd, "%d\n", (ns->ed.vesb.contr.nAnBdiscrtype));
                fprintf(fd, "%d\n", (ns->ed.vesb.contr.nAnBconvecdiscrtype));     
            }
            fclose(fd);
        }
        else
        {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}


// Loading the elastoviscoplastic parameters
void higflow_load_elastoviscoplastic_parameters(higflow_solver *ns, int myrank)
{
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.visceplpar", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL)
    {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%lf", &(ns->ed.vepl.par.De));
        ifd = fscanf(fd, "%lf", &(ns->ed.vepl.par.beta));
        ifd = fscanf(fd, "%lf", &(ns->ed.vepl.par.Bi));
        ifd = fscanf(fd, "%lf", &(ns->ed.vepl.par.zeta));
        ifd = fscanf(fd, "%lf", &(ns->ed.vepl.par.kernel_tol));
        ifd = fscanf(fd, "%lf", &(ns->ed.vepl.par.epsilon));
        ifd = fscanf(fd, "%lf", &(ns->ed.vepl.par.Np));
        fclose(fd);
        if (myrank == 0)
        {
            printf("=+=+=+= Deborah Number: %f =+=+=+=\n", ns->ed.vepl.par.De);
            printf("=+=+=+= Beta: %f =+=+=+=\n", ns->ed.vepl.par.beta);
            printf("=+=+=+= Bingham number: %f =+=+=+=\n", ns->ed.vepl.par.Bi);
            printf("=+=+=+= Zeta parameter (viscoelastic): %f =+=+=+=\n", ns->ed.vepl.par.zeta);
            if ((ns->ed.vepl.contr.model == 1)|| (ns->ed.vepl.contr.model == -1))
            {
                //Oldroyd-B-Herschel-Bulkley
                printf("=+=+=+= Power-law coefficient (OBHB): %f =+=+=+=\n", ns->ed.vepl.par.Np);
            }
            if ((ns->ed.vepl.contr.model == 2) || (ns->ed.vepl.contr.model == 3) ||
                (ns->ed.vepl.contr.model == -1))
            {
                //LPTT-Bingham and EPTT-Bingham
                printf("=+=+=+= Epsilon parameter: %f =+=+=+=\n", ns->ed.vepl.par.epsilon);
            }
        }
    }
    else
    {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the elastoviscoplastic parameters
void higflow_save_elastoviscoplastic_parameters(higflow_solver *ns, int myrank)
{
    if (myrank == 0)
    {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.visceplpar", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL)
        {
            // Saving the parameters
            fprintf(fd, "%lf\n", (ns->ed.vepl.par.De));
            fprintf(fd, "%lf\n", (ns->ed.vepl.par.beta));
            fprintf(fd, "%lf\n", (ns->ed.vepl.par.Bi));
            fprintf(fd, "%lf\n", (ns->ed.vepl.par.zeta));
            fprintf(fd, "%lf\n", (ns->ed.vepl.par.kernel_tol));
            fprintf(fd, "%lf\n", (ns->ed.vepl.par.epsilon));
            fprintf(fd, "%lf\n", (ns->ed.vepl.par.Np));
            fclose(fd);
        }
        else
        {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}

// Loading the elastoviscoplastic controllers
void higflow_load_elastoviscoplastic_controllers(higflow_solver *ns, int myrank)
{
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.visceplcontr", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL)
    {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%d", &(ns->ed.vepl.contr.model));
        ifd = fscanf(fd, "%d", &(ns->ed.vepl.contr.discrtype));
        ifd = fscanf(fd, "%d", &(ns->ed.vepl.contr.convecdiscrtype));
        fclose(fd);
        if (myrank == 0)
        {
            switch (ns->ed.vepl.contr.model)
            {
            case -1:
                printf("=+=+=+= Constitutive Equation Model: General Saramito Model (User Defined) =+=+=+=\n");
                break;
            case 0:
                printf("=+=+=+= Constitutive Equation Model: Oldroyd-B Bingham =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Constitutive Equation Model: Oldroyd-B-Herschel-Bulkley =+=+=+=\n");
                break;
            case 2:
                printf("=+=+=+= Constitutive Equation Model: LPTT-Bingham =+=+=+=\n");
                break;
            case 3:
                printf("=+=+=+= Constitutive Equation Model: EPTT-Bingham =+=+=+=\n");
                break;
            }
            switch (ns->ed.vepl.contr.discrtype)
            {
            case 0:
                printf("=+=+=+= Constitutive Equation Discretization: Explicit =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Constitutive Equation Discretization: Implicit =+=+=+=\n");
                break;
            }
            switch (ns->ed.vepl.contr.convecdiscrtype)
            {
            case 0:
                printf("=+=+=+= Constitutive Equation Convective Term: Upwind  =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Constitutive Equation Convective Term: CUBISTA =+=+=+=\n");
                break;
            }
        }
    }
    else
    {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the elastoviscoplastic controllers
void higflow_save_elastoviscoplastic_controllers(higflow_solver *ns, int myrank)
{
    if (myrank == 0)
    {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.visceplcontr", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL)
        {
            // Saving the parameters
            fprintf(fd, "%d\n", (ns->ed.vepl.contr.model));
            fprintf(fd, "%d\n", (ns->ed.vepl.contr.discrtype));
            fprintf(fd, "%d\n", (ns->ed.vepl.contr.convecdiscrtype));
            fclose(fd);
        }
        else
        {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}


// Loading the shear-thickening suspension parameters
void higflow_load_shear_thickening_suspension_parameters(higflow_solver *ns, int myrank)
{
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.stsppar", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL)
    {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%lf", &(ns->ed.stsp.par.alpha));
        ifd = fscanf(fd, "%lf", &(ns->ed.stsp.par.eta0));
        ifd = fscanf(fd, "%lf", &(ns->ed.stsp.par.beta));
        ifd = fscanf(fd, "%lf", &(ns->ed.stsp.par.chij1));
        ifd = fscanf(fd, "%lf", &(ns->ed.stsp.par.chij2));
        ifd = fscanf(fd, "%lf", &(ns->ed.stsp.par.X0));
        ifd = fscanf(fd, "%lf", &(ns->ed.stsp.par.phi));
        ifd = fscanf(fd, "%lf", &(ns->ed.stsp.par.Pic));            
        ifd = fscanf(fd, "%lf", &(ns->ed.stsp.par.gdrms));
        ifd = fscanf(fd, "%lf", &(ns->ed.stsp.par.apsize));           
        fclose(fd);
        /*
        if (myrank == 0)
        {
            printf("=+=+=+= Deborah Number: %f =+=+=+=\n", ns->ed.vevv.par.De);
            printf("=+=+=+= Beta: %f =+=+=+=\n", ns->ed.vevv.par.beta);
            if (ns->contr.rheotype == THIXOTROPIC)
            {
                printf("=+=+=+= Lambda: %f =+=+=+=\n", ns->ed.vevv.par.Lambda);
                printf("=+=+=+= Phi: %f =+=+=+=\n", ns->ed.vevv.par.Phi);
                printf("=+=+=+= Gamma: %f =+=+=+=\n", ns->ed.vevv.par.Gamma);
            }
        }
        */
        if (myrank == 0)
        {
            printf("=+=+=+= Alpha parameter: %f =+=+=+=\n", ns->ed.stsp.par.alpha);
            printf("=+=+=+= Viscosity of the suspension: %f =+=+=+=\n", ns->ed.stsp.par.eta0);
            printf("=+=+=+= Microstructure association rate (beta): %f =+=+=+=\n", ns->ed.stsp.par.beta);
            if ((ns->ed.stsp.contr.model == GW_WC)|| (ns->ed.stsp.contr.model == USERSET_SM))
            {
                printf("=+=+=+= Extremal jamming point 1: %f =+=+=+=\n", ns->ed.stsp.par.chij1);
                printf("=+=+=+= Extremal jamming point 2: %f =+=+=+=\n", ns->ed.stsp.par.chij2);
                printf("=+=+=+= X0 parameter: %f =+=+=+=\n", ns->ed.stsp.par.X0);
                printf("=+=+=+= Particle volume fraction parameter: %f =+=+=+=\n", ns->ed.stsp.par.phi);
                printf("=+=+=+= Critical particle pressure: %f =+=+=+=\n", ns->ed.stsp.par.Pic);
            }
            //only for the model with particle migration
            if ((ns->ed.stsp.contr.model == GW_WC_IF))
            {
                printf("=+=+=+= Extremal jamming point 1: %f =+=+=+=\n", ns->ed.stsp.par.chij1);
                printf("=+=+=+= Extremal jamming point 2: %f =+=+=+=\n", ns->ed.stsp.par.chij2);
                printf("=+=+=+= X0 parameter: %f =+=+=+=\n", ns->ed.stsp.par.X0);
                printf("=+=+=+= Particle volume fraction at random close packing: %f =+=+=+=\n", ns->ed.stsp.par.phi);
                printf("=+=+=+= Critical particle pressure: %f =+=+=+=\n", ns->ed.stsp.par.Pic);
                printf("=+=+=+= Standard deviation of the shear rate fluctuations: %f =+=+=+=\n", ns->ed.stsp.par.gdrms);
                printf("=+=+=+= Particle size: %f =+=+=+=\n", ns->ed.stsp.par.apsize);
            }                          

        }
    }
    else
    {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}


// Saving the shear-thickening suspension parameters
void higflow_save_shear_thickening_suspension_parameters(higflow_solver *ns, int myrank)
{
    if (myrank == 0)
    {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.stsppar", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL)
        {
            // Saving the parameters
            fprintf(fd, "%lf\n", (ns->ed.stsp.par.alpha));
            fprintf(fd, "%lf\n", (ns->ed.stsp.par.eta0));
            fprintf(fd, "%lf\n", (ns->ed.stsp.par.beta));
            fprintf(fd, "%lf\n", (ns->ed.stsp.par.chij1));
            fprintf(fd, "%lf\n", (ns->ed.stsp.par.chij2));
            fprintf(fd, "%lf\n", (ns->ed.stsp.par.X0));
            fprintf(fd, "%lf\n", (ns->ed.stsp.par.phi));
            fprintf(fd, "%lf\n", (ns->ed.stsp.par.Pic));
            fprintf(fd, "%lf\n", (ns->ed.stsp.par.gdrms));
            fprintf(fd, "%lf\n", (ns->ed.stsp.par.apsize));
            fclose(fd);
        }
        else
        {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}


// Loading the shear-thickening suspension controllers
void higflow_load_shear_thickening_suspension_controllers(higflow_solver *ns, int myrank)
{
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.stspcontr", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL)
    {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%d", &(ns->ed.stsp.contr.model));
        ifd = fscanf(fd, "%d", &(ns->ed.stsp.contr.discrtype));
        ifd = fscanf(fd, "%d", &(ns->ed.stsp.contr.convecdiscrtype));
        ifd = fscanf(fd, "%d", &(ns->ed.stsp.contr.volfracdiscrtype));
        ifd = fscanf(fd, "%d", &(ns->ed.stsp.contr.volfracconvecdiscrtype));
        fclose(fd);
        if (myrank == 0)
        {
            switch (ns->ed.stsp.contr.model)
            {
            case -1:
                printf("=+=+=+= Constitutive Equation Model: (User Defined) =+=+=+=\n");
                break;
            case 0:
                printf("=+=+=+= Constitutive Equation Model: GW model for suspensions =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Constitutive Equation Model: GW-WC model for shear-thickening suspensions =+=+=+=\n");
                break;
            case 2:
                printf("=+=+=+= Constitutive Equation Model: GW-WC model for shear-thickening suspensions (Inhomogeneous flows) =+=+=+=\n");
                break;
            }
            switch (ns->ed.stsp.contr.discrtype)
            {
            case 0:
                printf("=+=+=+= Constitutive Equation Discretization: Explicit =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Constitutive Equation Discretization: Implicit =+=+=+=\n");
                break;
            }
            switch (ns->ed.stsp.contr.convecdiscrtype)
            {
            case 0:
                printf("=+=+=+= Constitutive Equation Convective Term: Upwind  =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Constitutive Equation Convective Term: CUBISTA =+=+=+=\n");
                break;
            }
            if (ns->ed.stsp.contr.model == GW_WC_IF)
            {
                switch (ns->ed.stsp.contr.volfracdiscrtype)
                {
                    case 0:
                        printf("=+=+=+= Volume Fraction Evolution Equation Discretization: Explicit  =+=+=+=\n"); 
                        break;
                    case 1:
                        printf("=+=+=+= Volume Fraction Evolution Equation Discretization: Implicit  =+=+=+=\n");
                        break;
                }
                switch (ns->ed.stsp.contr.volfracconvecdiscrtype)
                {
                    case 0:
                        printf("=+=+=+= Volume Fraction Evolution Equation Convective Term: Upwind  =+=+=+=\n");
                        break;
                    case 1:
                        printf("=+=+=+= Volume Fraction Evolution Equation Convective Term: CUBISTA  =+=+=+=\n");
                        break;            
                }
            }
        }
    }
    else
    {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}


// Saving the shear-thickening suspension controllers
void higflow_save_shear_thickening_suspension_controllers(higflow_solver *ns, int myrank)
{
    if (myrank == 0)
    {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.stspcontr", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL)
        {
            // Saving the parameters
            fprintf(fd, "%d\n", (ns->ed.stsp.contr.model));
            fprintf(fd, "%d\n", (ns->ed.stsp.contr.discrtype));
            fprintf(fd, "%d\n", (ns->ed.stsp.contr.convecdiscrtype));
            fprintf(fd, "%d\n", (ns->ed.stsp.contr.volfracdiscrtype));
            fprintf(fd, "%d\n", (ns->ed.stsp.contr.volfracconvecdiscrtype));
            fclose(fd);
        }
        else
        {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}


// Loading the non-isothermal flow parameters
void higflow_load_non_isothermal_flow_parameters(higflow_solver *ns, int myrank)
{
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.nifpar", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL)
    {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%lf", &(ns->ed.nif.par.Pr));
        ifd = fscanf(fd, "%lf", &(ns->ed.nif.par.Br));
        ifd = fscanf(fd, "%lf", &(ns->ed.nif.par.Gr));
        ifd = fscanf(fd, "%lf", &(ns->ed.nif.par.Ga));
        ifd = fscanf(fd, "%lf", &(ns->ed.nif.par.k0));
        ifd = fscanf(fd, "%lf", &(ns->ed.nif.par.k1));
        ifd = fscanf(fd, "%lf", &(ns->ed.nif.par.k2));
        ifd = fscanf(fd, "%lf", &(ns->ed.nif.par.k3));
        ifd = fscanf(fd, "%lf", &(ns->ed.nif.par.T0));
        ifd = fscanf(fd, "%lf", &(ns->ed.nif.par.T1));
        ifd = fscanf(fd, "%lf", &(ns->ed.nif.par.alphaT));
        ifd = fscanf(fd, "%lf", &(ns->ed.nif.par.A0));
        ifd = fscanf(fd, "%lf", &(ns->ed.nif.par.B0));
        ifd = fscanf(fd, "%lf", &(ns->ed.nif.par.c1));
        ifd = fscanf(fd, "%lf", &(ns->ed.nif.par.c2));
        fclose(fd);
        if (myrank == 0)
        {
            printf("=+=+=+= Prandtl number: %f =+=+=+=\n", ns->ed.nif.par.Pr);
            printf("=+=+=+= Brinkman number: %f =+=+=+=\n", ns->ed.nif.par.Br);
            printf("=+=+=+= Grashof number: %f =+=+=+=\n", ns->ed.nif.par.Gr);
            printf("=+=+=+= Galilei number: %f =+=+=+=\n", ns->ed.nif.par.Ga);
            printf("=+=+=+= Temperature of reference T0 (in Kelvin degrees): %f =+=+=+=\n", ns->ed.nif.par.T0);
            printf("=+=+=+= Temperature of reference T1 (in Kelvin degrees): %f =+=+=+=\n", ns->ed.nif.par.T1);


            //Arrhenius or modified Arrhenius model
            if ((ns->ed.nif.contr.ftmodel == 0) || (ns->ed.nif.contr.ftmodel == 1))
            {
            printf("=+=+=+= Alpha number: %f =+=+=+=\n", ns->ed.nif.par.alphaT);
            }
            //VFT model
            if ((ns->ed.nif.contr.ftmodel == 2))
            {
            printf("=+=+=+= Parameter A0: %f =+=+=+=\n", ns->ed.nif.par.A0);
            printf("=+=+=+= Parameter A0: %f =+=+=+=\n", ns->ed.nif.par.B0);
            }
            //WLF model
            if ((ns->ed.nif.contr.ftmodel == 3))
            {
            printf("=+=+=+= Parameter c1: %f =+=+=+=\n", ns->ed.nif.par.c1);
            printf("=+=+=+= Parameter c2: %f =+=+=+=\n", ns->ed.nif.par.c2);
            }
        }
    }
    else
    {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the non-isothermal flow parameters
void higflow_save_non_isothermal_flow_parameters(higflow_solver *ns, int myrank)
{
    if (myrank == 0)
    {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.nifpar", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL)
        {
            // Saving the parameters
            fprintf(fd, "%lf\n", (ns->ed.nif.par.Pr));
            fprintf(fd, "%lf\n", (ns->ed.nif.par.Br));
            fprintf(fd, "%lf\n", (ns->ed.nif.par.Gr));
            fprintf(fd, "%lf\n", (ns->ed.nif.par.Ga));
            fprintf(fd, "%lf\n", (ns->ed.nif.par.k0));
            fprintf(fd, "%lf\n", (ns->ed.nif.par.k1));
            fprintf(fd, "%lf\n", (ns->ed.nif.par.k2));
            fprintf(fd, "%lf\n", (ns->ed.nif.par.k3));
            fprintf(fd, "%lf\n", (ns->ed.nif.par.T0));
            fprintf(fd, "%lf\n", (ns->ed.nif.par.T1));
            fprintf(fd, "%lf\n", (ns->ed.nif.par.alphaT));
            fprintf(fd, "%lf\n", (ns->ed.nif.par.A0));
            fprintf(fd, "%lf\n", (ns->ed.nif.par.B0));
            fprintf(fd, "%lf\n", (ns->ed.nif.par.c1));
            fprintf(fd, "%lf\n", (ns->ed.nif.par.c2));
            fclose(fd);
        }
        else
        {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}


// Loading the non-isothermal controllers
void higflow_load_non_isothermal_flow_controllers(higflow_solver *ns, int myrank)
{
    // Parameters file name
    char namefile[1024];
    snprintf(namefile, sizeof namefile, "%s.nifcontr", ns->par.nameload);
    FILE *fd = fopen(namefile, "r");
    if (fd != NULL)
    {
        // Loading the parameters
        int ifd;
        ifd = fscanf(fd, "%d", &(ns->ed.nif.contr.bousinessq));
        ifd = fscanf(fd, "%d", &(ns->ed.nif.contr.visdissipation));
        ifd = fscanf(fd, "%d", &(ns->ed.nif.contr.thermodiscrtype));
        ifd = fscanf(fd, "%d", &(ns->ed.nif.contr.thermoconvecdiscrtype));
        ifd = fscanf(fd, "%d", &(ns->ed.nif.contr.ftmodel));
        fclose(fd);
        if (myrank == 0)
        {
            switch (ns->ed.nif.contr.bousinessq)
            {
            case 0:
                printf("=+=+=+= Bousinessq approximation: Deactivated =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Bousinessq approximation: Activated =+=+=+=\n");
                break;
            }
            switch (ns->ed.nif.contr.visdissipation)
            {
            case 0:
                printf("=+=+=+= Viscous dissipation: Deactivated =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Viscous dissipation: Activated =+=+=+=\n");
                break;
            }
            switch (ns->ed.nif.contr.thermodiscrtype)
            {
            case 0:
                printf("=+=+=+= Energy Equation Discretization: Explicit =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Energy Equation Discretization: Implicit =+=+=+=\n");
                break;
            }
            switch (ns->ed.nif.contr.thermoconvecdiscrtype)
            {
            case 0:
                printf("=+=+=+= Energy Equation Convective Term: Upwind =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Energy Equation Convective Term: CUBISTA =+=+=+=\n");
                break;
            }
            switch (ns->ed.nif.contr.ftmodel)
            {
            case 0:
                printf("=+=+=+= Function f(T): Arrhenius =+=+=+=\n");
                break;
            case 1:
                printf("=+=+=+= Function f(T): Modified-Arrhenius =+=+=+=\n");
                break;
            case 2:
                printf("=+=+=+= Function f(T): VFT =+=+=+=\n");
                break;
            case 3:
                printf("=+=+=+= Function f(T): WLF =+=+=+=\n");
                break;
            case 4:
                printf("=+=+=+= Function f(T): Constant =+=+=+=\n");
                break;
            }


        }
    }
    else
    {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n", namefile);
        exit(1);
    }
}

// Saving the non-isothermal controllers
void higflow_save_non_isothermal_flow_controllers(higflow_solver *ns, int myrank)
{
    if (myrank == 0)
    {
        // Parameters file name
        char namefile[1024];
        snprintf(namefile, sizeof namefile, "%s.nifcontr", ns->par.namesave);
        FILE *fd = fopen(namefile, "w");
        if (fd != NULL)
        {
            // Saving the parameters
            fprintf(fd, "%d\n", (ns->ed.nif.contr.bousinessq));
            fprintf(fd, "%d\n", (ns->ed.nif.contr.visdissipation));
            fprintf(fd, "%d\n", (ns->ed.nif.contr.thermodiscrtype));
            fprintf(fd, "%d\n", (ns->ed.nif.contr.thermoconvecdiscrtype));
            fprintf(fd, "%d\n", (ns->ed.nif.contr.ftmodel));
            fclose(fd);
        }
        else
        {
            // Error in open the file
            printf("=+=+=+= Error saving file %s =+=+=+=\n", namefile);
            exit(1);
        }
    }
}
