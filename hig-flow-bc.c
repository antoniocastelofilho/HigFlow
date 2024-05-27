// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Boundary Condition - version 26/05/2021
// *******************************************************************
// *******************************************************************

#include "hig-flow-bc.h"
#include <string.h>
#include <libfyaml.h>

// Make the boundary condition
sim_boundary *higflow_make_bc(hig_cell *bcg, bc_type type, int id, bc_valuetype valuetype) {
    // Setting the map
    mp_mapper *bm = mp_create();
    higcit_celliterator *it;
    it = higcit_create_all_leaves(bcg);
    mp_assign_from_celliterator(bm, it, 0);
    higcit_destroy(it);
    // Create the bondary contition
    sim_boundary *bc = sb_create(bcg, type, bm);
    // Set the identifier boundary condition
    sb_set_userid(bc, id);
    // Set the valuetype of the boundary condition
    sb_set_valuetype(bc, valuetype);
    return bc;
}

// Creating and setting the boundary condition for the pressure
void higflow_set_boundary_condition_for_pressure(higflow_solver *ns, int numbcs, int id[numbcs], char bcfilenames[numbcs][1024], bc_type pbctypes[], bc_valuetype pbcvaluetype[]) {
    // Get the HigTree from amr file
    hig_cell *bcg[numbcs];
    for(int h = 0; h < numbcs; h++) {
        FILE *fd = fopen(bcfilenames[h], "r");
        bcg[h] = higio_read_from_amr(fd);
        fclose(fd);
    }
    // Loop for each boundary condition
    for(int h = 0; h < numbcs; h++) {
        // Get the local domain for cell center
        sim_domain *sd = psd_get_local_domain(ns->psdp);
        // Create the bounary condition
        sim_boundary *bc = higflow_make_bc(bcg[h], pbctypes[h], id[h], pbcvaluetype[h]);
        // Adding the boundary condition 
        sd_add_boundary(sd, bc);
        // Get the mapper for the boundary condition
        mp_mapper *bm = sb_get_mapper(bc);
        // Loop for the cells of the boundaries conditions 
        higcit_celliterator *it;
        for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell 
            hig_cell *bcell = higcit_getcell(it);
            // Get the cell center
            Point bccenter;
            hig_get_center(bcell, bccenter);
            // Get the id of the cell
            int bclid = mp_lookup(bm, hig_get_cid(bcell));
            // Set the time to get the pressure
            real t = ns->par.t + ns->par.dt;
            // Get the pressure defined by the user
            real val = ns->func.get_boundary_pressure(id[h], bccenter, ns->par.t);
            sb_set_value(bc, bclid, val);
        }
        // Destroy the iterator
        higcit_destroy(it);
    }
}

// Creating and setting the boundary condition for the velocity
void higflow_set_boundary_condition_for_velocities(higflow_solver *ns, int numbcs, int id[numbcs], char bcfilenames[numbcs][1024], bc_type bctypes[DIM][numbcs], bc_valuetype bcvaluetype[DIM][numbcs]) {
    // Get the HigTree from amr file
    hig_cell *bcg[numbcs];
    for(int h = 0; h < numbcs; h++) {
        FILE *fd = fopen(bcfilenames[h], "r");
        bcg[h] = higio_read_from_amr(fd);
        fclose(fd);
    }
    sim_facet_domain *sfd;
    // Loop for the dimension
    for(int dim = 0; dim < DIM; dim++) {
        // Loop for the boundaries conditions
        for(int h = 0; h < numbcs; h++) {
            // Get the local domain for the facet
            sfd = psfd_get_local_domain(ns->psfdu[dim]);
            // Create the bounary condition
            sim_boundary *bc = higflow_make_bc(bcg[h], bctypes[dim][h], id[h], bcvaluetype[dim][h]);
            // Adding the boundary condition 
            sfd_add_boundary(sfd, bc);
            // Get the map
            mp_mapper *bm = sb_get_mapper(bc);
            // loop for the facet of the cells of the boundaries conditions 
            higcit_celliterator *it;
            for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
                // Getting the cell 
                hig_cell *bcell = higcit_getcell(it);
                // Get the cell center
                Point bccenter;
                hig_get_center(bcell, bccenter);
                // Set the value bcval for the center of the facet cell 
                int  bclid = mp_lookup(bm, hig_get_cid(bcell));
                //real bcval = bcvalues[dim][h];
                // Get the velocity defined by the user
                real bcval = ns->func.get_boundary_velocity(id[h], bccenter, dim, ns->par.t);
                // Set the value pbcvalues for the center of the cell 
                sb_set_value(bc, bclid, bcval);
            }
            // Destroying the iterator 
            higcit_destroy(it);
        }
        /////////// Already runs on kernel.c
        // // Mapping the properties in the domain (velocities)
        // psfd_compute_sfbi(ns->psfdu[dim]);
        // // Sync mapper for velocities
        // psfd_synced_mapper(ns->psfdu[dim]);
    }
}

// Creating and setting the boundary condition for the electro-osmotic source term
void higflow_set_boundary_condition_for_electroosmotic_source_term(higflow_solver *ns, int numbcs, int id[numbcs], char bcfilenames[numbcs][1024], bc_type bctypes[DIM][numbcs], bc_valuetype bcvaluetype[DIM][numbcs]) {    
    if (ns->contr.eoflow == true || (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true)) {
        real bcval, fracvol;

        // Get the HigTree from amr file
        hig_cell *bcg[numbcs];
        for(int h = 0; h < numbcs; h++) {
            FILE *fd = fopen(bcfilenames[h], "r");
            bcg[h] = higio_read_from_amr(fd);
            fclose(fd);
        }
        sim_facet_domain *sfd;
        // Loop for the dimension
        for(int dim = 0; dim < DIM; dim++) {
            // Loop for the boundaries conditions
            for(int h = 0; h < numbcs; h++) {
                // Get the local domain for the facet
                sfd = psfd_get_local_domain(ns->ed.eo.psfdEOFeo[dim]);
                // Create the bounary condition
                sim_boundary *bc = higflow_make_bc(bcg[h], DIRICHLET, id[h], timedependent);
                // Adding the boundary condition 
                sfd_add_boundary(sfd, bc);
                // Get the map
                mp_mapper *bm = sb_get_mapper(bc);
                // loop for the facet of the cells of the boundaries conditions 
                higcit_celliterator *it;
                for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
                    // Getting the cell 
                    hig_cell *bcell = higcit_getcell(it);
                    // Get the cell center
                    Point bccenter;
                    hig_get_center(bcell, bccenter);
                    // Set the value bcval for the center of the facet cell 
                    int  bclid = mp_lookup(bm, hig_get_cid(bcell));
                    //real bcval = bcvalues[dim][h];
                    // Get the eo source term defined by the user
                    if (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true) {
                        fracvol = compute_value_at_point(ns->ed.mult.sdmult, bccenter, bccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
                        bcval = ns->ed.mult.eo.get_boundary_multiphase_electroosmotic_source_term(fracvol, id[h], bccenter, dim, ns->par.t);
                    } else
                        bcval = ns->ed.eo.get_boundary_electroosmotic_source_term(id[h], bccenter, dim, ns->par.t);

                    // Set the value pbcvalues for the center of the cell 
                    sb_set_value(bc, bclid, bcval);
                }
                // Destroying the iterator 
                higcit_destroy(it);
            }
            /////////// Already runs on kernel.c
            // // Mapping the properties in the domain (velocities)
            // psfd_compute_sfbi(ns->ed.eo.psfdEOFeo[dim]);
            // // Sync mapper for velocities
            // psfd_synced_mapper(ns->ed.eo.psfdEOFeo[dim]);
        }

    }
}

// Creating and setting the boundary condition for the electro-osmotic phi
void higflow_set_boundary_condition_for_electroosmotic_phi(higflow_solver *ns, int numbcs, int id[numbcs], char bcfilenames[numbcs][1024], bc_type phibctypes[], bc_valuetype phibcvaluetype[]) {  
    if (ns->contr.eoflow == true || (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true)) {
        real bcval, fracvol;

        // Get the HigTree from amr file
        hig_cell *bcg[numbcs];
        for(int h = 0; h < numbcs; h++) {
            FILE *fd = fopen(bcfilenames[h], "r");
            bcg[h] = higio_read_from_amr(fd);
            fclose(fd);
        }
        // Loop for each boundary condition
        for(int h = 0; h < numbcs; h++) {
            // Get the local domain for cell center
            sim_domain *sd = psd_get_local_domain(ns->ed.eo.psdEOphi);
            // Create the bounary condition
            sim_boundary *bc = higflow_make_bc(bcg[h], phibctypes[h], id[h], phibcvaluetype[h]);
            // Adding the boundary condition 
            sd_add_boundary(sd, bc);
            // Get the mapper for the boundary condition
            mp_mapper *bm = sb_get_mapper(bc);
            // Loop for the cells of the boundaries conditions 
            higcit_celliterator *it;
            for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
                // Get the cell 
                hig_cell *bcell = higcit_getcell(it);
                // Get the cell center
                Point bccenter;
                hig_get_center(bcell, bccenter);
                // Get the id of the cell
                int bclid = mp_lookup(bm, hig_get_cid(bcell));
                // Set the time to get the pressure
                real t = ns->par.t + ns->par.dt;
                // Get phi defined by the user
                if (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true) {
                    fracvol = compute_value_at_point(ns->ed.mult.sdmult, bccenter, bccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
                    bcval = ns->ed.mult.eo.get_boundary_multiphase_electroosmotic_phi(fracvol, id[h], bccenter, ns->par.t);
                } else
                    bcval = ns->ed.eo.get_boundary_electroosmotic_phi(id[h], bccenter, ns->par.t);
                // Set the value 
                sb_set_value(bc, bclid, bcval);
            }
            // Destroy the iterator
            higcit_destroy(it);
        }

    }
}

// Creating and setting the boundary condition for the electro-osmotic psi
void higflow_set_boundary_condition_for_electroosmotic_psi(higflow_solver *ns, int numbcs, int id[numbcs], char bcfilenames[numbcs][1024], bc_type psibctypes[], bc_valuetype psibcvaluetype[]) {
    if (ns->contr.eoflow == true || (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true)) {
        real bcval, fracvol;
    
        // Get the HigTree from amr file
        hig_cell *bcg[numbcs];
        for(int h = 0; h < numbcs; h++) {
            FILE *fd = fopen(bcfilenames[h], "r");
            bcg[h] = higio_read_from_amr(fd);
            fclose(fd);
        }
        // Loop for each boundary condition
        for(int h = 0; h < numbcs; h++) {
            // Get the local domain for cell center
            sim_domain *sd = psd_get_local_domain(ns->ed.eo.psdEOpsi);
            // Create the bounary condition
            sim_boundary *bc = higflow_make_bc(bcg[h], psibctypes[h], id[h], psibcvaluetype[h]);
            // Adding the boundary condition 
            sd_add_boundary(sd, bc);
            // Get the mapper for the boundary condition
            mp_mapper *bm = sb_get_mapper(bc);
            // Loop for the cells of the boundaries conditions 
            higcit_celliterator *it;
            for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
                // Get the cell 
                hig_cell *bcell = higcit_getcell(it);
                // Get the cell center
                Point bccenter;
                hig_get_center(bcell, bccenter);
                // Get the id of the cell
                int bclid = mp_lookup(bm, hig_get_cid(bcell));
                // Set the time to get the pressure
                real t = ns->par.t + ns->par.dt;
                // Get psi defined by the user
                if (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true) {
                    fracvol = compute_value_at_point(ns->ed.mult.sdmult, bccenter, bccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
                    bcval = ns->ed.mult.eo.get_boundary_multiphase_electroosmotic_psi(fracvol, id[h], bccenter, ns->par.t);
                } else
                    bcval = ns->ed.eo.get_boundary_electroosmotic_psi(id[h], bccenter, ns->par.t);
                // Set the value 
                sb_set_value(bc, bclid, bcval);
            }
            // Destroy the iterator
            higcit_destroy(it);
        }

    }
}

// Creating and setting the boundary condition for the electro-osmotic nplus
void higflow_set_boundary_condition_for_electroosmotic_nplus(higflow_solver *ns, int numbcs, int id[numbcs], char bcfilenames[numbcs][1024], bc_type nplusbctypes[], bc_valuetype nplusbcvaluetype[]) {
    if (ns->contr.eoflow == true || (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true)) {
        real bcval, fracvol;
    
        // Get the HigTree from amr file
        hig_cell *bcg[numbcs];
        for(int h = 0; h < numbcs; h++) {
            FILE *fd = fopen(bcfilenames[h], "r");
            bcg[h] = higio_read_from_amr(fd);
            fclose(fd);
        }
        // Loop for each boundary condition
        for(int h = 0; h < numbcs; h++) {
            // Get the local domain for cell center
            sim_domain *sd = psd_get_local_domain(ns->ed.eo.psdEOnplus);
            // Create the bounary condition
            sim_boundary *bc = higflow_make_bc(bcg[h], nplusbctypes[h], id[h], nplusbcvaluetype[h]);
            // Adding the boundary condition 
            sd_add_boundary(sd, bc);
            // Get the mapper for the boundary condition
            mp_mapper *bm = sb_get_mapper(bc);
            // Loop for the cells of the boundaries conditions 
            higcit_celliterator *it;
            for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
                // Get the cell 
                hig_cell *bcell = higcit_getcell(it);
                // Get the cell center
                Point bccenter;
                hig_get_center(bcell, bccenter);
                // Get the id of the cell
                int bclid = mp_lookup(bm, hig_get_cid(bcell));
                // Set the time to get the pressure
                real t = ns->par.t + ns->par.dt;
                // Get the pressure defined by the user
                if (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true) {
                    fracvol = compute_value_at_point(ns->ed.mult.sdmult, bccenter, bccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
                    bcval = ns->ed.mult.eo.get_boundary_multiphase_electroosmotic_nplus(fracvol, id[h], bccenter, ns->par.t);
                } else
                    bcval = ns->ed.eo.get_boundary_electroosmotic_nplus(id[h], bccenter, ns->par.t);
                // Set the value 
                sb_set_value(bc, bclid, bcval);
            }
            // Destroy the iterator
            higcit_destroy(it);
        }

    }
}

// Creating and setting the boundary condition for the electro-osmotic nminus
void higflow_set_boundary_condition_for_electroosmotic_nminus(higflow_solver *ns, int numbcs, int id[numbcs], char bcfilenames[numbcs][1024], bc_type nminusbctypes[], bc_valuetype nminusbcvaluetype[]) {
    if (ns->contr.eoflow == true || (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true)) {
        real bcval, fracvol;
    
        // Get the HigTree from amr file
        hig_cell *bcg[numbcs];
        for(int h = 0; h < numbcs; h++) {
            FILE *fd = fopen(bcfilenames[h], "r");
            bcg[h] = higio_read_from_amr(fd);
            fclose(fd);
        }
        // Loop for each boundary condition
        for(int h = 0; h < numbcs; h++) {
            // Get the local domain for cell center
            sim_domain *sd = psd_get_local_domain(ns->ed.eo.psdEOnminus);
            // Create the bounary condition
            sim_boundary *bc = higflow_make_bc(bcg[h], nminusbctypes[h], id[h], nminusbcvaluetype[h]);
            // Adding the boundary condition 
            sd_add_boundary(sd, bc);
            // Get the mapper for the boundary condition
            mp_mapper *bm = sb_get_mapper(bc);
            // Loop for the cells of the boundaries conditions 
            higcit_celliterator *it;
            for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
                // Get the cell 
                hig_cell *bcell = higcit_getcell(it);
                // Get the cell center
                Point bccenter;
                hig_get_center(bcell, bccenter);
                // Get the id of the cell
                int bclid = mp_lookup(bm, hig_get_cid(bcell));
                // Set the time to get the pressure
                real t = ns->par.t + ns->par.dt;
                // Get the pressure defined by the user
                if (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true) {
                    fracvol = compute_value_at_point(ns->ed.mult.sdmult, bccenter, bccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
                    bcval = ns->ed.mult.eo.get_boundary_multiphase_electroosmotic_nminus(fracvol, id[h], bccenter, ns->par.t);
                } else
                    bcval = ns->ed.eo.get_boundary_electroosmotic_nminus(id[h], bccenter, ns->par.t);
                // Set the value 
                sb_set_value(bc, bclid, bcval);
            }
            // Destroy the iterator
            higcit_destroy(it);
        }

    }
}

// Creating and setting the boundary condition for cell source term
void higflow_set_boundary_condition_for_cell_source_term(higflow_solver *ns, int numbcs, int id[numbcs], char bcfilenames[numbcs][1024], bc_type pbctypes[], bc_valuetype pbcvaluetype[]) {
    // Get the HigTree from amr file
    hig_cell *bcg[numbcs];
    for(int h = 0; h < numbcs; h++) {
        FILE *fd = fopen(bcfilenames[h], "r");
        bcg[h] = higio_read_from_amr(fd);
        fclose(fd);
    }
    // Loop for each boundary condition
    for(int h = 0; h < numbcs; h++) {
        // Get the local domain for cell center
        sim_domain *sd = psd_get_local_domain(ns->psdF);
        // Create the bounary condition
        sim_boundary *bc = higflow_make_bc(bcg[h], DIRICHLET, id[h], timedependent);
        // Adding the boundary condition 
        sd_add_boundary(sd, bc);
        // Get the mapper for the boundary condition
        mp_mapper *bm = sb_get_mapper(bc);
        // Loop for the cells of the boundaries conditions 
        higcit_celliterator *it;
        for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell 
            hig_cell *bcell = higcit_getcell(it);
            // Get the cell center
            Point bccenter;
            hig_get_center(bcell, bccenter);
            // Get the id of the cell
            int bclid = mp_lookup(bm, hig_get_cid(bcell));
            // Get the electro-osmotic nminus defined by the user
            real val = ns->func.get_boundary_source_term(id[h], bccenter, ns->par.t);
            // Set the value 
            sb_set_value(bc, bclid, val);
        }
        // Destroy the iterator
        higcit_destroy(it);
    }
}

// Creating and setting the boundary condition for facet source term
void higflow_set_boundary_condition_for_facet_source_term(higflow_solver *ns, int numbcs, int id[numbcs], char bcfilenames[numbcs][1024], bc_type bctypes[DIM][numbcs], bc_valuetype bcvaluetype[DIM][numbcs]) {
    // Get the HigTree from amr file
    hig_cell *bcg[numbcs];
    for(int h = 0; h < numbcs; h++) {
        FILE *fd = fopen(bcfilenames[h], "r");
        bcg[h] = higio_read_from_amr(fd);
        fclose(fd);
    }
    sim_facet_domain *sfdF;
    // Loop for the dimension
    for(int dim = 0; dim < DIM; dim++) {
        // Loop for the boundaries conditions
        for(int h = 0; h < numbcs; h++) {
            // Get the local domain for the facet
            sfdF = psfd_get_local_domain(ns->psfdF[dim]);
            // Create the bounary condition
            sim_boundary *bc = higflow_make_bc(bcg[h], DIRICHLET, id[h], timedependent);
            // Adding the boundary condition 
            sfd_add_boundary(sfdF, bc);
            // Get the map
            mp_mapper *bm = sb_get_mapper(bc);
            // loop for the facet of the cells of the boundaries conditions 
            higcit_celliterator *it;
            for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
                // Getting the cell 
                hig_cell *bcell = higcit_getcell(it);
                // Get the cell center
                Point bccenter;
                hig_get_center(bcell, bccenter);
                // Set the value bcval for the center of the facet cell 
                int  bclid = mp_lookup(bm, hig_get_cid(bcell));
                //real bcval = bcvalues[dim][h];
                // Get the velocity defined by the user
                real bcval = ns->func.get_boundary_facet_source_term(id[h], bccenter, dim, ns->par.t);
                // Set the value pbcvalues for the center of the cell 
                sb_set_value(bc, bclid, bcval);
            }
            // Destroying the iterator 
            higcit_destroy(it);
        }
        /////////// Already runs on kernel.c
        // // // Mapping the properties in the domain (velocities)
        // psfd_compute_sfbi(ns->psfdF[dim]);
        // // Sync mapper for velocities
        // psfd_synced_mapper(ns->psfdF[dim]);
    }
}

// Navier-Stokes initialize the domain and boudaries
void higflow_initialize_boundaries(higflow_solver *ns) {
    // Loading the boundary condition data
    char namefile[1024];
    sprintf(namefile,"%s.bc",ns->par.nameload);
    FILE *fbc = fopen(namefile, "r");
    if (fbc == NULL) {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n",namefile);
        exit(1);
    }
    // Number of boundaries
    int numbcs; 
    int ifd = fscanf(fbc,"%d\n",&numbcs);
    // Boudary condition data
    int           id[numbcs];
    char          amrBCfilename[numbcs][1024]; 
    bc_type       pbctypes[numbcs]; 
    bc_type       ubctypes[DIM][numbcs]; 
    bc_valuetype  pbcvaluetype[numbcs];
    bc_valuetype  ubcvaluetype[DIM][numbcs]; 
    // Setting the pressure desingularizadtion control
    ns->contr.desingpressure = true;
    for(int h = 0; h < numbcs; h++) {
        ifd = fscanf(fbc,"%d\n",&(id[h]));
        // HigTree Boundary condition file name
        __higflow_readstring(amrBCfilename[h],1024,fbc);
        // Pressure boundary condition type
        ifd = fscanf(fbc,"%d",(int *)&(pbctypes[h]));
        // Pressure boundary condition value
        ifd = fscanf(fbc,"%d",(int *)&(pbcvaluetype[h]));
        // Setting the pressure desingularizadtion control
        if (pbctypes[h] == DIRICHLET) {
            // Outflow
            ns->contr.desingpressure = false;
        }
        
        for (int dim = 0; dim < DIM; dim++) {
            // Velocity boundary condition type
            ifd = fscanf(fbc,"%d",(int *)&(ubctypes[dim][h]));
            // Velocity boundary condition valuetype
            ifd = fscanf(fbc,"%d",(int *)&(ubcvaluetype[dim][h]));
        }
    }
    fclose(fbc);
    // Setting the boundary conditions for the pressure
    higflow_set_boundary_condition_for_pressure(ns, numbcs, id, amrBCfilename, pbctypes, pbcvaluetype);
    // Setting the boundary conditions for the velocities
    higflow_set_boundary_condition_for_velocities(ns, numbcs, id, amrBCfilename, ubctypes, ubcvaluetype);
    // Setting the boundary conditions for the cell source term 
    higflow_set_boundary_condition_for_cell_source_term(ns, numbcs, id, amrBCfilename, pbctypes, pbcvaluetype);
    // Setting the boundary conditions for the facet source term
    higflow_set_boundary_condition_for_facet_source_term(ns, numbcs, id, amrBCfilename, ubctypes, ubcvaluetype);
    
    // Setting the boundary conditions for the electro-osmotic model
    if (ns->contr.eoflow == true || (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true)) {
        sprintf(namefile,"%s.bcelectroosmotic",ns->par.nameload);
        //printf("nome: %s", namefile);
        FILE *fbc = fopen(namefile, "r");
        if (fbc == NULL) {
            // Error in open the file
            printf("=+=+=+= Error loading file %s =+=+=+=\n",namefile);
            exit(1);
        }
        // Number of boundaries
        int numbcs; 
        int ifd = fscanf(fbc,"%d\n",&numbcs);
        // Boudary condition data
        int           phibc_timedependent = 0;
        int           psibc_timedependent = 0;
        int           id[numbcs];
        char          amrBCfilename[numbcs][1024]; 
        bc_type       phibctypes[numbcs]; 
        bc_valuetype  phibcvaluetype[numbcs];
        bc_type       psibctypes[numbcs]; 
        bc_valuetype  psibcvaluetype[numbcs];
        bc_type       nplusbctypes[numbcs]; 
        bc_valuetype  nplusbcvaluetype[numbcs];
        bc_type       nminusbctypes[numbcs]; 
        bc_valuetype  nminusbcvaluetype[numbcs];
        for(int h = 0; h < numbcs; h++) {
            ifd = fscanf(fbc,"%d\n",&(id[h]));
            // HigTree Boundary condition file name
            __higflow_readstring(amrBCfilename[h],1024,fbc);
            // Electro-osmotic potential boundary condition type
            ifd = fscanf(fbc,"%d",(int *)&(phibctypes[h]));
            // Electro-osmotic potential boundary condition value
            ifd = fscanf(fbc,"%d",(int *)&(phibcvaluetype[h]));
            if(phibcvaluetype[h] == timedependent) phibc_timedependent = true;

            // Electro-osmotic potential boundary condition type
            ifd = fscanf(fbc,"%d",(int *)&(psibctypes[h]));
            // Electro-osmotic potential boundary condition value
            ifd = fscanf(fbc,"%d",(int *)&(psibcvaluetype[h]));
            if(psibcvaluetype[h] == timedependent) psibc_timedependent = true;
    
            // Electro-osmotic concentration n+ boundary condition type
            ifd = fscanf(fbc,"%d",(int *)&(nplusbctypes[h]));
            // Electro-osmotic concentration n+ boundary condition value
            ifd = fscanf(fbc,"%d",(int *)&(nplusbcvaluetype[h]));
          
            // Electro-osmotic concentration n- boundary condition type
            ifd = fscanf(fbc,"%d",(int *)&(nminusbctypes[h]));
            // Electro-osmotic concentration n- boundary condition value
            ifd = fscanf(fbc,"%d",(int *)&(nminusbcvaluetype[h]));
        }
        fclose(fbc);

        // Setting the boundary conditions for the electro-osmotic phi
        higflow_set_boundary_condition_for_electroosmotic_phi(ns, numbcs, id, amrBCfilename, phibctypes, phibcvaluetype);
        // // Setting the boundary conditions for the electro-osmotic psi
        higflow_set_boundary_condition_for_electroosmotic_psi(ns, numbcs, id, amrBCfilename, psibctypes, psibcvaluetype);
        // // Setting the boundary conditions for the electro-osmotic nplus
        higflow_set_boundary_condition_for_electroosmotic_nplus(ns, numbcs, id, amrBCfilename, nplusbctypes, nplusbcvaluetype);
        // // Setting the boundary conditions for the electro-osmotic nminus
        higflow_set_boundary_condition_for_electroosmotic_nminus(ns, numbcs, id, amrBCfilename, nminusbctypes, nminusbcvaluetype);
        // Setting the boundary conditions for the electro-osmotic source term
        higflow_set_boundary_condition_for_electroosmotic_source_term(ns, numbcs, id, amrBCfilename, ubctypes, ubcvaluetype);

        eo_controllers *eo_contr;
        if (ns->contr.eoflow == true) eo_contr = &(ns->ed.eo.contr);
        else eo_contr = &(ns->ed.mult.eo.contr);
        
        eo_contr->is_phibc_timedependent = phibc_timedependent;
        if(eo_contr->is_phibc_timedependent == true){
            print0f("=+=+=+= Applied potential phi boundary conditions are time dependent =+=+=+=\n");
        } else {
            print0f("=+=+=+= Applied potential phi boundary conditions are time independent - Will only be calculated once in the Laplace equation=+=+=+=\n");
        }
        eo_contr->is_psibc_timedependent = psibc_timedependent;
        if(eo_contr->is_psibc_timedependent == true){
            print0f("=+=+=+= Electro-osmotic potential psi boundary conditions are time dependent =+=+=+=\n");
        } else {
            print0f("=+=+=+= Electro-osmotic potential psi boundary conditions are time independent ");
            if(eo_contr->eo_model == PB || eo_contr->eo_model == PBDH || eo_contr->eo_model == PBDH_ANALYTIC){
                print0f("- Will only be calculated once in the poisson equation");
            }
            print0f("=+=+=+=\n");
        }
    }
}
// Navier-Stokes initialize the domain and boudaries
void higflow_initialize_boundaries_yaml(higflow_solver *ns) {
    // Loading the boundary condition data
    char namefile[1024];
    sprintf(namefile,"%s.bc.yaml",ns->par.nameload);
    
    FILE *fbc = fopen(namefile, "r");
    struct fy_document *fyd = NULL;
    fyd = fy_document_build_from_file(NULL, namefile);
     
    if (fyd == NULL) {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n",namefile);
        exit(1);
    }
     
    // Number of boundaries
    int numbcs; 
    int ifd = fy_document_scanf(fyd,"/bc/number_bc %d",&numbcs);
 
    // Boudary condition data
    int           id[numbcs];
    char          amrBCfilename[numbcs][1024]; 
    bc_type       pbctypes[numbcs]; 
    bc_type       ubctypes[DIM][numbcs]; 
    bc_valuetype  pbcvaluetype[numbcs];
    bc_valuetype  ubcvaluetype[DIM][numbcs]; 
    // Setting the pressure desingularizadtion control
    ns->contr.desingpressure = true;
    
    for(int h = 0; h < numbcs; h++) {
        char atrib[1024];
        sprintf(atrib,"/bc/bc%d/id %%d",h);
        ifd = fy_document_scanf(fyd,atrib,&(id[h]));
        // HigTree Boundary condition file name
        sprintf(atrib,"/bc/bc%d/path %%s",h);
        ifd = fy_document_scanf(fyd,atrib,amrBCfilename[h]);
        //__higflow_readstring(amrBCfilename[h],1024,fbc);
        // Pressure boundary condition type
        char aux[1024];
        sprintf(atrib,"/bc/bc%d/pressure/type %%s",h);
        ifd = fy_document_scanf(fyd,atrib,aux);
        if (strcmp(aux,"dirichlet") == 0) {
           pbctypes[h] = DIRICHLET;
        } else if (strcmp(aux,"neumann") == 0) {
            pbctypes[h] = NEUMANN;
        } else {
            printf("=+=+=+= Error loading boundary condition type for the pressure in the boundary %d (may not be implemented yet) =+=+=+= \n",h);
            exit(1);
        }
        // Pressure boundary condition value
        sprintf(atrib,"/bc/bc%d/pressure/value_type %%s",h);
        ifd = fy_document_scanf(fyd,atrib,aux);
        if (strcmp(aux,"fixed_value") == 0) {
            pbcvaluetype[h] = fixedValue;
        } else if(strcmp(aux,"time_dependent") == 0) {
            pbcvaluetype[h] = timedependent;
        } 
           // Not implemented yet
        // else if(strcmp(aux,"zero_gradient") == 0) {
        //     pbcvaluetype[h] = zeroGradient;
        //  } else if(strcmp(aux,"freestream") == 0) {
        //     pbcvaluetype[h] = freestream;
        //  } else if(strcmp(aux,"empty") == 0) {
        //     pbcvaluetype[h] = empty;
        //  } else if(strcmp(aux,"outflow") == 0) {
        //      pbcvaluetype[h] = outflow;
        //  }
        else {
            printf("=+=+=+= Error loading boundary condition valuetype for the pressure in the boundary %d (may not be implemented yet) =+=+=+= \n",h);
            exit(1);
        }
        // Setting the pressure desingularizadtion control
        if (pbctypes[h] == DIRICHLET) {
            // Outflow
            //sprintf(atrib,"/boundary_conditions/bc%d/pressure/value_type/value %%lf",h);
            //ifd = fy_document_scanf(fyd,atrib,&(ns->contr.desingpressure));
            ns->contr.desingpressure = false;
        }
        // Velocity boundary condition value
        for (int dim = 0; dim < DIM; dim++) {
            sprintf(atrib,"/bc/bc%d/velocity_%d/type %%s",h,dim);
            ifd = fy_document_scanf(fyd,atrib,aux);
            if (strcmp(aux,"dirichlet") == 0) {
               ubctypes[dim][h] = DIRICHLET;
            } else if(strcmp(aux,"neumann") == 0) {
               ubctypes[dim][h] = NEUMANN;
            } else {
               printf("=+=+=+= Error loading boundary condition type for the velocity in the boundary %d (may not be implemented yet) =+=+=+= \n",h);
               exit(1);
            }
            sprintf(atrib,"/bc/bc%d/velocity_%d/value_type %%s",h,dim);
            ifd = fy_document_scanf(fyd,atrib,aux);
            // Velocity boundary condition valuetype
            if (strcmp(aux,"fixed_value") == 0) {
               ubcvaluetype[dim][h] = fixedValue;
            } else if(strcmp(aux,"time_dependent") == 0) {
               ubcvaluetype[dim][h] = timedependent;
            } else {
               printf("=+=+=+= Error loading boundary condition valuetype for the velocity in the boundary %d (may not be implemented yet) =+=+=+= \n",h);
               exit(1);
            }
        }
    }

    // Setting the boundary conditions for the pressure
    higflow_set_boundary_condition_for_pressure(ns, numbcs, id, amrBCfilename, pbctypes, pbcvaluetype);
    // Setting the boundary conditions for the velocities
    higflow_set_boundary_condition_for_velocities(ns, numbcs, id, amrBCfilename, ubctypes, ubcvaluetype);
    // Setting the boundary conditions for the cell source term 
    higflow_set_boundary_condition_for_cell_source_term(ns, numbcs, id, amrBCfilename, pbctypes, pbcvaluetype);
    // Setting the boundary conditions for the facet source term
    higflow_set_boundary_condition_for_facet_source_term(ns, numbcs, id, amrBCfilename, ubctypes, ubcvaluetype);

    // Setting the boundary conditions for the electro-osmotic model
    if (ns->contr.eoflow == true || (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true)) {
        int ifd = fy_document_scanf(fyd,"/bc_electroosmotic/number_bc %d",&numbcs);
 
        // Boudary condition data
        int           phibc_timedependent = 0;
        int           psibc_timedependent = 0;
        int           id[numbcs];
        char          amrBCfilename[numbcs][1024]; 
        bc_type       phibctypes[numbcs]; 
        bc_valuetype  phibcvaluetype[numbcs];
        bc_type       psibctypes[numbcs]; 
        bc_valuetype  psibcvaluetype[numbcs];
        bc_type       nplusbctypes[numbcs]; 
        bc_valuetype  nplusbcvaluetype[numbcs];
        bc_type       nminusbctypes[numbcs]; 
        bc_valuetype  nminusbcvaluetype[numbcs];
        
        for(int h = 0; h < numbcs; h++) {
            char atrib[1024];
            sprintf(atrib,"/bc_electroosmotic/bc%d/id %%d",h);
            ifd = fy_document_scanf(fyd,atrib,&(id[h]));
            // HigTree Boundary condition file name
            sprintf(atrib,"/bc_electroosmotic/bc%d/path %%s",h);
            ifd = fy_document_scanf(fyd,atrib,amrBCfilename[h]);
            // Applied Potential boundary condition type
            char aux[1024];
            sprintf(atrib,"/bc_electroosmotic/bc%d/phi/type %%s",h);
            ifd = fy_document_scanf(fyd,atrib,aux);
            if (strcmp(aux,"dirichlet") == 0) {
                phibctypes[h] = DIRICHLET;
            } else if (strcmp(aux,"neumann") == 0) {
                phibctypes[h] = NEUMANN;
            } else {
                printf("=+=+=+= Error loading boundary condition type for phi in the boundary %d (may not be implemented yet) =+=+=+= \n",h);
                exit(1);
            }
            // Applied Potential boundary condition value type
            sprintf(atrib,"/bc_electroosmotic/bc%d/phi/value_type %%s",h);
            ifd = fy_document_scanf(fyd,atrib,aux);
            if (strcmp(aux,"fixed_value") == 0) {
                phibcvaluetype[h] = fixedValue;
            } else if(strcmp(aux,"time_dependent") == 0) {
                phibcvaluetype[h] = timedependent;
            } else {
                printf("=+=+=+= Error loading boundary condition valuetype for phi in the boundary %d (may not be implemented yet) =+=+=+= \n",h);
                exit(1);
            }

            // Zeta Potential boundary condition type
            sprintf(atrib,"/bc_electroosmotic/bc%d/psi/type %%s",h);
            ifd = fy_document_scanf(fyd,atrib,aux);
            if (strcmp(aux,"dirichlet") == 0) {
                psibctypes[h] = DIRICHLET;
            } else if (strcmp(aux,"neumann") == 0) {
                psibctypes[h] = NEUMANN;
            } else {
                printf("=+=+=+= Error loading boundary condition type for psi in the boundary %d (may not be implemented yet) =+=+=+= \n",h);
                exit(1);
            }
            // Zeta Potential boundary condition value type
            sprintf(atrib,"/bc_electroosmotic/bc%d/psi/value_type %%s",h);
            ifd = fy_document_scanf(fyd,atrib,aux);
            if (strcmp(aux,"fixed_value") == 0) {
                psibcvaluetype[h] = fixedValue;
            } else if(strcmp(aux,"time_dependent") == 0) {
                psibcvaluetype[h] = timedependent;
            } else {
                printf("=+=+=+= Error loading boundary condition valuetype for psi in the boundary %d (may not be implemented yet) =+=+=+= \n",h);
                exit(1);
            }

            // Positive Charge Concentration boundary condition type
            sprintf(atrib,"/bc_electroosmotic/bc%d/nplus/type %%s",h);
            ifd = fy_document_scanf(fyd,atrib,aux);
            if (strcmp(aux,"dirichlet") == 0) {
                nplusbctypes[h] = DIRICHLET;
            } else if (strcmp(aux,"neumann") == 0) {
                nplusbctypes[h] = NEUMANN;
            } else {
                printf("=+=+=+= Error loading boundary condition type for nplus in the boundary %d (may not be implemented yet) =+=+=+= \n",h);
                exit(1);
            }
            // Positive Charge Concentration boundary condition value type
            sprintf(atrib,"/bc_electroosmotic/bc%d/nplus/value_type %%s",h);
            ifd = fy_document_scanf(fyd,atrib,aux);
            if (strcmp(aux,"fixed_value") == 0) {
                nplusbcvaluetype[h] = fixedValue;
            } else if(strcmp(aux,"time_dependent") == 0) {
                nplusbcvaluetype[h] = timedependent;
            } else {
                printf("=+=+=+= Error loading boundary condition valuetype for nplus in the boundary %d (may not be implemented yet) =+=+=+= \n",h);
                exit(1);
            }

            // Negative Charge Concentration boundary condition type
            sprintf(atrib,"/bc_electroosmotic/bc%d/nminus/type %%s",h);
            ifd = fy_document_scanf(fyd,atrib,aux);
            if (strcmp(aux,"dirichlet") == 0) {
                nminusbctypes[h] = DIRICHLET;
            } else if (strcmp(aux,"neumann") == 0) {
                nminusbctypes[h] = NEUMANN;
            } else {
                printf("=+=+=+= Error loading boundary condition type for mminus in the boundary %d (may not be implemented yet) =+=+=+= \n",h);
                exit(1);
            }
            // Negative Charge Concentration boundary condition value type
            sprintf(atrib,"/bc_electroosmotic/bc%d/nminus/value_type %%s",h);
            ifd = fy_document_scanf(fyd,atrib,aux);
            if (strcmp(aux,"fixed_value") == 0) {
                nminusbcvaluetype[h] = fixedValue;
            } else if(strcmp(aux,"time_dependent") == 0) {
                nminusbcvaluetype[h] = timedependent;
            } else {
                printf("=+=+=+= Error loading boundary condition valuetype for nminus in the boundary %d (may not be implemented yet) =+=+=+= \n",h);
                exit(1);
            }
          
        }
        
        // Setting the boundary conditions for the electro-osmotic phi
        higflow_set_boundary_condition_for_electroosmotic_phi(ns, numbcs, id, amrBCfilename, phibctypes, phibcvaluetype);
        // Setting the boundary conditions for the electro-osmotic psi
        higflow_set_boundary_condition_for_electroosmotic_psi(ns, numbcs, id, amrBCfilename, psibctypes, psibcvaluetype);
        // Setting the boundary conditions for the electro-osmotic nplus
        higflow_set_boundary_condition_for_electroosmotic_nplus(ns, numbcs, id, amrBCfilename, nplusbctypes, nplusbcvaluetype);
        // Setting the boundary conditions for the electro-osmotic nminus
        higflow_set_boundary_condition_for_electroosmotic_nminus(ns, numbcs, id, amrBCfilename, nminusbctypes, nminusbcvaluetype);
        // Setting the boundary conditions for the electro-osmotic source term
        higflow_set_boundary_condition_for_electroosmotic_source_term(ns, numbcs, id, amrBCfilename, ubctypes, ubcvaluetype);

        eo_controllers *eo_contr;
        if (ns->contr.eoflow == true) eo_contr = &(ns->ed.eo.contr);
        else eo_contr = &(ns->ed.mult.eo.contr);
        
        eo_contr->is_phibc_timedependent = phibc_timedependent;
        if(eo_contr->is_phibc_timedependent == true){
            print0f("=+=+=+= Applied potential phi boundary conditions are time dependent =+=+=+=\n");
        } else {
            print0f("=+=+=+= Applied potential phi boundary conditions are time independent - Will only be calculated once in the Laplace equation=+=+=+=\n");
        }
        eo_contr->is_psibc_timedependent = psibc_timedependent;
        if(eo_contr->is_psibc_timedependent == true){
            print0f("=+=+=+= Electro-osmotic potential psi boundary conditions are time dependent =+=+=+=\n");
        } else {
            print0f("=+=+=+= Electro-osmotic potential psi boundary conditions are time independent ");
            if(eo_contr->eo_model == PB || eo_contr->eo_model == PBDH || eo_contr->eo_model == PBDH_ANALYTIC){
                print0f("- Will only be calculated once in the poisson equation");
            }
            print0f("=+=+=+=\n");
        }

    }

    fy_document_destroy(fyd);
}
