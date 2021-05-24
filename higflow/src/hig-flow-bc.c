// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Boundary Condition - version 10/11/2016
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
        // Set the type bondary condition
        bc_type bc_t = ((pbctypes[h]==0)?DIRICHLET:NEUMANN);
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
            int bcgid = mp_lookup(bm, hig_get_cid(bcell));
            // Set the time to get the pressure
            real t = ns->par.t + ns->par.dt;
            // Get the pressure defined by the user
            real val;
            //if (pbctypes[h] == NEUMANN) {
            //     // Non incremental projection method
            //     val = ns->func.get_boundary_pressure(id[h], bccenter, ns->par.t);
            //} else {
            //     if (ns->contr.projtype == 0) {
            //          // Non incremental projection method
            //          val = ns->func.get_boundary_pressure(id[h], bccenter, ns->par.t);
            //     } else {
            //          // Incremental projection method
            //          val = ns->func.get_boundary_pressure(id[h], bccenter, t) - ns->func.get_boundary_pressure(id[h], bccenter, ns->par.t);
            //     }
            //}
            // Set the value 
            val = ns->func.get_boundary_pressure(id[h], bccenter, ns->par.t);
            sb_set_value(bc, bcgid, val);
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
            // Set the type bondary condition
            bc_type bc_t = ((bctypes[dim][h]==0)?DIRICHLET:NEUMANN);
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
                int  bcgid = mp_lookup(bm, hig_get_cid(bcell));
                //real bcval = bcvalues[dim][h];
                // Get the velocity defined by the user
                real bcval = ns->func.get_boundary_velocity(id[h], bccenter, dim, ns->par.t);
                // Set the value pbcvalues for the center of the cell 
                sb_set_value(bc, bcgid, bcval);
            }
            // Destroying the iterator 
            higcit_destroy(it);
        }
	
        // Mapping the properties in the domain (velocities)
        psfd_compute_sfbi(ns->psfdu[dim]);
        // Sync mapper for velocities
        psfd_synced_mapper(ns->psfdu[dim]);
    }
}

// Creating and setting the boundary condition for the electro-osmotic source term
void higflow_set_boundary_condition_for_electroosmotic_source_term(higflow_solver *ns, int numbcs, int id[numbcs], char bcfilenames[numbcs][1024], bc_type bctypes[DIM][numbcs], bc_valuetype bcvaluetype[DIM][numbcs]) {
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
            // Set the type bondary condition
            bc_type bc_t = DIRICHLET;
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
                int  bcgid = mp_lookup(bm, hig_get_cid(bcell));
                //real bcval = bcvalues[dim][h];
                // Get the velocity defined by the user
                real bcval = ns->ed.eo.get_boundary_electroosmotic_source_term(id[h], bccenter, dim, ns->par.t);
                // Set the value pbcvalues for the center of the cell 
                sb_set_value(bc, bcgid, bcval);
            }
            // Destroying the iterator 
            higcit_destroy(it);
        }
	
        // Mapping the properties in the domain (velocities)
        psfd_compute_sfbi(ns->ed.eo.psfdEOFeo[dim]);
        // Sync mapper for velocities
        psfd_synced_mapper(ns->ed.eo.psfdEOFeo[dim]);
    }
}

// Creating and setting the boundary condition for the electro-osmotic phi
void higflow_set_boundary_condition_for_electroosmotic_phi(higflow_solver *ns, int numbcs, int id[numbcs], char bcfilenames[numbcs][1024], bc_type phibctypes[], bc_valuetype phibcvaluetype[]) {
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
        // Set the type bondary condition
        bc_type bc_t = ((phibctypes[h]==0)?DIRICHLET:NEUMANN);
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
            int bcgid = mp_lookup(bm, hig_get_cid(bcell));
            // Set the time to get the pressure
            real t = ns->par.t + ns->par.dt;
            // Get the pressure defined by the user
            real val;
            if (phibctypes[h] == NEUMANN) {
                // Get the electro-osmotic phi defined by the user
                val = ns->ed.eo.get_boundary_electroosmotic_phi(id[h], bccenter, ns->par.t);
            } else {
                // Dirichlet boundary type
                val = ns->ed.eo.get_boundary_electroosmotic_phi(id[h], bccenter, ns->par.t);
            }
            // Set the value 
            sb_set_value(bc, bcgid, val);
        }
        // Destroy the iterator
        higcit_destroy(it);
    }
}

// Creating and setting the boundary condition for the electro-osmotic psi
void higflow_set_boundary_condition_for_electroosmotic_psi(higflow_solver *ns, int numbcs, int id[numbcs], char bcfilenames[numbcs][1024], bc_type psibctypes[], bc_valuetype psibcvaluetype[]) {
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
        // Set the type bondary condition
        bc_type bc_t = ((psibctypes[h]==0)?DIRICHLET:NEUMANN);
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
            int bcgid = mp_lookup(bm, hig_get_cid(bcell));
            // Set the time to get the pressure
            real t = ns->par.t + ns->par.dt;
            // Get the pressure defined by the user
            real val;
            if (psibctypes[h] == NEUMANN) {
                // Get the electro-osmotic psi defined by the user
                val = ns->ed.eo.get_boundary_electroosmotic_psi(id[h], bccenter, ns->par.t);
            } else {
                // Dirichlet boundary type
                val = ns->ed.eo.get_boundary_electroosmotic_psi(id[h], bccenter, ns->par.t);
            }
            // Set the value 
            sb_set_value(bc, bcgid, val);
        }
        // Destroy the iterator
        higcit_destroy(it);
    }
}

// Creating and setting the boundary condition for the electro-osmotic nplus
void higflow_set_boundary_condition_for_electroosmotic_nplus(higflow_solver *ns, int numbcs, int id[numbcs], char bcfilenames[numbcs][1024], bc_type nplusbctypes[], bc_valuetype nplusbcvaluetype[]) {
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
        // Set the type bondary condition
        bc_type bc_t = ((nplusbctypes[h]==0)?DIRICHLET:NEUMANN);
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
            int bcgid = mp_lookup(bm, hig_get_cid(bcell));
            // Set the time to get the pressure
            real t = ns->par.t + ns->par.dt;
            // Get the pressure defined by the user
            real val;
            if (nplusbctypes[h] == NEUMANN) {
                // Get the electro-osmotic nplus defined by the user
                val = ns->ed.eo.get_boundary_electroosmotic_nplus(id[h], bccenter, ns->par.t);
            } else {
                // Dirichlet boundary type
                val = ns->ed.eo.get_boundary_electroosmotic_nplus(id[h], bccenter, ns->par.t);
            }
            // Set the value 
            sb_set_value(bc, bcgid, val);
        }
        // Destroy the iterator
        higcit_destroy(it);
    }
}

// Creating and setting the boundary condition for the electro-osmotic nminus
void higflow_set_boundary_condition_for_electroosmotic_nminus(higflow_solver *ns, int numbcs, int id[numbcs], char bcfilenames[numbcs][1024], bc_type nminusbctypes[], bc_valuetype nminusbcvaluetype[]) {
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
        // Set the type bondary condition
        bc_type bc_t = ((nminusbctypes[h]==0)?DIRICHLET:NEUMANN);
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
            int bcgid = mp_lookup(bm, hig_get_cid(bcell));
            // Set the time to get the pressure
            real t = ns->par.t + ns->par.dt;
            // Get the pressure defined by the user
            real val;
            if (nminusbctypes[h] == NEUMANN) {
                // Get the electro-osmotic nplus defined by the user
                val = ns->ed.eo.get_boundary_electroosmotic_nminus(id[h], bccenter, ns->par.t);
            } else {
                // Dirichlet boundary type
                val = ns->ed.eo.get_boundary_electroosmotic_nminus(id[h], bccenter, ns->par.t);
            }
            // Set the value 
            sb_set_value(bc, bcgid, val);
        }
        // Destroy the iterator
        higcit_destroy(it);
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
        // Set the type bondary condition
        pbctypes[h] = DIRICHLET;
        pbcvaluetype[h] = timedependent;
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
            int bcgid = mp_lookup(bm, hig_get_cid(bcell));
            // Get the electro-osmotic nminus defined by the user
            real val = ns->func.get_boundary_source_term(id[h], bccenter, ns->par.t);
            // Set the value 
            sb_set_value(bc, bcgid, val);
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
            // Set the type bondary condition
            //bc_type bc_t = ((bctypes[dim][h]==0)?DIRICHLET:NEUMANN);
            // Create the bounary condition
            //sim_boundary *bc = higflow_make_bc(bcg[h], bctypes[dim][h], id[h], bcvaluetype[dim][h]);
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
                int  bcgid = mp_lookup(bm, hig_get_cid(bcell));
                //real bcval = bcvalues[dim][h];
                // Get the velocity defined by the user
                real bcval = ns->func.get_boundary_facet_source_term(id[h], bccenter, dim, ns->par.t);
                // Set the value pbcvalues for the center of the cell 
                sb_set_value(bc, bcgid, bcval);
            }
            // Destroying the iterator 
            higcit_destroy(it);
        }
	
        // Mapping the properties in the domain (velocities)
        psfd_compute_sfbi(ns->psfdF[dim]);
        // Sync mapper for velocities
        psfd_synced_mapper(ns->psfdF[dim]);
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
    ns->contr.desingpressure = 1;
    for(int h = 0; h < numbcs; h++) {
        ifd = fscanf(fbc,"%d\n",&(id[h]));
        // HigTree Boundary condition file name
        __higflow_readstring(amrBCfilename[h],1024,fbc);
        // Pressure boundary condition type
        int aux;
        ifd = fscanf(fbc,"%d",&aux);
        switch (aux) {
            case 0:
                pbctypes[h] = DIRICHLET;
                break;
            case 1:
                pbctypes[h] = NEUMANN;
                break;
        }
        // Pressure boundary condition value
        ifd = fscanf(fbc,"%d",&aux);
        switch (aux) {
            case 0:
                pbcvaluetype[h] = fixedValue;
                break;
            case 1:
                pbcvaluetype[h] = zeroGradient;
                break;
            case 2:
                pbcvaluetype[h] = freestream;
                break;
            case 3:
                pbcvaluetype[h] = empty;
                break;
            case 4:
                pbcvaluetype[h] = timedependent;
                break;
            case 5:
                pbcvaluetype[h] = outflow;
                break;
        }
        // Setting the pressure desingularizadtion control
        if (pbctypes[h] == DIRICHLET) {
            // Outflow
            ns->contr.desingpressure = 0;
        }
        for (int dim = 0; dim < DIM; dim++) {
            // Velocity boundary condition type
            ifd = fscanf(fbc,"%d",&aux);
            switch (aux) {
                case 0:
                    ubctypes[dim][h] = DIRICHLET;
                    break;
                case 1:
                    ubctypes[dim][h] = NEUMANN;
                    break;
            }
            // Velocity boundary condition valuetype
            ifd = fscanf(fbc,"%d",&aux);
            switch (aux) {
                case 0:
                    ubcvaluetype[dim][h] = fixedValue;
                    break;
                case 1:
                    ubcvaluetype[dim][h] = zeroGradient;
                    break;
                case 2:
                    ubcvaluetype[dim][h] = freestream;
                    break;
                case 3:
                    ubcvaluetype[dim][h] = empty;
                    break;
                case 4:
                    ubcvaluetype[dim][h] = timedependent;
                    break;
                case 5:
                    ubcvaluetype[dim][h] = outflow;
                    break;
            }
        }
    }
    fclose(fbc);
    // Setting the boundary conditions for the pressure
    higflow_set_boundary_condition_for_pressure(ns, numbcs, id, amrBCfilename, pbctypes, pbcvaluetype);
    // Setting the boundary conditions for the velocities
    higflow_set_boundary_condition_for_velocities(ns, numbcs, id, amrBCfilename, ubctypes, ubcvaluetype);
    // Setting the boundary conditions for the electro-osmotic model
    if (ns->contr.modelflowtype == 1) {
        sprintf(namefile,"%s.bcelectroosmotic",ns->par.nameload);
        printf("nome: %s", namefile);
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
            // Pressure boundary condition type
            int aux;
            ifd = fscanf(fbc,"%d",&aux);
            // Electro-osmotic potential boundary condition type
            switch (aux) {
                case 0:
                    phibctypes[h] = DIRICHLET;
                    break;
                case 1:
                    phibctypes[h] = NEUMANN;
                    break;
            }
            // Electro-osmotic potential boundary condition value
            ifd = fscanf(fbc,"%d",&aux);
            switch (aux) {
                case 0:
                    phibcvaluetype[h] = fixedValue;
                    break;
                case 1:
                    phibcvaluetype[h] = zeroGradient;
                    break;
                case 2:
                    phibcvaluetype[h] = freestream;
                    break;
                case 3:
                    phibcvaluetype[h] = empty;
                    break;
                case 4:
                    phibcvaluetype[h] = timedependent;
                    break;
                case 5:
                    phibcvaluetype[h] = outflow;
                    break;
            }
            // Electro-osmotic potential boundary condition type
            ifd = fscanf(fbc,"%d",&aux);
            switch (aux) {
                case 0:
                    psibctypes[h] = DIRICHLET;
                    break;
                case 1:
                    psibctypes[h] = NEUMANN;
                    break;
            }
            // Electro-osmotic potential boundary condition value
            ifd = fscanf(fbc,"%d",&aux);
            switch (aux) {
                case 0:
                    psibcvaluetype[h] = fixedValue;
                    break;
                case 1:
                    psibcvaluetype[h] = zeroGradient;
                    break;
                case 2:
                    psibcvaluetype[h] = freestream;
                    break;
                case 3:
                    psibcvaluetype[h] = empty;
                    break;
                case 4:
                    psibcvaluetype[h] = timedependent;
                    break;
                case 5:
                    psibcvaluetype[h] = outflow;
                    break;
            }
            // Electro-osmotic concentration n+ boundary condition type
            ifd = fscanf(fbc,"%d",&aux);
            switch (aux) {
                case 0:
                    nplusbctypes[h] = DIRICHLET;
                    break;
                case 1:
                    nplusbctypes[h] = NEUMANN;
                    break;
            }
            // Electro-osmotic concentration n+ boundary condition value
            ifd = fscanf(fbc,"%d",&aux);
            switch (aux) {
                case 0:
                    nplusbcvaluetype[h] = fixedValue;
                    break;
                case 1:
                    nplusbcvaluetype[h] = zeroGradient;
                    break;
                case 2:
                    nplusbcvaluetype[h] = freestream;
                    break;
                case 3:
                    nplusbcvaluetype[h] = empty;
                    break;
                case 4:
                    nplusbcvaluetype[h] = timedependent;
                    break;
                case 5:
                    nplusbcvaluetype[h] = outflow;
                    break;
            }
            // Electro-osmotic concentration n- boundary condition type
            ifd = fscanf(fbc,"%d",&aux);
            switch (aux) {
                case 0:
                    nminusbctypes[h] = DIRICHLET;
                    break;
                case 1:
                    nminusbctypes[h] = NEUMANN;
                    break;
            }
            // Electro-osmotic concentration n- boundary condition value
            ifd = fscanf(fbc,"%d",&aux);
            switch (aux) {
                case 0:
                    nminusbcvaluetype[h] = fixedValue;
                    break;
                case 1:
                    nminusbcvaluetype[h] = zeroGradient;
                    break;
                case 2:
                    nminusbcvaluetype[h] = freestream;
                    break;
                case 3:
                    nminusbcvaluetype[h] = empty;
                    break;
                case 4:
                    nminusbcvaluetype[h] = timedependent;
                    break;
                case 5:
                    nminusbcvaluetype[h] = outflow;
                    break;
            }
        }
        fclose(fbc);
        // Setting the boundary conditions for the electro-osmotic source term
        higflow_set_boundary_condition_for_electroosmotic_source_term(ns, numbcs, id, amrBCfilename, ubctypes, ubcvaluetype);
        // Setting the boundary conditions for the electro-osmotic phi
        higflow_set_boundary_condition_for_electroosmotic_phi(ns, numbcs, id, amrBCfilename, phibctypes, phibcvaluetype);
        // Setting the boundary conditions for the electro-osmotic psi
        higflow_set_boundary_condition_for_electroosmotic_psi(ns, numbcs, id, amrBCfilename, psibctypes, psibcvaluetype);
        // Setting the boundary conditions for the electro-osmotic nplus
        higflow_set_boundary_condition_for_electroosmotic_nplus(ns, numbcs, id, amrBCfilename, nplusbctypes, nplusbcvaluetype);
        // Setting the boundary conditions for the electro-osmotic nminus
        higflow_set_boundary_condition_for_electroosmotic_nminus(ns, numbcs, id, amrBCfilename, nminusbctypes, nminusbcvaluetype);
        // Setting the boundary conditions for the cell source term 
        higflow_set_boundary_condition_for_cell_source_term(ns, numbcs, id, amrBCfilename, pbctypes, pbcvaluetype);
        higflow_set_boundary_condition_for_facet_source_term(ns, numbcs, id, amrBCfilename, ubctypes, ubcvaluetype);
    }
}

// Navier-Stokes initialize the domain and boudaries
void higflow_initialize_boundaries_and_domain(higflow_solver *ns) {
    
    // Loading the boundary condition data
    char namefile[1024];
    sprintf(namefile,"%s.bc.domain.yaml",ns->par.nameload);
    
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
    int ifd = fy_document_scanf(fyd,"/boundary_condition/number_bc %d",&numbcs);
    
    // Boudary condition data
    int           id[numbcs];
    char          amrBCfilename[numbcs][1024]; 
    bc_type       pbctypes[numbcs]; 
    bc_type       ubctypes[DIM][numbcs]; 
    bc_valuetype  pbcvaluetype[numbcs];
    bc_valuetype  ubcvaluetype[DIM][numbcs]; 
    // Setting the pressure desingularizadtion control
    ns->contr.desingpressure = 1;
    
    for(int h = 0; h < numbcs; h++) {
        char atrib[1024];
        sprintf(atrib,"/boundary_condition/bc%d/id %%d",h);
        ifd = fy_document_scanf(fyd,atrib,&(id[h]));
        // HigTree Boundary condition file name
        sprintf(atrib,"/boundary_condition/bc%d/path %%s",h);
        ifd = fy_document_scanf(fyd,atrib,amrBCfilename[h]);
        //__higflow_readstring(amrBCfilename[h],1024,fbc);
        // Pressure boundary condition type
        char aux[1024];
        sprintf(atrib,"/boundary_condition/bc%d/pressure/type %%s",h);
        ifd = fy_document_scanf(fyd,atrib,aux);
        if (strcmp(aux,"dirichelet") == 0) {
           pbctypes[h] = DIRICHLET;
        } else if (strcmp(aux,"neumann") == 0) {
                pbctypes[h] = NEUMANN;
        } else {
            printf("=+=+=+= Error loading boundary condition type for the pressure in the boundary %d =+=+=+= \n",h);
            exit(1);
        }
        // Pressure boundary condition value
        sprintf(atrib,"/boundary_condition/bc%d/pressure/sub_type %%s",h);
        ifd = fy_document_scanf(fyd,atrib,aux);
        if (strcmp(aux,"fixed_value") == 0) {
            pbcvaluetype[h] = fixedValue;
         } else if(strcmp(aux,"zero_gradient") == 0) {
            pbcvaluetype[h] = zeroGradient;
         } else if(strcmp(aux,"freestream") == 0) {
            pbcvaluetype[h] = freestream;
         } else if(strcmp(aux,"empty") == 0) {
            pbcvaluetype[h] = empty;
         } else if(strcmp(aux,"time_independent") == 0) {
            pbcvaluetype[h] = timedependent;
         } else if(strcmp(aux,"outflow") == 0) {
             pbcvaluetype[h] = outflow;
         } else {
            printf("=+=+=+= Error loading boundary condition valuetype for the pressure in the boundary %d =+=+=+= \n",h);
            exit(1);
        }
        // Setting the pressure desingularizadtion control
        if (pbctypes[h] == DIRICHLET) {
            // Outflow
            //sprintf(atrib,"/boundary_condition/bc%d/pressure/sub_type/value %%lf",h);
            //ifd = fy_document_scanf(fyd,atrib,&(ns->contr.desingpressure));
            ns->contr.desingpressure = 0;
        }
        // Velocity boundary condition value
        for (int dim = 0; dim < DIM; dim++) {
            sprintf(atrib,"/boundary_condition/bc%d/velocity_%d/type %%s",h,dim);
            ifd = fy_document_scanf(fyd,atrib,aux);
            if (strcmp(aux,"dirichelet") == 0) {
               ubctypes[dim][h] = DIRICHLET;
            } else if(strcmp(aux,"neumann") == 0) {
               ubctypes[dim][h] = NEUMANN;
            } else {
               printf("=+=+=+= Error loading boundary condition type for the velocity in the boundary %d =+=+=+= \n",h);
               exit(1);
            }
            sprintf(atrib,"/boundary_condition/bc%d/velocity_%d/sub_type %%s",h,dim);
            ifd = fy_document_scanf(fyd,atrib,aux);
            // Velocity boundary condition valuetype
            if (strcmp(aux,"fixed_value") == 0) {
               ubcvaluetype[dim][h] = fixedValue;
            } else if(strcmp(aux,"zero_gradient") == 0) {
               ubcvaluetype[dim][h] = zeroGradient;
            } else if(strcmp(aux,"freestream") == 0) {
               ubcvaluetype[dim][h] = freestream;
            } else if(strcmp(aux,"empty") == 0) {
               ubcvaluetype[dim][h] = empty;
            } else if(strcmp(aux,"time_independent") == 0) {
               ubcvaluetype[dim][h] = timedependent;
            } else if(strcmp(aux,"outflow") == 0) {
               ubcvaluetype[dim][h] = outflow;
            } else {
               printf("=+=+=+= Error loading boundary condition valuetype for the velocity in the boundary %d =+=+=+= \n",h);
               exit(1);
            }
        }
    }
    
    fy_document_destroy(fyd);
    // Setting the boundary conditions for the pressure
    higflow_set_boundary_condition_for_pressure(ns, numbcs, id, amrBCfilename, pbctypes, pbcvaluetype);
    // Setting the boundary conditions for the velocities
    higflow_set_boundary_condition_for_velocities(ns, numbcs, id, amrBCfilename, ubctypes, ubcvaluetype);
    // Setting the boundary conditions for the electro-osmotic model
    
    if (ns->contr.modelflowtype == 1) {
         printf("=+=+=+= electroosmotic - NOT IMPLEMENTED READ YAML =+=+=+= \n");
        sprintf(namefile,"%s.bcelectroosmotic",ns->par.nameload);
        printf("nome: %s", namefile);
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
            // Pressure boundary condition type
            int aux;
            ifd = fscanf(fbc,"%d",&aux);
            // Electro-osmotic potential boundary condition type
            switch (aux) {
                case 0:
                    phibctypes[h] = DIRICHLET;
                    break;
                case 1:
                    phibctypes[h] = NEUMANN;
                    break;
            }
            // Electro-osmotic potential boundary condition value
            ifd = fscanf(fbc,"%d",&aux);
            switch (aux) {
                case 0:
                    phibcvaluetype[h] = fixedValue;
                    break;
                case 1:
                    phibcvaluetype[h] = zeroGradient;
                    break;
                case 2:
                    phibcvaluetype[h] = freestream;
                    break;
                case 3:
                    phibcvaluetype[h] = empty;
                    break;
                case 4:
                    phibcvaluetype[h] = timedependent;
                    break;
                case 5:
                    phibcvaluetype[h] = outflow;
                    break;
            }
            // Electro-osmotic potential boundary condition type
            ifd = fscanf(fbc,"%d",&aux);
            switch (aux) {
                case 0:
                    psibctypes[h] = DIRICHLET;
                    break;
                case 1:
                    psibctypes[h] = NEUMANN;
                    break;
            }
            // Electro-osmotic potential boundary condition value
            ifd = fscanf(fbc,"%d",&aux);
            switch (aux) {
                case 0:
                    psibcvaluetype[h] = fixedValue;
                    break;
                case 1:
                    psibcvaluetype[h] = zeroGradient;
                    break;
                case 2:
                    psibcvaluetype[h] = freestream;
                    break;
                case 3:
                    psibcvaluetype[h] = empty;
                    break;
                case 4:
                    psibcvaluetype[h] = timedependent;
                    break;
                case 5:
                    psibcvaluetype[h] = outflow;
                    break;
            }
            // Electro-osmotic concentration n+ boundary condition type
            ifd = fscanf(fbc,"%d",&aux);
            switch (aux) {
                case 0:
                    nplusbctypes[h] = DIRICHLET;
                    break;
                case 1:
                    nplusbctypes[h] = NEUMANN;
                    break;
            }
            // Electro-osmotic concentration n+ boundary condition value
            ifd = fscanf(fbc,"%d",&aux);
            switch (aux) {
                case 0:
                    nplusbcvaluetype[h] = fixedValue;
                    break;
                case 1:
                    nplusbcvaluetype[h] = zeroGradient;
                    break;
                case 2:
                    nplusbcvaluetype[h] = freestream;
                    break;
                case 3:
                    nplusbcvaluetype[h] = empty;
                    break;
                case 4:
                    nplusbcvaluetype[h] = timedependent;
                    break;
                case 5:
                    nplusbcvaluetype[h] = outflow;
                    break;
            }
            // Electro-osmotic concentration n- boundary condition type
            ifd = fscanf(fbc,"%d",&aux);
            switch (aux) {
                case 0:
                    nminusbctypes[h] = DIRICHLET;
                    break;
                case 1:
                    nminusbctypes[h] = NEUMANN;
                    break;
            }
            // Electro-osmotic concentration n- boundary condition value
            ifd = fscanf(fbc,"%d",&aux);
            switch (aux) {
                case 0:
                    nminusbcvaluetype[h] = fixedValue;
                    break;
                case 1:
                    nminusbcvaluetype[h] = zeroGradient;
                    break;
                case 2:
                    nminusbcvaluetype[h] = freestream;
                    break;
                case 3:
                    nminusbcvaluetype[h] = empty;
                    break;
                case 4:
                    nminusbcvaluetype[h] = timedependent;
                    break;
                case 5:
                    nminusbcvaluetype[h] = outflow;
                    break;
            }
        }
        fclose(fbc);
        // Setting the boundary conditions for the electro-osmotic source term
        higflow_set_boundary_condition_for_electroosmotic_source_term(ns, numbcs, id, amrBCfilename, ubctypes, ubcvaluetype);
        // Setting the boundary conditions for the electro-osmotic phi
        higflow_set_boundary_condition_for_electroosmotic_phi(ns, numbcs, id, amrBCfilename, phibctypes, phibcvaluetype);
        // Setting the boundary conditions for the electro-osmotic psi
        higflow_set_boundary_condition_for_electroosmotic_psi(ns, numbcs, id, amrBCfilename, psibctypes, psibcvaluetype);
        // Setting the boundary conditions for the electro-osmotic nplus
        higflow_set_boundary_condition_for_electroosmotic_nplus(ns, numbcs, id, amrBCfilename, nplusbctypes, nplusbcvaluetype);
        // Setting the boundary conditions for the electro-osmotic nminus
        higflow_set_boundary_condition_for_electroosmotic_nminus(ns, numbcs, id, amrBCfilename, nminusbctypes, nminusbcvaluetype);
        // Setting the boundary conditions for the cell source term 
        higflow_set_boundary_condition_for_cell_source_term(ns, numbcs, id, amrBCfilename, pbctypes, pbcvaluetype);
        higflow_set_boundary_condition_for_facet_source_term(ns, numbcs, id, amrBCfilename, ubctypes, ubcvaluetype);
    }
}

