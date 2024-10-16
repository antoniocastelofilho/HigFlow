// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Boundary Condition - version 10/11/2016
// *******************************************************************
// *******************************************************************

#include "hig-flow-bc.h"


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
            real val = ns->func.get_boundary_pressure(id[h], bccenter, ns->par.t);
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

// Creating and setting the boundary condition for the density number nA (viscoelastic flows with shear-banding)
void higflow_set_boundary_condition_for_viscoelastic_shear_banding_nA(higflow_solver *ns, int numbcs, int id[numbcs], char bcfilenames[numbcs][1024], bc_type nAbctypes[], bc_valuetype nAbcvaluetype[]) {
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
        sim_domain *sd = psd_get_local_domain(ns->ed.vesb.psdSBnA);
        // Set the type bondary condition
        bc_type bc_t = ((nAbctypes[h]==0)?DIRICHLET:NEUMANN);
        // Create the bounary condition
        sim_boundary *bc = higflow_make_bc(bcg[h], nAbctypes[h], id[h], nAbcvaluetype[h]);
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
            if (nAbctypes[h] == NEUMANN) {
                // Get the density number defined by the user
                val = ns->ed.vesb.get_boundary_nA(id[h], bccenter, ns->par.t);
            } else {
                // Dirichlet boundary type
                val = ns->ed.vesb.get_boundary_nA(id[h], bccenter, ns->par.t);
            }
            // Set the value 
            sb_set_value(bc, bcgid, val);
        }
        // Destroy the iterator
        higcit_destroy(it);
    }
}

// Creating and setting the boundary condition for the density number nB (viscoelastic flows with shear-banding)
void higflow_set_boundary_condition_for_viscoelastic_shear_banding_nB(higflow_solver *ns, int numbcs, int id[numbcs], char bcfilenames[numbcs][1024], bc_type nBbctypes[], bc_valuetype nBbcvaluetype[]) {
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
        sim_domain *sd = psd_get_local_domain(ns->ed.vesb.psdSBnB);
        // Set the type bondary condition
        bc_type bc_t = ((nBbctypes[h]==0)?DIRICHLET:NEUMANN);
        // Create the bounary condition
        sim_boundary *bc = higflow_make_bc(bcg[h], nBbctypes[h], id[h], nBbcvaluetype[h]);
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
            if (nBbctypes[h] == NEUMANN) {
                // Get the electro-osmotic nplus defined by the user
                val = ns->ed.vesb.get_boundary_nB(id[h], bccenter, ns->par.t);
            } else {
                // Dirichlet boundary type
                val = ns->ed.vesb.get_boundary_nB(id[h], bccenter, ns->par.t);
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
    // Setting the boundary conditions for the VCM model (viscoelastic with shear-banding)
    if (ns->contr.rheotype == VCM) {
        sprintf(namefile,"%s.bcvesb",ns->par.nameload);
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
        bc_type       nAbctypes[numbcs]; 
        bc_valuetype  nAbcvaluetype[numbcs];
        bc_type       nBbctypes[numbcs]; 
        bc_valuetype  nBbcvaluetype[numbcs];
        for(int h = 0; h < numbcs; h++) {
            ifd = fscanf(fbc,"%d\n",&(id[h]));
            // HigTree Boundary condition file name
            __higflow_readstring(amrBCfilename[h],1024,fbc);
            // Pressure boundary condition type
            int aux;
            // Density number of spece A nA boundary condition type
            ifd = fscanf(fbc,"%d",&aux);
            switch (aux) {
                case 0:
                    nAbctypes[h] = DIRICHLET;
                    break;
                case 1:
                    nAbctypes[h] = NEUMANN;
                    break;
            }
            // Density number of spece A nA boundary condition value
            ifd = fscanf(fbc,"%d",&aux);
            switch (aux) {
                case 0:
                    nAbcvaluetype[h] = fixedValue;
                    break;
                case 1:
                    nAbcvaluetype[h] = zeroGradient;
                    break;
                case 2:
                    nAbcvaluetype[h] = freestream;
                    break;
                case 3:
                    nAbcvaluetype[h] = empty;
                    break;
                case 4:
                    nAbcvaluetype[h] = timedependent;
                    break;
                case 5:
                    nAbcvaluetype[h] = outflow;
                    break;
            }
            // Density number of spece B nB boundary condition type
            ifd = fscanf(fbc,"%d",&aux);
            switch (aux) {
                case 0:
                    nBbctypes[h] = DIRICHLET;
                    break;
                case 1:
                    nBbctypes[h] = NEUMANN;
                    break;
            }
            // Density number of spece B nB boundary condition value
            ifd = fscanf(fbc,"%d",&aux);
            switch (aux) {
                case 0:
                    nBbcvaluetype[h] = fixedValue;
                    break;
                case 1:
                    nBbcvaluetype[h] = zeroGradient;
                    break;
                case 2:
                    nBbcvaluetype[h] = freestream;
                    break;
                case 3:
                    nBbcvaluetype[h] = empty;
                    break;
                case 4:
                    nBbcvaluetype[h] = timedependent;
                    break;
                case 5:
                    nBbcvaluetype[h] = outflow;
                    break;
            }
        }
        fclose(fbc);
        // Setting the boundary conditions for the density number nB of specie B
        higflow_set_boundary_condition_for_viscoelastic_shear_banding_nA(ns, numbcs, id, amrBCfilename, nAbctypes, nAbcvaluetype);
        // Setting the boundary conditions for the density number nA of specie A
        higflow_set_boundary_condition_for_viscoelastic_shear_banding_nB(ns, numbcs, id, amrBCfilename, nBbctypes, nBbcvaluetype);
    }
    }
}
