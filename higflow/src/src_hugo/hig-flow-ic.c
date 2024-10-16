// *******************************************************************
// *******************************************************************
#include "hig-flow-ic.h"

// *******************************************************************
// Navier-Stokes Initialize Properties
// *******************************************************************
// Navier-Stokes initialize the domain
void higflow_initialize_domain(higflow_solver *ns, int ntasks, int myrank, int order) {
    // Loading the domain data
    char namefile[1024];
    sprintf(namefile,"%s.domain",ns->par.nameload);
    FILE *fdomain = fopen(namefile, "r");
    if (fdomain == NULL) {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n",namefile);
        exit(1);
    }
    // Number of HigTrees
    int numhigs;
    int ifd = fscanf(fdomain,"%d\n",&numhigs);
    higio_amr_info *mi[numhigs];
    for(int h = 0; h < numhigs; h++) {
        // Name of the HigTree file
        //char *amrfilename = argv[1]; argv++;
        char amrfilename[1024];
        __higflow_readstring(amrfilename,1024,fdomain);
        // Open the AMR format file
        FILE *fd = fopen(amrfilename, "r");
        // Reading the higtree information from the file 
        mi[h] = higio_read_amr_info(fd);
        // Close the AMR format file
        fclose(fd);
    }
    fclose(fdomain);
    // Taking the sub-domain of the myrank process
    // Creates a partition table.
    partition_graph *pg = pg_create(MPI_COMM_WORLD);
    // Initializing partition table
    higflow_partition_domain(ns, pg, numhigs, mi, ntasks, myrank);
    // Creating the partitioned sub-domain to simulation
    higflow_create_partitioned_domain(ns, pg, order);
    // Creating the stencil for properties interpolation
    higflow_create_stencil(ns);
    // Creating the partitioned sub-domain for extra properties
    if (ns->contr.flowtype == 1) {
        // Initialize generalized newtonian domain
        higflow_create_partitioned_domain_generalized_newtonian(ns, pg, order);
        // Creating the stencil for properties interpolation
        higflow_create_stencil_for_extra_domain(ns);
    } else if (ns->contr.flowtype == 2) {
        // Initialize multiphase domain
        higflow_create_partitioned_domain_multiphase(ns, pg, order);
        // Creating the stencil for properties interpolation
        higflow_create_stencil_for_extra_domain(ns);
    } else if (ns->contr.flowtype == 3) {
        // Initialize visoelastic tensor domain
        higflow_create_partitioned_domain_viscoelastic(ns, pg, order);
        // Creating the stencil for properties interpolation
        higflow_create_stencil_for_extra_domain(ns);
    } else if (ns->contr.flowtype == 4) {
        // Initialize integral tensors domain
        higflow_create_partitioned_domain_viscoelastic_integral(ns, pg, order);
        // Creating the stencil for properties interpolation
        higflow_create_stencil_for_extra_domain(ns);
    } else if (ns->contr.flowtype == 6) {
        // Initialize viscoelastic flow with variable viscosity domain
        higflow_create_partitioned_domain_viscoelastic_variable_viscosity(ns, pg, order);
        // Creating the stencil for properties interpolation
        higflow_create_stencil_for_extra_domain(ns);
    } else if (ns->contr.flowtype == 7) {
        // Initialize viscoelastic flow with shear-banding domain
        higflow_create_partitioned_domain_viscoelastic_shear_banding(ns, pg, order);
        // Creating the stencil for properties interpolation
        higflow_create_stencil_for_extra_domain(ns);
    } else if (ns->contr.flowtype == 8) {
        // Initialize elastoviscoplastic tensor domain
        higflow_create_partitioned_domain_elastoviscoplastic(ns, pg, order);
        // Creating the stencil for properties interpolation
        higflow_create_stencil_for_extra_domain(ns);
    } else if (ns->contr.flowtype == 9) {
        // Initialize shear-thickening supension tensor domain
        higflow_create_partitioned_domain_shear_thickening_suspension(ns, pg, order);
        // Creating the stencil for properties interpolation
        higflow_create_stencil_for_extra_domain(ns);
    }
    if (ns->contr.modelflowtype == 1) {
        // Initialize electro-osmotic domain
        higflow_create_partitioned_domain_electroosmotic(ns, pg, order);
        // Creating the stencil for properties interpolation
        higflow_create_stencil_for_extra_domain(ns);
    }
    if (ns->contr.modelflowtype == 5) {
        // Initialize non-isothermal flows domain
        higflow_create_partitioned_domain_non_isothermal_flows(ns, pg, order);
        // Creating the stencil for properties interpolation
        higflow_create_stencil_for_extra_domain(ns);
    }
}

// Initialize the pressure
void higflow_initialize_pressure(higflow_solver *ns) {
    // Setting the cell iterator
    higcit_celliterator *it;
    // Getting the local domain
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    // Getting the Mapper for the local domain
    mp_mapper *m = sd_get_domain_mapper(sdp);
    // Traversing the cells of local domain
    for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Getting the cell
        hig_cell *c = higcit_getcell(it);
        // Get the cell identifier
        int cgid = mp_lookup(m, hig_get_cid(c));
        // Get the center of the cell
        Point center;
        hig_get_center(c, center);
        // Get the value for pressure in this cell
        real val = ns->func.get_pressure(center, ns->par.t);
        // Set the value for pressure distributed property
        dp_set_value(ns->dpp, cgid, val);
    }
    // Destroying the iterator
    higcit_destroy(it);
    // Sync initial values among processes
    dp_sync(ns->dpp);
}

// Initialize the viscosity
void higflow_initialize_viscosity_gn(higflow_solver *ns) {
    // Setting the cell iterator
    higcit_celliterator *it;
    // Getting the local domain 
    sim_domain *sdvisc = psd_get_local_domain(ns->ed.psdED);
    // Getting the Mapper for the local domain 
    mp_mapper *m = sd_get_domain_mapper(sdvisc);
    // Traversing the cells of local domain
    for(it = sd_get_domain_celliterator(sdvisc); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Getting the cell
        hig_cell *c = higcit_getcell(it);
        // Get the cell identifier
        int cgid = mp_lookup(m, hig_get_cid(c));
        // Get the center of the cell
        Point center;
        hig_get_center(c, center);
        // Get the value for pressure in this cell
        real val = ns->ed.gn.get_viscosity(center, 0.0, ns->par.t);
        // Set the value for pressure distributed property
        dp_set_value(ns->ed.gn.dpvisc, cgid, val);
    }
    // Destroying the iterator
    higcit_destroy(it);
    // Sync initial values among processes
    dp_sync(ns->ed.gn.dpvisc);
}

// Initialize the viscosity
void higflow_initialize_viscosity_mult(higflow_solver *ns) {
    // Setting the cell iterator
    higcit_celliterator *it;
    // Getting the local domain 
    sim_domain *sdvisc = psd_get_local_domain(ns->ed.psdED);
    // Getting the Mapper for the local domain 
    mp_mapper *m = sd_get_domain_mapper(sdvisc);
    // Traversing the cells of local domain
    for(it = sd_get_domain_celliterator(sdvisc); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Getting the cell
        hig_cell *c = higcit_getcell(it);
        // Get the cell identifier
        int cgid = mp_lookup(m, hig_get_cid(c));
        // Get the center of the cell
        Point center;
        hig_get_center(c, center);
        real fracvol = compute_value_at_point(sdvisc, center, center, 1.0, ns->ed.mult.dpfracvol, ns->stn);
        // Calculate the viscosity
        real visc0 = ns->ed.mult.get_viscosity0(center, ns->par.t);
        real visc1 = ns->ed.mult.get_viscosity1(center, ns->par.t);
	real visc  = (1.0-fracvol)*visc0 + fracvol*visc1;
        // Set the viscosity in the distributed viscosity property
        dp_set_value(ns->ed.mult.dpvisc, cgid, visc);
    }
    // Destroying the iterator
    higcit_destroy(it);
    // Sync initial values among processes
    dp_sync(ns->ed.mult.dpvisc);
}

// Initialize the volume fraction 
void higflow_initialize_fracvol(higflow_solver *ns) {
    // Setting the cell iterator
    higcit_celliterator *it;
    // Getting the local domain 
    sim_domain *sdfracvol = psd_get_local_domain(ns->ed.psdED);
    // Getting the Mapper for the local domain 
    mp_mapper *m = sd_get_domain_mapper(sdfracvol);
    // Traversing the cells of local domain
    for(it = sd_get_domain_celliterator(sdfracvol); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Getting the cell
        hig_cell *c = higcit_getcell(it);
        // Get the cell identifier
        int cgid = mp_lookup(m, hig_get_cid(c));
        // Get the center of the cell
        Point center, delta;
        hig_get_center(c, center);
        hig_get_delta(c, delta);
        // Get the value for pressure in this cell
        real val = ns->ed.mult.get_fracvol(center, delta, ns->par.t);
        // Set the value for pressure distributed property
        dp_set_value(ns->ed.mult.dpfracvol, cgid, val);
    }
    // Destroying the iterator
    higcit_destroy(it);
    // Sync initial values among processes
    dp_sync(ns->ed.mult.dpfracvol);
}

// Initialize the density
void higflow_initialize_density(higflow_solver *ns) {
    // Setting the cell iterator
    higcit_celliterator *it;
    // Getting the local domain
    sim_domain *sddens = psd_get_local_domain(ns->ed.psdED);
    // Getting the Mapper for the local domain
    mp_mapper *m = sd_get_domain_mapper(sddens);
    // Traversing the cells of local domain
    for(it = sd_get_domain_celliterator(sddens); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Getting the cell
        hig_cell *c = higcit_getcell(it);
        // Get the cell identifier
        int cgid = mp_lookup(m, hig_get_cid(c));
        // Get the center of the cell
        Point center;
        hig_get_center(c, center);
        // Calculate the density
        real fracvol  = compute_value_at_point(sddens, center, center, 1.0, ns->ed.mult.dpfracvol, ns->stn);
        // Calculate the density
        real dens0 = ns->ed.mult.get_density0(center, ns->par.t);
        real dens1 = ns->ed.mult.get_density1(center, ns->par.t);
	real dens  = (1.0-fracvol)*dens0 + fracvol*dens1;
        // Set the viscosity in the distributed viscosity property
        dp_set_value(ns->ed.mult.dpdens, cgid, dens);
    }
    // Destroying the iterator
    higcit_destroy(it);
    // Sync initial values among processes
    dp_sync(ns->ed.mult.dpdens);
}

// Initialize the cell source term
void higflow_initialize_cell_source_term(higflow_solver *ns) {
    // Setting the cell iterator
    higcit_celliterator *it;
    // Getting the local domain
    sim_domain *sdp = psd_get_local_domain(ns->psdF);
    // Getting the Mapper for the local domain
    mp_mapper *m = sd_get_domain_mapper(sdp);
    // Traversing the cells of local domain
    for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Getting the cell
        hig_cell *c = higcit_getcell(it);
        // Get the cell identifier
        int cgid = mp_lookup(m, hig_get_cid(c));
        // Get the center of the cell
        Point center;
        hig_get_center(c, center);
        // Get the value for the source term in this cell
        real val = ns->func.get_source_term(center, ns->par.t);
        // Set the value for source term distributed property
        dp_set_value(ns->dpF, cgid, val);
    }
    // Destroying the iterator
    higcit_destroy(it);
    // Sync initial values among processes
    dp_sync(ns->dpF);
}

// Initialize the Non-Newtonian Tensor
void higflow_initialize_viscoelastic_tensor(higflow_solver *ns) {
    // Non Newtonian flow
    if (ns->contr.flowtype == 3) {
        // Setting the cell iterator
        higcit_celliterator *it;
        // Getting the local domain
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Getting the Mapper for the local domain
        mp_mapper *m = sd_get_domain_mapper(sdp);
        // Traversing the cells of local domain
        for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Getting the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int cgid = mp_lookup(m, hig_get_cid(c));
            // Get the center of the cell
            Point center;
            hig_get_center(c, center);
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get the value for the tensor in this cell
                    real val = ns->ed.ve.get_tensor(center, i, j, ns->par.t);
                    // Set the value for tensor distributed property
                    dp_set_value(ns->ed.ve.dpS[i][j], cgid, val);
                }
            }
        }
        // Destroying the iterator
        higcit_destroy(it);
        // Sync initial values among processes
	for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.ve.dpS[i][j]);
	    }
	}
    }
}


// Initialize the viscoelastic integral Tensor
void higflow_initialize_viscoelastic_integral_tensor(higflow_solver *ns) {
    // Non Newtonian flow
    if (ns->contr.flowtype == 4) {
        // Setting the cell iterator
        higcit_celliterator *it;
        // Getting the local domain
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Getting the Mapper for the local domain
        mp_mapper *m = sd_get_domain_mapper(sdp);
        // Traversing the cells of local domain
        for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Getting the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int cgid = mp_lookup(m, hig_get_cid(c));
            // Get the center of the cell
            Point center;
            hig_get_center(c, center);
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get the value for the tensor in this cell
                    real val = ns->ed.im.get_tensor(center, i, j, ns->par.t);
                    // Set the value for tensor distributed property
                    dp_set_value(ns->ed.im.dpS[i][j], cgid, val);
                }
            }
        }
        // Destroying the iterator
        higcit_destroy(it);
        // Sync initial values among processes
	for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.im.dpS[i][j]);
	    }
	}
    }
}

// Initialize the viscoelastic integral finger Tensor
void higflow_initialize_viscoelastic_integral_finger_tensor(higflow_solver *ns) {
    // Non Newtonian flow
    if (ns->contr.flowtype == 4) {
        // Setting the cell iterator
        higcit_celliterator *it;
        // Getting the local domain
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Getting the Mapper for the local domain
        mp_mapper *m = sd_get_domain_mapper(sdp);
        // Traversing the cells of local domain
        for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Getting the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int cgid = mp_lookup(m, hig_get_cid(c));
            // Get the center of the cell
            Point center;
            hig_get_center(c, center);
            for (int k = 0; k <= NDT; k++) {
                for (int i = 0; i < DIM; i++) {
                    for (int j = 0; j < DIM; j++) {
                        // Define the value for the tensor in this cell
                        real val = 0.0;
                        if (i == j) val = 1.0;
                        // Set the value for tensor distributed property
                        dp_set_value(ns->ed.im.dpB[k][i][j], cgid, val);
                    }
                }
            }
        }
        // Destroying the iterator
        higcit_destroy(it);
        // Sync initial values among processes
        for (int k = 0; k <= NDT; k++) {
	    for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    dp_sync(ns->ed.im.dpB[k][i][j]);
	        }
	    }
	}
    }
}

// Initialize the electro-osmotic source term 
void higflow_initialize_electroosmotic_source_term(higflow_solver *ns) {
    // Setting the facet-cell iterator
    higfit_facetiterator *fit;
    // Setting the volocity values for the domain
    sim_facet_domain *sfdFeo[DIM];
    // Setting the velocity U(dim)
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local domain property
        sfdFeo[dim] = psfd_get_local_domain(ns->ed.eo.psfdEOFeo[dim]);
        // Get the Mapper for the local domain
        mp_mapper *m = sfd_get_domain_mapper(sfdFeo[dim]);
        for(fit = sfd_get_domain_facetiterator(sfdFeo[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Getting the cell
            hig_facet *f = higfit_getfacet(fit);
            // Get the cell identifier
            int flid = mp_lookup(m, hig_get_fid(f));
            // Get the center of the facet
            Point center;
            hig_get_facet_center(f, center);
            // Get the value for the velocity in this cell facet
            real val = ns->ed.eo.get_electroosmotic_source_term(center, dim, ns->par.t);
            // Set the velocity value for the velocity distributed property
            dp_set_value(ns->ed.eo.dpFeo[dim], flid, val);
        }
        // Destroying the iterator
        higfit_destroy(fit);
        // Sync initial values among processes
        dp_sync(ns->ed.eo.dpFeo[dim]);
    }
}

// Initialize the electro-osmotic phi
void higflow_initialize_electroosmotic_phi(higflow_solver *ns) {
    // Setting the cell iterator
    higcit_celliterator *it;
    // Getting the local domain
    sim_domain *sdp = psd_get_local_domain(ns->ed.eo.psdEOphi);
    // Getting the Mapper for the local domain
    mp_mapper *m = sd_get_domain_mapper(sdp);
    // Traversing the cells of local domain
    for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Getting the cell
        hig_cell *c = higcit_getcell(it);
        // Get the cell identifier
        int cgid = mp_lookup(m, hig_get_cid(c));
        // Get the center of the cell
        Point center;
        hig_get_center(c, center);
        // Get the value for electro-osmotic phi in this cell
        real val = ns->ed.eo.get_electroosmotic_phi(center, ns->par.t);
        // Set the value for pressure distributed property
        dp_set_value(ns->ed.eo.dpphi, cgid, val);
    }
    // Destroying the iterator
    higcit_destroy(it);
    // Sync initial values among processes
    dp_sync(ns->ed.eo.dpphi);
}

// Initialize the electro-osmotic psi
void higflow_initialize_electroosmotic_psi(higflow_solver *ns) {
    // Setting the cell iterator
    higcit_celliterator *it;
    // Getting the local domain
    sim_domain *sdp = psd_get_local_domain(ns->ed.eo.psdEOpsi);
    // Getting the Mapper for the local domain
    mp_mapper *m = sd_get_domain_mapper(sdp);
    // Traversing the cells of local domain
    for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Getting the cell
        hig_cell *c = higcit_getcell(it);
        // Get the cell identifier
        int cgid = mp_lookup(m, hig_get_cid(c));
        // Get the center of the cell
        Point center;
        hig_get_center(c, center);
        // Get the value for electro-osmotic psi in this cell
        real val = ns->ed.eo.get_electroosmotic_psi(center, ns->par.t);
        // Set the value for pressure distributed property
        dp_set_value(ns->ed.eo.dppsi, cgid, val);
    }
    // Destroying the iterator
    higcit_destroy(it);
    // Sync initial values among processes
    dp_sync(ns->ed.eo.dppsi);
}

// Initialize the electro-osmotic nplus
void higflow_initialize_electroosmotic_nplus(higflow_solver *ns) {
    // Setting the cell iterator
    higcit_celliterator *it;
    // Getting the local domain
    sim_domain *sdp = psd_get_local_domain(ns->ed.eo.psdEOnplus);
    // Getting the Mapper for the local domain
    mp_mapper *m = sd_get_domain_mapper(sdp);
    // Traversing the cells of local domain
    for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Getting the cell
        hig_cell *c = higcit_getcell(it);
        // Get the cell identifier
        int cgid = mp_lookup(m, hig_get_cid(c));
        // Get the center of the cell
        Point center;
        hig_get_center(c, center);
        // Get the value for electro-osmotic nplus in this cell
        real val = ns->ed.eo.get_electroosmotic_nplus(center, ns->par.t);
        // Set the value for pressure distributed property
        dp_set_value(ns->ed.eo.dpnplus, cgid, val);
    }
    // Destroying the iterator
    higcit_destroy(it);
    // Sync initial values among processes
    dp_sync(ns->ed.eo.dpnplus);
}

// Initialize the electro-osmotic nminus
void higflow_initialize_electroosmotic_nminus(higflow_solver *ns) {
    // Setting the cell iterator
    higcit_celliterator *it;
    // Getting the local domain
    sim_domain *sdp = psd_get_local_domain(ns->ed.eo.psdEOnminus);
    // Getting the Mapper for the local domain
    mp_mapper *m = sd_get_domain_mapper(sdp);
    // Traversing the cells of local domain
    for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Getting the cell
        hig_cell *c = higcit_getcell(it);
        // Get the cell identifier
        int cgid = mp_lookup(m, hig_get_cid(c));
        // Get the center of the cell
        Point center;
        hig_get_center(c, center);
        // Get the value for electro-osmotic nminus in this cell
        real val = ns->ed.eo.get_electroosmotic_nminus(center, ns->par.t);
        // Set the value for pressure distributed property
        dp_set_value(ns->ed.eo.dpnminus, cgid, val);
    }
    // Destroying the iterator
    higcit_destroy(it);
    // Sync initial values among processes
    dp_sync(ns->ed.eo.dpnminus);
}

// Initialize the viscosity for viscoelastic flows with variable viscosity
void higflow_initialize_viscosity_vevv(higflow_solver *ns) {
    // Setting the cell iterator
    higcit_celliterator *it;
    // Getting the local domain (viscosity)
    sim_domain *sdvisc = psd_get_local_domain(ns->ed.vevv.psdVisc);
    // Getting the Mapper for the local domain 
    mp_mapper *m = sd_get_domain_mapper(sdvisc);
    // Traversing the cells of local domain
    for(it = sd_get_domain_celliterator(sdvisc); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Getting the cell
        hig_cell *c = higcit_getcell(it);
        // Get the cell identifier
        int cgid = mp_lookup(m, hig_get_cid(c));
        // Get the center of the cell
        Point center;
        hig_get_center(c, center);
        real val;
       if (ns->contr.rheotype == PLM) {
           val = ns->ed.vevv.get_viscosity(center, 0.0, ns->par.t, ns->ed.vevv.par.beta, 0.0);
       }
       if (ns->contr.rheotype == THIXOTROPIC) {
           real valstruct;
           valstruct = ns->ed.vevv.get_structpar(center, 0.0, ns->par.t, ns->ed.vevv.par.beta, ns->ed.vevv.par.Phi, ns->ed.vevv.par.Lambda, ns->ed.vevv.par.Gamma);
           val = ns->ed.vevv.get_viscosity(center, 0.0, ns->par.t, ns->ed.vevv.par.beta, valstruct);
       }
        // Set the value for viscosity distributed property
        dp_set_value(ns->ed.vevv.dpvisc, cgid, val);
    }
    // Destroying the iterator
    higcit_destroy(it);
    // Sync initial values among processes
    dp_sync(ns->ed.vevv.dpvisc);
}

// Initialize the structural parameter for viscoelastic flows
void higflow_initialize_structural_parameter(higflow_solver *ns) {
    //Get the BMP model parameters
    real Lambda = ns->ed.vevv.par.Lambda;
    real Phi = ns->ed.vevv.par.Phi;
    real Gamma = ns->ed.vevv.par.Gamma;
    real beta = ns->ed.vevv.par.beta;
    // Setting the cell iterator
    higcit_celliterator *it;
    // Getting the local domain 
    sim_domain *sdstructpar = psd_get_local_domain(ns->ed.vevv.psdVisc);
    // Getting the Mapper for the local domain 
    mp_mapper *m = sd_get_domain_mapper(sdstructpar);
    // Traversing the cells of local domain
    for(it = sd_get_domain_celliterator(sdstructpar); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Getting the cell
        hig_cell *c = higcit_getcell(it);
        // Get the cell identifier
        int cgid = mp_lookup(m, hig_get_cid(c));
        // Get the center of the cell
        Point center;
        hig_get_center(c, center);
        // Get the value for structural parameter in this cell
        real val = ns->ed.vevv.get_structpar(center, 0.0, ns->par.t, beta, Phi, Lambda, Gamma);
        // Set the value for structural parameter distributed property
        dp_set_value(ns->ed.vevv.dpStructPar, cgid, val);
    }
    // Destroying the iterator
    higcit_destroy(it);
    // Sync initial values among processes
    dp_sync(ns->ed.vevv.dpStructPar);
}

// Initialize the Non-Newtonian Tensor for viscoelastic flows with variable viscosity
void higflow_initialize_viscoelastic_tensor_variable_viscosity(higflow_solver *ns) {
    // Non Newtonian flow
    if (ns->contr.flowtype == 6) {
        // Setting the cell iterator
        higcit_celliterator *it;
        // Getting the local domain
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Getting the Mapper for the local domain
        mp_mapper *m = sd_get_domain_mapper(sdp);
        // Traversing the cells of local domain
        for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Getting the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int cgid = mp_lookup(m, hig_get_cid(c));
            // Get the center of the cell
            Point center;
            hig_get_center(c, center);
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get the value for the tensor in this cell
                    real val = ns->ed.vevv.get_tensor(center, i, j, ns->par.t);
                    // Set the value for tensor distributed property
                    dp_set_value(ns->ed.vevv.dpS[i][j], cgid, val);
                }
            }
        }
        // Destroying the iterator
        higcit_destroy(it);
        // Sync initial values among processes
	for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.vevv.dpS[i][j]);
	    }
	}
    }
}

// Initialize the shear-banding density number nA
void higflow_initialize_shear_banding_nA(higflow_solver *ns) {
    // Setting the cell iterator
    higcit_celliterator *it;
    // Getting the local domain
    sim_domain *sdp = psd_get_local_domain(ns->ed.vesb.psdSBnA);
    // Getting the Mapper for the local domain
    mp_mapper *m = sd_get_domain_mapper(sdp);
    // Traversing the cells of local domain
    for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Getting the cell
        hig_cell *c = higcit_getcell(it);
        // Get the cell identifier
        int cgid = mp_lookup(m, hig_get_cid(c));
        // Get the center of the cell
        Point center;
        hig_get_center(c, center);
        // Get the value for nA in this cell
        real val = ns->ed.vesb.get_nA(center, ns->par.t);
        // Set the value for distributed property
        dp_set_value(ns->ed.vesb.dpnA, cgid, val);
    }
    // Destroying the iterator
    higcit_destroy(it);
    // Sync initial values among processes
    dp_sync(ns->ed.vesb.dpnA);
}

// Initialize the shear-banding density number nB
void higflow_initialize_shear_banding_nB(higflow_solver *ns) {
    // Setting the cell iterator
    higcit_celliterator *it;
    // Getting the local domain
    sim_domain *sdp = psd_get_local_domain(ns->ed.vesb.psdSBnB);
    // Getting the Mapper for the local domain
    mp_mapper *m = sd_get_domain_mapper(sdp);
    // Traversing the cells of local domain
    for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Getting the cell
        hig_cell *c = higcit_getcell(it);
        // Get the cell identifier
        int cgid = mp_lookup(m, hig_get_cid(c));
        // Get the center of the cell
        Point center;
        hig_get_center(c, center);
        // Get the value for nB in this cell
        real val = ns->ed.vesb.get_nB(center, ns->par.t);
        // Set the value for distributed property
        dp_set_value(ns->ed.vesb.dpnB, cgid, val);
    }
    // Destroying the iterator
    higcit_destroy(it);
    // Sync initial values among processes
    dp_sync(ns->ed.vesb.dpnB);
}

// Initialize the concentration of specie A (cA)
void higflow_initialize_shear_banding_cA(higflow_solver *ns) {
    real CAeq = ns->ed.vesb.par.CAeq;
    real chi = ns->ed.vesb.par.chi;
    // Setting the cell iterator
    higcit_celliterator *it;
    // Getting the local domain 
    sim_domain *sdca = psd_get_local_domain(ns->ed.psdED);
    // Getting the Mapper for the local domain 
    mp_mapper *m = sd_get_domain_mapper(sdca);
    // Traversing the cells of local domain
    for(it = sd_get_domain_celliterator(sdca); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Getting the cell
        hig_cell *c = higcit_getcell(it);
        // Get the cell identifier
        int cgid = mp_lookup(m, hig_get_cid(c));
        // Get the center of the cell
        Point center;
        hig_get_center(c, center);
        // Get the value for pressure in this cell
        real val = ns->ed.vesb.get_cA(center, ns->par.t, CAeq, chi, 0.0);;
        // Set the value for pressure distributed property
        dp_set_value(ns->ed.vesb.dpcA, cgid, val);
    }
    // Destroying the iterator
    higcit_destroy(it);
    // Sync initial values among processes
    dp_sync(ns->ed.vesb.dpcA);
}

// Initialize the concentration of specie B (cB)
void higflow_initialize_shear_banding_cB(higflow_solver *ns) {
    real CBeq = ns->ed.vesb.par.CBeq;
    real chi = ns->ed.vesb.par.chi;
    // Setting the cell iterator
    higcit_celliterator *it;
    // Getting the local domain 
    sim_domain *sdcb = psd_get_local_domain(ns->ed.psdED);
    // Getting the Mapper for the local domain 
    mp_mapper *m = sd_get_domain_mapper(sdcb);
    // Traversing the cells of local domain
    for(it = sd_get_domain_celliterator(sdcb); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Getting the cell
        hig_cell *c = higcit_getcell(it);
        // Get the cell identifier
        int cgid = mp_lookup(m, hig_get_cid(c));
        // Get the center of the cell
        Point center;
        hig_get_center(c, center);
        // Get the value for pressure in this cell
        real val = ns->ed.vesb.get_cB(center, ns->par.t, CBeq, chi, 0.0);;
        // Set the value for pressure distributed property
        dp_set_value(ns->ed.vesb.dpcB, cgid, val);
    }
    // Destroying the iterator
    higcit_destroy(it);
    // Sync initial values among processes
    dp_sync(ns->ed.vesb.dpcB);
}


// Initialize the Non-Newtonian Tensor for viscoelastic flows with shear-banding
void higflow_initialize_viscoelastic_tensor_shear_banding(higflow_solver *ns) {
    // Non Newtonian flow
    if (ns->contr.flowtype == 7) {
        // Setting the cell iterator
        higcit_celliterator *it;
        // Getting the local domain
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Getting the Mapper for the local domain
        mp_mapper *m = sd_get_domain_mapper(sdp);
        // Traversing the cells of local domain
        for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Getting the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int cgid = mp_lookup(m, hig_get_cid(c));
            // Get the center of the cell
            Point center;
            hig_get_center(c, center);
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get the value for the tensor in this cell
                    real val = ns->ed.vesb.get_tensor(center, i, j, ns->par.t);
                    // Set the value for tensor distributed property
                    dp_set_value(ns->ed.vesb.dpS[i][j], cgid, val);
                }
            }
        }
        // Destroying the iterator
        higcit_destroy(it);
        // Sync initial values among processes
	for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.vesb.dpS[i][j]);
	    }
	}
    }
}

// Initialize the conformation tensor of specie A for viscoelastic flows with shear-banding
void higflow_initialize_conformation_tensor_A_shear_banding(higflow_solver *ns) {
    // Non Newtonian flow
    if (ns->contr.rheotype == VCM) {
        // Setting the cell iterator
        higcit_celliterator *it;
        // Getting the local domain
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Getting the Mapper for the local domain
        mp_mapper *m = sd_get_domain_mapper(sdp);
        // Traversing the cells of local domain
        for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Getting the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int cgid = mp_lookup(m, hig_get_cid(c));
            // Get the center of the cell
            Point center;
            hig_get_center(c, center);
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get the value for the tensor in this cell
                    real val = ns->ed.vesb.get_tensor_A(center, i, j, ns->par.t);
                    // Set the value for tensor distributed property
                    dp_set_value(ns->ed.vesb.dpA[i][j], cgid, val);
                }
            }
        }
        // Destroying the iterator
        higcit_destroy(it);
        // Sync initial values among processes
	for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.vesb.dpA[i][j]);
	    }
	}
    }
}

// Initialize the conformation tensor of specie B for viscoelastic flows with shear-banding
void higflow_initialize_conformation_tensor_B_shear_banding(higflow_solver *ns) {
    // Non Newtonian flow
    if (ns->contr.rheotype == VCM) {
        // Setting the cell iterator
        higcit_celliterator *it;
        // Getting the local domain
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Getting the Mapper for the local domain
        mp_mapper *m = sd_get_domain_mapper(sdp);
        // Traversing the cells of local domain
        for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Getting the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int cgid = mp_lookup(m, hig_get_cid(c));
            // Get the center of the cell
            Point center;
            hig_get_center(c, center);
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get the value for the tensor in this cell
                    real val = ns->ed.vesb.get_tensor_B(center, i, j, ns->par.t);
                    // Set the value for tensor distributed property
                    dp_set_value(ns->ed.vesb.dpB[i][j], cgid, val);
                }
            }
        }
        // Destroying the iterator
        higcit_destroy(it);
        // Sync initial values among processes
	for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.vesb.dpB[i][j]);
	    }
	}
    }
}

// Initialize the Non-Newtonian Tensor
void higflow_initialize_elastoviscoplastic_tensor(higflow_solver *ns) {
    // Non Newtonian flow
    if (ns->contr.flowtype == 8) {
        // Setting the cell iterator
        higcit_celliterator *it;
        // Getting the local domain
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Getting the Mapper for the local domain
        mp_mapper *m = sd_get_domain_mapper(sdp);
        // Traversing the cells of local domain
        for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Getting the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int cgid = mp_lookup(m, hig_get_cid(c));
            // Get the center of the cell
            Point center;
            hig_get_center(c, center);
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get the value for the tensor in this cell
                    real val = ns->ed.vepl.get_tensor(center, i, j, ns->par.t);
                    // Set the value for tensor distributed property
                    dp_set_value(ns->ed.vepl.dpS[i][j], cgid, val);
                }
            }
        }
        // Destroying the iterator
        higcit_destroy(it);
        // Sync initial values among processes
	for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.vepl.dpS[i][j]);
	    }
	}
    }
}

// Initialize the Tensor S for shear-thickening suspension
void higflow_initialize_shear_thickening_suspension_tensor(higflow_solver *ns) {
    // Non Newtonian flow
    if (ns->contr.flowtype == 9) {
        // Setting the cell iterator
        higcit_celliterator *it;
        // Getting the local domain
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Getting the Mapper for the local domain
        mp_mapper *m = sd_get_domain_mapper(sdp);
        // Traversing the cells of local domain
        for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Getting the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int cgid = mp_lookup(m, hig_get_cid(c));
            // Get the center of the cell
            Point center;
            hig_get_center(c, center);
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get the value for the tensor in this cell
                    real val = ns->ed.stsp.get_tensor(center, i, j, ns->par.t);
                    // Set the value for tensor distributed property
                    dp_set_value(ns->ed.stsp.dpS[i][j], cgid, val);
                }
            }
        }
        // Destroying the iterator
        higcit_destroy(it);
        // Sync initial values among processes
	for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.stsp.dpS[i][j]);
	    }
	}
    }
}

// Initialize the microstructure tensor A for shear-thickening suspension
void higflow_initialize_shear_thickening_suspension_microstructure_tensor(higflow_solver *ns) {
    // Non Newtonian flow
    if (ns->contr.flowtype == 9) {
        // Setting the cell iterator
        higcit_celliterator *it;
        // Getting the local domain
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Getting the Mapper for the local domain
        mp_mapper *m = sd_get_domain_mapper(sdp);
        // Traversing the cells of local domain
        for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Getting the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int cgid = mp_lookup(m, hig_get_cid(c));
            // Get the center of the cell
            Point center;
            hig_get_center(c, center);
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get the value for the tensor in this cell
                    real val = ns->ed.stsp.get_tensor_A(center, i, j, ns->par.t);
                    // Set the value for tensor distributed property
                    dp_set_value(ns->ed.stsp.dpA[i][j], cgid, val);
                }
            }
        }
        // Destroying the iterator
        higcit_destroy(it);
        // Sync initial values among processes
	for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.stsp.dpA[i][j]);
	    }
	}
    }
}

// Initialize the volume fraction for shear thickening suspensions (with particle migration)
void higflow_initialize_volume_fraction(higflow_solver *ns) {
    // Setting the cell iterator
    higcit_celliterator *it;
    // Getting the local domain 
    sim_domain *sdphi = psd_get_local_domain(ns->ed.stsp.psdphi);
    // Getting the Mapper for the local domain 
    mp_mapper *m = sd_get_domain_mapper(sdphi);
    // Traversing the cells of local domain
    for(it = sd_get_domain_celliterator(sdphi); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Getting the cell
        hig_cell *c = higcit_getcell(it);
        // Get the cell identifier
        int cgid = mp_lookup(m, hig_get_cid(c));
        // Get the center of the cell
        Point center;
        hig_get_center(c, center);
        // Get the value for structural parameter in this cell
        real val = ns->ed.stsp.get_vol_frac(center, ns->par.t);
        //printf("===> volfrac = %lf <===\n", val);
        // Set the value for structural parameter distributed property
        dp_set_value(ns->ed.stsp.dpphi, cgid, val);
    }
    // Destroying the iterator
    higcit_destroy(it);
    // Sync initial values among processes
    dp_sync(ns->ed.stsp.dpphi);
    //printf("=+=+=+= WE ARE HERE AFTER initialising volfrac =+=+=+=\n");
}

// Initialize the velocities
void higflow_initialize_velocity(higflow_solver *ns) {
    // Setting the facet-cell iterator
    higfit_facetiterator *fit;
    // Setting the volocity values for the domain
    sim_facet_domain *sfdu[DIM];
    // Setting the velocity U(dim)
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local domain property
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the Mapper for the local domain
        mp_mapper *m = sfd_get_domain_mapper(sfdu[dim]);
        for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Getting the cell
            hig_facet *f = higfit_getfacet(fit);
            // Get the cell identifier
            int flid = mp_lookup(m, hig_get_fid(f));
            // Get the center of the facet
            Point center;
            hig_get_facet_center(f, center);
            // Get the value for the velocity in this cell facet
            real val = ns->func.get_velocity(center, dim, ns->par.t);
            // Set the velocity value for the velocity distributed property
            dp_set_value(ns->dpu[dim], flid, val);
        }
        // Destroying the iterator
        higfit_destroy(fit);
        // Sync initial values among processes
        dp_sync(ns->dpu[dim]);
    }
}

// Initialize the facet source term
void higflow_initialize_facet_source_term(higflow_solver *ns) {
    // Setting the facet-cell iterator
    higfit_facetiterator *fit;
    // Setting the volocity values for the domain
    sim_facet_domain *sfdu[DIM];
    // Setting the velocity U(dim)
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local domain property
        sfdu[dim] = psfd_get_local_domain(ns->psfdF[dim]);
        // Get the Mapper for the local domain 
        mp_mapper *m = sfd_get_domain_mapper(sfdu[dim]);
        for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Getting the cell
            hig_facet *f = higfit_getfacet(fit);
            // Get the cell identifier
            int flid = mp_lookup(m, hig_get_fid(f));
            // Get the center of the facet
            Point center;
            hig_get_facet_center(f, center);
            // Get the value for the facet source term in this cell facet
            real val = ns->func.get_facet_source_term(center, dim, ns->par.t);
            // Set the value for the facet source term distributed property
            dp_set_value(ns->dpFU[dim], flid, val);
        }
        // Destroying the iterator
        higfit_destroy(fit);
        // Sync initial values among processes
        dp_sync(ns->dpFU[dim]);
    }
}

// Initialize distributed properties
void higflow_initialize_distributed_properties(higflow_solver *ns) {
    // Initialize pressure distributed property
    higflow_initialize_pressure(ns);
    // Initialize velocity distributed property
    higflow_initialize_velocity(ns);
    // Initialize cell source term distributed property
    higflow_initialize_cell_source_term(ns);
    // Initialize facet source term distributed property
    higflow_initialize_facet_source_term(ns);
    if (ns->contr.flowtype == 1) {
        // Initialize viscosity distributed property
        higflow_initialize_viscosity_gn(ns);
    } else if (ns->contr.flowtype == 2) {
        // Initialize volume fraction distributed property
        higflow_initialize_fracvol(ns);
        // Initialize viscosity distributed property
        higflow_initialize_viscosity_mult(ns);
        // Initialize density distributed property
        higflow_initialize_density(ns);
    } else if (ns->contr.flowtype == 3) {
        // Initialize non Newtonian tensor distributed property
        higflow_initialize_viscoelastic_tensor(ns);
    } else if (ns->contr.flowtype == 4) {
        // Initialize non Newtonian integral tensor distributed property
        higflow_initialize_viscoelastic_integral_tensor(ns);
        // Initialize non Newtonian integral finger tensor distributed property
        higflow_initialize_viscoelastic_integral_finger_tensor(ns);
    } else if (ns->contr.flowtype == 6) {
        if (ns->contr.rheotype == THIXOTROPIC) {
            //Initialize structural parameter
            higflow_initialize_structural_parameter(ns);
        }
        //Initialize viscosity distributed property for viscoelastic flows with variable viscosity
        higflow_initialize_viscosity_vevv(ns);
        //Initialize non Newtonian and viscoelastic tensor distributed property
        higflow_initialize_viscoelastic_tensor_variable_viscosity(ns);
    } else if (ns->contr.flowtype == 7) {
       if (ns->contr.rheotype == VCM) {
           //Initialize the density number nA
           higflow_initialize_shear_banding_nA(ns);
           //Initialize the density number nB
           higflow_initialize_shear_banding_nB(ns);
           // Initialize the concentration of specie A (cA)
           higflow_initialize_shear_banding_cA(ns);
           // Initialize the concentration of specie B (cB)
           higflow_initialize_shear_banding_cB(ns);
           // Initialize the conformation tensor of specie A for viscoelastic flows with shear-banding
           higflow_initialize_conformation_tensor_A_shear_banding(ns);
           // Initialize the conformation tensor of specie B for viscoelastic flows with shear-banding
           higflow_initialize_conformation_tensor_B_shear_banding(ns);
       }
       //Initialize the viscoelastic tensor distributed property    
       higflow_initialize_viscoelastic_tensor_shear_banding(ns);
    } else if (ns->contr.flowtype == 8) {
       // Initialize non Newtonian tensor distributed property
       higflow_initialize_elastoviscoplastic_tensor(ns);
    } else if (ns->contr.flowtype == 9) {
        //only for the model that considers particle migration
        //printf("=+=+=+= WE ARE HERE IN INITIALIZE DP =+=+=+=\n");
        if (ns->ed.stsp.contr.model == GW_WC_IF) {
            //Initialize structural parameter
            higflow_initialize_volume_fraction(ns);
            //printf("=+=+=+= WE ARE HERE after =+=+=+=\n");
        }
        // Initialize shear-thickening suspension tensor distributed properties
        higflow_initialize_shear_thickening_suspension_tensor(ns);
        higflow_initialize_shear_thickening_suspension_microstructure_tensor(ns);
    }
    if (ns->contr.modelflowtype == 1) {
        // Initialize electro-osmotic phi distributed property
        higflow_initialize_electroosmotic_phi(ns);
        // Initialize electro-osmotic psi distributed property
        higflow_initialize_electroosmotic_psi(ns);
        // Initialize electro-osmotic nplus distributed property
        higflow_initialize_electroosmotic_nplus(ns);
        // Initialize electro-osmotic nminus distributed property
        higflow_initialize_electroosmotic_nminus(ns);
    }
    //Non-isothermal flows
    if (ns->contr.modelflowtype == 5) {
        // Initialize electro-osmotic phi distributed property
        higflow_initialize_temperature(ns);
        //Initialize the energy source term
        higflow_initialize_energy_source_term(ns);
    }
}

