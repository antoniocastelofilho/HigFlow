// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Initial Condition - version 10/11/2016
// *******************************************************************
// *******************************************************************

#include "hig-flow-ic.h"
#include <libfyaml.h>

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
    }
    if (ns->contr.modelflowtype == 1) {
        // Initialize electro-osmotic domain
        higflow_create_partitioned_domain_electroosmotic(ns, pg, order);
        // Creating the stencil for properties interpolation
        higflow_create_stencil_for_extra_domain(ns);
    }
}

// Navier-Stokes initialize the domain yaml
void higflow_initialize_domain_yaml(higflow_solver *ns, int ntasks, int myrank, int order) {
    
	 // Loading the boundary condition data
    char namefile[1024];
    sprintf(namefile,"%s.domain.yaml",ns->par.nameload);
    
    FILE *fdomain = fopen(namefile, "r");
    struct fy_document *fyd = NULL;
    fyd = fy_document_build_from_file(NULL, namefile);
     
    if (fyd == NULL) {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n",namefile);
        exit(1);
    }
	  
    // Number of HigTrees
    int numhigs;
    int ifd = fy_document_scanf(fyd,"/domain/number_domain %d",&numhigs);
    higio_amr_info *mi[numhigs];
    for(int h = 0; h < numhigs; h++) {
        char atrib[1024], amrfilename[1024];
		  sprintf(atrib,"/domain/domain%d/id %%d",h);
        //ifd = fy_document_scanf(fyd,atrib,&(id[h]));
		  
		  // Name of the HigTree file
        //char *amrfilename = argv[1]; argv++;
        
        sprintf(atrib,"/domain/domain%d/path %%s",h);
        ifd = fy_document_scanf(fyd,atrib,amrfilename);

        FILE *fd = fopen(amrfilename, "r");
        mi[h] = higio_read_amr_info(fd);
        fclose(fd);
    }
    fy_document_destroy(fyd);
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
    // the partitioned sub-domain for extra properties
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
    }
    if (ns->contr.modelflowtype == 1) {
        // Initialize electro-osmotic domain
        higflow_create_partitioned_domain_electroosmotic(ns, pg, order);
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

// Initialize the viscosity - Multiphase
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

// Initialize the Non-Newtonian Tensor - multiphase
void higflow_initialize_viscoelastic_mult_tensor(higflow_solver *ns) {
    // Non Newtonian flow - multiphase
    if (ns->contr.flowtype == 2) {
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
            // Calculate the density
            real fracvol  = compute_value_at_point(sddens, center, center, 1.0, ns->ed.mult.dpfracvol, ns->stn);
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get the value for the tensor in this cell
                    real S0 = ns->ed.mult.get_tensor0(center, i, j, ns->par.t);
                    real S1 = ns->ed.mult.get_tensor1(center, i, j, ns->par.t);
                    if (fracvol == 0.0){
                       S1 = 0.0;
                    }
                    if (fracvol == 1.0){
                       S0 = 0.0;
                    }
                    // Set the value for tensor distributed property - phase 0
                    dp_set_value(ns->ed.mult.dpS0[i][j], cgid, S0);
                    // Set the value for tensor distributed property - phase 1
                    dp_set_value(ns->ed.mult.dpS1[i][j], cgid, S1);
                }
            }
        }
        // Destroying the iterator
        higcit_destroy(it);
        // Sync initial values among processes
	     for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.mult.dpS0[i][j]);
                dp_sync(ns->ed.mult.dpS1[i][j]);
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
            int fgid = mp_lookup(m, hig_get_fid(f));
            // Get the center of the facet
            Point center;
            hig_get_facet_center(f, center);
            // Get the value for the velocity in this cell facet
            real val = ns->ed.eo.get_electroosmotic_source_term(center, dim, ns->par.t);
            // Set the velocity value for the velocity distributed property
            dp_set_value(ns->ed.eo.dpFeo[dim], fgid, val);
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
            int fgid = mp_lookup(m, hig_get_fid(f));
            // Get the center of the facet
            Point center;
            hig_get_facet_center(f, center);
            // Get the value for the velocity in this cell facet
            real val = ns->func.get_velocity(center, dim, ns->par.t);
            // Set the velocity value for the velocity distributed property
            dp_set_value(ns->dpu[dim], fgid, val);
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
            int fgid = mp_lookup(m, hig_get_fid(f));
            // Get the center of the facet
            Point center;
            hig_get_facet_center(f, center);
            // Get the value for the facet source term in this cell facet
            real val = ns->func.get_facet_source_term(center, dim, ns->par.t);
            // Set the value for the facet source term distributed property
            dp_set_value(ns->dpFU[dim], fgid, val);
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
        // Initialize Tensor distributed property
        higflow_initialize_viscoelastic_mult_tensor(ns);
    } else if (ns->contr.flowtype == 3) {
        // Initialize non Newtonian tensor distributed property
        higflow_initialize_viscoelastic_tensor(ns);
    } else if (ns->contr.flowtype == 4) {
        // Initialize non Newtonian integral tensor distributed property
        higflow_initialize_viscoelastic_integral_tensor(ns);
        // Initialize non Newtonian integral finger tensor distributed property
        higflow_initialize_viscoelastic_integral_finger_tensor(ns);
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
}

