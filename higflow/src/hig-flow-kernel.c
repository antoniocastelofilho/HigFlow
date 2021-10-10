// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Kernel - version 10/11/2016
// *******************************************************************
// *******************************************************************

#include "hig-flow-kernel.h"

// *******************************************************************
// Navier-Stokes Create and Destroy Object
// *******************************************************************

// Create the NS object
higflow_solver *higflow_create (void) {
    DECL_AND_ALLOC(higflow_solver, ns, 1);	
    return ns;
}

// Destroy the NS object
void higflow_destroy (higflow_solver *ns) {
    // Destroy the distributed properties
    // Pressure
    dp_destroy(ns->dpp);
    // Pressure difference
    dp_destroy(ns->ddeltap);
    for(int dim = 0; dim < DIM; dim++) {
        // Final velocity
        dp_destroy(ns->dpu[dim]);
        // Intermediate velocity
        dp_destroy(ns->dpustar[dim]);
        // Auxiliar intermediate velocity for Runge-Kutta method
        dp_destroy(ns->dpuaux[dim]);
        // Facet source term
        dp_destroy(ns->dpFU[dim]);
    }
    // Cell source term 
    dp_destroy(ns->dpF);
    // Non Newtonian tensor
    switch (ns->contr.flowtype) {
        case 1:
            // Generalized Newtonian
            dp_destroy(ns->ed.gn.dpvisc);
            for (int i = 0; i < DIM; i++)
                for (int j = 0; j < DIM; j++) {
                    // Tensor terms destroy
                    dp_destroy(ns->ed.gn.dpD[i][j]);
                }
            // Destroy the stencil for extra domains
            stn_destroy(ns->ed.stn);
            break;
        case 2:
            // Multiphase
            dp_destroy(ns->ed.mult.dpvisc);
            dp_destroy(ns->ed.mult.dpdens);
            dp_destroy(ns->ed.mult.dpfracvol);
            dp_destroy(ns->ed.mult.dpfracvolaux);
            dp_destroy(ns->ed.mult.dpcurvature);
            dp_destroy(ns->ed.mult.dpdistance);
            for (int i = 0; i < DIM; i++){
               dp_destroy(ns->ed.mult.dpIF[i]);
               dp_destroy(ns->ed.mult.dpnormal[i]);
            }
            // Destroy the beta viscoelastic
            dp_destroy(ns->ed.mult.dpbeta);
            for (int i = 0; i < DIM; i++)
                for (int j = 0; j < DIM; j++) {
                    // Tensor terms destroy
                    dp_destroy(ns->ed.mult.dpS0[i][j]);
                    dp_destroy(ns->ed.mult.dpS1[i][j]);
                    dp_destroy(ns->ed.mult.dpKernel0[i][j]);
                    dp_destroy(ns->ed.mult.dpKernel1[i][j]);
                }
            
            // Destroy the stencil for extra domains
            stn_destroy(ns->ed.stn);
            break;
        case 3:
            // Viscoelastic
            for (int i = 0; i < DIM; i++)
                for (int j = 0; j < DIM; j++) {
                    // Tensor terms destroy
                    dp_destroy(ns->ed.ve.dpD[i][j]);
                    dp_destroy(ns->ed.ve.dpS[i][j]);
                    dp_destroy(ns->ed.ve.dpKernel[i][j]);
                }
            // Destroy the stencil for extra domains
            stn_destroy(ns->ed.stn);
            break;
         case 4:
           // Viscoelastic integral model
           for (int i = 0; i < DIM; i++)
               for (int j = 0; j < DIM; j++) {
                   // Tensor terms destroy
                   dp_destroy(ns->ed.im.dpD[i][j]);
                   dp_destroy(ns->ed.im.dpS[i][j]);
                   for (int k = 0; k <= NDT; k++)
                       dp_destroy(ns->ed.im.dpB[k][i][j]);
               }
           // Destroy the stencil for extra domains
           stn_destroy(ns->ed.stn);
           break;   
    }
    if (ns->contr.modelflowtype == 1) {
        // Destroy the distributed properties for Electro-osmotic model
        for(int dim = 0; dim < DIM; dim++) {
            dp_destroy(ns->ed.eo.dpFeo[dim]);
        }
        dp_destroy(ns->ed.eo.dpphi);
        dp_destroy(ns->ed.eo.dppsi);
        dp_destroy(ns->ed.eo.dpnplus);
        dp_destroy(ns->ed.eo.dpnminus);
        // Destroy the stencil for extra domains
        stn_destroy(ns->ed.stn);
        // Destroy the solver for potential psi 
        slv_destroy(ns->ed.eo.slvpsi);
        slv_destroy(ns->ed.eo.slvphi);
        slv_destroy(ns->ed.eo.slvnplus);
        slv_destroy(ns->ed.eo.slvnminus);
    }
    // Destroy the stencil for pressure and velocities
    stn_destroy(ns->stn);
    // Destroy the stencil for source force
    stn_destroy(ns->stnF);
    // Destroy the solver for pressure
    slv_destroy(ns->slvp);
    // Destroy the solver for the velocity
    if (ns->contr.tempdiscrtype >= 3) {
        for (int i = 0; i < DIM; i++) {
            slv_destroy(ns->slvu[i]);
        }
    }
}

// *******************************************************************
// Navier-Stokes Initilize Linear System Solver and MPI
// *******************************************************************

// Initializing Navier-Stokes solver
void higflow_initialize(int *argc, char **argv[], int *myrank, int *ntasks) {
    // Initialize the Hig-Tree
    higtree_initialize(argc, argv);
    // Setting myrank and tasks
    MPI_Comm_rank(MPI_COMM_WORLD, myrank);
    MPI_Comm_size(MPI_COMM_WORLD, ntasks);
}

// *******************************************************************
// Navier-Stokes Create and Realloc the Linear System Solvers
// *******************************************************************

// Create the linear system solvers
void higflow_create_solver(higflow_solver *ns) {
    // Get the localdomainsize for cell center
    int localdomainsize = psd_get_local_domain_size(ns->psdp);
    // Creates a solver for pressure
    ns->slvp            = slv_create(SOLVER_ANY, psd_get_first_id(ns->psdp), localdomainsize);
    // Set the maximum of non zeros 
    slv_set_maxnonzeros(ns->slvp, 800);
    // Create the solver for the implicit methods
    if (ns->contr.tempdiscrtype >= 3) {
        for (int dim = 0; dim < DIM; dim++) {
            // Get the localdomainsize for facet center
            localdomainsize = psfd_get_local_domain_size(ns->psfdu[dim]);
            // Creates a solver for velocity
            ns->slvu[dim]   = slv_create(SOLVER_ANY, psfd_get_first_id(ns->psfdu[dim]), localdomainsize);
            // Set the maximum of non zeros 
            slv_set_maxnonzeros(ns->slvu[dim], 800);
        }
    }
    if (ns->contr.modelflowtype == 1) {
        // Get the localdomainsize for cell center
        int localdomainsizepsi = psd_get_local_domain_size(ns->ed.eo.psdEOpsi);
        // Creates a solver for pressure
        ns->ed.eo.slvpsi            = slv_create(SOLVER_ANY, psd_get_first_id(ns->ed.eo.psdEOpsi), localdomainsizepsi);
        // Set the maximum of non zeros 
        slv_set_maxnonzeros(ns->ed.eo.slvpsi, 800);
        // Get the localdomainsize for cell center
        int localdomainsizephi = psd_get_local_domain_size(ns->ed.eo.psdEOphi);
        // Creates a solver for pressure
        ns->ed.eo.slvphi            = slv_create(SOLVER_ANY, psd_get_first_id(ns->ed.eo.psdEOphi), localdomainsizephi);
        // Set the maximum of non zeros 
        slv_set_maxnonzeros(ns->ed.eo.slvphi, 800);
        // Get the localdomainsize for cell center
        int localdomainsizenplus = psd_get_local_domain_size(ns->ed.eo.psdEOnplus);
        // Creates a solver for pressure
        ns->ed.eo.slvnplus            = slv_create(SOLVER_ANY, psd_get_first_id(ns->ed.eo.psdEOnplus), localdomainsizenplus);
        // Set the maximum of non zeros 
        slv_set_maxnonzeros(ns->ed.eo.slvnplus, 800);
        // Get the localdomainsize for cell center
        int localdomainsizenminus = psd_get_local_domain_size(ns->ed.eo.psdEOnminus);
        // Creates a solver for pressure
        ns->ed.eo.slvnminus            = slv_create(SOLVER_ANY, psd_get_first_id(ns->ed.eo.psdEOnminus), localdomainsizenminus);
        // Set the maximum of non zeros 
        slv_set_maxnonzeros(ns->ed.eo.slvnminus, 800);
    }
}

// Realloc the linear system solvers
void higflow_realloc_solver(higflow_solver *ns) {
    // Destroy the solver for pressure
    slv_destroy(ns->slvp);
    // Get the localdomainsize for cell center
    int localdomainsize = psd_get_local_domain_size(ns->psdp);
    // Creates a solver for pressure
    ns->slvp            = slv_create(SOLVER_ANY, psd_get_first_id(ns->psdp), localdomainsize);
    // Set the maximum of non zeros 
    slv_set_maxnonzeros(ns->slvp, 800);
    // Realloc the solver for the implicit methods
    if (ns->contr.tempdiscrtype >= 3) {
        for (int dim = 0; dim < DIM; dim++) {
            // Destroy the solver for the velocity
            slv_destroy(ns->slvu[dim]);
            // Get the localdomainsize for facet center
            localdomainsize = psfd_get_local_domain_size(ns->psfdu[dim]);
            // Creates a solver for velocity
            ns->slvu[dim]   = slv_create(SOLVER_ANY, psfd_get_first_id(ns->psfdu[dim]), localdomainsize);
            // Set the maximum of non zeros 
            slv_set_maxnonzeros(ns->slvu[dim], 800);
        }
    }
    if (ns->contr.modelflowtype == 1) {
        // Destroy the solver for pressure
        slv_destroy(ns->ed.eo.slvpsi);
        // Get the localdomainsize for cell center
        int localdomainsizepsi = psd_get_local_domain_size(ns->ed.eo.psdEOpsi);
        // Creates a solver for pressure
        ns->ed.eo.slvpsi            = slv_create(SOLVER_ANY, psd_get_first_id(ns->ed.eo.psdEOpsi), localdomainsizepsi);
        // Set the maximum of non zeros 
        slv_set_maxnonzeros(ns->ed.eo.slvpsi, 800);
        // Destroy the solver for pressure
        slv_destroy(ns->ed.eo.slvphi);
        // Get the localdomainsize for cell center
        int localdomainsizephi = psd_get_local_domain_size(ns->ed.eo.psdEOphi);
        // Creates a solver for pressure
        ns->ed.eo.slvphi            = slv_create(SOLVER_ANY, psd_get_first_id(ns->ed.eo.psdEOpsi), localdomainsizephi);
        // Set the maximum of non zeros 
        slv_set_maxnonzeros(ns->ed.eo.slvphi, 800);
        // Destroy the solver for pressure
        slv_destroy(ns->ed.eo.slvnplus);
        // Get the localdomainsize for cell center
        int localdomainsizenplus = psd_get_local_domain_size(ns->ed.eo.psdEOnplus);
        // Creates a solver for pressure
        ns->ed.eo.slvnplus            = slv_create(SOLVER_ANY, psd_get_first_id(ns->ed.eo.psdEOnplus), localdomainsizenplus);
        // Set the maximum of non zeros 
        slv_set_maxnonzeros(ns->ed.eo.slvnplus, 800);
        // Destroy the solver for pressure
        slv_destroy(ns->ed.eo.slvnminus);
        // Get the localdomainsize for cell center
        int localdomainsizenminus = psd_get_local_domain_size(ns->ed.eo.psdEOnminus);
        // Creates a solver for pressure
        ns->ed.eo.slvnminus            = slv_create(SOLVER_ANY, psd_get_first_id(ns->ed.eo.psdEOnminus), localdomainsizenminus);
        // Set the maximum of non zeros 
        slv_set_maxnonzeros(ns->ed.eo.slvnminus, 800);
    }
}

// *******************************************************************
// Navier-Stokes Create Domain and Distributed Properties 
// *******************************************************************

// Create the simulation domain for NS object
void higflow_create_domain (higflow_solver *ns, int cache, int order) {
    // simulation domain SD (pressure)
    ns->sdp = sd_create(NULL);
    //reuse interpolation, 0 on, 1 off
    sd_use_cache(ns->sdp, cache);      
    //Sets the order of the interpolation to bhe used for the SD (pressure) 
    sd_set_interpolator_order(ns->sdp, order);
    // simulation domain SD (source force)
    ns->sdF = sd_create(NULL);
    //reuse interpolation, 0 on, 1 off
    sd_use_cache(ns->sdF, cache);      
    //Sets the order of the interpolation to bhe used for the SD (source force) 
    sd_set_interpolator_order(ns->sdF, order);
}

// Create the simulation domain for generalized newtonian flow
void higflow_create_domain_generalized_newtonian (higflow_solver *ns, int cache, int order, real (*get_viscosity)(Point center, real q, real t)) {
    if (ns->contr.flowtype == 1) {
       // simulation domain (SD) extra domain
       ns->ed.sdED = sd_create(NULL);
       //reuse interpolation, 0 on, 1 off
       sd_use_cache(ns->ed.sdED, cache);      
       //Sets the order of the interpolation to bhe used for the SD. 
       sd_set_interpolator_order(ns->ed.sdED, order);
       // function for the domain
       ns->ed.gn.get_viscosity = get_viscosity;
    }
}

// Create the simulation domain for multiphase flow
void higflow_create_domain_multiphase(higflow_solver *ns, int cache, int order, 
real (*get_viscosity0)(Point center, real t), 
real (*get_viscosity1)(Point center, real t), 
real (*get_density0)(Point center, real t),
real (*get_density1)(Point center, real t),
real (*get_fracvol)(Point center, Point delta, real t),
real (*get_tensor0)(Point center, int i, int j, real t),
real (*get_tensor1)(Point center, int i, int j, real t),
real (*get_kernel)(int dim, real lambda, real tol),
real (*get_kernel_inverse)(int dim, real lambda, real tol),
real (*get_kernel_jacobian)(int dim, real lambda, real tol)) {
    if (ns->contr.flowtype == 2) {
       // simulation domain (SD) extra domain
       ns->ed.sdED = sd_create(NULL);
       //reuse interpolation, 0 on, 1 off
       sd_use_cache(ns->ed.sdED, cache);      
       //Sets the order of the interpolation to bhe used for the SD. 
       sd_set_interpolator_order(ns->ed.sdED, order);
       // function for the domain
       ns->ed.mult.get_viscosity0      = get_viscosity0;
       // function for the domain
       ns->ed.mult.get_viscosity1      = get_viscosity1;
       // function for the domain
       ns->ed.mult.get_density0        = get_density0;
       // function for the domain
       ns->ed.mult.get_density1        = get_density1;
       // function for the domain
       ns->ed.mult.get_fracvol         = get_fracvol;
       // function for the domain
       ns->ed.mult.get_tensor0         = get_tensor0;
       // function for the domain
       ns->ed.mult.get_tensor1         = get_tensor1;
       // function for the kernel transformation
       ns->ed.mult.get_kernel          = get_kernel;
       // function for the inverse kernel transformation
       ns->ed.mult.get_kernel_inverse  = get_kernel_inverse;
       // function for the kernel transformation jacobian
       ns->ed.mult.get_kernel_jacobian = get_kernel_jacobian;
    }
}

// Create the simulation domain for viscoelastic flow
void higflow_create_domain_viscoelastic(higflow_solver *ns, int cache, int order,
real (*get_tensor)(Point center, int i, int j, real t),
real (*get_kernel)(int dim, real lambda, real tol),
real (*get_kernel_inverse)(int dim, real lambda, real tol),
real (*get_kernel_jacobian)(int dim, real lambda, real tol)) {
    if ((ns->contr.flowtype == 3) || (ns->contr.flowtype == 2)) {
       // simulation domain (SD) extra domain
       ns->ed.sdED = sd_create(NULL);
       //reuse interpolation, 0 on, 1 off
       sd_use_cache(ns->ed.sdED, cache);      
       //Sets the order of the interpolation to bhe used for the SD. 
       sd_set_interpolator_order(ns->ed.sdED, order);
       // function for the domain
       ns->ed.ve.get_tensor          = get_tensor;
       // function for the kernel transformation
       ns->ed.ve.get_kernel          = get_kernel;
       // function for the inverse kernel transformation
       ns->ed.ve.get_kernel_inverse  = get_kernel_inverse;
       // function for the kernel transformation jacobian
       ns->ed.ve.get_kernel_jacobian = get_kernel_jacobian;
    }
}

// Define the user function for viscoelastic flow
void higflow_define_user_function_viscoelastic(higflow_solver *ns, 
void (*calculate_m_user)(real Re, real De, real beta, real tr, real lambda[DIM],real R[DIM][DIM], real M[DIM][DIM], real M_aux[DIM][DIM], real tol)) {
    if (ns->contr.flowtype == 3) {
       // function for the user viscoelastic model
       ns->ed.ve.calculate_m_user = calculate_m_user;
    }
}


// Create the simulation domain for viscoelastic flow integral model
void higflow_create_domain_viscoelastic_integral(higflow_solver *ns, int cache, int order,
real (*get_tensor)(Point center, int i, int j, real t)) {
   // simulation domain (SD) extra domain
   if (ns->contr.flowtype == 4) {
       ns->ed.sdED = sd_create(NULL);
       //reuse interpolation, 0 on, 1 off
       sd_use_cache(ns->ed.sdED, cache);      
       //Sets the order of the interpolation to bhe used for the SD. 
       sd_set_interpolator_order(ns->ed.sdED, order);
       // function for the domain
       ns->ed.im.get_tensor          = get_tensor;
    }
}

// Create the simulation domain for electroosmotic flow
void higflow_create_domain_electroosmotic (higflow_solver *ns, int cache, int order,
real (*get_electroosmotic_source_term)(Point center, int dim, real t),
real (*get_electroosmotic_phi)(Point center, real t),
real (*get_electroosmotic_psi)(Point center, real t),
real (*get_electroosmotic_nplus)(Point center, real t),
real (*get_electroosmotic_nminus)(Point center, real t),
real (*get_boundary_electroosmotic_source_term)(int id, Point center, int dim, real t),
real (*get_boundary_electroosmotic_phi)(int id, Point center, real t),
real (*get_boundary_electroosmotic_psi)(int id, Point center, real t),
real (*get_boundary_electroosmotic_nplus)(int id, Point center, real t),
real (*get_boundary_electroosmotic_nminus)(int id, Point center, real t)) {
    if (ns->contr.modelflowtype == 1) {
       // simulation domain (EO) for phi
       ns->ed.eo.sdEOphi = sd_create(NULL);
       //reuse interpolation, 0 on, 1 off
       sd_use_cache(ns->ed.eo.sdEOphi, cache);      
       //Sets the order of the interpolation to bhe used for the SD. 
       sd_set_interpolator_order(ns->ed.eo.sdEOphi, order);
       // simulation domain (EO) for psi
       ns->ed.eo.sdEOpsi = sd_create(NULL);
       //reuse interpolation, 0 on, 1 off
       sd_use_cache(ns->ed.eo.sdEOpsi, cache);      
       //Sets the order of the interpolation to bhe used for the SD. 
       sd_set_interpolator_order(ns->ed.eo.sdEOpsi, order);
       // simulation domain (EO) for nplus
       ns->ed.eo.sdEOnplus = sd_create(NULL);
       //reuse interpolation, 0 on, 1 off
       sd_use_cache(ns->ed.eo.sdEOnplus, cache);      
       //Sets the order of the interpolation to bhe used for the SD. 
       sd_set_interpolator_order(ns->ed.eo.sdEOnplus, order);
       // simulation domain (EO) for nminus
       ns->ed.eo.sdEOnminus = sd_create(NULL);
       //reuse interpolation, 0 on, 1 off
       sd_use_cache(ns->ed.eo.sdEOnminus, cache);      
       //Sets the order of the interpolation to bhe used for the SD. 
       sd_set_interpolator_order(ns->ed.eo.sdEOnminus, order);
       // function for the eo source term 
       ns->ed.eo.get_electroosmotic_source_term = get_electroosmotic_source_term;
       // function for the eo phi
       ns->ed.eo.get_electroosmotic_phi         = get_electroosmotic_phi;
       // function for the eo psi
       ns->ed.eo.get_electroosmotic_psi         = get_electroosmotic_psi ;
       // function for the eo positive ion concentration 
       ns->ed.eo.get_electroosmotic_nplus       = get_electroosmotic_nplus;
       // function for the eo negative ion concentration 
       ns->ed.eo.get_electroosmotic_nminus      = get_electroosmotic_nminus;
       // function for the boundary source term
       ns->ed.eo.get_boundary_electroosmotic_source_term = get_boundary_electroosmotic_source_term;
       // function for the boundary potential phi
       ns->ed.eo.get_boundary_electroosmotic_phi         = get_boundary_electroosmotic_phi;
       // function for the boundary potential psi
       ns->ed.eo.get_boundary_electroosmotic_psi         = get_boundary_electroosmotic_psi ;
       // function for the boundary concentration of positive ion
       ns->ed.eo.get_boundary_electroosmotic_nplus       = get_boundary_electroosmotic_nplus;
       // function for the boundary concentration of negative ion
       ns->ed.eo.get_boundary_electroosmotic_nminus      = get_boundary_electroosmotic_nminus;
    }
}

// Create the partitioned simulation sub-domain for NS object
void higflow_create_partitioned_domain (higflow_solver *ns, partition_graph *pg, int order) {
    // Creating the partitioned sub-domain to simulation for pressure
    ns->psdp = psd_create(ns->sdp, pg);
    // Synced mapper for pressure
    psd_synced_mapper(ns->psdp);
    // Creating the partitioned sub-domain to simulation for source force
    ns->psdF = psd_create(ns->sdF, pg);
    // Synced mapper for source force
    psd_synced_mapper(ns->psdF);
    // Linking the property to the facets 
    for (int dim = 0; dim < DIM; dim++) {
        // Creating the list of cells center in the local domain (velocities)
        ns->sfdu[dim] = sfd_create(NULL, dim);
        // Setting the order for the properties interpolation (velocities)
        sfd_set_interpolator_order(ns->sfdu[dim], order);   
        // Copying the list of cells center in the domain for pressure 
        sfd_copy_higtrees_from_center_domain(ns->sfdu[dim], ns->sdp);
        sfd_adjust_facet_ids(ns->sfdu[dim]);
        // Creating property for the facets (velocities)
        ns->psfdu[dim] = psfd_create(ns->sfdu[dim], ns->psdp);
        // Creating the list of cells center in the local domain (source force)
        ns->sfdF[dim] = sfd_create(NULL, dim);
        // Setting the order for the properties interpolation (source force)
        sfd_set_interpolator_order(ns->sfdF[dim], order);   
        // Copying the list of cells center in the domain for pressure 
        sfd_copy_higtrees_from_center_domain(ns->sfdF[dim], ns->sdp);
        sfd_adjust_facet_ids(ns->sfdF[dim]);
        // Creating property for the facets (source force)
        ns->psfdF[dim] = psfd_create(ns->sfdF[dim], ns->psdp);
    }
    for(int dim = 0; dim < DIM; dim++) {
        // Mapping the properties in the domain (source force)
        psfd_compute_sfbi(ns->psfdF[dim]);
        // Sync mapper for source force
        psfd_synced_mapper(ns->psfdF[dim]);
    }
}

// Create the partitioned simulation sub-domain for generalized newtonian simulation
void higflow_create_partitioned_domain_generalized_newtonian (higflow_solver *ns, partition_graph *pg, int order) {
    if (ns->contr.flowtype == 1) {
        // Creating the partitioned sub-domain to simulation
        ns->ed.psdED = psd_create(ns->ed.sdED, pg);
        // Synced mapper
        psd_synced_mapper(ns->ed.psdED);
    }
}

// Create the partitioned simulation sub-domain for multiphase simulation
void higflow_create_partitioned_domain_multiphase (higflow_solver *ns, partition_graph *pg, int order) {
    if (ns->contr.flowtype == 2) {
        // Creating the partitioned sub-domain to simulation
        ns->ed.psdED = psd_create(ns->ed.sdED, pg);
        // Synced mapper
        psd_synced_mapper(ns->ed.psdED);
    }
}

// Create the partitioned simulation sub-domain for viscoelastic simulation
void higflow_create_partitioned_domain_viscoelastic (higflow_solver *ns, partition_graph *pg, int order) {
    if (ns->contr.flowtype == 3) {
        // Creating the partitioned sub-domain to simulation
        ns->ed.psdED = psd_create(ns->ed.sdED, pg);
        // Synced mapper
        psd_synced_mapper(ns->ed.psdED);
    }
}

// Create the partitioned simulation sub-domain for viscoelastic simulation integral model
void higflow_create_partitioned_domain_viscoelastic_integral (higflow_solver *ns, partition_graph *pg, int order) {
    if (ns->contr.flowtype == 4) {
        // Creating the partitioned sub-domain to simulation
        ns->ed.psdED = psd_create(ns->ed.sdED, pg);
        // Synced mapper
        psd_synced_mapper(ns->ed.psdED);
    }
}

// Create the partitioned simulation sub-domain for electroosmotic simulation
void higflow_create_partitioned_domain_electroosmotic (higflow_solver *ns, partition_graph *pg, int order) {
    if (ns->contr.modelflowtype == 1) {
        // Creating the partitioned sub-domain to simulation EO phi
        ns->ed.eo.psdEOphi = psd_create(ns->ed.eo.sdEOphi, pg);
        // Synced mapper
        psd_synced_mapper(ns->ed.eo.psdEOphi);
        // Creating the partitioned sub-domain to simulation EO psi
        ns->ed.eo.psdEOpsi = psd_create(ns->ed.eo.sdEOpsi, pg);
        // Synced mapper
        psd_synced_mapper(ns->ed.eo.psdEOpsi);
        // Creating the partitioned sub-domain to simulation EO nplus
        ns->ed.eo.psdEOnplus = psd_create(ns->ed.eo.sdEOnplus, pg);
        // Synced mapper
        psd_synced_mapper(ns->ed.eo.psdEOnplus);
        // Creating the partitioned sub-domain to simulation EO nminus
        ns->ed.eo.psdEOnminus = psd_create(ns->ed.eo.sdEOnminus, pg);
        // Synced mapper
        psd_synced_mapper(ns->ed.eo.psdEOnminus);
        // Linking the property to the facets 
        for (int dim = 0; dim < DIM; dim++) {
            // Creating the list of cells center in the local domain (source force)
            ns->ed.eo.sfdEOFeo[dim] = sfd_create(NULL, dim);
            // Setting the order for the properties interpolation (source force)
            sfd_set_interpolator_order(ns->ed.eo.sfdEOFeo[dim], order);   
            // Copying the list of cells center in the domain for pressure 
            sfd_copy_higtrees_from_center_domain(ns->ed.eo.sfdEOFeo[dim], ns->ed.eo.sdEOpsi);
            sfd_adjust_facet_ids(ns->ed.eo.sfdEOFeo[dim]);
            // Creating property for the facets (source force)
            ns->ed.eo.psfdEOFeo[dim] = psfd_create(ns->ed.eo.sfdEOFeo[dim], ns->ed.eo.psdEOpsi);
        }
        for(int dim = 0; dim < DIM; dim++) {
            // Mapping the properties in the domain (source force)
            psfd_compute_sfbi(ns->ed.eo.psfdEOFeo[dim]);
            // Sync mapper for source force
            psfd_synced_mapper(ns->ed.eo.psfdEOFeo[dim]);
        }
    }
}

// Create the distributed properties for NS object
void higflow_create_ditributed_properties(higflow_solver *ns) {
    // Distributed property for pressure 
    ns->dpp     = psd_create_property(ns->psdp);
    // Distributed property for pressure difference 
    ns->ddeltap = psd_create_property(ns->psdp);
    // Distributed property for source term 
    ns->dpF     = psd_create_property(ns->psdF);
    // Distributed property for facets
    for(int dim = 0; dim < DIM; dim++) {
        // Final velocity 
        ns->dpu[dim]     = psfd_create_property(ns->psfdu[dim]);
        // Intermediate velocity 
        ns->dpustar[dim] = psfd_create_property(ns->psfdu[dim]);
        // Auxiliary intermediate velocity for Runge-Kutta method
        ns->dpuaux[dim]  = psfd_create_property(ns->psfdu[dim]);
        // Facet source term
        ns->dpFU[dim]    = psfd_create_property(ns->psfdF[dim]);
    }
}

// Create the distributed properties for generalized newtonian simulation
void higflow_create_ditributed_properties_generalized_newtonian(higflow_solver *ns) {
    // Non Newtonian tensor
    if (ns->contr.flowtype == 1) {
         // Distributed property for viscosity 
         ns->ed.gn.dpvisc  = psd_create_property(ns->ed.psdED);
         // Distributed property for non newtonian tensor 
         for (int i = 0; i < DIM; i++) {
              for (int j = 0; j < DIM; j++) {
                  ns->ed.gn.dpD[i][j] = psd_create_property(ns->ed.psdED);
              }
         }
    }
}

// Create the distributed properties for multiphase simulation
void higflow_create_ditributed_properties_multiphase(higflow_solver *ns) {
	// Non Newtonian tensor
	if (ns->contr.flowtype == 2) {
		// Distributed property for viscosity
		ns->ed.mult.dpvisc = psd_create_property(ns->ed.psdED);
		// Distributed property for density
		ns->ed.mult.dpdens = psd_create_property(ns->ed.psdED);
		// Distributed property for volume fraction
		ns->ed.mult.dpfracvol = psd_create_property(ns->ed.psdED);
		// Distributed property for auxiliary volume fraction
		ns->ed.mult.dpfracvolaux = psd_create_property(ns->ed.psdED);
		// Distributed property for curvature
		ns->ed.mult.dpcurvature = psd_create_property(ns->ed.psdED);
		// Distributed property for curvature
		ns->ed.mult.dpdistance = psd_create_property(ns->ed.psdED);
		// Distributed property for interfacial force
		for (int i = 0; i < DIM; i++) {
			ns->ed.mult.dpIF[i]     = psd_create_property(ns->ed.psdED);
			ns->ed.mult.dpnormal[i] = psd_create_property(ns->ed.psdED);
		}
		// Distributed property for beta
		ns->ed.mult.dpbeta = psd_create_property(ns->ed.psdED);
		
		for (int i = 0; i < DIM; i++) {
			for (int j = 0; j < DIM; j++) {
				ns->ed.mult.dpS0[i][j]      = psd_create_property(ns->ed.psdED);
				ns->ed.mult.dpS1[i][j]      = psd_create_property(ns->ed.psdED);
            ns->ed.mult.dpKernel0[i][j] = psd_create_property(ns->ed.psdED);
            ns->ed.mult.dpKernel1[i][j] = psd_create_property(ns->ed.psdED);
			}
		}
		
	}
}

// Create the distributed properties for viscoelastic simulation
void higflow_create_ditributed_properties_viscoelastic(higflow_solver *ns) {
    // Non Newtonian tensor
    if ((ns->contr.flowtype == 3) || (ns->contr.flowtype == 2)) {
         // Distributed property for viscoelastic tensor 
         for (int i = 0; i < DIM; i++) {
              for (int j = 0; j < DIM; j++) {
                  ns->ed.ve.dpD[i][j]      = psd_create_property(ns->ed.psdED);
                  ns->ed.ve.dpS[i][j]      = psd_create_property(ns->ed.psdED);
                  ns->ed.ve.dpKernel[i][j] = psd_create_property(ns->ed.psdED);
              }
         }
    }
}


// Create the distributed properties for viscoelastic simulation integral model
void higflow_create_ditributed_properties_viscoelastic_integral(higflow_solver *ns) {
    // Non Newtonian tensor
    if (ns->contr.flowtype == 4) {
         // Distributed property for viscoelastic tensor 
         for (int i = 0; i < DIM; i++) {
              for (int j = 0; j < DIM; j++) {
                  ns->ed.im.dpD[i][j] = psd_create_property(ns->ed.psdED);
                  ns->ed.im.dpS[i][j] = psd_create_property(ns->ed.psdED);
                  for (int k = 0; k <= NDT; k++)
                      ns->ed.im.dpB[k][i][j] = psd_create_property(ns->ed.psdED);
              }
         }
    }
}

// Create the distributed properties for electroosmotic simulation
void higflow_create_ditributed_properties_electroosmotic(higflow_solver *ns) {
    if (ns->contr.modelflowtype == 1) {
        for(int dim = 0; dim < DIM; dim++) {
            ns->ed.eo.dpFeo[dim]      = psfd_create_property(ns->ed.eo.psfdEOFeo[dim]);
        }
        ns->ed.eo.dpphi      = psd_create_property(ns->ed.eo.psdEOphi);
        ns->ed.eo.dppsi      = psd_create_property(ns->ed.eo.psdEOpsi);
        ns->ed.eo.dpnplus    = psd_create_property(ns->ed.eo.psdEOnplus);
        ns->ed.eo.dpnminus   = psd_create_property(ns->ed.eo.psdEOnminus);
    }
}

// Create the stencil for the NS object
void higflow_create_stencil(higflow_solver *ns) {
    // Create the stencil for properties interpolation (pressure and velocities)
    ns->stn   = stn_create();
    // Create the stencil for properties interpolation (source force)
    ns->stnF  = stn_create();
}

// Create the stencil for the extra domain
void higflow_create_stencil_for_extra_domain(higflow_solver *ns) {
    // Create the stencil for properties interpolation 
    ns->ed.stn = stn_create();
}

// Partition table initialize
void higflow_partition_domain (higflow_solver *ns, partition_graph *pg, int numhigs, higio_amr_info *mi[numhigs], int ntasks, int myrank) {
    // Setting the fringe size of the sub-domain
    // The fringe is a buffer around the cells of a given node
    pg_set_fringe_size(pg, 5);
    /* Partitioning the grid from AMR information */
    load_balancer *lb = lb_create(MPI_COMM_WORLD, 1);
    for(unsigned i = myrank; i < numhigs; i += ntasks) {
        lb_add_input_tree(lb, higio_read_from_amr_info(mi[i]), true, 0);
    }
    lb_calc_partition(lb, pg);
    numhigs = lb_get_num_local_trees(lb);
    for(int h = 0; h < numhigs; h++) {
        /* Creating the distributed HigTree data structure */
        hig_cell *root = lb_get_local_tree(lb, h, NULL);
        // Add higtree for SD
        sd_add_higtree(ns->sdp, root);
        if (ns->contr.flowtype > 0) {
            sd_add_higtree(ns->ed.sdED, root);
        }
        if (ns->contr.modelflowtype == 1) {
            sd_add_higtree(ns->ed.eo.sdEOphi, root);
            sd_add_higtree(ns->ed.eo.sdEOpsi, root);
            sd_add_higtree(ns->ed.eo.sdEOnplus, root);
            sd_add_higtree(ns->ed.eo.sdEOnminus, root);
        }
    }
    lb_destroy(lb);
}

// Set the external functions in the NS solver
void higflow_set_external_functions(higflow_solver *ns,
real (*get_pressure)(Point center, real t),
real (*get_velocity)(Point center, int dim, real t),
real (*get_source_term)(Point center, real t),
real (*get_facet_source_term)(Point center, int dim, real t),
real (*get_boundary_pressure)(int id, Point center, real t),
real (*get_boundary_velocity)(int id, Point center, int dim, real t),
real (*get_boundary_source_term)(int id, Point center, real t),
real (*get_boundary_facet_source_term)(int id, Point center, int dim, real t)) {
    // functions for the domain
    ns->func.get_pressure = get_pressure;
    ns->func.get_velocity = get_velocity;
    ns->func.get_source_term = get_source_term;
    ns->func.get_facet_source_term = get_facet_source_term;
    // functions for the boundary
    ns->func.get_boundary_pressure = get_boundary_pressure;
    ns->func.get_boundary_velocity = get_boundary_velocity;
    ns->func.get_boundary_source_term = get_boundary_source_term;
    ns->func.get_boundary_facet_source_term = get_boundary_facet_source_term;
}

// Auxiliar
void __higflow_readstring(char s[], int max, FILE *file) {
    int  i;
    for (i = 0; i < (max - 1); i++) {
       char letra = fgetc(file);
       if ((letra == '\n') && (i == 0)) {
          printf("=+=+=+ Erro na leitura do string =+=+=+\n");
          exit(1);
       }
       if (letra == '\n') break;
       s[i] = letra;
    }
    s[i] = 0;
    return;
}

