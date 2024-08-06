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
        case GENERALIZED_NEWTONIAN:
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
        case MULTIPHASE:
            // Multiphase
            dp_destroy(ns->ed.mult.dpvisc);
            dp_destroy(ns->ed.mult.dpdens);
            dp_destroy(ns->ed.mult.dpfracvol);
            dp_destroy(ns->ed.mult.dpfracvolaux);
            dp_destroy(ns->ed.mult.dpcurvature);
            dp_destroy(ns->ed.mult.dpdistance);
            // dp_destroy(ns->ed.mult.dpbeta);
            // Interfacial Term and Normal
            for (int i = 0; i < DIM; i++) {
                dp_destroy(ns->ed.mult.dpIF[i]);
                dp_destroy(ns->ed.mult.dpnormal[i]);
                // // Tensor terms destroy
                // for (int j = 0; j < DIM; j++) {
                //     dp_destroy(ns->ed.mult.dpD[i][j]);
                //     dp_destroy(ns->ed.mult.dpS[i][j]);
                //     // dp_destroy(ns->ed.mult.dpS0[i][j]);
                //     // dp_destroy(ns->ed.mult.dpS1[i][j]);
                //     // dp_destroy(ns->ed.mult.dpKernel0[i][j]);
                //     // dp_destroy(ns->ed.mult.dpKernel1[i][j]);
                // }
            }
            if(ns->ed.mult.contr.viscoelastic_either == true) {
                // Viscoelastic
                for (int i = 0; i < DIM; i++) {
                    for (int j = 0; j < DIM; j++) {
                        // Tensor terms destroy
                        dp_destroy(ns->ed.ve.dpD[i][j]);
                        dp_destroy(ns->ed.ve.dpS[i][j]);
                        dp_destroy(ns->ed.ve.dpKernel[i][j]);
                    }
                }
                // Destroy the stencil for extra domains
                stn_destroy(ns->ed.stn);
            }
            // Destroy the stencil for multiphase
            stn_destroy(ns->ed.mult.stn);
            break;
        case VISCOELASTIC:
            // Viscoelastic
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Tensor terms destroy
                    dp_destroy(ns->ed.ve.dpD[i][j]);
                    // dp_destroy(ns->ed.ve.dpD_prev[i][j]);
                    dp_destroy(ns->ed.ve.dpS[i][j]);
                    dp_destroy(ns->ed.ve.dpKernel[i][j]);
                }
            }
            // Destroy the stencil for extra domains
            stn_destroy(ns->ed.stn);
            break;
         case VISCOELASTIC_INTEGRAL:
           // Viscoelastic integral model
           for (int i = 0; i < DIM; i++) {
               for (int j = 0; j < DIM; j++) {
                   // Tensor terms destroy
                   dp_destroy(ns->ed.im.dpD[i][j]);
                   dp_destroy(ns->ed.im.dpS[i][j]);
                   for (int k = 0; k <= NDT; k++) {
                       dp_destroy(ns->ed.im.dpB[k][i][j]);
                   }
               }
           }
           // Destroy the stencil for extra domains
           stn_destroy(ns->ed.stn);
           break;   
         case VISCOELASTIC_VAR_VISCOSITY:
           //Destroy the variable viscosity
            dp_destroy(ns->ed.vevv.dpvisc);
            if (ns->ed.nn_contr.rheotype == THIXOTROPIC) {
                dp_destroy(ns->ed.vevv.dpStructPar);
            }
            //Destroy the viscoelastic tensor
            for (int i = 0; i < DIM; i++)
                for (int j = 0; j < DIM; j++) {
                    // Tensor terms destroy
                    dp_destroy(ns->ed.vevv.dpD[i][j]);
                    dp_destroy(ns->ed.vevv.dpS[i][j]);
                    dp_destroy(ns->ed.vevv.dpKernel[i][j]);
                }
            // Destroy the stencil for extra domains
            stn_destroy(ns->ed.stn);
            break;
        case SHEAR_BANDING:
            //VCM Model
            if (ns->ed.nn_contr.rheotype == VCM) {
                //Destroy nA and nB
                dp_destroy(ns->ed.vesb.dpnA);
                dp_destroy(ns->ed.vesb.dpnB);
                //Destroy cA and cB
                dp_destroy(ns->ed.vesb.dpcA);
                dp_destroy(ns->ed.vesb.dpcB);
                //Destroy the conformation tensors A and B
                for (int i = 0; i < DIM; i++)
                    for (int j = 0; j < DIM; j++) {
                        // Tensor terms destroy
                        dp_destroy(ns->ed.vesb.dpA[i][j]);
                        dp_destroy(ns->ed.vesb.dpB[i][j]);
                    }
                //Destroy the solvers of nA and nB    
                slv_destroy(ns->ed.vesb.slvnA);
                slv_destroy(ns->ed.vesb.slvnB);    
            }
            //Destroy the deformation and polymeric tensors
            for (int i = 0; i < DIM; i++)
                for (int j = 0; j < DIM; j++) {
                    // Tensor terms destroy
                    dp_destroy(ns->ed.vesb.dpD[i][j]);
                    dp_destroy(ns->ed.vesb.dpS[i][j]);
                }
            // Destroy the stencil for extra domains
            stn_destroy(ns->ed.stn);
            break;
        case ELASTOVISCOPLASTIC:
            // Elastoviscoplastic
            for (int i = 0; i < DIM; i++)
                for (int j = 0; j < DIM; j++) {
                    // Tensor terms destroy
                    dp_destroy(ns->ed.vepl.dpD[i][j]);
                    dp_destroy(ns->ed.vepl.dpS[i][j]);
                    dp_destroy(ns->ed.vepl.dpKernel[i][j]);
                }
            // Destroy the stencil for extra domains
            stn_destroy(ns->ed.stn);
            break;
        case SUSPENSIONS:
            // Shear-thickening suspensions
            //Destroy this for the model that considers the particle migration equation
            if (ns->ed.stsp.contr.model == GW_WC_IF) {
                dp_destroy(ns->ed.stsp.dpphi);
            }
            for (int i = 0; i < DIM; i++)
                for (int j = 0; j < DIM; j++) {
                    // Tensor terms destroy
                    dp_destroy(ns->ed.stsp.dpD[i][j]);
                    dp_destroy(ns->ed.stsp.dpS[i][j]);
                    dp_destroy(ns->ed.stsp.dpA[i][j]);
                }
            // Destroy the stencil for extra domains
            stn_destroy(ns->ed.stn);
            break;
    }
    if (ns->contr.eoflow == true || (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true)) {
        // Destroy the distributed properties for Electro-osmotic model
        for(int dim = 0; dim < DIM; dim++) {
            dp_destroy(ns->ed.eo.dpFeo[dim]);
        }
        dp_destroy(ns->ed.eo.dpphi);
        dp_destroy(ns->ed.eo.dppsi);
        dp_destroy(ns->ed.eo.dpnplus);
        dp_destroy(ns->ed.eo.dpnminus);
        dp_destroy(ns->ed.eo.dpnplus_temp);
        dp_destroy(ns->ed.eo.dpnminus_temp);
        // Destroy the solver for potential psi 
        slv_destroy(ns->ed.eo.slvpsi);
        slv_destroy(ns->ed.eo.slvphi);
        slv_destroy(ns->ed.eo.slvnplus);
        slv_destroy(ns->ed.eo.slvnminus);
        // Destroy the stencils
        stn_destroy(ns->ed.eo.stnphi);
        stn_destroy(ns->ed.eo.stnpsi);
        stn_destroy(ns->ed.eo.stnnplus);
        stn_destroy(ns->ed.eo.stnnminus);
    }
    // Destroy the stencil for pressure and velocities
    stn_destroy(ns->stn);
    // Destroy the stencil for source force
    stn_destroy(ns->stnF);
    // Destroy the solver for pressure
    slv_destroy(ns->slvp);
    // Destroy the solver for the velocity
    if (ns->contr.tempdiscrtype == SEMI_IMPLICIT_EULER || ns->contr.tempdiscrtype == SEMI_IMPLICIT_CN || ns->contr.tempdiscrtype == SEMI_IMPLICIT_BDF2) {
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
    if (ns->contr.tempdiscrtype == SEMI_IMPLICIT_EULER || ns->contr.tempdiscrtype == SEMI_IMPLICIT_CN || ns->contr.tempdiscrtype == SEMI_IMPLICIT_BDF2) {
        for (int dim = 0; dim < DIM; dim++) {
            // Get the localdomainsize for facet center
            localdomainsize = psfd_get_local_domain_size(ns->psfdu[dim]);
            // Creates a solver for velocity
            ns->slvu[dim]   = slv_create(SOLVER_ANY, psfd_get_first_id(ns->psfdu[dim]), localdomainsize);
            // Set the maximum of non zeros 
            slv_set_maxnonzeros(ns->slvu[dim], 800);
        }
    }
    if (ns->contr.eoflow == true || (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true)) {

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
    if (ns->ed.nn_contr.rheotype == VCM) {
        // Get the localdomainsize for cell center
        int localdomainsizenA = psd_get_local_domain_size(ns->ed.vesb.psdSBnA);
        // Creates a solver for nA
        ns->ed.vesb.slvnA            = slv_create(SOLVER_ANY, psd_get_first_id(ns->ed.vesb.psdSBnA), localdomainsizenA);
        // Set the maximum of non zeros 
        slv_set_maxnonzeros(ns->ed.vesb.slvnA, 800);
        // Get the localdomainsize for cell center
        int localdomainsizenB = psd_get_local_domain_size(ns->ed.vesb.psdSBnB);
        // Creates a solver for nB
        ns->ed.vesb.slvnB            = slv_create(SOLVER_ANY, psd_get_first_id(ns->ed.vesb.psdSBnB), localdomainsizenB);
        // Set the maximum of non zeros 
        slv_set_maxnonzeros(ns->ed.vesb.slvnB, 800);
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
    if (ns->contr.tempdiscrtype == SEMI_IMPLICIT_EULER || ns->contr.tempdiscrtype == SEMI_IMPLICIT_CN || ns->contr.tempdiscrtype == SEMI_IMPLICIT_BDF2) {
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
    if (ns->contr.eoflow == true || (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true)) {
    
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
    if (ns->ed.nn_contr.rheotype == VCM) {
        // Get the localdomainsize for cell center
        int localdomainsizenA = psd_get_local_domain_size(ns->ed.vesb.psdSBnA);
        // Creates a solver for nA
        ns->ed.vesb.slvnA            = slv_create(SOLVER_ANY, psd_get_first_id(ns->ed.vesb.psdSBnA), localdomainsizenA);
        // Set the maximum of non zeros 
        slv_set_maxnonzeros(ns->ed.vesb.slvnA, 800);
        // Get the localdomainsize for cell center
        int localdomainsizenB = psd_get_local_domain_size(ns->ed.vesb.psdSBnB);
        // Creates a solver for nB
        ns->ed.vesb.slvnB            = slv_create(SOLVER_ANY, psd_get_first_id(ns->ed.vesb.psdSBnB), localdomainsizenB);
        // Set the maximum of non zeros 
        slv_set_maxnonzeros(ns->ed.vesb.slvnB, 800);
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
    if (ns->contr.flowtype == GENERALIZED_NEWTONIAN) {
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
real (*get_fracvol)(Point center, Point delta, real t)) {
    if (ns->contr.flowtype == MULTIPHASE) {
       // simulation domain (SD) extra domain
       ns->ed.mult.sdmult = sd_create(NULL);
       //reuse interpolation, 0 on, 1 off
       sd_use_cache(ns->ed.mult.sdmult, cache);      
       //Sets the order of the interpolation to bhe used for the SD. 
       sd_set_interpolator_order(ns->ed.mult.sdmult, order);

       if(ns->ed.mult.contr.viscoelastic_either == true){
            // simulation domain (SD) extra domain
            ns->ed.sdED = sd_create(NULL);
            //reuse interpolation, 0 on, 1 off
            sd_use_cache(ns->ed.sdED, cache);      
            //Sets the order of the interpolation to bhe used for the SD. 
            sd_set_interpolator_order(ns->ed.sdED, order);
       }

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

       ns->ed.mult.num_plic_lines = 0;
    }
}

void higflow_create_domain_multiphase_viscoelastic (higflow_solver *ns,
real (*get_tensor_multiphase)(real fracvol, Point center, int i, int j, real t),
real (*get_kernel)(int dim, real lambda, real tol),
real (*get_kernel_inverse)(int dim, real lambda, real tol),
real (*get_kernel_jacobian)(int dim, real lambda, real tol)) {
    if (ns->contr.flowtype == MULTIPHASE) {
        if(ns->ed.mult.contr.viscoelastic_either == true) {
            // function for the domain
            ns->ed.mult.ve.get_tensor_multiphase          = get_tensor_multiphase;
            // function for the kernel transformation
            ns->ed.mult.ve.get_kernel          = get_kernel;
            // function for the inverse kernel transformation
            ns->ed.mult.ve.get_kernel_inverse  = get_kernel_inverse;
            // function for the kernel transformation jacobian
            ns->ed.mult.ve.get_kernel_jacobian = get_kernel_jacobian;
        }
    }
}

// Define the user function for viscoelastic flow
void higflow_define_user_function_multiphase_viscoelastic(higflow_solver *ns, 
void (*calculate_m_user_multiphase)(real fracvol, real lambda[DIM], real jlambda[DIM],real M_aux[DIM][DIM], real Re, real trS, ve_parameters *par0, ve_parameters *par1)) {
    if (ns->contr.flowtype == MULTIPHASE) {
       // function for the user viscoelastic model
       ns->ed.mult.ve.calculate_m_user_multiphase = calculate_m_user_multiphase;
    }
}

// Create the simulation domain for viscoelastic flow
void higflow_create_domain_viscoelastic(higflow_solver *ns, int cache, int order,
real (*get_tensor)(Point center, int i, int j, real t),
real (*get_kernel)(int dim, real lambda, real tol),
real (*get_kernel_inverse)(int dim, real lambda, real tol),
real (*get_kernel_jacobian)(int dim, real lambda, real tol)) {
    if (ns->contr.flowtype == VISCOELASTIC) {
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
void (*calculate_m_user)(real lambda[DIM], real jlambda[DIM],real M_aux[DIM][DIM], real Re, real trS, ve_parameters *par)) {
    if (ns->contr.flowtype == VISCOELASTIC) {
       // function for the user viscoelastic model
       ns->ed.ve.calculate_m_user = calculate_m_user;
    }
}

// Create the simulation domain for viscoelastic flow integral model
void higflow_create_domain_viscoelastic_integral(higflow_solver *ns, int cache, int order,
real (*get_tensor)(Point center, int i, int j, real t)) {
   // simulation domain (SD) extra domain
   if (ns->contr.flowtype == VISCOELASTIC_INTEGRAL) {
       ns->ed.sdED = sd_create(NULL);
       //reuse interpolation, 0 on, 1 off
       sd_use_cache(ns->ed.sdED, cache);      
       //Sets the order of the interpolation to bhe used for the SD. 
       sd_set_interpolator_order(ns->ed.sdED, order);
       // function for the domain
       ns->ed.im.get_tensor          = get_tensor;
    }
}

// Create the simulation domain for multiphase electroosmotic flow
void higflow_create_domain_multiphase_electroosmotic(higflow_solver* ns, int cache, int order,
    real(*get_multiphase_electroosmotic_source_term)(real fracvol, Point center, int dim, real t),
    real(*get_multiphase_electroosmotic_phi)(real fracvol, Point center, real t),
    real(*get_multiphase_electroosmotic_psi)(real fracvol, Point center, real t),
    real(*get_multiphase_electroosmotic_nplus)(real fracvol, Point center, real t),
    real(*get_multiphase_electroosmotic_nminus)(real fracvol, Point center, real t),
    real(*get_boundary_multiphase_electroosmotic_source_term)(real fracvol, int id, Point center, int dim, real t),
    real(*get_boundary_multiphase_electroosmotic_phi)(real fracvol, int id, Point center, real t),
    real(*get_boundary_multiphase_electroosmotic_psi)(real fracvol, int id, Point center, real t),
    real(*get_boundary_multiphase_electroosmotic_nplus)(real fracvol, int id, Point center, real t),
    real(*get_boundary_multiphase_electroosmotic_nminus)(real fracvol, int id, Point center, real t),
    real(*get_permittivity)(real fracvol, Point center, real t)) {
    if (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true) {
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
        ns->ed.mult.eo.get_multiphase_electroosmotic_source_term = get_multiphase_electroosmotic_source_term;
        // function for the eo phi
        ns->ed.mult.eo.get_multiphase_electroosmotic_phi = get_multiphase_electroosmotic_phi;
        // function for the eo psi
        ns->ed.mult.eo.get_multiphase_electroosmotic_psi = get_multiphase_electroosmotic_psi;
        // function for the eo positive ion concentration 
        ns->ed.mult.eo.get_multiphase_electroosmotic_nplus = get_multiphase_electroosmotic_nplus;
        // function for the eo negative ion concentration 
        ns->ed.mult.eo.get_multiphase_electroosmotic_nminus = get_multiphase_electroosmotic_nminus;
        // function for the boundary source term
        ns->ed.mult.eo.get_boundary_multiphase_electroosmotic_source_term = get_boundary_multiphase_electroosmotic_source_term;
        // function for the boundary potential phi
        ns->ed.mult.eo.get_boundary_multiphase_electroosmotic_phi = get_boundary_multiphase_electroosmotic_phi;
        // function for the boundary potential psi
        ns->ed.mult.eo.get_boundary_multiphase_electroosmotic_psi = get_boundary_multiphase_electroosmotic_psi;
        // function for the boundary concentration of positive ion
        ns->ed.mult.eo.get_boundary_multiphase_electroosmotic_nplus = get_boundary_multiphase_electroosmotic_nplus;
        // function for the boundary concentration of negative ion
        ns->ed.mult.eo.get_boundary_multiphase_electroosmotic_nminus = get_boundary_multiphase_electroosmotic_nminus;
        // function for the permittivity
        ns->ed.mult.eo.get_permittivity = get_permittivity;
    }
}

// Create the simulation domain for electroosmotic flow
void higflow_create_domain_electroosmotic(higflow_solver* ns, int cache, int order,
    real(*get_electroosmotic_source_term)(Point center, int dim, real t),
    real(*get_electroosmotic_phi)(Point center, real t),
    real(*get_electroosmotic_psi)(Point center, real t),
    real(*get_electroosmotic_nplus)(Point center, real t),
    real(*get_electroosmotic_nminus)(Point center, real t),
    real(*get_boundary_electroosmotic_source_term)(int id, Point center, int dim, real t),
    real(*get_boundary_electroosmotic_phi)(int id, Point center, real t),
    real(*get_boundary_electroosmotic_psi)(int id, Point center, real t),
    real(*get_boundary_electroosmotic_nplus)(int id, Point center, real t),
    real(*get_boundary_electroosmotic_nminus)(int id, Point center, real t),
    real(*get_permittivity)(Point center, real t)) {
    if (ns->contr.eoflow == true) {
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
        ns->ed.eo.get_electroosmotic_phi = get_electroosmotic_phi;
        // function for the eo psi
        ns->ed.eo.get_electroosmotic_psi = get_electroosmotic_psi;
        // function for the eo positive ion concentration 
        ns->ed.eo.get_electroosmotic_nplus = get_electroosmotic_nplus;
        // function for the eo negative ion concentration 
        ns->ed.eo.get_electroosmotic_nminus = get_electroosmotic_nminus;
        // function for the boundary source term
        ns->ed.eo.get_boundary_electroosmotic_source_term = get_boundary_electroosmotic_source_term;
        // function for the boundary potential phi
        ns->ed.eo.get_boundary_electroosmotic_phi = get_boundary_electroosmotic_phi;
        // function for the boundary potential psi
        ns->ed.eo.get_boundary_electroosmotic_psi = get_boundary_electroosmotic_psi;
        // function for the boundary concentration of positive ion
        ns->ed.eo.get_boundary_electroosmotic_nplus = get_boundary_electroosmotic_nplus;
        // function for the boundary concentration of negative ion
        ns->ed.eo.get_boundary_electroosmotic_nminus = get_boundary_electroosmotic_nminus;
        // function for the permittivity
        ns->ed.eo.get_permittivity = get_permittivity;
    }
}

// Create the simulation domain for viscoelastic flow with variable viscosity
void higflow_create_domain_viscoelastic_variable_viscosity (higflow_solver *ns, int cache, int order,
real (*get_tensor)(Point center, int i, int j, real t),
real (*get_kernel)(int dim, real lambda, real tol),
real (*get_kernel_inverse)(int dim, real lambda, real tol),
real (*get_kernel_jacobian)(int dim, real lambda, real tol),
real (*get_viscosity)(Point center, real q, real t, real beta, real struct_par),
real (*get_structpar)(Point center, real q, real t, real beta, real Phi, real Lambda, real Gamma)) {
    if ((ns->contr.flowtype == 6)) {
       // simulation domain (SD) extra domain
       ns->ed.sdED = sd_create(NULL);
       //reuse interpolation, 0 on, 1 off
       sd_use_cache(ns->ed.sdED, cache);      
       //Sets the order of the interpolation to bhe used for the SD. 
       sd_set_interpolator_order(ns->ed.sdED, order);
       // function for the domain
       ns->ed.vevv.get_tensor          = get_tensor;
       // function for the kernel transformation
       ns->ed.vevv.get_kernel          = get_kernel;
       // function for the inverse kernel transformation
       ns->ed.vevv.get_kernel_inverse  = get_kernel_inverse;
       // function for the kernel transformation jacobian
       ns->ed.vevv.get_kernel_jacobian = get_kernel_jacobian;
       // simulation domain for viscosity
       ns->ed.vevv.sdVisc = sd_create(NULL);
       //reuse interpolation, 0 on, 1 off
       sd_use_cache(ns->ed.vevv.sdVisc, cache);      
       //Sets the order of the interpolation to bhe used for the SD. 
       sd_set_interpolator_order(ns->ed.vevv.sdVisc, order);
       // viscosity function for the domain
       ns->ed.vevv.get_viscosity = get_viscosity;
       if (ns->ed.nn_contr.rheotype == THIXOTROPIC) {
           // Structural parameter function for the domain (For the BMP model, structural parameter=fluidity)
           ns->ed.vevv.get_structpar = get_structpar;
       }
    }
}

// Create the simulation domain for viscoelastic flow with shear-banding
void higflow_create_domain_viscoelastic_shear_banding (higflow_solver *ns, int cache, int order,
    real (*get_tensor)(Point center, int i, int j, real t),
    real (*get_tensor_A)(Point center, int i, int j, real t),
    real (*get_tensor_B)(Point center, int i, int j, real t),
    real (*get_nA)(Point center, real t),
    real (*get_nB)(Point center, real t),
    real (*get_cA)(Point center, real t, real CAeq, real chi, real ANA),
    real (*get_cB)(Point center, real t, real CBeq, real chi, real ANA),
    real (*get_boundary_nA)(int id, Point center, real t),
    real (*get_boundary_nB)(int id, Point center, real t)) {
    if ((ns->contr.flowtype == 7)) {
       // simulation domain (SD) extra domain
       ns->ed.sdED = sd_create(NULL);
       //reuse interpolation, 0 on, 1 off
       sd_use_cache(ns->ed.sdED, cache);      
       //Sets the order of the interpolation to bhe used for the SD. 
       sd_set_interpolator_order(ns->ed.sdED, order);
       // function for the domain
       ns->ed.vesb.get_tensor          = get_tensor;
       if (ns->ed.nn_contr.rheotype == VCM) {
           //Get the function for tensor_A
           ns->ed.vesb.get_tensor_A = get_tensor_A;
           //Get the function for tensor_B
           ns->ed.vesb.get_tensor_B = get_tensor_B;
           //Get the function for nA
           ns->ed.vesb.get_nA = get_nA;
           //Get the function for nB
           ns->ed.vesb.get_nB = get_nB;
           //Get the function for cA
           ns->ed.vesb.get_cA = get_cA;
           //Get the function for cB
           ns->ed.vesb.get_cB = get_cB;
           // function for the boundary density number of specie A
           ns->ed.vesb.get_boundary_nA = get_boundary_nA;
           // function for the boundary density number of specie B
           ns->ed.vesb.get_boundary_nB = get_boundary_nB;
           // simulation domain (SB) for nA
           ns->ed.vesb.sdSBnA = sd_create(NULL);
           //reuse interpolation, 0 on, 1 off
           sd_use_cache(ns->ed.vesb.sdSBnA, cache);      
           //Sets the order of the interpolation to bhe used for the SD. 
           sd_set_interpolator_order(ns->ed.vesb.sdSBnA, order);
           // simulation domain (SB) for nB
           ns->ed.vesb.sdSBnB = sd_create(NULL);
           //reuse interpolation, 0 on, 1 off
           sd_use_cache(ns->ed.vesb.sdSBnB, cache);      
           //Sets the order of the interpolation to bhe used for the SD. 
           sd_set_interpolator_order(ns->ed.vesb.sdSBnB, order);
       }
    }
}


// Create the simulation domain for elastoviscoplastic flow
void higflow_create_domain_elastoviscoplastic (higflow_solver *ns, int cache, int order,
real (*get_tensor)(Point center, int i, int j, real t),
real (*get_kernel)(int dim, real lambda, real tol),
real (*get_kernel_inverse)(int dim, real lambda, real tol),
real (*get_kernel_jacobian)(int dim, real lambda, real tol)) {
    if ((ns->contr.flowtype == 8)) {
       // simulation domain (SD) extra domain
       ns->ed.sdED = sd_create(NULL);
       //reuse interpolation, 0 on, 1 off
       sd_use_cache(ns->ed.sdED, cache);      
       //Sets the order of the interpolation to bhe used for the SD. 
       sd_set_interpolator_order(ns->ed.sdED, order);
       // function for the domain
       ns->ed.vepl.get_tensor          = get_tensor;
       // function for the kernel transformation
       ns->ed.vepl.get_kernel          = get_kernel;
       // function for the inverse kernel transformation
       ns->ed.vepl.get_kernel_inverse  = get_kernel_inverse;
       // function for the kernel transformation jacobian
       ns->ed.vepl.get_kernel_jacobian = get_kernel_jacobian;
    }
}

// Define the user function for viscoelastic flow
void higflow_define_user_function_elastoviscoplastic (higflow_solver *ns, 
void (*calculate_m_user)(real Re, real De, real beta, real Bi, real zeta, real epsilon, real Np, real tr, real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real M_aux[DIM][DIM], real tol, real smallTD, real SD)) {
    if (ns->ed.vepl.contr.model == GENERAL_SARAMITO) {
       // function for the user viscoelastic model
       ns->ed.vepl.calculate_m_user = calculate_m_user;
    }
}


// Create the simulation domain for shear-thickening suspension flow
void higflow_create_domain_shear_thickening_suspension (higflow_solver *ns, int cache, int order,
real (*get_tensor)(Point center, int i, int j, real t),
real (*get_tensor_A)(Point center, int i, int j, real t),
real (*get_X)(Point center, real t, real X0, real chi, real chi_J),
real (*get_vol_frac)(Point center, real t),
real (*get_alpha)(Point center, real t, real alpha, real phi, real phircp)) {
    if ((ns->contr.flowtype == 9)) {
        // simulation domain (SD) extra domain
        ns->ed.sdED = sd_create(NULL);
        //reuse interpolation, 0 on, 1 off
        sd_use_cache(ns->ed.sdED, cache);      
        //Sets the order of the interpolation to bhe used for the SD. 
        sd_set_interpolator_order(ns->ed.sdED, order);
        // function for the domain
        ns->ed.stsp.get_tensor          = get_tensor;
        // function for the microstructure tensor
        ns->ed.stsp.get_tensor_A        = get_tensor_A;
        // function for the X-term direct interparticle forces
        ns->ed.stsp.get_X  = get_X;
        //printf("=+=+=+= WE ARE HERE AFTER creating the domain =+=+=+=\n");
        //Only for th model that considers particle migration
        //printf("=+=+=+= WE ARE HERE BEFORE creating the domain =+=+=+=\n");
        //printf("=+=+=+= WE ARE HERE AFTER creating the domain =+=+=+=\n");
        // simulation domain for viscosity
        ns->ed.stsp.sdphi = sd_create(NULL);
        //reuse interpolation, 0 on, 1 off
        sd_use_cache(ns->ed.stsp.sdphi, cache);      
        //Sets the order of the interpolation to bhe used for the SD. 
        sd_set_interpolator_order(ns->ed.stsp.sdphi, order);
        // volume fraction value for the domain
        ns->ed.stsp.get_vol_frac = get_vol_frac;
        // volume fraction value for the domain
        ns->ed.stsp.get_alpha    = get_alpha;
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
        sfd_copy_higtrees_from_center_domain(ns->sfdF[dim], ns->sdF);
        sfd_adjust_facet_ids(ns->sfdF[dim]);
        // Creating property for the facets (source force)
        ns->psfdF[dim] = psfd_create(ns->sfdF[dim], ns->psdF);
    }
    for(int dim = 0; dim < DIM; dim++) {
        // Mapping the properties in the domain (source force)
        psfd_compute_sfbi(ns->psfdu[dim]);
        // Sync mapper for source force
        psfd_synced_mapper(ns->psfdu[dim]);
        // Mapping the properties in the domain (source force)
        psfd_compute_sfbi(ns->psfdF[dim]);
        // Sync mapper for source force
        psfd_synced_mapper(ns->psfdF[dim]);
    }
}

// Create the partitioned simulation sub-domain for generalized newtonian simulation
void higflow_create_partitioned_domain_generalized_newtonian (higflow_solver *ns, partition_graph *pg, int order) {
    if (ns->contr.flowtype == GENERALIZED_NEWTONIAN) {
        // Creating the partitioned sub-domain to simulation
        ns->ed.psdED = psd_create(ns->ed.sdED, pg);
        // Synced mapper
        psd_synced_mapper(ns->ed.psdED);
    }
}

// Create the partitioned simulation sub-domain for multiphase simulation
void higflow_create_partitioned_domain_multiphase (higflow_solver *ns, partition_graph *pg, int order) {
    if (ns->contr.flowtype == MULTIPHASE) {
        // Creating the partitioned sub-domain to simulation
        ns->ed.mult.psdmult = psd_create(ns->ed.mult.sdmult, pg);
        // Synced mapper
        psd_synced_mapper(ns->ed.mult.psdmult);

        if(ns->ed.mult.contr.viscoelastic_either == true) {
            // Creating the partitioned sub-domain to simulation
            ns->ed.psdED = psd_create(ns->ed.sdED, pg);
            // Synced mapper
            psd_synced_mapper(ns->ed.psdED);
        }
    }
}

// Create the partitioned simulation sub-domain for viscoelastic simulation
void higflow_create_partitioned_domain_viscoelastic (higflow_solver *ns, partition_graph *pg, int order) {
    if (ns->contr.flowtype == VISCOELASTIC) {
        // Creating the partitioned sub-domain to simulation
        ns->ed.psdED = psd_create(ns->ed.sdED, pg);
        // Synced mapper
        psd_synced_mapper(ns->ed.psdED);
    }
}

// Create the partitioned simulation sub-domain for viscoelastic simulation integral model
void higflow_create_partitioned_domain_viscoelastic_integral (higflow_solver *ns, partition_graph *pg, int order) {
    if (ns->contr.flowtype == VISCOELASTIC_INTEGRAL) {
        // Creating the partitioned sub-domain to simulation
        ns->ed.psdED = psd_create(ns->ed.sdED, pg);
        // Synced mapper
        psd_synced_mapper(ns->ed.psdED);
    }
}

// Create the partitioned simulation sub-domain for electroosmotic simulation
void higflow_create_partitioned_domain_electroosmotic (higflow_solver *ns, partition_graph *pg, int order) {
    if (ns->contr.eoflow == true || (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true)) {
    
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

// Create the partitioned simulation sub-domain for viscoelastic flow with variable viscosity simulation
void higflow_create_partitioned_domain_viscoelastic_variable_viscosity (higflow_solver *ns, partition_graph *pg, int order) {
    if (ns->contr.flowtype == VISCOELASTIC_VAR_VISCOSITY) {
        // Creating the partitioned sub-domain to simulation (for viscoelastic tensors)
        ns->ed.psdED = psd_create(ns->ed.sdED, pg);
        // Synced mapper (for viscoelastic tensors)
        psd_synced_mapper(ns->ed.psdED);
        // Creating the partitioned sub-domain to simulation for viscosity (or for the fluidity, in case rheotype==THIXOTROPIC)
        ns->ed.vevv.psdVisc = psd_create(ns->ed.vevv.sdVisc, pg);
        // Synced mapper for viscosity
        psd_synced_mapper(ns->ed.vevv.psdVisc);
    }
}

// Create the partitioned simulation sub-domain for viscoelastic flow with shear-banding simulation
void higflow_create_partitioned_domain_viscoelastic_shear_banding (higflow_solver *ns, partition_graph *pg, int order) {
    if (ns->contr.flowtype == SHEAR_BANDING) {
        // Creating the partitioned sub-domain to simulation (for viscoelastic tensors)
        ns->ed.psdED = psd_create(ns->ed.sdED, pg);
        // Synced mapper (for viscoelastic tensors)
        psd_synced_mapper(ns->ed.psdED);
        if (ns->ed.nn_contr.rheotype == VCM) {
            // Creating the partitioned sub-domain to simulation SB nA
            ns->ed.vesb.psdSBnA = psd_create(ns->ed.vesb.sdSBnA, pg);
            // Synced mapper
            psd_synced_mapper(ns->ed.vesb.psdSBnA);
            // Creating the partitioned sub-domain to simulation SB nB
            ns->ed.vesb.psdSBnB = psd_create(ns->ed.vesb.sdSBnB, pg);
            // Synced mapper
            psd_synced_mapper(ns->ed.vesb.psdSBnB);
        }
    }
}

// Create the partitioned simulation sub-domain for elastoviscoplastic simulation
void higflow_create_partitioned_domain_elastoviscoplastic (higflow_solver *ns, partition_graph *pg, int order) {
    if (ns->contr.flowtype == ELASTOVISCOPLASTIC) {
        // Creating the partitioned sub-domain to simulation
        ns->ed.psdED = psd_create(ns->ed.sdED, pg);
        // Synced mapper
        psd_synced_mapper(ns->ed.psdED);
    }
}

// Create the partitioned simulation sub-domain for shear-thickening suspension simulation
void higflow_create_partitioned_domain_shear_thickening_suspension (higflow_solver *ns, partition_graph *pg, int order) {
    if (ns->contr.flowtype == SUSPENSIONS) {
        // Creating the partitioned sub-domain to simulation
        ns->ed.psdED = psd_create(ns->ed.sdED, pg);
        // Synced mapper
        psd_synced_mapper(ns->ed.psdED);
        //only for the model that considers particle migration
        // Creating the partitioned sub-domain to simulation for volume fraction
        ns->ed.stsp.psdphi = psd_create(ns->ed.stsp.sdphi, pg);
        // Synced mapper for volume fraction
        psd_synced_mapper(ns->ed.stsp.psdphi);
    }
}

// Create the distributed properties for generalized newtonian simulation
void higflow_create_distributed_properties_generalized_newtonian(higflow_solver *ns) {
    // Non Newtonian tensor
    if (ns->contr.flowtype == GENERALIZED_NEWTONIAN) {
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
void higflow_create_distributed_properties_multiphase(higflow_solver *ns) {
    // Non Newtonian tensor
    if (ns->contr.flowtype == MULTIPHASE) {
      // Distributed property for viscosity
      ns->ed.mult.dpvisc = psd_create_property(ns->ed.mult.psdmult);
      // Distributed property for density
      ns->ed.mult.dpdens = psd_create_property(ns->ed.mult.psdmult);
      // Distributed property for volume fraction
      ns->ed.mult.dpfracvol = psd_create_property(ns->ed.mult.psdmult);
      // Distributed property for auxiliary volume fraction
      ns->ed.mult.dpfracvolaux = psd_create_property(ns->ed.mult.psdmult);
      // Distributed property for curvature
      ns->ed.mult.dpcurvature = psd_create_property(ns->ed.mult.psdmult);
      // Distributed property for curvature
      ns->ed.mult.dpdistance = psd_create_property(ns->ed.mult.psdmult);
      // // Distributed property for beta
      // ns->ed.mult.dpbeta = psd_create_property(ns->ed.mult.psdmult);
      // Distributed property for interfacial force and Normal
      for (int i = 0; i < DIM; i++) {
         ns->ed.mult.dpIF[i]     = psd_create_property(ns->ed.mult.psdmult);
            ns->ed.mult.dpnormal[i] = psd_create_property(ns->ed.mult.psdmult);
        }
        // for (int i = 0; i < DIM; i++) {
        //     for (int j = 0; j < DIM; j++) {
        //         ns->ed.mult.dpD[i][j]       = psd_create_property(ns->ed.psdED);
        //         ns->ed.mult.dpS[i][j]       = psd_create_property(ns->ed.psdED);
        //         // ns->ed.mult.dpS0[i][j]      = psd_create_property(ns->ed.psdED);
        //         // ns->ed.mult.dpS1[i][j]      = psd_create_property(ns->ed.psdED);
        //         // ns->ed.mult.dpKernel0[i][j] = psd_create_property(ns->ed.psdED);
        //         // ns->ed.mult.dpKernel1[i][j] = psd_create_property(ns->ed.psdED);
        //     }
        // }
    }
}

// Create the distributed properties for multiphase viscoelastic simulation
void higflow_create_distributed_properties_multiphase_viscoelastic(higflow_solver *ns) {
    if (ns->contr.flowtype == MULTIPHASE) {
        if(ns->ed.mult.contr.viscoelastic_either == true) {
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    ns->ed.ve.dpD[i][j]       = psd_create_property(ns->ed.psdED);
                    ns->ed.ve.dpS[i][j]       = psd_create_property(ns->ed.psdED);
                    ns->ed.ve.dpKernel[i][j] = psd_create_property(ns->ed.psdED);
                }
            }
        }
    }
}

// Create the distributed properties for viscoelastic simulation
void higflow_create_distributed_properties_viscoelastic(higflow_solver *ns) {
    // Non Newtonian tensor
    if (ns->contr.flowtype == VISCOELASTIC) {
         // Distributed property for viscoelastic tensor 
         for (int i = 0; i < DIM; i++) {
              for (int j = 0; j < DIM; j++) {
                  ns->ed.ve.dpD[i][j]      = psd_create_property(ns->ed.psdED);
                //   ns->ed.ve.dpD_prev[i][j] = psd_create_property(ns->ed.psdED);
                  ns->ed.ve.dpS[i][j]      = psd_create_property(ns->ed.psdED);
                  ns->ed.ve.dpKernel[i][j] = psd_create_property(ns->ed.psdED);
              }
         }
    }
}

// Create the distributed properties for viscoelastic simulation integral model
void higflow_create_distributed_properties_viscoelastic_integral(higflow_solver *ns) {
    // Non Newtonian tensor
    if (ns->contr.flowtype == VISCOELASTIC_INTEGRAL) {
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
void higflow_create_distributed_properties_electroosmotic(higflow_solver *ns) {
    if (ns->contr.eoflow == true || (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true)) {
        for(int dim = 0; dim < DIM; dim++) {
            ns->ed.eo.dpFeo[dim]      = psfd_create_property(ns->ed.eo.psfdEOFeo[dim]);
        }
        ns->ed.eo.dpphi      = psd_create_property(ns->ed.eo.psdEOphi);
        ns->ed.eo.dppsi      = psd_create_property(ns->ed.eo.psdEOpsi);
        ns->ed.eo.dpnplus    = psd_create_property(ns->ed.eo.psdEOnplus);
        ns->ed.eo.dpnminus   = psd_create_property(ns->ed.eo.psdEOnminus);
        ns->ed.eo.dpnplus_temp = psd_create_property(ns->ed.eo.psdEOnplus);
        ns->ed.eo.dpnminus_temp = psd_create_property(ns->ed.eo.psdEOnminus);
    }
}

// Create the distributed properties for simulation of viscoelastic flows with variable viscosity
void higflow_create_distributed_properties_viscoelastic_variable_viscosity(higflow_solver *ns) {
    // Viscoelastic flow with variable viscosity distributed properties
    if (ns->contr.flowtype == VISCOELASTIC_VAR_VISCOSITY) {
         // Distributed property for variable viscosity
         ns->ed.vevv.dpvisc  = psd_create_property(ns->ed.vevv.psdVisc);
         //Distributed property for structural parameter
         if (ns->ed.nn_contr.rheotype == THIXOTROPIC) {
             ns->ed.vevv.dpStructPar  = psd_create_property(ns->ed.vevv.psdVisc);
         }
         // Distributed property for viscoelastic tensor 
         for (int i = 0; i < DIM; i++) {
              for (int j = 0; j < DIM; j++) {
                  ns->ed.vevv.dpD[i][j]      = psd_create_property(ns->ed.psdED);
                  ns->ed.vevv.dpS[i][j]      = psd_create_property(ns->ed.psdED);
                  ns->ed.vevv.dpKernel[i][j] = psd_create_property(ns->ed.psdED);
              }
         }
    }
}

// Create the distributed properties for simulation of viscoelastic flows with shear-banding
void higflow_create_distributed_properties_viscoelastic_shear_banding(higflow_solver *ns) {
    // Viscoelastic flow with shear-banding distributed properties
    if (ns->contr.flowtype == SHEAR_BANDING) {
        for (int i = 0; i < DIM; i++) {
              for (int j = 0; j < DIM; j++) {
                  ns->ed.vesb.dpD[i][j]      = psd_create_property(ns->ed.psdED);
                  ns->ed.vesb.dpS[i][j]      = psd_create_property(ns->ed.psdED);
              }
        }
         //Distributed property for density numbers nA and nB
         if (ns->ed.nn_contr.rheotype == VCM) {
             ns->ed.vesb.dpnA  = psd_create_property(ns->ed.vesb.psdSBnA);
             ns->ed.vesb.dpnB  = psd_create_property(ns->ed.vesb.psdSBnB);
             ns->ed.vesb.dpcA  = psd_create_property(ns->ed.psdED);
             ns->ed.vesb.dpcB  = psd_create_property(ns->ed.psdED);
             // Distributed property for viscoelastic tensor 
             for (int i = 0; i < DIM; i++) {
                 for (int j = 0; j < DIM; j++) {
                     ns->ed.vesb.dpA[i][j]      = psd_create_property(ns->ed.psdED);
                     ns->ed.vesb.dpB[i][j]      = psd_create_property(ns->ed.psdED);
                     }
            }
         }
    }
}

// Create the distributed properties for elastoviscoplastic simulation
void higflow_create_distributed_properties_elastoviscoplastic(higflow_solver *ns) {
    // Non Newtonian tensor
    if ((ns->contr.flowtype == ELASTOVISCOPLASTIC)) {
         // Distributed property for viscoelastic tensor 
         for (int i = 0; i < DIM; i++) {
              for (int j = 0; j < DIM; j++) {
                  ns->ed.vepl.dpD[i][j]      = psd_create_property(ns->ed.psdED);
                  ns->ed.vepl.dpS[i][j]      = psd_create_property(ns->ed.psdED);
                  ns->ed.vepl.dpKernel[i][j] = psd_create_property(ns->ed.psdED);
              }
         }
    }
}

// Create the distributed properties for simulation of shear-thickening suspensions
void higflow_create_distributed_properties_shear_thickening_suspension(higflow_solver *ns) {
    // Shear-thickening suspension distributed properties
    if (ns->contr.flowtype == SUSPENSIONS) {
        //printf("=+=+=+= WE ARE HERE =+=+=+=\n");
        //Distributed property for volumen fraction
        if (ns->ed.stsp.contr.model == GW_WC_IF) {
            //printf("=+=+=+= WE ARE HERE BEFORE CREATING DP =+=+=+=\n");
            ns->ed.stsp.dpphi  = psd_create_property(ns->ed.stsp.psdphi);
            //printf("=+=+=+= WE ARE HERE AFTER CREATING DP =+=+=+=\n");
        }
         // Distributed property for tensors
         for (int i = 0; i < DIM; i++) {
              for (int j = 0; j < DIM; j++) {
                  ns->ed.stsp.dpD[i][j]      = psd_create_property(ns->ed.psdED);
                  ns->ed.stsp.dpS[i][j]      = psd_create_property(ns->ed.psdED);
                  ns->ed.stsp.dpA[i][j]      = psd_create_property(ns->ed.psdED);
              }
         }
    }
}

// Create the distributed properties for NS object
void higflow_create_distributed_properties(higflow_solver *ns) {
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

    switch (ns->contr.flowtype) {
        case GENERALIZED_NEWTONIAN:
            higflow_create_distributed_properties_generalized_newtonian(ns);
            break;
        case MULTIPHASE:
            higflow_create_distributed_properties_multiphase(ns);
            if(ns->ed.mult.contr.viscoelastic_either == true) 
                higflow_create_distributed_properties_multiphase_viscoelastic(ns);
            break;
        case VISCOELASTIC:
            higflow_create_distributed_properties_viscoelastic(ns);
            break;
        case VISCOELASTIC_INTEGRAL:
            higflow_create_distributed_properties_viscoelastic_integral(ns);
            break;
        case VISCOELASTIC_VAR_VISCOSITY:
            higflow_create_distributed_properties_viscoelastic_variable_viscosity(ns);
            break;
        case SHEAR_BANDING:
            higflow_create_distributed_properties_viscoelastic_shear_banding(ns);
            break;
        case ELASTOVISCOPLASTIC:
            higflow_create_distributed_properties_elastoviscoplastic(ns);
            break;
        case SUSPENSIONS:
            higflow_create_distributed_properties_shear_thickening_suspension(ns);
            break;
    }

    if (ns->contr.eoflow == true || (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true))
        higflow_create_distributed_properties_electroosmotic(ns);
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

// Create the stencil for multiphase
void higflow_create_stencil_multiphase(higflow_solver *ns) {
    // Create the stencil for properties interpolation 
    ns->ed.mult.stn = stn_create();
}

// Create the stencils for electroosmotic flow
void higflow_create_stencils_electroosmotic(higflow_solver *ns) {
    // Create the stencils for properties interpolation 
    if (ns->contr.eoflow == true || (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true)) {

        ns->ed.eo.stnphi = stn_create();
        ns->ed.eo.stnpsi = stn_create();
        ns->ed.eo.stnnplus = stn_create();
        ns->ed.eo.stnnminus = stn_create();
    }
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
        // Add higtree for SDs
        sd_add_higtree(ns->sdp, root);
        sd_add_higtree(ns->sdF, root);
        if (ns->contr.flowtype != NEWTONIAN && ns->contr.flowtype != MULTIPHASE) {
            sd_add_higtree(ns->ed.sdED, root);
        }
        if (ns->contr.flowtype == MULTIPHASE) {
            if(ns->ed.mult.contr.viscoelastic_either == true) {
                sd_add_higtree(ns->ed.sdED, root);
            }
        }
        if (ns->contr.eoflow == true || (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true)) {
            sd_add_higtree(ns->ed.eo.sdEOphi, root);
            sd_add_higtree(ns->ed.eo.sdEOpsi, root);
            sd_add_higtree(ns->ed.eo.sdEOnplus, root);
            sd_add_higtree(ns->ed.eo.sdEOnminus, root);
        }

        //Viscoelastic flow with variable viscosity
        if ((ns->ed.nn_contr.rheotype == PLM) || (ns->ed.nn_contr.rheotype == THIXOTROPIC)) {
            sd_add_higtree(ns->ed.vevv.sdVisc, root);
        }
        //Viscoelastic flow with shear-banding
        if ((ns->ed.nn_contr.rheotype == VCM)) {
            sd_add_higtree(ns->ed.vesb.sdSBnA, root);
            sd_add_higtree(ns->ed.vesb.sdSBnB, root);
        }
        //Only for Particle migration in shear thickening suspensions
        if ((ns->contr.flowtype == SUSPENSIONS)) {
            sd_add_higtree(ns->ed.stsp.sdphi, root);
        }
    }
    lb_destroy(lb);
}

// Partition table multiphase initialize
void higflow_partition_domain_multiphase (higflow_solver *ns, partition_graph *pg, int numhigs, higio_amr_info *mi[numhigs], higio_amr_info *mi_mult[numhigs], int ntasks, int myrank) {
    if(ns->contr.flowtype == MULTIPHASE) {
        // The fringe is a buffer around the cells of a given node
        pg_set_fringe_size(pg, 5);
        /* Partitioning the grid from AMR information */
        load_balancer *lb = lb_create(MPI_COMM_WORLD, 2);
        for(unsigned i = myrank; i < numhigs; i += ntasks) {
            lb_add_input_tree(lb, higio_read_from_amr_info(mi[i]), true, 0);
        }
        for(unsigned i = myrank; i < numhigs; i += ntasks) {
            lb_add_input_tree(lb, higio_read_from_amr_info(mi_mult[i]), true, 1);
        }
        lb_calc_partition(lb, pg);
        int numhigs0 = lb_get_num_local_trees_in_group(lb, 0);
        for(int h = 0; h < numhigs0; h++) {
            /* Creating the distributed HigTree data structure */
            hig_cell *root = lb_get_local_tree_in_group(lb, h, 0);
            // Add higtree for SDs
            sd_add_higtree(ns->sdp, root);
            sd_add_higtree(ns->sdF, root);
            if (ns->contr.flowtype != NEWTONIAN && ns->contr.flowtype != MULTIPHASE) {
                sd_add_higtree(ns->ed.sdED, root);
            }
            if (ns->contr.flowtype == MULTIPHASE) {
                if(ns->ed.mult.contr.viscoelastic_either == true) {
                    sd_add_higtree(ns->ed.sdED, root);
                }
            }
            if (ns->contr.eoflow == true || (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true)) {
                sd_add_higtree(ns->ed.eo.sdEOphi, root);
                sd_add_higtree(ns->ed.eo.sdEOpsi, root);
                sd_add_higtree(ns->ed.eo.sdEOnplus, root);
                sd_add_higtree(ns->ed.eo.sdEOnminus, root);
            }
        }
        int numhigs1 = lb_get_num_local_trees_in_group(lb, 1);
        for(int h = 0; h < numhigs1; h++) {
            /* Creating the distributed HigTree data structure */
            hig_cell *root = lb_get_local_tree_in_group(lb, h, 1);
            // Add higtree for SDs
            sd_add_higtree(ns->ed.mult.sdmult, root);
        }
        lb_destroy(lb);
    }
}

// Use the given mesh information to create to create amr info about the multiphase domain
// Corresponds to a uniform higtree with the cell sizes of the last given level
higio_amr_info *higflow_create_amr_info_mult(higio_amr_info *mi) {
    DECL_AND_ALLOC(higio_amr_info, mi_mult, 1);
    for(int dim = 0; dim < DIM; dim++) {
		mi_mult->l[dim] = mi->l[dim];
		mi_mult->h[dim] = mi->h[dim];
	}
    mi_mult->numlevels = 1;
    ALLOC(__higio_level_info, mi_mult->levels, mi_mult->numlevels);
    for(int dim = 0; dim < DIM; dim++) {
        // this uniform mesh has the same size as the last level of the original mesh
        mi_mult->levels[0].delta[dim] = mi->levels[mi->numlevels - 1].delta[dim];
	}
    mi_mult->levels[0].numpatches = 1;
    ALLOC(__higio_patch_info, mi_mult->levels[0].patches, mi_mult->levels[0].numpatches);
    for(int dim = 0; dim < DIM; dim++) {
        mi_mult->levels[0].patches[0].initcell[dim] = 1;
        // calculate the number of cells in each direction
        int nc = (int)round((mi_mult->h[dim] - mi_mult->l[dim])/mi_mult->levels[0].delta[dim]);
        mi_mult->levels[0].patches[0].patchsize[dim] = nc;
	}

    return mi_mult;
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

