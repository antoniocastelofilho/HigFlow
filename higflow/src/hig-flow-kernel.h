// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Kernel - version 10/11/2016
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_KERNEL
#define HIG_FLOW_KERNEL

#include "higtree.h"
#include "higtree-io.h"
#include "higtree-iterator.h"
#include "pdomain.h"
#include "domain.h"
#include "lbal.h"
#include "solver.h"

#define DEBUG
#include "Debug-c.h"
#include <stdio.h>
#include <stdlib.h>

// *******************************************************************
// HiG-Flow Solver Data Structure
// *******************************************************************

typedef enum convecdiscr_type{
    CENTRAL = 0,
    FIRST_ORDER = 1,
    SECOND_ORDER = 2
} convecdiscr_type;

// Computational cell data structure
typedef struct higflow_compcell{
    // Convective Method
    convecdiscr_type     convec_type;
    // Central velocity at facet
    real    ufacet;
    // Central velocity at cell
    real    ucell;
    // Velocity for convective term
    real    v[DIM];
    // Cubista velocity for convective term
    real    vc[DIM];
    // Source term
    real    F;
    // Electroosmotic source term
    real    Feo;
    // Density
    real    dens;
    // Presure derivative
    real    dpdx;
    // Central first order velocity derivative
    real    dudxc[DIM];
    // Forward first order velocity derivative
    real    dudxr[DIM];
    // Backward first order velocity derivative
    real    dudxl[DIM];
    // Forward first order velocity derivative
    real    dudxrr[DIM];
    // Backward first order velocity derivative
    real    dudxll[DIM];
    // Second order velocity derivative
    real    du2dx2[DIM];
    // Tensor derivative
    real    dSdx[DIM];
    // Cental ionic concentration
    real    ncell;
    // Ionic potential derivative
    real    dpsidx;
    // Applied potential derivative
    real    dphidx;
    // Ionic potential second order derivative
    real    d2psidx2;
    // Applied potential second order derivative
    real    d2phidx2;
    // First order ionic derivative
    real    dndx;
    // Second order ionic derivative
    real    d2ndx2;
    // Psi at the center
    real    psicell;
    // Viscosity at left cell
    real    viscl;
    // Viscosity at right cell
    real    viscr;
    //  Curvature cell
    real    curv;
    // Interfacial Force
    real    IF;
    // Central density number
    real nABcell;
    // First order density number derivative
    real dnABdx;
    // Second order density number derivative
    real d2nABdx2;
    // Cental volume fraction for particle suspensions
    real phivfcell;
    // Derivative of the volume fraction
    real dphivfdx;
    // Stress tensor of particle suspension derivative
    real dTdx[DIM];
    // Second derivative of the stress tensor of particle suspension
    real d2Tdx2;
} higflow_compcell;

// Parameters for simulation data structure
typedef struct higflow_parameters{
    // // Density
    // real   rho;
    // // Viscosity
    // real   nu;
    // Reynolds number
    real   Re;
    // Froude number
    real   Fr;
    // Time
    real   t;
    // Difference time
    real   dt;
    // Initial step number
    int    initstep;
    // Current step number
    int    step;
    int    stepaux;
    // number of steps
    int    numsteps;
    // Last step number
    int    finalstep;
    // Save time
    real   ts;
    // Difference save time
    real   dts;
    // Print time
    real   tp;
    // Difference print time
    real   dtp;
    // Frame number for print
    int    frame;
    // File name for print 
    char   *nameprint;
    // File name for load
    char   *nameload;
    // File name for save
    char   *namesave;
    // Folder name for residual
    char   nameres[1024];
} higflow_parameters;

typedef enum projection_type{
    NON_INCREMENTAL = 0,
    INCREMENTAL = 1
} projection_type;

typedef enum flow_type{
    NEWTONIAN = 0,
    GENERALIZED_NEWTONIAN = 1,
    MULTIPHASE = 2,
    VISCOELASTIC = 3,
    VISCOELASTIC_INTEGRAL = 4,
    VISCOELASTIC_VAR_VISCOSITY = 6,
    SHEAR_BANDING = 7,
    ELASTOVISCOPLASTIC = 8,
    SUSPENSIONS = 9,
} flow_type;

typedef enum tempdiscr_type{
    EXPLICIT_EULER = 0,
    EXPLICIT_RK2 = 1,
    EXPLICIT_RK3 = 2,
    SEMI_IMPLICIT_EULER = 3,
    SEMI_IMPLICIT_CN = 4,
    SEMI_IMPLICIT_BDF2 = 5
} tempdiscr_type;

typedef enum rheo_type{
    PLM = 2,
    THIXOTROPIC = 3,
    VCM = 4
} rheo_type;

typedef enum spatialdiscr_type{
    ORDER2 = 0,
    ORDER4 = 1
} spatialdiscr_type;

typedef enum secondconvecdiscr_type{
    MCU = 0,
    CUBISTA = 1,
    QUICK = 2
} secondconvecdiscr_type;

// Methods controllers for simulation data structure
typedef struct higflow_controllers{
    // Projection type: 0 non incremental, 1 incremental
    projection_type   projtype;
    // Model type: 0 Newtonian, 1 Generalized Newtonian, 2 Multifase, 3 Viscoelastic
    flow_type   flowtype;
    // Electroosmotic flow? : 0 No, 1 Electro-osmotic, 
    int   eoflow;
    // Rheological Flow Type : 2 Power Law Model, 3 ThixoTropic Model, 4 VCM
    rheo_type rheotype;
    // Temporal discretization type: 0 Explicit Euler, 1 Explicit RK2, .....
    tempdiscr_type   tempdiscrtype;
    // Spatial discretization type: 0 second order, 1 forth order
    spatialdiscr_type   spatialdiscrtype;
    // Convective discretization type: 0 central, 1 first order, 2 second order
    convecdiscr_type   convecdiscrtype;
    // Second order discretization type: 0 Modiffied Coefficient Upwind, 1 CUBISTA, 2 Quick
    secondconvecdiscr_type   secondconvecdiscrtype;
    // Desingualrization controllers for the pressure: 0 no desing. and 1 desing.
    int   desingpressure; 
} higflow_controllers;

// Methods controllers for simulation data structure
typedef struct higflow_functions{
    // Function to get the pressure
    real (*get_pressure)(Point center, real t);
    // Function to get the pressure at boundary
    real (*get_boundary_pressure)(int id, Point center, real t);
    // Function to get the velocity
    real (*get_velocity)(Point center, int dim, real t);
    // Function to get the velocity at boundary
    real (*get_boundary_velocity)(int id, Point center, int dim, real t);
    // Function to get the source term
    real (*get_source_term)(Point center, real t);
    // Function to get the source term at boundary
    real (*get_boundary_source_term)(int id, Point center, real t);
    // Function to get the facet source term
    real (*get_facet_source_term)(Point center, int dim, real t);
    // Function to get the facet source term at boundary
    real (*get_boundary_facet_source_term)(int id, Point center, int dim, real t);
} higflow_functions;

// Domains and distributed properties for generalized newtonian
typedef struct higflow_gen_newtonian{
    // Distributed property for velocity derivative tensor 
    distributed_property *dpD[DIM][DIM];
    // Distributed property for viscosity 
    distributed_property *dpvisc;
    // Function to get the viscosity
    real (*get_viscosity)(Point center, real q, real t);
} higflow_gen_newtonian;



// Parameters for multiphase viscoelastic simulation
// typedef struct mult_ve_parameters{
//     // Parameters - phase 0
//     // Deborah Number
//     real   De0;
//     // Ratio of solvent to total viscosity
//     real   beta0;
//     // LPTT parameter
//     real   epsilon0;
//     // LPTT parameter
//     real   psi0;
//     // Giesekus parameter
//     real   alpha0;
//     // Kernel tolerance parameter
//     real   kernel_tol0;
//     // GPTT parameters
//     real   alpha_gptt0;
//     real   beta_gptt0;
//     real   gamma_gptt0;
//     // FENE parameters
//     real   L2_fene0;
//     real   lambda_fene0;
//     real   E_fene0;
//     // Parameters - phase 1
//     // Deborah Number
//     real   De1;
//     // Ratio of solvent to total viscosity
//     real   beta1;
//     // LPTT parameter
//     real   epsilon1;
//     // LPTT parameter
//     real   psi1;
//     // Giesekus parameter
//     real   alpha1;
//     // Kernel tolerance parameter
//     real   kernel_tol1;
//     // GPTT parameters
//     real   alpha_gptt1;
//     real   beta_gptt1;
//     real   gamma_gptt1;
//     // FENE parameters
//     real   L2_fene1;
//     real   lambda_fene1;
//     real   E_fene1;
// } mult_ve_parameters;


typedef enum visc_model_type{
    USERSET = -1,
    OLDROYD_B = 0,
    GIESEKUS = 1,
    LPTT = 2,
    GPTT = 3,
    FENE_P = 4,
    E_FENE = 5
} visc_model_type;

typedef enum inhomogenous_discr_type{
    EXPLICIT = 0,
    IMPLICIT = 1
} inhomogenous_discr_type;

typedef enum cell_convecdiscr_type{
    CELL_CENTRAL = 0,
    CELL_CUBISTA = 1
} cell_convecdiscr_type;

// Parameters for electro-osmotic simulation
typedef struct eo_parameters{
    // alpha = e*z*zeta0/K*T 
    real   alpha;
    // delta = n0*e*z*H*H/(epsilon*zeta0)
    real   delta;
    // Peclet number = U*H/D 
    real   Pe;
    // External electric field
    real   Ex;
} eo_parameters;

typedef enum eo_model_type{
    PNP = 0,
    PB = 1,
    PBDH = 2,
    PBDH_ANALYTIC = 3
} eo_model_type;

// Controllers for viscoelastic simulation
typedef struct eo_controllers{
    // EO model Nernst-Planck, Poisson-Boltzmann or Debye-Hückel 
    eo_model_type    eo_model;
    tempdiscr_type    tempdiscrtype;             // temporal discretization (relative to diffusive term)
    cell_convecdiscr_type    convecdiscrtype;       // (u + grad (phi + psi)) dot grad n term
    // cell_convecdiscr_type    electric_convecdiscrtype; // grad (phi + psi) dot grad n    term
    // inhomogenous_discr_type    electricdiv_discrtype; // n lap (phi + psi) term
    int    is_phibc_timedependent; // if boundary conditions are time dependent - determines if laplace equation should be solved again
    int    is_psibc_timedependent; // if boundary conditions are time dependent - determines if poisson equation should be solved again
                                  // only relevant for poisson-boltzmann and poisson-boltzmann-debye huckel
} eo_controllers;

// Domains and distributed properties for electro-osmotic
typedef struct higflow_electroosmotic{
    // EO controllers
    eo_controllers        contr;
    // EO parameters
    eo_parameters         par;
    // Sub-domain to simulation for facets (EO source force)
    sim_facet_domain     *sfdEOFeo[DIM];
    // Partitioned sub-domain to simulation for facets (EO source force)
    psim_facet_domain    *psfdEOFeo[DIM];
    // Sub-domain to simulation for electro-osmotic phi
    sim_domain           *sdEOphi;
    // Partitioned sub-domain to simulation for electro-osmotic phi
    psim_domain          *psdEOphi;
    // Sub-domain to simulation for electro-osmotic psi
    sim_domain           *sdEOpsi;
    // Partitioned sub-domain to simulation for electro-osmotic psi
    psim_domain          *psdEOpsi;
    // Sub-domain to simulation for electrocosmotic nplus
    sim_domain           *sdEOnplus;
    // Partitioned sub-domain to simulation for electrocosmotic nminus
    psim_domain          *psdEOnplus;
    // Sub-domain to simulation for electro-osmotic nminus
    sim_domain           *sdEOnminus;
    // Partitioned sub-domain to simulation for electro-osmotic nminus
    psim_domain          *psdEOnminus;
    // Distributed property for source term 
    distributed_property *dpFeo[DIM];
    // Distributed property for applied potential 
    distributed_property *dpphi;
    // Distributed property for induced potential 
    distributed_property *dppsi;
    // Distributed property for positive charge concentration 
    distributed_property *dpnplus;
    // Distributed property for negative charge concentration 
    distributed_property *dpnminus;
    // stencil for phi
    sim_stencil *stnphi;
    // stencil for psi and Feo
    sim_stencil *stnpsi;
    // stencil for nplus
    sim_stencil *stnnplus;
    // stencil for nminus
    sim_stencil *stnnminus;
    // Function to get source term
    real (*get_electroosmotic_source_term)(Point center, int dim, real t);
    // Function to get phi
    real (*get_electroosmotic_phi)(Point center, real t);
    // Function to get psi
    real (*get_electroosmotic_psi)(Point center, real t);
    // Function to get nplus
    real (*get_electroosmotic_nplus)(Point center, real t);
    // Function to get nminus
    real (*get_electroosmotic_nminus)(Point center, real t);
    // Function to get the source term at boundary
    real (*get_boundary_electroosmotic_source_term)(int id, Point center, int dim, real t);
    // Function to get phi at boundary
    real (*get_boundary_electroosmotic_phi)(int id, Point center, real t);
    // Function to get psi at boundary
    real (*get_boundary_electroosmotic_psi)(int id, Point center, real t);
    // Function to get nplus at boundary
    real (*get_boundary_electroosmotic_nplus)(int id, Point center, real t);
    // Function to get nminus at boundary
    real (*get_boundary_electroosmotic_nminus)(int id, Point center, real t);
    // Linear system solver for potential psi
    solver                     *slvpsi;
    solver                     *slvphi;
    solver                     *slvnplus;
    solver                     *slvnminus;
} higflow_electroosmotic;

// Parameters for viscoelastic simulation
typedef struct ve_parameters{
    // Deborah Number
    real   De;
    // Ratio of solvent to total viscosity
    real   beta;
    // LPTT parameter
    real   epsilon;
    // LPTT parameter
    real   xi;
    // Giesekus parameter
    real   alpha;
    // Kernel tolerance parameter
    real   kernel_tol;
    // GPTT parameters
    real   alpha_gptt;
    real   beta_gptt;
    real   gamma_gptt;
    // FENE parameters
    real   L2_fene;
    real   lambda_fene;
    real   E_fene;
} ve_parameters;

// Controllers for viscoelastic simulation
typedef struct ve_controllers{
    // Viscoelastic model
    visc_model_type    model;
    // Viscoelastic discretization type
    inhomogenous_discr_type    discrtype;
    // Viscoelastic convective discretization type
    cell_convecdiscr_type    convecdiscrtype;
} ve_controllers;

// Domains and distributed properties for viscoelastic
typedef struct higflow_viscoelastic{
    // Viscoelastic controllers
    ve_controllers       contr;
    // Viscoelastic parameters
    ve_parameters        par;
    // Distributed property for velocity derivative tensor 
    distributed_property *dpD[DIM][DIM];
    // // Distributed property for velocity derivative tensor previous (used for e-FENE)
    // distributed_property *dpD_prev[DIM][DIM];
    // Distributed property for polymeric tensor 
    distributed_property *dpS[DIM][DIM];
    // Distributed property for kernel tensor 
    distributed_property *dpKernel[DIM][DIM];
    // Function to get the tensor
    real (*get_tensor)(Point center, int i, int j, real t);
    // Function to get the kernel transformation
    real (*get_kernel)(int dim, real lambda, real tol);
    // Function to get the inverse kernel transformation
    real (*get_kernel_inverse)(int dim, real lambda, real tol);
    // Function to get the kernel jacobian
    real (*get_kernel_jacobian)(int dim, real lambda, real tol);
    // User function to define the viscoelastic model
    void (*calculate_m_user)(real lambda[DIM], real jlambda[DIM],  real B[DIM][DIM], real M_aux[DIM][DIM], real Re, real trS, ve_parameters *par);
} higflow_viscoelastic;


// Define the number of relaxation time for integral viscoelastic model
#define  NDT 80

// Define the number of parameters for integral viscoelastic model
#define  NRP 8

// Parameters for viscoelastic simulation integral model
typedef struct im_parameters{
    // Deborah Number
    real   De;
    // Number of relaxation time
    int    M;
    // Parameters of Damping function
    real   alpha;
    real   beta;
    // Relaxation parameters for KBKZ model
    real   a[NRP];
    real   lambda[NRP];
    //  Integral points
    real   eps_intpoints;
    real   scorte; 
    // 
    real rho;
    real v_ref;
    real l_ref;
    // Relaxation parameters for KBKZ-Fractional model
    real   Phi1;
    real   Phi2;
    // Order of fractinal derivative for KBKZ-Fractional model
    real   alpha_frac;
    real   beta_frac;
} im_parameters;

typedef enum visc_integral_model_type {
    KBKZ = 0,
    KBKZ_FRACTIONAL = 1
} visc_integral_model_type;

typedef enum visc_integral_relax_model_type {
    PSM = 0,
    UCM = 1
} visc_integral_relax_model_type;

// Controllers for viscoelastic simulation integral models
typedef struct im_controllers{
    // Viscoelastic model
    visc_integral_model_type    model;
    // Relaxation model
    visc_integral_relax_model_type    model_H;
    // Viscoelastic discretization type
    inhomogenous_discr_type    discrtype;
    // Viscoelastic convective discretization type
    cell_convecdiscr_type    convecdiscrtype;
} im_controllers;

// Domains and distributed properties for viscoelastic integral models
typedef struct higflow_viscoelastic_integral{
    // Viscoelastic controllers
    im_controllers       contr;
    // Viscoelastic parameters
    im_parameters        par;
    // Discretization time
    real                 s[NDT+1];
    // Discretization time (old step)
    real                 sold[NDT+1];
    // Distributed property for velocity derivative tensor 
    distributed_property *dpB[NDT+1][DIM][DIM];
    // Distributed property for velocity derivative tensor 
    distributed_property *dpD[DIM][DIM];
    // Distributed property for polymeric tensor 
    distributed_property *dpS[DIM][DIM];
    // Function to get the tensor
    real (*get_tensor)(Point center, int i, int j, real t);
} higflow_viscoelastic_integral;

// Parameters for simulation of viscoelastic flows with variable viscosity
typedef struct vevv_parameters
{
    // Deborah Number
    real De;
    // Ratio of solvent to total viscosity
    real beta;
    // LPTT parameter
    real epsilon;
    // LPTT parameter
    real psi;
    // Giesekus parameter
    real alpha;
    // Kernel tolerance parameter
    real kernel_tol;
    // GPTT parameters
    real alpha_gptt;
    real beta_gptt;
    real gamma_gptt;
    //BMP model parameters
    real Lambda;
    real Phi;
    real Gamma;
} vevv_parameters;

// Controllers for viscoelastic simulation
typedef struct vevv_controllers
{
    // Viscoelastic model
    int model;
    // Viscoelastic discretization type
    int discrtype;
    // Viscoelastic convective discretization type
    int convecdiscrtype;
    // Discretization type for the structural parameter evolution equation (only for the models that have a kinetic equation for the structural parameter)
    int structpdiscrtype;
    // Convective discretization type for the structural parameter evolution equation (only for the models that have a kinetic equation for the structural parameter)
    int structpconvecdiscrtype;
    // Kinetic equation to be used for the structural parameter (Original BMP model, BMP model with solvent, MBM, NM_taup or NM_T)
    int structparmodel;
} vevv_controllers;

//Domains and distributed properties for viscoelastic flows with variable viscosity
typedef struct higflow_viscoelastic_variable_viscosity
{
    // Controllers for viscoelastic flows with variable viscosity
    vevv_controllers contr;
    // Parameters for viscoelastic flows with variable viscosity
    vevv_parameters par;
    //Sub-domain to simulation for viscosity
    sim_domain  *sdVisc;
    // Partitioned sub-domain to simulation for viscosity
    psim_domain *psdVisc;
    // Distributed property for variable viscosity
    distributed_property *dpvisc;
    // Function to get the viscosity
    real (*get_viscosity)(Point center, real q, real t, real beta, real struct_par);
    // Distributed property for variable structural parameter
    distributed_property *dpStructPar;
    // Function to get the structural parameter
    real (*get_structpar)(Point center, real q, real t, real beta, real Phi, real Lambda, real Gamma); 
    // Distributed property for velocity derivative tensor
    distributed_property *dpD[DIM][DIM];
    // Distributed property for polymeric tensor
    distributed_property *dpS[DIM][DIM];
    // Distributed property for kernel tensor
    distributed_property *dpKernel[DIM][DIM];
    // Function to get the tensor
    real (*get_tensor)(Point center, int i, int j, real t);
    // Function to get the kernel transformation
    real (*get_kernel)(int dim, real lambda, real tol);
    // Function to get the inverse kernel transformation
    real (*get_kernel_inverse)(int dim, real lambda, real tol);
    // Function to get the kernel jacobian
    real (*get_kernel_jacobian)(int dim, real lambda, real tol);
    // User function to define the viscoelastic model
    void (*calculate_m_user)(real Re, real De, real beta, real tr, real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real M_aux[DIM][DIM], real tol);
    // Function to get viscosity at boundary
    //real (*get_boundary_viscosity)(int id, Point center, real q, real t, real beta, real struct_par);
    // Function to get structura parameter at boundary
    //real (*get_boundary_structural_parameter)(int id, Point center, real q, real t, real beta);
    // Linear system solver for structural parameter (?)
    solver *slvstructpar;
} higflow_viscoelastic_variable_viscosity;

// Parameters for simulation of viscoelastic flows with shear-banding
typedef struct vesb_parameters
{
    // Effective Deborah Number De_eff= lambda_eff *(U0/L)
    real De;
    // Ratio of solvent to total viscosity
    real beta;
    // Deborah Number of specie A DeA = lambda_A*(U0/L)
    real DeA;
    // Epsilon = lambda_B/lambda_A
    real epsilon;
    // Peclet parameter of species A
    real PeA;
    // Peclet parameter of species B
    real PeB;
    // Chi parameter Chi = (1+CAeq*lambda_A)*Xi*(U0/L)
    real chi;
    // Equilibrium Concentration of species A
    real CAeq;
    // Equilibrium Concentration of species B
    real CBeq;
} vesb_parameters;

// Controllers for viscoelastic simulation
typedef struct vesb_controllers
{
    // Shear-banding model
    int model;
    // Discretization type for the viscoelastic equations of both species
    int discrtype;
    // Viscoelastic convective discretization type for both species
    int convecdiscrtype;
    // Discretization type for the evolution equations of the density numbers of both species 
    int nAnBdiscrtype;
    // Convective discretization type for the  evolution equations of the density numbers of both species
    int nAnBconvecdiscrtype;
} vesb_controllers;

//Domains and distributed properties for viscoelastic flows with variable viscosity
typedef struct higflow_viscoelastic_shear_banding
{
    // Controllers for viscoelastic flows with variable viscosity
    vesb_controllers contr;
    // Parameters for viscoelastic flows with variable viscosity
    vesb_parameters par;
    // Sub-domain to simulation for density number of specie A nA
    sim_domain *sdSBnA;
    // Partitioned sub-domain for density number of specie A nA
    psim_domain *psdSBnA;
    // Sub-domain to simulation for density number of specie B nB
    sim_domain *sdSBnB;
    // Partitioned sub-domain for density number of specie B nB
    psim_domain *psdSBnB;
    // Distributed property for density number of specie A nA
    distributed_property *dpnA;
    // Distributed property for density number of specie B nB
    distributed_property *dpnB;
    // Distributed property for concentration of specie A cA
    distributed_property *dpcA;
    // Distributed property for concentration of specie B cB
    distributed_property *dpcB;
    // Distributed property for velocity derivative tensor
    distributed_property *dpD[DIM][DIM];
    // Distributed property for polymeric tensor
    distributed_property *dpS[DIM][DIM];
    // Distributed property for conformation tensor of specie A
    distributed_property *dpA[DIM][DIM];
    // Distributed property for conformation tensor of specie B
    distributed_property *dpB[DIM][DIM];
    // Function to get the tensor
    real (*get_tensor)(Point center, int i, int j, real t);
    // Function to get the A tensor (conformation tensor of specie A)
    real (*get_tensor_A)(Point center, int i, int j, real t);
    // Function to get the B tensor (conformation tensor of specie B)
    real (*get_tensor_B)(Point center, int i, int j, real t);
    // Function to get nA
    real (*get_nA)(Point center, real t);
    // Function to get nB
    real (*get_nB)(Point center, real t);
    // Function to get cA
    real (*get_cA)(Point center, real t, real CAeq, real chi, real ANA);
    // Function to get cB
    real (*get_cB)(Point center, real t, real CBeq, real chi, real ANA);
    // Function to get nA at boundary
    real (*get_boundary_nA)(int id, Point center, real t);
    // Function to get nB at boundary
    real (*get_boundary_nB)(int id, Point center, real t);
    // Linear system solver for density numbers
    solver *slvnA;
    solver *slvnB;
} higflow_viscoelastic_shear_banding;


// Parameters for elastoviscoplastic simulation
typedef struct vepl_parameters
{
    // Deborah Number
    real De;
    // Ratio of solvent to total viscosity
    real beta;
    // Bingham number
    real Bi;
    // zeta parameter (from the Gordon-Schowalter derivative)
    real zeta;
    // Kernel tolerance parameter
    real kernel_tol;
    // LPTT and EPTT parameter
    real epsilon;
    // Power-law coefficient (for the Oldroyd-B-Herschel-Bulkley)
    real Np;
    // GPTT parameters
    //real alpha_gptt;
    //real beta_gptt;
    //real gamma_gptt;
} vepl_parameters;

// Controllers for elastoviscoplastic simulation
typedef struct vepl_controllers
{
    // Elastoviscoplastic model
    int model;
    // Elastoviscoplastic discretization type
    int discrtype;
    // Elastoviscoplastic convective discretization type
    int convecdiscrtype;
} vepl_controllers;


// Domains and distributed properties for elastoviscoplastic simulation
typedef struct higflow_elastoviscoplastic
{
    // Elastoviscoplastic controllers
    vepl_controllers contr;
    // Elastoviscoplastic parameters
    vepl_parameters par;
    // Distributed property for velocity derivative tensor
    distributed_property *dpD[DIM][DIM];
    // Distributed property for polymeric tensor
    distributed_property *dpS[DIM][DIM];
    // Distributed property for kernel tensor
    distributed_property *dpKernel[DIM][DIM];
    // Function to get the tensor
    real (*get_tensor)(Point center, int i, int j, real t);
    // Function to get the kernel transformation
    real (*get_kernel)(int dim, real lambda, real tol);
    // Function to get the inverse kernel transformation
    real (*get_kernel_inverse)(int dim, real lambda, real tol);
    // Function to get the kernel jacobian
    real (*get_kernel_jacobian)(int dim, real lambda, real tol);
    // User function to define the viscoelastic model
    void (*calculate_m_user)(real Re, real De, real beta, real Bi, real zeta, real epsilon, real Np, real tr, real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real M_aux[DIM][DIM], real tol, real smallTD, real SD);
} higflow_elastoviscoplastic;


// Parameters for shear-thickening suspensions simulation
typedef struct stsp_parameters
{
    // Alpha parameter
    real alpha;
    // Solvent viscosity (with proper viscosity units)
    real eta0;
    // Creation and destruction of particle pairs parameter
    real beta;
    // Extremal jamming point 1
    real chij1;
    // Extremal jamming point 2
    real chij2;
    // X0 parameter 
    real X0;
    // Particle volume fraction parameter
    real phi;
    // Critical particle pressure
    real Pic;
    //Standard deviation of the shear rate fluctuations
    real gdrms;
    //Particle size
    real apsize;
} stsp_parameters;

typedef enum suspension_model_type {
    USERSET_SM = -1,
    GW = 0,
    GW_WC = 1,
    GW_WC_IF = 2,
} suspension_model_type;

// Controllers for shear-thickening suspensions simulation
typedef struct stsp_controllers
{
    // Suspension model
    suspension_model_type model;
    // Discretization type
    tempdiscr_type discrtype;
    // Discretization type of the convective term of the microstructure tensor evolution equation
    cell_convecdiscr_type convecdiscrtype;
    // Discretization type for the volume fraction (particle migration equation)
    tempdiscr_type volfracdiscrtype;
    // Convective discretization type for the volume fraction (particle migration equation)
    cell_convecdiscr_type volfracconvecdiscrtype;
} stsp_controllers;


// Domains and distributed properties for shear-thickening suspensions simulations
typedef struct higflow_shear_thickening_suspension
{
    // Elastoviscoplastic controllers
    stsp_controllers contr;
    // Elastoviscoplastic parameters
    stsp_parameters par;
    //Sub-domain to simulation for volume fraction
    sim_domain  *sdphi;
    // Partitioned sub-domain to simulation for volume fraction
    psim_domain *psdphi;
    // Distributed property for volume fraction
    distributed_property *dpphi;
    // Distributed property for velocity derivative tensor
    distributed_property *dpD[DIM][DIM];
    // Distributed property for polymeric tensor S = tau -2 (eta_0)D
    distributed_property *dpS[DIM][DIM];
    // Distributed property for microstructure tensor (A = <nn>)
    distributed_property *dpA[DIM][DIM];
    // Function to get the tensor S
    real (*get_tensor)(Point center, int i, int j, real t);
    // Function to get the microstructure tensor A
    real (*get_tensor_A)(Point center, int i, int j, real t);
    // Function to get X-term direct interparticle forces
    real (*get_X)(Point center, real t, real X0, real chi, real chi_J);
    // Function to get the volume fraction
    real (*get_vol_frac)(Point center, real t);
    // Function to get the value of alpha = alpha(phi)
    real (*get_alpha)(Point center, real t, real alpha, real phi, real phircp);
} higflow_shear_thickening_suspension;

typedef struct higflow_multiphase_electroosmotic{
     // Mult electroosmotic controllers
    eo_controllers   contr;
    // Mult viscoelastic parameters
    eo_parameters    par0;
    eo_parameters    par1;
    // // Sub-domain to simulation for facets (EO source force)
    // sim_facet_domain     *sfdEOFeo[DIM];
    // // Partitioned sub-domain to simulation for facets (EO source force)
    // psim_facet_domain    *psfdEOFeo[DIM];
    // // Sub-domain to simulation for electro-osmotic phi
    // sim_domain           *sdEOphi;
    // // Partitioned sub-domain to simulation for electro-osmotic phi
    // psim_domain          *psdEOphi;
    // // Sub-domain to simulation for electro-osmotic psi
    // sim_domain           *sdEOpsi;
    // // Partitioned sub-domain to simulation for electro-osmotic psi
    // psim_domain          *psdEOpsi;
    // // Sub-domain to simulation for electrocosmotic nplus
    // sim_domain           *sdEOnplus;
    // // Partitioned sub-domain to simulation for electrocosmotic nminus
    // psim_domain          *psdEOnplus;
    // // Sub-domain to simulation for electro-osmotic nminus
    // sim_domain           *sdEOnminus;
    // // Partitioned sub-domain to simulation for electro-osmotic nminus
    // psim_domain          *psdEOnminus;
    // // Distributed property for source term 
    // distributed_property *dpFeo[DIM];
    // // Distributed property for applied potential 
    // distributed_property *dpphi;
    // // Distributed property for induced potential 
    // distributed_property *dppsi;
    // // Distributed property for positive charge concentration 
    // distributed_property *dpnplus;
    // // Distributed property for negative charge concentration 
    // distributed_property *dpnminus;
    // // Distributed property for reference charge concentration in phases
    // distributed_property *dpnref;
    // // stencil for phi
    // sim_stencil *stnphi;
    // // stencil for psi, Feo and nref
    // sim_stencil *stnpsi;
    // // stencil for nplus
    // sim_stencil *stnnplus;
    // // stencil for nminus
    // sim_stencil *stnnminus;
    // // Function to get reference charge concentration in phases
    // real (*get_electroosmotic_nref)(real fracvol, Point center, real t);
    // Function to get source term
    real (*get_multiphase_electroosmotic_source_term)(real fracvol, Point center, int dim, real t);
    // Function to get phi
    real (*get_multiphase_electroosmotic_phi)(real fracvol, Point center, real t);
    // Function to get psi
    real (*get_multiphase_electroosmotic_psi)(real fracvol, Point center, real t);
    // Function to get nplus
    real (*get_multiphase_electroosmotic_nplus)(real fracvol, Point center, real t);
    // Function to get nminus
    real (*get_multiphase_electroosmotic_nminus)(real fracvol, Point center, real t);
    // Function to get the source term at boundary
    real (*get_boundary_multiphase_electroosmotic_source_term)(real fracvol, int id, Point center, int dim, real t);
    // Function to get phi at boundary
    real (*get_boundary_multiphase_electroosmotic_phi)(real fracvol, int id, Point center, real t);
    // Function to get psi at boundary
    real (*get_boundary_multiphase_electroosmotic_psi)(real fracvol, int id, Point center, real t);
    // Function to get nplus at boundary
    real (*get_boundary_multiphase_electroosmotic_nplus)(real fracvol, int id, Point center, real t);
    // Function to get nminus at boundary
    real (*get_boundary_multiphase_electroosmotic_nminus)(real fracvol, int id, Point center, real t);
    // // Linear system solver for potential psi
    // solver                     *slvpsi;
    // solver                     *slvphi;
    // solver                     *slvnplus;
    // solver                     *slvnminus;
} higflow_multiphase_electroosmotic;


// Controllers for multiphase viscoelastic simulation
typedef struct mult_ve_controllers{
    // Controlers - phase 0
    // Viscoelastic model
    visc_model_type    model0;
    // Controlers - phase 1
    // Viscoelastic model
    visc_model_type    model1;

    // Viscoelastic discretization type
    inhomogenous_discr_type    discrtype;
    // Viscoelastic convective discretization type
    cell_convecdiscr_type    convecdiscrtype;
} mult_ve_controllers;


typedef struct higflow_multiphase_viscoelastic{
     // Mult viscoelastic controllers
    mult_ve_controllers   contr;
    // Mult viscoelastic parameters
    ve_parameters    par0;
    ve_parameters    par1;
    // Distributed property for velocity derivative tensor 
    // distributed_property *dpD[DIM][DIM];
    // // Distributed property for beta viscoelastic
    // distributed_property *dpbeta;
    // Distributed property for S viscoelastic - Momento
    // distributed_property *dpS[DIM][DIM];
    // // Distributed property for S viscoelastic - phase 0
    // distributed_property *dpS0[DIM][DIM];
    // // Distributed property for S viscoelastic - phase 1
    // distributed_property *dpS1[DIM][DIM];
    // Distributed property for kernel tensor - phase 0
    // distributed_property *dpKernel0[DIM][DIM];
    // // Distributed property for kernel tensor - phase 1
    // distributed_property *dpKernel1[DIM][DIM];
    // Distributed property for kernel tensor
    // distributed_property *dpKernel[DIM][DIM];

    // Function to get the tensor - phase 0
    // real (*get_tensor0)(Point center, int i, int j, real t);
    // // Function to get the tensor - phase 1
    // real (*get_tensor1)(Point center, int i, int j, real t);


    // Function to get the tensor
    real (*get_tensor_multiphase)(real fracvol, Point center, int i, int j, real t);
    // Function to get the kernel transformation
    real (*get_kernel)(int dim, real lambda, real tol);
    // Function to get the inverse kernel transformation
    real (*get_kernel_inverse)(int dim, real lambda, real tol);
    // Function to get the kernel jacobian
    real (*get_kernel_jacobian)(int dim, real lambda, real tol);
    // User function to define the viscoelastic model
    void (*calculate_m_user_multiphase)(real fracvol, real lambda[DIM], real jlambda[DIM],  real B[DIM][DIM], real M_aux[DIM][DIM], real Re, real trS, ve_parameters *par0, ve_parameters *par1);
} higflow_multiphase_viscoelastic;


typedef struct mult_parameters{
    // Capillary number
    real Ca;
} mult_parameters;

typedef struct mult_controllers{
    // Controlers - phase 0
    flow_type flowtype0;
    int eoflow0;

    // Controlers - phase 1
    flow_type flowtype1;
    int eoflow1;

    // Controllers - set whenever either phase is
    // viscoelastic
    int viscoelastic_either;
    // Controllers - set whenever either phase is
    // electro-osmotic
    int eoflow_either;
} mult_controllers;

// Domains and distributed properties for multiphase flows 
typedef struct higflow_multiphase{
    // Sub-domain to simulation for multiphase
    sim_domain                     *sdmult;
    // Partitioned sub-domain to simulation for multiphase
    psim_domain                    *psdmult;
    // Stencil for properties interpolation for multiphase
    sim_stencil           *stn;
    // Multiphase viscoelastic
    higflow_multiphase_viscoelastic ve;
    // Multiphase electroosmotic
    higflow_multiphase_electroosmotic eo;
    // Mult controllers
    mult_controllers      contr;
    // Mult parameters
    mult_parameters       par;
    // Distributed property for viscosity
    distributed_property *dpvisc;
    // Distributed property for density 
    distributed_property *dpdens;
    // Distributed property for volume fraction 
    distributed_property *dpfracvol;
    // Distributed property for auxiliary volume fraction
    distributed_property *dpfracvolaux;
    // Distributed property for normal
    distributed_property *dpnormal[DIM];
    // Distributed property for curvature 
    distributed_property *dpcurvature;
    // Distributed property for distance (plic)
    distributed_property *dpdistance;
    // Distributed property for interfacial force
    distributed_property *dpIF[DIM];
    
    // Function to get the viscosity for phase 0
    real (*get_viscosity0)(Point center, real t);
    // Function to get the viscosity for phase 1
    real (*get_viscosity1)(Point center, real t);
    // Function to get the density for phase 0
    real (*get_density0)(Point center, real t);
    // Function to get the density for phase 1
    real (*get_density1)(Point center, real t);
    // Function to get the volume fraction 
    real (*get_fracvol)(Point center, Point delta, real t);
} higflow_multiphase;

// Domains and distributed properties for non newtonian and multiphase flows
typedef struct higflow_extra_domains{
    // Sub-domain to simulation for extra domain
    sim_domain                     *sdED;
    // Partitioned sub-domain to simulation for extra domain
    psim_domain                    *psdED;
    // Sub-domain to simulation for extra domain
    // Generalized newtonian flows
    higflow_gen_newtonian          gn;
    // Multifase flows
    higflow_multiphase             mult;
    // Viscoelastic flows
    higflow_viscoelastic           ve;
    // Viscoelastic flows integral models
    higflow_viscoelastic_integral  im;
    // Electro-osmotic flows
    higflow_electroosmotic         eo;
    //Viscoelastic flows with variable viscosity
    higflow_viscoelastic_variable_viscosity vevv;
    //Viscoelastic flows with shear-banding
    higflow_viscoelastic_shear_banding vesb;
    //Elastoviscoplastic flows
    higflow_elastoviscoplastic vepl;
    //Shear-thickening suspensions
    higflow_shear_thickening_suspension stsp;
    // Stencil for properties interpolation for extra domains
    sim_stencil           *stn;
} higflow_extra_domains;

// #define COMPUTE_RESIDUALS
typedef struct sim_residuals sim_residuals;

// Navier-Stokes Solver Data Structure
typedef struct higflow_solver {
    // Parameters
    higflow_parameters         par;
    // Controlers 
    higflow_controllers        contr;
    // Computational cell
    higflow_compcell           cc;
    // Extra domains
    higflow_extra_domains      ed;
    // External functions
    higflow_functions          func;
    // Sub-domain to simulation for cells (pressure)
    sim_domain                 *sdp;
    // Partitioned sub-domain to simulation for cells (pressure)
    psim_domain                *psdp;
    // Sub-domain to simulation for cells (source force)
    sim_domain                 *sdF;
    // Partitioned sub-domain to simulation for cells (source force)
    psim_domain                *psdF;
    // Sub-domain to simulation for facets (velocities)
    sim_facet_domain           *sfdu[DIM];
    // Partitioned sub-domain to simulation for facets (velocities)
    psim_facet_domain          *psfdu[DIM];
    // Sub-domain to simulation for facets (source force)
    sim_facet_domain           *sfdF[DIM];
    // Partitioned sub-domain to simulation for facets (source force)
    psim_facet_domain          *psfdF[DIM];
    // Stencil for properties interpolation for pressure and velocities
    sim_stencil                *stn;
    // Stencil for properties interpolation for source force
    sim_stencil                *stnF;
    // Linear system solver for pressure
    solver                     *slvp;
    // Linear system solver for velocity
    solver                     *slvu[DIM];
    // Distributed properties from cells in the domain 
    // Distributed property for pressure 
    distributed_property       *dpp;
    // Distributed property for pressure difference
    distributed_property       *ddeltap;
    // Disttibuted property for source term 
    distributed_property       *dpF;
    // Distributed properties from facets in the domain 
    // Final velocity 
    distributed_property       *dpu[DIM];
    // Intermediate velocity 
    distributed_property       *dpustar[DIM];
    distributed_property       *dpuaux[DIM];
    // Source term for the facets
    distributed_property       *dpFU[DIM];
    // structure responsible for tracking residuals of consecutive iterations
    sim_residuals              *residuals;
} higflow_solver;

// *******************************************************************
// Navier-Stokes Create and Destroy Object
// *******************************************************************

// Create the NS object
higflow_solver *higflow_create (void); 

// Destroy the NS object
void higflow_destroy (higflow_solver *ns); 

// *******************************************************************
// Navier-Stokes Initilize Linear System Solver and MPI 
// *******************************************************************

// Initializing Navier-Stokes solver
void higflow_initialize(int *argc, char **argv[], int *myrank, int *ntasks);

// *******************************************************************
// Navier-Stokes Create and Realloc the Linear System Solvers
// *******************************************************************

// Create the linear system solvers
void higflow_create_solver(higflow_solver *ns);

// Realloc the linear system solvers
void higflow_realloc_solver(higflow_solver *ns);

// *******************************************************************
// Navier-Stokes Create Domain and Distributed Properties 
// *******************************************************************

// Create the simulation domain for NS object
void higflow_create_domain (higflow_solver *ns, int cache, int order); 

// Create the simulation domain for generalized newtonian flow
void higflow_create_domain_generalized_newtonian (higflow_solver *ns, int cache, int order, real (*get_viscosity)(Point center, real q, real t)); 

// Create the simulation domain for multiphase flow
void higflow_create_domain_multiphase(higflow_solver *ns, int cache, int order,
real (*get_viscosity0)(Point center, real t),
real (*get_viscosity1)(Point center, real t),
real (*get_density0)(Point center, real t),
real (*get_density1)(Point center, real t),
real (*get_fracvol)(Point center, Point delta, real t));

// Create the simulation domain for multiphase viscoelastic flow
void higflow_create_domain_multiphase_viscoelastic (higflow_solver *ns,
real (*get_tensor_multiphase)(real fracvol, Point center, int i, int j, real t),
real (*get_kernel)(int dim, real lambda, real tol),
real (*get_kernel_inverse)(int dim, real lambda, real tol),
real (*get_kernel_jacobian)(int dim, real lambda, real tol));

// Define the user function for Multiphase viscoelastic flow
void higflow_define_user_function_multiphase_viscoelastic (higflow_solver *ns, 
void (*calculate_m_user_multiphase)(real fracvol, real lambda[DIM], real jlambda[DIM],  real B[DIM][DIM], real M_aux[DIM][DIM], real Re, real trS, ve_parameters *par0, ve_parameters *par1));


// Create the simulation domain for viscoelastic flow
void higflow_create_domain_viscoelastic (higflow_solver *ns, int cache, int order,
real (*get_tensor)(Point center, int i, int j, real t),
real (*get_kernel)(int dim, real lambda, real tol),
real (*get_kernel_inverse)(int dim, real lambda, real tol),
real (*get_kernel_jacobian)(int dim, real lambda, real tol)); 

// Define the user function for viscoelastic flow
void higflow_define_user_function_viscoelastic (higflow_solver *ns, 
void (*calculate_m_user)(real lambda[DIM], real jlambda[DIM],  real B[DIM][DIM], real M_aux[DIM][DIM], real Re, real trS, ve_parameters *par));


// Create the simulation domain for viscoelastic flow integral model
void higflow_create_domain_viscoelastic_integral (higflow_solver *ns, int cache, int order,
real (*get_tensor)(Point center, int i, int j, real t));

// Create the simulation domain for multiphase electroosmotic flow
void higflow_create_domain_multiphase_electroosmotic (higflow_solver *ns, int cache, int order,
real (*get_multiphase_electroosmotic_source_term)(real fracvol, Point center, int dim, real t),
real (*get_multiphase_electroosmotic_phi)(real fracvol, Point center, real t),
real (*get_multiphase_electroosmotic_psi)(real fracvol, Point center, real t),
real (*get_multiphase_electroosmotic_nplus)(real fracvol, Point center, real t),
real (*get_multiphase_electroosmotic_nminus)(real fracvol, Point center, real t),
real (*get_multiphase_boundary_electroosmotic_source_term)(real fracvol, int id, Point center, int dim, real t),
real (*get_multiphase_boundary_electroosmotic_phi)(real fracvol, int id, Point center, real t),
real (*get_multiphase_boundary_electroosmotic_psi)(real fracvol, int id, Point center, real t),
real (*get_multiphase_boundary_electroosmotic_nplus)(real fracvol, int id, Point center, real t),
real (*get_multiphase_boundary_electroosmotic_nminus)(real fracvol, int id, Point center, real t));

// Create the simulation domain for electro-osmotic flow
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
real (*get_boundary_electroosmotic_nminus)(int id, Point center, real t));

// Create the simulation domain for viscoelastic flow with variable viscosity
void higflow_create_domain_viscoelastic_variable_viscosity(higflow_solver *ns, int cache, int order,
                                                           real (*get_tensor)(Point center, int i, int j, real t),
                                                           real (*get_kernel)(int dim, real lambda, real tol),
                                                           real (*get_kernel_inverse)(int dim, real lambda, real tol),
                                                           real (*get_kernel_jacobian)(int dim, real lambda, real tol),
                                                           real (*get_viscosity)(Point center, real q, real t, real beta, real struct_par),
                                                           real (*get_structpar)(Point center, real q, real t, real beta, real Phi, real Lambda, real Gamma));

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
                                                       real (*get_boundary_nB)(int id, Point center, real t));

// Create the simulation domain for elastoviscoplastic flow
void higflow_create_domain_elastoviscoplastic (higflow_solver *ns, int cache, int order,
                                               real (*get_tensor)(Point center, int i, int j, real t),
                                               real (*get_kernel)(int dim, real lambda, real tol),
                                               real (*get_kernel_inverse)(int dim, real lambda, real tol),
                                               real (*get_kernel_jacobian)(int dim, real lambda, real tol));

// Define the user function for viscoelastic flow
void higflow_define_user_function_elastoviscoplastic (higflow_solver *ns, 
                                                      void (*calculate_m_user)(real Re, real De, real beta, real Bi, real zeta, real epsilon, real Np, real tr, real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real M_aux[DIM][DIM], real tol, real smallTD, real SD));

// Create the simulation domain for shear-thickening suspension flow
void higflow_create_domain_shear_thickening_suspension (higflow_solver *ns, int cache, int order,
                                                        real (*get_tensor)(Point center, int i, int j, real t),
                                                        real (*get_tensor_A)(Point center, int i, int j, real t),
                                                        real (*get_X)(Point center, real t, real X0, real chi, real chi_J),
                                                        real (*get_vol_frac)(Point center, real t),
                                                        real (*get_alpha)(Point center, real t, real alpha, real phi, real phircp));
// Create the partitioned simulation sub-domain for NS object
void higflow_create_partitioned_domain (higflow_solver *ns, partition_graph *pg, int order); 

// Create the partitioned simulation sub-domain for NS object Generalized Newtonian
void higflow_create_partitioned_domain_generalized_newtonian (higflow_solver *ns, partition_graph *pg, int order); 

// Create the partitioned simulation sub-domain for NS object Multifase
void higflow_create_partitioned_domain_multiphase (higflow_solver *ns, partition_graph *pg, int order); 

// Create the partitioned simulation sub-domain for viscoelastic simulation
void higflow_create_partitioned_domain_viscoelastic (higflow_solver *ns, partition_graph *pg, int order); 

// Create the partitioned simulation sub-domain for NS object Viscoelastic integral model
void higflow_create_partitioned_domain_viscoelastic_integral (higflow_solver *ns, partition_graph *pg, int order); 

// Create the partitioned simulation sub-domain for NS object Electro-osmotic
void higflow_create_partitioned_domain_electroosmotic (higflow_solver *ns, partition_graph *pg, int order); 

// Create the partitioned simulation sub-domain for viscoelastic flow with variable viscosity simulation
void higflow_create_partitioned_domain_viscoelastic_variable_viscosity(higflow_solver *ns, partition_graph *pg, int order);

// Create the partitioned simulation sub-domain for viscoelastic flow with shear-banding simulation
void higflow_create_partitioned_domain_viscoelastic_shear_banding (higflow_solver *ns, partition_graph *pg, int order);

// Create the partitioned simulation sub-domain for elastoviscoplastic simulation
void higflow_create_partitioned_domain_elastoviscoplastic (higflow_solver *ns, partition_graph *pg, int order);

// Create the partitioned simulation sub-domain for shear-thickening suspension simulation
void higflow_create_partitioned_domain_shear_thickening_suspension (higflow_solver *ns, partition_graph *pg, int order);

// Create the distributed properties for generalized newtonian simulation
//void higflow_create_distributed_properties_generalized_newtonian(higflow_solver *ns); 

// Create the distributed properties for multiphase simulation
//void higflow_create_distributed_properties_multiphase(higflow_solver *ns); 

// Create the distributed properties for multiphase viscoelastic simulation
// void higflow_create_distributed_properties_multiphase_viscoelastic(higflow_solver *ns); 

// Create the distributed properties for viscoelastic simulation
// void higflow_create_distributed_properties_viscoelastic(higflow_solver *ns); 

// Create the distributed properties for viscoelastic simulation integral model
// void higflow_create_distributed_properties_viscoelastic_integral(higflow_solver *ns); 

// Create the distributed properties for electro-osmotic simulation
// void higflow_create_distributed_properties_electroosmotic(higflow_solver *ns); 

// Create the distributed properties for simulation of viscoelastic flows with variable viscosity
// void higflow_create_distributed_properties_viscoelastic_variable_viscosity(higflow_solver *ns);

// Create the distributed properties for simulation of viscoelastic flows with shear-banding
// void higflow_create_distributed_properties_viscoelastic_shear_banding(higflow_solver *ns);

// Create the distributed properties for elastoviscoplastic simulation
// void higflow_create_distributed_properties_elastoviscoplastic(higflow_solver *ns);

// Create the distributed properties for simulation of shear-thickening suspensions
// void higflow_create_distributed_properties_shear_thickening_suspension(higflow_solver *ns);

// Create the distributed properties for NS object
void higflow_create_distributed_properties(higflow_solver *ns); 

// // Create the distributed properties for the non newtonian simulation
// void higflow_create_distributed_properties_non_newtonian(higflow_solver *ns);

// Create the stencil for the NS object
void higflow_create_stencil(higflow_solver *ns); 

// Create the stencil for the extra domain
void higflow_create_stencil_for_extra_domain(higflow_solver *ns); 

// Create the stencil for multiphase
void higflow_create_stencil_multiphase(higflow_solver *ns);

// Create the stencils for electroosmotic flow
void higflow_create_stencils_electroosmotic(higflow_solver *ns);

// Partition table initialize
void higflow_partition_domain (higflow_solver *ns, partition_graph *pg, int numhigs, higio_amr_info *mi[numhigs], int ntasks, int myrank);

// Partition table multiphase initialize
void higflow_partition_domain_multiphase (higflow_solver *ns, partition_graph *pg, int numhigs, higio_amr_info *mi[numhigs], higio_amr_info *mi_mult[numhigs], int ntasks, int myrank);

// Use the given mesh information to create to create amr info about the multiphase domain
// Corresponds to a uniform higtree with the cell sizes of the last given level
higio_amr_info *higflow_create_amr_info_mult(higio_amr_info *mi);

// Set the external functions in the NS solver
void higflow_set_external_functions(higflow_solver *ns,
real (*get_pressure)(Point center, real t),
real (*get_velocity)(Point center, int dim, real t),
real (*get_source_term)(Point center, real t),
real (*get_facet_source_term)(Point center, int dim, real t),
real (*get_boundary_pressure)(int id, Point center, real t),
real (*get_boundary_velocity)(int id, Point center, int dim, real t),
real (*get_boundary_source_term)(int id, Point center, real t),
real (*get_boundary_facet_source_term)(int id, Point center, int dim, real t));

// Auxiliar
void __higflow_readstring(char s[], int max, FILE *file);

#endif
