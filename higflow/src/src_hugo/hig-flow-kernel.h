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
    int convec_type;
    // Cental velocity
    real ucell;
    // Velocity for convective term
    real v[DIM];
    // Cubista velocity for convective term
    real vc[DIM];
    // Source term
    real F;
    // Electroosmotic source term
    real Feo;
    // Density
    real dens;
    // Presure derivative
    real dpdx;
    // Central first order velocity derivative
    real dudxc[DIM];
    // Forward first order velocity derivative
    real dudxr[DIM];
    // Backward first order velocity derivative
    real dudxl[DIM];
    // Forward first order velocity derivative
    real dudxrr[DIM];
    // Backward first order velocity derivative
    real dudxll[DIM];
    // Second order velocity derivative
    real du2dx2[DIM];
    // Tensor derivative
    real dSdx[DIM];
    // Cental ionic concentration
    real ncell;
    // Ionic potential derivative
    real dpsidx;
    // Applied potential derivative
    real dphidx;
    // Ionic potential second order derivative
    real d2psidx2;
    // Applied potential second order derivative
    real d2phidx2;
    // First order ionic derivative
    real dndx;
    // Second order ionic derivative
    real d2ndx2;
    // Viscosity at left cell
    real viscl;
    // Viscosity at right cell
    real viscr;
    //  Curvature cell
    real curv;
    // Interfacial Force
    real IF;
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
    // Density
    real rho;
    // Viscosity
    real nu;
    // Reynolds number
    real Re;
    // Time
    real t;
    // Difference time
    real dt;
    // Current step number
    int step;
    int stepaux;
    // Maximal step number
    int finalstep;
    // Save time
    real ts;
    // Difference save time
    real dts;
    // Print time
    real tp;
    // Difference print time
    real dtp;
    // Frame number for print
    int frame;
    // File name for print
    char *nameprint;
    // File name for load
    char *nameload;
    // File name for save
    char *namesave;
} higflow_parameters;

// Methods controllers for simulation data structure
typedef struct higflow_controllers{
    // Projection type: 0 non incremental, 1 incremental
    int projtype;
    // Model type: 0 Newtonian, 1 Generalized Newtonian, 2 Multifase, 3 Viscoelastic
    int flowtype;
    // Diferencial Model type: 0 Single, 1 Electro-osmotic,
    int modelflowtype;
    // Temporal discretization type: 0 Explicit Euler, 1 Explicit RK2, .....
    int tempdiscrtype;
    // Spatial discretization type: 0 second order, 1 forth order
    int spatialdiscrtype;
    // Convective discretization type: 0 central, 1 first order, 2 second order
    int convecdiscrtype;
    // Second order discretization type: 0 Cubista, 1 Quick
    int secondconvecdiscrtype;
    // Desingualrization controllers for the pressure: 0 no desing. and 1 desing.
    int desingpressure;
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

// Domains and distributed properties for multiphase flows
typedef struct higflow_multiphase
{
    // Distributed property for velocity derivative tensor
    distributed_property *dpD[DIM][DIM];
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
    // Distributed property for beta viscoelastic
    distributed_property *dpbeta;
    // Distributed property for S viscoelastic
    distributed_property *dpS[DIM][DIM];
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

// Parameters for electro-osmotic simulation
typedef struct eo_parameters
{
    // alpha = e*z*zeta0/K*T
    real alpha;
    // delta = n0*e*z*H*H/(epsilon*zeta0)
    real delta;
    // Peclet number = U*H/D
    real Pe;
    // External electric field
    real dphidx;
} eo_parameters;

// Controllers for viscoelastic simulation
typedef struct eo_controllers
{
    // EO model Nernst-Planck, Poisson-Boltzmann or Debye-Hückel
    int eo_model;
    int convecdiscrtype;
} eo_controllers;

// Domains and distributed properties for electro-osmotic
typedef struct higflow_electroosmotic{
    // EO controllers
    eo_controllers contr;
    // EO parameters
    eo_parameters par;
    // Sub-domain to simulation for facets (EO source force)
    sim_facet_domain *sfdEOFeo[DIM];
    // Partitioned sub-domain to simulation for facets (EO source force)
    psim_facet_domain *psfdEOFeo[DIM];
    // Sub-domain to simulation for electro-osmotic phi
    sim_domain *sdEOphi;
    // Partitioned sub-domain to simulation for electro-osmotic phi
    psim_domain *psdEOphi;
    // Sub-domain to simulation for electro-osmotic psi
    sim_domain *sdEOpsi;
    // Partitioned sub-domain to simulation for electro-osmotic psi
    psim_domain *psdEOpsi;
    // Sub-domain to simulation for electrocosmotic nplus
    sim_domain *sdEOnplus;
    // Partitioned sub-domain to simulation for electrocosmotic nminus
    psim_domain *psdEOnplus;
    // Sub-domain to simulation for electro-osmotic nminus
    sim_domain *sdEOnminus;
    // Partitioned sub-domain to simulation for electro-osmotic nminus
    psim_domain *psdEOnminus;
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
    solver *slvpsi;
    solver *slvphi;
    solver *slvnplus;
    solver *slvnminus;
} higflow_electroosmotic;

// Parameters for viscoelastic simulation
typedef struct ve_parameters{
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
} ve_parameters;

// Controllers for viscoelastic simulation
typedef struct ve_controllers{
    // Viscoelastic model
    int model;
    // Viscoelastic discretization type
    int discrtype;
    // Viscoelastic convective discretization type
    int convecdiscrtype;
} ve_controllers;

// Domains and distributed properties for viscoelastic
typedef struct higflow_viscoelastic{
    // Viscoelastic controllers
    ve_controllers contr;
    // Viscoelastic parameters
    ve_parameters par;
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
} higflow_viscoelastic;

// Define the number of relaxation time for integral viscoelastic model
#define NDT 80

// Define the number of parameters for integral viscoelastic model
#define NRP 8

// Parameters for viscoelastic simulation integral model
typedef struct im_parameters
{
    // Deborah Number
    real De;
    // Number of relaxation time
    int M;
    // Parameters of Damping function
    real alpha;
    real beta;
    // Relaxation parameters for KBKZ model
    real a[NRP];
    real lambda[NRP];
    //  Integral points
    real eps_intpoints;
    real scorte;
    //
    real rho;
    real v_ref;
    real l_ref;
    // Relaxation parameters for KBKZ-Fractional model
    real Phi1;
    real Phi2;
    // Order of fractinal derivative for KBKZ-Fractional model
    real alpha_frac;
    real beta_frac;
} im_parameters;

// Controllers for viscoelastic simulation integral models
typedef struct im_controllers{
    // Viscoelastic model
    int model;
    // Relaxation model
    int model_H;
    // Viscoelastic discretization type
    int discrtype;
    // Viscoelastic convective discretization type
    int convecdiscrtype;
} im_controllers;

// Domains and distributed properties for viscoelastic integral models
typedef struct higflow_viscoelastic_integral{
    // Viscoelastic controllers
    im_controllers contr;
    // Viscoelastic parameters
    im_parameters par;
    // Discretization time
    real s[NDT + 1];
    // Discretization time (old step)
    real sold[NDT + 1];
    // Distributed property for velocity derivative tensor
    distributed_property *dpB[NDT + 1][DIM][DIM];
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

// Controllers for shear-thickening suspensions simulation
typedef struct stsp_controllers
{
    // Suspension model
    int model;
    // Discretization type
    int discrtype;
    // Discretization type of the convective term of the microstructure tensor evolution equation
    int convecdiscrtype;
    // Discretization type for the volume fraction (particle migration equation)
    int volfracdiscrtype;
    // Convective discretization type for the volume fraction (particle migration equation)
    int volfracconvecdiscrtype;
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


// Domains and distributed properties for non newtonian and multiphase flows
typedef struct higflow_extra_domains
{
    // Sub-domain to simulation for extra domain
    sim_domain *sdED;
    // Partitioned sub-domain to simulation for extra domain
    psim_domain *psdED;
    // Sub-domain to simulation for extra domain
    // Generalized newtonian flows
    higflow_gen_newtonian gn;
    // Multifase flows
    higflow_multiphase mult;
    // Viscoelastic flows
    higflow_viscoelastic ve;
    // Viscoelastic flows integral models
    higflow_viscoelastic_integral im;
    // Electro-osmotic flows
    higflow_electroosmotic eo;
    //Viscoelastic flows with variable viscosity
    higflow_viscoelastic_variable_viscosity vevv;
    //Viscoelastic flows with shear-banding
    higflow_viscoelastic_shear_banding vesb;
    //Elastoviscoplastic flows
    higflow_elastoviscoplastic vepl;
    //Shear-thickening suspensions
    higflow_shear_thickening_suspension stsp;
    // Stencil for properties interpolation for extra domains
    sim_stencil *stn;
} higflow_extra_domains;

// Navier-Stokes Solver Data Structure
typedef struct higflow_solver{
    // Parameters
    higflow_parameters par;
    // Controlers
    higflow_controllers contr;
    // Computational cell
    higflow_compcell cc;
    // Extra domains
    higflow_extra_domains ed;
    // External functions
    higflow_functions func;
    // Sub-domain to simulation for cells (pressure)
    sim_domain *sdp;
    // Partitioned sub-domain to simulation for cells (pressure)
    psim_domain *psdp;
    // Sub-domain to simulation for cells (source force)
    sim_domain *sdF;
    // Partitioned sub-domain to simulation for cells (source force)
    psim_domain *psdF;
    // Sub-domain to simulation for facets (velocities)
    sim_facet_domain *sfdu[DIM];
    // Partitioned sub-domain to simulation for facets (velocities)
    psim_facet_domain *psfdu[DIM];
    // Sub-domain to simulation for facets (source force)
    sim_facet_domain *sfdF[DIM];
    // Partitioned sub-domain to simulation for facets (source force)
    psim_facet_domain *psfdF[DIM];
    // Stencil for properties interpolation for pressure and velocities
    sim_stencil *stn;
    // Stencil for properties interpolation for source force
    sim_stencil *stnF;
    // Linear system solver for pressure
    solver *slvp;
    // Linear system solver for velocity
    solver *slvu[DIM];
    // Distributed properties from cells in the domain
    // Distributed property for pressure
    distributed_property *dpp;
    // Distributed property for pressure difference
    distributed_property *ddeltap;
    // Disttibuted property for source term
    distributed_property *dpF;
    // Distributed properties from facets in the domain
    // Final velocity
    distributed_property *dpu[DIM];
    // Intermediate velocity
    distributed_property *dpustar[DIM];
    distributed_property *dpuaux[DIM];
    // Source term for the facets
    distributed_property *dpFU[DIM];
} higflow_solver;

// *******************************************************************
// Navier-Stokes Create and Destroy Object
// *******************************************************************

// Create the NS object
higflow_solver *higflow_create(void);

// Destroy the NS object
void higflow_destroy(higflow_solver *ns);

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
void higflow_create_domain(higflow_solver *ns, int cache, int order);

// Create the simulation domain for generalized newtonian flow
void higflow_create_domain_generalized_newtonian(higflow_solver *ns, int cache, int order, real (*get_viscosity)(Point center, real q, real t));

// Create the simulation domain for multiphase flow
void higflow_create_domain_multiphase(higflow_solver *ns, int cache, int order,
                                      real (*get_viscosity0)(Point center, real t),
                                      real (*get_viscosity1)(Point center, real t),
                                      real (*get_density0)(Point center, real t),
                                      real (*get_density1)(Point center, real t),
                                      real (*get_fracvol)(Point center, Point delta, real t));

// Create the simulation domain for viscoelastic flow
void higflow_create_domain_viscoelastic(higflow_solver *ns, int cache, int order,
                                        real (*get_tensor)(Point center, int i, int j, real t),
                                        real (*get_kernel)(int dim, real lambda, real tol),
                                        real (*get_kernel_inverse)(int dim, real lambda, real tol),
                                        real (*get_kernel_jacobian)(int dim, real lambda, real tol));

// Define the user function for viscoelastic flow
void higflow_define_user_function_viscoelastic(higflow_solver *ns,
                                               void (*calculate_m_and_b_user)(real Re, real De, real beta, real tr, real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real M_aux[DIM][DIM], real tol));

// Create the simulation domain for viscoelastic flow integral model
void higflow_create_domain_viscoelastic_integral(higflow_solver *ns, int cache, int order,
                                                 real (*get_tensor)(Point center, int i, int j, real t));

// Create the simulation domain for electro-osmotic flow
void higflow_create_domain_electroosmotic(higflow_solver *ns, int cache, int order,
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
void higflow_create_partitioned_domain(higflow_solver *ns, partition_graph *pg, int order);

// Create the partitioned simulation sub-domain for NS object Generalized Newtonian
void higflow_create_partitioned_domain_generalized_newtonian(higflow_solver *ns, partition_graph *pg, int order);

// Create the partitioned simulation sub-domain for NS object Multifase
void higflow_create_partitioned_domain_multiphase(higflow_solver *ns, partition_graph *pg, int order);

// Create the partitioned simulation sub-domain for viscoelastic simulation
void higflow_create_partitioned_domain_viscoelastic(higflow_solver *ns, partition_graph *pg, int order);

// Create the partitioned simulation sub-domain for NS object Viscoelastic integral model
void higflow_create_partitioned_domain_viscoelastic_integral(higflow_solver *ns, partition_graph *pg, int order);

// Create the partitioned simulation sub-domain for NS object Electro-osmotic
void higflow_create_partitioned_domain_electroosmotic(higflow_solver *ns, partition_graph *pg, int order);

// Create the partitioned simulation sub-domain for viscoelastic flow with variable viscosity simulation
void higflow_create_partitioned_domain_viscoelastic_variable_viscosity(higflow_solver *ns, partition_graph *pg, int order);

// Create the partitioned simulation sub-domain for viscoelastic flow with shear-banding simulation
void higflow_create_partitioned_domain_viscoelastic_shear_banding (higflow_solver *ns, partition_graph *pg, int order);

// Create the partitioned simulation sub-domain for elastoviscoplastic simulation
void higflow_create_partitioned_domain_elastoviscoplastic (higflow_solver *ns, partition_graph *pg, int order);

// Create the partitioned simulation sub-domain for shear-thickening suspension simulation
void higflow_create_partitioned_domain_shear_thickening_suspension (higflow_solver *ns, partition_graph *pg, int order);

// Create the distributed properties for generalized newtonian simulation
void higflow_create_ditributed_properties_generalized_newtonian(higflow_solver *ns);

// Create the distributed properties for multiphase simulation
void higflow_create_ditributed_properties_multiphase(higflow_solver *ns);

// Create the distributed properties for viscoelastic simulation
void higflow_create_ditributed_properties_viscoelastic(higflow_solver *ns);

// Create the distributed properties for viscoelastic simulation integral model
void higflow_create_ditributed_properties_viscoelastic_integral(higflow_solver *ns);

// Create the distributed properties for electro-osmotic simulation
void higflow_create_ditributed_properties_electroosmotic(higflow_solver *ns);

// Create the distributed properties for simulation of viscoelastic flows with variable viscosity
void higflow_create_ditributed_properties_viscoelastic_variable_viscosity(higflow_solver *ns);

// Create the distributed properties for simulation of viscoelastic flows with shear-banding
void higflow_create_ditributed_properties_viscoelastic_shear_banding(higflow_solver *ns);

// Create the distributed properties for elastoviscoplastic simulation
void higflow_create_ditributed_properties_elastoviscoplastic(higflow_solver *ns);

// Create the distributed properties for simulation of shear-thickening suspensions
void higflow_create_ditributed_properties_shear_thickening_suspension(higflow_solver *ns);
// Create the distributed properties for NS object
void higflow_create_ditributed_properties(higflow_solver *ns);

// Create the distributed properties for the non newtonian simulation
void higflow_create_ditributed_properties_non_newtonian(higflow_solver *ns);

// Create the stencil for the NS object
void higflow_create_stencil(higflow_solver *ns);

// Create the stencil for the extra domain
void higflow_create_stencil_for_extra_domain(higflow_solver *ns);

// Partition table initialize
void higflow_partition_domain(higflow_solver *ns, partition_graph *pg, int numhigs, higio_amr_info *mi[numhigs], int ntasks, int myrank);

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
