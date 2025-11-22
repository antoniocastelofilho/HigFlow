// *******************************************************************
// *******************************************************************
//  Utility functions for directory - version 03/2023
// *******************************************************************
// *******************************************************************

#include "../src/hig-flow-step.h"
#include "../src/hig-flow-io.h"
#include "../src/hig-flow-ic.h"
#include "../src/hig-flow-bc.h"
#define DEBUG
#include "Debug-c.h"
#include <stdio.h>
#include <stdlib.h>
#include <petsc.h>
#include "mpi.h"
#include <sys/resource.h>

typedef struct dimensional_unit{
    char kg;
    char m;
    char s;
    char K;
    char mol;
    char A;
    char cd;
}dimensional_unit;

typedef struct physical_quantity{
    real val;
    dimensional_unit unit;
}physical_quantity;

typedef struct physical_parameters{
    //////////////////// fluid parameters ////////////////////
    // uniform fluid density
    physical_quantity rho_ref;
    // reference length
    physical_quantity H;
    // characteristic velocity
    physical_quantity U;
    // reference viscosity
    physical_quantity mu_ref;
    //dynamic viscosity
    physical_quantity nu;
    // Reynolds number
    physical_quantity Re;
    // timestep
    physical_quantity dt;
    //reference time
    physical_quantity tau_time;
    // real timestep
    physical_quantity dt_physical;
    // gravitational acceleration
    physical_quantity g;
    // Froude number
    physical_quantity Fr;
    ////////////////////////// length parameters /////////////////////////
    // ### may need to be changed if the domain is not a cuboid
    // lengths
    physical_quantity L[DIM];
    //physical lengths
    physical_quantity L_physical[DIM];
    //center
    physical_quantity center[DIM];
    ////////////////////////// electro-osmotic parameters /////////////////////////
    // elementary charge
    physical_quantity e;
    // ion valence number
    physical_quantity Z;
    // reference potential
    physical_quantity zeta_ref;
    // vaccum permittivity
    physical_quantity epsilon_0;
    // relative permittivity of fluid
    physical_quantity epsilon_r;
    // permittivity of fluid
    physical_quantity epsilon_e;
    // boltzmann constant
    physical_quantity k_B;
    // temperature of flow
    physical_quantity T;
    // reference concentration
    physical_quantity n_ref;
    // ion diffusivity
    physical_quantity D;
    // electro-osmotic alpha
    physical_quantity alpha_eo;
    // electro-osmotic delta
    physical_quantity delta_eo;
    // Péclet number
    physical_quantity Pe;
    // Debye length
    physical_quantity lambda_D;
    // Debye number
    physical_quantity kappa_eo;
    // potential difference per unit of length
    physical_quantity Ex_physical;
    // potential difference relative to reference potential for every unit of reference length
    physical_quantity Ex;
    // Helmholtz-Smoluchowski velocity
    physical_quantity u_hs_physical;
    // // Helmholtz-Smoluchowski velocity relative to characteristic velocity
    // physical_quantity u_hs;
    // Conversion factor between eletro-osmotic and momentum systems
    physical_quantity G_x;
    ////////////////////////// viscoelastic parameters /////////////////////////
    // polymeric viscosity
    physical_quantity mu_p;
    // solvent viscosity
    physical_quantity mu_s;
    // Hooke's spring constant of polymer
    physical_quantity K_s;
    // polymer relaxation time
    physical_quantity lambda_time;
    // Deborah number
    physical_quantity De;
    // Ratio of solvent to total viscosity
    physical_quantity beta;
    // Giesekus parameter
    physical_quantity alpha_giesekus;
    // LPTT parameter
    physical_quantity epsilon_ptt;
    // LPTT parameter (wrongly called psi in the ve parameters) - slipperyness of polymer over continuum
    physical_quantity xi_ptt;
    // Kernel tolerance parameter (for kernel conformation if necessary)
    physical_quantity kernel_tol;
    // GPTT alpha
    physical_quantity alpha_gptt;
    // GPTT beta
    physical_quantity beta_gptt;
    // FENE L^2
    physical_quantity L2_fene;
    // e-FENE lambda
    physical_quantity lambda_fene;
    // e-FENE E
    physical_quantity E_fene;
    // steady state fully developed newtonian velocity normalized by flow rate
    physical_quantity Un_U;
    ////////////////////////// multiphase parameters /////////////////////////
    // surface tension coefficient
    physical_quantity sigma;
    // Bond (or Eotvos) number
    physical_quantity Bo;
    // Weber number
    physical_quantity We;
    // Capillary number
    physical_quantity Ca;
    //////// phase 0 ////////
    // density
    physical_quantity rho0;
    // viscosity
    physical_quantity mu0;
    //////// phase 1 ////////
    // density
    physical_quantity rho1;
    // viscosity
    physical_quantity mu1;
    ////////////////////////// multiphase viscoelastic parameters /////////////////////////
    //////// phase 0 ////////
    // polymeric viscosity
    physical_quantity mu_p0;
    // solvent viscosity
    physical_quantity mu_s0;
    // Hooke's spring constant of polymer
    physical_quantity K_s0;
    // polymer relaxation time
    physical_quantity lambda_time0;
    // Deborah number
    physical_quantity De0;
    // Ratio of solvent to total viscosity
    physical_quantity beta0;
    // Giesekus parameter
    physical_quantity alpha_giesekus0;
    // LPTT parameter
    physical_quantity epsilon_ptt0;
    // LPTT parameter (wrongly called psi in the ve parameters) - slipperyness of polymer over continuum
    physical_quantity xi_ptt0;
    // Kernel tolerance parameter (for kernel conformation if necessary)
    physical_quantity kernel_tol0;
    // GPTT alpha
    physical_quantity alpha_gptt0;
    // GPTT beta
    physical_quantity beta_gptt0;
    // FENE L^2
    physical_quantity L2_fene0;
    // e-FENE lambda
    physical_quantity lambda_fene0;
    // e-FENE E
    physical_quantity E_fene0;
    //////// phase 1 ////////
    // polymeric viscosity
    physical_quantity mu_p1;
    // solvent viscosity
    physical_quantity mu_s1;
    // Hooke's spring constant of polymer
    physical_quantity K_s1;
    // polymer relaxation time
    physical_quantity lambda_time1;
    // Deborah number
    physical_quantity De1;
    // Ratio of solvent to total viscosity
    physical_quantity beta1;
    // Giesekus parameter
    physical_quantity alpha_giesekus1;
    // LPTT parameter
    physical_quantity epsilon_ptt1;
    // LPTT parameter (wrongly called psi in the ve parameters) - slipperyness of polymer over continuum
    physical_quantity xi_ptt1;
    // Kernel tolerance parameter (for kernel conformation if necessary)
    physical_quantity kernel_tol1;
    // GPTT alpha
    physical_quantity alpha_gptt1;
    // GPTT beta
    physical_quantity beta_gptt1;
    // FENE L^2
    physical_quantity L2_fene1;
    // e-FENE lambda
    physical_quantity lambda_fene1;
    // e-FENE E
    physical_quantity E_fene1;
    ////////////////////////// multiphase electroosmotic parameters /////////////////////////
    // phase 0
    // reference potential
    physical_quantity zeta_ref0;
    // reference concentration
    physical_quantity n_ref0;
    // adimensional relative permittivity of fluid given by user
    physical_quantity perm0;
    // relative permittivity of fluid
    physical_quantity epsilon_r0;
    // permittivity of fluid
    physical_quantity epsilon_e0;
    // ion diffusivity
    physical_quantity D0;
    // electro-osmotic alpha
    physical_quantity alpha_eo0;
    // electro-osmotic delta
    physical_quantity delta_eo0;
    // Péclet number
    physical_quantity Pe0;
    // Debye length
    physical_quantity lambda_D0;
    // Debye number
    physical_quantity kappa_eo0;
    // potential difference per unit of length
    physical_quantity Ex_physical0;
    // potential difference relative to reference potential for every unit of reference length
    physical_quantity Ex0;
    // Helmholtz-Smoluchowski velocity
    physical_quantity u_hs_physical0;
    // Conversion factor between eletro-osmotic and momentum systems
    physical_quantity G_x0;
    // phase 1
    // reference potential
    physical_quantity zeta_ref1;
    // reference concentration
    physical_quantity n_ref1;
    // adimensional relative permittivity of fluid given by user
    physical_quantity perm1;
    // relative permittivity of fluid
    physical_quantity epsilon_r1;
    // permittivity of fluid
    physical_quantity epsilon_e1;
    // ion diffusivity
    physical_quantity D1;
    // electro-osmotic alpha
    physical_quantity alpha_eo1;
    // electro-osmotic delta
    physical_quantity delta_eo1;
    // Péclet number
    physical_quantity Pe1;
    // Debye length
    physical_quantity lambda_D1;
    // Debye number
    physical_quantity kappa_eo1;
    // potential difference per unit of length
    physical_quantity Ex_physical1;
    // potential difference relative to reference potential for every unit of reference length
    physical_quantity Ex1;
    // Helmholtz-Smoluchowski velocity
    physical_quantity u_hs_physical1;
    // Conversion factor between eletro-osmotic and momentum systems
    physical_quantity G_x1;
}physical_parameters;

/*#define RES_STORE_NUM_BASE 5
#define RES_STORE_NUM_POW 3
#define RES_STORE_NUM 125 //access to the last num residuals

#define NUM_TIMEINFOS 3
#define NUM_TIMESTATS 4

#define COMPUTE_MIDRANGE // channel only
#define COMPUTE_MIDLINE // channel only

typedef struct residual_normgroup_buffer{
    real res_max;
    real res_2;
    real res_1;
    real area;
} residual_normgroup_buffer;

typedef struct residual_buffer{
    #ifdef COMPUTE_MIDLINE
        residual_normgroup_buffer midline; // for channels only
        real midlinex;
    #endif // COMPUTE_MIDLINE
    #ifdef COMPUTE_MIDRANGE
        residual_normgroup_buffer midrange; // for channels only
        real leftx;
        real rightx;
    #endif // COMPUTE_MIDRANGE
    residual_normgroup_buffer all;
}residual_buffer;

typedef struct residual_timeinfo{
    real stored[RES_STORE_NUM]; // last residual, last base residuals, last base^2 residuals and so on;
    real avg[RES_STORE_NUM_POW]; // average of last base residuals, average of last base^2 residuals and so on;
    real geoavg[RES_STORE_NUM_POW]; // geometric average of last base residuals, geometric average of last 100 residuals and so on;
    //real med[RES_STORE_NUM_POW]; // median of last base residuals, median of last base^2 residuals and so on;
    real max[RES_STORE_NUM_POW]; // max of last base residuals, max of last base^2 residuals and so on;
    real min[RES_STORE_NUM_POW]; // min of last base residuals, min of last base^2 residuals and so on;
    real last; // last residual
} residual_timeinfo;

typedef struct residual_normgroup{
    residual_timeinfo *res_max;
    residual_timeinfo *res_2;
    residual_timeinfo *res_1;
} residual_normgroup;

typedef struct dp_residuals{
    #ifdef COMPUTE_MIDLINE
        residual_normgroup *midline; // for channels only
    #endif // COMPUTE_MIDLINE
    #ifdef COMPUTE_MIDRANGE
        residual_normgroup *midrange; // for channels only
    #endif // COMPUTE_MIDRANGE
    residual_normgroup *all;
}dp_residuals;

typedef struct sim_residuals{
    // dp_residuals *u_star[DIM];
    dp_residuals *u[DIM];
    // dp_residuals *p;

    // dp_residuals *D[DIM][DIM];
    // dp_residuals *S[DIM][DIM];
    dp_residuals *Kernel[DIM][DIM];

    // dp_residuals *phi;
    dp_residuals *psi;
    // dp_residuals *nplus;
    // dp_residuals *nminus;
    // dp_residuals *Feo[DIM];

    dp_residuals *fracvol;
    // dp_residuals *curvature;
    // dp_residuals *normal[DIM];
    // dp_residuals *IF[DIM];

    residual_buffer *res_buffer;
}sim_residuals;*/

typedef struct dp_copies{
    // distributed_property *u_star[DIM];
    distributed_property *u[DIM];
    // distributed_property *p;

    // distributed_property *D[DIM][DIM];
    // distributed_property *S[DIM][DIM];
    distributed_property *Kernel[DIM][DIM];

    // distributed_property *phi;
    distributed_property *psi;
    // distributed_property *nplus;
    // distributed_property *nminus;
    // distributed_property *Feo[DIM];

    distributed_property *fracvol;
    // distributed_property *curvature;
    // distributed_property *normal[DIM];
    // distributed_property *IF[DIM];
} dp_copies;

physical_parameters *create_initialize_physical_parameters(higflow_solver *ns, int myrank);
void free_physical_parameters(physical_parameters *p_par);

distributed_property  **create_velocity_copy(higflow_solver *ns);
void free_velocity_copy(higflow_solver *ns, distributed_property *u_copy[DIM]);
void copy_velocity_facet(higflow_solver *ns, distributed_property *u_copy[DIM], distributed_property *dpu[DIM]);

real get_fdp_value_at_point(higflow_solver *ns, distributed_property *dp, psim_facet_domain *psfd, Point p);
real max_dif_fdp(higflow_solver *ns, distributed_property *fdp1, distributed_property *fdp2, psim_facet_domain *psfd);
real max_dif_fdp_midline(higflow_solver *ns, distributed_property *fdp1, distributed_property *fdp2, psim_facet_domain *psfd, real midlinex);
real max_dif_fdp_midrange(higflow_solver *ns, distributed_property *fdp1, distributed_property *fdp2, psim_facet_domain *psfd, real midlinex, real Lx);
real max_dp(higflow_solver *ns, distributed_property *dp);
real max_fdp(higflow_solver *ns, distributed_property *fdp);

dp_copies *create_dp_copies(higflow_solver *ns);
void free_dp_copies(dp_copies *copies, higflow_solver *ns);
void copy_dps(higflow_solver *ns, dp_copies *copies);

real compute_mid_end_err_x(higflow_solver *ns, real midlinex, real Lx);
real compute_inlet_mid_err_x(higflow_solver *ns, real midlinex, real Lx);

/*sim_residuals *create_initialize_sim_residuals(higflow_solver *ns);
void free_sim_residuals(sim_residuals *sim_res, higflow_controllers hig_contr, int mrank);

void compute_residuals(higflow_solver *ns, sim_residuals *sim_res, dp_copies *copies, physical_parameters* p_par, int myrank);
void get_nameres(char **source, char *destination);
void write_residuals(sim_residuals *sim_res, char *nameres, real dt, int step, int myrank);*/

long int get_current_mem_usage();
void write_mem_usage(higflow_solver *ns, char *brief_desc);

real solve_un_u(visc_model_type visc_model, physical_parameters *p_par);
real calc_u_ptt_fene(visc_model_type visc_model, physical_parameters *p_par, real y);

void solve_psi_in(real *psi_in, physical_parameters *p_par, int npoints);
real get_psi_val_from_array(real *psi_in, real y, int npoints);

real get_psi_channel_sol(real alpha, real kappa, real psi_up, real psi_down, real y);
real get_nplus_channel_sol(real alpha, real kappa, real psi_up, real psi_down, real y);
real get_nminus_channel_sol(real alpha, real kappa, real psi_up, real psi_down, real y);

// writing to file functions
void write_init(higflow_solver *ns);
void write_xdmf(higflow_solver *ns);