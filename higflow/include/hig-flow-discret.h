// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Discretization - version 10/11/2016
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_DISCRET
#define HIG_FLOW_DISCRET

#include "hig-flow-kernel.h"
#include "hig-flow-eval.h"

// *******************************************************************
// Navier-Stokes Discretization
// *******************************************************************

// Compute the first derivative using the central method
real compute_facet_dudxc(Point delta, int dim, real alpha, real uc, real ul, real ur); 

// Compute the first derivative using the forward method
real compute_facet_dudxr(Point delta, int dim, real alpha, real uc, real ul, real ur); 

// Compute the first derivative using the backward method
real compute_facet_dudxl(Point delta, int dim, real alpha, real uc, real ul, real ur); 

// Compute the second derivative using the central method
real compute_facet_du2dx2(Point delta, int dim, real alpha, real u, real ul, real ur); 

// Compute the first derivative using the central method
real compute_dpdx_at_point(Point delta, int dim, real alpha, real valuel, real valueh); 

// Compute the first derivative using the forward method
real compute_dpdxr_at_point(Point delta, int dim, real alpha, real valuec, real valuer); 

// Compute the first derivative using the backward method
real compute_dpdxl_at_point(Point delta, int dim, real alpha, real valuel, real valuec); 

// Compute the mid point value
real compute_value_at_mid_point(real valuel, real valuer); 

// *******************************************************************
// Navier-Stokes Computational Cell
// *******************************************************************

// Computing the necessary term for the Navier-Stokes equation
void higflow_computational_cell(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int flid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]); 

// Computing the necessary term for the Generalized Newtonian Navier-Stokes equation
void higflow_computational_cell_gen_newt(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int flid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]); 

// Computing the necessary term for the Multifase Navier-Stokes equation
void higflow_computational_cell_multiphase(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int flid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]); 

// Computing the necessary term for the Viscoelastic Navier-Stokes equation
void higflow_computational_cell_viscoelastic(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int flid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]); 

// Computing the necessary term for the Navier-Stokes equation
void higflow_computational_cell_electroosmotic(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int flid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]); 

// Computing the necessary term for the ionic
void higflow_computational_cell_electroosmotic_ionic(higflow_solver *ns, sim_domain *sdn, int clid, Point ccenter, Point cdelta, int dim, distributed_property *dpn, sim_stencil *stnn);

// Computing the necessary term for the ionic multiphase
void higflow_computational_cell_multiphase_electroosmotic_ionic(higflow_solver *ns, sim_domain *sdn, int clid, Point ccenter, Point cdelta, int dim, distributed_property *dpn, sim_stencil *stnn,
                                                                real alphaeo, real alphaeol, real alphaeor, real Pe, real Pel, real Per);

// Computing the necessary term for the Viscoelastic Navier-Stokes equation
void higflow_computational_cell_viscoelastic_integral(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int flid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]);

//Computing the necessary terms for the Navier-Stokes equations of the viscoelastic flows with variable viscosity
void higflow_computational_cell_viscoelastic_variable_viscosity(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]);

// Computing the necessary term for the viscoelastic Navier-Stokes equation
void higflow_computational_cell_viscoelastic_shear_banding(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]);

//Computing the terms needed for the simulation of viscoelastic flows that exhibit shear-banding behaviour using the VCM model
void higflow_computational_cell_shear_banding_VCM_model(higflow_solver *ns, sim_domain *sdn, int clid, Point ccenter, Point cdelta, int dim, distributed_property *dpn);

// Computing the necessary term for the elastoviscoplastic Navier-Stokes equation
void higflow_computational_cell_elastoviscoplastic(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]);

// Computing the necessary term for the shear-thickening suspenions Navier-Stokes equation
void higflow_computational_cell_shear_thickening_suspensions(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]);

//Computational cell used to calculate the terms of the particle migration equation (evolution of the fraction volume)
void higflow_computational_cell_volume_fraction_suspensions(higflow_solver *ns, sim_domain *sdphi, int clid, Point ccenter, Point cdelta, int dim,  distributed_property *dpphi);

//Computing and storing the temperature value to be included in the Bousinessq term of the momentum equation
void higflow_computational_cell_bousinessq_term_momentum_equation(higflow_solver *ns, int clid, Point ccenter, distributed_property *dpn);

//Computational cell used to calculate the terms needed to solve the energy equation
void higflow_computational_cell_energy_equation(higflow_solver *ns, sim_domain *sdT, int clid, Point ccenter, Point cdelta, int dim, distributed_property *dpT);

#endif
