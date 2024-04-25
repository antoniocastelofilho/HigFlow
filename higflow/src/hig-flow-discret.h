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
void higflow_computational_cell(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]); 

// Computing the necessary term for the Generalized Newtonian Navier-Stokes equation
void higflow_computational_cell_gen_newt(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]); 

// Computing the necessary term for the Multifase Navier-Stokes equation
void higflow_computational_cell_multiphase(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]); 

// Computing the necessary term for the Viscoelastic Navier-Stokes equation
void higflow_computational_cell_viscoelastic(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]); 

// Computing the necessary term for the Navier-Stokes equation
void higflow_computational_cell_electroosmotic(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]); 

void higflow_computational_cell_electroosmotic_ionic(higflow_solver *ns, sim_domain *sdn, sim_domain *sdpsi, sim_domain *sdphi, int clid, Point ccenter, Point cdelta, int dim, distributed_property *dpn, distributed_property *dppsi, distributed_property *dpphi);

// Computing the necessary term for the Viscoelastic Navier-Stokes equation
void higflow_computational_cell_viscoelastic_integral(higflow_solver *ns, sim_domain *sdp, sim_facet_domain *sfdu[DIM], int fgid, Point fcenter, Point fdelta, int dim, distributed_property *dpu[DIM]);
#endif
