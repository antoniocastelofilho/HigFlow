// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Eval - version 10/11/2016
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_EVAL
#define HIG_FLOW_EVAL

#include "hig-flow-kernel.h"

// *******************************************************************
// Navier-Stokes Computing Value Property
// *******************************************************************

// Get the value of a facet property
real compute_facet_u_4_left(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn);

// Get the value of a facet property
real compute_facet_u_4_right(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn);

// Get the value of a facet property
real compute_facet_u_4_corners_left(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn, int *infacet);

// Get the value of a facet property
real compute_facet_u_4_corners_right(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn, int *infacet);

// Get the value of a facet property
real compute_facet_value_at_point(sim_facet_domain *sfdu, Point center, Point p, real weight, distributed_property *dpu, sim_stencil *stn); 

// Get the left facet value of a facet property
real compute_facet_u_2(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn); 

// Get the left facet value of a facet property
real compute_facet_u_left(sim_facet_domain *sfdu, Point center, Point delta, int dim, real alpha, distributed_property *dpu, sim_stencil *stn, int *infacet); 

// Get the right facet value of a facet property
real compute_facet_u_right(sim_facet_domain *sfdu, Point center, Point delta, int dim, real alpha, distributed_property *dpu, sim_stencil *stn, int *infacet); 

// Get the left facet value of a facet property
real compute_facet_u_left_2(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn); 

// Get the left facet value of a facet property
real compute_facet_u_left_22(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn); 

// Get the right facet value of a facet property
real compute_facet_u_right_2(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn); 

// Get the right facet value of a facet property
real compute_facet_u_right_22(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn); 

// Get the value of a facet property
real compute_facet_u_2_left(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn); 

// Get the value of a facet property
real compute_facet_u_2_right(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn);

// Get the value of a cell property
real compute_value_at_point(sim_domain *sd, Point center, Point p, real weight, distributed_property *dp, sim_stencil *stn); 

// Get the left facet value of a cell property
real compute_center_p_left(sim_domain *sdp, Point center, Point delta, int dim, real alpha, distributed_property *dpp, sim_stencil *stn); 

// Get the right facet value of a cell property
real compute_center_p_right(sim_domain *sdp, Point center, Point delta, int dim, real alpha, distributed_property *dpp, sim_stencil *stn); 

// Get the left facet value of a cell property
real compute_center_p_left_2(sim_domain *sdp, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpp, sim_stencil *stn); 

// Get the right facet value of a cell property
real compute_center_p_right_2(sim_domain *sdp, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpp, sim_stencil *stn); 

// Get the left facet value of a cell property
real compute_center_p_left_22(sim_domain *sdp, Point center, Point delta, int dim, real alpha, distributed_property *dpp, sim_stencil *stn, int *incell); 

// Get the right facet value of a cell property
real compute_center_p_right_22(sim_domain *sdp, Point center, Point delta, int dim, real alpha, distributed_property *dpp, sim_stencil *stn, int *incell); 

// Get the right facet value of a facet property for CUBISTA term
real compute_facet_u_bar_1(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn);

// Get the left facet value of a facet property for velocity CUBISTA term
real compute_facet_u_bar_2(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn);

// Get the right facet value of a facet property for tensor CUBISTA term
real compute_cell_u_bar_1(sim_domain *sdp, Point center, Point delta, int dim, real alpha, distributed_property *dpu, sim_stencil *stn, int *incell);

// Get the left facet value of a facet property for tensor CUBISTA term
real compute_cell_u_bar_2(sim_domain *sdp, Point center, Point delta, int dim, real alpha, distributed_property *dpu, sim_stencil *stn, int *incell);


#endif
