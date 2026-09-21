// Reading a velocity or a property at a point that is not a degree of freedom --
// every one of these goes through the stencil machinery in domain.c.
//
// THE NAMES ARE THE ONLY THING DISTINGUISHING NEAR-IDENTICAL VARIANTS, and the
// suffixes are positional, not descriptive: _left/_right pick the side, a trailing
// _2 or _22 selects how far and along which of dim/dim2, and the _4_ family works on
// the four-point stencils used by the higher-order convective schemes.  They all
// have the same signature shape, so calling the wrong one compiles cleanly and
// returns a value from the wrong place.  Check the definition, not the name.
//
// The `infacet` out-parameter on some variants reports whether the point landed on a
// facet of the domain; callers use it to decide whether the value is a degree of
// freedom or an interpolation.  Ignoring it is how a boundary value gets treated as
// an interior one.

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

//These functions are used to calculate the second deritatives of the stress tensor (crossed derivatives)

// Get the left facet value of a cell property
real compute_center_p_left_2_1(sim_domain *sdp, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpp, sim_stencil *stn);

// Get the right facet value of a cell property
real compute_center_p_right_2_1(sim_domain *sdp, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpp, sim_stencil *stn);

// Get the left facet value of a cell property
real compute_center_p_left_2_r(sim_domain *sdp, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpp, sim_stencil *stn);

// Get the right facet value of a cell property
real compute_center_p_right_2_r(sim_domain *sdp, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpp, sim_stencil *stn);

#endif
