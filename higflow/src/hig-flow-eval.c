// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver - version 10/11/2016
// *******************************************************************
// *******************************************************************

#include "hig-flow-eval.h"

// *******************************************************************
// Navier-Stokes Computing Value Property
// *******************************************************************

real interpolate_from_stencil(distributed_property *dp, sim_stencil *stn) {
    /* Computing property at a point using the Moving Least Square method */
    /* Evaluating the property at the point using the stencil using the center of the cells */
    real res = 0.0;
    /* Geting value from boundary (Diriclet boundary condition) */
    res += -stn_get_rhs(stn);
    DEBUG_INSPECT(-stn_get_rhs(stn),%10.7f); 
    
    /* Getting the number of elements of internal centers cells */
    int numelems = stn_get_numelems(stn);
    DEBUG_INSPECT(stn_get_numelems(stn),%d);
    
    /* Getting the center and weights (MLS calculation) from neighborhood */
    for (int i = 0; i < numelems; i++) {
        /* Getting the identifier */
        int id = stn_get_id(stn, i);
        DEBUG_INSPECT(stn_get_id(stn,i),%d);
        /* Getting the weight value */
        real w = stn_get_val(stn, i);
        DEBUG_INSPECT(stn_get_val(stn,i),%10.7f);
        
        /* Getting the property value at the point with identifier id */
        real value = dp_get_value(dp, id);
        DEBUG_INSPECT(dp_get_value(dp,id),%10.7f);
        /* Computing the value at the given point */
        res += value * w;
    }
    return res;
}
// Get the value of a facet property
real compute_facet_value_at_point(sim_facet_domain *sfdu, Point center, Point p, real weight, distributed_property *dpu, sim_stencil *stn) {
    // Reset the stencil
//    stn_reset(stn);
    // Set the interpolater center
 //   sfd_set_interpolator_center(sfdu, center);
    // Get the stencil parameters
  //  sfd_get_stencil(sfdu, center, p, weight, stn);
    // Get the value of property at facet
   // real value1 = interpolate_from_stencil(dpu, stn);
    stn_reset(stn);
    // Get the stencil parameters
    sfd_get_stencil(sfdu, center, p, weight, stn);
    real value = dp_interpolate_from_stencil(dpu, stn);
    //DEBUG_INSPECT(p[0], %lf);
    //DEBUG_INSPECT(p[1], %lf);
    //DEBUG_INSPECT(stn_get_numelems(stn), %d);
    //DEBUG_INSPECT(value, %lf);
    return value;
}

// Get the value of a facet property
real compute_facet_u_2(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn) {
    // Set the left point
    Point p;
    POINT_ASSIGN(p, center);
    // Get the value at the right right point
    p[dim]       = center[dim]  + 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] + 0.5 * alpha * delta[dim2];
    real valuerr = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    //printf("u++ = %15.10f\n",valuerr);
    // Get the value at the left right point
    p[dim]       = center[dim]  + 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] - 0.5 * alpha * delta[dim2];
    real valuelr = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    //printf("u+- = %15.10f\n",valuelr);
    // Get the value at the right left point
    p[dim]       = center[dim]  - 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] + 0.5 * alpha * delta[dim2];
    real valuerl = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    //printf("u-+ = %15.10f\n",valuerl);
    // Get the value at the left left point
    p[dim]       = center[dim]  - 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] - 0.5 * alpha * delta[dim2];
    real valuell = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    //printf("u-- = %15.10f\n",valuell);

    real value  = 0.25 * (valuell + valuelr + valuerl + valuerr);
    return value;
}

// Get the value of a facet property
real compute_facet_u_4_right(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn) {
    // Set the left point
    Point p;
    POINT_ASSIGN(p, center);
    // Get the value at the right right point
    p[dim]       = center[dim]  + 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] + alpha * delta[dim2];
    real valuerr = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    //printf("u++ = %15.10f\n",valuerr);
    // Get the value at the left right point
    p[dim]       = center[dim]  + 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] ;
    real valuer = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    //printf("u+- = %15.10f\n",valuelr);
    // Get the value at the right left point
    p[dim]       = center[dim]  - 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] ;
    real valuel = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    //printf("u-+ = %15.10f\n",valuerl);
    // Get the value at the left left point
    p[dim]       = center[dim]  - 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] + alpha * delta[dim2];
    real valuell = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    //printf("u-- = %15.10f\n",valuell);
    real value  = 0.25 * (valuell + valuel + valuer + valuerr);
    return value;
}

// Get the value of a facet property
real compute_facet_u_4_corners_right(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn, int *infacet) {
    hig_facet f;
    real      value;
    // Set the left point
    Point p;
    POINT_ASSIGN(p, center);
    // Get the value at the right right point
    p[dim]       = center[dim]  + 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] + alpha * delta[dim2];
    *infacet     = sfd_get_facet_with_point(sfdu, p, &f);
    real valuerr = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    //printf("u++ = %15.10f\n",valuerr);
    // Get the value at the left right point
    p[dim]       = center[dim]  + 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] ;
    real valuer = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    //printf("u+- = %15.10f\n",valuelr);
    // Get the value at the right left point
    p[dim]       = center[dim]  - 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] ;
    real valuel = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    //printf("u-+ = %15.10f\n",valuerl);
    // Get the value at the left left point
    p[dim]       = center[dim]  - 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] + alpha * delta[dim2];
    real valuell = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    //printf("u-- = %15.10f\n",valuell);
    if (*infacet == 0) {
        value  = 0.5 * (valuell + valuel);
    } else {
        value  = 0.25 * (valuell + valuel + valuer + valuerr);
    }
    return value;
}

// Get the value of a facet property
real compute_facet_u_4_left(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn) {
    // Set the left point
    Point p;
    POINT_ASSIGN(p, center);
    // Get the value at the right right point
    p[dim]       = center[dim]  + 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] - alpha * delta[dim2];
    real valuerr = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    //printf("u++ = %15.10f\n",valuerr);
    // Get the value at the left right point
    p[dim]       = center[dim]  + 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] ;
    real valuer = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    //printf("u+- = %15.10f\n",valuelr);
    // Get the value at the right left point
    p[dim]       = center[dim]  - 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] ;
    real valuel = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    //printf("u-+ = %15.10f\n",valuerl);
    // Get the value at the left left point
    p[dim]       = center[dim]  - 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] - alpha * delta[dim2];
    real valuell = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    //printf("u-- = %15.10f\n",valuell);
    real value  = 0.25 * (valuell + valuel + valuer + valuerr);
    return value;
}

// Get the value of a facet property
real compute_facet_u_4_corners_left(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn, int *infacet) {
    hig_facet f;
    real      value;
    // Set the left point
    Point p;
    POINT_ASSIGN(p, center);
    // Get the value at the right right point
    p[dim]       = center[dim]  + 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] - alpha * delta[dim2];
    *infacet     = sfd_get_facet_with_point(sfdu, p, &f);
    real valuerr = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    //printf("u++ = %15.10f\n",valuerr);
    // Get the value at the left right point
    p[dim]       = center[dim]  + 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] ;
    real valuer = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    //printf("u+- = %15.10f\n",valuelr);
    // Get the value at the right left point
    p[dim]       = center[dim]  - 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] ;
    real valuel = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    //printf("u-+ = %15.10f\n",valuerl);
    // Get the value at the left left point
    p[dim]       = center[dim]  - 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] - alpha * delta[dim2];
    real valuell = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    //printf("u-- = %15.10f\n",valuell);
    if (*infacet == 0) {
        value  = 0.5 * (valuerr + valuer);
    } else {
        value  = 0.25 * (valuell + valuel + valuer + valuerr);
    }
    return value;
}

// Get the left facet value of a facet property
real compute_facet_u_left(sim_facet_domain *sfdu, Point center, Point delta, int dim, real alpha, distributed_property *dpu, sim_stencil *stn, int *infacet) {
    // Set the left point
    Point p;
    hig_facet f;
    POINT_ASSIGN(p, center);
    p[dim]      = center[dim] - alpha * delta[dim];
    *infacet    = sfd_get_facet_with_point(sfdu, p, &f);
    // Get the value at the left point
    real valuel = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    return valuel;
}

// Get the left facet value of a facet property
real compute_facet_u_left_2(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn) {
    // Set the left point
    Point p;
    POINT_ASSIGN(p, center);
    p[dim2]      = center[dim2] - alpha * delta[dim2];
    p[dim]       = center[dim] - 0.5*alpha * delta[dim];
    // Get the value at the left right point
    real valuell = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    p[dim2]      = center[dim2] - alpha * delta[dim2];
    p[dim]       = center[dim] + 0.5*alpha * delta[dim];
    // Get the value at the left right point
    real valuelr = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    // Get the value at the left point
    real valuel  = 0.5*(valuell + valuelr);
    return valuel;
}

// Get the left facet value of a facet property
real compute_facet_u_left_22(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn) {
    // Set the left point
    Point p;
    POINT_ASSIGN(p, center);
    //p[dim2]      = center[dim2] - alpha * delta[dim2];
    p[dim]       = center[dim] - 0.5*alpha * delta[dim];
    // Get the value at the left right point
    real valuell = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    //p[dim2]      = center[dim2] - alpha * delta[dim2];
    p[dim]       = center[dim] + 0.5*alpha * delta[dim];
    // Get the value at the left right point
    real valuelr = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    // Get the value at the left point
    real valuel  = 0.5*(valuell + valuelr);
    return valuel;
}

// Get the right facet value of a facet property
real compute_facet_u_right(sim_facet_domain *sfdu, Point center, Point delta, int dim, real alpha, distributed_property *dpu, sim_stencil *stn, int *infacet) {
    // Set the right point
    Point p;
    hig_facet f;
    POINT_ASSIGN(p, center);
    p[dim]      = center[dim] + alpha * delta[dim];
    *infacet    = sfd_get_facet_with_point(sfdu, p, &f);
    // Get the value at the right point
    real valueh = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    return valueh;
}

// Get the right facet value of a facet property
real compute_facet_u_right_2(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn) {
    // Set the right point
    Point p;
    POINT_ASSIGN(p, center);
    p[dim2]     = center[dim2] + alpha * delta[dim2];
    p[dim]      = center[dim] - 0.5*alpha * delta[dim];
    // Get the value at the right left point
    real valuerl = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    p[dim2]     = center[dim2] + alpha * delta[dim2];
    p[dim]      = center[dim] + 0.5*alpha * delta[dim];
    // Get the value at the right left point
    real valuerr = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    // Get the value at the right point
    real valuer  = 0.5*(valuerl + valuerr);
    return valuer;
}

// Get the value of a facet property
real compute_facet_u_2_left(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn) {
    // Set the left point
    Point p;
    POINT_ASSIGN(p, center);
    // Get the value at the left right point
    p[dim]       = center[dim]  + 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] - 0.5 * alpha * delta[dim2];
    real valuer = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    // Get the value at the right left point
    p[dim]       = center[dim]  - 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] - 0.5 * alpha * delta[dim2];
    real valuel = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    real value  = 0.5 * (valuel + valuer);
    return value;
}

// Get the value of a facet property
real compute_facet_u_2_right(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn) {
    // Set the left point
    Point p;
    POINT_ASSIGN(p, center);
    // Get the value at the left right point
    p[dim]       = center[dim]  + 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] + 0.5 * alpha * delta[dim2];
    real valuer = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    // Get the value at the right left point
    p[dim]       = center[dim]  - 0.5 * alpha * delta[dim];
    p[dim2]      = center[dim2] + 0.5 * alpha * delta[dim2];
    real valuel = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    real value  = 0.5 * (valuel + valuer);
    return value;
}

real compute_cell_u_bar_1(sim_domain *sdp, Point center, Point delta, int dim, real alpha, distributed_property *dpu, sim_stencil *stn, int *incell) {
    // Set the right point
    Point p;
    hig_cell *c;
    POINT_ASSIGN(p, center);
    p[dim]     = center[dim] + 0.5*alpha * delta[dim];
    c  = sd_get_cell_with_point(sdp, p);
    if (c == NULL) *incell = 0; else *incell = 1;
    // Get the value at the right point (i+1/2,j)
    real valuer = compute_value_at_point(sdp, center, p, 1.0, dpu, stn);
    // Get the value at the right point
    return valuer;
}

real compute_cell_u_bar_2(sim_domain *sdp, Point center, Point delta, int dim, real alpha, distributed_property *dpu, sim_stencil *stn, int *incell) {
    // Set the right point
    Point p;
    hig_cell *c;
    POINT_ASSIGN(p, center);
    p[dim]     = center[dim] - 0.5*alpha * delta[dim];
    c  = sd_get_cell_with_point(sdp, p);
    if (c == NULL) *incell = 0; else *incell = 1;
    // Get the value at the right point (i-1/2,j)
    real valuel = compute_value_at_point(sdp, center, p, 1.0, dpu, stn);
    // Get the value at the right point
    return valuel;
}


real compute_facet_u_bar_1(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn) {
    // Set the right point
    Point p;
    POINT_ASSIGN(p, center);
    p[dim]     = center[dim] + 0.5*alpha * delta[dim];
    p[dim2]    = center[dim2]  - 0.5*alpha * delta[dim2];
    // Get the value at the right left point
    real valuerl = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    p[dim]     = center[dim] + 0.5*alpha * delta[dim];
    p[dim2]    = center[dim2]  + 0.5*alpha * delta[dim2];
    // Get the value at the right left point
    real valuerr = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    // Get the value at the right point
    real valuer  = 0.5*(valuerl + valuerr);
    return valuer;
}

real compute_facet_u_bar_2(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn) {
    // Set the right point
    Point p;
    POINT_ASSIGN(p, center);
    p[dim]     = center[dim] - 0.5*alpha * delta[dim];
    p[dim2]    = center[dim2]  - 0.5*alpha * delta[dim2];
    // Get the value at the right left point
    real valuerl = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    p[dim]     = center[dim] - 0.5*alpha * delta[dim];
    p[dim2]    = center[dim2]  + 0.5*alpha * delta[dim2];
    // Get the value at the right left point
    real valuerr = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    // Get the value at the right point
    real valuer  = 0.5*(valuerl + valuerr);
    return valuer;
}


// Get the right facet value of a facet property
real compute_facet_u_right_22(sim_facet_domain *sfdu, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpu, sim_stencil *stn) {
    // Set the right point
    Point p;
    POINT_ASSIGN(p, center);
    //p[dim2]     = center[dim2] + alpha * delta[dim2];
    p[dim]      = center[dim] - 0.5*alpha * delta[dim];
    // Get the value at the right left point
    real valuerl = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    //p[dim2]     = center[dim2] + alpha * delta[dim2];
    p[dim]      = center[dim] + 0.5*alpha * delta[dim];
    // Get the value at the right left point
    real valuerr = compute_facet_value_at_point(sfdu, center, p, 1.0, dpu, stn);
    // Get the value at the right point
    real valuer  = 0.5*(valuerl + valuerr);
    return valuer;
}

// Get the value of a cell property
real compute_value_at_point(sim_domain *sd, Point center, Point p, real weight, distributed_property *dp, sim_stencil *stn) {
    // Reset the stencil
    stn_reset(stn);
    // Get the stencil parameters
    sd_get_stencil(sd, center, p, weight, stn);
    // Get the value of property at cell
    real value = dp_interpolate_from_stencil(dp, stn);
    return value;
}

// Get the left facet value of a cell property
real compute_center_p_left(sim_domain *sdp, Point center, Point delta, int dim, real alpha, distributed_property *dpp, sim_stencil *stn) {
    // Set the left point
    Point p;
    POINT_ASSIGN(p, center);
    p[dim] = center[dim] - alpha * delta[dim];
    // Get the value at the left point
    real valuel = compute_value_at_point(sdp, center, p, 1.0, dpp, stn);
    return valuel;
}

// Get the right facet value of a cell property
real compute_center_p_right(sim_domain *sdp, Point center, Point delta, int dim, real alpha, distributed_property *dpp, sim_stencil *stn) {
    // Set the right point
    Point p;
    POINT_ASSIGN(p, center);
    p[dim] = center[dim] + alpha * delta[dim];
    // Get the value at the right point
    real valueh = compute_value_at_point(sdp, center, p, 1.0, dpp, stn);
    return valueh;
}

// Get the left facet value of a cell property
real compute_center_p_left_2(sim_domain *sdp, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpp, sim_stencil *stn) {
    // Set the left point
    Point p;
    POINT_ASSIGN(p, center);
    // Get the value at the left left point
    p[dim2] = center[dim2] - alpha * delta[dim2];
    p[dim]  = center[dim] - 0.5*alpha * delta[dim];
    real valuel = compute_value_at_point(sdp, center, p, 1.0, dpp, stn);
    // Get the value at the left right point
    p[dim2] = center[dim2] - alpha * delta[dim2];
    p[dim]  = center[dim] + 0.5*alpha * delta[dim];
    real valuer = compute_value_at_point(sdp, center, p, 1.0, dpp, stn);
    // Eval at left point
    real value = 0.5*(valuel + valuer);
    return value;
}

// Get the right facet value of a cell property
real compute_center_p_right_2(sim_domain *sdp, Point center, Point delta, int dim, int dim2, real alpha, distributed_property *dpp, sim_stencil *stn) {
    // Set the right point
    Point p;
    POINT_ASSIGN(p, center);
    // Get the value at the right left point
    p[dim2] = center[dim2] + alpha * delta[dim2];
    p[dim]  = center[dim] - 0.5*alpha * delta[dim];
    real valuel = compute_value_at_point(sdp, center, p, 1.0, dpp, stn);
    // Get the value at the right right point
    p[dim2] = center[dim2] + alpha * delta[dim2];
    p[dim]  = center[dim] + 0.5*alpha * delta[dim];
    real valuer = compute_value_at_point(sdp, center, p, 1.0, dpp, stn);
    // Eval at left point
    real value = 0.5*(valuel + valuer);
    return value;
}

real compute_center_p_left_22(sim_domain *sdp, Point center, Point delta, int dim, real alpha, distributed_property *dpp, sim_stencil *stn, int *incell) {
    // Set the left point
    Point p;
    hig_cell *c;
    POINT_ASSIGN(p, center);
    p[dim] = center[dim] - alpha * delta[dim];
    c  = sd_get_cell_with_point(sdp, p);
    if (c == NULL) *incell = 0; else *incell = 1;
    // Get the value at the left point
    real valuel = compute_value_at_point(sdp, center, p, 1.0, dpp, stn);
    return valuel;
}

// Get the right facet value of a cell property
real compute_center_p_right_22(sim_domain *sdp, Point center, Point delta, int dim, real alpha, distributed_property *dpp, sim_stencil *stn, int *incell) {
    // Set the right point
    Point p;
    hig_cell *c;
    POINT_ASSIGN(p, center);
    p[dim] = center[dim] + alpha * delta[dim];
    c = sd_get_cell_with_point(sdp, p);
    if (c == NULL) *incell = 0; else *incell = 1;
    // Get the value at the right point
    real valueh = compute_value_at_point(sdp, center, p, 1.0, dpp, stn);
    return valueh;
}
