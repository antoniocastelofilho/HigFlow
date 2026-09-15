#include "ns-example-2d.h"


/******************************************************************************************/
/******************************************************************************************/
/***************************** newtonian user functions ***********************************/
/******************************************************************************************/
/******************************************************************************************/

// initial pressure
real get_pressure(Point center, real t) {
    return 0.0;
}
// initial velocity
real get_velocity(Point center, int dim, real t) {
    real value = 0.0;
    real y = center[1];
    real x = center[0];
    if(flowtype == MULTIPHASE && eoflow_either == false)
        if(dim==0) value = y;

    if(flowtype == MULTIPHASE && eoflow_either == true) {
        real fracvol = compute_value_at_point(*sdmult_ptr, center, center, 1.0, *dpfracvol_ptr, *stnmult_ptr);
        real alpha0 = p_par->alpha_eo0.val;
        real kappa0 = p_par->kappa_eo0.val;
        real alpha1 = p_par->alpha_eo1.val;
        real kappa1 = p_par->kappa_eo1.val;
        real alpha = (1.0 - fracvol) * alpha0 + fracvol * alpha1;
        real kappa = (1.0 - fracvol) * kappa0 + fracvol * kappa1;
        real psi_pb = get_psi_channel_sol(alpha, kappa, psi_up, psi_down, y);
        if(dim==0)
            value = psi_pb - 0.5*(psi_up-psi_down)*y - 0.5*(psi_up+psi_down);
    }

    return value;
}
// initial source term
real get_source_term(Point center, real t) {
    return 0.0;
}
// initial facet source term
real get_facet_source_term(Point center, int dim, real t) {
    return 0.0;
}

real get_boundary_pressure(int id, Point center, real t) {
    real value;
    switch (id) {
    case 0:
        value = 0.0;
        break;
    case 1:
        value = 0.0;
        break;
    case 2:
        value = 0.0;
        break;
    case 3:
        value = 0.0;
        break;
    }
    return value;
}

real get_boundary_velocity(int id, Point center, int dim, real t) {
    real value = 0.0;
    real y = center[1];
    real x = center[0];
    real psi_pb, fracvol;

    switch (id) {
    case 0:
        switch (dim) {
        case 0:;
            real alpha, kappa;

            if(u_inlet == DIRICHLET) {
                if (flowtype == VISCOELASTIC && eoflow == false && (visc_model == LPTT || visc_model == FENE_P)) {
                    ///////// steady state fully developed solution
                    value = calc_u_ptt_fene(visc_model, p_par, center[1]);
                }
                else if (flowtype == NEWTONIAN && eoflow == true) {
                    alpha = p_par->alpha_eo.val;
                    kappa = p_par->kappa_eo.val;
                    psi_pb = get_psi_channel_sol(alpha, kappa, psi_up, psi_down, y);
                    value = psi_pb - 0.5*(psi_up-psi_down)*y - 0.5*(psi_up+psi_down);
                }
                else if(flowtype == MULTIPHASE && eoflow_either == false) value = 0.0;
                else if(flowtype == MULTIPHASE && eoflow_either == true) {
                    fracvol = compute_value_at_point(*sdmult_ptr, center, center, 1.0, *dpfracvol_ptr, *stnmult_ptr);
                    alpha = (1.0 - fracvol) * p_par->alpha_eo0.val + fracvol * p_par->alpha_eo1.val;
                    kappa = (1.0 - fracvol) * p_par->kappa_eo0.val + fracvol * p_par->kappa_eo1.val;
                    psi_pb = get_psi_channel_sol(alpha, kappa, psi_up, psi_down, y);
                    value = psi_pb - 0.5*(psi_up-psi_down)*y - 0.5*(psi_up+psi_down);
                }
                else {
                    value = 1.5 * (1.0 - center[1] * center[1]);
                }
            }
            else value = 0.0;

            break;
        case 1:
            value = 0.0;
            break;
        }
        break;
    case 1:
        switch (dim) {
        case 0:
            if(flowtype == MULTIPHASE && eoflow_either == false)
                value = 1.0;
            else
                value = 0.0;
            break;
        case 1:
            if(flowtype == MULTIPHASE && eoflow_either == false)
                value = 0.1;
            else
                value = 0.0;
            break;
        }
        break;
    case 2:
        switch (dim) {
        case 0:
            value = 0.0;
            break;
        case 1:
            value = 0.0;
            break;
        }
        break;
    case 3:
        switch (dim) {
        case 0:
            if(flowtype == MULTIPHASE && eoflow_either == false)
                value = 0.0;
            else
                value = 0.0;
            break;
        case 1:
            value = 0.0;
            break;
        }
        break;
    }
    return value;
}

real get_boundary_source_term(int id, Point center, real t) {
    return 0.0;
}

real get_boundary_facet_source_term(int id, Point center, int dim, real t) {
    return 0.0;
}


/******************************************************************************************/
/******************************************************************************************/
/******************** generalized newtonian user functions ********************************/
/******************************************************************************************/
/******************************************************************************************/

// initial viscosity
real get_viscosity_gn(Point center, real q, real t) {
    return 1.0;
}
