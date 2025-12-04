#include "ns-example-2d.h"

/******************************************************************************************/
/******************************************************************************************/
/*************************** electroosmotic user functions ********************************/
/******************************************************************************************/
/******************************************************************************************/

// initial electroosmotic source term
real get_electroosmotic_source_term(Point center, int dim, real t) {
    return 0.0;
}
// initial applied potential phi
real get_electroosmotic_phi(Point center, real t) {
    return 0.0;
}
// initial induced potential psi
real get_electroosmotic_psi(Point center, real t) {
    real y = center[1];
    real kappa = p_par->kappa_eo.val;
    real alpha = p_par->alpha_eo.val;
    real value = 0.0;

    //value = get_psi_channel_sol(alpha, kappa, psi_up, psi_down, y);
    
    return value;
}
// initial positive charge concentration n^+
real get_electroosmotic_nplus(Point center, real t) {
    real y = center[1];
    real kappa = p_par->kappa_eo.val;
    real alpha = p_par->alpha_eo.val;
    real value = 1.0;

    //value = get_nplus_channel_sol(alpha, kappa, psi_up, psi_down, y);
    
    return value;
}
// initial negative charge concentration n^-
real get_electroosmotic_nminus(Point center, real t) {
    real y = center[1];
    real kappa = p_par->kappa_eo.val;
    real alpha = p_par->alpha_eo.val;
    real value = 1.0;

    //value = get_nminus_channel_sol(alpha, kappa, psi_up, psi_down, y);

    return value;
}

real get_boundary_electroosmotic_source_term(int id, Point center, int dim, real t) {
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

real get_boundary_electroosmotic_phi(int id, Point center, real t) {
    real value;
    switch (id) {
    case 0:
        value = p_par->Ex.val * p_par->L[0].val;
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

real get_boundary_electroosmotic_psi(int id, Point center, real t) {
    real value;
    real y = center[1];
    real alpha = p_par->alpha_eo.val;
    real kappa = p_par->kappa_eo.val;
    //real psi_pb  = -cosh(kappa*y)/(cosh(kappa));
    //real psi_pb = get_psi_val_from_array(psi_in, y, NLPB_SIZE);
    real psi_pb;
    
    switch (id) {
    case 0:
        if(psi_inlet == DIRICHLET) {
            psi_pb = get_psi_channel_sol(alpha, kappa, psi_up, psi_down, y);
            value = psi_pb;
        }
        else value = 0.0;
        break;
    case 1:
        value = psi_up;
        break;
    case 2:
        value = 0.0;
        break;
    case 3:
        value = psi_down;
        break;
    }
    return value;
}

real get_boundary_electroosmotic_nplus(int id, Point center, real t) {
    real value;
    real y = center[1];
    real kappa = p_par->kappa_eo.val;
    //real psi_pb  = -cosh(kappa*y)/(cosh(kappa));
    //real psi_pb = get_psi_val_from_array(psi_in, y, NLPB_SIZE);
    real np_pb;
    real alpha = p_par->alpha_eo.val;
    switch (id) {
    case 0:
        value = get_nplus_channel_sol(alpha, kappa, psi_up, psi_down, y);
        //value = 1.0;
        //value = 0.0;
        break;
    case 1:
        value = exp(-alpha * psi_up);
        break;
    case 2:
        value = 0.0;
        //value = exp(-alpha*psi_pb);
        break;
    case 3:
        value = exp(-alpha * psi_down);
        break;
    }
    return value;
}

// Function to get nminus at boundary
real get_boundary_electroosmotic_nminus(int id, Point center, real t) {
    real value;
    real y = center[1];
    real kappa = p_par->kappa_eo.val;
    //real psi_pb  = -cosh(kappa*y)/(cosh(kappa));
    //real psi_pb = get_psi_val_from_array(psi_in, y, NLPB_SIZE);
    real nm_pb;
    real alpha = p_par->alpha_eo.val;
    switch (id) {
    case 0:
        value = get_nminus_channel_sol(alpha, kappa, psi_up, psi_down, y);
        //value = 1.0;
        //value = 0.0;
        break;
    case 1:
        value = exp(alpha * psi_up);
        break;
    case 2:
        value = 0.0;
        //value = exp(-alpha*psi_pb);
        break;
    case 3:
        value = exp(alpha * psi_down);
        break;
    }
    return value;
}

real get_electroosmotic_permittivity(Point center, real t) {
    return 1.0;
}


/******************************************************************************************/
/******************************************************************************************/
/*********************** multiphase electroosmotic user functions ***************************/
/******************************************************************************************/
/******************************************************************************************/


// initial electroosmotic source term
real get_multiphase_electroosmotic_source_term(real fracvol, Point center, int dim, real t) {
    return 0.0;
}
// initial applied potential phi
real get_multiphase_electroosmotic_phi(real fracvol, Point center, real t) {
    return 0.0;
}
// initial induced potential psi
real get_multiphase_electroosmotic_psi(real fracvol, Point center, real t) {
    return 0.0;
}
// initial positive charge concentration n^+
real get_multiphase_electroosmotic_nplus(real fracvol, Point center, real t) {
    return 1.0;
}
// initial negative charge concentration n^-
real get_multiphase_electroosmotic_nminus(real fracvol, Point center, real t) {
    return 1.0;
}

real get_boundary_multiphase_electroosmotic_source_term(real fracvol, int id, Point center, int dim, real t) {
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

real get_boundary_multiphase_electroosmotic_phi(real fracvol, int id, Point center, real t) {
    real value, Ex;
    switch (id) {
    case 0:
        Ex = (1.0 - fracvol) * p_par->Ex0.val + fracvol * p_par->Ex1.val;
        value = Ex * p_par->L[0].val;
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

real get_boundary_multiphase_electroosmotic_psi(real fracvol, int id, Point center, real t) {
    real value;
    real y = center[1];
    real alpha = (1.0 - fracvol) * p_par->alpha_eo0.val + fracvol * p_par->alpha_eo1.val;
    real kappa = (1.0 - fracvol) * p_par->kappa_eo0.val + fracvol * p_par->kappa_eo1.val;
    real psi_pb;
    switch (id) {
    case 0:
        psi_pb = get_psi_channel_sol(alpha, kappa, psi_up, psi_down, y);
        //value = psi_pb;
        value = 0.0;
        break;
    case 1:
        value = psi_up;
        break;
    case 2:
        value = 0.0;
        //value = psi_pb;
        break;
    case 3:
        value = psi_down;
        break;
    }
    return value;
}

real get_boundary_multiphase_electroosmotic_nplus(real fracvol, int id, Point center, real t) {
    real value;
    real y = center[1];
    real kappa = (1.0 - fracvol) * p_par->kappa_eo0.val + fracvol * p_par->kappa_eo1.val;
    real alpha = (1.0 - fracvol) * p_par->alpha_eo0.val + fracvol * p_par->alpha_eo1.val;
    switch (id) {
    case 0:
        value = get_nplus_channel_sol(alpha, kappa, psi_up, psi_down, y);
        //value = 1.0;
        //value = 0.0;
        break;
    case 1:
        value = exp(-alpha * psi_up);
        break;
    case 2:
        value = 0.0;
        //value = exp(-alpha*psi_pb);
        break;
    case 3:
        value = exp(-alpha * psi_down);
        break;
    }
    return value;
}

// Function to get nminus at boundary
real get_boundary_multiphase_electroosmotic_nminus(real fracvol, int id, Point center, real t) {
    real value;
    real y = center[1];
    real kappa = (1.0 - fracvol) * p_par->kappa_eo0.val + fracvol * p_par->kappa_eo1.val;
    real alpha = (1.0 - fracvol) * p_par->alpha_eo0.val + fracvol * p_par->alpha_eo1.val;
    switch (id) {
    case 0:
        value = get_nminus_channel_sol(alpha, kappa, psi_up, psi_down, y);
        //value = 1.0;
        //value = 0.0;
        break;
    case 1:
        value = exp(alpha * psi_up);
        break;
    case 2:
        value = 0.0;
        //value = exp(-alpha*psi_pb);
        break;
    case 3:
        value = exp(alpha * psi_down);
        break;
    }
    return value;
}

real get_multiphase_electroosmotic_permittivity(real fracvol, Point center, real t) {
    return p_par->perm0.val * (1.0 - fracvol) + p_par->perm1.val * fracvol;
}