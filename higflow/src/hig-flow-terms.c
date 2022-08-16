// *******************************************************************
//  HiG-Flow Solver - version 25/01/2022
// *******************************************************************
#include "hig-flow-terms.h"
// *******************************************************************
// Navier-Stokes Equation Components
// *******************************************************************
// Pressure term contribution for the Navier-Stokes equation
real higflow_pressure_term(higflow_solver *ns) {
    // Set the pressure term
    real value = 0.0;
    if (ns->contr.projtype == 1) {
       value = ns->cc.dpdx;
       if (ns->contr.flowtype == 2) {
           value /= ns->cc.dens;
       }
    }
    return value;
}

// Cell source term contribution for the Navier-Stokes equation
real higflow_source_term(higflow_solver *ns) {
    // Set the source term
    real value = ns->cc.F;
    if (ns->contr.flowtype == 2) {
        value /= ns->cc.dens;
    }
    return value;
}

// Cell term contribution for the interfacial tension
real higflow_interfacial_tension_term(higflow_solver *ns) {
    real Bo, value;
    /*Bo = 10.0;
    value = ns->cc.IF;
    value = value/(Bo*ns->cc.dens);*/
    value = 0.0;
    return value;
}

// Cell term contribution for the gravity
real higflow_gravity_term(higflow_solver *ns) {
    //real Fr = 1.0;
    // real value = ns->cc.curv;
    //real value = 1.0/pow(Fr,2.0);
    //real value = 0.98;
    real value = 0.0;
    return value;
}

// Tensor term contribution for the Navier-Stokes equation
real higflow_tensor_term(higflow_solver *ns) {
    // Set the tensor term
    real value = 0.0;
    // Two-phase contribution
    if (ns->contr.flowtype == 2) {
        for (int dim2 = 0; dim2 < DIM; dim2++) {
            value += ns->cc.dSdx[dim2];
        }
        value /= ns->cc.dens;
    }
    // Non Newtonian contribution
    if (ns->contr.flowtype == 3) {
        for (int dim2 = 0; dim2 < DIM; dim2++) {
            value += ns->cc.dSdx[dim2];
        }
    }
    // Non Newtonian contribution
    if (ns->contr.flowtype == 4) {
        for (int dim2 = 0; dim2 < DIM; dim2++) {
            value += ns->cc.dSdx[dim2];
        }
    }
    return value;
}

// Convective term: Quick scheme
real higflow_convective_term_quick(higflow_solver *ns, Point delta, int dim) {
    // Set the convective term
    real value = 0.0;
    for (int dim2 = 0; dim2 < DIM; dim2++) {
        value += (ns->cc.v[dim2]+fabs(ns->cc.v[dim2]))*0.1875*ns->cc.dudxl[dim2]; 
        value += (ns->cc.v[dim2]+fabs(ns->cc.v[dim2]))*0.3750*ns->cc.dudxr[dim2];        
        value -= (ns->cc.v[dim2]+fabs(ns->cc.v[dim2]))*0.0625*ns->cc.dudxrr[dim2];    
        value += (ns->cc.v[dim2]-fabs(ns->cc.v[dim2]))*0.3750*ns->cc.dudxl[dim2]; 
        value += (ns->cc.v[dim2]-fabs(ns->cc.v[dim2]))*0.1875*ns->cc.dudxr[dim2];    
        value -= (ns->cc.v[dim2]-fabs(ns->cc.v[dim2]))*0.0625*ns->cc.dudxll[dim2];
    }
    return value;
}

// Convective term: Cubista scheme
real higflow_convective_term_cubista(higflow_solver *ns, Point delta, int dim) {
    // Set the convective term
    real value = 0.0;
    for (int dim2 = 0; dim2 < DIM; dim2++) {
        value += ns->cc.vc[dim2];
    }
    return value;
}

// Convective term: Modified Coefficient Upwind scheme
real higflow_convective_term_mcupwind(higflow_solver *ns, Point delta, int dim) {
    // Set the convective term
    real value = 0.0;
    for (int dim2 = 0; dim2 < DIM; dim2++) {
        value += 0.5*(ns->cc.v[dim2]+fabs(ns->cc.v[dim2]))*ns->cc.dudxl[dim2];
        value += 0.5*(ns->cc.v[dim2]-fabs(ns->cc.v[dim2]))*ns->cc.dudxr[dim2];
        if (ns->cc.v[dim2] > 0.0) {
            if (fabs(ns->cc.dudxl[dim2]) > 1.0e-8) {
                value += 0.5*ns->cc.v[dim2]*delta[dim2]*ns->cc.du2dx2[dim2]/ns->cc.dudxl[dim2];
            } 
        } else {
            if (fabs(ns->cc.dudxr[dim2]) > 1.0e-8) {
                value -= 0.5*ns->cc.v[dim2]*delta[dim2]*ns->cc.du2dx2[dim2]/ns->cc.dudxr[dim2];
            } 
        }
    }
    return value;
}

// Convective term: second order scheme
real higflow_convective_term_second_order(higflow_solver *ns, Point delta, int dim) {
    // Set the second order convective term
    real value;
    switch (ns->contr.secondconvecdiscrtype) {
        // Modified Coefficient Upwind scheme
        case 0:
           value = higflow_convective_term_mcupwind(ns, delta, dim);
           break;
        // Cubista scheme
        case 1:
           value = higflow_convective_term_cubista(ns, delta, dim);
           break;
        // Quick scheme
        case 2:
           value = higflow_convective_term_quick(ns, delta, dim);
           break;
    }
    return value;
}

// Convective term: first order scheme
real higflow_convective_term_first_order(higflow_solver *ns, Point delta, int dim) {
    // Set the convective term
    real value = 0.0;
    for (int dim2 = 0; dim2 < DIM; dim2++) {
        value += 0.5*(ns->cc.v[dim2]+fabs(ns->cc.v[dim2]))*ns->cc.dudxl[dim2];
        value += 0.5*(ns->cc.v[dim2]-fabs(ns->cc.v[dim2]))*ns->cc.dudxr[dim2];
    }
    return value;
}

// Convective term: central scheme
real higflow_convective_term_central(higflow_solver *ns, Point delta, int dim) {
    // Set the convective term
    real value = 0.0;
    for (int dim2 = 0; dim2 < DIM; dim2++) {
        value += ns->cc.v[dim2]*ns->cc.dudxc[dim2];
    }
    return value;
}

// Convective term contribution for the Navier-Stokes equation
real higflow_convective_term(higflow_solver *ns, Point delta, int dim) {
    // Set the convective term according choice
    real value;
    switch (ns->cc.convec_type) {
        // Cental differening scheme
        case 0:
           value = higflow_convective_term_central(ns, delta, dim);
           break;
        // First order upwind scheme
        case 1:
           value = higflow_convective_term_first_order(ns, delta, dim);
           break;
        // Second order upwind scheme
        case 2:
           value = higflow_convective_term_second_order(ns, delta, dim);
           break;
    }
    return value;
}

// Difusive term contribution for the Navier-Stokes equation
real higflow_difusive_term(higflow_solver *ns, Point delta) {
    real value = 0.0;
    switch (ns->contr.flowtype) {
        // Newtonian
        case 0:
           for (int dim2 = 0; dim2 < DIM; dim2++) {
               value += ns->cc.du2dx2[dim2];
           }
           value /= ns->par.Re;
           break;
        // Generalized Newtonian
        case 1:
           for (int dim2 = 0; dim2 < DIM; dim2++) {
               value += ns->cc.du2dx2[dim2];
           }
           value /= ns->par.Re;
           break;
        // Multiphase
        case 2:
           for (int dim2 = 0; dim2 < DIM; dim2++) {
               value += ns->cc.du2dx2[dim2];
           }
           value /= ns->par.Re;
           value /= ns->cc.dens;
           break;
        // Viscoelastic
        case 3:
           for (int dim2 = 0; dim2 < DIM; dim2++) {
               value += ns->cc.du2dx2[dim2];
           }
           value /= ns->par.Re;
           break;
       // Viscoelastic
        case 4:
           for (int dim2 = 0; dim2 < DIM; dim2++) {
               value += ns->cc.du2dx2[dim2];
           }
           value /= ns->par.Re;
           break;    
    }
    return value;
}

// Cell electroosmotic source term contribution for the Navier-Stokes equation
real higflow_electroosmotic_source_term(higflow_solver *ns) {
    // Set the source term
    real value = ns->cc.Feo;
    return value;
}

// Cell electroosmotic diffusive ionic term contribution for the Navier-Stokes equation
real higflow_diffusive_ionic_term(higflow_solver *ns) {
    // Set the diffusive ionic term
    real value = ns->cc.d2ndx2/ns->ed.eo.par.Pe;
    return value;
}

// Cell electroosmotic potential ionic term contribution for the Navier-Stokes equation
real higflow_potential_ionic_term(higflow_solver *ns) {
    // Set the potential ionic term
    real value;
    value   = ns->cc.dndx*(ns->cc.dphidx + ns->cc.dpsidx) + ns->cc.ncell*(ns->cc.d2psidx2 + ns->cc.d2phidx2);
    value   *= ns->ed.eo.par.alpha/ns->ed.eo.par.Pe;
    return value;
}
