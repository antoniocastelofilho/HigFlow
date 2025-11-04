// *******************************************************************
//  HiG-Flow Solver Step Electro-osmotic - version 03/2024
// *******************************************************************

#include "hig-flow-step-electroosmotic-viscoelastic.h"

// One step of the Navier-Stokes the projection method for viscoelastic flow
void higflow_solver_step_electroosmotic_viscoelastic(higflow_solver *ns) {

    // Compute deformation rate tensor
    higflow_compute_velocity_derivative_tensor(ns);
    // Constitutive Equation Step for viscoelastic flow
    switch (ns->ed.ve.contr.discrtype) {
    case EXPLICIT:
        higflow_explicit_euler_constitutive_equation(ns);
        break;
    case IMPLICIT:
        higflow_implicit_euler_constitutive_equation(ns);
        break;
    }
    // Calculate the elastic tensor to be used in the momentum equation
    higflow_compute_polymeric_tensor(ns);

    // Calculate the electroosmotic flow
    higflow_solver_step_electroosmotic(ns);
}

