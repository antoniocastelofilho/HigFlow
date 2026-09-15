// *******************************************************************
//  HiG-Flow Solver Step Electro-osmotic - version 12/05/2024
// *******************************************************************

#include "hig-flow-step-multiphase-electroosmotic-viscoelastic.h"

void higflow_solver_step_multiphase_electroosmotic_viscoelastic(higflow_solver *ns) {
    // Calculate the deformation rate tensor
    higflow_compute_velocity_derivative_tensor(ns);
    // Constitutive Equation Step for viscoelastic flow
    switch (ns->ed.mult.ve.contr.discrtype) {
        case EXPLICIT:
            higflow_explicit_euler_constitutive_equation_multiphase_viscoelastic(ns);
            break;
        case IMPLICIT:
            higflow_implicit_euler_constitutive_equation_multiphase_viscoelastic(ns);
            break;
    }
    // Calculate the elastic tensor to be used in the momentum equation
    higflow_compute_polymeric_tensor_multiphase_viscoelastic(ns);

    // Calculate the multiphase electroosmotic flow
    higflow_solver_step_multiphase_electroosmotic(ns);
}
