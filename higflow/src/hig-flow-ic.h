// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Initial Condition - version 10/11/2016
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_IC
#define HIG_FLOW_IC

#include "hig-flow-kernel.h"
#include "hig-flow-eval.h"

// *******************************************************************
// Navier-Stokes Initialize Domain and Initial Conditions
// *******************************************************************

// Navier-Stokes initialize the domain
void higflow_initialize_domain(higflow_solver *ns, int ntasks, int myrank, int order); 

// *******************************************************************
// Navier-Stokes Initialize Properties
// *******************************************************************

// Initialize the pressure
void higflow_initialize_pressure(higflow_solver *ns); 

// Initialize the viscosity
void higflow_initialize_viscosity_gn(higflow_solver *ns); 

// Initialize the viscosity
void higflow_initialize_viscosity_mult(higflow_solver *ns); 

// Initialize the volume fraction
void higflow_initialize_fracvol(higflow_solver *ns); 

// Initialize the density
void higflow_initialize_density(higflow_solver *ns); 

// Initialize the cell source term
void higflow_initialize_cell_source_term(higflow_solver *ns); 

// Initialize the viscoelastic Tensor
void higflow_initialize_viscoelastic_tensor(higflow_solver *ns);

// Initialize the viscoelastic integral Tensor
void higflow_initialize_viscoelastic_integral_tensor(higflow_solver *ns); 

// Initialize the viscoelastic integral Tensor
void higflow_initialize_viscoelastic_integral_finger_tensor(higflow_solver *ns); 

// Initialize the electro-osmotic phi
void higflow_initialize_electroosmotic_phi(higflow_solver *ns); 

// Initialize the electro-osmotic psi
void higflow_initialize_electroosmotic_psi(higflow_solver *ns); 

// Initialize the electro-osmotic nplus
void higflow_initialize_electroosmotic_nplus(higflow_solver *ns); 

// Initialize the electro-osmotic nminus
void higflow_initialize_electroosmotic_nminus(higflow_solver *ns);

// Initialize the velocities
void higflow_initialize_velocity(higflow_solver *ns); 

// Initialize the facet source term
void higflow_initialize_facet_source_term(higflow_solver *ns); 

// Initialize distributed properties
void higflow_initialize_distributed_properties(higflow_solver *ns); 

// Load string from file
void __higflow_readstring(char s[], int max, FILE *file);

#endif
