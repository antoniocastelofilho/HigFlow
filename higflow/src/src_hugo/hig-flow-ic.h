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

// Initialize the viscosity for viscoelastic flows
void higflow_initialize_viscosity_vevv(higflow_solver *ns);

// Initialize the structural parameter for viscoelastic flows
void higflow_initialize_structural_parameter(higflow_solver *ns);

// Initialize the Non-Newtonian Tensor for viscoelastic flows with variable viscosity
void higflow_initialize_viscoelastic_tensor_variable_viscosity(higflow_solver *ns);

// Initialize the shear-banding density number nA
void higflow_initialize_shear_banding_nA(higflow_solver *ns);

// Initialize the shear-banding density number nB
void higflow_initialize_shear_banding_nB(higflow_solver *ns);

// Initialize the concentration of specie A (cA)
void higflow_initialize_shear_banding_cA(higflow_solver *ns);

// Initialize the concentration of specie B (cB)
void higflow_initialize_shear_banding_cB(higflow_solver *ns);

// Initialize the Non-Newtonian Tensor for viscoelastic flows with shear-banding
void higflow_initialize_viscoelastic_tensor_shear_banding(higflow_solver *ns);

// Initialize the conformation tensor of specie A for viscoelastic flows with shear-banding
void higflow_initialize_conformation_tensor_A_shear_banding(higflow_solver *ns);

// Initialize the conformation tensor of specie B for viscoelastic flows with shear-banding
void higflow_initialize_conformation_tensor_B_shear_banding(higflow_solver *ns);

// Initialize the Non-Newtonian Tensor
void higflow_initialize_elastoviscoplastic_tensor(higflow_solver *ns);

// Initialize the Tensor S for shear-thickening suspension
void higflow_initialize_shear_thickening_suspension_tensor(higflow_solver *ns);

// Initialize the microstructure tensor A for shear-thickening suspension
void higflow_initialize_shear_thickening_suspension_microstructure_tensor(higflow_solver *ns);

// Initialize the volume fraction for shear thickening suspensions (with particle migration)
void higflow_initialize_volume_fraction(higflow_solver *ns);

// Initialize the energy source term for non-isothermal flows
void higflow_initialize_energy_source_term(higflow_solver *ns);

// Initialize the velocities
void higflow_initialize_velocity(higflow_solver *ns); 

// Initialize the facet source term
void higflow_initialize_facet_source_term(higflow_solver *ns); 

// Initialize distributed properties
void higflow_initialize_distributed_properties(higflow_solver *ns); 

// Load string from file
void __higflow_readstring(char s[], int max, FILE *file);

#endif
