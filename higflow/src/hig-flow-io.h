// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver IO - version 10/11/2016
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_IO
#define HIG_FLOW_IO

#include "hig-flow-kernel.h"
#include "hig-flow-eval.h"

// *******************************************************************
// Navier-Stokes Print for Visualize
// *******************************************************************

// Print the electro-osmotic proporties
void higflow_print_electroosmotic_properties(higflow_solver *ns, FILE *data, int dimprint, real pprint);
 
// Print the Polymeric Tensor
void higflow_print_polymeric_tensor(higflow_solver *ns, FILE *data, int dimprint, real pprint);

// Print the velocity
void higflow_print_velocity(higflow_solver *ns, FILE *data, int dimprint, real pprint);

// Print the VTK file for vilusalize
void higflow_print_vtk(higflow_solver *ns, int rank); 

// Print the kinetic energy
void higflow_kinetic_energy(higflow_solver *ns, FILE *data);

// Print the VTK 2D file for vilusalize
void higflow_print_vtk2D(higflow_solver *ns, int rank); 

// Print a single vtk file for vilusalize 2d
void higflow_print_vtk2D_parallel_single(higflow_solver *ns, int rank, int nprocs);

// Print the VTK 3D file for vilusalize
void higflow_print_vtk3D(higflow_solver *ns, int rank); 

// Print the VTK file for vilusalize 3D
void higflow_print_vtk3D_viscoelastic(higflow_solver *ns, int rank); 

// *******************************************************************
// Navier-Stokes Load and Save
// *******************************************************************

// Loading the data files
void higflow_load_data_files(int argc, char *argv[], higflow_solver *ns); 

// Loading the controllers and Parameters
void higflow_load_controllers_and_parameters_yaml(higflow_solver *ns, int myrank);

// Loading the controllers
void higflow_load_controllers(higflow_solver *ns, int myrank); 

// Saving the controllers
void higflow_save_controllers(higflow_solver *ns, int myrank); 

// Loading the parameters
void higflow_load_parameters(higflow_solver *ns, int myrank); 

// Saving the parameters
void higflow_save_parameters(higflow_solver *ns, int myrank); 

// Loading the properties
void higflow_load_properties(higflow_solver *ns, int myrank, int ntasks);

// Saving the properties
void higflow_save_properties(higflow_solver *ns, int myrank, int ntasks); 

// Save the domain amr info
void higflow_save_domain(higflow_solver *ns, int myrank, int ntasks);

// save the boundary amr info
void higflow_save_boundaries(higflow_solver *ns, int myrank, int ntasks);

// Save the boundary amr info for electroosmotic
void higflow_save_boundaries_electroosmotic(higflow_solver *ns, int myrank, int ntasks);

// Loading the multiphase parameters
void higflow_load_multiphase_parameters(higflow_solver *ns, int myrank);

// Saving the multiphase parameters
void higflow_save_multiphase_parameters(higflow_solver *ns, int myrank);

// Loading the multiphase controllers
void higflow_load_multiphase_controllers(higflow_solver *ns, int myrank);

// Saving the multiphase controllers
void higflow_save_multiphase_controllers(higflow_solver *ns, int myrank);

// Loading the multiphase viscoelastic parameters
void higflow_load_multiphase_viscoelastic_parameters(higflow_solver *ns, int myrank);

// Saving the multiphase viscoelastic parameters
void higflow_save_multiphase_viscoelastic_parameters(higflow_solver *ns, int myrank);

// Loading the multiphase viscoelastic controllers
void higflow_load_multiphase_viscoelastic_controllers(higflow_solver *ns, int myrank);

// Saving the multiphase viscoelastic controllers
void higflow_save_multiphase_viscoelastic_controllers(higflow_solver *ns, int myrank);

// Loading the viscoelastic parameters
void higflow_load_viscoelastic_parameters(higflow_solver *ns, int myrank);

// Saving the viscoelastic parameters
void higflow_save_viscoelastic_parameters(higflow_solver *ns, int myrank); 

// Loading the viscoelastic controllers
void higflow_load_viscoelastic_controllers(higflow_solver *ns, int myrank);

// Saving the viscoelastic controllers
void higflow_save_viscoelastic_controllers(higflow_solver *ns, int myrank);

// Loading the electroosmotic parameters
void higflow_load_electroosmotic_parameters(higflow_solver *ns, int myrank); 

// Saving the electroosmotic parameters
void higflow_save_electroosmotic_parameters(higflow_solver *ns, int myrank);

// Loading the electroosmotic controllers
void higflow_load_electroosmotic_controllers(higflow_solver *ns, int myrank);

// Saving the electroosmotic controllers
void higflow_save_electroosmotic_controllers(higflow_solver *ns, int myrank);

// Loading the viscoelastic parameters
void higflow_load_viscoelastic_integral_parameters(higflow_solver *ns, int myrank); 

// Saving the viscoelastic integral parameters
void higflow_save_viscoelastic_integral_parameters(higflow_solver *ns, int myrank); 

// Loading the viscoelastic integral controllers
void higflow_load_viscoelastic_integral_controllers(higflow_solver *ns, int myrank); 

// Saving the viscoelastic integral controllers
void higflow_save_viscoelastic_integral_controllers(higflow_solver *ns, int myrank); 

#endif
