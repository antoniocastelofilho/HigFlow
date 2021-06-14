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
 
// Print the velocity
void higflow_print_velocity(higflow_solver *ns, FILE *data, int dimprint, real pprint);

// Print the VTK file for vilusalize
void higflow_print_vtk(higflow_solver *ns, int rank); 

// Print the VTK 2D file for vilusalize
void higflow_print_vtk2D(higflow_solver *ns, int rank); 

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

// Loading the viscoelastic parameters
void higflow_load_viscoelastic_parameters(higflow_solver *ns, int myrank); 

// Loading the electroosmotic parameters
void higflow_load_electroosmotic_parameters(higflow_solver *ns, int myrank); 

// Saving the viscoelastic parameters
void higflow_save_viscoelastic_parameters(higflow_solver *ns, int myrank); 
// Saving the electroosmotic parameters
void higflow_save_electroosmotic_parameters(higflow_solver *ns, int myrank);

// Loading the viscoelastic controllers
void higflow_load_viscoelastic_controllers(higflow_solver *ns, int myrank); 

// Loading the electroosmotic controllers
void higflow_load_electroosmotic_controllers(higflow_solver *ns, int myrank);

// Saving the viscoelastic controllers
void higflow_save_viscoelastic_controllers(higflow_solver *ns, int myrank); 

// Saving the electroosmotic controllers
void higflow_save_electroosmotic_controllers(higflow_solver *ns, int myrank);

// Print the Polymeric Tensor
void higflow_print_polymeric_tensor(higflow_solver *ns, FILE *data, int dimprint, real pprint);

// Loading the viscoelastic parameters
void higflow_load_viscoelastic_integral_parameters(higflow_solver *ns, int myrank); 

// Saving the viscoelastic integral parameters
void higflow_save_viscoelastic_integral_parameters(higflow_solver *ns, int myrank); 

// Loading the viscoelastic integral controllers
void higflow_load_viscoelastic_integral_controllers(higflow_solver *ns, int myrank); 

// Saving the viscoelastic integral controllers
void higflow_save_viscoelastic_integral_controllers(higflow_solver *ns, int myrank); 

#endif
