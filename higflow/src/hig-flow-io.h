// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver IO - version 10/11/2016
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_IO
#define HIG_FLOW_IO

#include "hig-flow-kernel.h"
#include "hig-flow-eval.h"
#include <string.h>
#include <libfyaml.h>

// *******************************************************************
// Navier-Stokes Print for Visualize
// *******************************************************************

// Print the electro-osmotic proporties
void higflow_print_electroosmotic_properties(higflow_solver *ns, FILE *data, int dimprint, real pprint);
 
// Print the Polymeric Tensor
void higflow_print_polymeric_tensor(higflow_solver *ns, FILE *data, int dimprint, real pprint);

// Print the velocity
void higflow_print_velocity(higflow_solver *ns, FILE *data, int dimprint, real pprint);

// Print the kinetic energy
void higflow_kinetic_energy(higflow_solver *ns, FILE *data);

// Print the VTK file for vilusalize
void higflow_print_vtk(higflow_solver *ns, int rank); 

// Print the VTK 2D file for vilusalize
void higflow_print_vtk2D(higflow_solver *ns, int rank); 

void higflow_print_vtk2D_multiphase(higflow_solver *ns, int rank);

// Print a single vtk file for vilusalize 2d
void higflow_print_vtk2D_parallel_single(higflow_solver *ns, int rank, int nprocs);

void higflow_print_vtk2D_multiphase_parallel_single(higflow_solver *ns, int rank, int nprocs);

// Print the VTK 3D file for vilusalize
void higflow_print_vtk3D(higflow_solver *ns, int rank); 

// Print the VTK file for vilusalize 3D
void higflow_print_vtk3D_viscoelastic(higflow_solver *ns, int rank); 

// Print the VTK file for vilusalize 3D viscoelastic flows with variable viscosity
void higflow_print_vtk3D_viscoelastic_variable_viscosity(higflow_solver *ns, int rank);

// Print the VTK file for vilusalize 3D viscoelastic flows with shear-banding
void higflow_print_vtk3D_viscoelastic_shear_banding(higflow_solver *ns, int rank);

// Print the VTK file for vilusalize 3D elastoviscoplastic flows 
void higflow_print_vtk3D_elastoviscoplastic(higflow_solver *ns, int rank);

// Print the VTK file for vilusalize 3D
void higflow_print_vtk3D_shear_thickening_suspensions(higflow_solver *ns, int rank);
// *******************************************************************
// Navier-Stokes Load and Save
// *******************************************************************

// Loading the data files
void higflow_load_data_file_names(int argc, char *argv[], higflow_solver *ns); 

// Save the domain amr info
void higflow_save_domain(higflow_solver *ns, int myrank, int ntasks);

// save the boundary amr info
void higflow_save_boundaries(higflow_solver *ns, int myrank, int ntasks);

// Save the boundary amr info for electroosmotic
void higflow_save_boundaries_electroosmotic(higflow_solver *ns, int myrank, int ntasks);

// Save all the boundary amr info
void higflow_save_all_boundaries(higflow_solver *ns, int myrank, int ntasks);

// Save the yaml domain amr info
void higflow_save_domain_yaml(higflow_solver *ns, int myrank, int ntasks);

// Save the yaml boundary amr info
void higflow_save_all_boundaries_yaml(higflow_solver *ns, int myrank, int ntasks);

// Loading the properties
void higflow_load_properties(higflow_solver *ns, int myrank, int ntasks);

// Saving the properties
void higflow_save_properties(higflow_solver *ns, int myrank, int ntasks); 

/////////////// Controllers and Parameters //////////////////////////

// Loading the controllers and Parameters from yaml
void higflow_load_all_controllers_and_parameters_yaml(higflow_solver *ns, int myrank);

// Saving the controllers and Parameters to yaml
void higflow_save_all_controllers_and_parameters_yaml(higflow_solver* ns, int myrank);

// Loading all the controllers
void higflow_load_all_controllers(higflow_solver *ns, int myrank);

// Saving all the controllers
void higflow_save_all_controllers(higflow_solver *ns, int myrank);

// Loading all the parameters
void higflow_load_all_parameters(higflow_solver *ns, int myrank);

// Saving all the parameters
void higflow_save_all_parameters(higflow_solver *ns, int myrank);

// Loading the controllers
void higflow_load_controllers(higflow_solver *ns, int myrank); 

// Saving the controllers
void higflow_save_controllers(higflow_solver *ns, int myrank); 

// Loading the parameters
void higflow_load_parameters(higflow_solver *ns, int myrank); 

// Saving the parameters
void higflow_save_parameters(higflow_solver *ns, int myrank); 

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

// Loading the multiphase electroosmotic parameters
void higflow_load_multiphase_electroosmotic_parameters(higflow_solver *ns, int myrank);

// Saving the multiphase electroosmotic parameters
void higflow_save_multiphase_electroosmotic_parameters(higflow_solver *ns, int myrank);

// Loading the multiphase electroosmotic controllers
void higflow_load_multiphase_electroosmotic_controllers(higflow_solver *ns, int myrank);

// Saving the multiphase electroosmotic controllers
void higflow_save_multiphase_electroosmotic_controllers(higflow_solver *ns, int myrank);

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

// Loading the parameters of viscoelastic flows with variable viscosity 
void higflow_load_viscoelastic_variable_viscosity_parameters(higflow_solver *ns, int myrank);

// Saving the parameters of the viscoelastic flows with variable viscosity 
void higflow_save_viscoelastic_variable_viscosity_parameters(higflow_solver *ns, int myrank);

// Loading the controlles for viscoelastic flows with variable viscosity
void higflow_load_viscoelastic_variable_viscosity_controllers(higflow_solver *ns, int myrank);

// Saving the controllers for viscoelastic flows with variable viscosity
void higflow_save_viscoelastic_variable_viscosity_controllers(higflow_solver *ns, int myrank);

// Loading the parameters of viscoelastic flows with shear-banding
void higflow_load_viscoelastic_shear_banding_parameters(higflow_solver *ns, int myrank);

// Saving the parameters of the viscoelastic flows with shear-banding
void higflow_save_viscoelastic_shear_banding_parameters(higflow_solver *ns, int myrank);

// Loading the controllers for viscoelastic flows with shear-banding
void higflow_load_viscoelastic_shear_banding_controllers(higflow_solver *ns, int myrank);

// Saving the controllers for viscoelastic flows with shear_banding
void higflow_save_viscoelastic_shear_banding_controllers(higflow_solver *ns, int myrank);

// Loading the elastoviscoplastic parameters
void higflow_load_elastoviscoplastic_parameters(higflow_solver *ns, int myrank);

// Saving the elastoviscoplastic parameters
void higflow_save_elastoviscoplastic_parameters(higflow_solver *ns, int myrank);

// Loading the elastoviscoplastic controllers
void higflow_load_elastoviscoplastic_controllers(higflow_solver *ns, int myrank);

// Saving the elastoviscoplastic controllers
void higflow_save_elastoviscoplastic_controllers(higflow_solver *ns, int myrank);

// Loading the shear-thickening suspension parameters
void higflow_load_shear_thickening_suspension_parameters(higflow_solver *ns, int myrank);

// Saving the shear-thickening suspension parameters
void higflow_save_shear_thickening_suspension_parameters(higflow_solver *ns, int myrank);

// Loading the shear-thickening suspension controllers
void higflow_load_shear_thickening_suspension_controllers(higflow_solver *ns, int myrank);

// Saving the shear-thickening suspension controllers
void higflow_save_shear_thickening_suspension_controllers(higflow_solver *ns, int myrank);
#endif
