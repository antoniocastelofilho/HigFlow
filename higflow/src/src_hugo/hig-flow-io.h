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
void higflow_load_data_files(int argc, char *argv[], higflow_solver *ns); 

// Loading the controllers
void higflow_load_controllers(higflow_solver *ns, int myrank); 

// Saving the controllers
void higflow_save_controllers(higflow_solver *ns, int myrank); 

// Loading the parameters
void higflow_load_parameters(higflow_solver *ns, int myrank); 

// Saving the parameters
void higflow_save_parameters(higflow_solver *ns, int myrank); 

// Loading the properties
void higflow_load_properties(higflow_solver *ns); 

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

// Loading the parameters of viscoelastic flows with variable viscosity 
void higflow_load_viscoelastic_variable_viscosity_parameters(higflow_solver *ns, int myrank);

// Saving the parameters of the viscoelastic flows with variable viscosity 
void higflow_save_viscoelastic_variable_viscosity_parameters(higflow_solver *ns, int myrank);

// Loading the controlles for viscoelastic flows with variable viscosity
void higflow_load_viscoelastic_variable_viscosity_controllers(higflow_solver *ns, int myrank);

// Saving the controllers for viscoelastic flows with variable viscosity
void higflow_save_viscoelastic_variable_viscosity_controllers(higflow_solver *ns, int myrank);

// Loading the parameters of viscoelastic flows with shear-banding
void higflow_load_viscoelastic_variable_shear_banding_parameters(higflow_solver *ns, int myrank);

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

// Loading the non-isothermal flow parameters
void higflow_load_non_isothermal_flow_parameters(higflow_solver *ns, int myrank);

// Saving the non-isothermal flow parameters
void higflow_save_non_isothermal_flow_parameters(higflow_solver *ns, int myrank);

// Loading the non-isothermal controllers
void higflow_load_non_isothermal_flow_controllers(higflow_solver *ns, int myrank);

// Saving the non-isothermal controllers
void higflow_save_non_isothermal_flow_controllers(higflow_solver *ns, int myrank);

#endif
