// *******************************************************************
// *******************************************************************
//  Hig-Flow Solver Boundary Condition - version 10/11/2016
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_BC
#define HIG_FLOW_BC

#include "hig-flow-kernel.h"

// Make the boundary condition
sim_boundary *higflow_make_bc(hig_cell *bcg, bc_type type, int id, bc_valuetype valuetype); 

// Creating and setting the boundary condition for the pressure
void higflow_set_boundary_condition_for_pressure(higflow_solver *ns, int numbcs, int id[numbcs], char bcfilenames[numbcs][1024], bc_type pbctypes[], bc_valuetype pbcvaluetype[]); 

// Creating and setting the boundary condition for the velocity
void higflow_set_boundary_condition_for_velocities(higflow_solver *ns, int numbcs, int id[numbcs], char bcfilenames[numbcs][1024], bc_type bctypes[DIM][numbcs], bc_valuetype bcvaluetype[DIM][numbcs]); 

// Creating and setting the boundary condition for the electro-osmotic source term
void higflow_set_boundary_condition_for_electroosmotic_source_term(higflow_solver *ns, int numbcs, int id[numbcs], char bcfilenames[numbcs][1024], bc_type bctypes[DIM][numbcs], bc_valuetype bcvaluetype[DIM][numbcs]); 

// Creating and setting the boundary condition for the electro-osmotic phi
void higflow_set_boundary_condition_for_electroosmotic_phi(higflow_solver *ns, int numbcs, int id[numbcs], char bcfilenames[numbcs][1024], bc_type pbctypes[], bc_valuetype pbcvaluetype[]); 

// Creating and setting the boundary condition for the electro-osmotic psi
void higflow_set_boundary_condition_for_electroosmotic_psi(higflow_solver *ns, int numbcs, int id[numbcs], char bcfilenames[numbcs][1024], bc_type pbctypes[], bc_valuetype pbcvaluetype[]); 

// Creating and setting the boundary condition for the electro-osmotic nplus
void higflow_set_boundary_condition_for_electroosmotic_nplus(higflow_solver *ns, int numbcs, int id[numbcs], char bcfilenames[numbcs][1024], bc_type pbctypes[], bc_valuetype pbcvaluetype[]); 

// Creating and setting the boundary condition for the electro-osmotic nminus
void higflow_set_boundary_condition_for_electroosmotic_nminus(higflow_solver *ns, int numbcs, int id[numbcs], char bcfilenames[numbcs][1024], bc_type pbctypes[], bc_valuetype pbcvaluetype[]); 

// Navier-Stokes initialize the domain and boudaries
void higflow_initialize_boundaries(higflow_solver *ns); 

// Navier-Stokes initialize the domain and boudaries
void higflow_initialize_boundaries_and_domain(higflow_solver *ns);

#endif
