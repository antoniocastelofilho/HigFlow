// Registering boundary conditions on the solver: one setter per field, each taking
// the whole set at once.
//
// EVERY SETTER TAKES FOUR PARALLEL ARRAYS -- id[], bcfilenames[][1024], the type
// array and the value-type array -- all indexed together and all trusted to have
// numbcs entries.  Nothing checks that, and a short array is read past its end.
//
// A FACE WITH NO ENTRY IS NOT AN ERROR AND IS NOT REPORTED.  There is no
// completeness check against the mesh: an id listed in the YAML whose .amr patch is
// missing, or a wall nobody registered, simply leaves that stretch of the domain
// with no condition, and the stencil closure falls through to an interior formula.
// The run then produces a plausible field with a hole in the wall.  This has
// happened here, and it was found by a prediction failing, not by a diagnostic.
//
// higflow_set_bc_refine_hook() lets the caller refine a boundary patch as it is
// read, which is how a BC mesh is matched to a refined region of the domain.

// *******************************************************************
// *******************************************************************
//  Hig-Flow Solver Boundary Condition - version 10/11/2016
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_BC
#define HIG_FLOW_BC

#include "hig-flow-kernel.h"
#include "hig-flow-eval.h"

// Make the boundary condition
sim_boundary *higflow_make_bc(hig_cell *bcg, bc_type type, int id, bc_valuetype valuetype);

// Register a callback invoked on each boundary higtree immediately after it is
// read from disk, before the sim_boundary is created.  Use this to refine the
// higtree in place to match the adjacent internal mesh.  Pass NULL to disable.
void higflow_set_bc_refine_hook(void (*hook)(hig_cell *bc_root, int bc_id));

// Creating and setting the boundary condition for the pressure
void higflow_set_boundary_condition_for_pressure(higflow_solver *ns, int numbcs, int id[], char bcfilenames[][1024], bc_type pbctypes[], bc_valuetype pbcvaluetype[]); 

// Creating and setting the boundary condition for the velocity
void higflow_set_boundary_condition_for_velocities(higflow_solver *ns, int
    numbcs, int id[], char bcfilenames[][1024], bc_type
    bctypes[][DIM], bc_valuetype bcvaluetype[][DIM]); 

// Creating and setting the boundary condition for the electro-osmotic source term
void higflow_set_boundary_condition_for_electroosmotic_source_term(higflow_solver *ns, int numbcs, int id[], char bcfilenames[][1024], bc_type bctypes[][DIM], bc_valuetype bcvaluetype[][DIM]); 

// Creating and setting the boundary condition for the electro-osmotic phi
void higflow_set_boundary_condition_for_electroosmotic_phi(higflow_solver *ns, int numbcs, int id[], char bcfilenames[][1024], bc_type pbctypes[], bc_valuetype pbcvaluetype[]); 

// Creating and setting the boundary condition for the electro-osmotic psi
void higflow_set_boundary_condition_for_electroosmotic_psi(higflow_solver *ns, int numbcs, int id[], char bcfilenames[][1024], bc_type pbctypes[], bc_valuetype pbcvaluetype[]); 

// Creating and setting the boundary condition for the electro-osmotic nplus
void higflow_set_boundary_condition_for_electroosmotic_nplus(higflow_solver *ns, int numbcs, int id[], char bcfilenames[][1024], bc_type pbctypes[], bc_valuetype pbcvaluetype[]); 

// Creating and setting the boundary condition for the electro-osmotic nminus
void higflow_set_boundary_condition_for_electroosmotic_nminus(higflow_solver *ns, int numbcs, int id[], char bcfilenames[][1024], bc_type pbctypes[], bc_valuetype pbcvaluetype[]); 

// Creating and setting the boundary condition for volume fraction and other multiphase properties
void higflow_set_boundary_condition_for_mult(higflow_solver *ns, int numbcs, int id[], char bcfilenames[][1024], bc_type bctypes[], bc_valuetype bcvaluetype[]);

// Creating and setting the boundary condition for extra domains domain (for viscoelastic tensors and such)
void higflow_set_boundary_condition_for_tensors(higflow_solver *ns, int numbcs, int id[], char bcfilenames[][1024], bc_type bctypes[], bc_valuetype bcvaluetype[]);

// Navier-Stokes initialize the domain and boudaries
void higflow_initialize_boundaries(higflow_solver *ns); 

// Creating and setting the boundary condition for the density number nA (viscoelastic flows with shear-banding)
void higflow_set_boundary_condition_for_viscoelastic_shear_banding_nA(higflow_solver *ns, int numbcs, int id[], char bcfilenames[][1024], bc_type nAbctypes[], bc_valuetype nAbcvaluetype[]);

// Creating and setting the boundary condition for the density number nB (viscoelastic flows with shear-banding)
void higflow_set_boundary_condition_for_viscoelastic_shear_banding_nB(higflow_solver *ns, int numbcs, int id[], char bcfilenames[][1024], bc_type nBbctypes[], bc_valuetype nBbcvaluetype[]);

// Navier-Stokes initialize the domain and boudaries
void higflow_initialize_boundaries_yaml(higflow_solver *ns);

#endif
