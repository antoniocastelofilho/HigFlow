// *******************************************************************
// *******************************************************************
//  Example for HiG-Flow Solver - version 10/11/206
// *******************************************************************
// *******************************************************************

//#include "../src/hig-flow-step-viscoelastic.h"
#include "../src/hig-flow-step-generalized-newtonian.h"
#include "../src/hig-flow-io.h"
#include "../src/hig-flow-ic.h"
#include "../src/hig-flow-bc.h"
#define DEBUG
#include "Debug-c.h"
#include <stdio.h>
#include <stdlib.h>

// Define the time constants count
DECL_CLOCK(total)
DECL_CLOCK(firstiter)

// *******************************************************************
// Extern functions for the Navier-Stokes program
// *******************************************************************

// Value of the pressure
real get_pressure(Point center, real t); 

// Value of the velocity
real get_velocity(Point center, int dim, real t); 

// Value of the cell source term
real get_source_term(Point center, real t); 

// Value of the facet source term
real get_facet_source_term(Point center, int dim, real t); 

// Value of the pressure at boundary
real get_boundary_pressure(int id, Point center, real t); 

// Value of the velocity at boundary
real get_boundary_velocity(int id, Point center, int dim, real t); 

// Value of the cell source term at boundary
real get_boundary_source_term(int id, Point center, real t); 

// Value of the facet source term at boundary
real get_boundary_facet_source_term(int id, Point center, int dim, real t); 

// Value of the Tensor
real get_tensor(Point center, int i, int j, real t);

// Value of the Kernel
real get_kernel(int dim, real lambda, real tol);

// Value of the Kernel inverse
real get_kernel_inverse(int dim, real lambda, real tol);

// Value of the Kernel Jacobian
real get_kernel_jacobian(int dim, real lambda, real tol);

// Print the velocity
void print_velocity (higflow_solver *ns, int myrank, int dim, int dimprint1, real pprint1, int dimprint2, real pprint2); 

// Print the Polymeric Tensor
void print_tensor(higflow_solver *ns, int myrank, int i, int j, int dimprint1, real pprint1, int dimprint2, real pprint2); 
