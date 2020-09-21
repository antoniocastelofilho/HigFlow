// *******************************************************************
// *******************************************************************
//  Example for HiG-Flow Solver - version 10/11/2016
// *******************************************************************
// *******************************************************************

#include "../src/hig-flow-step-multiphase.h"
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

// Value of the density 
real get_density0(Point center, real t);
real get_density1(Point center, real t);

// Value of the density 
real get_viscosity0(Point center, real t);
real get_viscosity1(Point center, real t);

real func (Point p); 

int square_case (real f0, real f1, real f10, real f3); 


void Intersec(Point p0, Point p1, real f0, real f1, Point p); 

// Volume fraction
real get_fracvol(Point center, Point delta, real t); 
	
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

real analyticP(Point pt, real t);
real analyticU(int dim, Point pt, real t);
void compute_error_norm_velocity(int dim, higflow_solver *ns, real *norm2, real *norminf);
void compute_error_norm_pressure(higflow_solver *ns, real *norm2, real *norminf);
void compute_and_print_error_norm_velocity(higflow_solver *ns, int myrank);
void compute_and_print_error_norm_pressure(higflow_solver *ns, int myrank);
