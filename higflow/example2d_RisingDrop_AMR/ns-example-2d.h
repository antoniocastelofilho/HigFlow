#include "../src/hig-flow-step-multiphase.h"
#include "../src/hig-flow-io.h"
#include "../src/hig-flow-ic.h"
#include "../src/hig-flow-bc.h"
#define DEBUG
#include "Debug-c.h"
#include <stdio.h>
#include <stdlib.h>

// Adaptive mesh refinement around the interface.
// NOTE: AMR currently only works serial (NP=1) — the parallel rebuild does
// not partition the adapted tree. Run with NP=1.
#define ADAPT_ENABLED 1
#define ADAPT_FREQ 5

DECL_CLOCK(total)
DECL_CLOCK(firstiter)

real get_pressure(Point center, real t);
real get_velocity(Point center, int dim, real t);
real get_source_term(Point center, real t);
real get_density0(Point center, real t);
real get_density1(Point center, real t);
real get_viscosity0(Point center, real t);
real get_viscosity1(Point center, real t);
real func(Point p);
int square_case(real f0, real f1, real f10, real f3);
void Intersec(Point p0, Point p1, real f0, real f1, Point p);
real get_fracvol(Point center, Point delta, real t);
real get_facet_source_term(Point center, int dim, real t);
real get_boundary_pressure(int id, Point center, real t);
real get_boundary_velocity(int id, Point center, int dim, real t);
real get_boundary_source_term(int id, Point center, real t);
real get_boundary_facet_source_term(int id, Point center, int dim, real t);
hig_cell *higflow_make_adapted_tree_params(higflow_solver *ns, real *thresholds);
