// *******************************************************************
// *******************************************************************
//  Example for HiG-Flow Solver - version 03/2023
// *******************************************************************
// *******************************************************************

#ifndef NS_EXAMPLE_2D_H
#define NS_EXAMPLE_2D_H

#include "../src/hig-flow-step-multiphase-electroosmotic-viscoelastic.h"
#include "../src/hig-flow-step-electroosmotic-viscoelastic.h"
#include "../src/hig-flow-step-viscoelastic-integral.h"
#include "../src/hig-flow-io.h"
#include "../src/hig-flow-ic.h"
#include "../src/hig-flow-bc.h"
#include "utilities-2d.h"
#define DEBUG
#include "Debug-c.h"
#include <stdio.h>
#include <stdlib.h>

#define pi 3.1415926535897932384626434

// parameters for access in external functions
physical_parameters* p_par;
flow_type flowtype; 
flow_type flowtype0; flow_type flowtype1;
int viscoelastic_either;
int eoflow;
int eoflow0; int eoflow1;
int eoflow_either;
bc_type u_inlet;
bc_type psi_inlet;
sim_domain **sdmult_ptr;
distributed_property **dpfracvol_ptr;
sim_stencil **stnmult_ptr;

sim_domain **sdp_ptr;
distributed_property **dpp_ptr;
sim_stencil **stn_ptr;
sim_facet_domain **sfdv_ptr;
distributed_property **dpvstar_ptr;

sim_facet_domain **sfdFeoy_ptr;
distributed_property **dpFeoy_ptr;
sim_stencil **stnFeoy_ptr;

// induced potentials at the walls
real psi_up = -1.0, psi_down = 1.0; 


visc_model_type visc_model;
// checking for the positivity of the conformation tensor
real max_neg_lambda;
int num_neg_lambda;

// Define the time constants count
DECL_CLOCK(total)
DECL_CLOCK(firstiter)
DECL_CLOCK(currentiter)
DECL_CLOCK(iter_total)

real func (Point p);
real get_kernel(int dim, real lambda, real tol);

#endif