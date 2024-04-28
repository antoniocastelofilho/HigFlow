// *******************************************************************
//  Example for HiG-Flow Solver - version 10/11/2016
// *******************************************************************

#include "ns-example-2d.h"

// *******************************************************************
// Extern functions for the Navier-Stokes program
// *******************************************************************

// Value of the Tensor
real get_tensor(Point center, int i, int j, real t) {
    real value = 0.0;
    return value; 
}

// Value of the Tensor
real get_boundary_tensor(int id, Point center, int i, int j, real t) {
    real value = 0.0;
    return value; 
}

// Value of the Kernel
real get_kernel(int dim, real lambda, real tol) {
    real value;
    //if (lambda < tol)
    //   value = log(tol);
    //else
    //   value = log(lambda);
    //if (lambda < tol)
    //   value = sqrt(tol);
    //else
    //   value = sqrt(lambda);
    value = lambda;
    return value; 
}

// Value of the Kernel inverse
real get_kernel_inverse(int dim, real lambda, real tol) {
    real value;
    //value = exp(lambda);
    //value = lambda*lambda;
    value = lambda;
    return value; 
}

// Value of the Kernel Jacobian
real get_kernel_jacobian(int dim, real lambda, real tol) {
    real value;
    //if (lambda < tol)
    //   value = 1.0/tol;
    //else
    //   value = 1.0/lambda;
    //if (lambda < tol)
    //   value = 0.5/sqrt(tol);
    //else
    //   value = 0.5/sqrt(lambda);
    value = 1.0;
    return value; 
}

real dpdx = 3.0;
real L = 8.0;
//higflow_solver *nsaux;

// Value of the pressure
real get_pressure(Point center, real t) {
    real value = 0.0;
    return value; 
}

// Value of the velocity
real get_velocity(Point center, int dim, real t) {
    real value = 0.0;
    return value; 
}

// Value of the cell source term
real get_source_term(Point center, real t) {
    real value = 0.0;
    return value; 
}

// Value of the facet source term
real get_facet_source_term(Point center, int dim, real t) {
    real value = 0.0;
    return value; 
}

// Value of the viscosity
real get_viscosity(Point center, real q, real t) {
    real value = 1.0;
    return value; 
}

// Value of the pressure at boundary
real get_boundary_pressure(int id, Point center, real t) {
    real value;
    switch (id) {
        case 0:
            value = 0.0;      
            break;
        case 1:
            value = 0.0;        
            break;
        case 2:
            value = 0.0;        
            break;
        case 3:
            value = 0.0;        
            break;
    }
    return value; 
}

// Value of the velocity at boundary
real get_boundary_velocity(int id, Point center, int dim, real t) {
    real value;
    switch (id) {
        case 0:
            switch (dim) {
                case 0: ;
                    //set max velocity = 8.0e-4 
                    //value = 3.2e-3*(-center[1]*center[1] + center[1]);
                    //value = 4.0*(-center[1]*center[1] + 0.25);
                    //set max velocity = 1.5 
                    //value = -4.0*center[1]*(center[1] - 1.0);
                    value = 1.5*(1.0 - center[1]*center[1]);
                    //value = 1.25*(1.0 - center[1]*center[1]*center[1]*center[1]);
                    //value = 1.0*(1.0 - fabs(center[1]));
                    //value = 2.0*(1.0 - sqrt(fabs(center[1])));
                    //value = 1.0;   
                    //value = 0.0;
                                        // if(t==0.0) value = 0.0;
                    // else{ // set stagnation pressure
                    //     hig_cell *c = sd_get_cell_with_point(nsaux->sdp, center);
                    //     Point ccenter;
                    //     hig_get_center(c, ccenter);
                    //     sim_stencil *stn = stn_create();
                    //     real p = compute_value_at_point(nsaux->sdp, ccenter, center, 1.0, nsaux->dpp, stn);
                    //     real p_0 = dpdx * L;
                    //     value = sqrt(2.0*fabs(p_0-p));
                    //     printf("p = %f, p_0 = %f, value = %f\n", p, p_0, value);
                    //     stn_destroy(stn);
                    // }
                    break;
                case 1:
                    value = 0.0;
                                        break;
            }
            break;
        case 1:
            switch (dim) {
                case 0:
                    value = 0.0;
                    //value = 1.0;
                                        break;
                case 1:
                    value = 0.0;
                                        break;
            }
            break;
        case 2:
            switch (dim) {
                case 0:
                    value = 0.0;
                                        break;
                case 1:
                    value = 0.0;
                                        break;
            }
            break;
        case 3:
            switch (dim) {
                case 0:
                    value = 0.0;
                                        break;
                case 1:
                    value = 0.0;
                                        break;
            }
            break;
    }
    return value; 
}

// Value of the cell source term at boundary
real get_boundary_source_term(int id, Point center, real t) {
    real value = 0.0;
    return value; 
}

// Value of the facet source term at boundary
real get_boundary_facet_source_term(int id, Point center, int dim, real t) {
    real value = 0.0;
    return value; 
}

// Value of the boundary viscosity
real get_boundary_viscosity(int id, Point center, real q, real t) {
    real value = 1.0;
    return value; 
}

// Define the user function for viscoelastic flow
void calculate_m_user(real Re, real De, real beta, real tr, real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real M_aux[DIM][DIM], real tol) {
    // Calculate the matrix MM and BB for Oldroyd-B model
    for (int i = 0; i < DIM; i++) {
        for (int j = i+1; j < DIM; j++) {
            M_aux[i][j] = 0.0;
            M_aux[j][i] = 0.0;
        }
        real jlambda = get_kernel_jacobian(i, lambda[i], tol);
        M_aux[i][i]  = (1.0-lambda[i])*jlambda;
    }
}

// *******************************************************************
// Navier-Stokes main program
// *******************************************************************

// Main program for the Navier-Stokes simulation 
int main (int argc, char *argv[]) {
    // Initialize the total time counting
    START_CLOCK(total);
    // Number of tasks
    int ntasks;
    // Identifier of the process
    int myrank;
    // Initializing Navier-Stokes solver
    higflow_initialize(&argc, &argv, &myrank, &ntasks);
    // Create Navier-Stokes solver
    higflow_solver *ns = higflow_create();
    // Load the data files
    higflow_load_data_files(argc, argv, ns); 
	print0f("=+=+=+= Load Controllers and Parameters =+=+=+=+=+=+=+=+=+=+=+=+=\n");
    higflow_load_controllers_and_parameters_yaml(ns, myrank);
        // set the external functions
    higflow_set_external_functions(ns, get_pressure, get_velocity, 
        get_source_term, get_facet_source_term,
        get_boundary_pressure, get_boundary_velocity,
        get_boundary_source_term, get_boundary_facet_source_term); 
    // Set the order of the interpolation to be used in the SD. 
    int order_center = 3;
    int order_facet = 3;
    // Set the cache: Reuse interpolation, 0 on, 1 off
    int cache = 1;

    // Create the simulation domain
    higflow_create_domain(ns, cache, order_center); 
    
    // Initialize the domain
    print0f("=+=+=+= Load Domain =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=\n");
    //higflow_initialize_domain(ns, ntasks, myrank, order_facet); 
    higflow_initialize_domain_yaml(ns, ntasks, myrank, order_facet); 

    // Initialize the boundaries
    print0f("=+=+=+= Load Bondary Condtions =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=\n");
    //higflow_initialize_boundaries(ns);
    higflow_initialize_boundary_conditions_yaml(ns);

    // Creating distributed property  
    higflow_create_distributed_properties(ns);
    // Initialize distributed properties
    if (ns->par.step == 0) higflow_initialize_distributed_properties(ns);
    // Create the linear system solvers
    higflow_create_solver(ns);

    // Load the properties form 
    if (ns->par.step > 0) {
        // Loading the velocities 
        print0f("===> Loading t = %f <===\n",ns->par.t);
        higflow_load_properties(ns, myrank, ntasks);
    }
    // Printing the properties to visualize: first step
    //higflow_load_properties(ns, myrank, ntasks);
    if (ns->par.step == 0) {
        print0f("===> Printing frame: %4d <====> tp = %15.10lf <===\n",ns->par.frame, ns->par.tp);
        higflow_print_vtk(ns, myrank);
        // frame update
        ns->par.frame++;
    }
    
    // ********************************************************
    // Begin Loop for the Navier-Stokes equations integration
    // ********************************************************
    // Descomente para carregar uma simulacao anterior

    for (int step0 = ns->par.initstep; ns->par.step < ns->par.finalstep; ns->par.step++) {
        // Print the step
        print0f("===> Step:        %7d <====> t  = %15.10lf <===\n", ns->par.step, ns->par.t);
        // Start the first step time
        if (ns->par.step == step0)  START_CLOCK(firstiter); 
        // Update velocities and pressure using the projection method 
        higflow_solver_step(ns);
        // Time update 
        ns->par.t += ns->par.dt;
        // Stop the first step time
        if (ns->par.step == step0) STOP_CLOCK(firstiter); 
        // Printing
        if (ns->par.t >= ns->par.tp + ns->par.dtp) {
            // tprint update
            ns->par.tp += ns->par.dtp;
            // Printing the properties to visualize 
            if (myrank == 0) 
                printf("===> Printing frame: %4d <====> tp = %15.10lf <===\n",ns->par.frame, ns->par.tp);
            higflow_print_vtk(ns, myrank);
            // frame update
            ns->par.frame++;
        }
        // Saving the properties
        if (ns->par.t >= ns->par.ts + ns->par.dts) {
            // tsave update
            ns->par.ts += ns->par.dts;
            if (myrank == 0) 
                printf("===> Saving         <====> ts = %15.10lf <===\n",ns->par.ts);
            // Saving the properties
            higflow_save_properties(ns, myrank, ntasks);
            higflow_save_parameters(ns, myrank);
            higflow_save_controllers(ns, myrank);
        }
    }
    // ********************************************************
    // End Loop for the Navier-Stokes equations integration
    // ********************************************************

    // Saving
    if (myrank == 0) 
        print0f("===> Saving               <====> t  = %15.10lf <===\n",ns->par.t);
    // Saving the properties
    higflow_save_properties(ns, myrank, ntasks);
    higflow_save_parameters(ns, myrank);
    higflow_save_controllers(ns, myrank);
    // Destroy the Navier-Stokes object
    higflow_destroy(ns);
    // Stop the total time
    STOP_CLOCK(total);
    // Getting the execution time 
    if(myrank == 0) {
        DEBUG_INSPECT(GET_NSEC_CLOCK(total)/1.0e9, %lf);
        DEBUG_INSPECT(GET_NSEC_CLOCK(firstiter)/1.0e9, %lf);
    }
}
