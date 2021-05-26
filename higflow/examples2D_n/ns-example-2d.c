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

// Value of the pressure
real get_pressure(Point center, real t) {
    real value      = 0.0;
    return value; 
}

// Value of the velocity
real get_velocity(Point center, int dim, real t) {
    real value      = 0.0;
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
                case 0:
                    //set max velocity = 8.0e-4 
                    //value = 3.2e-3*(-center[1]*center[1] + center[1]);
                    //value = 4.0*(-center[1]*center[1] + 0.25);
                    //set max velocity = 1.5 
                    value = -4.0*center[1]*(center[1] - 1.0);
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
    // Calculate the matrix MM and BB for Oldroyd model
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
    /* NAO MAIS NECESSARIO
	    // Load the parameters data for Navier-Stokes simulation
       //higflow_load_parameters(ns, myrank);
       // Load the controllers data for Navier-Stokes simulation
       //higflow_load_controllers(ns, myrank);
       // Load the controllers and parameters data for Navier-Stokes simulation
    */
	 printf("=+=+=+= Load Controllers and Parameters =+=+=+=+=+=+=+=+=+=+=+=+=\n");
    higflow_load_controllers_and_parameters(ns, myrank);
    // set the external functions
    higflow_set_external_functions(ns, get_pressure, get_velocity, 
        get_source_term, get_facet_source_term,
        get_boundary_pressure,  get_boundary_velocity,
        get_boundary_source_term, get_boundary_facet_source_term); 
    // Set the order of the interpolation to be used in the SD. 
    int order_center = 2;
    int order_facet = 2;
    // Set the cache: Reuse interpolation, 0 on, 1 off
    int cache = 1;
    // Create the simulation domain
    higflow_create_domain(ns, cache, order_center); 
    //higflow_create_domain_generalized_newtonian(ns, cache, order_center, get_viscosity); 
    //higflow_create_domain_viscoelastic(ns, cache, order_center, get_tensor, get_boundary_tensor, get_kernel, get_kernel_inverse, get_kernel_jacobian); 

    // Initialize the domain
    higflow_initialize_domain(ns, ntasks, myrank, order_facet); 
    // Load the controllers data for viscoelastic simulation
    //higflow_load_viscoelastic_controllers(ns, myrank);
    // Load the parameters data for viscoelastic simulation
    //higflow_load_viscoelastic_parameters(ns, myrank);
    // Set the user model
    //higflow_define_user_function_viscoelastic(ns, calculate_m_user);

    // Initialize the boundaries
    //higflow_initialize_boundaries(ns);
    printf("=+=+=+= Load Domain and Bondary Condtions =+=+=+=+=+=+=+=+=+=+=+=+=\n");
    higflow_initialize_boundaries_conditions(ns);

    // Creating distributed property  
    higflow_create_ditributed_properties(ns);
    // Creating distributed property for generalized newtonian simulation
    // higflow_create_ditributed_properties_generalized_newtonian(ns);
    // higflow_create_ditributed_properties_viscoelastic(ns);
    // Initialize distributed properties
    higflow_initialize_distributed_properties(ns);
    // Create the linear system solvers
    higflow_create_solver(ns);
    // Load the properties form 
    if (ns->par.step > 0) {
        // Loading the velocities 
        printf("===> Loading t = %f <===\n",ns->par.t);
        higflow_load_properties(ns, myrank, ntasks);
    }
    // Printing the properties to visualize: first step
    //higflow_load_properties(ns, myrank, ntasks);
    if (ns->par.step == 0) {
        if (myrank == 0) 
            printf("===> Printing frame: %4d <====> tp = %15.10lf <===\n",ns->par.frame, ns->par.tp);
        higflow_print_vtk(ns, myrank);
        // frame update
        ns->par.frame++;
    }
    // ********************************************************
    // Begin Loop for the Navier-Stokes equations integration
    // ********************************************************
    // Descomente para carregar uma simulacao anterior
    int step0 = ns->par.step;
    for (int step = ns->par.step; step < ns->par.finalstep; step++) {
        // Print the step
        if (myrank == 0) 
            printf("===> Step:        %7d <====> t  = %15.10lf <===\n", step, ns->par.t);
        // Start the first step time
        if (step == step0)  START_CLOCK(firstiter); 
        // Update velocities and pressure using the projection method 
        higflow_solver_step(ns);
        // higflow_solver_step_gen_newt(ns);
        // higflow_solver_step_viscoelastic(ns);
        // Time update 
        ns->par.t += ns->par.dt;
        //--- compute norm ---
        // compute_and_print_error_norm_velocity(ns, myrank);
        // compute_and_print_error_norm_pressure(ns, myrank);
        //--- compute norm ---
        // Stop the first step time
        if (step == step0) STOP_CLOCK(firstiter); 
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
            // Saving the parameters
            //higflow_save_parameters(ns, myrank);
            // Saving the controllers
            //higflow_save_controllers(ns, myrank);
        }
        // Step update
        ns->par.step = step;
    }
    (ns->par.step)++;
    // ********************************************************
    // End Loop for the Navier-Stokes equations integration
    // ********************************************************

    // Saving
    if (myrank == 0) 
        printf("===> Saving               <====> t  = %15.10lf <===\n",ns->par.t);
    // Saving the properties
    higflow_save_properties(ns, myrank, ntasks);
    // Saving the parameters
    //higflow_save_parameters(ns, myrank);
    // Saving the controllers
    //higflow_save_controllers(ns, myrank);
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
