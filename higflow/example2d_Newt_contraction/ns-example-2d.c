// *******************************************************************
//  Example for HiG-Flow Solver - version 10/11/2016
// *******************************************************************

#include "ns-example-2d.h"

// *******************************************************************
// Extern functions for the Navier-Stokes program
// *******************************************************************

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
        case 5:
            value = 0.0;        
            break;
        case 6:
            value = 0.0;        
            break;
        case 7:
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
                                        // if(t == 0.0) value = 0.0;
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
        case 4:
            switch (dim) {
                case 0:
                    value = 0.0;
                                        break;
                case 1:
                    value = 0.0;
                                        break;
            }
            break;
        case 5:
            switch (dim) {
                case 0:
                    value = 0.0;
                                        break;
                case 1:
                    value = 0.0;
                                        break;
            }
            break;
        case 6:
            switch (dim) {
                case 0:
                    value = 0.0;
                                        break;
                case 1:
                    value = 0.0;
                                        break;
            }
            break;
        case 7:
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
    higflow_load_data_file_names(argc, argv, ns); 
    print0f("=+=+=+= Load Controllers and Parameters =+=+=+=+=+=+=+=+=+=+=+=+=\n");
    higflow_load_all_controllers_and_parameters_yaml(ns, myrank);
        // set the external functions
    higflow_set_external_functions(ns, get_pressure, get_velocity, 
        get_source_term, get_facet_source_term,
        get_boundary_pressure, get_boundary_velocity,
        get_boundary_source_term, get_boundary_facet_source_term); 
    // Set the order of the interpolation to be used in the SD. 
    int order_center = 2;
    int order_facet = 2;
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
    higflow_initialize_boundaries_yaml(ns);

    // Creating distributed property  
    higflow_create_distributed_properties(ns);
    // Initialize distributed properties
    if (ns->par.step == 0) higflow_initialize_distributed_properties(ns);
    // Create the linear system solvers
    higflow_create_solver(ns);
    // Load the properties form 
    if (ns->par.step > 0) {
        // Loading the velocities 
        if (myrank == 0) {
            printf("*********************************************************************************\n");
            printf("*********************************************************************************\n");
            printf("===> Reloading properties from previous simulation <====> step = %d <====> t = %15.10lf <===\n", ns->par.step, ns->par.t);
            printf("*********************************************************************************\n");
            printf("*********************************************************************************\n");
        }
        higflow_load_properties(ns, myrank, ntasks);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    print0f("=+=+ Saving Domain and Boundary Properties =+=+\n");
    higflow_save_domain_yaml(ns, myrank, ntasks);
    higflow_save_all_boundaries_yaml(ns, myrank, ntasks);
    higflow_save_all_controllers_and_parameters_yaml(ns, myrank); //copying necessary yamls

    // Printing the properties to visualize: first step
    if (ns->par.step == 0) {
        print0f("===> Printing frame: %4d <====> tp = %15.10lf <===\n",ns->par.frame, ns->par.tp);
        higflow_print_vtk(ns, myrank);
        //higflow_print_vtk2D_parallel_single(ns, myrank, ntasks);
        ns->par.tp += ns->par.dtp;
        ns->par.frame++;
        print0f("===> Saving               <====> ts = %15.10lf <===\n", ns->par.ts);
        higflow_save_all_controllers_and_parameters_yaml(ns, myrank);
        higflow_save_properties(ns, myrank, ntasks);
        ns->par.ts += ns->par.dts;
    }
    
    // ********************************************************
    // Begin Loop for the Navier-Stokes equations integration
    // ********************************************************

    for (int step0 = ns->par.initstep; ns->par.step <= ns->par.finalstep; ns->par.step++) {
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
        if (ns->par.t >= ns->par.tp) {
            print0f("===> Printing frame: %4d <====> tp = %15.10lf <===\n",ns->par.frame, ns->par.tp);
            higflow_print_vtk(ns, myrank);
            //higflow_print_vtk2D_parallel_single(ns, myrank, ntasks);
            ns->par.tp += ns->par.dtp;
            ns->par.frame++;
        }
        // Saving the properties
        if (ns->par.t >= ns->par.ts) {
            print0f("===> Saving               <====> ts = %15.10lf <===\n", ns->par.ts);
            higflow_save_all_controllers_and_parameters_yaml(ns, myrank);
            higflow_save_properties(ns, myrank, ntasks);
            ns->par.ts += ns->par.dts;
        }
    }
    // ********************************************************
    // End Loop for the Navier-Stokes equations integration
    // ********************************************************

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
