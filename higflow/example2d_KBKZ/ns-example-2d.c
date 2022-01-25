// *******************************************************************
// *******************************************************************
//  Example for HiG-Flow Solver - version 03/05/2019
// *******************************************************************
// *******************************************************************

#include "ns-example-2d.h"
real tensormax00 = 0.0;
real tensormax10 = 0.0;
real tensormax01 = 0.0;
real tensormax11 = 0.0;

// *******************************************************************
// Extern functions for the Navier-Stokes program
// *******************************************************************

// Value of the Tensor
real get_tensor(Point center, int i, int j, real t) {
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


// Print the velocity
void print_velocity (higflow_solver *ns, int myrank, int dim, real x, real yf, real yl, real dy) {
    char filename[1024];
    snprintf(filename,sizeof filename,"Dados/Velocidade/velocity_%d_%d_%d.dat",dim,myrank,ns->par.frame);
    FILE *fd = fopen(filename, "w");

    if (fd != NULL) {
        // Get the local sub-domain
        sim_domain *sdp = psd_get_local_domain(ns->psdp);
        // Get the local partitioned domain for facets
        sim_facet_domain *sfdu[DIM];
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        Point P;
	P[0] = x;
	for (real y = yf; y <= yl; y += dy) {
           P[1] = y;
	   real u = compute_facet_value_at_point(sfdu[dim],P,P,1.0,ns->dpu[dim],ns->stn); 
	   fprintf(fd,"%16.10lf %16.10lf %16.10lf \n",x,y,u);
	}
        // Destroy the iterator
        fclose(fd);
    } else {
        printf("Arquivo %s nao aberto\n",filename);
        exit(1);
    }
}

// Print the Polymeric Tensor
real print_tensor(higflow_solver *ns, int myrank, int i, int j, real x, real yf, real yl, real dy) {
    char filename[1024];
    snprintf(filename,sizeof filename,"Dados/Tensores/tensor_%d_%d_%d_%d.dat",i,j,myrank,ns->par.frame);
    FILE *fd = fopen(filename, "w");
    real max=0.0;
     if (fd != NULL) {
        // Get the cosntants
        real Re   = ns->par.Re;
        //real De   = ns->ed.im.par.De;
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        Point P;
	P[0] = x;
	for (real y = yf; y <= yl; y += dy) {
            P[1] = y;
            // Get the velocity derivative tensor Du and the Kernel tensor
            real Du[DIM][DIM];
            // Get Du
            Du[i][j] = compute_value_at_point(ns->ed.sdED, P, P, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
            Du[j][i] = compute_value_at_point(ns->ed.sdED, P, P, 1.0, ns->ed.im.dpD[j][i], ns->ed.stn);
            // Get T tensor
            real D  = (Du[i][j]+Du[j][i]);
            real S = compute_value_at_point(ns->ed.sdED, P, P, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
            real T = S + D/Re; 
	    fprintf(fd,"%16.10lf %16.10lf %16.10lf %16.10lf \n",x,y,T,S);
	}
        // Loop for each cell
        fclose(fd);
    } else {
        printf("Arquivo %s nao aberto\n",filename);
        exit(1);
    }
   return max; 
}



// Print the Pressure
void print_pressure(higflow_solver *ns, int myrank, int i, int j, real x, real yf, real yl, real dy,int p) {
    char filename[1024];
    snprintf(filename,sizeof filename,"Dados/Pressao/pressure_%d_%d_%d.dat",myrank,p,ns->par.frame);
    FILE *fd = fopen(filename, "w");

    if (fd != NULL) {
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        Point P;
	    P[0] = x;
     	for (real y = yf; y <= yl; y += dy) {
            P[1] = y;
            real p = compute_value_at_point(ns->sdp, P, P, 1.0, ns->dpp, ns->ed.stn);
	        fprintf(fd,"%16.10lf %16.10lf %16.10lf \n",x,y,p);
	    }
        // Loop for each cell
        fclose(fd);
    } else {
        printf("Arquivo %s nao aberto\n",filename);
        exit(1);
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
    // Load the parameters data for Navier-Stokes simulation
    higflow_load_parameters(ns, myrank);
    // Load the controllers data for Navier-Stokes simulation
    higflow_load_controllers(ns, myrank);
    // set the external functions
    higflow_set_external_functions(ns, get_pressure, get_velocity, 
        get_source_term, get_facet_source_term,
        get_boundary_pressure,  get_boundary_velocity,
        get_boundary_source_term, get_boundary_facet_source_term); 
    // Set the order of the interpolation to be used in the SD. 
    int order_center = 3;
    int order_facet = 2;
    // Set the cache: Reuse interpolation, 0 on, 1 off
    int cache = 1;
    // Create the simulation domain
    higflow_create_domain(ns, cache, order_center); 
    // Create the simulation domain for non newtonian simulation
    //higflow_create_domain_generalized_newtonian(ns, cache, order_center,
    //    get_viscosity, get_boundary_viscosity); 
    //higflow_create_domain_viscoelastic(ns, cache, order_center,
    //    get_tensor, get_kernel, get_kernel_inverse, get_kernel_jacobian); 
    higflow_create_domain_viscoelastic_integral(ns, cache, order_center, get_tensor); 
    // Initialize the domain
    higflow_initialize_domain(ns, ntasks, myrank, order_facet); 
    // Load the controllers data for viscoelastic simulation
   // higflow_load_viscoelastic_controllers(ns, myrank);
    // Load the parameters data for viscoelastic simulation
   // higflow_load_viscoelastic_parameters(ns, myrank);
    // Load the controllers data for viscoelastic simulation
    higflow_load_viscoelastic_integral_controllers(ns, myrank);
    // Load the parameters data for viscoelastic simulation
    higflow_load_viscoelastic_integral_parameters(ns, myrank);
    // Set the user model
    //higflow_define_user_function_viscoelastic(ns, calculate_m_user);
    // Initialize the boundaries
    higflow_initialize_boundaries(ns); 
    // Creating distributed property  
    higflow_create_ditributed_properties(ns);
    // Creating distributed property for generalized newtonian simulation
    //higflow_create_ditributed_properties_generalized_newtonian(ns);
    //higflow_create_ditributed_properties_viscoelastic(ns);
    higflow_create_ditributed_properties_viscoelastic_integral(ns);
    // Initialize distributed properties
    higflow_initialize_distributed_properties(ns);
    // Create the linear system solvers
    higflow_create_solver(ns);
    // Load the properties form 
    if (ns->par.step > 0) {
        // Loading the velocities 
        printf("===> Loading t = %f <===\n",ns->par.t);
        //higflow_load_properties(ns);
    }
    // Printing the properties to visualize: first step
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
    int step0 = ns->par.step;
    for (int step = ns->par.step; step < ns->par.finalstep; step++) {
        // Print the step
        if (myrank == 0) 
            printf("===> Step:        %7d <====> t     = %15.10lf <===\n", step, ns->par.t);
            printf("===> Step:        %7d <====> tdim  = %15.10lf seg <===\n", step, (ns->par.t)*0.06);
        // Start the first step time
        if (step == step0)  START_CLOCK(firstiter); 
        // Update velocities and pressure using the projection method 
        //higflow_solver_step(ns);
        //higflow_solver_step_gen_newt(ns);
        //higflow_solver_step_viscoelastic(ns);
          higflow_solver_step_viscoelastic_integral(ns);
        
        // Time update 
        ns->par.t += ns->par.dt;
        //--- compute norm ---
        //compute_and_print_error_norm_velocity(ns, myrank);
        //compute_and_print_error_norm_pressure(ns, myrank);
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
        ////////////////////////////////////////////////////////////////////////////////////    
            print_velocity(ns, myrank, 0, 5.0, 0.0, 1.0, 0.125);
	    print_pressure(ns, myrank, 0, 0, 5.0-0.125/2.0,  0.0, 1.0, 0.125,1);
            print_pressure(ns, myrank, 0, 0, 5.0+0.125/2.0,  0.0, 1.0, 0.125,2);
	    real max00, max01, max10, max11, tol;
	    tol=1.0e-8;
            max00 = print_tensor(ns, myrank, 0, 0, 5.0,  0.0, 1.0, 0.125);
            max01 = print_tensor(ns, myrank, 0, 1, 5.0,  0.0, 1.0, 0.125);
            max10 = print_tensor(ns, myrank, 1, 0, 5.0,  0.0, 1.0, 0.125);
            max11 = print_tensor(ns, myrank, 1, 1, 5.0,  0.0, 1.0, 0.125);
	    printf("max00 = %lf, max01 = %lf, max10 = %lf, max11 = %lf\n",max00,max01,max10,max11);
           // printf("ResiduoT00 = %lf, ResiduoT01 = %lf, ResiduoT10 = %lf, ResiduoT11 = %lf\n",tensormax00-max00, tensormax01-max01, tensormax10-max10, tensormax11-max11);
            //if ((fabs(tensormax00 - max00) < tol)&&(fabs(tensormax01 - max01) < tol) &&(fabs(tensormax10 - max10) < tol) &&(fabs(tensormax11 - max11) < tol) ) {
               //step = ns->par.finalstep; 
 	    //}
            tensormax00 = max00;
            tensormax01 = max01;
            tensormax10 = max10;
            tensormax11 = max11;
       /////////////////////////////////////////////////////////////////////////////     
            
            // frame update
            ns->par.frame++;
        }
        // Saving the properties
        if (ns->par.t >= ns->par.ts + ns->par.dts) {
            // tsave update
            ns->par.ts += ns->par.dts;
            if (myrank == 0) 
                printf("===> Saving               <====> ts = %15.10lf <===\n",ns->par.ts);
            // Saving the properties
            //higflow_save_properties(ns, myrank, ntasks);
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
    //higflow_save_properties(ns, myrank, ntasks);
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
