// *******************************************************************
// *******************************************************************
//  Example for HiG-Flow Solver - version 03/05/2019
// *******************************************************************
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
                    value = 1.5*(1.0 - center[1]*center[1]);
                    //value = 1.0;
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
void calculate_m_user(real lambda[DIM], real jlambda[DIM],  real B[DIM][DIM], real M_aux[DIM][DIM], real Re, real trS, ve_parameters *par) {
    // Calculate the matrix MM and BB for Oldroyd-B model
    for (int i = 0; i < DIM; i++) {
        for (int j = i+1; j < DIM; j++) {
            M_aux[i][j] = 0.0;
            M_aux[j][i] = 0.0;
        }
        M_aux[i][i]  = (1.0-lambda[i])*jlambda[i];
    }
}

/*

real Calculaerro(higflow_solver *ns, int myrank, real new, real old) {

    real erro= 0.0;
    erro = fabs((old-new)/new)*100.0;

    return erro;
}

// calcula u the velocity
real calc_u(higflow_solver *ns, int myrank, int dim, real Px, real Py) {

        // Get the local sub-domain
        sim_domain *sdp = psd_get_local_domain(ns->psdp);
        // Get the local partitioned domain for facets
        sim_facet_domain *sfdu[DIM];
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        Point P;
	P[0] = Px;
        P[1] = Py;
        real u = compute_facet_value_at_point(sfdu[dim],P,P,1.0,ns->dpu[dim],ns->stn); 
        
    return u;
}

// calcula u the velocity
real calc_tau(higflow_solver *ns, int myrank, int i, int j, real Px, real Py) {
        // Get the cosntants
        real Re   = ns->par.Re;
        real De   = ns->ed.im.par.De;
        real beta = ns->ed.ve.par.beta;
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        Point P;
	P[0] = Px;
        P[1] = Py;
        // Get the velocity derivative tensor Du and the Kernel tensor
        real Du[DIM][DIM];
        // Get Du
        Du[i][j] = compute_value_at_point(ns->ed.sdED, P, P, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
        Du[j][i] = compute_value_at_point(ns->ed.sdED, P, P, 1.0, ns->ed.ve.dpD[j][i], ns->ed.stn);
        // Get T tensor
        real D  = 0.5*(Du[i][j]+Du[j][i]);
        real S = compute_value_at_point(ns->ed.sdED, P, P, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
        real tau = S + 2.0*(1.0-beta)*(D/Re);  
        
    return tau;
}
*/

/*
// Print the velocity
void print_velocity(higflow_solver *ns, int myrank, int dim, real x, real yf, real yl, real dy) {
    char filename[1024];
    snprintf(filename,sizeof filename,"Dados/Velocidadex5/velocity_%d_%d_%d.dat",dim,myrank,ns->par.frame);
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
    snprintf(filename,sizeof filename,"Dados/Tensoresx5/tensor_%d_%d_%d_%d.dat",i,j,myrank,ns->par.frame);
    FILE *fd = fopen(filename, "w");
    real max=0.0;
     if (fd != NULL) {
        // Get the cosntants
        real Re   = ns->par.Re;
        real De   = ns->ed.im.par.De;
        real beta = ns->ed.ve.par.beta;
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        Point P;
	P[0] = x;
	for (real y = yf; y <= yl; y += dy) {
            P[1] = y;
            // Get the velocity derivative tensor Du and the Kernel tensor
            real Du[DIM][DIM];
            // Get Du
            Du[i][j] = compute_value_at_point(ns->ed.sdED, P, P, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
            Du[j][i] = compute_value_at_point(ns->ed.sdED, P, P, 1.0, ns->ed.ve.dpD[j][i], ns->ed.stn);
            // Get T tensor
            real D  = 0.5*(Du[i][j]+Du[j][i]);
            real S = compute_value_at_point(ns->ed.sdED, P, P, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
            real T = S + 2.0*(1.0-beta)*(D/Re); 
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

// Print the velocity
void print_velocity_2 (higflow_solver *ns, int myrank, int dim, real y, real xf, real xl, real dx) {
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
	P[1] = y;
	for (real x = xf; x <= xl; x += dx) {
           P[0] = x;
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
real print_tensor_2(higflow_solver *ns, int myrank, int i, int j, real y, real xf, real xl, real dx) {
    char filename[1024];
    snprintf(filename,sizeof filename,"Dados/Tensores/tensor_%d_%d_%d_%d.dat",i,j,myrank,ns->par.frame);
    FILE *fd = fopen(filename, "w");
    real max=0.0;
     if (fd != NULL) {
        // Get the cosntants
        real Re   = ns->par.Re;
        real De   = ns->ed.im.par.De;
        real beta = ns->ed.ve.par.beta;
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        Point P;
	P[1] = y;
	for (real x = xf; x <= xl; x += dx) {
            P[0] = x;
            // Get the velocity derivative tensor Du and the Kernel tensor
            real Du[DIM][DIM];
            // Get Du
            Du[i][j] = compute_value_at_point(ns->ed.sdED, P, P, 1.0, ns->ed.ve.dpD[i][j], ns->ed.stn);
            Du[j][i] = compute_value_at_point(ns->ed.sdED, P, P, 1.0, ns->ed.ve.dpD[j][i], ns->ed.stn);
            // Get T tensor
            real D  = 0.5*(Du[i][j]+Du[j][i]);
            real S = compute_value_at_point(ns->ed.sdED, P, P, 1.0, ns->ed.ve.dpS[i][j], ns->ed.stn);
            real T = S + 2.0*(1.0-beta)*(D/Re); 
	    fprintf(fd,"%16.10lf %16.10lf %16.10lf %16.10lf \n",x,y,T,S);
	}
        // Loop for each cell
        fclose(fd);
    } else {
        printf("Arquivo %s nao aberto\n",filename);
        exit(1);
    }
   return max; 
} */

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
    // Create the simulation domain for non newtonian simulation
    //higflow_create_domain_generalized_newtonian(ns, cache, order_center,
    //    get_viscosity, get_boundary_viscosity); 
    higflow_create_domain_viscoelastic(ns, cache, order_center,get_tensor, 
		    get_kernel, get_kernel_inverse, get_kernel_jacobian); 
    // Initialize the domain
    higflow_initialize_domain_yaml(ns, ntasks, myrank, order_facet); 
    // define the user function for viscoelastic flow in case model is user set
    //higflow_define_user_function_viscoelastic(ns, calculate_m_user);

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
        higflow_solver_step_viscoelastic(ns);
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
