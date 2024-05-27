// *******************************************************************
// *******************************************************************
//  Example for HiG-Flow Solver - version 10/11/2016
// *******************************************************************
// *******************************************************************

#include "ns-example-2d.h"
#define pi 3.1415926535897932384626434

// *******************************************************************
// Extern functions for the Navier-Stokes program
// *******************************************************************

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

// initial multiphase kernel conformation tensor
real get_tensor_multiphase(real fracvol, Point center, int i, int j, real t) {
    return get_kernel(i, 1.0, 0.0) * (i == j);
}

void calculate_m_user_multiphase(real fracvol, real lambda[DIM], real jlambda[DIM],  real B[DIM][DIM], real M_aux[DIM][DIM], real Re, real trS, ve_parameters *par0, ve_parameters *par1) {
    // // Calculate the matrix MM and BB for Oldroyd-B model
    // for (int i = 0; i < DIM; i++) {
    //     for (int j = i+1; j < DIM; j++) {
    //         M_aux[i][j] = 0.0;
    //         M_aux[j][i] = 0.0;
    //     }
    //     real jlambda = get_kernel_jacobian(i, lambda[i], tol);
    //     M_aux[i][i]  = (1.0-lambda[i])*jlambda;
    // }
}

real func (Point p) {
	real value;
	Point c, pc, r, cav;

	c[0] = 0.5; c[1] = 0.5;
	r[0] = 1.0/6.0; r[1] = 1.0/6.0;
	pc[0] = p[0] - c[0]; pc[1] = p[1] - c[1];
	real outcirc = 1.0 - (pc[0]*pc[0]/(r[0]*r[0]) + pc[1]*pc[1]/(r[1]*r[1]));
	value = outcirc;
	return value;
}

void Intersec(Point p0, Point p1, real f0, real f1, Point p) {
	real lambda = -f0/(f1-f0);
	p[0] = p0[0] + lambda*(p1[0]-p0[0]);
	p[1] = p0[1] + lambda*(p1[1]-p0[1]);
	return;
}

real volume_convex_facet(int n, real P[n][DIM]) {
	real area   = 0.0;
	real altura = 0.0;
	for (int i = 0; i< n-1; i++) {
		area   += P[i+1][0] - P[i][0];
		altura += P[i][1];
	}
	altura += P[n-1][1];

	return altura*area/n;
}

real volume_convex_13(Point p0, Point p01, Point p02){
	real facet0[2][DIM], facet1[2][DIM], facet2[2][DIM];
	real volume = 0.0;

	facet0[0][0] = p0[0];
	facet0[0][1] = p0[1];
	facet0[1][0] = p01[0];
	facet0[1][1] = p01[1];

	volume += volume_convex_facet(2,facet0);

	facet1[0][0] = p01[0];
	facet1[0][1] = p01[1];
	facet1[1][0] = p02[0];
	facet1[1][1] = p02[1];

	volume += volume_convex_facet(2,facet1);

	facet2[0][0] = p02[0];
	facet2[0][1] = p02[1];
	facet2[1][0] = p0[0];
	facet2[1][1] = p0[1];

	volume += volume_convex_facet(2,facet2);

	return volume;
}

real volume_convex_22(Point p0, Point p1, Point p13, Point p02){
	real facet0[2][DIM], facet1[2][DIM], facet2[2][DIM], facet3[2][DIM];
	real volume = 0.0;

	facet0[0][0] = p0[0];
	facet0[0][1] = p0[1];
	facet0[1][0] = p1[0];
	facet0[1][1] = p1[1];

	volume += volume_convex_facet(2,facet0);

	facet1[0][0] = p1[0];
	facet1[0][1] = p1[1];
	facet1[1][0] = p13[0];
	facet1[1][1] = p13[1];

	volume += volume_convex_facet(2,facet1);

	facet2[0][0] = p13[0];
	facet2[0][1] = p13[1];
	facet2[1][0] = p02[0];
	facet2[1][1] = p02[1];

	volume += volume_convex_facet(2,facet2);

	facet3[0][0] = p02[0];
	facet3[0][1] = p02[1];
	facet3[1][0] = p0[0];
	facet3[1][1] = p0[1];

	volume += volume_convex_facet(2,facet3);

	return volume;
}

real square_case_13(Point center, Point delta, Point p0, Point p1, Point p2, Point p3, real f0, real f1, real f2, real f3, int var){
	Point p01, p02, p13, p23;
	real value = 0.0;
	int caso;
	if (var == 1) {
		if (f0 > 0.0) {
			//caso = 2;
			Intersec(p0, p1, f0, f1, p01);
			Intersec(p0, p2, f0, f2, p02);
			value = fabs(volume_convex_13(p0, p01, p02));
			//printf("Caso %d, volume = %.18lf\n", caso, value);
			return value;
		} else if (f1 > 0.0) {
			//caso = 3;
			Intersec(p0, p1, f0, f1, p01);
			Intersec(p1, p3, f1, f3, p13);
			value = fabs(volume_convex_13(p1, p13, p01));
			//printf("Caso %d, volume = %.18lf\n", caso, value);
			return value;
		} else if (f2 > 0.0) {
			//caso = 4;
			Intersec(p0, p2, f0, f2, p02);
			Intersec(p2, p3, f2, f3, p23);
			value = fabs(volume_convex_13(p2, p02, p23));
			//printf("Caso %d, volume = %.18lf\n", caso, value);
			return value;
		} else if (f3 > 0.0) {
			//caso = 5;
			Intersec(p1, p3, f1, f3, p13);
			Intersec(p2, p3, f2, f3, p23);
			value = fabs(volume_convex_13(p3, p23, p13));
			//printf("Caso %d, volume = %.18lf\n", caso, value);
			return value;
		}
	} else {
		if (f0 <= 0.0) {
			//caso = 6;
			Intersec(p0, p1, f0, f1, p01);
			Intersec(p0, p2, f0, f2, p02);
			value = delta[0] * delta[1] - fabs(volume_convex_13(p0, p01, p02));
			//printf("Caso %d, volume = %.18lf\n", caso, value);
			return value;
		} else if (f1 <= 0.0) {
			//caso = 7;
			Intersec(p0, p1, f0, f1, p01);
			Intersec(p1, p3, f1, f3, p13);
			value = delta[0] * delta[1] - fabs(volume_convex_13(p1, p13, p01));
			//printf("Caso %d, volume = %.18lf\n", caso, value);
			return value;
		} else if (f2 <= 0.0) {
			//caso = 8;
			Intersec(p0, p2, f0, f2, p02);
			Intersec(p2, p3, f2, f3, p23);
			value = delta[0] * delta[1] - fabs(volume_convex_13(p2, p02, p23));
			//printf("Caso %d, volume = %.18lf\n", caso, value);
			return value;
		} else if (f3 <= 0.0) {
			//caso = 9;
			Intersec(p1, p3, f1, f3, p13);
			Intersec(p2, p3, f2, f3, p23);
			value = delta[0] * delta[1] - fabs(volume_convex_13(p3, p13, p23));
			//printf("Caso %d, volume = %.18lf\n", caso, volume);
			return value;
		}
	}
}

real square_case_22(Point center, Point delta, Point p0, Point p1, Point p2, Point p3, real f0, real f1, real f2, real f3){
	Point p01, p02, p13, p23;
	real value  = 0.0;
	int caso;
	if ((f0 > 0.0) && (f1 > 0.0)){
		//caso = 10;
		Intersec(p0,p2,f0,f2,p02);
		Intersec(p1,p3,f1,f3,p13);
		value = fabs(volume_convex_22(p0, p1, p13, p02));
		//printf("Caso %d, volume = %.18lf\n", caso, value);
		return value;
	} else if ((f0 > 0.0) && (f2 > 0.0)){
		//caso = 11;
		Intersec(p0,p1,f0,f1,p01);
		Intersec(p2,p3,f2,f3,p23);
		value = fabs(volume_convex_22(p2, p0, p01, p23));
		//printf("Caso %d, volume = %.18lf\n", caso, value);
		return value;
	} else if ((f1 > 0.0) && (f3 > 0.0)){
		//caso = 12;
		Intersec(p0,p1,f0,f1,p01);
		Intersec(p2,p3,f2,f3,p23);
		value = fabs(volume_convex_22(p1, p3, p23, p01));
		//printf("Caso %d, volume = %.18lf\n", caso, value);
		return value;
	} else if ((f2 > 0.0) && (f3 > 0.0)){
		//caso = 13;
		Intersec(p0,p2,f0,f2,p02);
		Intersec(p1,p3,f1,f3,p13);
		value = fabs(volume_convex_22(p3, p2, p02, p13));
		//printf("Caso %d, volume = %.18lf", caso, value);
		return value;
	} else {
		// caso = 14 e caso = 15
		printf("######## Malha pouco refinada para este problema ########\n");
		exit(1);
	}
}

// Volume fraction
real get_fracvolN(Point center, Point delta, real t) {
	Point p0, p1, p2, p3;
	real  f0, f1, f2, f3;
	real  value;
	int var = 0;

	// Canto inferior esquerdo
	p0[0] = center[0] - 0.5*delta[0];
	p0[1] = center[1] - 0.5*delta[1];
	f0    = func(p0);
	if (f0 > 0.0) var += 1;

	// Canto inferior direito
	p1[0] = center[0] + 0.5*delta[0];
	p1[1] = center[1] - 0.5*delta[1];
	f1    = func(p1);
	if (f1 > 0.0) var += 1;

	// Canto superior esquerdo
	p2[0] = center[0] - 0.5*delta[0];
	p2[1] = center[1] + 0.5*delta[1];
	f2    = func(p2);
	if (f2 > 0.0) var += 1;

	// Canto superior direito
	p3[0] = center[0] + 0.5*delta[0];
	p3[1] = center[1] + 0.5*delta[1];
	f3    = func(p3);
	if (f3 > 0.0) var += 1;

	if (var == 0){
		value = 0.0;
	} else if (var == 4){
		value = delta[0]*delta[1];
	} else if ((var == 1)||(var == 3)){
		value = square_case_13(center, delta, p0, p1, p2, p3, f0, f1, f2, f3, var);
	} else {
		value = square_case_22(center, delta, p0, p1, p2, p3, f0, f1, f2, f3);
	}

	return value;
}

real get_fracvol(Point center, Point delta, real t) {
	Point p0, p1, p2, p3;
	real  f0, f1, f2, f3;
	real  value;
	int var = 0;

	// Canto inferior esquerdo
	p0[0] = center[0] - 0.5*delta[0];
	p0[1] = center[1] - 0.5*delta[1];
	f0    = func(p0);
	if (f0 > 0.0) var += 1;

	// Canto inferior direito
	p1[0] = center[0] + 0.5*delta[0];
	p1[1] = center[1] - 0.5*delta[1];
	f1    = func(p1);
	if (f1 > 0.0) var += 1;

	// Canto superior esquerdo
	p2[0] = center[0] - 0.5*delta[0];
	p2[1] = center[1] + 0.5*delta[1];
	f2    = func(p2);
	if (f2 > 0.0) var += 1;

	// Canto superior direito
	p3[0] = center[0] + 0.5*delta[0];
	p3[1] = center[1] + 0.5*delta[1];
	f3    = func(p3);
	if (f3 > 0.0) var += 1;

	if (var == 0){
		value = 0.0;
	} else if (var == 4){
		value = delta[0]*delta[1];
	}
	else{
		value =0.0;
		int N=16;
		real delta_new[2];
		delta_new[0]=delta[0]/N;delta_new[1]=delta[1]/N;
		real center_new[2];
		for(int i=0;i<N;i++)
		{	
			for (int j=0;j<N;j++)
			{
				center_new[0]=p0[0]+(i+0.5)*delta_new[0];
				center_new[1]=p0[1]+(j+0.5)*delta_new[1];
				value = value + get_fracvolN(center_new,delta_new, t);
			}
		}
	}	

	value = value/delta[0]/delta[1];
	return value;
}


// Value of the pressure
real get_pressure(Point center, real t) {
    real value = 0.0;
	return value;
}

// Value of the velocity
real get_velocity(Point center, int dim, real t) {
	real value;
	real x = center[0];
	real y = center[1];
	value = 0.0;
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
real get_viscosity0(Point center, real t) {
	real value;
	value = 1.0;
	return value;
}

// Value of the viscosity 
real get_viscosity1(Point center, real t) {
	real value;
	value = 1.0;
	return value;
}

// Value of the density 
real get_density0(Point center, real t) {
	real value;
	value = 1.0;
	return value;
}

// Value of the density 
real get_density1(Point center, real t) {
	real value;
	value = 1.0;
	//value = 1.0;
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
	real x = center[0];
	real y = center[1];
	switch (id) {
	case 0:
		switch (dim) {
		case 0:
			// value = 1.5*(1.0 - center[1]*center[1]);
			value = 0.0;
			break;
		case 1:
			value = 0.0;
			break;
		}
		break;
	case 1:
		switch (dim) {
		case 0:
			value = 8.0*(1.0 + tanh(8.0*t - 4.0))*x*x*(1.0 - x)*(1.0 - x);
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
	higflow_create_domain_multiphase(ns, cache, order_center, get_viscosity0, get_viscosity1, get_density0, get_density1, get_fracvol);
	higflow_create_domain_multiphase_viscoelastic(ns, get_tensor_multiphase, get_kernel,
                                                  get_kernel_inverse, get_kernel_jacobian);
    higflow_define_user_function_multiphase_viscoelastic(ns, calculate_m_user_multiphase);
	
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
		higflow_solver_step_multiphase_viscoelastic(ns);
		ns->par.stepaux=ns->par.stepaux+1;
		
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
