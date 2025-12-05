#include "ns-example-2d.h"

/******************************************************************************************/
/******************************************************************************************/
/*************************** multiphase user functions ************************************/
/******************************************************************************************/
/******************************************************************************************/


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

real func (Point p) {
    real value;
	Point c, pc, r, cav;
    ////////////// Zalesak's disk ////////////
	// c[0]     = 0.5; c[1]     = 0.75;
	// r[0]     = 0.15; r[1]     = 0.15;
	// cav[0]   = 0.05; cav[1]   = 0.2479;
	// pc[0]    = p[0] - c[0]; pc[1]    = p[1] - c[1];
	// // negative when in the circle
	// real outcirc = 1.0 - (pc[0]*pc[0]/(r[0]*r[0]) + pc[1]*pc[1]/(r[1]*r[1]));
	// real incav1  = (fabs(pc[0]) - 0.5*cav[0])/r[0]; ;
	// real incav2  = (pc[1] - (cav[1] - r[1]))/r[1];
	// // negative when in the cavity - both incav1 and incav2 must be negative
	// real incav   = max(incav1, incav2);
	// value = min(outcirc, incav);

    /////// Cavity ////////////
	// c[0] = 0.5; c[1] = 0.5;
	// r[0] = 1.0/6.0; r[1] = 1.0/6.0;

    /////// Shearing Droplet ////////////
    c[0] = p_par->center[0].val; c[1] = p_par->center[1].val;
    r[0] = 0.25; r[1] = 0.25;

    pc[0] = p[0] - c[0]; pc[1] = p[1] - c[1];
	real outcirc = 1.0 - (pc[0]*pc[0]/(r[0]*r[0]) + pc[1]*pc[1]/(r[1]*r[1]));
	value = outcirc;
	return value;
}

// normalized viscosity at phase 0 over all time
real get_viscosity0(Point center, real t) {
	return p_par->mu0.val;
}
// normalized viscosity at phase 1 over all time
real get_viscosity1(Point center, real t) {
	return p_par->mu1.val;
}
// normalized density at phase 0 over all time
real get_density0(Point center, real t) {
	return p_par->rho0.val;
}
// normalized density at phase 1 over all time
real get_density1(Point center, real t) {
	return p_par->rho1.val;
}

real compute_deformation_parameter(higflow_solver *ns) {
    int num_plic = ns->ed.mult.num_plic_lines;
    int num_plic_global;
    Point center, center_global, midseg;
    
    center[0] = 0.0; center[1] = 0.0;
    for(int i=0; i<num_plic; i++) {
        midseg[0] = 0.5*(ns->ed.mult.plic_lines[i][0][0] + ns->ed.mult.plic_lines[i][1][0]);
        midseg[1] = 0.5*(ns->ed.mult.plic_lines[i][0][1] + ns->ed.mult.plic_lines[i][1][1]);
        center[0] += midseg[0];
        center[1] += midseg[1];
    }
    MPI_Allreduce(center, center_global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&num_plic, &num_plic_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    center_global[0] /= num_plic_global;
    center_global[1] /= num_plic_global;
    if(psi_down == 0.0 && ns->ed.mult.contr.eoflow_either == true) {
        center_global[1] = 0.0;
    }

    real Rmin=INFINITY, Rmax=0.0, dist2;
    real Rmin_global, Rmax_global;
    for(int i=0; i<num_plic; i++) {
        midseg[0] = 0.5*(ns->ed.mult.plic_lines[i][0][0] + ns->ed.mult.plic_lines[i][1][0]);
        midseg[1] = 0.5*(ns->ed.mult.plic_lines[i][0][1] + ns->ed.mult.plic_lines[i][1][1]);
        dist2 = (midseg[0] - center_global[0])*(midseg[0] - center_global[0]) 
              + (midseg[1] - center_global[1])*(midseg[1] - center_global[1]);
        if(dist2 < Rmin) Rmin = dist2;
        if(dist2 > Rmax) Rmax = dist2;
    }
    Rmin = sqrt(Rmin); Rmax = sqrt(Rmax);
    MPI_Allreduce(&Rmin, &Rmin_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&Rmax, &Rmax_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    return (Rmax_global-Rmin_global)/(Rmax_global+Rmin_global);
}

void save_deformation_parameter(higflow_solver *ns, int myrank) {
    real defpar = compute_deformation_parameter(ns);
    if(myrank==0) {
        char filename[1024];
        sprintf(filename, "%s.deformation_parameter.txt", ns->par.namesave);
        FILE *f;
        if(ns->par.step==0) f = fopen(filename, "w");
        else f = fopen(filename, "a");
        fprintf(f, "%d %.10lf %.10lf\n", ns->par.step, ns->par.t, defpar);
        fclose(f);
    }
}
