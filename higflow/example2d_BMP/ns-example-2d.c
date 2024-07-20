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

// Value of the velocity of the BMP model at the boundary
real get_velocity_BMP(Point center, real Re, real GMI, real LA, real Phi, real GP) {
    real y1 = center[1] -0.5;
    real lh2 = 0.25;
    real lh4 = 0.0625;
    real y2 = pow(y1,2.0);
    real y4 = pow(y1,4.0);
    real GM = GMI/Phi;
    real b0 = pow(GM*LA,2.0)*pow(Re*GP,4.0);
    real c0 = 2.0*GM*LA*pow(Re*GP,2.0)*(1.0-2.0*Phi);
    real sqb0 = sqrt(b0);
    real U0 = 0.25*Re*GP*(lh2-y2);
    real U1 = sqrt(b0*y4-c0*y2+1.0)-sqrt(b0*lh4-c0*lh2+1.0);
    real U2a = 2.0*sqb0*sqrt(b0*lh4-c0*lh2+1.0) + 2.0*b0*lh2 - c0;
    real U2b = 2.0*sqb0*sqrt(b0*y4-c0*y2+1.0) + 2.0*b0*y2 - c0;
    real U2 = 0.5*(c0/sqb0)*log(U2a/U2b);
    real U3a = 2.0-c0*lh2 + 2.0*sqrt(b0*lh4-c0*lh2+1.0);
    real U3b = 2.0-c0*y2 + 2.0*sqrt(b0*y4-c0*y2+1.0);
    real U3 = log(U3a/U3b);
    real Ux = U0 - (1.0/(4.0*GM*LA*Re*GP))*(U1+U2+U3);
    real value      = (1/Phi)*Ux;
    return value; 
}

// Value of the velocity of the BMP model at the boundary
real get_velocity_BMP2(Point center, real Re, real GMI, real LA, real Phi, real GP) {
    real y1 = center[1] -0.5;
    real lh2 = 0.25;
    real lh4 = 0.0625;
    real y2 = pow(y1,2.0);
    real y4 = pow(y1,4.0);
    real GM = GMI;
    real b0 = pow(GM*LA,2.0)*pow(Re*GP,4.0);
    real c0 = 2.0*GM*LA*pow(Re*GP,2.0)*(1.0-2.0*Phi);
    real sqb0 = sqrt(b0);
    real U0 = 0.25*Re*GP*(lh2-y2);
    real U1 = sqrt(b0*y4-c0*y2+1.0)-sqrt(b0*lh4-c0*lh2+1.0);
    real U2a = 2.0*sqb0*sqrt(b0*lh4-c0*lh2+1.0) + 2.0*b0*lh2 - c0;
    real U2b = 2.0*sqb0*sqrt(b0*y4-c0*y2+1.0) + 2.0*b0*y2 - c0;
    real U2 = 0.5*(c0/sqb0)*log(U2a/U2b);
    real U3a = 2.0-c0*lh2 + 2.0*sqrt(b0*lh4-c0*lh2+1.0);
    real U3b = 2.0-c0*y2 + 2.0*sqrt(b0*y4-c0*y2+1.0);
    real U3 = log(U3a/U3b);
    real Ux = U0 - (1.0/(4.0*GM*LA*Re*GP))*(U1+U2+U3);
    real value      = Ux;
    return value; 
}


// Value of the velocity of the BMP model at the boundary
real get_velocity_BMP_proof(real yc, real Re, real GMI, real LA, real Phi, real GP) {
    real y1 = yc -0.5;
    real lh2 = 0.25;
    real lh4 = 0.0625;
    real y2 = pow(y1,2.0);
    real y4 = pow(y1,4.0);
    real GM = GMI/Phi;
    real b0 = pow(GM*LA,2.0)*pow(Re*GP,4.0);
    real c0 = 2.0*GM*LA*pow(Re*GP,2.0)*(1.0-2.0*Phi);
    real sqb0 = sqrt(b0);
    real U0 = 0.25*Re*GP*(lh2-y2);
    real U1 = sqrt(b0*y4-c0*y2+1.0)-sqrt(b0*lh4-c0*lh2+1.0);
    real U2a = 2.0*sqb0*sqrt(b0*lh4-c0*lh2+1.0) + 2.0*b0*lh2 - c0;
    real U2b = 2.0*sqb0*sqrt(b0*y4-c0*y2+1.0) + 2.0*b0*y2 - c0;
    real U2 = 0.5*(c0/sqb0)*log(U2a/U2b);
    real U3a = 2.0-c0*lh2 + 2.0*sqrt(b0*lh4-c0*lh2+1.0);
    real U3b = 2.0-c0*y2 + 2.0*sqrt(b0*y4-c0*y2+1.0);
    real U3 = log(U3a/U3b);
    real Ux = U0 - (1.0/(4.0*GM*LA*Re*GP))*(U1+U2+U3);
    real value      = (1/Phi)*Ux;
    return value; 
}

// Value of the velocity of the BMP model at the boundary
real get_velocity_BMP_proof2(real yc, real Re, real GMI, real LA, real Phi, real GP) {
    real y1 = yc -0.5;
    real lh2 = 0.25;
    real lh4 = 0.0625;
    real y2 = pow(y1,2.0);
    real y4 = pow(y1,4.0);
    real GM = GMI;
    real b0 = pow(GM*LA,2.0)*pow(Re*GP,4.0);
    real c0 = 2.0*GM*LA*pow(Re*GP,2.0)*(1.0-2.0*Phi);
    real sqb0 = sqrt(b0);
    real U0 = 0.25*Re*GP*(lh2-y2);
    real U1 = sqrt(b0*y4-c0*y2+1.0)-sqrt(b0*lh4-c0*lh2+1.0);
    real U2a = 2.0*sqb0*sqrt(b0*lh4-c0*lh2+1.0) + 2.0*b0*lh2 - c0;
    real U2b = 2.0*sqb0*sqrt(b0*y4-c0*y2+1.0) + 2.0*b0*y2 - c0;
    real U2 = 0.5*(c0/sqb0)*log(U2a/U2b);
    real U3a = 2.0-c0*lh2 + 2.0*sqrt(b0*lh4-c0*lh2+1.0);
    real U3b = 2.0-c0*y2 + 2.0*sqrt(b0*y4-c0*y2+1.0);
    real U3 = log(U3a/U3b);
    real Ux = U0 - (1.0/(4.0*GM*LA*Re*GP))*(U1+U2+U3);
    real value      = Ux;
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
//real get_viscosity(Point center, real q, real t) {
//    real value = 1.0;
//    return value; 
//}

// Value of the non-Newtonian viscosity
real get_viscosity(Point center, real q, real t, real beta, real struct_par) {
    //For the BMP model, viscosity = 1/fluidity (fluidity is the structural parameter)
    real value = 1.0/struct_par;
    //real value = 1.0;
    return value; 
}

// Value of the viscosity
real get_viscosity2(Point center, real q, real t, real beta, real struct_par) {
    real ap = 2.0;
    real Kcy = 1.5;
    real n = 0.1;
    real m = n-1.0;
    real ma = m/ap;
    //real thixratio = 0.01;
    real glam = q*Kcy;
    real vlam = 1.0 + pow(glam,ap);
    //real etaP = 0.99;
    real value;
    //value = thixratio + (1.0-thixratio)*(pow(vlam,ma));
    value = pow(vlam,ma);
    return value; 
}

// Value of the viscosity
real get_viscosity4(Point center, real q, real t, real beta, real struct_par) {
    real qmin = 0.01;
    real n = 0.50;
    real m = n-1.0;
    real etamax = pow(qmin,m);
    real value;
    if (q == 0 || q <= qmin) {
        value = 1.0;
    }
    else {
        value = pow(q,m)/etamax;
    }
    return value; 
}

// Value of the viscosity
real get_viscosity1(Point center, real q, real t, real beta, real struct_par) {
    real qmin = 1E-5;
    real n = 0.20;
    real m = n-1.0;
    real etamax = pow(qmin,m);
    real value;
    if (q == 0 || q <= qmin) {
        value = etamax;
    }
    else {
        value = pow(q,m);
    }
    return value; 
}

// Value of the viscosity
real get_viscosity3(Point center, real q, real t, real beta, real struct_par) {
    real qmin = 0.01;
    real qmax = 1.0;
    real n = 0.20;
    real m = n-1.0;
    real etamax = pow(qmin,m);
    real value;
    if (q == 0 || q <= qmin) {
        value = 1.0;
    }
    else if (q >= qmax) {
        value = pow(qmax,m)/etamax;
    }
    else {
        value = pow(q,m)/etamax;
    }
    return value; 
}

// Value of the structural parameter
real get_structpar(Point center, real q, real t, real beta, real Phi, real Lambda, real Gamma) {
    //For the BMP model, the structural parameter has a value of Phi when the fluid is at rest
    //real value = Phi/(1.0 - (beta/(1.0-beta))*Phi);
    //real value = 1.0/(1.0 - beta/(1.0-beta));
    real value = 1.0;
    //real value = Phi;
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
                    //value = 0.5*(-4.0*center[1]*(center[1] - 1.0));
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

// Value of the velocity at boundary
real get_boundary_velocity1(int id, Point center, int dim, real t) {
    real value;
    // m = (n+1)/n
    real n = 0.30;	
    real m = (n+1)/n;
    real vx = 1.0-2.0*center[1];
    real vmax = 1.0;
    switch (id) {
        case 0:
            switch (dim) {
                case 0:
		     //set max velocity = 8.0e-4 
                    //value = 3.2e-3*(-center[1]*center[1] + center[1]);
                    //value = 4.0*(-center[1]*center[1] + 0.25);
                    //set max velocity = 1.5 
                   //value = -4.0*center[1]*(center[1] - 1.0);
		   //value = 1.0-pow(1.0-2*fabs(center[1]),3.0);
		            if (center[1] <= 0.5 ) {
                        value = vmax*(1.0-pow(vx,m));
                        }
		            else {
                        value = vmax*(1.0-pow((-vx),m));
                        }
		   //value = -4.0*(pow(fabs(center[1]),2.0)-0.25);
		   // value =  1.0 - pow(center[1],2.0);
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

// Value of the velocity of the BMP model at the boundary
real get_boundary_velocity2(int id, Point center, int dim, real t) {
    real value;
    //BMP model parameter values
    real Re = 0.1;
    real GM = 0.5;
    real LA = 0.2;
    real Phi = 0.1;
    //real GP = 13.13806962;
    real GP = 26.8324;
    switch (id) {
        case 0:
            switch (dim) {
                case 0:
                    //set max velocity = 8.0e-4 
                    //value = 3.2e-3*(-center[1]*center[1] + center[1]);
                    //value = 4.0*(-center[1]*center[1] + 0.25);
                    //set max velocity = 1.5 
                    value = get_velocity_BMP(center, Re, GM, LA, Phi, GP);
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

//Calculates the error
real Calculaerro(higflow_solver *ns, int myrank, real new, real old) {

    real erro= 0.0;
    erro = fabs((old-new)/new)*100.0;

    return erro;
}

// calculates u the velocity
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

//calculates the viscosity
real calc_eta(higflow_solver *ns, int myrank, int dim, real Px, real Py) {
        // Get the local sub-domain
        sim_domain *sdp = psd_get_local_domain(ns->ed.vevv.psdVisc);
        // Get the local partitioned domain for facets
        sim_facet_domain *sfdu[DIM];
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        Point P;
	    P[0] = Px;
        P[1] = Py;
        real eta = compute_value_at_point(ns->ed.vevv.sdVisc, P, P, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
    return eta;
}

// calculates the elastic stress tensor
real calc_tau(higflow_solver *ns, int myrank, int i, int j, real Px, real Py) {
        // Get the constants
        real Re    = ns->par.Re;
        real De    = ns->ed.vevv.par.De;
        real beta  = ns->ed.vevv.par.beta;
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        Point P;
	    P[0] = Px;
        P[1] = Py;
        // Get the velocity derivative tensor Du and the Kernel tensor
        real Du[DIM][DIM];
        // Get Du
        Du[i][j] = compute_value_at_point(ns->ed.sdED, P, P, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
        Du[j][i] = compute_value_at_point(ns->ed.sdED, P, P, 1.0, ns->ed.vevv.dpD[j][i], ns->ed.stn);
        // Get T tensor
        real D  = 0.5*(Du[i][j]+Du[j][i]);
        real S = compute_value_at_point(ns->ed.sdED, P, P, 1.0, ns->ed.vevv.dpS[i][j], ns->ed.stn);
        //Get the viscosity value
        //real eta = compute_value_at_point(ns->ed.vevv.sdVisc, P, P, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
        //real tau = S + 2.0*(1.0-beta)*eta*D/Re -2.0*beta/Re*D;
        real tau = S + 2.0*(1.0-beta)*(D/Re);
        
    return tau;
}

// Print the velocity
void print_velocity(higflow_solver *ns, int myrank, int dim, real x, int np, real yf, real yl, real dy) {
    char filename[1024];
    //snprintf(filename,sizeof filename,"Profiles/Velocities/velocity_%d_%d_%d.dat",dim,myrank,ns->par.frame);
    snprintf(filename,sizeof filename,"Profiles/Velocities/velocity_%d_%d_%d_%d.dat",dim,myrank,np,ns->par.frame);
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
            fprintf(fd,"%16.10lf %16.10lf %16.10lf %16.10lf \n",y,x,ns->par.t,u);
        }
        // Destroy the iterator
        fclose(fd);
    } else {
        printf("File %s could not be opened \n",filename);
        exit(1);
    }
}

// Print the velocity
void print_viscosity(higflow_solver *ns, int myrank, real x, int np, real yf, real yl, real dy) {
    char filename[1024];
    snprintf(filename,sizeof filename,"Profiles/Viscosities/viscosity_%d_%d_%d.dat",myrank,np,ns->par.frame);
    FILE *fd = fopen(filename, "w");
    if (fd != NULL) {
        // Get the local sub-domain
        sim_domain *sdp = psd_get_local_domain(ns->ed.vevv.psdVisc);
        Point P;
	    P[0] = x;
	    for (real y = yf; y <= yl; y += dy) {
            P[1] = y;
            real eta = compute_value_at_point(ns->ed.vevv.sdVisc, P, P, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            fprintf(fd,"%16.10lf %16.10lf %16.10lf %16.10lf \n",y,x,ns->par.t,eta);
        }
        // Destroy the iterator
        fclose(fd);
    } else {
        printf("File %s could not be opened \n",filename);
        exit(1);
    }
}


// Print the Polymeric Tensor
real print_tensor(higflow_solver *ns, int myrank, int i, int j, real x, int np, real yf, real yl, real dy) {
    char filename[1024];
    //snprintf(filename,sizeof filename,"Profiles/Tensors/T_%d_%d_%d_%d.dat",i,j,myrank,ns->par.frame);
    snprintf(filename,sizeof filename,"Profiles/Tensors/T_%d_%d_%d_%d_%d.dat",i,j,myrank,np,ns->par.frame);
    FILE *fd = fopen(filename, "w");
    //real max=0.0;
     if (fd != NULL) {
        // Get the constants
        real Re    = ns->par.Re;
        real De    = ns->ed.vevv.par.De;
        real beta  = ns->ed.vevv.par.beta;
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        Point P;
	    P[0] = x;
        for (real y = yf; y <= yl; y += dy) {
            P[1] = y;
            // Get the velocity derivative tensor Du and the Kernel tensor
            real Du[DIM][DIM];
            // Get Du
            Du[i][j] = compute_value_at_point(ns->ed.sdED, P, P, 1.0, ns->ed.vevv.dpD[i][j], ns->ed.stn);
            Du[j][i] = compute_value_at_point(ns->ed.sdED, P, P, 1.0, ns->ed.vevv.dpD[j][i], ns->ed.stn);
            // Get T tensor
            real D  = 0.5*(Du[i][j]+Du[j][i]);
            real S = compute_value_at_point(ns->ed.sdED, P, P, 1.0, ns->ed.vevv.dpS[i][j], ns->ed.stn);
            //Get the viscosity value
            //real eta = compute_value_at_point(ns->ed.vevv.sdVisc, P, P, 1.0, ns->ed.vevv.dpvisc, ns->ed.stn);
            //real T = S + 2.0*(1.0-beta)*eta*D/Re -(2.0*beta/Re)*D; 
            real T = S + 2.0*(1.0-beta)*D/Re;
	        fprintf(fd,"%16.10lf %16.10lf %16.10lf %16.10lf \n",y,x,ns->par.t,T);
	    }
        // Loop for each cell
        fclose(fd);
    } else {
        printf("File %s could not be opened\n",filename);
        exit(1);
    }
   //return max; 
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
    // Load the controllers data for Navier-Stokes simulation
    higflow_load_all_controllers(ns, myrank);
    // Load the parameters data for Navier-Stokes simulation
    higflow_load_all_parameters(ns, myrank);
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
    //higflow_create_domain_viscoelastic(ns, cache, order_center,get_tensor, 
	//	    get_kernel, get_kernel_inverse, get_kernel_jacobian); 
    higflow_create_domain_viscoelastic_variable_viscosity(ns,cache,order_center,get_tensor, 
		    get_kernel, get_kernel_inverse, get_kernel_jacobian,get_viscosity,get_structpar);
    //higflow_create_domain_viscoelastic_integral(ns, cache, order_center, get_tensor); 
    // Initialize the domain
    higflow_initialize_domain(ns, ntasks, myrank, order_facet); 
    // Load the controllers data for viscoelastic simulation
    // higflow_load_non_newtonian_controllers(ns, myrank);
    // //higflow_load_viscoelastic_controllers(ns, myrank);
    // higflow_load_viscoelastic_variable_viscosity_controllers(ns,myrank);
    // // Load the parameters data for viscoelastic simulation
    // //higflow_load_viscoelastic_parameters(ns, myrank);
    // higflow_load_viscoelastic_variable_viscosity_parameters(ns,myrank);
    // // Load the controllers data for viscoelastic simulation
    // //higflow_load_viscoelastic_integral_controllers(ns, myrank);
    // // Load the parameters data for viscoelastic simulation
    // //higflow_load_viscoelastic_integral_parameters(ns, myrank);
    // // Set the user model
    // //higflow_define_user_function_viscoelastic(ns, calculate_m_user);
    // Initialize the boundaries
    higflow_initialize_boundaries(ns); 
    // Creating distributed property  
    higflow_create_distributed_properties(ns);
    // Creating distributed property for generalized newtonian simulation
    //higflow_create_distributed_properties_generalized_newtonian(ns);
    //higflow_create_distributed_properties_viscoelastic(ns);
    //higflow_create_distributed_properties_viscoelastic_variable_viscosity(ns);
    //higflow_create_distributed_properties_viscoelastic_integral(ns);
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
        //real Ux = get_velocity_BMP_proof(0.4, 0.1, 0.5, 0.2, 0.1, 26.8324);
        //printf("===> Ux BMP model =%lf <===\n", Ux);
        // Print the step
        if (myrank == 0) 
            printf("===> Step:        %7d <====> t     = %15.10lf <===\n", step, ns->par.t);
            //printf("===> Step:        %7d <====> tdim  = %15.10lf seg <===\n", step, (ns->par.t)*0.06);
        // Start the first step time
        if (step == step0)  START_CLOCK(firstiter); 

            real Celc1_uold=calc_u(ns, myrank, 0, 1.0, 0.8);
            real Celc2_uold=calc_u(ns, myrank, 0, 9.5, 0.5);
            real Celc3_uold=calc_u(ns, myrank, 0, 9.0, 0.2);

            real Celc1_etaold=calc_eta(ns, myrank, 0, 1.0, 0.8);
            real Celc2_etaold=calc_eta(ns, myrank, 0, 9.5, 0.5);
            real Celc3_etaold=calc_eta(ns, myrank, 0, 9.0, 0.2);

            real Celc1_Txxold=calc_tau(ns, myrank, 0, 0, 0.1, 0.8);
            real Celc2_Txxold=calc_tau(ns, myrank, 0, 0, 9.5, 0.5);
            real Celc3_Txxold=calc_tau(ns, myrank, 0, 0, 9.0, 0.2);
         
            real Celc1_Txyold=calc_tau(ns, myrank, 1, 0, 0.1, 0.8);
            real Celc2_Txyold=calc_tau(ns, myrank, 1, 0, 9.5, 0.5);
            real Celc3_Txyold=calc_tau(ns, myrank, 1, 0, 9.0, 0.2);

        higflow_solver_step_viscoelastic_variable_viscosity(ns);
        //higflow_solver_step_viscoelastic_integral(ns);

                 ////////////////////////////////////////////////////////////////////////////////
         //Calculation of the errors of the velocity and the elastic stress tensor

         real Celc1_u=calc_u(ns, myrank, 0, 1.0, 0.8);
         real Celc2_u=calc_u(ns, myrank, 0, 9.5, 0.5);
         real Celc3_u=calc_u(ns, myrank, 0, 9.0, 0.2);

         real Celc1_eta=calc_eta(ns, myrank, 0, 1.0, 0.8);
         real Celc2_eta=calc_eta(ns, myrank, 0, 9.5, 0.5);
         real Celc3_eta=calc_eta(ns, myrank, 0, 9.0, 0.2);
         
         real Celc1_Txx=calc_tau(ns, myrank, 0, 0, 1.0, 0.8);
         real Celc2_Txx=calc_tau(ns, myrank, 0, 0, 9.5, 0.5);
         real Celc3_Txx=calc_tau(ns, myrank, 0, 0, 9.0, 0.2);
         
         real Celc1_Txy=calc_tau(ns, myrank, 1, 0, 1.0, 0.8);
         real Celc2_Txy=calc_tau(ns, myrank, 1, 0, 9.5, 0.5);
         real Celc3_Txy=calc_tau(ns, myrank, 1, 0, 9.0, 0.2);

        //////////////////////////////////////////////////////
         real erro_Celc1_u= Calculaerro(ns, myrank, Celc1_u, Celc1_uold);
         real erro_Celc2_u= Calculaerro(ns, myrank, Celc2_u, Celc2_uold);
         real erro_Celc3_u= Calculaerro(ns, myrank, Celc3_u, Celc3_uold);

         real erro_Celc1_eta= Calculaerro(ns, myrank, Celc1_eta, Celc1_etaold);
         real erro_Celc2_eta= Calculaerro(ns, myrank, Celc2_eta, Celc2_etaold);
         real erro_Celc3_eta= Calculaerro(ns, myrank, Celc3_eta, Celc3_etaold);
         
         real erro_Celc1_Txx= Calculaerro(ns, myrank, Celc1_Txx, Celc1_Txxold);
         real erro_Celc2_Txx= Calculaerro(ns, myrank, Celc2_Txx, Celc2_Txxold);
         real erro_Celc3_Txx= Calculaerro(ns, myrank, Celc3_Txx, Celc3_Txxold);
         
         real erro_Celc1_Txy= Calculaerro(ns, myrank, Celc1_Txy, Celc1_Txyold);
         real erro_Celc2_Txy= Calculaerro(ns, myrank, Celc2_Txy, Celc2_Txyold);
         real erro_Celc3_Txy= Calculaerro(ns, myrank, Celc3_Txy, Celc3_Txyold);
          
        //Calculation of the errors
         printf("===> Errors ===>   Ux  : Cel1 = %15.10lf <=> Cel2 = %15.10lf <=> Cel3 = %15.10lf  <===\n",erro_Celc1_u, erro_Celc2_u, erro_Celc3_u);
         printf("===> Errors ===>   eta  : Cel1 = %15.10lf <=> Cel2 = %15.10lf <=> Cel3 = %15.10lf  <===\n",erro_Celc1_eta, erro_Celc2_eta, erro_Celc3_eta);
         printf("===> Errors ===> Txx : Cel1 = %15.10lf <=> Cel2 = %15.10lf <=> Cel3 = %15.10lf  <===\n",erro_Celc1_Txx, erro_Celc2_Txx, erro_Celc3_Txx);
         printf("===> Errors ===> Txy : Cel1 = %15.10lf <=> Cel2 = %15.10lf <=> Cel3 = %15.10lf  <===\n",erro_Celc1_Txy, erro_Celc2_Txy, erro_Celc3_Txy);

        
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

            ///////////////////////////////////////////////////////////////////////////////////    
              real dx, dy;
              dx = 0.0625;
              dy = 0.01;
            ///////////////////////////////////////////////////////////////////////////////////    
            // Printing em y=0.5 para todo x
              //print_velocity_2(ns, myrank, 0, 0.5, 0.0, 10.0, dx);
              
              //print_tensor_2(ns, myrank, 0, 0, 0.5, 0.0, 10.0, dx);
              //print_tensor_2(ns, myrank, 0, 1, 0.5, 0.0, 10.0, dx);
              //print_tensor_2(ns, myrank, 1, 0, 0.5, 0.0, 10.0, dx);
              //print_tensor_2(ns, myrank, 1, 1, 0.5, 0.0, 10.0, dx);
          ///////////////////////////////////////////////////////////////////////////////////    
              // Printing in x=5 for every y
              real ymin = 0.0;
              real ymax = 1.0 + dy;
              real xfixed = 1.0;
              int np1 = 1;
              //print_velocity(ns, myrank, 0, xfixed, np1, ymin, ymax, dy);
              //print_viscosity(ns, myrank, xfixed, np1, ymin, ymax, dy);
              
              //print_tensor(ns, myrank, 0, 0, xfixed, np1, ymin, ymax, dy);
              //print_tensor(ns, myrank, 0, 1, xfixed, np1, ymin, ymax, dy);
              //print_tensor(ns, myrank, 1, 0, xfixed, np1, ymin, ymax, dy);
              //print_tensor(ns, myrank, 1, 1, xfixed, np1, ymin, ymax, dy);

              real x2 = 9.0;
              //real x4 = 0.4;  
              int np2 = 2;
              //int np4 = 4;
              print_velocity(ns, myrank, 0, x2, np2, ymin, ymax, dy);
              print_viscosity(ns, myrank, x2, np2, ymin, ymax, dy);

              //print_velocity(ns, myrank, 0, x4, np4, ymin, ymax, dy);
              //print_viscosity(ns, myrank, x4, np4, ymin, ymax, dy);
              
              print_tensor(ns, myrank, 0, 0, x2, np2, ymin, ymax, dy);
              print_tensor(ns, myrank, 0, 1, x2, np2, ymin, ymax, dy);
              //print_tensor(ns, myrank, 1, 0, x2, np2, ymin, ymax, dy);
              //print_tensor(ns, myrank, 1, 1, x2, np2, ymin, ymax, dy);

              //print_tensor(ns, myrank, 0, 0, x4, np4, ymin, ymax, dy);
              //print_tensor(ns, myrank, 0, 1, x4, np4, ymin, ymax, dy);
              //print_tensor(ns, myrank, 1, 0, x4, np4, ymin, ymax, dy);
              //print_tensor(ns, myrank, 1, 1, x4, np4, ymin, ymax, dy);
        ///////////////////////////////////////////////////////////////////////////////////  
            
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
    //printf("We are here");
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
