// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Step Viscoelastic Integral - version 20/11/2017
// *******************************************************************
// *******************************************************************

#include "hig-flow-step-viscoelastic-integral.h"

//*******************************************************
// Estruture for Mittag-Leffler function
//*******************************************************

#define pi 3.1415926535897932384626434
typedef struct{
  double real;
  double imag;
} numc;

numc kk(double r, double alfa, double beta1, numc z){
 double two = 2, um = 1;
 double A, B, C, D, E, norz;
 numc w1;
 A = pow(r,(um-beta1)/alfa)*exp(-pow(r,um/alfa))/(pi*alfa);
 B = r*sin(pi*(um-beta1));
 C = sin(pi*(um-beta1+alfa));
 D = pow(r, two);
 E = two*r*cos(pi*alfa);
 norz =  pow(z.real,2) + pow(z.imag,2); // norz= |z|^2
 
 w1.real = A*(B*D - (B*E + C*D)*z.real + pow(z.real,2)*(B+C*E-C*z.real) + pow(z.imag,2)*(-B + C*E - C*z.real))/(pow(D-E*z.real + pow(z.real,2) - pow(z.imag,2),2) + pow(-E*z.imag+2*z.real*z.imag,2));
 w1.imag = A*z.imag*(B*E - C*D - 2*B*z.real + C*norz)/(pow(D-E*z.real + pow(z.real,2) - pow(z.imag,2),2) + pow(-E*z.imag+2*z.real*z.imag,2));
return w1;
}

numc pp(double r, double alphaf, double betaf, numc z, double epsn){
 double two = 2, um = 1;
 double w, A;
 numc p;
 
 A = pow(epsn, um+(um-betaf)/alphaf)*exp(pow(epsn,(um/alphaf))*cos(r/alphaf))/(two*pi*alphaf);
 w= (pow(epsn,um/alphaf))*sin(r/alphaf)+r*(um+(um-betaf)/alphaf);

 p.real = A*(epsn*cos(w-r) - z.real*cos(w) -z.imag*sin(w))/(pow(epsn,2) -2*epsn*(z.real*cos(r) + z.imag*sin(r)) +z.real*z.real + z.imag*z.imag );
 p.imag = A*(epsn*sin(w-r) - z.real*sin(w) +z.imag*cos(w))/(pow((epsn*cos(r)-z.real),2) + pow((epsn*sin(r)-z.imag),2)); 
 // printf("%f\n",w);
return p;
} 

/* Romberg Integration*/

numc rombint(char funfcn, double a, double b, int order, double v1, double v2, numc v3, double v4){
	int iorder, ipower,i,j,k;
	double hh;
	double zero = 0.0;
	double two = 2, um = 1, six = 6, pcin = 0.5;
	numc a1,a2,auxsum;
	numc rom[3][9];
	
	iorder=order;

	if(funfcn == 'K') { 
		iorder = six;
	  //printf("iorder = %d\n", iorder); //cout << "iorder = 6!" << endl;
	}
	if(order > 8) {
		//printf("need to increase size of matrix r which is 8, order is %d\n", iorder);
  		//cout << "need to increase size of matrix r which is 8, order is" <<  iorder << endl;
  		/*break;*/
	}
	for (i=1; i<=2; i++) {
  		for (j=1; j<=iorder; j++) {
    		rom[i][j].real=zero;
    		rom[i][j].imag=zero;
  		}
	}
	hh=b-a;
	//printf("hh = %f\n ", hh);  //cout << "hh=" << hh << endl;
	if(funfcn == 'K') {
  		a1=kk(a,v1,v2,v3);
  		a2=kk(b,v1,v2,v3);
	}else{
  		a1=pp(a,v1,v2,v3,v4);
  		a2=pp(b,v1,v2,v3,v4);
	}
        
	rom[1][1].real=hh*(a1.real+a2.real)/two;
	rom[1][1].imag=hh*(a1.imag+a2.imag)/two;
	
	ipower=1;
	for (i=2; i<=iorder; i++) {
  		auxsum.real=zero;
		auxsum.imag=zero;
  		for (j=1;j<=ipower;j++) {

    			if(funfcn == 'K') {
      				a1=kk((a+hh*(j-pcin)),v1,v2,v3);
    			}else{
      				a1=pp((a+hh*(j-pcin)),v1,v2,v3,v4);
    			}
    			auxsum.real=auxsum.real+a1.real;
			auxsum.imag=auxsum.imag+a1.imag;
  		}
  		rom[2][1].real=(rom[1][1].real+hh*auxsum.real)/two;
		rom[2][1].imag=(rom[1][1].imag+hh*auxsum.imag)/two;
  		for (k=1; k<=(i-1); k++) {
    			rom[2][k+1].real=((pow(4,k))*rom[2][k].real-rom[1][k].real)/((pow(4,k))-1);
			rom[2][k+1].imag=((pow(4,k))*rom[2][k].imag-rom[1][k].imag)/((pow(4,k))-1);
  		}
  		for (j=0; j<=(i-1); j++) {
    			rom[1][j+1].real=rom[2][j+1].real;
			rom[1][j+1].imag=rom[2][j+1].imag;
  		}
  		ipower=ipower*2;
  		hh=hh/two;
	}
	//printf("rom[1][%d]=%f + i %f \n",iorder,rom[1][iorder].real, rom[1][iorder].imag );
	//cout << "rom[1][iorder]="<< rom[1][iorder] << endl;	
 
return rom[1][iorder];
}

// *******************************************************************
// Calculate the value of the Mittag-Leffler function
// *******************************************************************


numc mlfv(double alpha, double beta, numc z, int fi){
	int k,m,h;
	//double pi = 3.1415926535897932384626434;
	double two = 2, um = 1, pcin = 0.5, ten = 10, cinc=5, six = 6, zero = 0.0, dum=2.1, ccinc=5.5, vint=20, cinq=50;
	double r0, rc, angz, eps, aaz, norz, pr1, pr2, pq1, pq2;
	int l1,l2,l3,l4;
	numc newsum, term, aux, aux1, a1, a2, oldsum, zn, ic, res, raiz1z, raiz2z;
	char PN = 'P';
	char KN = 'K';
	ic.real = 0;
	ic.imag = 1;

	//printf("z=%f +i %f\n",z.real, z.imag); //cout << "z=" << z << endl;
	//printf("alpha = %f\n",alpha); //cout << "alpha=" << alpha << endl;

	a1.real=zero;
	a1.imag=zero;
	a2.real=zero;
	a2.imag=zero;
	newsum.real=zero; 
	newsum.imag=zero; 
	oldsum.real=zero; 
	oldsum.imag=zero; 
	res.real=zero;
	res.imag=zero;

	if (alpha < 0) {
		printf("alpha must be greater than zero\n"); //cout << "alpha must be greater than zero" << endl;
	} else if (alpha > 5) {
	printf("alpha must be smaller than five\n"); //cout << "alpha must be smaller than five" << endl;
	}
	if (beta > 5) {
	printf("beta must be smaller than zero\n");//cout << "beta must be smaller than zero" << endl;
	} else {
	 // printf("beta is Ok!\n"); //cout << "beta is OK!" << endl;
	}

	rc=pow((-two*log((pow(ten,-fi))*pi/six)),alpha);
	norz = sqrt(pow(z.real,2) + pow(z.imag,2));
	r0=fmax (um, two*norz);
	r0=fmax (r0, rc);
	// printf("two*abs(z)=%f\n", two*norz); //cout << "two*abs(z)=" << two*abs(z) << endl;
	// printf("rc = %f\n", rc);//cout << "rc=" << rc << endl;
	// printf("r0 = %f\n", r0);//cout << "r0=" << r0 << endl;
	angz=atan2 (z.imag, z.real);  //angz=atan2 (imag(z), real(z));
	aaz=fabs(angz);
	
	// printf("aaz= %f\n", aaz); //cout << "aaz=" << aaz << endl;
	if (alpha == 1 && beta == 1) {
		res.real=exp(z.real)*cos(z.imag);
		res.imag=exp(z.real)*sin(z.imag);
	}

	l1=(alpha < 1 && norz <= 1);
	// printf("abs(z) = %f\n", norz); //cout << "abs(z)=" << abs(z) << endl;
	l2=(1 <= alpha && alpha < two);
	l3=(norz <= floor (vint/pow(dum-alpha, ccinc-two*alpha)));
	// printf("teste: %f\n", floor (vint/pow(dum-alpha, ccinc-two*alpha))); //cout << "teste" << floor (vint/pow(dum-alpha, ccinc-two*alpha)) << endl;
	l4=(alpha >= two && norz <= cinq);
	// printf("l1= %d\n", l1); //cout << "l1=" << l1 << endl;
	// printf("l2= %d\n", l2); // cout << "l2=" << l2 << endl;
	// printf("l3= %d\n", l3); //cout << "l3=" << l3 << endl;
	// printf("l4= %d\n", l4); //cout << "l4=" << l4 << endl;

	if (l1 || l2 && l3 || l4) {
  		 /*complex<double> oldsum(0.0, 0.0);
  		oldsum=dcmplx(0d0,0d0);*/
		oldsum.real=zero;
		oldsum.imag=zero; 	
		// printf("oldsum = %f+ i %f\n", oldsum.real, oldsum.imag); //cout << "oldsum=" << oldsum << endl;
		k=0;
  		while ((alpha*k+beta) <= zero) {
    			k=k+1;
  		}
  		/*newsum=pow(real(z),k/tgamma(alpha*k+beta));*/
		// printf("k= %d\n", k);//cout << "k=" << k << endl;
		// printf("z= %f+ i %f\n",z.real, z.imag);//cout << "z=" << z << endl;
		// printf("pow(z,k)= %f + i %f\n", pow(norz,k)*cos(angz*k), pow(norz,k)*sin(angz*k));//cout << "pow(z,k)=" << pow(z,k) << endl;
  		newsum.real=pow(norz,k)*cos(angz*k)/tgamma(alpha*k+beta); //pow(z,k)/tgamma(alpha*k+beta);
  		newsum.imag=pow(norz,k)*sin(angz*k)/tgamma(alpha*k+beta);//pow(z,k)/tgamma(alpha*k+beta);

		/* double summation because z can be negative */
  		while (newsum.real != oldsum.real || newsum.imag != oldsum.imag) {//while (newsum != oldsum) {
    			oldsum.real=newsum.real;
			oldsum.imag=newsum.imag;
   			k=k+1;
    			term.real=pow(norz,k)*cos(angz*k)/tgamma(alpha*k+beta); //term=pow(z,k)/tgamma(alpha*k+beta);
			term.imag=pow(norz,k)*sin(angz*k)/tgamma(alpha*k+beta);	
    			newsum.real=newsum.real+term.real;//newsum=newsum+term;
			newsum.imag=newsum.imag+term.imag;
    			k=k+1;
    			term.real=pow(norz,k)*cos(angz*k)/tgamma(alpha*k+beta); //term=pow(z,k)/tgamma(alpha*k+beta);
			term.imag=pow(norz,k)*sin(angz*k)/tgamma(alpha*k+beta);
    			newsum.real=newsum.real+term.real;//newsum=newsum+term;
			newsum.imag=newsum.imag+term.imag;
  		}
  		res.real=newsum.real;
		res.imag=newsum.imag;
		// printf("newsummmmm= %f + i %f\n",newsum.real, newsum.imag); //cout << "newsummmmm=" << newsum << endl;
	return res;
	}

/*! the matlab function fix rounds toward zero, can  use floor since alpha is positive*/
	if (alpha <= 1 && norz <= floor(cinc*alpha+ten)) {
   		if ((aaz > pi*alpha) && (fabs(aaz-(pi*alpha)) > pow(ten,(-fi)))) {
    			if (beta <= 1) {
    				res=rombint(KN,zero,r0,fi,alpha,beta,z,zero);
				//printf("res1=%f+ i%f\n",res.real, res.imag);//cout << "res1=" << res<< endl;
   			 } else {

      				eps=1;
      				aux=rombint(PN,-pi*alpha,pi*alpha,fi,alpha,beta,z,eps);
      				res=rombint(KN,eps,r0,fi,alpha,beta,z,zero);//+aux;
				res.real = res.real + aux.real;
				res.imag = res.imag + aux.imag;
				// printf("res2=%f+ i%f\n",res.real, res.imag);//cout << "res2=" << res<< endl;
			}

  		} else if (aaz < pi*alpha && fabs(aaz-(pi*alpha)) > pow(ten,(-fi))) {
    				if (beta <= 1) {
      					//aux=(pow(z,((1-beta)/alpha)))*(exp(pow(z,1/alpha))/alpha);
					aux.real = (pow(norz,((1-beta)/alpha))/alpha)*exp(pow(norz,(1/alpha))*cos(angz/alpha))*( cos(angz*(1-beta)/alpha)*cos(pow(norz,(1/alpha))*sin(angz/alpha)) - sin(angz*(1-beta)/alpha)*sin(pow(norz,(1/alpha))*sin(angz/alpha))); //cos(angz*(1-beta)/alpha);
					aux.imag = (pow(norz,((1-beta)/alpha))/alpha)*exp(pow(norz,(1/alpha))*cos(angz/alpha))*( sin(angz*(1-beta)/alpha)*cos(pow(norz,(1/alpha))*sin(angz/alpha)) - cos(angz*(1-beta)/alpha)*sin(pow(norz,(1/alpha))*sin(angz/alpha)));
					//cout << "pow(z,((1-beta)/alpha))=" << pow(z,((1-beta)/alpha)) << endl;
					pr1 = pow(norz,((1-beta)/alpha))*cos(angz*(1-beta)/alpha);
					pr2 = pow(norz,((1-beta)/alpha))*sin(angz*(1-beta)/alpha);
					// printf("pow(z,((1-beta)/alpha))= %f+ i %f\n",pr1, pr2);
					//cout << "exp(pow(z,1/alpha))/alpha=" << exp(pow(z,1/alpha))/alpha << endl;
					pq1= exp(pow(norz,(1/alpha))*cos(angz/alpha))*cos(pow(norz,(1/alpha))*sin(angz/alpha));
					pq2 = exp(pow(norz,(1/alpha))*cos(angz/alpha))*sin(pow(norz,(1/alpha))*sin(angz/alpha));
					// printf("exp(pow(z,1/alpha))/alpha= %f+ i %f\n",pq1, pq2);
					//cout << "aux=" << aux << endl;
					//printf("aux = %f+ i%f\n",aux.real,aux.imag);
     					res=rombint(KN,zero,r0,fi,alpha,beta,z,zero);//+aux;
					res.real = res.real + aux.real;
					res.imag = res.imag + aux.imag;
					//cout << "res3=" << res<< endl;
					// printf("res3=%f+i %f\n", res.real, res.imag);
				} else {
      					eps=norz/two;//eps=abs(z)/two;
     					aux1=rombint(KN,eps,r0,fi,alpha,beta,z,zero);
					aux=rombint(PN,-pi*alpha,pi*alpha,fi,alpha,beta,z,eps);//+aux1;
					aux.real = aux.real+ aux1.real;
					aux.imag = aux.imag+ aux1.imag;
      					//res=aux+(pow(z,((um-beta)/alpha)))*(exp (pow(z,um/alpha))/alpha);
					res.real=aux.real+ (pow(norz,((1-beta)/alpha))/alpha)*exp(pow(norz,(1/alpha))*cos(angz/alpha))*( cos(angz*(1-beta)/alpha)*cos(pow(norz,(1/alpha))*sin(angz/alpha)) - sin(angz*(1-beta)/alpha)*sin(pow(norz,(1/alpha))*sin(angz/alpha)));
					res.imag = aux.imag + (pow(norz,((1-beta)/alpha))/alpha)*exp(pow(norz,(1/alpha))*cos(angz/alpha))*( sin(angz*(1-beta)/alpha)*cos(pow(norz,(1/alpha))*sin(angz/alpha)) - cos(angz*(1-beta)/alpha)*sin(pow(norz,(1/alpha))*sin(angz/alpha)));
					//cout << "res4=" << res<< endl;
					// printf("res4=%f+i %f\n", res.real, res.imag);
    				}
  		} else {
    			eps=norz+pcin;
    			aux=rombint(PN,-pi*alpha,pi*alpha,fi,alpha,beta,z,eps);
    			res=rombint(KN,eps,r0,fi,alpha,beta,z,zero);//+aux;
			res.real = res.real + aux.real;
			res.imag = res.imag + aux.imag;
			//cout << "res5=" << res<< endl;
			//printf("res4=%f+i %f\n", res.real, res.imag);
  		}
	//cout << "res6=" << res<< endl;
	// printf("res6=%f+i %f\n", res.real, res.imag);
	return res;

	}

	if (alpha <= 1) { 

 	 	if (aaz<(pi*alpha/two+fmin(pi,pi*alpha))/two) { 
    			//newsum=(pow(z,((um-beta)/alpha)))*exp(pow(z,(um/alpha)))/alpha;
			newsum.real = (pow(norz,((1-beta)/alpha))/alpha)*exp(pow(norz,(1/alpha))*cos(angz/alpha))*( cos(angz*(1-beta)/alpha)*cos(pow(norz,(1/alpha))*sin(angz/alpha)) - sin(angz*(1-beta)/alpha)*sin(pow(norz,(1/alpha))*sin(angz/alpha)));
			newsum.imag = (pow(norz,((1-beta)/alpha))/alpha)*exp(pow(norz,(1/alpha))*cos(angz/alpha))*( sin(angz*(1-beta)/alpha)*cos(pow(norz,(1/alpha))*sin(angz/alpha)) - cos(angz*(1-beta)/alpha)*sin(pow(norz,(1/alpha))*sin(angz/alpha)));

    			for (k=1; k<=floor(fi/log10(norz));k++) {
       				if (ceil(beta-alpha*k)!=floor(beta-alpha*k)) {
         			//newsum=newsum-((pow(z,(-k)))/tgamma(beta-alpha*k));
				newsum.real = newsum.real - ((pow(norz,-k)*cos(angz*k))/tgamma(beta-alpha*k));
				newsum.imag = newsum.imag - (-(pow(norz,-k)*sin(angz*k))/tgamma(beta-alpha*k));
				//cout << "passei aqui 1  " << gamma <<endl;
				// printf("passei aqui 1: \n ");
				//cout << "passei aqui 1  " << k <<endl;
				// printf("passei aqui 1: %d\n", k);
       				}
    			}
    			//res=newsum;
			res.real=newsum.real;
			res.imag=newsum.imag;
			//cout << "res6=" << res<< endl;
			// printf("res6=%f+i %f\n", res.real, res.imag);
  		}else{ 
    			//newsum=zero;
			newsum.real =zero;
			newsum.imag =zero;
    			for (k=1;k<=floor(fi/log10(norz));k++) {
       				if (ceil(beta-alpha*k) != floor(beta-alpha*k)) {
          			//newsum=newsum-((pow(z,(-k)))/tgamma(beta-alpha*k));
				newsum.real = newsum.real - ((pow(norz,-k)*cos(angz*k))/tgamma(beta-alpha*k));
				newsum.imag = newsum.imag - (-(pow(norz,-k)*sin(angz*k))/tgamma(beta-alpha*k));
				//cout << "passei aqui 2" << gamma(beta-alpha*k) <<endl;
				// printf("passei aqui 2: %f \n ", tgamma(beta-alpha*k));
       				}
    			}
    			//res=newsum;
			res.real=newsum.real;
			res.imag=newsum.imag;
			//cout << "res7=" << res<< endl;
			// printf("res7=%f+i %f\n", res.real, res.imag);
  		}
	}else{
 		 if (alpha >= two) {
   			m=floor(alpha/two);
    		 	//aux=zero;
			aux.real = zero;
			aux.imag = zero;
   			for (h=0;h<=m;h++) {
      				//zn=(pow(z,um/(m+um)))*exp((h*two*pi*ic)/(m+um));
				zn.real = pow(norz,um/(m+um))*cos((2*pi*h+angz)/(m+um));
				zn.imag = pow(norz,um/(m+um))*sin((2*pi*h+angz)/(m+um));
      				// aux=aux+mlfv(alpha/(m+1),beta,zn,fi);
				aux1=mlfv(alpha/(m+1),beta,zn,fi);
				aux.real= aux.real + aux1.real;
				aux.imag = aux.imag + aux1.imag;
				/*cout << "conta=" << (2.0*pi*h)/(m+1) << endl;*/     
   			}
    			//res=(um/(m+um))*aux;
			res.real = (um/(m+um))*aux.real;
			res.imag = (um/(m+um))*aux.imag;
			//cout << "res8=" << res<< endl;
			// printf("res8=%f+i %f\n", res.real, res.imag);
  		}else{
			raiz1z.real = sqrt(norz)*cos(angz/2);
			raiz1z.imag = sqrt(norz)*sin(angz/2);
			raiz2z.real = -sqrt(norz)*cos(angz/2);
			raiz2z.imag = -sqrt(norz)*sin(angz/2);
			a1=mlfv(alpha/two,beta,raiz1z,fi);
    			a2=mlfv(alpha/two,beta,raiz2z,fi);
    			//res=(a1+a2)/two;
			res.real = (a1.real+a2.real)/two;
			res.imag = (a1.imag+a2.imag)/two;
			//cout << "res9=" << res<< endl;
			// printf("res9=%f+i %f\n", res.real, res.imag);
  		}
	return res;
	//cout << "res10=" << res<< endl;
	// printf("res10=%f+i %f\n", res.real, res.imag);
	}

return res;
//cout << "res11=" << res<< endl;
// printf("res10=%f+i %f\n", res.real, res.imag);

}



// *******************************************************************
// Constitutive Equations
// *******************************************************************


// *******************************************************************
// Constitutive Equation Step for the Explicit Euler Method
// *******************************************************************
void higflow_explicit_euler_constitutive_equation_integral(higflow_solver *ns) {
    if (ns->contr.flowtype == 4) {
        // Get the cosntants
        real Re    = ns->par.Re;
        real De    = ns->ed.im.par.De;
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int i = 0; i < DIM; i++) {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdp);
        // Looping for each Finger Tensor
        for (int k = 0; k <= NDT; k++) {
            // Loop for each cell
            higcit_celliterator *it;
            for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                // Get the cell
                hig_cell *c = higcit_getcell(it);
                // Get the cell identifier
                int clid    = mp_lookup(mp, hig_get_cid(c));
                // Get the center of the cell
                Point ccenter;
                hig_get_center(c, ccenter);
                // Get the delta of the cell
                Point cdelta;
                hig_get_delta(c, cdelta);
                // Get the velocity derivative tensor Du, S and B tensor
                real Du[DIM][DIM], B[DIM][DIM];
                for (int i = 0; i < DIM; i++) {
                    for (int j = 0; j < DIM; j++) {
                        // Get Du
                        Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
                        B[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.im.dpB[k][i][j], ns->ed.stn);
                    }
                }
                // Calculate RHS = Omega Kernel - Kernel Omega + 2BB + MM/De
                // Get the velocity at cell center 
                real u[DIM], dBdx[DIM], RHS[DIM][DIM];
                hig_flow_velocity_at_center_cell_integral(ns, ccenter, cdelta, u);
                // Solving the Constitutive Equation using the Euler Method
                hig_flow_b_rhs(B, Du, RHS);
                for (int i = 0; i < DIM; i++) {
                    for (int j = i; j < DIM; j++) {
                        // Right hand side equation
                        real rhs = 0.0;
                        switch (ns->ed.im.contr.convecdiscrtype) {
                            case 0: 
                                // Kernel derivative at cell center
                                hig_flow_derivative_b_at_center_cell(ns, ccenter, cdelta, k, i, j, B[i][j], dBdx);
                                for (int dim = 0; dim < DIM; dim++) {
                                    //Compute convective tensor term in rhs
                                    rhs -= u[dim]*dBdx[dim];
                                }
                                break;
                            case 1: 
                                //Compute convective tensor term CUBISTA in rhs
                                for (int dim = 0; dim < DIM; dim++) {
                                    rhs -= hig_flow_convective_tensor_term_b_cubista(ns, ns->dpu[dim], ns->ed.sdED, ns->ed.stn, B, ccenter, cdelta, dim, k, i, j);
                                }
                                break;
                        }
                        // Compute the final rhs
                        rhs    += RHS[i][j];
                        // Compute the Kernel at next time
                        real b  = B[i][j] + ns->par.dt * rhs;
                        // Store Kernel in S
                        dp_set_value(ns->ed.im.dpS[i][j], clid, b);
                        if (i != j) {
                             dp_set_value(ns->ed.im.dpS[j][i], clid, b);
                        }
                    }
                }
            }
            // Destroy the iterator
            higcit_destroy(it);
            // Sync the ditributed pressure property
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    dp_sync(ns->ed.im.dpS[i][j]);
                }
            }
            // Store the Kernel Tensor
            for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
                // Get the cell
                hig_cell *c = higcit_getcell(it);
                // Get the cell identifier
                int clid    = mp_lookup(mp, hig_get_cid(c));
                // Get the center of the cell
                Point ccenter;
                hig_get_center(c, ccenter);
                // Get the delta of the cell
                Point cdelta;
                hig_get_delta(c, cdelta);
                real B[DIM][DIM];
                for (int i = 0; i < DIM; i++) {
                    for (int j = 0; j < DIM; j++) {
                        // Get S
                        B[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
                    }
                }
                // Get the B tensor and store in Finger Tensor
                for (int i = 0; i < DIM; i++) {
                    for (int j = 0; j < DIM; j++) {
                        // Get S
                        B[i][j]  = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
                         //if ( B[i][j] >= 5.0)
                        // printf("B[%d][%d][%d] = %f \n",i,j,k, B[i][j]);
                        // Store Kernel
                        dp_set_value(ns->ed.im.dpB[k][i][j], clid, B[i][j]);
                    }
                }
            }
            // Destroy the iterator
            higcit_destroy(it);
            // Sync the ditributed pressure property
            for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    dp_sync(ns->ed.im.dpB[k][i][j]);
                }
            }
        }
     }
}

// Calculate the integral equation 
void hig_flow_integral_equation (higflow_solver *ns) {
    if (ns->ed.im.contr.model == 0) {
        hig_flow_integral_equation_KBKZ(ns);
    } else if (ns->ed.im.contr.model == 1) {
        hig_flow_integral_equation_KBKZ_Fractional(ns);
    }
}


// Calculate the integral equation for KBKZ model
void hig_flow_integral_equation_KBKZ (higflow_solver *ns) {
    real tensao[DIM][DIM], Bold[NDT+1][DIM][DIM], B[NDT+1][DIM][DIM];
    real I1a, I2a, I1b, I2b,I1c, I2c, welambda, alpha, beta;
    real A0, A1, A2, b0, b1, b2, ds, c1, c2;
    real De, Re;
    real a[8], lambda[8];
    real STensor; 
    int  M, Nnow;
    
    real s[NDT+1];
    real sold[NDT+1];
    real t, dt;
    real Q, ds1;

    real Smax[DIM][DIM], Smin[DIM][DIM];
    for (int i = 0; i < DIM; i++) {
       for (int j = 0; j < DIM; j++) {
           Smax[i][j]   = -1.0e16;
           Smin[i][j]   =  1.0e16;
           tensao[i][j] =  0.0;
       }
    }
    
    dt    = ns->par.dt;
    t     = ns->par.t;
    beta  = ns->ed.im.par.beta;
    alpha = ns->ed.im.par.alpha;
    M     = ns->ed.im.par.M;
    De    = ns->ed.im.par.De;
    Re    = ns->par.Re;
    for (int i = 0; i < M; i++) {
       a[i]      = ns->ed.im.par.a[i];
       lambda[i] = ns->ed.im.par.lambda[i];
    }  
    
    //divide the interval [0,t] into (NDT+1) points  (vectors s[k], sold[k])
    Nnow = NDT;
    //real translad;
    if (t > 0.0) { 
     if ( ns->par.step+1 <= Nnow){
       Nnow = ns->par.step+1;
       s[0] = 0.0;
       s[Nnow] = t;
       for (int k = 1; k <= Nnow; k++) {
          s[k]    = s[k-1] + dt;
          sold[k] = s[k];
      //printf(" s[%d] = %f, t= %f, Nnow = %d \n", k, s[k], t, Nnow); 
     }
    } else{
       for (int k = 0; k <= Nnow; k++) {
          sold[k] = ns->ed.im.s[k];
         //printf(" sold[%d] = %f ", k, sold[k]);
       }
       s[0]     = 0.0;
       s[Nnow]  = t;
       Q        = pow((t/dt),1.0/(NDT));  
     //  ds1      = (dt*pow(Q,1));
     //  translad = 0.98*(ds1);
       ds       = 0.0;
       for (int k = 1; k < NDT; k++) {
          ds        = (dt*pow(Q,k));
          s[NDT-k] = t-ds;
       //   s[NDT-k] = t-ds+translad;
          if (ds <= 0.0) {
             printf("ERROR: ds= s[m][%d] - s[m][%d] = %f <= 0!\n",k, k-1, ds);
             exit(1);
          }
          if (s[NDT-k] < 0.0) {
             printf("ERROR: s[%d] = %f ", NDT-k, s[NDT-k]);
             exit(1);
          }
          //printf(" s[%d] = %f ", NDT-k, s[NDT-k]);  
       }  
    }
    for (int k = 0; k <= Nnow; k++) {
      //printf("s[%d] = %f \n",k,s[k]);
     ns->ed.im.s[k] = s[k];
    }
   }
   
    // Get the local sub-domain for the cells
    sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
    // Get the map for the domain properties
    mp_mapper *mp = sd_get_domain_mapper(sdp);
    // Loop for each cell
    higcit_celliterator *it;
    if (ns->par.step+1 >= NDT) {
      for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Get the cell
        hig_cell *c = higcit_getcell(it);
        // Get the cell identifier
        int clid    = mp_lookup(mp, hig_get_cid(c));
        // Get the center of the cell
        Point ccenter;
        hig_get_center(c, ccenter);
        // Get the delta of the cell
        Point cdelta;
        hig_get_delta(c, cdelta);
        // Loop for each Finger Tensor
        for (int k = 0; k <= Nnow; k++) {
           for (int i = 0; i < DIM; i++) {
              for (int j = 0; j < DIM; j++) {
                  Bold[k][i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.im.dpB[k][i][j], ns->ed.stn);
                  B[k][i][j]    = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.im.dpB[k][i][j], ns->ed.stn);
              }
           } 
        } 
    
        ////////// interpolate the finger tensor at the new points s[k] ////////////////// 
        int  cont=0;
        real aux[DIM][DIM];
        real s0, s1, s2, sms0, sms0sms1, sms0sms12;
        for (int i = 0; i < DIM; i++) {
           for (int j = 0; j < DIM; j++) {
              if (i == j){
                 aux[i][j] = 1.0;
              } else {
                 aux[i][j] = 0.0;
              }
           }
        }  
        sold[NDT+1]=t; 
        for (int k = 0; k < NDT; k++) {
           for(int l = 1;  l <= NDT  ; l++){        
              s0 =sold[l-1];
              s1 =sold[l];
              s2 =sold[l+1];
              if (s[k]>=s0 && s[k]<=s2) {
                 cont = l;
                 break;
              }
           }
           sms0      = (s[k]-s0)/(s1-s0);
           sms0sms1  = ((s[k]-s0)*(s[k]-s1))/((s2-s0)*(s2-s1));
           sms0sms12 = ((s[k]-s0)*(s[k]-s1))/((s1-s0)*(s2-s0)); 
 
           if (cont == (NDT)){
              for (int i = 0; i < DIM; i++) {
                 for (int j = 0; j < DIM; j++) {  	     
 	            B[k][i][j] = Bold[cont-1][i][j]+sms0*(Bold[cont][i][j]-Bold[cont-1][i][j])+
		          sms0sms1*(aux[i][j]-Bold[cont][i][j]) -
		   sms0sms12*(Bold[cont][i][j]-Bold[cont-1][i][j]);
                 }
              }  
           } else {
              for (int i = 0; i < DIM; i++) {
                 for (int j = 0; j < DIM; j++) {  	     
                B[k][i][j] = Bold[cont-1][i][j]+sms0*(Bold[cont][i][j]-Bold[cont-1][i][j])+
		         sms0sms1*(Bold[cont+1][i][j]-Bold[cont][i][j]) -
		  sms0sms12*(Bold[cont][i][j]-Bold[cont-1][i][j]);
                 }
              }       
           } 

        }   
        // condition for B at time t  
        for (int i = 0; i < DIM; i++) {
          for (int j = 0; j < DIM; j++) {
            if (i == j) {
              B[NDT][i][j] = 1.0;
            } else {
              B[NDT][i][j] = 0.0;
            }
          }
        }
          
        /////////      Save the finger tensor  ////////////
        for (int k = 0; k <= NDT; k++) {
           for (int i = 0; i < DIM; i++) {
              for (int j = 0; j < DIM; j++) {
                 dp_set_value(ns->ed.im.dpB[k][i][j], clid, B[k][i][j]);
              }  
           }   
        } 
      }
    // Destroy the iterator
    higcit_destroy(it);
    
    // Sync the finger tensor   
    for (int k = 0; k <= NDT; k++) {
       for (int i = 0; i < DIM; i++) {
          for (int j = 0; j < DIM; j++) {
             dp_sync(ns->ed.im.dpB[k][i][j]);
          }
       }
    } 
   }
    // Loop for each cell
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
      // Get the cell
      hig_cell *c = higcit_getcell(it);
      // Get the cell identifier
      int clid    = mp_lookup(mp, hig_get_cid(c));
      // Get the center of the cell
      Point ccenter;
      hig_get_center(c, ccenter);
      // Get the delta of the cell
      Point cdelta;
      hig_get_delta(c, cdelta);

       // calclute of integral equation (Tau = integrate (-inf, t) M(t,s[k]) H B ds[k])  
    
       // for t>-inf and t<=0   
      if (ns->par.step >= 3) {
       for (int k = 0; k <= Nnow; k++) {
          for (int i = 0; i < DIM; i++) {
             for (int j = 0; j < DIM; j++) {
//                  Bold[k][i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.im.dpB[k][i][j], ns->ed.stn);
                 B[k][i][j]    = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.im.dpB[k][i][j], ns->ed.stn);
             }
          } 
       }        
              
       for (int i = 0; i < DIM; i++) {
         for (int j = 0; j < DIM; j++){
           tensao[i][j] = 0.0;
           I1a          = 0.0;
           I2a          = 0.0;
         }    
       }
       for (int i = 0; i < DIM; i++) {
         for (int j = 0; j < DIM; j++) {
           if (i == j) {
              I1a += B[0][i][j];   
           }    
         }
       }
       I1a += 1.0;
       I2a = B[0][0][0]*B[0][1][1] +  B[0][0][0] +  B[0][1][1] - B[0][0][1]*B[0][0][1]; 
       for (int k = 0; k < M; k++) {
          welambda = lambda[k]*De;  
          for (int i = 0; i < DIM; i++) {
             for (int j = 0; j < DIM; j++) {
                tensao[i][j] += a[k]*(alpha/(alpha-3.0 + beta*I1a + (1.0-beta)*I2a))*exp(-t/welambda)*B[0][i][j];
             }
          }
       }  
       // for t>=0 and t<=tmax   
       for (int k = 1; k <= Nnow-1; k+=2) {
          I1a          = 0.0;
          I2a          = 0.0;
          I1b          = 0.0;
          I2b          = 0.0;
          I1c          = 0.0;
          I2c          = 0.0;
          for (int i = 0; i < DIM; i++) {
             for (int j = 0; j < DIM; j++) {
                if (i == j) {
                    I1a += B[k-1][i][j];   
                    I1b += B[k  ][i][j];
                    I1c += B[k+1][i][j];
                }    
             }    
          }
          I1a += 1.0;   
          I1b += 1.0;
          I1c += 1.0;
          I2a = B[k-1][0][0]*B[k-1][1][1] + B[k-1][0][0] + B[k-1][1][1] - B[k-1][0][1]*B[k-1][0][1]; 
          I2b = B[k  ][0][0]*B[k  ][1][1] + B[k  ][0][0] + B[k  ][1][1] - B[k  ][0][1]*B[k  ][0][1]; 
          I2c = B[k+1][0][0]*B[k+1][1][1] + B[k+1][0][0] + B[k+1][1][1] - B[k+1][0][1]*B[k+1][0][1]; 
          b0 = s[k+1]-s[k-1];
          b1 = (s[k+1]*s[k+1])/2.0 -(s[k-1]*s[k-1])/2.0; 
          b2 = (s[k+1]*s[k+1]*s[k+1])/3.0 -(s[k-1]*s[k-1]*s[k-1])/3.0;
        
          A0 = (b2-(b1*s[k])+(b0*s[k]-b1)*s[k+1])/((s[k]-s[k-1])*s[k+1]-(s[k]*s[k-1])+(s[k-1]*s[k-1]));
          A1 = -((b2-(b1*s[k-1])+(b0*s[k-1]-b1)*s[k+1])/((s[k]-s[k-1])*s[k+1]-(s[k]*s[k])+(s[k-1]*s[k])));
          A2 =  b0 - A0 - A1;
          for (int r = 0; r < M; r++) {
             welambda = lambda[r]*De;  
             for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                   tensao[i][j] +=A0 * (a[r]/welambda)*(alpha/(alpha-3.0+ beta*I1a+((1.0-beta)*I2a)))*exp((s[k-1]-t)/welambda)*B[k-1][i][j]+ 
                                A1 * (a[r]/welambda)*(alpha/(alpha-3.0+ beta*I1b+((1.0-beta)*I2b)))*exp((s[k  ]-t)/welambda)*B[k  ][i][j]+
                                A2 * (a[r]/welambda)*(alpha/(alpha-3.0+ beta*I1c+((1.0-beta)*I2c)))*exp((s[k+1]-t)/welambda)*B[k+1][i][j];
               }
             }
          }
        }

        if ((Nnow%2) != 0) {
           I1a          = 0.0;
           I2a          = 0.0;
           I1b          = 0.0;
           I2b          = 0.0;
           ds    = s[Nnow]-s[Nnow-1];
           for (int i = 0; i < DIM; i++) {
              for (int j = 0; j < DIM; j++) {
                 if (i == j) {
                    I1a += B[Nnow-1][i][j];   
                    I1b += B[Nnow  ][i][j];
                 }    
              }    
           }
           I1a += 1.0;   
           I1b += 1.0;
           I2a = B[Nnow-1][0][0]*B[Nnow-1][1][1] + B[Nnow-1][0][0] + B[Nnow-1][1][1] - B[Nnow-1][0][1]*B[Nnow-1][0][1]; 
           I2b = B[Nnow  ][0][0]*B[Nnow  ][1][1] + B[Nnow  ][0][0] + B[Nnow  ][1][1] - B[Nnow  ][0][1]*B[Nnow  ][0][1];
           c1 = 0.0;
           c2 = 0.0;
           for (int r = 0; r < M; r++) {
              welambda = lambda[r]*De;
              c1 += (a[r]/welambda)*exp((s[Nnow-1]-t)/welambda);
              c2 += (a[r]/welambda)*exp((s[Nnow  ]-t)/welambda);
           }
           for (int i = 0; i < DIM; i++) {
              for (int j = 0; j < DIM; j++) {
                 tensao[i][j] += (ds/2.0)*(c1*((alpha/(alpha-3.0+beta*I1a+(1.0-beta)*I2a))*B[Nnow-1][i][j])+c2*((alpha/(alpha-3.0+beta*I1b+(1.0-beta)*I2b))*B[Nnow][i][j]));
              }
           }
       }
      }
       
       //  Compute S= tau - (1/Re)*(nablaU + nablaUT)  
       real Du[DIM][DIM];
       for (int i = 0; i < DIM; i++) {
          for (int j = 0; j < DIM; j++) {
             // Get Du
             Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
          }
       }
       //  Save S tensor  
       for (int i = 0; i < DIM; i++) {
          for (int j = 0; j < DIM; j++) {
              STensor = tensao[i][j]- ((Du[i][j]+Du[j][i])/(Re));
              // Store Kernel
              dp_set_value(ns->ed.im.dpS[i][j], clid, STensor);
              if (STensor > Smax[i][j]) Smax[i][j] = STensor;
              if (STensor < Smin[i][j]) Smin[i][j] = STensor;
           }  
       }
   }
   // Destroy the iterator
   higcit_destroy(it);
   //   Sync S Tensor 
   for (int i = 0; i < DIM; i++) {
      for (int j = 0; j < DIM; j++) {
         dp_sync(ns->ed.im.dpS[i][j]);
      }
   }
   for (int i = 0; i < DIM; i++) {
      for (int j = 0; j < DIM; j++) {
         // Printing the min and max tensor
         printf("===> %d %d: Tmin = %lf <===> Tmax = %lf <===\n",i,j,Smin[i][j],Smax[i][j]);
      }
   }
   
}
/*
// Calculate the integral equation for KBKZ model
void hig_flow_integral_equation_KBKZ (higflow_solver *ns) {
    real tensao[DIM][DIM], Bold[NDT+1][DIM][DIM], B[NDT+1][DIM][DIM];
    real I1a, I2a, I1b, I2b,I1c, I2c, welambda, alpha, beta;
    real A0, A1, A2, b0, b1, b2, ds, c1, c2;
    real De, Re;
    real a[8], lambda[8];
    real STensor; 
    int  M, Nnow;
    
    real s[NDT+1];
    real sold[NDT+1];
    real t, dt;
    real Q, ds1;

    real Smax[DIM][DIM], Smin[DIM][DIM];
    for (int i = 0; i < DIM; i++) {
       for (int j = 0; j < DIM; j++) {
           Smax[i][j]   = -1.0e16;
           Smin[i][j]   =  1.0e16;
           tensao[i][j] =  0.0;
       }
    }
    
    dt    = ns->par.dt;
    t     = ns->par.t;
    beta  = ns->ed.im.par.beta;
    alpha = ns->ed.im.par.alpha;
    M     = ns->ed.im.par.M;
    De    = ns->ed.im.par.De;
    Re    = ns->par.Re;
    for (int i = 0; i < M; i++) {
       a[i]      = ns->ed.im.par.a[i];
       lambda[i] = ns->ed.im.par.lambda[i];
    }  
    
    //divide the interval [0,t] into (NDT+1) points  (vectors s[k], sold[k])
    Nnow = NDT;
    //real translad;
    if (t > 0.0) { 
     if ( ns->par.step+1 <= Nnow){
       Nnow = ns->par.step+1;
       s[0] = 0.0;
       s[Nnow] = t;
       for (int k = 1; k <= Nnow; k++) {
          s[k]    = s[k-1] + dt;
          sold[k] = s[k];
      //printf(" s[%d] = %f, t= %f, Nnow = %d \n", k, s[k], t, Nnow); 
     }
    } else{
       for (int k = 0; k <= Nnow; k++) {
          sold[k] = ns->ed.im.s[k];
         //printf(" sold[%d] = %f ", k, sold[k]);
       }
       s[0]     = 0.0;
       s[Nnow]  = t;
       Q        = pow((t/dt),1.0/(NDT));  
     //  ds1      = (dt*pow(Q,1));
     //  translad = 0.98*(ds1);
       ds       = 0.0;
       for (int k = 1; k < NDT; k++) {
          ds        = (dt*pow(Q,k));
          s[NDT-k] = t-ds;
       //   s[NDT-k] = t-ds+translad;
          if (ds <= 0.0) {
             printf("ERROR: ds= s[m][%d] - s[m][%d] = %f <= 0!\n",k, k-1, ds);
             exit(1);
          }
          if (s[k] < 0.0) {
             printf("ERROR: s[%d] = %f ", k, s[k]);
             exit(1);
          }
          //printf(" s[%d] = %f ", NDT-k, s[NDT-k]);  
       }  
    }
    for (int k = 0; k <= Nnow; k++) {
      //printf("s[%d] = %f \n",k,s[k]);
     ns->ed.im.s[k] = s[k];
    }
   }
   
    // Get the local sub-domain for the cells
    sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
    // Get the map for the domain properties
    mp_mapper *mp = sd_get_domain_mapper(sdp);
    // Loop for each cell
    higcit_celliterator *it;
    if (ns->par.step+1 >= NDT) { 
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
       // Get the cell
       hig_cell *c = higcit_getcell(it);
       // Get the cell identifier
       int clid    = mp_lookup(mp, hig_get_cid(c));
       // Get the center of the cell
       Point ccenter;
       hig_get_center(c, ccenter);
       // Get the delta of the cell
       Point cdelta;
       hig_get_delta(c, cdelta);
       for (int k = 0; k <= Nnow; k++) {
          for (int i = 0; i < DIM; i++) {
             for (int j = 0; j < DIM; j++) {
                 Bold[k][i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.im.dpB[k][i][j], ns->ed.stn);
                 B[k][i][j]    = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.im.dpB[k][i][j], ns->ed.stn);
             }
          } 
       } 
    
       ////////// interpolate the finger tensor at the new points s[k] ////////////////// 
      int  cont=0;
       real aux[DIM][DIM];
       real s0, s1, s2, sms0, sms0sms1, sms0sms12;
       for (int i = 0; i < DIM; i++) {
          for (int j = 0; j < DIM; j++) {
             if (i == j){
                aux[i][j] = 1.0;
             } else {
                aux[i][j] = 0.0;
             }
          }
       }  
      sold[NDT+1]=t; 
//       if (ns->par.step+1 >= NDT) {
         for (int k = 0; k < NDT; k++) {
             for(int l = 1;  l <= NDT  ; l++){        
                s0 =sold[l-1];
                s1 =sold[l];
                s2 =sold[l+1];
                if (s[k]>=s0 && s[k]<=s2) {
                   cont = l;
                   break;
                }
             }
             sms0      = (s[k]-s0)/(s1-s0);
             sms0sms1  = ((s[k]-s0)*(s[k]-s1))/((s2-s0)*(s2-s1));
             sms0sms12 = ((s[k]-s0)*(s[k]-s1))/((s1-s0)*(s2-s0)); 
 
             if (cont == (NDT)){
                for (int i = 0; i < DIM; i++) {
                   for (int j = 0; j < DIM; j++) {  	     
 	              B[k][i][j] = Bold[cont-1][i][j]+sms0*(Bold[cont][i][j]-Bold[cont-1][i][j])+sms0sms1*(aux[i][j]-Bold[cont][i][j])			    -sms0sms12*(Bold[cont][i][j]-Bold[cont-1][i][j]);
                   }
                }  
             } else {
                for (int i = 0; i < DIM; i++) {
                   for (int j = 0; j < DIM; j++) {  	     
 	              B[k][i][j] = Bold[cont-1][i][j]+sms0*(Bold[cont][i][j]-Bold[cont-1][i][j])+sms0sms1*(Bold[cont+1][i][j]-Bold[cont][i][j])			    -sms0sms12*(Bold[cont][i][j]-Bold[cont-1][i][j]);
                   }
                }       
             } 

          }   
        // condition for B at time t  
        for (int i = 0; i < DIM; i++) {
          for (int j = 0; j < DIM; j++) {
             if (i == j){
                B[NDT][i][j] = 1.0;
             } else {
                B[NDT][i][j] = 0.0;
             }
          }
       }
           
      /////////      Save the finger tensor  ////////////
       for (int k = 0; k <= NDT; k++) {
          for (int i = 0; i < DIM; i++) {
             for (int j = 0; j < DIM; j++) {
                dp_set_value(ns->ed.im.dpB[k][i][j], clid, B[k][i][j]);
             }  
          }   
       } 
    } // end iterator
     // Destroy the iterator
     higcit_destroy(it);
       ///////      Sync the finger tensor    //////////////  
       for (int k = 0; k <= NDT; k++) {
          for (int i = 0; i < DIM; i++) {
             for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.im.dpB[k][i][j]);
             }
          }
       } 
     } // end    if (ns->par.step+1 >= NDT)
     
     
      for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
       // Get the cell
       hig_cell *c = higcit_getcell(it);
       // Get the cell identifier
       int clid    = mp_lookup(mp, hig_get_cid(c));
       // Get the center of the cell
       Point ccenter;
       hig_get_center(c, ccenter);
       // Get the delta of the cell
       Point cdelta;
       hig_get_delta(c, cdelta);
     
       /////  calclute of integral equation ( Tau= integrate (-inf, t) M(t,s[k]) H B ds[k])    //////////////// 
      //////      for t>-inf and t<=0   //////////////////////////////////////////////////////////////////////
       if (ns->par.step+1 >=3){ 
              for (int k = 0; k <= Nnow; k++) {
          for (int i = 0; i < DIM; i++) {
             for (int j = 0; j < DIM; j++) {
//                  Bold[k][i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.im.dpB[k][i][j], ns->ed.stn);
                 B[k][i][j]    = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.im.dpB[k][i][j], ns->ed.stn);
             }
          } 
       }    
           
           
        for (int i = 0; i < DIM; i++) {
          for (int j = 0; j < DIM; j++){
             tensao[i][j] = 0.0;
             I1a          = 0.0;
             I2a          = 0.0;
          }    
       }
       for (int i = 0; i < DIM; i++) {
          for (int j = 0; j < DIM; j++) {
             if (i == j) {
                I1a += B[0][i][j];   
             }    
          }
       }
        I1a += 1.0;
        I2a = B[0][0][0]*B[0][1][1] +  B[0][0][0] +  B[0][1][1] - B[0][0][1]*B[0][0][1]; 
    
       for (int k = 0; k < M; k++) {
          welambda = lambda[k]*De;  
          for (int i = 0; i < DIM; i++) {
             for (int j = 0; j < DIM; j++) {
                tensao[i][j] += a[k]*(alpha/(alpha-3.0 + beta*I1a + (1.0-beta)*I2a))*exp(-t/welambda)*B[0][i][j];
             }
          }
       }  
       //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       /////      for t>=0 and t<=tmax      /////////////////////////////////////////////////////////////////////////////
       for (int k = 1; k <= Nnow-1; k+=2) {
                I1a          = 0.0;
                I2a          = 0.0;
                I1b          = 0.0;
                I2b          = 0.0;
                I1c          = 0.0;
                I2c          = 0.0;
        for (int i = 0; i < DIM; i++) {
           for (int j = 0; j < DIM; j++) {
              if (i == j) {
                  I1a += B[k-1][i][j];   
                  I1b += B[k  ][i][j];
                  I1c += B[k+1][i][j];
              }    
           }    
        }
          I1a += 1.0;   
          I1b += 1.0;
          I1c += 1.0;
          I2a = B[k-1][0][0]*B[k-1][1][1] + B[k-1][0][0] + B[k-1][1][1] - B[k-1][0][1]*B[k-1][0][1]; 
          I2b = B[k  ][0][0]*B[k  ][1][1] + B[k  ][0][0] + B[k  ][1][1] - B[k  ][0][1]*B[k  ][0][1]; 
          I2c = B[k+1][0][0]*B[k+1][1][1] + B[k+1][0][0] + B[k+1][1][1] - B[k+1][0][1]*B[k+1][0][1]; 
      
          b0 = s[k+1]-s[k-1];
          b1 = (s[k+1]*s[k+1])/2.0 -(s[k-1]*s[k-1])/2.0; 
          b2 = (s[k+1]*s[k+1]*s[k+1])/3.0 -(s[k-1]*s[k-1]*s[k-1])/3.0;
          
          A0 = (b2-(b1*s[k])+(b0*s[k]-b1)*s[k+1])/((s[k]-s[k-1])*s[k+1]-(s[k]*s[k-1])+(s[k-1]*s[k-1]));
          A1 = -((b2-(b1*s[k-1])+(b0*s[k-1]-b1)*s[k+1])/((s[k]-s[k-1])*s[k+1]-(s[k]*s[k])+(s[k-1]*s[k])));
          A2 =  b0 - A0 - A1;

          for (int r = 0; r < M; r++) {
             welambda = lambda[r]*De;  
             for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                   tensao[i][j] +=A0 * (a[r]/welambda)*(alpha/(alpha-3.0+ beta*I1a+((1.0-beta)*I2a)))*exp((s[k-1]-t)/welambda)*B[k-1][i][j]+ 
                                  A1 * (a[r]/welambda)*(alpha/(alpha-3.0+ beta*I1b+((1.0-beta)*I2b)))*exp((s[k  ]-t)/welambda)*B[k  ][i][j]+
                                  A2 * (a[r]/welambda)*(alpha/(alpha-3.0+ beta*I1c+((1.0-beta)*I2c)))*exp((s[k+1]-t)/welambda)*B[k+1][i][j];
               }
             }
          } 
       }  
       
        if ((Nnow%2) != 0) {
          I1a          = 0.0;
          I2a          = 0.0;
          I1b          = 0.0;
          I2b          = 0.0;
          int k = Nnow;
          ds    = s[k]-s[k-1];
          for (int i = 0; i < DIM; i++) {
             for (int j = 0; j < DIM; j++) {
                if (i == j) {
                   I1a += B[k-1][i][j];   
                   I1b += B[k  ][i][j];
                }    
             }    
          }
          I1a += 1.0;   
          I1b += 1.0;
          I2a = B[k-1][0][0]*B[k-1][1][1] + B[k-1][0][0] + B[k-1][1][1] - B[k-1][0][1]*B[k-1][0][1]; 
          I2b = B[k  ][0][0]*B[k  ][1][1] + B[k  ][0][0] + B[k  ][1][1] - B[k  ][0][1]*B[k  ][0][1];
          c1 = 0.0;
          c2 = 0.0;
          for (int r = 0; r < M; r++) {
             welambda = lambda[r]*De;
             c1 += (a[r]/welambda)*exp((s[k-1]-t)/welambda);
             c2 += (a[r]/welambda)*exp((s[k  ]-t)/welambda);
          }
          for (int i = 0; i < DIM; i++) {
             for (int j = 0; j < DIM; j++) {
                tensao[i][j] += (ds/2.0)*(c1*((alpha/(alpha-3.0+beta*I1a+(1.0-beta)*I2a))*B[k-1][i][j])+c2*((alpha/(alpha-3.0+beta*I1b+(1.0-beta)*I2b))*B[k][i][j]));
             }
          }
       }
       }
       
       //////////  Compute S= tau - (1/Re)*(nablaU + nablaUT)  ///////////////
       real Du[DIM][DIM];
       for (int i = 0; i < DIM; i++) {
          for (int j = 0; j < DIM; j++) {
             // Get Du
             Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
          }
       }
       /////////  Save S tensor  //////////////////////////////////
       for (int i = 0; i < DIM; i++) {
          for (int j = 0; j < DIM; j++) {
              STensor = tensao[i][j]- ((Du[i][j]+Du[j][i])/(Re));
              // Store Kernel
              dp_set_value(ns->ed.im.dpS[i][j], clid, STensor);
              if (STensor > Smax[i][j]) Smax[i][j] = STensor;
              if (STensor < Smin[i][j]) Smin[i][j] = STensor;
           }  
       }
   }
   // Destroy the iterator
   higcit_destroy(it);
   ///////   Sync S Tensor ///////////////////////////////
   for (int i = 0; i < DIM; i++) {
      for (int j = 0; j < DIM; j++) {
         dp_sync(ns->ed.im.dpS[i][j]);
      }
   }
   for (int i = 0; i < DIM; i++) {
      for (int j = 0; j < DIM; j++) {
         // Printing the min and max tensor
         printf("===> %d %d: Tmin = %lf <===> Tmax = %lf <===\n",i,j,Smin[i][j],Smax[i][j]);
      }
   }
}

*/

// Calculate the integral equation for KBKZ-Fractional model
void hig_flow_integral_equation_KBKZ_Fractional (higflow_solver *ns) {
    real tensao[DIM][DIM], Bold[NDT+1][DIM][DIM], B[NDT+1][DIM][DIM];
    real I1a, I2a, I1b, I2b,I1c, I2c, welambda, alpha, beta;
    real A0, A1, A2, b0, b1, b2, ds, c1, c2;
    real De, Re;
    real alphafr, betafr, V, G; // for Mittag-Leffler function
    real a[8], lambda[8];
    real STensor; 
    int  M, Nnow;
    numc tm, mitt, tm0, mitt0, tm1, mitt1, tm2, mitt2;
    
    real s[NDT+1];
    real sold[NDT+1];
    real t, dt;
    real Q, ds1;

    real Smax[DIM][DIM], Smin[DIM][DIM];
    for (int i = 0; i < DIM; i++) {
       for (int j = 0; j < DIM; j++) {
           Smax[i][j]   = -1.0e16;
           Smin[i][j]   =  1.0e16;
           tensao[i][j] =  0.0;
       }
    }
    
    dt    = ns->par.dt;
    t     = ns->par.t;
    beta  = ns->ed.im.par.beta;
    alpha = ns->ed.im.par.alpha;
    De    = ns->ed.im.par.De;
    Re    = ns->par.Re;
    alphafr = ns->ed.im.par.alpha_frac;
    betafr  = ns->ed.im.par.beta_frac;
    V       = ns->ed.im.par.Phi1;
    G       = ns->ed.im.par.Phi2; 
    
    //divide the interval [0,t] into (NDT+1) points  (vectors s[k], sold[k])
    Nnow = NDT;
    real translad;
    //dt = t/NDT;
    if (t > 0.0) { 
     if ( ns->par.step+1 <= Nnow){
       Nnow = ns->par.step+1;
       //dt = t/Nnow;
       s[0] = 0.0;
       s[Nnow] = t;
       for (int k = 1; k <= Nnow; k++) {
          s[k]    = s[k-1] + dt;
          sold[k] = s[k];
      //printf(" s[%d] = %f, t= %f, Nnow = %d \n", k, s[k], t, Nnow); 
     }
    } else{
       for (int k = 0; k <= Nnow; k++) {
          sold[k] = ns->ed.im.s[k];
         //printf(" sold[%d] = %f ", k, sold[k]);
       }
       s[0]     = 0.0;
       s[Nnow]  = t;
       Q        = pow((t/dt),1.0/(NDT));  
       ds1      = (dt*pow(Q,1));
       translad = 0.98*(ds1);
       ds       = 0.0;
       for (int k = 1; k < NDT; k++) {
          ds        = (dt*pow(Q,k));
          s[NDT-k] = t-ds+translad;
          if (ds <= 0.0) {
             printf("ERROR: ds= s[m][%d] - s[m][%d] = %f <= 0!\n",k, k-1, ds);
             exit(1);
          }
          if (s[k] < 0.0) {
             printf("ERROR: s[%d] = %f ", k, s[k]);
             exit(1);
          }
          //printf(" s[%d] = %f ", NDT-k, s[NDT-k]);  
       }  
    }
    for (int k = 0; k <= Nnow; k++) {
      //printf("s[%d] = %f \n",k,s[k]);
     ns->ed.im.s[k] = s[k];
    }
   }
   
    // Get the local sub-domain for the cells
    sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
    // Get the map for the domain properties
    mp_mapper *mp = sd_get_domain_mapper(sdp);
    // Loop for each cell
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
       // Get the cell
       hig_cell *c = higcit_getcell(it);
       // Get the cell identifier
       int clid    = mp_lookup(mp, hig_get_cid(c));
       // Get the center of the cell
       Point ccenter;
       hig_get_center(c, ccenter);
       // Get the delta of the cell
       Point cdelta;
       hig_get_delta(c, cdelta);
       for (int k = 0; k <= Nnow; k++) {
          for (int i = 0; i < DIM; i++) {
             for (int j = 0; j < DIM; j++) {
                 Bold[k][i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.im.dpB[k][i][j], ns->ed.stn);
                 B[k][i][j]    = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.im.dpB[k][i][j], ns->ed.stn);
                //  Bold[k][i][j] = ns->ed.im.dpB[k][i][j];
             }
          } 
       } 
    
       ////////// interpolate the finger tensor at the new points s[k] ////////////////// 
       int  cont=0;
       real aux[DIM][DIM];
       real s0, s1, s2, sms0, sms0sms1, sms0sms12;
       for (int i = 0; i < DIM; i++) {
          for (int j = 0; j < DIM; j++) {
             if (i == j){
                aux[i][j] = 1.0;
             } else {
                aux[i][j] = 0.0;
             }
          }
       }  
      sold[NDT+1]=t; 
      if (ns->par.step+1 >= NDT) {
         // printf("Entreiiiiiiii ........ \n");
          for (int k = 0; k < NDT; k++) {
             for(int l = 1;  l <= NDT  ; l++){        
                s0 =sold[l-1];
                s1 =sold[l];
                s2 =sold[l+1];
                if (s[k]>=s0 && s[k]<=s2) {
                   cont = l;
                   break;
                }
             }
             sms0      = (s[k]-s0)/(s1-s0);
             sms0sms1  = ((s[k]-s0)*(s[k]-s1))/((s2-s0)*(s2-s1));
             sms0sms12 = ((s[k]-s0)*(s[k]-s1))/((s1-s0)*(s2-s0)); 
             //if (sms0 >=1.0) {
               //  printf("sms0 = %f   sms0sms1 = %f  e  sms0sms12 = %f   k=%d  cont = %d\n", sms0, sms0sms1,sms0sms12,k,cont);
             //}  
             if (cont == (NDT)){
                for (int i = 0; i < DIM; i++) {
                   for (int j = 0; j < DIM; j++) {  	     
 	              B[k][i][j] = Bold[cont-1][i][j]+sms0*(Bold[cont][i][j]-Bold[cont-1][i][j])+sms0sms1*(aux[i][j]-Bold[cont][i][j])			    -sms0sms12*(Bold[cont][i][j]-Bold[cont-1][i][j]);
                   }
                }  
             } else {
                for (int i = 0; i < DIM; i++) {
                   for (int j = 0; j < DIM; j++) {  	     
 	              B[k][i][j] = Bold[cont-1][i][j]+sms0*(Bold[cont][i][j]-Bold[cont-1][i][j])+sms0sms1*(Bold[cont+1][i][j]-Bold[cont][i][j])			    -sms0sms12*(Bold[cont][i][j]-Bold[cont-1][i][j]);
                      //  if ( B[k][i][j] >= 5.0)
                       // printf("B[%d][%d][%d] = %f - apos interpolar \n",i,j,k, B[k][i][j]);                   
                       
                }
                }       
             } 

          }   
        // condition for B at time t  
        for (int i = 0; i < DIM; i++) {
          for (int j = 0; j < DIM; j++) {
             if (i == j){
                B[NDT][i][j] = 1.0;
             } else {
                B[NDT][i][j] = 0.0;
             }
          }
       }
           
      /////////      Save the finger tensor  ////////////
       for (int k = 0; k <= NDT; k++) {
          for (int i = 0; i < DIM; i++) {
             for (int j = 0; j < DIM; j++) {
                dp_set_value(ns->ed.im.dpB[k][i][j], clid, B[k][i][j]);
             }  
          }   
       } 
       ///////      Sync the finger tensor    //////////////  
       for (int k = 0; k <= NDT; k++) {
          for (int i = 0; i < DIM; i++) {
             for (int j = 0; j < DIM; j++) {
                dp_sync(ns->ed.im.dpB[k][i][j]);
             }
          }
       } 
      } 
      // calclute of integral equation ( Tau= integrate (-inf, t) M(t,s[k]) H B ds[k])
      // for t>-inf and t<=0   
       if (ns->par.step+1 >=3){       
        for (int i = 0; i < DIM; i++) {
          for (int j = 0; j < DIM; j++){
             tensao[i][j] = 0.0;
             I1a          = 0.0;
             I2a          = 0.0;
          }    
       }
       for (int i = 0; i < DIM; i++) {
          for (int j = 0; j < DIM; j++) {
             if (i == j) {
                I1a += B[0][i][j];   
             }    
            //printf(" I1a = %f   ",I1a);  
          }
       }
        I1a += 1.0;
        I2a = B[0][0][0]*B[0][1][1] +  B[0][0][0] +  B[0][1][1] - B[0][0][1]*B[0][0][1];
                  //  printf(" \n" );    
       // for (int k = 0; k < M; k++) {
       //   welambda = lambda[k]*De;  
	  tm.real = -G*pow(t,alphafr-betafr)/V;
	  tm.imag = 0;
	  mitt =  mlfv(alphafr-betafr,-betafr, tm, 6);
          for (int i = 0; i < DIM; i++) {
             for (int j = 0; j < DIM; j++) {
		tensao[i][j] += G*pow(t,-1-betafr)*mitt.real*(alpha/(alpha-3.0 + beta*I1a + (1.0-beta)*I2a))*B[0][i][j]; //a[k]*(alpha/(alpha-3.0 + beta*I1a + (1.0-beta)*I2a))*exp(-t/welambda)*B[0][i][j];
                //if (tensao[i][j] > 0.0)
                   //printf("I1a = %f, I2a = %f, Tensao[%d][%d] = %f \n",I1a, I2a, i,j,tensao[i][j]);  
             }
          }
      // }  
      // for t>=0 and t<=tmax  
       for (int k = 1; k <= Nnow-1; k+=2) {
          //for (int i = 0; i < DIM; i++) {
           //  for (int j = 0; j < DIM; j++) {
                //tensao[i][j] = 0.0;
                I1a          = 0.0;
                I2a          = 0.0;
                I1b          = 0.0;
                I2b          = 0.0;
                I1c          = 0.0;
                I2c          = 0.0;
            // }
         // }
          for (int i = 0; i < DIM; i++) {
             for (int j = 0; j < DIM; j++) {
                if (i == j) {
                   I1a += B[k-1][i][j];   
                   I1b += B[k  ][i][j];
                   I1c += B[k+1][i][j];
                }    
             }    
          }
          I1a += 1.0;   
          I1b += 1.0;
          I1c += 1.0;
          I2a = B[k-1][0][0]*B[k-1][1][1] + B[k-1][0][0] + B[k-1][1][1] - B[k-1][0][1]*B[k-1][0][1]; 
          I2b = B[k  ][0][0]*B[k  ][1][1] + B[k  ][0][0] + B[k  ][1][1] - B[k  ][0][1]*B[k  ][0][1]; 
          I2c = B[k+1][0][0]*B[k+1][1][1] + B[k+1][0][0] + B[k+1][1][1] - B[k+1][0][1]*B[k+1][0][1]; 
      
          b0 = s[k+1]-s[k-1];
          b1 = (s[k+1]*s[k+1])/2.0 -(s[k-1]*s[k-1])/2.0; 
          b2 = (s[k+1]*s[k+1]*s[k+1])/3.0 -(s[k-1]*s[k-1]*s[k-1])/3.0;
          
          A0 = (b2-(b1*s[k])+(b0*s[k]-b1)*s[k+1])/((s[k]-s[k-1])*s[k+1]-(s[k]*s[k-1])+(s[k-1]*s[k-1]));
          A1 = -((b2-(b1*s[k-1])+(b0*s[k-1]-b1)*s[k+1])/((s[k]-s[k-1])*s[k+1]-(s[k]*s[k])+(s[k-1]*s[k])));
          A2 =  b0 - A0 - A1;
          
          //printf(" k = %d  A0= %f,   A1 = %f,  A2 = %f\n",k, A0, A1, A2);
          //printf("M = %d ", M);
          //for (int r = 0; r < M; r++) {
             //welambda = lambda[r]*De;  
            // printf("welambda = %f ",welambda );
	  tm0.real = -G*pow((t-s[k-1]),alphafr-betafr)/V;
	  tm0.imag = 0;
	  mitt0 =  mlfv(alphafr-betafr,-betafr, tm0, 6);
	  tm1.real = -G*pow((t-s[k  ]),alphafr-betafr)/V;
	  tm1.imag = 0;
	  mitt1 =  mlfv(alphafr-betafr,-betafr, tm1, 6);
	  tm2.real = -G*pow((t-s[k+1]),alphafr-betafr)/V;
	  tm2.imag = 0;
	  mitt2 =  mlfv(alphafr-betafr,-betafr, tm2, 6);
          
             for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    //if (alpha-3.0+beta*I1a+(1.0-beta)*I2a == 0.0) printf(" Divisão = %f  \n ",(alpha-3.0+beta*I1a+(1.0-beta)*I2a));
                   tensao[i][j] +=A0 * G*pow(t-(s[k-1]),-1-betafr)*mitt0.real*(alpha/(alpha-3.0+ beta*I1a+((1.0-beta)*I2a)))*B[k-1][i][j]+ 
                                  A1 * G*pow((t-s[k  ]),-1-betafr)*mitt1.real*(alpha/(alpha-3.0+ beta*I1b+((1.0-beta)*I2b)))*B[k  ][i][j]+
                                  A2 * G*pow((t-s[k+1]),-1-betafr)*mitt2.real*(alpha/(alpha-3.0+ beta*I1c+((1.0-beta)*I2c)))*B[k+1][i][j];
		   //printf("tm0 = %f, mitt0 = %f, tm1 = %f, mitt1 = %f, tm2 = %f, mitt2 = %f\n", tm0.real, mitt0.real, tm1.real, mitt1.real, tm2.real, mitt2.real);
                   //printf("I1a = %f, I1b = %f, I1c = %f, I2a = %f, I2b = %f, I2c = %f, Tensao[%d][%d] = %f \n",I1a, I1b, I1c, I2a, I2b, I2c, i,j,tensao[i][j]);  
		   //exit(1);
                }
             }
          //} 
       }  
       
        if ((Nnow%2) != 0) {
          I1a          = 0.0;
          I2a          = 0.0;
          I1b          = 0.0;
          I2b          = 0.0;
          int k = Nnow;
          ds    = s[k]-s[k-1];
          for (int i = 0; i < DIM; i++) {
             for (int j = 0; j < DIM; j++) {
                if (i == j) {
                   I1a += B[k-1][i][j];   
                   I1b += B[k  ][i][j];
                }    
             }    
          }
          I1a += 1.0;   
          I1b += 1.0;
          I2a = B[k-1][0][0]*B[k-1][1][1] + B[k-1][0][0] + B[k-1][1][1] - B[k-1][0][1]*B[k-1][0][1]; 
          I2b = B[k  ][0][0]*B[k  ][1][1] + B[k  ][0][0] + B[k  ][1][1] - B[k  ][0][1]*B[k  ][0][1]; 
          //c1 = 0.0;
          //c2 = 0.0;
          //for (int r = 0; r < M; r++) {
            // welambda = lambda[r]*De;
             //c1 += (a[r]/welambda)*exp((s[k-1]-t)/welambda);
             //c2 += (a[r]/welambda)*exp((s[k  ]-t)/welambda);
          //}
          tm0.real = -G*pow((t-s[k-1]),alphafr-betafr)/V;
	  tm0.imag = 0;
	  mitt0 =  mlfv(alphafr-betafr,-betafr, tm0, 6);
	  tm1.real = -G*pow((t-s[k  ]),alphafr-betafr)/V;
	  tm1.imag = 0;
	  mitt1 =  mlfv(alphafr-betafr,-betafr, tm1, 6);
          for (int i = 0; i < DIM; i++) {
             for (int j = 0; j < DIM; j++) {
                tensao[i][j] += (ds/2.0)*(G*pow((t-s[k-1]),-1-betafr)*mitt0.real*((alpha/(alpha-3.0+beta*I1a+(1.0-beta)*I2a))*B[k-1][i][j])+G*pow((t-s[k  ]),-1-betafr)*mitt1.real*((alpha/(alpha-3.0+beta*I1b+(1.0-beta)*I2b))*B[k][i][j]));
                //if ( tensao[i][j] != 0.0)
                   //printf("I1a = %f, I1b = %f, I2a = %f, I2b = %f, Tensao[%d][%d] = %f \n",I1a, I1b, I2a, I2b, i,j,tensao[i][j]);  
             }
          }
       }
       }
       
       //////////  Compute S= tau - (1/Re)*(nablaU + nablaUT)  ///////////////
       real Du[DIM][DIM];
       for (int i = 0; i < DIM; i++) {
          for (int j = 0; j < DIM; j++) {
             // Get Du
             Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
          }
       }
       /////////  Save S tensor  //////////////////////////////////
       for (int i = 0; i < DIM; i++) {
          for (int j = 0; j < DIM; j++) {
              STensor = tensao[i][j]- ((Du[i][j]+Du[j][i])/(Re));
              // Store Kernel
              dp_set_value(ns->ed.im.dpS[i][j], clid, STensor);
              //printf("Tensao[%d][%d] = %lf \n",i,j,tensao[i][j]);
              if (tensao[i][j] > Smax[i][j]) Smax[i][j] = tensao[i][j];
              if (tensao[i][j] < Smin[i][j]) Smin[i][j] = tensao[i][j];
           }  
       }
   }
   // Destroy the iterator
   higcit_destroy(it);
   ///////   Sync S Tensor ///////////////////////////////
   for (int i = 0; i < DIM; i++) {
      for (int j = 0; j < DIM; j++) {
         dp_sync(ns->ed.im.dpS[i][j]);
      }
   }
   for (int i = 0; i < DIM; i++) {
      for (int j = 0; j < DIM; j++) {
         // Printing the min and max tensor
         printf("===> %d %d: Tmin = %lf <===> Tmax = %lf <===\n",i,j,Smin[i][j],Smax[i][j]);
      }
   }
}

// Calculate RHS = Du^t B + B Du
void hig_flow_b_rhs (real B[DIM][DIM], real Du[DIM][DIM], real RHS[DIM][DIM]) {
    real DutB[DIM][DIM], BDu[DIM][DIM];
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            DutB[i][j] = 0.0;
            BDu[i][j] = 0.0;
            for (int k = 0; k < DIM; k++) {
                DutB[i][j] += Du[k][i]*B[k][j];
                BDu[i][j]  += B[i][k]*Du[k][j];
            }
        }
    }
    // Calculate RHS = Du^t B + B Du
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            RHS[i][j] = DutB[i][j] + BDu[i][j];
        }
    }
}


// Get the velocity at cell center 
void hig_flow_velocity_at_center_cell_integral (higflow_solver *ns, Point ccenter, Point cdelta, real u[DIM]) {
    for (int dim = 0; dim < DIM; dim++) {
        // Verity if is in facet
        int infacet;
        // Get the velocity in the left facet center
        real ul = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
        // Get the velocity in the right facet center
        real ur = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
        // Setting the velocity at cell center
        u[dim]  = 0.5*(ul + ur);
    }
}

// Get the derivative of B tensor 
void hig_flow_derivative_b_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, int k, int i, int j, real Bcenter, real dBdx[DIM]) {
    for (int dim = 0; dim < DIM; dim++) {
        int incell_left, incell_right;
        // Get the Kernel in the left cell
        real Bleft = compute_center_p_left_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.im.dpB[k][i][j], ns->ed.stn, &incell_left);
        // Get the Kernel in the right cell
        real Bright = compute_center_p_right_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.im.dpB[k][i][j], ns->ed.stn, &incell_right);
        // Compute the Kernel derivative
        if ((incell_left == 1) && (incell_right == 1)) { 
           dBdx[dim] = compute_dpdx_at_point(cdelta, dim, 1.0, Bleft, Bright);
        } else if (incell_right == 1) { 
           dBdx[dim] = compute_dpdxr_at_point(cdelta, dim, 1.0, Bcenter, Bright);
        } else {
           dBdx[dim] = compute_dpdxl_at_point(cdelta, dim, 1.0, Bleft, Bcenter);
        }
    }
}


// *******************************************************************
// Calculate convective tensor term CUBISTA
// *******************************************************************
real hig_flow_convective_tensor_term_b_cubista(higflow_solver *ns, distributed_property *dpu, sim_domain *sdp, sim_stencil *stn, real B[DIM][DIM], Point ccenter, Point cdelta, int dim, int k, int i, int j) {
    real  vbar[DIM], dBdx[dim], kr, krr, kl, kll, kc, a, b, c, d, e, tol, fi, conv1,conv2;
    a     = 1.7500;
    b     = 0.3750;
    c     = 0.7500;
    d     = 0.1250;
    e     = 0.2500;
    tol   = 1.0e-14;
    conv1 = 0.0;
    conv2 = 0.0;
    int   incell_r, incell_l, incell_ll, incell_rr, infacet;
    // Get the kernel at center cell
    kc  = B[i][j];
    // Get the low, high, lowlow, highhigh component kernel at center cell
    kl  = compute_center_p_left_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.im.dpB[k][i][j], ns->ed.stn, &incell_l); 
    kr  = compute_center_p_right_22(ns->ed.sdED, ccenter, cdelta, dim, 1.0, ns->ed.im.dpB[k][i][j], ns->ed.stn, &incell_r); 
    kll = compute_center_p_left_22(ns->ed.sdED, ccenter, cdelta, dim, 2.0, ns->ed.im.dpB[k][i][j], ns->ed.stn, &incell_ll);
    krr = compute_center_p_right_22(ns->ed.sdED, ccenter, cdelta, dim, 2.0, ns->ed.im.dpB[k][i][j], ns->ed.stn, &incell_rr);
    // Get the velocity  v1bar(i+1/2,j) in the facet center
    vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
    if (vbar[dim] > 0.0){
        if (fabs(kr - kl) <= tol){
            conv1 = vbar[dim]*kc;
        }else {
            fi = (kc - kl)/(kr - kl);
            if ((fi <= 0.0) || (fi >= 1.0)) {
                conv1 = vbar[dim]*kc;
            }else {
                if (fi < b){ 
                    if (incell_l == 1)                    conv1 = vbar[dim]*(a*kc - c*kl);
                    else                                  conv1 = vbar[dim]*kc;
                }
	        if ((fi >= b) && (fi <= c)){
                    if ((incell_l == 1)&&(incell_r == 1)) conv1 = vbar[dim]*(c*kc + b*kr -d*kl);
                    else                                  conv1 = vbar[dim]*kc;
                }
	        if (fi > c){ 
                    if (incell_r == 1)                    conv1 = vbar[dim]*(e*kc + c*kr);
                    else                                  conv1 = vbar[dim]*kc;
                }
                    
            }    
        }
    //v1bar < 0.0
    }else {
        if ((incell_r == 1) && (incell_rr == 1)){
            if (fabs(kc - krr) <= tol){
                conv1 = vbar[dim]*kr;
            }else {
                fi = (kr- krr)/(kc - krr);
                if ((fi <= 0.0) || (fi >= 1.0)) {
                    conv1 = vbar[dim]*kr;
                }else {
		    if (fi < b) 
                        conv1 = vbar[dim]*(a*kr - c*krr);
                    if ((fi >= b) && (fi <= c))
                        conv1 = vbar[dim]*(c*kr + b*kc -d*krr);
	            if (fi > c) 
                        conv1 = vbar[dim]*(c*kc + e*kr);
                }
            }
        //Return upwind value at boundary
        }else if ((incell_r == 1) && (incell_rr == 0)){
            if (fabs(kc - krr) <= tol){
                conv1 = vbar[dim]*kr;
            }else {
                fi = (kr- krr)/(kc - krr);
                if ((fi <= 0.0) || (fi >= 1.0)) {
                    conv1 = vbar[dim]*kr;
                }else {
		    if (fi <= c) 
                        conv1 = vbar[dim]*kr;
	            if (fi > c) 
                        conv1 = vbar[dim]*(c*kc + e*kr);
                }
            }/*
            vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            if (vbar[dim] > 0.0) conv1 = vbar[dim]*kc;
            else                 conv1 = vbar[dim]*kr;
            vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            if (vbar[dim] > 0.0) conv2 = vbar[dim]*kl;
            else                 conv2 = vbar[dim]*kc;
            return ((conv1 - conv2)/cdelta[dim]); */
        }else {
                vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
                if (vbar[dim] > 0.0) conv1 = vbar[dim]*kc;
                else                 conv1 = vbar[dim]*kc;
                vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
                if (vbar[dim] > 0.0) conv2 = vbar[dim]*kl;
                else                 conv2 = vbar[dim]*kc;
                return ((conv1 - conv2)/cdelta[dim]); 
        }
        
    }
    // Get the velocity  v2bar(i-1/2,j) in the facet center
    vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
    if (vbar[dim] > 0.0){
        if ((incell_l == 1) && (incell_ll == 1)){
            if (fabs(kc-kll) <= tol) {
	        conv2 = vbar[dim]*kl;
            }else {
	        fi = (kl - kll)/(kc - kll);
	        if ((fi <= 0.0) || (fi >= 1.0)) {
	            conv2 = vbar[dim]*kl;
	        }else {
	            if (fi < b)
	                conv2 = vbar[dim]*(a*kl - c*kll);
	            if ((fi >= b) && (fi <= c))
	                conv2 = vbar[dim]*(b*kc + c*kl - d*kll);
	            if (fi > c)  
	                conv2 = vbar[dim]*(c*kc + e*kl);
	        }
	    }
        }else if ((incell_l == 1) && (incell_ll == 0)){
            if (fabs(kc-kll) <= tol) {
	        conv2 = vbar[dim]*kl;
            }else {
	        fi = (kl - kll)/(kc - kll);
	        if ((fi <= 0.0) || (fi >= 1.0)) {
	            conv2 = vbar[dim]*kl;
	        }else {
	            if (fi <= c)
	                conv2 = vbar[dim]*kl;
	            if (fi > c)  
	                conv2 = vbar[dim]*(c*kc + e*kl);
	        }
	    }/*
            vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            if (vbar[dim] > 0.0) conv1 = vbar[dim]*kc;
            else                 conv1 = vbar[dim]*kr;
            vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
            if (vbar[dim] > 0.0) conv2 = vbar[dim]*kl;
            else                 conv2 = vbar[dim]*kc;
            return ((conv1 - conv2)/cdelta[dim]); */
       }else {
                vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
                if (vbar[dim] > 0.0) conv1 = vbar[dim]*kc;
                else                 conv1 = vbar[dim]*kr;
                vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
                if (vbar[dim] > 0.0) conv2 = vbar[dim]*kc;
                else                 conv2 = vbar[dim]*kc;
                return ((conv1 - conv2)/cdelta[dim]); 
        } 
    }else {
    //v2bar < 0.0 
        if (fabs(kl - kr) <= tol) {
            conv2 = vbar[dim]*kc;
        }else {
            fi = (kc - kr)/(kl - kr);
            if ((fi <= 0.0) || (fi >= 1.0)) {
                conv2 = vbar[dim]*kc;
            }else {
	        if (fi < b){
                    if (incell_r == 1)                    conv2 = vbar[dim]*(a*kc - c*kr);
                    else                                  conv2 = vbar[dim]*kc;
                }
	        if ((fi >= b) && (fi <= c)){
                    if ((incell_l == 1)&&(incell_r == 1)) conv2 = vbar[dim]*(c*kc + b*kl -d*kr);
                    else                                  conv2 = vbar[dim]*kc;
                }
	        if (fi > c){ 
                    if (incell_l == 1)                    conv2 = vbar[dim]*(e*kc + c*kl);
                    else                                  conv2 = vbar[dim]*kc;
                }
	    }
        }
    }
    return ((conv1-conv2)/cdelta[dim]);
}


// *******************************************************************
// Navier-Stokes Step for the Explicit Euler Method
// *******************************************************************
void higflow_explicit_euler_intermediate_velocity_viscoelastic_integral(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]) {
    // Get the facet iterator
    higfit_facetiterator *fit;
    // Get the local domain for cell
    sim_domain *sdp    = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++) {
        // Get the local domain for facet cell
        for (int dim2 = 0; dim2 < DIM; dim2++) {
            sfdu[dim2] = psfd_get_local_domain(ns->psfdu[dim2]);
        }
        // Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Set the computational cell
            higflow_computational_cell_viscoelastic_integral(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Difusive term contribution
            rhs += higflow_difusive_term(ns, fdelta);
            // Compute the intermediate velocity
            real ustar = ns->cc.ucell + ns->par.dt * rhs;
            // Update the distributed property intermediate velocity
            dp_set_value(dpustar[dim], flid, ustar);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Syncing the intermediate velocity
        dp_sync(dpustar[dim]);
    }
}


// *******************************************************************
// Navier-Stokes Step for the Second Order Explicit Runge-Kutta Method
// *******************************************************************
void higflow_explicit_runge_kutta_2_intermediate_velocity_viscoelastic_integral(higflow_solver *ns) {
    // Calculate the auxiliar velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_viscoelastic_integral(ns, ns->dpu, ns->dpuaux);
    // Calculate the star velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_viscoelastic_integral(ns, ns->dpuaux, ns->dpustar);
    // Calculate the order 2 Runge-Kutta method using the euler method
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Loop for each dimension
    higfit_facetiterator *fit;
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local partitioned domain for facets
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Get the intermediate velocity
            real u     = dp_get_value(ns->dpu[dim], flid);
            real ustar = dp_get_value(ns->dpustar[dim], flid);
            // Compute the intermediate velocity
            real urk2  = 0.5*(u + ustar);
            // Set the final velocity in the distributed velocity property
            dp_set_value(ns->dpustar[dim], flid, urk2);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the ditributed velocity property
        dp_sync(ns->dpustar[dim]);
    }
}

// *******************************************************************
// Navier-Stokes Step for third Order Explicit Runge-Kutta Method
// *******************************************************************
void higflow_explicit_runge_kutta_3_intermediate_velocity_viscoelastic_integral(higflow_solver *ns) {
    // Calculate the auxiliar velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_viscoelastic_integral(ns, ns->dpu, ns->dpuaux);
    // Calculate the second stage velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_viscoelastic_integral(ns, ns->dpuaux, ns->dpustar);
    // Calculate the order 2 Runge-Kutta method using the euler method
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Loop for each dimension
    higfit_facetiterator *fit;
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local partitioned domain for facets
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Get the intermediate velocity
            real u     = dp_get_value(ns->dpu[dim], flid);
            real ustar = dp_get_value(ns->dpustar[dim], flid);
            // Compute the intermediate velocity
            real urk3  = 0.75*u + 0.25*ustar;
            // Set the final velocity in the distributed velocity property
            dp_set_value(ns->dpuaux[dim], flid, urk3);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the ditributed velocity property
        dp_sync(ns->dpuaux[dim]);
    }
    // Calculate the order 2 Runge-Kutta method using the euler method
    higflow_explicit_euler_intermediate_velocity_viscoelastic_integral(ns, ns->dpuaux, ns->dpustar);
    // Loop for each dimension
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local partitioned domain for facets
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
    // Calculate the third stage velocity by the explicit euler method
        for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Get the intermediate velocity
            real u     = dp_get_value(ns->dpu[dim], flid);
            real ustar = dp_get_value(ns->dpustar[dim], flid);
            // Compute the intermediate velocity
            real urk3  = u/3.0 + 2.0*ustar/3.0;
            // Set the final velocity in the distributed velocity property
            dp_set_value(ns->dpustar[dim], flid, urk3);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the ditributed velocity property
        dp_sync(ns->dpustar[dim]);
    }
}

// *******************************************************************
// Navier-Stokes Step for the Implicit Euler Method
// *******************************************************************
void higflow_semi_implicit_euler_intermediate_velocity_viscoelastic_integral(higflow_solver *ns) {
    // Get the facet iterator
    higfit_facetiterator *fit;
    // Get the local domain for cell
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Get the local domain for facet cell
    for (int dim2 = 0; dim2 < DIM; dim2++) {
        sfdu[dim2] = psfd_get_local_domain(ns->psfdu[dim2]);
    }
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++) {
        // Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Set the computational cell
            higflow_computational_cell_viscoelastic_integral(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Total contribuition terms by delta t
            rhs *= ns->par.dt;
            // Velocity term contribution
            rhs += ns->cc.ucell;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w = - ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
                alpha -= 2.0 * w ;
                Point p;
                POINT_ASSIGN(p, fcenter);
                // Stencil point update: right point
                p[dim2] = fcenter[dim2] + fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
                // Stencil point update: left point
                p[dim2] = fcenter[dim2] - fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
            }
            alpha = 1.0 + alpha;
            // Get the stencil
            sfd_get_stencil(sfdu[dim], fcenter, fcenter, alpha, ns->stn);
            // Get the index of the stencil
            int *ids   = psfd_stn_get_gids(ns->psfdu[dim],ns->stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->stn);
            // Get the cell identifier of the cell
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Set the right side of solver linear system
            slv_set_bi(ns->slvu[dim], fgid, stn_get_rhs(ns->stn));
            // Set the line of matrix of the solver linear system
            slv_set_Ai(ns->slvu[dim], fgid, numelems, ids, vals);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Assemble the solver
        slv_assemble(ns->slvu[dim]);
        // Solve the linear system
        slv_solve(ns->slvu[dim]);
        // Get the solution of linear system

        // Gets the values of the solution
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Get the value of ustar
            real ustar = slv_get_xi(ns->slvu[dim], fgid);
            // Set the value of ustar
            dp_set_value(ns->dpustar[dim], flid, ustar);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Syncing the intermediate velocity
        dp_sync(ns->dpustar[dim]);
    }
}


// *******************************************************************
// Navier-Stokes Step for the Semi-Implicit Crank-Nicolson Method
// *******************************************************************
void higflow_semi_implicit_crank_nicolson_intermediate_velocity_viscoelastic_integral(higflow_solver *ns) {
    // Get the facet iterator
    higfit_facetiterator *fit;
    // Get the local domain for cell
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Get the local domain for facet cell
    for (int dim2 = 0; dim2 < DIM; dim2++) {
        sfdu[dim2] = psfd_get_local_domain(ns->psfdu[dim2]);
    }
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++) {
        // Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Set the computational cell
            higflow_computational_cell_viscoelastic_integral(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Diffusive term term contribution
            rhs += 0.5 * higflow_difusive_term(ns, fdelta);
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Total contribuition terms times delta t
            rhs *= ns->par.dt;
            // Velocity term contribution
            rhs += ns->cc.ucell;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w = - 0.5 * ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
                alpha -= 2.0 * w ;
                Point p;
                POINT_ASSIGN(p, fcenter);
                // Stencil point update: right point
                p[dim2] = fcenter[dim2] + fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
                // Stencil point update: left point
                p[dim2] = fcenter[dim2] - fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
            }
            alpha = 1.0 + alpha;
            // Get the stencil
            sfd_get_stencil(sfdu[dim], fcenter, fcenter, alpha, ns->stn);
            // Get the index of the stencil
            int *ids   = psfd_stn_get_gids(ns->psfdu[dim],ns->stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->stn);
            // Get the cell identifier of the cell
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Set the right side of solver linear system
            slv_set_bi(ns->slvu[dim], fgid, stn_get_rhs(ns->stn));
            // Set the line of matrix of the solver linear system
            slv_set_Ai(ns->slvu[dim], fgid, numelems, ids, vals);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Assemble the solver
        slv_assemble(ns->slvu[dim]);
        // Solve the linear system
        slv_solve(ns->slvu[dim]);
        // Gets the values of the solution
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Get the value of ustar
            real ustar = slv_get_xi(ns->slvu[dim], fgid);
            // Set the value of ustar
            dp_set_value(ns->dpustar[dim], flid, ustar);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Syncing the intermediate velocity
        dp_sync(ns->dpustar[dim]);
    }
}

// *******************************************************************
// Navier-Stokes Step for the Implicit BDF2 Method
// *******************************************************************
void higflow_semi_implicit_bdf2_intermediate_velocity_viscoelastic_integral(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]) {
    // Firt stage of Tr-BDF2 method
    // Get the facet iterator
    higfit_facetiterator *fit;
    // Get the local domain for cell
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Get the local domain for facet cell
    for (int dim2 = 0; dim2 < DIM; dim2++) {
        sfdu[dim2] = psfd_get_local_domain(ns->psfdu[dim2]);
    }
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++) {
        // Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Set the computational cell
            higflow_computational_cell_viscoelastic_integral(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Total contribuition terms times delta t
            rhs *= 0.25*ns->par.dt;
            // Velocity term contribution
            rhs += ns->cc.ucell;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w = - 0.25*ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
                alpha -= 2.0 * w ; //divide po 4 para usar regra trapezio em t(n+1/2)
                Point p;
                POINT_ASSIGN(p, fcenter);
                // Stencil point update: right point
                p[dim2] = fcenter[dim2] + fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
                // Stencil point update: left point
                p[dim2] = fcenter[dim2] - fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
            }
            alpha = 1.0 + alpha;
            // Get the stencil
            sfd_get_stencil(sfdu[dim], fcenter, fcenter, alpha, ns->stn);
            // Get the index of the stencil
            int *ids   = psfd_stn_get_gids(ns->psfdu[dim],ns->stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->stn);
            // Get the cell identifier of the cell
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Set the right side of solver linear system
            slv_set_bi(ns->slvu[dim], fgid, stn_get_rhs(ns->stn));
            // Set the line of matrix of the solver linear system
            slv_set_Ai(ns->slvu[dim], fgid, numelems, ids, vals);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Assemble the solver
        slv_assemble(ns->slvu[dim]);
        // Solve the linear system
        slv_solve(ns->slvu[dim]);
        // Gets the values of the solution
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Get the value of ustar
            real uaux = slv_get_xi(ns->slvu[dim], fgid);
            // Set the value of ustar
            dp_set_value(ns->dpuaux[dim], flid, uaux);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Syncing the intermediate velocity
        dp_sync(ns->dpuaux[dim]);
    }
    //Second Stage of Tr-BDF2
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++) {
        // Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Set the computational cell
            higflow_computational_cell_viscoelastic_integral(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            //Get the uaux
            real uaux = dp_get_value(ns->dpuaux[dim], flid);
            // Right hand side equation
            real rhs = 0.0;
            rhs = (4.0*uaux - ns->cc.ucell)/3.0;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w = - 1.0/3.0*ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
                alpha -=  2.0 * w ;
                Point p;
                POINT_ASSIGN(p, fcenter);
                // Stencil point update: right point
                p[dim2] = fcenter[dim2] + fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
                // Stencil point update: left point
                p[dim2] = fcenter[dim2] - fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
            }
            alpha = 1.0 + alpha;
            // Get the stencil
            sfd_get_stencil(sfdu[dim], fcenter, fcenter, alpha, ns->stn);
            // Get the index of the stencil
            int *ids   = psfd_stn_get_gids(ns->psfdu[dim],ns->stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->stn);
            // Get the cell identifier of the cell
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Set the right side of solver linear system
            slv_set_bi(ns->slvu[dim], fgid, stn_get_rhs(ns->stn));
            // Set the line of matrix of the solver linear system
            slv_set_Ai(ns->slvu[dim], fgid, numelems, ids, vals);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Assemble the solver
        slv_assemble(ns->slvu[dim]);
        // Solve the linear system
        slv_solve(ns->slvu[dim]);
        // Get the solution of linear system
        //Vec *vecu = slv_get_solution_vec(ns->slvu[dim]);
        // Gets the values of the solution
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Get the value of ustar
            real ustar = slv_get_xi(ns->slvu[dim], fgid);
            // Set the value of ustar
            dp_set_value(ns->dpustar[dim], flid, ustar);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Syncing the intermediate velocity
        dp_sync(ns->dpustar[dim]);
    }
}

// Print the Polymeric Tensor
void higflow_print_polymeric_tensor_integral(higflow_solver *ns) {
    if (ns->contr.flowtype == 4) {
        FILE *data1;
        //FILE *data2;
        //FILE *data3;
        //FILE *data4;
        //FILE *data5;
        //FILE *data6;
        data1 = fopen("Taus.dat", "w");
        //data2 = fopen("Sl.dat", "w");
        //data3 = fopen("Tr.dat", "w");
        //data4 = fopen("Sr.dat", "w");
        //data5 = fopen("T8.dat", "w");
        //data6 = fopen("S8.dat", "w");
        real Du[DIM][DIM], S[DIM][DIM], D[DIM][DIM];
        real Re   = ns->par.Re;
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdp);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid    = mp_lookup(mp, hig_get_cid(c));
            // Get the inside/outside inflow point cell
            int inflowcell, outflowcell;
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Get the velocity derivative tensor Du and the Kernel tensor
                        for (int i = 0; i < DIM; i++) {
                for (int j = 0; j < DIM; j++) {
                    // Get Du
                    Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.im.dpD[i][j], ns->ed.stn);
                    D[i][j]  = (Du[i][j]+Du[j][i]);
                    // Get S tensor
                    S[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.im.dpS[i][j], ns->ed.stn);
                }
            }
            if (ccenter[0]<4.99 && ccenter[0]>4.94){       //if (ccenter[0]<0.06e-0 && ccenter[0]>0.04e-0){
             fprintf(data1, "%lf  %lf  %15.12lf  %15.12lf  %15.12lf  %15.12lf\n", ccenter[0], ccenter[1], S[0][0]+ D[0][0]/Re, S[0][1]+ D[0][1]/Re, S[1][0] + D[1][0]/Re, S[1][1]+ D[1][1]/Re); 
            }
            //if (ccenter[0]<0.504 && ccenter[0]>0.501){
            //if (ccenter[0]>3.94e-0 && ccenter[0]<3.96e-0){
                   //fprintf(data3, "%lf  %lf  %15.12lf  %15.12lf  %15.12lf  %15.12lf\n", ccenter[0], ccenter[1], S[0][0]+ 2.0*(1-beta)*D[0][0]/Re, S[0][1]+ 2.0*(1-beta)*D[0][1]/Re, S[1][0]+ 2.0*(1-beta)*D[1][0]/Re, S[1][1]+ 2.0*(1-beta)*D[1][1]/Re); 
                  // fprintf(data4, "%15.08lf  %15.08lf  %15.12lf  %15.12lf  %15.12lf  %15.12lf\n", ccenter[0], ccenter[1], S[0][0]/100.0, S[0][1], S[1][0], S[1][1]); 
            //}
        }
           // Destroy the iterator
        higcit_destroy(it);
        fclose(data1);
        //fclose(data2);
        //fclose(data3);
        //fclose(data4);
        //fclose(data5);
        //fclose(data6);
    }
}



// One step of the Navier-Stokes the projection method
void higflow_solver_step_viscoelastic_integral(higflow_solver *ns) {
    // Boundary condition for velocity
    higflow_boundary_condition_for_velocity(ns);
    // Calculate the source term
    higflow_calculate_source_term(ns);
    // Calculate the facet source term
    higflow_calculate_facet_source_term(ns);
    // Calculate the intermediated velocity
    switch (ns->contr.tempdiscrtype) {
        case 0:
           // Explicit Euler method
           higflow_explicit_euler_intermediate_velocity_viscoelastic_integral(ns, ns->dpu, ns->dpustar);
           break;
        case 1: 
           // Explicit RK2 method
           higflow_explicit_runge_kutta_2_intermediate_velocity_viscoelastic_integral(ns);
           break;
        case 2: 
           // Explicit RK3 method
           higflow_explicit_runge_kutta_3_intermediate_velocity_viscoelastic_integral(ns);
           break;
        case 3: 
           // Semi-Implicit Euler Method
           higflow_semi_implicit_euler_intermediate_velocity_viscoelastic_integral(ns);
           break;
        case 4: 
           // Semi-Implicit Crank-Nicolson Method
           higflow_semi_implicit_crank_nicolson_intermediate_velocity_viscoelastic_integral(ns);
           break;
        case 5: 
           // Semi-Implicit Crank-Nicolson Method
           higflow_semi_implicit_bdf2_intermediate_velocity_viscoelastic_integral(ns, ns->dpu, ns->dpustar);
           break;
    }
    // Set outflow for ustar velocity 
    //higflow_outflow_ustar_step(ns);
    // Boundary condition for pressure
    higflow_boundary_condition_for_pressure(ns);
    // Calculate the pressure
    higflow_pressure(ns);
    // Calculate the final velocity
    higflow_final_velocity(ns);
    // Boundary condition for velocity
    higflow_boundary_condition_for_velocity(ns);
    // Set outflow for velocity
    //higflow_outflow_u_step(ns);
    // Calculate the final pressure
    higflow_final_pressure(ns);
    // Calculate the velocity derivative tensor
    higflow_compute_velocity_derivative_tensor(ns);
    // Constitutive Equation Step for the Explicit Euler Method
    higflow_explicit_euler_constitutive_equation_integral(ns);
    // Computing the Polymeric and Extra-Stress Tensors
    hig_flow_integral_equation(ns);
    // Print the Polymeric Tensor
    //higflow_print_polymeric_tensor_integral(ns);
}

