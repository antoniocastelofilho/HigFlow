//*******************************************************
// Estruture for Mittag-Leffler function
//*******************************************************

#include <math.h>
#include <stdio.h>

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


numc mlfv(double alpha, double beta, numc z, int fi) {
   int k,m,h;
   double two = 2, um = 1, pcin = 0.5, ten = 10, cinc=5, six = 6;
   double  zero = 0.0, dum=2.1, ccinc=5.5, vint=20, cinq=50;
   double r0, rc, angz, eps, aaz, norz, pr1, pr2, pq1, pq2;
   int l1,l2,l3,l4;
   numc   newsum, term, aux, aux1, a1, a2, oldsum, zn, ic, res, raiz1z, raiz2z;
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
      printf("beta must be smaller than five\n");//cout << "beta must be smaller than zero" << endl;
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
      // complex<double> oldsum(0.0, 0.0);
      // oldsum=dcmplx(0d0,0d0);
      oldsum.real=zero;
      oldsum.imag=zero;    
      // printf("oldsum = %f+ i %f\n", oldsum.real, oldsum.imag); //cout << "oldsum=" << oldsum << endl;
      k=0;
      while ((alpha*k+beta) <= zero) {
         k=k+1;
      }
      // newsum=pow(real(z),k/tgamma(alpha*k+beta));
      // printf("k= %d\n", k);//cout << "k=" << k << endl;
      // printf("z= %f+ i %f\n",z.real, z.imag);//cout << "z=" << z << endl;
      // printf("pow(z,k)= %f + i %f\n", pow(norz,k)*cos(angz*k), pow(norz,k)*sin(angz*k));//cout << "pow(z,k)=" << pow(z,k) << endl;
      newsum.real=pow(norz,k)*cos(angz*k)/tgamma(alpha*k+beta); //pow(z,k)/tgamma(alpha*k+beta);
      newsum.imag=pow(norz,k)*sin(angz*k)/tgamma(alpha*k+beta);//pow(z,k)/tgamma(alpha*k+beta);

      // double summation because z can be negative 
      while (newsum.real != oldsum.real || newsum.imag != oldsum.imag) {
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

   //! the matlab function fix rounds toward zero, can  use floor since alpha is positive
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
            //cout << "conta=" << (2.0*pi*h)/(m+1) << endl;    
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


