// The Mittag-Leffler function, needed by the fractional viscoelastic models where it
// plays the role the exponential plays for an ordinary relaxation.
//
// Self-contained numerics over a local complex type (`numc`) -- it does not use the
// solver or the mesh, and can be tested on its own.  NOTE THE `pi` MACRO: this
// header defines the bare identifier `pi`, so any translation unit that includes it
// cannot use `pi` as a name of its own.

//*******************************************************
// Estruture for Mittag-Leffler function
//*******************************************************

#ifndef HIG_FLOW_MITTAG_LEFFLER
#define HIG_FLOW_MITTAG_LEFFLER

#define pi 3.1415926535897932384626434
typedef struct{
  double real;
  double imag;
} numc;

numc kk(double r, double alfa, double beta1, numc z);

numc pp(double r, double alphaf, double betaf, numc z, double epsn) ;

/* Romberg Integration*/

numc rombint(char funfcn, double a, double b, int order, double v1, double v2, numc v3, double v4);

numc mlfv(double alpha, double beta, numc z, int fi);

#endif