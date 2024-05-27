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