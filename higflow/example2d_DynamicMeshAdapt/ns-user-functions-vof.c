#include "ns-example-2d.h"

/******************************************************************************************/
/******************************************************************************************/
/*************************** multiphase user functions ************************************/
/******************************************************************************************/
/******************************************************************************************/

// Geometria 2-D comum aos exemplos VOF: as seis funcoes de recorte de poligono
// eram identicas em cinco exemplos e passaram a viver num lugar so'.
#include "../examples-common/vof-geometry-2d.c"

// Volume fraction


real func (Point p) {
  real value;
  Point c, pc, r, cav;

  ////////////// Zalesak's disk ////////////
  c[0]   = 0.5;  c[1]   = 0.75;
  r[0]   = 0.15; r[1]   = 0.15;
  cav[0] = 0.05; cav[1] = 0.2479;
  pc[0]  = p[0] - c[0]; pc[1] = p[1] - c[1];
  real outcirc = 1.0 - (pc[0]*pc[0]/(r[0]*r[0]) + pc[1]*pc[1]/(r[1]*r[1]));
  real incav1  = (fabs(pc[0]) - 0.5*cav[0]) / r[0];
  real incav2  = (pc[1] - (cav[1] - r[1]))  / r[1];
  real incav   = fmax(incav1, incav2);
  value = fmin(outcirc, incav);

  /////// Shearing Droplet (disabled) ////////////
  c[0] = p_par->center[0].val; c[1] = p_par->center[1].val;
  r[0] = 0.25; r[1] = 0.25;
  pc[0] = p[0] - c[0]; pc[1] = p[1] - c[1];
  value = 1.0 - (pc[0]*pc[0]/(r[0]*r[0]) + pc[1]*pc[1]/(r[1]*r[1]));

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


