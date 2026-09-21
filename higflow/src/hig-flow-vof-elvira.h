// ELVIRA normal estimation: build candidate interface lines from the column and row
// sums of the 3x3 neighbourhood, then keep the one that best reproduces the
// neighbours' fractions.
//
// It is a SELECTION among candidates, not a formula, which makes it robust where a
// finite-difference normal is noisy and more expensive than one.
//
// PART OF THE VOF FAMILY: a volume fraction per cell says how much of each phase is
// there, and the interface is RECONSTRUCTED from it rather than tracked.  The steps
// are: estimate a normal, place an interface from the normal and the fraction, then
// advect the fraction.  These files are the alternative methods for those steps, not
// a pipeline -- several are alternatives to each other.

#ifndef HIG_FLOW_VOF_ELVIRA
#define HIG_FLOW_VOF_ELVIRA

#include "hig-flow-discret.h"
#include "hig-flow-vof-plic.h"
#include "hig-flow-vof-elvira.h"
#include "hig-flow-vof-adap-hf.h"
#include "hig-flow-step-multiphase.h"
#include "hig-flow-vof-finite-difference-normal-curvature.h"


void fraction_correction_at_get(real *fracvol);

void elvira_vertical_collumn(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *vertical, int *aux);

void elvira_horizontal_row(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *horizontal, int *aux);

void elvira_vertical_collumn_adap(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *vertical, int *aux, real *orig);

void elvira_horizontal_row_adap(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *horizontal, int *aux, real *orig);

void higflow_compute_normal_multiphase_2D_elvira(higflow_solver *ns, sim_domain *sdp, mp_mapper *mp, higcit_celliterator *it, hig_cell *c,int clid, Point center, Point delta, Point p);

void higflow_compute_normal_multiphase_2D_elvira_adap(higflow_solver *ns, sim_domain *sdp, mp_mapper *mp, higcit_celliterator *it, hig_cell *c,int clid, Point center, Point delta, Point p);

#endif
