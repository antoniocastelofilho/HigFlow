#ifndef HIG_FLOW_VOF_9_CELLS
#define HIG_FLOW_VOF_9_CELLS

#include "hig-flow-discret.h"

void shirani_9_cells(sim_domain *sdp, higflow_solver *ns, int clid, Point center, Point p, Point delta);

#endif
