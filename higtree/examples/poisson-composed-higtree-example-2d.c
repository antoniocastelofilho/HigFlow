
#include<stdio.h>
#include<stdlib.h>
#define DEBUG
#include "Debug-c.h"

#include "utils.h"
#include "higtree.h"
#include "higtree-io.h"
#include "mapper.h"
#include "wls.h"
#include "glib.h"
#include "domain.h"
#include "solver.h"

int func = 0;
real f(real x, real y) {
	switch(func) {
		case 0:	return sin(x)*cos(y);
		case 1: return 4*x*x*y+2*x*x-4*y*x-5*x;
		case 2: return cos(x*y)+sin(x+y);
		case 3: return cos(x-y)*sin(x+y);
		case 4: return cos(x+sin(y));
		case 5: return x*x+y*y-x*y-1.1234;
	}
	printf("undefined function! Should be in 0..4\n");
	exit(-1);
}

real g(real x, real y) {
	switch(func) {
		case 0: return -2.0*sin(x)*cos(y);
		case 1: return 4 + 8*y;
		case 2: return (-x*x-y*y)*cos(x*y)-2*sin(x+y);
		case 3: return -4*cos(x-y)*sin(x+y);
		case 4: return -(1+cos(y)*cos(y))*cos(x+sin(y))+sin(y)*(sin(x+sin(y)));
		case 5: return 4.0;
	}
}

real dfdy(real x, real y) {
	switch(func) {
		case 0: return -sin(x)*sin(y);
		case 1: return -4*x+4*x*x;
		case 2: return cos(x+y) - x*sin(x*y);
		case 3: return cos(x-y)*cos(x+y) + sin(x-y)*sin(x+y);
		case 4: return -cos(y)*sin(x+sin(y));
		case 5: return 2.0*y-x;
	}
}

real dfdx(real x, real y) {
	switch(func) {
		case 0: return cos(x)*cos(y);
		case 1: return -5 + 4*x - 4*y + 8*x*y;
		case 2: return cos(x+y) - y*sin(x*y);
		case 3: return cos(x-y)*cos(x+y) - sin(x-y)*sin(x+y);
		case 4: return -sin(x+sin(y));
		case 5: return 2.0*x-y;
	}
}

const char * f_str() {
	switch(func) {
		case 0:	return "sin(x)*cos(y)";
		case 1: return "4*x*x*y+2*x*x-4*y*x-5*x";
		case 2: return "cos(x*y)+sin(x+y)";
		case 3: return "cos(x-y)*sin(x+y)";
		case 4: return "cos(x+sin(y))";
		case 5: return "x*x+y*y-x*y-1.1234";
	}
	return "undefined function!";
}

void __refine_non_uniform(hig_cell *root[], int numdoms, int nc[]) {
	for(int i = 0; i < numdoms; i++) {
		Point bl, bh;
		POINT_ASSIGN_REALS(bl, -2.59, -0.59);
		POINT_ASSIGN_REALS(bh, 0.59, 0.59);
		higcit_celliterator *it;
		for(it = higcit_create_bounding_box(root[i], bl, bh); !higcit_isfinished(it); higcit_nextcell(it)) {
			hig_cell *c = higcit_getcell(it);
			hig_refine_uniform(c, nc);
		}
		higcit_destroy(it);
		POINT_ASSIGN_REALS(bl, -2.19, -0.19);
		POINT_ASSIGN_REALS(bh, 0.19, 0.19);
		for(it = higcit_create_bounding_box(root[i], bl, bh); !higcit_isfinished(it); higcit_nextcell(it)) {
			hig_cell *c = higcit_getcell(it);
			hig_refine_uniform(c, nc);
		}
		higcit_destroy(it);
	}
}

void __usage() {
	printf("./poisson-example \n\t[--side s (default 10)] \n\t[--refinement-level r (default 0)]\n\t[--uniform u (default: 1)]\n\t[--bctypes (default: DDDD)]\n\t[--writevtk vtkfilename | --writebc bcfilename]\n\n");
	exit(-1);
}

void __get_opts(int argc, char *argv[], int *side, int *refinementLevel, int *isUniform, int *bctypes, int *refineFactor, int *writevtk, int *writebc, char vtkfilename[], int maxlen) {
	while(argc > 2) {
		if(strcmp(argv[1], "--side") == 0) {
			argv++;
			argc--;
			*side = atoi(argv[1]);
			argv++;
			argc--;
		} else if(strcmp(argv[1], "--refinement-level") == 0) {
			argv++;
			argc--;
			*refinementLevel = atoi(argv[1]);
			argv++;
			argc--;
		} else if(strcmp(argv[1], "--uniform") == 0) {
			argv++;
			argc--;
			*isUniform = atoi(argv[1]);
			argv++;
			argc--;
		} else if(strcmp(argv[1], "--refinement-factor") == 0) {
			argv++;
			argc--;
			*refineFactor = atoi(argv[1]);
			argv++;
			argc--;
		} else if(strcmp(argv[1], "--bctypes") == 0) {
			argv++;
			argc--;
			*bctypes = 0;
			for(int i = 0; i < 4; i++) {
				if (argv[1][i] == 'N' || argv[1][i] == 'n') {
					*bctypes |= (1 << i);
				} else if (argv[1][i] != 'D' && argv[1][i] != 'd') {
					__usage();
				}
			}
			argv++;
			argc--;
		} else if(strcmp(argv[1], "--writevtk") == 0) {
			argv++;
			argc--;
			*writevtk = 1;
			strncpy(vtkfilename, argv[1], maxlen);
			argv++;
			argc--;
		} else if(strcmp(argv[1], "--writebc") == 0) {
			argv++;
			argc--;
			*writebc = 1;
			strncpy(vtkfilename, argv[1], maxlen);
			argv++;
			argc--;
		} else if(strcmp(argv[1], "--func") == 0) {
			argv++;
			argc--;
			func = atoi(argv[1]);
			argv++;
			argc--;
		} else {
			__usage();
		}
	}
}

void __print_opts(int side, int size, int refinementLevel, int isUniform, int bctypes, int refineFactor, int writevtk, int writebc, char vtkfilename[]) {
	printf("func = %s\n", f_str());
	printf("side = %d\n", side);
	printf("number of cells = %d\n", size);
	printf("refinement levels = %d\n", refinementLevel);
	printf("%s\n", (isUniform?"uniform":"non-uniform"));
	printf("refinement factor = %d\n", refineFactor);
	printf("bc types:\n\t  Left = %s\n\t Rigth = %s\n\tBottom = %s\n\t   Top = %s\n",
		((bctypes & (1<<0))?"NEUMANN":"DIRICHLET"),
		((bctypes & (1<<1))?"NEUMANN":"DIRICHLET"),
		((bctypes & (1<<2))?"NEUMANN":"DIRICHLET"),
		((bctypes & (1<<3))?"NEUMANN":"DIRICHLET"));
	printf("writevtk = %d\n", writevtk);
	printf("writebc = %d\n", writebc);
	if(writebc || writevtk) {
		printf("vtkfilename = %s\n", vtkfilename);
	}
}

sim_domain * __build_domain(int side, int *size, int refinementLevel, int isUniform, int bctypes, int refineFactor, int writevtk, int writebc, char vtkfilename[]) {
	Point l, h;
	POINT_ASSIGN_REALS(l, -0.2, -1.0);
	POINT_ASSIGN_REALS(h, 0.2, -0.6);
	hig_cell *root[3];
	root[0] = hig_create_root(l, h);
	int nc[DIM];
	int part = 8;
	POINT_ASSIGN_INTS(nc, 2*part, 2*part);
	hig_refine_uniform(root[0], nc);

	POINT_ASSIGN_REALS(l, -1.0, -0.6);
	POINT_ASSIGN_REALS(h, 1.0, 0.6);
	root[1] = hig_create_root(l, h);
	POINT_ASSIGN_INTS(nc, 10*part, 6*part);
	hig_refine_uniform(root[1], nc);

	POINT_ASSIGN_REALS(l, 0.2, 0.6);
	POINT_ASSIGN_REALS(h, 1.0, 1.0);
	root[2] = hig_create_root(l, h);
	POINT_ASSIGN_INTS(nc, 4*part, 2*part);
	hig_refine_uniform(root[2], nc);

	POINT_ASSIGN_INTS(nc, refineFactor, refineFactor);
	higcit_celliterator *it;
	if (!isUniform) {
		__refine_non_uniform(root, 3, nc);
	}
	POINT_ASSIGN_INTS(nc, 2, 2);
	for (int i = 0; i < refinementLevel; i++) {
		for (int j = 0; j < 3; j++) {
			for(it = higcit_create_all_leaves(root[j]); !higcit_isfinished(it); higcit_nextcell(it)) {
				hig_cell *c = higcit_getcell(it);
				hig_refine_uniform(c, nc);
			}
			higcit_destroy(it);
		}
		side *= 2;
	}
	mp_mapper *m = mp_create();
	sim_domain *sd = sd_create(m);
	*size = 0;
	for(int j = 0; j < 3; j++) {
		it = higcit_create_all_leaves(root[j]);
		*size = mp_assign_from_celliterator(m, it, *size);
		higcit_destroy(it);

		sd_add_higtree(sd, root[j]);
	}

	if (writevtk) {
		FILE *fd = fopen(vtkfilename, "w");
		it = sd_get_domain_celliterator(sd);
		higio_print_celliterator_in_vtk2d(fd, it);
		higcit_destroy(it);
		fclose(fd);
	}

	Point ls[10], hs[10];
	int ncs[10][DIM];
	hig_cell *bcs[10];
	mp_mapper *ms[10];
	sim_boundary *sbs[10];
	bc_type types[10];

	POINT_ASSIGN_REALS(ls[0], -1.0 - EPSDELTA, -0.6 - EPSDELTA);
	POINT_ASSIGN_REALS(hs[0], -1.0 + EPSDELTA, 0.6  + EPSDELTA);
	POINT_ASSIGN_INTS(ncs[0], 1, 6*part);
	types[0] = DIRICHLET;

	POINT_ASSIGN_REALS(ls[1], -1.0 - EPSDELTA, 0.6 - EPSDELTA);
	POINT_ASSIGN_REALS(hs[1], 0.2  + EPSDELTA, 0.6  + EPSDELTA);
	POINT_ASSIGN_INTS(ncs[1], 6*part, 1);
	types[1] = DIRICHLET;

	POINT_ASSIGN_REALS(ls[2], 0.2  - EPSDELTA, 0.6 - EPSDELTA);
	POINT_ASSIGN_REALS(hs[2], 0.2  + EPSDELTA, 1.0 + EPSDELTA);
	POINT_ASSIGN_INTS(ncs[2], 1, 2*part);
	types[2] = DIRICHLET;

	POINT_ASSIGN_REALS(ls[3], 0.2  - EPSDELTA, 1.0 - EPSDELTA);
	POINT_ASSIGN_REALS(hs[3], 1.0  + EPSDELTA, 1.0 + EPSDELTA);
	POINT_ASSIGN_INTS(ncs[3], 4*part, 1);
	types[3] = DIRICHLET;

	POINT_ASSIGN_REALS(ls[4], 1.0  - EPSDELTA, -0.6 - EPSDELTA);
	POINT_ASSIGN_REALS(hs[4], 1.0  + EPSDELTA, 1.0 + EPSDELTA);
	POINT_ASSIGN_INTS(ncs[4], 1, 8*part);
	types[4] = DIRICHLET;

	POINT_ASSIGN_REALS(ls[5], 0.2  - EPSDELTA, -0.6 - EPSDELTA);
	POINT_ASSIGN_REALS(hs[5], 1.0  + EPSDELTA, -0.6 + EPSDELTA);
	POINT_ASSIGN_INTS(ncs[5], 4*part, 1);
	types[5] = DIRICHLET;

	POINT_ASSIGN_REALS(ls[6], 0.2  - EPSDELTA, -1.0 - EPSDELTA);
	POINT_ASSIGN_REALS(hs[6], 0.2  + EPSDELTA, -0.6 + EPSDELTA);
	POINT_ASSIGN_INTS(ncs[6], 1, 2*part);
	types[6] = DIRICHLET;

	POINT_ASSIGN_REALS(ls[7], -0.2 - EPSDELTA, -1.0 - EPSDELTA);
	POINT_ASSIGN_REALS(hs[7], 0.2  + EPSDELTA, -1.0 + EPSDELTA);
	POINT_ASSIGN_INTS(ncs[7], 2*part, 1);
	types[7] = DIRICHLET;

	POINT_ASSIGN_REALS(ls[8], -0.2 - EPSDELTA, -1.0 - EPSDELTA);
	POINT_ASSIGN_REALS(hs[8], -0.2 + EPSDELTA, -0.6 + EPSDELTA);
	POINT_ASSIGN_INTS(ncs[8], 1, 2*part);
	types[8] = DIRICHLET;

	POINT_ASSIGN_REALS(ls[9], -1.0 - EPSDELTA, -0.6 - EPSDELTA);
	POINT_ASSIGN_REALS(hs[9], -0.2 + EPSDELTA, -0.6 + EPSDELTA);
	POINT_ASSIGN_INTS(ncs[9], 4*part, 1);
	types[9] = DIRICHLET;

	for(int i = 0; i < 10; i++) {
		bcs[i] = hig_create_root(ls[i], hs[i]);
		hig_refine_uniform(bcs[i], ncs[i]);
		ms[i] = mp_create();
		it = higcit_create_all_leaves(bcs[i]);
		mp_assign_from_celliterator(ms[i], it, 0);
		higcit_destroy(it);
		sbs[i] = sb_create(bcs[i], types[i], ms[i]);
		sd_add_boundary(sd, sbs[i]);
	}

	for(int i = 0; i < 10; i++) {
		for(it = higcit_create_all_leaves(bcs[i]); !higcit_isfinished(it); higcit_nextcell(it)) {
			hig_cell *c = higcit_getcell(it);
			Point ccenter;
			hig_get_center(c, ccenter);
			int cgid = mp_lookup(ms[i], hig_get_cid(c));
			if(types[i] == DIRICHLET) {
				sb_set_value(sbs[i], cgid, f(ccenter[0], ccenter[1]));
			} else {
				switch (i) {
					case 0: sb_set_value(sbs[i], cgid, -dfdx(ccenter[0], ccenter[1]));
						break;
					case 1: sb_set_value(sbs[i], cgid, dfdx(ccenter[0], ccenter[1]));
						break;
					case 2: sb_set_value(sbs[i], cgid, -dfdy(ccenter[0], ccenter[1]));
						break;
					case 3: sb_set_value(sbs[i], cgid, dfdy(ccenter[0], ccenter[1]));
						break;
					default : exit(-1);
				}
			}
		}
		higcit_destroy(it);
	}
	if (writebc) {
		FILE *fd = fopen(vtkfilename, "w");
		it = sd_get_bcs_celliterator(sd);
		higio_print_celliterator_in_vtk2d(fd, it);
		higcit_destroy(it);
		fclose(fd);
	}

	return sd;
}


void __calc_poisson_stencil(sim_domain *sd, real cx, real cy, real dx, real dy, sim_stencil *stn) {
}

int main(int argc, char *argv[]) {
	int nc[DIM];
	int side = 10;
	int refinementLevel = 0;
	int isUniform = 1;
	int bctypes = 0;
	int writevtk = 0;
	int writebc = 0;
	int refineFactor = 2;
	int size;
	char vtkfilename[100];

	higtree_initialize(&argc, &argv);
	DEBUG_DIFF_TIME;
	__get_opts(argc, argv, &side, &refinementLevel, &isUniform, &bctypes, &refineFactor, &writevtk, &writebc, vtkfilename, 100);

	sim_domain *sd = __build_domain(side, &size, refinementLevel, isUniform, bctypes, refineFactor, writevtk, writebc, vtkfilename);

	__print_opts(side, size, refinementLevel, isUniform, bctypes, refineFactor, writevtk, writebc, vtkfilename);

	sim_stencil *stn = stn_create();

	printf("generating matrix...\n");
	DEBUG_DIFF_TIME;
	solver *slv = slv_create(SOLVER_ANY, 0, size);
	slv_set_maxnonzeros(slv, 300);

	higcit_celliterator *it;
	mp_mapper * m = sd_get_domain_mapper(sd);
	sd_set_interpolator_order(sd, 3);
	for(it = sd_get_domain_celliterator(sd); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		Point ccenter;
		Point cdelta;
		hig_get_center(c, ccenter);
		real cx = ccenter[0];
		real cy = ccenter[1];
		hig_get_delta(c, cdelta);
		real dx = cdelta[0];
		real dy = cdelta[1];
		dy *= 1.0;

		Point center = {cx     , cy     };
		Point left   = {cx - dx, cy     };
		Point right  = {cx + dx, cy     };
		Point top    = {cx,      cy + dy};
		Point bottom = {cx,      cy - dy};

		stn_reset(stn);
		stn_set_rhs(stn, g(cx, cy));
		sd_get_stencil(sd, center, center, (-2.0/(dx*dx)-2.0/(dy*dy)), stn);
		sd_get_stencil(sd, center, left,    (1.0/(dx*dx)),              stn);
		sd_get_stencil(sd, center, right,   (1.0/(dx*dx)),              stn);
		sd_get_stencil(sd, center, top,     (1.0/(dy*dy)),              stn);
		sd_get_stencil(sd, center, bottom,  (1.0/(dy*dy)),              stn);

		int *ids = stn_get_ids(stn);
		real *vals = stn_get_vals(stn);
		int numelems = stn_get_numelems(stn);

		int cgid = mp_lookup(m, hig_get_cid(c));
		slv_set_Ai(slv, cgid, numelems, ids, vals);
		slv_set_bi(slv, cgid, stn_get_rhs(stn));
	}
	higcit_destroy(it);
	DEBUG_DIFF_TIME;
	printf("assembling...\n");
	slv_assemble(slv);
	DEBUG_DIFF_TIME;

	printf("solving...\n");
	slv_solve(slv);
	DEBUG_DIFF_TIME;

	real norm_inf = -1.0;
	real norm_2 = 0.0;

	for(int i = 0; i < 3; i++) {
		hig_cell *root = sd_get_higtree(sd, i);
		mp_mapper * m = sd_get_domain_mapper(sd);
		for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
			hig_cell *c = higcit_getcell(it);
			int cgid = mp_lookup(m, hig_get_cid(c));
			Point ccenter;
			hig_get_center(c, ccenter);
			real cx = ccenter[0];
			real cy = ccenter[1];

			real error = fabs(f(cx, cy) - slv_get_xi(slv, cgid));
			norm_2 += error*error;
			norm_inf = ((norm_inf<error) ? error : norm_inf);
		}
	}
	norm_2 /= (real) size;
	norm_2 = sqrt(norm_2);

	printf("norm_inf = %f\n", norm_inf);
	printf("norm_2 = %f\n", norm_2);
	printf("{%e, %e, %e}\n", 2.0/(side*pow(2.0, refinementLevel)), norm_inf, norm_2);
	higcit_destroy(it);
	return 0;
}

/*
		Point center = {cx         , cy     };
		Point px1    = {cx - 1.0*dx, cy};
		Point px2    = {cx + 0.5*dx, cy};
		Point px3    = {cx + 2.0*dx, cy};
		Point py1    = {cx,          cy - 1.0*dy};
		Point py2    = {cx,          cy + 0.5*dy};
		Point py3    = {cx,          cy + 2.0*dy};

		stn_reset(stn);
		stn_set_rhs(stn, g(cx, cy));
		real wc = -27.0/9.0, w1 = 10.0/9.0, w2 = 16.0/9.0, w3 = 1.0/9.0;

		sd_get_stencil(sd, center, center,  (wc/(dx*dx)), stn);
		sd_get_stencil(sd, center, px1,     (w1/(dx*dx)),  stn);
		sd_get_stencil(sd, center, px2,     (w2/(dx*dx)),  stn);
		sd_get_stencil(sd, center, px3,     (w3/(dx*dx)), stn);
		sd_get_stencil(sd, center, center,  (wc/(dy*dy)), stn);
		sd_get_stencil(sd, center, py1,     (w1/(dy*dy)),  stn);
		sd_get_stencil(sd, center, py2,     (w2/(dy*dy)),  stn);
		sd_get_stencil(sd, center, py3,     (w3/(dy*dy)), stn);
*/

/*
		Point center = {cx     , cy     };
		Point left   = {cx - dx, cy     };
		Point right  = {cx + dx, cy     };
		Point p1    = {cx,       cy - 1.0*dy};
		Point p2    = {cx,       cy + 2.0*dy};
		Point p3    = {cx,       cy + 3.0*dy};

		stn_reset(stn);
		stn_set_rhs(stn, g(cx, cy));
		sd_get_stencil(sd, center, center,  (-2.0/(dx*dx)), stn);
		sd_get_stencil(sd, center, left,    (1.0/(dx*dx)),  stn);
		sd_get_stencil(sd, center, right,   (1.0/(dx*dx)),  stn);
		sd_get_stencil(sd, center, center,  ((-8.0/6.0)/(dy*dy)), stn);
		sd_get_stencil(sd, center, p1,  ((5.0/6.0)/(dy*dy)),  stn);
		sd_get_stencil(sd, center, p2,  ((4.0/6.0)/(dy*dy)),  stn);
		sd_get_stencil(sd, center, p3,  ((-1.0/6.0)/(dy*dy)),  stn);
*/

/*
		Point center = {cx     , cy     };
		Point left   = {cx - dx, cy     };
		Point right  = {cx + dx, cy     };
		Point top    = {cx,      cy + dy};
		Point bottom = {cx,      cy - dy};

		stn_reset(stn);
		stn_set_rhs(stn, g(cx, cy));
		sd_get_stencil(sd, center, center, (-2.0/(dx*dx)-2.0/(dy*dy)), stn);
		sd_get_stencil(sd, center, left,    (1.0/(dx*dx)),              stn);
		sd_get_stencil(sd, center, right,   (1.0/(dx*dx)),              stn);
		sd_get_stencil(sd, center, top,     (1.0/(dy*dy)),              stn);
		sd_get_stencil(sd, center, bottom,  (1.0/(dy*dy)),              stn);
*/

/*
		Point center = {cx     , cy     };
		Point left   = {cx - dx, cy     };
		Point right  = {cx + dx, cy     };
		Point top    = {cx,      cy + dy};
		Point bottom = {cx,      cy - dy};

		stn_reset(stn);
		stn_set_rhs(stn, g(cx, cy));
		sd_get_stencil(sd, center, center,  (-2.0/(dx*dx)), stn);
		sd_get_stencil(sd, center, left,    (1.0/(dx*dx)),  stn);
		sd_get_stencil(sd, center, right,   (1.0/(dx*dx)),  stn);
		sd_get_stencil(sd, center, center,  (-2.0/(dy*dy)), stn);
		sd_get_stencil(sd, center, top,     (1.0/(dy*dy)),  stn);
		sd_get_stencil(sd, center, bottom,  (1.0/(dy*dy)),  stn);
*/

/*
		Point center = {cx     , cy     };
		Point left   = {cx - dx, cy     };
		Point right  = {cx + dx, cy     };
		Point top    = {cx,      cy + 2.0*dy};
		Point bottom = {cx,      cy - 2.0*dy};

		stn_reset(stn);
		stn_set_rhs(stn, g(cx, cy));
		sd_get_stencil(sd, center, center,  (-2.0/(dx*dx)), stn);
		sd_get_stencil(sd, center, left,    (1.0/(dx*dx)),  stn);
		sd_get_stencil(sd, center, right,   (1.0/(dx*dx)),  stn);
		sd_get_stencil(sd, center, center,  (-2.0/(4.0*dy*dy)), stn);
		sd_get_stencil(sd, center, top,     (1.0/(4.0*dy*dy)),  stn);
		sd_get_stencil(sd, center, bottom,  (1.0/(4.0*dy*dy)),  stn);
*/
