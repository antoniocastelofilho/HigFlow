
#include<stdio.h>
#include<stdlib.h>

#include "utils.h"
#include "higtree.h"
#include "higtree-io.h"
#include "mapper.h"
#include "wls.h"
#include "glib.h"
#include "domain.h"
#include "solver.h"

#define DEBUG
#include "Debug-c.h"

int func = 0;
#if DIM == 2
real f(Point p) {
	real x = p[0];
	real y = p[1];
	switch(func) {
		case 0:	return sin(x)*cos(y);
		case 1: return 4*x*x*y+2*x*x-4*y*x-5*x;
		case 2: return cos(x*y)+sin(x+y);
		case 3: return cos(x-y)*sin(x+y);
		case 4: return cos(x+sin(y));
		case 5: return x*x+y*y-x*y-1.1234;
	}
	printf("undefined function! Should be in 0..4\n");
	exit(0);
}

real g(Point p) {
	real x = p[0];
	real y = p[1];
	switch(func) {
		case 0: return -2.0*sin(x)*cos(y);
		case 1: return 4 + 8*y;
		case 2: return (-x*x-y*y)*cos(x*y)-2*sin(x+y);
		case 3: return -4*cos(x-y)*sin(x+y);
		case 4: return -(1+cos(y)*cos(y))*cos(x+sin(y))+sin(y)*(sin(x+sin(y)));
		case 5: return 4.0;
	}
}

real dfdy(Point p) {
	real x = p[0];
	real y = p[1];
	switch(func) {
		case 0: return -sin(x)*sin(y);
		case 1: return -4*x+4*x*x;
		case 2: return cos(x+y) - x*sin(x*y);
		case 3: return cos(x-y)*cos(x+y) + sin(x-y)*sin(x+y);
		case 4: return -cos(y)*sin(x+sin(y));
		case 5: return 2.0*y-x;
	}
}

real dfdx(Point p) {
	real x = p[0];
	real y = p[1];
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
#elif DIM == 3
real f(Point p) {
	real x = p[0];
	real y = p[1];
	real z = p[2];
	switch(func) {
		case 0:	return cos(x+z)*sin(y+z);
		case 1:	return x+y+z;
	}
	printf("undefined function! Should be in 0..4\n");
	exit(0);
}

real g(Point p) {
	real x = p[0];
	real y = p[1];
	real z = p[2];
	switch(func) {
		case 0: return -2.0*cos(y + z)*sin(x + z) - 4.0*cos(x + z)*sin(y + z);
		case 1: return 0.0;
	}
}

real dfdx(Point p) {
	real x = p[0];
	real y = p[1];
	real z = p[2];
	switch(func) {
		case 0: return -(sin(x + z)*sin(y + z));
	}
}

real dfdy(Point p) {
	real x = p[0];
	real y = p[1];
	real z = p[2];
	switch(func) {
		case 0: return cos(x + z)*cos(y + z);
	}
}

real dfdz(Point p) {
	real x = p[0];
	real y = p[1];
	real z = p[2];
	switch(func) {
		case 0: return cos(x + z)*cos(y + z) - sin(x + z)*sin(y + z);
	}
}

const char * f_str() {
	switch(func) {
		case 0:	return "cos(x + z)*sin(y + z)";
	}
	return "undefined function!";
}
#endif

void __refine_non_uniform(hig_cell *root, int nc[]) {
		Point bl, bh;
		POINT_ASSIGN_SCALAR(bl, -0.59);
		bl[0] = -2.59;
		POINT_ASSIGN_SCALAR(bh, 0.59);
		higcit_celliterator *it;
		for(it = higcit_create_bounding_box(root, bl, bh); !higcit_isfinished(it); higcit_nextcell(it)) {
			hig_cell *c = higcit_getcell(it);
			hig_refine_uniform(c, nc);
		}
		higcit_destroy(it);
		POINT_ASSIGN_SCALAR(bl, -0.19);
		bl[0] = -2.19;
		POINT_ASSIGN_SCALAR(bh, 0.19);
		for(it = higcit_create_bounding_box(root, bl, bh); !higcit_isfinished(it); higcit_nextcell(it)) {
			hig_cell *c = higcit_getcell(it);
			hig_refine_uniform(c, nc);
		}
		higcit_destroy(it);
}

void __usage() {
	printf("./poisson-example \n\t[--side s (default 10)] \n\t[--refinement-level r (default 0)]\n\t[--uniform u (default: 1)]\n\t[--bctypes (default: DDDD)]\n\t[--writevtk vtkfilename]\n\n");
	//exit(0);
}

void __get_opts(int argc, char *argv[], int *side, int *refinementLevel, int *isUniform, int *bctypes, int *refineFactor, int *writevtk, char vtkfilename[], int maxlen) {
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

void __print_opts(int side, int size, int refinementLevel, int isUniform, int bctypes, int refineFactor, int writevtk, char vtkfilename[]) {
	printf("func = %s\n", f_str());
	printf("side = %d\n", side);
	printf("number of cells = %d\n", size);
	printf("refinement levels = %d\n", refinementLevel);
	printf("%s\n", (isUniform?"uniform":"non-uniform"));
	printf("refinement factor = %d\n", refineFactor);
#if DIM == 2
	printf("bc types:\n\t  Left = %s\n\t Rigth = %s\n\tBottom = %s\n\t   Top = %s\n",
		((bctypes & (1<<0))?"NEUMANN":"DIRICHLET"),
		((bctypes & (1<<1))?"NEUMANN":"DIRICHLET"),
		((bctypes & (1<<2))?"NEUMANN":"DIRICHLET"),
		((bctypes & (1<<3))?"NEUMANN":"DIRICHLET"));
#elif DIM == 3
	printf("bc types:\n\t  Left = %s\n\t Rigth = %s\n\tBottom = %s\n\t   Top = %s\n\tBack = %s\n\tFront = %s\n",
		((bctypes & (1<<0))?"NEUMANN":"DIRICHLET"),
		((bctypes & (1<<1))?"NEUMANN":"DIRICHLET"),
		((bctypes & (1<<2))?"NEUMANN":"DIRICHLET"),
		((bctypes & (1<<3))?"NEUMANN":"DIRICHLET"),
		((bctypes & (1<<4))?"NEUMANN":"DIRICHLET"),
		((bctypes & (1<<5))?"NEUMANN":"DIRICHLET"));
#endif
	if(writevtk) {
		printf("vtkfilename = %s\n", vtkfilename);
	} else {
		printf("no vtk\n");
	}
}

sim_domain * __build_domain(int side, int *size, int refinementLevel, int isUniform, int bctypes, int refineFactor, int writevtk, char vtkfilename[]) {
	DEBUG_INSPECT(DIM, %d);
	Point l, h;
	POINT_ASSIGN_SCALAR(l, -1.0);
	POINT_ASSIGN_SCALAR(h, 1.0);
	hig_cell *root = hig_create_root(l, h);
	DEBUG_PASS;
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, side);
	hig_refine_uniform(root, nc);
	DEBUG_PASS;
	POINT_ASSIGN_SCALAR(nc, refineFactor);
	higcit_celliterator *it;
	if (!isUniform) {
		__refine_non_uniform(root, nc);
	}
	POINT_ASSIGN_SCALAR(nc, 2);
	for (int i = 0; i < refinementLevel; i++) {
		for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
			hig_cell *c = higcit_getcell(it);
			hig_refine_uniform(c, nc);
		}
		higcit_destroy(it);
		side *= 2;
	}
	if (writevtk) {
		FILE *fd = fopen(vtkfilename, "w");
#if DIM == 2
		higio_print_in_vtk2d(fd, root);
#elif DIM == 3
		higio_print_in_vtk3d(fd, root);
#endif
		fclose(fd);
	}

	mp_mapper *m = mp_create();
	it = higcit_create_all_leaves(root);
	*size = mp_assign_from_celliterator(m, it, 0);
	higcit_destroy(it);

	sim_domain *sd = sd_create(m);
	sd_add_domain(sd, root);

#define TWODIM (2*DIM)
	Point ls[TWODIM], hs[TWODIM];
	int ncs[TWODIM][DIM];
	hig_cell *bcs[TWODIM];
	mp_mapper *ms[TWODIM];
	sim_boundary *sbs[TWODIM];
	bc_type types[TWODIM];

	int k = 0;
	for(int i = 0; i < DIM; i++) {
		for(int j = 0; j < 2; j++) {
			POINT_ASSIGN(ls[k], l);
			ls[k][i] = ((j == 0)?l[i]:h[i]) - EPSDELTA;
			POINT_ASSIGN(hs[k], h);
			hs[k][i] = ((j == 0)?l[i]:h[i]) + EPSDELTA;
			POINT_ASSIGN_SCALAR(ncs[k], side);
			ncs[k][i] = 1;
			types[k] = ((bctypes&(1<<k))?NEUMANN:DIRICHLET);
			k++;
		}
	}

	for(int i = 0; i < TWODIM; i++) {
		bcs[i] = hig_create_root(ls[i], hs[i]);
		hig_refine_uniform(bcs[i], ncs[i]);
		POINT_ASSIGN_SCALAR(nc, refineFactor);
		nc[i/2]=1;
		char filename[1024];
		sprintf(filename, "/tmp/bc%d.vtk", i);
			FILE *fd = fopen(filename, "w");
#if DIM == 2
		higio_print_in_vtk2d(fd, bcs[i]);
#elif DIM == 3
		higio_print_in_vtk3d(fd, bcs[i]);
#endif
			fclose(fd);
		if (!isUniform && i == 0) {
			__refine_non_uniform(bcs[i], nc);
		}
		ms[i] = mp_create();
		it = higcit_create_all_leaves(bcs[i]);
		mp_assign_from_celliterator(ms[i], it, 0);
		higcit_destroy(it);
		sbs[i] = sb_create(bcs[i], types[i], ms[i]);
		sd_add_boundary(sd, sbs[i]);
	}

	for(int i = 0; i < TWODIM; i++) {
		for(it = higcit_create_all_leaves(bcs[i]); !higcit_isfinished(it); higcit_nextcell(it)) {
			hig_cell *c = higcit_getcell(it);
			Point ccenter;
			hig_get_center(c, ccenter);
			int cgid = mp_lookup(ms[i], hig_get_id(c, 0));
			if(types[i] == DIRICHLET) {
				sb_set_value(sbs[i], cgid, f(ccenter));
			} else {
				switch (i) {
					case 0: sb_set_value(sbs[i], cgid, -dfdx(ccenter));
						break;
					case 1: sb_set_value(sbs[i], cgid, dfdx(ccenter));
						break;
					case 2: sb_set_value(sbs[i], cgid, -dfdy(ccenter));
						break;
					case 3: sb_set_value(sbs[i], cgid, dfdy(ccenter));
						break;
#if DIM == 3
					case 4: sb_set_value(sbs[i], cgid, -dfdz(ccenter));
						break;
					case 5: sb_set_value(sbs[i], cgid, dfdz(ccenter));
#endif
					default : exit(0);
				}
			}
		}
		higcit_destroy(it);
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
	int refineFactor = 2;
	int size;
	char vtkfilename[100];

	higtree_initialize(&argc, &argv);
	DEBUG_DIFF_TIME;
	__get_opts(argc, argv, &side, &refinementLevel, &isUniform, &bctypes, &refineFactor, &writevtk, vtkfilename, 100);

	sim_domain *sd = __build_domain(side, &size, refinementLevel, isUniform, bctypes, refineFactor, writevtk, vtkfilename);

	__print_opts(side, size, refinementLevel, isUniform, bctypes, refineFactor, writevtk, vtkfilename);

	sim_stencil *stn = stn_create();

	printf("generating matrix...\n");
	DEBUG_DIFF_TIME;
	solver *slv = slv_create(SOLVER_ANY, 0, size);
	slv_set_maxnonzeros(slv, 20);

	higcit_celliterator *it;
	hig_cell *root = sd_get_domain(sd, 0);
	mp_mapper * m = sd_get_domain_mapper(sd);
	sd_set_interpolator_order(sd, 3);
	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		Point ccenter;
		Point cdelta;
		hig_get_center(c, ccenter);
		hig_get_delta(c, cdelta);
		Point center;
		POINT_ASSIGN(center, ccenter);

		stn_reset(stn);
		stn_set_rhs(stn, g(center));
		sd_set_interpolator_center(sd, center);
		real alphacenter = 0.0;
		for(int i = 0; i < DIM; i++) {
			for(int j = -1; j <= 1; j +=2) {
				Point p;
				POINT_ASSIGN(p, center);
				p[i] += ((real)j)*cdelta[i];
				real alphap = 1.0/(cdelta[i]*cdelta[i]);
				alphacenter -= alphap;
				sd_get_stencil(sd, p, alphap, stn);
			}
		}

		sd_get_stencil(sd, center, alphacenter, stn);

		int *ids = stn_get_ids(stn);
		real *vals = stn_get_vals(stn);
		int numelems = stn_get_numelems(stn);

		int cgid = mp_lookup(m, hig_get_id(c, 0));
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

	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		int cgid = mp_lookup(m, hig_get_id(c, 0));
		Point ccenter;
		hig_get_center(c, ccenter);
		real cx = ccenter[0];
		real cy = ccenter[1];

		real error = fabs(f(ccenter) - slv_get_xi(slv, cgid));
		norm_2 += error*error;
		norm_inf = ((norm_inf<error) ? error : norm_inf);
	}
	norm_2 /= (real) size;
	norm_2 = sqrt(norm_2);

	printf("norm_inf = %f\n", norm_inf);
	printf("norm_2 = %f\n", norm_2);
	printf("{%e, %e, %e}\n", 2.0/(side*pow(2.0, refinementLevel)), norm_inf, norm_2);
	higcit_destroy(it);
	slv_destroy(slv);
	sd_destroy(sd);
}

