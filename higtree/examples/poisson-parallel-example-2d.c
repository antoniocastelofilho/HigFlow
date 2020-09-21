
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "utils.h"
#include "higtree.h"
#include "higtree-io.h"
#include "mapper.h"
#include "wls.h"
#include "glib.h"
#include "pdomain.h"
#include "solver.h"
#include "lbal.h"

#define DEBUG
#include "Debug-c.h"

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
	exit(0);
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
	exit(0);
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

void __build_trees(int side, int refinementLevel, int isUniform, int bctypes, int refineFactor, int writevtk, int writebc, char vtkfilename[], load_balancer *lb) {
	Point l, h;
	POINT_ASSIGN_REALS(l, -1.0, -1.0);
	POINT_ASSIGN_REALS(h, 1.0, 1.0);
	hig_cell *root[1];
	int numdoms = 1;
	root[0] = hig_create_root(l, h);
	int nc[DIM];
	int part = 100;
	POINT_ASSIGN_INTS(nc, part, part);
	hig_refine_uniform(root[0], nc);

	POINT_ASSIGN_INTS(nc, refineFactor, refineFactor);
	higcit_celliterator *it;
	if (!isUniform) {
		__refine_non_uniform(root, numdoms, nc);
	}
	POINT_ASSIGN_INTS(nc, 2, 2);
	for (int i = 0; i < refinementLevel; i++) {
		for (int j = 0; j < numdoms; j++) {
			for(it = higcit_create_all_leaves(root[j]); !higcit_isfinished(it); higcit_nextcell(it)) {
				hig_cell *c = higcit_getcell(it);
				hig_refine_uniform(c, nc);
			}
			higcit_destroy(it);
		}
		side *= 2;
	}

	for(int j = 0; j < numdoms; j++) {
		it = higcit_create_all_leaves(root[j]);
		higcit_destroy(it);

		lb_add_input_tree(lb, root[j], true, 0);
	}

	if (writevtk) {
		higcit_celliterator *its[numdoms];
		for(unsigned j = 0; j < numdoms; ++j) {
			its[j] = higcit_create_all_leaves(root[j]);
		}
		it = higcit_create_concat(its, numdoms);
		FILE *fd = fopen(vtkfilename, "w");
		higio_print_celliterator_in_vtk2d(fd, it);
		higcit_destroy(it);
		fclose(fd);
	}

	Point ls[4], hs[4];
	int ncs[4][DIM];
	hig_cell *bcs[4];
	bc_type types[4];

	POINT_ASSIGN_REALS(ls[0], -1.0 - EPSDELTA, -1.0 - EPSDELTA);
	POINT_ASSIGN_REALS(hs[0],  1.0 + EPSDELTA, -1.0 + EPSDELTA);
	POINT_ASSIGN_INTS(ncs[0], part, 1);
	types[0] = DIRICHLET;

	POINT_ASSIGN_REALS(ls[1], -1.0 - EPSDELTA, -1.0 - EPSDELTA);
	POINT_ASSIGN_REALS(hs[1], -1.0  + EPSDELTA, 1.0 + EPSDELTA);
	POINT_ASSIGN_INTS(ncs[1], 1, part);
	types[1] = DIRICHLET;

	POINT_ASSIGN_REALS(ls[2], -1.0 - EPSDELTA, 1.0 - EPSDELTA);
	POINT_ASSIGN_REALS(hs[2],  1.0 + EPSDELTA, 1.0 + EPSDELTA);
	POINT_ASSIGN_INTS(ncs[2], part, 1);
	types[2] = DIRICHLET;

	POINT_ASSIGN_REALS(ls[3], 1.0 - EPSDELTA, -1.0 - EPSDELTA);
	POINT_ASSIGN_REALS(hs[3], 1.0 + EPSDELTA,  1.0 + EPSDELTA);
	POINT_ASSIGN_INTS(ncs[3], 1, part);
	types[3] = DIRICHLET;

	for(int i = 0; i < 4; i++) {
		bcs[i] = hig_create_root(ls[i], hs[i]);
		hig_refine_uniform(bcs[i], ncs[i]);
		lb_add_input_tree(lb, bcs[i], true, i+1);
	}
	if (writebc) {
		higcit_celliterator *its[4];
		for(unsigned j = 0; j < 4; ++j) {
			its[j] = higcit_create_all_leaves(bcs[j]);
		}
		it = higcit_create_concat(its, 4);
		FILE *fd = fopen(vtkfilename, "w");
		higio_print_celliterator_in_vtk2d(fd, it);
		higcit_destroy(it);
		fclose(fd);
	}
}

void __setup_bc(load_balancer *lb, sim_domain *sd)
{
	higcit_celliterator *it;
	for(int i = 0; i < 4; i++) {
		unsigned num_bcs = lb_get_num_local_trees_in_group(lb, i);
		for(unsigned j = 0; j < num_bcs; ++j) {
			hig_cell *bc = lb_get_local_tree_in_group(lb, j, i);
			sim_boundary *sb = sb_create(bc, DIRICHLET, NULL);
			sd_add_boundary(sd, sb);
			mp_mapper *m = sb_get_mapper(sb);

			for(it = higcit_create_all_leaves(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
				hig_cell *c = higcit_getcell(it);
				Point ccenter;
				hig_get_center(c, ccenter);
				int cgid = mp_lookup(m, hig_get_cid(c));
				sb_set_value(sb, cgid, f(ccenter[0], ccenter[1]));
			}
			higcit_destroy(it);
		}
	}

}

void __calc_poisson_stencil(sim_domain *sd, real cx, real cy, real dx, real dy, sim_stencil *stn) {
}

int main(int argc, char *argv[]) {
	int nc[DIM];
	int side = 30;
	int refinementLevel = 0;
	int isUniform = 1;
	int bctypes = 0;
	int writevtk = 0;
	int writebc = 0;
	int refineFactor = 2;
	char vtkfilename[100];
	DEBUG_DIFF_TIME;
	int myrank;
	int ntasks;

	higtree_initialize(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

	partition_graph *pg = pg_create(MPI_COMM_WORLD);
	pg_set_fringe_size(pg, 4);

	load_balancer *lb = lb_create(MPI_COMM_WORLD, 5);

	if (myrank == 0) {
		__get_opts(argc, argv, &side, &refinementLevel, &isUniform, &bctypes, &refineFactor, &writevtk, &writebc, vtkfilename, 100);
		__build_trees(side, refinementLevel, isUniform, bctypes, refineFactor, writevtk, writebc, vtkfilename, lb);

	}
	lb_calc_partition(lb, pg);
	sim_domain *sd = sd_create(NULL);
	{
		unsigned numlocal_dtrees = lb_get_num_local_trees_in_group(lb, 0);
		for(unsigned i = 0; i < numlocal_dtrees; ++i) {
			sd_add_higtree(sd, lb_get_local_tree_in_group(lb, i, 0));
		}
	}
	__setup_bc(lb, sd);
	lb_destroy(lb);

	psim_domain *psd = psd_create(sd, pg);
	psd_synced_mapper(psd);
	int localdomainsize = psd_get_local_domain_size(psd);

	/*if (myrank == 0) {
		char filename[1024];
		sprintf(filename, "/tmp/d.vtk");
		FILE *fd = fopen(filename, "w");
		higcit_celliterator *it = sd_get_domain_celliterator(psd_get_global_domain(psd));
		higio_print_celliterator_in_vtk2d(fd, it);
		fclose(fd);
		higcit_destroy(it);
	}*/
	char filename[1024];
	sprintf(filename, "/tmp/sbd-%d.vtk", myrank);
	FILE *fd = fopen(filename, "w");

	higcit_celliterator *it = sd_get_domain_celliterator(sd);
	higio_print_celliterator_in_vtk2d(fd, it);
	fclose(fd);
	higcit_destroy(it);
	sprintf(filename, "/tmp/sbfd-%d.vtk", myrank);
	fd = fopen(filename, "w");
	it = higcit_create_all_leaves(sd_get_higtree(sd, 0));
	higio_print_celliterator_in_vtk2d(fd, it);
	fclose(fd);
	higcit_destroy(it);

	sprintf(filename, "/tmp/bc-%d.vtk", myrank);
	fd = fopen(filename, "w");
	it = sd_get_bcs_celliterator(sd);
	higio_print_celliterator_in_vtk2d(fd, it);
	fclose(fd);
	higcit_destroy(it);

	sim_stencil *stn = stn_create();

	printf("generating matrix...\n");
	DEBUG_DIFF_TIME;
	solver *slv = slv_create(SOLVER_ANY, psd_get_first_id(psd), localdomainsize);
	slv_set_maxnonzeros(slv, 300);

	mp_mapper *m = sd->m;
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

		int *ids = psd_stn_get_gids(psd, stn);
		real *vals = stn_get_vals(stn);
		int numelems = stn_get_numelems(stn);

		int cgid = psd_get_global_id(psd, c);
		slv_set_Ai(slv, cgid, numelems, ids, vals);
		slv_set_bi(slv, cgid, stn_get_rhs(stn));
	}
	higcit_destroy(it);
	DEBUG_DIFF_TIME;
	printf("assembling...\n");
//	slv_assemble(slv);
//	DEBUG_DIFF_TIME;
//
//	printf("solving...\n");
//	int psize[2];
//	slv_solve(slv);
//
//	Vec fval;
//	VecDuplicate(slv->x, &fval);
//	MatView(slv->A, PETSC_VIEWER_DRAW_WORLD);
//	getchar();
//
//	real norm_inf = -1.0;
//	real norm_2 = 0.0;
//
//	for(int i = 0; i < 3; i++) {
//		hig_cell *root = sd_get_higtree(sd, i);
//		mp_mapper * m = sd_get_domain_mapper(sd);
//	for(it = sd_get_domain_celliterator(sd); !higcit_isfinished(it); higcit_nextcell(it)) {
//		hig_cell *c = higcit_getcell(it);
//		int cgid = mp_lookup(m, hig_get_id(c, 0));
//		Point ccenter;
//		hig_get_center(c, ccenter);
//		real cx = ccenter[0];
//		real cy = ccenter[1];
//		VecSetValue(fval, cgid, f(cx, cy), INSERT_VALUES);
//	}
//	VecAssemblyBegin(fval);
//	VecAssemblyEnd(fval);
//	VecAXPY(fval, -1.0, slv->x);
//	real norm_inf, norm_1, norm_2;
//	VecNorm(fval, NORM_1, &norm_1);
//	VecNorm(fval, NORM_2, &norm_2);
//	VecNorm(fval, NORM_INFINITY, &norm_inf);
//	if (myrank == 0) {
//		DEBUG_INSPECT(norm_1, %lg);
//		DEBUG_INSPECT(norm_2, %lg);
//		DEBUG_INSPECT(norm_inf, %lg);
//	}

	//VecView(fval, PETSC_VIEWER_STDOUT_WORLD);
//	}
//	norm_2 /= (real) size;
//	norm_2 = sqrt(norm_2);
//
//	printf("norm_inf = %f\n", norm_inf);
//	printf("norm_2 = %f\n", norm_2);
//	printf("{%e, %e, %e}\n", 2.0/(side*pow(2.0, refinementLevel)), norm_inf, norm_2);
//	higcit_destroy(it);
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
