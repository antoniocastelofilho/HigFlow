
#include<stdio.h>
#include<stdlib.h>

#include "utils.h"
#include "higtree.h"
#include "mapper.h"
#include "wls.h"
#include "glib.h"
#include "higtree-io.h"

#define DEBUG
#include "Debug-c.h"

int is_boundary(hig_cell *root, hig_cell *c) {
	Point rlow, rhigh;
	hig_get_lowpoint(root, rlow);
	hig_get_highpoint(root, rhigh);

	Point clow, chigh;
	hig_get_lowpoint(c, clow);
	hig_get_highpoint(c, chigh);
	for(int i = 0; i < DIM; i++) {
		if (FLT_EQ(clow[i], rlow[i])) {
			return 1;
		}
		if (FLT_EQ(chigh[i], rhigh[i])) {
			return 1;
		}
	}
	return 0;
}

#if DIM == 2
real func(Point x) {
	return x[0]*x[0] + x[1]*x[1] - 2.0*x[0]*x[1];
	//return cos(x[0])*cos(x[1]);
}
#elif DIM == 3
real func(Point x) {
//	return x[0]*x[0] + x[1]*x[1] - 2.0*x[0]*x[1] + 2.0*x[2]*x[2];
	return cos(x[0])*cos(x[1])*cos(x[2]);
}
#endif

int main(int argc, char *argv[]) {
	if (argc <= 1) {
		printf("usage: %s amr2d.file\n", argv[0]);
		exit(0);
	}
	FILE * fdin = fopen(argv[1], "r");
#if DIM == 2
	hig_cell * root = higio_read_from_amr2d(fdin);
#elif DIM == 3
	hig_cell * root = higio_read_from_amr3d(fdin);
#endif
	fclose(fdin);

	higcit_celliterator *it;
	mp_mapper *m = mp_create();

	it = higcit_create_all_leaves(root);
	mp_assign_from_celliterator(m, it, 0);
	higcit_destroy(it);

	int numleaves = hig_get_number_of_leaves(root);

	real val[numleaves];
	real error[numleaves];
	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		Point center;
		hig_get_center(c, center);
		uniqueid cid = hig_get_cid(c);
		mp_value_t cgid = mp_lookup(m, cid);
		val[cgid] = func(center);
	}

	Point low, high;

	const int maxpts = 2 * DIM * wls_num_min_points(2);
	Point pts[maxpts];
	real w[maxpts];
	uniqueid gids[maxpts];
	real norm_inf = -10000.0;
	real minerror = 10000.0;

	DEBUG_DIFF_TIME;

	wls_interpolator *wls = wls_create(2, maxpts);

	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		uniqueid cid = hig_get_cid(c);
		int cgid = mp_lookup(m, cid);
		if (is_boundary(root, c)) {
			error[cgid] = 0.0;
			continue;
		}
		Point ccenter;
		hig_get_center(c, ccenter);
		higcit_celliterator *itn;
		int numpts = 0;
		for(itn = higcit_create_neighbours(c); !higcit_isfinished(itn) && numpts < maxpts; higcit_nextcell(itn), numpts++) {
			hig_cell *n = higcit_getcell(itn);
			Point ncenter;
			hig_get_center(n, ncenter);
			uniqueid nid = hig_get_cid(n);
			int ngid = mp_lookup(m, nid);
			gids[numpts] = ngid;
			POINT_ASSIGN(pts[numpts], ncenter);
		}
		gids[numpts] = cgid;
		POINT_ASSIGN(pts[numpts], ccenter);
		numpts++;
		higcit_destroy(itn);
		Point highpoint;
		hig_get_highpoint(c, highpoint);
		wls_set_points(wls, numpts, pts, highpoint);
		wls_calc(wls, highpoint, numpts, w);
		real fval = 0.0;
		for(int i = 0; i < numpts; i++) {
			fval += w[i] * val[gids[i]];
		}
		real err = fabs(fval - func(highpoint));
		if(norm_inf < err) {
			norm_inf = err;
		}
		if(minerror > err) {
			minerror = err;
		}
		error[cgid] = err;
	}

	DEBUG_DIFF_TIME;

	higcit_destroy(it);

	printf("%g %g\n", minerror, norm_inf);

	printf("Infity Norm: %g\n\n", norm_inf);

	hig_destroy(root);
}
