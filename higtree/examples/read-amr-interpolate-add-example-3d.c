
#include<stdio.h>
#include<stdlib.h>

#include "utils.h"
#include "higtree.h"
#include "mapper.h"
#include "wls.h"
#include "glib.h"
#include "higtree-io.h"

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
	//return x[0]*x[0] + x[1]*x[1] - 2.0*x[0]*x[1];
	return cos(x[0])*cos(x[1]);
}
#elif DIM == 3
real func(Point x) {
	//return x[0]*x[0] + x[1]*x[1] - 2.0*x[0]*x[1] + 2.0*x[2]*x[2];
	return cos(x[0])*cos(x[1])*cos(x[2]);
}
#endif

int main(int argc, char *argv[]) {
	if (argc <= 2) {
		printf("usage: %s amr2d.file numptsadd\n", argv[0]);
		exit(-1);
	}
	FILE * fdin = fopen(argv[1], "r");
	int nadd = atoi(argv[2]);
	printf("nadd = %d\n",nadd);
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
	real erroradd[numleaves];
	real diag[numleaves];
	real diagadd[numleaves];
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
	Point ptsadd[maxpts];
	real w[maxpts];
	uniqueid gids[maxpts];


        wls_interpolator *wls = wls_create(2, maxpts);
	int nbadpts = 0;
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
		higcit_destroy(itn);
                int numptsaux = wls->numpoly+nadd;
		int numptsadd = numpts-numptsaux;
                for (int i = 0; i < numptsadd; i++) {
			for (int j = 0; j < DIM; j++) {
				ptsadd[i][j] = pts[i+numptsaux][j];
			}
		}
		diag[cgid] = wls_set_points_and_calc(wls, numptsaux, pts, ccenter, w);
		real fval = 0.0;
		for(int i = 0; i < numptsaux; i++) {
			fval += w[i] * val[gids[i]];
		}
		error[cgid] = fabs(fval - val[cgid]);
		diagadd[cgid] = wls_add_points_and_calc(wls, numptsaux, numptsadd, ptsadd, ccenter, w);
		fval = 0.0;
		for(int i = 0; i < numpts; i++) {
			fval += w[i] * val[gids[i]];
		}
		erroradd[cgid] = fabs(fval - val[cgid]);
		if (diag[cgid] < 1.0e-12) {
			nbadpts++;
			//printf("Erro = %15.10g  ErroAdd = %15.10g\n",error[cgid],erroradd[cgid]);
		}
	}
	higcit_destroy(it);
	wls_destroy(wls);

	real norm_inf = error[0];
	real diag_norm_inf = diag[0];
	int  nerro = 0;
	for(int i = 1; i < numleaves; i++) {
		if(norm_inf < error[i]) {
			norm_inf = error[i];
			diag_norm_inf = diag[i];
			nerro = i;
		}
	}
	printf("Before %d %g %g\n", nerro, norm_inf, diag_norm_inf);
	norm_inf = error[0];
	real norm_inf_add = erroradd[0];
	diag_norm_inf = diag[0];
	real diag_norm_inf_add = diagadd[0];
	nerro = 0;
	for(int i = 1; i < numleaves; i++) {
		if(norm_inf < erroradd[i]) {
			norm_inf = error[i];
			norm_inf_add = erroradd[i];
			diag_norm_inf = diag[i];
			diag_norm_inf_add = diagadd[i];
			nerro = i;
		}
	}
	printf("After %d %g %g %g %g\n", nerro, norm_inf, diag_norm_inf, norm_inf_add, diag_norm_inf_add);
	printf("Nbadpoints = %d\n",nbadpts);

	hig_destroy(root);
	return 0;
}
