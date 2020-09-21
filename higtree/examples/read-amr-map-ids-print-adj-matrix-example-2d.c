
#include<stdio.h>
#include<stdlib.h>

#include "utils.h"
#include "types.h"
#include "higtree.h"
#include "mapper.h"
#include "glib.h"
#include "higtree-io.h"

int hig_get_level(hig_cell *c) {
	int level = 0;
	hig_cell *p = c;
	while(p != NULL) {
		p = hig_get_parent(p);
		level++;
	}
	return level;
}

bool filter_by_level(hig_cell *c, void *d) {
	int *level = (int *)d;
	if (hig_get_level(c) == *level) {
		return 1;
	} else {
		return 0;
	}
}

bool ffilter(hig_cell *c, void *d) {
	double *deltas = (double *) d;
	double dmin = deltas[0];
	double dmax = deltas[1];
	Point delta;
	hig_get_delta(c, delta);
	if (FLT_LE(dmin, delta[0]) && FLT_GE(dmax, delta[1])) {
		return 1;
	}
	return 0;
}

mp_mapper *calc_good_enumeration2(hig_cell *root) {
	mp_mapper *m = mp_create();
	higcit_celliterator *it;
	int firstid = 0;
	double dmax = -10000.0;
	int nlevels = 0;
	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		Point delta;
		hig_get_delta(c, delta);
		if (FLT_LT(dmax, delta[0])) {
			dmax = delta[0];
			nlevels++;
		}
	}
	higcit_destroy(it);

	for(int i = 0; i < nlevels; i++) {
		double pars[2];
		pars[0] = dmax - EPSDELTA;
		pars[1] = dmax + EPSDELTA;
		it = higcit_create_all_leaves(root);
		higcit_celliterator *itf = higcit_create_filter(it, ffilter, pars);
		firstid = mp_assign_from_celliterator(m, itf, firstid);
		higcit_destroy(itf);
		dmax /= 2.0;
	}
	it = higcit_create_all_leaves(root);
	firstid = mp_assign_from_celliterator(m, it, firstid);
	higcit_destroy(it);
	return m;
}

mp_mapper *calc_good_enumeration3(hig_cell *root) {
	mp_mapper *m = mp_create();
	higcit_celliterator *it;
	int firstid = 0;

	for(int level = 15; level >= 1; level--) {
		it = higcit_create_all_leaves(root);
		higcit_celliterator *itf = higcit_create_filter(it, filter_by_level, &level);
		firstid = mp_assign_from_celliterator(m, itf, firstid);
		higcit_destroy(itf);
	}
	it = higcit_create_all_leaves(root);
	firstid = mp_assign_from_celliterator(m, it, firstid);
	higcit_destroy(it);
	return m;
}

mp_mapper *calc_good_enumeration(hig_cell *root) {
	mp_mapper *m = mp_create();
	higcit_celliterator *it;
	int firstid = 0;
	Point p;
	hig_cell *c;
	uniqueid cid;

	POINT_ASSIGN_REALS(p, 0.85, 0.75);
	c = hig_get_cell_with_point(root, p);
	cid = hig_get_cid(c);
	mp_assign(m, cid, firstid);
	firstid++;

	POINT_ASSIGN_REALS(p, 0.85, 0.65);
	c = hig_get_cell_with_point(root, p);
	cid = hig_get_cid(c);
	mp_assign(m, cid, firstid);
	firstid++;

	Point l, h;
	POINT_ASSIGN_REALS(l, 0.55, 0.45);
	POINT_ASSIGN_REALS(h, 0.85, 0.75);
	it = higcit_create_bounding_box(root, l, h);
	firstid = mp_assign_from_celliterator(m, it, firstid);
	higcit_destroy(it);

	POINT_ASSIGN_REALS(l, 0.15, 0.25);
	POINT_ASSIGN_REALS(h, 0.65, 0.55);
	it = higcit_create_bounding_box(root, l, h);
	firstid = mp_assign_from_celliterator(m, it, firstid);
	higcit_destroy(it);

	it = higcit_create_all_leaves(root);
	firstid = mp_assign_from_celliterator(m, it, firstid);
	higcit_destroy(it);
	return m;
}

int main(int argc, char *argv[]) {
	if (argc <= 1) {
		printf("usage: %s amr2d.file\n", argv[0]);
		exit(0);
	}
	FILE * fdin = fopen(argv[1], "r");
	hig_cell * root = higio_read_from_amr2d(fdin);
	fclose(fdin);

	mp_mapper *m = calc_good_enumeration3(root);

	int numcells = hig_get_number_of_leaves(root);
	int mat[numcells][numcells];
	for (int i = 0; i < numcells; i++) {
		for (int j = 0; j < numcells; j++) {
			mat[i][j] = 0;
		}
	}
	higcit_celliterator *it;
	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		uniqueid cid = hig_get_cid(c);
		int cgid = mp_lookup(m, cid);
		higcit_celliterator *itn;
		for(itn = higcit_create_neighbours(c); !higcit_isfinished(itn); higcit_nextcell(itn)) {
			hig_cell *n = higcit_getcell(itn);
			uniqueid nid = hig_get_cid(n);
			int ngid = mp_lookup(m, nid);
			mat[cgid][ngid] = ngid;
		}
		higcit_destroy(itn);
	}
	higcit_destroy(it);
	mp_destroy(m);

	for (int i = 0; i < numcells; i++) {
//		printf("%d: ", i);
		for (int j = 0; j < numcells; j++) {
//			if (mat[i][j] != 0) {
				printf("%d ", mat[i][j]);
//			}
		}
		printf("\n");
	}

	hig_destroy(root);
	return 0;
}
