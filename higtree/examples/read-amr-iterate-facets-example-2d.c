
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

int main(int argc, char *argv[]) {
	FILE *fd = fopen(argv[1], "r");
	hig_cell *root = higio_read_from_amr2d(fd);
	fclose(fd);
	fd = fopen("/tmp/t.vtk", "w");
	higio_print_in_vtk2d(fd, root);
	fclose(fd);
	mp_mapper *m[DIM];
	for(int dim = 0; dim < DIM; dim++) {
		m[dim] = mp_create();
	}
	DEBUG_PASS;
	int cnt[DIM];
	POINT_ASSIGN_SCALAR(cnt, 0);
	higcit_celliterator *it;
	it = higcit_create_all_leaves(root);
	int dimofinterest[DIM];
	POINT_ASSIGN_SCALAR(dimofinterest, 1);
	higfit_facetiterator *fit;
	DEBUG_PASS;
	for(fit = higfit_create_allfacets(it, dimofinterest);!higfit_isfinished(fit); higfit_nextfacet(fit)) {
		hig_facet *f = higfit_getfacet(fit);
		hig_cell *c = hig_get_facet_cell(f);
		int dim = hig_get_facet_dim(f);
		int dir = hig_get_facet_dir(f);
		int cid = hig_get_fid_with_dim_dir(c, dim, dir);
		mp_assign(m[dim], cid, cnt[dim]++);
	}
	DEBUG_PASS;
	higfit_destroy(fit);
	DEBUG_PASS;
	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		for(int dim = 0; dim < DIM; dim++) {
			printf(" DIM = %d", dim);
			for(int dir = 0; dir < 2; dir++) {
				int cid = hig_get_fid_with_dim_dir(c, dim, dir);
				int gid = -1;
				if (cid > 0) {
					gid = mp_lookup(m[dim], cid);
				} else if (cid < 0) {
					gid = mp_lookup(m[dim], -cid);
				}
				printf(" %5d(%5d)", gid, cid);
			}
		}
		printf("\n");
	}
	higcit_destroy(it);
	printf("***********\n");
	Point l, h;
	POINT_ASSIGN_SCALAR(l, 0.1);
	POINT_ASSIGN_SCALAR(h, 0.9);
	it = higcit_create_bounding_box(root, l, h);
	POINT_ASSIGN_SCALAR(dimofinterest, 0);
	for(fit = higfit_create_allfacets(it, dimofinterest); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
		hig_facet *f = higfit_getfacet(fit);
		Point center;
		hig_get_facet_center(f, center);
		printf("(%f, %f)\n", center[0], center[1]);
	}
	higfit_destroy(fit);
	printf("***********\n");
	it = higcit_create_bounding_box(root, l, h);
	POINT_ASSIGN_SCALAR(dimofinterest, 0);
	int dim = 1;
	dimofinterest[dim] = 1;
	for(fit = higfit_create_allfacets(it, dimofinterest);!higfit_isfinished(fit); higfit_nextfacet(fit)) {
		hig_facet *f = higfit_getfacet(fit);
		Point fl, fh;
		hig_get_facet_lowpoint(f, fl);
		hig_get_facet_highpoint(f, fh);
		printf("(%f, %f)(%f, %f)\n", fl[0], fl[1], fh[0], fh[1]);
		if (hig_facet_intersect_bounding_box(f, l, h)) {
			Point center;
			hig_get_facet_center(f, center);
			printf("[%f, %f]@\n", center[0], center[1]);
			hig_facet facet;
			int hasf = hig_get_facet_with_point(root, dim, center, &facet);
			int badfacet = 0;
			if (!hasf) {
				badfacet = 1;
			} else {
				Point center2;
				hig_get_facet_center(&facet, center2);
				if (!Pointequal(center, center2)) {
					badfacet = 1;
				}
			}
			if (badfacet) {
				printf("!!!!!!!!!!!!\n");
				exit(0);
			}
		}
	}
	higfit_destroy(fit);
	dimofinterest[0] = 1;
	POINT_ASSIGN_REALS(l, 0.05, 0.15);
	POINT_ASSIGN_REALS(h, 0.05, 0.46);
	for(fit = higfit_create_bounding_box_facets(root, dimofinterest, l, h); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
		hig_facet *f = higfit_getfacet(fit);
		Point center;
		hig_get_facet_center(f, center);
		printf("[%f, %f]*\n", center[0], center[1]);
	}
	higfit_destroy(fit);
	for(it = higcit_create_bounding_box(root, l, h); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		for(int dim = 0; dim < DIM; dim++) {
			printf(" DIM = %d", dim);
			for(int dir = 0; dir < 2; dir++) {
				int cid = hig_get_fid_with_dim_dir(c, dim, dir);
				int gid = -1;
				if (cid > 0) {
					gid = mp_lookup(m[dim], cid);
				} else if (cid < 0) {
					gid = mp_lookup(m[dim], -cid);
				}
				printf(" %5d(%5d)", gid, cid);
			}
		}
		printf("\n");
	}
	higcit_destroy(it);
	int numpts = 8;
	Point pts[numpts];
	POINT_ASSIGN_REALS(pts[0], 0.1, 0.05);
	POINT_ASSIGN_REALS(pts[1], 0.1, 0.1);
	POINT_ASSIGN_REALS(pts[2], 0.12, 0.12);
	POINT_ASSIGN_REALS(pts[3], 0.2, 0.12);
	POINT_ASSIGN_REALS(pts[4], 0.2, 0.125);
	POINT_ASSIGN_REALS(pts[5], 0.0, 0.125);
	POINT_ASSIGN_REALS(pts[6], 0.9, 0.125);
	POINT_ASSIGN_REALS(pts[7], 0.3, 0.225);
	for(int i = 0; i < numpts; i++) {
		hig_facet facet;
		int hasf = hig_get_facet_with_point(root, 0, pts[i], &facet);
		printf("[%f, %f]#\n", pts[i][0], pts[i][1]);
		if (hasf) {
			Point fl, fh;
			hig_get_facet_lowpoint(&facet, fl);
			hig_get_facet_highpoint(&facet, fh);
			int dir = hig_get_facet_dir(&facet);
			printf("(%f, %f)(%f, %f) %d\n", fl[0], fl[1], fh[0], fh[1], dir);
		} else {
			printf("no facet\n");
		}
	}
	printf("sizeof(hig_cell) = %ld\n", sizeof(hig_cell));
	return 0;
}
