
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
	Point l, h;
	int nc[DIM];
#if DIM == 2
	POINT_ASSIGN_INTS(nc, 2, 2);
	POINT_ASSIGN_REALS(l, 0.0, 0.0);
	POINT_ASSIGN_REALS(h, 8.0, 8.0);
	Point ps[] = {{4.0, 4.0}, {2.0, 2.0}, {2.0, 6.0}, {1.0, 5.0}, {3.0, 5.0}};
	int numps = 5;
#elif DIM == 3
	POINT_ASSIGN_INTS(nc, 2, 2, 2);
	POINT_ASSIGN_REALS(l, 0.0, 0.0, 0.0);
	POINT_ASSIGN_REALS(h, 8.0, 8.0, 8.0);
	Point ps[] = {{4.0, 4.0, 4.0}, {2.0, 2.0, 2.0}};
	int numps = 2;
#endif

	hig_cell *root = hig_create_root(l, h);
	hig_cell *c;
	for (int i = 0; i < numps; i++) {
		c = hig_get_cell_with_point(root, ps[i]);
		hig_refine_uniform(c, nc);
	}
	DEBUG_PASS;

	FILE *fd = fopen("/tmp/t.vtk", "w");
#if DIM == 2
	higio_print_in_vtk2d(fd, root);
#elif DIM == 3
	higio_print_in_vtk3d(fd, root);
#endif
	fclose(fd);

	higcit_celliterator *it;
	int cnt = 0;
	int numDIM = 1;
	for(int i = 0; i < DIM; i++) {
		numDIM *= 2;
	}
	int eqCoord = atoi(argv[1]);
	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		Point cl, ch, ccenter;
		hig_get_lowpoint(c, cl);
		hig_get_highpoint(c, ch);
		hig_get_center(c, ccenter);
		for(int i = 0; i < numDIM; i++) {
			for(int j = i; j < numDIM; j++) {
				int numEq = 0;
				int numEqOne = 0;
				int valid = 1;
				int print = 0;
				for(int k = 0; k < DIM; k++) {
					int bi = i&(1<<k);
					int bj = j&(1<<k);
					if (bi > bj) {
						valid = 0;
						break;
					} else if (bi == bj) {
						numEq++;
						if (bi) {
							numEqOne++;
						}
					}
				}
				if (!valid || numEq != eqCoord) {
					continue;
				}
				Point fl, fh;
				Point fcenter;
				Point fbl;
				Point fbh;
				for(int k = 0; k < DIM; k++) {
					int bi = i&(1<<k);
					int bj = j&(1<<k);
					fl[k] = (bi?ch[k]:cl[k]);
					fh[k] = (bj?ch[k]:cl[k]);
					if(!print) {
						fcenter[k] = (fl[k] + fh[k])/2.0;
						if(bi == bj) {
							fbl[k] = fcenter[k] - EPSDELTA;
							fbh[k] = fcenter[k] + EPSDELTA;
						} else {
							fbl[k] = fl[k] + EPSDELTA;
							fbh[k] = fh[k] - EPSDELTA;
						}
					}
				}
				if (!print) {
					higcit_celliterator *it;
					print = 1;
					for (it = higcit_create_bounding_box(root, fbl, fbh); !higcit_isfinished(it) && print; higcit_nextcell(it)) {
						hig_cell *ofc = higcit_getcell(it);
						if (ofc == c) {
							continue;
						}
						Point ofcl;
						Point ofch;
						hig_get_lowpoint(ofc, ofcl);
						hig_get_highpoint(ofc, ofch);
						int isBelow = 0;
						for(int k = 0; k < DIM; k++) {
							if (FLT_EQ(ofch[k], cl[k])) {
								isBelow = 1;
								break;
							}
						}
						if (isBelow) {
							int coincide = 1;
							for(int k = 0; k < DIM && print; k++) {
								int bi = i&(1<<k);
								int bj = j&(1<<k);
								if (bi == bj) {
									if (FLT_NE(ofch[k], fcenter[k]) && FLT_NE(ofcl[k], fcenter[k])) {
										coincide = 0;
										break;
									}
								} else {
									if (FLT_NE(ofch[k], ch[k]) || FLT_NE(ofcl[k], cl[k])) {
										coincide = 0;
										break;
									}
								}
							}
							if (coincide) {
								print = 0;
							}
						}
					}
					higcit_destroy(it);
				}
				if (print) {
#if DIM == 2
					printf("%d = ((%f, %f), (%f, %f))\n", cnt++, fl[0], fl[1], fh[0], fh[1]);
#elif DIM == 3
					printf("%d = ((%f, %f, %f), (%f, %f, %f))\n", cnt++, fl[0], fl[1], fl[2], fh[0], fh[1], fh[2]);
#endif
				}
			}
		}
		printf("\n");
	}
	higcit_destroy(it);
	hig_destroy(root);
	return 0;
}
