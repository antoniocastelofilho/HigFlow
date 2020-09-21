
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
	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		Point cl, ch, ccenter;
		hig_get_lowpoint(c, cl);
		hig_get_highpoint(c, ch);
		hig_get_center(c, ccenter);
		for(int i = 0; i < DIM; i++) {
			Point fcenter;
			Point ofcenter;
			POINT_ASSIGN(fcenter, ccenter);
			POINT_ASSIGN(ofcenter, ccenter);
			for(int j = 0; j < 2; j++) {
				fcenter[i] = ((j==0)?cl[i]:ch[i]);
				ofcenter[i] = ((j==0)?cl[i]-EPSDELTA:ch[i]+EPSDELTA);
				int print = 0;
				if (j == 1) {
					print = 1;
				} else {
					hig_cell *ofc = hig_get_cell_with_point(root, ofcenter);
					if (ofc == NULL) {
						print = 1;
					} else {
						Point ofl;
						hig_get_lowpoint(ofc, ofl);
						for(int k = 0; k < DIM; k++) {
							if (k != i && FLT_NE(cl[k], ofl[k])) {
								//print = 1;
								break;
							}
						}
						if (!print) {
							Point ofh;
							hig_get_highpoint(ofc, ofh);
							for(int k = 0; k < DIM; k++) {
								if (k != i && FLT_NE(ch[k], ofh[k])) {
									//print = 1;
									break;
								}
							}
						}
					}
				}
				if (print) {
					Point fl, fh;
					POINT_ASSIGN(fl, cl);
					POINT_ASSIGN(fh, ch);
					if (j == 0) {
						fh[i] = cl[i];
					} else {
						fl[i] = ch[i];
					}
#if DIM == 2
					printf("%d = ((%f, %f), (%f, %f))\n", cnt++, fl[0], fl[1], fh[0], fh[1]);
#elif DIM == 3
					printf("%d = ((%f, %f, %f), (%f, %f, %f))\n", cnt++, fl[0], fl[1], fl[2], fh[0], fh[1], fh[2]);
#endif
				}
			}
		}
	}
	higcit_destroy(it);
	return 0;
}
