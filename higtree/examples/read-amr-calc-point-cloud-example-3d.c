
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

Point allpts[30000][100];
unsigned long hashcods[30000];
int numallpts[30000];
int numhcs = 0;
int numconflictsperhc[30000];
int maxconflictperhc = 0;
int numchecks;

int main(int argc, char *argv[]) {
	if (argc <= 1) {
		printf("usage: %s amr3d.file\n", argv[0]);
		exit(-1);
	}
	FILE * fdin = fopen(argv[1], "r");
	hig_cell * root = higio_read_from_amr3d(fdin);
	fclose(fdin);
	const int maxpts = 100;

	const int numirs = 301;
	real irrationals[numirs];
	for(int i = 0; i < numirs; i++) {
		irrationals[i] = sqrt(i+1.0);
	}
	int numconflicts = 0;
	int numcached = 0;

	higcit_celliterator *it;
	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		Point ccenter;
		hig_get_center(c, ccenter);
		Point diff;
		Point pts[maxpts];
		higcit_celliterator *itn;
		int numpts = 0;
		real val = 0.0;
		real ir = 1.0;
		real sq = 1.0;
		real scale = -1.0;
		int posir = 0;
		for(itn = higcit_create_neighbours(c); !higcit_isfinished(itn) && numpts < maxpts; higcit_nextcell(itn), numpts++) {
			hig_cell *n = higcit_getcell(itn);
			Point ncenter;
			hig_get_center(n, ncenter);
			POINT_SUB(diff, ncenter, ccenter);
			if (scale < 0) {
				for (int i = 0; i < DIM; i++) {
					//DEBUG_INSPECT(diff[i], %lg);
					scale = fabs(diff[i]);
					if (scale > EPSDELTA) {
						break;
					}
				}
			}
			//POINT_DIV_SCALAR(diff, diff, scale);
			POINT_ASSIGN(pts[numpts], diff);
			for(int i = 0; i < DIM; i++) {
				val += diff[i] * sqrt(i+5.0) * irrationals[posir];
		//		printf("%0.10g ", diff[i]);
				//ir = (ir + 1) % numirs;
				ir = (ir + sq);
			}
			posir = (posir + 1) % numirs;
		}
		higcit_destroy(itn);
		//printf("\n");
		//int hc = (int)((val + (0.5*EPSDELTA))/EPSDELTA);
		unsigned long hc = (((*(unsigned long *)&val)<<11)>>11);
		//int hc = bit;
		//printf("%d\n", hc);
		int iscached = 0;
		for(int j = 0; j < numhcs; j++) {
			if (hc == hashcods[j]) {
				int isequal = 1;
				if (numpts != numallpts[j]) {
					isequal = 0;
				}
				for(int k = 0; isequal && k < numpts; k++) {
					for(int l = 0; isequal && l < DIM; l++) {
						if(FLT_NE(allpts[j][k][l], pts[k][l])) {
							isequal = 0;
						}
					}
				}
				if (!isequal) {
					numconflictsperhc[j]++;
					if(maxconflictperhc < numconflictsperhc[j]) {
						maxconflictperhc = numconflictsperhc[j];
					}
					numchecks += numconflictsperhc[j];
					numconflicts++;
					for(int k = 0; k < numpts; k++) {
						printf("(");
						for(int l = 0; l < DIM; l++) {
							printf("[%0.10g %0.10g] ", allpts[j][k][l], pts[k][l]);
						}
						printf(") ");
					}
					printf("\n");
				}
				iscached = 1;
				break;
			}
		}
		if (iscached) {
			numcached++;
		} else {
			hashcods[numhcs] = hc;
			numallpts[numhcs] = numpts;
			for(int k = 0; k < numpts; k++) {
				for(int l = 0; l < DIM; l++) {
					allpts[numhcs][k][l] = pts[k][l];
				}
			}
			numconflictsperhc[numhcs] = 1;
			numhcs++;
		}
		printf("%d %d %d %d %d\n", numhcs, numcached, numconflicts, maxconflictperhc, numchecks);
	}
	higcit_destroy(it);


	hig_destroy(root);
	return 0;
}
