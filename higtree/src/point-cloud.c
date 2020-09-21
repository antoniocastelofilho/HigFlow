#include <math.h>
#include <stdlib.h>
#include "Debug-c.h"
#include "utils.h"
#include "point-cloud.h"

point_cloud *pc_create() {
	DECL_AND_ALLOC(point_cloud, pc, 1);
	pc->numpts = 0;
	int maxpts = 20;
	pc->maxpts = maxpts;
	DECL_AND_ALLOC(Point, pts, maxpts);
	pc->pts = pts;
	return pc;
}

void pc_reset(point_cloud *pc) {
	pc->numpts = 0;
}

void pc_add_point(point_cloud *pc, Point pt) {
	if (pc->numpts >= pc->maxpts) {
		int maxpts = pc->maxpts * 2;
		DECL_AND_ALLOC(Point, pts, maxpts);
		for(int i = 0; i < pc->maxpts; i++) {
			POINT_ASSIGN(pts[i], pc->pts[i]);
		}
		free(pc->pts);
		pc->pts = pts;
		pc->maxpts = maxpts;
	}
	POINT_ASSIGN(pc->pts[pc->numpts], pt);
	pc->numpts++;
}

void pc_add_points(point_cloud *pc, int numpts, Point *pts) {
	for(int i = 0; i < numpts; i++) {
		pc_add_point(pc, pts[i]);
	}
}

int pc_get_numpts(point_cloud *pc) {
	return pc->numpts;
}

void pc_copy_point(point_cloud *pc, int numpt, Point pt) {
	assert(numpt < pc->numpts);
	POINT_ASSIGN(pt, pc->pts[numpt]);
}

PPoint pc_get_point(point_cloud *pc, int numpt) {
	assert(numpt < pc->numpts);
	return pc->pts[numpt];
}

void pc_calc_diff(point_cloud *pcin, Point x, point_cloud *pcout) {
	pc_reset(pcout);
	for(int i = 0; i < pcin->numpts; i++) {
		Point diff;
		POINT_SUB(diff, pcin->pts[i], x);
		pc_add_point(pcout, diff);
	}
}

int pc_hash_key(point_cloud *pc) {
	real val = 0.0;
	static real factor = 0.70710678118;
	for(int i = 0; i < pc->numpts; i++) {
		for(int j = 0; j < DIM; j++) {
			val = val * factor + pc->pts[i][j];
		}
	}
	return (int) (sqrt(val) / EPSDELTA);
}

int pc_equal(point_cloud *pc1, point_cloud *pc2) {
	if (pc1->numpts != pc2->numpts) {
		return 0;
	} else {
		for(int i = 0; i < pc1->numpts; i++) {
			for(int j = 0; j < DIM; j++) {
				if (FLT_NE(pc1->pts[i][j], pc2->pts[i][j])) {
					return 0;
				}
			}
		}
	}
	return 1;
}

void pc_destroy(point_cloud *pc) {
	free(pc->pts);
	free(pc);
}

