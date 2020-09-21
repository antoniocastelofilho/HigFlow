
#include<stdio.h>
#include<stdlib.h>

#include "mtree.h"
#include "mtree-io.h"

int main(int argc, char *argv[]) {
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	mt_cell *root = mt_create_root(lp, hp);
	int num_cell[DIM];
	POINT_ASSIGN_SCALAR(num_cell, 2);
	mt_refine_uniform(root, num_cell);

	FILE *fd = fopen("tut01.vtk", "w");
	mtio_print_in_vtk2d(fd, root);
	fclose(fd);

	mt_destroy(root);
}
