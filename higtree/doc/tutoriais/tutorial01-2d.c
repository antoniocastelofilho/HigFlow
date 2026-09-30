
#include<stdio.h>
#include<stdlib.h>

#include "higtree.h"
#include "higtree-io.h"

int main (int argc, char *argv[]) {
    higtree_initialize(&argc, &argv);
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int num_cell[DIM];
	POINT_ASSIGN_SCALAR(num_cell, 2);
	hig_refine_uniform(root, num_cell);

	FILE *fd = fopen("VTKS/tut01.vtk", "w");
	higio_print_in_vtk2d(fd, root);
	fclose(fd);

	hig_destroy(root);
}
