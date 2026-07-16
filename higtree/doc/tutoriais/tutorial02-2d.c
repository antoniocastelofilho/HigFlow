
#include<stdio.h>
#include<stdlib.h>

#include "higtree.h"
#include "higtree-io.h"

int main(int argc, char *argv[]) {
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int num_cells[DIM];
        num_cells[0] = 4;
        num_cells[1] = 5;
	hig_refine_uniform(root, num_cells);

        Point p;
        p[0] = 0.5;
        p[1] = 0.5;
        hig_cell *c = hig_get_cell_with_point(root,p);
        num_cells[0] = 2;
        num_cells[1] = 2;
        hig_refine_uniform(c,num_cells);

	FILE *fd = fopen("VTKS/tut02.vtk", "w");
	higio_print_in_vtk2d(fd, root);
	fclose(fd);

	hig_destroy(root);
}
