
#include<stdio.h>
#include<stdlib.h>

#include "mtree.h"
#include "mtree-io.h"

int main(int argc, char *argv[]) {
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	mt_cell *root = mt_create_root(lp, hp);
	int num_cells[DIM];
        num_cells[0] = 4;
        num_cells[1] = 5;
	mt_refine_uniform(root, num_cells);

        Point p;
        p[0] = 0.5;
        p[1] = 0.5;
        mt_cell *c = mt_get_cell_with_point(root,p);
        num_cells[0] = 2;
        num_cells[1] = 2;
        mt_refine_uniform(c,num_cells);

	FILE *fd = fopen("tut02.vtk", "w");
	mtio_print_in_vtk2d(fd, root);
	fclose(fd);

	mt_destroy(root);
}
