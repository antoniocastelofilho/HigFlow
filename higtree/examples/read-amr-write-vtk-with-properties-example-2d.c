
#include<stdio.h>
#include<stdlib.h>

#include "higtree.h"
#include "higtree-io.h"

real f(real x, real y) {
	return cos(x)*sin(y);
}

void write_property_cell(FILE *fd, hig_cell *root) {
	higcit_celliterator *it;
	it = higcit_create_all_leaves(root);
	int numcells = higcit_count(it);
	higcit_destroy(it);
	fprintf(fd, "CELL_DATA %d\nSCALARS u float 1\nLOOKUP_TABLE default\n", numcells);
	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		Point center;
		hig_get_center(c, center);
		real value = f(center[0], center[1]);
		fprintf(fd, "%f\n", value);
	}
	higcit_destroy(it);
}

int main(int argc, char *argv[]) {
	if (argc <= 2) {
		printf("usage: %s amr2d.file vtk.file\n", argv[0]);
		exit(-1);
	}
	FILE * fdin = fopen(argv[1], "r");
	FILE * fdout = fopen(argv[2], "w");
	hig_cell * root = higio_read_from_amr2d(fdin);
	fclose(fdin);
	higio_print_in_vtk2d(fdout, root);
	write_property_cell(fdout, root);
	fclose(fdout);
	hig_destroy(root);
	return 0;
}
