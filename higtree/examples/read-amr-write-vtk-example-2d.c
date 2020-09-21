
#include<stdio.h>
#include<stdlib.h>

#include "higtree.h"
#include "higtree-io.h"

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
	fclose(fdout);
	hig_destroy(root);
	return 0;
}
