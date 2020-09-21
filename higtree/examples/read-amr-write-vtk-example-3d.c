
#include<stdio.h>
#include<stdlib.h>

#include "higtree.h"
#include "higtree-io.h"

int main(int argc, char *argv[]) {
	if (argc <= 2) {
		printf("usage: %s amr3d.file vtk.file\n", argv[0]);
		exit(0);
	}
	FILE * fdin = fopen(argv[1], "r");
	FILE * fdout = fopen(argv[2], "w");
	hig_cell * root = higio_read_from_amr3d(fdin);
	fclose(fdin);
        printf("Nodes %ld\n", hig_get_number_of_leaves(root));
	higio_print_in_vtk3d(fdout, root);
	fclose(fdout);
	hig_destroy(root);
}
