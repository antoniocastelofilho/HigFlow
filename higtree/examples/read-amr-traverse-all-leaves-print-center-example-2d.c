
#include<stdio.h>
#include<stdlib.h>

#include "higtree.h"
#include "higtree-io.h"

int main(int argc, char *argv[]) {
	if (argc <= 1) {
		printf("usage: %s amr2d.file\n", argv[0]);
		exit(-1);
	}
	FILE * fdin = fopen(argv[1], "r");
	hig_cell * root = higio_read_from_amr2d(fdin);
	fclose(fdin);

	higcit_celliterator *it;
	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		Point center;
		hig_get_center(c, center);
		printf("center of cell %d: (%lg, %lg)\n", hig_get_cid(c), center[0], center[1]);
	}
	higcit_destroy(it);

	hig_destroy(root);
	return 0;
}
