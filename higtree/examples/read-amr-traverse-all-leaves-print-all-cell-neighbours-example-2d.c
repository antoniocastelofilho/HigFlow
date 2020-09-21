
#include<stdio.h>
#include<stdlib.h>

#include "higtree.h"
#include "higtree-io.h"

int main(int argc, char *argv[]) {
	if (argc <= 1) {
		printf("usage: %s amr2d.file\n", argv[0]);
		exit(0);
	}
	FILE * fdin = fopen(argv[1], "r");
	hig_cell * root = higio_read_from_amr2d(fdin);
	fclose(fdin);

	higcit_celliterator *it;
	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		printf("neighbours of cell %d:", hig_get_cid(c));
		higcit_celliterator *itn;
		for(itn = higcit_create_neighbours(c); !higcit_isfinished(itn); higcit_nextcell(itn)) {
			hig_cell *n = higcit_getcell(itn);
			printf(" %d", hig_get_cid(n));
		}
		printf("\n");
		higcit_destroy(itn);
	}
	higcit_destroy(it);

	hig_destroy(root);
	return 0;
}
