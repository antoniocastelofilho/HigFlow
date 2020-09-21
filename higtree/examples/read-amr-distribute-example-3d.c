
#include<stdio.h>
#include<stdlib.h>

#include "utils.h"
#include "higtree.h"
#include "higtree-io.h"
#include "mapper.h"
#include "wls.h"
#include "glib.h"
#include "pdomain.h"
#include "solver.h"
#include "lbal.h"

#define DEBUG
#include "Debug-c.h"
#define PI 3.1415926535897932385


int main(int argc, char *argv[]) {
	int myrank;
	int ntasks;
	higtree_initialize(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

	partition_graph *pg = pg_create(MPI_COMM_WORLD);
	psim_domain *psd;

	DEBUG_DIFF_TIME;
	FILE *fd = fopen(argv[1], "r");
	higio_amr_info *mi = higio_read_amr_info(fd);
	fclose(fd);
	load_balancer *lb = higio_partition_from_amr_info(pg, mi);
	unsigned numtrees = lb_get_num_local_trees(lb);
	hig_cell *roots[numtrees];
	for(unsigned i = 0; i < numtrees; ++i) {
		roots[i] = lb_get_local_tree(lb, i, NULL);
	}
	lb_destroy(lb);
	DEBUG_DIFF_TIME;
	return 0;
}

