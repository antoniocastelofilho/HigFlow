
#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>

#include "higtree.h"
#include "higtree-io.h"
#include "lbal.h"

int main(int argc, char *argv[]) {
	if (argc <= 2) {
		printf("usage: %s amr3d.file xdmf_prefix\n", argv[0]);
		exit(0);
	}

	// MPI stuff
	int myrank, ntasks;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

	// Load AMR
	FILE * fdin = fopen(argv[1], "r");
	higio_amr_info *mi = higio_read_amr_info(fdin);
	fclose(fdin);

	// Partition and create tree
	partition_graph *pg = pg_create(MPI_COMM_WORLD);
	load_balancer *lb = higio_partition_from_amr_info(pg, mi);

	// Use local tree in sim-domain
	sim_domain *sd = sd_create(NULL);
	{
		unsigned numtrees = lb_get_num_local_trees(lb);
		for(unsigned i = 0; i < numtrees; ++i) {
			hig_cell *root = lb_get_local_tree(lb, i, NULL);
			sd_add_higtree(sd, root);
		}
	}
	lb_destroy(lb);

	// Writing stuff
	xdmf_output *output_ctx = xdmf_init(argv[2], sd);
	xdmf_write_timestep(output_ctx, 0, false);
	xdmf_destroy(output_ctx);

	sd_destroy(sd);

	MPI_Finalize();
}
