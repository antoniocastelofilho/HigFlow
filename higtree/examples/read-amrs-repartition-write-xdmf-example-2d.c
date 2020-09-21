#include <string.h>
#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<assert.h>

#include "higtree.h"
#include "higtree-io.h"
#include "lbal.h"

bool manage = true;

int main(int argc, char *argv[]) {
	if (argc < 3) {
		printf("usage: %s file1.amr [file2.amr ...] [-b bc_file1.amr [bc_file2.amr...]] xdmf_prefix\n", argv[0]);
		exit(0);
	}
	const char *xdmf_prefix = argv[argc-1];
	size_t xdmf_prefix_len = strlen(xdmf_prefix);

	// Parse the command line parameters
	DECL_AND_ALLOC(char *, domain_trees, argc-2);
	size_t domain_trees_count = 0;

	DECL_AND_ALLOC(char *, bc_trees, argc-2);
	size_t bc_trees_count = 0;
	{
		int i = 1;
		for(; i < argc-1; ++i) {
			if(strcmp("-b", argv[i]) == 0)
				break;
			domain_trees[domain_trees_count++] = argv[i];
		}

		for(++i; i < argc-1; ++i) {
			bc_trees[bc_trees_count++] = argv[i];
		}
	}

	// MPI stuff
	int myrank, ntasks;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

	// Simulation domain
	sim_domain *sd = sd_create(NULL);

	// Distributed partition context
	load_balancer *part_ctx = lb_create(MPI_COMM_WORLD, 2);
	lb_set_group_weight(part_ctx, 1, 0);

	// Randomly shift the start of the Round-Robin among processes
	int rand_rank;
	if(myrank == 0) {
		puts("Scrambling processes ids");
		int ranks[ntasks];
		// Using /dev/urandom because rand seems biased
		FILE *frand = fopen("/dev/urandom", "r");

		ranks[0] = 0;
		for(int i = 1; i < ntasks; ++i) {
			unsigned j;
			int ret = fread(&j, sizeof j, 1, frand);
			assert(ret == 1);
			j %= (i+1);

			ranks[i] = ranks[j];
			ranks[j] = i;
		}
		fclose(frand);
		rand_rank = ranks[0];
		for(int i = 1; i < ntasks; ++i) {
			MPI_Send(&ranks[i], 1, MPI_INT, i, 33, MPI_COMM_WORLD);
		}
	} else {
		MPI_Recv(&rand_rank, 1, MPI_INT, 0, 33, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank == 0)
		puts("Loading trees from AMR");
	size_t max_trees = domain_trees_count > bc_trees_count ? domain_trees_count : bc_trees_count;
	for(unsigned i = rand_rank; i < max_trees; i += ntasks) {
		if(i < domain_trees_count) {
			FILE * fdin = fopen(domain_trees[i], "r");
			hig_cell *tree = higio_read_from_amr2d(fdin);
			fclose(fdin);

			// Add original trees to simulation domain for printing
			// If original tree will be managed by the partition,
			// we must copy it for the domain.
			sd_add_higtree(sd, manage ? hig_clone(tree) : tree);

			// Set the tree for partitioning at domain group
			lb_add_input_tree(part_ctx, tree, manage, 0);
		}
		if(i < bc_trees_count) {
			FILE * fdin = fopen(bc_trees[i], "r");
			hig_cell *tree = higio_read_from_amr2d(fdin);
			fclose(fdin);

			sd_add_higtree(sd, manage ? hig_clone(tree) : tree);

			// Set the tree for partitioning at boundary condition group
			lb_add_input_tree(part_ctx, tree, manage, 1);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank == 0)
		puts("Writing out original partition");
	{
		char prefix[xdmf_prefix_len+10];
		snprintf(prefix, sizeof prefix, "%s_original", xdmf_prefix);
		xdmf_output *output_ctx = xdmf_init(prefix, sd);
		xdmf_write_timestep(output_ctx, 0, false);
		xdmf_destroy(output_ctx);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank == 0)
		puts("Partitioning the domain");
	partition_graph pg;
	lb_calc_partition(part_ctx, &pg);

	sd_destroy(sd);

	// Get the partitioned trees
	sd = sd_create(NULL);
	printf("Num local trees: %zu\n", lb_get_num_local_trees(part_ctx));
	for(int i = lb_get_num_local_trees(part_ctx) - 1; i >= 0; --i) {
		unsigned group;
		hig_cell *tree = lb_get_local_tree(part_ctx, i, &group);
		sd_add_higtree(sd, tree);
	}
	lb_destroy(part_ctx);

	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank == 0)
		puts("Writing the newly partitioned domain");
	{
		char prefix[xdmf_prefix_len+20];
		snprintf(prefix, sizeof prefix, "%s_partitioned", xdmf_prefix);
		xdmf_output *output_ctx = xdmf_init(prefix, sd);
		xdmf_write_timestep(output_ctx, 0, false);
		xdmf_destroy(output_ctx);
	}

	sd_destroy(sd);

	MPI_Finalize();
}
