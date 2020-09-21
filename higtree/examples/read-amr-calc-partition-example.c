
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

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

void writevtk(psim_domain *psd, int myrank) {
	higcit_celliterator *it;

	sim_domain *local_domain = psd_get_local_domain(psd);
	char filename[1024];
	sprintf(filename, "/tmp/sbd-%d.vtk", myrank);
	FILE *fd = fopen(filename, "w");
	it = sd_get_domain_celliterator(local_domain);
	printf("part[%d] = %ld\n", myrank, higcit_count_without_advancing(it));
#if DIM == 3
		higio_print_celliterator_in_vtk3d(fd, it);
#elif DIM == 2
		higio_print_celliterator_in_vtk2d(fd, it);
#endif
	fclose(fd);
	higcit_destroy(it);

	sprintf(filename, "/tmp/sbfd-%d.vtk", myrank);
	fd = fopen(filename, "w");
	it = higcit_create_all_leaves(sd_get_higtree(local_domain, 0));
#if DIM == 3
		higio_print_celliterator_in_vtk3d(fd, it);
#elif DIM == 2
		higio_print_celliterator_in_vtk2d(fd, it);
#endif
	fclose(fd);
	higcit_destroy(it);

	sprintf(filename, "/tmp/bc-%d.vtk", myrank);
	fd = fopen(filename, "w");
	it = sd_get_bcs_celliterator(local_domain);
#if DIM == 3
		higio_print_celliterator_in_vtk3d(fd, it);
#elif DIM == 2
		higio_print_celliterator_in_vtk2d(fd, it);
#endif
	fclose(fd);
	higcit_destroy(it);
}

int main(int argc, char *argv[]) {
	int nc[DIM];
	int size;
	DEBUG_DIFF_TIME;
	int myrank;
	int ntasks;

	higtree_initialize(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

	partition_graph *pg = pg_create(MPI_COMM_WORLD);
	psim_domain *psd;

	load_balancer *lb = lb_create(MPI_COMM_WORLD, 1);
	lb_set_group_has_fringe(lb, 0, false);

	if (myrank == 0) {
		// Reading domain
		FILE *fd = fopen(argv[1], "r");
#if DIM == 3
		hig_cell *root = higio_read_from_amr3d(fd);
#elif DIM == 2
		hig_cell *root = higio_read_from_amr2d(fd);
#endif
		fclose(fd);
		lb_add_input_tree(lb, root, true, 0);
	}

	lb_calc_partition(lb, pg);

	sim_domain *sd = sd_create(NULL);

	{
		unsigned num_trees = lb_get_num_local_trees(lb);
		assert(num_trees <= 1);
		for(unsigned i = 0; i < num_trees; ++i) {
			hig_cell *root = lb_get_local_tree(lb, i, NULL);
			sd_add_higtree(sd, root);

			Point lg;
			Point hg;
			hig_get_lowpoint(root, lg);
			hig_get_highpoint(root, hg);

			char out[100];
			int ccount = snprintf(out, 100, "%g %g | %g %g\n", lg[0], hg[0], lg[1], hg[1]);
#if DIM == 3
			--ccount;
			ccount += snprintf(&out[ccount], 100 - ccount, " | %g %g\n", lg[2], hg[2]);
#endif
			// Use write because it is as close to an atomic operation as we can get...
			write(STDOUT_FILENO, out, ccount);
		}
	}

	unsigned num_elems;
	{
		higcit_celliterator *it;
		it = sd_get_domain_celliterator(sd);
		num_elems = higcit_count(it);
		higcit_destroy(it);
	}

	if(myrank == 0) {
		unsigned int total;
		unsigned int maximum;
		MPI_Reduce(&num_elems, &total, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&num_elems, &maximum, 1, MPI_UNSIGNED, MPI_MAX, 0, MPI_COMM_WORLD);
		printf("Partition imbalance: %g\n", ((double)total) / (ntasks * maximum));
	} else {
		MPI_Reduce(&num_elems, NULL, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&num_elems, NULL, 1, MPI_UNSIGNED, MPI_MAX, 0, MPI_COMM_WORLD);
 	}

	psd = psd_create(sd, pg);

	DEBUG_DIFF_TIME;
	psd_synced_mapper(psd);

	writevtk(psd, myrank);

	return 0;
}

