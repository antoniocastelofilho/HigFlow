/*********************************************************************************************/
/***    Including files: Begin                                                             ***/
/*********************************************************************************************/

#include "higtree.h"
#include "higtree-io.h"
#include "higtree-iterator.h"
#include "lbal.h"
#define DEBUG
#include "Debug-c.h"
#include <stdbool.h>
#include <assert.h>
#include <string.h>
#include <mpi.h>

/*********************************************************************************************/
/***    Including files: End                                                               ***/
/*********************************************************************************************/

/* Test if BFS yelds the same number of nodes as the other iterator implemented. */
static void test_bfs(hig_cell *root)
{
	higcit_celliterator *it;

	size_t all_count = 0;
	DEBUG_DIFF_TIME_START;
	for (it = higcit_create_all_higtree(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		++all_count;
	}
	higcit_destroy(it);
	DEBUG_DIFF_TIME_FINISH("All iteration duration");

	printf("All higtree iterator count: %lu\n", all_count);

	size_t bfs_count = 0;
	DEBUG_DIFF_TIME_START;
	for (it = higcit_create_breadth_first(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		++bfs_count;
	}
	higcit_destroy(it);
	DEBUG_DIFF_TIME_FINISH("BFS iteration duration");

	printf("BFS iterator count: %lu\n", bfs_count);
	assert(all_count == bfs_count);
}

static void compare_trees(hig_cell *a, hig_cell *b)
{
	higcit_celliterator *orig_it = higcit_create_all_higtree(a);
	higcit_celliterator *rest_it = higcit_create_all_higtree(b);

	while(!higcit_isfinished(rest_it)) {
		assert(!higcit_isfinished(orig_it));

		hig_cell *orig = higcit_getcell(orig_it);
		hig_cell *rest = higcit_getcell(rest_it);

		assert(hig_get_number_of_children(orig) == hig_get_number_of_children(rest));
		Point a, b;

		hig_get_center(orig, a);
		hig_get_center(rest, b);
		assert(memcmp(a, b, sizeof a) == 0);

		hig_get_lowpoint(orig, a);
		hig_get_lowpoint(rest, b);
		assert(memcmp(a, b, sizeof a) == 0);

		hig_get_highpoint(orig, a);
		hig_get_highpoint(rest, b);
		assert(memcmp(a, b, sizeof a) == 0);

		higcit_nextcell(orig_it);
		higcit_nextcell(rest_it);
	};
	assert(higcit_isfinished(orig_it));

	higcit_destroy(orig_it);
	higcit_destroy(rest_it);
}

static void finalize()
{
	MPI_Finalize();
}

int main (int argc, char *argv[]) {
        /* Initializing MPI */
	int myrank;
	int ntasks;

	MPI_Init(&argc, &argv);
	atexit(finalize);
	int mpierr = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

        /* Seeing the success or not of initializing */
	DEBUG_INSPECT(mpierr, %d);
	DEBUG_INSPECT(MPI_SUCCESS, %d);
	DEBUG_INSPECT(myrank, %d);

        /* Partitioning the domain */
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
	DEBUG_INSPECT(myrank, %d);

        /* Taking the sub-domain of the myrank process */
	partition_graph *pg = pg_create(MPI_COMM_WORLD);

        /* Setting the fringe size of the sub-domain */
	pg_set_fringe_size(pg, 4);

        /* Getting the time */
	DEBUG_DIFF_TIME;

	if(argc < 2) {
		fputs("Missing parameter: amr_file\n", stderr);
		return 1;
	}

        /* Adapting the data input by the command line */
	char *amrfilename = argv[1];

	printf("AMR file: %s\n", amrfilename);

        /* Open the AMR format file */
	FILE *fd = fopen(amrfilename, "r");

        /* Reanding the IO information of the file */
	higio_amr_info *mi = higio_read_amr_info(fd);

        /* Closing the AMR format file */
	fclose(fd);

	/* Partitioning the grid from AMR information */
	load_balancer *lb = higio_partition_from_amr_info(pg, mi);

	/* We don't need AMR info anymore... */
	higio_amr_info_destroy(mi);

	/* Creating the simulation domain */
	sim_domain *sdp = sd_create(NULL);

        /* Setting the order of the polynomial interpolation */
	sd_set_interpolator_order(sdp, 2);

	/* Get the number of trees in the local process. */
	unsigned numtrees = lb_get_num_local_trees(lb);

        /* Linking the HigTree grid to the simulation domain */
	hig_cell *roots[numtrees];
	for(unsigned i = 0; i < numtrees; ++i) {
		roots[i] = lb_get_local_tree(lb, i, NULL);
		sd_add_higtree(sdp, roots[i]);
		test_bfs(roots[i]);
		hig_adjust_facet_ids(roots[i]);
	}

	/* Destroy partition context.*/
	lb_destroy(lb);

        /* Creating the partitioned sub-domain */
	psim_domain *psdp = psd_create(sdp, pg);

	/* Creating the distributed properties */
	psd_synced_mapper(psdp);
	distributed_property * dp1 = psd_create_property(psdp);
	distributed_property * dp2 = psd_create_property(psdp);

	/* Fill the cell property with random values */
	puts("Filling cell properties with random values...");
	fflush(stdout);
	{
		DEBUG_DIFF_TIME_START;
		srand48(time(NULL));

		mp_mapper *m = sd_get_domain_mapper(sdp);
		higcit_celliterator *it;
		for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
			hig_cell *cell = higcit_getcell(it);
			dp_set_value(dp1, mp_lookup(m, hig_get_cid(cell)), mrand48() * drand48());
			dp_set_value(dp2, mp_lookup(m, hig_get_cid(cell)), mrand48() * drand48());
		}
		higcit_destroy(it);
		DEBUG_DIFF_TIME_FINISH("  ...duration");
	}

	/* Creating the vector properties from facets */

        /* Facets list */
	sim_facet_domain *sfdu[DIM];

        /* Property facet list */
	psim_facet_domain *psfdu[DIM];

        /* Distributed properties from facets in the domain.
	 * The property in the first direction U will have 2 components.
	 * The others V and W will have one component each. */

        /* Bidimensional facet property on one direction. */
	distributed_property *u[2];

	/* Scalar facect property of other directions. */
	distributed_property *vw[DIM-1];

        /* Linking the property to the facets */
	for(int dim = 0; dim < DIM; dim++) {
                /* Creating the facet domain */
		sfdu[dim] = sfd_create(NULL, dim);

                /* Setting the order 2 for the properties interpolation */
		sfd_set_interpolator_order(sfdu[dim], 2);

                /* Copying the list of cells center in the domain */
		sfd_copy_higtrees_from_center_domain(sfdu[dim], sdp);

		sfd_adjust_facet_ids(sfdu[dim]);

                /* Creating property for the facets */
		psfdu[dim] = psfd_create(sfdu[dim], psdp);

		psfd_compute_sfbi(psfdu[dim]);
                /* Mapping the properties in the domain */
		psfd_synced_mapper(psfdu[dim]);

		/* Create the property itself. */
		if(dim == 0) {
			for(size_t i = 0; i < 2; ++i)
				u[i] = psfd_create_property(psfdu[dim]);
		} else {
			vw[dim-1] = psfd_create_property(psfdu[dim]);
		}
	}

	puts("Filling facet properties with random values...");
	fflush(stdout);
	{
		DEBUG_DIFF_TIME_START;

		// First dimension, with 2 components.
		{
			mp_mapper *m = sfd_get_domain_mapper(sfdu[0]);
			higfit_facetiterator *it;
			for(it = sfd_get_domain_facetiterator(sfdu[0]); !higfit_isfinished(it); higfit_nextfacet(it))
			{
				hig_facet *facet = higfit_getfacet(it);
				int id = mp_lookup(m, hig_get_fid(facet));
				if(id < 0)
					continue;
				for(size_t i = 0; i < 2; ++i)
					dp_set_value(u[i], id, mrand48() * drand48());
			}
			higfit_destroy(it);
		}

		// Second and third dimensions, with one component each
		for(int dim = 1; dim < DIM; dim++) {
			mp_mapper *m = sfd_get_domain_mapper(sfdu[dim]);
			higfit_facetiterator *it;
			for(it = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(it); higfit_nextfacet(it))
			{
				hig_facet *facet = higfit_getfacet(it);
				int id = mp_lookup(m, hig_get_fid(facet));
				if(id < 0)
					continue;
				dp_set_value(vw[dim-1], id, mrand48() * drand48());
			}
			higfit_destroy(it);
		}
		DEBUG_DIFF_TIME_FINISH("  ...duration");
	}

	/* Write the data in HDF5 format. */
	{
		distributed_property *v[2] = {dp1, dp2};
		higio_hdf5 *fd;
		{
			char path[30];
			snprintf(path, 30, "tree.%d.h5", myrank);
			fd = higio_create_hdf5(path);
		}

		for(unsigned i = 0; i < numtrees; ++i) {
			char name[30];
			snprintf(name, 30, "/tree%05u", i);
			printf("Writing tree %u (%s) to file...\n", i, name);
			fflush(stdout);

			DEBUG_DIFF_TIME_START;
			higio_write_tree_hdf5(fd, name, roots[i]);
			DEBUG_DIFF_TIME_FINISH("  ...duration");
		}

		puts("Writing 1d_cell_property to file...");
		fflush(stdout);
		{
			DEBUG_DIFF_TIME_START;
			higio_write_cell_property_hdf5(fd, "/1d_cell_property", sdp, v, 1);
			DEBUG_DIFF_TIME_FINISH("  ...duration");
		}

		puts("Writing 2d_cell_property to file...");
		fflush(stdout);
		{
			DEBUG_DIFF_TIME_START;
			higio_write_cell_property_hdf5(fd, "/2d_cell_property", sdp, v, 2);
			DEBUG_DIFF_TIME_FINISH("  ...duration");
		}

		puts("Writing 2d_facet_property_U to file...");
		fflush(stdout);
		{
			DEBUG_DIFF_TIME_START;
			higio_write_facet_property_hdf5(fd, "/2d_facet_property_U",
				sfdu[0], u, 2);
			DEBUG_DIFF_TIME_FINISH("  ...duration");
		}

		puts("Writing 1d_facet_property_V to file...");
		fflush(stdout);
		{
			DEBUG_DIFF_TIME_START;
			higio_write_facet_property_hdf5(fd, "/1d_facet_property_V",
				sfdu[1], &vw[0], 1);
			DEBUG_DIFF_TIME_FINISH("  ...duration");
		}
#if DIM > 2
		puts("Writing 1d_facet_property_W to file...");
		fflush(stdout);
		{
			DEBUG_DIFF_TIME_START;
			higio_write_facet_property_hdf5(fd, "/1d_facet_property_W",
				sfdu[2], &vw[1], 1);
			DEBUG_DIFF_TIME_FINISH("  ...duration");
		}
#endif
		higio_close_hdf5(fd);
	}

	/* Restore the data from file. */
	puts("\nRestoring data from file:");
	sim_domain *rest_dp = sd_create(NULL);
	sd_set_interpolator_order(rest_dp, 2);

	hig_cell *restored[numtrees];
	distributed_property *dp1d_rest;
	distributed_property *dp2d_rest[2];
	psim_domain *rest_psdp;

	sim_facet_domain *r_sfdu[DIM];
	for(size_t dim = 0; dim < DIM; ++dim) {
		r_sfdu[dim] = sfd_create(NULL, dim);
		sfd_set_interpolator_order(r_sfdu[dim], 2);
	}

	psim_facet_domain *r_psfdu[DIM];
	distributed_property *r_u[2];
	distributed_property *r_vw[DIM-1];

	{
		higio_hdf5 *fd;
		{
			char path[30];
			snprintf(path, 30, "tree.%d.h5", myrank);
			fd = higio_open_hdf5(path);
		}

		for(unsigned i = 0; i < numtrees; ++i) {
			char name[30];
			snprintf(name, 30, "/tree%05u", i);
			printf("Rebuilding tree %u (%s) from file...\n", i, name);
			fflush(stdout);

			DEBUG_DIFF_TIME_START;
			restored[i] = higio_load_tree_hdf5(fd, name);
			DEBUG_DIFF_TIME_FINISH("  ...duration");
			sd_add_higtree(rest_dp, restored[i]);
		}

		rest_psdp = psd_create(rest_dp, pg);

		for(size_t dim = 0; dim < DIM; ++dim) {
			sfd_copy_higtrees_from_center_domain(r_sfdu[dim], rest_dp);
			sfd_adjust_facet_ids(r_sfdu[dim]);
			r_psfdu[dim] = psfd_create(r_sfdu[dim], rest_psdp);

			psfd_compute_sfbi(r_psfdu[dim]);
			psfd_synced_mapper(r_psfdu[dim]);
			if(dim == 0) {
				for(size_t i = 0; i < 2; ++i) {
					r_u[i] = psfd_create_property(r_psfdu[dim]);
				}
			} else {
				r_vw[dim-1] = psfd_create_property(r_psfdu[dim]);
			}
		}

		psd_synced_mapper(rest_psdp);

		dp1d_rest = psd_create_property(rest_psdp);
		dp2d_rest[0] = psd_create_property(rest_psdp);
		dp2d_rest[1] = psd_create_property(rest_psdp);

		puts("Reading 1d_cell_property from file...");
		{
			DEBUG_DIFF_TIME_START;
			higio_read_cell_property_hdf5(fd, "/1d_cell_property", rest_dp, &dp1d_rest, 1);
			dp_sync(dp1d_rest);
			DEBUG_DIFF_TIME_FINISH("  ...duration");
		}

		puts("Reading 2d_cell_property to file...");
		{
			DEBUG_DIFF_TIME_START;
			higio_read_cell_property_hdf5(fd, "/2d_cell_property", rest_dp, dp2d_rest, 2);
			dp_sync(dp2d_rest[0]);
			dp_sync(dp2d_rest[1]);
			DEBUG_DIFF_TIME_FINISH("  ...duration");
		}

		puts("Reading 2d_facet_property_U from file...");
		{
			DEBUG_DIFF_TIME_START;
			higio_read_facet_property_hdf5(fd, "/2d_facet_property_U", r_sfdu[0],
					r_u, 2);
			dp_sync(r_u[0]);
			dp_sync(r_u[1]);
			DEBUG_DIFF_TIME_FINISH(" ...duration");
		}

		puts("Reading 1d_facet_property_V from file...");
		{
			DEBUG_DIFF_TIME_START;
			higio_read_facet_property_hdf5(fd, "/1d_facet_property_V", r_sfdu[1],
					&r_vw[0], 1);
			dp_sync(r_vw[0]);
			DEBUG_DIFF_TIME_FINISH(" ...duration");
		}
#if DIM > 2
		puts("Reading 1d_facet_property_W from file...");
		{
			DEBUG_DIFF_TIME_START;
			higio_read_facet_property_hdf5(fd, "/1d_facet_property_W", r_sfdu[2],
					&r_vw[1], 1);
			dp_sync(r_vw[1]);
			DEBUG_DIFF_TIME_FINISH(" ...duration");
		}
#endif
		higio_close_hdf5(fd);
	}

	/* Check if the restored trees has the same properties as the original. */
	for(unsigned i = 0; i < numtrees; ++i) {
		compare_trees(roots[i], restored[i]);
	}

	/* Check if the restored cell data is identical to the original */
	{
		mp_mapper *orig_m = sd_get_domain_mapper(sdp);
		mp_mapper *rest_m = sd_get_domain_mapper(rest_dp);

		higcit_celliterator *orig_it = sd_get_domain_celliterator(sdp);
		higcit_celliterator *rest_it = sd_get_domain_celliterator(rest_dp);
		{
			long co, ca;
			co = higcit_count(orig_it);
			ca = higcit_count(rest_it);
			printf("Original cell property count: %ld\n", co);
			printf("Restored cell property count: %ld\n", ca);
			assert(co == ca);
		}
		higcit_destroy(orig_it);
		higcit_destroy(rest_it);

		orig_it = sd_get_domain_celliterator(sdp);
		rest_it = sd_get_domain_celliterator(rest_dp);

		while(!higcit_isfinished(rest_it)) {
			assert(!higcit_isfinished(orig_it));

			hig_cell *orig = higcit_getcell(orig_it);
			hig_cell *rest = higcit_getcell(rest_it);

			assert(memcmp(orig->ids, rest->ids, sizeof(orig->ids)) == 0);

			// The tests bellow are redundant...
			assert(memcmp(orig->lowpoint, rest->lowpoint, sizeof(orig->lowpoint)) == 0);
			assert(memcmp(orig->highpoint, rest->highpoint, sizeof(orig->highpoint)) == 0);

			assert(memcmp(orig->numcells, rest->numcells, sizeof(orig->numcells)) == 0);
			assert(hig_get_number_of_children(orig) == hig_get_number_of_children(rest));
			// ...end of redundant tests.

			hig_cell *oparent = hig_get_parent(orig);
			hig_cell *rparent = hig_get_parent(rest);
			assert((oparent && rparent) || (!oparent && !rparent));
			if(rparent) {
				assert(hig_get_cid(rparent) == hig_get_cid(oparent));
			}

			assert(hig_get_number_of_children(rest) == 0);
			{
				/* Check restored properties. */
				size_t orig_id = mp_lookup(orig_m, hig_get_cid(orig));
				size_t rest_id = mp_lookup(rest_m, hig_get_cid(rest));

				assert(dp_get_value(dp1d_rest, rest_id) == dp_get_value(dp1, orig_id));
				assert(dp_get_value(dp2d_rest[0], rest_id) == dp_get_value(dp1, orig_id));
				assert(dp_get_value(dp2d_rest[1], rest_id) == dp_get_value(dp2, orig_id));
			}

			higcit_nextcell(orig_it);
			higcit_nextcell(rest_it);
		}
		assert(higcit_isfinished(orig_it));

		higcit_destroy(orig_it);
		higcit_destroy(rest_it);
	}


	/* Check restored facet properties. */
	for(size_t dim = 0; dim < DIM; ++dim) {
		mp_mapper *orig_fm = sfd_get_domain_mapper(sfdu[dim]);
		mp_mapper *rest_fm = sfd_get_domain_mapper(r_sfdu[dim]);

		higfit_facetiterator *orig_it = sfd_get_domain_facetiterator(sfdu[dim]);
		higfit_facetiterator *rest_it = sfd_get_domain_facetiterator(r_sfdu[dim]);
		{
			long co, ca;
			co = higfit_count(orig_it);
			ca = higfit_count(rest_it);
			printf("Original facet property count for dim %zd: %ld\n", dim, co);
			printf("Restored facet property count for dim %zd: %ld\n", dim, ca);
			assert(co == ca);
		}
		higfit_destroy(orig_it);
		higfit_destroy(rest_it);

		orig_it = sfd_get_domain_facetiterator(sfdu[dim]);
		rest_it = sfd_get_domain_facetiterator(r_sfdu[dim]);

		while(!higfit_isfinished(rest_it)) {
			assert(!higfit_isfinished(orig_it));

			hig_facet *orig = higfit_getfacet(orig_it);
			hig_facet *rest = higfit_getfacet(rest_it);

			assert(hig_get_facet_dir(orig) == hig_get_facet_dir(rest));
			assert(hig_get_facet_dim(orig) == hig_get_facet_dim(rest));
			{
				Point oc, rc;
				hig_get_center(hig_get_facet_cell(orig), oc);
				hig_get_center(hig_get_facet_cell(rest), rc);
				assert(memcmp(oc, rc, sizeof oc) == 0);
			}

			size_t orig_id = mp_lookup(orig_fm, hig_get_fid(orig));
			size_t rest_id = mp_lookup(rest_fm, hig_get_fid(rest));

			if(dim == 0) {
				assert(dp_get_value(r_u[0], rest_id) == dp_get_value(u[0], orig_id));
				assert(dp_get_value(r_u[1], rest_id) == dp_get_value(u[1], orig_id));
			} else {
				assert(dp_get_value(r_vw[dim-1], rest_id) == dp_get_value(vw[dim-1], orig_id));
			}

			higfit_nextfacet(orig_it);
			higfit_nextfacet(rest_it);
		}
		assert(higfit_isfinished(orig_it));

		higfit_destroy(orig_it);
		higfit_destroy(rest_it);
	}

#ifdef NDEBUG
	puts("Tests not enabled.");
#else
	puts("All tests passed!");
#endif
}

/*********************************************************************************************/
/***    Main Program: End                                                                  ***/
/*********************************************************************************************/

