#include<string.h>
#include<assert.h>
#include<stdlib.h>
#include<stdio.h>
#include<stdarg.h>
#include<hdf5.h>

#include "Debug-c.h"
#include "utils.h"
#include "higtree.h"
#include "coord.h"
#include "higtree-iterator.h"
#include "higtree-parallel.h"
#include "higtree-serialize.h"
#include "lbal.h"

#include "higtree-io.h"

void higio_print(FILE *f, hig_cell *cell) {
	if (cell->children == NULL) {
		fprintf(f, "(");
		co_print_coord(f, cell->lowpoint, ", ");
		fprintf(f, ", ");
		co_print_coord(f, cell->highpoint, ", ");
		fprintf(f, " posincell: ");
		fprintf(f, " %d ", cell->posinparent);
		fprintf(f, ")\n");
	} else {
		int numChildren = hig_get_number_of_children(cell);
		for(int j = 0; j < numChildren; j++) {
			higio_print(f, hig_get_child(cell, j));
		}
	}
}

static void __higio_write_to_file(FILE *f, hig_cell *p) {
	if(p->children == NULL) {
		fprintf(f, " 0\n");
	} else {
		for(int i = 0; i < DIM; i++) {
			fprintf(f, " %d", p->numcells[i]);
		}
		fprintf(f, "\n");
		int numChildren = hig_get_number_of_children(p);
		for(int j = 0; j < numChildren; j++) {
			__higio_write_to_file(f, hig_get_child(p, j));
		}
	}
}

void higio_write_to_file(FILE * f, hig_cell *root) {
	assert(root != NULL);
	fprintf(f, "%d", DIM);
	for(int i = 0; i < DIM; i++) {
		fprintf(f, " %lf", root->lowpoint[i]);
	}
	for(int i = 0; i < DIM; i++) {
		fprintf(f, " %lf", root->highpoint[i]);
	}
	fprintf(f, "\n");
	__higio_write_to_file(f, root);
}

static void __higio_read_from_file(FILE *f, hig_cell *p) {
	int numcells[DIM];
	fscanf(f, "%d", &numcells[0]);
	if (numcells[0] != 0) {
		for(int i = 1; i < DIM; i++) {
			fscanf(f, "%d", &numcells[i]);
			assert(numcells[i] > 0);
		}
		hig_refine_uniform(p, numcells);
		int numChildren = hig_get_number_of_children(p);
		for(int j = 0; j < numChildren; j++) {
			__higio_read_from_file(f, hig_get_child(p, j));
		}
	}
}

hig_cell * higio_read_from_file(FILE *f) {
	int dim;
	fscanf(f, "%d", &dim);
	assert(dim == DIM);
	Point lowpoint;
	Point highpoint;
	for(int i = 0; i < DIM; i++) {
		fscanf(f, "%lf", &lowpoint[i]);
	}
	for(int i = 0; i < DIM; i++) {
		fscanf(f, "%lf", &highpoint[i]);
	}
	hig_cell * root = hig_create_root(lowpoint, highpoint);
	__higio_read_from_file(f, root);
	return root;
}

hig_cell * higio_read_from_amr_info(higio_amr_info *mi)
{
	hig_cell *root = hig_create_root(mi->l, mi->h);

	for(int i = 0; i < DIM; i++) {
		assert(FLT_EQ(mi->levels[0].delta[i],
			(mi->h[i] - mi->l[i])/mi->levels[0].patches[0].patchsize[i]));
	}

	hig_refine_uniform(root, mi->levels[0].patches[0].patchsize);

	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 2);
	for(int l = 1; l < mi->numlevels; l++) {
		__higio_level_info *level = mi->levels + l;
		for(int p = 0; p < level->numpatches; p++) {
			__higio_patch_info *patch = level->patches + p;

			Point initpatch;
			POINT_MULT(initpatch, level->delta, patch->initcell);
			POINT_ADD(initpatch, initpatch, mi->l);

			int k[DIM];
			POINT_ASSIGN_SCALAR(k, 0);
			while(1) {
				Point offset;
				POINT_MULT(offset, level->delta, k);
				Point point;
				POINT_ADD(point, initpatch, offset);
				hig_cell *c = hig_get_cell_with_point(root, point);
				Point cdelta;
				hig_get_delta(c, cdelta);
				int isrefinedenough = 1;
				for(int i = 0; i < DIM; i++) {
					if(FLT_LT(level->delta[i], cdelta[i])) {
						isrefinedenough = 0;
						break;
					}
				}
				if(!isrefinedenough) {
					hig_refine_uniform(c, nc);
				}
				int i = 0;
				while(i < DIM && k[i] + 2 >= patch->patchsize[i] - 1) {
					k[i] = 0;
					i++;
				}
				if (i < DIM) {
					k[i]+=2;
				} else {
					break;
				}
			}
		}
	}
	return root;
}

hig_cell * higio_read_from_amr(FILE *fd) {
	higio_amr_info *mi = higio_read_amr_info(fd);
	hig_cell *ret = higio_read_from_amr_info(mi);
	higio_amr_info_destroy(mi);

	return ret;
}

#if DIM == 2

void higio_print_celliterator_in_vtk(FILE * f, higcit_celliterator *it) {
	long numleafs = higcit_count_without_advancing(it);
	fprintf(f, "# vtk DataFile Version 3.0\n");
	fprintf(f, "higtree\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET POLYDATA\n");
	fprintf(f, "POINTS %ld float\n", 4*numleafs);
	hig_cell *c;
	for (; !higcit_isfinished(it); higcit_nextcell(it)) {
		c = higcit_getcell(it);
		fprintf(f, "%lf %lf 0\n", c->lowpoint[0], c->lowpoint[1]);
		fprintf(f, "%lf %lf 0\n", c->lowpoint[0], c->highpoint[1]);
		fprintf(f, "%lf %lf 0\n", c->highpoint[0], c->highpoint[1]);
		fprintf(f, "%lf %lf 0\n", c->highpoint[0], c->lowpoint[1]);
	}
	fprintf(f, "POLYGONS %ld %ld\n", numleafs, 5*numleafs);
	for(int i = 0; i < numleafs; i++) {
		fprintf(f, "%d %d %d %d %d\n", 4, 4*i, 4*i+1, 4*i+2, 4*i+3);
	}
}

void higio_print_in_vtk(FILE * f, hig_cell *cell) {
	higcit_celliterator *it = higcit_create_all_leaves(cell);
	higio_print_celliterator_in_vtk2d(f, it);
	higcit_destroy(it);
}

hig_cell * higio_read_from_amr2d(FILE *fd) {
	return higio_read_from_amr(fd);
}

void higio_print_celliterator_in_vtk2d(FILE * f, higcit_celliterator *it) {
	higio_print_celliterator_in_vtk(f, it);
}

void higio_print_in_vtk2d(FILE * f, hig_cell *cell) {
	higio_print_in_vtk(f, cell);
}

#elif DIM == 3

void higio_print_celliterator_in_vtk(FILE * f, higcit_celliterator *it) {
	long numleafs = higcit_count_without_advancing(it);
	fprintf(f, "# vtk DataFile Version 3.0\n");
	fprintf(f, "higtree\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET POLYDATA\n");
	fprintf(f, "POINTS %ld float\n", 8*numleafs);
	hig_cell *c;
	for (; !higcit_isfinished(it); higcit_nextcell(it)) {
		c = higcit_getcell(it);
		fprintf(f, "%lf %lf %lf\n", c->lowpoint [0], c->lowpoint [1], c->lowpoint [2]);
		fprintf(f, "%lf %lf %lf\n", c->highpoint[0], c->lowpoint [1], c->lowpoint [2]);
		fprintf(f, "%lf %lf %lf\n", c->highpoint[0], c->highpoint[1], c->lowpoint [2]);
		fprintf(f, "%lf %lf %lf\n", c->lowpoint [0], c->highpoint[1], c->lowpoint [2]);
		fprintf(f, "%lf %lf %lf\n", c->lowpoint [0], c->lowpoint [1], c->highpoint[2]);
		fprintf(f, "%lf %lf %lf\n", c->highpoint[0], c->lowpoint [1], c->highpoint[2]);
		fprintf(f, "%lf %lf %lf\n", c->highpoint[0], c->highpoint[1], c->highpoint[2]);
		fprintf(f, "%lf %lf %lf\n", c->lowpoint [0], c->highpoint[1], c->highpoint[2]);
	}
	fprintf(f, "POLYGONS %ld %ld\n", 6*numleafs, 5*6*numleafs);
	for(int i = 0; i < numleafs; i++) {
		fprintf(f, "%d %d %d %d %d\n", 4, 8*i+0, 8*i+3, 8*i+7, 8*i+4); //x == 0
		fprintf(f, "%d %d %d %d %d\n", 4, 8*i+1, 8*i+2, 8*i+6, 8*i+5); //x == 1
		fprintf(f, "%d %d %d %d %d\n", 4, 8*i+0, 8*i+1, 8*i+5, 8*i+4); //y == 0
		fprintf(f, "%d %d %d %d %d\n", 4, 8*i+3, 8*i+2, 8*i+6, 8*i+7); //y == 1
		fprintf(f, "%d %d %d %d %d\n", 4, 8*i+0, 8*i+1, 8*i+2, 8*i+3); //z == 0
		fprintf(f, "%d %d %d %d %d\n", 4, 8*i+4, 8*i+5, 8*i+6, 8*i+7); //z == 1
	}
}

void higio_print_in_vtk(FILE * f, hig_cell *cell) {
	higcit_celliterator *it = higcit_create_all_leaves(cell);
	higio_print_celliterator_in_vtk3d(f, it);
	higcit_destroy(it);
}

hig_cell * higio_read_from_amr3d(FILE *fd) {
	return higio_read_from_amr(fd);
}

void higio_print_celliterator_in_vtk3d(FILE * f, higcit_celliterator *it) {
	higio_print_celliterator_in_vtk(f, it);
}

void higio_print_in_vtk3d(FILE * f, hig_cell *cell) {
	higio_print_in_vtk(f, cell);
}
#endif

void higio_write_leaves(FILE *f, hig_cell *root) {
	higcit_celliterator *it;
	hig_cell *c;
	for (it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		c = higcit_getcell(it);
		if (DIM >= 2) {
			fprintf(f, "%2.5f %2.5f %2.5f %2.5f\n", c->lowpoint[0], c->lowpoint[1], c->highpoint[0], c->highpoint[1]);
		} else {
			fprintf(f, "%2.5f %2.5f %2.5f %2.5f\n", c->lowpoint[0], 0.0, c->highpoint[0], 0.0);
		}
	}
	higcit_destroy(it);
}

higio_amr_info *higio_read_amr_info(FILE *fd) {
	DECL_AND_ALLOC(higio_amr_info, mi, 1);
	for(int dim = 0; dim < DIM; dim++) {
		fscanf(fd, "%lf", &mi->l[dim]);
		fscanf(fd, "%lf", &mi->h[dim]);
	}
	fscanf(fd, "%d", &mi->numlevels);
	ALLOC(__higio_level_info, mi->levels, mi->numlevels);
	for(int l = 0; l < mi->numlevels; l++) {
		for(int dim = 0; dim < DIM; dim++) {
			fscanf(fd, "%lf", &mi->levels[l].delta[dim]);
		}
		fscanf(fd, "%d", &mi->levels[l].numpatches);
		ALLOC(__higio_patch_info, mi->levels[l].patches, mi->levels[l].numpatches);
		for(int p = 0; p < mi->levels[l].numpatches; p++) {
			for(int dim = 0; dim < DIM; dim++) {
				fscanf(fd, "%d", &mi->levels[l].patches[p].initcell[dim]);
			}
			for(int dim = 0; dim < DIM; dim++) {
				fscanf(fd, "%d", &mi->levels[l].patches[p].patchsize[dim]);
			}
		}
	}
	return mi;
}

higio_amr_info *higio_clone_amr_info(higio_amr_info *orig) {
	DECL_AND_ALLOC(higio_amr_info, mi, 1);

	POINT_ASSIGN(mi->l, orig->l);
	mi->numlevels = orig->numlevels;

	ALLOC(__higio_level_info, mi->levels, mi->numlevels);
	for(int l = 0; l < mi->numlevels; l++) {
		POINT_ASSIGN(mi->levels[l].delta, orig->levels[l].delta);

		mi->levels[l].numpatches = orig->levels[l].numpatches;
		ALLOC(__higio_patch_info, mi->levels[l].patches, mi->levels[l].numpatches);
		for(int p = 0; p < mi->levels[l].numpatches; p++) {
			for(int dim = 0; dim < DIM; dim++) {
				mi->levels[l].patches[p].initcell[dim] = orig->levels[l].patches[p].initcell[dim];
				mi->levels[l].patches[p].patchsize[dim] = orig->levels[l].patches[p].patchsize[dim];
			}
		}
	}
	return mi;
}

void higio_amr_info_destroy(higio_amr_info *mi) {
	for(int l = 0; l < mi->numlevels; l++) {
		free(mi->levels[l].patches);
	}
	free(mi->levels);
	free(mi);
}

void __higio_amr_info_print(higio_amr_info *mi) {
	printf("%f %f %f %f\n", mi->l[0], mi->h[0], mi->l[1], mi->h[1]);
	printf("%d\n", mi->numlevels);
	for(int l = 0; l < mi->numlevels; l++) {
		printf("%d\n", mi->levels[l].numpatches);
		for(int p = 0; p < mi->levels[l].numpatches; p++) {
			printf("%d %d ", mi->levels[l].patches[p].initcell[0], mi->levels[l].patches[p].initcell[1]);
			printf("%d %d\n", mi->levels[l].patches[p].patchsize[0], mi->levels[l].patches[p].patchsize[1]);
		}
	}
}

void higio_compute_cells_per_coarser_block_from_amr_info(higio_amr_info *mi, int numcells[], int numblocks) {
	for(int c = 0; c < numblocks; c++) {
		numcells[c] = 1;
	}
	int div = 1;
	int nc[DIM];
	POINT_ASSIGN(nc, mi->levels[0].patches[0].patchsize);
	for(int l = 1; l < mi->numlevels; l++) {
		for(int p = 0; p < mi->levels[l].numpatches; p++) {
			assert(mi->levels[l].patches[p].initcell[0] % 2 == 1);
			int cl[DIM];
			int ch[DIM];
			for(int d = 0; d < DIM; d++) {
				cl[d] = mi->levels[l].patches[p].initcell[d] - 1;
				ch[d] = cl[d] + mi->levels[l].patches[p].patchsize[d];
				cl[d] >>= 1;
				ch[d] >>= 1;
			}
			int ind[DIM];
			POINT_ASSIGN(ind, cl);
			while(1) {
				int pos = 0;
				for(int dim = 0; dim < DIM; dim++) {
					pos *= nc[dim];
					pos += (ind[dim] >> (l-1));
				}
				numcells[pos] += (1<<DIM)-1;
				int dim = 0;
				while(dim < DIM && ind[dim] >= ch[dim] - 1) {
					ind[dim] = cl[dim];
					dim++;
				}
				if (dim < DIM) {
					ind[dim]++;
				} else {
					break;
				}
			}
		}
	}
}

load_balancer* higio_partition_from_amr_info(partition_graph *pg, higio_amr_info *mi) {
	MPI_Comm comm = pg_get_MPI_comm(pg);
	int rank;
	MPI_Comm_rank(comm, &rank);

	load_balancer *lb = lb_create(comm, 1);
	if(rank == 0) {
		hig_cell *root = higio_read_from_amr_info(mi);
		lb_add_input_tree(lb, root, true, 0);
	}
	lb_calc_partition(lb, pg);

	return lb;
}

hig_cell *higio_read_bc_from_amr(higio_amr_info *mi, int bcdim, int bcdir) {
	Point lowpoint, highpoint;
	POINT_ASSIGN(lowpoint, mi->l);
	POINT_ASSIGN(highpoint, mi->h);
	if (bcdir == 0) {
		lowpoint[bcdim]  = mi->l[bcdim] - EPSDELTA;
		highpoint[bcdim] = mi->l[bcdim] + EPSDELTA;
	} else {
		lowpoint[bcdim]  = mi->h[bcdim] - EPSDELTA;
		highpoint[bcdim] = mi->h[bcdim] + EPSDELTA;
	}
	hig_cell *root = hig_create_root(lowpoint, highpoint);
	int nc[DIM];
	POINT_ASSIGN(nc, mi->levels[0].patches[0].patchsize);
	nc[bcdim] = 1;
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	nc[bcdim] = 1;
	for(int l = 1; l < mi->numlevels; l++) {
		Point delta;
		POINT_ASSIGN(delta, mi->levels[l].delta);
		int numpatches = mi->levels[l].numpatches;
		for(int p = 0; p < numpatches; p++) {
			int initcell[DIM];
			int patchsize[DIM];
			POINT_ASSIGN(initcell, mi->levels[l].patches[p].initcell);
			int istorefine = 0;
			if(bcdir == 0) {
				if(initcell[bcdim] == 1) {
					istorefine = 1;
				}
			} else {
				int lastcell = initcell[bcdim] + mi->levels[l].patches[p].patchsize[bcdim];
				//DEBUG_INSPECT(lastcell, %d);
				//DEBUG_INSPECT(mi->levels[0].patches[0].patchsize[bcdim] * (1 << l) + 1, %d);
				if (lastcell == mi->levels[0].patches[0].patchsize[bcdim] * (1 << l) + 1) {
					istorefine = 1;
				}
			}
			if (istorefine) {
				POINT_ASSIGN(patchsize, mi->levels[l].patches[p].patchsize);
				patchsize[bcdim] = 1;
				Point initpatch;
				POINT_MULT(initpatch, delta, initcell);
				POINT_ADD(initpatch, initpatch, mi->l);
				int k[DIM];
				int firstid[DIM];
				POINT_ASSIGN_SCALAR(firstid, 0);
				int lastid[DIM];
				POINT_ASSIGN(lastid, patchsize);
				POINT_SUB_SCALAR(lastid, lastid, 1);
				if (bcdir == 0) {
					lastid[bcdim] = firstid[bcdim];
				} else {
					firstid[bcdim] = lastid[bcdim];
				}
				POINT_ASSIGN(k, firstid);
				while(1) {
					Point offset;
					POINT_MULT(offset, delta, k);
					Point point;
					POINT_ADD(point, initpatch, offset);
					if(bcdir == 0) {
						point[bcdim] = lowpoint[bcdim];
					} else {
						point[bcdim] = highpoint[bcdim];
					}
					hig_cell *c = hig_get_cell_with_point(root, point);
					if (c != NULL) {
						Point cdelta;
						hig_get_delta(c, cdelta);
						int isrefinedenough = 1;
						for(int dim = 0; dim < DIM; dim++) {
							if(FLT_LT(delta[dim], cdelta[dim])) {
								isrefinedenough = 0;
								break;
							}
						}
						if(!isrefinedenough) {
							hig_refine_uniform(c, nc);
						}
					}
					int dim = 0;
					while(dim < DIM && k[dim] + 2 >= lastid[dim]) {
						k[dim] = firstid[dim];
						dim++;
					}
					if (dim < DIM) {
						k[dim]+=2;
					} else {
						break;
					}
				}
			}
		}
	}
	return root;
}


void higio_get_domain_lowpoint_from_mesh_info(higio_amr_info *mi, Point l) {
	POINT_ASSIGN(l, mi->l);
}

void higio_get_domain_delta_first_level_from_mesh_info(higio_amr_info *mi, Point delta) {
	POINT_ASSIGN(delta, mi->levels[0].delta);
}

struct higio_hdf5 {
	hid_t fd;
};

//! Returns the equivalent HDF5 type for real
static hid_t _get_hdf5_real() {
	if(sizeof(real) == sizeof(float))
		return H5T_NATIVE_FLOAT;
	else if(sizeof(real) == sizeof(double))
		return H5T_NATIVE_DOUBLE;
	else {
		assert(sizeof(real) == sizeof(long double));
		return H5T_NATIVE_LDOUBLE;
	}
}

//! Returns the equivalent HDF5 type for uniqueid
static inline hid_t _get_hdf5_uniqueid() {
	if(sizeof(uniqueid) == sizeof(int8_t))
		return H5T_NATIVE_B8;
	else if(sizeof(uniqueid) == sizeof(int16_t))
		return H5T_NATIVE_B16;
	else if(sizeof(uniqueid) == sizeof(int32_t))
		return H5T_NATIVE_B32;
	else {
		assert(sizeof(uniqueid) == sizeof(int64_t));
		return H5T_NATIVE_B64;
	}
}

// TODO: write points as a separated dataset; verify points dataset interface for this...

higio_hdf5* higio_create_hdf5(const char *filename)
{
	DECL_AND_ALLOC(higio_hdf5, ctx, 1);

	/* HDF5 library compulsivelly call this function everywhere,
	 * so it seems wise to call it here, too. */
	H5open();
	ctx->fd = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	// TODO: write attributes to the file with dimmension, version and identification
	// as a HigTree file, so to it be specific.

	return ctx;
}

higio_hdf5* higio_open_hdf5(const char *filename)
{
	DECL_AND_ALLOC(higio_hdf5, ctx, 1);

	/* HDF5 library compulsivelly call this function everywhere,
	 * so it seems wise to call it here, too. */
	H5open();
	ctx->fd = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

	// TODO: check if the file is a valid HigTree file, check if dimmension and version
	// matches.

	return ctx;
}

void higio_close_hdf5(higio_hdf5* ctx)
{
	H5Fclose(ctx->fd);

	free(ctx);
}

static const hsize_t chunk_size = 64*1024;

void higio_write_tree_hdf5(higio_hdf5* ctx, const char* path, hig_cell* root)
{
	static hid_t dt_higcell = 0;
	if(!dt_higcell) {
		// Create the datatypes and dataspaces used to store the tree

		// Point and children cell count
		hid_t dt_point;
		hid_t dt_cellcount;
		{
			const hsize_t dim = DIM;
			dt_point = H5Tarray_create(_get_hdf5_real(), 1, &dim);
			dt_cellcount = H5Tarray_create(H5T_NATIVE_INT, 1, &dim);
		}

		// Cell
		{
			// TODO: find a way to deallocate this (maybe not necessary)
			// can be done with atexit()

			dt_higcell = H5Tcreate(H5T_COMPOUND, sizeof(hig_serial_tree));
			H5Tinsert(dt_higcell, "lowpoint",
				HOFFSET(hig_serial_tree, lowpoint), dt_point);
			H5Tinsert(dt_higcell, "highpoint",
				HOFFSET(hig_serial_tree, highpoint), dt_point);
			H5Tinsert(dt_higcell, "numcells",
				HOFFSET(hig_serial_tree, numcells), dt_cellcount);
		}

		H5Tclose(dt_cellcount);
	}

	// The dataspace (memory accessor and the on-file layout)
	hid_t ds_mem, ds_file;
	{
		hsize_t unlimited = H5S_UNLIMITED;
		ds_mem = H5Screate_simple(1, &chunk_size, &unlimited);
		ds_file = H5Scopy(ds_mem);
	}

	// Prepare the dataset to receive data
	hid_t dset;
	{
		hid_t dset_plist = H5Pcreate(H5P_DATASET_CREATE);
		H5Pset_chunk(dset_plist, 1, &chunk_size);

		unsigned int compression_level = 1; // Fastest compression
		H5Pset_filter(dset_plist, H5Z_FILTER_DEFLATE, H5Z_FLAG_OPTIONAL, 1, &compression_level);

		/* Create the dataspace and the chunked dataset */
		dset = H5Dcreate(ctx->fd, path, dt_higcell, ds_mem, H5P_DEFAULT, dset_plist, H5P_DEFAULT);
		H5Pclose(dset_plist);
	}

	// Model the data to the correct format
	{
		hig_serial_tree buff[chunk_size];

		higcit_celliterator *iter;
		hsize_t i = 0;
		hsize_t written_size = 0;

		for (iter = higcit_create_breadth_first(root); !higcit_isfinished(iter); higcit_nextcell(iter))
		{
			hig_cell *c = higcit_getcell(iter);
			hig_get_lowpoint(c, buff[i].lowpoint);
			hig_get_highpoint(c, buff[i].highpoint);
			memcpy(buff[i].numcells, c->numcells, sizeof(buff[i].numcells));

			hig_cell *parent = hig_get_parent(c);

			if(++i == chunk_size) {
				hsize_t offset = written_size;
				i = 0;

				// Set the new size of the dataset
				written_size += chunk_size;
				H5Dset_extent(dset, &written_size);

				// Configure the hyperslab of where to write
				H5Sclose(ds_file);
				ds_file = H5Dget_space(dset);
				H5Sselect_hyperslab(ds_file, H5S_SELECT_SET, &offset, NULL, &chunk_size, NULL);

				H5Dwrite(dset, dt_higcell, ds_mem, ds_file, H5P_DEFAULT, buff);
			}
		}

		// Last write, incomplente chunk...
		if(i) {
			hsize_t offset = written_size;

			// Reconfigure memory dataspace to fit the number of elements in buffer
			hsize_t unlimited = H5S_UNLIMITED;
			H5Sset_extent_simple(ds_mem, 1, &i, &unlimited);

			// Set the new size of the dataset
			written_size += i;
			H5Dset_extent(dset, &written_size);

			// Configure the hyperslab of where to write
			H5Sclose(ds_file);
			ds_file = H5Dget_space(dset);
			H5Sselect_hyperslab(ds_file, H5S_SELECT_SET, &offset, NULL, &i, NULL);

			H5Dwrite(dset, dt_higcell, ds_mem, ds_file, H5P_DEFAULT, buff);
		}
		higcit_destroy(iter);

	}
	H5Dclose(dset);
	H5Sclose(ds_mem);
	H5Sclose(ds_file);
}

hig_cell *higio_load_tree_hdf5(higio_hdf5* ctx, const char* path)
{
	hig_cell *tree = NULL;

	// Open tree dataset
	hid_t dset = H5Dopen(ctx->fd, path, H5P_DEFAULT);
	hid_t ds_file = H5Dget_space(dset);

	// Get the chunk size
	hsize_t chunk_size;
	{
		hid_t dset_plist = H5Dget_create_plist(dset);
		H5Pget_chunk(dset_plist, 1, &chunk_size);
		H5Pclose(dset_plist);
	}

	// Get the total size of the tree
	hsize_t total_size;
	H5Sget_simple_extent_dims(ds_file, &total_size, NULL);

	if(total_size) {
		hig_serial_tree buff[chunk_size];
		hid_t dtype = H5Dget_type(dset);
		hid_t ds_mem = H5Screate_simple(1, &chunk_size, &chunk_size);

		// Calculate the number of chunks the dataset is split into
		size_t num_chunks = total_size / chunk_size;
		if(total_size % chunk_size)
			++num_chunks;

		higcit_celliterator *iter = NULL;
		hig_cell *parent = NULL;
		size_t numchildren;
		size_t child_iter;

		hsize_t rem_size = total_size;
		for(size_t chunk = 0; chunk < num_chunks; ++chunk) {
			hsize_t offset = total_size - rem_size;
			hsize_t read_size;
			if(rem_size >= chunk_size) {
				// Reading a whole chunk
				read_size = chunk_size;
				rem_size -= chunk_size;
			} else {
				// Final reading with a smaller chunk
				read_size = rem_size;
				H5Sset_extent_simple(ds_mem, 1, &rem_size, &rem_size);
				rem_size = 0;
			}

			H5Sselect_hyperslab(ds_file, H5S_SELECT_SET, &offset, NULL, &read_size, NULL);
			H5Dread(dset, dtype, ds_mem, ds_file, H5P_DEFAULT, buff);

			// Build the root node if this is the first read.
			size_t i = 0;
			if(!tree) {
				tree = hig_create_root(buff[0].lowpoint, buff[0].highpoint);
				hig_refine_empty(tree, buff[0].numcells);

				iter = higcit_create_breadth_first(tree);
				parent = tree; //higcit_getcell(iter);

				numchildren = hig_get_number_of_children(tree);
				child_iter = 0;

				++i;
			}

			// Use the loaded nodes to build the tree.
			for(; i < read_size; ++i) {
				// Get the node that will receive the loaded data
				if(child_iter >= numchildren) {
					do {
						// Sanity check.
						if(higcit_isfinished(iter)) {
							char filename[100];

							H5Fget_name(ctx->fd, filename, 100);
							fprintf(stderr, "Error: Tree stored in %s:%s is inconsistent: stored nodes does not fit in the claimed tree size.\n", filename, path);

							hig_destroy(tree);
							tree = NULL;
							goto clean_and_exit;
						}

						higcit_nextcell(iter);
						parent = higcit_getcell(iter);
						numchildren = hig_get_number_of_children(parent);
					} while(!numchildren);

					child_iter = 0;
				}

				hig_cell* target = hig_create_leaf(buff[i].lowpoint, buff[i].highpoint, parent, child_iter++);
				hig_refine_empty(target, buff[i].numcells);
			}
		}

		// Another sanity check.
		if(rem_size) {
			char filename[100];

			H5Fget_name(ctx->fd, filename, 100);
			fprintf(stderr, "Error: Tree stored in %s:%s is inconsistent: there are not enought stored nodes to fill-in the claimed tree size.\n", filename, path);

			hig_destroy(tree);
			tree = NULL;
		}

clean_and_exit:
		higcit_destroy(iter);
		H5Sclose(ds_mem);
		H5Tclose(dtype);
	}

	H5Sclose(ds_file);
	H5Dclose(dset);

	return tree;
}

struct _prop_writer {
	int rank;
	hid_t ds_mem, ds_file;
	hid_t dset;
	hid_t dtype;
	hsize_t size[2];
	hsize_t max_size[2];
	hsize_t dimension;

	// Stuff to iterate through the leafs of all trees
	hsize_t offset[2];
	hid_t datatype;
};

static void _prop_writer_init_dataset(higio_hdf5 *ctx, const char *dataset_path,
		hsize_t max_size, size_t prop_dimension, hid_t datatype,
		bool use_var_env_compression, struct _prop_writer* wprops)
{
	assert(prop_dimension > 0);
	assert(max_size > 0);

	enum compression_type {
		DEFAULT,
		NONE,
		GZIP,
		MAFISC,
		UNKNOWN
	};
	static enum compression_type env_comp_type = UNKNOWN;

	wprops->rank = prop_dimension == 1 ? 1 : 2;

	wprops->size[0] = 0;
	wprops->size[1] = prop_dimension;

	wprops->max_size[0] = max_size;
	wprops->max_size[1] = prop_dimension;

	memset(wprops->offset, 0, sizeof wprops->offset);
	wprops->datatype = datatype;

	wprops->ds_file = H5Screate_simple(wprops->rank, wprops->size, wprops->max_size);

	hsize_t csize[2] = {chunk_size, prop_dimension};
	wprops->ds_mem = H5Screate_simple(wprops->rank, csize, csize);

	hid_t dset_plist = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_chunk(dset_plist, wprops->rank, csize);

	enum compression_type used = use_var_env_compression ? env_comp_type : DEFAULT;

	if(used == UNKNOWN) {
		// Find what compression to use, per environment variable
		char *choice = getenv("HIGTREE_OUTPUT_COMPRESSION");
		const char *const choices[] = {
			"DEFAULT",
			"NONE",
			"GZIP",
			"MAFISC"
		};

		env_comp_type = DEFAULT;
		if(choice) {
			for(size_t i = 1; i < (sizeof(choices)/sizeof(char*)); ++i) {
				if(strcmp(choices[i], choice) == 0) {
					env_comp_type = i;
					break;
				}
			}
		}
		used = env_comp_type;
	}

	unsigned int filter_option;
	const unsigned int GZIP_LEVEL = 4;
	switch(used) {
	case DEFAULT:
		{
			static enum compression_type available = UNKNOWN;
			switch(available) {
			case GZIP:
				H5Pset_deflate(dset_plist, GZIP_LEVEL);
				break;
			case MAFISC:
				filter_option = 0;
				H5Pset_filter(dset_plist, (H5Z_filter_t)32002,
					H5Z_FLAG_OPTIONAL, 1, &filter_option);
				break;
			case UNKNOWN:
			default:
				// Disable automatic error message to not alarm user...
				H5Eset_auto(H5E_DEFAULT, NULL, NULL);

				// Try to use MAFISC filter, then fallback to DEFLATE
				filter_option = 0;
				if(H5Pset_filter(dset_plist, (H5Z_filter_t)32002,
					H5Z_FLAG_OPTIONAL, 1, &filter_option) < 0) {
					int myrank;
					MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
					// Use our own, more friendly warning message:
					if(myrank == 0)
						fputs("Warning: MAFISC compression filter not found, falling back to GZIP.\n", stderr);
					H5Pset_deflate(dset_plist, GZIP_LEVEL);
					available = GZIP;
				} else {
					available = MAFISC;
				}

				H5Eset_auto(H5E_DEFAULT, (H5E_auto_t) H5Eprint, stderr);
			}
		}
		break;
	case NONE:
		break;
	case GZIP:
		H5Pset_deflate(dset_plist, GZIP_LEVEL);
		break;
	case MAFISC:
		filter_option = 0;
		H5Pset_filter(dset_plist, (H5Z_filter_t)32002,
			H5Z_FLAG_OPTIONAL, 1, &filter_option);
		break;
	case UNKNOWN:
		; // Can't happen.
	}

	wprops->dset = H5Dcreate(ctx->fd, dataset_path, datatype,
			wprops->ds_file, H5P_DEFAULT, dset_plist, H5P_DEFAULT);
	H5Pclose(dset_plist);
}

static void _prop_write_chunk(size_t size, void* buff, struct _prop_writer *wprops)
{
	hsize_t count[2] = {size, wprops->max_size[1]};

	// Set the size to fit this chunk
	wprops->size[0] += size;
	H5Dset_extent(wprops->dset, wprops->size);

	// Set the slice on where to write the chunk
	H5Sset_extent_simple(wprops->ds_file, wprops->rank, wprops->size, wprops->max_size);
	H5Sselect_hyperslab(wprops->ds_file, H5S_SELECT_SET,
			wprops->offset, NULL, count, NULL);
	wprops->offset[0] += size;

	// Resize the dataspace if last chunk
	if(size != chunk_size) {
		H5Sset_extent_simple(wprops->ds_mem, wprops->rank, count, count);
	}

	// Do the writing
	H5Dwrite(wprops->dset, wprops->datatype, wprops->ds_mem,
			wprops->ds_file, H5P_DEFAULT, buff);
}

static void _prop_writer_destroy(struct _prop_writer *wprops)
{
	H5Sclose(wprops->ds_mem);
	H5Sclose(wprops->ds_file);
	H5Dclose(wprops->dset);
}

static void _higio_write_cell_property_hdf5(higio_hdf5 *ctx, const char *dataset_path,
		sim_domain *sd,	distributed_property *const *dp_component,
		size_t num_components, bool use_var_env_compression)
{
	// Prepare the dataspaces (memory accessor and the on-file layout) and dataset
	struct _prop_writer wprops;
	_prop_writer_init_dataset(ctx, dataset_path, dp_component[0]->pdata->local_count,
		num_components, _get_hdf5_real(), use_var_env_compression, &wprops);

	// The buffer to store the ordered data before write
	real buff[chunk_size][num_components];
	size_t i = 0;

	// Build iterator
	higcit_celliterator *iter = sd_get_domain_celliterator(sd);

	mp_mapper *m = sd_get_domain_mapper(sd);
	for(; !higcit_isfinished(iter); higcit_nextcell(iter))
	{
		hig_cell *cell = higcit_getcell(iter);
		int id = mp_lookup(m, hig_get_cid(cell));

		for(size_t j = 0; j < num_components; ++j) {
			buff[i][j] = dp_get_value(dp_component[j], id);
		}
		if(++i == chunk_size) {
			_prop_write_chunk(i, buff, &wprops);
			i = 0;
		}
	}
	higcit_destroy(iter);

	// Finalize write with data in buff
	if(i) {
		_prop_write_chunk(i, buff, &wprops);
	}
	_prop_writer_destroy(&wprops);
}

void higio_write_cell_property_hdf5(higio_hdf5 *ctx, const char *dataset_path,
	sim_domain *sd,	distributed_property *const *dp_component,
	size_t num_components)
{
	_higio_write_cell_property_hdf5(ctx, dataset_path, sd, dp_component,
		num_components, false);
}

/* Used only inside an assert, compiler will warn when building optimized. */
static bool _components_matches(distributed_property * const* dp_component, size_t num_components) {
	for(size_t i = 1; i < num_components; ++i) {
		if(dp_component[0]->pdata->local_count != dp_component[i]->pdata->local_count) {
			return false;
		}
	}
	return true;
}

static int _create_read_prop_dataset(higio_hdf5 *ctx, const char *dataset_path,
		size_t prop_dimension, size_t rank,
		size_t* num_full_chunks, hsize_t *chunk_size,
		hid_t *ds_file, hid_t *dset, size_t *num_elems)
{
	// Open tree dataset
	*dset = H5Dopen(ctx->fd, dataset_path, H5P_DEFAULT);
	*ds_file = H5Dget_space(*dset);

	// Get the chunk size
	{
		hid_t dset_plist = H5Dget_create_plist(*dset);
		H5Pget_chunk(dset_plist, 1, chunk_size);
		H5Pclose(dset_plist);
	}

	// Size and dimension sanity checks
	{
		int read_rank = H5Sget_simple_extent_ndims(*ds_file);
		hsize_t total_size[read_rank];
		H5Sget_simple_extent_dims(*ds_file, total_size, NULL);

		*num_full_chunks = total_size[0] / *chunk_size;

		// Check if dimensions matches
		if(rank != read_rank || (read_rank > 1 && total_size[1] != prop_dimension)) {
			ssize_t fname_size = 1 + H5Fget_name(ctx->fd, NULL, 0);
			char fname[fname_size];
			H5Fget_name(ctx->fd, fname, fname_size);
			fprintf(stderr, "Error: Property stored in %s:%s has rank %d, "
				"but rank %d was expected.\n", fname, dataset_path,
				(int)(read_rank == 1 ? read_rank : total_size[1]),
				(int)prop_dimension);

			H5Sclose(*ds_file);
			H5Dclose(*dset);
			return -1;
		}
		*num_elems = total_size[0];
	}
	return 0;
}

int higio_read_cell_property_hdf5(higio_hdf5* ctx, const char* dataset_path, sim_domain *sd, distributed_property * const* dp_component, size_t num_components)
{
	assert(num_components > 0);

	int rank;
	if(num_components == 1) {
		rank = 1;
	} else {
		rank = 2;

		// Assert all property components have the same size
		assert(_components_matches(dp_component, num_components));
	}

	size_t num_full_chunks, total_size;
	hsize_t chunk_size;
	hid_t ds_file;
	hid_t dset;

	// Open file, dataset, read sizes and dataspace and check for inconsistencies.
	if(_create_read_prop_dataset(ctx, dataset_path, num_components, rank,
		&num_full_chunks, &chunk_size, &ds_file, &dset, &total_size) < 0)
	{
		return -1;
	}

	// Build iterator
	higcit_celliterator *iter = sd_get_domain_celliterator(sd);

	// The buffer to store the ordered data before write
	real buff[chunk_size][num_components];
	size_t i = 0;

	// Stuff to iterate through the leafs of all trees
	hsize_t offset[2] = {0};
	hsize_t count[2] = {chunk_size, num_components};

	hid_t ds_mem;
	ds_mem = H5Screate_simple(rank, count, count);

	// Load the full chunks of data
	mp_mapper *m = sd_get_domain_mapper(sd);
	for(size_t i = 0; i < num_full_chunks; ++i)
	{
		H5Sselect_hyperslab(ds_file, H5S_SELECT_SET, offset, NULL, count, NULL);
		H5Dread(dset, _get_hdf5_real(), ds_mem, ds_file, H5P_DEFAULT, buff);
		offset[0] += chunk_size;

		for(size_t j = 0; j < chunk_size; ++j) {
			hig_cell *cell = higcit_getcell(iter);
			assert(cell);
			higcit_nextcell(iter);

			size_t id = mp_lookup(m, hig_get_cid(cell));
			for(size_t k = 0; k < num_components; ++k) {
				dp_set_value(dp_component[k], id, buff[j][k]);
			}
		}
	}

	// Load the last partial chunk
	size_t last_partial_chunk = total_size % chunk_size;
	if(last_partial_chunk) {
		count[0] = last_partial_chunk;
		H5Sselect_hyperslab(ds_file, H5S_SELECT_SET, offset, NULL, count, NULL);
		H5Sset_extent_simple(ds_mem, rank, count, count);
		H5Dread(dset, _get_hdf5_real(), ds_mem, ds_file, H5P_DEFAULT, buff);

		for(size_t j = 0; j < last_partial_chunk; ++j) {
			hig_cell *cell = higcit_getcell(iter);
			assert(cell);
			higcit_nextcell(iter);

			size_t id = mp_lookup(m, hig_get_cid(cell));
			for(size_t k = 0; k < num_components; ++k) {
				dp_set_value(dp_component[k], id, buff[j][k]);
			}
		}
	}

	H5Sclose(ds_mem);

	assert(higcit_isfinished(iter));
	higcit_destroy(iter);

	H5Sclose(ds_file);
	H5Dclose(dset);

	return 0;
}

void higio_write_facet_property_hdf5(higio_hdf5* ctx, const char* dataset_path,
		sim_facet_domain *sfd, distributed_property * const* dp_component,
		size_t num_components)
{
	// Prepare the dataspaces (memory accessor and the on-file layout) and dataset
	struct _prop_writer wprops;
	_prop_writer_init_dataset(ctx, dataset_path, dp_component[0]->pdata->local_count,
		num_components, _get_hdf5_real(), false, &wprops);

	// The buffer to store the ordered data before write
	real buff[chunk_size][num_components];
	size_t i = 0;

	// Build iterator
	higfit_facetiterator *iter = sfd_get_domain_facetiterator(sfd);

	mp_mapper *m = sfd_get_domain_mapper(sfd);
	for(; !higfit_isfinished(iter); higfit_nextfacet(iter))
	{
		hig_facet *facet = higfit_getfacet(iter);
		int fid = mp_lookup(m, hig_get_fid(facet));

		// Fill one element in the write out buffer.
		for(size_t j = 0; j < num_components; ++j) {
			buff[i][j] = dp_get_value(dp_component[j], fid);
		}

		// Buffer is full, flush to file
		if(++i == chunk_size) {
			_prop_write_chunk(i, buff, &wprops);
			i = 0;
		}
	}
	higfit_destroy(iter);

	// Finalize write with data in buff
	if(i) {
		_prop_write_chunk(i, buff, &wprops);
	}
	_prop_writer_destroy(&wprops);
}

int higio_read_facet_property_hdf5(higio_hdf5* ctx, const char* dataset_path, sim_facet_domain *sfd, distributed_property *const * dp, size_t num_components)
{
	assert(num_components > 0);

	int rank;
	if(num_components == 1) {
		rank = 1;
	} else {
		rank = 2;

		// Assert all property components have the same size
		assert(_components_matches(dp, num_components));
	}

	size_t num_full_chunks, total_size;
	hsize_t chunk_size;
	hid_t ds_file;
	hid_t dset;

	// Open file, dataset, read sizes and dataspace and check for inconsistencies.
	if(_create_read_prop_dataset(ctx, dataset_path, num_components, rank,
		&num_full_chunks, &chunk_size, &ds_file, &dset, &total_size) < 0)
	{
		return -1;
	}

	// Build iterator
	higfit_facetiterator *iter = sfd_get_domain_facetiterator(sfd);

	// The buffer to store the ordered data before write
	real buff[chunk_size][num_components];
	size_t i = 0;

	// Stuff to iterate through the leafs of all trees
	hsize_t offset[2] = {0};
	hsize_t count[2] = {chunk_size, num_components};

	hid_t ds_mem;
	ds_mem = H5Screate_simple(rank, count, count);

	// Load the full chunks of data
	mp_mapper *m = sfd_get_domain_mapper(sfd);
	for(size_t i = 0; i < num_full_chunks; ++i)
	{
		H5Sselect_hyperslab(ds_file, H5S_SELECT_SET, offset, NULL, count, NULL);
		H5Dread(dset, _get_hdf5_real(), ds_mem, ds_file, H5P_DEFAULT, buff);
		offset[0] += chunk_size;

		for(size_t j = 0; j < chunk_size; ++j) {
			assert(!higfit_isfinished(iter));
			hig_facet *facet = higfit_getfacet(iter);
			assert(facet);

			size_t id = mp_lookup(m, hig_get_fid(facet));
			for(size_t k = 0; k < num_components; ++k) {
				dp_set_value(dp[k], id, buff[j][k]);
			}
			higfit_nextfacet(iter);
		}
	}

	// Load the last partial chunk
	size_t last_partial_chunk = total_size % chunk_size;
	if(last_partial_chunk) {
		count[0] = last_partial_chunk;
		H5Sselect_hyperslab(ds_file, H5S_SELECT_SET, offset, NULL, count, NULL);
		H5Sset_extent_simple(ds_mem, rank, count, count);
		H5Dread(dset, _get_hdf5_real(), ds_mem, ds_file, H5P_DEFAULT, buff);

		for(size_t j = 0; j < last_partial_chunk; ++j) {
			assert(!higfit_isfinished(iter));
			hig_facet *facet = higfit_getfacet(iter);
			assert(facet);

			size_t id = mp_lookup(m, hig_get_fid(facet));
			for(size_t k = 0; k < num_components; ++k) {
				dp_set_value(dp[k], id, buff[j][k]);
			}
			higfit_nextfacet(iter);
		}
	}

	H5Sclose(ds_mem);

	assert(higfit_isfinished(iter));
	higfit_destroy(iter);

	H5Sclose(ds_file);
	H5Dclose(dset);

	return 0;
}

/* Starting XDMF write stuff. */

static gboolean _point_equal(gconstpointer a, gconstpointer b)
{
	PPoint p1 = (PPoint)a, p2 = (PPoint)b;

	for(size_t i = 0; i < DIM; ++i) {
		if(p1[i] != p2[i])
			return FALSE;
	}
	return TRUE;
}

static guint _point_hash(gconstpointer key)
{
	guint *p = (guint*)key;
	const size_t n = DIM * sizeof(double) / sizeof(guint);
	guint hash = 0;
	for(size_t i = 0; i < n; ++i) {
		hash ^= p[i] + 0x9e3779b9 + (hash << 6) + (hash >> 2);
	}
	return hash;
}

static void _xdmf_write_h5_grid(higio_hdf5 *ctx, sim_domain *sd,
	unsigned int *num_points, unsigned int *num_cells)
{
	// Corner permutations
	const size_t PERMS = 1u << DIM;

	struct _prop_writer writer_points, writer_cells;
	_prop_writer_init_dataset(ctx, "/grid/points", H5S_UNLIMITED, DIM,
		_get_hdf5_real(), true, &writer_points);
	/* TODO: use portable point index size */
	_prop_writer_init_dataset(ctx, "/grid/cells", H5S_UNLIMITED, PERMS,
		H5T_NATIVE_ULONG, true, &writer_cells);

	*num_points = *num_cells = 0;

	size_t idx_points = 0;
	size_t idx_cells = 0;

	unsigned long (*buff_cells)[PERMS];
	ALLOC_INFER(buff_cells, chunk_size);
	real (*buff_points)[DIM];
	ALLOC_INFER(buff_points, chunk_size);

	PPoint p = malloc(sizeof(Point));
	GHashTable *points;
	points = g_hash_table_new_full(_point_hash, _point_equal, free, NULL);
	gsize points_ids = 0;

	/* The corners are enumerated in the following order:
	 *  2 +--------+ 3
	 *    |        |
	 *    |        |
	 *    |        |
	 *  0 +--------+ 1
	 *
	 *  But XDMF expect them in the following order (3 swapped with 2):
	 *
	 *  3 +--------+ 2
	 *    |        |
	 *    |        |
	 *    |        |
	 *  0 +--------+ 1
	 *
	 * The following index map is used to correct the placing:
	 */
	static const size_t idxmap[] = {0, 1, 3, 2};

	higcit_celliterator *iter;
	for(iter = sd_get_domain_celliterator(sd);
		!higcit_isfinished(iter); higcit_nextcell(iter))
	{
		Point lims[2];
		hig_cell *c = higcit_getcell(iter);
		hig_get_highpoint(c, lims[0]);
		hig_get_lowpoint(c, lims[1]);

		// Find all corner permutations
		for(size_t i = 0; i < PERMS; ++i) {
			for(size_t j = 0; j < DIM; ++j) {
				p[j] = lims[(i>>j) & 1u][j];
			}

			gsize curr_id;

			gboolean found = g_hash_table_lookup_extended(points,
				p, NULL, (gpointer *)&curr_id);
			if(!found) {
				POINT_ASSIGN(buff_points[idx_points], p);
				++idx_points;
				if(idx_points == chunk_size) {
					_prop_write_chunk(idx_points, buff_points, &writer_points);
					*num_points += idx_points;
					idx_points = 0;
				}
				curr_id = points_ids++;
				g_hash_table_insert(points, p, GSIZE_TO_POINTER(curr_id));
				p = malloc(sizeof(Point));
			}
			buff_cells[idx_cells][idxmap[i%4] + 4 * (i/4)] = curr_id;
		}
		++idx_cells;
		if(idx_cells == chunk_size) {
			_prop_write_chunk(idx_cells, buff_cells, &writer_cells);
			*num_cells += idx_cells;
			idx_cells = 0;
		}
	}
	free(p);
	higcit_destroy(iter);
	g_hash_table_destroy(points);

	// Flush buffer
	if(idx_points) {
		_prop_write_chunk(idx_points, buff_points, &writer_points);
		*num_points += idx_points;
	}
	if(idx_cells) {
		_prop_write_chunk(idx_cells, buff_cells, &writer_cells);
		*num_cells += idx_cells;
	}

	free(buff_cells);
	free(buff_points);

	_prop_writer_destroy(&writer_points);
	_prop_writer_destroy(&writer_cells);
}

struct _cell_prop
{
	struct _cell_prop* next;
	char* prop_name;
	size_t num_components;
	distributed_property* dps[];
};

struct _facet_prop_component
{
	sim_facet_domain* sfd;
	distributed_property *dp;
};

struct _facet_prop
{
	struct _facet_prop* next;
	char* prop_name;
	size_t num_components;
	struct _facet_prop_component components[];
};

struct xdmf_output
{
	size_t last_grid_time_index;
	size_t str_buf_size;
	unsigned int local_sizes[2];
	char* file_prefix;
	char* filename_prefix;
	sim_domain* sd;

	struct _cell_prop *cell_props;
	struct _facet_prop *facet_props;
};

xdmf_output *xdmf_init(const char* file_prefix, sim_domain* sd)
{
	// TODO: treat no local file names
	assert((DIM == 2 || DIM == 3) && "XDMF write is implemented only in 'normal' dimensions!");

	xdmf_output *ret = calloc(1, sizeof *ret);
	ret->str_buf_size = strlen(file_prefix) + 100;
	ret->file_prefix = strdup(file_prefix);
	char *last_slash = rindex(ret->file_prefix, '/');
	ret->filename_prefix = last_slash ? last_slash + 1 : ret->file_prefix;
	ret->sd = sd;
	return ret;
}

void xdmf_destroy(xdmf_output *ctx)
{
	free(ctx->file_prefix);

	while(ctx->cell_props != NULL) {
		struct _cell_prop* next = ctx->cell_props->next;
		free(ctx->cell_props->prop_name);
		free(ctx->cell_props);
		ctx->cell_props = next;
	}

	while(ctx->facet_props != NULL) {
		struct _facet_prop* next = ctx->facet_props->next;
		free(ctx->facet_props->prop_name);
		free(ctx->facet_props);
		ctx->facet_props = next;
	}

	free(ctx);
}

const char* xdmf_get_file_prefix(xdmf_output *ctx)
{
	return ctx->file_prefix;
}

void xdmf_register_cell_property(xdmf_output *ctx, const char* prop_name,
	size_t num_components, distributed_property *dp, ...)
{
	assert(num_components >= 1);
	struct _cell_prop *new_prop = malloc(sizeof(*new_prop) + num_components * sizeof(distributed_property *));

	new_prop->num_components = num_components;
	new_prop->prop_name = strdup(prop_name);

	new_prop->dps[0] = dp;

	va_list ap;
	va_start(ap, dp);
	for(size_t i = 1; i < num_components; ++i) {
		new_prop->dps[i] = va_arg(ap, distributed_property *);
	}
	va_end(ap);

	new_prop->next = ctx->cell_props;
	ctx->cell_props = new_prop;
}

static struct _facet_prop *
_xdmf_create_facet_property(xdmf_output *ctx, const char* prop_name,
	size_t num_components)
{
	assert(num_components >= 1);
	struct _facet_prop *new_prop = malloc(sizeof(*new_prop) + num_components * sizeof(struct _facet_prop_component));

	new_prop->num_components = num_components;
	new_prop->prop_name = strdup(prop_name);

	new_prop->next = ctx->facet_props;
	ctx->facet_props = new_prop;

	return new_prop;
}

void xdmf_register_facet_property(xdmf_output *ctx, const char* prop_name,
	size_t num_components, sim_facet_domain *sfd, distributed_property *dp, ...)
{
	struct _facet_prop *new_prop =
		_xdmf_create_facet_property(ctx, prop_name, num_components);

	new_prop->components[0].sfd = sfd;
	new_prop->components[0].dp = dp;

	va_list ap;
	va_start(ap, dp);
	for(size_t i = 1; i < num_components; ++i) {
		new_prop->components[i].sfd = va_arg(ap, sim_facet_domain *);
		new_prop->components[i].dp = va_arg(ap, distributed_property *);
	}
	va_end(ap);
}

void xdmf_register_facet_property_array(xdmf_output *ctx,
	const char* prop_name, size_t num_components,
	sim_facet_domain *sfd[], distributed_property *dp[])
{
	struct _facet_prop *new_prop =
		_xdmf_create_facet_property(ctx, prop_name, num_components);

	for(size_t i = 0; i < num_components; ++i) {
		new_prop->components[i].sfd = sfd[i];
		new_prop->components[i].dp = dp[i];
	}
}

static void _output_xdmf_prop(FILE *xdmf, const char* prop_name,
	const char* h5name, const char* gridh5name, size_t num_components,
	unsigned int grid_size)
{
	// Evil hack for 2D vectors to work on Paraview...
	// Split the 2D vector data into 2 scalar arrays, and reassemble into a
	// 3D vector array, with 3rd arguent being zero.
	if(DIM == 2 && num_components == 2) {
		fprintf(xdmf, ""
"        <Attribute Name=\"%1$s\" AttributeType=\"Vector\" Center=\"Cell\">\n"
"	   <DataItem Dimensions=\"%2$u 3\" Function=\"JOIN($0, $1, $2)\" ItemType=\"Function\">\n"
"		<DataItem ItemType=\"HyperSlab\" Dimensions=\"%2$u\" Type=\"HyperSlab\">\n"
"			<DataItem Dimensions=\"6\" Format=\"XML\">0 0 1 1 %2$u 1</DataItem>\n"
"			<DataItem Dimensions=\"%2$u 2\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\">\n"
"				%3$s:/properties/%1$s\n"
"			</DataItem>\n"
"		</DataItem>\n"
"		<DataItem ItemType=\"HyperSlab\" Dimensions=\"%2$u\" Type=\"HyperSlab\">\n"
"			<DataItem Dimensions=\"6\" Format=\"XML\">0 1 1 1 %2$u 1</DataItem>\n"
"			<DataItem Dimensions=\"%2$u 2\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\">\n"
"				%3$s:/properties/%1$s\n"
"			</DataItem>\n"
"		</DataItem>\n"
"		<DataItem Dimensions=\"%2$u\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\">\n"
"			%4$s:/workaround/zero\n"
"		</DataItem>\n"
"	   </DataItem>\n"
"        </Attribute>\n",
		prop_name, grid_size, h5name, gridh5name);
	} else {
		const char* dtype;
		if(num_components == 1) {
			dtype = "Scalar";
		} else {
			dtype = "Vector";
		}

		fprintf(xdmf, ""
"        <Attribute Name=\"%s\" AttributeType=\"%s\" Center=\"Cell\">\n"
"          <DataItem Format=\"HDF\" Dimensions=\"%u",
			prop_name, dtype, grid_size);
		if(num_components != 1) {
			fprintf(xdmf, " %zu", num_components);
		}
		fprintf(xdmf, "\" NumberType=\"Float\" Precision=\"8\">\n%s:/properties/%s\n"
"          </DataItem>\n"
"        </Attribute>\n",
			h5name, prop_name);
	}
}

void xdmf_write_timestep(xdmf_output* ctx, size_t time_index, bool reuse_last_grid)
{
	higio_hdf5 *hdf5;

	FILE *xdmf = NULL;
	//Number of points, then number of cells for this process
	char fname[ctx->str_buf_size];

	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	if(myrank == 0) {
		snprintf(fname, sizeof fname, "%s.t%05zu.xdmf", ctx->file_prefix, time_index);
		xdmf = fopen(fname, "w");
	}

	snprintf(fname, sizeof fname, "%s.t%05zu.p%05d.h5", ctx->file_prefix, time_index, myrank);
	hdf5 = higio_create_hdf5(fname);

	// Maybe write current grid
	if(!reuse_last_grid) {
		hid_t group = H5Gcreate(hdf5->fd, "grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Gclose(group);
		_xdmf_write_h5_grid(hdf5, ctx->sd, &ctx->local_sizes[0], &ctx->local_sizes[1]);
		ctx->last_grid_time_index = time_index;

		// Workaround for 2D output in Paraview: create a zero-filled dataset
		if(DIM == 2) {
			group = H5Gcreate(hdf5->fd, "workaround", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			H5Gclose(group);

			hsize_t size = ctx->local_sizes[1];
			hid_t ds = H5Screate_simple(1, &size, &size);

			hid_t dset_plist = H5Pcreate(H5P_DATASET_CREATE);
			H5Pset_chunk(dset_plist, 1, &size);

			// The dataset is created with sparse storage, no chunks
			// is allocated, and reading from it returns 0 for every element,
			// what is precisely what we need.
			hid_t dset = H5Dcreate(hdf5->fd, "/workaround/zero", _get_hdf5_real(),
				ds, H5P_DEFAULT, dset_plist, H5P_DEFAULT);

			H5Dclose(dset);
			H5Pclose(dset_plist);
			H5Sclose(ds);
		}
	}

	// Create properties group
	{
		hid_t group = H5Gcreate(hdf5->fd, "properties", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Gclose(group);
	}

	// Write cell centered properties
	for(struct _cell_prop *prop = ctx->cell_props; prop; prop = prop->next)
	{
		char dset_name[15 + strlen(prop->prop_name)];
		snprintf(dset_name, sizeof dset_name, "/properties/%s", prop->prop_name);

		_higio_write_cell_property_hdf5(hdf5, dset_name, ctx->sd,
			prop->dps, prop->num_components, true);
	}

	// Write facet properties
	sim_stencil *stn = stn_create();
	for(struct _facet_prop *prop = ctx->facet_props; prop; prop = prop->next)
	{
		struct _prop_writer pwriter;
		{
			char dset_name[15 + strlen(prop->prop_name)];
			snprintf(dset_name, sizeof dset_name,
				"/properties/%s", prop->prop_name);
			_prop_writer_init_dataset(hdf5, dset_name, ctx->local_sizes[1],
				prop->num_components, _get_hdf5_real(), true, &pwriter);
		}

		real buff[chunk_size][prop->num_components];
		higcit_celliterator *it;
		size_t idx = 0;
		for(it = sd_get_domain_celliterator(ctx->sd);
				!higcit_isfinished(it); higcit_nextcell(it))
		{
			hig_cell *c = higcit_getcell(it);

			Point ccenter;
			Point cdelta;
			hig_get_center(c, ccenter);
			hig_get_delta(c, cdelta);
			POINT_DIV_SCALAR(cdelta, cdelta, 2.0);

			// Interpolate and fill-in value for each component
			for(size_t j = 0; j < prop->num_components; ++j) {
				size_t dim = sfd_get_dim(prop->components[j].sfd);
				Point dcenter;
				POINT_ASSIGN(dcenter, ccenter);
				dcenter[dim] += cdelta[dim];
				buff[idx][j] = sfd_dp_interpolate(prop->components[j].sfd, prop->components[j].dp, dcenter, ccenter, stn);
			}
			if(++idx == chunk_size) {
				_prop_write_chunk(idx, buff, &pwriter);
				idx = 0;
			}
		}
		if(idx)
			_prop_write_chunk(idx, buff, &pwriter);
		higcit_destroy(it);
		_prop_writer_destroy(&pwriter);
	}
	stn_destroy(stn);

	higio_close_hdf5(hdf5);

	if(xdmf) {
		int dsize;
		MPI_Comm_size(MPI_COMM_WORLD, &dsize);
		unsigned int total_sizes[dsize][2];

		MPI_Gather(ctx->local_sizes, 2, MPI_UNSIGNED, total_sizes, 2, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

		fputs(""
		"<?xml version=\"1.0\"?>\n"
		"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
		"<Xdmf Version=\"2.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\">\n"
		"  <Domain>\n"
		"    <Grid Name=\"Domain\" GridType=\"Collection\">\n",
		xdmf);

		const char* topo_type = NULL;
		const char* geom_type = NULL;
		switch(DIM) {
			case 2:
				topo_type = "Quadrilateral";
				geom_type = "XY";
				break;
			case 3:
				topo_type = "Hexahedron";
				geom_type = "XYZ";
				break;
		};

		for(int i = 0; i < dsize; ++i) {
			size_t grid_datac = 0;

			char grid_fname[ctx->str_buf_size];
			snprintf(grid_fname, sizeof fname, "%s.t%05zu.p%05d.h5", ctx->filename_prefix, ctx->last_grid_time_index, i);
			fprintf(xdmf, ""
"      <Grid Name=\"Proc %d\" GridType=\"Uniform\">\n"
"        <Topology NumberOfElements=\"%u\" TopologyType=\"%s\">\n"
"          <DataItem Format=\"HDF\" NumberType=\"UInt\" Precision=\"8\" Dimensions=\"%u %u\">\n"
"%s:/grid/cells\n"
"          </DataItem>\n"
"        </Topology>\n"
"        <Geometry GeometryType=\"%s\">\n"
"          <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\"%u %u\">\n"
"%s:/grid/points\n"
"          </DataItem>\n"
"        </Geometry>\n",
				i, total_sizes[i][1], topo_type,
				total_sizes[i][1], 1u << DIM, grid_fname,
				geom_type,
				total_sizes[i][0], DIM, grid_fname);

			snprintf(fname, sizeof fname, "%s.t%05zu.p%05d.h5", ctx->filename_prefix, time_index, i);

			for(struct _cell_prop *prop = ctx->cell_props; prop; prop = prop->next)
				_output_xdmf_prop(xdmf, prop->prop_name, fname, grid_fname, prop->num_components, total_sizes[i][1]);

			for(struct _facet_prop *prop = ctx->facet_props; prop; prop = prop->next)
				_output_xdmf_prop(xdmf, prop->prop_name, fname, grid_fname, prop->num_components, total_sizes[i][1]);

			fputs("      </Grid>\n", xdmf);
		}

		fputs(""
		"    </Grid>\n"
		"  </Domain>\n"
		"</Xdmf>\n",
		xdmf);

		fclose(xdmf);
	} else {
		MPI_Gather(ctx->local_sizes, 2, MPI_UNSIGNED, NULL, 2, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	}
}
