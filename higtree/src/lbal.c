#include <assert.h>
#include <zoltan.h>
#include <string.h>
#include "allocator.h"
#include "isend-pool.h"
#include "rect.h"
#include "build-fringe.h"
#include "term-det.h"
#include "higtree-parallel.h"
#include "higtree-serialize.h"
#include "Debug-c.h"

#include "lbal.h"

struct _part_params {
	struct load_balancer *ctx;
	struct _lb_input_tree **tree_index;
};

load_balancer * lb_create(MPI_Comm group, unsigned num_groups)
{
	assert(num_groups >= 1);
	DECL_AND_ALLOC(load_balancer, ret, 1);

	ret->num_groups = num_groups;
	ALLOC_INFER(ret->group_weight, num_groups);
	ret->group_weight[0] = 1.0;
	for(size_t i = 1; i < num_groups; ++i) {
		ret->group_weight[i] = 0.0001;
	}

	// By default, all tree groups will generate fringe.
	ALLOC_INFER(ret->group_fringe, num_groups);
	for(size_t i = 0; i < num_groups; ++i) {
		ret->group_fringe[i] = true;
	}

	ret->group = group;
	ret->trees = NULL;
	ret->output = NULL;
	ret->portals = NULL;
	ret->out_higs_group_start = NULL;
	ret->num_out_higs = 0;
	ret->first_tree_id = 0;
	ret->num_in_higs = 0;
	ret->num_portals = 0;
	return ret;
}

void lb_set_group_weight(load_balancer *ctx, unsigned group, float weight)
{
	assert(group < ctx->num_groups);
	ctx->group_weight[group] = weight;
}

void lb_set_group_has_fringe(load_balancer *ctx, unsigned group, bool has_fringe)
{
	assert(group < ctx->num_groups);
	ctx->group_fringe[group] = has_fringe;
}

void lb_add_input_tree(load_balancer *ctx, hig_cell *tree, bool managed,
		unsigned group)
{
	DECL_AND_ALLOC(struct _lb_input_tree, new, 1);
	new->tree = tree;
	new->managed = managed;
	new->next = ctx->trees;
	assert(group < ctx->num_groups);
	new->group = group;
	ctx->trees = new;
	++ctx->num_in_higs;
}

void lb_add_portal(load_balancer *lb, const Rect *a, const Rect *b,
		unsigned group)
{
	assert(rect_degenerate_count(a) == 1);
	assert(rect_degenerate_count(b) == 1);
	assert(rect_has_same_size(a, b));

	DECL_AND_ALLOC(struct _lb_portal, new, 1);
	new->r[0] = *a;
	new->r[1] = *b;
	new->group = group;
	new->next = lb->portals;
	lb->portals = new;

	++lb->num_portals;
}

static int _zt_cb_num_obj(size_t *total_elems, int *ierr)
{
	return *total_elems;
}

static void _zt_cb_obj_list(load_balancer *ctx, int num_gid_entries, int num_lid_entries,
	ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int wgt_dim, float *obj_wgts, int *ierr)
{
	assert(num_gid_entries == 2);

	// Iterate through all elements and calculate weight
	size_t idx = 0;
	size_t tree_id = ctx->first_tree_id;
	for(struct _lb_input_tree *ptr = ctx->trees; ptr; ptr = ptr->next)
	{
		size_t num_children = hig_get_number_of_children(ptr->tree);
		for(size_t i = 0; i < num_children; ++i) {
			// Fill in weight
                	hig_cell *b = hig_get_child(ptr->tree, i);
			higcit_celliterator *it = higcit_create_all_leaves(b);
        		obj_wgts[idx] = higcit_count(it) * ctx->group_weight[ptr->group];
                	higcit_destroy(it);

			// Fill in ID
			global_ids[2*idx] = tree_id;
			global_ids[2*idx + 1] = i;

			++idx;
		}
		++tree_id;
	}
}

static int _zt_cb_num_geom(void *ignore, int *ierr)
{
	return DIM;
}

/* Converge numbers close to each other to the same value. */
static real _normalize_float(real n)
{
	return EPSDELTA * round(n / EPSDELTA);
}

static void _zt_cb_geom_multi(struct _part_params *params, int num_gid_entries, int num_lid_entries,
	int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int num_dim, double *geom_vec,
	int *ierr)
{
	assert(num_gid_entries == 2);
	assert(num_dim == DIM);

	// Iterate through all elementes and calculate center
	for(unsigned i = 0; i < num_obj; ++i)
	{
		unsigned tree = global_ids[2 * i] - params->ctx->first_tree_id;
		unsigned elem = global_ids[2 * i + 1];

		hig_cell *c = hig_get_child(params->tree_index[tree]->tree, elem);

		// Fill in center
		{
			Point p;
			hig_get_center(c, p);
			for(size_t j = 0; j < DIM; ++j) {
				// Not the exact center value, but a value within
				// acceptable error (EPSDELTA) difference
				geom_vec[j] = _normalize_float(p[j]);
			}
		}

		geom_vec += DIM;
	}
}

static int _id_cmp(const unsigned int *a, const unsigned int *b)
{
	// Assuming that converting unsigned to int will keep the same binary
	// representation, and 2-complement is used...
	unsigned int ret = a[0] - b[0];
	if(ret)
		return ret;
	return a[1] - b[1];
}

struct _to_export_elem
{
	int proc; // proc must be the first element to be used by _create_proc_list.
	unsigned int fid[2];
};

static int _proc_id_cmp(struct _to_export_elem *a, struct _to_export_elem *b)
{
	int val = a->proc - b->proc;
	if(val)
		return val;

	return _id_cmp(a->fid, b->fid);
}

static int _just_id_cmp(struct _to_export_elem *a, struct _to_export_elem *b)
{
	return _id_cmp(a->fid, b->fid);
}

static int _out_group_cmp(struct _lb_output_tree *a, struct _lb_output_tree *b)
{
	return a->group - b->group;
}

static unsigned _create_proc_list(void *procs, unsigned count, size_t size, size_t total_procs, int **ret_procs)
{
	uint8_t *ptr = procs;

	int last_p = -1;
	unsigned num_procs = 0;

	// Create overestimated output vector
	size_t size_est = count < total_procs ? count : total_procs;
	if(ret_procs) {
		ALLOC(int, (*ret_procs), size_est);
	}

	// Count how many processes we have.
	for(unsigned i = 0; i < count; ++i) {
		int proc = *(int*)(ptr + i * size);
		if(proc != last_p) {
			if(ret_procs)
				(*ret_procs)[num_procs] = proc;
			last_p = proc;
			++num_procs;
		}
	}

	// Free unused part of the vector
	if(ret_procs)
		*ret_procs = realloc(*ret_procs, num_procs * sizeof(int));

	return num_procs;
}

// Alertas de nao usados !!!
/*static unsigned 
_calc_chunk_numcells(hig_cell *original_tree, int initial_node_idx, int final_node_idx, unsigned *out_numcells)
{
	// Find the dimensional size of the new tree
	int new_low_pp[DIM];
	hig_tobase(initial_node_idx, original_tree->numcells, new_low_pp);

	int new_high_pp[DIM];
	hig_tobase(final_node_idx, original_tree->numcells, new_high_pp);

	unsigned total_size = 1;
	for(unsigned i = 0; i < DIM; ++i) {
		out_numcells[i] = new_high_pp[i] - new_low_pp[i] + 1;
		total_size *= out_numcells[i];
	}

	return total_size;
}
*/

static void _calc_boundingbox_from_children(hig_cell *tree)
{
	unsigned numoch = hig_get_number_of_children(tree);

	if(!numoch) {
		memset(tree->lowpoint, 0, sizeof tree->lowpoint);
		memset(tree->highpoint, 0, sizeof tree->highpoint);
		return;
	}

	hig_cell *child = hig_get_child(tree, 0);
	hig_get_highpoint(child, tree->highpoint);
	hig_get_lowpoint(child, tree->lowpoint);

	for(unsigned i = 1; i < numoch; ++i) {
		child = hig_get_child(tree, i);
		Point lowpoint, highpoint;

		hig_get_highpoint(child, highpoint);
		hig_get_lowpoint(child, lowpoint);

		for(unsigned d = 0; d < DIM; ++d) {
			if(highpoint[d] > tree->highpoint[d])
				tree->highpoint[d] = highpoint[d];
			if(lowpoint[d] < tree->lowpoint[d])
				tree->lowpoint[d] = lowpoint[d];
		}
	}
}

static bool _try_merge_chunk_into_tree(Point lowpoint, Point highpoint, hig_cell **rcells,
	unsigned *numcells, hig_cell *tree)
{
	int matchdim = -1;

	// For each dimension, try to find one where this chunk matches with the tree
	for(unsigned dim = 0; dim < DIM; ++dim)
	{
		if(POS_EQ(tree->lowpoint[dim], highpoint[dim])
			|| POS_EQ(tree->highpoint[dim], lowpoint[dim]))
		{
			// Found a candidate, assert that all other sizes and dimensions matches
			for(unsigned odim = 0; odim < DIM; ++odim) {
				if(dim == odim)
					continue;

				if(POS_NE(highpoint[odim], tree->highpoint[odim])
					|| POS_NE(lowpoint[odim], tree->lowpoint[odim])
					|| numcells[odim] != tree->numcells[odim])
					return false;
			}
			matchdim = dim;

			// Inneficient goto substitute:
			break;
		}
	}

	// Inneficient goto label substitute:
	if(matchdim == -1)
		return false;

	// Calculate the new dimension
	unsigned newncells[DIM];
	memcpy(newncells, numcells, sizeof newncells);
	newncells[matchdim] += tree->numcells[matchdim];

	unsigned total_size = 1;
	for(unsigned dim = 0; dim < DIM; ++dim) {
		total_size *= newncells[dim];
	}

	// Order the chunks...
	hig_cell **lower;
	int *lower_numcells;

	hig_cell **higher;
	int *higher_numcells;

	unsigned middle;
	if(POS_EQ(tree->lowpoint[matchdim], highpoint[matchdim])) {
		// The tree comes after the new chunk
		lower = rcells;
		lower_numcells = (int*)numcells;
		higher = tree->children;
		higher_numcells = tree->numcells;
		middle = numcells[matchdim];

		// Adjust tree's lowpoint
		tree->lowpoint[matchdim] = lowpoint[matchdim];
	} else {
		// The new chunk comes after the tree
		lower = tree->children;
		lower_numcells = tree->numcells;
		higher = rcells;
		higher_numcells = (int*)numcells;
		middle = tree->numcells[matchdim];

		// Adjust tree's highpoint
		tree->highpoint[matchdim] = highpoint[matchdim];
	}

	// Merge the nodes
	DECL_AND_ALLOC(hig_cell *, children, total_size);
	int pp[DIM] = {0};
	for(unsigned i = 0; i < total_size; ++i) {
		// Find the destination idx
		int idx = hig_toposition((int *)newncells, pp);

		// Assign the element to its new position
		if(pp[matchdim] < middle) {
			int s_idx = hig_toposition(lower_numcells, pp);
			children[idx] = lower[s_idx];
		} else {
			int npp[DIM];
			memcpy(npp, pp, sizeof pp);
			npp[matchdim] -= middle;

			int s_idx = hig_toposition(higher_numcells, npp);
			children[idx] = higher[s_idx];
		}

		// Adjust children to tree
		children[idx]->posinparent = idx;
		children[idx]->parent = tree;

		// Increment element's position
		pp[0]++;
		for(unsigned dim = 0; dim < DIM && pp[dim] == newncells[dim];) {
			pp[dim] = 0;
			++pp[++dim];
		}
	}

	// Adjust the tree's root node
	free(tree->children);
	tree->children = children;

	tree->numcells[matchdim] = newncells[matchdim];

	return true;
}

// If merged, ob will be merged into oa
static bool
_try_merge_trees(struct _lb_output_tree *oa, struct _lb_output_tree *ob)
{
	if(oa->group != ob->group)
		return false;

	hig_cell *a = oa->tree;
	hig_cell *b = ob->tree;
	bool merged = _try_merge_chunk_into_tree(b->lowpoint, b->highpoint, b->children,
		(unsigned *)b->numcells, a);

	if(merged) {
		free(b);
		return true;
	}

	return false;
}

static struct _lb_output_tree *
_try_merge_tree_list(
		size_t num_merged, struct _lb_output_tree *merged,
		size_t num_unmerged, struct _lb_output_tree *unmerged,
		size_t *ret_count
	)
{
	const size_t maximum_count = num_merged + num_unmerged;
	DECL_AND_ALLOC(struct _lb_output_tree, next_merged, maximum_count);
	size_t next_merged_count = 0;

	size_t next_unmerged_count = 0;

	// Try to match each merged with each unmerged
	for(size_t i = 0; i < num_merged; ++i) {
		for(size_t j = next_unmerged_count; j < num_unmerged; ++j) {
			bool merge_res = _try_merge_trees(&merged[i], &unmerged[j]);
			if(merge_res) {
				// Put all the newly created trees at the begining of the
				// unmerged vector, swapping places the just merged tree
				unmerged[j] = unmerged[next_unmerged_count];
				unmerged[next_unmerged_count] = merged[i];
				++next_unmerged_count;

				// Fill the hole with the last element in vector
				merged[i] = merged[--num_merged];

				break;
			}
		}
	}

	// Try to match each unmerged with itself
	while(next_unmerged_count < num_unmerged) {
		bool has_merged = false;
		for(size_t j = next_unmerged_count+1; j < num_unmerged; ++j) {
			has_merged = _try_merge_trees(&unmerged[next_unmerged_count],
					&unmerged[j]);
			if(has_merged) {
				++next_unmerged_count;
				unmerged[j] = unmerged[--num_unmerged];

				//goto continue_while;
				break;
			}
		}

		// Poor goto replacement
		if(!has_merged) {
			// The tree has been tried agains all input trees
			next_merged[next_merged_count++] = unmerged[next_unmerged_count];
			unmerged[next_unmerged_count] = unmerged[--num_unmerged];
		}

		//continue_while:;
	}

	// Put all merged in a single list
	memcpy(&next_merged[next_merged_count], merged, num_merged * sizeof *merged);
	next_merged_count += num_merged;
	next_merged = realloc(next_merged, next_merged_count * sizeof *next_merged);

	// Test if recursion end
	if(num_unmerged) {
		struct _lb_output_tree *ret = _try_merge_tree_list(next_merged_count, next_merged,
			next_unmerged_count, unmerged, ret_count);
		free(next_merged);
		return ret;
	} else {
		*ret_count = next_merged_count;
		return next_merged;
	}
}

static void
_parallel_partition_and_distribute_trees(load_balancer *ctx)
{
	int group_size, rank;
	MPI_Comm_rank(ctx->group, &rank);
	MPI_Comm_size(ctx->group, &group_size);

	// Uniquely enumerate all trees taking part in the partitioning
	{
		MPI_Request reqs[group_size - rank - 1];

		for(int i = rank+1; i < group_size; ++i) {
			MPI_Isend(&ctx->num_in_higs, 1, MPI_UNSIGNED, i, 73, ctx->group, &reqs[i - rank - 1]);
		}
		for(int i = 0; i < rank; ++i) {
			unsigned int received;
			MPI_Recv(&received, 1, MPI_UNSIGNED, i, 73, ctx->group, MPI_STATUS_IGNORE);
			ctx->first_tree_id += received;
		}
		MPI_Waitall(group_size - rank - 1, reqs, MPI_STATUSES_IGNORE);
	}

	// Calculate the number of elements to be distributed and
	// build an array from the linked input list, for random access.
	struct _part_params params = {.ctx = ctx};
	ALLOC_INFER(params.tree_index, ctx->num_in_higs);
	size_t total_elems = 0;
	{
		size_t i = 0;
		for(struct _lb_input_tree *ptr = ctx->trees; ptr; ptr = ptr->next)
		{
			size_t num_children = hig_get_number_of_children(ptr->tree);
			total_elems += num_children;

			params.tree_index[i++] = ptr;
		}
		assert(i == ctx->num_in_higs);
	}

	// Create Zoltan context and pass the data input functions
	struct Zoltan_Struct *zz;
	zz = Zoltan_Create(ctx->group);

	/* Register query functions. */
	Zoltan_Set_Fn(zz, ZOLTAN_NUM_OBJ_FN_TYPE, (void (*)()) _zt_cb_num_obj, (void *)&total_elems);
	Zoltan_Set_Fn(zz, ZOLTAN_OBJ_LIST_FN_TYPE, (void (*)()) _zt_cb_obj_list, (void *)ctx);
	Zoltan_Set_Fn(zz, ZOLTAN_NUM_GEOM_FN_TYPE, (void (*)()) _zt_cb_num_geom, NULL);
	Zoltan_Set_Fn(zz, ZOLTAN_GEOM_MULTI_FN_TYPE, (void (*)()) _zt_cb_geom_multi, (void *)&params);

	/* Set some Zoltan parameters. */
	Zoltan_Set_Param(zz, "debug_level", "0");
	Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "0");
	Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "2");
	Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1");
	Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
	Zoltan_Set_Param(zz, "LB_APPROACH", "REPARTITION");
	Zoltan_Set_Param(zz, "REMAP", "0");
	Zoltan_Set_Param(zz, "RCB_RECTILINEAR_BLOCKS", "1");

	// Make the partition
	int is_new, num_imp, num_exp, *imp_procs, *exp_procs;
	int *imp_to_part, *exp_to_part;
	int num_gid_entries, num_lid_entries;
	ZOLTAN_ID_PTR imp_global_ids, exp_global_ids;
	ZOLTAN_ID_PTR imp_local_ids, exp_local_ids;

	/* Perform the partitioning */
	Zoltan_LB_Partition(zz, &is_new, &num_gid_entries, &num_lid_entries,
		&num_imp,&imp_global_ids,&imp_local_ids,&imp_procs,&imp_to_part,
		&num_exp,&exp_global_ids,&exp_local_ids,&exp_procs,&exp_to_part);

	if(!is_new) {
		Zoltan_LB_Free_Part(&exp_global_ids, &exp_local_ids, &exp_procs, &exp_to_part);
		Zoltan_LB_Free_Part(&imp_global_ids, &imp_local_ids, &imp_procs, &imp_to_part);
		free(params.tree_index);

		// The partitioning was unaltered.
		// Just use the same input as output (or clone it, if not managed)
		ctx->num_out_higs = ctx->num_in_higs;
		ALLOC_INFER(ctx->output, ctx->num_out_higs);
		size_t i = 0;

		while(ctx->trees)
		{
			struct _lb_input_tree *ptr = ctx->trees;

			if(ptr->managed) {
				ctx->output[i].tree = ptr->tree;
			} else {
				ctx->output[i].tree = hig_clone(ptr->tree);
			}
			ctx->output[i].group = ptr->group;
			++i;

			ctx->trees = ptr->next;
			free(ptr);
		}
	} else {
		// Group equal imp process together:
		qsort(imp_procs, num_imp, sizeof *imp_procs, (int (*)(const void *, const void *))int_cmp);

		// Count how many imp processes
		unsigned num_imp_procs = _create_proc_list(imp_procs, num_imp, sizeof *imp_procs,
			group_size, NULL);

		Zoltan_LB_Free_Part(&imp_global_ids, &imp_local_ids, &imp_procs, &imp_to_part);

		// Assemble ids for export with its respective destination processes.
		DECL_AND_ALLOC(struct _to_export_elem, to_export, num_exp);
		for(size_t i = 0; i < num_exp; ++i) {
			{
				unsigned noc = hig_get_number_of_children(params.tree_index[exp_global_ids[2*i] - ctx->first_tree_id]->tree);
				unsigned id = exp_global_ids[2*i + 1];
				assert(id <= noc);;
			}

			// Copy the tree id and the child index:
			memcpy(&to_export[i].fid, &exp_global_ids[2*i], 2*sizeof(unsigned int));
			to_export[i].proc = exp_procs[i];
		}
		Zoltan_LB_Free_Part(&exp_global_ids, &exp_local_ids, &exp_procs, &exp_to_part);

		// Place the returned lists in a defined order (Zoltan documentation
		// doesn't guarantees any particular order).
		qsort(to_export, num_exp, sizeof *to_export, (int (*)(const void *, const void *))_proc_id_cmp);

		// Get a list of the exp processes
		int *all_exp_procs;
		unsigned num_exp_procs = _create_proc_list(to_export, num_exp, sizeof *to_export,
			group_size, &all_exp_procs);

		// Count how many trees and tree nodes will be exported to each process
		size_t num_nodes_to_send[num_exp_procs];
		size_t num_trees_to_send[num_exp_procs];
		memset(num_nodes_to_send, 0, sizeof num_nodes_to_send);
		memset(num_trees_to_send, 0, sizeof num_trees_to_send);
		{
			size_t trees_counted = 0;
			int last_tree = -1;
			int last_p = -1;
			int pidx = -1;
			for(size_t i = 0; i < num_exp; ++i) {
				size_t hig_idx = to_export[i].fid[0] - ctx->first_tree_id;
				int p = to_export[i].proc;

				if(hig_idx != last_tree || p != last_p) {
					if(p != last_p) {
						++pidx;
						last_p = p;
					}
					++num_trees_to_send[pidx];
					last_tree = hig_idx;
				}

				hig_cell *ptree = params.tree_index[hig_idx]->tree;
				hig_cell *elem = hig_get_child(ptree, to_export[i].fid[1]);
				higcit_celliterator *it = higcit_create_all_higtree(elem);
				num_nodes_to_send[pidx] += higcit_count(it);

				higcit_destroy(it);
			}
		}

		// Allocate the preliminary output vector.
		ALLOC_INFER(ctx->output, ctx->num_in_higs);

		// For each process we must send to, the following code will:
		// a) count and send the number of elements per tree the remote process
		//    will receive;
		// b) serialize and send each individual element.
		// c) store the rectangular chunks of nodes removed from this tree
		hig_serial_tree *tree_buffs[num_exp_procs];
		unsigned *sizes_buffs[num_exp_procs];
		unsigned *groups_buffs[num_exp_procs];
		MPI_Request send_reqs[3 * num_exp_procs];
		{
			unsigned k = 0; // Elems iterator.

			// Fill in and send the buffers
			MPI_Request *size_send_reqs = &send_reqs[0];
			MPI_Request *tree_send_reqs = &send_reqs[num_exp_procs];
			MPI_Request *groups_send_reqs = &send_reqs[2 * num_exp_procs];

			for(unsigned pidx = 0; pidx < num_exp_procs; ++pidx) {
				// Allocate the buffers to hold the data to send to pidx:
				DECL_AND_ALLOC(hig_serial_tree, buf, num_nodes_to_send[pidx]);
				tree_buffs[pidx] = buf;

				ALLOC_INFER(sizes_buffs[pidx], num_trees_to_send[pidx] * DIM);
				memset(sizes_buffs[pidx], 0, num_trees_to_send[pidx] * DIM * sizeof *sizes_buffs[pidx]);
				ALLOC_INFER(groups_buffs[pidx], num_trees_to_send[pidx]);

				// Serialize each tree node outgoing to pidx
				int last_tree = -1;
				int first_tree_child = -1;
				int p = all_exp_procs[pidx];
				unsigned otree_it = 0;
				unsigned tree_child_count;
				while(buf - tree_buffs[pidx] < num_nodes_to_send[pidx]) {
					// Get the current element to export
					struct _to_export_elem* to_exp = &to_export[k++];
					assert(to_exp->proc == p);

					// Get its parent tree
					size_t ptree_idx = to_exp->fid[0] - ctx->first_tree_id;
					hig_cell *ptree = params.tree_index[ptree_idx]->tree;

					// Get the subtree from the parent
					hig_cell *elem = hig_get_child(ptree, to_exp->fid[1]);

					// Serialize it
					buf += hs_serialize(buf, elem);
					assert(num_nodes_to_send[pidx] >= buf - tree_buffs[pidx]);

					// For each different tree chunk, increment the counter and
					// find its bounding box
					if(ptree_idx != last_tree) {
						if(first_tree_child != -1) {
                     // Alertas de nao usados !!!
							/*unsigned count;
							count = _calc_chunk_numcells(
								params.tree_index[last_tree]->tree,
								first_tree_child, to_export[k-2].fid[1],
								&sizes_buffs[pidx][DIM * otree_it++]);
							assert(count == tree_child_count);*/
						}
						tree_child_count = 1;
						first_tree_child = to_exp->fid[1];

						// Fill-in the group of this tree:
						groups_buffs[pidx][otree_it] =
							params.tree_index[ptree_idx]->group;

						last_tree = ptree_idx;
					} else {
						++tree_child_count;
					}
				}
				if(first_tree_child != -1) {
					size_t last_tree_child = to_export[k-1].fid[1];
					last_tree = to_export[k-1].fid[0] - ctx->first_tree_id;
               // Alertas de nao usados !!!
					/*unsigned count;
					count = _calc_chunk_numcells(params.tree_index[last_tree]->tree,
						first_tree_child, last_tree_child,
						&sizes_buffs[pidx][DIM * otree_it++]);
					assert(count == tree_child_count);*/
				}
				assert(otree_it == num_trees_to_send[pidx]);

				// Send the actual tree nodes
				MPI_Isend(tree_buffs[pidx], num_nodes_to_send[pidx] * sizeof *buf,
					MPI_BYTE, p, 465, ctx->group, &tree_send_reqs[pidx]);

				// Send the group and the size of each tree
				MPI_Isend(sizes_buffs[pidx], num_trees_to_send[pidx] * DIM,
					MPI_UNSIGNED, p, 546, ctx->group, &size_send_reqs[pidx]);
				MPI_Isend(groups_buffs[pidx], num_trees_to_send[pidx], MPI_UNSIGNED,
					p, 654, ctx->group, &groups_send_reqs[pidx]);

				assert(num_nodes_to_send[pidx] == buf - tree_buffs[pidx]);
			}
			free(all_exp_procs);

			// We now iterate through input trees and build the output trees vector
			// for this process, removing sent elements, and taking into consideration
			// if trees must be reused or copyied.

			// Reorder the elements to being exported, this time, disregarding
			// destination process
			qsort(to_export, num_exp, sizeof *to_export,
				(int (*)(const void *, const void *))_just_id_cmp);

			k = 0; // Reseting iterator

			// Handle each input tree
			for(unsigned i = 0; i < ctx->num_in_higs; ++i)
			{
				hig_cell *tree = params.tree_index[i]->tree;

				// Get the number of original elements
				unsigned numoch = hig_get_number_of_children(tree);
				bool managed = params.tree_index[i]->managed;
				hig_cell **nodes;
				if(managed)
					nodes = tree->children;
				else {
					ALLOC(hig_cell *, nodes, numoch);
				}

				// Iterate through each exported element, checking each against
				// the corresponding node in tree, to remove the exported ones
				struct _to_export_elem* to_exp = &to_export[k];
				unsigned new_count = 0;
				for(unsigned j = 0; j < numoch; ++j) {
					assert(k <= num_exp);

					// If we are not over the exported elements,
					// then we aren't ahead of its tree
					assert(k == num_exp
						|| to_exp->fid[0] - ctx->first_tree_id >= i);

					// It we are not over the exported elements,
					// and if the exported element is in this tree,
					// then j hasn't passed it yet:
					assert(k == num_exp
						|| to_exp->fid[0] - ctx->first_tree_id != i
						|| j <= to_exp->fid[1]);

					if(k >= num_exp
						|| to_exp->fid[0] - ctx->first_tree_id != i
						|| j != to_exp->fid[1])
					{
						// Node remains in tree...
						hig_cell *node = hig_get_child(tree, j);

						if(!managed) {
							hig_cell *tmp = hig_clone(node);
							tmp->posinparent = node->posinparent;
							node = tmp;
						}

						nodes[new_count++] = node;
					} else {
						// Node is removed from the tree...
						if(managed) {
							hig_cell *node = hig_get_child(tree, j);
							hig_destroy(node);
						}

						// Advance exported node iterator
						to_exp = &to_export[++k];
					}
				}
				nodes = realloc(nodes, new_count * sizeof *nodes);

				// Adapt the actual hig_cell structure to fit the changes
				if(new_count == 0) {
					// No elements left, no tree should be output for this one
					if(managed)
						// Can't use hig_destroy here because the tree->children
						// is now invalid, it has been freed by the realloc
						// just above, and all children have been destroyed in
						// the loop.
						free(tree);
					else
						// If not managed and has no tree to output, do nothing.
						;
				} else {
					// There are some remaining elements on this tree,
					// build the output tree.
					hig_cell *out;
					if(managed) {
						out = tree;
						// The pointer has been realloc'ed and might
						// have changed
						tree->children = nodes;
					} else {
						Point p = {0};
						out = hig_create_root(p, p);
						out->children = nodes;
						for(unsigned i = 0; i < new_count; ++i)
						{
							nodes[i]->parent = out;
						}
					}

					// Find the dimensional size of the new tree,
					// relies on the information of the old tree
               // Alertas de nao usados !!!
					/*unsigned calculated_size;
					calculated_size = _calc_chunk_numcells(tree, nodes[0]->posinparent, nodes[new_count - 1]->posinparent, (unsigned int*)out->numcells);
					assert(calculated_size == new_count);*/

					// Adjust posinparent for children:
					for(unsigned i = 0; i < new_count; ++i) {
						nodes[i]->posinparent = i;
					}

					// Adjust new bounding box
					_calc_boundingbox_from_children(out);

					// Add the new tree in the output list
					ctx->output[ctx->num_out_higs].group = params.tree_index[i]->group;
					ctx->output[ctx->num_out_higs].tree  = out;

					ctx->num_out_higs++;
				}
			} // end of "for each input tree"
			free(to_export);
		}

		// Try to coalesce the trees already set to output
		size_t output_size;
		{
			struct _lb_output_tree *new_output;
			new_output = _try_merge_tree_list(0, NULL, ctx->num_out_higs, ctx->output,
				&output_size);
			free(ctx->output);
			ctx->output = new_output;
			ctx->num_out_higs = output_size;
		}

		// Input tree data is no longer needed from now on
		free(params.tree_index);
		ctx->num_in_higs = 0;
		while(ctx->trees) {
			void *tmp = ctx->trees;
			ctx->trees = ctx->trees->next;
			free(tmp);
		}

		// Creates a list to be filled with trees that may still be merged.
		// Use an optimistically allocated size.
		size_t mergeable_size = 2 * output_size + 1;
		size_t mergeable_count = 0;
		DECL_AND_ALLOC(struct _lb_output_tree, mergeable, mergeable_size);

		// Receive the elements from the remote processes
		size_t recv_elem_buff_size = 0;
		size_t recv_sizes_buff_size = 0;
		unsigned *sizes_buf = NULL;
		unsigned *groups_buf = NULL;
		hig_serial_tree *elem_buf = NULL;
		for(int pidx = 0; pidx < num_imp_procs; ++pidx) {
			// Receive the size of the tree_buf
			MPI_Status status;
			MPI_Probe(MPI_ANY_SOURCE, 465, ctx->group, &status);
			int elem_count = 0;
			MPI_Get_count(&status, MPI_BYTE, &elem_count);

			// Adjust tree_buf size
			if(elem_count > recv_elem_buff_size) {
				recv_elem_buff_size = elem_count;

				// won't use realloc here, because I am not
				// interested in keeping old data in the buffer
				free(elem_buf);
				elem_buf = malloc(elem_count);
			}

			// Now we asynchronously receive the trees' elements,
			// while synchronously receiving trees sizes.

			// The receiving of the trees' nodes:
			MPI_Request elems_req;
			MPI_Irecv(elem_buf, elem_count, MPI_BYTE, status.MPI_SOURCE, 465,
				ctx->group, &elems_req);

			// Adjust the tree_count: from bytes to number of elements:
			assert((elem_count % sizeof *elem_buf) == 0);
			elem_count /= sizeof *elem_buf;

			// Sync probe for the message of number of elements per tree:
			MPI_Probe(status.MPI_SOURCE, 546, ctx->group, &status);
			int sizes_count;
			MPI_Get_count(&status, MPI_UNSIGNED, &sizes_count);

			// Adjust sizes_buf to fit the data:
			if(sizes_count > recv_sizes_buff_size) {
				recv_sizes_buff_size = sizes_count;

				free(sizes_buf);
				ALLOC_INFER(sizes_buf, sizes_count);

				free(groups_buf);
				ALLOC_INFER(groups_buf, sizes_count / DIM);
			}

			// Async receive the groups each incoming tree belongs to:
			MPI_Request groups_req;
			MPI_Irecv(groups_buf, sizes_count / DIM, MPI_UNSIGNED,
				status.MPI_SOURCE, 654, ctx->group, &groups_req);

			// Sync receive the tree chunk dimmensions.
			MPI_Recv(sizes_buf, sizes_count, MPI_UNSIGNED, status.MPI_SOURCE,
				546, ctx->group, MPI_STATUS_IGNORE);

			// Adjust sizes_count:
			assert(sizes_count % DIM == 0);
			sizes_count /= DIM;

			MPI_Wait(&elems_req, MPI_STATUS_IGNORE);
			MPI_Wait(&groups_req, MPI_STATUS_IGNORE);

			// Create a buffer to hold temporarily the deserialized
			// tree nodes.
			hig_cell **rcells = NULL;
			size_t rcells_size = 0;

			// Deserialize one tree at a time
			hig_serial_tree *elemptr = elem_buf;
			for(size_t tree = 0; tree < sizes_count; ++tree) {
				size_t tree_size = 1;
				for(unsigned j = 0; j < DIM; ++j) {
					tree_size *= sizes_buf[tree * DIM + j];
				}
				if(rcells_size < tree_size) {
					free(rcells);
					rcells_size = tree_size;
					ALLOC(hig_cell *, rcells, tree_size);
				}

				elemptr += hs_deserialize(&rcells[0], elemptr);
				Point highpoint, lowpoint;
				hig_get_highpoint(rcells[0], highpoint);
				hig_get_lowpoint(rcells[0], lowpoint);

				for(size_t i = 1; i < tree_size; ++i) {
					elemptr += hs_deserialize(&rcells[i], elemptr);
					Point chigh, clow;
					hig_get_highpoint(rcells[i], chigh);
					hig_get_lowpoint(rcells[i], clow);

					for(size_t dim = 0; dim < DIM; ++dim) {
						if(chigh[dim] > highpoint[dim])
							highpoint[dim] = chigh[dim];
						if(clow[dim] < lowpoint[dim])
							lowpoint[dim] = clow[dim];
					}
				}

				// Try to integrate the new tree into each of the already
				// present trees:
				bool merged = false;
				for(unsigned i = 0; i < ctx->num_out_higs; ++i) {
					merged = groups_buf[tree] == ctx->output[i].group
						&& _try_merge_chunk_into_tree(lowpoint, highpoint, rcells,
							&sizes_buf[tree * DIM], ctx->output[i].tree);
					if(merged) {
						// The new merged tree may now be merged again...
						// Add to the corresponding list
						assert(mergeable_count <= mergeable_size);
						if(mergeable_count >= mergeable_size) {
							mergeable_size = 1 + mergeable_size * 2;
							mergeable = realloc(mergeable,
								mergeable_size * sizeof *mergeable);
						}
						mergeable[mergeable_count++] = ctx->output[i];
						ctx->output[i] = ctx->output[--ctx->num_out_higs];
						break;
					}
				}

				if(!merged) {
					// If couldn't merge, create a new output tree:
					hig_cell *ntree = hig_create_root(lowpoint, highpoint);
					hig_refine_empty(ntree, (int *)&sizes_buf[tree * DIM]);

					for(unsigned i = 0; i < tree_size; ++i) {
						rcells[i]->posinparent = i;
						rcells[i]->parent = ntree;
						ntree->children[i] = rcells[i];
					}

					// Add the new tree to the output list
					if(ctx->num_out_higs >= output_size) {
						output_size = 1 + output_size * 2;
						ctx->output = realloc(ctx->output,
							output_size * sizeof *ctx->output);
					}
					ctx->output[ctx->num_out_higs].tree = ntree;
					ctx->output[ctx->num_out_higs].group = groups_buf[tree];
					ctx->num_out_higs++;
				}
			}
			free(rcells);
		}
		free(sizes_buf);
		free(elem_buf);
		free(groups_buf);

		// Final merge of all output trees
		{
			struct _lb_output_tree *new_output;
			new_output = _try_merge_tree_list(ctx->num_out_higs, ctx->output,
				mergeable_count, mergeable, &output_size);
			free(ctx->output);
			free(mergeable);

			ctx->output = new_output;
			ctx->num_out_higs = output_size;
		}

		// Free the sending buffers
		MPI_Waitall((sizeof send_reqs) / (sizeof send_reqs[0]), send_reqs,
			MPI_STATUSES_IGNORE);

		for(int p = 0; p < num_exp_procs; ++p) {
			free(sizes_buffs[p]);
			free(tree_buffs[p]);
			free(groups_buffs[p]);
		}

		assert(output_size >= ctx->num_out_higs);
		ctx->output = realloc(ctx->output, ctx->num_out_higs * sizeof *ctx->output);
	}
	Zoltan_Destroy(&zz);
}

static void
_proc_conn_add_touching(struct _proc_conectivity *pc, int neighbor_rank, Rect *intersection)
{
	if(pc->num_neighbors >= pc->max_neighbors) {
		assert(pc->num_neighbors == pc->max_neighbors);
		pc->max_neighbors = pc->max_neighbors * 2 + 5;
		pc->neighbors = realloc(pc->neighbors, pc->max_neighbors);
	}

	unsigned n = pc->num_neighbors++;

	pc->neighbors[n].intersection = *intersection;
	pc->neighbors[n].proc = neighbor_rank;
}

static void
_proc_conn_destroy(struct _proc_conectivity *pc)
{
	free(pc->neighbors);
	free(pc);
}

static struct _proc_conectivity*
_proc_conn_parallel_find(load_balancer *ctx)
{
	DECL_AND_ALLOC(struct _proc_conectivity, pc, 1);

	if(DIM <= 5)
		// relativelly small number of possible neighbors
		pc->max_neighbors = 2 * upow(3, DIM);
	else
		// number would be too big, just use some arbitrary value
		pc->max_neighbors = 600;

	ALLOC_INFER(pc->neighbors, pc->max_neighbors);
	pc->num_neighbors = 0;

	MPI_Comm comm;

	// Take those procs who have no tress out of computation
	if(!ctx->num_out_higs) {
		MPI_Comm_split(ctx->group, MPI_UNDEFINED, 0, &comm);
		return pc;
	}

	int outer_rank;
	MPI_Comm_rank(ctx->group, &outer_rank);
	MPI_Comm_split(ctx->group, 0, outer_rank, &comm);

	int ntasks;
	MPI_Comm_size(comm, &ntasks);

	// Find local bounding box
	Rect lbounds;
	hig_get_bounding_box(ctx->output[0].tree, &lbounds);
	for(unsigned i = 1; i < ctx->num_out_higs; ++i) {
		Rect tbounds;
		hig_get_bounding_box(ctx->output[i].tree, &tbounds);

		for(unsigned dim = 0; dim < DIM; ++dim) {
			if(lbounds.lo[dim] > tbounds.lo[dim])
				lbounds.lo[dim] = tbounds.lo[dim];

			if(lbounds.hi[dim] < tbounds.hi[dim])
				lbounds.hi[dim] = tbounds.hi[dim];
		}
	}

	// Gather all procs outer ranks
	DECL_AND_ALLOC(int, proc_ranks, ntasks);
	MPI_Allgather(&outer_rank, 1, MPI_INT, proc_ranks, 1, MPI_INT, comm);

	// Gather all procs bounding boxes
	DECL_AND_ALLOC(Rect, rbounds, ntasks);
	MPI_Allgather(&lbounds, sizeof(Rect), MPI_BYTE, rbounds,
			sizeof(Rect), MPI_BYTE, comm);

	// Find the intersection between this process and every other
	for(int i = 0; i < ntasks; ++i) {
		if(outer_rank == proc_ranks[i])
			// Skip if the current process...
			continue;

		Rect isect;
		if(rect_intersect(&rbounds[i], &lbounds, &isect))
			_proc_conn_add_touching(pc, proc_ranks[i], &isect);
	}

	free(rbounds);
	free(proc_ranks);
	MPI_Comm_free(&comm);

	return pc;
}

static void _local_connectivity_append_tree_neighbor(struct _local_connectivity *lc,
	size_t tidx_orig, size_t tidx_dest, Rect *isect_orig, Point displacement)
{
	// Augment local conectivity with relevant
	// tree neighbors.
	const size_t nbs_count =
		++lc->tree_neighbors[tidx_orig].num_neighbors;

	REALLOC_INFER(lc->tree_neighbors[tidx_orig].neighbors, nbs_count);
	struct _local_neighbor *new =
		&lc->tree_neighbors[tidx_orig].neighbors[nbs_count-1];

	new->intersection = *isect_orig;
	new->neighbor_tree = tidx_dest;
	POINT_ASSIGN(new->displacement, displacement);
}

struct _local_connectivity *
_local_connectivity_build(load_balancer *ctx, struct _local_portals *lps)
{
	DECL_AND_ALLOC(struct _local_connectivity, lc, 1);
	ALLOC_INFER(lc->tree_neighbors, ctx->num_out_higs);
	memset(lc->tree_neighbors, 0, ctx->num_out_higs * sizeof *lc->tree_neighbors);
	lc->num_trees = ctx->num_out_higs;

	// Find geometric contact between trees
	{
		struct _local_neighbor curr[ctx->num_out_higs];
		unsigned curr_size = 0;
		for(unsigned i = 0; i < ctx->num_out_higs; ++i) {
			// If "i" doesn't have fringe, skip it
			if(!ctx->group_fringe[ctx->output[i].group])
				continue;

			for(unsigned j = 0; j < ctx->num_out_higs; ++j) {
				// If "i" and "j" have different groups, skip "j".
				if(ctx->output[i].group != ctx->output[j].group || i == j)
					continue;

				Rect a, b;
				hig_get_bounding_box(ctx->output[i].tree, &a);
				hig_get_bounding_box(ctx->output[j].tree, &b);
				if(rect_intersect(&a, &b, &curr[curr_size].intersection)) {
					curr[curr_size].neighbor_tree = j;
					POINT_ASSIGN_SCALAR(curr[curr_size].displacement, 0);
					++curr_size;
				}
			}
			if(curr_size) {
				ALLOC_INFER(lc->tree_neighbors[i].neighbors, curr_size);
				memcpy(lc->tree_neighbors[i].neighbors, curr,
					curr_size * sizeof *curr);
				lc->tree_neighbors[i].num_neighbors = curr_size;

				curr_size = 0;
			}
		}
		assert(curr_size == 0);
	}

	// Find portal contact between trees.
	for(struct _local_portal *iter = lps->start; iter; iter = iter->next) {
		for(struct _portal_tree_isect *i0 = iter->side_isect[0]; i0; i0 = i0->next) {
			for(struct _portal_tree_isect *i1 = iter->side_isect[1]; i1; i1 = i1->next) {
				// Move side 1 intersection to postion of side 0:
				// Step 1: calculate the displacement
				Point delta;
				POINT_SUB(delta, iter->p->r[0].lo, iter->p->r[1].lo);
				// Step 2: apply the displacement to the side 1 intersection
				Rect s1;
				rect_translate(delta, &i1->isect, &s1);

				// Test for intersection:
				Rect isect_s0;
				if(rect_intersect(&i0->isect, &s1, &isect_s0)) {
					// Warp inside the same process found!

					Point opposite_delta;
					POINT_MULT_SCALAR(opposite_delta, delta, -1);

					// Calculate the contact area on the side 1:
					Rect isect_s1;
					rect_translate(opposite_delta, &isect_s0, &isect_s1);

					// Augment local conectivity with relevant
					// tree neighbors.
					// On side 0:
					_local_connectivity_append_tree_neighbor(lc,
						i0->tree_idx, i1->tree_idx, &isect_s0,
						opposite_delta);
					// On side 1:
					_local_connectivity_append_tree_neighbor(lc,
						i1->tree_idx, i0->tree_idx, &isect_s1,
						delta);
				}
			}
		}
	}

	return lc;
}

struct _local_neighbor *
_local_connectivity_get_neighbors(struct _local_connectivity *lc,
	unsigned tree_id, unsigned *num_neighbors)
{
	*num_neighbors = lc->tree_neighbors[tree_id].num_neighbors;
	return lc->tree_neighbors[tree_id].neighbors;
}

void
_local_connectivity_destroy(struct _local_connectivity *lc)
{
	for(unsigned i = 0; i < lc->num_trees; ++i) {
		free(lc->tree_neighbors[i].neighbors);
	}
	free(lc->tree_neighbors);
	free(lc);
}


struct _search_work {
	int answer_to_proc;
	struct _search_work_start sstarts[];
};

void
_dispatch_remote_fringe_search(TermDetection *td, int to_proc_rank, int ref_proc_rank,
	unsigned num_starts, struct _search_work_start *starts)
{
	for(unsigned i = 0; i < num_starts; ++i) {
		assert(starts[i].search_from_dist > 0);
	}

	struct _search_work header = {
		.answer_to_proc = ref_proc_rank,
	};

	const void *bufs[] = {
		&header,
		starts
	};

	size_t buf_sizes[] = {
		sizeof header,
		num_starts * sizeof *starts,
		0
	};
	term_det_send_asyncv(td, bufs, buf_sizes, to_proc_rank);
}

static void destroy_portal_tree_isect(struct _portal_tree_isect *l)
{
	while(l) {
		struct _portal_tree_isect *tmp = l;
		l = l->next;
		free(tmp);
	}
}

struct portal_tree_isect_list
{
	struct _portal_tree_isect *start;
};

static void destroy_portal_tree_isect_list(struct portal_tree_isect_list *l)
{
	destroy_portal_tree_isect(l->start);
	free(l);
}

struct _local_portals* _local_portals_create(load_balancer *ctx)
{
	DECL_AND_ALLOC(struct _local_portals, ret, 1);
	ret->local_portal_count = 0;
	ret->start = NULL;

	struct _local_portal new_lp;
	memset(&new_lp, 0, sizeof new_lp);

	// Empty struct to be filled
	struct _portal_tree_isect new_isect;
	memset(&new_isect, 0, sizeof new_isect);

	// For each portal, test every local tree.
	size_t portal_gid = 0;
	for(struct _lb_portal *cur = ctx->portals; cur; cur = cur->next) {

		bool isect_found = false;
		for(size_t tidx = 0; tidx < ctx->num_out_higs; ++tidx) {
			if(ctx->output[tidx].group != cur->group) {
				continue;
			}

			// Get tree bounding box.
			hig_cell *t = ctx->output[tidx].tree;
			Rect bbox;
			hig_get_bounding_box(t, &bbox);

			// Test if each side of the portal intersects the tree.
			for(uint8_t side = 0; side < 2; ++side) {
				if(rect_intersect(&bbox, &cur->r[side], &new_isect.isect)) {
					isect_found = true;

					new_isect.tree_idx = tidx;
					new_isect.next = new_lp.side_isect[side];

					ALLOC_INFER(new_lp.side_isect[side], 1);
				       	*new_lp.side_isect[side] = new_isect;
					++new_lp.side_count[side];

					memset(&new_isect, 0, sizeof new_isect);
				}
			}
		}

		// If there is intersection with a tree, add to return.
		if(isect_found) {
			++ret->local_portal_count;

			new_lp.p = cur;
			new_lp.next = ret->start;
			new_lp.portal_gid = portal_gid;

			ALLOC_INFER(ret->start, 1);
			*ret->start = new_lp;

			memset(&new_lp, 0, sizeof new_lp);
		}
		++portal_gid;
	}

	ret->by_proc = g_hash_table_new_full(g_direct_hash, g_direct_equal,
		NULL, (GDestroyNotify)destroy_portal_tree_isect_list);

	return ret;
}

static struct _portal_tree_isect *
_local_portals_get_by_proc(struct _local_portals *lps, int prank)
{
	struct portal_tree_isect_list *l =
	       g_hash_table_lookup(lps->by_proc, GINT_TO_POINTER(prank));
	if(l) {
		return l->start;
	}
	return NULL;
}

void _local_portals_destroy(struct _local_portals* lps)
{
	struct _local_portal *lp = lps->start;
	while(lp) {
		for(uint8_t side = 0; side < 2; ++side) {
			destroy_portal_tree_isect(lp->side_isect[side]);
		}

		struct _local_portal *next = lp->next;
		free(lp);
		lp = next;
	}
	if(lps->by_proc) {
		g_hash_table_destroy(lps->by_proc);
	}
	free(lps);
}

struct serialized_portal_isection
{
	size_t portal_gid;
	unsigned group;
	uint8_t side;
	struct {
		Rect r;
		unsigned group;
	} isects[];
};

const int portal_search_tag = 8632345;

static void send_portal_side(load_balancer *ctx, allocator *alc,
	isend_pool *isp, struct _local_portal *lp, uint8_t side, int dest)
{
	struct serialized_portal_isection *buf;
	const size_t bufsize = (sizeof *buf)
		+ lp->side_count[side] * sizeof buf->isects[0];
	buf = allocator_alloc(alc, bufsize);

	buf->portal_gid = lp->portal_gid;
	buf->side = side;

	struct _portal_tree_isect *it = lp->side_isect[side];
	size_t count;
	for(count = 0; it; ++count) {
		assert(count < lp->side_count[side]);

		buf->isects[count].r = it->isect;
		buf->isects[count].group = ctx->output[it->tree_idx].group;

		it = it->next;
	}
	assert(count == lp->side_count[side]);

	isend_pool_send(isp, buf, bufsize, MPI_BYTE, dest,
		portal_search_tag, ctx->group, false);
}

static void
dispatch_search_from_portals(load_balancer *ctx, _FringeBuilder *fb,
	TermDetection *td, struct _local_portals *lps, unsigned fringesize)
{
	int nprocs;
	MPI_Comm_size(ctx->group, &nprocs);

	int myrank;
	MPI_Comm_rank(ctx->group, &myrank);

	// Number of portals each process have.
	unsigned proc_portal_count[nprocs];
	MPI_Allgather(&lps->local_portal_count, 1, MPI_UNSIGNED,
		proc_portal_count, 1, MPI_UNSIGNED, ctx->group);

	// Portal information sent to other processes
	typedef struct {
		size_t portal_gid;
		uint8_t side0: 1;
		uint8_t side1: 1;
	} portal_info;

	// Global portal information vector
	size_t total_portal_instances;
	portal_info *pinfos;

	// Portal data index, for direct reference
	struct _local_portal *pindex[ctx->num_portals];
	memset(pindex, 0, sizeof pindex);

	// Gather portal info from all processes.
	{
		// Assemble my local portal information:
		portal_info locpinfo[lps->local_portal_count];
		struct _local_portal *lp = lps->start;
		for(unsigned i = 0; i < lps->local_portal_count; ++i) {
			locpinfo[i].portal_gid = lp->portal_gid;
			locpinfo[i].side0 = lp->side_isect[0] != NULL;
			locpinfo[i].side1 = lp->side_isect[1] != NULL;

			assert(!pindex[lp->portal_gid]);
			pindex[lp->portal_gid] = lp;

			lp = lp->next;
		}

		// Organize receive buffer parameters
		int displs[nprocs];
		int counts[nprocs];
		total_portal_instances = 0;
		for(unsigned i = 0; i < nprocs; ++i) {
			displs[i] = total_portal_instances
				* sizeof(portal_info);
			total_portal_instances += proc_portal_count[i];
			counts[i] = proc_portal_count[i] * sizeof(portal_info);
		}
		ALLOC_INFER(pinfos, total_portal_instances);

		MPI_Allgatherv(locpinfo, sizeof locpinfo, MPI_BYTE, pinfos,
			counts, displs, MPI_BYTE, ctx->group);
	}

	allocator *al = allocator_init();
	isend_pool *isp = isend_pool_create();

	// Look through every received portal_info, and see if there is a match.
	unsigned to_recv_count = 0;
	unsigned i = 0;
	for(int proc = 0; proc < nprocs; ++proc) {
		for(unsigned pc = 0; pc < proc_portal_count[proc]; ++pc) {
			// Get local instance of given portal:
			struct _local_portal *lp = pindex[pinfos[i].portal_gid];
			if(lp) {
				// Check possible match and send local intersections:
				if(lp->side_isect[0] && pinfos[i].side1) {
					send_portal_side(ctx, al, isp, lp, 0, proc);
					to_recv_count++;
				}
				if(lp->side_isect[1] && pinfos[i].side0) {
					send_portal_side(ctx, al, isp, lp, 1, proc);
					to_recv_count++;
				}
			}

			++i;
		}
	}

	free(pinfos);

	// Receive each "to_recv_count" and compare if we have a neighbor
	// through a portal.
	for(i = 0; i < to_recv_count; ++i) {
		MPI_Status s;
		MPI_Probe(MPI_ANY_SOURCE, portal_search_tag, ctx->group, &s);

		int byte_count;
		MPI_Get_count(&s, MPI_BYTE, &byte_count);

		struct serialized_portal_isection *rcbuf;
		size_t isect_count =
			(byte_count - sizeof *rcbuf) / sizeof rcbuf->isects[0];

		ALLOC_INFER_FLEX(rcbuf, isects, isect_count);
		MPI_Recv(rcbuf, byte_count, MPI_BYTE, s.MPI_SOURCE,
			portal_search_tag, ctx->group, MPI_STATUS_IGNORE);

		struct _local_portal *lp = pindex[rcbuf->portal_gid];
		assert(lp);

		uint8_t lside = !rcbuf->side;

		// For each intersection between local side and remote side,
		// store it, add a remote tree neighbor and dispatch a remote
		// fringe search.
		struct _portal_tree_isect *it = lp->side_isect[lside];
		assert(it);

		Point loc2rem;
		POINT_SUB(loc2rem,
			lp->p->r[rcbuf->side].lo, lp->p->r[lside].lo);

		Point rem2loc;
		POINT_SUB(rem2loc,
			lp->p->r[lside].lo, lp->p->r[rcbuf->side].lo);

		struct _search_work_start sworks[lp->side_count[lside] * isect_count];
		unsigned sworks_count = 0;

		// Get contact list for process, or allocate if doesn't exist.
		struct portal_tree_isect_list *proc_contacts =
			g_hash_table_lookup(lps->by_proc,
				GINT_TO_POINTER(s.MPI_SOURCE));
		if(!proc_contacts) {
			ALLOC_INFER(proc_contacts, 1);
			proc_contacts->start = NULL;
			g_hash_table_insert(lps->by_proc,
				GINT_TO_POINTER(s.MPI_SOURCE), proc_contacts);
		}

		while(it) {
			Rect disp_lisect;
			rect_translate(loc2rem, &it->isect, &disp_lisect);
			for(size_t i = 0; i < isect_count; ++i) {
				// Skip if groups are different.
				if(rcbuf->isects[i].group
					!= ctx->output[it->tree_idx].group)
				{
					continue;
				}

				// For each intersection between local side and
				// remote side, store it, add a remote tree
				// neighbor and dispatch a remote fringe search.
				if(rect_intersect(&disp_lisect, &rcbuf->isects[i].r,
					&sworks[sworks_count].isect))
				{
					struct _search_work_start *sw =
						&sworks[sworks_count++];
					sw->tree_group = rcbuf->isects[i].group;
					sw->search_from_dist = fringesize;
					POINT_ASSIGN(sw->disp, loc2rem);

					DECL_AND_ALLOC(struct _portal_tree_isect, new, 1);
					new->next = proc_contacts->start;
					proc_contacts->start = new;
					new->tree_idx = it->tree_idx;

					rect_translate(rem2loc, &sw->isect,
						&new->isect);

					_fringe_builder_add_remote_tree_neighbor(fb,
						s.MPI_SOURCE, it->tree_idx, loc2rem,
						&new->isect);
				}
			}

			it = it->next;
		}

		if(sworks_count) {
			_dispatch_remote_fringe_search(td, s.MPI_SOURCE,
				myrank, sworks_count, sworks);
		}

		free(rcbuf);
	}

	isend_pool_destroy(isp);
	allocator_destroy(al);
}

static void
_parallel_extend_fringe(load_balancer *ctx, partition_graph *pg)
{
	// Find which half-portals belong to this process:
	struct _local_portals *lps = _local_portals_create(ctx);

	// Find the process bounding box intersection
	struct _proc_conectivity *pc = _proc_conn_parallel_find(ctx);

	DECL_AND_ALLOC(struct _search_work_start, isects, ctx->num_out_higs * pc->num_neighbors);
	DECL_AND_ALLOC(unsigned, touching_trees, ctx->num_out_higs * pc->num_neighbors);
	DECL_AND_ALLOC(unsigned, isects_per_proc_count, pc->num_neighbors);
	DECL_AND_ALLOC(unsigned, proc_start, pc->num_neighbors);
	memset(isects_per_proc_count, 0, pc->num_neighbors * sizeof *isects_per_proc_count);

	TermDetection *td = term_det_init(ctx->group);

	// Creates the structure holding the in construction fringes for all relevant
	// neighbouring processes.
	_FringeBuilder *fb = _fringe_builder_create(pc, lps, ctx, td);

	// Find where each tree touches other processes.
	{
		int myrank;
		MPI_Comm_rank(ctx->group, &myrank);
		Point zero_disp = {};

		unsigned next_isect = 0;
		for(size_t j = 0; j < pc->num_neighbors; ++j) {
			proc_start[j] = next_isect;

			for(size_t i = 0; i < ctx->num_out_higs; ++i) {
				// If doesn't have fringe, skip it
				if(!ctx->group_fringe[ctx->output[i].group])
					continue;

				Rect tree_box;
				hig_get_bounding_box(ctx->output[i].tree, &tree_box);

				if(rect_intersect(&pc->neighbors[j].intersection, &tree_box,
					&isects[next_isect].isect))
				{
					isects[next_isect].search_from_dist = pg->fringesize;
					isects[next_isect].tree_group = ctx->output[i].group;
					POINT_ASSIGN_SCALAR(isects[next_isect].disp, 0);

					touching_trees[next_isect] = i;
					const int nb_rank = pc->neighbors[j].proc;
					_fringe_builder_add_remote_tree_neighbor(fb, nb_rank,
						i, zero_disp, &isects[next_isect].isect);

					// Handle the special case a intersecting flat tree is
					// coplanar with a process boundary, in which case
					// a fringe copy of the local tree must be available
					// in the remote process. This is needed to properly
					// handle boundary conditions that are separated from
					// their domain by the partition.
					if(rect_are_parallel(&tree_box,
						&pc->neighbors[j].intersection))
					{
						Point disp = {0};
						_fringe_builder_search_from(fb,
							&isects[next_isect].isect, 1,
							pc->neighbors[j].proc, disp, i);
					}

					++next_isect;
					++isects_per_proc_count[j];
				}
			}

			// Send the found intersections to the other processes
			if(isects_per_proc_count[j]) {
				_dispatch_remote_fringe_search(td, pc->neighbors[j].proc,
					myrank,	isects_per_proc_count[j],
					&isects[proc_start[j]]);
			}
		}
	}

	// Jumpstart seach from portals, both inside and outside current
	// process. Don't forget to _fringe_builder_add_remote_tree_neighbor for
	// portal neighbors.
	dispatch_search_from_portals(ctx, fb, td, lps, pg->fringesize);

	// Receive tree intersections from neighbors, compute their
	// fringes and send the result back.
	{
		int from;
		int msg_size;
		struct _search_work *work;
		while((work = term_det_wait_for_msg(td, &from, &msg_size, true))) {
			size_t count = (msg_size - sizeof *work) / sizeof work->sstarts[0];

			// Get proc idx, if geometrically connected:
			unsigned geom_proc;
			bool geom_conn = _fringe_builder_get_neighbor_idx(fb, from, &geom_proc);

			// Get per process portal connection info:
			struct _portal_tree_isect *ppc =
				_local_portals_get_by_proc(lps, from);

			// We already have the local intersections with that process,
			// which are the same rects sent to it. Check each new received
			// rect with the rects sent to that process.
			for(size_t i = 0; i < count; ++i) {
				struct _search_work_start *ss = &work->sstarts[i];

				// Check the geometrical connections:
				if(geom_conn) {
					for(size_t k = 0; k < isects_per_proc_count[geom_proc];
						++k)
					{
						size_t idx = k + proc_start[geom_proc];

						Rect isect;
						const unsigned tid = touching_trees[idx];
						if(ss->tree_group != ctx->output[tid].group
							|| !rect_intersect(&ss->isect,
								&isects[idx].isect, &isect))
						{
							continue;
						}

						// The areas received and sent touch each
						// other at isect, and the trees are of the
						// same group.
						// Build the remote fringe from the local
						// tree, and span the search to other
						// neighboring (possibly remote) trees.
						_fringe_builder_search_from(
							fb, &isect,
							work->sstarts[i].search_from_dist,
							work->answer_to_proc,
							work->sstarts[i].disp,
							tid
						);
					}
				}

				// Check portal connections:
				struct _portal_tree_isect *it =
					_local_portals_get_by_proc(lps, from);
				for(; it; it = it->next) {
					if(ctx->output[it->tree_idx].group != ss->tree_group) {
						continue;
					}

					Rect isect;
					if(!rect_intersect(&ss->isect, &it->isect, &isect)) {
						continue;
					}

					_fringe_builder_search_from(fb, &isect,
						work->sstarts[i].search_from_dist,
						work->answer_to_proc,
						work->sstarts[i].disp,
						it->tree_idx
					);
				}
			}
		}
	}
	term_det_destroy(td);

	free(isects_per_proc_count);
	free(proc_start);
	free(isects);
	free(touching_trees);

	// By now, this process have in its FringeBuilder the fringe it must provide to
	// all remote processes (not just neighbors' fringes, but any process' that might
	// span as far as here).

	// Build the tree slices of the fringe's cells and sync them with the neighboring
	// processed.
	{
		unsigned count;
		struct _lb_output_tree *fringes = _fringe_builder_assemble_result(fb, &count, pg);

		// Add the fringe trees to the output.
		REALLOC_INFER(ctx->output, ctx->num_out_higs + count);
		memcpy(ctx->output + ctx->num_out_higs, fringes, count * sizeof *fringes);
		ctx->num_out_higs += count;
	}

	_fringe_builder_destroy(fb);
	_proc_conn_destroy(pc);
	_local_portals_destroy(lps);
}

void lb_calc_partition(load_balancer *ctx, partition_graph *pg)
{
	// First, distribute the trees
	_parallel_partition_and_distribute_trees(ctx);
	assert(ctx->trees == NULL);

	// Index the base trees properties in the partition_graph
	for(unsigned i = 0; i < ctx->num_out_higs; ++i)
	{
		DECL_AND_ALLOC(struct tree_properties, prop, 1);
		prop->is_fringe = false;
		prop->partition_group = ctx->output[i].group;

		g_hash_table_insert(pg->tree_props, ctx->output[i].tree, prop);
	}

	// Then extend the output to include fringe sizes
	if(pg_get_fringe_size(pg) > 0)
		_parallel_extend_fringe(ctx, pg);

	// Sort output trees according to its group
	qsort(ctx->output, ctx->num_out_higs, sizeof *ctx->output, (int(*)(const void *, const void *))_out_group_cmp);

	// Count how many trees are there in each group, and calculate the
	// stating index for the groups:
	assert(ctx->out_higs_group_start == NULL);
	ALLOC_INFER(ctx->out_higs_group_start, ctx->num_groups + 1);
	memset(ctx->out_higs_group_start, 0,
		(ctx->num_groups + 1) * sizeof *ctx->out_higs_group_start);

	// First just count how many trees in each group
	for(size_t i = 0; i < ctx->num_out_higs; ++i) {
		ctx->out_higs_group_start[ctx->output[i].group + 1]++;
	}

	// Then calculate the starting index.
	for(size_t i = 1; i <= ctx->num_groups; ++i) {
		ctx->out_higs_group_start[i] += ctx->out_higs_group_start[i-1];
	}
}

size_t lb_get_num_local_trees(load_balancer *ctx)
{
	return ctx->num_out_higs;
}

hig_cell *lb_get_local_tree(load_balancer *ctx, size_t idx, unsigned *group)
{
	assert(idx < ctx->num_out_higs);
	if(group)
		*group = ctx->output[idx].group;
	return ctx->output[idx].tree;
}

unsigned lb_get_num_local_trees_in_group(load_balancer *ctx, unsigned group)
{
	assert(group < ctx->num_groups);
	return ctx->out_higs_group_start[group+1]
		- ctx->out_higs_group_start[group];
}

hig_cell *
lb_get_local_tree_in_group(load_balancer *ctx, size_t group_idx, unsigned group)
{
	assert(group < ctx->num_groups);
	size_t idx = group_idx + ctx->out_higs_group_start[group];
	assert(idx < ctx->out_higs_group_start[group+1]);
	return ctx->output[idx].tree;
}

void lb_destroy(load_balancer *ctx)
{
	free(ctx->output);
	for(struct _lb_input_tree *ptr = ctx->trees; ptr;) {
		void *del = ptr;
		ptr = ptr->next;
		free(del);
	}
	for(struct _lb_portal *ptr = ctx->portals; ptr;) {
		void *del = ptr;
		ptr = ptr->next;
		free(del);
	}
	free(ctx->group_weight);
	free(ctx->group_fringe);
	free(ctx->out_higs_group_start);
	free(ctx);
}
