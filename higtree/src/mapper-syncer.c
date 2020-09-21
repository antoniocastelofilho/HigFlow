#include <string.h>

#include "mapper-syncer.h"

static const int _ms_tag = 3850298;

static int *
_ms_organize_send(mapper_syncer_info *ms, sim_domain *sd, struct neighbor_proc *np,
	unsigned *elem_count)
{
	// First calculate total buffer size...
	*elem_count = 0;
	for(unsigned i = 0; i < np->to_send_count; ++i) {
		hig_cell *tree = np->to_send[i].local_tree;
		unsigned tree_idx = GPOINTER_TO_UINT(
			g_hash_table_lookup(ms->psd->tree_map, tree)
		);

		int *lo = np->to_send[i].lo_idx;
		int *hi = np->to_send[i].hi_idx;

		int pp[DIM];
		memcpy(pp, lo, sizeof(pp));
		do {
			unsigned idx = hig_toposition(tree->numcells, pp);
			*elem_count += ms->elem_counter(ms->data, idx, tree_idx, tree);
		} while(subgrid_coord_advance(pp, lo, hi));
	}

	// Allocate and fill te buffer
	DECL_AND_ALLOC(int, buff, *elem_count);
	int *ret = buff;
	mp_mapper *m = sd_get_domain_mapper(sd);
	for(unsigned i = 0; i < np->to_send_count; ++i) {
		hig_cell *tree = np->to_send[i].local_tree;
		unsigned tree_idx = GPOINTER_TO_UINT(
			g_hash_table_lookup(ms->psd->tree_map, tree)
		);

		int *lo = np->to_send[i].lo_idx;
		int *hi = np->to_send[i].hi_idx;

		int pp[DIM];
		memcpy(pp, lo, sizeof(pp));
		do {
			unsigned idx = hig_toposition(tree->numcells, pp);
			buff += ms->send_fill(ms->data, idx, tree_idx, tree, buff);
		} while(subgrid_coord_advance(pp, lo, hi));
	}

	return ret;
}

static void _ms_receive(mapper_syncer_info *ms, sim_domain *sd, MPI_Status *status)
{
	psim_domain *psd = ms->psd;

	int count;
	MPI_Get_count(status, MPI_INT, &count);

	partition_graph *pg = psd_get_partition_graph(psd);

	int buff[count];
	unsigned buf_idx = 0;
	MPI_Recv(buff, count, MPI_INT, status->MPI_SOURCE, _ms_tag,
		pg_get_MPI_comm(pg), MPI_STATUS_IGNORE);

	mp_mapper *m = sd_get_domain_mapper(sd);

	struct neighbor_proc *np;
	np = g_hash_table_lookup(psd->filtered_neighbors,
		GINT_TO_POINTER(status->MPI_SOURCE));

	for(unsigned i = 0; i < np->to_recv_count; ++i) {
		hig_cell *tree = np->to_recv_trees[i];
		unsigned tree_idx = GPOINTER_TO_UINT(
			g_hash_table_lookup(ms->psd->tree_map, tree)
		);

		unsigned numchildren = hig_get_number_of_children(tree);
		for(unsigned j = 0; j < numchildren; ++j) {
			buf_idx += ms->recv_parse(ms->data, j,
				tree_idx, tree, &buff[buf_idx]);
		}
	}
}

void ms_syncmapper(mapper_syncer_info *ms) {
	psim_domain *psd = ms->psd;
	partition_graph *pg = psd_get_partition_graph(psd);
	sim_domain *sd = psd_get_local_domain(psd);
	MPI_Comm comm = pg_get_MPI_comm(pg);
	unsigned nbs_count = g_hash_table_size(psd->filtered_neighbors);

	MPI_Request send_reqs[nbs_count];
	int *send_buff[nbs_count];
	unsigned send_count = 0;
	unsigned recv_count = 0;

	// First, the asynchronous send loop,
	GHashTableIter iter;
	g_hash_table_iter_init (&iter, psd->filtered_neighbors);

	gpointer key;
	struct neighbor_proc *np;
	while (g_hash_table_iter_next (&iter, &key, (gpointer *)&np)) {
		int proc = GPOINTER_TO_INT(key);
		// Organize and send data to neighbor process.
		unsigned count;
		send_buff[send_count] = _ms_organize_send(ms, sd, np, &count);
		MPI_Isend(send_buff[send_count], count, MPI_INT, proc,
			_ms_tag, comm, &send_reqs[send_count]);
		++send_count;

		// Before going on, we try to receive as much pending messages
		// as we can, so to not leave the network buffers full.
		// TODO: test if this makes any difference in runtime or not
		for(;;) {
			MPI_Status status;
			int recv_ready = 0;
			MPI_Iprobe(MPI_ANY_SOURCE, _ms_tag, comm, &recv_ready, &status);
			if(recv_ready) {
				_ms_receive(ms, sd, &status);
				++recv_count;
			} else {
				break;
			}
		}
	}

	// Then receive what hasn't been received yet
	for(; recv_count < nbs_count; ++recv_count) {
		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE, _ms_tag, comm, &status);
		_ms_receive(ms, sd, &status);
	}

	// Sending cleanup
	assert(send_count == nbs_count);
	MPI_Waitall(send_count, send_reqs, MPI_STATUSES_IGNORE);
	for(unsigned i = 0; i < send_count; ++i) {
		free(send_buff[i]);
	}
}
