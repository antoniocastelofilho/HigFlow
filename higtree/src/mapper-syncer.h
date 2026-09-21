// Generic exchange of per-cell data across a partition, parameterized by three
// callbacks: how many elements a tree contributes, how to pack them, how to unpack.
//
// The pairing rule from higtree-parallel.h applies with full force here: send_fill
// and recv_parse must agree on the layout, and elem_counter must return exactly the
// count send_fill will write.  A counter that disagrees with the filler is a buffer
// overrun on the sender or a misparse on the receiver, with nothing in between to
// notice.

#ifndef MAPPER_SYNCER_H
#define MAPPER_SYNCER_H

#include "pdomain.h"

typedef unsigned (*ms_node_handler) (void *data, unsigned child_idx, unsigned tree_idx,
	hig_cell *tree, int *buff);

typedef unsigned (*ms_node_elem_counter) (void *data, unsigned child_idx, unsigned tree_idx,
	hig_cell *tree);

struct mapper_syncer_info {
	ms_node_elem_counter elem_counter;
	ms_node_handler send_fill;
	ms_node_handler recv_parse;
	psim_domain *psd;
	void *data;
};
typedef struct mapper_syncer_info mapper_syncer_info;

void ms_syncmapper(mapper_syncer_info *ms);

#endif
