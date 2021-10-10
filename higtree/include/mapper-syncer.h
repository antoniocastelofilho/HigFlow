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
