#ifndef BUILD_FRINGE
#define BUILD_FRINGE

#include "higtree.h"
#include "term-det.h"
#include "rect.h"
#include "lbal.h"

struct _FringeBuilder;
typedef struct _FringeBuilder _FringeBuilder;

_FringeBuilder *_fringe_builder_create(struct _proc_conectivity *pcon,
	struct _local_portals *lps, load_balancer *ctx, TermDetection *td);

void _fringe_builder_add_remote_tree_neighbor(_FringeBuilder *fb, int nb_rank,
	unsigned tree_idx, CPPoint displacement, Rect *isect);

bool _fringe_builder_get_neighbor_idx(_FringeBuilder *fb, int mpi_rank, unsigned *idx);

void _fringe_builder_search_from(_FringeBuilder *fb, Rect *isect, unsigned start,
	int dest_proc_rank, CPPoint disp, unsigned tree_idx);

struct _lb_output_tree*
_fringe_builder_assemble_result(_FringeBuilder *fb, unsigned *fringe_tree_count,
	partition_graph *pg);

void _fringe_builder_destroy(_FringeBuilder *fb);

#endif
