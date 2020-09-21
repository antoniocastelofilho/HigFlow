
#include "mpi.h"
#include "Debug-c.h"
#include "higtree-parallel.h"


enum {CHILDREN_TAG = 0, COORD_TAG = 1, ID_TAG = 2};


static void __higp_send_uniform_refinement_to_node_aux(hig_cell *root, int torank) {
  int sendhaschildren;
  int numchildren = hig_get_number_of_children(root);
  if (numchildren == 0) {
    sendhaschildren = 0;
  } else {
    sendhaschildren = 1;
  }
  (void)MPI_Send(&sendhaschildren, 1, MPI_INT, torank, CHILDREN_TAG, MPI_COMM_WORLD);
  if (numchildren != 0) {
    int nc[DIM];
    hig_get_cells_per_dim(root, nc);
    (void)MPI_Send(nc, DIM, MPI_INT, torank, CHILDREN_TAG, MPI_COMM_WORLD);
    for (int i = 0; i < numchildren; i++) {
      __higp_send_uniform_refinement_to_node_aux(hig_get_child(root, i), torank);
    }
  }
}

void higp_send_uniform_refinement_to_node(hig_cell *root, int torank) {
  assert(root != NULL);
  __higp_send_uniform_refinement_to_node_aux(root, torank);
}

void higp_send_uniform_tree_to_node(hig_cell *root, int torank) {
  assert(root != NULL);
  (void)MPI_Send(root->lowpoint, DIM, MPI_HIGREAL, torank, COORD_TAG, MPI_COMM_WORLD);
  (void)MPI_Send(root->highpoint, DIM, MPI_HIGREAL, torank, COORD_TAG, MPI_COMM_WORLD);
  higp_send_uniform_refinement_to_node(root, torank);
}

static void __higp_receive_uniform_refinement_from_node_aux(hig_cell *r, int fromrank) {
  MPI_Status status;
  int haschildren;
  (void)MPI_Recv(&haschildren, 1, MPI_INT, fromrank, CHILDREN_TAG, MPI_COMM_WORLD, &status);
  if (haschildren != 0) {
     int nc[DIM];
     (void)MPI_Recv(nc, DIM, MPI_INT, fromrank, CHILDREN_TAG, MPI_COMM_WORLD, &status);
     hig_refine_uniform(r, nc);
     int numchildren = hig_get_number_of_children(r);
     for (int i = 0; i < numchildren; i++) {
       __higp_receive_uniform_refinement_from_node_aux(hig_get_child(r, i), fromrank);
     }
  }
}

void higp_receive_uniform_refinement_from_node(hig_cell *r, int fromrank) {
  __higp_receive_uniform_refinement_from_node_aux(r, fromrank);
}

hig_cell * higp_receive_uniform_tree_from_node(int fromrank) {
  MPI_Status status;
  Point l, h;
  (void)MPI_Recv(l, DIM, MPI_HIGREAL, fromrank, COORD_TAG, MPI_COMM_WORLD, &status);
  (void)MPI_Recv(h, DIM, MPI_HIGREAL, fromrank, COORD_TAG, MPI_COMM_WORLD, &status);
  hig_cell *root = hig_create_root(l, h);
  higp_receive_uniform_refinement_from_node(root, fromrank);
  return root;
}

void higp_sync_mapper_to_node(hig_cell *root, mp_mapper *m, int torank) {
	higcit_celliterator *it;
	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		int cid = hig_get_cid(c);
		int gid = mp_lookup(m, cid);
		(void)MPI_Send(&gid, 1, MPI_INT, torank, ID_TAG, MPI_COMM_WORLD);
	}
	higcit_destroy(it);
}

void higp_sync_mapper_from_node(hig_cell *root, mp_mapper *m, int fromrank) {
	MPI_Status status;
	higcit_celliterator *it;
	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		int cid = hig_get_cid(c);
		int gid;
		(void)MPI_Recv(&gid, 1, MPI_INT, fromrank, ID_TAG, MPI_COMM_WORLD, &status);
		mp_assign(m, cid, gid);
	}
	higcit_destroy(it);
}

unsigned
higp_sync_mapper_fill_send_buf(hig_cell *root, mp_mapper *m, int *gid_buff)
{
	unsigned count = 0;
	higcit_celliterator *it;
	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it))
	{
		hig_cell *c = higcit_getcell(it);
		int cid = hig_get_cid(c);
		gid_buff[count++] = mp_lookup(m, cid);
	}
	higcit_destroy(it);

	return count;
}

unsigned higp_sync_mapper_from_recv_buf(hig_cell *root, mp_mapper *m, int *gid_buff)
{
	unsigned count = 0;
	higcit_celliterator *it;
	for(it = higcit_create_all_leaves(root);
		!higcit_isfinished(it); higcit_nextcell(it))
	{
		hig_cell *c = higcit_getcell(it);
		int cid = hig_get_cid(c);
		mp_assign(m, cid, gid_buff[count++]);
	}
	higcit_destroy(it);

	return count;
}

unsigned higp_sync_facet_mapper_from_facetiterator(higfit_facetiterator *fit,
	int is_receive, mp_mapper *m, int *gid_buff)
{
	unsigned count;
	for(count = 0; !higfit_isfinished(fit); ++count, higfit_nextfacet(fit))
	{
		hig_facet *f = higfit_getfacet(fit);
		uniqueid id = hig_get_fid(f);
		if (is_receive) {
			mp_assign(m, id, gid_buff[count]);
		} else {
			gid_buff[count] = mp_lookup(m, id);
		}
	}
	return count;
}
