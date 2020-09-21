#ifndef HIGTREE_PARALLEL
#define HIGTREE_PARALLEL

#include "higtree.h"
#include "mapper.h"
#include "pdomain.h"

//! \brief Sends a hig-tree to the node torank. The process in torank should call higp_receive_uniform_tree_from_node.
//! The hig-tree is assumed to be uniform.
void higp_send_uniform_tree_to_node(hig_cell *root, int torank);

//! \brief Receives a hig-tree from the node fromrank. The process in fromrank should call higp_send_uniform_tree_to_node.
//! The hig-tree is assumed to be uniform.
hig_cell *higp_receive_uniform_tree_from_node(int fromrank);

//! \brief Sends the information about the refinement of the hig-tree root to the note torank. The process in torank should call higp_receive_uniform_refinement_from_node.
void higp_send_uniform_refinement_to_node(hig_cell *root, int torank);

//! \brief Receives the information about the refinement of the hig-tree root from the note fromrank. The process in fromrank should call higp_send_uniform_refinement_to_node.
void higp_receive_uniform_refinement_from_node(hig_cell *root, int fromrank);

//! \brief Syncs the information of the mapper to the node torank. The process in torank should call higp_sync_mapper_from_node.
void higp_sync_mapper_to_node(hig_cell *root, mp_mapper *m, int torank);

//! \brief Syncs the information of the mapper from the node fromrank. The process in fromrank should call higp_sync_mapper_to_node.
void higp_sync_mapper_from_node(hig_cell *root, mp_mapper *m, int fromrank);

//! Fill the sending buffer with information for remote process to sync with our mapper. Pairs with higp_sync_mapper_from_recv_buf().
unsigned higp_sync_mapper_fill_send_buf(hig_cell *root, mp_mapper *m, int *gid_buff);

//! Syncs mapper with information received from remote process. Pairs with higp_sync_mapper_fill_send_buf().
unsigned higp_sync_mapper_from_recv_buf(hig_cell *root, mp_mapper *m, int *gid_buff);

//! \brief Syncs the information of a mapper associated with a facet.
unsigned higp_sync_facet_mapper_from_facetiterator(higfit_facetiterator *fit, int is_receive, mp_mapper *m, int *gid_buff);

#endif
