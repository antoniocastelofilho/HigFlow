#ifndef __LBAL__
#define __LBAL__

#include <mpi.h>
#include <stdbool.h>
#include "pdomain.h"
#include "higtree-io.h"

struct _lb_input_tree {
	struct _lb_input_tree *next;
	hig_cell *tree;
	unsigned group;
	bool managed;
};

struct _lb_output_tree {
	hig_cell *tree;
	unsigned group;
};

struct _lb_portal {
	struct _lb_portal *next;

	Rect r[2];
	unsigned group;
};

//! Connection between a portal and local trees.
struct _local_portal {
	struct _local_portal *next;
	struct _lb_portal *p;

	//! Global unique identifier for the portal.
	size_t portal_gid;

	//! The intersection between a portal and a local output tree.
	size_t side_count[2];
	struct _portal_tree_isect {
		size_t tree_idx;
		Rect isect;
		struct _portal_tree_isect *next;
	} *side_isect[2];
};

struct _local_portals {
	unsigned local_portal_count;
	struct _local_portal *start;
	GHashTable *by_proc;
};

//! Filter which portals intersect with local trees.
struct _local_portals* _local_portals_create(load_balancer *ctx);

void _local_portals_destroy(struct _local_portals* lp);


struct load_balancer {
	MPI_Comm group;
	struct _lb_input_tree *trees;
	struct _lb_output_tree *output;
	struct _lb_portal *portals;
	float *group_weight;
	bool *group_fringe;

	unsigned int first_tree_id;
	unsigned int num_in_higs;
	unsigned int num_portals;
	unsigned int num_groups;
	size_t num_out_higs;
	size_t *out_higs_group_start;
};
typedef struct load_balancer load_balancer;

/** Creates the context for partitioning HigTrees among processes.
 *
 * Creates a context for distributed load-balance partitioning of trees among
 * all process in MPI communication group @p comm.
 *
 * This is a collective call, and must be performed by all process in @p comm.
 *
 * @param comm An MPI communication group containing all process that will
 * take part on the partition.
 *
 * @param num_groups
 * @parblock
 * Trees taking part in partition are organized in groups. Each
 * group has its own weight, and each tree will be kept in its own group after
 * the partitioning.
 *
 * Groups are useful in differentiating trees of the boundary conditions from
 * trees of the internal domain, for instance.
 *
 * This parameter @p num_groups gives the number of groups that will be used on this
 * partitioning, and must be given the same for all processes in MPI group @p comm. Must
 * be 1 or more. Tree groups are identified from 0 to (@p num_groups - 1).
 * @endparblock
 *
 * @return The parallel load-balance context to be passed to other parallel load-balance
 * related functions. Must be freed with lb_destroy().
 */
load_balancer * lb_create(MPI_Comm comm, unsigned num_groups);

/** Sets the weight the leaves of the trees in the group will have in the partitioning.
 *
 * Is *not* a collective call, and will affect only the group's trees supplied by
 * the calling process.
 *
 * By default, group 0 have @p weight set to 1.0, and other groups have it
 * set to 0.0001.
 *
 * This function must be called before lb_calc_partition().
 *
 * @param lb Partition context created with lb_create().
 * @param group Group's id whose weight will be set.
 * @param weight The weight of each tree leaf in the partition, default is 1. Set
 * to 0 if you don't want the trees of this group to affect the partitioning.
 */
void lb_set_group_weight(load_balancer *lb, unsigned group, float weight);

/** Sets a flag telling if the trees of a particular group will have remote fringes.
 *
 * Is *not* a collective call, but fringe will be generated only between process
 * that set this true to the that group.
 *
 * By default, all groups have @p has_fringe flag set to true.
 *
 * This function must be called before lb_calc_partition().
 *
 * @param lb Partition context created with lb_create().
 * @param group Group's id whose flag be set.
 * @param has_fringe Set to true if this group's trees must be extended with a
 * fringe, and false if they must not.
 */
void lb_set_group_has_fringe(load_balancer *lb, unsigned group, bool has_fringe);

/** Adds a local tree to be partitioned and distributed.
 *
 * Is *not* a collective call. Each process may add zero or more trees to be
 * redistributed. Each tree is entitled to a group, that will define the weight
 * of its leaves.
 *
 * This function must be called before lb_calc_partition().
 *
 * @param lb Partition context created with lb_create().
 * @param tree Pointer to the tree to have its nodes redistributed.
 * @param managed Tells if the partitioning can change the trees in-place. If
 * @p managed is true, the memory management of the tree is handled by the
 * load-balance procedure, and may be modified or even freed by it, and may be part
 * of the output trees. If @p managed is false, the original tree will be left
 * unchanged after the procedure, and the corresponding (possibly modified) output
 * tree will be a clone.
 * @param group The id of the group the tree will belong to.
 */
void lb_add_input_tree(load_balancer *lb, hig_cell *tree, bool managed,
		unsigned group);

/** Adds a space warping portal connecting two remote areas of the domain.
 *
 * This is a collective call. All portals must be specified in the same
 * order for all processes. Used for periodic domains.
 *
 * Rects a and b must have exactly the same size, and must have exactly one
 * dimension with size 0. These preconditions are enforced with assert().
 * The caller is responsible for sanitizing user input.
 *
 * The behavior of portals is well defined if, and only if:
 * 1) they are placed on the limit of a tree and free space;
 * 2) there is enough free space behind the portal, so a fringe fits there;
 * 3) a-side wall must face the opposite direction of b-side wall.
 */
void lb_add_portal(load_balancer *lb, const Rect *a, const Rect *b,
		unsigned group);

/** Perform the actual partitioning and load-balancing.
 *
 * This is a collective call.
 *
 * This function must be called only once per partition context @p lb. After
 * that, the partitioned local trees may be retrieved with
 * lb_get_num_local_trees_per_group() and
 * lb_get_local_tree().
 *
 * @param lb The partitining context, containing the trees contributed by this
 * process.
 * @param[out] pt The to-be-filled partition table.
 */
void lb_calc_partition(load_balancer *lb, partition_graph *pg);

/** Retrieve the number of trees output by the partitioning.
 *
 * Is *not* a collective call. Must be called after lb_calc_partition().
 *
 * @return The number of trees entitled to this process.
 */
size_t lb_get_num_local_trees(load_balancer *lb);

/** Retrieve a tree after the load balancing.
 *
 * Is *not* a collective call. Must be called after lb_calc_partition().
 *
 * Each call to this function will return one of the trees in the new, balanced,
 * partition, indexed by @p idx. Indexes starts from 0, and must be less than
 * the total number of trees, as given by lb_get_num_local_trees()
 * function.
 *
 * Each tree that can be returned by this function must
 * be managed by the caller, and it is an error to not retrieve all the available
 * trees.
 *
 * @param idx The index of the output tree. The output trees are ordered and
 * grouped by their group's ids, so trees from group 0 have the lower @p idx values,
 * from group 1 have the next ones, and so forth.
 * @param[out] group A pointer to be filled by the tree's group's id. May be NULL, if
 * you don't need that information.
 */
hig_cell *lb_get_local_tree(load_balancer *lb, size_t idx, unsigned *group);

/** Retrieve the number of trees in group output by the partitioning.
 *
 * Is *not* a collective call. Must be called after lb_calc_partition().
 *
 * @param group The group you want to know how many trees you have.
 * @return The number of local trees in group.
 */
unsigned lb_get_num_local_trees_in_group(load_balancer *lb, unsigned group);

/** Retrieve a tree in specific group after the load balancing.
 *
 * Is *not* a collective call. Must be called after lb_calc_partition().
 *
 * Each call to this function will return one of the trees in the new, balanced,
 * partition, indexed by @p group_idx. Group index starts from 0, and must be less
 * than the total number of trees in that group, as given by
 * lb_get_num_local_trees_in_group() function.
 *
 * Each tree that can be returned by this function must be managed by the caller,
 * and it is an error to not retrieve all the available trees in the partition.
 *
 * @param group_idx The index of the output tree inside its group.
 * @param group The group from where to retrieve the tree.
 */
hig_cell *
lb_get_local_tree_in_group(load_balancer *lb, size_t group_idx, unsigned group);


/** Destroys a parallel load balance context.
 *
 * Destroys a parallel load balance context created with lb_create().
 */
void lb_destroy(load_balancer *lb);

/** Wire format of a finge search activity initiated by a remote process. */
struct _search_work_start {
	unsigned search_from_dist;
	unsigned tree_group;
	Point disp;
	Rect isect;
};
struct TermDetection;

/** Initiate an asynchronous remote fringe search using a termination detection. */
void _dispatch_remote_fringe_search(struct TermDetection *td, int to_proc_rank,
	int ref_proc_rank, unsigned num_starts,
	struct _search_work_start *starts);

/** Struct used internally during partition. */
struct _local_neighbor {
	unsigned neighbor_tree;
	Rect intersection;

	/** Displacement needed to reach neighbor_tree.
	 *
	 * Different from zero only when the connection is via portal.
	 */
	Point displacement;
};

/** Describes what trees interfaces with each other in the local process.
 * Used temporarily during the building of the fringes.
 */
struct _local_connectivity
{
	struct {
		struct _local_neighbor *neighbors;
		size_t num_neighbors;
	} *tree_neighbors;
	unsigned num_trees;
};

struct _local_connectivity *_local_connectivity_build(load_balancer *lb,
	struct _local_portals *lps);

struct _local_neighbor *_local_connectivity_get_neighbors(struct _local_connectivity *lc,
	unsigned tree_id, unsigned *num_neighbors);

void _local_connectivity_destroy(struct _local_connectivity *lc);

/** Rough connectivity among processes. Used internally during fringe building.
 *
 * Holds the touching area between local process' bounding box and neighboring
 * processes' bounding boxes.
 */
struct _proc_conectivity {
	struct {
		Rect intersection;
		int proc;
	} *neighbors;
	unsigned max_neighbors;
	unsigned num_neighbors;
};

#endif
