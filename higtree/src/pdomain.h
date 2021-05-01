#ifndef __PDOMAIN_H
#define __PDOMAIN_H

#include <stdbool.h>
#include <mpi.h>

#include "domain.h"
#include "solver.h"

#define MAXPARTITIONS 256

struct neighbor_proc
{
	//! Index of this entry inside the hash table.
	unsigned idx;

	//! List of all remote fringes that mirrors this process
	//
	// This list of fringe rects must be strictly on the same order as the
	// corresponding list on the remote process.
	struct to_send_fringe {
		hig_cell* local_tree;

		/** Where this fringe starts, in tree coordinates. */
		int lo_idx[DIM];

		/** One after the position this fringe ends, in tree coordinates. */
		int hi_idx[DIM];
	} *to_send;
	unsigned to_send_count;

	//! List of all local fringes, organized by what process they mirror
	//
	// The list of rects per proc must be strictly on the same order as the
	// corresponding list on the remote process.
	hig_cell** to_recv_trees;
	unsigned to_recv_count;
};

struct tree_properties {
	//! Group of the tree inside the partition procedure.
	unsigned partition_group;

	bool is_fringe;
};

//! Defines a structure which keeps the information of the partition.
struct partition_graph {
	unsigned fringesize;

	//! The MPI communicator holding all the nodes of the full domain graph
	MPI_Comm comm;

	//! Maps neighbor's MPI rank to struct neighbor_proc
	GHashTable *neighbors;

	//! Maps tree pointers to partition related struct tree_properties
	GHashTable *tree_props;
};
typedef struct partition_graph partition_graph;

partition_graph *pg_create(MPI_Comm comm);
void pg_set_fringe_size(partition_graph *pg, unsigned fs);
unsigned pg_get_fringe_size(partition_graph *pg);

//! Retruns current MPI process rank relative to this partition MPI communicator group.
int pg_get_rank(partition_graph *pg);

//! Returns the number of MPI process in the whole partitioned domain.
int pg_get_numnodes(partition_graph *pg);

//! Returns the number of neighbors of the local process.
unsigned pg_get_num_neighbors(partition_graph *pg);

//! Returng MPI communicator for the partition.
MPI_Comm pg_get_MPI_comm(partition_graph *pg);

void pg_destroy(partition_graph *pg);

typedef struct {
	/*! Cell property send datatype.
	 *
	 * Use to send all cell centered properties in one go.
	 * May be MPI_DATATYPE_NULL if psd_create_property was not called.
	 */
	MPI_Datatype send;

	/*! Cell property receive datatype.
	 *
	 * Use to receive all cell centered properties in one go.
	 */
	MPI_Datatype recv;
} property_sync_datatype;

//! Struct containing all domain data a distributed property shares with its
//! siblings and its psim(_facet)_domain struct.
typedef struct {
	/*! Array of datatypes used in properties synchronization.
	 *
	 * One element per neighbor, indexed by nb->idx. May be NULL if
	 * psd_create_property was not called.
	 */
	property_sync_datatype *psync;

	/*! Base psim_domain.
	 *
	 * Needed for access to MPI communicator and filtered neighbors.
	 */
	struct psim_domain *bpsd;

	/*! Mapping of local ids to corresponding global ids. */
	int *gid_map;

	/*! The first global id of this process. */
	int firstid;

	/*! Number of elements owned by this process */
	int local_count;

	/*! The total count of elements locally stored. Defined as local_count
	 *! plus the number of fringe elements stored locally. */
	int total_count;
} _dp_shared;

//! Defines a partitioned simulation domain (PSD).
struct psim_domain{
	_dp_shared dp_data;

	int globaldomainsize;

	sim_domain *localdomain;
	partition_graph *pg;

	/*!
	 * Neighbors from the partition graph, but filtered
	 * with only the trees that are part of this psim domain.
	 * This is done because not every tree for the partition
	 * needs to be used on the sim domain.
	 */
	GHashTable *filtered_neighbors;

	/*! Maps the tree pointers to its index inside the sim_domain. */
	GHashTable *tree_map;
};
typedef struct psim_domain psim_domain;

//! Defines a partitioned simulation facet domain (PSFD).
typedef struct {
	_dp_shared dp_data;

	int globaldomainsize;

	sim_facet_domain *localdomain;

	/*! Keep a reference to psim_domain to reuse its partition_graph
	 *! and neighborhood data.
	 */
	psim_domain *psd;
	bool psd_managed;
} psim_facet_domain;

//! Creates a partitioned simulation domain.
psim_domain *psd_create(sim_domain *subd, partition_graph *pg);

//! Gets the local subdomain of calling process.
sim_domain *psd_get_local_domain(psim_domain *psd);

//! Gets the partition graph.
partition_graph *psd_get_partition_graph(psim_domain *psd);

//! Syncs the mapper of all subdomains, so that the ids of all cells shared by different subdomains are assigned to the respective global ids.
void psd_synced_mapper(psim_domain *psd);

//! Syncs the mapper of all subdomains, so that the ids of all facets shared by different subdomains are assigned to the respective global ids.
void psd_facet_synced_mapper(partition_graph *pt, hig_cell *root, int dim, mp_mapper *m, Point lp, Point hp, int rank, int size, int *firstid, int *numlocalelems, int *numglobalelems);

//! Sets the local subdomain.
void psd_set_local_domain(psim_domain *psd, sim_domain *d);

//! Destroys the PSD.
void psd_destroy(psim_domain *psd);

//! Gets the first global id owned by a cell of the subdomain.
int psd_get_first_id(psim_domain *psd);

//! Get global id of higcell inside the domain.
//
// This global id can be used to identify uniquely the cell inside
// the solver. The higcell must be contained in the domain.
int psd_get_global_id(psim_domain *psd, hig_cell *c);

//! Get global id from a local id of a higcell inside the domain.
//
// Prefer to use this function instead of psd_get_global_id(), when possible.
int psd_lid_to_gid(psim_domain *psd, int localid);

//! Gets the number of cells in the local subdomain.
int psd_get_local_domain_size(psim_domain *psd);

//! Gets the number of cells in the whole domain.
int psd_get_global_domain_size(psim_domain *psd);

//! Gets the global ids from an stencil.
//! Modifies the stencil state so it can't be use for interpolation.
int* psd_stn_get_gids(psim_domain *psd, sim_stencil *stn);

/*! Calculates and synchronizes joint unique global ID for all elements of all
 * psim_(facet_)domain given.
 *
 * The global ids inside a single process are contiguous. To be used when you
 * need to place multiple proprerties inside a single linear system.
 *
 * @returns The first global id of this process.
 */
int psd_set_joint_gid(size_t psd_count, psim_domain *psd[],
		size_t psfd_count, psim_facet_domain *psfd[]);

//! Creates a partitioned simulation facet domain (PSFD) sharing an existing cell psim_domain.
psim_facet_domain *psfd_create(sim_facet_domain *subd, psim_domain *psd);

//! Creates a partitioned simulation facet domain (PSFD) directly from a partition_graph.
psim_facet_domain *psfd_create_from_pg(sim_facet_domain *subd, partition_graph *pg);

//! Gets the partition graph of a psfd.
partition_graph *psfd_get_partition_graph(psim_facet_domain *psfd);

//! Distributes the facet domain from rootnode to all other nodes.
void psfd_distribute_domain(psim_facet_domain *psfd, int myrank, int rootnode, int size);

//! Syncs the mapper of all subdomains, so that the ids of all facets shared by different subdomains are assigned to the respective global ids.
void psfd_synced_mapper(psim_facet_domain *psfd);

//! Gets the local facet subdomain.
sim_facet_domain *psfd_get_local_domain(psim_facet_domain *psfd);

//! Gets the number of facets in the whole domain.
int psfd_get_global_domain_size(psim_facet_domain *psfd);

//! Destroys the PSFD.
void psfd_destroy(psim_facet_domain *psfd);

//! Gets the first global id owned by a facet of the subdomain.
int psfd_get_first_id(psim_facet_domain *psd);

//! Gets the number of facets in the local facet subdomain.
int psfd_get_local_domain_size(psim_facet_domain *psd);

void psfd_compute_sfbi(psim_facet_domain *psfd);

//! Gets the global ids from an stencil.
//! Modifies the stencil state so it can't be use for interpolation.
int* psfd_stn_get_gids(psim_facet_domain *psfd, sim_stencil *stn);

//! \brief Defines a distributed property. Which process contains only part of the values
//! of the property in the cells owned by the that process. psd_sync_distributed_property and
//! psfd_sync_distributed_property should be called to sync the values shared by
//! different nodes.
typedef struct distributed_property {
	_dp_shared *pdata;
	real *values;
} distributed_property;

//! Gets the value of the property at id.
real dp_get_value(distributed_property *dp, int id);

//! Sets the value of the property at it to val.
void dp_set_value(distributed_property *dp, int id, real val);

//! Adds val the value of the property at it.
void dp_add_value(distributed_property *dp, int id, real val);

//! Gets local size.
int dp_local_size(distributed_property *dp);

/*! Copy all values from one distributed_property to another.
 *
 * Copy both local and fringe values. Both distributed properties must be from
 * the same psim_(facet_)domain.
 */
void dp_copy_values(distributed_property *to, distributed_property *from);

/*! Sets all local values of a distributed_property to a single scalar. */
void dp_set_all_values(distributed_property *dp, real val);

/*! Loads the solution vector of the solver into the distributed property.
 *
 * There is no need to sync the property after this call.
 */
void dp_slv_load_from_solver(distributed_property *dp, solver *s);

//! Interpolates a value of the distributed property from a fully set stencil.
real dp_interpolate_from_stencil(distributed_property *dp, sim_stencil *stn);

//! Interpolates the value of a distributed property at a given sim_domain point.
//
// The distributed_property \p dp must have been created from a psim_domain
// created from sim_domain \p sd. I.e., \p sd and \p dp must be related.
//
// @param stn The stencil that will be filled and used for interpolation. Must
// be provided only for performance, so there is no need to allocate it
// internally.
real sd_dp_interpolate(sim_domain *sd, distributed_property *dp,
	const Point center, const Point p, sim_stencil *stn);

//! Interpolates the value of a distributed property at a given sim_facet_domain point.
//
// The distributed_property \p dp must have been created from a
// psim_facet_domain created from sim_facet_domain \p sfd. I.e., \p sfd and
// \p dp must be related.
//
// @param stn The stencil that will be filled and used for interpolation. Must
// be provided only for performance, so there is no need to allocate it
// internally.
real sfd_dp_interpolate(sim_facet_domain *sfd, distributed_property *dp,
	const Point center, const Point p, sim_stencil *stn);

//! Destroys a distributed property.
void dp_destroy(distributed_property *dp);

//! Creates a distributed property for the partitioned simulation domain.
distributed_property *psd_create_property(psim_domain *psd);

//! Creates a distributed property for the partitioned simulation facet domain.
distributed_property *psfd_create_property(psim_facet_domain *psfd);

//! Prints the values of a distributed property of the calling node.
void dp_print(const char *propname, int rank, distributed_property *dp);

//! \brief Syncs the distributed property over all nodes.
void dp_sync(distributed_property *dp);

//! Get global id of hig_facet inside the domain.
//
// This global id can be used to identify uniquely the facet inside
// the solver. The higfacet must be contained in the domain.
int psfd_get_global_id(psim_facet_domain *pfsd, hig_facet *f);

//! Get global id from a local id of a hig_facet inside the domain.
//
// Prefer to use this function instead of psfd_get_global_id(), when possible.<F8>
int psfd_lid_to_gid(psim_facet_domain *psfd, int localid);

// Restart simulation ---> 20_01_26 ---> Properties
/**
 * I removed the definition of this function from the code because it relied on
 * a single process having the complete global data of a distributed_property.
 * If you want to restore and fix this function to work with **actually** distributed
 * properties, the code is in commit 5a4c19ae3c315b876848107e7008513372c977a5.
 */
void loadUV(psim_facet_domain *psfdu[DIM], distributed_property *dpu[DIM], FILE *fd, int myrank, int ntasks);
void loadOneU(psim_facet_domain *psfdu, distributed_property *dpu, FILE *fd, int myrank, int ntasks);

void loadP(psim_domain *psdu, distributed_property *dpu, FILE *fd, int myrank, int ntasks);
void loadOneP(psim_domain *psdu, distributed_property *dpu, FILE *fd, int myrank, int ntasks);

void loadVF(psim_domain *psdu, distributed_property *dpu, FILE *fd, int myrank, int ntasks);
void loadOneVF(psim_domain *psdu, distributed_property *dpu, FILE *fd, int myrank, int ntasks);

void loadS(psim_domain *psdu, distributed_property *dpu[DIM][DIM], FILE *fd, int myrank, int ntasks);
void loadOneS(psim_domain *psdu, distributed_property *dpu, FILE *fd, int myrank, int ntasks);

/**
 * I removed the definition of this function from the code because it relied on
 * a single process having the complete global data of a distributed_property.
 * If you want to restore and fix this function to work with **actually** distributed
 * properties, the code is in commit 5a4c19ae3c315b876848107e7008513372c977a5.
 */
void saveUV(psim_facet_domain *psfdu[DIM], distributed_property *dpu[DIM], FILE *fd, int myrank, int ntasks);
void saveOneU(psim_facet_domain *psfdu, distributed_property *dpu, FILE *fd, int myrank, int ntasks);

void saveP(psim_domain *psdu, distributed_property *dpu, FILE *fd, int myrank, int ntasks);
void saveOneP(psim_domain *psdu, distributed_property *dpu, FILE *fd, int myrank, int ntasks);

void saveVF(psim_domain *psdu, distributed_property *dpu, FILE *fd, int myrank, int ntasks);
void saveOneVF(psim_domain *psdu, distributed_property *dpu, FILE *fd, int myrank, int ntasks);

void saveS(psim_domain *psdu, distributed_property *dpu[DIM][DIM], FILE *fd, int myrank, int ntasks);
void saveOneS(psim_domain *psdu, distributed_property *dpu, FILE *fd, int myrank, int ntasks);

// Restart simulation ---> 20_01_26 ---> Properties Cells
void load_property_cell(psim_domain *psd, distributed_property *dp, FILE *fd);

void save_property_cell(psim_domain *psd, distributed_property *dp, FILE *fd, int myrank, int ntasks);

#endif
