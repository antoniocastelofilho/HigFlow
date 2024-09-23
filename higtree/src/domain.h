#ifndef DOMAIN_H
#define DOMAIN_H

#include "higtree.h"
#include "higtree-iterator.h"
#include "mapper.h"
#include "point-mapper.h"
#include "wls.h"

//! Type of boundary conditions.
typedef enum bc_type {
	DIRICHLET = 0,
	NEUMANN = 1,
} bc_type;

//! Value type of boundary contition
typedef enum bc_valuetype {
	fixedValue = 0,
	zeroGradient = 1,
        freestream = 2,
        empty = 3,
        timedependent = 4,
	outflow = 5,
} bc_valuetype;

//! \brief A boundary condition (BC).
//! A BC contains a HiG-Tree which represents the points where the values are defined (at the cell center). The HiG-tree may be shared by another BC; thus, the is_own_hig attribute indicates whether this BC is the owner of the tree and, hence, should destroy it when not using anymore.
//! The BC also has a mapping from cell id to a number in the interval 0..numvals - 1. The attribute vals contains the value of the BC at the cell center.
typedef struct sim_boundary {
	hig_cell *bc;
	int is_own_hig;
	bc_type type;
	bc_valuetype valuetype;
        int userid;
	int numvals;
	mp_mapper *m;
	real *vals;
} sim_boundary;

//! Creates a BC of a given type, given HiG-Tree bc and a mapper m.
//! m may be NULL, in which case it will be created from the tree
sim_boundary *sb_create(hig_cell *bc, bc_type t, mp_mapper *m);

//! Sets the values of all BC's cell centers.
void sb_set_values(sim_boundary *sb, real *vals);

//! Sets the value of a given position in the interval 0..numvals - 1.
void sb_set_value(sim_boundary *sb, int i, real val);

//! Sets the user id of a given position in the interval 0..numvals - 1.
void sb_set_userid(sim_boundary *sb, int userid);

//! Sets the valuetype of a given position in the interval 0..numvals - 1.
void sb_set_valuetype(sim_boundary *sb, bc_valuetype valuetype);

//! Gets the value of a given position in the interval 0..numvals - 1.
real sb_get_value(sim_boundary *sb, int i);

//! Gets the user id of a given position in the interval 0..numvals - 1.
int sb_get_userid(sim_boundary *sb);

//! Gets the valuetype of a given position in the interval 0..numvals - 1.
bc_valuetype sb_get_valuetype(sim_boundary *sb);

//! Destroys the BC.
void sb_destroy(sim_boundary *sb);

//! Gets the mapper of the BC.
mp_mapper *sb_get_mapper(sim_boundary *bc);

//! Gets a cell iterator which traverses all cells of the BC.
higcit_celliterator *sb_get_celliterator(sim_boundary *bc);

//! Gets the HiG-Tree associated with the BC
hig_cell *sb_get_higtree(sim_boundary *bc);

#define MAXHIGTREESPERDOMAIN 40
#define MAXBCSPERDOMAIN 200

//! \brief A simulation domain (SD) contains a set of HiG-Trees which composes the domain; each domain has a low point lp and a high point hp. It has a mapper shared by all domains.
//! An SD has a set of dirichlet BCs and a set of neumann BCs.
//! An SD has an interpolator, which can be used multiple times without needing creating/destroying one each time.
//! The cwls keeps the interpolation within a cell (it is a cached interpolation).
typedef struct sim_domain {
	unsigned max_numhigtrees; //!< Currently allocated size of higtrees vector.
	int numhigtrees; //!< Actual number of higtrees.
	hig_cell **higtrees;

	//! First tree of the fringe. Trees after this index will not take part in the
	//! domain cell iterator.
	unsigned first_fringe;

	int is_own_hig;
	mp_mapper *m;
	int numdirichlet_bcs;
	sim_boundary **dirichlet_bcs;
	int numneumann_bcs;
	sim_boundary **neumann_bcs;
	struct interpolator {
		wls_interpolator *wls;
		int dim;
		int order;
		int numminpoints;
		unsigned maxpts;
	} inter, bc_inter;
	int use_cache;
	point_mapper *cwls;
} sim_domain;

typedef struct sim_facet_block_info {
	unsigned char l: 1;
	unsigned char h: 1;
} sim_facet_block_info;

//! \brief A simulation facet domain (SFD) contains a simulation domain SD, thus, SFD is a specialization of SD.
//! An SFD has a valid bounding box and a dimension of the facets which are used.
typedef struct sim_facet_domain {
	mp_mapper *__cm;
	mp_mapper *fm;
	sim_domain *cdom;
	sim_facet_block_info *sfbi[MAXHIGTREESPERDOMAIN];
	int dimofinterest[DIM];
	int dim;
} sim_facet_domain;

//! Creates a simulation domain.
sim_domain *sd_create(mp_mapper *m);

//! Adds a local domain HiG-tree to the SD.
int sd_add_higtree(sim_domain *d, hig_cell *c);

//! Adds a fringe to remote domain Higtree to the SD.
int sd_add_fringe_higtree(sim_domain *d, hig_cell *c);

//! Marks domain trees as fringe.
/*! Trees will be internally reordered, so the index of
 * the trees may change. Input array will be sorted.
 * @param trees_idx Array of tree indexes to be marked as fringe.
 * @param trees_idx_size Number of elements in trees_idx array.
 */
void sd_set_trees_as_fringe(sim_domain *d, unsigned* trees_idx,
	unsigned trees_idx_size);

//! Marks domains trees as not-fringe.
/*! Trees will be internally reordered, so the index of
 * the trees may change. Input array will be sorted.
 * @param trees_idx Array of tree indexes to be marked as local domain.
 * @param trees_idx_size Number of elements in trees_idx array.
 */
void sd_set_trees_as_local_domain(sim_domain *d, unsigned* trees_idx,
	unsigned trees_idx_size);

//! Sets whether the SD is the owner of the HiG-Trees
void sd_set_higtrees_ownership(sim_domain *d, int is_own_hig);

//! Gets the bounding box of the locally controlled domain.
void sd_get_domain_bounding_box(sim_domain *d, Rect *bbox);

//! Adds a BC.
void sd_add_boundary(sim_domain *d, sim_boundary *bc);

//! Destroys the SD.
void sd_destroy(sim_domain *d);

//! Gets the total number of HiG-Trees of the SD.
int sd_get_num_higtrees(sim_domain *d);

//! Gets the number of HiG-Trees of the SD that are fringe.
unsigned sd_get_num_fringe_higtrees(sim_domain *d);

//! Gets the number of HiG-Trees of the SD that are local domain.
unsigned sd_get_num_local_higtrees(sim_domain *d);

//! Gets the number of boundary conditions of the SD.
int sd_get_num_bcs(sim_domain *d, int type);

//! Creates the boundary conditions around the hig-tree rooted at root. The boundary conditions will be of the given type.
void sd_create_boundary(hig_cell *root, sim_domain *sd, int type);

//! Gets the i-th BC of a given type.
sim_boundary *sd_get_bc(sim_domain *d, int type, int i);

//! Gets a cell iterator for the whole domain.
higcit_celliterator *sd_get_domain_celliterator(sim_domain *d);

//! Gets a cell iterator for all BCs.
higcit_celliterator *sd_get_bcs_celliterator(sim_domain *d);

//! Gets the i-th HiG-Tree of the SD.
hig_cell *sd_get_higtree(sim_domain *d, int i);

//! Gets the i-th fringe HiG-Tree of the SD.
hig_cell *sd_get_fringe_higtree(sim_domain *d, unsigned i);

//! Gets the i-th local domain HiG-Tree of the SD.
hig_cell *sd_get_local_higtree(sim_domain *d, unsigned i);

//! Get the mapper of the SD.
mp_mapper *sd_get_domain_mapper(sim_domain *d);

//! Get local id of higcell inside the domain.
//
// This local id can be used to index local cell centered properties.
// The higcell must be contained in the domain.
int sd_get_local_id(sim_domain *sd, hig_cell *c);

//! Get global id of higcell inside the domain.
//
// This global id can be used to identify uniquely the cell inside
// the solver. The higcell must be contained in the domain.
int sd_get_global_id(sim_domain *sd, hig_cell *c);

hig_cell *sd_get_cell_with_point(sim_domain *d, const Point x);

void sd_use_cache(sim_domain *d, int v);

//! Creates a simulation facet domain (SFD).
sim_facet_domain *sfd_create(mp_mapper *m, int dim);

//! Uses the HiG-Trees of a SD (center domain) as HiG-Trees for the SFD.
void sfd_copy_higtrees_from_center_domain(sim_facet_domain *sfd, sim_domain *sd);

//! Adds a HiG-Tree to the SFD.
int sfd_add_higtree(sim_facet_domain *d, hig_cell *c);

//! Sets whether the SFD is the owner of the HiG-Trees.
void sfd_set_higtrees_ownership(sim_facet_domain *d, int is_own_hig);

//! Adds a boundary condition to the SFD.
void sfd_add_boundary(sim_facet_domain *d, sim_boundary *bc);

//! Destroys the SFD.
void sfd_destroy(sim_facet_domain *d);

//! Gets the number of HiG-Trees of the SFD.
int sfd_get_num_higtrees(sim_facet_domain *d);

//! Gets the number of HiG-Trees of the SFD that are fringe.
unsigned sfd_get_num_fringe_higtrees(sim_facet_domain *d);

//! Gets the number of HiG-Trees of the SFD that are local domain.
unsigned sfd_get_num_local_higtrees(sim_facet_domain *d);

//! Gets the number of boundary conditions of the SFD.
int sfd_get_num_bcs(sim_facet_domain *d, int type);

//! Creates the boundary conditions around the hig-tree rooted at root. The boundary conditions will be of the given type.
void sfd_create_boundary(hig_cell *root, sim_facet_domain *sd, int type);

//! Gets in which dimension the facet are at.
int sfd_get_dim(sim_facet_domain *d);

//! Gets the i-th BC of a given type.
sim_boundary *sfd_get_bc(sim_facet_domain *d, int type, int i);

//! Gets a facet iterator of all locally owned facets of the SFD.
higfit_facetiterator *sfd_get_domain_facetiterator(sim_facet_domain *d);

//! Gets a facet iterator of all fringe facets of the SFD.
higfit_facetiterator *sfd_get_fringe_facetiterator(sim_facet_domain *d);

higfit_facetiterator *sfd_get_facetiterator_for_block(sim_facet_domain *d, int hig, int pos);

//! Gets a cell iterator of all BCs of the SFD.
higcit_celliterator *sfd_get_bcs_celliterator(sim_facet_domain *d);

//! Gets the i-th HiG-Tree of the SFD.
hig_cell *sfd_get_higtree(sim_facet_domain *d, int i);

//! Gets the i-th fringe HiG-Tree of the SFD.
hig_cell *sfd_get_fringe_higtree(sim_facet_domain *d, unsigned i);

//! Gets the i-th local domain HiG-Tree of the SFD.
hig_cell *sfd_get_local_higtree(sim_facet_domain *d, unsigned i);

//! Get local id of hig_facet inside the domain.
//
// This local id can be used to index local facet properties.
// The hig_facet must be contained in the domain.
int sfd_get_local_id(sim_facet_domain *sfd, hig_facet *f);

//! Get the SFD mapper.
mp_mapper *sfd_get_domain_mapper(sim_facet_domain *d);

//! Adjusts the ids of facets, using the following principles:<BR>
//! 1) if the same facet which is shared by two cells, the cell with the lowest coordinate keeps the id and the other cell sets it to -id. <BR>
//! 2) if a facet f1 of cell c1 includes the facet f2 of cell c2 (i.e., f2 is within the bounding box of f1), then c2 keeps the id of f2 and c1 sets the id of f1 to 0 (an invalid id) <BR>
//! It assumes that a facet area is can partitioned into the areas of a set of facets.
void sfd_adjust_facet_ids(sim_facet_domain *d);

void sfd_set_sfbi(sim_facet_domain *sfd, int i, sim_facet_block_info *sfbi);

hig_cell *sfd_get_cell_with_point(sim_facet_domain *d, Point x);

int sfd_get_facet_with_point(sim_facet_domain *d, Point x, hig_facet *facet);

void sfd_use_cache(sim_facet_domain *d, int v);

//! \brief Defines the structure that keeps a stencil. A stencil has a set of
//! (id, val) pairs, where val is the weight of element represents by id.
//! It also has a rhs (right hand side) value, which is a value not associated with any element,
//! usually used to keep the dirichlet values.
typedef struct sim_stencil {
	int numelems; //! Defines how many elements are in the stencil
	int maxelems; //! Defines the maximum number of elements the stencil can keep.
	int *ids; //! Keeps the ids of the elements of the stencil
	real *vals; //! Keeps the values associated with each element, usually the weight
	real rhs; //! Keeps the values which is not associated with any element, usually the dirichlet value
} sim_stencil;

//! Sets the order of the interpolation to be used for the SD
void sd_set_interpolator_order(sim_domain *d, int order);

//! \brief Computes the stencil associated with point x, scaled by alpha, using the center of cells in SD. The results
//! are accumulates in stn, i.e., if stn is not a resetted stencil, the result of this interpolation will be
//! added to the current result of stn. Boundary conditions are considered.
void sd_get_stencil(sim_domain *d, const Point center, const Point x, real alpha, sim_stencil *stn);

//! \brief Computes the stencil associated with point x, scaled by alpha, using the center of cells in SD. The results
//! are accumulates in stn, i.e., if stn is not a resetted stencil, the result of this interpolation will be
//! added to the current result of stn. Boundary conditions are not considered.
void sd_get_stencil_without_bc(sim_domain *d, Point x, real delta, real alpha, sim_stencil *stn);

//! \brief Sets the order of the interpolation to bhe used for the SFD
void sfd_set_interpolator_order(sim_facet_domain *d, int order);

//! \brief Computes the stencil associated with point x, scaled by alpha, using the center of facets in SFD. The results
//! are accumulates in stn, i.e., if stn is not a resetted stencil, the result of this interpolation will be
//! added to the current result of stn. Boundary conditions are considered.
void sfd_get_stencil(sim_facet_domain *d, const Point center, const Point x, real alpha, sim_stencil *stn);

//! \brief Computes the stencil associated with point x, scaled by alpha, using the center of facets in SFD. The results
//! are accumulates in stn, i.e., if stn is not a resetted stencil, the result of this interpolation will be
//! added to the current result of stn. Boundary conditions are not considered.
void sfd_get_stencil_without_bc(sim_facet_domain *d, Point x, real delta, real alpha, sim_stencil *stn);

//! Creates a stencil
sim_stencil *stn_create();

//! Resets a stencil. After calling this function, the stencil will be empty.
void stn_reset(sim_stencil *stn);

//! Adds v to the value associated with id. If no element exists, creates one and set to v.
void stn_add_to_element(sim_stencil *stn, int id, real v);

//! Sets the value associated with id to v.
void stn_set_element(sim_stencil *stn, int id, real v);

//! Adds v to the right hand side.
void stn_add_to_rhs(sim_stencil *stn, real v);

//! Sets the value of the right hand side to v.
void stn_set_rhs(sim_stencil *stn, real v);

//! Gets the i-th id.
int stn_get_id(sim_stencil *stn, int i);

//! Gets the i-th value.
real stn_get_val(sim_stencil *stn, int i);

//! Gets a pointer the vector of ids.
int *stn_get_ids(sim_stencil *stn);

//! Gets a pointer the vector of values.
real *stn_get_vals(sim_stencil *stn);

//! Gets the right hand side value.
real stn_get_rhs(sim_stencil *stn);

//! Gets the number of elements in the stencil.
int stn_get_numelems(sim_stencil *stn);

//! Adds the values associated with elements in the stencil stnfrom to the stencil associated with the stencil stnto.
void stn_add_stencil(sim_stencil *stnto, sim_stencil *stnfrom);

/*! Adds s times the values associated with elements in the stencil stnfrom
 * to the stencil associated with the stencil stnto.
 *
 * I.e.: stnto += s * stnfrom
 */
void stn_mult_scalar_add(sim_stencil *stnto, sim_stencil *stnfrom, real s);

//! Destroys a stencil.
void stn_destroy(sim_stencil *stn);

//! Computes the scalar product of the stencil with the values in v.
real stn_mult_vector(sim_stencil *stn, real *v);

//! \brief Caches the result of an interpolation. It should not be used directly.
typedef struct _cache_wls_item {
	int numpts;
	int numbcpts;
	int *ids;
	sim_boundary **bcpoints;
	real *w;
} _cache_wls_item;


/*! Caches two interpolation stencils of a point. One with and one without dirichlet BCs. */
typedef struct _cache_wls_choice {
	struct _cache_wls_item *with_bc;
	struct _cache_wls_item *without_bc;
} _cache_wls_choice;

#endif
