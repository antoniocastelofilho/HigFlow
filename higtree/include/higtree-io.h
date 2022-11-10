#ifndef HIGTREE_IO_H
#define HIGTREE_IO_H

#include<stdbool.h>

#include "types.h"
#include "higtree.h"
#include "higtree-iterator.h"
#include "pdomain.h"

//! Prints in file f the information about the tree rooted in *parent*. <BR>
//! For debuggin purpose.
void higio_print(FILE *f, hig_cell *parent);

//! Prints in file f the 2D (dim = 2) cells returned by iterator *it*, in vtk format.
void higio_print_celliterator_in_vtk2d(FILE *f, higcit_celliterator *it);

//! Prints in file f the 3D (dim = 3) cells returned by iterator *it*, in vtk format.
void higio_print_celliterator_in_vtk3d(FILE *f, higcit_celliterator *it);

//! Prints in file f the DIM-D cells returned by iterator *it*, in vtk format.
void higio_print_celliterator_in_vtk(FILE *f, higcit_celliterator *it);

//! Prints in file f the 2D (dim = 2) cells in the tree rooted at *parent*, in vtk format.
void higio_print_in_vtk2d(FILE *f, hig_cell *parent);

//! Prints in file f the 3D (dim = 3) cells in the tree rooted at *parent*, in vtk format.
void higio_print_in_vtk3d(FILE *f, hig_cell *parent);

//! Prints in file f the DIM-D cells in the tree rooted at *parent*, in vtk format.
void higio_print_in_vtk(FILE *f, hig_cell *parent);

//! Prints in file f the coordinates of the leafes of the tree rooted at *root*.
void higio_write_leaves(FILE *f, hig_cell *root);

//! Reads the information about a tree, constructs the tree and returns it. <BR>
//! The file should have been written by higio_write_to_file.
hig_cell *higio_read_from_file(FILE *f);

//! Writes the information the tree rooted at *root* in file f.<BR>
//! Use higio_read_from_file to read it.
void higio_write_to_file(FILE * f, hig_cell *root);

//! Reads the information of an AMR 2D and builds the corresponding tree.
hig_cell * higio_read_from_amr2d(FILE *fd);

//! Reads the information of an AMR 3D and builds the corresponding tree.
hig_cell * higio_read_from_amr3d(FILE *fd);

//! Reads the information of an AMR DIM-D and builds the corresponding tree.
hig_cell * higio_read_from_amr(FILE *fd);

//! \brief Represents one patch in an AMR description.
typedef struct __higio_patch_info {
	int initcell[DIM];
	int patchsize[DIM];
} __higio_patch_info;

//! \brief Represents one level in an AMR description, i.e. a set of patches.
typedef struct __higio_level_info {
	int level;
	int numpatches;
	Point delta;
	__higio_patch_info *patches;
} __higio_level_info;

//! \brief Represents an AMR description, i.e. a set of levels, which has a set of patches.
typedef struct higio_amr_info {
	Point l, h;
	int numlevels;
	__higio_level_info *levels;
} higio_amr_info;

//! Gets the low point of a tree represented by the AMR description.
void higio_get_domain_lowpoint_from_mesh_info(higio_amr_info *mi, Point l);

//! Gets the delta in the coarest refined level.
void higio_get_domain_delta_first_level_from_mesh_info(higio_amr_info *mi, Point delta);

//! Reads the AMR description.
higio_amr_info *higio_read_amr_info(FILE *fd);

//! Duplicates original AMR description.
higio_amr_info *higio_clone_amr_info(higio_amr_info *orig);

//! Free memory used by an instance of higio_amr_info returned from cloning or reading a file.
void higio_amr_info_destroy(higio_amr_info *mi);

//! Reads the information of an AMR info and builds the corresponding tree.
hig_cell * higio_read_from_amr_info(higio_amr_info *mi);

//! Computes how many cells there are in the tree rooted in each cell of the first refinement leval.
//! The complete tree is not built.
void higio_compute_cells_per_coarser_block_from_amr_info(higio_amr_info *mi, int numcells[], int numblocks);

struct load_balancer;
typedef struct load_balancer load_balancer;

//! Builts the trees corresponding to node myrank and fills the partition graph.
//! It should be called by each process node in the gaph. Parameter @p mi is only
//! used by rank 0 process.
load_balancer* higio_partition_from_amr_info(partition_graph *pg, higio_amr_info *mi);

//! Builts the tree corresponding to boundary condition a dimension bcdim and direction bcdir at node myrank. It should be called by each process node.
hig_cell *higio_read_bc_from_amr(higio_amr_info *mi, int bcdim, int bcdir);

//! \brief An HDF5 file-context for reading/writing trees and domains
struct higio_hdf5;
typedef struct higio_hdf5 higio_hdf5;

//! Creates an HDF5 file for writing.
higio_hdf5 *higio_create_hdf5(const char *filename);

//! Opens an HDF5 file for reading.
higio_hdf5 *higio_open_hdf5(const char *filename);

//! Closes an HDF5 file, and free the resources used by it.
void higio_close_hdf5(higio_hdf5* ctx);

//! Writes down a HigTree in serialized format.
void higio_write_tree_hdf5(higio_hdf5* ctx, const char* dataset_path, hig_cell* root);

//! Loads up a HigTree from a dataset in HDF5. The returned tree must be freed with hig_destroy().
hig_cell *higio_load_tree_hdf5(higio_hdf5* ctx, const char* dataset_path);

//! Writes down a cell centered distributed property in HDF5.
//
// The values are laid down in the same sequence as they appears in the domain's trees
// breadth-first searchs.
void higio_write_cell_property_hdf5(higio_hdf5* ctx, const char* dataset_path, sim_domain *sd, distributed_property * const* dp, size_t prop_dimension);

//! Loads a cell centered distributed property written down by higio_write_cell_property_hdf5() .
//
// @returns 0 on success, -1 on error
int higio_read_cell_property_hdf5(higio_hdf5* ctx, const char* dataset_path, sim_domain *sd, distributed_property * const* dp, size_t prop_dimension);

//! Writes down a facet centered distributed property in HDF5.
//
// The values are laid down in the same sequence as they appears in the domain's trees
// breadth-first searchs.
//
// The property may have multiple dimensions, but all dimensions must share the same
// facet domain, i.e. all dimensions must be colocalized at the same facets centers.
//
// For instance, if every facet has 2 components U and V, you may use a single call
// to write them all with prop_dimension == 2. Otherwise, if facets containing U (facing
// sideways) are different from facets containing V (facing upwards), then you must make
// one call for each component, each with prop_dimension == 1.
void higio_write_facet_property_hdf5(higio_hdf5* ctx, const char* dataset_path, sim_facet_domain *sfd, distributed_property *const * dp, size_t prop_dimension);

//! Loads a facet centered distributed property written down by higio_write_face_property_hdf5() .
// @returns 0 on success, -1 on error
int higio_read_facet_property_hdf5(higio_hdf5* ctx, const char* dataset_path, sim_facet_domain *sfd, distributed_property *const * dp, size_t prop_dimension);

//! \brief An HDF5 file-context for reading/writing trees and domains
struct xdmf_output;
typedef struct xdmf_output xdmf_output;

//! Creates a context for writing visualization data to a XDMF+HDF5 file;
//
// The file_prefix parameter will prefix the file names generated when writing
// the data.
//
// The sd paramater is the domain containing the tree where all written
// properties belongs to.
//
// The pointer returned by this function must be
// freed with xdmf_destroy()
xdmf_output *xdmf_init(const char* file_prefix, sim_domain* sd);

//! Destroy a XDMF context created with xdmf_init()
void xdmf_destroy(xdmf_output *ctx);

//! Returns XDMF files prefix, set when creating the context.
const char* xdmf_get_file_prefix(xdmf_output *ctx);

//! Register a cell centered property to be written by the XDMF context.
//
// Every property registered with this function will be output at each
// call to the function xdmf_write_timestep(). This property must
// belong to the domain passed to xdmf_init().
//
// Argument prop_name is the meaningful name referring to the property
// inside the XDMF file (i.e. "Pressure", "Velocity", etc). It will be
// used by the visualization tools.
//
// The property may be a scalar or vector field, depending on the value
// of num_componencts argument. Following it, one distributed_property
// for each component must be passed. Recommended values of num_components
// best accepeted by visualization tools are either 1, for scalar fields,
// or DIM, for vector fields.
void xdmf_register_cell_property(xdmf_output *ctx, const char* prop_name,
	size_t num_components, distributed_property *dp, ...);

//! Register a facet centered property to be written by the XDMF context.
//
// Every property registered with this function will be output at each
// call to the function xdmf_write_timestep(). This property must
// belong to the domain passed to xdmf_init().
//
// Argument prop_name is the meaningful name referring to the property
// inside the XDMF file (i.e. "Pressure", "Velocity", etc). It will be
// used by the visualization tools.
//
// The property may be a scalar or vector field, depending on the value
// of num_componencts argument. Following it, one sim_facet_domain and
// one distributed_property, for each component, must be passed.
// Recommended values of num_components best accepeted by visualization
// tools are either 1, for scalar fields, or DIM, for vector fields.
void xdmf_register_facet_property(xdmf_output *ctx, const char* prop_name,
	size_t num_components, sim_facet_domain *sfd, distributed_property *dp, ...);

//! Register a facet centered property to be written by the XDMF context.
//
// Just like xdmf_register_facet_property(), but with array inputs instead
// of variable argument list.
void xdmf_register_facet_property_array(xdmf_output *ctx,
	const char* prop_name, size_t num_components,
	sim_facet_domain *sfd[], distributed_property *dp[]);

//! Outputs to XDMF+HDF5 the current state of the registered properties.
//
// Parameter time_index is the ordering number used in file names.
//
// Parameter reuse_last_grid may be set to true case the domain
// grid hasn't change since the last call to this same function. If
// so, a previous timestep containing the last written grid will
// be referenced by the newly written timestep.
//
// Either at the first call of the function or when grid has changed,
// reuse_last_grid must necessarily be given as false.
void xdmf_write_timestep(xdmf_output* ctx, size_t time_index, bool reuse_last_grid);

#endif
