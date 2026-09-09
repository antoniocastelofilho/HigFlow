# example2d_DynamicMeshAdapt

Generic benchmark case for 2-D multiphase VOF flows with adaptive mesh
refinement.  The level-set function `func()` in `ns-user-functions-vof.c`
defines the initial interface geometry, and the boundary/initial conditions
in `ns-user-functions-newtonian-gn.c` define the flow type.  Change these
functions to simulate different problems without modifying the AMR kernel
or the solver loop.

## Running

### Zalesak's disk (default)

A rotating notched circle in a uniform flow field — standard test for
VOF interface reconstruction.

```bash
make build_run NP=1 MESH=square_61_br                             \
  IN="mult newt-newt Re=0.005 dt=0.001 numsteps=1001              \
      dts=1.0 dtp=0.01 zalesak_disk adapt=true"                   \
  ENAME=$(date +"%Y%m%d_%H%M")
```

### Shearing droplet

Circular droplet deformed by a lid-driven cavity flow.

```bash
make build_run NP=1 MESH=square_61_br                             \
  IN="mult newt-newt De1=1.0 beta1=0.5 Re=0.005 dt=0.001          \
      numsteps=1001 dts=1.0 dtp=0.01 shearing_droplet adapt=true" \
  ENAME=$(date +"%Y%m%d_%H%M")
```

To switch between configurations, comment/uncomment the appropriate
`func()` in `ns-user-functions-vof.c` and adjust the boundary
conditions in `ns-user-functions-newtonian-gn.c` as needed.

## Fixed issues

### 1. Initial adaptation (step 0) with inconsistent partitioned domains

**Symptom:** The refined mesh at step 0 was not applied correctly — the
simulation used the original unrefined mesh or showed inconsistent behavior.

**Cause:** The step-0 adaptation block called `sd_add_higtree(ns->sdp, root)`
to swap the mesh tree *after* `create_initialize_all_domains()` had already
built the partitioned domains (`psdp`, `psfdu`, `psdmult`) on top of the
original mesh. The partitioned domains kept referencing the old tree; the
refined tree was orphaned in the serial domain while the partitioned domains
operated on the wrong topology.

**Solution:** Use the same approach as dynamic adaptation (step % 5): create a
new solver `ns2`, build the adapted tree via `higflow_make_adapted_tree_params`
(clones and refines), add it to the serial domains of `ns2`, partition `ns2`
with an empty `partition_graph` (valid in serial only), create stencils, and
swap `ns = ns2`. This way DPs, stencils and partitioned domains are created
from scratch on the already-refined tree.

```c
higflow_solver *ns2 = higflow_create();
higflow_load_data_file_names(argc, argv, ns2);
higflow_load_all_controllers_and_parameters_yaml(ns2, myrank);
higflow_set_external_functions(ns2, ...);
higflow_create_domain(ns2, cache, order_center);
higflow_create_domain_multiphase(ns2, cache, order_center, ...);

hig_cell *root = higflow_make_adapted_tree_params(ns, REFINE_THRESHOLDS);
sd_add_higtree(ns2->sdp, root);
sd_add_higtree(ns2->sdF, root);
sd_add_higtree(ns2->ed.mult.sdmult, root);

partition_graph *pg = pg_create(MPI_COMM_WORLD);
pg_set_fringe_size(pg, 5);
load_balancer *lb = lb_create(MPI_COMM_WORLD, 1);
lb_destroy(lb);
higflow_create_partitioned_domain(ns2, pg, order_center);
higflow_create_partitioned_domain_multiphase(ns2, pg, order_center);
higflow_create_stencil(ns2);
higflow_create_stencil_multiphase(ns2);

ns2->par = ns->par;
ns2->contr = ns->contr;
ns = ns2;
```

**File:** `ns-example-2d.c`, block `INITIAL ADAPTATION BASED ON ANALYTICAL
INTERFACE`.

### 2. Manual `num_levels` and duplicated `thresholds` array

**Symptom:** The refinement level count was set manually (`num_levels = 2`)
while the `thresholds` vector had 3 values — dead code and an inconsistency
risk.

**Cause:** Two call sites of `higflow_make_adapted_tree_params` defined
different `thresholds` vectors, and `num_levels` was passed as a separate
parameter with no relation to the actual vector size.

**Solution:** Thresholds were unified into a single global array
`REFINE_THRESHOLDS[]` in `ns-example-2d.c`, terminated by a sentinel `-1.0`.
The function `higflow_make_adapted_tree_params` now computes `num_levels`
internally by walking the array to the sentinel, eliminating the redundant
parameter.

**Files:** `ns-example-2d.c` (global array), `mesh_adapt_function.c` (adapted
function).

### 3. Empty `partition_graph` in adaptations (blocks parallelism)

**Symptom:** The simulation only worked in serial (1 MPI rank). In parallel,
results were
incorrect.

**Cause:** Both adaptation blocks (step 0 and step % 5) created a new, empty
`partition_graph` with `pg_create(MPI_COMM_WORLD)` and discarded the
`load_balancer`
without using it. This empty pg had no MPI neighbor information. When
`dp_sync()` was
called it consulted `pg->filtered_neighbors` (empty) and sent/received nothing —
ghost
(fringe) cells were never synchronized.

**Attempt 1 — Reuse the old `partition_graph`:**
`psd_get_partition_graph(ns->psdp)`. Failed because pg stores pointers to the
old trees in `pg->tree_props` and `pg->neighbors`. After `hig_clone` +
`hig_refine_uniform` the new tree pointers are different; queries with the new
pointers find nothing and `_create_filtered_neighbors` produces an empty
neighbor list.

**Attempt 2 — In-place tree refinement:** `higflow_refine_tree_inplace`
modifies the tree directly without cloning, preserving pointers. Works for step
0 but breaks at step 5 because the in-place-refined tree, when cloned by
`higflow_make_adapted_tree_params`, produces a tree whose interpolation
corrupts memory (heap corruption).

**Attempt 3 — Load balancer with cloned trees:** Use `lb_add_input_tree` +
`lb_calc_partition` to re-partition the adapted trees. Failed with
`MPI_ERR_TRUNCATE` because the load balancer expects complete trees, not
partial ones (each rank holds only its local portion).

**Current solution:** Empty `partition_graph` (`pg_create` +
`pg_set_fringe_size` without `lb_calc_partition`). Valid in serial (1 MPI rank)
only. In parallel, `dp_sync` does not work because no neighbors are configured.
Parallel support for dynamic mesh adaptation remains a future task.

```c
partition_graph *pg = pg_create(MPI_COMM_WORLD);
pg_set_fringe_size(pg, 5);
load_balancer *lb = lb_create(MPI_COMM_WORLD, 1);
lb_destroy(lb);
```

**File:** `ns-example-2d.c`, both adaptation blocks.

### 4. `fracvol` filter with per-cell `printf`

**Symptom:** Thousands of "Step 5" lines in the simulation log, polluting output
and
degrading performance.

**Cause:** The empirical `fracvol` filter in `higflow_interpolate_viscosity`
printed
`printf("Step %d\n", ns->par.step)` for every cell in the new mesh during
interpolation.

**Solution:** Kept as-is per user decision (the filter is needed for
volume-fraction
sharpening, but the `printf` can be removed if output volume becomes an issue).

**File:** `ns-example-2d.c`, function `higflow_interpolate_viscosity`.

### 5. Lid-driven cavity velocity

**Change:** Top wall velocity (bc1) was reduced from `u=2.0` to `u=1.0` and a
small
perpendicular component `v=0.1` was added. The velocity initial condition was
updated
from `u=2.0*y` to `u=y` for consistency with the boundary condition.

Physics note: Couette flow (u=y, v=0) is an exact Navier-Stokes solution with
uniform
pressure, fully compatible with the Neumann pressure BCs (dp/dn=0) in the YAML.

**Files:** `ns-user-functions-newtonian-gn.c`, functions `get_boundary_velocity`
(case 1)
and `get_velocity`.

### 10. Boundary higtree refinement (BC AMR)

**Problem:** Boundary higtrees (files `mesh/__using/bc/ch-bc-N.amr`) have a
fixed
resolution defined on disk. When the internal mesh is refined by AMR (interface
cells
advance from level 2 to level 3, for example), the boundary cells keep their
original
coarse resolution — even when the interface touches the domain boundary, that
region of
the boundary is not refined.

**Cause:** `higflow_initialize_boundaries_yaml` reads each `.amr` file, builds
the
boundary higtree and creates the `sim_boundary` objects directly, without
comparing
against the resolution of the adjacent internal cells.

**Solution:** A post-read hook in `higflow/src/hig-flow-bc.c` lets any caller
inject a
refinement callback:

```c
// Register a callback invoked on each boundary higtree immediately after it
// is read from disk, before the sim_boundary is created.  Pass NULL to disable.
void higflow_set_bc_refine_hook(void (*hook)(hig_cell *bc_root, int bc_id));
```

The refinement logic itself lives in `ns-example-2d.c` (case-specific code stays
in the
case folder):

```c
static hig_cell *_bc_domain_root = NULL;

static void _refine_bc_tree(hig_cell *bc_root, int bc_id) {
    // Walk bc_root leaves; for each, probe the adjacent internal cell.
    // If the internal cell is deeper (finer), split the BC leaf in the
    // direction parallel to the boundary and restart the walk.
}
```

Usage around `higflow_initialize_boundaries_yaml`:

```c
_bc_domain_root = root;               // adapted internal domain root
higflow_set_bc_refine_hook(_refine_bc_tree);
higflow_initialize_boundaries_yaml(ns2);
higflow_set_bc_refine_hook(NULL);
_bc_domain_root = NULL;
```

**Algorithm (`_refine_bc_tree`):**

1. Iterate leaves of the boundary higtree.
2. For each leaf, find the adjacent internal cell by probing a point just
   inside the domain (offset ε in the inward-normal direction) via
   `hig_get_cell_with_point`.
3. Compare levels: if the internal cell is deeper than the boundary leaf,
   call `hig_refine_uniform` with `{1,2}` for left/right boundaries or
   `{2,1}` for top/bottom boundaries, then restart the iteration.
4. Repeat until the boundary matches the internal mesh resolution everywhere.

**Boundary interpolation correctness:** After BC higtree refinement,
`higflow_interpolate_all_bcs` re-evaluates velocity values via
`get_boundary_velocity` at each refined cell center — correct for any
analytical BC.  Neumann and fixedValue pressure BCs are skipped by the
interpolator.

**Validation:** `print_bc_cell_counts` reports the number of cells per
boundary before and after each adaptation.  With a 2-level adapted internal
mesh (61 → 122 → 244 cells per side), boundaries adjacent to the interface
should report extra cells in that region.

**Files:** `higflow/src/hig-flow-bc.c` (hook mechanism only),
`higflow/src/hig-flow-bc.h` (declaration), `ns-example-2d.c`
(refinement logic + `#if ADAPT_ENABLED` call sites).

## Current status

- **Serial (NP=1):** Working.  Adaptation at step 0 and every
  `ADAPT_FREQ` steps with
  mass conservation ~0.05%.
- **Parallel (NP>1):** Not working.  The empty `partition_graph` has no MPI neighbors
  configured, so `dp_sync` cannot communicate ghost cells between ranks. Pending
  future
  implementation.

## Future optimisations (pending)

### 6. Optimise seed search in coarsening (O(N²) → O(N))

**Problem:** In `mesh_adapt_function.c`, coarsening (lines 169-251) checks each
parent
cell against ALL interface seeds:

```c
for(int s=0; s<seed_count; s++) {
    real d2 = dist_sq_c(center, seeds[s].center);
    if (d2 < min_d2) min_d2 = d2;
}
```

For 1000 seeds × 500 parents × 4 children × 4 passes = **8M distance
computations**.
Complexity: O(seed_count × parent_count) = O(N²) near the interface.

**Proposed solution:** Spatial hashing (2D bin grid):

```c
// Divide domain into bins of size max(thresholds) + coarsen_hys
// Each bin stores indices of seeds in that region
// For each child, query only the same bin + 8 neighbors (9 bins max)
// Reduces from O(N²) to O(N) with O(seed_count) memory overhead
```

**File:** `mesh_adapt_function.c`, function `higflow_make_adapted_tree_params`,
coarsening step.

### 7. Merge the 8 AMR interpolations into a single loop

**Problem:** In the AMR block (`ns-example-2d.c`), 8 interpolation functions
iterate the
full mesh sequentially:

1. `higflow_interpolate_velocity` — all facets
2. `higflow_interpolate_viscosity` — all cells (2× `compute_value_at_point`)
3. `higflow_interpolate_density` — all cells
4. `higflow_interpolate_pressure` — all cells
5. `higflow_interpolate_bc_for_velocity` — boundary cells
6. `higflow_interpolate_bc_for_pressure` — boundary cells

Each function opens its own `celliterator`/`facetiterator`, does `mp_lookup` per
cell,
and closes the iterator. Iterator creation/destruction overhead and poor cache
locality
make this ~30% slower than a single pass.

**Proposed solution:** A `higflow_interpolate_all` function that performs a
single
iteration over the new mesh:

```c
void higflow_interpolate_all(higflow_solver *ns, higflow_solver *ns2) {
    higcit_celliterator *it = sd_get_domain_celliterator(sdm2);
    while (!higcit_isfinished(it)) {
        hig_cell *c = higcit_getcell(it);
        int clid = mp_lookup(mp2, hig_get_cid(c));
        Point ccenter;
        hig_get_center(c, ccenter);
        real fracvol = compute_value_at_point(...);
        real visc    = compute_value_at_point(...);
        dp_set_value(ns2->ed.mult.dpfracvol, clid, fracvol);
        dp_set_value(ns2->ed.mult.dpvisc,    clid, visc);
        // ...
    }
}
```

**Expected gain:** ~30% reduction in AMR block time.

**Files:** `ns-example-2d.c`, `higflow_interpolate_*` functions, AMR block.

### 8. Incremental AMR (Basilisk style)

**Problem:** Currently AMR **destroys and recreates the entire solver** on every
adaptation:

- `higflow_destroy(ns)` — frees domains, stencils, DPs, PETSc KSP
- `higflow_create()` — allocates new solver
- `higflow_create_domain()` + `higflow_create_domain_multiphase()` — new domains
- `higflow_create_partitioned_domain()` — partitioning
- `higflow_create_stencil()` — stencils
- `higflow_create_distributed_properties()` — new DPs
- `higflow_create_solver()` — new PETSc KSP
- 8 full-mesh interpolations to transfer data

Total cost ≈ **5-10× a normal solver step**. In Basilisk, AMR typically costs
<5% of
one integration step.

**Proposed Basilisk-style solution:**

1. **Refine/coarsen in-place** — use `hig_refine_uniform` / `hig_merge_children`
   directly on the active tree without cloning.  `higflow_refine_tree_inplace` in
   `mesh_adapt_function.c` is partially implemented but unused.

2. **Update only affected cells** — after refine/coarsen, only new and removed cells
   need initialisation.  Transfer values by local interpolation
   (edge-by-edge), not
   full-mesh.

3. **AMR every step but cheap** — instead of adapting 100% of the
   interface every N
   steps, adapt 1-2% of cells per step (those that crossed the refinement
   threshold).
   The amortised cost is negligible and the mesh stays optimal at all times.

4. **Reuse the PETSc KSP** — `KSPSetOperators` can reinitialise an
   existing KSP instead
   of destroying and recreating it. PETSc matrices can be refilled with new
   topology
   without reallocating the KSP object.

5. **Local wavelet criterion** (like Basilisk) — instead of a global
   seed search, use
   the volume-fraction gradient as a per-cell error estimator:

   ```c
   // If |∇f| > threshold and level < max_level  → refine
   // If |∇f| < threshold/coarsen_factor and level > min_level → coarsen
   ```

**Expected gain:** ~80% reduction in AMR cost (from 5-10× steps to <1× step).

**Files:** `mesh_adapt_function.c`, `ns-example-2d.c`.

### 9. Fix `create_velocity_copy` — uninitialised pointer

**Problem:** In `utilities-2d.c:1031`, `create_velocity_copy` declares
`distributed_property **u_copy` without `malloc`, then writes to `u_copy[dim]`.
This is undefined behaviour (stack corruption).

```c
distributed_property **create_velocity_copy(higflow_solver *ns){
    distributed_property **u_copy;  // NEVER allocated!
    for(int dim=0; dim<DIM; dim++){
        u_copy[dim] = psfd_create_property(ns->psfdu[dim]);  // writes to garbage
    }
    return u_copy;
}
```

**Solution:** Add `malloc`:

```c
distributed_property **u_copy = malloc(DIM * sizeof(distributed_property*));
```

**File:** `utilities-2d.c`, function `create_velocity_copy`.

## Post-processing

### `snapshot_last.py`

A ParaView Python script (run with `pvpython`) that opens the last
VTK timestep of the most recent simulation, displays the `FracVol`
field with a blue-green-orange colormap, and saves a PNG screenshot.

```bash
pvpython snapshot_last.py                    # latest VTK in ./output/
pvpython snapshot_last.py output/MyRun       # or a specific output dir
pvpython snapshot_last.py --step 10 --size 1600x1200
```

Arguments:

- `output_dir` — root folder to scan for `*/vtk/` directories (default: `./output/`)
- `--step N` — render timestep N instead of the last one
- `--size WxH` — image resolution (default: `1200x1200`)

The PNG is saved next to the script with the simulation name and
step number, e.g. `mySimulation__step36.fracvol.png`.
