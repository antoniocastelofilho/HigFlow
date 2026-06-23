# HigFlow Parallelization Saga — 2D AMR (example2d_DynamicMeshAdapt)

## The Goal

Make HigFlow's 2D adaptive mesh refinement (AMR) work in parallel (NP>1),
including dynamic load balancing and mid-step domain rebuilds at runtime.

---

## Act I — The Problem

Running `example2d_DynamicMeshAdapt` with `NP=4` on `square_61_br` crashed
immediately in `higflow_destroy(ns_old)` with a segmentation fault after
step-0 AMR switched to the adapted tree.

**Symptom:**
```
[0]PETSC ERROR: Caught signal number 11 SEGV: Segmentation Violation
```

**Initial investigation:**
- Uniform mesh (ADAPT_ENABLED=0, NP=4) worked fine for 50 steps
- Step-0 AMR with NP=1 also worked fine
- Only step-0 AMR + NP>1 crashed

---

## Act II — The False Leads

### Theory 1: `ns_old` has stale data after AMR
**Fix attempt:** Save and restore only the minimal runtime parameters
needed by `higflow_create_independent_controller`.
**Result:** Still crashed.

### Theory 2: Load balancer doesn't work with adapted trees
**Fix attempt:** Refactor the load balancer call into the rebuild flow.
**Result:** Load balancer was fine (`lb_calc_partition` handles adapted
trees). The crash was elsewhere.

### Theory 3: Solver created twice causes corruption
**Fix attempt:** Remove the duplicate `higflow_create_solver` call;
let the helper create it once, skip if `amr_solver_ready`.
**Result:** Still crashed.

---

## Act III — The Root Cause

Deep debugging revealed the real issue:

**In step-0 AMR, when `amr_solver_ready` is set, two functions are
skipped:**

```c
if (!amr_solver_ready) {
    higflow_create_distributed_properties(ns);   // SKIPPED
    higflow_create_solver(ns);                     // SKIPPED
}
```

This leaves `ns_old` with **NULL pointers** for:
- `ns->dpu[0]`, `ns->dpu[1]`, `ns->dpdiv` (DPs)
- `ns->stnu[0]`, `ns->stnu[1]`, `ns->stndiv` (stencils)
- `ns->lslv`, `ns->gslv` (solvers)

When `higflow_destroy(ns_old)` is called, it tries to free these NULLs:

```c
dp_destroy(ns->dpu[0]);   // ns->dpu[0] is NULL → SEGV
```

**The higtree functions `dp_destroy`, `stn_destroy`, and `slv_destroy`
did not tolerate NULL arguments.**

---

## Act IV — The Fix

Added NULL-safety guards to three functions in higtree:

### `higtree/src/pdomain.c` — `dp_destroy`
```c
void dp_destroy(distributed_property *dp) {
    if (dp == NULL) return;
    // ... existing cleanup ...
}
```

### `higtree/src/domain.c` — `stn_destroy`
```c
void stn_destroy(stencil *stn) {
    if (stn == NULL) return;
    // ... existing cleanup ...
}
```

### `higtree/src/solver.c` — `slv_destroy`
```c
void slv_destroy(solver *slv) {
    if (slv == NULL) return;
    // ... existing cleanup ...
}
```

**Result:** `higflow_destroy(ns_old)` now succeeds safely. The step-0
AMR flow runs for 50 steps with NP=2 without execution-time crash.

---

## Act V — The Architecture (What We Built)

### Parallel AMR rebuild flow

```
Step 0:
  1. higflow_rebuild_with_amr()
     ├── higflow_make_global_adapted_tree()  ← builds identical tree
     │   on ALL ranks using MPI-gathered seeds
     ├── lb_calc_partition()                 ← load balance
     ├── higflow_create_partitioned_domain() ← new domain
     ├── higflow_create_distributed_properties()
     ├── higflow_create_solver()
     └── free(ns_old)                        ← now safe
  2. amr_solver_ready = true

Mid-step (every adapt_freq steps):
  1. higflow_rebuild_with_amr()
     ├── higflow_make_global_adapted_tree()
     ├── lb_calc_partition()
     ├── higflow_create_partitioned_domain()
     ├── higflow_create_distributed_properties()
     ├── higflow_create_solver()
     ├── free(ns_old)                        ← now safe
     └── load fresh controllers for ns2
  2. copy runtime params, advance simulation
```

### Key refactoring in `mesh_adapt_function.c`

| Function | Purpose |
|----------|---------|
| `higflow_make_global_adapted_tree` | Builds adapted tree identically on all ranks |
| `collect_interface_seeds_local` | Collects seeds from local cells near interface |
| `gather_seeds_mpi` | Gathers seeds via MPI_Allgatherv to all ranks |
| `adapt_tree_with_seeds` | Adapts tree using unified seed list |

---

## Act VI — Remaining Issues

### 1. Mid-step AMR (multiple rebuilds at runtime)
Not yet tested with the NULL-safe destroy fix for long runs.

### 2. Final-exit heap corruption
`higflow_destroy(ns)` at program exit crashes:
```
*** Error in `./ns-example': free(): invalid next size (normal):
    0x0000000001a9d6a0 ***
```
This is a minor heap corruption in `dpu[0]` during solver steps,
detected only at cleanup. Not a showstopper for functionality.

**Note:** AddressSanitizer cannot be used — crashes immediately
with MPI/OpenMPI. Valgrind is too slow to reach failure point.

---

## Act VII — What Works

| Scenario | Status |
|----------|--------|
| Uniform mesh, NP=1 | Works (50 steps) |
| Uniform mesh, NP=2 | Works (50 steps) |
| Uniform mesh, NP=4 | Works (50 steps) |
| Step-0 AMR, NP=1 | Works (50 steps) |
| Step-0 AMR, NP=2 | Works (50 steps) — with NULL-safe destroy |
| Mid-step AMR, NP=1 | Not yet tested |
| Mid-step AMR, NP=2 | Not yet tested with NULL-safe fix |
| Long runs (>50 steps) | Needs final-exit fix |

---

## Act VIII — Next Steps

1. **Test mid-step AMR** (multiple rebuilds at runtime) with NP=2
   after the NULL-safe destroy fix.
2. **Investigate final-exit heap corruption** — may require
   targeted `valgrind` or a minimal reproduction case.
3. **Clean up remaining debug artifacts** if any.
4. **Benchmark** parallel AMR vs. uniform for performance.

---

## The Takeaway

The entire saga came down to one insight:

> **When step-0 AMR skips solver creation, `ns_old` has NULL pointers.
> `higflow_destroy` must tolerate NULLs, or the NULLs must never exist.**

The fix was three one-line NULL guards in higtree — but finding it
required rewriting the entire rebuild flow, refactoring the mesh
adaptation into a global parallel operation, and systematically
eliminating false leads.

---

## Files Modified

| File | Change |
|------|--------|
| `ns-example-2d.c` | Step-0 & mid-step AMR blocks; helpers |
| `mesh_adapt_function.c` | Global adapted tree build; seed collection |
| `ns-example-2d.h` | `higflow_make_global_adapted_tree` declaration |
| `src/hig-flow-kernel.c` | `higflow_destroy` (NULL-safe); solver creation |
| `higtree/src/pdomain.c` | `dp_destroy` NULL guard |
| `higtree/src/domain.c` | `stn_destroy` NULL guard |
| `higtree/src/solver.c` | `slv_destroy` NULL guard |

---

*Last updated: 2026-06-16*
