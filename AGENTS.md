# HigFlow — AGENTS.md

## Environment setup (must source before build/run)

- `source varsrc` — sets `HIGTREE_DIR`, `HIGFLOW_DIR`, `PETSC_DIR`, `PETSC_ARCH`

## Build (traditional Make)

- **Order matters:** `higtree` → `higflow` → example
- `make DIM=2` or `make DIM=3` (required variable; targets
  `libhig$(DIM)d.a` / `libhigflow$(DIM)d.a`)
- Debug: `make DEBUG=1 DIM=2`
- Sanitize: `make SANITIZE=1 DIM=2` (address sanitizer)
- Compiler: `mpicc`, C standard `-std=gnu99`, C++ standard `-std=gnu++11/14`
- Release: `-O3 -flto=8 -fno-fat-lto-objects -march=native -mtune=native`
- From the example dir, rebuild deps with:
  `make rebuild_higtree`, `make rebuild_higflow`, `make rebuild_all`

## Build (CMake)

```bash
mkdir build && cd build
cmake .. -Wno-dev -Ddim=2           # or -Ddim=3
make && make install                # installs to <repo>/install/ by default
```

CMake debug: `-Ddebug=1`. Install prefix: `-Dprefix=PATH`.

## Run pattern (reference: example2d_DynamicMeshAdapt)

```bash
cd higflow/example2d_DynamicMeshAdapt
make build_run NP=4 IN="mult newt-newt De1=1.0 Re=0.005 dt=0.001 numsteps=1001 shearing_droplet adapt=true" MESH=square_61_br ENAME=$(date +"%Y%m%d_%H%M")
```

- `build_run` = `make clean` + `make` + `make run`
- `IN` — YAML config overrides parsed by `setup.py`; defaults should live
  in `input/load.par.contr.yaml`
- `NP` — MPI processes (serial=1)
- `MESH` — selects mesh from `mesh/$(MESH)/`; runs `setup_mesh`
  to populate `mesh/__using/`
- `ENAME` — date suffix appended to `OUTNAME` for unique output folders
- `OUTNAME = $(shell python3 setup.py $(INPUT))__$(MESHNAME)`
  → `output/$(OUTNAME)/{vtk,save,res}/`
- `adapt=true` enables AMR (rewrites `ADAPT_ENABLED` macro in header)
- `adapt_freq=N` controls AMR interval (default 5)
- Restart: `RESTART=output/previous_outname`
- Thread binding: `THREAD=1`

## 3D cases

```bash
cd higflow/example3d_RisingDropCap
make NP=1 run                    # uniform mesh, 20×40×20 cells
```

- AMR: toggle `ADAPT_ENABLED 1/0` in `ns-exemple-3d.h`
- Mesh files: `amrs/ch-{d,bc-N}.amr` (6 BCs for 3D)
- ⚠️ **Gravity:** the multiphase solver hardcodes gravity on `dim == 1`
  (Y direction).  Domains must be oriented with **Y as the tall dimension**
  so the bubble can rise the full height.  The OF reference (gravity in Z)
  is adapted by swapping Y ↔ Z in the mesh.
- Bubble: center (0.5, 0.5, 0.5), radius 0.25, rising in +Y (gravity in -Y)
- Setup mesh (first time only):

  ```bash
  mkdir -p mesh/__using/{domain,bc}
  cp amrs/ch-d.amr mesh/__using/domain/ch-d.amr
  for i in 0 1 2 3 4 5; do cp "amrs/ch-bc-$i.amr" "mesh/__using/bc/ch-bc-$i.amr";
  done
  ```

- Parameters: Re=35, Ca=3.57, Fr=1.0, ρ_bubble/ρ_amb=0.001, μ_bubble/μ_amb=0.01
- 3D AMR solver: `gmres(200) + hypre boomeramg` (set in `_AMR/Makefile`)

## Dam break cases

```bash
# AMR version (recommended)
cd higflow/example2d_DamBreak_AMR
make build_run NP=1 MESH=dambreak_160x100 \
  IN="mult newt-newt Re=3130000 Ca=0.045 dt=0.001 numsteps=8000 dts=0.05 dtp=0.005" \
  ENAME=$(date +"%Y%m%d_%H%M")

# Uniform mesh version
cd higflow/example2d_DamBreak
make build_run NP=1 MESH=dambreak_160x100 \
  IN="mult newt-newt Re=3130000 Ca=0.045 dt=0.001 numsteps=8000 dts=0.05 dtp=0.005" \
  ENAME=$(date +"%Y%m%d_%H%M")
```

## Rising drop / bubble cases

```bash
# 2D AMR
cd higflow/example2d_RisingDropCap_AMR
make build_run NP=1 MESH=bubble_50x100 \
  IN="mult newt-newt Re=35 Ca=3.57 dt=0.001 numsteps=10000 dts=0.1 dtp=0.05" \
  ENAME=$(date +"%Y%m%d_%H%M")

# 2D uniform
cd higflow/example2d_RisingDropCap
make build_run NP=1 MESH=bubble_50x100 \
  IN="mult newt-newt Re=35 Ca=3.57 dt=0.001 numsteps=10000 dts=0.1 dtp=0.05" \
  ENAME=$(date +"%Y%m%d_%H%M")

# 3D uniform (AMR via #define in header)
cd higflow/example3d_RisingDropCap
make NP=1 run
```

## Multi-instance job submission

```bash
cd higflow/example2d_DynamicMeshAdapt
python3 master_prog/master.py     # queue manager for batch runs
```

- Monitors `master_prog/shared/submission/` for `.sh` scripts
- Executes them one at a time with cooldown

## Mesh setup

- Generate: `python3 mesh/gen_mesh_channel.py`
  (creates AMR-format meshes)
- Activate: `make setup_mesh MESH=name` (copies `mesh/name/` → `mesh/__using/`)
- Format (2D): `xmin xmax ymin ymax` / `levels` / `dx dy npatch` / `xs ys w h`
- Format (3D):
  `xmin xmax ymin ymax zmin zmax` / `levels` / `dx dy dz npatch`
  / `xs ys zs w h d`

## Dependencies

- **Required:** PETSc 3.14.0 (tarball in `bibliotecas/`), MPI, HDF5,
  GLib (`glib-2.0`), libfyaml, Trilinos/Zoltan, Hypre
- **pkg-config:** glib-2.0, hdf5, libfyaml — all needed for both
  `--cflags` and `--libs`
- **Python:** `ruamel.yaml` (`pip3 install ruamel.yaml`) — required
  by `setup.py`
- Install scripts: `install_higflow_ubuntu22.sh` or `install_higflow_arch.sh`

## Docker

```bash
# PETSc + deps (slow build)
docker build . -t ubuntu_petsc3.14:v01 -f conteiner/Dockerfile.petsc
# compile higtree + higflow libraries
docker build . -t higflow:v01         -f conteiner/Dockerfile
```

## Apptainer / Singularity

```bash
apptainer build higflow.sif docker-daemon://higflow:v01
```

- Use `--no-home` to preserve container's `/home/hig_user/` files (env setup):

  ```bash
  apptainer exec --no-home --cleanenv higflow.sif <command>
  apptainer shell --no-home --cleanenv higflow.sif
  ```

- Env vars are stored in `/.singularity.d/env/99-higflow.sh` inside the image
- `--cleanenv` strips host env; combine with `--no-home`

## Testing

- `./TUM_2D` — build everything 2D + run Newt and Newt_contraction examples
- `./TUM_3D` — build 3D + run example3d_complex
- `./TUM_VISC` — build 2D + run Gptt viscoelastic
- `./TUM_TEST_ALL_CASES` — run all 8 2D examples sequentially
- `./TUM_CLEAR_ALL` — clean all example build artifacts
- Benchmark: `./benchmark_speed.sh` (in `example2d_DynamicMeshAdapt/`)

### The two suites, and what each one covers

```bash
set -a; . ./varsrc; set +a          # from the REPOSITORY ROOT; varsrc uses $(pwd)
python3 ci/run_suite.py             # 33 runs: the HiGFlow examples, 2D and 3D
python3 ci/run_higtree_tests.py     # 104 cases: the HiGTree library itself
```

`run_suite.py` **builds** higtree and never tests it — the library is a
prerequisite of the examples. Every `domain.c` defect fixed between 15 and
18/09/2026 lived there and was found by a 3D run aborting, not by a test. That
is what `run_higtree_tests.py` exists for, and it asserts **value, not form**.

It also enforces the **Mesh contract** (`higtree/src/hig-mesh-contract.h`): the
driver reads the clauses from the header itself and requires each to have a
passing case. A clause whose test disappeared prints `SEM TESTE QUE RODE`; one
whose test failed prints `REPROVADA`. Both count as suite failures — a clause
without a test is a comment, not a clause.

### Second mesh backend (t8code), optional

```bash
python3 ci/run_higtree_tests.py --t8code bibliotecas/t8code/install
```

With the flag, the 18 contract clauses are verified for **two** mesh
implementations instead of one (142 cases). Without it nothing of t8code is
built and the suite is identical to before.

t8code is not packaged anywhere and must be built from source — the how, the
three non-obvious integration constraints, and what is still missing are in
**`higtree/tests/t8code/README.md`**. Read it before touching that directory.

## README per case

Every new test case must have a `README.md` describing:

- What problem it solves / what physics it models
- Reference paper or benchmark
- Run command (`make build_run …`)
- Key parameters (Re, Ca, ρ ratios, mesh)
- Output structure

## Documentation language

All documentation must be written in English — README files,
code comments, commit messages, function docstrings.  Portuguese is not
used in any written artefact checked into the repository.

## Python scripts

- All Python scripts must produce **zero pylint warnings or errors**.
- Structure: top-level functions + `if __name__ == "__main__":` guard.
- Late imports (e.g. `paraview.simple`) go inside functions with
  `# pylint: disable=import-outside-toplevel`.
- Every function needs a docstring.

## Saving instructions ("lembre")

When the user says **"lembre"** in a request, save the instruction to the
auto-memory system (`~/.claude/projects/.../memory/`) **and** add it to this
file if it is a standing code or workflow convention.

## Code style guide

- **Language:** C99 (`-std=gnu99`), single `return` per function where possible
- **Line length:** max 80 characters
- **Indentation:** 4 spaces, no tabs
- **Braces:** Allman style (`{` on next line) for function definitions;
  K&R style for control flow (`if () {`)
- **Documentation:** Doxygen `/*! \brief … */` on every public function
- **Comments:** English only.  Explain *why*, not *what* (the code is the what)
- **Variables:** descriptive names.  When they can't be (e.g. `p_par`),
  add a comment at the declaration
- **Headers:** include guard `#ifndef NS_EXAMPLE_2D_H` / `#define …`
- **Globals:** group in `/*! @{ */` / `/*! @} */` blocks, document each
- **Linting:** run `clang-format --dry-run -Werror --style=gnu *.c` to check
  C formatting; MD files are also linted — use `markdownlint-cli2` with the
  config at `.markdownlint.jsonc`.  Run before every commit:

  ```bash
  markdownlint-cli2 AGENTS.md **/README.md
  ```

  Table style: compact (no spaces around pipe characters) — use
  format like:

  ```markdown
  |Cell|Content|
  |-----|-------|
  |val1|val2|
  ```

  not spaced like `| Cell | Content |`.

## Committing

- **Before every commit** run `git status` and inspect the list of files.
  Watch for:
  - `*.o`, `*.a`, `*.so`, executables (`ns-example`, `ns-exemple-3d`) —
    these are build artifacts, never commit them.
  - `core.*`, `vgcore.*` — core dumps, never commit.
  - `*.mp4`, `*.avi`, `*.log`, `*.vtk`, `*.png` — simulation outputs,
    add to `.gitignore` if they appear.
  - Large generated files (check with `git ls-files --cached -z
    \| xargs -0 ...`).
- When adding new directories, use `git add <dir>/<specific-file>` instead of
  `git add <dir>/` to avoid picking up build artifacts and media files.
- **Commit messages** must be written in English with a detailed description
  of the root cause, what was done, and verification results.  Follow the
  structure: summary line, blank line, `Root cause`, `What was done`,
  `Verification` sections.

- **Dimension:** `DIM=2`/`DIM=3` via make variable or `-Ddim=` CMake
  variable. Source files are DIM-agnostic; compiled with `-DDIM=$(DIM)`.
- **PETSc env:** `PETSC_DIR` + `PETSC_ARCH` must be set (from `varsrc`
  or the container env)
- **Static libs** (`.a`) are gitignored — rebuild from source after clone
- **`bibliotecas/`** gitignored — contains PETSc tarball + libfyaml source
- **`**/output/**`** gitignored — simulation results
- **Setup scripts** use `ruamel.yaml` to rewrite
  `input/load.par.contr.yaml`.  The `ruamel` package must be installed
  system-wide (not user-local) for Apptainer compatibility.
- **HDF5 headers** on Ubuntu: at `/usr/include/hdf5/openmpi/`,
  resolved by `pkg-config --cflags hdf5`
- **Viscoelastic models:** oldroyd_b, giesekus, lptt, gptt, fene_p, e_fene
- **Electroosmotic models:** pnp, pb, pbdh
- **BC types:** cavity, channel
- **VTK output** viewable in ParaView

## Markdown formatting: tables

Use **compact style** for all tables — no spaces around pipe characters:

```markdown
|Cell|Content|
|-----|-------|
|val1|val2|
```

NOT:

```markdown
| Cell | Content |
|------|---------|
| val1 | val2   |
```

This avoids MD060/table-column-style linter warnings.

## Reading traps in this tree

Three tools answer a question *close to* the one you meant. Each cost real time on
2026-09-16, and each produced a plausible wrong answer rather than an error.

### `git diff` does not tell you what is left to commit

Several sessions share one working tree, so they share the index. A `git add`
without a pathspec leaves a stale snapshot there, and `git status --short` then
marks files `MM`. From that point on, plain `git diff` compares the tree against
that stale index and shows *already committed* work as if it were pending.

|Question|Command|
|-----|-------|
|What is left to commit?|`git diff HEAD`|
|Which files?|`git diff HEAD --name-only`|
|Commit without touching the shared index|`git commit -- <paths>`|

Never commit by staging into the shared index: it can **revert** another session's
committed work, which is worse than picking it up.

### A pipeline hides the exit status you care about

`cmd | tr -d '\0'` makes `$?` the status of `tr`. A run killed by `timeout` reports
success. Use `set -o pipefail` or `${PIPESTATUS[0]}`.

### `pgrep -f` and `pkill -f` match their own command line

`pgrep -f run_suite.py` matches the shell that runs it, so an
`until pgrep …; do sleep; done` loop never ends. Worse, `pkill -f 'while pgrep'`
kills that shell (exit 144) — once taking a nearly finished suite with it.

The bracket form, `run_suite[.]py`, does stop the pattern from matching itself —
measured, in isolation. It is still not enough, for two reasons that bit us both:

- **A second, unbracketed occurrence on the same command line.** In a single
  `bash -c` that waits for X and then runs X, the literal from the "then runs"
  half sits in the same cmdline, and that is what matches. Brackets protect the
  pattern, not the command that contains it.
- **Another session's process.** With several sessions on one machine, "is any
  `run_suite.py` alive" does not distinguish yours from theirs, and
  `pgrep … || break` keeps waiting while *any* one exists.

So do not wait by process name at all. Capture the PID at launch, chain with `&&`
so no wait is needed, or poll for a marker the job itself writes.

### Every rank opens the same file and they clobber each other

Instrumentation that writes a dump from inside the solver runs on *every* rank. With
a fixed filename, the ranks open the same path with `"w"` and overwrite one another,
so a two-rank run leaves roughly one rank's worth of data:

    ustar_pre_1.txt   np=1  24.5 MB
    ustar_pre_1.txt   np=2  12.2 MB     <- not half the work, half the data

Comparing that against a complete serial dump finds "divergence" across every row the
surviving rank never wrote — confirming whatever you hoped to confirm, for the wrong
reason. Put the rank in the name (`dump_%d_r%d.txt`, from `MPI_Comm_rank`) and
concatenate afterwards.

The habit that caught it: **state the check before looking at the result.** Having
already committed to "verify the matched-key count before reading any difference"
made a half-sized file impossible to wave away. A control decided after seeing the
number is worth much less — by then you know which answer you want.

### A reference recorded from a broken configuration hides the break forever

`example-3d.load.bc.yaml` pointed `bc10` at `mesh-channel3D-bc-101.amr` — one digit
too many. One patch was loaded twice, another never, and a one-cell hole was left in
an obstacle wall. **97% of the inflow came through that hole instead of the channel
inlet**: the case was solving a different problem, and it *passed the suite*, because
the reference had been generated with the hole in place.

That is the worst failure mode here, because it is self-sealing: once the reference
encodes the defect, the defect becomes indistinguishable from correct behaviour for
every test that consults that reference. The antidote is not to distrust references —
it is to keep at least one criterion that does not depend on one:

|Criterion|What it caught|
|-----|-------|
|Flux through cross-sections|the hole, and a mass leak in example2d_Newt_contraction|
|Symmetry in a symmetric case|the same contraction leak, for free, in the diff|
|Agreement across process counts|the whole partitioned-stencil investigation|

Before trusting such an instrument, plant a defect and check it is found *in the
right place*. The planted defect needs an independently known signature — position,
magnitude, or both — otherwise "it detected something" is not evidence. The flux
measurement above was run twice, the second leg deliberately restoring the missing
patch at a known x, and only then did the clean leg's 1.4e-13 mean anything.

### One number can hide the profile that explains it

That same flux measurement, summarised as a mean over all 130 sections, reported
"100% deviation" for *both* legs and would have refuted a working instrument. The
mean included sections the flow had not reached in two steps. The profile showed the
structure at once: flat at the correct value, a factor-32 step exactly at the plane
of the missing patch, then the transient front. Look at the profile before reducing
it to a number.

### VTK output is partition-dependent; aggregates hide it

In a multi-block domain the VTK writer interpolates velocity at cell *corners*
(`compute_facet_value_at_point` → `sfd_get_stencil`), and the least-squares support
is capped at `maxpts` (120 for order 2 in 3D). At a block interface there are more
candidates than that at the same distance, so **which ones make the cut depends on
the order the trees are visited — which depends on the partition.**

Measured on `example3d_complex`, np=1 against np=2, same fields, two instruments:

|Path|Points compared|Max difference|
|-----|-------|-------|
|`ns->dpu` per facet, coordinate as key|458,400|1.0e-10 (the dump's own floor)|
|VTK nodal values|190,729|1.78e-02 — **0.56% of scale**|

The solution agrees to the limit of the instrument; the *output* does not. At the
worst node the two supports share 114 of 120 points, and the six that differ are all
3.08 cells away — well inside the 5-cell fringe. Nothing is missing: the tie at the
cutoff is broken differently.

**The suite does not catch this, and cannot.** It compares min/max/mean per field.
Those aggregates are blind to a defect localised in a few hundred corner nodes:

|Component|np=1 min / max / mean|np=2 min / max / mean|Relative difference|
|-----|-------|-------|-------|
|u|-3.180224 / 1.615838 / -0.085414|-3.180224 / 1.615838 / -0.085416|0, 0, 1.9e-05|
|v|-1.617355 / 1.617351 / -0.014271|-1.617355 / 1.617351 / -0.014271|0, 0, 4.6e-06|

Min and max agree *exactly* — the 236 divergent nodes never reach the extremes (their
u spans [-2.586, 1.528] against a global [-3.180, 1.616]) — and the rest dilutes in
the mean. A green suite at tolerance 1e-5 means three statistics agree, **not** that
the fields do.

Consequences:

- **Do not compare VTK between different `np`.** To compare decompositions, dump
  `ns->dpu` per facet with the coordinate as key.
- **A reference recorded from VTK in a multi-block domain is only valid for the
  decomposition that recorded it.** `ci/run_suite.py` generates at np=1 by design and
  checks np=2,3,4 against it; that check passes on aggregates, not on fields.

This is inherent to the writer, not a bug introduced by a change: the support cut is
order-dependent by construction. It was left unfixed deliberately — making the
tie-break deterministic would alter the same interpolation the *solver* uses and force
every reference to be regenerated, a certain cost in shared machinery for a benefit
confined to pictures.

### `varsrc` uses `$(pwd)`, so a worktree gets a `PETSC_DIR` that does not exist

`varsrc` sets every path with `$(pwd)`:

    export HIGTREE_DIR=$(pwd)/higtree
    export PETSC_DIR=$(pwd)/bibliotecas/petsc-3.25.4/x86_64

Sourcing it from the main tree is right. Sourcing it from a **git worktree** is not:
`HIGTREE_DIR` and `HIGFLOW_DIR` correctly follow the worktree, but `PETSC_DIR` points
at a `bibliotecas/` that is *gitignored and therefore absent there* — the worktree
holds only the tarballs, never the built PETSc. The build then loses
`PETSC_CC_INCLUDES` entirely and dies on `petsc.h: No such file or directory`, a
hundred lines down in the output.

What makes it cost a whole run rather than a minute: `ci/run_suite.py` reports the
failure as `build failed: <first 80 chars of the command>`, truncated before the
compiler ever gets to say what went wrong. Every case shows the same truncated line,
so it reads like a problem with the tree, not with one environment variable.

Run a suite in a worktree with PETSc taken from the main tree — it is an external,
read-only dependency and sharing it is correct:

    HIGTREE_DIR=$W/higtree HIGFLOW_DIR=$W/higflow \
    PETSC_DIR=$MAIN/bibliotecas/petsc-3.25.4/x86_64 PETSC_ARCH= \
    python3 ci/run_suite.py

Note also that `build_for_dim` uses `env.setdefault`, so a `HIGTREE_DIR` already
exported in the shell **wins over the worktree**. Sourcing `varsrc` in the main tree
and then running the suite from a worktree silently builds and links the main tree.

### The physics library compiled with no warnings at all until 2026-09-18

`higtree/Makefile` has carried `-Wall -Wno-unused-result -Wno-unused-variable` all
along. `higflow/Makefile` had **no `-W` flag of any kind** — the entire physics
library, solvers included, compiled with warnings off. Two of the bugs found on
2026-09-18 are ones `-Wall` reports at every optimization level, `-O0` included.

`higflow/Makefile` now mirrors `higtree`. The baseline it exposes is large and mostly
benign, so read it by category rather than by count:

|DIM|total|what dominates|
|---|---|---|
|2|544|176 `array-bounds` in the 3D branch that `switch (DIM)` makes unreachable|
|3|389|154 `unused-but-set-variable`, 133 `switch`|

The categories worth reading are the small ones: `maybe-uninitialized` (21 at DIM=2,
24 at DIM=3 — this is the one that found both bugs), `return-type`, `dangling-else`,
`parentheses`. Note `maybe-uninitialized` needs `-O1` or higher: a `-fsyntax-only`
pass reports none of it.

The example Makefiles are unaffected. They carry their own `CFLAGS` with `-Werror`,
but their `%.o: %.c` rule only reaches sources in the example's own directory; the
shared library arrives as `libhigflow$(DIM)d.a`, built by `higflow/Makefile`.


### `git push <branch>` answers about the branch you named, not the work you did

`git push origin migracao-cpp` run from a checkout of a *different* branch replies
`Everything up-to-date`. That is true about `migracao-cpp` and says nothing about
the commit you just made, which is still only on disk. Read as confirmation, it
reports work as sent that was never sent.

**Check `git branch -vv` before pushing, not the branch name you remember.** It
shows where you are, each branch's upstream, and the divergence, so
`Everything up-to-date` can no longer be read as "my work went out".

Two more in the same family, both from 2026-09-19:

- **`git diff master...branch` (THREE dots) compares against the merge base, not
  against the tip of `master`.** On a branch whose base is from 2024, it presents
  as that branch's novelty everything `master` also gained since then. It reported
  `+638` lines and a set of functions as belonging to a branch that in fact has
  `domain.c` byte-identical to `master`, and inflated a file-overlap count from 26
  to 79. For "what does this branch change", use TWO dots.
- **`make` reports a missing *library* as `No rule to make target 'X-3d'`**, which
  reads as a missing *source*. In `higtree/tests` it meant one of two things:
  `libhig2d.a` did not exist because the tree was built for `DIM=3`, or
  `INSPATH`/`HIGPATH` were empty because `HIGTREE_DIR` was not exported — `varsrc`
  does not set it, only the test driver does. That Makefile now checks the library
  and fails with its own message.


### The shared pattern

None of these failed loudly. A step fails or lies, later steps run on stale
state, and the number that comes out looks reasonable. When a result is suspiciously
clean — an empty log read as "zero occurrences", a `rc=0` from a run that should have
taken ten minutes — verify the step produced what you assumed before trusting it.
