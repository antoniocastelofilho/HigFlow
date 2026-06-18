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
