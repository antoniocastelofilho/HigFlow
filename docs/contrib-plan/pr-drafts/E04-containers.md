# PR draft - E04

**Branch:** `juniormar/04-containers` para `antoniocastelofilho/HigFlow:master`
**Title:** `Add reproducible container images for desktop and HPC`

---

## Summary

A single `docker build` produces a working HigFlow - PETSc, OpenMPI, HDF5,
Zoltan, libfyaml and the compiled libraries - with nothing installed on the host.
The same image converts to an Apptainer `.sif` for clusters. Adds a container
guide written for someone who has not used containers before.

Builds on `conteiner/Dockerfile` and `conteiner/Dockerfile.petsc` by Pedro
Coimbra on `PC_ImproveDocumentation`, and on
`stacks/singularity/higflow_image.def`.

## Motivation

The dependency stack is the hardest part of using HigFlow, and the hardest part
of that is making the pieces agree. PETSc, HDF5 and the solver all have to be
built against the same MPI. When they are not, the failure does not announce
itself: a run hangs, or produces different numbers on more than one rank, and it
looks like a solver bug.

`install_higflow_ubuntu22` installs `openmpi-bin`, `libopenmpi-dev` **and**
`mpich`, then asks PETSc to `--download-openmpi` on top of both - three MPI
stacks whose `mpicc`, `mpirun` and `libmpi.so` compete on `PATH` and in the
linker.

A container settles this once, and is also what makes HigFlow usable on Windows,
where a native build is not viable.

## Changes

```
containers/
├── Dockerfile              four stages: base, petsc, fyaml, higflow
├── Dockerfile.dev          + gdb, valgrind, clangd; sources mounted, not copied
├── docker-compose.yml      build, shell, info, examples, case, dev
├── apptainer.def           HPC image, bootstrapped from the Docker image
├── entrypoint.sh           command dispatch
├── .gitattributes          LF on everything copied into an image
└── README.md               file map and the decisions behind them
.dockerignore               build-context exclusions, at the repository root
docs/install/containers.md  the guide
```

### Corrections to the existing container material

| | Was | Now |
|---|---|---|
| PETSc configure | `--with-debubbing=yes` | `--with-debugging=0`. PETSc rejects unrecognised options, so the typo aborts the build |
| `.bashrc` | `echo '. $HOME/.varsrc"' >> $HOME/.bashrc` | Environment set with `ENV`. The unbalanced quote wrote a broken line into `.bashrc`, breaking every later shell in the container |
| MPI | OpenMPI and MPICH from apt, plus `--download-openmpi` | One MPI; PETSc pointed at the system OpenMPI |
| Stages | `FROM ubuntu_petsc3.14:v01` - a tag to build and name by hand from a second Dockerfile, in an undocumented order | One `docker build` |
| Context | no `.dockerignore` | present, at the repository root |
| Permissions | `chmod 777` | a normal user owning what it needs |
| Header | "Ubuntu22.04x64 + OpenFOAM-9 + Python 3.10" | describes this image |
| Boost | `libboost-all-dev` | `libboost-dev` - only `boost/geometry` and `boost/numeric/ublas` are used, both header-only |
| `$HOME` | `"/home/hig_user/"` - trailing slash producing `//` in derived paths | no trailing slash |
| PETSc source | 37 MB archive copied from the context | downloaded during the build, verified against its SHA-256 |

`stacks/singularity/higflow_image.def` could not build as written either: `%post`
invokes `./scripts/pip.sh` by a relative path that does not resolve there, and
the PETSc script copies `configure.log` into `/pacotes`, a directory nothing
creates.

### Four defects found while getting this to run, worked around rather than fixed

All four belong to the build system, so the image reproduces existing behaviour
rather than silently diverging from what a host build produces. Each was found
because the image would not build or would not run a case.

**`higflow/Makefile:83` omits the dimension from object file names.** It uses
`hig-flow-%.o`, while `higtree/Makefile:90` correctly uses `%-$(DIM)d.o`. So
`make DIM=2 && make DIM=3` in `higflow/` finds the `DIM=2` objects up to date and
links them into `libhigflow3d.a` - an archive carrying the wrong dimension with
nothing to indicate it. `make clean` is no escape: it deletes
`$(HIGFLOW_LIBPATH)/*.a`, taking the other library with it. The image clears the
objects by hand between dimensions and keeps one archive aside.

**The examples cannot be relocated.** All eleven include `../src/hig-flow-*.h`
and link `../src/hig-flow-*.o` directly rather than `libhigflow<dim>d.a`. A case
copied anywhere outside `higflow/` stops compiling. Combined with the previous
defect, it also means the nine 2D and two 3D examples cannot both be buildable at
once, since `higflow/src` holds one set of objects. The image builds 3D first and
2D last, records the dimension in `.built-dim`, runs cases where they live, and
symlinks only `VTKS`, `DATA` and `output` into the mounted directory.

**The examples compile with `gcc`, not `mpicc`.** All eleven set `CC = gcc`,
while both library Makefiles set `CC = mpicc`. An example's only source of the
MPI include path is therefore `$(PETSC_CC_INCLUDES)`.

This explains something in `install_higflow_ubuntu22` that had looked purely
wrong. `--download-openmpi` puts `mpi.h` inside the PETSc prefix, so
`PETSC_CC_INCLUDES` happens to cover it and the examples compile - by accident of
layout. Point PETSc at a system MPI and the accident stops working:
`pdomain.h:5:10: fatal error: mpi.h: No such file or directory`. On Debian and
Ubuntu there is no single MPI root either: `--with-mpi-dir=/usr` records an
include path with no `mpi.h`, and the OpenMPI prefix makes configure fail with
`Fortran error! mpi_init() could not be located!` because the Fortran bindings
live elsewhere. The image configures PETSc with the wrappers and sets `CPATH` and
`LIBRARY_PATH` so plain `gcc` sees the same directories.

**`-lHYPRE_krylov`.** HYPRE was split into `libHYPRE_krylov`, `libHYPRE_IJ_mv`
and others until roughly version 2.11 and has shipped as a single `libHYPRE`
since. The Makefiles still link the old name; the install script papers over it
with a symlink, and the image does the same. The fix is to link `-lHYPRE`.

### One further fix, needed to land this change at all

`.gitignore` used `**build**` and `**install**`. Outside the three special forms
(`**/`, `/**`, `/**/`), consecutive asterisks collapse to a single `*`, so those
patterns match any path containing the substring anywhere:

```
$ git check-ignore --no-index -v docs/install/containers.md
.gitignore:74:**install**  docs/install/containers.md
```

`git add` of the guide silently did nothing. The same rules match
`higtree/src/build-fringe.cpp`, a file both build systems compile; it survives
only because it was committed before the rules existed, and deleting and
re-adding it would drop it.

Both are now anchored to the repository root and matched as directories, which
is what they were meant to express: the CMake build output, and the install
prefix `CMakeLists.txt` places at `<source>/install`.

## Verification

Built and run end to end on Ubuntu 22.04 under WSL2, Docker Engine 29.7.2,
sixteen cores.

**Build.** 6 min 30 s from a cold cache. Per stage:

| Step | Measured |
|---|---|
| system packages | 27 s |
| PETSc configure + build + install | 5 min 3 s |
| libfyaml | 12 s |
| HigTree, both dimensions | 11 s |
| HigFlow, both dimensions | 6 s |
| build context transfer | 14 s |

Context sent to the daemon: **11.67 MB**, from a 96 MB tree - the rest excluded
by `.dockerignore`. Image: 299 MB of content, 1.34 GB on disk.

**Libraries.** All four archives produced, and the two HigFlow ones differ,
which is what shows the dimension workaround did its job rather than producing
two copies of the same objects:

```
libhig2d.a       10393132
libhig3d.a       10498110
libhigflow2d.a    2662630   sha256 29f5ec14…
libhigflow3d.a    2669766   sha256 6ab832ef…
```

**Cases.** Each run writes to a host directory through the bind mount:

| Case | Ranks | Output |
|---|---|---|
| `example2d_Newt` | 1 | 101 VTK files, 166 MB |
| `example2d_Newt` | 2 | 202 VTK files, 166 MB - one per rank per frame |
| `example3d_lid_driven` | 1 | 880 VTK files, 509 MB, stopped once verified |

The 3D case exercised the object rebuild path, reporting
`higflow/src objects are compiled for DIM=2; this case needs DIM=3` before
building. Output is well-formed:

```
# vtk DataFile Version 3.0
higtree
ASCII
DATASET UNSTRUCTURED_GRID

POINTS 25600 float
```

25 600 points is the 160 × 40 channel at four corners per cell, which is the
mesh `amrs/domain/ch-d-0.amr` describes.

**Failure paths.** A read-only mount is rejected before anything runs, naming
the uid on each side and both commands that fix it. An unknown case name exits 2
with a usable message.

**Static checks.** `docker compose config` resolves all seven services with
`shm_size` applied; `bash -n` passes on the entrypoint; the Apptainer definition
is syntactically valid.

## Notes for reviewers

**1. The Singularity definition is superseded, not deleted.**
`stacks/singularity/` is left in place; `containers/apptainer.def` is the
replacement. Removing the old one is the maintainers' call.

**2. `apptainer.def` bootstraps from the Docker image rather than repeating its
build steps.** This is deliberate. `CMakeLists.txt` and the Makefiles already
list different source files, so neither builds the whole project and two people
using two build methods get binaries with different features. A second
hand-maintained container recipe would drift the same way.

**3. Publishing the image is not included.** A GitHub Actions workflow that
builds and pushes to `ghcr.io` would let users skip the 30-50 minute PETSc
compile entirely, and belongs with the CI change rather than here.

**4. The base image is pinned by tag, not digest.** A digest pins exactly and is
the stronger choice for a release. The Dockerfile records how to resolve and
substitute one; it is not done here because the right digest is a decision about
release policy.
