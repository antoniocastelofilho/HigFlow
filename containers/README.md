# containers/

Container definitions for HigFlow. The user-facing guide is
[`docs/install/containers.md`](../docs/install/containers.md); this file
describes the files themselves and the decisions behind them.

| File | Purpose |
|---|---|
| `Dockerfile` | The image. Four stages: system packages, PETSc, libfyaml, HigFlow |
| `Dockerfile.dev` | Adds gdb, valgrind, clangd; mounts sources instead of copying |
| `docker-compose.yml` | Named services for build, shell, info, examples, case, dev |
| `apptainer.def` | HPC image, bootstrapped from the Docker image |
| `entrypoint.sh` | Command dispatch inside the container |
| `.gitattributes` | Forces LF on everything copied into an image |
| `../.dockerignore` | Build-context exclusions - lives at the repository root |

Build from the repository root, never from this directory:

```bash
docker build -f containers/Dockerfile -t higflow:latest .
```

---

## Relationship to the existing container work

This replaces `conteiner/Dockerfile` and `conteiner/Dockerfile.petsc`, written by
**Pedro Coimbra** on the `PC_ImproveDocumentation` branch, and
`stacks/singularity/higflow_image.def`. Those were the starting point; the
structure and the corrections below are what changed.

### Corrections

| | Was | Now |
|---|---|---|
| PETSc configure | `--with-debubbing=yes` | `--with-debugging=0`. PETSc rejects unrecognised options, so the typo aborts the build |
| `.bashrc` | `echo '. $HOME/.varsrc"' >> $HOME/.bashrc` | Environment set with `ENV`. The unbalanced quote wrote a broken line into `.bashrc`, so every later shell in the container failed to parse it |
| MPI | OpenMPI *and* MPICH from apt, then `--download-openmpi` in PETSc | One MPI. PETSc is pointed at the system OpenMPI with `--with-mpi-dir=/usr` |
| Stages | `FROM ubuntu_petsc3.14:v01`, a tag you had to build and name by hand from a second Dockerfile, in an order nothing documented | One `docker build`, four stages |
| Build context | no `.dockerignore` | `.dockerignore` at the repository root |
| Permissions | `chmod 777` | A normal user owning what it needs |
| Header comment | "Ubuntu22.04x64 + OpenFOAM-9 + Python 3.10" | Describes this image |
| Boost | `libboost-all-dev` | `libboost-dev`. The code uses `boost/geometry` and `boost/numeric/ublas`, both header-only |
| `$HOME` | `"/home/hig_user/"`, trailing slash, producing `//` in every derived path | No trailing slash |
| PETSc source | 37 MB archive copied from the context | Downloaded during the build and checked against its SHA-256 |

The Singularity definition could not build as written either: `%post` invokes
`./scripts/pip.sh` by a relative path that does not resolve in that context, and
the PETSc script copies `configure.log` into `/pacotes`, a directory nothing
creates.

---

## Why the stages are ordered this way

Docker caches each stage and invalidates everything after the first one that
changes. PETSc takes five minutes to compile on sixteen cores and considerably
longer on fewer, and changes approximately never; the project sources change
every commit. Putting them in that order means editing a `.c` file rebuilds the
last stage only, which is seconds rather than minutes.

This is also why PETSc is downloaded inside the build rather than copied from
the context: `COPY` of anything from the repository would tie the PETSc layer to
the repository's contents, and every source edit would trigger a full rebuild of
it.

libfyaml goes the other way - it *is* copied from the context. The committed
`libfyaml-master.zip` is a snapshot of upstream `master`, not a tagged release,
so cloning `master` at build time would give a different revision than the one
the project is known to build against. 451 KB is a fair price for that
certainty, and libfyaml compiles in under a minute so the cache cost is
negligible.

## One recipe, not two

`apptainer.def` bootstraps from the Docker image rather than repeating its build
steps.

This repository is a live demonstration of why. `CMakeLists.txt` and the
Makefiles each list a different set of source files - `hig-flow-timestep.c` is
in one, `hig-flow-vof-elvira.c` in the other - so neither builds the whole
project, and two people using two build methods end up with binaries that have
different features. A second, hand-maintained container recipe would drift the
same way, and the drift would be invisible until someone's cluster run behaved
differently from their laptop run.

## Known defects worked around, not fixed

Three things in the image reproduce existing behaviour rather than correcting
it, because the correction belongs to the build system rather than to a
container. All three were found while getting this image to run a case.

**`-lHYPRE_krylov`.** The Makefiles link `-lHYPRE_krylov`. HYPRE was split into
`libHYPRE_krylov`, `libHYPRE_IJ_mv` and others until roughly version 2.11, and
has shipped as a single `libHYPRE` ever since. The image creates the same
symlink the repository's install script creates. The fix is to link `-lHYPRE`.

**`higflow/Makefile:83` omits the dimension from object names.** It uses
`hig-flow-%.o`, while `higtree/Makefile:90` correctly uses `%-$(DIM)d.o`. So
`make DIM=2 && make DIM=3` in `higflow/` finds the `DIM=2` objects up to date and
links them into `libhigflow3d.a` - an archive carrying the wrong dimension with
nothing to indicate it. `make clean` is no way out either: it deletes
`$(HIGFLOW_LIBPATH)/*.a`, taking the other library with it.

The Dockerfile clears the objects by hand between dimensions and keeps one
archive aside. The fix is to put the dimension in the object name, as HigTree
already does.

**The examples cannot be moved, and share one set of objects.** Every one of the
eleven examples includes `../src/hig-flow-*.h` and links `../src/hig-flow-*.o`
directly rather than `libhigflow<dim>d.a`. Two consequences follow.

A case copied anywhere outside `higflow/` stops compiling, because `../src` no
longer resolves. And since `higflow/src/` holds a single set of objects carrying
whichever dimension was built last, the nine 2D examples and the two 3D examples
cannot both be buildable at the same time.

The image handles this by building 3D first and 2D last, so the objects left in
place suit the nine 2D examples, and recording the dimension in
`higflow/.built-dim`. The entrypoint reads that marker and rebuilds the objects -
a few seconds - when a 3D case is requested. Cases are built and run where they
live, with only `VTKS`, `DATA` and `output` symlinked into the mounted `/work`,
so output reaches the host without the case having to move.

The fix is for examples to link the library rather than raw objects, and to
include headers through an include path rather than a relative path. That is
what the case-runner work is for.

**The examples compile with `gcc`, not `mpicc`.** All eleven set `CC = gcc`,
while `higflow/Makefile` and `higtree/Makefile` both set `CC = mpicc`. An example
therefore has exactly one source for the MPI include path: `$(PETSC_CC_INCLUDES)`,
inherited from PETSc's `petscvariables`.

This turns out to explain a choice in the original install script that looks
simply wrong. `--download-openmpi` makes PETSc build its own MPI into the PETSc
prefix, so `mpi.h` lands in `/opt/petsc/include` and `PETSC_CC_INCLUDES` picks it
up. The examples then compile - by accident of layout rather than by design.
Point PETSc at a system OpenMPI instead, and the accident stops working:

```
pdomain.h:5:10: fatal error: mpi.h: No such file or directory
```

on Debian and Ubuntu, where the OpenMPI headers live under
`/usr/lib/x86_64-linux-gnu/openmpi/include` and only the `mpicc` wrapper knows to
add that path.

Pointing `--with-mpi-dir` at that prefix instead is not a fix either: Debian
splits OpenMPI's Fortran bindings and modules out of it, so PETSc's configure
fails with `Fortran error! mpi_init() could not be located!`. There is no single
directory that is the MPI root on these distributions.

The image therefore does two things. PETSc is configured with the wrappers -
`--with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90` - which know every path.
And `CPATH` and `LIBRARY_PATH` are set in the image so that plain `gcc` and `ld`
see the same directories the wrapper would have added, without any Makefile being
modified. The underlying fix is for the examples to compile with `mpicc`, as the
libraries already do.

## Verifying a build

```bash
docker run --rm higflow:latest info
docker run --rm higflow:latest examples
mkdir -p cases && docker run --rm -v "$PWD/cases:/work" higflow:latest case example2d_Newt
```

The third produces VTK files under `cases/example2d_Newt/VTKS/` on the host.
