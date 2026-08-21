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
| `../.dockerignore` | Build-context exclusions — lives at the repository root |

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
changes. PETSc takes 30 to 50 minutes to compile and changes approximately
never; the project sources change every commit. Putting them in that order means
editing a `.c` file rebuilds the last stage only.

This is also why PETSc is downloaded inside the build rather than copied from
the context: `COPY` of anything from the repository would tie the PETSc layer to
the repository's contents, and every source edit would trigger an hour-long
rebuild.

libfyaml goes the other way — it *is* copied from the context. The committed
`libfyaml-master.zip` is a snapshot of upstream `master`, not a tagged release,
so cloning `master` at build time would give a different revision than the one
the project is known to build against. 451 KB is a fair price for that
certainty, and libfyaml compiles in under a minute so the cache cost is
negligible.

## One recipe, not two

`apptainer.def` bootstraps from the Docker image rather than repeating its build
steps.

This repository is a live demonstration of why. `CMakeLists.txt` and the
Makefiles each list a different set of source files — `hig-flow-timestep.c` is
in one, `hig-flow-vof-elvira.c` in the other — so neither builds the whole
project, and two people using two build methods end up with binaries that have
different features. A second, hand-maintained container recipe would drift the
same way, and the drift would be invisible until someone's cluster run behaved
differently from their laptop run.

## Known workarounds carried forward

Two things in the image reproduce existing behaviour rather than fixing it,
because fixing them belongs to the build system rather than to a container.

**`libHYPRE_krylov.so`.** The Makefiles link `-lHYPRE_krylov`. HYPRE was split
into `libHYPRE_krylov`, `libHYPRE_IJ_mv` and others until roughly version 2.11,
and has shipped as a single `libHYPRE` ever since. The image creates the same
symlink the repository's install script creates. The real fix is to change the
flag to `-lHYPRE` in `higflow/Makefile` and the example Makefiles.

**Building both dimensions of `libhigflow`.** `higflow/Makefile:83` names its
objects `hig-flow-%.o`, with no dimension in the name, while
`higtree/Makefile:90` correctly uses `%-$(DIM)d.o`. So `make DIM=2 && make DIM=3`
in `higflow/` finds the `DIM=2` objects up to date and links them into
`libhigflow3d.a` — an archive carrying the wrong dimension with nothing to
indicate it. `make clean` is no way out either, since it deletes
`$(HIGFLOW_LIBPATH)/*.a`, taking the 2D library with it.

The Dockerfile clears the objects by hand between dimensions and keeps the 2D
archive aside. The real fix is to put the dimension in the object name, as
HigTree already does.

## Verifying a build

```bash
docker run --rm higflow:latest info
docker run --rm higflow:latest examples
mkdir -p cases && docker run --rm -v "$PWD/cases:/work" higflow:latest case example2d_Newt
```

The third produces VTK files under `cases/example2d_Newt/VTKS/` on the host.
