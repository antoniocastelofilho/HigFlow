# Running HigFlow in a container

This is the shortest path to a working HigFlow. It needs one tool installed and
nothing else - no PETSc build, no MPI configuration, no environment variables.

```bash
git clone https://github.com/antoniocastelofilho/HigFlow.git
cd HigFlow
docker build -f containers/Dockerfile -t higflow:latest .
mkdir cases
docker run --rm -v "$PWD/cases:/work" higflow:latest case example2d_Newt
```

The first build takes a few minutes, almost all of it compiling PETSc - six and
a half minutes on sixteen cores, proportionally longer on fewer. It happens
once. After that, starting a simulation takes seconds.

---

## Contents

- [Why a container, for this project specifically](#why-a-container-for-this-project-specifically)
- [The three concepts you need](#the-three-concepts-you-need)
- [Docker or Apptainer](#docker-or-apptainer)
- [Installing Docker](#installing-docker)
- [Building the image](#building-the-image)
- [Running a case](#running-a-case)
- [Running in parallel](#running-in-parallel)
- [Developing inside the container](#developing-inside-the-container)
- [Using compose](#using-compose)
- [On a cluster with Apptainer](#on-a-cluster-with-apptainer)
- [Troubleshooting](#troubleshooting)
- [What the image contains](#what-the-image-contains)

---

## Why a container, for this project specifically

The generic argument for containers is reproducibility. The specific argument
here is that HigFlow's dependency stack is the hardest part of using it.

A native build needs PETSc compiled from source, one MPI implementation
consistently used by every component, HDF5 built against that same MPI, Zoltan
from Trilinos, glib, and libfyaml - which has no distribution package and must
be built from a source snapshot. These have to agree with each other. A binary
compiled against one MPI and launched by another does not report "wrong MPI"; it
hangs, or produces wrong numbers on more than one rank, and the failure looks
like a bug in the solver.

The repository's own install script shows how easy this is to get wrong: it
installs OpenMPI and MPICH side by side, then asks PETSc to download a third
OpenMPI on top of both.

A container settles all of it once. The image is built with exactly one MPI, and
every component in it was compiled against that one.

The second reason is that it makes HigFlow usable on Windows, where it otherwise
is not. A native Windows build is not viable - OpenMPI and libfyaml have no
supported Windows port, and `CMakeLists.txt` requires `libnuma`, which exists
only on Linux - so a container, or WSL2, is the whole story there. See
[the Windows guide](windows.md).

## The three concepts you need

You can use containers productively knowing only three things.

**An image** is a filesystem plus a command to run, built once and then frozen.
`higflow:latest` is an image: Ubuntu, PETSc, the compiled HigFlow libraries and
the example cases. Nothing in an image ever changes.

**A container** is one running instance of an image. Start it, it does something,
it exits, and everything it wrote is discarded. That last part surprises people:
by default a simulation's output dies with the container.

**A volume, or bind mount,** is a directory on your machine made visible inside
the container. This is how results survive. The flag is `-v host_path:container_path`:

```bash
docker run --rm -v "$PWD/cases:/work" higflow:latest case example2d_Newt
```

`$PWD/cases` on your machine appears as `/work` inside. The simulation writes to
`/work`, which is your `cases/` directory, so when the container exits the VTK
files are sitting on your disk.

The `--rm` deletes the container when it exits. Use it always. Without it,
stopped containers accumulate silently until `docker system df` surprises you.

The mental model that makes this click: **the image is read-only and disposable;
your data lives on the host and is mounted in.** Nothing important should ever
be inside a container.

## Docker or Apptainer

Both are supported here. Which one you use is usually decided for you.

| | Docker | Apptainer |
|---|---|---|
| Your own machine | yes | yes |
| Shared cluster | almost never permitted | the standard choice |
| Runs as | a root daemon | your own user |
| Image format | layers in a daemon-managed store | a single `.sif` file |
| Mounts your directory by default | no, `-v` required | yes |

Clusters refuse Docker for a concrete reason: anyone who can talk to the Docker
daemon can start a container that mounts the host filesystem as root, which
makes daemon access equivalent to root on the node. Apptainer runs unprivileged
and keeps your identity inside the container, so a shared machine can allow it.

Use Docker on your laptop and workstation. Use Apptainer on the cluster. The
Apptainer image is built from the Docker one, so they are the same software.

## Installing Docker

### Windows

Two options.

**Docker Desktop** is the usual choice. Install it from
[docker.com](https://www.docker.com/products/docker-desktop/), make sure the
WSL2 backend is enabled in its settings, and `docker` then works from
PowerShell. Note that Docker Desktop requires a paid subscription for large
companies; it is free for personal use, education and small businesses.

**Docker Engine inside WSL2** avoids that entirely and is lighter. Install a
WSL2 distribution first:

```powershell
wsl --install -d Ubuntu-22.04
```

then make sure systemd is enabled in it, which Docker's service needs:

```bash
# inside WSL, as root
cat /etc/wsl.conf     # should contain [boot] with systemd=true
```

If it does not, add it and run `wsl --shutdown` from PowerShell, then reopen the
distribution. Then install Docker's official packages:

```bash
sudo install -m 0755 -d /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg \
  | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg
sudo chmod a+r /etc/apt/keyrings/docker.gpg
echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] \
https://download.docker.com/linux/ubuntu $(. /etc/os-release && echo "$VERSION_CODENAME") stable" \
  | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io \
                        docker-buildx-plugin docker-compose-plugin
sudo systemctl enable --now docker
docker --version
```

`docker` then works inside WSL, and from PowerShell as `wsl docker ...`.

### Linux

Use your distribution's Docker packages, or Docker's own repository as above.
Then add yourself to the `docker` group so you do not need `sudo` for every
command:

```bash
sudo usermod -aG docker "$USER"
newgrp docker
```

Be aware of what that grants: membership in the `docker` group is equivalent to
root on that machine, because you can start a container that mounts `/`. That is
fine on your own workstation and not fine on a shared one.

### macOS

Docker Desktop, or [Colima](https://github.com/abiosoft/colima) if you prefer an
open alternative. On Apple Silicon the image builds natively for `arm64`; the
build takes the same time and the result is a native `arm64` image.

## Building the image

From the repository root, not from `containers/`:

```bash
docker build -f containers/Dockerfile -t higflow:latest .
```

The trailing `.` is the build context - the directory Docker sends to the
daemon. The root `.dockerignore` keeps out what the image does not need: the
PETSc archive, compiled binaries, simulation output and the LaTeX manuals.

### What takes the time

The build has four stages, ordered by how often each one changes.

Measured on sixteen cores; the PETSc stage is the one that scales with core
count, so expect it to dominate more on a smaller machine.

| Stage | What it does | Measured |
|---|---|---|
| `base` | Ubuntu 22.04 plus system packages | 27 s |
| `petsc` | downloads and compiles PETSc 3.14.0 | 5 min 3 s |
| `fyaml` | compiles libfyaml | 12 s |
| `higflow` | compiles HigTree and HigFlow, 2D and 3D | 17 s |
| | sending the build context | 14 s |
| | **total** | **6 min 30 s** |

The ordering matters. Docker caches each stage, and editing a `.c` file only
invalidates the last one - so a rebuild after a source change takes minutes, not
an hour. This is why the Dockerfile downloads PETSc rather than copying it from
the repository: a change to any file in the context would otherwise invalidate
the PETSc layer.

The PETSc archive is verified against a SHA-256 recorded in the Dockerfile. If
the download is corrupted or the mirror serves something else, the build fails
immediately rather than producing a subtly different PETSc.

### Rebuilding

```bash
docker build -f containers/Dockerfile -t higflow:latest .          # uses cache
docker build --no-cache -f containers/Dockerfile -t higflow:latest .   # from scratch
```

Use `--no-cache` when you have changed something the cache cannot see, such as
an apt package that has been updated upstream.

## Running a case

```bash
mkdir -p cases
docker run --rm -v "$PWD/cases:/work" higflow:latest case example2d_Newt
```

This copies `example2d_Newt` out of the image into `cases/`, builds it, and runs
it. Everything - the binary, the VTK output, the restart files - ends up in
`cases/example2d_Newt/` on your machine.

List what is available:

```bash
docker run --rm higflow:latest examples
```

Inspect the image's toolchain:

```bash
docker run --rm higflow:latest info
```

Get a shell and work by hand:

```bash
docker run --rm -it -v "$PWD/cases:/work" higflow:latest shell
```

Anything the entrypoint does not recognise is run verbatim, so the image is also
usable as a plain toolbox:

```bash
docker run --rm -v "$PWD/cases:/work" higflow:latest \
    bash -c 'cd example2d_Newt && make && make run NP=2'
```

### Viewing results

The VTK files are on your host, so open them with whatever you already use.
ParaView and VisIt both read them directly. On Windows, files written inside WSL
are reachable from Explorer at `\\wsl$\Ubuntu-22.04\home\...`, so a ParaView
installed on Windows can open output produced in the container without copying
anything.

## Running in parallel

```bash
docker run --rm --shm-size=1g -v "$PWD/cases:/work" \
    higflow:latest case example2d_Newt 4
```

`--shm-size=1g` matters. OpenMPI's shared-memory transport allocates from
`/dev/shm`, and Docker gives a container 64 MB of it by default. A parallel run
that exhausts that fails with a message about `vader` or `shm` that gives no
hint the cause is a container flag.

The container has as many cores as Docker gives it. On Linux that is all of
them; on Docker Desktop it is whatever the settings allow, which is worth
checking before concluding that HigFlow scales badly.

Two limits worth knowing. There is no point running more ranks than you have
physical cores - MPI ranks spin while waiting, so oversubscribing makes a
simulation slower, not faster. And a single container is a single machine: to
run across several nodes you need Apptainer and the cluster's scheduler, covered
below.

## Developing inside the container

The development image adds a debugger and a memory checker, and mounts the
repository live rather than copying it in, so an edit on your host is
immediately visible inside:

```bash
docker build -f containers/Dockerfile.dev -t higflow:dev .
docker run --rm -it \
    --cap-add=SYS_PTRACE --security-opt seccomp=unconfined \
    --shm-size=1g \
    -v "$PWD:/src" -v "$PWD/cases:/work" \
    higflow:dev
```

`--cap-add=SYS_PTRACE` and `--security-opt seccomp=unconfined` are what gdb and
valgrind need in order to attach to a process. Docker blocks that by default.

Inside:

```bash
cd /src/higflow/example2d_Newt
make DEBUG=1
gdb --args ./ns-example input/example-2d.load output/x.save VTKS/x.print
```

To generate a `compile_commands.json` for clangd from the Makefile build:

```bash
cd /src/higtree && bear -- make DIM=2
```

## Using compose

Compose gives the same operations shorter names.

```bash
docker compose -f containers/docker-compose.yml build
docker compose -f containers/docker-compose.yml run --rm info
docker compose -f containers/docker-compose.yml run --rm examples
docker compose -f containers/docker-compose.yml run --rm shell
docker compose -f containers/docker-compose.yml run --rm case example2d_Newt 4
docker compose -f containers/docker-compose.yml run --rm dev
```

`shm_size: 1gb` and the volume mount are set for every service, so parallel runs
work without remembering the flags.

## On a cluster with Apptainer

Build the `.sif` on a machine where you have Docker, then copy that one file to
the cluster:

```bash
docker build -f containers/Dockerfile -t higflow:latest .
apptainer build higflow.sif containers/apptainer.def
scp higflow.sif user@cluster:~/
```

Apptainer mounts your current directory and `$HOME` automatically, so no volume
flag is needed:

```bash
cd my-case
apptainer run higflow.sif case example2d_Newt 4
```

Under Slurm, let the scheduler place the ranks and hand each one the container,
rather than launching `mpirun` inside a single container:

```bash
#!/bin/bash
#SBATCH --job-name=higflow
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --time=04:00:00

module load apptainer          # or singularity, depending on the site

srun apptainer exec higflow.sif \
    /opt/higflow/higflow/example2d_Newt/ns-example \
    input/example-2d.load output/run.save VTKS/run.print \
    -ksp_type bcgs -pc_type bjacobi -ksp_atol 1e-10 -ksp_rtol 1e-10
```

**The one thing that will catch you on a multi-node run is MPI compatibility.**
The image carries OpenMPI from Ubuntu 22.04. Within a single node that is all
that matters. Across nodes, the container's MPI has to work with the host's
interconnect drivers - InfiniBand or Omni-Path - and if it cannot, it silently
falls back to TCP or hangs. Ask the site's support what they expect; some bind
the host MPI into the container, others provide a matching module. Every cluster
answers this differently, and the symptom is a job that works on one node and
hangs on two.

## Troubleshooting

**`permission denied` writing to the mounted directory.** The container runs as
a user with UID 1000. If your host UID differs, files it writes may not be
yours. Check with `id -u`; if you are not 1000, run as yourself:

```bash
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD/cases:/work" higflow:latest ...
```

**Files owned by root after a run.** You ran the container as root. Fix the
ownership with `sudo chown -R "$USER:$USER" cases/` and use `--user` as above
next time.

**A parallel run dies mentioning `vader`, `shm` or `/dev/shm`.** Shared memory
is exhausted. Add `--shm-size=1g`.

**`mpirun` refuses to run as root.** Do not run the container as root. The image
already provides a normal user; this happens when `--user 0:0` is passed or a
derived image switches back to root.

**The build is killed partway through PETSc.** Out of memory. Docker Desktop
defaults to a fraction of your RAM - raise it to at least 4 GB in Settings →
Resources. On WSL2, create `%UserProfile%\.wslconfig`:

```ini
[wsl2]
memory=8GB
```

then `wsl --shutdown` and reopen.

**`no space left on device`.** Images and build cache accumulate. See what is
using space and reclaim it:

```bash
docker system df
docker system prune -a      # deletes unused images - you will rebuild
```

Building HigFlow needs roughly 10 GB free while it runs, and the finished image
is a few GB.

**The build fails at `COPY bibliotecas/libfyaml-master.zip`.** You are building
from `containers/` instead of the repository root. The context must be the root:
`docker build -f containers/Dockerfile -t higflow:latest .`

**`exec /usr/local/bin/entrypoint.sh: no such file or directory`, but the file
is clearly there.** The script has Windows line endings, so the kernel is
looking for an interpreter named `/usr/bin/env bash\r`. The Dockerfile strips
carriage returns during the build, so this means something reintroduced them -
check `git config core.autocrlf` and that `containers/.gitattributes` survived
your checkout.

**The simulation runs but `cases/` stays empty.** The `-v` flag is missing, so
everything was written inside the container and discarded when it exited.

## What the image contains

| | |
|---|---|
| Base | Ubuntu 22.04 |
| MPI | OpenMPI, from Ubuntu - the only MPI present, deliberately |
| PETSc | 3.14.0, `--with-debugging=0`, shared libraries, built with the `mpicc`/`mpicxx`/`mpif90` wrappers; HYPRE and fblaslapack downloaded by PETSc |
| HDF5 | Ubuntu's OpenMPI build |
| Zoltan | Trilinos Zoltan, from Ubuntu |
| libfyaml | built from the snapshot committed to the repository |
| Boost | headers only - the code uses `geometry` and `ublas`, both header-only |
| HigFlow | HigTree and HigFlow compiled for both 2D and 3D, plus every example case |

Paths inside:

```
/opt/higflow      the repository, built
/opt/petsc        PETSc install prefix
/opt/libfyaml     libfyaml install prefix
/work             working directory - mount your host directory here
```

Environment already set: `PETSC_DIR`, `PETSC_ARCH`, `PETSC_EXTRA_LIB`,
`HIGTREE_DIR`, `HIGFLOW_DIR`, `PKG_CONFIG_PATH`, `LD_LIBRARY_PATH`, `CPATH` and
`LIBRARY_PATH`. There is no `source varsrc` step - the image does not need one.

`CPATH` and `LIBRARY_PATH` carry the OpenMPI paths, because the example
Makefiles compile with `gcc` rather than `mpicc` and would otherwise not find
`mpi.h`. See [`containers/README.md`](../../containers/README.md) for why.
