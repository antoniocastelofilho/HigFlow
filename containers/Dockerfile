# HigFlow — reproducible build image
#
# Ubuntu 22.04, OpenMPI, PETSc 3.14.0, libfyaml, Zoltan, HDF5.
#
# Build from the REPOSITORY ROOT, not from this directory:
#     docker build -f containers/Dockerfile -t higflow:latest .
#
# Based on conteiner/Dockerfile and conteiner/Dockerfile.petsc by Pedro Coimbra
# on the PC_ImproveDocumentation branch. Restructured into a real multi-stage
# build and corrected — see containers/README.md for what changed and why.
#
# ─────────────────────────────────────────────────────────────────────────────
# Stage layout, ordered by how often each layer changes. PETSc takes about five
# minutes to compile on sixteen cores, proportionally longer on fewer, and almost
# never changes; the project sources change every commit. Keeping them in
# separate stages means editing a .c file rebuilds seconds of work rather than
# the whole chain.
#
#   base     → system packages
#   petsc    → PETSc, built from a checksum-verified release archive
#   fyaml    → libfyaml
#   higflow  → the project
# ─────────────────────────────────────────────────────────────────────────────

# Pinned by digest, not only by tag. The ubuntu:22.04 tag moves with every
# security refresh, so two builds a month apart would otherwise start from
# different filesystems. This is the digest the tag resolved to when the image
# was last built and verified end to end.
#
# To move to a newer base, resolve the current digest and replace it here:
#     docker buildx imagetools inspect ubuntu:22.04
FROM ubuntu:22.04@sha256:2edbbc5dc405e9612ba3584ce95480277e3eb374407b5505fe26f17df77c7dbc AS base

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Etc/UTC \
    LANG=C.UTF-8

# One MPI implementation, deliberately.
#
# The repository's install script installs openmpi-bin, libopenmpi-dev AND
# mpich, and then asks PETSc to --download-openmpi on top, leaving three MPI
# stacks whose mpicc, mpirun and libmpi.so compete on PATH and in the linker.
# A binary compiled against one and launched by another fails in ways that look
# like application bugs. This image installs OpenMPI only, and PETSc is pointed
# at it with --with-mpi-dir.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        gfortran \
        make \
        cmake \
        pkg-config \
        git \
        ca-certificates \
        curl \
        unzip \
        python3 \
        autoconf \
        automake \
        libtool \
        libltdl-dev \
        openmpi-bin \
        libopenmpi-dev \
        libhdf5-openmpi-dev \
        libtrilinos-zoltan-dev \
        libglib2.0-dev \
        libboost-dev \
        libyaml-dev \
        libnuma-dev \
        zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


# ─────────────────────────────────────────────────────────────────────────────
# PETSc
# ─────────────────────────────────────────────────────────────────────────────
FROM base AS petsc

ARG PETSC_VERSION=3.14.0
# SHA-256 of the archive this project has always built against.
ARG PETSC_SHA256=bf37688f68889e82b7609f6b36d535899c892bbe7335e93721d42e90ea8b03c3

WORKDIR /tmp/build

# Downloaded here rather than copied from the build context: the archive is
# 37 MB, and pulling it in would put it in a layer for no benefit. The checksum
# makes a corrupted or substituted download a hard failure.
RUN curl -fL --retry 3 -o petsc.tar.gz \
      "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-${PETSC_VERSION}.tar.gz" \
    && echo "${PETSC_SHA256}  petsc.tar.gz" | sha256sum -c - \
    && tar xzf petsc.tar.gz \
    && rm petsc.tar.gz

# --with-debugging=0 is the corrected form of --with-debubbing=yes, which
# appears in install_higflow_ubuntu22 and in both existing Dockerfiles. PETSc's
# configure rejects unrecognised options, so that typo aborts the build.
#
# The MPI compiler wrappers are given directly rather than a --with-mpi-dir.
# Debian and Ubuntu split OpenMPI across directories — headers under
# /usr/lib/x86_64-linux-gnu/openmpi/include, Fortran bindings and modules
# elsewhere — so neither /usr nor the openmpi prefix is a valid single MPI root.
# Pointing at /usr makes PETSc record MPI_INCLUDE = -I/usr/include, which does
# not contain mpi.h; pointing at the openmpi prefix makes configure fail with
# "Fortran error! mpi_init() could not be located!". The wrappers know every
# path, which is what they exist for.
RUN cd "petsc-${PETSC_VERSION}" \
    && ./configure \
        --prefix=/opt/petsc \
        --with-debugging=0 \
        --with-shared-libraries=1 \
        --with-cc=mpicc \
        --with-cxx=mpicxx \
        --with-fc=mpif90 \
        --with-hdf5-dir=/usr/lib/x86_64-linux-gnu/hdf5/openmpi \
        --download-hypre \
        --download-fblaslapack \
        COPTFLAGS=-O3 CXXOPTFLAGS=-O3 FOPTFLAGS=-O3 \
    && make PETSC_DIR="$(pwd)" PETSC_ARCH=arch-linux-c-opt all \
    && make PETSC_DIR="$(pwd)" PETSC_ARCH=arch-linux-c-opt install \
    && cd / && rm -rf /tmp/build

# Compatibility shim for a stale link flag.
#
# The Makefiles link -lHYPRE_krylov. HYPRE was split into libHYPRE_krylov,
# libHYPRE_IJ_mv and others until roughly 2.11; it has shipped as a single
# libHYPRE ever since. The repository's install script papers over this with the
# same symlink, in the opposite direction.
#
# The real fix is to change the flag to -lHYPRE in higflow/Makefile and the
# example Makefiles. That belongs to the build-system change, not here, so the
# image reproduces the existing behaviour rather than silently diverging from
# what a host build produces.
RUN if [ -f /opt/petsc/lib/libHYPRE.so ] && [ ! -e /opt/petsc/lib/libHYPRE_krylov.so ]; then \
        ln -s libHYPRE.so /opt/petsc/lib/libHYPRE_krylov.so ; \
    fi


# ─────────────────────────────────────────────────────────────────────────────
# libfyaml
# ─────────────────────────────────────────────────────────────────────────────
FROM base AS fyaml

WORKDIR /tmp/build

# Taken from the build context rather than cloned. The committed archive is a
# snapshot of upstream master, not a tagged release, so cloning master today
# would give a different revision than the one the project is known to build
# against. 451 KB is a fair price for that certainty.
COPY bibliotecas/libfyaml-master.zip .

RUN unzip -q libfyaml-master.zip \
    && cd libfyaml-master \
    && ./bootstrap.sh \
    && ./configure --prefix=/opt/libfyaml \
    && make -j"$(nproc)" \
    && make install \
    && cd / && rm -rf /tmp/build


# ─────────────────────────────────────────────────────────────────────────────
# HigFlow
# ─────────────────────────────────────────────────────────────────────────────
FROM base AS higflow

COPY --from=petsc /opt/petsc /opt/petsc
COPY --from=fyaml /opt/libfyaml /opt/libfyaml

ENV PETSC_DIR=/opt/petsc \
    PETSC_ARCH= \
    PETSC_EXTRA_LIB=-lmpi \
    PKG_CONFIG_PATH=/opt/libfyaml/lib/pkgconfig \
    LD_LIBRARY_PATH=/opt/petsc/lib:/opt/libfyaml/lib

# CPATH and LIBRARY_PATH give plain gcc and ld the OpenMPI paths that the mpicc
# wrapper would have added.
#
# This is not gratuitous: all eleven example Makefiles set `CC = gcc`, while
# higflow/Makefile and higtree/Makefile set `CC = mpicc`. An example's only
# source of the MPI include path is therefore $(PETSC_CC_INCLUDES), inherited
# from PETSc.
#
# That explains a choice in the repository's install script that otherwise looks
# simply wrong. --download-openmpi puts mpi.h inside the PETSc prefix, so
# PETSC_CC_INCLUDES happens to cover it and the examples compile — by accident of
# layout rather than by design. Using a system MPI removes the accident, so the
# image supplies the paths explicitly instead. The underlying fix is for the
# examples to compile with mpicc, as the libraries already do.
ENV CPATH=/usr/lib/x86_64-linux-gnu/openmpi/include \
    LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/openmpi/lib

# A normal user, because OpenMPI refuses to run as root without
# --allow-run-as-root, and because a container process writing into a mounted
# host directory as root leaves files the host user cannot delete.
RUN useradd --create-home --shell /bin/bash higflow \
    && mkdir -p /work \
    && chown higflow:higflow /work

ENV HIGFLOW_ROOT=/opt/higflow \
    HIGTREE_DIR=/opt/higflow/higtree \
    HIGFLOW_DIR=/opt/higflow/higflow

COPY containers/entrypoint.sh /usr/local/bin/entrypoint.sh
# The carriage returns are stripped as well as the mode being set. containers/
# has a .gitattributes forcing LF, but that only helps a clone that has it;
# an archive download, a copy through a Windows share, or an editor saving with
# CRLF would still leave "#!/usr/bin/env bash\r" here, and the container would
# fail to start with "bad interpreter: /usr/bin/env bash^M".
RUN sed -i 's/\r$//' /usr/local/bin/entrypoint.sh \
    && chmod 0755 /usr/local/bin/entrypoint.sh

COPY --chown=higflow:higflow . ${HIGFLOW_ROOT}

USER higflow
WORKDIR ${HIGFLOW_ROOT}

# DIM is a compile-time macro, so 2D and 3D are genuinely separate archives.
#
# higtree names its objects with the dimension — `%-$(DIM)d.o`, Makefile:90 —
# so building one after the other is correct.
RUN cd higtree && make DIM=2 && make DIM=3

# higflow does not. Its objects are `hig-flow-%.o`, Makefile:83, with no
# dimension in the name, so a second `make DIM=3` finds the DIM=2 objects up to
# date and links them into libhigflow3d.a. The archive would carry the wrong
# dimension with nothing to indicate it. `make clean` is no escape either: it
# deletes $(HIGFLOW_LIBPATH)/*.a, taking the 2D library with it.
#
# So the objects are cleared by hand between dimensions and one archive is kept
# aside. The proper fix is to put the dimension in the object name, as higtree
# already does; that belongs to the build-system change, not here.
#
# The order — 3D first, 2D last — is not arbitrary either. The example Makefiles
# link ../src/hig-flow-*.o directly rather than libhigflow<dim>d.a, so whichever
# dimension is built last decides which examples can be compiled. Nine of the
# eleven examples are 2D, so 2D is what is left in place, and .built-dim records
# it so the entrypoint can rebuild when a 3D case is asked for.
RUN cd higflow \
    && make DIM=3 \
    && cp lib/libhigflow3d.a /tmp/libhigflow3d.a \
    && rm -f src/*.o \
    && make DIM=2 \
    && cp /tmp/libhigflow3d.a lib/libhigflow3d.a \
    && rm -f /tmp/libhigflow3d.a \
    && echo 2 > .built-dim \
    && ls -l lib/

# /work is where a mounted host directory lands, so cases and their output stay
# on the host and the container stays disposable.
#
# Deliberately not declared as a VOLUME: that makes Docker create an anonymous
# volume whenever nothing is mounted there, and those accumulate silently for
# anyone who forgets --rm. A bind mount works without it.
WORKDIR /work

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["help"]
