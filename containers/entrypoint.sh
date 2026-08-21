#!/usr/bin/env bash
# Entry point for the HigFlow container.
#
# Dispatches a small set of verbs and otherwise gets out of the way: anything
# unrecognised is executed as-is, so `docker run ... higflow bash` and
# `docker run ... higflow mpirun -n 4 ./ns-example ...` both work.

set -euo pipefail

: "${HIGFLOW_ROOT:=/opt/higflow}"

usage() {
    cat <<'EOF'
HigFlow container

  docker run --rm -v "$PWD:/work" higflow <command>

Commands

  help                 this message
  info                 versions and paths of the toolchain in this image
  shell                interactive bash in /work
  examples             list the example cases built into the image
  case <name> [np]     copy example <name> into /work and run it on <np> ranks
                       (default 1), leaving all output on the host

Anything else is run verbatim, so the image is usable as a plain toolbox:

  docker run --rm -v "$PWD:/work" higflow bash -c 'mpirun -n 2 ./ns-example ...'

Notes

  /work is the working directory and the place to mount a host directory.
  Nothing written outside it survives the container.

  MPI inside a container needs shared memory. If a parallel run dies with an
  error mentioning shm or vader, raise it:  --shm-size=1g
EOF
}

info() {
    echo "HigFlow container"
    echo
    echo "  root          ${HIGFLOW_ROOT}"
    echo "  HIGTREE_DIR   ${HIGTREE_DIR:-unset}"
    echo "  HIGFLOW_DIR   ${HIGFLOW_DIR:-unset}"
    echo "  PETSC_DIR     ${PETSC_DIR:-unset}"
    echo
    echo "  gcc           $(gcc -dumpversion 2>/dev/null || echo '-')"
    echo "  mpicc         $(mpicc -showme:version 2>&1 | head -1 || echo '-')"
    echo "  mpirun        $(mpirun --version 2>/dev/null | head -1 || echo '-')"
    echo
    echo "  libraries built:"
    ls -1 "${HIGFLOW_ROOT}/higtree/lib" 2>/dev/null | sed 's/^/    /' || echo "    none"
    ls -1 "${HIGFLOW_ROOT}/higflow/lib" 2>/dev/null | sed 's/^/    /' || true
}

examples() {
    echo "Example cases in this image:"
    for d in "${HIGFLOW_ROOT}"/higflow/example*/; do
        [ -d "$d" ] || continue
        printf '  %s\n' "$(basename "$d")"
    done
    echo
    echo "Run one with:  case <name> [np]"
}

run_case() {
    local name="${1:-}"
    local np="${2:-1}"

    if [ -z "$name" ]; then
        echo "error: which case? run 'examples' to list them" >&2
        exit 2
    fi

    local src="${HIGFLOW_ROOT}/higflow/${name}"
    if [ ! -d "$src" ]; then
        echo "error: no example named '${name}'" >&2
        echo "run 'examples' to list them" >&2
        exit 2
    fi

    # Copied into /work rather than run in place, so every artefact — the
    # rebuilt binary, the VTKs, the restart files — lands on the host where the
    # user can see it, and the image stays untouched.
    local dst="/work/${name}"
    if [ -d "$dst" ]; then
        echo "==> ${dst} already exists, reusing it"
    else
        echo "==> copying ${name} into /work"
        cp -r "$src" "$dst"
    fi

    cd "$dst"
    echo "==> building"
    make
    echo "==> running on ${np} rank(s)"
    make run NP="$np"
    echo
    echo "==> done. Output is on the host under ./${name}/"
}

cmd="${1:-help}"
case "$cmd" in
    help|--help|-h) usage ;;
    info)           info ;;
    examples)       examples ;;
    shell)          exec bash ;;
    case)           shift; run_case "$@" ;;
    *)              exec "$@" ;;
esac
