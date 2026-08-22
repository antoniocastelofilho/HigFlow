#!/usr/bin/env bash
# Entry point for the HigFlow container.
#
# Dispatches a small set of verbs and otherwise gets out of the way: anything
# unrecognised is executed as-is, so `docker run ... higflow bash` and
# `docker run ... higflow mpirun -n 4 ./ns-example ...` both work.

set -uo pipefail

: "${HIGFLOW_ROOT:=/opt/higflow}"
: "${HIGFLOW_DIR:=${HIGFLOW_ROOT}/higflow}"

DIM_MARKER="${HIGFLOW_DIR}/.built-dim"

usage() {
    cat <<'EOF'
HigFlow container

  docker run --rm -v "$PWD/cases:/work" higflow:latest <command>

Commands

  help                 this message
  info                 versions and paths of the toolchain in this image
  shell                interactive bash in /work
  examples             list the example cases built into the image
  case <name> [np]     build and run example <name> on <np> ranks (default 1),
                       writing all output to /work/<name>/ on the host
  extract <name>       copy example <name> to /work/<name>/ so you can read or
                       edit it (see the note it prints about building it)

Anything else is run verbatim, so the image is usable as a plain toolbox:

  docker run --rm -v "$PWD/cases:/work" higflow:latest \
      bash -c 'cd /opt/higflow/higflow/example2d_Newt && make && make run NP=2'

Notes

  /work is the working directory and the place to mount a host directory.
  Nothing written outside it survives the container.

  The container runs as uid 1000. If the directory you mount is owned by a
  different user, it cannot write there - run with --user "$(id -u):$(id -g)"
  or chown the directory.

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
    ls -1 "${HIGFLOW_DIR}/lib" 2>/dev/null | sed 's/^/    /' || true
    echo
    echo "  higflow/src objects currently compiled for DIM=$(cat "$DIM_MARKER" 2>/dev/null || echo '?')"
}

examples() {
    echo "Example cases in this image:"
    for d in "${HIGFLOW_DIR}"/example*/; do
        [ -d "$d" ] || continue
        local n dim
        n="$(basename "$d")"
        dim="$(case_dim "$n")"
        printf '  %-28s %sD\n' "$n" "$dim"
    done
    echo
    echo "Run one with:  case <name> [np]"
}

# The dimension an example needs, taken from its name. Each example Makefile
# hardcodes DIM, and the names are consistent with it.
case_dim() {
    case "$1" in
        *3d*|*3D*) echo 3 ;;
        *)         echo 2 ;;
    esac
}

require_writable_work() {
    if [ ! -d /work ]; then
        echo "error: /work does not exist inside the container" >&2
        echo "  mount a host directory there:  -v \"\$PWD/cases:/work\"" >&2
        exit 2
    fi
    if [ ! -w /work ]; then
        echo "error: /work is not writable by this container" >&2
        echo >&2
        echo "  the container runs as uid $(id -u), and the mounted directory is owned by" >&2
        echo "  uid $(stat -c %u /work 2>/dev/null || echo '?')." >&2
        echo >&2
        echo "  either give the directory to that uid on the host:" >&2
        echo "      chown $(id -u):$(id -g) <the directory you mounted>" >&2
        echo >&2
        echo "  or run the container as yourself:" >&2
        echo "      docker run --rm --user \"\$(id -u):\$(id -g)\" -v ... higflow:latest ..." >&2
        exit 2
    fi
}

# The example Makefiles link ../src/*.o directly rather than libhigflow<dim>d.a,
# and higflow/Makefile does not put the dimension in its object file names. So
# there is one shared set of objects and it carries whatever dimension was built
# last. A 3D example therefore needs those objects rebuilt before it will link.
ensure_dim() {
    local want="$1"
    local have
    have="$(cat "$DIM_MARKER" 2>/dev/null || echo '')"
    if [ "$have" = "$want" ]; then
        return 0
    fi
    echo "==> higflow/src objects are compiled for DIM=${have:-unknown}; this case needs DIM=${want}"
    echo "    rebuilding them (a few seconds)"
    ( cd "$HIGFLOW_DIR" && rm -f src/*.o && make "DIM=${want}" >/dev/null ) || {
        echo "error: rebuilding higflow objects for DIM=${want} failed" >&2
        exit 1
    }
    echo "$want" > "$DIM_MARKER" 2>/dev/null || true
}

extract() {
    local name="${1:-}"
    [ -n "$name" ] || { echo "error: which case? run 'examples' to list them" >&2; exit 2; }
    local src="${HIGFLOW_DIR}/${name}"
    [ -d "$src" ] || { echo "error: no example named '${name}'" >&2; exit 2; }
    require_writable_work

    cp -r "$src" "/work/${name}"
    echo "copied ${name} to /work/${name}/"
    echo
    echo "Note: this copy cannot be built where it now sits. Every example"
    echo "includes ../src/hig-flow-*.h and links ../src/hig-flow-*.o, so it only"
    echo "compiles inside the source tree. To edit and build, use the development"
    echo "image with the repository mounted:"
    echo
    echo "    docker run --rm -it -v \"\$PWD:/src\" higflow:dev"
    echo "    cd /src/higflow/${name} && make && make run"
}

run_case() {
    local name="${1:-}"
    local np="${2:-1}"

    [ -n "$name" ] || { echo "error: which case? run 'examples' to list them" >&2; exit 2; }
    local src="${HIGFLOW_DIR}/${name}"
    [ -d "$src" ] || { echo "error: no example named '${name}'" >&2; exit 2; }

    require_writable_work
    ensure_dim "$(case_dim "$name")"

    # The case is built and run where it lives, because it cannot be moved: the
    # header includes and the Makefile both reach into ../src by relative path.
    # Only the output directories are redirected, by making them symlinks into
    # the mounted /work, so every file the run produces lands on the host while
    # the case itself stays in the image.
    local out="/work/${name}"
    mkdir -p "${out}/VTKS" "${out}/DATA" "${out}/output"

    cd "$src" || exit 1
    for d in VTKS DATA output; do
        rm -rf "./$d"
        ln -s "${out}/${d}" "./$d"
    done

    echo "==> building ${name}"
    if ! make; then
        echo "error: build failed" >&2
        exit 1
    fi

    echo "==> running on ${np} rank(s)"
    if ! make run NP="$np"; then
        echo "error: run failed" >&2
        exit 1
    fi

    echo
    echo "==> done. Output is on the host under ./${name}/"
    find "$out" -type f -name '*.vtk' | wc -l | sed 's/^/    vtk files: /'
}

cmd="${1:-help}"
case "$cmd" in
    help|--help|-h) usage ;;
    info)           info ;;
    examples)       examples ;;
    shell)          exec bash ;;
    case)           shift; run_case "$@" ;;
    extract)        shift; extract "$@" ;;
    *)              exec "$@" ;;
esac
