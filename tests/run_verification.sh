#!/usr/bin/env bash
# Run the plane channel case and check it against its exact solution.
#
# Exits non-zero when a check fails, so it works as a CTest test or a CI step.
#
# Two ways to run the case, chosen by environment:
#
#   HIGFLOW_IMAGE=higflow:latest   run it in the container (needs docker)
#   HIGFLOW_CASE_DIR=/path/to/case use a case directory that is already built
#
# With neither set it looks for the container image, since that is the route
# the documentation recommends and the one CI uses.

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
CASE="${HIGFLOW_TEST_CASE:-example2d_Newt}"
WORK="${HIGFLOW_TEST_WORK:-$(mktemp -d)}"
IMAGE="${HIGFLOW_IMAGE:-higflow:latest}"

log() { printf '%s\n' "$*"; }
fail() { printf 'error: %s\n' "$*" >&2; exit 1; }

command -v python3 >/dev/null || fail "python3 is required"
python3 -c "import numpy" 2>/dev/null || fail "numpy is required"

log "verification: $CASE"
log "  work directory: $WORK"

if [ -n "${HIGFLOW_CASE_DIR:-}" ]; then
    # A case directory that has already been built and run.
    VTKS="$HIGFLOW_CASE_DIR/VTKS"
    [ -d "$VTKS" ] || fail "no VTKS directory under $HIGFLOW_CASE_DIR"
    log "  using existing output in $VTKS"
else
    command -v docker >/dev/null || fail "docker is required, or set HIGFLOW_CASE_DIR"
    docker image inspect "$IMAGE" >/dev/null 2>&1 \
        || fail "image $IMAGE not found; build it with: docker build -f containers/Dockerfile -t $IMAGE ."

    mkdir -p "$WORK"
    # The container runs as uid 1000 and writes here.
    chmod 0777 "$WORK" 2>/dev/null || true

    log "  running the case in $IMAGE"
    if ! docker run --rm --shm-size=1g -v "$WORK:/work" "$IMAGE" case "$CASE" > "$WORK/run.log" 2>&1; then
        log "  the case failed to run; last lines of its log:"
        tail -20 "$WORK/run.log" | sed 's/^/    /'
        exit 1
    fi
    VTKS="$WORK/$CASE/VTKS"
    [ -d "$VTKS" ] || fail "the case produced no VTKS directory"
fi

log ""
python3 "$HERE/verification/verify.py" check --vtks "$VTKS"
status=$?

log ""
if [ $status -eq 0 ]; then
    log "verification passed"
else
    log "verification failed"
fi
exit $status
