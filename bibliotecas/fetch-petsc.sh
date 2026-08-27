#!/usr/bin/env bash
# Download the PETSc source tarball this project builds against.
#
# The tarball used to be committed to the repository, which made every clone
# 37 MB heavier and put a third-party release under this project's version
# control. It is fetched on demand instead, and verified against the checksum
# of the archive that was previously committed, so the build gets byte-for-byte
# the same source as before.
#
# Idempotent: a present, valid tarball is left alone.
#
# Usage:  bash bibliotecas/fetch-petsc.sh

set -euo pipefail

VERSION="3.14.0"
ARCHIVE="petsc-${VERSION}.tar.gz"
SHA256="bf37688f68889e82b7609f6b36d535899c892bbe7335e93721d42e90ea8b03c3"

# Official Argonne release snapshots. Both mirrors serve the same archive.
MIRRORS=(
    "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/${ARCHIVE}"
    "https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/${ARCHIVE}"
)

cd "$(dirname "${BASH_SOURCE[0]}")"

checksum() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "$1" | awk '{print $1}'
    else
        echo "error: neither sha256sum nor shasum is available" >&2
        exit 1
    fi
}

if [ -f "$ARCHIVE" ]; then
    if [ "$(checksum "$ARCHIVE")" = "$SHA256" ]; then
        echo "${ARCHIVE}: already present and verified"
        exit 0
    fi
    echo "${ARCHIVE}: present but checksum does not match, downloading again" >&2
    mv -f "$ARCHIVE" "${ARCHIVE}.bad"
fi

for url in "${MIRRORS[@]}"; do
    echo "downloading ${url}"
    if command -v curl >/dev/null 2>&1; then
        curl -fL --retry 3 --connect-timeout 30 -o "${ARCHIVE}.part" "$url" && break
    elif command -v wget >/dev/null 2>&1; then
        wget -O "${ARCHIVE}.part" "$url" && break
    else
        echo "error: neither curl nor wget is available" >&2
        exit 1
    fi
    echo "mirror failed, trying the next one" >&2
done

if [ ! -f "${ARCHIVE}.part" ]; then
    echo "error: every mirror failed" >&2
    exit 1
fi

got="$(checksum "${ARCHIVE}.part")"
if [ "$got" != "$SHA256" ]; then
    echo "error: checksum mismatch for ${ARCHIVE}" >&2
    echo "  expected ${SHA256}" >&2
    echo "  got      ${got}" >&2
    echo "The downloaded file is not the archive this project was built against." >&2
    echo "It has been left as ${ARCHIVE}.part for inspection." >&2
    exit 1
fi

mv -f "${ARCHIVE}.part" "$ARCHIVE"
echo "${ARCHIVE}: downloaded and verified"
