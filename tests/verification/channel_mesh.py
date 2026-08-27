#!/usr/bin/env python3
"""Generate the mesh files for the plane channel at a chosen resolution.

A grid convergence study needs the same case on several meshes, and the meshes
that ship with `example2d_Newt` are a single fixed 160 x 40. This writes the
domain and the four boundary patches for any resolution, in the same format.

The .amr format is positional and has no header. Four lines:

    xmin xmax ymin ymax     the bounding box
    nlevels                 refinement levels
    dx dy 1                 cell size
    i0 j0 nx ny             first cell index and the count per direction

A boundary patch is the same thing with the box collapsed in one direction, so
one of dx, dy is zero and the matching count is 1.
"""

from __future__ import annotations

import argparse
import os
import sys

# The channel example2d_Newt solves: length 8, half-height 1.
X0, X1 = 0.0, 8.0
Y0, Y1 = -1.0, 1.0


def _amr(box, dx, dy, nx, ny):
    xmin, xmax, ymin, ymax = box
    return (
        f"{xmin} {xmax} {ymin} {ymax}\n"
        f"1\n"
        f"{dx} {dy} 1\n"
        f"1 1 {nx} {ny}\n"
    )


def write_channel(out_dir: str, nx: int, ny: int) -> dict:
    """Write domain and boundary meshes for an nx by ny channel.

    Returns the cell size, which is what a convergence study plots against.
    """
    dx = (X1 - X0) / nx
    dy = (Y1 - Y0) / ny

    os.makedirs(os.path.join(out_dir, "domain"), exist_ok=True)
    os.makedirs(os.path.join(out_dir, "bc"), exist_ok=True)

    def put(path, text):
        with open(path, "w", newline="\n") as fh:
            fh.write(text)

    put(os.path.join(out_dir, "domain", "ch-d-0.amr"),
        _amr((X0, X1, Y0, Y1), dx, dy, nx, ny))

    # The boundary identifiers match what example-2d.load.bc.yaml expects, and
    # the driver's get_boundary_velocity switches on them: 0 is the inlet with
    # the parabolic profile, 2 the outlet, 1 and 3 the walls.
    put(os.path.join(out_dir, "bc", "ch-bc-0.amr"),
        _amr((X0, X0, Y0, Y1), 0.0, dy, 1, ny))          # inlet,  x = 0
    put(os.path.join(out_dir, "bc", "ch-bc-1.amr"),
        _amr((X0, X1, Y1, Y1), dx, 0.0, nx, 1))          # top wall,   y = +1
    put(os.path.join(out_dir, "bc", "ch-bc-2.amr"),
        _amr((X1, X1, Y0, Y1), 0.0, dy, 1, ny))          # outlet, x = 8
    put(os.path.join(out_dir, "bc", "ch-bc-3.amr"),
        _amr((X0, X1, Y0, Y0), dx, 0.0, nx, 1))          # bottom wall, y = -1

    return {"nx": nx, "ny": ny, "dx": dx, "dy": dy, "cells": nx * ny}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True, help="directory to write amrs/ into")
    ap.add_argument("--nx", type=int, required=True)
    ap.add_argument("--ny", type=int, required=True)
    args = ap.parse_args()

    info = write_channel(args.out, args.nx, args.ny)
    print(f"  {info['nx']} x {info['ny']}  dx = {info['dx']:g}  dy = {info['dy']:g}  "
          f"{info['cells']} cells  ->  {args.out}")


if __name__ == "__main__":
    sys.exit(main())
