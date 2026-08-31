#!/usr/bin/env python3
"""Reader for the VTK files HigFlow writes.

Shared by the gallery renderer and the verification suite, so the two cannot
disagree about what a field means.

HigFlow writes ASCII UNSTRUCTURED_GRID. Two things about the format matter:
velocity is POINT_DATA and pressure is CELL_DATA in the same file, so the two
have different lengths, and the points are not shared between cells, so a
"point" field is really a value per cell corner.
"""

from __future__ import annotations

import os
import re

import numpy as np


# VTK reader
class Grid:
    """A 2D unstructured grid of quadrilaterals with its point and cell data.

    HigFlow writes velocity as POINT_DATA and pressure as CELL_DATA in the same
    file, so the two have different lengths - 25 600 against 6 400 for the
    channel case. The points are not shared between cells either: each cell
    carries its own four corners, so a "point" field is really a value per cell
    corner. `as_cell` collapses one to a per-cell value by averaging the four.
    """

    def __init__(self, points, cells, point_data, cell_data):
        self.points = points        # (npoints, 3)
        self.cells = cells          # (ncells, 4) indices
        self.point_data = point_data
        self.cell_data = cell_data

    @property
    def ncells(self):
        return self.cells.shape[0]

    @property
    def polygons(self):
        return self.points[self.cells][:, :, :2]

    @property
    def centroids(self):
        return self.polygons.mean(axis=1)

    def as_cell(self, name):
        """A per-cell array for `name`, wherever the field happens to live."""
        if name in self.cell_data:
            return self.cell_data[name]
        if name in self.point_data:
            v = self.point_data[name]
            return v[self.cells].mean(axis=1)   # works for scalars, vectors and tensors
        raise KeyError(
            f"no field named {name!r}; "
            f"point data {sorted(self.point_data)}, cell data {sorted(self.cell_data)}"
        )

    def speed(self):
        v = self.as_cell("vel")
        return np.hypot(v[..., 0], v[..., 1])

    def tensor_field(self):
        """The polymer stress tensor, whatever it is called in this file.

        The writer names it with subscript characters, so it is found by shape
        rather than by matching a literal name.
        """
        for store in (self.cell_data, self.point_data):
            for name, arr in store.items():
                if getattr(arr, "ndim", 0) == 3 and arr.shape[1:] == (3, 3):
                    return name, self.as_cell(name)
        raise KeyError("no 3x3 tensor field in this file - is it a viscoelastic run?")

    def first_normal_stress_difference(self):
        """N1 = tau_xx - tau_yy.

        The reason to plot this rather than the velocity: N1 is identically zero
        for a Newtonian fluid and non-zero for a viscoelastic one, so it shows
        the thing the model was added for. A velocity field alone looks much the
        same either way.
        """
        _, t = self.tensor_field()
        return t[:, 0, 0] - t[:, 1, 1]


def _numbers_after(text, start, count):
    return np.fromstring(" ".join(text[start:].split()[:count]), sep=" ")


def read_vtk(path: str) -> Grid:
    with open(path, "r", errors="replace") as fh:
        text = fh.read()

    m = re.search(r"^POINTS\s+(\d+)\s+\w+\s*$", text, re.M)
    if not m:
        raise ValueError(f"{path}: no POINTS section")
    npoints = int(m.group(1))
    points = _numbers_after(text, m.end(), npoints * 3).reshape(npoints, 3)

    m = re.search(r"^CELLS\s+(\d+)\s+(\d+)\s*$", text, re.M)
    if not m:
        raise ValueError(f"{path}: no CELLS section")
    ncells, total = int(m.group(1)), int(m.group(2))
    raw = _numbers_after(text, m.end(), total).astype(int)
    stride = total // ncells          # "count i0 i1 …" per cell; quads here
    cells = raw.reshape(ncells, stride)[:, 1:5]

    # Walk the attribute sections in order, so each field is sized by whichever
    # of POINT_DATA / CELL_DATA it falls under.
    point_data, cell_data = {}, {}
    marks = [
        (mm.start(), mm.group(1), int(mm.group(2)))
        for mm in re.finditer(r"^(POINT_DATA|CELL_DATA)\s+(\d+)\s*$", text, re.M)
    ]

    def owner(pos):
        current = (None, 0)
        for start, kind, n in marks:
            if start < pos:
                current = (kind, n)
        return current

    for sm in re.finditer(
        r"^SCALARS\s+(\S+)\s+\w+(?:\s+\d+)?\s*$\s*^LOOKUP_TABLE\s+\S+\s*$", text, re.M
    ):
        kind, n = owner(sm.start())
        if kind is None:
            continue
        vals = _numbers_after(text, sm.end(), n)
        if vals.size == n:
            (point_data if kind == "POINT_DATA" else cell_data)[sm.group(1)] = vals

    for vm in re.finditer(r"^VECTORS\s+(\S+)\s+\w+\s*$", text, re.M):
        kind, n = owner(vm.start())
        if kind is None:
            continue
        vals = _numbers_after(text, vm.end(), n * 3)
        if vals.size == n * 3:
            (point_data if kind == "POINT_DATA" else cell_data)[vm.group(1)] = vals.reshape(n, 3)

    # A viscoelastic run adds the polymer stress as a 3x3 tensor per point.
    # The field name carries subscripts, so the pattern cannot assume ASCII.
    for tm in re.finditer(r"^TENSORS\s+(\S+)\s+\w+\s*$", text, re.M):
        kind, n = owner(tm.start())
        if kind is None:
            continue
        vals = _numbers_after(text, tm.end(), n * 9)
        if vals.size == n * 9:
            (point_data if kind == "POINT_DATA" else cell_data)[tm.group(1)] = vals.reshape(n, 3, 3)

    return Grid(points, cells, point_data, cell_data)


def latest_vtk(directory: str) -> str:
    """The highest-numbered .vtk in a directory, which is the last frame."""
    files = [f for f in os.listdir(directory) if f.endswith(".vtk")]
    if not files:
        raise SystemExit(f"no .vtk files in {directory}")

    def frame_of(name):
        m = re.search(r"-(\d+)\.vtk$", name)
        return int(m.group(1)) if m else -1

    return os.path.join(directory, max(files, key=frame_of))


def frames(directory: str):
    """Every .vtk in a directory, in frame order, as (frame number, path)."""
    out = []
    for f in os.listdir(directory):
        if not f.endswith(".vtk"):
            continue
        m = re.search(r"-(\d+)\.vtk$", f)
        if m:
            out.append((int(m.group(1)), os.path.join(directory, f)))
    return sorted(out)
