#!/usr/bin/env python3
"""Checks on the verification suite itself, none of which need the solver.

The suite measures the solver against exact solutions. This measures the suite
against things that are true by construction: the exact solution against its
own definition, the order estimator against synthetic data whose order is
known, the mesh generator against the mesh that ships with the example, and the
VTK reader against a file written here.

It exists because the expensive checks cannot run on every push. A container
build takes tens of minutes; this takes under a second, and it still catches
the mistakes that matter, since a wrong error measure would report a wrong
convergence order just as confidently as a right one.

    python3 tests/verification/selftest.py
"""

from __future__ import annotations

import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tools"))
sys.path.insert(0, HERE)

import channel                      # noqa: E402
import channel_mesh                 # noqa: E402
from higflow_vtk import read_vtk    # noqa: E402


class SkipTest(Exception):
    pass


CASES = []


def case(fn):
    CASES.append(fn)
    return fn


# The exact solution against its own definition

@case
def flow_rate_is_the_integral_of_the_profile():
    """4 u_max h / 3 is a closed form; check it against the integral it stands for."""
    y = np.linspace(-channel.HALF_HEIGHT, channel.HALF_HEIGHT, 200001)
    numeric = float(np.trapz(channel.exact_u(y), y))
    assert abs(numeric - channel.exact_flow_rate()) < 1e-9, numeric


@case
def pressure_gradient_follows_from_the_profile():
    """For fully developed flow mu u_yy = dp/dx. Differentiate and compare."""
    y = np.linspace(-1.0, 1.0, 20001)
    d2u = np.gradient(np.gradient(channel.exact_u(y), y), y)
    interior = d2u[10:-10]
    assert float(np.ptp(interior)) < 1e-6, "the second derivative should be constant"
    got = channel.MU * float(interior.mean())
    assert abs(got - channel.exact_pressure_gradient()) < 1e-6, got


@case
def the_profile_vanishes_at_the_walls():
    assert channel.exact_u(channel.HALF_HEIGHT) == 0.0
    assert channel.exact_u(-channel.HALF_HEIGHT) == 0.0
    assert channel.exact_u(0.0) == channel.U_MAX


# The order estimator against data whose order is known

@case
def order_recovers_a_second_order_sequence():
    h = np.array([0.2, 0.1, 0.05, 0.025])
    res = channel.observed_order(h, 0.3 * h**2)
    for p in res["orders"]:
        assert abs(p - 2.0) < 1e-9, res["orders"]


@case
def order_recovers_a_first_order_sequence():
    h = np.array([0.2, 0.1, 0.05])
    res = channel.observed_order(h, 7.0 * h)
    for p in res["orders"]:
        assert abs(p - 1.0) < 1e-9, res["orders"]


@case
def order_is_insensitive_to_the_order_of_its_input():
    h = [0.05, 0.2, 0.1]
    e = [0.3 * x**2 for x in h]
    for p in channel.observed_order(h, e)["orders"]:
        assert abs(p - 2.0) < 1e-9


# The mesh generator against the mesh that ships with the example

@case
def generated_mesh_reproduces_the_shipped_one():
    shipped = os.path.join(ROOT, "higflow", "example2d_Newt", "amrs")
    if not os.path.isdir(shipped):
        raise SkipTest("the example is not present in this checkout")
    with tempfile.TemporaryDirectory() as tmp:
        channel_mesh.write_channel(tmp, nx=160, ny=40)
        names = [("domain", "ch-d-0.amr")] + [("bc", "ch-bc-%d.amr" % i) for i in range(4)]
        for sub, name in names:
            with open(os.path.join(tmp, sub, name)) as fh:
                got = fh.read().split()
            with open(os.path.join(shipped, sub, name)) as fh:
                want = fh.read().split()
            assert got == want, "%s/%s\n  generated %s\n  shipped   %s" % (sub, name, got, want)


@case
def mesh_cell_size_halves_when_the_count_doubles():
    with tempfile.TemporaryDirectory() as tmp:
        coarse = channel_mesh.write_channel(tmp, nx=40, ny=10)
        fine = channel_mesh.write_channel(tmp, nx=80, ny=20)
    assert abs(coarse["dx"] / fine["dx"] - 2.0) < 1e-12
    assert fine["cells"] == 4 * coarse["cells"]


# The reader and the error measures against a field built here

def _synthetic_vtk(path, nx, ny, u_of_y, pressure_slope):
    """Write a channel-shaped VTK the way HigFlow writes one.

    Velocity as POINT_DATA at four unshared corners per cell, pressure as
    CELL_DATA in the same file. That layout is the thing the reader has to cope
    with, so a synthetic file is only useful if it reproduces it.
    """
    dx = channel.LENGTH / nx
    dy = 2.0 * channel.HALF_HEIGHT / ny
    pts, cells, vel = [], [], []
    for i in range(nx):
        for j in range(ny):
            x0 = i * dx
            y0 = -channel.HALF_HEIGHT + j * dy
            corners = [(x0, y0), (x0 + dx, y0), (x0 + dx, y0 + dy), (x0, y0 + dy)]
            cells.append([len(pts) + k for k in range(4)])
            for cx, cy in corners:
                pts.append((cx, cy))
                vel.append((u_of_y(cy), 0.0))

    with open(path, "w", newline="\n") as f:
        f.write("# vtk DataFile Version 3.0\nsynthetic\nASCII\nDATASET UNSTRUCTURED_GRID\n")
        f.write("POINTS %d float\n" % len(pts))
        for x, y in pts:
            f.write("%r %r 0\n" % (x, y))
        f.write("\nCELLS %d %d\n" % (len(cells), 5 * len(cells)))
        for c in cells:
            f.write("4 " + " ".join(str(k) for k in c) + "\n")
        f.write("\nCELL_TYPES %d\n" % len(cells))
        f.write("9\n" * len(cells))
        f.write("\nPOINT_DATA %d\nVECTORS vel FLOAT\n" % len(pts))
        for ux, uy in vel:
            f.write("%r %r 0\n" % (ux, uy))
        f.write("\nCELL_DATA %d\nSCALARS p FLOAT\nLOOKUP_TABLE default\n" % len(cells))
        for c in cells:
            xc = sum(pts[k][0] for k in c) / 4.0
            f.write("%r\n" % (pressure_slope * xc))


def _grid(nx=16, ny=8, u_of_y=None, slope=-3.0):
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "synthetic.vtk")
        _synthetic_vtk(path, nx, ny, u_of_y or channel.exact_u, slope)
        return read_vtk(path)


@case
def reader_splits_point_data_from_cell_data():
    g = _grid(nx=16, ny=8)
    assert g.ncells == 128
    assert g.point_data["vel"].shape == (512, 3)
    assert g.cell_data["p"].shape == (128,)
    assert g.as_cell("p").shape == (128,)
    assert g.as_cell("vel").shape == (128, 3)


@case
def an_exact_field_has_almost_no_error():
    """The measures have to report near zero on the solution they compare against."""
    g = _grid(nx=16, ny=8)
    e = channel.field_error(g)
    assert e["cells"] == 128
    # The corner average of a parabola is not the parabola at the centre. The
    # gap is u_max dy^2 / (4 h^2), which here is 0.023, so this bound checks
    # that nothing larger than that discretisation has crept in.
    assert e["l2"] < 0.03, e


@case
def flow_rate_of_an_exact_field_is_the_exact_one():
    g = _grid(nx=16, ny=32)
    q = channel.flow_rates(g)
    assert q["relative_error"] < 5e-3, q
    assert q["relative_spread"] < 1e-12, q["relative_spread"]


@case
def pressure_gradient_recovers_the_slope_it_was_given():
    g = _grid(nx=16, ny=8, slope=-3.0)
    dp = channel.pressure_gradient(g)
    assert abs(dp["measured"] + 3.0) < 1e-6, dp


@case
def skipping_edge_columns_drops_exactly_those_columns():
    g = _grid(nx=16, ny=8)
    full = channel.field_error(g, skip_edge_columns=0)
    inner = channel.field_error(g, skip_edge_columns=1)
    assert full["cells"] == 128
    assert inner["cells"] == 128 - 2 * 8, inner
    assert inner["skipped"] == 16


@case
def a_wrong_first_column_is_caught_and_then_excluded():
    """The reason skip_edge_columns exists, on a field whose defect we placed."""
    g = _grid(nx=16, ny=8)
    cx = g.polygons[:, :, 0].mean(axis=1)
    first = np.isclose(cx, cx.min())
    g.point_data["vel"][g.cells[first].ravel(), 0] = 2.7

    whole = channel.field_error(g, skip_edge_columns=0)
    inner = channel.field_error(g, skip_edge_columns=1)
    assert whole["l2"] > 10 * inner["l2"], (whole, inner)

    b = channel.boundary_corner_error(g)
    assert b["error_on_boundary"] > 1.0, b


def main():
    width = max(len(c.__name__) for c in CASES)
    failed = skipped = 0
    for c in CASES:
        try:
            c()
        except SkipTest as exc:
            print("  SKIP  %s  %s" % (c.__name__.ljust(width), exc))
            skipped += 1
        except AssertionError as exc:
            print("  FAIL  %s  %s" % (c.__name__.ljust(width), exc))
            failed += 1
        else:
            print("  ok    %s" % c.__name__)

    total = len(CASES)
    print()
    print("  %d passed, %d failed, %d skipped, of %d"
          % (total - failed - skipped, failed, skipped, total))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
