#!/usr/bin/env python3
"""Exact solution and error measures for plane Poiseuille flow.

The case is `higflow/example2d_Newt`: a channel of length L = 8 and half-height
h = 1, driven by a parabolic inlet profile, at Re = 1.

At steady state the solution is known in closed form, which is what makes this
case worth verifying against rather than merely looking at:

    u(y) = u_max (1 - (y/h)^2),     v = 0,     dp/dx = -2 mu u_max / h^2

The pressure gradient follows from the momentum balance: for fully developed
flow the inertial terms vanish and mu d2u/dy2 = dp/dx, and d2u/dy2 = -2u_max/h^2
is constant. With mu = 1, u_max = 1.5 and h = 1 that is dp/dx = -3.

The exact solution does not depend on x, so it holds at every station once the
flow is developed. That is what lets the error be measured over the whole field
rather than on one profile.
"""

from __future__ import annotations

import numpy as np

U_MAX = 1.5
HALF_HEIGHT = 1.0
LENGTH = 8.0
MU = 1.0


def exact_u(y, u_max=U_MAX, h=HALF_HEIGHT):
    return u_max * (1.0 - (y / h) ** 2)


def exact_pressure_gradient(u_max=U_MAX, h=HALF_HEIGHT, mu=MU):
    return -2.0 * mu * u_max / h**2


def exact_flow_rate(u_max=U_MAX, h=HALF_HEIGHT):
    """Integral of u over the cross-section: 4 u_max h / 3."""
    return 4.0 * u_max * h / 3.0


def _cell_geometry(grid):
    polys = grid.polygons
    cx = polys[:, :, 0].mean(axis=1)
    cy = polys[:, :, 1].mean(axis=1)
    dy = polys[:, :, 1].max(axis=1) - polys[:, :, 1].min(axis=1)
    dx = polys[:, :, 0].max(axis=1) - polys[:, :, 0].min(axis=1)
    return cx, cy, dx, dy


def field_error(grid):
    """Error of u against the exact profile, over every cell in the domain.

    The norm is area-weighted, so it does not change meaning when the mesh is
    refined, which is what a convergence study needs.
    """
    cx, cy, dx, dy = _cell_geometry(grid)
    u = grid.as_cell("vel")[:, 0]
    err = u - exact_u(cy)

    area = dx * dy
    total = area.sum()
    l2 = float(np.sqrt(np.sum(err**2 * area) / total))
    linf = float(np.max(np.abs(err)))
    return {
        "l2": l2,
        "linf": linf,
        "relative_l2": l2 / U_MAX,
        "cells": int(err.size),
    }


def profile_error(grid, x_fraction=0.75):
    """Same comparison, on one column of cells, for plotting."""
    cx, cy, _, _ = _cell_geometry(grid)
    x_target = cx.min() + x_fraction * (cx.max() - cx.min())
    column = np.isclose(cx, cx[np.argmin(np.abs(cx - x_target))])
    y = cy[column]
    order = np.argsort(y)
    y = y[order]
    u = grid.as_cell("vel")[:, 0][column][order]
    err = u - exact_u(y)
    return {
        "x": float(cx[column][0]),
        "y": y,
        "u": u,
        "exact": exact_u(y),
        "l2": float(np.sqrt(np.mean(err**2))),
        "linf": float(np.max(np.abs(err))),
    }


def flow_rates(grid, stations=7):
    """Volumetric flow rate at several stations along the channel.

    In an incompressible flow this is the same everywhere, so the spread across
    stations measures how well discrete mass conservation is holding.
    """
    cx, cy, _, dy = _cell_geometry(grid)
    u = grid.as_cell("vel")[:, 0]
    xs = np.unique(np.round(cx, 9))
    picks = xs[np.linspace(0, xs.size - 1, stations).astype(int)]

    out = []
    for x in picks:
        m = np.isclose(cx, x)
        out.append({"x": float(x), "Q": float(np.sum(u[m] * dy[m])), "cells": int(m.sum())})
    q = np.array([r["Q"] for r in out])
    return {
        "stations": out,
        "mean": float(q.mean()),
        "spread": float(q.max() - q.min()),
        "relative_spread": float((q.max() - q.min()) / abs(q.mean())),
        "relative_error": float(abs(q.mean() - exact_flow_rate()) / exact_flow_rate()),
    }


def pressure_gradient(grid):
    """Least-squares slope of the cell-centre pressure against x."""
    cx, _, _, _ = _cell_geometry(grid)
    p = grid.as_cell("p")
    slope, _ = np.polyfit(cx, p, 1)
    exact = exact_pressure_gradient()
    return {
        "measured": float(slope),
        "exact": float(exact),
        "relative_error": float(abs(slope - exact) / abs(exact)),
    }


def observed_order(h_values, errors):
    """Convergence order between consecutive mesh refinements.

    p = log(e1/e2) / log(h1/h2). Reported per pair rather than as a single
    fitted number, because a pair that disagrees with the others says something
    the fit would hide.
    """
    h = np.asarray(h_values, dtype=float)
    e = np.asarray(errors, dtype=float)
    idx = np.argsort(-h)                       # coarse to fine
    h, e = h[idx], e[idx]
    pairs = []
    for i in range(len(h) - 1):
        if e[i + 1] <= 0 or e[i] <= 0:
            pairs.append(None)
            continue
        pairs.append(float(np.log(e[i] / e[i + 1]) / np.log(h[i] / h[i + 1])))
    return {"h": h.tolist(), "errors": e.tolist(), "orders": pairs}
