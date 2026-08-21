#!/usr/bin/env python3
"""Render HigFlow VTK output into the figures used by the gallery.

Every figure in docs/gallery.md is produced by this script, so anyone can
regenerate them from a run rather than take them on trust.

    python3 tools/gallery/render.py --input /path/to/case/VTKS --out docs/images/gallery

It reads the ASCII UNSTRUCTURED_GRID files HigFlow writes directly and draws the
quadrilateral cells as they are, with no interpolation or resampling, so what the
figure shows is what the solver computed.

Dependencies: numpy and matplotlib. No VTK library, no ParaView, no display.
"""

from __future__ import annotations

import argparse
import os
import re
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")  # headless: no display needed, and none is available in CI
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.ticker import MaxNLocator

# ── design tokens ────────────────────────────────────────────────────────────
# The two categorical hues were checked with the palette validator: chroma above
# the gray floor, adjacent-pair separation ΔE 23.0 under protanopia and 29.2 for
# normal vision, contrast above 3:1 against the surface.
INK = "#17222E"
MUTED = "#6B7A8B"
FAINT = "#C8D2DB"
SURFACE = "#FCFCFB"
BLUE = "#0F70B7"
AMBER = "#C56A18"

# Sequential ramps are single-hue, light to dark. Deliberately not a rainbow:
# a rainbow map invents boundaries where the field is smooth and hides them
# where it is not, which is why it is the wrong default for a scalar field.
SEQ_BLUE = LinearSegmentedColormap.from_list(
    "seq_blue", ["#F2F7FB", "#CFE2F0", "#8FBEDD", "#3E8CC4", BLUE, "#0A4F80", "#072F4D"]
)
SEQ_AMBER = LinearSegmentedColormap.from_list(
    "seq_amber", ["#FDF6EE", "#F6E1C6", "#E8BC85", "#D89546", AMBER, "#8E4A11", "#5C2F09"]
)
# Diverging carries a neutral midpoint, so zero is visibly neutral rather than
# being another colour competing with the poles.
DIV_BA = LinearSegmentedColormap.from_list(
    "div_ba", ["#072F4D", BLUE, "#8FBEDD", "#F0F0EE", "#E8BC85", AMBER, "#5C2F09"]
)

plt.rcParams.update(
    {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.edgecolor": FAINT,
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.linewidth": 0.8,
        "axes.grid": False,
        "legend.frameon": False,
    }
)


# ── VTK reader ───────────────────────────────────────────────────────────────
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


# ── figures ──────────────────────────────────────────────────────────────────
def _style_axes(ax, xlabel, ylabel):
    ax.set_xlabel(xlabel, color=MUTED, fontsize=8.5)
    ax.set_ylabel(ylabel, color=MUTED, fontsize=8.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(FAINT)


def _titles(ax, title, subtitle=""):
    """Title and subtitle stacked above the axes.

    Offsets are in points rather than axes fractions. A field plot with
    `set_aspect("equal")` on a long thin domain has a very short axes box, and a
    fractional offset there collapses to a couple of pixels - which is how the
    first version of these figures ended up printing the subtitle on top of the
    title.
    """
    ax.annotate(title, xy=(0, 1), xycoords="axes fraction",
                xytext=(0, 24 if subtitle else 8), textcoords="offset points",
                fontsize=11.5, fontweight="600", color=INK, va="bottom", ha="left")
    if subtitle:
        ax.annotate(subtitle, xy=(0, 1), xycoords="axes fraction",
                    xytext=(0, 8), textcoords="offset points",
                    fontsize=8.5, color=MUTED, va="bottom", ha="left")


def field_figure(grid, values, title, subtitle, bar_label, cmap, out,
                 diverging=False, width=8.4, height=None):
    polys = grid.polygons
    xmin, ymin = polys.reshape(-1, 2).min(axis=0)
    xmax, ymax = polys.reshape(-1, 2).max(axis=0)
    aspect = (ymax - ymin) / (xmax - xmin)
    if height is None:
        height = max(2.0, width * aspect + 1.15)

    fig, ax = plt.subplots(figsize=(width, height))

    if diverging:
        lim = float(np.max(np.abs(values)))
        norm = TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim) if lim > 0 else None
    else:
        norm = None

    pc = PolyCollection(polys, array=values, cmap=cmap, norm=norm,
                        edgecolors="none", linewidths=0)
    ax.add_collection(pc)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal")
    _style_axes(ax, "x", "y")

    _titles(ax, title, subtitle)

    cb = fig.colorbar(pc, ax=ax, fraction=0.028, pad=0.015)
    cb.set_label(bar_label, color=MUTED, fontsize=8.5)
    cb.ax.tick_params(labelsize=7.5, colors=MUTED, length=2)
    cb.outline.set_edgecolor(FAINT)
    cb.outline.set_linewidth(0.6)

    fig.tight_layout()
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return out


def poiseuille_validation(grid, out, u_max=1.5, half_height=1.0):
    """Numerical profile against the exact plane-Poiseuille solution.

    With no-slip at y = ±h and the parabolic inlet the driver imposes, the
    fully developed solution is u(y) = u_max (1 - (y/h)^2). Sampling a column
    near the outlet, where the flow has developed, is therefore a real check
    and not a restatement of the boundary condition.
    """
    cx, cy = grid.centroids[:, 0], grid.centroids[:, 1]
    u = grid.as_cell("vel")[:, 0]

    # A column three quarters of the way downstream.
    x_target = cx.min() + 0.75 * (cx.max() - cx.min())
    col = np.isclose(cx, cx[np.argmin(np.abs(cx - x_target))], atol=1e-9)
    y = cy[col]
    un = u[col]
    order = np.argsort(y)
    y, un = y[order], un[order]

    exact = u_max * (1.0 - (y / half_height) ** 2)
    err = un - exact
    l2 = float(np.sqrt(np.mean(err ** 2)))
    linf = float(np.max(np.abs(err)))
    rel = l2 / u_max

    fig, (ax, axe) = plt.subplots(
        1, 2, figsize=(9.2, 3.9), gridspec_kw={"width_ratios": [1.5, 1]}
    )

    ys = np.linspace(-half_height, half_height, 400)
    ax.plot(u_max * (1 - (ys / half_height) ** 2), ys,
            color=AMBER, lw=2.0, zorder=2, label="Exact  $u = u_{max}(1-y^2)$")
    ax.plot(un, y, linestyle="none", marker="o", ms=4.2,
            markerfacecolor=BLUE, markeredgecolor=SURFACE, markeredgewidth=0.7,
            zorder=3, label="HigFlow, cell centres")
    ax.axhline(0, color=FAINT, lw=0.7, zorder=1)
    _style_axes(ax, "streamwise velocity $u$", "$y$")
    _titles(ax, "Plane Poiseuille flow", f"profile at x = {x_target:.2f}, Re = 1")
    # Upper right, because the profile puts its maximum at y = 0 and runs to
    # zero at both walls, leaving that corner empty. Upper left sits on the curve.
    leg = ax.legend(frameon=False, fontsize=8.5, loc="upper right",
                    handletextpad=0.6, borderaxespad=0.4)
    for t in leg.get_texts():
        t.set_color(INK)

    axe.plot(err, y, color=BLUE, lw=1.8)
    axe.axvline(0, color=FAINT, lw=0.7)
    _style_axes(axe, "numerical − exact", "$y$")
    _titles(axe, "Pointwise error", "against the exact profile")

    # Four ticks, not eight: at 1e-3 magnitudes the default locator packs this
    # narrow panel with labels that overlap each other.
    axe.xaxis.set_major_locator(MaxNLocator(4))
    axe.ticklabel_format(axis="x", style="sci", scilimits=(-2, 2), useMathText=True)
    axe.xaxis.get_offset_text().set(color=MUTED, size=7.5)

    # The norms go inside the panel. As a subtitle they ran off its right edge.
    axe.text(0.05, 0.05,
             f"$L_2$ = {l2:.2e}\n$L_\\infty$ = {linf:.2e}\nrelative $L_2$ = {rel:.2e}",
             transform=axe.transAxes, fontsize=8, color=MUTED,
             va="bottom", ha="left", linespacing=1.6)

    fig.tight_layout()
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return out, {"l2": l2, "linf": linf, "relative_l2": rel,
                 "x": float(x_target), "samples": int(y.size)}


def mesh_figure(grid, title, subtitle, out, width=8.4):
    polys = grid.polygons
    xmin, ymin = polys.reshape(-1, 2).min(axis=0)
    xmax, ymax = polys.reshape(-1, 2).max(axis=0)
    aspect = (ymax - ymin) / (xmax - xmin)
    fig, ax = plt.subplots(figsize=(width, max(2.0, width * aspect + 0.9)))
    pc = PolyCollection(polys, facecolors="none", edgecolors=BLUE, linewidths=0.35)
    ax.add_collection(pc)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal")
    _style_axes(ax, "x", "y")
    _titles(ax, title, subtitle)
    fig.tight_layout()
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return out


def latest_vtk(directory: str) -> str:
    files = [f for f in os.listdir(directory) if f.endswith(".vtk")]
    if not files:
        raise SystemExit(f"no .vtk files in {directory}")

    def frame_of(name):
        m = re.search(r"-(\d+)\.vtk$", name)
        return int(m.group(1)) if m else -1

    return os.path.join(directory, max(files, key=frame_of))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, help="a VTKS directory, or a single .vtk file")
    ap.add_argument("--out", required=True, help="directory to write figures into")
    ap.add_argument("--name", default="field", help="basename for the output files")
    ap.add_argument("--kind", default="speed",
                    choices=["speed", "pressure", "fracvol", "n1", "mesh", "poiseuille"])
    ap.add_argument("--title", default="")
    ap.add_argument("--subtitle", default="")
    args = ap.parse_args()

    path = args.input if args.input.endswith(".vtk") else latest_vtk(args.input)
    grid = read_vtk(path)
    os.makedirs(args.out, exist_ok=True)
    dest = os.path.join(args.out, f"{args.name}.png")

    print(f"  source   {path}")
    print(f"  cells    {grid.cells.shape[0]}   points {grid.points.shape[0]}")
    print(f"  point data {sorted(grid.point_data)}")
    print(f"  cell data  {sorted(grid.cell_data)}")

    if args.kind == "speed":
        v = grid.speed()
        field_figure(grid, v, args.title or "Velocity magnitude", args.subtitle,
                     "$|u|$", SEQ_BLUE, dest)
    elif args.kind == "pressure":
        p = grid.as_cell("p")
        field_figure(grid, p - p.mean(), args.title or "Pressure", args.subtitle,
                     "$p - \\bar{p}$", DIV_BA, dest, diverging=True)
    elif args.kind == "fracvol":
        field_figure(grid, grid.as_cell("FracVol"), args.title or "Volume fraction",
                     args.subtitle, "$\\phi$", SEQ_AMBER, dest)
    elif args.kind == "n1":
        n1 = grid.first_normal_stress_difference()
        name, _ = grid.tensor_field()
        print(f"  tensor   {name!r}   N1 range [{n1.min():.4g}, {n1.max():.4g}]")
        field_figure(grid, n1,
                     args.title or "First normal stress difference",
                     args.subtitle, r"$N_1 = \tau_{xx} - \tau_{yy}$",
                     DIV_BA, dest, diverging=True)
    elif args.kind == "mesh":
        mesh_figure(grid, args.title or "Computational grid", args.subtitle, dest)
    elif args.kind == "poiseuille":
        _, stats = poiseuille_validation(grid, dest)
        print(f"  L2 {stats['l2']:.6e}  Linf {stats['linf']:.6e} "
              f"relative L2 {stats['relative_l2']:.6e} over {stats['samples']} cells")

    print(f"  wrote    {dest}")


if __name__ == "__main__":
    sys.exit(main())
