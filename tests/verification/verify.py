#!/usr/bin/env python3
"""Verification checks for the plane channel case.

Three subcommands, meant to be used in this order:

    transient   error against the exact solution at every output frame, to
                see whether the run has reached steady state
    check       the checks on one run, with tolerances; exits non-zero on
                failure, so it can be a CTest test
    order       convergence order across several runs at different resolutions

`order` is only meaningful once `transient` shows the error has stopped
falling. If a run is still developing in time, refining the mesh will not
reduce its error, and the measured order will say more about the time
integration than about the spatial scheme.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from higflow_vtk import read_vtk, latest_vtk, frames   # noqa: E402
import channel                                          # noqa: E402


def cmd_transient(args):
    fs = frames(args.vtks)
    if not fs:
        raise SystemExit(f"no .vtk files in {args.vtks}")
    if args.every > 1:
        fs = fs[:: args.every] + ([fs[-1]] if fs[-1] not in fs[:: args.every] else [])

    print(f"  frames: {len(fs)} of {len(frames(args.vtks))}")
    print()
    print("   frame        L2          Linf     change vs previous")
    prev = None
    rows = []
    for n, path in fs:
        e = channel.field_error(read_vtk(path))
        change = "" if prev is None else f"{(e['l2'] - prev) / prev:+.2%}"
        print(f"  {n:6d}   {e['l2']:.6e}  {e['linf']:.6e}   {change:>10}")
        rows.append({"frame": n, **e})
        prev = e["l2"]

    if len(rows) >= 3:
        last = [r["l2"] for r in rows[-3:]]
        drift = (max(last) - min(last)) / max(last)
        print()
        print(f"  spread over the last three frames: {drift:.2%}")
        if drift < args.steady_tol:
            print("  the error has stopped changing: the run is at steady state")
        else:
            print("  the error is still moving: this run is NOT at steady state,")
            print("  so a convergence study on it would measure the transient")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(rows, fh, indent=2)
    return 0


def cmd_check(args):
    grid = read_vtk(latest_vtk(args.vtks))
    whole = channel.field_error(grid, skip_edge_columns=0)
    inner = channel.field_error(grid, skip_edge_columns=args.skip_edges)
    q = channel.flow_rates(grid, skip_edge_columns=args.skip_edges)
    q_all = channel.flow_rates(grid, skip_edge_columns=0)
    dp = channel.pressure_gradient(grid)

    print(f"  source           {os.path.basename(latest_vtk(args.vtks))}")
    print(f"  cells            {whole['cells']}, of which {inner['cells']} are interior")
    print()
    print("  velocity against u = u_max (1 - y^2)")
    print(f"    interior L2        {inner['l2']:.6e}")
    print(f"    interior relative  {inner['relative_l2']:.6e}   tolerance {args.tol_velocity:.1e}")
    print(f"    whole field L2     {whole['l2']:.6e}   (not gated, see below)")
    print()
    print("  flow rate, which incompressibility makes constant along the channel")
    print(f"    interior mean      {q['mean']:.6f}   exact {channel.exact_flow_rate():.6f}")
    print(f"    interior error     {q['relative_error']:.6e}   tolerance {args.tol_flow:.1e}")
    print(f"    interior spread    {q['relative_spread']:.6e}   tolerance {args.tol_spread:.1e}")
    print(f"    spread with edges  {q_all['relative_spread']:.6e}   (not gated)")
    print()
    print("  pressure gradient, from the momentum balance")
    print(f"    measured           {dp['measured']:.6f}   exact {dp['exact']:.6f}")
    print(f"    relative error     {dp['relative_error']:.6e}   tolerance {args.tol_dpdx:.1e}")

    failures = []
    if inner["relative_l2"] > args.tol_velocity:
        failures.append(f"interior velocity relative L2 {inner['relative_l2']:.3e} > {args.tol_velocity:.1e}")
    if q["relative_error"] > args.tol_flow:
        failures.append(f"interior flow rate error {q['relative_error']:.3e} > {args.tol_flow:.1e}")
    if q["relative_spread"] > args.tol_spread:
        failures.append(f"interior flow rate spread {q['relative_spread']:.3e} > {args.tol_spread:.1e}")
    if dp["relative_error"] > args.tol_dpdx:
        failures.append(f"pressure gradient relative error {dp['relative_error']:.3e} > {args.tol_dpdx:.1e}")

    print()
    if failures:
        print("  FAIL")
        for f in failures:
            print(f"    {f}")
        return 1

    print("  PASS")
    if whole["l2"] > inner["l2"] * 5:
        b = channel.boundary_corner_error(grid)
        print()
        print("  Note, not a failure: the whole-field L2 above is "
              f"{whole['l2'] / inner['l2']:.0f} times the interior")
        print(f"  one, and the flow rate spread reaches {q_all['relative_spread']:.1%} once the column against")
        print("  the inlet is counted. That column is an artefact of how the VTK file is")
        print("  written, not of the solution:")
        print()
        print(f"    corners at x = {b['x_one_cell_in']:.4f}, one cell in    "
              f"error {b['error_one_cell_in']:.2e}")
        print(f"    corners at x = {b['x_on_boundary']:.4f}, on the inlet   "
              f"error {b['error_on_boundary']:.2e}")
        print()
        print("  Both belong to the same cells and come from the same facet values. Only")
        print("  the pair lying on the boundary plane is wrong, and it reaches "
              f"{b['worst_value']:.2f}")
        print(f"  against a u_max of {b['u_max']:.1f}, which no solution of this problem attains.")
        print("  The checks above are gated on the interior for that reason.")
    return 0


def cmd_order(args):
    runs = []
    for spec in args.run:
        if "=" not in spec:
            raise SystemExit(f"expected NX=DIR, got {spec!r}")
        nx, directory = spec.split("=", 1)
        grid = read_vtk(latest_vtk(directory))
        e = channel.field_error(grid, skip_edge_columns=args.skip_edges)
        h = channel.LENGTH / int(nx)
        runs.append({"nx": int(nx), "h": h, **e})

    runs.sort(key=lambda r: -r["h"])
    res = channel.observed_order([r["h"] for r in runs], [r["l2"] for r in runs])

    print("     nx        h         cells        L2          order vs coarser")
    for i, r in enumerate(runs):
        p = res["orders"][i - 1] if i > 0 and res["orders"][i - 1] is not None else None
        ptxt = f"{p:.2f}" if p is not None else ""
        print(f"  {r['nx']:5d}  {r['h']:.5f}  {r['cells']:8d}  {r['l2']:.6e}   {ptxt:>8}")

    good = [p for p in res["orders"] if p is not None]
    if not good:
        print("\n  no order could be computed")
        return 1
    mean_p = float(np.mean(good))
    print(f"\n  mean observed order: {mean_p:.2f}")
    print(f"  formal order of the scheme in use: {args.expect:.1f}")

    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"runs": runs, "orders": res["orders"], "mean": mean_p}, fh, indent=2)

    lo, hi = args.expect - args.order_tol, args.expect + args.order_tol
    if not (lo <= mean_p <= hi):
        print(f"\n  FAIL: observed order outside [{lo:.1f}, {hi:.1f}]")
        return 1
    print("\n  PASS")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    t = sub.add_parser("transient", help="error at every frame, to test for steady state")
    t.add_argument("--vtks", required=True)
    t.add_argument("--every", type=int, default=1, help="sample every Nth frame")
    t.add_argument("--steady-tol", type=float, default=1e-3,
                   help="relative spread over the last three frames that counts as steady")
    t.add_argument("--json")
    t.set_defaults(func=cmd_transient)

    c = sub.add_parser("check", help="the checks on one run, with tolerances")
    c.add_argument("--vtks", required=True)
    c.add_argument("--skip-edges", type=int, default=1,
                   help="cell columns to drop from each end; the inlet column is anomalous")
    c.add_argument("--tol-velocity", type=float, default=5e-3)
    c.add_argument("--tol-flow", type=float, default=5e-3)
    c.add_argument("--tol-spread", type=float, default=1e-3)
    c.add_argument("--tol-dpdx", type=float, default=5e-2)
    c.set_defaults(func=cmd_check)

    o = sub.add_parser("order", help="convergence order across resolutions")
    o.add_argument("--run", action="append", required=True, metavar="NX=DIR",
                   help="repeat once per resolution")
    o.add_argument("--skip-edges", type=int, default=0,
                   help="drop this many cell columns from each end before measuring")
    o.add_argument("--expect", type=float, default=2.0)
    o.add_argument("--order-tol", type=float, default=0.4)
    o.add_argument("--json")
    o.set_defaults(func=cmd_order)

    args = ap.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
