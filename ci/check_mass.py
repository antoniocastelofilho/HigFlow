#!/usr/bin/env python3
"""Check mass conservation from VTK output files.

Parses SCALARS FracVol from CELL_DATA section of legacy VTK files
and verifies the sum is constant across timesteps (within tolerance).
"""

import sys
import os
import glob
import re


def extract_fracvol(filename):
    """Return sum of FracVol values from a VTK file, or None if absent."""
    found = False
    in_fracvol = False
    total = 0.0
    count = 0
    with open(filename) as fp:
        for line in fp:
            if line.startswith("SCALARS FracVol"):
                in_fracvol = True
                found = True
                next(fp, None)  # skip LOOKUP_TABLE line
                continue
            if in_fracvol:
                if line.startswith(("SCALARS ", "POINT_DATA ", "CELL_DATA ")):
                    break
                try:
                    val = float(line.strip())
                    total += val
                    count += 1
                except ValueError:
                    break
    if not found:
        return None
    return total, count


def parse_step(filename):
    """Extract step number from filename like ns.print_0-3.vtk"""
    m = re.search(r'print[_-]\d+-(\d+)', filename)
    return int(m.group(1)) if m else 0


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Check FracVol mass conservation in VTK files")
    parser.add_argument("vtk_dir", help="Directory containing VTK files")
    parser.add_argument("--tolerance", type=float, default=0.001,
                        help="Max allowed fractional change (default 0.001 = 0.1%%)")
    parser.add_argument("--single-phase", action="store_true",
                        help="Single-phase case (skip mass check)")
    args = parser.parse_args()

    vtk_files = sorted(glob.glob(os.path.join(args.vtk_dir, "*.vtk")),
                       key=parse_step)

    if not vtk_files:
        print(f"No VTK files found in {args.vtk_dir}")
        sys.exit(1)

    if args.single_phase:
        print(f"Single-phase case, mass conservation N/A")
        sys.exit(0)

    fracvol_by_step = {}
    for f in vtk_files:
        step = parse_step(f)
        result = extract_fracvol(f)
        if result is None:
            print(f"{os.path.basename(f)}: no FracVol field")
            continue
        total, count = result
        fracvol_by_step[step] = total
        print(f"{os.path.basename(f)}: sum_fracvol={total:.6f}, cells={count}")

    if not fracvol_by_step:
        print("No FracVol data found in any VTK file, skipping")
        sys.exit(0)

    steps_sorted = sorted(fracvol_by_step.keys())
    ref_sum = fracvol_by_step[steps_sorted[0]]
    max_change = 0.0

    for step in steps_sorted[1:]:
        s = fracvol_by_step[step]
        change = abs(s - ref_sum) / ref_sum if ref_sum != 0 else abs(s - ref_sum)
        max_change = max(max_change, change)
        print(f"  step {step}: sum={s:.6f}, change={change:.4e} ({change*100:.4f}%)")

    print(f"\nSummary: ref_sum={ref_sum:.6f}, max_change={max_change:.4e}")
    if max_change > args.tolerance:
        print(f"FAIL: max change {max_change*100:.4f}% > tolerance {args.tolerance*100:.2f}%")
        sys.exit(1)
    print(f"PASS: mass conserved within {args.tolerance*100:.2f}%")
    sys.exit(0)


if __name__ == "__main__":
    main()
