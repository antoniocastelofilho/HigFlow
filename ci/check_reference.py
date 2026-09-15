#!/usr/bin/env python3
"""Check simulation output against stored reference metrics.

Extracts min/max/mean of velocity components (POINT_DATA) and
pressure (CELL_DATA) from legacy VTK files and compares against
a YAML reference file.

Usage:
  ./ci/check_reference.py vtk_dir/ ci/reference/mycasename.yaml
  ./ci/check_reference.py --generate vtk_dir/ -o ref.yaml
"""

import sys
import os
import glob
import re
import math
import yaml


def parse_vtk(filename):
    """Parse a legacy VTK file and return dict of field metrics.

    Returns:
        {"vel": {"u": {"min":, "max":, "mean":},
                 "v": {...}, "w": {...}},
         "p": {"min":, "max":, "mean":},
         "fracvol_sum": float or None}
    """
    data = {
        "vel": {"u": {"min": None, "max": None, "sum": 0.0, "count": 0},
                "v": {"min": None, "max": None, "sum": 0.0, "count": 0},
                "w": {"min": None, "max": None, "sum": 0.0, "count": 0}},
        "p": {"min": None, "max": None, "sum": 0.0, "count": 0},
        "fracvol_sum": None
    }

    with open(filename) as fp:
        lines = fp.readlines()

    i = 0
    nlines = len(lines)

    # Determine num points / num cells from POINTS and CELL_DATA headers
    num_cells = 0
    while i < nlines:
        line = lines[i]
        if line.startswith("CELL_DATA"):
            num_cells = int(line.split()[1])
            break
        i += 1

    # Parse fields sequentially
    i = 0
    section = None
    field_name = None
    field_type = None
    expecting_lookup = False

    while i < nlines:
        line = lines[i]

        if line.startswith("SCALARS "):
            parts = line.split()
            field_name = parts[1]
            field_type = parts[2] if len(parts) > 2 else "FLOAT"
            section = "scalar"
            expecting_lookup = True
            i += 1
            continue

        if expecting_lookup and line.startswith("LOOKUP_TABLE"):
            expecting_lookup = False
            i += 1
            # Read scalar data
            count = 0
            if field_name == "p":
                data_p = data["p"]
            while i < nlines and count < num_cells:
                val_line = lines[i].strip()
                if not val_line or val_line.startswith(("SCALARS ", "VECTORS ",
                                                        "TENSORS ", "POINT_DATA ",
                                                        "CELL_DATA ", "LOOKUP_TABLE")):
                    break
                try:
                    val = float(val_line)
                except ValueError:
                    break
                if field_name == "p":
                    if data_p["min"] is None or val < data_p["min"]:
                        data_p["min"] = val
                    if data_p["max"] is None or val > data_p["max"]:
                        data_p["max"] = val
                    data_p["sum"] += val
                    data_p["count"] += 1
                elif field_name == "FracVol":
                    if data["fracvol_sum"] is None:
                        data["fracvol_sum"] = 0.0
                    data["fracvol_sum"] += val
                count += 1
                i += 1
            continue

        if line.startswith("VECTORS "):
            parts = line.split()
            field_name = parts[1]
            field_type = parts[2] if len(parts) > 2 else "FLOAT"
            section = "vectors"
            i += 1
            continue

        if section == "vectors":
            # vel vectors: read until next header or EOF
            if line.startswith(("SCALARS ", "VECTORS ", "TENSORS ",
                                "POINT_DATA ", "CELL_DATA ", "LOOKUP_TABLE")):
                section = None
                continue

            if field_name == "vel":
                vals = line.strip().split()
                if len(vals) >= 3 and field_name == "vel":
                    for k, comp in enumerate(["u", "v", "w"]):
                        try:
                            v = float(vals[k])
                        except ValueError:
                            continue
                        d = data["vel"][comp]
                        if d["min"] is None or v < d["min"]:
                            d["min"] = v
                        if d["max"] is None or v > d["max"]:
                            d["max"] = v
                        d["sum"] += v
                        d["count"] += 1
            i += 1
            continue

        i += 1

    # Convert sums to means
    for comp in ["u", "v", "w"]:
        d = data["vel"][comp]
        if d["count"] > 0:
            d["mean"] = d["sum"] / d["count"]
        del d["sum"]
        del d["count"]

    if data["p"]["count"] > 0:
        data["p"]["mean"] = data["p"]["sum"] / data["p"]["count"]
    del data["p"]["sum"]
    del data["p"]["count"]

    return data


def metrics_close(a, b, tol=0.001):
    """Return True if a and b are within relative tolerance."""
    if a is None and b is None:
        return True
    if a is None or b is None:
        return False
    denom = max(abs(a), abs(b), 1e-16)
    return abs(a - b) / denom < tol


def compare_data(current, reference, case_name, tol=0.001):
    """Compare current metrics against reference, return (pass, report_lines)."""
    report = []
    passed = True

    # Compare vel components
    for comp in ["u", "v", "w"]:
        for metric in ["min", "max", "mean"]:
            v_cur = current.get("vel", {}).get(comp, {}).get(metric)
            v_ref = reference.get("vel", {}).get(comp, {}).get(metric)
            if not metrics_close(v_cur, v_ref, tol):
                passed = False
                report.append(
                    f"  vel.{comp}.{metric}: current={v_cur}, ref={v_ref}")

    # Compare pressure
    for metric in ["min", "max", "mean"]:
        v_cur = current.get("p", {}).get(metric)
        v_ref = reference.get("p", {}).get(metric)
        if not metrics_close(v_cur, v_ref, tol):
            passed = False
            report.append(f"  p.{metric}: current={v_cur}, ref={v_ref}")

    # Compare FracVol sum if present
    fc_cur = current.get("fracvol_sum")
    fc_ref = reference.get("fracvol_sum")
    if fc_cur is not None and fc_ref is not None:
        if not metrics_close(fc_cur, fc_ref, tol):
            passed = False
            report.append(
                f"  fracvol_sum: current={fc_cur}, ref={fc_ref}")

    return passed, report


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Check VTK output against reference")
    parser.add_argument("vtk_dir", help="VTK directory")
    parser.add_argument("reference", nargs="?", help="Reference YAML file")
    parser.add_argument("--generate", action="store_true",
                        help="Generate reference YAML from current output")
    parser.add_argument("-o", "--output",
                        help="Output file for --generate (default: stdout)")
    parser.add_argument("--tolerance", type=float, default=0.001,
                        help="Relative tolerance (default 0.001)")
    args = parser.parse_args()

    vtk_files = sorted(glob.glob(os.path.join(args.vtk_dir, "*.vtk")))
    if not vtk_files:
        print(f"No VTK files found in {args.vtk_dir}")
        sys.exit(1)

    # Parse all VTK files and average per-step metrics
    all_steps = {}
    for f in vtk_files:
        step_match = re.search(r'print[_-]\d+-(\d+)', os.path.basename(f))
        step = int(step_match.group(1)) if step_match else 0
        if step not in all_steps:
            all_steps[step] = []
        all_steps[step].append(parse_vtk(f))

    # Average across ranks for each step, then take last step for reference
    last_step = max(all_steps.keys())
    step_data_list = all_steps[last_step]

    # Aggregate across ranks
    current = {
        "vel": {"u": [], "v": [], "w": []},
        "p": [],
        "fracvol_sum": 0.0
    }
    for d in step_data_list:
        for comp in ["u", "v", "w"]:
            current["vel"][comp].append(d["vel"][comp])
        current["p"].append(d["p"])
        if d["fracvol_sum"] is not None:
            current["fracvol_sum"] += d["fracvol_sum"]
        else:
            current["fracvol_sum"] = None

    # Compute global min/max/mean from per-rank data
    def agg_min(arr, key):
        return min(d[key] for d in arr if d[key] is not None)
    def agg_max(arr, key):
        return max(d[key] for d in arr if d[key] is not None)
    def agg_mean(arr, key):
        vals = [d[key] for d in arr if d[key] is not None]
        if not vals:
            return None
        return sum(vals) / len(vals)

    result = {}
    result["vel"] = {}
    for comp in ["u", "v", "w"]:
        result["vel"][comp] = {
            "min": agg_min(current["vel"][comp], "min"),
            "max": agg_max(current["vel"][comp], "max"),
            "mean": agg_mean(current["vel"][comp], "mean"),
        }
    result["p"] = {
        "min": agg_min(current["p"], "min"),
        "max": agg_max(current["p"], "max"),
        "mean": agg_mean(current["p"], "mean"),
    }
    if current["fracvol_sum"] is not None:
        result["fracvol_sum"] = current["fracvol_sum"]

    if args.generate:
        output = args.output
        if not output:
            print(yaml.dump({"last_step": last_step, "metrics": result},
                            default_flow_style=False, sort_keys=False))
        else:
            with open(output, "w") as fp:
                yaml.dump({"last_step": last_step, "metrics": result},
                          fp, default_flow_style=False, sort_keys=False)
            print(f"Reference written to {output}")
        sys.exit(0)

    # Comparison mode
    if not args.reference:
        print("Either specify reference file or use --generate")
        sys.exit(1)

    with open(args.reference) as fp:
        ref_data = yaml.safe_load(fp)

    ref_metrics = ref_data.get("metrics", {})
    case_name = os.path.splitext(os.path.basename(args.reference))[0]

    print(f"Comparing {case_name} against {os.path.basename(args.reference)}")
    print(f"  Last step: current={last_step}, ref={ref_data.get('last_step')}")
    print(f"  Velocity:")
    for comp in ["u", "v", "w"]:
        print(f"    {comp}: min={result['vel'][comp]['min']:.6e}, "
              f"max={result['vel'][comp]['max']:.6e}, "
              f"mean={result['vel'][comp]['mean']:.6e}")
    print(f"  Pressure:")
    print(f"    min={result['p']['min']:.6e}, "
          f"max={result['p']['max']:.6e}, "
          f"mean={result['p']['mean']:.6e}")
    if result.get("fracvol_sum") is not None:
        print(f"  FracVol sum={result['fracvol_sum']:.6f}")

    passed, diffs = compare_data(result, ref_metrics, case_name, args.tolerance)
    if diffs:
        for d in diffs:
            print(f"  DIFF: {d}")

    if passed:
        print(f"PASS: all metrics within {args.tolerance*100:.2f}% tolerance")
        sys.exit(0)
    else:
        print(f"FAIL: metrics exceed {args.tolerance*100:.2f}% tolerance")
        sys.exit(1)


if __name__ == "__main__":
    main()
