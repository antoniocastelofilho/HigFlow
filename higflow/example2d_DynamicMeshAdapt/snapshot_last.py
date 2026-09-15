"""
Render the last VTK timestep of the most recent simulation as a PNG.

Shows the FracVol field with Surface With Edges and the Blue-Green-Orange
colormap, matching the interactive ParaView view.

Usage:
    pvpython snapshot_last.py [output_dir] [--step N] [--size WxH]

Arguments:
    output_dir   Root output directory to search (default: ./output/)
    --step N     Render timestep N instead of the last one
    --size WxH   Image resolution, e.g. --size 1600x1600 (default: 1200x1200)

The PNG is written in the same directory as this script:
    <sim_name>__step<N>.fracvol.png

Note: ParaView 6.1 exits with code 134 after a successful save — this is a
known upstream bug and can be ignored.
"""

import os
import sys
import re
import glob
import subprocess

_SENTINEL = "--_offscreen_done"
_OFFSCREEN_FLAG = "--force-offscreen-rendering"


def _reexec_offscreen():
    """Re-launch under pvpython --force-offscreen-rendering (once)."""
    if _SENTINEL not in sys.argv:
        clean = [a for a in sys.argv if a != _SENTINEL]
        sys.exit(subprocess.call(
            ["pvpython", _OFFSCREEN_FLAG] + clean + [_SENTINEL]
        ))
    sys.argv = [a for a in sys.argv if a != _SENTINEL]


def _parse_args():
    """Return (output_root, target_step, img_w, img_h) from sys.argv."""
    tokens = sys.argv[1:]
    output_root = None
    target_step = None
    img_w, img_h = 1200, 1200

    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok == "--step":
            i += 1
            target_step = int(tokens[i])
        elif tok == "--size":
            i += 1
            w, h = tokens[i].split("x")
            img_w, img_h = int(w), int(h)
        elif not tok.startswith("--"):
            output_root = tok
        i += 1

    return output_root, target_step, img_w, img_h


def _dir_mtime(path):
    """Return the mtime of the newest file inside a directory."""
    try:
        files = [os.path.join(path, f) for f in os.listdir(path)]
        return max(os.path.getmtime(f) for f in files) if files else 0.0
    except OSError:
        return 0.0


def _find_vtk(output_root, target_step):
    """Return (sim_name, vtk_dir, step, filename) for the requested timestep."""
    vtk_dirs = glob.glob(os.path.join(output_root, "*", "vtk"))
    if not vtk_dirs:
        sys.exit("ERROR: no vtk/ directories found under " + output_root)

    vtk_dir = max(vtk_dirs, key=_dir_mtime)
    sim_name = os.path.basename(os.path.dirname(vtk_dir))
    print(f"Simulation : {sim_name}")

    pattern = re.compile(r"^(.+)\.print_0-(\d+)\.vtk$")
    candidates = []
    for fname in os.listdir(vtk_dir):
        match = pattern.match(fname)
        if match:
            candidates.append((int(match.group(2)), fname))

    if not candidates:
        sys.exit("ERROR: no *.print_0-N.vtk files found in " + vtk_dir)

    candidates.sort()

    if target_step is not None:
        hits = [(n, f) for n, f in candidates if n == target_step]
        if not hits:
            available = [n for n, _ in candidates]
            sys.exit(f"ERROR: timestep {target_step} not found; available: {available}")
        step, fname = hits[0]
    else:
        step, fname = candidates[-1]

    return sim_name, vtk_dir, step, fname


def _world_to_viewport_x(world_x, focal_x, parallel_scale, aspect):
    """Map a world-space X coordinate to normalized viewport [0, 1]."""
    return 0.5 + (world_x - focal_x) / (2.0 * parallel_scale * aspect)


def _render(vtk_path, png_path, img_w, img_h):
    """Load VTK, set up FracVol display, and save a PNG screenshot."""
    # pylint: disable=import-outside-toplevel
    from paraview.simple import (
        LegacyVTKReader, GetActiveViewOrCreate, Show, Render, ResetCamera,
        ColorBy, GetColorTransferFunction, GetScalarBar, SaveScreenshot,
    )

    view = GetActiveViewOrCreate("RenderView")
    view.Background = [0.0, 0.0, 0.0]
    view.CameraParallelProjection = 1

    rw = view.GetRenderWindow()
    rw.SetOffScreenRendering(1)
    rw.SetSize(img_w, img_h)

    reader = LegacyVTKReader(FileNames=[vtk_path])
    display = Show(reader, view)
    display.Representation = "Surface With Edges"
    display.EdgeColor = [0.0, 0.0, 0.0]

    ColorBy(display, ("CELLS", "FracVol"))

    lut = GetColorTransferFunction("FracVol")
    lut.ApplyPreset("Blue - Green - Orange", True)
    lut.ScalarRangeInitialized = 1

    sb = GetScalarBar(lut, view)
    sb.Title = "FracVol"
    sb.ComponentTitle = ""
    sb.Visibility = 1
    sb.Orientation = "Horizontal"
    sb.TitleFontSize = 14
    sb.LabelFontSize = 11

    ResetCamera()
    view.CameraParallelScale = 0.59
    Render()

    # Align the colorbar with the geometry's x-extent in viewport coordinates.
    # The scalar bar Position/Length are in normalized viewport space [0, 1],
    # so we project the data bounds through the current camera parameters.
    focal = view.CameraFocalPoint
    scale = view.CameraParallelScale
    aspect = img_w / img_h

    bounds = reader.GetDataInformation().GetBounds()
    vp_left = _world_to_viewport_x(bounds[0], focal[0], scale, aspect)
    vp_right = _world_to_viewport_x(bounds[1], focal[0], scale, aspect)

    margin = 0.015
    sb.Position = [vp_left + margin, 0.02]
    sb.ScalarBarLength = max(0.0, vp_right - vp_left - 2.0 * margin)

    Render()
    SaveScreenshot(png_path, view)


def main():
    """Parse arguments, find the VTK file, and render the snapshot."""
    output_root, target_step, img_w, img_h = _parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    if output_root is None:
        output_root = os.path.join(script_dir, "output")

    if not os.path.isdir(output_root):
        sys.exit(f"ERROR: output directory not found: {output_root}")

    sim_name, vtk_dir, step, fname = _find_vtk(output_root, target_step)
    vtk_path = os.path.join(vtk_dir, fname)
    png_path = os.path.join(script_dir, f"{sim_name}__step{step}.fracvol.png")

    print(f"Main VTK   : {fname}  (timestep {step})")
    print(f"Output PNG : {png_path}")

    _render(vtk_path, png_path, img_w, img_h)
    print(f"Screenshot saved: {png_path}")


if __name__ == "__main__":
    _reexec_offscreen()
    main()
