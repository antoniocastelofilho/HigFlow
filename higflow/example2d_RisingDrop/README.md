# example2d_RisingDrop — 2D rising drop, cap regime (uniform)

2D analogue of `example3d_RisingDropCap`, built on the mature 2D multiphase VOF
(the same machinery as `example2d_DamBreak`).  A light circular drop rises through
a denser ambient fluid in a `[0,1] × [0,2]` box.  Initial drop: circle of radius
0.25 centred at `(0.5, 0.5)`.  All four walls are no-slip.

**Gravity direction:** the solver applies gravity on `dim == 1` (Y), so the domain
is oriented with **Y as the tall dimension** (height 2.0); the drop rises in +Y.
This matches the OpenFOAM reference `~/foam/pedro-5.0/run/risingDrop_cap`
(gravity in Z, mesh 50×50×100), nondimensionalised by the ambient fluid.

## Phase indexing

Mixing rule is `x = (1 − FracVol)·x₀ + FracVol·x₁`, and `get_fracvol` returns 1
inside the circle (`func() ≤ 0`).  So the **drop is phase 1** (light:
rho 0.001, mu 0.01) and the **ambient is phase 0** (heavy: rho 1.0, mu 1.0).

## Key parameters

| Parameter | Value | Notes |
| --- | --- | --- |
| Re | 35 | U = √(g·D), D = 0.5 |
| Ca | 3.57 | σ = 1.96 |
| Fr | 1.0 | |
| ρ₀ / ρ₁ | 1.0 / 0.001 | ambient / drop |
| μ₀ / μ₁ | 1.0 / 0.01 | ambient / drop |
| Drop | centre (0.5, 0.5), r = 0.25 | |

## Meshes

| Name | Cells | Δ | drop radius |
| --- | --- | --- | --- |
| `drop_32x64` | 32 × 64 (2 048) | 0.03125 | ~8 cells |
| `drop_64x128` | 64 × 128 (8 192) | 0.015625 | ~16 cells |

`.amr` cell sizes are written as exact `(h−l)/patchsize` (required within 1e-12 by
`higio_read_from_amr_info`; rounding them corrupts the run).

## Running

```bash
cd higflow/example2d_RisingDrop
make build_run MESH=drop_64x128 NP=4 \
  IN="mult newt-newt Re=35 Ca=3.57 dt=0.001 numsteps=3000 dts=0.5 dtp=0.05" \
  ENAME=run1
```

`numsteps=3000`, `dt=0.001` → t = 3.0 s (same physical time as the OpenFOAM case);
`dtp=0.05` writes 60 VTK frames.  Output in `output/<OUTNAME>/{vtk,save,res}`.

## Why a 2D version

The 3D VOF (`example3d_RisingDropCap`) produces a physically correct, upward
velocity field but the interface does not translate with it (FracVol erodes at the
top instead of the drop rising) — a 3D-VOF advection issue still under
investigation.  The 2D VOF is mature (and has an adaptive-mesh variant), so this
case isolates the physics on the trusted 2D path.
