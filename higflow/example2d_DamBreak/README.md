# example2d_DamBreak — 2D Dam break (uniform mesh)

**Problem:** A rectangular water column (1.2 m × 0.6 m) collapses under
gravity inside a closed domain (3.22 m × 2.0 m).  Right, left, and
bottom walls are no-slip; the top wall is open (Neumann).  Gravity:
(0, -9.81) m/s².

**Physics:** Newtonian two-phase (water + air), VOF with piecewise-linear
interface reconstruction (PLIC), surface tension, gravity.

| Property | Water (phase 0) | Air (phase 1) |
|----------|:---------------:|:-------------:|
| ρ (kg/m³) | 998.2 | 1.225 |
| μ (Pa·s)  | 1.0e-3 | 1.81e-5 |
| ν (m²/s)  | 1.0e-6 | 1.48e-5 |
| σ (N/m)   | 0.07 | — |

**Characteristic scales:** H = 1.0 m, U = √(g·H) ≈ 3.13 m/s.
**Dimensionless numbers:** Re ≈ 3.13×10⁶, Ca = 0.045, Fr = 1.0.

**Reference:** OpenFOAM tutorial `dambreakValidation` (foam-extend 5.0);
experimental data in Zhainakov & Kurbanaliev (2013).

## Running

```bash
cd higflow/example2d_DamBreak
make MESH=dambreak_160x100 build_run NP=1 \
  IN="mult newt-newt Re=3130000 Ca=0.045 dt=0.001 numsteps=8000 dts=0.05 dtp=0.005" \
  ENAME=$(date +"%Y%m%d_%H%M")
```

Output: `output/OUTNAME__MESHNAME__ENAME/{vtk,save,res}/`.

## Mesh

160 × 100 cells (dx = 0.020125, dy = 0.02) over [0, 3.22] × [0, 2.0].

## Water column

Rectangle [0, 1.2] × [0, 0.602] (0.602 avoids alignment with cell faces
on the 160×100 mesh).
