# example3d_RisingDropCap — 3D rising bubble, cap regime (uniform)

**Problem:** A light bubble (ρ = 1 kg/m³, μ = 0.1 Pa·s) rises through
a denser ambient fluid (ρ = 1000 kg/m³, μ = 10 Pa·s) in a
[0,1] × [0,2] × [0,1] domain.  Initial bubble: sphere of radius 0.25
centred at (0.5, 0.5, 0.5).  All walls: no-slip.

**⚠️ Gravity direction:** The HigFlow multiphase solver applies gravity
on `dim == 1` (Y axis).  The domain is oriented with **Y as the tall
dimension** (height = 2.0) so the bubble rises from y = 0.5 toward
y = 2.0, matching the OpenFOAM reference where gravity acts in Z.

**Reference:** Silva et al. (2023), Mathematics 11, 3900, Fig. 9 (case 2).
Comparison data:
`~/foam/pedro-5.0/run/risingDrop_ellipsoidal/fig9_caso1_dados.csv`

## Running

```bash
cd higflow/example3d_RisingDropCap
make MESH=bubble_50x50x100 build_run NP=1 \
  IN="mult newt-newt Re=35 Ca=3.57 dt=0.001 numsteps=30000 dts=0.1 dtp=0.05" \
  ENAME=$(date +"%Y%m%d_%H%M")
```

## Mesh

| Direction | Range | Cells | Size |
| ----------- | :----- | :-----: | :----: |
| X | [0, 1] | 20 | 0.05 |
| Y (tall) | [0, 2] | 40 | 0.05 |
| Z | [0, 1] | 20 | 0.05 |

Total: 16 000 cells.  AMR (up to 2 levels) available by toggling
`ADAPT_ENABLED` in `ns-exemple-3d.h`.

## Key parameters

| Parameter | Value | Notes |
| ----------- | ------- | ------- |
| Re | 35 | U = √(g·D) |
| Ca | 3.57 | σ = 1.96 |
| Fr | 1.0 | |
| Eötvös | ~125 | Cap regime |
| Morton | ~1.3 | |
| ρ₀/ρ₁ | 1.0 / 0.001 | ambient / drop |
| μ₀/μ₁ | 1.0 / 0.01 | ambient / drop |
| Bubble | centre (0.5, 0.5, 0.5), r = 0.25 | |

Phase indexing follows the solver mixing rule `x = (1−FracVol)·x₀ + FracVol·x₁`:
the drop is the sphere where `FracVol = 1`, so it is **phase 1** (light), and
the surrounding ambient is **phase 0** (heavy).  The light drop rises in +Y.

## Solver

- Uniform: `-ksp_type fbcgsr -pc_type hypre -pc_hypre_type boomeramg`
- AMR: `-ksp_type gmres -ksp_gmres_restart 200 -pc_type hypre -pc_hypre_type boomeramg`
