# example3d_RisingDropCap_AMR — 3D rising bubble, cap regime, AMR

**Problem:** A light bubble (ρ = 1 kg/m³, μ = 0.1 Pa·s) rises through
a denser ambient fluid (ρ = 1000 kg/m³, μ = 10 Pa·s) in a
[0,1] × [0,2] × [0,1] domain.  Initial bubble: sphere of radius 0.25
centred at (0.5, 0.5, 0.5).

**Gravity direction:** Y axis (`dim == 1`), because the HigFlow
multiphase solver hardcodes the gravity term on `dim == 1`.  The
domain is oriented with Y as the tall dimension (height = 2.0) so
the bubble rises from y = 0.5 toward y = 2.0, matching the
OpenFOAM reference where gravity acts in Z on a [0,1] × [0,1] × [0,2]
domain.

**Reference:** Silva et al. (2023), Mathematics 11, 3900, Fig. 9 (case 2).

## Running

```bash
cd higflow/example3d_RisingDropCap_AMR
make MESH=bubble_50x50x100 build_run NP=1 \
  IN="mult newt-newt Re=35 Ca=3.57 dt=0.001 numsteps=30000 dts=0.1 dtp=0.05" \
  ENAME=$(date +"%Y%m%d_%H%M")
```

| Parameter | Value | Notes |
|-----------|-------|-------|
| Re        | 35    | Based on U = √(g·D) |
| Ca        | 3.57  | σ = 1.96 |
| Fr        | 1.0   | U / √(g·D) |
| Eötvös    | ~125  | Cap regime (Eo > 40) |
| Morton    | ~1.3  | |
| ρ₀/ρ₁    | 0.001 / 1.0 | bubble / ambient |
| μ₀/μ₁    | 0.01 / 1.0 | bubble / ambient |
| Mesh      | 20 × 40 × 20 | Y is tall = gravity dir |
| AMR       | Up to 2 levels | `ADAPT_ENABLED 1` in header |
| Solver    | `gmres(200) + hypre boomeramg` | 3D-optimised |

## Mesh setup (first time)

```bash
mkdir -p mesh/__using/{domain,bc}
cp mesh/bubble_50x50x100/domain/bubble_50x50x100-d.amr mesh/__using/domain/ch-d.amr
for i in 0 1 2 3 4 5; do
  cp "mesh/bubble_50x50x100/bc/bubble_50x50x100-bc-$i.amr" "mesh/__using/bc/ch-bc-$i.amr"
done
```

## Output

- VTK (initial adapted mesh + time steps): `output/*/vtk/`
- Restart checkpoints: `output/*/save/`
- Plot z_c(t) and v_c(t) with reference data: use
  `../../foam/pedro-5.0/run/risingDrop_ellipsoidal/postprocess.py`
