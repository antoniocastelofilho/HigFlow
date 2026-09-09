# example2d_DamBreak_AMR — 2D Dam break with adaptive mesh refinement

**Problem:** A rectangular water column (1.2 m × 0.6 m) collapses under
gravity inside a closed domain (3.22 m × 2.0 m). The right, left, and
bottom walls are no-slip; the top wall is open (Neumann).

**Physics:** Newtonian two-phase (water + air), VOF, surface tension,
gravity.  Water: ρ = 1000 kg/m³, ν = 1.0e-6 m²/s.  Air: ρ = 1.225 kg/m³,
ν = 1.48e-5 m²/s.

**Reference:** OpenFOAM tutorial `dambreakValidation` (foam-extend 5.0);
experimental data from Zhainakov & Kurbanaliev (2013).

## Running

```bash
cd higflow/example2d_DamBreak_AMR
make MESH=dambreak_160x100 build_run NP=1 \
  IN="mult newt-newt Re=3130000 Ca=0.045 dt=0.001 numsteps=8000 dts=0.05 dtp=0.005" \
  ENAME=$(date +"%Y%m%d_%H%M")
```

Output goes to `output/OUTNAME__MESHNAME__ENAME/{vtk,save,res}/`.

| Parameter | Value | Notes |
|-----------|-------|-------|
| Re        | 3.13e6 | Based on U = √(g·H), water viscosity |
| Ca        | 0.045  | σ = 0.07 N/m |
| Fr        | 1.0    | U / √(g·H) |
| ρ₀/ρ₁    | 815:1 | water / air |
| μ₀/μ₁    | 0.018 | air dynamic-viscosity ratio |

## Mesh

Base mesh: 160 × 100 cells (dx = dy = 0.02).  AMR adds up to 2
refinement levels near the interface (cells ¼ the base size).

## Output

- VTK files (ParaView): `output/*/vtk/ns.print_*.vtk`
- Restart checkpoints: `output/*/save/`
- Simulation log: `output/*/save/ns.log`
