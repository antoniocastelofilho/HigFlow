# Verification

Checks the solver against solutions that are known in closed form, so that
"it runs" and "it is right" stop being the same claim.

```bash
docker build -f containers/Dockerfile -t higflow:latest .
bash tests/run_verification.sh
```

Exits non-zero when a check fails.

## What is checked

The case is `higflow/example2d_Newt`: a channel of length 8 and half-height 1,
driven by a parabolic inlet profile, at Re = 1. At steady state its solution is
known exactly:

```
u(y) = u_max (1 - y^2)        u_max = 1.5
dp/dx = -2 mu u_max / h^2 = -3
Q     = 4 u_max h / 3 = 2
```

Three independent things follow from that, and all three are checked.

| Check | Measured at 160 x 40 | Tolerance |
|---|---|---|
| Velocity against the exact profile, relative L2 | 6.85e-04 | 5e-03 |
| Flow rate against the exact 2.0 | 6.47e-04 | 5e-03 |
| Flow rate spread along the channel | 1.46e-05 | 1e-03 |
| Pressure gradient against -3 | 9.65e-04 | 5e-02 |

They are independent in a useful way: the velocity check compares a field to a
formula, the flow rate is a conservation law the discretisation should satisfy
whatever the profile looks like, and the pressure gradient comes from the
momentum balance rather than from the boundary condition.

## Steady state comes first

```bash
python3 tests/verification/verify.py transient --vtks <VTKS dir> --every 10
```

The other checks mean nothing without this one. If a run has not settled, its
error is dominated by the transient, refining the mesh will not reduce it, and
a convergence study would measure the time integration rather than the spatial
scheme.

At the shipped settings the channel does settle. The error falls by 93 % over
the first ten frames and then stops moving:

```
   frame        L2          change vs previous
       0   1.082657e+00
      10   6.919085e-02      -93.61%
      20   5.339843e-02      -22.82%
      30   5.253928e-02       -1.61%
      40   5.245338e-02       -0.16%
      ...
     100   5.243904e-02       -0.00%
```

## Convergence order

```bash
python3 tests/verification/verify.py order \
    --run 40=<dir> --run 80=<dir> --run 160=<dir> --run 320=<dir> --skip-edges 1
```

`tests/verification/channel_mesh.py` writes the meshes, so the same case can be
run at any resolution. At 160 x 40 it reproduces the meshes that ship with the
example byte for byte.

Measured over four meshes, each a halving of the last:

```
   nx        h        cells        L2         order vs coarser
   40    0.20000       380    1.678394e-02
   80    0.10000      1560    4.140728e-03       2.02
  160    0.05000      6320    1.028249e-03       2.01
  320    0.02500     25440    2.562075e-04       2.00
```

**Observed order 2.01 against a formal order of 2.** Three consecutive
refinements agree to two decimal places, which is about as clean as this
measurement gets.

## The inlet column

The numbers above exclude one cell column at each end. That is not a
convenience, and it is worth stating plainly.

**The cell column against the inlet does not carry the interior solution.** It
reports a nearly flat profile close to u_max rather than the parabola:

```
  80 x 20, first column at x = 0.05
    y = -0.95   u = 0.676   exact 0.146
    y = -0.85   u = 1.415   exact 0.416
    y = -0.75   u = 1.479   exact 0.656
    y = -0.65   u = 1.479   exact 0.866
```

Two consequences follow.

Its flow rate is wrong by a wide margin, so the spread along the channel reads
42.7 % when that column is counted and 1.5e-05 when it is not.

And it dominates the norm badly enough to invert the convergence study. Over
the whole field the measured order is 0.19 and the error rises between the two
coarsest meshes, because the anomaly grows with refinement while everything
else shrinks:

```
   nx      whole field L2    order       interior L2      order
   40      5.243904e-02                  1.678394e-02
   80      5.747045e-02      -0.13       4.140728e-03      2.02
  160      4.737065e-02       0.28       1.028249e-03      2.01
  320      3.561326e-02       0.41       2.562075e-04      2.00
```

Read only the left-hand columns and the scheme looks first order at best. Read
the right-hand ones and it is exactly second order. The difference is one
column of cells out of 320.

What has not been established is whether the solver computes a wrong value
there or writes a wrong one to the VTK file. Distinguishing those needs a look
at the field before it is written, which the checks here do not do. Both checks
and convergence study are therefore gated on the interior, and the whole-field
number is reported next to it rather than hidden.

## Files

| | |
|---|---|
| `run_verification.sh` | runs the case and checks it; the entry point for CI |
| `verification/channel.py` | the exact solution and the error measures |
| `verification/verify.py` | `transient`, `check` and `order` |
| `verification/channel_mesh.py` | writes channel meshes at any resolution |

The VTK reader is `tools/higflow_vtk.py`, shared with the gallery renderer so
the two cannot disagree about what a field means.

## What is not here yet

**The method of manufactured solutions.** Every check here uses a case whose
exact solution happens to be known. MMS would let any model be verified, by
choosing a solution, substituting it into the equations and taking the residual
as a source term. It needs the driver to accept a source term, so it is a
larger change than this.

**Anything but the Newtonian channel.** The viscoelastic, multiphase and
electro-osmotic solvers have no check here. Two of the eleven shipped cases do
not run at all, which is the more pressing problem.
