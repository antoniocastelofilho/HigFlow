# Verification

Checks the solver against solutions that are known in closed form, so that
"it runs" and "it is right" stop being the same claim.

```bash
docker build -f containers/Dockerfile -t higflow:latest .
bash tests/run_verification.sh
```

Exits non-zero when a check fails.

## Through CTest

`tests/CMakeLists.txt` registers the same checks as CTest tests. It is a
standalone configuration, so it needs nothing from the project's own build:

```bash
cmake -S tests -B build-tests
ctest --test-dir build-tests --output-on-failure
```

Configured that way, one test runs: it builds the container, runs the case and
checks it. Point it at a run that already exists and it registers two faster
ones instead, and a third if the convergence meshes are also on disk:

```bash
cmake -S tests -B build-tests \
    -DHIGFLOW_VTKS=<VTKS dir> \
    -DHIGFLOW_ORDER_RUNS="40=<dir>;80=<dir>;160=<dir>;320=<dir>"
```

```
    Start 1: channel_steady_state
1/3 Test #1: channel_steady_state .............   Passed    1.20 sec
    Start 2: channel_exact_solution
2/3 Test #2: channel_exact_solution ...........   Passed    0.49 sec
    Start 3: channel_convergence_order
3/3 Test #3: channel_convergence_order ........   Passed    0.75 sec

100% tests passed, 0 tests failed out of 3
```

The slow ones carry the `slow` label, so `ctest -L verification -LE slow`
selects the checks that run in a couple of seconds.

The configuration is deliberately separate from the two build systems the
project already has, which compile different source sets. Choosing between them
is its own change; once there is one, `add_subdirectory(tests)` from the root
picks these up unchanged.

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

## What runs without the solver

`tests/verification/selftest.py` checks the suite against things that are true
by construction, and needs nothing but Python and numpy:

```bash
python3 tests/verification/selftest.py
```

```
  14 passed, 0 failed, 0 skipped, of 14
```

It exists because a wrong error measure would report a wrong convergence order
just as confidently as a right one, and nothing else here would notice. So the
closed-form flow rate is checked against the integral it stands for, the
pressure gradient against the second derivative of the profile, the order
estimator against synthetic sequences whose order is known, the mesh generator
against the mesh that ships with the example, and the reader and the error
measures against a VTK file written by the test itself.

The mesh check is worth singling out: the generator has to reproduce the
shipped 160 x 40 meshes exactly, otherwise a convergence study would be
comparing the solver against itself on a different problem.

It is also the only part of the suite that can run on every push, which is what
the continuous integration below is built around.

## Continuous integration

`.github/workflows/verification.yml` has two jobs, on different triggers
because they cost very different amounts.

| | Runs on | Takes |
|---|---|---|
| `selftest` | every push and pull request | seconds |
| `channel` | pushes to the main branches, pull requests, on demand | tens of minutes |

`channel` builds the container and runs the solver. PETSc alone accounts for
most of that, so the image build is cached by layer; the first three of the
four stages change rarely. The run output is uploaded as an artifact whether it
passes or fails, since the failing run is the one worth looking at. A pull
request labelled `docs-only` skips it.

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
convenience, and it is worth stating plainly what is being excluded and why.

**The column against the inlet is an artefact of the VTK output, not of the
solution.** HigFlow writes velocity as POINT_DATA, so `hig-flow-io.c`
interpolates the staggered facet values to the four corners of every cell with
`compute_facet_value_at_point`. For a cell against the inlet, two of those four
corners lie exactly on the boundary plane x = 0, where that interpolation has
only a one-sided stencil to work with.

Measuring the two sides of those same cells separately settles it:

```
  corners at x = 0.0500, one cell in    error 9.42e-05
  corners at x = 0.0000, on the inlet   error 2.56e+00
```

The two pairs belong to the same cells and are computed from the same facet
values. If the solution in the first cell were wrong, both pairs would be
wrong. Only the pair on the boundary is, and it reaches 2.70 against a u_max of
1.5, which no solution of this problem attains:

```
  y_corner    u at x = 0    u at x = 0.05    exact
   -0.9500       2.7043           0.1463    0.1463
   -0.9000       2.7043           0.2850    0.2850
   -0.8500       2.5731           0.4162    0.4163
```

Every other column is clean on both sides, at 7.4e-04 in mid-channel.

Two consequences follow for anything that reads these files. The cell average
of four corners mixes the bad pair with the good one, so the first column's
flow rate is wrong by a wide margin: the spread along the channel reads 42.7 %
with that column and 1.5e-05 without it.

And it dominates the norm badly enough to invert the convergence study. Over
the whole field the measured order is 0.19 and the error rises between the two
coarsest meshes, because the artefact grows with refinement while everything
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
column of corner values that the solver never uses.

This is worth fixing at the source, by giving the corner interpolation the
boundary condition instead of letting it extrapolate. That is a change to the
output path of every case rather than to the tests, so it is not made here.
Until it is, the checks and the convergence study are gated on the interior,
`verify.py check` prints the two-sided measurement above so the reason is
visible in the output, and the whole-field number is reported next to the
interior one rather than hidden.

## Files

| | |
|---|---|
| `run_verification.sh` | runs the case and checks it; the entry point for CI |
| `CMakeLists.txt` | registers the checks as CTest tests |
| `verification/selftest.py` | checks on the suite itself, no solver needed |
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
