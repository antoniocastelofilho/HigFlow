# PR draft - E02

**Branch:** `juniormar/02-readme` para `antoniocastelofilho/HigFlow:master`
**Title:** `Rewrite README in English with project overview, model reference and quick start`

---

## Summary

Replaces the README with an English document that describes what the solver does, which
models and methods it implements, how to install and run it on Linux and on Windows
through WSL2, and how the code is organised. The Portuguese text is preserved as
`README.pt-BR.md` with a corrupted passage repaired.

## Motivation

The previous README was 92 lines, in Portuguese only, with no description of the physics
the code solves, no images, and no reference to the models it implements. Someone
arriving at the repository could not tell whether it was relevant to their problem.

It also contained a corrupted passage where a PETSc `configure` line had been pasted
into the middle of the word *realizar*, splitting a sentence across a shell command and
a `sleep` call:

```
Após re./configure --prefix=/opt/petsc-3.14.0-openmnpi-hypre-hdf5 ... ;
sleep 5alizar um dos passos anterior você já pode utilizar o sistema!
```

## Changes

- **`README.md`** - rewritten in English, 12 sections
- **`README.pt-BR.md`** - the Portuguese text, with the corrupted sentence repaired
- **`docs/images/architecture.svg`** - a diagram of the three software layers and the
  six stages of one projection-method time step

Existing figures in `higtree/doc/figuras/` are referenced in place, so no binary is
duplicated into the tree.

### Content

The model tables list every model reachable from the configuration file, with the exact
key that selects it. Keys and model sets were read from the enumerations in
`higflow/src/hig-flow-kernel.h` rather than from the prose documentation, which in
places describes as unimplemented models that have since been implemented.

Standard references are given where the formulation is unambiguous. Where a model
follows a variant whose authoritative source is better supplied by its author, it is
listed by name and the gap is stated, rather than a reference being guessed. A note
invites contributions completing the table.

### Two behaviours documented that were previously unstated

- **The second-order convective stencil degrades to first order at boundaries.**
  `hig-flow-discret.c:131-138` demotes `SECOND_ORDER` to `FIRST_ORDER` where the stencil
  crosses a boundary. This looks deliberate, and it matters to anyone measuring
  convergence order, so the README says so.
- **`forth_order` is accepted but never applied.** `ORDER4` appears only in the YAML
  read and write paths (`hig-flow-io.c:6229` and `:8015`) and in the enum declaration.
  No discretisation routine consults it, so selecting it silently yields second order.

### Three statements that are uncomfortable but accurate

Included because a user needs them, and phrased factually:

1. The bundled installer passes `--with-debubbing=yes` to PETSc's `configure`, which is
   not a recognised option; it installs three MPI implementations side by side; and it
   sets `PKG_CONFIG_PATH` to a file rather than a directory.
2. The committed `varsrc` points `PETSC_DIR` at
   `bibliotecas/petsc-3.14.0/x86_64` with `PETSC_ARCH=arch-linux-c-debug`, while the
   installer installs to `/opt/petsc-3.14.0-openmnpi-hypre-hdf5` with
   `PETSC_ARCH=x86_64`. Following the instructions verbatim does not produce a working
   environment.
3. There is no `LICENSE` file, so default copyright applies and the terms of use are
   undefined.

The first two are corrected in a follow-up change to the installer; the README flags
them so users are not stranded meanwhile. The third is stated because users need to
know where they stand - the choice of licence is left entirely to the owners.

### On the absence of a native Windows build

The README explains it by its actual cause rather than leaving it unmentioned: OpenMPI
and libfyaml have no supported Windows port, and `CMakeLists.txt` requires `libnuma`,
which exists only on Linux, so configuration fails before anything compiles. WSL2 is
documented as the supported route.

## Verification

- Every image path, file link and internal anchor resolves
- All markdown tables well-formed
- The SVG validates as XML and no element exceeds its viewBox
- Model keys, method keys and enum values cross-checked against
  `higflow/src/hig-flow-kernel.h`
- Contributor names reproduced exactly as `git log` records them

## Notes for reviewers

- **Badges** are static and truthful (language, parallelism, solvers, dimensions). No
  CI badge is included because there is no CI yet.
- **The gallery holds one existing figure** and a note that it will be expanded with
  cases that ship with the repository, each with the command that reproduces it. Adding
  reproducible result figures needs a working build environment and is a separate
  change.
- **`CITATION.cff` and `CONTRIBUTING.md` are named as absent** in their sections rather
  than being invented here.
