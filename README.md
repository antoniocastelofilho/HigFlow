<h1 align="center">HigFlow</h1>

<p align="center">
  <strong>A finite-difference solver for incompressible Newtonian, viscoelastic
  and multiphase flows on hierarchical adaptive grids.</strong>
</p>

<p align="center">
  <img alt="Language: C" src="https://img.shields.io/badge/language-C99-4A8199">
  <img alt="Parallelism: MPI" src="https://img.shields.io/badge/parallel-MPI-6E9488">
  <img alt="Linear algebra: PETSc" src="https://img.shields.io/badge/solvers-PETSc%20%7C%20HYPRE-B09760">
  <img alt="Dimensions: 2D and 3D" src="https://img.shields.io/badge/dimensions-2D%20%7C%203D-BE7430">
</p>

<p align="center">
  <a href="README.pt-BR.md">Leia em português</a>
</p>

---

## What HigFlow is

HigFlow solves the incompressible Navier–Stokes equations by a projection method,
discretised with finite differences on a **staggered, hierarchically refined grid**.
Pressure and scalar properties live at cell centres; velocity components live on cell
faces.

What sets it apart from a general-purpose CFD code is its **rheology**. HigFlow was
built to simulate fluids whose stress does not follow a constant viscosity: polymer
solutions and melts, wormlike micellar solutions, thixotropic and elastoviscoplastic
materials, dense suspensions, and electrolytes driven by electro-osmosis. It carries a
broad library of differential and integral constitutive models, coupled to a
volume-of-fluid method for free surfaces and interfaces.

The code is organised in two layers:

| Layer | Responsibility |
|---|---|
| **HigTree** | The hierarchical adaptive grid: cell trees, domain decomposition across MPI ranks, moving least-squares interpolation between refinement levels, and the interface to linear solvers (PETSc, HYPRE, ViennaCL, SOR) |
| **HigFlow** | The physics: the Navier–Stokes projection loop, constitutive models, volume-of-fluid interface tracking, electro-osmotic coupling, and I/O |

HigFlow is developed at the **Institute of Mathematical and Computer Sciences (ICMC),
University of São Paulo**.

## Problems it is built for

- **Viscoelastic flow** through contractions, expansions and channels — corner
  vortices, stress boundary layers, high Weissenberg number behaviour
- **Free-surface and two-phase flow** — bubbles, droplets, dam break, jet breakup
- **Shear banding** in wormlike micellar solutions
- **Thixotropic and elastoviscoplastic** materials with evolving microstructure
- **Electro-osmotic flow** in microchannels, including electrolyte transport coupled to
  a viscoelastic solvent
- **Dense suspensions** exhibiting shear thickening

## Table of contents

- [Physical models](#physical-models)
- [Numerical methods](#numerical-methods)
- [Architecture](#architecture)
- [Gallery](#gallery)
- [Installation](#installation)
- [Running your first case](#running-your-first-case)
- [Repository layout](#repository-layout)
- [Documentation](#documentation)
- [Citing HigFlow](#citing-higflow)
- [Contributing](#contributing)
- [Authors](#authors)
- [License](#license)
