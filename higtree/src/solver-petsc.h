// PETSc-specific extension to the solver interface: naming index ranges as fields,
// which is what lets PETSc's block preconditioners see the velocity/pressure split.
//
// THE PETSc VERSION IS PINNED, and the pin is not visible from this file: the build
// expects 3.25.4 (see varsrc).  Note that bibliotecas/petsc-3.14.0.tar.gz IS
// VERSIONED IN THIS REPOSITORY and is a DIFFERENT, incompatible version -- a fresh
// clone therefore contains a PETSc tarball that must not be used to satisfy this
// dependency.  instalar.sh says so explicitly when PETSc is missing.

#ifndef SOLVER_PETSC_H
#define SOLVER_PETSC_H

#include "solver.h"
#include <petsc.h>
#include <petscksp.h>

void slv_PETSc_set_field(solver *cs, const char *name, int start, size_t size);

#endif
