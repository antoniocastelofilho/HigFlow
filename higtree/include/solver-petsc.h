#ifndef SOLVER_PETSC_H
#define SOLVER_PETSC_H

#include "solver.h"

void slv_PETSc_set_field(solver *cs, const char *name, int start, size_t size);

#endif
