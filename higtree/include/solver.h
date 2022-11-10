#ifndef SOLVER_H
#define SOLVER_H

#include "utils.h"
#include "higtree.h"

#include<stdbool.h>
#include<stdio.h>
#include<math.h>

struct solver;

// Keep the same order on both following lists:
static const char* engine_names[] = {
	"hypre",
	"petsc",
	"hollow",
#ifdef USE_SOR
	"sor",
#endif
#ifdef USE_VIENNACL
	"viennacl",
#endif
	"debug-write"
};

typedef enum {
    SOLVER_ANY,
    SOLVER_HYPRE,
    SOLVER_PETSC,
    SOLVER_HOLLOW,
#ifdef USE_SOR
    SOLVER_SOR,
#endif
#ifdef USE_VIENNACL
    SOLVER_VIENNACL,
#endif
    SOLVER_DEBUG_WRITE
} SolverEngine;

//! Internal struct to hold the implementation dependant pointer functions
typedef struct _solver_vtable {
    void (*set_maxnonzeros)(struct solver *cs);
    void (*set_Ai)(struct solver *cs, int i, int numjs, const int *j, const real *v);
    void (*set_bi)(struct solver *cs, int i, real v);
    void (*add_bi)(struct solver *cs, int i, real v);
    void (*set_xi)(struct solver *cs, int i, real v);
    void (*set_x)(struct solver *cs, const real *x);
    real (*get_xi)(struct solver *cs, int i);
    void (*get_x)(struct solver *cs, real *x);
    void (*get_x_scatter)(struct solver *cs, int size, const int *idx, real *x);

    void (*assemble_matrix)(struct solver *cs);
    void (*assemble_rhs)(struct solver *cs);
    void (*assemble_output)(struct solver *cs);

    void (*solve)(struct solver *cs);

    void (*destroy)(struct solver *cs);
} _solver_vtable;

//! Defines a solver for a linear system.
typedef struct solver {
    const _solver_vtable *vtable;

    // Define constants
    int  MaxIteration;
    int  MaxNonZeros;
    real relative_tolerance;
    real absolute_tolerance;

    // Cache stuff
    int lasti;
    int lastj;
    int *j;
    real *vs;
    int ntasks;

    // Post-solver varibles
    int  SolverIterations;
    real SolverResidual;
    bool converged;

    size_t size; //!< Process-local size of the vector

    int imposed_line;
} solver;

//! Creates a solver context
solver* slv_create(SolverEngine engine, int first_gid, size_t local_size);

//! Returns the local size of the domain worked by this solver.
size_t slv_get_local_size(solver *s);

//! Sets the maximum number of iteration of a iterative method.
void slv_set_maxiteration( solver *s, int _maxit );

//! Sets the number of non zero values per line.
void slv_set_maxnonzeros( solver *s, int _maxnz );

//! Sets the tolerance error of the iterative method.
//
// Uses the same value for relative and absolute tolerance.
void slv_set_eps( solver *s, real _eps );

//! Sets the absolute tolerance error of the iterative method.
void slv_set_abs_tol( solver *s, real tol);

//! Sets the relative tolerance error of the iterative method.
void slv_set_rel_tol( solver *s, real tol);

//! Sets the value of the j-th element of the i-th row.
void slv_set_Aij(solver *s, int i, int j, real v);

//! Sets the values of the i-th row.
void slv_set_Ai(solver *s, int i, int numjs, int *j, real *v);

//! Sets the right hand side value of the i-th equation.
void slv_set_bi(solver *s, int i, real v);

//! Adds v to the right hand side value of the i-th equation.
void slv_add_bi(solver *s, int i, real v);

//! Sets the value of the i-th variable of the solution.
void slv_set_xi(solver *s, int i, real v);

//! Sets the solution.
void slv_set_x(solver *s, real *x);

//! Gets the value of the i-th variable of the solution.
real slv_get_xi(solver *s, int i);

//! \brief Gets the values of the solution. x should be long enough to receive all local elements of the solution.
void slv_get_x(solver *s, real *x);

//! \brief Gets 'size' values from the solution, in the specified order. The order of the elements returned in 'x' must be specified by 'idx'.
void slv_get_x_scatter(solver *s, int size, const int *idx, real *x);

//! Impose value for one unknown. Use before assembly.
void slv_impose_value(solver *s, int line, real value);

/** Assemble matrix for use in solver.
 *
 * The matrix can be assembled only once while
 * the tree structure is unchanged.
 */
void slv_assemble_matrix(solver *s);

/** Assemble the right hand side vector.
 *
 * Will possibly change at every time step.
 */
void slv_assemble_rhs(solver *s);

/** Assemble the right hand side vector.
 *
 * Will possibly change at every time step.
 */
void slv_assemble_output(solver *s);

//! Assembles the matrixes and vectors prior solving.
void slv_assemble(solver *s);

//! Solves the system.
void slv_solve(solver *s);

//! Returns true if last solver execution has converged.
bool slv_has_converged(solver *s);

//! Returns the number of iterations performed by the solver
int slv_get_iteration_count(solver *s);

//! Returns the achieved relative residual norm
// The relative residual norm is given by the residual norm divided by
// vectorial norm of the right-hand side values.
real slv_get_relative_residual_norm(solver *s);

//! Destroys the solver.
void slv_destroy(solver *s);

#endif
