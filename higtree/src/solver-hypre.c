#include <assert.h>
#include <stdbool.h>
#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_mv.h>
#include <HYPRE_parcsr_ls.h>
#include <HYPRE_krylov.h>

#include "solver.h"

// This solver uses the generic inteface of HYPRE for solving matrices...

typedef struct solver_impl {
    solver s;

    // HYPRE specifics:

    int iupper, ilower;

    HYPRE_IJMatrix matrix;
    HYPRE_ParCSRMatrix assembled_matrix;
    bool matrix_assembled;
    bool matrix_new_assembly;

    HYPRE_IJVector rhs;
    HYPRE_ParVector assembled_rhs;
    bool rhs_assembled;

    HYPRE_IJVector sol;
    HYPRE_ParVector assembled_sol;
    bool sol_assembled;

    HYPRE_Solver ksp;
    HYPRE_Solver pc;
} solver_impl;

static inline solver_impl* scast(solver *cs){
	return (solver_impl*)(cs);
}

static void mat_reset(solver_impl *s)
{
    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, s->ilower, s->iupper, s->ilower, s->iupper, &s->matrix);
    HYPRE_IJMatrixSetObjectType(s->matrix, HYPRE_PARCSR);

    s->matrix_assembled = false;
}

// Assumes it will be called once, and at maximum once, before matrix elements are set
static void set_max_nonzeros(solver *cs)
{
    solver_impl *s = scast(cs);
    {
        size_t elem_count = cs->size;
        int vals[elem_count];
        for(size_t i = 0; i < elem_count; ++i) {
            vals[i] = cs->MaxNonZeros;
        }
        HYPRE_IJMatrixSetRowSizes(s->matrix, vals);
    }
    HYPRE_IJMatrixInitialize(s->matrix);
}

static void set_Ai(solver *cs, int i, int numjs, const int *j, const real *v)
{
    solver_impl *s = scast(cs);
    if(s->matrix_assembled) {
        // TODO: find a way to know when sparcity will not change,
        // so the matrix won't need to be recreated...
        HYPRE_IJMatrixDestroy(s->matrix);
        mat_reset(s);
        set_max_nonzeros(cs);
    }

    if(sizeof(real) != sizeof(double)) {
        double vals[numjs];
        for(int i = 0; i < numjs; ++i) {
            vals[i] = v[i];
        }
        HYPRE_IJMatrixSetValues(s->matrix, 1, &numjs, &i, j, vals);
    } else {
        HYPRE_IJMatrixSetValues(s->matrix, 1, &numjs, &i, j, v);
    }
}

static void set_bi(struct solver *cs, int i, real v)
{
    solver_impl *s = scast(cs);
    if(s->rhs_assembled) {
        HYPRE_IJVectorInitialize(s->rhs);
        s->rhs_assembled = false;
    }

    if(sizeof(double) == sizeof(real)) {
        HYPRE_IJVectorSetValues(s->rhs, 1, &i, &v);
    } else {
        double val = v;
        HYPRE_IJVectorSetValues(s->rhs, 1, &i, &val);
    }
}

static void add_bi(struct solver *cs, int i, real v)
{
    solver_impl *s = scast(cs);
    if(s->rhs_assembled) {
        HYPRE_IJVectorInitialize(s->rhs);
        s->rhs_assembled = false;
    }


    if(sizeof(double) == sizeof(real)) {
        HYPRE_IJVectorAddToValues(s->rhs, 1, &i, &v);
    } else {
        double val = v;
        HYPRE_IJVectorAddToValues(s->rhs, 1, &i, &val);
    }
}

static void set_xi(struct solver *cs, int i, real v)
{
    solver_impl *s = scast(cs);
    if(s->sol_assembled) {
        HYPRE_IJVectorInitialize(s->sol);
        s->sol_assembled = false;
    }

    if(sizeof(double) == sizeof(real)) {
        HYPRE_IJVectorSetValues(s->sol, 1, &i, &v);
    } else {
        double val = v;
        HYPRE_IJVectorSetValues(s->sol, 1, &i, &val);
    }
}

static void set_x(struct solver *cs, const real *x)
{
    solver_impl *s = scast(cs);
    if(s->sol_assembled) {
        HYPRE_IJVectorInitialize(s->sol);
        s->sol_assembled = false;
    }

    size_t size = cs->size;
    int indices[size];

    for(int i = 0; i < size; ++i) {
        indices[i] = i + s->ilower;
    }

    if(sizeof(double) == sizeof(real)) {
        HYPRE_IJVectorSetValues(s->sol, size, indices, x);
    } else {
        double val[size];
        for(int i = 0; i < size; ++i) {
            val[i] = x[i];
        }
        HYPRE_IJVectorSetValues(s->sol, size, indices, val);
    }
}

static real get_xi(struct solver *cs, int i)
{
    solver_impl *s = scast(cs);
    double val;

    HYPRE_IJVectorGetValues(s->sol, 1, &i, &val);

    return (real)val;
}

static void get_x(struct solver *cs, real *x)
{
    solver_impl *s = scast(cs);
    size_t size = cs->size;
    int indices[size];

    for(int i = 0; i < size; ++i) {
        indices[i] = i + s->ilower;
    }

    if(sizeof(double) == sizeof(real)) {
        HYPRE_IJVectorGetValues(s->sol, size, indices, x);
    } else {
        double vals[size];
        HYPRE_IJVectorGetValues(s->sol, size, indices, vals);

        for(int i = 0; i < size; ++i) {
            x[i] = vals[i];
        }
    }
}

static void get_x_scatter(struct solver *cs, int size, const int *idx, real *x)
{
    solver_impl *s = scast(cs);

    if(sizeof(double) == sizeof(real)) {
        HYPRE_IJVectorGetValues(s->sol, size, idx, x);
    } else {
        double vals[size];
        HYPRE_IJVectorGetValues(s->sol, size, idx, vals);

        for(int i = 0; i < size; ++i) {
            x[i] = vals[i];
        }
    }
}

static void assemble_matrix(struct solver *cs)
{
    solver_impl *s = scast(cs);
    HYPRE_IJMatrixAssemble(s->matrix);
    HYPRE_IJMatrixGetObject(s->matrix, (void **) &s->assembled_matrix);
    s->matrix_assembled = true;
    s->matrix_new_assembly = true;

    HYPRE_BiCGSTABSetTol(s->ksp, cs->relative_tolerance);
    HYPRE_BiCGSTABSetAbsoluteTol(s->ksp, cs->absolute_tolerance);
    HYPRE_BiCGSTABSetMaxIter(s->ksp, cs->MaxIteration);
}

static void assemble_rhs(struct solver *cs)
{
    solver_impl *s = scast(cs);
    HYPRE_IJVectorAssemble(s->rhs);
    HYPRE_IJVectorGetObject(s->rhs, (void **) &s->assembled_rhs);
    s->rhs_assembled = true;
}

static void assemble_output(struct solver *cs)
{
    solver_impl *s = scast(cs);
    HYPRE_IJVectorAssemble(s->sol);
    HYPRE_IJVectorGetObject(s->sol, (void **) &s->assembled_sol);
    s->sol_assembled = true;
}

static void solve(struct solver *cs)
{
    solver_impl *s = scast(cs);

    if(s->matrix_new_assembly) {
	    HYPRE_BiCGSTABSetup(s->ksp, (HYPRE_Matrix)s->assembled_matrix, (HYPRE_Vector)s->assembled_rhs, (HYPRE_Vector)s->assembled_sol);
	    s->matrix_new_assembly = false;
    }

    HYPRE_ClearError(0);
    HYPRE_Int ret = HYPRE_BiCGSTABSolve(s->ksp, (HYPRE_Matrix)s->assembled_matrix, (HYPRE_Vector)s->assembled_rhs, (HYPRE_Vector)s->assembled_sol);

    cs->converged = !HYPRE_CheckError(ret, HYPRE_ERROR_CONV);

    HYPRE_BiCGSTABGetNumIterations(s->ksp, &cs->SolverIterations);
    HYPRE_BiCGSTABGetFinalRelativeResidualNorm(s->ksp, &cs->SolverResidual);
}

static void destroy(struct solver *cs)
{
    solver_impl *s = scast(cs);
    HYPRE_IJMatrixDestroy(s->matrix);

    HYPRE_IJVectorDestroy(s->rhs);
    HYPRE_IJVectorDestroy(s->sol);

    HYPRE_ParCSRBiCGSTABDestroy(s->ksp);
    HYPRE_BoomerAMGDestroy(s->pc);
}

static const _solver_vtable vtable = {
	.set_maxnonzeros = set_max_nonzeros,
	.set_Ai = set_Ai,
	.set_bi = set_bi,
	.add_bi = add_bi,
	.set_xi = set_xi,
	.set_x = set_x,
	.get_xi = get_xi,
	.get_x = get_x,
	.get_x_scatter = get_x_scatter,
	.assemble_matrix = assemble_matrix,
	.assemble_rhs = assemble_rhs,
	.assemble_output = assemble_output,
	.solve = solve,
	.destroy = destroy,
};

solver * _slv_create_hypre(int ilower, int iupper)
{
    // Initialize solver
    DECL_AND_ALLOC(solver_impl, s, 1);

    s->s.vtable = &vtable;

    s->iupper = iupper;
    s->ilower = ilower;
    mat_reset(s);

    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &s->rhs);
    HYPRE_IJVectorSetObjectType(s->rhs, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(s->rhs);

    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &s->sol);
    HYPRE_IJVectorSetObjectType(s->sol, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(s->sol);

    HYPRE_BoomerAMGCreate(&s->pc);
    HYPRE_BoomerAMGSetNumSweeps(s->pc, 1);
    HYPRE_BoomerAMGSetTol(s->pc, 0.0);
    HYPRE_BoomerAMGSetMaxIter(s->pc, 1);

    HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, &s->ksp);
    HYPRE_BiCGSTABSetPrecond(s->ksp, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
		     (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, s->pc);
    s->rhs_assembled = s->sol_assembled = s->matrix_assembled = false;
    s->matrix_new_assembly = false;

    return &s->s;
}
