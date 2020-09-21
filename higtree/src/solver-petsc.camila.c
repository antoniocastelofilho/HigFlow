#include "solver.h"
#include "Debug-c.h"

#include <petsc.h>
#include <petscksp.h>

// At the first time the solver is set up,
// this structure is filled with the matrix
// coefficients in order to create an space
// efficient PETSc matrix.
typedef struct mat_prealloc
{
    PetscInt **col_ids;
    PetscReal **vals;

    PetscInt *d_nnz;
    PetscInt *o_nnz;

    int size;
    int max_gid;
} mat_prealloc;

static mat_prealloc* mat_prealloc_create(int start, int size)
{
    assert(size >= 0);

    DECL_AND_ALLOC(mat_prealloc, mp, 1);
    mp->size = size;
    mp->max_gid = start + size - 1;

    ALLOC_INFER(mp->col_ids, size);
    memset(mp->col_ids, 0, size * sizeof *mp->col_ids);

    ALLOC_INFER(mp->vals, size);
    memset(mp->vals, 0, size * sizeof *mp->vals);

    ALLOC_INFER(mp->d_nnz, size);
    memset(mp->d_nnz, 0, size * sizeof *mp->d_nnz);

    ALLOC_INFER(mp->o_nnz, size);
    memset(mp->o_nnz, 0, size * sizeof *mp->o_nnz);

    return mp;
}

static void mat_prealloc_destroy(mat_prealloc *mp)
{
    for(int i = 0; i < mp->size; ++i) {
	free(mp->col_ids[i]);
	free(mp->vals[i]);
    }
    free(mp->col_ids);
    free(mp->vals);
    free(mp->d_nnz);
    free(mp->o_nnz);

    free(mp);
}

// Petsc solver internal struct.
typedef struct _solver_impl_petsc {
    solver s;

    // Petsc stuff
    Mat A;
    Vec x;
    Vec b;

    // Used to calculate residual
    Vec tmp;

    KSP ksp;
    PC  pc;

    int first_gid;

    void (*pre_solve_func)(struct _solver_impl_petsc *s);

    //! Non-NULL only in first iteration.
    mat_prealloc *mp;
} _solver_impl_petsc;

static void mat_prealloc_set_Ai(solver *cs, int gid, int numjs, const int *j, const real *v)
{
    _solver_impl_petsc *s = (_solver_impl_petsc*)cs;
    assert(s->mp);

    int lid = gid - s->first_gid;

    // We do not support setting the same row twice!
    // If this is ever needed, don't forget to adjust
    // o_nnz and d_nnz when removing previously placed row.
    assert(!s->mp->col_ids[lid] && !s->mp->vals[lid]);

    DECL_AND_ALLOC(PetscReal, vals, numjs);
    memcpy(vals, v, numjs * sizeof *vals);
    s->mp->vals[lid] = vals;

    DECL_AND_ALLOC(PetscInt, cols, numjs);
    for(int i = 0; i < numjs; ++i) {
	cols[i] = j[i];
	if(j[i] < s->first_gid || j[i] > s->mp->max_gid) {
	    ++s->mp->o_nnz[lid];
	} else {
	    ++s->mp->d_nnz[lid];
	}
    }
    s->mp->col_ids[lid] = cols;
}

static void mat_prealloc_assemble_matrix(solver *cs)
{
    _solver_impl_petsc *s = (_solver_impl_petsc*)cs;
    assert(s->mp);

    // Preallocate matrix inside PETSc:
    if (cs->ntasks == 1) {
        MatSeqAIJSetPreallocation(s->A, cs->MaxNonZeros, s->mp->d_nnz);
    } else {
        MatMPIAIJSetPreallocation(s->A, cs->MaxNonZeros, s->mp->d_nnz, cs->MaxNonZeros, s->mp->o_nnz);
    }

    // Set matrix values:
    for(int i = 0; i < s->mp->size; ++i) {
	assert(s->mp->col_ids[i] && s->mp->vals[i]);

	int row_id = i + s->first_gid;
	int num_cols = s->mp->d_nnz[i] + s->mp->o_nnz[i];
	MatSetValues(s->A, 1, &row_id, num_cols, s->mp->col_ids[i], s->mp->vals[i], INSERT_VALUES);
    }

    // Destroy prealloc and set vtable to not use it anymore.
    mat_prealloc_destroy(s->mp);
    s->mp = NULL;

    extern const _solver_vtable _petsc_vtable;
    cs->vtable = &_petsc_vtable;

    // Finish setting up the matrix:
    cs->vtable->assemble_matrix(cs);
}

static void petsc_set_max_nonzeros(solver *cs)
{
    /* Don't use this value anymore. */
}

static void petsc_slv_set_Ai(solver *cs, int i, int numjs, const int *j, const real *v)
{
    _solver_impl_petsc *s = (_solver_impl_petsc*)cs;
    // MatSetValues(s->A, 1, &i, numjs, j, v, INSERT_VALUES);
    // camila
    if(numjs==-1)
    {
        MatSetOption(s->A,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);
        VecAssemblyBegin(s->b);
        VecAssemblyEnd(s->b);
        MatZeroRows(s->A, 1, &i, 1.0, PETSC_NULL, PETSC_NULL);
        printf("singularidade\n");
    }
    else
    {
        MatSetValues(s->A, 1, &i, numjs, j, v, INSERT_VALUES);
    }
}

static void petsc_slv_set_bi(solver *cs, int i, real v)
{
    VecSetValue(((_solver_impl_petsc*)cs)->b, i, v, INSERT_VALUES);
}

static void petsc_slv_add_bi(solver *cs, int i, real v)
{
    VecSetValue(((_solver_impl_petsc*)cs)->b, i, v, ADD_VALUES);
}

static void petsc_slv_set_xi(solver *cs, int i, real v)
{
    VecSetValue(((_solver_impl_petsc*)cs)->x, i, v, INSERT_VALUES);
}

static void petsc_slv_set_x(solver *cs, const real *x)
{
    _solver_impl_petsc *s = (_solver_impl_petsc*)cs;
    int ind[cs->size];
    for (int c = 0; c < cs->size; c++) ind[c] = s->first_gid + c;

    VecSetValues(s->x, cs->size, ind, x, INSERT_VALUES);
}

static real petsc_slv_get_xi(solver *cs, int i)
{
    real v;

    VecGetValues(((_solver_impl_petsc*)cs)->x, 1, &i, &v);

    return v;
}

static void petsc_slv_get_x_scatter(solver *cs, int size, const int *idx, real *x) {
    VecGetValues(((_solver_impl_petsc*)cs)->x, size, idx, x);
}

static void petsc_slv_get_x(solver *cs, real *x) {
    _solver_impl_petsc *s = (_solver_impl_petsc*)cs;
    int ind[cs->size];
    for (int c = 0; c < cs->size; c++) ind[c] = s->first_gid + c;

    VecGetValues(s->x, cs->size, ind, x);
}

static void petsc_slv_assemble_matrix(solver *cs)
{
    _solver_impl_petsc *s = (_solver_impl_petsc*)cs;

    MatAssemblyBegin(s->A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd  (s->A, MAT_FINAL_ASSEMBLY);

#if ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR > 4)) || (PETSC_VERSION_MAJOR > 3)
    KSPSetOperators(s->ksp, s->A, s->A);
#else
    KSPSetOperators(s->ksp, s->A, s->A, DIFFERENT_NONZERO_PATTERN);
#endif

    KSPSetTolerances(s->ksp, cs->relative_tolerance, cs->absolute_tolerance,
	PETSC_DEFAULT, cs->MaxIteration);
}

static void petsc_slv_assemble_rhs(solver *cs) {
    _solver_impl_petsc *s = (_solver_impl_petsc*)cs;
    VecAssemblyBegin(s->b);
    VecAssemblyEnd  (s->b);
}

static void petsc_slv_assemble_output(solver *cs) {
    _solver_impl_petsc *s = (_solver_impl_petsc*)cs;
    VecAssemblyBegin(s->x);
    VecAssemblyEnd  (s->x);
}

static real calc_res_norm(_solver_impl_petsc *s)
{
	MatMult(s->A, s->x, s->tmp);
	VecAXPY(s->tmp, -1.0, s->b);

	real abs_res_norm;
	VecNorm(s->tmp, NORM_INFINITY, &abs_res_norm);

	return abs_res_norm;
}

static PetscErrorCode abs_convergence_test(KSP ksp, PetscInt it, PetscReal rnorm, KSPConvergedReason *reason, void *cctx)
{
	_solver_impl_petsc *s = (_solver_impl_petsc*)cctx;

	real abs_res_norm = calc_res_norm(s);

	if(abs_res_norm <= s->s.absolute_tolerance) {
		*reason = KSP_CONVERGED_ATOL;
	} else {
		*reason = KSP_CONVERGED_ITERATING;
	}

	return 0;
}

static void petsc_slv_solve(solver *cs) {
    _solver_impl_petsc *s = (_solver_impl_petsc*)cs;

    if(s->pre_solve_func) {
	    s->pre_solve_func(s);
    }

    // Solve system
    KSPSolve(s->ksp, s->b, s->x);

    KSPGetIterationNumber(s->ksp, &cs->SolverIterations);

    real abs_res_norm;
    KSPGetResidualNorm(s->ksp, &abs_res_norm);
    real rhs_norm;
    VecNorm(s->b, NORM_2, &rhs_norm);
    cs->SolverResidual = abs_res_norm / rhs_norm;

    KSPConvergedReason r;
    KSPGetConvergedReason(s->ksp, &r);
    cs->converged = r > 0;

    DEBUG_EXEC({
	KSPConvergedReason r;
	KSPGetConvergedReason(s->ksp, &r);
	print0f("PETSc DEBUG: Converged reason: %d, real abs norm: %g\n",
		r, calc_res_norm(s));
    });
}

static void petsc_slv_destroy(solver *cs) {
    _solver_impl_petsc *s = (_solver_impl_petsc*)cs;

    // Destroy Petsc stuff
    MatDestroy(&s->A);
    VecDestroy(&s->x);
    VecDestroy(&s->b);

    VecDestroy(&s->tmp);

    KSPDestroy(&s->ksp); // PC is destroyed with KSP11
}

static const _solver_vtable prealloc_vtable = {
	.set_maxnonzeros = petsc_set_max_nonzeros,
	.set_Ai = mat_prealloc_set_Ai,
	.set_bi = petsc_slv_set_bi,
	.add_bi = petsc_slv_add_bi,
	.set_xi = petsc_slv_set_xi,
	.set_x = petsc_slv_set_x,
	.get_xi = petsc_slv_get_xi,
	.get_x = petsc_slv_get_x,
	.get_x_scatter = petsc_slv_get_x_scatter,
	.assemble_matrix = mat_prealloc_assemble_matrix,
	.assemble_rhs = petsc_slv_assemble_rhs,
	.assemble_output = petsc_slv_assemble_output,
	.solve = petsc_slv_solve,
	.destroy = petsc_slv_destroy,
};

const _solver_vtable _petsc_vtable = {
	.set_maxnonzeros = petsc_set_max_nonzeros,
	.set_Ai = petsc_slv_set_Ai,
	.set_bi = petsc_slv_set_bi,
	.add_bi = petsc_slv_add_bi,
	.set_xi = petsc_slv_set_xi,
	.set_x = petsc_slv_set_x,
	.get_xi = petsc_slv_get_xi,
	.get_x = petsc_slv_get_x,
	.get_x_scatter = petsc_slv_get_x_scatter,
	.assemble_matrix = petsc_slv_assemble_matrix,
	.assemble_rhs = petsc_slv_assemble_rhs,
	.assemble_output = petsc_slv_assemble_output,
	.solve = petsc_slv_solve,
	.destroy = petsc_slv_destroy,
};


solver * _slv_create_petsc(int first_gid, int size)
{
    // Try to initialize PETSc
    void _try_initialize_petsc(void);
    _try_initialize_petsc();

    // Allocate solver structure
    DECL_AND_ALLOC(_solver_impl_petsc, s, 1);

    s->first_gid = first_gid;

    MPI_Comm_size(MPI_COMM_WORLD, &s->s.ntasks);

    //MatCreateSeqAIJ(PETSC_COMM_WORLD,s->size,s->size,s->MaxNonZeros,PETSC_NULL,&s->A);
    MatCreate(PETSC_COMM_WORLD,&s->A);

    MatSetSizes(s->A, size, size, PETSC_DECIDE, PETSC_DECIDE);
    MatSetFromOptions(s->A);

    VecCreate(PETSC_COMM_WORLD, &s->x);
    VecSetSizes(s->x, size, PETSC_DECIDE);
    VecSetFromOptions(s->x);

    VecCreate(PETSC_COMM_WORLD, &s->b);
    VecSetSizes(s->b, size, PETSC_DECIDE);
    VecSetFromOptions(s->b);

    VecCreate(MPI_COMM_WORLD, &s->tmp);
    VecSetSizes(s->tmp, size, PETSC_DECIDE);
    VecSetFromOptions(s->tmp);

    KSPCreate(PETSC_COMM_WORLD, &s->ksp);
    KSPGetPC(s->ksp, &s->pc);

    //KSPSetType(s->ksp, KSPGMRES);
    KSPSetType(s->ksp, KSPBCGS);

    // camila comentei
    KSPSetNormType(s->ksp, KSP_NORM_UNPRECONDITIONED);

    //KSPSetInitialGuessNonzero(s->ksp, PETSC_FALSE);

    PCSetType(s->pc, PCHYPRE);

    KSPSetFromOptions(s->ksp);
    PCSetFromOptions(s->pc);

    KSPSetConvergenceTest(s->ksp, abs_convergence_test, s, NULL);

    s->pre_solve_func = NULL;

    // Initially, setup preallocation:
    s->mp = mat_prealloc_create(first_gid, size);
    s->s.vtable = &prealloc_vtable;

    return &s->s;
}

typedef struct {
	void *default_ctx;
	real tolerance;
} singular_converged_ctx;

static PetscErrorCode singular_converged_ctx_destroy(void *ctx)
{
	singular_converged_ctx *c = ctx;
	PetscErrorCode ret = KSPConvergedDefaultDestroy(c->default_ctx);
	free(c);
	return ret;
}

static PetscErrorCode singular_converged(KSP ksp, PetscInt n,
	PetscReal rnorm, KSPConvergedReason *reason, void *ctx)
{
	singular_converged_ctx *c = ctx;
	PetscErrorCode ret = KSPConvergedDefault(ksp, n, rnorm, reason, c->default_ctx);

	if(*reason == 0) {
		PetscReal emax, emin;
		KSPComputeExtremeSingularValues(ksp, &emax, &emin);
		if(fabs(emax/emin) > c->tolerance) {
			*reason = KSP_DIVERGED_BREAKDOWN;
		}
	}

	return ret;
}

static void set_sub_tolerances(_solver_impl_petsc *s)
{
	int n_subksp;
	KSP *subksp = NULL;

	PCFieldSplitGetSubKSP(s->pc, &n_subksp, &subksp);
	assert(n_subksp == 2);

	const real subEPS = s->s.absolute_tolerance * 10;
	for(int i = 0; i < n_subksp; ++i) {
		KSPSetTolerances(subksp[i], s->s.relative_tolerance, subEPS,
			PETSC_DEFAULT, s->s.MaxIteration);
	}

	PetscFree(subksp);
}

static void set_singular_convergence_test(_solver_impl_petsc *s)
{
	// Setup custom convergence test.
	KSPSetUp(s->ksp);
	PCSetUp(s->pc);

	int n_subksp;
	KSP *subksp = NULL;

	PCFieldSplitGetSubKSP(s->pc, &n_subksp, &subksp);
	assert(n_subksp == 2);

	KSPSetComputeSingularValues(subksp[1], PETSC_TRUE);

	DECL_AND_ALLOC(singular_converged_ctx, cctx, 1);
	KSPConvergedDefaultCreate(&cctx->default_ctx);
	cctx->tolerance = 1e4;

	KSPSetConvergenceTest(subksp[1], singular_converged, cctx,
		singular_converged_ctx_destroy);

	PetscFree(subksp);

	set_sub_tolerances(s);
	s->pre_solve_func = set_sub_tolerances;
}

solver * _slv_create_hollow_setup_petsc(int first_gid, int _size)
{
	// TODO: find a way to get and restore original options
	//char *orig_opts = NULL;
	//PetscOptionsGetAll(NULL, &orig_opts);

	// Set the options to use Schur complement:
	PetscOptionsInsertString(
#if ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 7)) || (PETSC_VERSION_MAJOR > 3)
		NULL,
#endif
		"-ksp_type fgmres "
		"-pc_type fieldsplit "
		"-pc_fieldsplit_type schur "
		"-pc_fieldsplit_schur_fact_type full "
		"-pc_fieldsplit_schur_precondition self "
		"-pc_fieldsplit_detect_saddle_point "
		"-fieldsplit_0_ksp_type bcgs "
		"-fieldsplit_0_pc_type hypre "
		"-fieldsplit_1_ksp_type gmres "
		"-fieldsplit_1_pc_type lsc "
		"-fieldsplit_1_lsc_pc_type hypre "
		"-fieldsplit_1_lsc_pc_hypre_boomeramg_cycle_type w");

	_solver_impl_petsc *s =
		(_solver_impl_petsc*)_slv_create_petsc(first_gid, _size);

	// Restore the original options:
	//PetscOptionsClear(NULL);
	//PetscOptionsInsertString(NULL, orig_opts);

	//PetscFree(orig_opts);

	// Configure to use custom convergence test.
	s->pre_solve_func = set_singular_convergence_test;

	return &s->s;
}


void slv_PETSc_set_field(solver *cs, const char *name, int start, size_t size)
{
	// Ignore if this is not a PETSc solver.
	if(cs->vtable->solve != petsc_slv_solve) {
		return;
	}

	_solver_impl_petsc *s = (_solver_impl_petsc*)cs;

	// This function is only relevant on PCFIELDSPLIT
	// so exit if using different preconditioner.
	PCType type;
	PCGetType(s->pc, &type);
	if(strcmp(type, PCFIELDSPLIT) != 0) {
		return;
	}

	// Build the index set
	IS is;
	{
		PetscInt *idx;
		PetscMalloc(size * sizeof *idx, &idx);

		for(size_t i = 0; i < size; ++i) {
			idx[i] = i + start;
		}

		ISCreateGeneral(MPI_COMM_WORLD, size, idx, PETSC_OWN_POINTER, &is);
	}

	// Set the field.
	PCFieldSplitSetIS(s->pc, name, is);

	ISDestroy(&is);
}
