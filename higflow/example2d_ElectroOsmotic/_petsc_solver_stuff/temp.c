#include "solver-petsc.c"

/*
create
set_maxnonzeros

set_bi
set_Ai

assemble
solve
load_from_solver || get_xi
*/


void initial () {
    int first_gid;
    int size;

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

    KSPSetNormType(s->ksp, KSP_NORM_UNPRECONDITIONED);

    //KSPSetInitialGuessNonzero(s->ksp, PETSC_FALSE);

    PCSetType(s->pc, PCHYPRE);

    KSPSetFromOptions(s->ksp);
    PCSetFromOptions(s->pc);

    real rtol, atol, dtol;
    int maxits;
    KSPGetTolerances(s->ksp, &rtol, &atol, &dtol, &maxits);
    print0f("PETSc: atol=%g, rtol=%g, dtol=%g, maxits=%d\n", atol, rtol, dtol, maxits);

    slv_set_rel_tol(&s->s, rtol);
    slv_set_abs_tol(&s->s, atol);
    slv_set_maxiteration(&s->s, maxits);

    //KSPSetConvergenceTest(s->ksp, rel_convergence_test, s, NULL);

    s->pre_solve_func = NULL;

    // Initially, setup preallocation:
    
    {
        int start = first_gid;

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

        s->mp = mp;
    }
    s->s.vtable = &prealloc_vtable;


    {
        // Generic stuff
        s->s.size         = size;
        s->s.MaxNonZeros  = 0;
        if(!SOLVER_PETSC) { // already set on petsc
            s->s.MaxIteration = 10000;
            s->s.absolute_tolerance = 1e-8;
            s->s.relative_tolerance = 1e-8;
        }

        // Cache stuff
        s->s.lasti = -1;
        s->s.lastj = 0;
        s->s.j = NULL;
        s->s.vs = NULL;

        // No imposed value:
        s->s.imposed_line = -1;
    }

    {int _maxnz = 800;
        s->s.MaxNonZeros = _maxnz;
        s->s.vtable->set_maxnonzeros(&s->s); // no use
    }

    

}

void during_loop (_solver_impl_petsc *s) {

    {int i; real v;
        assert(i >= 0);

        if(s->s.imposed_line != i)
            VecSetValue(s->b, i, v, INSERT_VALUES);
    }

    {int i; int numjs; const int *j; const real *v;

        assert(i >= 0);
        for(unsigned k = 0; k < numjs; ++k) {
            assert(j[k] >= 0);
        }

        // Preallocate before setting values
        if (s->s.MaxNonZeros == 0) {
            slv_set_maxnonzeros(&s->s, 1000);
        }

        MatSetValues(s->A, 1, &i, numjs, j, v, INSERT_VALUES);
    }

}

void first_iteration_cache (_solver_impl_petsc *s) {

    { int gid; int numjs; const int *j; const real *v; // called in the for loop
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

    {
        assert(s->mp);

        // Preallocate matrix inside PETSc:
        if (s->s.ntasks == 1) {
            MatSeqAIJSetPreallocation(s->A, s->s.MaxNonZeros, s->mp->d_nnz);
        } else {
            MatMPIAIJSetPreallocation(s->A, s->s.MaxNonZeros, s->mp->d_nnz, s->s.MaxNonZeros, s->mp->o_nnz);
        }

        // Set matrix values:
        for(int i = 0; i < s->mp->size; ++i) {
            assert(s->mp->col_ids[i] && s->mp->vals[i]);

            int row_id = i + s->first_gid;
            int num_cols = s->mp->d_nnz[i] + s->mp->o_nnz[i];
            MatSetValues(s->A, 1, &row_id, num_cols, s->mp->col_ids[i], s->mp->vals[i], INSERT_VALUES);
        }

        // Destroy prealloc and set vtable to not use it anymore.
        {mat_prealloc *mp = s->mp;
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
        s->mp = NULL;

        extern const _solver_vtable _petsc_vtable;
        s->s.vtable = &_petsc_vtable;

        // Finish setting up the matrix:
        {
            // Store cached last entries
            {int i = s->s.lasti; int numjs = s->s.lastj; const int *j = s->s.j; const real *v = s->s.vs;

                assert(i >= 0);
                for(unsigned k = 0; k < numjs; ++k) {
                    assert(j[k] >= 0);
                }

                // Preallocate before setting values
                if (s->s.MaxNonZeros == 0) {
                    slv_set_maxnonzeros(&s->s, 1000);
                }

                MatSetValues(s->A, 1, &i, numjs, j, v, INSERT_VALUES);
            }

            s->s.lasti = -1;

            // (s->lasti != -1) if, and only if, the cache vectors were allocated
            // Freeing them for while no matrix is being assembled.
            free(s->s.j);
            free(s->s.vs);
            s->s.j = NULL;
            s->s.vs = NULL;
            

            MatAssemblyBegin(s->A, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd  (s->A, MAT_FINAL_ASSEMBLY);


            KSPSetOperators(s->ksp, s->A, s->A);


            KSPSetTolerances(s->ksp, s->s.relative_tolerance, s->s.absolute_tolerance,
            PETSC_DEFAULT, s->s.MaxIteration);
        }
    }

}

void every_iteration (_solver_impl_petsc *s) {

    //KSPSetGuess
    //KSPSetReusePreconditioner
    //KSPSetDiagonalScale
    //SAME_NONZERO_PATTERN 
    //Reuse same matrix
    VecAssemblyBegin(s->b);
    VecAssemblyEnd  (s->b);

    VecAssemblyBegin(s->x);
    VecAssemblyEnd  (s->x);

    {
        MatAssemblyBegin(s->A, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd  (s->A, MAT_FINAL_ASSEMBLY);


        KSPSetOperators(s->ksp, s->A, s->A);


        KSPSetTolerances(s->ksp, s->s.relative_tolerance, s->s.absolute_tolerance,
        PETSC_DEFAULT, s->s.MaxIteration);
    }

    {
        KSPSolve(s->ksp, s->b, s->x);

        KSPGetIterationNumber(s->ksp, &s->s.SolverIterations);

        real abs_res_norm;
        KSPGetResidualNorm(s->ksp, &abs_res_norm);
        real rhs_norm;
        VecNorm(s->b, NORM_2, &rhs_norm);
        s->s.SolverResidual = abs_res_norm / rhs_norm;

        KSPConvergedReason r;
        KSPGetConvergedReason(s->ksp, &r);
        s->s.converged = r > 0;

        DEBUG_EXEC({
        KSPConvergedReason r;
        KSPGetConvergedReason(s->ksp, &r);
        print0f("PETSc DEBUG: Converged reason: %d, real abs norm: %g\n",
            r, calc_res_norm(s));
        });
    }



   
    {int size; const int *idx; real *x;
        VecGetValues(s->x, size, idx, x);
    }

   

}