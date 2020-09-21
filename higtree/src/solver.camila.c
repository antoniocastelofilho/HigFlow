#include <stdlib.h>
#include <string.h>
#include "solver.h"
#include "Debug-c.h"

size_t slv_get_local_size(solver *s)
{
    return s->size;
}

void slv_set_maxiteration( solver *s, int _maxit )
{
    s->MaxIteration = _maxit;
}

void slv_set_maxnonzeros( solver *s, int _maxnz ) {
    s->MaxNonZeros = _maxnz;
    s->vtable->set_maxnonzeros(s);
}

void slv_set_eps( solver *s, real _eps )
{
    s->absolute_tolerance = _eps;
    s->relative_tolerance = _eps;
}

void slv_set_abs_tol( solver *s, real tol)
{
    s->absolute_tolerance = tol;
}

void slv_set_rel_tol( solver *s, real tol)
{
    s->relative_tolerance = tol;
}

/** Add lines directly to the matrix structure */
void slv_set_Ai(solver *s, int i, int numjs, int *j, real *v)
{
    assert(i >= 0);
    for(unsigned k = 0; k < numjs; ++k) {
	    assert(j[k] >= 0);
    }

    // Preallocate before setting values
    if (s->MaxNonZeros == 0) {
	slv_set_maxnonzeros(s, 1000);
    }

    s->vtable->set_Ai(s, i, numjs, j, v);
}

/** Add values to the matrix, caching lines */
void slv_set_Aij(solver *s, int i, int j, real v) {

    // Preallocate before setting values
    if (s->MaxNonZeros == 0) {
	slv_set_maxnonzeros(s, 1000);
    }

    if (i == s->lasti) {

        // Cache values to be added at once in the matrix
        s->j [s->lastj] = j;
        s->vs[s->lastj] = v;
        s->lastj++;
    }

    else {

        if (s->lasti != -1) {
	    // Refrain from setting an imposed value
            if(s->imposed_line != s->lasti) {
                // Add cached values to the matrix
                s->vtable->set_Ai(s, s->lasti, s->lastj, s->j, s->vs);
            }

            // Reset pointers
            for (int c = 0; c < s->lastj; c++) {
                s->j[c]  = 0;
                s->vs[c] = 0.0;
            }
        } else {
            // The cache is not allocated
            // allocate it!
            ALLOC(int, s->j, s->size);
            ALLOC(real, s->vs, s->size);
        }

        // Now cache new value
        s->lasti = i;
        s->j [0] = j;
        s->vs[0] = v;
        s->lastj = 1;
    }
}

/** Set RHS */
void slv_set_bi(solver *s, int i, real v)
{
    assert(i >= 0);

    if(s->imposed_line != i)
        s->vtable->set_bi(s, i, v);
}

/** Add values to RHS */
void slv_add_bi(solver *s, int i, real v)
{
    assert(i >= 0);

    if(s->imposed_line != i)
        s->vtable->add_bi(s, i, v);
}

/** Set ith entry of the initial guess */
void slv_set_xi(solver *s, int i, real v)
{
    s->vtable->set_xi(s, i, v);
}

/** Set entire initial guess */
void slv_set_x(solver *s, real *x)
{
    s->vtable->set_x(s, x);
}

/** Get ith entry of the soluction vector */
real slv_get_xi(solver *s, int i)
{
    return s->vtable->get_xi(s, i);
}

/** Get entire solution vector */
void slv_get_x(solver *s, real *x) {
    s->vtable->get_x(s, x);
}

/** Impose value for one unknown. Use before assembling matrix rows and right hand side.
 * TODO: Maybe implements slv_unimpose_value ?
 */
void slv_impose_value(solver *s, int line, real value) {
    s->imposed_line = line;
    s->vtable->set_bi(s, line, value);
    const real one = 1.0;
    s->vtable->set_Ai(s, line, -1, &line, &one);

}

/** Get and arbitrary chunk from the solution vector */
void slv_get_x_scatter(solver *s, int size, const int *idx, real *x) {
    s->vtable->get_x_scatter(s, size, idx, x);
}

/** Assemble matrix for use in solver.
 *
 * The matrix can be assembled only once while
 * the tree structure is unchanged.
 */
void slv_assemble_matrix(solver *s) {
    if (s->lasti != -1) {
        // Store cached last entries
        s->vtable->set_Ai(s, s->lasti, s->lastj, s->j, s->vs);
        s->lasti = -1;

	// (s->lasti != -1) if, and only if, the cache vectors were allocated
	// Freeing them for while no matrix is being assembled.
	free(s->j);
	free(s->vs);
	s->j = NULL;
       	s->vs = NULL;
    }

    s->vtable->assemble_matrix(s);
}


/** Assemble the right hand side vector.
 *
 * Will possibly change at every time step.
 */
void slv_assemble_rhs(solver *s) {
    s->vtable->assemble_rhs(s);
}


/** Assemble the right hand side vector.
 *
 * Will possibly change at every time step.
 */
void slv_assemble_output(solver *s) {
    s->vtable->assemble_output(s);
}


/** Assemble matrix and vectors before solving  */
void slv_assemble(solver *s) {
    slv_assemble_rhs(s);
    slv_assemble_output(s);
    slv_assemble_matrix(s);
}


/** Solve current linear system */
void slv_solve(solver *s) {
    s->vtable->solve(s);
}

bool slv_has_converged(solver *s)
{
	return s->converged;
}

int slv_get_iteration_count(solver *s)
{
    return s->SolverIterations;
}

real slv_get_relative_residual_norm(solver *s)
{
    return s->SolverResidual;
}

/** Destroy structure */
void slv_destroy(solver *s) {

    s->vtable->destroy(s);

    //// free cache
    free (s->j);
    free (s->vs);

    free(s);
}

solver * slv_create(SolverEngine engine, int first_gid, size_t size)
{
    extern solver* _slv_create_petsc(int first_gid, int _size);
    extern solver * _slv_create_hollow_setup_petsc(int first_gid, int _size);
    extern solver* _slv_create_viennacl(size_t start, size_t size);
    extern solver* _slv_create_hypre(int ilower, int iupper);
    extern solver* _slv_create_sor(size_t start, size_t size);
    extern solver* _slv_create_debug_write(size_t start, size_t size);

    // First try to use solver specified in environment variable
    if(engine == SOLVER_ANY) {
	const char* env_val = getenv("HIGTREE_SOLVER_ENGINE");
	if(env_val) {
		for(int i = 0; i < (sizeof(engine_names)/sizeof(char*)); ++i) {
			if(strcmp(env_val, engine_names[i]) == 0) {
				engine = i+1;
				break;
			}
		}
	}
    }

    // If still undefined, choose according some heuristic
    if(engine == SOLVER_ANY) {
        /*if(MPI_Comm_size() == 1) {
            // TODO: Find a good reason to choose ViennaCL
            engine = SOLVER_VIENNACL;
        } else if(sizeof(double) == sizeof(real)) {
	    // Very fast, but unstable...
	    // TODO: Is speed more important than stability?
            engine = SOLVER_HYPRE;
        } else */ {
            // Assuming linked PETSC has same precision as HigTree
            engine = SOLVER_PETSC;
        }
    }

    print0f("Solver engine used: %s\n", engine_names[engine-1]);

    solver *s = NULL;

    switch(engine) {
        case SOLVER_HYPRE:
            if(sizeof(double) != sizeof(real)) {
                print0f("Warning: using linear sover HYPRE with double precision,\n"
			"  but HigTree was built with a different precision.\n"
			"  This is suboptimal: an extra data conversion step is used.");
            }
            s = _slv_create_hypre(first_gid, first_gid + size - 1);
            break;
        case SOLVER_PETSC:
            s = _slv_create_petsc(first_gid, size);
            break;
	case SOLVER_HOLLOW:
	    s = _slv_create_hollow_setup_petsc(first_gid, size);
	    break;
#ifdef USE_SOR
	case SOLVER_SOR:
	    s = _slv_create_sor(first_gid, size);
	    break;
#endif
	case SOLVER_DEBUG_WRITE:
	    s = _slv_create_debug_write(first_gid, size);
	    break;
#ifdef USE_VIENNACL
        case SOLVER_VIENNACL:
            s = _slv_create_viennacl(first_gid, size);
            break;
#endif
	default:
	    ; // Must never happen.
    }

    // Generic stuff
    s->size         = size;
    s->MaxIteration = 10000;
    s->MaxNonZeros  = 0;
    s->absolute_tolerance = 1e-8;
    s->relative_tolerance = 1e-8;

    // Cache stuff
    s->lasti = -1;
    s->lastj = 0;
    s->j = NULL;
    s->vs = NULL;

    // No imposed value:
    s->imposed_line = -1;

    return s;
}
