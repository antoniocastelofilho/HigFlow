#ifndef SOLVER_CENTRALIZED_H
#define SOLVER_CENTRALIZED_H

#include <stdbool.h>
#include "solver.h"

/* This file contains the internal interface used by solver engines that
 * doesn't work via MPI, so the it cares to centralize the data in a single
 * node, solve, and redistribute the result. */

struct _solver_centralized;

typedef struct _solver_centralized_vtable {
	void (*set_full_b)(struct _solver_centralized *cs, real *rhs);
	void (*set_full_A)(struct _solver_centralized *cs,
			unsigned int *row_jumper, unsigned int *col_buffer,
			real *coefs, size_t coefs_size);
	void (*solve)(struct _solver_centralized *cs, real *solution);
	void (*destroy)(struct _solver_centralized *cs);
} _solver_centralized_vtable;

typedef union _solver_centralized_split {
	struct{
		unsigned int start;
		unsigned int size;
	};
	unsigned long serial[2];
} _solver_centralized_split;

typedef struct _solver_centralized_master {
	unsigned int total_size;
	int *procs_order;
	real *solution;
	_solver_centralized_split all[];
} _solver_centralized_master;

typedef struct _solver_centralized {
	solver s;

	_solver_centralized_split local;

	size_t *A_sizes;
	unsigned int **A_cols_cpu;
	real **A_vals_cpu;
	real *B_cpu;
	real *X_cpu;

	_solver_centralized_master *master;

	const _solver_centralized_vtable *cvtable;
} _solver_centralized;

void _solver_centralized_init_master(_solver_centralized *slv, size_t start, size_t size, const _solver_centralized_vtable *internal_vtable);

_solver_centralized* _solver_centralized_create_slave(size_t start, size_t size);

#endif // Include guard
