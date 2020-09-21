#include <mpi.h>
#include <stdio.h>
#include "solver-centralized.h"

static void
set_full_A(_solver_centralized *cs, unsigned int *row_jumper,
	unsigned int *col_buffer, real *coefs, size_t coefs_size)
{
	FILE *f = fopen("matrix-row-refs.dat", "wb");
	fwrite(row_jumper, sizeof(*row_jumper), cs->master->total_size+1, f);
	fclose(f);

	f = fopen("matrix-col-ids.dat", "wb");
	fwrite(col_buffer, sizeof(*col_buffer), coefs_size, f);
	fclose(f);

	f = fopen("matrix-elems.dat", "wb");
	fwrite(coefs, sizeof(*coefs), coefs_size, f);
	fclose(f);
}

static void
set_full_b(_solver_centralized *cs, real *rhs)
{
	FILE *f = fopen("matrix-rhs.dat", "wb");
	fwrite(rhs, sizeof(*rhs), cs->master->total_size, f);
	fclose(f);
}

static void solve(_solver_centralized *cs, real *solution)
{
	puts("Linear system printed to files. Exiting...");
	MPI_Abort(MPI_COMM_WORLD, 0);
}

static void destroy(_solver_centralized *cs)
{
}

static const _solver_centralized_vtable vtable = {
	.set_full_b = set_full_b,
	.set_full_A = set_full_A,
	.solve = solve,
	.destroy = destroy,
};

solver * _slv_create_debug_write(size_t start, size_t size)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	_solver_centralized *s;

	if(rank == 0) {
		s = malloc(sizeof(_solver_centralized));
		_solver_centralized_init_master(s, start, size, &vtable);
	} else {
		s = _solver_centralized_create_slave(start, size);
	}

	return &s->s;
}
