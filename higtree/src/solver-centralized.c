/* This is not an implementation of a solver per-se, instead, this file is a
 * generic framework for implementing solvers that needs data do be
 * centralized on a single node (MPI node 0) before being solved. */

#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <assert.h>
#include <time.h>
#include "solver-centralized.h"

// TODO: refactor to use MPI_Scatterv, MPI_Gatherv...

enum LocallyUsedMPITags {
	INT_EXCHANGE_TAG = 42,
	REAL_EXCHANGE_TAG
};

static inline _solver_centralized* scast(solver *cs){
	return (_solver_centralized*)(cs);
}

// Assumes it will be called once, and at maximum once, before matrix elements are set
static void set_max_nonzeros(solver *cs)
{
	// Empty, there is no benefit in preallocating in this structure implemented here
}

static void set_Ai(solver *cs, int i, int numjs, const int *j, const real *v)
{
	_solver_centralized *s = scast(cs);
	i -= s->local.start;

	if(numjs > s->A_sizes[i]) {
		free(s->A_cols_cpu[i]);
		s->A_cols_cpu[i] = malloc(numjs * sizeof(unsigned int));

		free(s->A_vals_cpu[i]);
		s->A_vals_cpu[i] = malloc(numjs * sizeof(real));
	}
	s->A_sizes[i] = numjs;

	memcpy(s->A_cols_cpu[i], j, numjs * sizeof(unsigned int));
	memcpy(s->A_vals_cpu[i], v, numjs * sizeof(real));
}

static void set_bi(struct solver *cs, int i, real v)
{
	_solver_centralized *s = scast(cs);
	s->B_cpu[i - s->local.start] = v;
}

static void add_bi(struct solver *cs, int i, real v)
{
	_solver_centralized *s = scast(cs);
	s->B_cpu[i - s->local.start] += v;
}

static void set_xi(struct solver *cs, int i, real v)
{
	_solver_centralized *s = scast(cs);
	s->X_cpu[i - s->local.start] = v;
}

static void set_x(struct solver *cs, const real *x)
{
	_solver_centralized *s = scast(cs);
	memcpy(s->X_cpu, x, s->local.size * sizeof(real));
}

static real get_xi(struct solver *cs, int i)
{
	_solver_centralized *s = scast(cs);
	return s->X_cpu[i - s->local.start];
}

static void get_x(struct solver *cs, real *x)
{
	_solver_centralized *s = scast(cs);
	memcpy(x, s->X_cpu, s->local.size * sizeof(real));
}

static void get_x_scatter(struct solver *cs, int size, const int *idx, real *x)
{
	_solver_centralized *s = scast(cs);
	for(int i = 0; i < size; ++i)
	{
		x[i] = s->X_cpu[idx[i] - s->local.start];
	}
}

static void local_row_size_assemble(const _solver_centralized *s, unsigned int *row_jumper)
{
	int i;

	unsigned int counter = 0;
	for(i = 0; i < s->local.size; ++i) {
		row_jumper[i] = (counter += s->A_sizes[i]);
	}
}

static void local_row_assemble(const _solver_centralized *s, unsigned int *col_buffer, real *coefs)
{
	size_t counter = 0;
	for(int i = 0; i < s->local.size; ++i) {
		size_t size = s->A_sizes[i];
		memcpy(&col_buffer[counter], s->A_cols_cpu[i], size * sizeof(unsigned int));
		memcpy(&coefs[counter], s->A_vals_cpu[i], size * sizeof(real));
		counter += size;
	}
}

static inline MPI_Datatype get_mpi_real()
{
	assert(sizeof(real) == sizeof(float) || sizeof(real) == sizeof(double));
	return sizeof(real) == sizeof(double) ? MPI_DOUBLE : MPI_FLOAT;
}

//// Master Implementation
// The functions bellow will be called by the process that
// runs the actual solver

static void master_assemble_matrix(struct solver *cs)
{
	_solver_centralized *s = scast(cs);
	_solver_centralized_split *psplits = s->master->all;

	int num_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	unsigned int *row_jumper;
	unsigned int *col_buffer;
	real *coefs;

	row_jumper = malloc((s->master->total_size+1) * sizeof(unsigned int));
	row_jumper[0] = 0;

	// Get the start point for remote rows
	MPI_Request *reqs = malloc(2 * (num_procs-1) * sizeof(MPI_Request));
	for(int i = 1; i < num_procs; ++i) {
		_solver_centralized_split *psplit = &psplits[i-1];
		MPI_Irecv(&row_jumper[psplit->start+1], psplit->size, MPI_UNSIGNED, i, INT_EXCHANGE_TAG, MPI_COMM_WORLD, &reqs[i-1]);
	}

	// Calculate start point for local rows
	local_row_size_assemble(s, &row_jumper[s->local.start+1]);

	MPI_Waitall(num_procs-1, reqs, MPI_STATUSES_IGNORE);

	_solver_centralized_split proc_elems[num_procs];

	// Normalize starting points with global index
	size_t total_elems;
	{
		unsigned offset = 0;
		unsigned count = 0;
		const unsigned int num_jumpers = s->master->total_size + 1;
		int* procs_order = s->master->procs_order;

		proc_elems[procs_order[0]].start = 0;
		for(unsigned int i = 1; i < num_jumpers; ++i) {
			row_jumper[i] += offset;
			if(row_jumper[i] < row_jumper[i-1]) {
				unsigned int size = row_jumper[i-1] - offset;
				row_jumper[i] += size;
				proc_elems[procs_order[count]].size = size;
				offset = row_jumper[i-1];
				proc_elems[procs_order[++count]].start = offset;
			}

		}
		total_elems = row_jumper[num_jumpers-1];
		proc_elems[procs_order[count]].size = total_elems - offset;
		assert(count == num_procs-1);
	}

	// Receive the bulk data (coefficients and colums)
	col_buffer = malloc(total_elems * sizeof(unsigned int));
	coefs = malloc(total_elems * sizeof(real));
	for(int i = 1; i < num_procs; ++i) {
		size_t reqi = 2*(i-1);
		_solver_centralized_split *elems = &proc_elems[i];
		MPI_Irecv(&col_buffer[elems->start], elems->size, MPI_UNSIGNED, i, INT_EXCHANGE_TAG, MPI_COMM_WORLD, &reqs[reqi]);
		MPI_Irecv(&coefs[elems->start], elems->size, get_mpi_real(), i, REAL_EXCHANGE_TAG, MPI_COMM_WORLD, &reqs[reqi+1]);

	}

	// Assemble the local bulk data
	local_row_assemble(s, &col_buffer[proc_elems[0].start], &coefs[proc_elems[0].start]);

	MPI_Waitall(2 * (num_procs-1), reqs, MPI_STATUSES_IGNORE);
	free(reqs);

	s->cvtable->set_full_A(s, row_jumper, col_buffer, coefs, total_elems);

	free(coefs);
	free(col_buffer);
	free(row_jumper);
}

static void master_vector_assemble(_solver_centralized *s, real *vec, real *local_buf)
{
	int num_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	// Get remote
	MPI_Request reqs[num_procs-1];
	for(int i = 1; i < num_procs; ++i) {
		_solver_centralized_split *psplit = &s->master->all[i-1];
		MPI_Irecv(&vec[psplit->start], psplit->size, get_mpi_real(), i, REAL_EXCHANGE_TAG, MPI_COMM_WORLD, &reqs[i-1]);
	}

	// Copy local
	memcpy(&vec[s->local.start], local_buf, s->local.size * sizeof(real));

	// Wait remote copy to finish
	MPI_Waitall(num_procs-1, reqs, MPI_STATUSES_IGNORE);
}

static void master_assemble_rhs(struct solver *cs)
{
	_solver_centralized *s = scast(cs);

	real* rhs = malloc(s->master->total_size * sizeof(real));

	master_vector_assemble(s, rhs, s->B_cpu);
	s->cvtable->set_full_b(s, rhs);

	free(rhs);
}

static void master_assemble_output(struct solver *cs)
{
	_solver_centralized *s = scast(cs);

	master_vector_assemble(s, s->master->solution, s->X_cpu);
}

static void master_solve(struct solver *cs)
{
	_solver_centralized *s = scast(cs);
	real *sol = s->master->solution;

	// Solver the previously configured system
	s->cvtable->solve(s, sol);

	// Send the result back to the processes
	int num_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	MPI_Request reqs[num_procs];
	for(int i = 1; i < num_procs; ++i) {
		_solver_centralized_split *psplit = &s->master->all[i-1];
		MPI_Isend(&sol[psplit->start], psplit->size, get_mpi_real(), i, REAL_EXCHANGE_TAG, MPI_COMM_WORLD, &reqs[i-1]);
	}

	// Copy local
	memcpy(s->X_cpu, &sol[s->local.start], s->local.size * sizeof(real));

	// Wait remote copy to finish
	MPI_Waitall(num_procs-1, reqs, MPI_STATUSES_IGNORE);
}

//// Slave Implementation
// The functions bellow will be called by all other processes,
// sending the data to the master process, and waiting for the
// solution

static void slave_assemble_matrix(struct solver *cs)
{
	_solver_centralized *s = scast(cs);
	unsigned int *row_jumper = malloc(s->local.size * sizeof(unsigned int));

	local_row_size_assemble(s, row_jumper);
	MPI_Request req[3];
	MPI_Isend(row_jumper, s->local.size, MPI_UNSIGNED, 0, INT_EXCHANGE_TAG, MPI_COMM_WORLD, &req[0]);

	size_t elem_count = row_jumper[s->local.size-1];
	unsigned int *col_buffer = malloc(elem_count * sizeof(unsigned int));
	real *coefs = malloc(elem_count * sizeof(real));
	local_row_assemble(s, col_buffer, coefs);

	MPI_Isend(col_buffer, elem_count, MPI_UNSIGNED, 0, INT_EXCHANGE_TAG, MPI_COMM_WORLD, &req[1]);
	MPI_Isend(coefs, elem_count, get_mpi_real(), 0, REAL_EXCHANGE_TAG, MPI_COMM_WORLD, &req[2]);

	MPI_Waitall(3, req, MPI_STATUSES_IGNORE);

	free(coefs);
	free(col_buffer);
	free(row_jumper);
}

static void slave_assemble_rhs(struct solver *cs)
{
	_solver_centralized *s = scast(cs);
	MPI_Send(s->B_cpu, s->local.size, get_mpi_real(), 0, REAL_EXCHANGE_TAG, MPI_COMM_WORLD);
}

static void slave_assemble_output(struct solver *cs)
{
	_solver_centralized *s = scast(cs);
	MPI_Send(s->X_cpu, s->local.size, get_mpi_real(), 0, REAL_EXCHANGE_TAG, MPI_COMM_WORLD);
}

static void slave_solve(struct solver *cs)
{
	_solver_centralized *s = scast(cs);
	MPI_Request req;
	MPI_Irecv(s->X_cpu, s->local.size, get_mpi_real(), 0, REAL_EXCHANGE_TAG, MPI_COMM_WORLD, &req);

        const struct timespec sleep = {0, 100000000};
	int done;
	do {
		nanosleep(&sleep, NULL);
		MPI_Test(&req, &done, MPI_STATUS_IGNORE);
	} while(!done);
}

/// End of differing implementations

static void destroy(struct solver *cs)
{
	_solver_centralized *s = scast(cs);
	for(size_t i = 0; i < s->local.size; ++i) {
		free(s->A_cols_cpu[i]);
		free(s->A_vals_cpu[i]);
	}

	free(s->A_sizes);
	free(s->A_cols_cpu);
	free(s->A_vals_cpu);
	free(s->B_cpu);
	free(s->X_cpu);

	if(s->master) {
		free(s->master->procs_order);
		free(s->master->solution);
		free(s->master);
	}

	s->cvtable->destroy(s);
}

static const _solver_vtable master_vtable = {
	.set_maxnonzeros = set_max_nonzeros,
	.set_Ai = set_Ai,
	.set_bi = set_bi,
	.add_bi = add_bi,
	.set_xi = set_xi,
	.set_x = set_x,
	.get_xi = get_xi,
	.get_x = get_x,
	.get_x_scatter = get_x_scatter,
	.assemble_matrix = master_assemble_matrix,
	.assemble_rhs = master_assemble_rhs,
	.assemble_output = master_assemble_output,
	.solve = master_solve,
	.destroy = destroy,
};

static const _solver_vtable slave_vtable = {
	.set_maxnonzeros = set_max_nonzeros,
	.set_Ai = set_Ai,
	.set_bi = set_bi,
	.add_bi = add_bi,
	.set_xi = set_xi,
	.set_x = set_x,
	.get_xi = get_xi,
	.get_x = get_x,
	.get_x_scatter = get_x_scatter,
	.assemble_matrix = slave_assemble_matrix,
	.assemble_rhs = slave_assemble_rhs,
	.assemble_output = slave_assemble_output,
	.solve = slave_solve,
	.destroy = destroy,
};

static void _solver_centralized_init(_solver_centralized *slv, size_t start, size_t size, const _solver_centralized_vtable *internal_vtable)
{
	slv->A_sizes = calloc(size, sizeof(size_t));
	slv->A_cols_cpu = calloc(size, sizeof(unsigned int*));
	slv->A_vals_cpu = calloc(size, sizeof(real*));
	slv->B_cpu = calloc(size, sizeof(real));
	slv->X_cpu = calloc(size, sizeof(real));

	slv->local.size = size;
	slv->local.start = start;

	slv->cvtable = internal_vtable;
}

typedef struct ProcOrderer {
	unsigned int start;
	int proc_id;
} ProcOrderer;

int proc_compare(const ProcOrderer *a, const ProcOrderer *b)
{
	return (int)(a->start) - (int)(b->start);
}

void _solver_centralized_init_master(_solver_centralized *slv, size_t start, size_t size, const _solver_centralized_vtable *internal_vtable) {
	_solver_centralized_init(slv, start, size, internal_vtable);

	slv->s.vtable = &master_vtable;

	// Uses MPI to get to know other processes split:

	int num_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	_solver_centralized_master *master_data;
	master_data = slv->master = malloc(sizeof(_solver_centralized_master)
			+ (num_procs-1) * sizeof(_solver_centralized_split));

	master_data->procs_order = malloc(num_procs * sizeof(int));

	MPI_Request *reqs = malloc((num_procs-1) * sizeof(MPI_Request));
	for(int i = 1; i < num_procs; ++i) {
		MPI_Irecv(master_data->all[i-1].serial, 2, MPI_UNSIGNED,
				i, INT_EXCHANGE_TAG, MPI_COMM_WORLD, &reqs[i-1]);
	}
	MPI_Waitall(num_procs-1, reqs, MPI_STATUSES_IGNORE);
	free(reqs);

	unsigned long total_size = size;
	ProcOrderer order[num_procs];
	for(int i = 1; i < num_procs; ++i) {
		total_size += master_data->all[i-1].size;
		order[i].start = master_data->all[i-1].start;
		order[i].proc_id = i;
	}
	master_data->total_size = total_size;
	master_data->solution = malloc(total_size * sizeof(real));

	order[0].start = slv->local.start;
	order[0].proc_id = 0;
	qsort(order, num_procs, sizeof(ProcOrderer),
                  (int (*)(const void *, const void *))proc_compare);

	for(int i = 0; i < num_procs; ++i)
	{
		master_data->procs_order[i] = order[i].proc_id;
	}
}

static const _solver_centralized_vtable internal_slave_vtable = {
	.set_full_b = NULL,
	.set_full_A = NULL,
	.solve = NULL,
	.destroy = (void (*)(struct _solver_centralized *))(free),
};

_solver_centralized* _solver_centralized_create_slave(size_t start, size_t size)
{
	_solver_centralized* slv = malloc(sizeof(_solver_centralized));
	_solver_centralized_init(slv, start, size, &internal_slave_vtable);

	slv->s.vtable = &slave_vtable;
	slv->master = NULL;

	// Uses MPI to send metadata about this process' split to master
	MPI_Send(slv->local.serial, 2, MPI_UNSIGNED, 0, INT_EXCHANGE_TAG, MPI_COMM_WORLD);

	return slv;
}
