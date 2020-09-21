#include <omp.h>
#include <numa.h>
#include <numaif.h>
#include <errno.h>
#include <assert.h>
#include <unistd.h>
#include <stdlib.h>
#include <mpi.h>

#include "solver-centralized.h"

typedef struct Mat {
	size_t size;

	unsigned int *col_buffer;
	unsigned int *row_jumper;
	real *elements;
} Mat;

typedef struct ImplSOR {
	_solver_centralized s;
	Mat A;
	real rhs_norm;
	real *rhs;
} ImplSOR;

static int get_numa_count()
{
	// Get the number of threads the program is set to use
	int omp_threads = omp_get_max_threads();

	// Get the number of NUMA banks the process can use
	int numa_nodes = numa_num_task_nodes();

	// Return the minimum of the the two
	if(omp_threads < numa_nodes)
		return omp_threads;
	else
		return numa_nodes;
}

// In a space of "total_size" elements, split it in "total_parts" parts, store in "split_size"
// the size of the chunk get by "me" part, and returns where this chunk starts.
static unsigned int
spread_evenly(unsigned int me, unsigned int total_parts, unsigned int total_size, unsigned int *split_size)
{
	*split_size = total_size / total_parts;
	unsigned int start = me * (*split_size);
	unsigned int rem = total_size % total_parts;
	if(me >= total_parts - rem) {
		++(*split_size);
		start += me - (total_parts - rem);
	}

	return start;
}

static void insist_moving_pages(void* first_page, size_t num_pages, size_t pagesize, int node, int threads)
{
	uint8_t *bytes = (uint8_t*) first_page;
	int *nodes = malloc(num_pages * sizeof(int) * 2);
	int *status = nodes + num_pages;
	void** pages = malloc(num_pages * sizeof(void*));

	#pragma omp parallel for num_threads(threads) schedule(static, 16)
	for(int i = 0; i < num_pages; ++i) {
		nodes[i] = node;
		pages[i] = bytes + (i * pagesize);
	}

	unsigned long pages_to_move = num_pages;
	int retries = 10;
	while(retries-- && pages_to_move) {
		move_pages(0, pages_to_move, (void **)pages, nodes, status, MPOL_MF_MOVE);

		unsigned long pending = 0;
		for(int i = 0; i < pages_to_move; ++i) {
			if(status[i] == -EBUSY) {
				pages[pending++] = pages[i];
			}
		}
		pages_to_move = pending;
	}

	free(nodes);
	free(pages);
}

inline static real
solver_sweep(
	const unsigned int *rows,
	const unsigned int *cols,
	const real *coefs,
	const real *rhs,
	real *abs_sol,
	real sor_omega,
	unsigned int size,
	unsigned int start)
{
	real compl_omega = 1.0 - sor_omega;
	real sum_residual = 0.0;

	#pragma omp for nowait schedule(dynamic, 920)
	for(int i = 0; i < size; ++i) {
		real new_val = rhs[i];
		real sol_i;
		real diagonal;
		real residual;
		unsigned int end = rows[i+1];
		unsigned int abs_i = i + start;

		for(int j = rows[i]; j < end; ++j) {
			unsigned int abs_col = cols[j]; // Absolute reference
			if(abs_col != abs_i) {
				real tmp;
				#pragma omp atomic read
				tmp = abs_sol[abs_col];
				new_val -= coefs[j] * tmp;
			} else {
				diagonal = coefs[j];
			}
		}

		sol_i = abs_sol[abs_i];
		residual = new_val - diagonal * sol_i;
		sum_residual += residual * residual;

		new_val = sor_omega * new_val / diagonal + compl_omega * sol_i;
		#pragma omp atomic write
		abs_sol[abs_i] = new_val;
	}

	return sum_residual;
}

// Assumes there are no null elements on main diagonal
static unsigned int solve(const Mat *matrix, const real *rhs, real *sol, real rhs_norm, real sor_omega, unsigned int max_iters, real abs_tol, real rel_tol, bool *converged)
{
	const size_t pagesize = sysconf(_SC_PAGESIZE);
	int numa_nodes = get_numa_count();
	unsigned int count;
	real *orig_sol = NULL;
	real norm_residual = 0.0;
	real abs_residual = 0.0;

	assert(pagesize % sizeof(real) == 0);
	const size_t elems_per_page = pagesize / sizeof(real);

	if(numa_nodes > 1) {
		// Test if sol is aligned to page size
		if(((uintptr_t)(const void *)sol) % pagesize != 0)
		{
			orig_sol = sol;
			posix_memalign((void **)&sol, pagesize, matrix->size * sizeof(real));
			memcpy(sol, orig_sol, matrix->size * sizeof(real));
		}

		omp_set_nested(true);
		// Split the data among the NUMA nodes
		#pragma omp parallel num_threads(numa_nodes) shared(norm_residual, count)
		{
			int node = omp_get_thread_num();

			// Calculate how many threads this NUMA node will use
			unsigned int node_threads;
		       	spread_evenly(node, numa_nodes, omp_get_max_threads(), &node_threads);

			// Split the work evenly among the nodes, so that the
			// biggest work size difference between them is 1 memory page of reals:

			// The number of pages the solution array uses
			size_t pages = 1 + ((matrix->size - 1) / elems_per_page);

			// Split the pages among the nodes
			unsigned int split_size, split_start;
			split_start = spread_evenly(node, numa_nodes, pages, &split_size);

			// Bind the thread to the node
			numa_run_on_node(node);

			// We have the work size in pages, use it to move the relevant part of
			// the sol vector to this node
			insist_moving_pages(sol + (split_start * elems_per_page),
					split_size, elems_per_page, node, node_threads);

			// Transform pages to number of elements
			split_start *= elems_per_page;
			split_size *= elems_per_page;
			printf("numa: %d, split_size: %u, split_start: %u\n", node, split_size, split_start);

			// Last page may be only partially filled
			if(node + 1 == numa_nodes) {
				size_t remaining = split_start + split_size - matrix->size;
				split_size -= remaining;
			}

			// Place the rest of the data on the relevant node:

			real *node_rhs = numa_alloc_local(split_size * sizeof(real));
			memcpy(node_rhs, rhs + split_start, split_size * sizeof(real));

			unsigned int *rows = numa_alloc_local((split_size+1) * sizeof(unsigned int));
			memcpy(rows, matrix->row_jumper + split_start, (split_size+1) * sizeof(unsigned int));

			unsigned int elems_count = rows[split_size] - rows[0];
			real *coefs = numa_alloc_local(elems_count * sizeof(real));
			memcpy(coefs, matrix->elements + rows[0], elems_count * sizeof(real));

			unsigned int *cols = numa_alloc_local(elems_count * sizeof(unsigned int));
			memcpy(cols, matrix->col_buffer + rows[0], elems_count * sizeof(unsigned int));

			// Adjust the row jumper to the new referential
			if(rows[0]) {
				#pragma omp parallel for num_threads(node_threads) schedule(static, 16)
				for(int i = 1; i <= split_size; ++i)
				{
					rows[i] -= rows[0];
				}
				rows[0] = 0;
			}

			// The solver loop
			unsigned int node_count = 0;
			real node_residual;
			do {
				node_residual = 0.0;
				#pragma omp barrier
				#pragma omp single
				norm_residual = 0.0;
				#pragma omp parallel reduction(+:node_residual) num_threads(node_threads)
				{
					numa_run_on_node(node);
			       		node_residual = solver_sweep(
						rows, cols, coefs, node_rhs,
						sol, sor_omega,	split_size,
						split_start
					);
				}

				#pragma omp atomic
				norm_residual += node_residual;

				#pragma omp barrier
				abs_residual = sqrt(norm_residual);
				node_residual = abs_residual / rhs_norm;
			} while(++node_count < max_iters && node_residual > rel_tol && abs_residual > abs_tol);

			numa_free(cols, elems_count * sizeof(unsigned int));
			numa_free(coefs, elems_count * sizeof(real));
			numa_free(rows, (split_size+1) * sizeof(unsigned int));
			numa_free(node_rhs, split_size * sizeof(real));

			#pragma omp single nowait
			{
				count = node_count;
				norm_residual = node_residual;
			}
		}

		if(orig_sol)
		{
			memcpy(orig_sol, sol, matrix->size * sizeof(real));
			free(sol);
		}
	} else {
		real *coefs = matrix->elements;
		unsigned int *cols = matrix->col_buffer;
		unsigned int *rows = matrix->row_jumper;
		int size = matrix->size;
		int num_threads = omp_get_max_threads();
		count = 0;
		real abs_residual;

		do {
			norm_residual = 0.0;
			#pragma omp parallel reduction(+:norm_residual)
			{
				norm_residual = solver_sweep(
					rows, cols, coefs, rhs,
					sol, sor_omega, size, 0
				);
			}

			abs_residual = sqrt(norm_residual);
			norm_residual = abs_residual / rhs_norm;
		} while(++count < max_iters && norm_residual > rel_tol && abs_residual > abs_tol);
	}

	*converged = norm_residual < rel_tol || abs_residual < abs_tol;
	return count;
}

static inline ImplSOR * tosor(_solver_centralized *cs)
{
	return (ImplSOR *)cs;
}

static void
set_full_A(_solver_centralized *cs, unsigned int *row_jumper,
	unsigned int *col_buffer, real *coefs, size_t coefs_size)
{
	ImplSOR *s = tosor(cs);

	size_t old_size = s->A.row_jumper[s->A.size];
	size_t new_size = row_jumper[s->A.size];

	memcpy(s->A.row_jumper, row_jumper, (s->A.size+1) * sizeof(unsigned int));
	if(new_size > old_size) {
		free(s->A.col_buffer);
		s->A.col_buffer = NULL;
		free(s->A.elements);
		s->A.elements = NULL;
	}

	s->A.col_buffer = realloc(s->A.col_buffer, new_size * sizeof(unsigned int));
	s->A.elements = realloc(s->A.elements, new_size * sizeof(real));

	memcpy(s->A.col_buffer, col_buffer, new_size * sizeof(unsigned int));
	memcpy(s->A.elements, coefs, new_size * sizeof(real));
}

static void
set_full_b(_solver_centralized *cs, real *rhs)
{
	ImplSOR *s = tosor(cs);
	memcpy(s->rhs, rhs, s->A.size);

	real sq_sum = 0.0;
	#pragma omp parallel for reduction(+:sq_sum) schedule(static, 16)
	for(size_t i = 0; i < s->A.size; ++i) {
		sq_sum += rhs[i] * rhs[i];
	}
	s->rhs_norm = sqrt(sq_sum);
}

static void
solve_interface(_solver_centralized *cs, real *solution)
{
	ImplSOR *s = tosor(cs);
	cs->s.SolverIterations = solve(&s->A, s->rhs, solution, s->rhs_norm, 1.99, cs->s.MaxIteration,
		cs->s.absolute_tolerance, cs->s.relative_tolerance, &cs->s.converged);
}

static void destroy(_solver_centralized *cs)
{
	ImplSOR *s = tosor(cs);
	free(s->rhs);
	free(s->A.row_jumper);
	free(s->A.col_buffer);
	free(s->A.elements);
	free(s);
}

const _solver_centralized_vtable vtable = {
	.set_full_b = set_full_b,
	.set_full_A = set_full_A,
	.solve = solve_interface,
	.destroy = destroy,
};

solver *_slv_create_sor(size_t start, size_t size)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	solver *ret;

	if(rank == 0) {
		ImplSOR *s = malloc(sizeof(ImplSOR));
		_solver_centralized_init_master(&s->s, start, size, &vtable);

		s->A.size = s->s.master->total_size;
		s->rhs = malloc(s->A.size * sizeof(real));
		s->A.row_jumper = calloc(s->A.size + 1,  sizeof(unsigned int));
		s->A.col_buffer = NULL;
		s->A.elements = NULL;

		ret = &s->s.s;
	} else {
		_solver_centralized *s = _solver_centralized_create_slave(start, size);
		ret = &s->s;
	}

	return ret;
}
