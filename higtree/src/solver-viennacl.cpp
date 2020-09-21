#define VIENNACL_WITH_OPENCL

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/vector.hpp>
#include <viennacl/linalg/bicgstab.hpp>
#include <viennacl/linalg/amg.hpp>
#include <utility>
#include <climits>
#include <algorithm>
#include <iostream>
#include <mpi.h>

extern "C" {
#include "solver-centralized.h"
solver *_slv_create_viennacl(size_t start, size_t size);
}

typedef viennacl::compressed_matrix<real> GPUMatrix;
typedef viennacl::vector<real> GPUVector;
typedef viennacl::linalg::amg_precond<GPUMatrix> PrecondType;

struct ImplVcl: public _solver_centralized {
	ImplVcl(size_t _start, size_t _size)
	{
		_solver_centralized_init_master(this, _start, _size, &master_vtable);
		b.resize(master->total_size);
	}

	GPUVector b;
	GPUMatrix A;

	std::unique_ptr<PrecondType> preconditioner;
	std::unique_ptr<viennacl::linalg::bicgstab_tag> solver_tag;

	static const _solver_centralized_vtable master_vtable;
};

static ImplVcl* to_vcl(_solver_centralized* cs)
{
	return static_cast<ImplVcl*>(cs);
}

static void set_full_A(_solver_centralized *cs, unsigned int *row_jumper,
			unsigned int *col_buffer, real *coefs, size_t coefs_size)
{
	auto *s = to_vcl(cs);

	viennacl::vcl_size_t size = s->master->total_size;
	s->A.set(row_jumper, col_buffer, coefs, size, size, coefs_size);

	/*std::ofstream file("saida.txt");
	for(size_t i = 0; i < size; ++i) {
		for(size_t j = row_jumper[i]; j < row_jumper[i+1]; ++j) {
			file << col_buffer[j] << ' ' << coefs[j] << ' ';
		}
		file << '\n';
	}
	file << '\n';
	for(size_t i = 0; i < size; ++i) {
		char col[size];
		memset(col, '.', size);
		for(size_t j = row_jumper[i]; j < row_jumper[i+1]; ++j) {
			col[col_buffer[j]] = '*';
		}
		for(size_t j = 0; j < size; ++j) {
			file << col[j] << ' ';
		}
		file << '\n';
	}*/


	// Setup the solver settings.
	/*s->preconditioner.reset(
		new PrecondType(s->A, viennacl::linalg::amg_tag())
	);*/

	s->solver_tag.reset(
		new viennacl::linalg::bicgstab_tag(s->s.relative_tolerance, s->s.MaxIteration)
	);
	s->solver_tag->abs_tolerance(s->s.absolute_tolerance);
}

static void set_full_b(_solver_centralized *cs, real *rhs)
{
	auto *s = to_vcl(cs);
	viennacl::fast_copy(rhs, rhs + s->master->total_size, s->b.begin());
}

static void solve(_solver_centralized *cs, real *solution)
{
	auto *s = to_vcl(cs);

	// TODO: use solution as initial guess
	// see: https://github.com/viennacl/viennacl-dev/issues/97

	//GPUVector x(s->master->total_size);
	//viennacl::fast_copy(solution, solution + s->master->total_size, x.begin());
	std::cout << "Solving with " << s->b.size() << " elements." << std::endl;

	GPUVector x = viennacl::linalg::solve(s->A, s->b, *s->solver_tag);
	std::cout << "Iters: " << s->solver_tag->iters()
		<< "\nError: " << s->solver_tag->error() << '\n';

	viennacl::fast_copy(x.begin(), x.end(), solution);

	cs->s.SolverIterations = s->solver_tag->iters();
	cs->s.SolverResidual = viennacl::linalg::norm_2(
			viennacl::linalg::prod(s->A, x) - s->b
		) / viennacl::linalg::norm_2(s->b);
	cs->s.converged = s->solver_tag->iters() >= cs->s.MaxIteration;

	// TODO: investigate wrong residual
	/*x = s->b - viennacl::linalg::prod(s->A, x);
	std::vector<real> res(x.size());
	viennacl::fast_copy(x.begin(), x.end(), res.begin());
	real mean_res = 0.0;
	for(auto r : res) {
		mean_res += fabs(r);
	}
	mean_res /= res.size();
	std::cout << "Calculated mean residual: " << mean_res << std::endl;*/
}

static void destroy(_solver_centralized *cs)
{
	auto *s = to_vcl(cs);
	delete s;
}

const _solver_centralized_vtable ImplVcl::master_vtable = {
	.set_full_b = set_full_b,
	.set_full_A = set_full_A,
	.solve = solve,
	.destroy = destroy,
};

solver *_slv_create_viennacl(size_t start, size_t size)
{
	static_assert(sizeof(unsigned int) == 4,
		"unsigned int is not 32 bits, code will be incompatible with OpenCL");

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	solver *ret;

	if(rank == 0) {
		ImplVcl *s = new ImplVcl(start, size);
		ret = &s->s;
	} else {
		_solver_centralized *s = _solver_centralized_create_slave(start, size);
		ret = &s->s;
	}

	return ret;
}
