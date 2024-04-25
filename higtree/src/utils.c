#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdbool.h>
#include <stdarg.h>

#define __USE_POSIX199309 1
#include <time.h>
#include <sys/resource.h>
#include <petsc.h>
#include <mpi.h>

#include "Debug-c.h"
#include "utils.h"

// Poiters to program arguments, used
// to control PETSc initialization.
static int _argc;
static char **_argv;

static bool is_petsc_initialized = false;

/* Must be called every time PETSc is used and it
 * is not know if it was already called.
 * Currently, this is called when creating a PETSc solver object. */
void _try_initialize_petsc()
{
	if(is_petsc_initialized)
		return;

	is_petsc_initialized = true;

	char **orig_argv = _argv;
	char *orig_argv_buf = *_argv;

	PetscInitialize(&_argc, &_argv, NULL, NULL);

	_argv = orig_argv;
	*_argv = orig_argv_buf;
}

static void _higtree_finalize()
{
	if(is_petsc_initialized) {
		PetscFinalize();
	}

	free(*_argv);
	free(_argv);

	MPI_Finalize();
}

void higtree_initialize(int *argc, char **argv[])
{
	MPI_Init(argc, argv);

	// Copy args contents because user may modify
	// it before we need to use it.
	_argc = *argc;
	ALLOC_INFER(_argv, _argc);
	size_t argv_size = 0;
	for(int i = 0; i < _argc; ++i) {
		argv_size += strlen((*argv)[i]) + 1;
	}
	ALLOC_INFER(*_argv, argv_size);
	for(int i = 0;;) {
		char *next = stpcpy(_argv[i], (*argv)[i]);

		if(++i >= _argc)
			break;

		_argv[i] = next + 1;
	}

	// Tries setting the stack size to unlimited.
	struct rlimit rlim =
	{
		.rlim_cur = RLIM_INFINITY,
		.rlim_max = RLIM_INFINITY,
	};
	if(setrlimit(RLIMIT_STACK, &rlim) < 0) {
		char *err = strerror(errno);
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		printf("Warning! Could not raise the stack limit on process %d:\n"
			"  %s.\n"
			"  Either raise the limit yourself, or expect segmentation fault.\n",
			rank, err);
	}

	atexit(_higtree_finalize);
}

bool subgrid_coord_advance(int pp[DIM], const int lo[DIM], const int hi[DIM])
{
	for(int dim = DIM-1; dim >= 0; --dim) {
		if(++pp[dim] < hi[dim]) {
			return true;
		} else {
			pp[dim] = lo[dim];
		}
	}
	return false;
}

unsigned upow(unsigned base, unsigned exp)
{
	/* Algorithm inspired by: http://stackoverflow.com/a/101613/578749 */
	unsigned result = 1;
	for(;;)
	{
		if (exp & 1)
			result *= base;
		exp >>= 1;

		if(exp) {
			base *= base;
		} else {
			return result;
		}
	}
}

real max(real a, real b)
{
	return a > b ? a : b;
}

real min(real a, real b)
{
	return a < b ? a : b;
}

int is_within(int pp[DIM], int l[DIM], int h[DIM]) {
	for(int i = 0; i < DIM; i++) {
		if (pp[i] < l[i] || pp[i] >= h[i]) {
			return 0;
		}
	}
	return 1;
}

int int_cmp(const int *a, const int *b)
{
	return *a - *b;
}

#ifdef __MACH__
#include <time.h>
//clock_gettime is not implemented on OSX
int clock_gettime(int clk_id, struct timespec* t) {
    struct timeval now;
    int rv = gettimeofday(&now, NULL);
    if (rv) return rv;
    t->tv_sec  = now.tv_sec;
    t->tv_nsec = now.tv_usec * 1000;
    return 0;
}
#endif

int print0f(const char *format, ...)
{
	va_list ap;
	int ret = 0;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(rank == 0) {
		va_start(ap, format);
		ret = vprintf(format, ap);
		fflush(stdout);
		va_end(ap);
	}
	return ret;
}
