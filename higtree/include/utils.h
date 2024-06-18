#ifndef __UTILS_H
#define __UTILS_H

#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<stdbool.h>
#include<math.h>
#define __USE_POSIX199309 1
#include<time.h>

#include "coord.h"

// Smaller than the smallest delta
//! Defines the smallest delta of a cell.
#define EPSDELTA 1.0e-8
//TODO: There is a mismatch between parent boundaries and children's which is of order e-7.

#define POS_EQ(x, y) (fabs((x) - (y)) < EPSDELTA)
#define POS_NE(x, y) (fabs((x) - (y)) >= EPSDELTA)
#define POS_GE(x, y) ((x) > ((y) - EPSDELTA))
#define POS_LE(x, y) ((x) < ((y) + EPSDELTA))

//! Defines the smallest difference between two reals; anything smaller than this is considered as rounding error.
#define EPSMACH 1.0e-12

//! Checkes whether x is greater than y, considering the rounding error.
#define FLT_GT(x, y) ((x) > ((y) + EPSMACH))

//! Checkes whether x is greater than or equal to y, considering the rounding error.
#define FLT_GE(x, y) ((x) > ((y) - EPSMACH))

//! Checkes whether x is lesser than y, considering the rounding error.
#define FLT_LT(x, y) ((x) < ((y) - EPSMACH))

//! Checkes whether x is lesser than or equal to y, considering the rounding error.
#define FLT_LE(x, y) ((x) < ((y) + EPSMACH))

//! Checkes whether x is equals to y, considering the rounding error.
//#define FLT_EQ(x, y) ((FLT_GE(x,y))&&(FLT_LE(x,y)))
#define FLT_EQ(x, y) (fabs((x) - (y)) < EPSMACH)

//! Checkes whether x is different from y, considering the rounding error.
#define FLT_NE(x, y) (!FLT_EQ(x,y))

#define ALLOCATION_ERROR_CHECK(type, v, size) do { \
	if (v == NULL) { \
		fprintf(stderr, \
			"Error allocating at %s:%d, variable %s * %s (size %ld)", \
			__FILE__, __LINE__, #type, #v, (long)size); \
		abort(); \
	} } while(0)

//! Allocs size elements of type in variable v.
// Aborts the program if the allocation fails.
#define ALLOC(type, v, size) do { v = (type *) malloc(sizeof(type)*(size)); \
	ALLOCATION_ERROR_CHECK(type, v, size); } while(0)


//! Allocate an element of type in variable v, with count elements in
// the flexible array member. Aborts the program if the allocation fails.
#define ALLOC_INFER_FLEX(v, flex_member, count) do { \
		v = malloc(count * (sizeof v->flex_member[0]) + sizeof *v); \
		ALLOCATION_ERROR_CHECK(type, v, 1); \
	} while(0)

//! Declares the variable v of type*, and allocs size elements.
// Aborts the program if the allocation fails.
#define DECL_AND_ALLOC(type, v, size) type * v; \
	ALLOC(type, v, size)

//! Declares the variable v of type*, and allocates it with count elements in
// the flexible array member. Aborts the program if the allocation fails.
#define DECL_AND_ALLOC_FLEX(type, v, flex_member, count) type * v; \
	ALLOC_INFER_FLEX(v, flex_member, count)

//! \brief Like ALLOC, but uses type of v
#define ALLOC_INFER(v, size) do { (v) = malloc((size) * sizeof *(v)); if (v == NULL) {fprintf(stderr, "Error allocating at %s:%d, variable %s (size %ld)", __FILE__, __LINE__, #v, (long)size); abort();} } while(0)

//! \brief Realloc to multiples of the type and check for error
#define REALLOC_INFER(v, size) do { (v) = realloc((v), (size) * sizeof *(v)); if ((size) && !(v)) {fprintf(stderr, "Error reallocating at %s:%d, variable %s (new size %ld)", __FILE__, __LINE__, #v, (long)(size)); abort();} } while(0)


//! \brief Declares the variable v of type, and initializes it with val.
#define DECL_AND_INIT(type, v, val) type v = (type) val;

/** Advance pp one step inside a subgrid so that lo <= pp < hi.
 * All parameters pp, lo and hi are assumed to be arrays of
 * dimension DIM. The positions are advanced in the same order
 * as stored in hig_cell.children.
 */
bool subgrid_coord_advance(int pp[DIM], const int lo[DIM], const int hi[DIM]);

/** Integer exponentiation. */
unsigned upow(unsigned base, unsigned exp) __attribute__ ((const));

int is_within(int pp[DIM], int l[DIM], int h[DIM]);

/** Compare function suitable for qsort on (unsigned) int. */
int int_cmp(const int *a, const int *b);

real max(real a, real b) __attribute__ ((const));
real min(real a, real b) __attribute__ ((const));

//clock_gettime is not implemented on OSX
#ifdef __MACH__
#include <time.h>
#define CLOCK_MONOTONIC 0
#define CLOCK_REALTIME 0
int clock_gettime(int clk_id, struct timespec* t);
#endif

#define DECL_CLOCK(cl) struct timespec (__start_clock_ ## cl) , (__end_clock_ ## cl) ; double __time_past_ ## cl = 0.0;

#define START_CLOCK(cl) { clock_gettime(CLOCK_MONOTONIC, &(__start_clock_ ## cl)); }

#define STOP_CLOCK(cl) { clock_gettime(CLOCK_MONOTONIC, &(__end_clock_ ## cl)); (__time_past_ ## cl) += ((__end_clock_ ## cl).tv_sec - (__start_clock_ ## cl).tv_sec) * 1000000000.0 + ((__end_clock_ ## cl).tv_nsec - (__start_clock_ ## cl).tv_nsec); }

#define GET_NSEC_CLOCK(cl) (__time_past_ ## cl)

/*! printf restricted to MPI process 0. */
int print0f(const char *format, ...);

/*! Initialize underlining libraries, like MPI and PETSc.
 *
 * Must be called first thing in the program.
 */
void higtree_initialize(int *argc, char **argv[]);

//#define DECL_CLOCK(cl)
//#define START_CLOCK(cl)
//#define STOP_CLOCK(cl)
//#define GET_NSEC_CLOCK(cl)



#endif
