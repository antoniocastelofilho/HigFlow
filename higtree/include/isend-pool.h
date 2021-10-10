#ifndef ISEND_POOL
#define ISEND_POOL

#include <mpi.h>
#include <stdbool.h>

/* This convenience module allows for making multiple MPI_Isend's and waiting
 * for all to complete at once.
 */

struct isend_pool;
typedef struct isend_pool isend_pool;

/*! Create an Isend pool. */
isend_pool *isend_pool_create();

/*! Send a MPI_Isend() menssage, managed by this pool.
 *
 * If parameter @p copy_buf is false, then @p buf must be available
 * until isend_pool_wait_all() or isend_pool_destroy() is called.
 */
void isend_pool_send(isend_pool *isp, const void *buf, size_t count,
	MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
	bool copy_buf);

void isend_pool_wait_all(isend_pool *isp);

/*! Destroys the pool, waiting for all, if not waited already. */
void isend_pool_destroy(isend_pool *isp);

#endif
