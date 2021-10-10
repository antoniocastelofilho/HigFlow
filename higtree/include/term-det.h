#ifndef TERM_DET
#define TERM_DET

#include <stdbool.h>
#include <mpi.h>

/*
 * Implements weight throwing termination detection algorithm over MPI.
 * If function term_det_wait_for_msg() is always called with may_terminate == true,
 * the algorithm should never deadlock.
 */

struct TermDetection;
typedef struct TermDetection TermDetection;

TermDetection *term_det_init(MPI_Comm comm);
void term_det_destroy(TermDetection *td);

void * term_det_wait_for_msg(TermDetection *td, int *from, int *msg_size,
		bool may_terminate);

void term_det_send_async(TermDetection *td, const void *buf, size_t size, int to);
void term_det_send_asyncv(TermDetection *td, void const * const * bufs, const size_t *buf_sizes,
		int to);

#endif
