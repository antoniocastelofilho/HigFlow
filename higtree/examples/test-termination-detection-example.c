#include <stdlib.h>
#include <stdio.h>
#include "term-det.h"

// This test will stress out the termination detection algorithm. By sending
// hundreds of randomly trigered messages to random targets.

/** Generate a random numper in the range [from, to) */
static int
rand_range(int from, int to)
{
	int r = rand();
	r = from + (int)((((double)r) / RAND_MAX) * (to - from));
	return r == to ? r - 1 : r;
}

unsigned sent_msgs = 0;

static void
send_random_message(TermDetection *td, int nprocs, int msg)
{
	int dest = rand_range(0, nprocs);
	term_det_send_async(td, &msg, sizeof msg, dest);
	++sent_msgs;
}

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);

	int rank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	srand(rank);

	if(rank == 0) {
		puts("Running termination detection test. If wrong, this may deadlock.");
	}

	TermDetection *td = term_det_init(MPI_COMM_WORLD);

	for(unsigned i = 0; i < 500; ++i) {
		send_random_message(td, nprocs, 1);
	}

	int max_round = 0;
	int *ret;
	int from, msg_size;
	unsigned recved_msgs = 0;
	while((ret = term_det_wait_for_msg(td, &from, &msg_size, true))) {
		++recved_msgs;
		if(max_round < *ret) {
			max_round = *ret;
		}

		/* 50% of chance of sending yet another message. */
		if(rand() & 1u)
			send_random_message(td, nprocs, 1 + *ret);
	}
	printf("Proc %d reached round %d.\n", rank, max_round);

	term_det_destroy(td);

	if(rank == 0) {
		unsigned total_recved_msgs;
		MPI_Reduce(&recved_msgs, &total_recved_msgs, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

		unsigned total_sent_msgs;
		MPI_Reduce(&sent_msgs, &total_sent_msgs, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

		MPI_Finalize();

		printf("Total sent messages: %u\nTotal received messages: %u\n", total_sent_msgs, total_recved_msgs);
		return total_sent_msgs != total_recved_msgs;
	} else {
		MPI_Reduce(&recved_msgs, NULL, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&sent_msgs, NULL, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

		MPI_Finalize();
	}
}
