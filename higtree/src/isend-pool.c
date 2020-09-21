#include <string.h>
#include "utils.h"

#include "isend-pool.h"

// Request data will be stored in a linked list.
typedef struct request_node
{
	struct request_node *next;
	MPI_Request req;
	char buf[];
} request_node;

struct isend_pool
{
	request_node *reqs;
};

isend_pool *isend_pool_create()
{
	DECL_AND_ALLOC(isend_pool, ret, 1);
	ret->reqs = NULL;
	return ret;
}

void isend_pool_send(isend_pool *isp, const void *buf, size_t count,
	MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
	bool copy_buf)
{
	size_t extra_size = 0;
	if(copy_buf) {
		int dsize;
		MPI_Type_size(datatype, &dsize);
		extra_size = count * dsize;
	}
	DECL_AND_ALLOC_FLEX(request_node, req, buf, extra_size);

	req->next = isp->reqs;
	isp->reqs = req;

	if(copy_buf) {
		memcpy(req->buf, buf, extra_size);
		buf = req->buf;
	}

	MPI_Isend(buf, count, datatype, dest, tag, comm, &req->req);
}

void isend_pool_wait_all(isend_pool *isp)
{
	request_node *it = isp->reqs;
	while(it) {
		MPI_Wait(&it->req, MPI_STATUS_IGNORE);
		request_node *tmp = it;
		it = it->next;
		free(tmp);
	}
	isp->reqs = NULL;
}

void isend_pool_destroy(isend_pool *isp)
{
	isend_pool_wait_all(isp);
	free(isp);
}
