#include <stddef.h>
#include "utils.h"

#include "allocator.h"

// Implemented as circular double linked list.

struct allocator
{
	allocator *next;
	allocator *prev;
	uint8_t buf[];
};

allocator *allocator_init()
{
	DECL_AND_ALLOC(allocator, ret, 1);
	ret->next = ret;
	ret->prev = ret;

	return ret;
}

void *allocator_alloc(allocator *al, size_t bytes_size)
{
	DECL_AND_ALLOC_FLEX(allocator, ret, buf, bytes_size);

	ret->next = al->next;
	al->next = ret;

	ret->prev = al;
	ret->next->prev = ret;

	return ret->buf;
}

void allocator_dealloc(allocator *al, void *ptr)
{
	// We are trusting the ptr was allocated with allocator_alloc() on al.
	// al doesn't really serves anything, but I was afraid of making this
	// function taking a single void* parameter and being confused with
	// allocator_destroy().
	allocator *node = (allocator *)((char*)ptr - offsetof(allocator, buf));
	assert(node != al);

	// Detach itself from the linked list:
	node->prev->next = node->next;
	node->next->prev = node->prev;

	// Deallocate:
	free(node);
}

void allocator_destroy(allocator *al)
{
	// Break the cycle
	al->prev->next = NULL;

	// Forward iterate while freeing
	do {
		allocator *tmp = al;
		al = al->next;
		free(tmp);
	} while(al);
}
