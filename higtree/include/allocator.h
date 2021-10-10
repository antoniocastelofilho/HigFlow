#ifndef ALLOCATOR_H
#define ALLOCATOR_H

#include <stdlib.h>

// This is a simple allocator module that allows bulk deallocation of all
// allocated buffers.

struct allocator;
typedef struct allocator allocator;

allocator *allocator_init();

/*! Allocates a buffer managed by this allocator.
 *
 * Do not write outside the boundary of allocated region!
 * If you do, code will break spectacularly.
 */
void *allocator_alloc(allocator *al, size_t bytes_size);

/*! Deallocates a pointer returned by allocator_alloc().
 *
 * If ptr was not returned by a call of allocator_alloc() on @p al
 * you will probably (hopefully?) have a segmentation fault.
 */
void allocator_dealloc(allocator *al, void *ptr);

void allocator_destroy(allocator *al);

#endif
