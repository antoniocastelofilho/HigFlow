// A bump-style pool whose point is BULK RELEASE: allocator_destroy() frees
// everything at once, so callers with many short-lived buffers of the same lifetime
// need not track them individually.
//
// The two warnings in the function comments are literal, not decorative -- this
// allocator keeps its bookkeeping next to the returned blocks, so a write past the
// end of one allocation corrupts the pool rather than the neighbouring data, and the
// failure surfaces far from its cause.

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
