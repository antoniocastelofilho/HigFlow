#ifndef RECT_H
#define RECT_H

#include <stdbool.h>
#include "coord.h"

typedef struct Rect {
	Point lo;
	Point hi;
} Rect;

/** Calculates the intersection of rects a and b, returning true if exists.
 *
 * If there is an intersection (even if degenerated) between a and b, sets
 * it to the third argument, and returns true. Returns false, otherwise.
 * Parameter @p intersection may be NULL.
 */
bool rect_intersect(const Rect *a, const Rect *b, Rect *intersection);

/** Sets the first parameter as the axis-aligned bounding box of both rects. */
void rect_expand(Rect *into, const Rect *other);

/** Checks if later rect is fully contained in first rect. */
bool rect_contains(const Rect *container, const Rect *contained);

/** Checks point x is contained in container. */
bool rect_contains_point(const Rect *container, CPPoint x);

/** Calculates the distance of a rect to a point. */
real rect_distance_to_point(const Rect *r, CPPoint x);

/** Check if two degenerate rects are parallel. */
bool rect_are_parallel(const Rect *a, const Rect *b);

/** Check if two rects have the same dimensions. */
bool rect_has_same_size(const Rect *a, const Rect *b);

/** Count the number of degenete dimensions in a rect. */
unsigned rect_degenerate_count(const Rect *r);

/** Translate rect by displacement. */
void rect_translate(const Point displacement, const Rect *from, Rect *to);

#endif
