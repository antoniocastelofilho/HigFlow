#include "utils.h"
#include "rect.h"

bool
rect_intersect(const Rect *a, const Rect *b, Rect *intersection)
{
	Rect isect;

	for(unsigned dim = 0; dim < DIM; ++dim) {
		isect.hi[dim] = min(a->hi[dim], b->hi[dim]);
		isect.lo[dim] = max(a->lo[dim], b->lo[dim]);
		if(!POS_LE(isect.lo[dim], isect.hi[dim]))
			return false;
	}

	if(intersection) {
		*intersection = isect;
	}

	return true;
}


void
rect_expand(Rect *into, const Rect *other)
{
	for(unsigned dim = 0; dim < DIM; ++dim) {
		into->hi[dim] = max(into->hi[dim], other->hi[dim]);
		into->lo[dim] = min(into->lo[dim], other->lo[dim]);
	}
}

bool rect_contains(const Rect *container, const Rect *contained)
{
	for(unsigned dim = 0; dim < DIM; ++dim) {
		if(!POS_LE(contained->hi[dim], container->hi[dim]))
			return false;
		if(!POS_GE(contained->lo[dim], container->lo[dim]))
			return false;
	}
	return true;
}

bool rect_contains_point(const Rect *container, CPPoint x)
{
	for(unsigned dim = 0; dim < DIM; ++dim) {
		if(!POS_LE(x[dim], container->hi[dim]))
			return false;
		if(!POS_GE(x[dim], container->lo[dim]))
			return false;
	}
	return true;
}

real rect_distance_to_point(const Rect *r, CPPoint x)
{
	real accum = 0;
	for(unsigned dim = 0; dim < DIM; ++dim) {
		const real delta = max(max(r->lo[dim] - x[dim], 0.0), x[dim] - r->hi[dim]);
		accum += delta * delta;
	}
	return sqrt(accum);
}

bool rect_are_parallel(const Rect *a, const Rect *b)
{
	for(unsigned dim = 0; dim < DIM; ++dim) {
		if(
			POS_EQ(a->lo[dim], a->hi[dim]) &&
			POS_EQ(b->lo[dim], b->hi[dim])
		) {
			return true;
		}
	}
	return false;
}

bool rect_has_same_size(const Rect *a, const Rect *b)
{
	for(unsigned i = 0; i < DIM; ++i) {
		if(!POS_EQ(a->hi[i] - a->lo[i], b->hi[i] - b->lo[i])) {
			return false;
		}
	}
	return true;
}

unsigned rect_degenerate_count(const Rect *r)
{
	unsigned count = 0;
	for(unsigned i = 0; i < DIM; ++i) {
		count += !!POS_EQ(r->hi[i], r->lo[i]);
	}
	return count;
}

void rect_translate(const Point displacement, const Rect *from, Rect *to)
{
	POINT_ADD(to->lo, from->lo, displacement);
	POINT_ADD(to->hi, from->hi, displacement);
}
