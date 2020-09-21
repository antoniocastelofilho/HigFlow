
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "coord.h"
#include "utils.h"

void co_print_coord(FILE *f, Point coord, const char sep[]) {
	fprintf(f, "%2.15f", coord[0]);
	for(int i = 1; i < DIM; i++) {
		fprintf(f, "%s%2.15f", sep, coord[i]);
	}
}

real co_distance(const Point l, const Point h) {
	real sum = 0;
	for(unsigned dim = 0; dim < DIM; ++dim) {
		const real delta = l[dim] - h[dim];
		sum += delta * delta;
	}
	return sqrt(sum);
}

real co_dot(CPPoint a, CPPoint b)
{
	real sum = 0;
	for(unsigned dim = 0; dim < DIM; ++dim) {
		sum += a[dim] * b[dim];
	}
	return sum;
}

bool point_equal_tol(const Point p1, const Point p2, real tol)
{
	for (int i = 0; i < DIM; i++) {
		if (fabs(p1[i] - p2[i]) >= tol) {
			return false;
		}
	}
	return true;
}

bool co_equal(const Point p1, const Point p2)
{
	return point_equal_tol(p1, p2, EPSDELTA);
}

bool Pointequal(const Point p1, const Point p2)
{
	return point_equal_tol(p1, p2, EPSMACH);
}

void randpoint(Point p1, Point p2, Point point)
{
	Point r, delta;
	POINT_ASSIGN_SCALAR(r, ((real)rand()/(real)RAND_MAX));
	POINT_SUB(delta, p2, p1);
	POINT_MULT(delta, delta, r);
	POINT_ADD(point, p1, delta);
}

int co_is_within_bounding_box(Point l, Point h, Point x) {
	for(int dim = 0; dim < DIM; dim++) {
		if (FLT_LT(x[dim], l[dim]) || FLT_GT(x[dim], h[dim])) {
			return 0;
		}
	}
	return 1;
}
