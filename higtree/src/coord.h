#ifndef COORD_H
#define COORD_H

#include "dim.h"

#include <stdio.h>
#include <assert.h>
#include <stdbool.h>
#include "types.h"
#include "utils.h"

//! Defines the type of a coordinate
typedef real coordtype;

//! Defines a point with DIM coordinates
typedef coordtype Point[DIM];

//! Defines a pointer to a point
typedef coordtype *PPoint;

//! Defines a pointer to const point
typedef const coordtype *CPPoint;

//! Assigns the contents of a point to another.
//! It can be used to assign any vector with DIM elements.
// x[0] = y[0]; x[1] = y[1]; ...
#define POINT_ASSIGN(x, y) \
{\
	int POINT__i;\
	for(POINT__i = 0; POINT__i < DIM; POINT__i++) {\
		x[POINT__i] = y[POINT__i];\
	}\
}

//! Assigns the scalar v to each coordinate of the point x.
//! It can be used to assign v to each element of a vector with DIM elments.
// x[0] = v; x[1] = v; ...
#define POINT_ASSIGN_SCALAR(x, v) \
{\
	int POINT__i; \
	for(POINT__i = 0; POINT__i < DIM; POINT__i++) {\
		x[POINT__i] = v;\
	}\
}

//! Assigns the scalars v_0, v_1... of type to each coordinate of the point x.
// x[0] = v_0; x[1] = v_1; ...
#define POINT_ASSIGN_SCALARS(x, type, vs...) \
{\
	type POINT__vals[] = {vs}; \
	POINT_ASSIGN(x, POINT__vals); \
}

//! Assigns the integer values v_0, v_1... to each coordinate of the point x.
#define POINT_ASSIGN_INTS(x, vs...) POINT_ASSIGN_SCALARS(x, int, vs)

//! Assigns the real values v_0, v_1... to each coordinate of the point x.
#define POINT_ASSIGN_REALS(x, vs...) POINT_ASSIGN_SCALARS(x, real, vs)

//! Assigns to each coordinate of the point x the result of operation op applied to each element of the points y and z.
// x = y op z;
#define POINT_OPERATOR(x, y, op, z) \
{\
	int POINT__i; \
	for(POINT__i = 0; POINT__i < DIM; POINT__i++) {\
		x[POINT__i] = (y[POINT__i]) op (z[POINT__i]);\
	}\
}

//! Assigns to each coordinate of the point x the result of operation op applied to each element of the point y and the scalar v.
#define POINT_OPERATOR_SCALAR(x, y, op, v) \
{\
	int POINT__i; \
	for(POINT__i = 0; POINT__i < DIM; POINT__i++) {\
		x[POINT__i] = y[POINT__i] op (v);\
	}\
}

//! Assigns to each coordinate of the point x the result of operation op applied to each element of the point y and the scalars v_0, v_1 ...
#define POINT_OPERATOR_SCALARS(x, y, op, type, vs...) \
{\
	type POINT__vals[] = {vs}; \
	POINT_OPERATOR(x, y, op, POINT__vals);\
}

//! Assigns to each coordinate of the point x the result of operation op applied to each element of the point y and the integer values v_0, v_1 ...
#define POINT_OPERATOR_INTS(x, y, op, vs...) POINT_OPERATOR_SCALARS(x, y, op, int, vs)
//! Assigns to each coordinate of the point x the result of operation op applied to each element of the point y and the real values v_0, v_1 ...
#define POINT_OPERATOR_REALS(x, y, op, vs...) POINT_OPERATOR_SCALARS(x, y, op, real, vs)

//! Assigns to each coordinate of the point x the sum of the respective coordinate of points y and z.
#define POINT_ADD(x, y, z)  POINT_OPERATOR(x, y, +, z)

//! Assigns to each coordinate of the point x the difference of the respective coordinate of points y and z.
#define POINT_SUB(x, y, z)  POINT_OPERATOR(x, y, -, z)

//! Assigns to each coordinate of the point x the ratio of the respective coordinate of points y and z.
#define POINT_DIV(x, y, z)  POINT_OPERATOR(x, y, /, z)

//! Assigns to each coordinate of the point x the product of the respective coordinate of points y and z.
#define POINT_MULT(x, y, z) POINT_OPERATOR(x, y, *, z)

//! Assigns to each coordinate of the point x the sum of the coordinate of point y and scalar v.
#define POINT_ADD_SCALAR(x, y, v)  POINT_OPERATOR_SCALAR(x, y, +, v)

//! Assigns to each coordinate of the point x the difference of the coordinate of point y and scalar v.
#define POINT_SUB_SCALAR(x, y, v)  POINT_OPERATOR_SCALAR(x, y, -, v)

//! Assigns to each coordinate of the point x the ratio of the coordinate of point y and scalar v.
#define POINT_DIV_SCALAR(x, y, v)  POINT_OPERATOR_SCALAR(x, y, /, v)

//! Assigns to each coordinate of the point x the product of the coordinate of point y and scalar v.
#define POINT_MULT_SCALAR(x, y, v) POINT_OPERATOR_SCALAR(x, y, *, v)

//! Assigns to each coordinate of the point x the sum of the coordinates of point y and integer values v_0, v_1....
#define POINT_ADD_INTS(x, y, vs...)  POINT_OPERATOR_INTS(x, y, +, vs)

//! Assigns to each coordinate of the point x the difference of the coordinates of point y and integer values v_0, v_1....
#define POINT_SUB_INTS(x, y, vs...)  POINT_OPERATOR_INTS(x, y, -, vs)

//! Assigns to each coordinate of the point x the ratio of the coordinates of point y and integer values v_0, v_1....
#define POINT_DIV_INTS(x, y, vs...)  POINT_OPERATOR_INTS(x, y, /, vs)

//! Assigns to each coordinate of the point x the ratio of the coordinates of point y and integer values v_0, v_1....
#define POINT_MULT_INTS(x, y, vs...)  POINT_OPERATOR_INTS(x, y, *, vs)

//! Assigns to each coordinate of the point x the sum of the coordinates of point y and real values v_0, v_1....
#define POINT_ADD_REALS(x, y, vs...)  POINT_OPERATOR_REALS(x, y, +, vs)

//! Assigns to each coordinate of the point x the difference of the coordinates of point y and real values v_0, v_1....
#define POINT_SUB_REALS(x, y, vs...)  POINT_OPERATOR_REALS(x, y, -, vs)

//! Assigns to each coordinate of the point x the ratio of the coordinates of point y and real values v_0, v_1....
#define POINT_DIV_REALS(x, y, vs...)  POINT_OPERATOR_REALS(x, y, /, vs)

//! Assigns to each coordinate of the point x the product of the coordinates of point y and real values v_0, v_1....
#define POINT_MULT_REALS(x, y, vs...)  POINT_OPERATOR_REALS(x, y, *, vs)

//! Computes in POINT__res the iterated calculation of ((init op x[0]) op x[1]) op ...)
#define POINT_FOLD(POINT__res, x, init, op) \
do {\
	POINT__res = (init); \
	int POINT__i; \
	for(POINT__i = 0; POINT__i < DIM; POINT__i++) {\
		POINT__res = ((POINT__res) op (x[POINT__i]));\
	}\
} while(0)

//! Determines whether x[0] comp && x[1] comp && ... is true. comp is of the form "rel-op value"
#define POINT_ASSERT_ALL_TRUE(x, comp) \
do {\
	int POINT__i; \
	for(POINT__i = 0; POINT__i < DIM; POINT__i++) {\
		assert(x[POINT__i] comp);\
	}\
} while(0)

//! Prints the coordinates of point coord, separeated by sep
void co_print_coord(FILE *f, Point coord, const char sep[]);

//! Computes the euclidean distance between points l and h
real co_distance(const Point l, const Point h);

//! Computes the scalar product (dot product) between a and b
real co_dot(CPPoint a, CPPoint b);

//! Computes a point with random coordinates within the bounding box determined by points p1 and p2
void randpoint(Point p1, Point p2, Point point);

//! Determines whether point x is within the bounding box determined by points l and h
int co_is_within_bounding_box(Point l, Point h, Point x);

//! Determines whether two points are equal, i.e., respective coordinates are
// the same within tolerance of EPSDELTA.
bool co_equal(const Point p1, const Point p2);

//! Determines whether two points are equal, i.e., respective coordinates are
// the same within tolerance of EPSMACH.
bool Pointequal(const Point p1, const Point p2);

//! Determines whether two points are equal within some specified tolerance.
bool point_equal_tol(const Point p1, const Point p2, real tol);

#endif
