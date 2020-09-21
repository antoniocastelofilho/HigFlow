#ifndef __WLS_H
#define __WLS_H

/* ****************************************************** */
/*                  incluindo bibliotecas                 */
/* ****************************************************** */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "coord.h"
#include "types.h"


/* ****************************************************** */
/*   declarando todas as rotinas em wls.c                 */
/* ****************************************************** */

//! \brief. Determines the minimum number of points required to have a
//! good approximation of a given dimension and order.
int wls_num_min_points (unsigned dim, int ord);

//! Type of least square method.
typedef enum {WLSMOVING, WLSSTATIC} wls_type;

//! Defines a single sample used for interpolation.
typedef struct wls_item {
	Point x;
	real dist;
	struct sim_boundary *bcpoint;
	int id;
	unsigned char type;
} wls_item;

//! Defines an interpolator.
typedef struct wls_interpolator {
        mreal  A, Q;
        vreal  b;
        vreal  D;
	vreal  w;
	int ord;
	wls_type type;
	int maxnumpts;
	int numpoly;
	unsigned dim;

	vreal poly_coefs;
} wls_interpolator;

//! Creates an interpolator, given the order of interpolation and the maximum number of points which the interpolator may be used.
wls_interpolator * wls_create(unsigned dim, int ord, int maxnumpts);

//! Set the interpolator type. The default is Moving.
void wls_set_type(wls_interpolator *wls, wls_type t);

//! Sets the which points to use in the interpolation and the point where to interpolate.
real wls_set_points(wls_interpolator *wls, int numpts, Point pts[], Point x);

//! Adds points to be used in the interpolation.
real wls_add_points(wls_interpolator *wls, int numpts, int addnumpts, Point pts[], Point x);

//! Computes the interpolation. Computes the w weights of numpts points.
void wls_calc(wls_interpolator *wls, Point x, int numpts, real w[]);

//! Sets the which points to use in the interpolation and the point where to interpolate. Then, computes the interpolation.
real wls_set_points_and_calc(wls_interpolator *wls, int numpts, Point pts[], Point x, real w[]);

//! Verifies if order was achieved.
bool wls_is_good_enough(wls_interpolator *wls, int numusedpts, wls_item items[], const Point x, real w[]);

//! Sets the which samples to use in the interpolation and the point where to interpolate. Then, computes the interpolation.
//
// This function is similar to wls_set_points_and_calc(), but takes the points from an array
// of struct wls_item
real wls_set_samples_and_calc(wls_interpolator *wls, int numpts, wls_item items[], const Point x, real w[]);

//! Adds points to be used in the interpolation. Then, computes the interpolation.
real wls_add_points_and_calc(wls_interpolator *wls, int numpts, int addnumpts, Point pts[], Point x, real w[]);

//! Destroys the interpolator.
void wls_destroy(wls_interpolator *wls);



#endif
