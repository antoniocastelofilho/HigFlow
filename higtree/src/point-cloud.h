#ifndef __POINT_CLOUD_H
#define __POINT_CLOUD_H

#include "coord.h"

//! Defines a set of points.
typedef struct point_cloud {
	int numpts;
	int maxpts;
	Point *pts;
} point_cloud;

//! Creates a point cloud.
point_cloud *pc_create() ;

//! Resets a point cloud, i.e., removes all points.
void pc_reset(point_cloud *pc) ;

//! Adds a point to the cloud.
void pc_add_point(point_cloud *pc, Point pt) ;

//! Adds numpts points to the cloud.
void pc_add_points(point_cloud *pc, int numpts, Point *pts) ;

//! Gets the number of points in the cloud.
int pc_get_numpts(point_cloud *pc) ;

//! Copies the i-th point
void pc_copy_point(point_cloud *pc, int i, Point pt) ;

//! Gets the point to the i-th point.
PPoint pc_get_point(point_cloud *pc, int i) ;

//! Adds the point x to all points in pcin.
void pc_calc_diff(point_cloud *pcin, Point x, point_cloud *pcout);

//! Computes a hash code for the point cloud.
int pc_hash_key(point_cloud *pc) ;

//! Checks whether two point clouds are equal.
int pc_equal(point_cloud *pc1, point_cloud *pc2) ;

//! Destroys the point cloud.
void pc_destroy(point_cloud *pc) ;


#endif
