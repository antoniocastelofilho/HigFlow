#include <string.h>
#include <atf-c.h>
#include "coord.h"
#include "wls.h"
#include "utils.h"
#define DEBUG
#include "Debug-c.h"


ATF_TC(wls_1_degree);
ATF_TC_HEAD(wls_1_degree, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(wls_1_degree, tc)
{
	Point pts[6];
	Point x;
	real w[6];
	POINT_ASSIGN_REALS(pts[0], -1.0, -1.0);
	POINT_ASSIGN_REALS(pts[1], -1.0, 0.0);
	POINT_ASSIGN_REALS(pts[2], -1.0, 1.0);
	POINT_ASSIGN_REALS(pts[3], 1.0, -1.0);
	POINT_ASSIGN_REALS(pts[4], 1.0, 0.0);
	POINT_ASSIGN_REALS(pts[5], 1.0, 1.0);

	wls_interpolator *wls = wls_create(1, 10);
	POINT_ASSIGN_REALS(x, 0.0, 0.0);
	wls_set_points(wls, 6, pts, x);
	real v[6] = {
		-5,
		-1,
		3,
		1,
		5,
		9,
	};

	wls_calc(wls, x, 6, w);
	real val = 0.0;
	for(int i = 0; i < 6; i++) {
		val += w[i] * v[i];
	}
	DEBUG_INSPECT(val, %g);
	ATF_CHECK(FLT_EQ(val, 2.0));

	POINT_ASSIGN_REALS(x, 0.5, 0.5);
	wls_set_points(wls, 6, pts, x);
	wls_calc(wls, x, 6, w);
	val = 0.0;
	for(int i = 0; i < 6; i++) {
		val += w[i] * v[i];
	}
	DEBUG_INSPECT(val, %g);
	ATF_CHECK(FLT_EQ(val, 5.5));

	POINT_ASSIGN_REALS(x, 1.0, 1.0);
	wls_set_points(wls, 6, pts, x);
	wls_calc(wls, x, 6, w);
	val = 0.0;
	for(int i = 0; i < 6; i++) {
		val += w[i] * v[i];
	}
	DEBUG_INSPECT(val, %g);
	ATF_CHECK(FLT_EQ(val, 9.0));

	Point pts2[6];
	Point x2;
	real w2[6];
	POINT_ASSIGN_REALS(pts2[0], -1.0, -1.0);
	POINT_ASSIGN_REALS(pts2[1], 0.0, -1.0);
	POINT_ASSIGN_REALS(pts2[2], 0.0, 1.0);
	POINT_ASSIGN_REALS(pts2[3], 1.0, -1.0);
	POINT_ASSIGN_REALS(pts2[4], 1.0, 0.0);
	POINT_ASSIGN_REALS(pts2[5], 1.0, 1.0);

	real v2[6] = {
		-5,
		-2,
		6,
		1,
		5,
		9,
	};

	POINT_ASSIGN_REALS(x2, 1.0, 1.0);
	wls_set_points(wls, 6, pts2, x2);
	wls_calc(wls, x2, 6, w2);
	val = 0.0;
	for(int i = 0; i < 6; i++) {
		val += w2[i] * v2[i];
	}
	DEBUG_INSPECT(val, %g);
	ATF_CHECK(FLT_EQ(val, 9.0));

	wls_destroy(wls);
}


ATF_TC(wls_2_degree);
ATF_TC_HEAD(wls_2_degree, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(wls_2_degree, tc)
{
	Point pts[6];
	Point x;
	real w[6];
	POINT_ASSIGN_REALS(pts[0], -1.0, -1.0);
	POINT_ASSIGN_REALS(pts[1], -1.0, 0.0);
	POINT_ASSIGN_REALS(pts[2], -1.0, 1.0);
	POINT_ASSIGN_REALS(pts[3], 1.0, -1.0);
	POINT_ASSIGN_REALS(pts[4], 1.0, 0.0);
	POINT_ASSIGN_REALS(pts[5], 1.0, 1.0);
	POINT_ASSIGN_REALS(pts[6], 0.0, 0.0);

	POINT_ASSIGN_REALS(x, 0.5, 0.5);

	wls_interpolator *wls = wls_create(2, 10);
	wls_set_points_and_calc(wls, 7, pts, x, w);

	real v[7] = {
		-1,
		1,
		5,
		5,
		1,
		-1,
		0,
	};
	real val = 0.0;
	for(int i = 0; i < 6; i++) {
		val += w[i] * v[i];
	}
	DEBUG_INSPECT(val, %g);
	ATF_CHECK(FLT_EQ(val, -0.25));
	wls_destroy(wls);
}

ATF_TC(wls_3_and_4_degree);
ATF_TC_HEAD(wls_3_and_4_degree, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(wls_3_and_4_degree, tc)
{
	Point pts[15];
	Point x;
	real w[15];
	POINT_ASSIGN_REALS(pts[0], -1.0, -1.0);
	POINT_ASSIGN_REALS(pts[1], -1.0, 0.0);
	POINT_ASSIGN_REALS(pts[2], -1.0, 1.0);
	POINT_ASSIGN_REALS(pts[3], 0.0, -1.0);
	POINT_ASSIGN_REALS(pts[4], 0.0, 0.0);
	POINT_ASSIGN_REALS(pts[5], 0.0, 1.0);
	POINT_ASSIGN_REALS(pts[6], 1.0, -1.0);
	POINT_ASSIGN_REALS(pts[7], 1.0, 0.0);
	POINT_ASSIGN_REALS(pts[8], 1.0, 1.0);
	POINT_ASSIGN_REALS(pts[9], -2.0, 0.0);
	POINT_ASSIGN_REALS(pts[10], 2.0, 0.0);
	POINT_ASSIGN_REALS(pts[11], 0.0, -2.0);
	POINT_ASSIGN_REALS(pts[12], 0.0, 2.0);
	POINT_ASSIGN_REALS(pts[13], -2.0, -2.0);
	POINT_ASSIGN_REALS(pts[14], 2.0, 2.0);

	wls_interpolator *wls = wls_create(3, 15);
	real v[15] = {
		-1 -1 -1 -1,
		-1 -1 -1,
		-1 -1 -1 -1,
		-1 -1,
		-1,
		-1 -1,
		1 -1 +1 -1,
		1 +1 -1,
		1 -1 +1 -1,
		-8 -2 -1,
		8 +2 -1,
		-4 -1,
		-4 -1,
		-8 -4 -2 -1,
		8 -4 +2 -1,
	};

	POINT_ASSIGN_REALS(x, 0.5, 0.5);
	wls_set_points(wls, 15, pts, x);
	wls_calc(wls, x, 15, w);
	real val = 0.0;
	for(int i = 0; i < 15; i++) {
		val += w[i] * v[i];
	}
	ATF_CHECK(FLT_EQ(val, -0.625));
	wls_destroy(wls);

	wls = wls_create(4, 15);
	POINT_ASSIGN_REALS(x, 0.5, 0.5);
	wls_set_points(wls, 15, pts, x);
	wls_calc(wls, x, 15, w);
	val = 0.0;
	for(int i = 0; i < 15; i++) {
		val += w[i] * v[i];
	}
	ATF_CHECK(FLT_EQ(val, -0.625));
	wls_destroy(wls);
}

ATF_TC(wls_colinears);
ATF_TC_HEAD(wls_colinears, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(wls_colinears, tc)
{
	Point pts[5], ptsadd[1];
	Point x;
	real w[5];
	POINT_ASSIGN_REALS(pts[0], -2.0, 0.0);
	POINT_ASSIGN_REALS(pts[1], -1.0, 0.0);
	POINT_ASSIGN_REALS(pts[2], 0.0, 0.0);
	POINT_ASSIGN_REALS(pts[3], 1.0, 0.0);
	POINT_ASSIGN_REALS(pts[4], 2.0, 0.0);
	POINT_ASSIGN_REALS(ptsadd[0], 2.0, 0.0);

	wls_interpolator *wls = wls_create(1, 5);
	real v[5] = {
		-1,
		0,
		1,
		2,
		3,
	};

	POINT_ASSIGN_REALS(x, 0.5, 0.5);
	wls_set_points(wls, 4, pts, x);
	wls_add_points(wls, 4, 1, ptsadd, x);
	wls_calc(wls, x, 5, w);
	real val = 0.0;
	for(int i = 0; i < 5; i++) {
		val += w[i] * v[i];
	}
	ATF_CHECK(FLT_EQ(val, 1.5));
	wls_destroy(wls);

	wls = wls_create(1, 5);
	wls_set_type(wls, WLSSTATIC);
	wls_set_points(wls, 4, pts, x);
	wls_add_points(wls, 4, 1, ptsadd, x);
	wls_calc(wls, x, 5, w);
	val = 0.0;
	for(int i = 0; i < 5; i++) {
		val += w[i] * v[i];
	}
	ATF_CHECK(FLT_EQ(val, 1.5));
	wls_destroy(wls);

	wls = wls_create(1, 5);
	wls_set_type(wls, WLSSTATIC);
        wls_set_points_and_calc(wls,5,pts,x,w);
	val = 0.0;
	for(int i = 0; i < 5; i++) {
		val += w[i] * v[i];
	}
	ATF_CHECK(FLT_EQ(val, 1.5));
	wls_destroy(wls);

	wls = wls_create(1, 5);
        wls_set_points_and_calc(wls,5,pts,x,w);
	val = 0.0;
	for(int i = 0; i < 5; i++) {
		val += w[i] * v[i];
	}
	ATF_CHECK(FLT_EQ(val, 1.5));
	wls_destroy(wls);
}

ATF_TC(wls_1_degree_add);
ATF_TC_HEAD(wls_1_degree_add, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(wls_1_degree_add, tc)
{
	Point pts[4], ptsadd[2];
	Point x;
	real w[6];
	POINT_ASSIGN_REALS(pts[0], -1.0, -1.0);
	POINT_ASSIGN_REALS(pts[1], -1.0, 0.0);
	POINT_ASSIGN_REALS(pts[2], -1.0, 1.0);
	POINT_ASSIGN_REALS(pts[3], 1.0, -1.0);
	POINT_ASSIGN_REALS(ptsadd[0], 1.0, 0.0);
	POINT_ASSIGN_REALS(ptsadd[1], 1.0, 1.0);

	wls_interpolator *wls = wls_create(1, 10);
	real v[6] = {
		-5,
		-1,
		3,
		1,
		5,
		9,
	};

	POINT_ASSIGN_REALS(x, 0.0, 0.0);
        wls_set_points_and_calc(wls,4,pts,x,w);
        wls_add_points_and_calc(wls,4,2,ptsadd,x,w);
	real val = 0.0;
	for(int i = 0; i < 6; i++) {
		val += w[i] * v[i];
	}
	DEBUG_INSPECT(val, %g);
	ATF_CHECK(FLT_EQ(val, 2.0));
	wls_destroy(wls);

	wls = wls_create(1, 10);
	wls_set_type(wls, WLSSTATIC);
        wls_set_points_and_calc(wls,4,pts,x,w);
        wls_add_points_and_calc(wls,4,2,ptsadd,x,w);
	val = 0.0;
	for(int i = 0; i < 6; i++) {
		val += w[i] * v[i];
	}
	DEBUG_INSPECT(val, %g);
	ATF_CHECK(FLT_EQ(val, 2.0));
	wls_destroy(wls);
}

ATF_TC(wls_3_and_4_degree_add);
ATF_TC_HEAD(wls_3_and_4_degree_add, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(wls_3_and_4_degree_add, tc)
{
	Point pts[10], ptsadd[5];
	Point x;
	real w[15];
	POINT_ASSIGN_REALS(pts[0], -1.0, -1.0);
	POINT_ASSIGN_REALS(pts[1], -1.0, 0.0);
	POINT_ASSIGN_REALS(pts[2], -1.0, 1.0);
	POINT_ASSIGN_REALS(pts[3], 0.0, -1.0);
	POINT_ASSIGN_REALS(pts[4], 0.0, 0.0);
	POINT_ASSIGN_REALS(pts[5], 0.0, 1.0);
	POINT_ASSIGN_REALS(pts[6], 1.0, -1.0);
	POINT_ASSIGN_REALS(pts[7], 1.0, 0.0);
	POINT_ASSIGN_REALS(pts[8], 1.0, 1.0);
	POINT_ASSIGN_REALS(pts[9], -2.0, 0.0);
	POINT_ASSIGN_REALS(ptsadd[0], 2.0, 0.0);
	POINT_ASSIGN_REALS(ptsadd[1], 0.0, -2.0);
	POINT_ASSIGN_REALS(ptsadd[2], 0.0, 2.0);
	POINT_ASSIGN_REALS(ptsadd[3], -2.0, -2.0);
	POINT_ASSIGN_REALS(ptsadd[4], 2.0, 2.0);

	wls_interpolator *wls = wls_create(3, 15);
	real v[15] = {
		-1 -1 -1 -1,
		-1 -1 -1,
		-1 -1 -1 -1,
		-1 -1,
		-1,
		-1 -1,
		1 -1 +1 -1,
		1 +1 -1,
		1 -1 +1 -1,
		-8 -2 -1,
		8 +2 -1,
		-4 -1,
		-4 -1,
		-8 -4 -2 -1,
		8 -4 +2 -1,
	};

	POINT_ASSIGN_REALS(x, 0.5, 0.5);
        wls_set_points_and_calc(wls,10,pts,x,w);
        wls_add_points_and_calc(wls,10,5,ptsadd,x,w);
	real val = 0.0;
	for(int i = 0; i < 15; i++) {
		val += w[i] * v[i];
	}
	ATF_CHECK(FLT_EQ(val, -0.625));
	wls_destroy(wls);

	wls = wls_create(4, 15);
	POINT_ASSIGN_REALS(x, 0.5, 0.5);
        wls_set_points_and_calc(wls,10,pts,x,w);
        wls_add_points_and_calc(wls,10,5,ptsadd,x,w);
	val = 0.0;
	for(int i = 0; i < 15; i++) {
		val += w[i] * v[i];
	}
	ATF_CHECK(FLT_EQ(val, -0.625));
	wls_destroy(wls);
}

ATF_TC(wls_3_and_4_degree_static_add);
ATF_TC_HEAD(wls_3_and_4_degree_static_add, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(wls_3_and_4_degree_static_add, tc)
{
	Point pts[10], ptsadd[5];
	Point x;
	real w[15];
	POINT_ASSIGN_REALS(pts[0], -1.0, -1.0);
	POINT_ASSIGN_REALS(pts[1], -1.0, 0.0);
	POINT_ASSIGN_REALS(pts[2], -1.0, 1.0);
	POINT_ASSIGN_REALS(pts[3], 0.0, -1.0);
	POINT_ASSIGN_REALS(pts[4], 0.0, 0.0);
	POINT_ASSIGN_REALS(pts[5], 0.0, 1.0);
	POINT_ASSIGN_REALS(pts[6], 1.0, -1.0);
	POINT_ASSIGN_REALS(pts[7], 1.0, 0.0);
	POINT_ASSIGN_REALS(pts[8], 1.0, 1.0);
	POINT_ASSIGN_REALS(pts[9], -2.0, 0.0);
	POINT_ASSIGN_REALS(ptsadd[0], 2.0, 0.0);
	POINT_ASSIGN_REALS(ptsadd[1], 0.0, -2.0);
	POINT_ASSIGN_REALS(ptsadd[2], 0.0, 2.0);
	POINT_ASSIGN_REALS(ptsadd[3], -2.0, -2.0);
	POINT_ASSIGN_REALS(ptsadd[4], 2.0, 2.0);

	wls_interpolator *wls = wls_create(3, 15);
	wls_set_type(wls, WLSSTATIC);
	real v[15] = {
		-1 -1 -1 -1,
		-1 -1 -1,
		-1 -1 -1 -1,
		-1 -1,
		-1,
		-1 -1,
		1 -1 +1 -1,
		1 +1 -1,
		1 -1 +1 -1,
		-8 -2 -1,
		8 +2 -1,
		-4 -1,
		-4 -1,
		-8 -4 -2 -1,
		8 -4 +2 -1,
	};

	POINT_ASSIGN_REALS(x, 0.5, 0.5);
        wls_set_points(wls,10,pts,x);
        wls_add_points(wls,10,5,ptsadd,x);
	wls_calc(wls, x, 15, w);
	real val = 0.0;
	for(int i = 0; i < 15; i++) {
		val += w[i] * v[i];
	}
	ATF_CHECK(FLT_EQ(val, -0.625));
	wls_destroy(wls);





	wls = wls_create(4, 15);
	POINT_ASSIGN_REALS(x, 0.5, 0.5);
        wls_set_points_and_calc(wls,10,pts,x,w);
        wls_add_points_and_calc(wls,10,5,ptsadd,x,w);
	val = 0.0;
	for(int i = 0; i < 15; i++) {
		val += w[i] * v[i];
	}
	ATF_CHECK(FLT_EQ(val, -0.625));
	wls_destroy(wls);
}

ATF_TP_ADD_TCS(tp)
{
    ATF_TP_ADD_TC(tp, wls_1_degree);
//    ATF_TP_ADD_TC(tp, wls_2_degree); // Fix it!
    ATF_TP_ADD_TC(tp, wls_3_and_4_degree);
    ATF_TP_ADD_TC(tp, wls_colinears);
    ATF_TP_ADD_TC(tp, wls_1_degree_add);
    ATF_TP_ADD_TC(tp, wls_3_and_4_degree_add);
    ATF_TP_ADD_TC(tp, wls_3_and_4_degree_static_add);
}
