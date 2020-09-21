#include <string.h>
#include <atf-c.h>
#include "coord.h"
#include "point-cloud.h"
#include "utils.h"
#define DEBUG
#include "Debug-c.h"

ATF_TC(point_cloud_basic);
ATF_TC_HEAD(point_cloud_basic, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(point_cloud_basic, tc)
{
	point_cloud *pc = pc_create();

	for(int i = 1; i <= 15; i++) {
		Point p;
		POINT_ASSIGN_REALS(p, i, -i);
		pc_add_point(pc, p);
	}

	ATF_CHECK(pc_get_numpts(pc) == 15);

	for(int i = 0; i < 15; i++) {
		PPoint p;
		Point x;
		p = pc_get_point(pc, i);
		pc_copy_point(pc, i, x);
		ATF_CHECK(FLT_EQ(p[0], (i + 1)));
		ATF_CHECK(FLT_EQ(p[1], -(i + 1)));
		ATF_CHECK(FLT_EQ(x[0], (i + 1)));
		ATF_CHECK(FLT_EQ(x[1], -(i + 1)));
	}

	pc_destroy(pc);
}


ATF_TC(point_cloud_calc_diff);
ATF_TC_HEAD(point_cloud_calc_diff, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(point_cloud_calc_diff, tc)
{
	point_cloud *pcin = pc_create();

	for(int i = 1; i <= 30; i++) {
		Point p;
		POINT_ASSIGN_REALS(p, 2*i, -i);
		pc_add_point(pcin, p);
	}

	point_cloud *pcout = pc_create();

	Point x;
	POINT_ASSIGN_REALS(x, -1.0, 2.0);
	pc_calc_diff(pcin, x, pcout);

	for(int i = 0; i < 30; i++) {
		PPoint p;
		p = pc_get_point(pcout, i);
		ATF_CHECK(FLT_EQ(p[0], (2*(i+1) + 1)));
		ATF_CHECK(FLT_EQ(p[1], (-(i + 1) - 2)));
	}

	pc_destroy(pcin);
	pc_destroy(pcout);
}

ATF_TC(point_cloud_equal);
ATF_TC_HEAD(point_cloud_equal, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(point_cloud_equal, tc)
{
	point_cloud *pcin = pc_create();

	for(int i = 1; i <= 30; i++) {
		Point p;
		POINT_ASSIGN_REALS(p, 2*i, -i);
		pc_add_point(pcin, p);
	}

	point_cloud *pcout = pc_create();

	Point x;
	POINT_ASSIGN_REALS(x, 0.0, 0.0);
	pc_calc_diff(pcin, x, pcout);

	ATF_CHECK(pc_equal(pcin, pcout));

	ATF_CHECK(pc_hash_key(pcin) == pc_hash_key(pcout));

	Point pts[2];
	POINT_ASSIGN_REALS(pts[0], -1.0, 1.0);
	POINT_ASSIGN_REALS(pts[1], 2.0, -1.0);
	pc_add_points(pcin, 2, pts);
	DEBUG_INSPECT(pc_hash_key(pcin), %d);
	DEBUG_INSPECT(pc_hash_key(pcout), %d);
	ATF_CHECK(pc_hash_key(pcin) != pc_hash_key(pcout));
	ATF_CHECK(!pc_equal(pcin, pcout));

	pc_destroy(pcin);
	pc_destroy(pcout);
}

ATF_TC(point_cloud_not_equal);
ATF_TC_HEAD(point_cloud_not_equal, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(point_cloud_not_equal, tc)
{
	point_cloud *pc1 = pc_create();
	point_cloud *pc2 = pc_create();
	Point p;
	POINT_ASSIGN_REALS(p, 0.5, 0.5);
	pc_add_point(pc1, p);
	POINT_ASSIGN_REALS(p, 0.5, -0.5);
	pc_add_point(pc2, p);

	ATF_CHECK(!pc_equal(pc1, pc2));

	pc_destroy(pc1);
	pc_destroy(pc2);
}

ATF_TP_ADD_TCS(tp)
{
    ATF_TP_ADD_TC(tp, point_cloud_basic);
    ATF_TP_ADD_TC(tp, point_cloud_calc_diff);
    ATF_TP_ADD_TC(tp, point_cloud_equal);
    ATF_TP_ADD_TC(tp, point_cloud_not_equal);
}
