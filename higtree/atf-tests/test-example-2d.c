#include <string.h>
#include <atf-c.h>
#include "higtree.h"
#include "utils.h"
#define DEBUG
#include "Debug-c.h"

ATF_TC(testname);
ATF_TC_HEAD(testname, tc)
{
    atf_tc_set_md_var(tc, "descr", "Some description...");
}
ATF_TC_BODY(testname, tc)
{
	// Initialize
	int x = 1;
	// Do something;
	x++;
	// Check if everything is ok
	ATF_CHECK(x == 2);
}

ATF_TP_ADD_TCS(tp)
{
    ATF_TP_ADD_TC(tp, testname);
}
