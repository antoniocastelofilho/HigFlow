#include <string.h>
#include <atf-c.h>
#include "higtree.h"
#include "higtree-iterator.h"
#include "higtree-io.h"
#include "utils.h"
#define DEBUG
#include "Debug-c.h"

static void __refine_grid(hig_cell *root, int nc[DIM]) {
	int numChildren = hig_get_number_of_children(root);
	for(int i = 0; i < numChildren; i++) {
		int pp[DIM];
		hig_tobase(i, root->numcells, pp);
		if ((pp[0] + pp[1])%2==0) {
		  hig_refine_uniform(hig_get_child(root, i), nc);
		}
	}
}

void __check_cells_n_levels(higcit_celliterator *it, hig_cell *cell, int pos[], int num, int levels) {
	for(int i = 0, ind = 0; i < num; i++, ind+=levels) {
		hig_cell *c = cell;
		for(int j = 0; j < levels; j++) {
			if (pos[ind+j] >= 0) {
				c = hig_get_child(c, pos[ind+j]);
			} else {
				break;
			}
		}
		hig_cell *ic = higcit_getcell(it);
		ATF_CHECK(ic == c);
		higcit_nextcell(it);
	}
	ATF_CHECK(higcit_isfinished(it));
}

ATF_TC(higio_print);
ATF_TC_HEAD(higio_print, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(higio_print, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 2);
	hig_refine_uniform(root, nc);
	
	char tmpfilename[L_tmpnam];
	tmpnam(tmpfilename);
	FILE *fd = fopen(tmpfilename, "w");
	higio_print(fd, root);
	fclose(fd);

	fd = fopen(tmpfilename, "r");
	char buf[1024];
	fread(buf, 1024, sizeof(char), fd);
	fclose(fd);

	char result[] = 
"(-1.000000000000000, -1.000000000000000, 0.000000000000000, 0.000000000000000 posincell:  0 )\n(-1.000000000000000, 0.000000000000000, 0.000000000000000, 1.000000000000000 posincell:  1 )\n(0.000000000000000, -1.000000000000000, 1.000000000000000, 0.000000000000000 posincell:  2 )\n(0.000000000000000, 0.000000000000000, 1.000000000000000, 1.000000000000000 posincell:  3 )";
	printf("%s", buf);
	ATF_CHECK(strncmp(buf, result,strlen(result)) == 0);

	hig_destroy(root);
}

ATF_TC(higio_write_leaves);
ATF_TC_HEAD(higio_write_leaves, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(higio_write_leaves, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 2);
	hig_refine_uniform(root, nc);
	
	char tmpfilename[L_tmpnam];
	tmpnam(tmpfilename);
	FILE *fd = fopen(tmpfilename, "w");
	higio_write_leaves(fd, root);
	fclose(fd);

	fd = fopen(tmpfilename, "r");
	char buf[1024];
	fread(buf, 1024, sizeof(char), fd);
	fclose(fd);

	char result[] = 
"-1.00000 -1.00000 0.00000 0.00000\n"
"-1.00000 0.00000 0.00000 1.00000\n"
"0.00000 -1.00000 1.00000 0.00000\n"
"0.00000 0.00000 1.00000 1.00000";
	printf("%s", buf);
	ATF_CHECK(strncmp(buf, result,strlen(result)) == 0);

	hig_destroy(root);
}

ATF_TC(higio_read_from_amr2d);
ATF_TC_HEAD(higio_read_from_amr2d, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(higio_read_from_amr2d, tc)
{
	char tmpfilename[L_tmpnam];
	tmpnam(tmpfilename);
	FILE *fd = fopen(tmpfilename, "w");
	char inp[] =
"0 0.9 0 0.8\n"
"3\n"
"0.1 0.1 1\n"
"1 1 9 8\n"
"0.05 0.05 2\n"
"5 3 4 10\n"
"9 5 2 6\n"
"0.025 0.025 1 \n"
"13 13 6 4\n";
	fwrite(inp, strlen(inp), sizeof(char), fd);
	fclose(fd);
	fd = fopen(tmpfilename, "r");
	hig_cell *root = higio_read_from_amr2d(fd);
	fclose(fd);

	higcit_celliterator *it = higcit_create_all_leaves(root);
	int pos[] = {
		0,-1,-1,
		1,-1,-1,
		2,-1,-1,
		3,-1,-1,
		4,-1,-1,
		5,-1,-1,
		6,-1,-1,
		7,-1,-1,
		8,-1,-1,
		9,-1,-1,
		10,-1,-1,
		11,-1,-1,
		12,-1,-1,
		13,-1,-1,
		14,-1,-1,
		15,-1,-1,
		16,-1,-1,
		17,0,-1,
		17,1,-1,
		17,2,-1,
		17,3,-1,
		18,0,-1,
		18,1,-1,
		18,2,-1,
		18,3,-1,
		19,0,-1,
		19,1,-1,
		19,2,-1,
		19,3,-1,
		20,0,-1,
		20,1,-1,
		20,2,-1,
		20,3,-1,
		21,0,-1,
		21,1,-1,
		21,2,-1,
		21,3,-1,
		22,-1,-1,
		23,-1,-1,
		24,-1,-1,
		25,0,-1,
		25,1,-1,
		25,2,-1,
		25,3,-1,
		26,0,-1,
		26,1,-1,
		26,2,-1,
		26,3,-1,
		27,0,0,
		27,0,1,
		27,0,2,
		27,0,3,
		27,1,0,
		27,1,1,
		27,1,2,
		27,1,3,
		27,2,0,
		27,2,1,
		27,2,2,
		27,2,3,
		27,3,0,
		27,3,1,
		27,3,2,
		27,3,3,
		28,0,-1,
		28,1,-1,
		28,2,-1,
		28,3,-1,
		29,0,-1,
		29,1,-1,
		29,2,-1,
		29,3,-1,
		30,-1,-1,
		31,-1,-1,
		32,-1,-1,
		33,-1,-1,
		34,0,-1,
		34,1,-1,
		34,2,-1,
		34,3,-1,
		35,0,0,
		35,0,1,
		35,0,2,
		35,0,3,
		35,1,0,
		35,1,1,
		35,1,2,
		35,1,3,
		35,2,-1,
		35,3,-1,
		36,0,-1,
		36,1,-1,
		36,2,-1,
		36,3,-1,
		37,-1,-1,
		38,-1,-1,
		39,-1,-1,
		40,-1,-1,
		41,-1,-1,
		42,-1,-1,
		43,-1,-1,
		44,-1,-1,
		45,-1,-1,
		46,-1,-1,
		47,-1,-1,
		48,-1,-1,
		49,-1,-1,
		50,-1,-1,
		51,-1,-1,
		52,-1,-1,
		53,-1,-1,
		54,-1,-1,
		55,-1,-1,
		56,-1,-1,
		57,-1,-1,
		58,-1,-1,
		59,-1,-1,
		60,-1,-1,
		61,-1,-1,
		62,-1,-1,
		63,-1,-1,
		64,-1,-1,
		65,-1,-1,
		66,-1,-1,
		67,-1,-1,
		68,-1,-1,
		69,-1,-1,
		70,-1,-1,
		71,-1,-1,
	};

	__check_cells_n_levels(it, root, pos, 129, 3);
	higcit_destroy(it);

	hig_destroy(root);
}


ATF_TC(higio_print_in_vtk);
ATF_TC_HEAD(higio_print_in_vtk, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(higio_print_in_vtk, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 2);
	hig_refine_uniform(root, nc);
	
	char tmpfilename[L_tmpnam];
	tmpnam(tmpfilename);
	FILE *fd = fopen(tmpfilename, "w");
	higio_print_in_vtk2d(fd, root);
	fclose(fd);

	fd = fopen(tmpfilename, "r");
	char buf[1024];
	fread(buf, 1024, sizeof(char), fd);
	fclose(fd);

	char result[] = 
"# vtk DataFile Version 3.0\n"
"higtree\n"
"ASCII\n"
"DATASET POLYDATA\n"
"POINTS 16 float\n"
"-1.000000 -1.000000 0\n"
"-1.000000 0.000000 0\n"
"0.000000 0.000000 0\n"
"0.000000 -1.000000 0\n"
"-1.000000 0.000000 0\n"
"-1.000000 1.000000 0\n"
"0.000000 1.000000 0\n"
"0.000000 0.000000 0\n"
"0.000000 -1.000000 0\n"
"0.000000 0.000000 0\n"
"1.000000 0.000000 0\n"
"1.000000 -1.000000 0\n"
"0.000000 0.000000 0\n"
"0.000000 1.000000 0\n"
"1.000000 1.000000 0\n"
"1.000000 0.000000 0\n"
"POLYGONS 4 20\n"
"4 0 1 2 3\n"
"4 4 5 6 7\n"
"4 8 9 10 11\n"
"4 12 13 14 15\n";

	printf("%s", buf);
	ATF_CHECK(strncmp(buf, result,strlen(result)) == 0);

	hig_destroy(root);
}

ATF_TC(write_read_file);
ATF_TC_HEAD(write_read_file, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(write_read_file, tc)
{
	Point lp, hp;
	POINT_ASSIGN_SCALAR(lp, -1.0);
	POINT_ASSIGN_SCALAR(hp, 1.0);
	hig_cell *root = hig_create_root(lp, hp);
	int nc[DIM];
	POINT_ASSIGN_SCALAR(nc, 4);
	hig_refine_uniform(root, nc);
	POINT_ASSIGN_SCALAR(nc, 2);
	__refine_grid(root, nc);
	
	char tmpfilename[L_tmpnam];
	tmpnam(tmpfilename);
	FILE *fd = fopen(tmpfilename, "w");
	higio_write_to_file(fd, root);
	fclose(fd);

	fd = fopen(tmpfilename, "r");
	hig_cell *r = higio_read_from_file(fd);
	fclose(fd);

	higcit_celliterator *it1, *it2;
	for(it1 = higcit_create_all_leaves(root), it2 = higcit_create_all_leaves(r); !higcit_isfinished(it1) && !higcit_isfinished(it2); higcit_nextcell(it1), higcit_nextcell(it2)) {
		hig_cell *c1 = higcit_getcell(it1);
		hig_cell *c2 = higcit_getcell(it2);
		Point p1, p2;
		hig_get_lowpoint(c1, p1);
		hig_get_lowpoint(c2, p2);
		ATF_CHECK(p1[0] == p2[0]);
		ATF_CHECK(p1[1] == p2[1]);
		hig_get_highpoint(c1, p1);
		hig_get_highpoint(c2, p2);
		ATF_CHECK(p1[0] == p2[0]);
		ATF_CHECK(p1[1] == p2[1]);
	}
	higcit_destroy(it1);
	higcit_destroy(it2);
	hig_destroy(root);
	hig_destroy(r);
}


ATF_TP_ADD_TCS(tp)
{
    ATF_TP_ADD_TC(tp, write_read_file);
    ATF_TP_ADD_TC(tp, higio_print);
    ATF_TP_ADD_TC(tp, higio_print_in_vtk);
    ATF_TP_ADD_TC(tp, higio_write_leaves);
    ATF_TP_ADD_TC(tp, higio_read_from_amr2d);
}

