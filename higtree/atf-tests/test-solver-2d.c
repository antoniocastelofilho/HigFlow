#include <string.h>
#include <atf-c.h>
#include "coord.h"
#include "solver.h"
#include "utils.h"
#define DEBUG
#include "Debug-c.h"

ATF_TC(solver_3x3_system);
ATF_TC_HEAD(solver_3x3_system, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(solver_3x3_system, tc)
{
    int procs = 1;
    
    PetscInitializeNoArguments();
    //PetscInitialize(NULL, NULL, (char *)0, "");
    MPI_Comm_size(PETSC_COMM_WORLD, &procs);

    real m[3][3] = {{-1.0,3.0,5.0},{3.0,-5.0,7.0},{-1.0,4.0,-10.0}};
    real b[3] = {1.0,-31.0,35.0};
    int  j[3] = {0,1,2};
    real x[3];
    
    solver *s = slv_create(3,0);
    
    slv_set_maxnonzeros(s, 3);
    slv_set_eps(s, 1e-8);
    
    for (int i = 0; i < 3; i++) {
        
        slv_set_Ai(s, i, 3, j, m[i]);
        slv_set_bi(s, i, b[i]);
    }
    
    slv_assemble(s);
    slv_solve(s);
    
    slv_get_x(s, x);
    
    ATF_CHECK(FLT_EQ(x[0],  1.0));
    ATF_CHECK(FLT_EQ(x[1],  4.0));
    ATF_CHECK(FLT_EQ(x[2], -2.0));
    
	slv_destroy(s);
    
    PetscFinalize();
}

ATF_TC(set_get_mat_vec);
ATF_TC_HEAD(set_get_mat_vec, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(set_get_mat_vec, tc)
{
    int procs = 1;
    
    PetscInitializeNoArguments();
    MPI_Comm_size(PETSC_COMM_WORLD, &procs);

    real m[4][4] = {{-1.0,3.0,5.0,8.0},{-2.0,3.0,-5.0,7.0},{-1.0,4.0,-10.0,-2.0},{2.3,-1.4,2.2,10.2}};
    real b[4] = {-1.0, 2.0, -4.0, 1.0};
    real r[4];
    int  j[4] = {0,1,2,3};
    real v;
    //real x[3];
    
    // **************** Testing slv_set and slv_add
    
    
    solver *s1 = slv_create(4,0);
    
    slv_set_maxiteration(s1, 1234);
    slv_set_maxnonzeros(s1, 4);
    slv_set_eps(s1, 1e-6);

    ATF_CHECK(s1->MaxIteration == 1234 );
    ATF_CHECK(s1->MaxNonZeros  == 4 );
    ATF_CHECK(FLT_EQ(s1->EPS, 1e-6));
    
    // Set A matrix and b vector
    for (int i = 0; i < 4; i++) {
        
        slv_set_Ai(s1, i, 4, j, m[i]);
        slv_set_bi(s1, i, b[i]);
        slv_set_xi(s1, i, b[i]);
    }

    slv_assemble(s1);
    // Compare
    for (int i=0; i < 4; i++) {
        for (int k=0; k < 4; k++) {
            
            MatGetValues(s1->A, 1, &i, 1, &k, &v);
            ATF_CHECK(FLT_EQ(v, m[i][k]));
        }
        
        VecGetValues(s1->b, 1, &i, &v);
        ATF_CHECK(FLT_EQ(v, b[i]));

        VecGetValues(s1->x, 1, &i, &v);
        ATF_CHECK(FLT_EQ(v, b[i]));
    }
    
    solver *s2 = slv_create(4,0);

    slv_set_maxiteration(s2, 4321);
    slv_set_maxnonzeros(s2, 4);
    slv_set_eps(s2, 1e-6);

    ATF_CHECK(s2->MaxIteration == 4321 );
    ATF_CHECK(s2->MaxNonZeros  == 4 );
    ATF_CHECK(FLT_EQ(s2->EPS, 1e-6));

    // Set A matrix and b vector
    for (int i=0; i < 4; i++) {
        
        for (int k=0; k < 4; k++) {
            
            slv_set_Aij(s2, i, k, m[i][k]);
        }
        
        slv_set_bi(s2, i, 0.0);
        slv_add_bi(s2, i, b[i]);
    }
    
    slv_set_x(s2, b);
    
    slv_assemble(s2);
    

    // Compare
    for (int i=0; i < 4; i++) {
        for (int k=0; k < 4; k++) {
            
            MatGetValues(s2->A, 1, &i, 1, &k, &v);
    
            //DEBUG_INSPECT(v,%e);
            //DEBUG_INSPECT(m[i][k],%e);

            ATF_CHECK(FLT_EQ(v, m[i][k]));
        }
        
        VecGetValues(s2->b, 1, &i, &v);
        ATF_CHECK(FLT_EQ(v, b[i]));

        VecGetValues(s2->x, 1, &i, &v);
        ATF_CHECK(FLT_EQ(v, b[i]));
    }
    
    // **************** Testing slv_get
    
    for (int i = 0; i < 4; i++) {
       
        v = slv_get_xi(s2,i);
        ATF_CHECK(FLT_EQ(v, b[i]));
    }
    
    slv_get_x(s2, r);
 
    for (int i = 0; i < 4; i++) {
       
        ATF_CHECK(FLT_EQ(r[i], b[i]));
    }

    
    DEBUG_PASS;
    solver *s3 = slv_create(100,0);
    
    int k = 10;
    slv_set_Ai(s3, k, 4, j, m[0]);
    
	DEBUG_PASS;
	slv_assemble(s3);
    
	DEBUG_PASS;
    MatGetValues(s3->A, 1, &k, 4, j, r);
        
    for (int i = 0; i < 4; i++) {
        
        ATF_CHECK(FLT_EQ(r[i],m[0][i]));
    }
    
    //s3->MaxNonZeros = 0;

    //MatDestroy(&s3->A);
    MatZeroEntries(s3->A);

    for (int i = 0; i < 4; i++) {
        
        slv_set_Aij(s3, k, i, m[0][i]);
    }
    
    DEBUG_PASS;
    slv_assemble(s3);

    MatGetValues(s3->A, 1, &k, 4, j, r);

    for (int i = 0; i < 4; i++) {
        
        ATF_CHECK(FLT_EQ(r[i],m[0][i]));
    }
    
    ATF_CHECK(FLT_EQ(s3->MaxNonZeros, 1000));    
    
    slv_set_maxnonzeros(s3,100);
    
    DEBUG_PASS;
    slv_destroy(s1);
    slv_destroy(s2);
    slv_destroy(s3);
    
    PetscFinalize();
}

ATF_TC(solver_sparse_lu);
ATF_TC_HEAD(solver_sparse_lu, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(solver_sparse_lu, tc)
{
    int procs = 1;
    int argc = 5;
    char **argv;
    
    argv = malloc(sizeof(char *)*10);
    
    for (int i = 0; i < 5; i++)
        argv[i] = malloc(sizeof(char) * 20);
    
    strcpy(argv[0],"./test-solver");
    strcpy(argv[1],"-ksp_type");
    strcpy(argv[2],"preonly");
    strcpy(argv[3],"-pc_type");
    strcpy(argv[4],"lu");
    
    PetscInitialize(&argc, &argv, PETSC_NULL, "HELP");
    
    MPI_Comm_size(PETSC_COMM_WORLD, &procs);

    int N = 100;
    int col[3];
    real x[N];
    real val[3] = {-1.0, 2.0, -1.0};
    int ncols = 3;
    
    solver *s = slv_create(N,0);

    slv_set_maxiteration(s, 1000);
    slv_set_maxnonzeros(s, 3);
    slv_set_eps(s, 1e-6);

    for (int i = 0; i < N; i++) {
        
        if (i==0) {
            col[0] = i;
            col[1] = i+1;
            ncols = 2;
            slv_set_bi(s, i, 1.0);
        }
        else if (i == N-1) {
            col[0] = i-1;
            col[1] = i;
            ncols = 2;
            slv_set_bi(s, i, 1.0);
        }
        else {
            col[0] = i-1;
            col[1] = i;
            col[2] = i+1;
            ncols = 3;
            slv_set_bi(s, i, 0.0);
        }
        
        slv_set_Ai(s, i, ncols, col, val);
    }
    
    slv_assemble(s);
    slv_solve(s);
    
    slv_get_x(s, x);
    
    for (int i = 0; i < N; i++) {
        ATF_CHECK(FLT_EQ(x[i],  1.0));
        DEBUG_INSPECT(x[i],%e);
    }
    
	slv_destroy(s);
    
    PetscFinalize();
}

ATF_TC(solver_sparse_gmres);
ATF_TC_HEAD(solver_sparse_gmres, tc)
{
    atf_tc_set_md_var(tc, "descr", "...");
}
ATF_TC_BODY(solver_sparse_gmres, tc)
{
    int procs = 1;
    
    int argc = 5;
    // not working:
    //char *argv[] = {"./test-solver", "-ksp_type", "gmres", "-pc_type", "ilu"};
    char **argv;
    
    argv = malloc(sizeof(char *)*10);
    
    for (int i = 0; i < 5; i++)
        argv[i] = malloc(sizeof(char) * 20);
    
    strcpy(argv[0],"./test-solver");
    strcpy(argv[1],"-ksp_type");
    strcpy(argv[2],"gmres");
    strcpy(argv[3],"-pc_type");
    strcpy(argv[4],"ilu");
    
    
    PetscInitialize(&argc, &argv, PETSC_NULL, "HELP");
    
    //PetscInitialize(NULL, NULL, PETSC_NULL, "");
    MPI_Comm_size(PETSC_COMM_WORLD, &procs);

    int N = 100;
    int col[3];
    real x[N];
    real val[3] = {-1.0, 2.0, -1.0};
    int ncols = 3;
    
    solver *s = slv_create(N,0);
    
    slv_set_maxiteration(s, 1000);
    slv_set_maxnonzeros(s, 3);
    slv_set_eps(s, 1e-6);

    for (int i = 0; i < N; i++) {
        
        if (i==0) {
            col[0] = i;
            col[1] = i+1;
            ncols = 2;
            slv_set_bi(s, i, 1.0);
        }
        else if (i == N-1) {
            col[0] = i-1;
            col[1] = i;
            ncols = 2;
            slv_set_bi(s, i, 1.0);
        }
        else {
            col[0] = i-1;
            col[1] = i;
            col[2] = i+1;
            ncols = 3;
            slv_set_bi(s, i, 0.0);
        }
        
        slv_set_Ai(s, i, ncols, col, val);
    }
    
    slv_assemble(s);
    slv_solve(s);
    
    slv_get_x(s, x);
    
    for (int i = 0; i < N; i++) {
        ATF_CHECK(x[i] < 1.0 + 1e-3);
        ATF_CHECK(x[i] > 1.0 - 1e-3);
        //DEBUG_INSPECT(x[i],%e);
    }
    
	slv_destroy(s);
    
    PetscFinalize();
}

ATF_TP_ADD_TCS(tp)
{
    ATF_TP_ADD_TC(tp, solver_3x3_system);
    ATF_TP_ADD_TC(tp, set_get_mat_vec);
    ATF_TP_ADD_TC(tp, solver_sparse_lu);
//    ATF_TP_ADD_TC(tp, solver_sparse_gmres); // Fix it!
}
