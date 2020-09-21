#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include <mpi.h>

#include "solver.h"

#define DEBUG
#include "Debug-c.h"

int main(int argc, char *argv[]) {

    clock_t start, stop;
    real cpu_time;
    int procs = 1;

    DEBUG_INSPECT(argc,%d);

    for (int i = 0; i < argc; i++) {

        DEBUG_INSPECT(argv[i],%s);

    }

    higtree_initialize(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

    int N = 100;
    int col[3];
    real x[N];
    real val[3] = {-1.0, 2.0, -1.0};
    int ncols = 3;

    solver *s = slv_create(SOLVER_ANY, 0, N);

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
        DEBUG_INSPECT(x[i],%e);
    }

	slv_destroy(s);









/*
    int N = 1e+6;
    int col[100];
    real val[100];
    real v;

    for (int j = 0; j < 100; j++) {
        col[j] = random() % N;
    }

    solver *s1 = slv_create_petsc(N,100,1000,1e-8);
    solver *s2 = slv_create_petsc(N,100,1000,1e-8);
    solver *s3 = slv_create_petsc(N,100,1000,1e-8);

    start = clock();

    for (int i = 0; i < N; i++) {

        for (int j = 0; j < 100; j++) {

            v = (real) (random() % 6);

            MatSetValues(s1->A, 1, &i, 1, &col[j], &v, INSERT_VALUES);
        }
    }

    MatAssemblyBegin(s1->A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd  (s1->A, MAT_FINAL_ASSEMBLY);

    stop = clock();

    cpu_time = (stop - start)/CLOCKS_PER_SEC;

    DEBUG_INSPECT("Usando MatSetValues",%s);
    DEBUG_INSPECT(cpu_time, %f);

    start = clock();

    for (int i = 0; i < N; i++) {

        for (int j = 0; j < 100; j++) {

            v = (real) (random() % 6);

            slv_set_Aij(s2, i, col[j], v);
        }
    }
    slv_set_Aij(s2,0,0,0); // ************* Flush


    MatAssemblyBegin(s2->A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd  (s2->A, MAT_FINAL_ASSEMBLY);

    stop = clock();

    cpu_time = (stop - start)/CLOCKS_PER_SEC;

    DEBUG_INSPECT("Usando slv_set_Aij",%s);
    DEBUG_INSPECT(cpu_time, %f);


    start = clock();

    for (int i = 0; i < N; i++) {

        for (int j = 0; j < 100; j++) {

            val[j] = (real) (random() % 6);
        }

        slv_set_Ai(s3, i, 100, col, val);
    }

    MatAssemblyBegin(s3->A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd  (s3->A, MAT_FINAL_ASSEMBLY);

    stop = clock();

    cpu_time = (stop - start)/CLOCKS_PER_SEC;

    DEBUG_INSPECT("Usando slv_set_Ai",%s);
    DEBUG_INSPECT(cpu_time, %f);

    slv_destroy(s1);
    slv_destroy(s2);
    slv_destroy(s3);
*/

    return 0;

}
