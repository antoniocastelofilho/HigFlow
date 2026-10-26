// Carga de trabalho FIXO, limitada por ULA, para medir a maquina e nao o codigo.
// Trabalho identico em toda invocacao: a duracao E' a velocidade.
#include <mpi.h>
#include <stdio.h>
#include <math.h>
int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double d[4096];
    for (int i = 0; i < 4096; i++) d[i] = i * 1e-6;
    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime(), s = 0.0;
    for (int r = 0; r < 300000; r++)
        for (int i = 0; i < 4096; i++) s += sin(d[i] + r);
    MPI_Barrier(MPI_COMM_WORLD);
    double t = MPI_Wtime() - t0;
    if (rank == 0) printf("%.3f %.1f\n", t, s);
    MPI_Finalize();
    return 0;
}
