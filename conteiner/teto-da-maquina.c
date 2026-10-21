// CONTROLE DA MAQUINA, nao do codigo.
//
// Sem isto nao da' para dizer se 2,64x em oito nucleos e' o solver que nao
// escala ou a estacao que nao entrega.  Sao duas conclusoes opostas e o numero
// do solver, sozinho, nao distingue.
//
// Dois regimes de proposito:
//   TRIADE  a[i] = b[i] + s*c[i] sobre vetores grandes -- limitado por BANDA DE
//           MEMORIA.  Num notebook de dois canais isto satura muito antes de
//           oito nucleos, e esse limite e' da maquina.
//   NUCLEO  soma de senos num vetor que cabe no cache -- limitado por ULA.
//           Escala quase linear ate' o numero de nucleos FISICOS.
//
// Os dois juntos dao a faixa: o solver nao pode escalar melhor que a triade se
// for limitado por memoria, nem pior que ela sem haver outra causa.
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, np;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    const size_t N = 24u * 1024u * 1024u;     // ~576 MB nos tres vetores
    double *a = malloc(N * sizeof *a);
    double *b = malloc(N * sizeof *b);
    double *c = malloc(N * sizeof *c);
    if (!a || !b || !c) { fprintf(stderr, "sem memoria\n"); MPI_Abort(MPI_COMM_WORLD, 1); }
    for (size_t i = 0; i < N; i++) { a[i] = 0.0; b[i] = 1.0; c[i] = 2.0; }

    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();
    for (int r = 0; r < 10; r++)
        for (size_t i = 0; i < N; i++) a[i] = b[i] + 3.0 * c[i];
    MPI_Barrier(MPI_COMM_WORLD);
    double t_triade = MPI_Wtime() - t0;

    const size_t M = 4096;                     // 32 KB: cabe no L1
    double *d = malloc(M * sizeof *d);
    for (size_t i = 0; i < M; i++) d[i] = i * 1e-6;
    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();
    double s = 0.0;
    for (int r = 0; r < 20000; r++)
        for (size_t i = 0; i < M; i++) s += sin(d[i] + r);
    MPI_Barrier(MPI_COMM_WORLD);
    double t_nucleo = MPI_Wtime() - t0;

    if (rank == 0) {
        // Vazao AGREGADA: o trabalho total cresce com np, entao o tempo
        // constante e' que significa escala perfeita.
        printf("np=%d  triade=%.3f s  nucleo=%.3f s  (soma %.1f)\n",
               np, t_triade, t_nucleo, s);
    }
    free(a); free(b); free(c); free(d);
    MPI_Finalize();
    return 0;
}
