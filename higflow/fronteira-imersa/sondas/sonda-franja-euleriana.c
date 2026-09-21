// Sonda: quantas celulas para FORA do que o rank possui a localizacao de ponto
// ainda resolve?  E' a largura de franja EULERIANA util, e e' ela que limita o
// suporte do delta regularizado -- nao a sobreposicao do DMPlex.
//
// O HiGFlow declara pg_set_fringe_size(pg, 5) em hig-flow-kernel.c:1465.  Isto
// mede o ALCANCE REAL, que e' outra coisa: constante declarada nao e' geometria
// entregue.
//
//   Roma  (3 pontos)   suporte 1,5 celulas para cada lado
//   Peskin(4 pontos)   suporte 2,0 celulas para cada lado
#include "higtree.h"
#include "higtree-io.h"
#include "higtree-iterator.h"
#include "pdomain.h"
#include "lbal.h"
#include <mpi.h>
#include <stdio.h>
#include <string.h>

#define NC       24      // celulas por direcao
#define MAXPASSO  8      // ate' onde tentar alcancar

int main(int argc, char *argv[]) {
    higtree_initialize(&argc, &argv);
    int rank, ntasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    const unsigned franjas[] = {1, 2, 5};
    for (int fi = 0; fi < 3; fi++) {
        const unsigned FRANJA = franjas[fi];

        partition_graph *pg = pg_create(MPI_COMM_WORLD);
        pg_set_fringe_size(pg, FRANJA);
        load_balancer *lb = lb_create(MPI_COMM_WORLD, 1);
        if (rank == 0) {
            Point lo, hi;
            POINT_ASSIGN_SCALAR(lo, 0.0);
            POINT_ASSIGN_SCALAR(hi, 1.0);
            hig_cell *raiz = hig_create_root(lo, hi);
            int nc[DIM];
            for (int d = 0; d < DIM; d++) nc[d] = NC;
            hig_refine_uniform(raiz, nc);
            lb_add_input_tree(lb, raiz, true, 0);
        }
        lb_calc_partition(lb, pg);

        sim_domain *sd = sd_create(NULL);
        sd_set_interpolator_order(sd, 2);
        for (unsigned i = 0; i < lb_get_num_local_trees(lb); ++i)
            sd_add_higtree(sd, lb_get_local_tree(lb, i, NULL));
        lb_destroy(lb);

        psim_domain *psd = psd_create(sd, pg);
        psd_synced_mapper(psd);

        // De cada celula PROPRIA, andar para fora em cada direcao ate' a
        // localizacao falhar.  So' contam os passos que continuam DENTRO do
        // dominio global -- fora dele o limite e' o dominio, nao a franja.
        const double h = 1.0 / NC;
        int alcance_min = MAXPASSO, alcance_max = 0, amostras = 0;

        const unsigned nproprias = sd_get_num_local_higtrees(sd);
        for (unsigned t = 0; t < nproprias; t++) {
            higcit_celliterator *it;
            for (it = higcit_create_all_leaves(sd_get_higtree(sd, t));
                 !higcit_isfinished(it); higcit_nextcell(it)) {
                hig_cell *c = higcit_getcell(it);
                Point cc;
                hig_get_center(c, cc);
                for (int d = 0; d < DIM; d++) {
                    for (int s = -1; s <= 1; s += 2) {
                        int k = 0;
                        for (; k < MAXPASSO; k++) {
                            Point p;
                            memcpy(p, cc, sizeof(Point));
                            p[d] += s * (k + 1) * h;
                            if (p[d] <= 0.0 || p[d] >= 1.0) { k = -1; break; }
                            if (sd_get_cell_with_point(sd, p) == NULL) break;
                        }
                        if (k < 0) continue;              // barrou no dominio
                        if (k < alcance_min) alcance_min = k;
                        if (k > alcance_max) alcance_max = k;
                        amostras++;
                    }
                }
            }
            higcit_destroy(it);
        }

        printf("  franja %u | rank %d | alcance para fora: min %d, max %d celulas "
               "(%d direcoes medidas)  Roma 1,5 %s  Peskin 2,0 %s\n",
               FRANJA, rank, alcance_min, alcance_max, amostras,
               alcance_min >= 2 ? "cabe" : "NAO CABE",
               alcance_min >= 2 ? "cabe" : "NAO CABE");
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);

        psd_destroy(psd);
        pg_destroy(pg);
    }

    // Sem MPI_Finalize de proposito: o higtree_initialize inicializa o PETSc
    // (utils.c:44) e o encerramento dele roda DEPOIS, com o MPI ja' derrubado --
    // "Local abort after MPI_FINALIZE started".  Os testes da suite tambem nao o
    // chamam, e e' por isso.
    return 0;
}
