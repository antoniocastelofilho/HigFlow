// A PARTICAO FINAL SENDO A DO T8CODE, e o conteudo da franja como oraculo.
//
// Ate' aqui o t8code produzia a malha e o `lbal` a reparticionava.  Aqui nao ha'
// `lb_calc_partition`: cada rank fica com o que a curva de preenchimento lhe deu,
// e o `partition_graph` e' montado a partir das caixas.
//
// O ORACULO E' O DO test-fringe-sync, e de proposito: a construcao do grafo falha
// EM SILENCIO se o pareamento entre quem envia e quem recebe sair trocado -- o
// dominio se monta, o sync roda, e a celula de franja recebe o valor de outra
// posicao.  Contagem de celulas e volume nao pegam isso; so' o conteudo pega.
//
//   dominio_esta_inteiro          a soma global das celulas locais e' a malha.
//   ha_franja_para_conferir       com np>1 ha' franja; sem ela o caso seguinte
//                                 nao afirmaria nada.
//   franja_recebe_o_valor_da_propria_posicao
//                                 depois do `dp_sync`, cada celula de franja tem
//                                 de conter o campo no PROPRIO centro.  Os
//                                 coeficientes diferem por direcao para que o
//                                 valor identifique a posicao.
//   estencil_atravessa_a_franja   o suporte de uma celula na borda alcanca celula
//                                 de franja -- a franja serve ao que existe para
//                                 servir, e nao so' existe.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "higtree.h"
#include "higtree-iterator.h"
#include "pdomain.h"
#include "utils.h"
#include "testing.h"
#include "t8code/t8-particao-grafo.h"

#define NX 32
#define NY 16
#define FRINGE 2

static real campo(const Point p) {
    return 1.0 + 2.0 * p[0] + 30.0 * p[1] + (DIM > 2 ? 400.0 * p[DIM - 1] : 0.0);
}

int main(int argc, char *argv[]) {
    higtree_initialize(&argc, &argv);

    int rank, ntasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    Point lo, hi;
    lo[0] = 0.0; lo[1] = -1.0;
    hi[0] = 4.0; hi[1] =  1.0;
    int nb[DIM];
    nb[0] = NX; nb[1] = NY;
#if DIM == 3
    lo[2] = 0.0; hi[2] = 2.0; nb[2] = 8;
#endif

    partition_graph *pg = pg_create(MPI_COMM_WORLD);
    pg_set_fringe_size(pg, FRINGE);

    sim_domain *sd = sd_create(NULL);
    sd_set_interpolator_order(sd, 2);

    const long ok = t8_monta_dominio_particionado(lo, hi, nb, sd, pg) ? 1 : 0;

    psim_domain *psd = psd_create(sd, pg);
    psd_synced_mapper(psd);
    mp_mapper *m = sd_get_domain_mapper(sd);
    distributed_property *dp = psd_create_property(psd);

    long n_local = 0;
    higcit_celliterator *cit;
    for (cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit);
         higcit_nextcell(cit)) {
        hig_cell *c = higcit_getcell(cit);
        Point ce;
        hig_get_center(c, ce);
        dp_set_value(dp, mp_lookup(m, hig_get_cid(c)), campo(ce));
        n_local++;
    }
    higcit_destroy(cit);

    dp_sync(dp);

    long n_franja = 0, erradas = 0;
    real pior = 0.0;
    char primeira[256]; primeira[0] = '\0';
    for (unsigned k = 0; k < sd_get_num_fringe_higtrees(sd); ++k) {
        higcit_celliterator *fit = higcit_create_all_leaves(sd_get_fringe_higtree(sd, k));
        for (; !higcit_isfinished(fit); higcit_nextcell(fit)) {
            hig_cell *c = higcit_getcell(fit);
            const int lid = mp_lookup(m, hig_get_cid(c));
            if (lid < 0) continue;
            Point ce;
            hig_get_center(c, ce);
            const real esp = campo(ce), obt = dp_get_value(dp, lid);
            const real e = fabs(obt - esp);
            n_franja++;
            if (e > 1.0e-12) {
                if (!erradas)
                    snprintf(primeira, sizeof primeira,
                        "franja em (%.4f, %.4f): esperado %.10g, obtido %.10g",
                        (double) ce[0], (double) ce[1], (double) esp, (double) obt);
                erradas++;
                if (e > pior) pior = e;
            }
        }
        higcit_destroy(fit);
    }

    // O estencil de uma celula alcanca a franja?
    long atravessa = 0;
    sim_stencil *stn = stn_create();
    for (cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit) && !atravessa;
         higcit_nextcell(cit)) {
        hig_cell *c = higcit_getcell(cit);
        // Ponto DESLOCADO do centro: no centro exato o `cell_find_in_center`
        // atalha e o estencil sai com a propria celula so' -- o caso passaria a
        // afirmar nada.
        Point ce, alvo, de;
        hig_get_center(c, ce);
        hig_get_delta(c, de);
        POINT_ASSIGN(alvo, ce);
        alvo[0] = ce[0] + 0.25 * de[0];
        stn_reset(stn);
        sd_get_stencil(sd, alvo, alvo, 1.0, stn);
        const int n = stn_get_numelems(stn);
        const int *ids = stn_get_ids(stn);
        for (int i = 0; i < n; i++)
            if (ids[i] >= (int) n_local) { atravessa = 1; break; }
    }
    higcit_destroy(cit);
    stn_destroy(stn);

    long n_local_g, n_franja_g, erradas_g, ok_min, atravessa_min;
    real pior_g;
    MPI_Allreduce(&n_local,   &n_local_g,   1, MPI_LONG,   MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&n_franja,  &n_franja_g,  1, MPI_LONG,   MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&erradas,   &erradas_g,   1, MPI_LONG,   MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&pior,      &pior_g,      1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&ok,        &ok_min,      1, MPI_LONG,   MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&atravessa, &atravessa_min,1, MPI_LONG,  MPI_MIN, MPI_COMM_WORLD);

    if (rank != 0) return 0;

    long esperado = 1;
    for (int d = 0; d < DIM; d++) esperado *= nb[d];

    t_case("montagem_nao_falha");
    T_CHECK_MSG(ok_min == 1, "t8_monta_dominio_particionado falhou em algum rank");

    t_case("dominio_esta_inteiro");
    T_CHECK_MSG(n_local_g == esperado,
        "soma global das celulas locais = %ld, esperado %ld (np=%d)",
        n_local_g, esperado, ntasks);

    t_case("ha_franja_para_conferir");
    if (ntasks == 1) {
        T_CHECK_MSG(n_franja_g == 0, "com np=1 nao deveria haver franja");
    } else {
        T_CHECK_MSG(n_franja_g > 0,
            "com np=%d nao ha' franja nenhuma -- o caso seguinte nao afirmaria "
            "nada", ntasks);
    }

    t_case("franja_recebe_o_valor_da_propria_posicao");
    T_CHECK_MSG(erradas_g == 0,
        "%ld de %ld celula(s) de franja com valor de outra posicao (pior %.6g).  "
        "%s", erradas_g, n_franja_g, (double) pior_g, primeira);

    t_case("estencil_atravessa_a_franja");
    if (ntasks == 1) {
        T_CHECK_MSG(atravessa_min == 0, "com np=1 nao ha' franja a atravessar");
    } else {
        T_CHECK_MSG(atravessa_min == 1,
            "com np=%d todo rank deveria ter celula cujo suporte alcanca a "
            "franja", ntasks);
    }

    return t_end();
}
