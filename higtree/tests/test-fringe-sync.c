// A TROCA DE FRANJA ENTREGA O VALOR CERTO?
//
// A franja ja' tinha dois testes: o test-fringe-support afirma o CONTRATO dela
// (o suporte a alcanca) e o test-fringe-parallel afirma que ela EXISTE sob
// particionamento real.  Nenhum dos dois olha o CONTEUDO: se o pareamento entre
// quem envia e quem recebe estivesse trocado, a celula de franja receberia o
// valor de OUTRA celula, e os dois continuariam verdes.
//
// POR QUE ISSO IMPORTA AGORA.  O `partition_graph` guarda, por vizinho, faixas
// retangulares a enviar e arvores a receber, e o cabecalho e' explicito: "must be
// strictly on the same order as the corresponding list on the remote process".
// Ordem trocada nao quebra nada visivelmente -- entrega dado errado.  Qualquer
// tentativa de construir esse grafo fora do `lbal` (particionar pelo t8code, por
// exemplo) precisa deste oraculo ANTES de comecar, nao depois.
//
// O ORACULO E' ANALITICO e nao depende de fisica: cada celula LOCAL recebe
// f(centro) de um campo conhecido; depois do `dp_sync`, cada celula de FRANJA tem
// de conter f(centro DELA).  Se o pareamento estiver trocado, o valor lido e' o de
// outra posicao e a diferenca e' grande -- nao e' arredondamento.
//
// VALIDADO NO SENTIDO INVERSO: escrevendo f(centro) + 1 nas celulas locais de um
// rank so', o caso REPROVA nos vizinhos dele -- o que confirma que ele le' o dado
// que atravessou, e nao o que o proprio rank escreveu.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "higtree.h"
#include "higtree-iterator.h"
#include "pdomain.h"
#include "lbal.h"
#include "utils.h"
#include "testing.h"

#define NC     8
#define FRINGE 2

static real campo(const Point p) {
    // Linear em cada direcao, com coeficientes distintos: assim o valor
    // identifica a POSICAO, e trocar duas celulas muda o numero.
    return 1.0 + 2.0 * p[0] + 30.0 * p[1] + (DIM > 2 ? 400.0 * p[DIM - 1] : 0.0);
}

int main(int argc, char *argv[]) {
    higtree_initialize(&argc, &argv);

    int rank, ntasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    partition_graph *pg = pg_create(MPI_COMM_WORLD);
    pg_set_fringe_size(pg, FRINGE);

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
    mp_mapper *m = sd_get_domain_mapper(sd);

    distributed_property *dp = psd_create_property(psd);

    // So' as celulas LOCAIS sao escritas.  A franja fica por conta do sync.
    long n_local = 0;
    higcit_celliterator *cit;
    for (cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit);
         higcit_nextcell(cit)) {
        hig_cell *c = higcit_getcell(cit);
        const int lid = mp_lookup(m, hig_get_cid(c));
        Point ce;
        hig_get_center(c, ce);
        dp_set_value(dp, lid, campo(ce));
        n_local++;
    }
    higcit_destroy(cit);

    dp_sync(dp);

    // Agora a franja: cada celula tem de conter o campo no PROPRIO centro.
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
            const real esperado = campo(ce);
            const real obtido = dp_get_value(dp, lid);
            const real e = fabs(obtido - esperado);
            n_franja++;
            if (e > 1.0e-12) {
                if (!erradas) {
                    snprintf(primeira, sizeof primeira,
                        "celula de franja em (%.6f, %.6f): esperado %.12g, obtido "
                        "%.12g", (double) ce[0], (double) ce[DIM > 1 ? 1 : 0],
                        (double) esperado, (double) obtido);
                }
                erradas++;
                if (e > pior) pior = e;
            }
        }
        higcit_destroy(fit);
    }

    long n_franja_g, erradas_g, n_local_g;
    real pior_g;
    MPI_Allreduce(&n_franja, &n_franja_g, 1, MPI_LONG,   MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&erradas,  &erradas_g,  1, MPI_LONG,   MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&n_local,  &n_local_g,  1, MPI_LONG,   MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&pior,     &pior_g,     1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (rank != 0) return 0;

    long esperado_total = 1;
    for (int d = 0; d < DIM; d++) esperado_total *= NC;

    t_case("dominio_esta_inteiro");
    T_CHECK_MSG(n_local_g == esperado_total,
        "soma global das celulas locais = %ld, esperado %ld", n_local_g,
        esperado_total);

    t_case("ha_franja_para_conferir");
    if (ntasks == 1) {
        T_CHECK_MSG(n_franja_g == 0, "com np=1 nao deveria haver franja");
    } else {
        T_CHECK_MSG(n_franja_g > 0,
            "com np=%d nao ha' celula de franja nenhuma -- o caso abaixo nao "
            "afirmaria nada", ntasks);
    }

    t_case("franja_recebe_o_valor_da_propria_posicao");
    T_CHECK_MSG(erradas_g == 0,
        "%ld de %ld celula(s) de franja com valor de outra posicao (pior erro "
        "%.6g).  %s", erradas_g, n_franja_g, (double) pior_g, primeira);

    return t_end();
}
