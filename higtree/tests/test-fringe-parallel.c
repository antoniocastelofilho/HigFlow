// A franja sob particionamento real: np = 1, 2 e 3.
//
// POR QUE ESTE TESTE EXISTE, e por que o serial nao bastava.  O
// test-fringe-support monta a franja A MAO, em serie, e afirma o contrato dela.
// Isso cobre o que o estencil pede, mas nao cobre quem PRODUZ a franja: o
// build-fringe.cpp e o balanceador, que so' rodam sob MPI.  E e' exatamente ai' que
// o t8code vai colidir, porque ele traz a propria camada de ghost.
//
// O QUE E' AFIRMADO, tudo por REDUCAO GLOBAL e nao por rank.  Um teste paralelo em
// que cada rank imprime o seu veredicto esconde a falha de um rank no meio das
// linhas dos outros; aqui os DADOS sao reduzidos e so' o rank 0 conclui.
//
//   particao_cobre_o_dominio_uma_vez   soma global dos volumes locais = volume do
//                                      dominio, e soma global das celulas locais =
//                                      total.  Pega celula perdida E celula contada
//                                      duas vezes, que contagem por rank nao ve.
//
//   franja_existe_quando_ha_vizinho    com np>1 TODO rank tem franja; com np=1
//                                      NENHUM tem.  Reduzido por MIN e por MAX, o
//                                      que e' mais forte que "alguem tem".
//
//   estencil_alcanca_a_franja          com np>1, em todo rank existe celula local
//                                      cujo suporte inclui celula de franja.  E'
//                                      o mesmo criterio do teste serial -- o que a
//                                      franja entrega, nao o tamanho dela.
//
//   valor_nao_depende_da_particao      um campo linear interpolado num ponto FIXO
//                                      da' o mesmo valor em qualquer np.  E' o
//                                      oraculo de INVARIANCIA, o unico capaz de
//                                      distinguir "particiona diferente" de
//                                      "particiona errado" -- e o t8code vai
//                                      particionar diferente de proposito.
//
// Ids: `psd_synced_mapper` da' aos locais [0, n_local) e a' franja os seguintes,
// entao `lid >= n_local` identifica celula de franja sem precisar de API nova.
//
// VALIDADO NO SENTIDO INVERSO com `pg_set_fringe_size(pg, 0)`:
//
//     np=1   4 de 4 passam        (nao se espera franja)
//     np=2   franja_existe e estencil_alcanca_a_franja FALHAM
//     np=3   idem
//     particao_cobre_o_dominio_uma_vez passa sempre -- a particao continua sendo
//     uma particao, e e' bom que esse caso NAO se mova por causa da franja
//
// O QUE `valor_nao_depende_da_particao` NAO PEGA, e fica escrito para ninguem
// confiar demais nele: com a franja zerada ele continua VERDE.  Campo linear e'
// reproduzido exato por minimos quadrados mesmo com suporte so' de um lado, a mesma
// limitacao registrada no test-fringe-support.  Ele guarda contra franja com
// GEOMETRIA errada -- peso aplicado a celula na posicao errada muda a soma -- e nao
// contra franja AUSENTE.  Quem guarda contra ausencia e' o caso anterior.

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

#define NC       8        // celulas por direcao na raiz
#define FRINGE   2

static real campo(const Point p) { return 1.0 + 2.0 * p[0] + 3.0 * p[1]; }

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
    for (unsigned i = 0; i < lb_get_num_local_trees(lb); ++i) {
        sd_add_higtree(sd, lb_get_local_tree(lb, i, NULL));
    }
    lb_destroy(lb);

    psim_domain *psd = psd_create(sd, pg);
    psd_synced_mapper(psd);
    mp_mapper *m = sd_get_domain_mapper(sd);

    // ------------------------------------------------ o que este rank enxerga
    long   n_local = 0;
    real   vol_local = 0.0;
    higcit_celliterator *cit;
    for (cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit);
         higcit_nextcell(cit)) {
        Point de;
        hig_get_delta(higcit_getcell(cit), de);
        real v = 1.0;
        for (int d = 0; d < DIM; d++) v *= de[d];
        vol_local += v;
        n_local++;
    }
    higcit_destroy(cit);

    long tem_franja = (sd_get_num_fringe_higtrees(sd) > 0) ? 1 : 0;

    // Alguma celula local cujo suporte inclua celula de franja?
    long atravessa = 0;
    for (cit = sd_get_domain_celliterator(sd);
         !higcit_isfinished(cit) && !atravessa; higcit_nextcell(cit)) {
        // A consulta tem de cair FORA do centro.  No centro exato o interpolador
        // devolve a propria celula com peso 1 (ne=1) e nunca alcanca vizinho
        // nenhum -- medido: ne=1, maxid=0 em todo rank.  Usa-se a face alta em x,
        // que obriga o suporte a pegar os dois lados.
        hig_cell *c = higcit_getcell(cit);
        Point ce, de;
        hig_get_center(c, ce);
        hig_get_delta(c, de);
        ce[0] += 0.5 * de[0];
        sim_stencil *stn = stn_create();
        sd_get_stencil(sd, ce, ce, 1.0, stn);
        int ne = stn_get_numelems(stn);
        int *ids = stn_get_ids(stn);
        for (int i = 0; i < ne; i++) {
            if (ids[i] >= n_local) { atravessa = 1; break; }
        }
        stn_destroy(stn);
    }
    higcit_destroy(cit);

    // Campo linear interpolado num ponto FIXO do dominio.  Só o rank que o possui
    // contribui; os demais mandam 0 e a reducao por soma junta.
    // O alvo tem de cair ESTRITAMENTE dentro de uma celula, em todas as direcoes.
    // Com 0,5 numa coordenada ele cai sobre uma FACE da malha de NC=8: dois ranks
    // reivindicam a posse e a reducao por soma soma duas contribuicoes -- medido em
    // np=3, 2 donos e valor dobrado.  E nao pode ser o CENTRO tampouco: la' o
    // interpolador devolve a propria celula com peso 1 e a invariancia fica trivial.
    // 0,59375 = centro (0,5625) mais um quarto de celula, dentro de [0,5; 0,625).
    Point alvo;
    POINT_ASSIGN_SCALAR(alvo, 0.59375);
    real interp_local = 0.0;
    long possui = 0;
    {
        sim_stencil *stn = stn_create();
        sd_get_stencil(sd, alvo, alvo, 1.0, stn);
        int ne = stn_get_numelems(stn);
        int *ids = stn_get_ids(stn);
        real *w = stn_get_vals(stn);
        if (ne > 0) {
            // casa id -> centro percorrendo locais e franja
            for (int t = 0; t < sd_get_num_higtrees(sd); t++) {
                higcit_celliterator *ct =
                    higcit_create_all_leaves(sd_get_higtree(sd, t));
                for (; !higcit_isfinished(ct); higcit_nextcell(ct)) {
                    hig_cell *c = higcit_getcell(ct);
                    int lid = mp_lookup(m, hig_get_cid(c));
                    Point ce;
                    hig_get_center(c, ce);
                    for (int i = 0; i < ne; i++) {
                        if (ids[i] == lid) interp_local += w[i] * campo(ce);
                    }
                }
                higcit_destroy(ct);
            }
            hig_cell *dono = NULL;
            for (unsigned i = 0; i < sd_get_num_local_higtrees(sd) && !dono; i++) {
                dono = hig_get_cell_with_point(sd_get_local_higtree(sd, i), alvo);
            }
            possui = dono ? 1 : 0;
        }
        if (!possui) interp_local = 0.0;
        stn_destroy(stn);
    }

    // ------------------------------------------------------- reducoes globais
    long  n_global, franja_min, franja_max, atravessa_min, donos;
    real  vol_global, interp_global;
    MPI_Allreduce(&n_local,     &n_global,      1, MPI_LONG,   MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&vol_local,   &vol_global,    1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&tem_franja,  &franja_min,    1, MPI_LONG,   MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&tem_franja,  &franja_max,    1, MPI_LONG,   MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&atravessa,   &atravessa_min, 1, MPI_LONG,   MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&possui,      &donos,         1, MPI_LONG,   MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&interp_local,&interp_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    long esperado = 1;
    for (int d = 0; d < DIM; d++) esperado *= NC;

    // So' o rank 0 conclui (ver o cabecalho).  E ninguem chama MPI_Finalize: o
    // `higtree_initialize` registra o PetscFinalize no atexit, e ele ja' finaliza o
    // MPI -- chamar de novo aborta com "MPI_Finalize called after MPI_FINALIZE".
    if (rank != 0) return 0;

    t_case("particao_cobre_o_dominio_uma_vez");
    T_CHECK_MSG(n_global == esperado,
        "soma global das celulas locais = %ld, esperado %ld (np=%d): celula "
        "perdida ou contada duas vezes", n_global, esperado, ntasks);
    T_NEAR(vol_global, 1.0, 1.0e-12, "soma global dos volumes locais");

    t_case("franja_existe_quando_ha_vizinho");
    if (ntasks == 1) {
        T_CHECK_MSG(franja_max == 0,
            "com np=1 nenhum rank deveria ter franja, e algum tem");
    } else {
        T_CHECK_MSG(franja_min == 1,
            "com np=%d TODO rank deveria ter franja, e ao menos um nao tem",
            ntasks);
    }

    t_case("estencil_alcanca_a_franja");
    if (ntasks == 1) {
        T_CHECK_MSG(atravessa_min == 0,
            "com np=1 nenhum suporte deveria sair do dominio local");
    } else {
        T_CHECK_MSG(atravessa_min == 1,
            "com np=%d todo rank deveria ter alguma celula cujo suporte alcanca "
            "a franja, e ao menos um nao tem", ntasks);
    }

    t_case("valor_nao_depende_da_particao");
    // O ponto pertence a exatamente um rank, em qualquer particao.
    T_CHECK_MSG(donos == 1,
        "o ponto alvo deveria pertencer a exatamente um rank, e pertence a %ld",
        donos);
    // Reproducao polinomial: com o suporte completo -- local mais franja -- o
    // linear sai exato, e o valor e' o mesmo em np=1, 2 ou 3.  E' a invariancia
    // que importa para trocar o particionador.
    T_NEAR(interp_global, campo(alvo), 1.0e-12,
        "campo linear no ponto fixo, independente da particao");

    return t_end();
}
