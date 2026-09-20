// O T8CODE PRODUZINDO A MALHA QUE O DOMINIO USA.
//
// Ate' aqui o t8code respondia ao contrato: dado um ponto, ele localiza; dado um
// dominio, ele diz quais celulas existem.  Aqui ele PRODUZ -- a arvore que entra
// no `lb_add_input_tree` tem a estrutura decidida por ele, e tudo o que vem depois
// e' o caminho de producao de verdade: `lb_calc_partition`, `psd_create`,
// `psd_synced_mapper`, instantaneo, estencil.
//
// O QUE ISTO NAO E'.  Nao e' o t8code particionando -- quem reparte aqui e' o
// `lbal` com Zoltan, e o particionador do t8code ja' tem o seu proprio caso
// (test-partition-t8code, clausulas C13 e P1..P4).  E nao e' um exemplo rodando
// sobre malha do t8code: os exemplos leem a malha de arquivo AMR e sao comparados
// com referencia gravada, entao trocar a malha mudaria os numeros por construcao.
// O que se afirma e' o CAMINHO: malha do t8code atravessa a producao inteira e
// cumpre o contrato do outro lado.
//
//   contagem_bate_com_a_do_t8code   a soma global das celulas locais e' exatamente
//       o numero de folhas que o t8code reportou ao produzir.  O numero vem DELE,
//       nao de uma formula repetida aqui -- repetir a formula testaria a minha
//       aritmetica, nao a malha.
//
//   particao_cobre_o_dominio_uma_vez   soma global dos volumes = volume da caixa.
//       Pega celula perdida E celula contada duas vezes.
//
//   malha_nao_e_uniforme            a malha produzida TEM refino: se o alvo nao
//       pegasse, tudo abaixo passaria sobre uma malha uniforme e nao diria nada
//       sobre producao de malha adaptada.  Este caso existe para que o teste nao
//       se auto-esvazie em silencio.
//
//   franja_existe_quando_ha_vizinho com np>1 todo rank tem franja; com np=1
//       nenhum.  Reduzido por MIN e por MAX.
//
//   valor_nao_depende_da_particao   campo linear num ponto FIXO da' o mesmo valor
//       em qualquer np.  E' a P4 aplicada a uma malha que o t8code produziu.
//
//   instantaneo_aprova_o_oraculo    o `sd_snapshot_verify` aprova o dominio
//       montado sobre a malha do t8code.  Fecha o circuito: producao do t8code ->
//       caminho de producao da HiGTree -> fronteira das consultas.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "higtree.h"
#include "higtree-iterator.h"
#include "pdomain.h"
#include "lbal.h"
#include "hig-mesh-snapshot.h"
#include "utils.h"
#include "testing.h"
#include "t8code/t8-mesh-producer.h"

#define NIVEL_BASE 3         // 2^3 = 8 celulas por direcao
#define REFINOS    2
#define ALVO       0.3125
#define FRINGE     2

static real campo(const Point p) { return 1.0 + 2.0 * p[0] + 3.0 * p[1]; }

int main(int argc, char *argv[]) {
    higtree_initialize(&argc, &argv);

    int rank, ntasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    partition_graph *pg = pg_create(MPI_COMM_WORLD);
    pg_set_fringe_size(pg, FRINGE);

    // ------------------------------------------------ o t8code produz a malha
    long folhas_t8 = 0;
    int niveis_vistos = 0;
    load_balancer *lb = lb_create(MPI_COMM_WORLD, 1);
    if (rank == 0) {
        hig_cell *raiz = t8_produz_malha_para_dominio(NIVEL_BASE, ALVO, REFINOS,
                                                      &folhas_t8);
        if (raiz != NULL) {
            // Quantos tamanhos de celula distintos a malha tem?  Um so' significa
            // uniforme, e o refino nao pegou.
            real menor = 1.0, maior = 0.0;
            higcit_celliterator *cit;
            for (cit = higcit_create_all_leaves(raiz); !higcit_isfinished(cit);
                 higcit_nextcell(cit)) {
                Point d;
                hig_get_delta(higcit_getcell(cit), d);
                if (d[0] < menor) menor = d[0];
                if (d[0] > maior) maior = d[0];
            }
            higcit_destroy(cit);
            niveis_vistos = (maior > menor * 1.5) ? 2 : 1;
            lb_add_input_tree(lb, raiz, true, 0);
        }
    }
    MPI_Bcast(&folhas_t8, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&niveis_vistos, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // ------------------------------------------------ o caminho de producao
    lb_calc_partition(lb, pg);

    sim_domain *sd = sd_create(NULL);
    sd_set_interpolator_order(sd, 2);
    for (unsigned i = 0; i < lb_get_num_local_trees(lb); ++i) {
        sd_add_higtree(sd, lb_get_local_tree(lb, i, NULL));
    }
    lb_destroy(lb);

    psim_domain *psd = psd_create(sd, pg);
    psd_synced_mapper(psd);

    // ------------------------------------------------ o que este rank enxerga
    long n_local = 0;
    real vol_local = 0.0;
    higcit_celliterator *cit;
    for (cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit);
         higcit_nextcell(cit)) {
        Point d;
        hig_get_delta(higcit_getcell(cit), d);
        real v = 1.0;
        for (int k = 0; k < DIM; k++) v *= d[k];
        vol_local += v;
        n_local++;
    }
    higcit_destroy(cit);

    const long tem_franja = (sd_get_num_fringe_higtrees(sd) > 0) ? 1 : 0;

    char detalhe[256];
    const long oraculo_ok = (sd_snapshot_verify(sd, detalhe, sizeof detalhe) == 0) ? 1 : 0;

    // Ponto fixo, estritamente dentro de uma celula em toda direcao.
    Point alvo;
    POINT_ASSIGN_SCALAR(alvo, 0.6796875);
    real v_local = 0.0;
    long possui = 0;
    {
        // `sd_get_cell_with_point` varre TAMBEM as arvores de franja, e a franja
        // de um vizinho contem o mesmo ponto -- sem o filtro, dois ranks o
        // reivindicam e a reducao por soma soma duas vezes.  O criterio de local
        // e' o mesmo do test-fringe-parallel: `lid < n_local`.
        hig_cell *c = sd_get_cell_with_point(sd, alvo);
        if (c != NULL) {
            mp_mapper *m = sd_get_domain_mapper(sd);
            const int lid = mp_lookup(m, hig_get_cid(c));
            if (lid >= 0 && lid < (int) n_local) {
                Point ce;
                hig_get_center(c, ce);
                possui = 1;
                v_local = campo(ce);
            }
        }
    }

    long n_global, franja_min, franja_max, oraculo_min, donos;
    real vol_global, v_global;
    MPI_Allreduce(&n_local,     &n_global,    1, MPI_LONG,   MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&vol_local,   &vol_global,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&tem_franja,  &franja_min,  1, MPI_LONG,   MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&tem_franja,  &franja_max,  1, MPI_LONG,   MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&oraculo_ok,  &oraculo_min, 1, MPI_LONG,   MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&possui,      &donos,       1, MPI_LONG,   MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&v_local,     &v_global,    1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (rank != 0) return 0;

    t_case("malha_nao_e_uniforme");
    T_CHECK_MSG(niveis_vistos == 2,
        "o t8code produziu uma malha UNIFORME (um tamanho de celula so'): o alvo "
        "de refino nao pegou, e todos os casos abaixo passariam sem dizer nada "
        "sobre producao de malha adaptada");

    t_case("contagem_bate_com_a_do_t8code");
    T_CHECK_MSG(folhas_t8 > 0, "o produtor do t8code devolveu NULL");
    T_CHECK_MSG(n_global == folhas_t8,
        "o t8code produziu %ld folha(s) e o dominio montado tem %ld celula(s) "
        "locais somadas (np=%d)", folhas_t8, n_global, ntasks);

    t_case("particao_cobre_o_dominio_uma_vez");
    T_NEAR(vol_global, 1.0, 1.0e-12, "soma global dos volumes locais");

    t_case("franja_existe_quando_ha_vizinho");
    if (ntasks == 1) {
        T_CHECK_MSG(franja_max == 0, "com np=1 nenhum rank deveria ter franja");
    } else {
        T_CHECK_MSG(franja_min == 1,
            "com np=%d todo rank deveria ter franja, e ao menos um nao tem", ntasks);
    }

    t_case("valor_nao_depende_da_particao");
    T_CHECK_MSG(donos == 1,
        "o ponto alvo deveria pertencer a exatamente um rank, e pertence a %ld",
        donos);
    {
        // O valor esperado vem da MALHA, nao da particao: 0,6796875 cai na celula
        // de indice 5 no nivel base (lado 0,125), cujo centro e' 0,6875 em toda
        // direcao -- longe do alvo de refino, entao ela nao foi subdividida.
        // Afirmar o valor EXATO e' mais forte que afirmar "nao mudou": pega
        // tambem a particao que entrega a celula errada.
        const real c = 0.6875;
        const real esperado = 1.0 + 2.0 * c + 3.0 * c;
        T_NEAR(v_global, esperado, 1.0e-12,
               "campo linear no ponto fixo, independente da particao");
    }

    t_case("instantaneo_aprova_o_oraculo");
    T_CHECK_MSG(oraculo_min == 1,
        "o oraculo diferencial reprovou o instantaneo em ao menos um rank, sobre "
        "malha produzida pelo t8code");

    return t_end();
}
