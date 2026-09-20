// PRODUCAO POR RANK, em CAIXAS COMPLETAS.
//
// A floresta nasce em COMM_WORLD -- o t8code ja' a reparte por curva de
// preenchimento -- e cada rank materializa so' o que e' dele.  Isso tira o
// gargalo de um processo ter de caber a malha inteira antes de existir particao.
//
// POR QUE CAIXAS, E NAO UMA ARVORE COM BURACOS.  A tentativa com buracos nao
// sobrevive: `hig_get_cell_coords_of_point` desreferencia os filhos para ler as
// caixas deles, e `hig_get_cell_with_point` ainda faz `cell = cell->children[p]`
// sem testar nulo.  O `lbal` tambem nao produz buraco -- ele preenche todos os
// compartimentos, e por isso um dominio tem VARIAS arvores por rank.
//
//   nenhuma_familia_dividida     o t8code reparte por ELEMENTO, entao ele pode
//       deixar parte dos filhos de uma celula base com o vizinho -- e uma caixa
//       completa nao representa isso.  `set_for_coarsening=1` pede que nao
//       divida; este caso CONFERE, em vez de confiar.  E' o unico caso que pode
//       falhar por motivo do t8code e nao meu, e por isso vem primeiro.
//
//   cobertura_global_bate_com_o_t8code  a soma das folhas materializadas em todos
//       os ranks e' o total que o t8code reporta.  O numero vem DELE.
//
//   volume_global_e_o_da_caixa   soma global dos volumes = volume da caixa.  Pega
//       celula perdida E celula materializada em dois ranks.
//
//   caixas_sao_navegaveis        `hig_get_cell_with_point` no centro de CADA
//       folha devolve aquela folha.  E' o caso que a versao com buracos nao
//       passaria -- ela morria em SEGV aqui --, e e' o que autoriza usar estas
//       arvores como arvore de dominio.
//
//   varias_caixas_quando_reparte com np>1 os ranks nao recebem o dominio inteiro:
//       a soma das caixas locais e' menor que o total.  Sem isto, um produtor que
//       ignorasse a particao passaria em tudo acima.
//
//   franja_existe_quando_ha_vizinho

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "higtree.h"
#include "higtree-iterator.h"
#include "utils.h"
#include "testing.h"
#include "t8code/t8-mesh-rank.h"

#define NIVEL_BASE 3
#define REFINOS    2

static void caixa(Point lo, Point hi) {
    const real LO[3] = { 0.1,  0.2, -0.3 };
    const real HI[3] = { 0.9,  1.4,  0.5 };
    for (int d = 0; d < DIM; d++) { lo[d] = LO[d]; hi[d] = HI[d]; }
}

static void para_dominio(real u, Point p) {
    Point lo, hi;
    caixa(lo, hi);
    for (int d = 0; d < DIM; d++) p[d] = lo[d] + u * (hi[d] - lo[d]);
}

int main(int argc, char *argv[]) {
    higtree_initialize(&argc, &argv);

    int rank, ntasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    Point lo, hi, alvo;
    caixa(lo, hi);
    para_dominio(0.3125, alvo);

    t8_producao_rank p;
    const long ok = t8_produz_por_rank(lo, hi, NIVEL_BASE, alvo, REFINOS, &p) ? 1 : 0;

    long folhas = 0, nav_ruins = 0;
    real vol = 0.0;
    for (int i = 0; i < p.n_locais; i++) {
        higcit_celliterator *it;
        for (it = higcit_create_all_leaves(p.locais[i]); !higcit_isfinished(it);
             higcit_nextcell(it)) {
            hig_cell *c = higcit_getcell(it);
            Point d, ce;
            hig_get_delta(c, d);
            hig_get_center(c, ce);
            real v = 1.0;
            for (int k = 0; k < DIM; k++) v *= d[k];
            vol += v; folhas++;
            // Navegabilidade: o centro da folha tem de voltar na propria folha.
            if (hig_get_cell_with_point(p.locais[i], ce) != c) nav_ruins++;
        }
        higcit_destroy(it);
    }

    const long tem_franja = (p.n_franja > 0) ? 1 : 0;
    const long dividida = p.base_dividida;
    const long n_global_t8 = p.n_global;

    long folhas_soma, div_soma, nav_max, franja_min, franja_max, ok_min;
    real vol_soma;
    MPI_Allreduce(&folhas,    &folhas_soma, 1, MPI_LONG,   MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&vol,       &vol_soma,    1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&dividida,  &div_soma,    1, MPI_LONG,   MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nav_ruins, &nav_max,     1, MPI_LONG,   MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&tem_franja,&franja_min,  1, MPI_LONG,   MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&tem_franja,&franja_max,  1, MPI_LONG,   MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&ok,        &ok_min,      1, MPI_LONG,   MPI_MIN, MPI_COMM_WORLD);

    const long folhas_locais = folhas;
    t8_producao_rank_destroi(&p);
    if (rank != 0) return 0;

    t_case("producao_por_rank_nao_falha");
    T_CHECK_MSG(ok_min == 1, "t8_produz_por_rank falhou em ao menos um rank");

    t_case("nenhuma_familia_dividida");
    T_CHECK_MSG(div_soma == 0,
        "%ld celula(s) base com familia dividida entre ranks (np=%d).  A caixa "
        "completa nao representa isso; a saida e' pedir ao t8code que nao divida "
        "familia, nao remendar a materializacao", div_soma, ntasks);

    t_case("cobertura_global_bate_com_o_t8code");
    T_CHECK_MSG(n_global_t8 > 0, "o t8code reportou floresta vazia");
    T_CHECK_MSG(folhas_soma == n_global_t8,
        "os ranks materializaram %ld folha(s) somadas e o t8code diz que a "
        "floresta tem %ld (np=%d)", folhas_soma, n_global_t8, ntasks);

    t_case("volume_global_e_o_da_caixa");
    {
        Point l, h;
        caixa(l, h);
        real v = 1.0;
        for (int d = 0; d < DIM; d++) v *= (h[d] - l[d]);
        T_NEAR(vol_soma, v, 1.0e-12, "soma global dos volumes materializados");
    }

    t_case("caixas_sao_navegaveis");
    T_CHECK_MSG(nav_max == 0,
        "ate' %ld folha(s) num mesmo rank nao sao encontradas por "
        "hig_get_cell_with_point a partir do proprio centro -- as caixas nao "
        "servem como arvore de dominio", nav_max);

    t_case("varias_caixas_quando_reparte");
    if (ntasks == 1) {
        T_CHECK_MSG(folhas_locais == n_global_t8,
            "com np=1 o rank unico deveria ter a malha inteira: %ld de %ld",
            folhas_locais, n_global_t8);
    } else {
        T_CHECK_MSG(folhas_locais > 0 && folhas_locais < n_global_t8,
            "com np=%d o rank 0 deveria ter uma PARTE da malha, e tem %ld de %ld",
            ntasks, folhas_locais, n_global_t8);
    }

    t_case("franja_existe_quando_ha_vizinho");
    if (ntasks == 1) {
        T_CHECK_MSG(franja_max == 0, "com np=1 nao deveria haver ghost");
    } else {
        T_CHECK_MSG(franja_min == 1, "com np=%d todo rank deveria ter franja", ntasks);
    }

    return t_end();
}
