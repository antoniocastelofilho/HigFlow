// O ORACULO DIFERENCIAL, e a prova de que ele acusa.
//
// O `sd_snapshot_verify` existe para AUTORIZAR a migracao dos lacos: um laco que
// deixa de navegar a arvore e passa a ler o arranjo esta' certo se o oraculo
// aprova, sem depender do resultado fisico para saber disso.  Uma autorizacao que
// aprova tudo nao autoriza nada, entao metade deste arquivo e' a validacao inversa.
//
// POR QUE ELE LOCALIZA POR PONTO EM VEZ DE PERCORRER.  O instantaneo foi
// preenchido percorrendo o iterador de dominio e indexando pelo mapeador; e o
// mapeador foi atribuido a partir DESSE MESMO iterador, no `_psd_setmapper`.
// Conferir um contra o outro e' medir a aritmetica contra ela propria.  Isso nao e'
// conjectura: MEDIDO, trocar `mp_lookup` por um contador de percurso dentro do
// `hms_from_domain` deixa a suite INTEIRA verde, nos dois DIM e nos tres np.  A
// busca por ponto e' o unico caminho ate' a celula que nao passa por onde o
// instantaneo passou.
//
//   aprova_instantaneo_correto     o recem-produzido passa com zero divergencias.
//
//   acusa_linha_trocada            trocando DUAS linhas do arranjo entre si, o
//       oraculo acusa.  E' a falha que o `hms_from_domain` nao poderia cometer
//       hoje, mas que um segundo backend preenchendo o arranjo direto da propria
//       estrutura pode -- e e' exatamente para esse caso que o oraculo existe.
//
//   acusa_caixa_deslocada          mexendo na caixa de UMA linha, o oraculo
//       acusa.  Separado do anterior de proposito: troca de linha e' erro de
//       INDICE, perturbacao e' erro de VALOR, e um oraculo pode pegar um e nao o
//       outro.
//
//   acusa_caixa_alargada           idem para o tamanho.  O centro certo com o
//       delta errado passaria por qualquer conferencia que so' olhasse posicao --
//       e o delta e' o que 144 sitios de higflow/src leem.
//
// O teste escreve no instantaneo por um ponteiro NAO-const obtido com cast.  A
// API o devolve como const de proposito -- ninguem em producao deve escrever nele
// --, e corromper de proposito e' justamente o que um teste de oraculo faz.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "higtree.h"
#include "higtree-iterator.h"
#include "domain.h"
#include "hig-mesh-snapshot.h"
#include "utils.h"
#include "testing.h"

#define NC 4

static sim_domain *monta(void) {
    Point lo, hi;
    POINT_ASSIGN_SCALAR(lo, 0.0);
    POINT_ASSIGN_SCALAR(hi, 1.0);
    hig_cell *raiz = hig_create_root(lo, hi);
    int nc[DIM];
    for (int d = 0; d < DIM; d++) nc[d] = NC;
    hig_refine_uniform(raiz, nc);

    // Um nivel a mais em uma celula: malha nao uniforme, para que a busca por
    // ponto tenha de descer a arvore e nao acerte por simetria.
    Point p;
    POINT_ASSIGN_SCALAR(p, 0.375);
    hig_cell *c = hig_get_cell_with_point(raiz, p);
    if (c != NULL) {
        int n2[DIM];
        for (int d = 0; d < DIM; d++) n2[d] = 2;
        hig_refine_uniform(c, n2);
    }

    sim_domain *sd = sd_create(NULL);
    sd_add_higtree(sd, raiz);
    mp_mapper *m = sd_get_domain_mapper(sd);
    higcit_celliterator *it = sd_get_domain_celliterator(sd);
    mp_assign_from_celliterator(m, it, 0);
    higcit_destroy(it);
    sd_compute_snapshot(sd);
    return sd;
}

// O instantaneo, escrevivel.  Ver o cabecalho.
static hig_mesh_snapshot *corrompivel(sim_domain *sd) {
    return (hig_mesh_snapshot *) sd_get_snapshot(sd);
}

int main(int argc, char *argv[]) {
    higtree_initialize(&argc, &argv);
    char detalhe[256];

    // ------------------------------------------------------------------
    t_case("aprova_instantaneo_correto");
    {
        sim_domain *sd = monta();
        const int ruins = sd_snapshot_verify(sd, detalhe, sizeof detalhe);
        T_CHECK_MSG(ruins == 0,
            "o oraculo reprovou um instantaneo recem-produzido: %d linha(s).  %s",
            ruins, detalhe);
        sd_destroy(sd);
    }

    // ------------------------------------------------------------------
    t_case("acusa_linha_trocada");
    {
        sim_domain *sd = monta();
        hig_mesh_snapshot *s = corrompivel(sd);
        T_CHECK_MSG(s != NULL && s->n >= 2, "malha pequena demais para trocar");
        if (s != NULL && s->n >= 2) {
            // Duas linhas de celulas DIFERENTES.  Trocar linhas identicas nao
            // seria corrupcao nenhuma.
            int a = 0, b = s->n - 1;
            for (int d = 0; d < DIM; d++) {
                real t = s->low[a * DIM + d];
                s->low[a * DIM + d] = s->low[b * DIM + d];
                s->low[b * DIM + d] = t;
                t = s->high[a * DIM + d];
                s->high[a * DIM + d] = s->high[b * DIM + d];
                s->high[b * DIM + d] = t;
            }
            const int ruins = sd_snapshot_verify(sd, detalhe, sizeof detalhe);
            T_CHECK_MSG(ruins >= 2,
                "troquei as linhas %d e %d e o oraculo acusou %d divergencia(s) "
                "-- deveria acusar as duas", a, b, ruins);
        }
        sd_destroy(sd);
    }

    // ------------------------------------------------------------------
    t_case("acusa_caixa_deslocada");
    {
        sim_domain *sd = monta();
        hig_mesh_snapshot *s = corrompivel(sd);
        if (s != NULL && s->n > 0) {
            // Deslocamento pequeno, dentro da MESMA celula: a busca por ponto
            // ainda devolve a celula certa, com o id certo.  Quem tem de pegar e'
            // a comparacao de geometria, e nao a de indice.
            const int i = s->n / 2;
            // Desloca a CAIXA inteira: o centro anda 1e-6/2, pouco para sair da
            // celula, entao a busca por ponto ainda devolve a celula certa e
            // quem tem de pegar e' a comparacao de geometria.
            s->low[i * DIM]  += 1.0e-6;
            s->high[i * DIM] += 1.0e-6;
            const int ruins = sd_snapshot_verify(sd, detalhe, sizeof detalhe);
            T_CHECK_MSG(ruins == 1,
                "desloquei o centro da linha %d em 1e-6 e o oraculo acusou %d "
                "divergencia(s), esperava 1.  %s", i, ruins, detalhe);
        }
        sd_destroy(sd);
    }

    // ------------------------------------------------------------------
    t_case("acusa_caixa_alargada");
    {
        sim_domain *sd = monta();
        hig_mesh_snapshot *s = corrompivel(sd);
        if (s != NULL && s->n > 0) {
            // Centro intacto, tamanho errado.  Passaria por qualquer conferencia
            // que so' olhasse posicao -- e o delta e' o que 144 sitios leem.
            const int i = s->n / 3;
            // Canto inferior intacto, superior dobrado: o tamanho muda e a
            // posicao do canto de baixo nao.  Passaria por qualquer conferencia
            // que so' olhasse o centro de uma direcao.
            s->high[i * DIM] += (s->high[i * DIM] - s->low[i * DIM]);
            const int ruins = sd_snapshot_verify(sd, detalhe, sizeof detalhe);
            T_CHECK_MSG(ruins == 1,
                "dobrei o delta da linha %d e o oraculo acusou %d divergencia(s), "
                "esperava 1.  %s", i, ruins, detalhe);
        }
        sd_destroy(sd);
    }

    return t_end();
}
