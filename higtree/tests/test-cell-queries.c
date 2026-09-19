// As quatro consultas de celula que a aplicacao mais usa, afirmadas por VALOR.
//
// POR QUE ESTE TESTE EXISTE.  `hig_get_center`, `hig_get_delta`, `hig_get_cid` e o
// iterador de celulas locais somam mais de 850 usos em higflow/src e sao quatro das
// seis consultas que o contrato de Mesh vai congelar.  Ate' aqui tinham apenas
// cobertura INDIRETA: passavam porque os outros testes dependem delas.  Cobertura
// indireta mede que a funcao nao explodiu, nao que ela devolveu o numero certo --
// e congelar interface sobre comportamento nao verificado e' exatamente o erro que
// esta suite existe para impedir.
//
// Todos os oraculos aqui sao ANALITICOS: saem da malha construida, nao de referencia
// gravada.  A malha tem DOIS niveis de propósito, porque acessor que acerta num
// tamanho de celula pode errar noutro.
//
//     raiz [0,1]^DIM refinada 4 por direcao      -> h = 0,25
//     uma celula refinada 2 por direcao          -> h = 0,125
//
// O que cada caso afirma:
//   centro_e_delta_coerentes_com_a_caixa   centro = (low+high)/2 e delta = high-low,
//                                          o que amarra os tres acessores entre si
//   centro_na_grade_analitica              cada centro cai em (k+0,5)*h exato
//   delta_assume_os_dois_tamanhos          os dois niveis aparecem, e nada alem
//   cid_e_bijecao                          ids distintos, e mp_lookup cobre [0,n)
//   iterador_cobre_o_dominio_uma_vez       soma dos volumes = volume do dominio

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "higtree.h"
#include "higtree-iterator.h"
#include "domain.h"
#include "testing.h"

#define NC   4
#define TOL  1.0e-14

static bool quase(real a, real b) { return fabs(a - b) <= TOL; }

// Refina, em 2 por direcao, a celula que contem `p`.
static void refina_em(hig_cell *raiz, const Point p) {
    hig_cell *c = hig_get_cell_with_point(raiz, p);
    T_CHECK_MSG(c != NULL, "nenhuma celula contem o ponto de refino");
    int nc[DIM];
    for (int d = 0; d < DIM; d++) nc[d] = 2;
    hig_refine_uniform(c, nc);
}

int main(void) {
    Point lo, hi;
    POINT_ASSIGN_SCALAR(lo, 0.0);
    POINT_ASSIGN_SCALAR(hi, 1.0);
    hig_cell *raiz = hig_create_root(lo, hi);
    int nc[DIM];
    for (int d = 0; d < DIM; d++) nc[d] = NC;
    hig_refine_uniform(raiz, nc);

    const real h_grosso = 1.0 / NC;          // 0,25
    const real h_fino   = h_grosso / 2.0;    // 0,125
    Point alvo;
    POINT_ASSIGN_SCALAR(alvo, 0.375);        // uma celula do meio
    refina_em(raiz, alvo);

    sim_domain *sd = sd_create(NULL);
    sd_add_higtree(sd, raiz);
    mp_mapper *m = sd_get_domain_mapper(sd);
    higcit_celliterator *cit = sd_get_domain_celliterator(sd);
    int n = mp_assign_from_celliterator(m, cit, 0);
    higcit_destroy(cit);

    const long folhas = hig_get_number_of_leaves(raiz);
    T_CHECK_MSG(n == folhas,
        "o mapeador registrou %d celulas e a arvore tem %ld folhas", n, folhas);

    // ---------------------------------------------------------------- acessores
    t_case("centro_e_delta_coerentes_com_a_caixa");
    // Reporta o PRIMEIRO contraexemplo e para.  O invariante e' uniforme: uma
    // linha por celula enterraria o sinal em centenas de repeticoes iguais.
    bool achou = false;
    for (cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit) && !achou;
         higcit_nextcell(cit)) {
        hig_cell *c = higcit_getcell(cit);
        Point ce, de, p0, p1;
        hig_get_center(c, ce);
        hig_get_delta(c, de);
        hig_get_lowpoint(c, p0);
        hig_get_highpoint(c, p1);
        for (int d = 0; d < DIM && !achou; d++) {
            if (!quase(ce[d], 0.5 * (p0[d] + p1[d]))) {
                _t_fail("centro[%d] = %.17g mas (low+high)/2 = %.17g",
                        d, ce[d], 0.5 * (p0[d] + p1[d]));
                achou = true;
            } else if (!quase(de[d], p1[d] - p0[d])) {
                _t_fail("delta[%d] = %.17g mas high-low = %.17g",
                        d, de[d], p1[d] - p0[d]);
                achou = true;
            }
        }
    }
    higcit_destroy(cit);

    t_case("centro_na_grade_analitica");
    achou = false;
    for (cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit) && !achou;
         higcit_nextcell(cit)) {
        hig_cell *c = higcit_getcell(cit);
        Point ce, de;
        hig_get_center(c, ce);
        hig_get_delta(c, de);
        for (int d = 0; d < DIM && !achou; d++) {
            // numa malha uniforme de passo de[d], o centro e' (k + 1/2)*de[d]
            real k = ce[d] / de[d] - 0.5;
            if (fabs(k - round(k)) > 1.0e-12) {
                _t_fail("centro[%d] = %.17g nao cai na grade de passo %.17g "
                        "(indice %.17g, nao inteiro)", d, ce[d], de[d], k);
                achou = true;
            }
        }
    }
    higcit_destroy(cit);

    t_case("delta_assume_os_dois_tamanhos");
    int n_grosso = 0, n_fino = 0, n_outro = 0;
    for (cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit);
         higcit_nextcell(cit)) {
        Point de;
        hig_get_delta(higcit_getcell(cit), de);
        bool g = true, f = true;
        for (int d = 0; d < DIM; d++) {
            if (!quase(de[d], h_grosso)) g = false;
            if (!quase(de[d], h_fino))   f = false;
        }
        if (g) n_grosso++; else if (f) n_fino++; else n_outro++;
    }
    higcit_destroy(cit);
    // a celula refinada saiu e entraram 2^DIM filhas
    int esperado_fino = 1;
    for (int d = 0; d < DIM; d++) esperado_fino *= 2;
    int esperado_grosso = 1;
    for (int d = 0; d < DIM; d++) esperado_grosso *= NC;
    esperado_grosso -= 1;
    T_CHECK_MSG(n_outro == 0, "%d celulas com delta que nao e' nenhum dos dois "
        "niveis (%.17g ou %.17g)", n_outro, h_grosso, h_fino);
    T_CHECK_MSG(n_grosso == esperado_grosso && n_fino == esperado_fino,
        "celulas por nivel: grossas %d (esperado %d), finas %d (esperado %d)",
        n_grosso, esperado_grosso, n_fino, esperado_fino);

    // --------------------------------------------------------------------- cid
    t_case("cid_e_bijecao");
    {
        int *visto = (int *) calloc(n, sizeof *visto);
        int fora = 0, repetido = 0, total = 0;
        for (cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit);
             higcit_nextcell(cit)) {
            uniqueid id = hig_get_cid(higcit_getcell(cit));
            int lid = mp_lookup(m, id);
            total++;
            if (lid < 0 || lid >= n) { fora++; continue; }
            if (visto[lid]++) repetido++;
        }
        higcit_destroy(cit);
        int nao_cobertos = 0;
        for (int i = 0; i < n; i++) if (visto[i] == 0) nao_cobertos++;
        free(visto);
        T_CHECK_MSG(total == n && fora == 0 && repetido == 0 && nao_cobertos == 0,
            "mp_lookup(hig_get_cid) deveria ser bijecao sobre [0,%d): "
            "visitadas %d, fora da faixa %d, repetidas %d, nao cobertas %d",
            n, total, fora, repetido, nao_cobertos);
    }

    t_case("cid_estavel_entre_percursos");
    {
        // o mesmo ponto geometrico tem de dar o mesmo id em dois percursos
        uniqueid *ids = (uniqueid *) malloc(n * sizeof *ids);
        Point *ces = (Point *) malloc(n * sizeof *ces);
        int k = 0;
        for (cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit);
             higcit_nextcell(cit)) {
            hig_cell *c = higcit_getcell(cit);
            ids[k] = hig_get_cid(c);
            hig_get_center(c, ces[k]);
            k++;
        }
        higcit_destroy(cit);
        int divergentes = 0;
        for (int i = 0; i < k; i++) {
            hig_cell *c = hig_get_cell_with_point(raiz, ces[i]);
            if (c == NULL || hig_get_cid(c) != ids[i]) divergentes++;
        }
        free(ids); free(ces);
        T_CHECK_MSG(divergentes == 0,
            "%d celulas mudaram de cid entre o percurso e a busca por ponto",
            divergentes);
    }

    // --------------------------------------------------------------- iterador
    t_case("iterador_cobre_o_dominio_uma_vez");
    {
        // soma dos volumes = volume do dominio.  Pega celula faltando E celula
        // repetida, que contagem sozinha nao distingue.
        real vol = 0.0;
        int visitadas = 0;
        for (cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit);
             higcit_nextcell(cit)) {
            Point de;
            hig_get_delta(higcit_getcell(cit), de);
            real v = 1.0;
            for (int d = 0; d < DIM; d++) v *= de[d];
            vol += v;
            visitadas++;
        }
        higcit_destroy(cit);
        T_CHECK_MSG(visitadas == n,
            "o iterador visitou %d celulas, o mapeador registrou %d", visitadas, n);
        T_NEAR(vol, 1.0, 1.0e-12, "soma dos volumes das celulas (dominio [0,1]^DIM)");
    }

    return t_end();
}
