// Estencil ATRAVESSANDO salto de nivel, inclusive salto MAIOR que 2:1.
//
// POR QUE ESTE TESTE PODE MUDAR A ESTRATEGIA, e nao apenas verifica-la.  A
// contribuicao de Sousa et al. (2019) que a HiGTree implementa e' minimos
// quadrados moveis em arvore NAO GRADUADA -- isto e', sem exigir que celulas
// vizinhas difiram no maximo um nivel.  O balanceamento 2:1 do t8code e' opcional,
// mas as rotinas de vizinhanca e de ghost dele foram construidas em torno do caso
// balanceado.
//
// Se a interpolacao atravessa um salto 4:1 aqui e reproduz um campo linear
// exatamente, entao a capacidade existe e precisa ser preservada por qualquer
// segunda implementacao de malha -- e' requisito, nao detalhe.  Se NAO atravessa,
// a premissa de que t8code e MTree podem ser pares intercambiaveis precisa ser
// revista ANTES da integracao.  Em qualquer dos dois casos a resposta e' mais
// barata agora do que depois.
//
// MALHA (2D; em 3D a construcao e' a mesma com a terceira direcao junto):
//
//     +--------+--------+        raiz [0,1]^DIM refinada 4x4  -> h = 0,25
//     |        |        |
//     |        |  celula escolhida, refinada 2x2  -> h = 0,125
//     +--------+---+----+        e uma NETA refinada de novo  -> h = 0,0625
//     |        | . |    |
//     |        +---+----+        vizinha grossa continua em 0,25:
//     +--------+--------+        salto de QUATRO para um.

#include <stdlib.h>
#include "testing.h"
#include "domain.h"
#include "utils.h"

static real campo(const Point p) {
    real v = 0.5;
    const real coef[3] = { 2.0, -3.0, 1.5 };
    for(int d = 0; d < DIM; d++) v += coef[d] * p[d];
    return v;
}

static void parede_dirichlet(sim_domain *sd, const Point lo, const Point hi,
                             const int nc[DIM]) {
    hig_cell *bc = hig_create_root((real *) lo, (real *) hi);
    hig_refine_uniform(bc, (int *) nc);
    mp_mapper *bm = mp_create();
    higcit_celliterator *it = higcit_create_all_leaves(bc);
    mp_assign_from_celliterator(bm, it, 0);
    higcit_destroy(it);
    sim_boundary *sb = sb_create(bc, DIRICHLET, bm);
    for(it = higcit_create_all_leaves(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        Point ct; hig_get_center(c, ct);
        sb_set_value(sb, mp_lookup(bm, hig_get_id(c, 0)), campo(ct));
    }
    higcit_destroy(it);
    sd_add_boundary(sd, sb);
}

// Refina a celula que contem `p` em 2 por direcao.  Devolve a celula refinada.
static hig_cell *refina_em(hig_cell *raiz, const Point p) {
    hig_cell *c = hig_get_cell_with_point(raiz, (real *) p);
    if(c == NULL) return NULL;
    int nc[DIM];
    for(int d = 0; d < DIM; d++) nc[d] = 2;
    hig_refine_uniform(c, nc);
    return c;
}

int main(void) {
    sim_stencil *stn = stn_create();

    Point l, h;
    for(int d = 0; d < DIM; d++) { l[d] = 0.0; h[d] = 1.0; }
    hig_cell *raiz = hig_create_root(l, h);
    int nc[DIM];
    for(int d = 0; d < DIM; d++) nc[d] = 4;
    hig_refine_uniform(raiz, nc);            // h = 0,25

    // Um passo de refino: 2:1 contra a vizinha grossa.
    Point p2;
    for(int d = 0; d < DIM; d++) p2[d] = 0.375;
    hig_cell *meia = refina_em(raiz, p2);    // h = 0,125

    // Segundo passo, numa neta: 4:1 contra a vizinha grossa.
    Point p4;
    for(int d = 0; d < DIM; d++) p4[d] = 0.3125;
    hig_cell *quarto = refina_em(raiz, p4);  // h = 0,0625

    t_case("malha_nao_graduada_foi_construida");
    T_CHECK_MSG(meia != NULL && quarto != NULL,
        "nao consegui refinar: meia=%p quarto=%p", (void *) meia, (void *) quarto);
    {   // O salto so' e' 4:1 se a vizinha imediata continuou grossa.
        Point viz;  for(int d = 0; d < DIM; d++) viz[d] = 0.625;
        hig_cell *cv = hig_get_cell_with_point(raiz, viz);
        Point dv, dq;
        hig_get_delta(cv, dv);
        hig_get_delta(hig_get_cell_with_point(raiz, p4), dq);
        T_NEAR(dv[0] / dq[0], 4.0, 1e-9, "razao de tamanho entre vizinha grossa e fina");
    }

    mp_mapper *m = mp_create();
    higcit_celliterator *it = higcit_create_all_leaves(raiz);
    const int n = mp_assign_from_celliterator(m, it, 0);
    higcit_destroy(it);

    sim_domain *sd = sd_create(m);
    sd_add_higtree(sd, raiz);
    for(int d = 0; d < DIM; d++) {
        for(int lado = 0; lado < 2; lado++) {
            Point lo, hi; int bnc[DIM];
            for(int k = 0; k < DIM; k++) { lo[k] = 0.0; hi[k] = 1.0; bnc[k] = 4; }
            const real plano = lado ? 1.0 : 0.0;
            lo[d] = plano - EPSDELTA; hi[d] = plano + EPSDELTA; bnc[d] = 1;
            parede_dirichlet(sd, lo, hi, bnc);
        }
    }

    real *fval = (real *) calloc(n, sizeof *fval);
    real *fdel = (real *) calloc(n, sizeof *fdel);   // tamanho de cada celula, por id
    for(it = higcit_create_all_leaves(raiz); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        Point ct, dt; hig_get_center(c, ct); hig_get_delta(c, dt);
        const int id = mp_lookup(m, hig_get_id(c, 0));
        fval[id] = campo(ct);
        fdel[id] = dt[0];
    }
    higcit_destroy(it);

    // O oraculo: minimos quadrados moveis de ordem >= 1 reproduz campo linear
    // EXATAMENTE, e isso nao deveria depender da graduacao da malha.  Se o suporte
    // atravessa o salto de nivel corretamente, o valor bate; se ele se confunde com
    // celulas de tamanhos diferentes, nao bate.
    t_case("reproduz_campo_linear_atravessando_salto_4_para_1");
    const real eps = 1e-10;
    real tamanhos_vistos = 1.0;      // maior razao de tamanhos dentro de um suporte
    // Pontos de um lado e do outro da interface entre a regiao fina e a grossa,
    // e na propria interface.
    const real xs[] = { 0.28125, 0.3125, 0.34375, 0.375, 0.40625, 0.4375, 0.5, 0.5625 };
    for(unsigned i = 0; i < sizeof xs / sizeof *xs; i++) {
        Point org, x;
        for(int d = 0; d < DIM; d++) { org[d] = 0.3125; x[d] = 0.3125; }
        x[0]   = xs[i];
        org[0] = xs[i] - 0.03125;
        stn_reset(stn);
        sd_get_stencil(sd, org, x, 1.0, stn);
        const real v = stn_mult_vector(stn, fval) - stn_get_rhs(stn);
        T_NEAR(v, campo(x), eps, "campo linear atravessando o salto");

        // PRE-CONDICAO por ponto: o suporte tem de conter celulas de MAIS DE UM
        // tamanho, senao ele nao atravessou salto nenhum e o acerto acima nao diz
        // nada sobre malha nao graduada.
        real menor = 1e30, maior = 0.0;
        for(int k = 0; k < stn_get_numelems(stn); k++) {
            const real d = fdel[stn_get_id(stn, k)];
            if(d < menor) menor = d;
            if(d > maior) maior = d;
        }
        if(maior > 0.0) {
            tamanhos_vistos = (maior / menor > tamanhos_vistos)
                            ? (maior / menor) : tamanhos_vistos;
        }
    }

    // Se nenhum suporte misturou tamanhos, o teste passou sem exercitar o salto.
    T_CHECK_MSG(tamanhos_vistos >= 3.9,
        "nenhum suporte atravessou o salto: maior razao de tamanhos dentro de um "
        "estencil foi %.2f, esperado 4 -- o teste nao exercitou malha nao graduada",
        tamanhos_vistos);

    free(fval);
    free(fdel);
    stn_destroy(stn);
    sd_destroy(sd);
    return t_end();
}
