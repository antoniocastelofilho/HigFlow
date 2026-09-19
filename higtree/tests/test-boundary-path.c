// O ramo ON_BOUNDARY do despacho, que os outros testes nao alcancam.
//
// LACUNA QUE ESTE ARQUIVO FECHA.  `get_stencil` (domain.c, por volta de 2161)
// despacha em tres ramos: ponto no centro de uma celula, ponto SOBRE o contorno
// (ON_BOUNDARY), e ponto FORA do dominio (OUTSIDE_DOMAIN).  Os cinco testes
// escritos antes deste entram todos por OUTSIDE_DOMAIN -- medido com contador:
// zero chamadas aos fechamentos `*_boundary` em 20 casos.
//
// E' metade do despacho sem cobertura, e nao a metade menos usada: uma sonda na
// suite do HiGFlow contou 312.480 chamadas a get_stencil_neumann_boundary_any_order
// e 468.222 a get_stencil_dirichlet_boundary numa unica passada de 33 casos.  O
// escritor VTK e' o grande consumidor, porque interpola nos CANTOS das celulas e
// esses caem sobre o contorno.
//
// Foi nesse ramo que sobreviveram tres dos sete sitios do defeito de compactacao
// corrigido em 5145c6e -- eles ficaram verdes por ausencia de teste, nao por
// estarem certos.
//
// SEGUNDO CASO, o invariante de faixa -- e uma ressalva medida sobre ele.
//
// Interpolacao nao produz valor muito fora da faixa dos dados que usa.  Na suite do
// HiGFlow esse criterio e' forte: o defeito de compactacao fazia o VTK reportar
// u.min = -159,283 para um campo que ia de -36,800 a +1,778, e vel.u.max = 205,517
// onde o maximo era 1,778 -- 115x fora.  Nao consulta referencia nenhuma.
//
// AQUI ELE NAO PEGA ESSE DEFEITO, e isso foi medido, nao suposto: revertendo os tres
// sitios, o primeiro caso falha em 14 consultas e este PASSA.  A razao e' que o valor
// errado (constante ao longo da parede) cai DENTRO da faixa, e a faixa inclui a
// contribuicao da propria BC -- que e' justamente o que esta' errado.  Excluir a BC
// da faixa nao resolve: num contorno de Dirichlet o valor legitimamente sai da faixa
// interior.
//
// Fica como detector de erro GROSSEIRO de pesos, que e' barato e cobre outras
// familias de defeito.  A forma forte do invariante precisa comparar a saida com a
// SOLUCAO, e a solucao nao esta' no VTK -- exige instrumentar o solver.

#include <stdlib.h>
#include "testing.h"
#include "domain.h"
#include "utils.h"

#define LADO 8

static real campo(const Point p) {
    // Varia em TODAS as direcoes, senao o defeito de compactacao fica invisivel:
    // com campo constante ao longo da parede, compactar errado nao muda o valor.
    real v = 0.5;
    const real coef[3] = { 2.0, -3.0, 1.5 };
    for(int d = 0; d < DIM; d++) v += coef[d] * p[d];
    return v;
}

int main(void) {
    const real h = 1.0 / LADO;
    sim_stencil *stn = stn_create();

    Point l, hi; int nc[DIM];
    for(int d = 0; d < DIM; d++) { l[d] = 0.0; hi[d] = 1.0; nc[d] = LADO; }
    hig_cell *raiz = hig_create_root(l, hi);
    hig_refine_uniform(raiz, nc);

    mp_mapper *m = mp_create();
    higcit_celliterator *it = higcit_create_all_leaves(raiz);
    const int n = mp_assign_from_celliterator(m, it, 0);
    higcit_destroy(it);

    sim_domain *sd = sd_create(m);
    sd_add_higtree(sd, raiz);

    for(int d = 0; d < DIM; d++) {
        for(int lado = 0; lado < 2; lado++) {
            Point bl, bh; int bnc[DIM];
            for(int k = 0; k < DIM; k++) { bl[k] = 0.0; bh[k] = 1.0; bnc[k] = LADO; }
            const real plano = lado ? 1.0 : 0.0;
            bl[d] = plano - EPSDELTA; bh[d] = plano + EPSDELTA; bnc[d] = 1;
            hig_cell *bc = hig_create_root(bl, bh);
            hig_refine_uniform(bc, bnc);
            mp_mapper *bm = mp_create();
            higcit_celliterator *bit = higcit_create_all_leaves(bc);
            mp_assign_from_celliterator(bm, bit, 0);
            higcit_destroy(bit);
            sim_boundary *sb = sb_create(bc, DIRICHLET, bm);
            for(bit = higcit_create_all_leaves(bc); !higcit_isfinished(bit); higcit_nextcell(bit)) {
                hig_cell *c = higcit_getcell(bit);
                Point ct; hig_get_center(c, ct);
                sb_set_value(sb, mp_lookup(bm, hig_get_id(c, 0)), campo(ct));
            }
            higcit_destroy(bit);
            sd_add_boundary(sd, sb);
        }
    }

    real *fval = (real *) calloc(n, sizeof *fval);
    for(it = higcit_create_all_leaves(raiz); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        Point ct; hig_get_center(c, ct);
        fval[mp_lookup(m, hig_get_id(c, 0))] = campo(ct);
    }
    higcit_destroy(it);

    // ------------------------------------------------ o ramo ON_BOUNDARY
    // Pontos EXATAMENTE sobre o plano do contorno, que e' o que dispara
    // in_domain == ON_BOUNDARY.  Sao tambem os pontos que o escritor VTK consulta.
    t_case("reproduz_campo_linear_SOBRE_o_contorno");
    for(int d = 0; d < DIM; d++) {
        for(int lado = 0; lado < 2; lado++) {
            for(int j = 1; j < LADO; j++) {
                Point org, x;
                for(int k = 0; k < DIM; k++) { org[k] = 0.5; x[k] = 0.5; }
                x[d]   = lado ? 1.0 : 0.0;          // SOBRE o plano
                x[(d + 1) % DIM] = j * h;           // varre ao longo da parede
                POINT_ASSIGN(org, x);
                org[d] = lado ? 1.0 - 0.5 * h : 0.5 * h;   // centro da celula vizinha
                stn_reset(stn);
                sd_get_stencil(sd, org, x, 1.0, stn);
                const real v = stn_mult_vector(stn, fval) - stn_get_rhs(stn);
                T_NEAR(v, campo(x), 1e-10, "campo linear sobre o contorno");
            }
        }
    }

    // --------------------------------------------- invariante de faixa
    // Interpolacao nao produz valor muito fora da faixa dos dados que usa.  O
    // limite e' generoso de proposito: serve para pegar erro GROSSEIRO de pesos,
    // que e' como o defeito de compactacao se manifestava (4,3x a 115x fora).
    t_case("nao_extrapola_grosseiramente_alem_do_suporte");
    real pior_fator = 1.0;
    Point pior_x;  POINT_ASSIGN_SCALAR(pior_x, 0.0);
    for(int d = 0; d < DIM; d++) {
        for(int lado = 0; lado < 2; lado++) {
            for(int j = 1; j < LADO; j++) {
                Point org, x;
                for(int k = 0; k < DIM; k++) { org[k] = 0.5; x[k] = 0.5; }
                x[d] = lado ? 1.0 : 0.0;
                x[(d + 1) % DIM] = j * h;
                POINT_ASSIGN(org, x);
                org[d] = lado ? 1.0 - 0.5 * h : 0.5 * h;
                stn_reset(stn);
                sd_get_stencil(sd, org, x, 1.0, stn);
                const real v = stn_mult_vector(stn, fval) - stn_get_rhs(stn);

                // faixa dos dados que o estencil de fato usou
                real lo = 1e30, hi2 = -1e30;
                for(int k = 0; k < stn_get_numelems(stn); k++) {
                    const real f = fval[stn_get_id(stn, k)];
                    if(f < lo) lo = f;
                    if(f > hi2) hi2 = f;
                }
                const real rhs = stn_get_rhs(stn);
                if(-rhs < lo) lo = -rhs;
                if(-rhs > hi2) hi2 = -rhs;
                if(hi2 <= lo) continue;

                const real amplitude = hi2 - lo;
                const real excesso = (v > hi2) ? (v - hi2) : (v < lo ? lo - v : 0.0);
                const real fator = 1.0 + excesso / amplitude;
                if(fator > pior_fator) { pior_fator = fator; POINT_ASSIGN(pior_x, x); }
            }
        }
    }
    T_BELOW(pior_fator, 2.0,
        "extrapolacao alem do suporte (1,0 = dentro da faixa; o defeito de "
        "compactacao dava 4,3 a 115)");
    if(pior_fator > 2.0) {
        _t_fail("pior ponto: (%.4f, %.4f)", pior_x[0], pior_x[1]);
    }

    free(fval);
    stn_destroy(stn);
    sd_destroy(sd);
    return t_end();
}
