// O VALOR que o estencil produz, nao a forma dele.
//
// Por que este arquivo existe: a suite ATF de 2020 ja' chamava `sd_get_stencil`
// em pontos fora do dominio, com contornos montados -- e verificava
//
//     ATF_CHECK(stn_get_numelems(stn) >= 3);
//     ATF_CHECK(FLT_NE(stn_get_rhs(stn), 0.0));
//
// isto e', "tem pelo menos tres elementos" e "o lado direito nao e' zero".  Os sete
// defeitos de fechamento corrigidos em 2026-09-15/16 e os dois de 2026-09-17
// passariam por essas duas assercoes sem excecao: todos produziam estencil com a
// FORMA certa e os NUMEROS errados.  Uma assercao so' vale se o defeito a violaria.
//
// A ORACULO usado aqui e' a reproducao polinomial: uma interpolacao de ordem k
// reproduz exatamente qualquer polinomio de grau <= k.  Preenchendo o campo e os
// valores de contorno com o mesmo polinomio analitico, o valor interpolado tem de
// bater com o analitico a precisao de maquina, em qualquer ponto -- inclusive fora
// do dominio, onde quem responde e' o fechamento por condicao de contorno.
//
// SUTILEZA QUE DECIDE O DESENHO, e que custou uma iteracao para perceber: exatidao
// polinomial NAO pega erro de SELECAO de parede.  Se o campo e' linear e os valores
// de contorno saem do mesmo campo, fechar contra a parede errada ainda reproduz o
// campo exatamente -- a interpolacao de Lagrange ao longo de qualquer eixo acerta um
// campo linear.  Por isso ha' dois campos aqui: o de grau <= ordem, para exatidao, e
// um de grau ACIMA, em que o erro da parede certa e' O(h^3) e o de uma parede
// distante e' ordens de grandeza maior.  Ver `campo_cubico`.

#include <stdlib.h>
#include <string.h>
#include "testing.h"
#include "domain.h"
#include "utils.h"

#define LADO   8               // celulas por direcao
#define LO    (-1.0)
#define HI    ( 1.0)

// ---------------------------------------------------------------- campos

// Grau 1: reproduzido EXATAMENTE por qualquer interpolacao de ordem >= 1.
static real campo_linear(const Point p) {
    real v = 0.5;
    const real coef[3] = { 2.0, -3.0, 1.5 };
    for(int d = 0; d < DIM; d++) v += coef[d] * p[d];
    return v;
}

// Grau 3: NAO reproduzido.  O erro de um fechamento correto e' O(h^3); o de um
// fechamento que escolhe parede errada e' muito maior, e e' isso que se mede.
static real campo_cubico(const Point p) {
    real v = 0.0;
    for(int d = 0; d < DIM; d++) v += p[d]*p[d]*p[d] - 0.7*p[d]*p[d];
    return v;
}

// ---------------------------------------------------------------- montagem

typedef struct {
    hig_cell    *root;
    mp_mapper   *m;
    sim_domain  *sd;
    real        *fval;        // campo nos centros das celulas, indexado pelo id local
    int          n;
} caixa;

// Dominio [LO,HI]^DIM com LADO celulas por direcao e Dirichlet nas 2*DIM faces,
// com os valores de contorno tirados do MESMO campo analitico.
static caixa monta(real (*campo)(const Point)) {
    caixa c;
    Point l, h;
    for(int d = 0; d < DIM; d++) { l[d] = LO; h[d] = HI; }
    c.root = hig_create_root(l, h);

    int nc[DIM];
    for(int d = 0; d < DIM; d++) nc[d] = LADO;
    hig_refine_uniform(c.root, nc);

    c.m = mp_create();
    higcit_celliterator *it = higcit_create_all_leaves(c.root);
    c.n = mp_assign_from_celliterator(c.m, it, 0);
    higcit_destroy(it);

    c.sd = sd_create(c.m);
    sd_add_higtree(c.sd, c.root);

    // Campo nos centros
    c.fval = (real *) calloc(c.n, sizeof *c.fval);
    for(it = higcit_create_all_leaves(c.root); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *cell = higcit_getcell(it);
        Point ct; hig_get_center(cell, ct);
        c.fval[mp_lookup(c.m, hig_get_id(cell, 0))] = campo(ct);
    }
    higcit_destroy(it);

    // Uma face Dirichlet por lado, em cada direcao
    for(int d = 0; d < DIM; d++) {
        for(int lado = 0; lado < 2; lado++) {
            Point bl, bh; int bnc[DIM];
            for(int k = 0; k < DIM; k++) {
                bl[k] = LO;   bh[k] = HI;   bnc[k] = LADO;
            }
            const real plano = lado ? HI : LO;
            bl[d] = plano - EPSDELTA;
            bh[d] = plano + EPSDELTA;
            bnc[d] = 1;

            hig_cell *bc = hig_create_root(bl, bh);
            hig_refine_uniform(bc, bnc);
            mp_mapper *bm = mp_create();
            higcit_celliterator *bit = higcit_create_all_leaves(bc);
            mp_assign_from_celliterator(bm, bit, 0);
            higcit_destroy(bit);

            sim_boundary *sb = sb_create(bc, DIRICHLET, bm);
            for(bit = higcit_create_all_leaves(bc); !higcit_isfinished(bit); higcit_nextcell(bit)) {
                hig_cell *bcell = higcit_getcell(bit);
                Point bct; hig_get_center(bcell, bct);
                sb_set_value(sb, mp_lookup(bm, hig_get_id(bcell, 0)), campo(bct));
            }
            higcit_destroy(bit);
            sd_add_boundary(c.sd, sb);
        }
    }
    return c;
}

static void desmonta(caixa *c) {
    sd_destroy(c->sd);
    free(c->fval);
}

// Valor interpolado no ponto x, com o estencil construido a partir de `origem`.
// Convencao de dp_interpolate_from_stencil (pdomain.c): -rhs + soma(w_i * f_i).
static real interpola(caixa *c, const Point origem, const Point x, sim_stencil *stn) {
    stn_reset(stn);
    sd_get_stencil(c->sd, origem, x, 1.0, stn);
    return stn_mult_vector(stn, c->fval) - stn_get_rhs(stn);
}

int main(void) {
    sim_stencil *stn = stn_create();
    const real h = (HI - LO) / LADO;

    // ---------------------------------------------------------- exatidao
    // Campo linear reproduzido a precisao de maquina em pontos INTERNOS que nao
    // sao centro de celula (centro de celula sai por acerto exato e nao exercita
    // a interpolacao).
    {
        caixa c = monta(campo_linear);
        t_case("linear_interior");
        for(int i = 1; i < LADO; i++) {
            Point org, x;
            for(int d = 0; d < DIM; d++) { org[d] = 0.0; x[d] = 0.0; }
            org[0] = LO + (i - 0.5) * h;
            x[0]   = LO + i * h;              // face entre duas celulas
            const real obtido = interpola(&c, org, x, stn);
            T_NEAR(obtido, campo_linear(x), 1e-10, "campo linear no interior");
        }
        desmonta(&c);
    }

    // Campo linear reproduzido tambem FORA do dominio, onde quem responde e' o
    // fechamento por condicao de contorno.  E' o caminho dos sete defeitos.
    {
        caixa c = monta(campo_linear);
        t_case("linear_fora_do_dominio");
        for(int d = 0; d < DIM; d++) {
            for(int lado = 0; lado < 2; lado++) {
                Point org, x;
                for(int k = 0; k < DIM; k++) { org[k] = 0.0; x[k] = 0.0; }
                const real plano = lado ? HI : LO;
                const real passo = lado ? h : -h;
                org[d] = plano - 0.5 * passo;     // ultima celula, dentro
                x[d]   = plano + 0.5 * passo;     // meia celula FORA
                const real obtido = interpola(&c, org, x, stn);
                T_NEAR(obtido, campo_linear(x), 1e-10, "campo linear fora do dominio");
            }
        }
        desmonta(&c);
    }

    // ------------------------------------------------- selecao da parede
    // Campo cubico: nao e' reproduzido, e o erro revela QUAL parede fechou.  Um
    // fechamento correto erra O(h^3); fechar contra parede de outra direcao, ou
    // contra uma distante, erra muito mais.  O limite abaixo e' generoso de
    // proposito -- serve para pegar erro grosseiro de selecao, nao para medir ordem.
    {
        caixa c = monta(campo_cubico);
        t_case("cubico_fora_do_dominio_pega_parede_errada");
        const real limite = 40.0 * h * h * h;
        for(int d = 0; d < DIM; d++) {
            for(int lado = 0; lado < 2; lado++) {
                Point org, x;
                for(int k = 0; k < DIM; k++) { org[k] = 0.0; x[k] = 0.0; }
                const real plano = lado ? HI : LO;
                const real passo = lado ? h : -h;
                org[d] = plano - 0.5 * passo;
                x[d]   = plano + 0.5 * passo;
                const real erro = fabs(interpola(&c, org, x, stn) - campo_cubico(x));
                T_BELOW(erro, limite, "erro do campo cubico fora do dominio");
            }
        }
        desmonta(&c);
    }

    stn_destroy(stn);
    return t_end();
}
