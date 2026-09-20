// C10 pela SEPARACAO certa: a malha fornece o suporte, o ajuste e' o mesmo.
//
// POR QUE ESTE TESTE EXISTE, E POR QUE ELE NAO CHAMA sd_get_stencil.
//
// A clausula C10 diz que interpolacao de ordem k reproduz polinomio de grau <= k.
// Quem reproduz e' o ajuste de minimos quadrados moveis, que pertence a
// Discretization e e' o MESMO para qualquer malha (`wls_set_points_and_calc`).  O
// que pertence a Mesh -- e o que uma segunda implementacao precisa entregar -- e'
// o SUPORTE: quais celulas ficam perto do ponto de consulta.
//
// Comparar `sd_get_stencil` do MTree com uma montagem propria do t8code mediria
// as duas coisas juntas e nao diria qual delas falhou.  Aqui os dois backends
// entregam apenas centros, os centros vao para o MESMO `wls`, e a reproducao e'
// verificada.  Se falhar com o suporte de um e acertar com o do outro, o defeito
// esta' na malha -- que e' exatamente o que se quer poder afirmar.
//
// O suporte do t8code e' colhido por busca em largura sobre VIZINHANCA DE FACE,
// nao por filtro de caixa sobre todas as folhas: filtro de caixa serviria igual
// para qualquer backend e nao exercitaria topologia nenhuma.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "higtree.h"
#include "higtree-iterator.h"
#include "domain.h"
#include "wls.h"
#include "utils.h"
#include "testing.h"

#ifdef HIGTREE_COM_T8CODE
#include "t8code/t8-stencil-support.h"
#endif

#define LADO   8
#define ORDEM  2
#define MAXPTS 64

// Campos de teste: um linear e um quadratico.  Ordem 2 reproduz os dois
// exatamente; e' a definicao da clausula.
static real linear(const Point p)    { return 1.0 + 2.0*p[0] + 3.0*p[1]; }
static real quadratico(const Point p) {
    return 1.0 + 2.0*p[0] - 0.5*p[1] + 0.25*p[0]*p[0] + 0.75*p[0]*p[1];
}

typedef struct {
    const char *nome;
    //! Centros do suporte em torno de `x`.  Devolve quantos escreveu.
    int (*suporte)(const Point x, int minpts, int maxpts, Point pts[]);
} FonteDeSuporte;

// ------------------------------------------------------------------ MTree
static hig_cell   *g_raiz = NULL;
static sim_domain *g_sd   = NULL;

static void monta_mtree(void) {
    if(g_sd != NULL) return;
    Point lo, hi; int nc[DIM];
    for(int d = 0; d < DIM; d++) { lo[d] = 0.0; hi[d] = 1.0; nc[d] = LADO; }
    g_raiz = hig_create_root(lo, hi);
    hig_refine_uniform(g_raiz, nc);
    g_sd = sd_create(NULL);
    sd_add_higtree(g_sd, g_raiz);
    mp_mapper *m = sd_get_domain_mapper(g_sd);
    higcit_celliterator *cit = sd_get_domain_celliterator(g_sd);
    mp_assign_from_celliterator(m, cit, 0);
    higcit_destroy(cit);
}

// Vizinhanca por caixa crescente em torno de x, que e' como a HiGTree colhe
// suporte (ver search_cells_in_tree_box em domain.c).
static int mtree_suporte(const Point x, int minpts, int maxpts, Point pts[]) {
    monta_mtree();
    const real h = 1.0 / LADO;
    for(int anel = 1; anel <= LADO; anel++) {
        int n = 0;
        const real r = anel * h;
        higcit_celliterator *cit;
        for(cit = sd_get_domain_celliterator(g_sd);
            !higcit_isfinished(cit) && n < maxpts; higcit_nextcell(cit)) {
            Point c; hig_get_center(higcit_getcell(cit), c);
            int dentro = 1;
            for(int d = 0; d < DIM; d++) if(fabs(c[d] - x[d]) > r + 1e-12) dentro = 0;
            if(dentro) { for(int d = 0; d < DIM; d++) pts[n][d] = c[d]; n++; }
        }
        higcit_destroy(cit);
        if(n >= minpts) return n;
    }
    return 0;
}

// ----------------------------------------------------------------- t8code
#ifdef HIGTREE_COM_T8CODE
static int t8_suporte(const Point x, int minpts, int maxpts, Point pts[]) {
    return t8_monta_suporte(1, x, minpts, maxpts, pts);
}
#endif

static const FonteDeSuporte FONTES[] = {
    { "mtree", mtree_suporte },
#ifdef HIGTREE_COM_T8CODE
    { "t8code", t8_suporte },
#endif
};

// Interpola `campo` em `x` usando o suporte da fonte e o wls compartilhado.
// Devolve 0 se a fonte nao entregou pontos bastante.
static int interpola(const FonteDeSuporte *f, wls_interpolator *wls,
                     real (*campo)(const Point), const Point x, real *saida,
                     int *npts) {
    Point pts[MAXPTS];
    const int minpts = wls_num_min_points(DIM, ORDEM);
    const int n = f->suporte(x, minpts, MAXPTS, pts);
    *npts = n;
    if(n < minpts) return 0;

    real w[MAXPTS];
    wls_set_points_and_calc(wls, n, pts, (real *) x, w);
    real v = 0.0;
    for(int i = 0; i < n; i++) v += w[i] * campo(pts[i]);
    *saida = v;
    return 1;
}

static void verifica(const FonteDeSuporte *f, wls_interpolator *wls) {
    const real h = 1.0 / LADO;

    t_case("suporte_tem_pontos_bastante");
    {
        Point x; for(int d = 0; d < DIM; d++) x[d] = 0.5 + 0.25 * h;
        Point pts[MAXPTS];
        const int minpts = wls_num_min_points(DIM, ORDEM);
        const int n = f->suporte(x, minpts, MAXPTS, pts);
        T_CHECK_MSG(n >= minpts,
            "[%s] o suporte trouxe %d pontos e a ordem %d em %dD precisa de %d",
            f->nome, n, ORDEM, DIM, minpts);
    }

    t_case("suporte_e_local_e_cai_na_grade");
    {
        // A REPRODUCAO POLINOMIAL NAO PEGA ISTO, e foi medido: deslocando todos
        // os centroides meia celula em x, o caso seguinte continua VERDE.  A
        // razao e' que o ajuste reproduz o polinomio onde quer que os pontos
        // estejam -- o campo e' avaliado nos proprios pontos, entao deslocar
        // amostra e coordenada juntas nao muda nada.  Ele verifica o AJUSTE, nao
        // se o suporte sao as celulas certas.
        //
        // Duas afirmacoes fecham o buraco: todo ponto de suporte cai na grade de
        // centros da malha, e todo ponto esta' PERTO do ponto de consulta.  A
        // segunda e' o que distingue estencil de "a malha inteira".
        Point x; for(int d = 0; d < DIM; d++) x[d] = 0.5 + 0.25 * h;
        Point pts[MAXPTS];
        const int minpts = wls_num_min_points(DIM, ORDEM);
        const int n = f->suporte(x, minpts, MAXPTS, pts);
        int fora_da_grade = 0, longe = 0;
        real pior = 0.0;
        for(int i = 0; i < n; i++) {
            for(int d = 0; d < DIM; d++) {
                const real k = pts[i][d] / h - 0.5;     // centro = (k + 1/2)*h
                if(fabs(k - round(k)) > 1e-9) { fora_da_grade++; break; }
            }
            real dist = 0.0;
            for(int d = 0; d < DIM; d++) {
                const real dd = fabs(pts[i][d] - x[d]);
                if(dd > dist) dist = dd;
            }
            if(dist > pior) pior = dist;
            if(dist > 3.0 * h) longe++;
        }
        T_CHECK_MSG(fora_da_grade == 0,
            "[%s] %d de %d pontos de suporte nao caem na grade de centros da "
            "malha (passo %.6f)", f->nome, fora_da_grade, n, h);
        T_CHECK_MSG(longe == 0,
            "[%s] %d de %d pontos de suporte estao a mais de 3 celulas do ponto "
            "de consulta; o mais distante a %.6f (= %.1f celulas).  Suporte que "
            "vira a malha inteira nao exercita topologia nenhuma",
            f->nome, longe, n, pior, pior / h);
    }

    t_case("reproduz_polinomio_de_grau_ate_a_ordem");
    {
        // Pontos no INTERIOR, longe do contorno: C10 e' sobre o interior; o
        // fechamento por condicao de contorno e' C12 e C14.
        int falhas = 0, sem_pontos = 0;
        char primeira[256]; primeira[0] = '\0';
        for(int i = 3; i <= LADO - 3; i++) {
            for(int j = 3; j <= LADO - 3; j++) {
                Point x;
                for(int d = 0; d < DIM; d++) x[d] = 0.5;
                x[0] = (i + 0.25) * h;      // fora do centro, para exigir ajuste
                x[1] = (j + 0.25) * h;
                real v; int n;
                if(!interpola(f, wls, linear, x, &v, &n)) { sem_pontos++; continue; }
                if(fabs(v - linear(x)) > 1e-10) {
                    if(!falhas) snprintf(primeira, sizeof primeira,
                        "linear em (%.4f, %.4f): obtido %.12f, esperado %.12f, "
                        "com %d pontos de suporte", x[0], x[1], v, linear(x), n);
                    falhas++;
                    continue;
                }
                if(!interpola(f, wls, quadratico, x, &v, &n)) continue;
                if(fabs(v - quadratico(x)) > 1e-10) {
                    if(!falhas) snprintf(primeira, sizeof primeira,
                        "quadratico em (%.4f, %.4f): obtido %.12f, esperado %.12f, "
                        "com %d pontos", x[0], x[1], v, quadratico(x), n);
                    falhas++;
                }
            }
        }
        T_CHECK_MSG(sem_pontos == 0,
            "[%s] %d pontos ficaram sem suporte suficiente", f->nome, sem_pontos);
        T_CHECK_MSG(falhas == 0,
            "[%s] %d ponto(s) nao reproduzem o polinomio.  %s",
            f->nome, falhas, primeira);
    }
}

int main(int argc, char *argv[]) {
    higtree_initialize(&argc, &argv);
    wls_interpolator *wls = wls_create(DIM, ORDEM, MAXPTS);
    for(unsigned i = 0; i < sizeof FONTES / sizeof *FONTES; i++) {
        verifica(&FONTES[i], wls);
    }
    wls_destroy(wls);
    return t_end();
}
