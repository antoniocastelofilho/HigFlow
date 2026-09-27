// Adaptador entre o solver e o front-tracking -- o analogo de
// examples-common/fronteira-imersa.c, para a INTERFACE (tensao superficial), nao
// para o corpo rigido.  Compilado POR EXEMPLO (2D), fora da lib: o nucleo
// hig-flow-front-tracking.{c,h} nao toca o solver, e este arquivo faz a ponte.
//
// FASE B2 (gota estatica, lei de Laplace): a forca de tensao superficial
// F = sigma*kappa*n dos marcadores e' espalhada em dpFU (o campo de fonte que a
// equacao de momento le'), ANTES do preditor -- mesmo gancho que o corpo rigido
// usa.  Reusa a busca de faceta escalonada do corpo rigido (fi_suporte_facetas),
// para nao duplicar a parte delicada.
//
// SERIAL POR ENQUANTO: o reduce de franja (dono<-franja) NAO esta' aqui.  Em np=1
// nao ha' franja e o espalhamento e' completo; o paralelo e' trabalho da fase de
// escala (o mesmo padrao do baseline VOF, que ficou serial primeiro).

#include "hig-flow-kernel.h"
#include "hig-flow-fronteira-imersa.h"   // fi_suporte_facetas, fi_suporte_capacidade
#include "hig-flow-front-tracking.h"
#include "hig-mesh-snapshot.h"           // hig_facet_snapshot, hfs_center (para o oraculo)
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct {
    ft_frente *frente;
    real       sigma;
    int        advecta;      // 0 = frente fixa (Laplace puro); 1 = advecta
    int        passo;        // contador proprio, para o diagnostico periodico
    Point      centro_ini;   // centroide inicial, para medir deriva
    real       area_ini;     // area inicial, para medir conservacao
} ft_ctx;

// Tamanho de celula na posicao X (h do nucleo ali).  Malha uniforme: delta[0].
static real _h_em(sim_facet_domain *sfd, const Point X)
{
    hig_cell *c = sfd_get_cell_with_point(sfd, (real *) X);
    if (c == NULL) return 0.0;
    Point d; hig_get_delta(c, d);
    return d[0];
}

// Espalha F = sigma*kappa*n na malha escalonada (dpFU), ACUMULANDO.  Serial.
//
//   contribuicao a' faceta = F_dim * ds * (produto phi) / h^DIM
//
// Conservacao (o mesmo oraculo do standalone): SUM_faceta valor*h^DIM = SUM_k
// F_k*ds_k, porque o produto phi soma 1 no suporte.  Para curva fechada isso
// tende a zero (Laplace: forca liquida nula).
static void ft_espalha_tensao_solver(ft_frente *frente, real sigma,
                                     sim_facet_domain *sfd[DIM],
                                     distributed_property *dpF[DIM])
{
    int n = ft_num(frente);
    Point *pos = (Point *) malloc((size_t) n * sizeof(Point));
    Point *F   = (Point *) malloc((size_t) n * sizeof(Point));
    real  *ds  = (real  *) malloc((size_t) n * sizeof(real));
    ft_forcas_tensao(frente, sigma, pos, F, ds);

    int capac = fi_suporte_capacidade();
    int  *lids  = (int  *) malloc((size_t) capac * sizeof(int));
    real *pesos = (real *) malloc((size_t) capac * sizeof(real));

    for (int k = 0; k < n; k++) {
        real h = _h_em(sfd[0], pos[k]);
        if (h <= 0.0) continue;               // marcador fora do dominio local
        real hd = 1.0;
        for (int d = 0; d < DIM; d++) hd *= h;

        for (int dim = 0; dim < DIM; dim++) {
            int m = fi_suporte_facetas(sfd[dim], dim, pos[k], h, lids, pesos, capac);
            real esc = F[k][dim] * ds[k] / hd;
            for (int i = 0; i < m; i++)
                dp_add_value(dpF[dim], lids[i], esc * pesos[i]);
        }
    }

    // ORACULO DE CONSERVACAO, uma vez (FT_DIAG_CONSERVA): a soma da forca na
    // grade tem de igualar a dos marcadores -- SUM_faceta valor*h^DIM == SUM_k
    // F_k*ds_k -- porque o nucleo de Roma soma 1.  E' o teste do espalhamento
    // escalonado, independente das condicoes de contorno.  So' vale em serial
    // (sem franja) e antes do dp_sync.
    static int ja_conferiu = 0;
    if (!ja_conferiu && getenv("FT_DIAG_CONSERVA") != NULL) {
        ja_conferiu = 1;
        // Malha uniforme: h^DIM constante, tirado de um marcador dentro do dominio.
        real h = 0.0;
        for (int k = 0; k < n && h <= 0.0; k++) h = _h_em(sfd[0], pos[k]);
        real hd = 1.0; for (int d = 0; d < DIM; d++) hd *= h;
        for (int dim = 0; dim < DIM; dim++) {
            real lado_marc = 0.0;
            for (int k = 0; k < n; k++) lado_marc += F[k][dim] * ds[k];
            real lado_grade = 0.0;
            const hig_facet_snapshot *hfs = sfd_get_snapshot(sfd[dim]);
            for (int flid = 0; flid < hfs->n; flid++)
                lado_grade += dp_get_value(dpF[dim], flid) * hd;
            fprintf(stderr, "FT conserva dim=%d: grade=%.6e  marcadores=%.6e  "
                    "erro=%.3e\n", dim, lado_grade, lado_marc,
                    fabs(lado_grade - lado_marc));
        }
    }

    free(lids); free(pesos);
    free(pos); free(F); free(ds);

    // Serial: sem reduce de franja.  As facetas de franja em np>1 precisariam do
    // PetscSFReduce como em fi_espalha -- deixado para o paralelo.
    for (int dim = 0; dim < DIM; dim++) dp_sync(dpF[dim]);
}

// --------------------------------------------------------------------------
// ADVECCAO ACOPLADA: os marcadores andam com o fluido.
//
// E' a mudanca 1 das tres do projeto, e o que separa o front-tracking da
// fronteira imersa rigida: onde o corpo rigido INTERPOLA u para calcular a forca
// que impoe u=0, aqui interpola u para MOVER o marcador.  Mesma interpolacao,
// destino oposto.
// --------------------------------------------------------------------------

// Contexto do callback: os dominios e os campos de velocidade da malha.
typedef struct {
    sim_facet_domain     **sfd;   // ns->sfdu
    distributed_property **dpu;   // ns->dpu
} interp_ctx;

// Diagnostico da PARTICAO DA UNIDADE (FT_DIAG_INTERP).  Os pesos de Roma tem de
// somar 1 em cada marcador; se o suporte estiver incompleto -- marcador perto da
// fronteira, ou suporte atravessando nivel de refino -- a soma cai abaixo de 1 e
// a velocidade interpolada sai PEQUENA DEMAIS, em silencio.  Foi esta checagem
// que pegou a forca "exatamente pela metade" no corpo rigido: ela compara com um
// valor absoluto conhecido, coisa que a conservacao nao faz.
static real _pior_desvio_unidade = 0.0;

// Maior |u| interpolado nos marcadores.  NAO e' enfeite: uma gota estatica que
// fica parada e' exatamente o que se veria se a interpolacao devolvesse ZERO em
// silencio -- area perfeita, deriva nula, circularidade constante.  Este numero
// separa "adveccao funcionando, velocidade genuinamente minuscula" de "adveccao
// que nao faz nada".  Comparar com o Vmax do campo euleriano.
static real _maior_u_marcador = 0.0;

// u(X) interpolado da malha escalonada.  O campo fica CONGELADO em u^n durante o
// passo (nao existe u em t+dt/2), entao `t` e' ignorado: o RK2 do ft_advecta
// vira avaliacao de ponto medio no ESPACO, e o esquema segue de primeira ordem
// no tempo -- coerente com o resto do acoplamento explicito.
static void _campo_da_malha(const Point x, real t, void *vctx, real u[DIM])
{
    interp_ctx *ic = (interp_ctx *) vctx;
    (void) t;

    for (int d = 0; d < DIM; d++) u[d] = 0.0;

    real h = _h_em(ic->sfd[0], x);
    if (h <= 0.0) return;                 // fora do dominio: marcador nao anda

    int capac = fi_suporte_capacidade();
    int  *lids  = (int  *) malloc((size_t) capac * sizeof(int));
    real *pesos = (real *) malloc((size_t) capac * sizeof(real));

    for (int dim = 0; dim < DIM; dim++) {
        int m = fi_suporte_facetas(ic->sfd[dim], dim, x, h, lids, pesos, capac);
        real soma = 0.0, val = 0.0;
        // u(X) = SUM u(x) d_h(x-X) h^DIM, e o h^DIM cancela o 1/h^DIM do nucleo:
        // sobra a soma ponderada pelo produto de phi -- a MESMA de fi_interpola.
        for (int i = 0; i < m; i++) {
            val  += dp_get_value(ic->dpu[dim], lids[i]) * pesos[i];
            soma += pesos[i];
        }
        u[dim] = val;
        real desvio = fabs(soma - 1.0);
        if (desvio > _pior_desvio_unidade) _pior_desvio_unidade = desvio;
        if (fabs(val) > _maior_u_marcador) _maior_u_marcador = fabs(val);
    }

    free(lids); free(pesos);
}

// Metricas de forma da gota, para o oraculo do B2 dinamico: area, deriva do
// centroide e circularidade (4*pi*A/P^2, que vale 1 para o circulo).  Uma gota
// estatica em equilibrio tem de manter as tres.
static void _metricas(const ft_frente *f, const Point centro_ini,
                      real *area, real *deriva, real *circ)
{
    int n = ft_num(f);
    Point *p = (Point *) malloc((size_t) n * sizeof(Point));
    ft_posicoes(f, p);
    real cx = 0.0, cy = 0.0;
    for (int i = 0; i < n; i++) { cx += p[i][0]; cy += p[i][1]; }
    cx /= n; cy /= n;
    *deriva = sqrt((cx - centro_ini[0]) * (cx - centro_ini[0])
                 + (cy - centro_ini[1]) * (cy - centro_ini[1]));
    *area = ft_area(f);
    real P = ft_perimetro(f);
    *circ = (P > 0.0) ? 4.0 * M_PI * (*area) / (P * P) : 0.0;
    free(p);
}

// Gancho: roda antes do preditor.  Nesse ponto `ns->dpu` e' u^n -- a velocidade
// final, JA' PROJETADA (discretamente livre de divergencia), do passo anterior.
// Entao a ordem e': move a frente com u^n, faz a cirurgia, e espalha a forca de
// tensao superficial nas posicoes NOVAS, que e' o que o preditor vai ver.
//
// Adveccao com campo livre de divergencia preserva a area fechada -- e' o
// oraculo mais barato de que a interpolacao esta' certa.
static void _aplica_tensao(higflow_solver *ns, void *vctx)
{
    ft_ctx *ctx = (ft_ctx *) vctx;

    if (ctx->advecta) {
        interp_ctx ic = { ns->sfdu, ns->dpu };
        ft_advecta(ctx->frente, _campo_da_malha, &ic, ns->par.t, ns->par.dt);
        ft_cirurgia(ctx->frente);

        if (getenv("FT_DIAG_INTERP") != NULL && ctx->passo % 50 == 0) {
            real area, deriva, circ;
            _metricas(ctx->frente, ctx->centro_ini, &area, &deriva, &circ);
            fprintf(stderr, "FT passo %5d: n=%4d  area=%.8f (dA/A=%.2e)  "
                    "deriva=%.3e  circ=%.6f  unidade_pior=%.2e  "
                    "max|u_marc|=%.3e\n",
                    ctx->passo, ft_num(ctx->frente), (double) area,
                    (double) fabs(area - ctx->area_ini) / ctx->area_ini,
                    (double) deriva, (double) circ,
                    (double) _pior_desvio_unidade,
                    (double) _maior_u_marcador);
        }
        ctx->passo++;
    }

    ft_espalha_tensao_solver(ctx->frente, ctx->sigma, ns->sfdF, ns->dpFU);
}

//! Instala o front-tracking (tensao superficial) no solver.  `sigma` e' o
//! coeficiente; a frente e' o corpo ja' criado (ft_cria_circulo etc.).
//! extern "C": compilado como C++, mas o exemplo o declara com ligacao C.
extern "C" void front_tracking_instala(higflow_solver *ns, ft_frente *frente,
                                       real sigma)
{
    ft_ctx *ctx = (ft_ctx *) malloc(sizeof *ctx);
    ctx->frente = frente;
    ctx->sigma  = sigma;
    const char *sa = getenv("FT_ADVECTA");
    ctx->advecta = (sa != NULL) ? atoi(sa) : 0;   // padrao: frente fixa (Laplace)
    ctx->passo   = 0;

    // Referencias para as metricas do oraculo: area e centroide iniciais.
    {
        int n = ft_num(frente);
        Point *p = (Point *) malloc((size_t) n * sizeof(Point));
        ft_posicoes(frente, p);
        real cx = 0.0, cy = 0.0;
        for (int i = 0; i < n; i++) { cx += p[i][0]; cy += p[i][1]; }
        ctx->centro_ini[0] = cx / n;
        ctx->centro_ini[1] = cy / n;
        for (int d = 2; d < DIM; d++) ctx->centro_ini[d] = 0.0;
        ctx->area_ini = ft_area(frente);
        free(p);
    }

    higflow_set_fronteira_imersa(ns, _aplica_tensao, ctx);
}
