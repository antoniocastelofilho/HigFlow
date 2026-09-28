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
#include "hig-flow-eval.h"               // compute_center_p_left/right
#include "hig-flow-discret.h"            // compute_dpdx_at_point
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct {
    ft_frente *frente;
    real       sigma;
    int        advecta;      // 0 = frente fixa (Laplace puro); 1 = advecta
    int        balanceado;   // 1 = forca como sigma*kappa*grad(H) balanceado
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
                      real *area, real *deriva, real *circ, real *defor)
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
    // DEFORMACAO COM SINAL, pelos semieixos equivalentes (momentos de area).
    // Troca de sinal a cada meio periodo -- a circularidade, sempre <= 1, nao
    // trocaria, e oscilaria no dobro da frequencia.
    real sa, sb;
    ft_semieixos(f, &sa, &sb);
    *defor = (sa + sb > 0.0) ? (sa - sb) / (sa + sb) : 0.0;
    free(p);
}

// ---------------------------------------------------------------------------
// FORCA BALANCEADA (FT_BALANCEADO=1) -- o teste da hipotese do balanco.
//
// A forma padrao espalha sigma*kappa*n pelo nucleo de Roma direto nas facetas.
// Esse operador NAO e' o gradiente discreto de escalar nenhum que a projecao
// produza, e o descasamento entre ele e o gradiente de pressao aparece como
// corrente parasita -- e' a hipotese levantada pela comparacao com o VOF.
//
// Aqui a forca e' montada como
//
//     F = sigma * kappa * grad(H),
//
// com H a fracao de area da gota por celula (da GEOMETRIA da frente, via
// ft_area_na_caixa) e `grad` o MESMO operador que higflow_final_velocity usa
// para corrigir a velocidade:
//
//     compute_center_p_left/right  +  compute_dpdx_at_point
//
// Com kappa constante -- que e' o caso do circulo -- isso e' exatamente
// grad(sigma*kappa*H), isto e', a forca PERTENCE a' imagem do gradiente
// discreto, e a pressao pode cancela-la termo a termo.  Se a hipotese estiver
// certa, as correntes parasitas caem sem que o salto de pressao piore.
// Ver Francois et al. (2006).
// ---------------------------------------------------------------------------

// Campo H por celula, preso ao dominio corrente.  Criado sob demanda.
static distributed_property *_dpH = NULL;

static void _preenche_indicadora(higflow_solver *ns, ft_frente *frente)
{
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    mp_mapper  *mp  = sd_get_domain_mapper(sdp);
    if (_dpH == NULL) _dpH = psd_create_property(ns->psdp);

    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it);
         higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        int clid = mp_lookup(mp, hig_get_cid(c));
        if (clid < 0) continue;
        Point lo, hi;
        hig_get_lowpoint(c, lo);
        hig_get_highpoint(c, hi);
        real vol = 1.0;
        for (int d = 0; d < DIM; d++) vol *= (hi[d] - lo[d]);
        real a = ft_area_na_caixa(frente, lo, hi);
        dp_set_value(_dpH, clid, (vol > 0.0) ? a / vol : 0.0);
    }
    higcit_destroy(it);
    dp_sync(_dpH);
}

// Curvatura interpolada dos marcadores em `x`, por MEDIA PONDERADA com o nucleo
// de Roma de largura `h`.  So' e' consultada onde grad(H) e' nao nulo, a menos
// de uma celula da frente.
//
// A VERSAO ANTERIOR USAVA O MARCADOR MAIS PROXIMO, e isso custou uma corrida.
// No circulo era exato (kappa constante) e o teste passou espetacularmente.  Na
// ELIPSE, com kappa variavel, o campo do vizinho mais proximo e' descontinuo:
// salta quando a atribuicao troca de marcador.  Enquanto a gota ainda estava
// longe do circulo a forca fisica dominava e nada aparecia; assim que ela ficou
// quase circular (passo ~900, circ 0,9887) o artefato passou a dominar, a
// circularidade REVERTEU e a corrida foi a instabilidade.  A descontinuidade
// alimenta a adveccao, que move a frente, que muda a atribuicao -- realimentacao.
//
// Fora do suporte de qualquer marcador devolve 0; ali grad(H) tambem e' nulo.
static real _kappa_interp(const Point *pos, const Point *forca, int n,
                          real sigma, const Point x, real h)
{
    real num = 0.0, den = 0.0;
    for (int k = 0; k < n; k++) {
        real w = ft_delta_roma((pos[k][0] - x[0]) / h)
               * ft_delta_roma((pos[k][1] - x[1]) / h);
        if (w == 0.0) continue;
        // |forca| = sigma*kappa, por construcao de ft_forcas_tensao.
        real fx = forca[k][0], fy = forca[k][1];
        num += w * sqrt(fx * fx + fy * fy) / sigma;
        den += w;
    }
    return (den > 0.0) ? num / den : 0.0;
}

static void ft_espalha_tensao_balanceada(higflow_solver *ns, ft_frente *frente,
                                         real sigma,
                                         sim_facet_domain *sfd[DIM],
                                         distributed_property *dpF[DIM])
{
    _preenche_indicadora(ns, frente);

    int n = ft_num(frente);
    Point *pos = (Point *) malloc((size_t) n * sizeof(Point));
    Point *F   = (Point *) malloc((size_t) n * sizeof(Point));
    real  *ds  = (real  *) malloc((size_t) n * sizeof(real));
    ft_forcas_tensao(frente, sigma, pos, F, ds);

    sim_domain *sdp = psd_get_local_domain(ns->psdp);

    for (int dim = 0; dim < DIM; dim++) {
        mp_mapper *mu = sfd_get_domain_mapper(sfd[dim]);
        (void) mu;
        const hig_facet_snapshot *hfs = sfd_get_snapshot(sfd[dim]);
        for (int flid = 0; flid < hfs->n; flid++) {
            Point fcenter, fdelta;
            hfs_center(hfs, flid, fcenter);
            hfs_delta(hfs, flid, fdelta);

            // O MESMO par de amostragem que higflow_final_velocity usa para a
            // pressao -- e' o que torna a forca cancelavel pelo gradiente.
            real Hl = compute_center_p_left (sdp, fcenter, fdelta, dim, 0.5,
                                             _dpH, ns->stn);
            real Hr = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5,
                                             _dpH, ns->stn);
            real dHdx = compute_dpdx_at_point(fdelta, dim, 0.5, Hl, Hr);
            if (dHdx == 0.0) continue;          // longe da interface

            real hcel = _h_em(sfd[dim], fcenter);
            if (hcel <= 0.0) hcel = fdelta[dim];
            real kappa = _kappa_interp(pos, F, n, sigma, fcenter, hcel);
            dp_add_value(dpF[dim], flid, sigma * kappa * dHdx);
        }
    }

    free(pos); free(F); free(ds);
    for (int dim = 0; dim < DIM; dim++) dp_sync(dpF[dim]);
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

        static int cada = 0;
        if (cada == 0) { const char *c = getenv("FT_DIAG_CADA");
                         cada = (c != NULL) ? atoi(c) : 50; if (cada < 1) cada = 50; }
        if (getenv("FT_DIAG_INTERP") != NULL && ctx->passo % cada == 0) {
            real area, deriva, circ, defor;
            _metricas(ctx->frente, ctx->centro_ini, &area, &deriva, &circ, &defor);
            fprintf(stderr, "FT passo %5d: n=%4d  area=%.8f (dA/A=%.2e)  "
                    "deriva=%.3e  circ=%.6f  D=%+.6e  unidade_pior=%.2e  "
                    "max|u_marc|=%.3e\n",
                    ctx->passo, ft_num(ctx->frente), (double) area,
                    (double) fabs(area - ctx->area_ini) / ctx->area_ini,
                    (double) deriva, (double) circ, (double) defor,
                    (double) _pior_desvio_unidade,
                    (double) _maior_u_marcador);
        }
        // Despejo das POSICOES da frente (FT_DUMP_FRENTE=prefixo): um arquivo
        // por quadro com x y por marcador, fechando no primeiro.  E' o que as
        // figuras do relatorio desenham -- a forma, nao so' o numero.
        { const char *pre = getenv("FT_DUMP_FRENTE");
          static int cadaq = 0, quadro = 0;
          if (cadaq == 0) { const char *c = getenv("FT_DUMP_CADA");
                            cadaq = (c != NULL) ? atoi(c) : 100; if (cadaq < 1) cadaq = 100; }
          if (pre != NULL && ctx->passo % cadaq == 0) {
              char nome[512];
              snprintf(nome, sizeof nome, "%s_%04d.dat", pre, quadro++);
              FILE *fp = fopen(nome, "w");
              if (fp != NULL) {
                  int nn = ft_num(ctx->frente);
                  Point *pp = (Point *) malloc((size_t) nn * sizeof(Point));
                  ft_posicoes(ctx->frente, pp);
                  fprintf(fp, "x y\n");
                  for (int i = 0; i <= nn; i++)
                      fprintf(fp, "%.8f %.8f\n", pp[i % nn][0], pp[i % nn][1]);
                  fclose(fp); free(pp);
              }
          } }
        ctx->passo++;
    }

    if (ctx->balanceado)
        ft_espalha_tensao_balanceada(ns, ctx->frente, ctx->sigma,
                                     ns->sfdF, ns->dpFU);
    else
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
    { const char *s = getenv("FT_BALANCEADO"); ctx->balanceado = (s != NULL) ? atoi(s) : 0; }

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
