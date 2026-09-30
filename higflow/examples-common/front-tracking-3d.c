// Adaptador entre o solver e o front-tracking 3D -- o analogo de
// examples-common/front-tracking.c, para a SUPERFICIE triangulada.  Compilado
// POR EXEMPLO (DIM=3), fora da lib: o nucleo hig-flow-front-tracking-3d.{c,h}
// nao toca o solver, e este arquivo faz a ponte.
//
// FASE B2 3D (gota estatica, lei de Laplace): a forca de tensao dos VERTICES e'
// espalhada em dpFU, ANTES do preditor -- mesmo gancho do corpo rigido.  A busca
// de faceta escalonada vem de fi_suporte_facetas, que ja' e' 3D e esta'
// verificada pelo example3d_SchaeferTurek: nao se duplica a parte delicada.
//
// UMA DIFERENCA DE UNIDADES QUE DARIA Dp ERRADO EM SILENCIO.  No adaptador 2D,
// ft_forcas_tensao devolve sigma*kappa -- forca por unidade de AREA --, e o
// espalhamento multiplica por `ds` (comprimento) para virar forca.  Aqui,
// ft3_forcas_tensao devolve a integral de linha sigma * contorno (t x n) ds, que
// JA' E' UMA FORCA em newtons.  Multiplicar pelo peso do vertice repetiria a
// integracao e erraria por um fator de area.
//
//     2D:  contribuicao = F_dim * ds   / h^2      (F em N/m^2)
//     3D:  contribuicao = F_dim        / h^3      (F em N)
//
// O oraculo de conservacao (FT3_DIAG_CONSERVA) existe exatamente para prender
// esse fator: SUM_faceta valor*h^3 tem de dar SUM_k F_k, e um erro de unidade
// aparece la' antes de aparecer no salto de pressao.
//
// SERIAL POR ENQUANTO, como o 2D comecou: em np=1 nao ha' franja e o
// espalhamento e' completo.  A regra de posse que o 2D acabou adotando
// (acumular so' em faceta propria, com a frente replicada) se transporta sem
// mudanca quando o paralelo entrar.

#include "hig-flow-kernel.h"
#include "hig-flow-fronteira-imersa.h"   // fi_suporte_facetas, fi_suporte_capacidade
#include "hig-flow-front-tracking-3d.h"
#include "hig-mesh-snapshot.h"           // hig_facet_snapshot (para o oraculo)
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct {
    ft3_superficie *sup;
    real            sigma;
    int             advecta;   // 0 = superficie fixa (Laplace puro)
    int             passo;
} ft3_ctx;

// Tamanho de celula na posicao X.  Malha uniforme: delta[0].
static real _h_em(sim_facet_domain *sfd, const Point X)
{
    hig_cell *c = sfd_get_cell_with_point(sfd, (real *) X);
    if (c == NULL) return 0.0;
    Point d; hig_get_delta(c, d);
    return d[0];
}

static int _posse_regra(void)
{
    static int r = -1;
    if (r < 0) { const char *s = getenv("FT3_POSSE"); r = (s != NULL) ? atoi(s) : 0; }
    return r;
}

static void ft3_espalha_tensao_solver(ft3_superficie *sup, real sigma,
                                      sim_facet_domain *sfd[DIM],
                                      distributed_property *dpF[DIM])
{
    const int n = ft3_num_vertices(sup);
    Point *pos  = (Point *) malloc((size_t) n * sizeof(Point));
    Point *F    = (Point *) malloc((size_t) n * sizeof(Point));
    real  *peso = (real  *) malloc((size_t) n * sizeof(real));
    ft3_forcas_tensao(sup, sigma, pos, F, peso);

    int capac = fi_suporte_capacidade();
    int  *lids  = (int  *) malloc((size_t) capac * sizeof(int));
    real *pesos = (real *) malloc((size_t) capac * sizeof(real));

    for (int k = 0; k < n; k++) {
        real h = _h_em(sfd[0], pos[k]);
        if (h <= 0.0) continue;               // vertice fora do dominio local
        real hd = h * h * h;                  // DIM=3

        for (int dim = 0; dim < DIM; dim++) {
            int m = fi_suporte_facetas(sfd[dim], dim, pos[k], h, lids, pesos, capac);
            // SEM multiplicar por peso[k]: F ja' e' forca.  Ver o cabecalho.
            real esc = F[k][dim] / hd;
            for (int i = 0; i < m; i++) {
                if (_posse_regra() == 2 && lids[i] >= fi_dp_num_proprias(dpF[dim]))
                    continue;                 // faceta de franja: nao e' minha
                dp_add_value(dpF[dim], lids[i], esc * pesos[i]);
            }
        }
    }

    // ORACULO DE CONSERVACAO (FT3_DIAG_CONSERVA), uma vez.  SUM_faceta valor*h^3
    // tem de igualar SUM_k F_k, porque o nucleo de Roma soma 1 no suporte.  E'
    // o que prende o fator de unidade.  So' vale em serial (sem franja) e antes
    // do dp_sync.  Numa superficie FECHADA os dois lados tendem a zero, entao o
    // que se compara e' o erro contra a ESCALA das contribuicoes, e nao contra
    // um total que e' zero por simetria.
    static int ja_conferiu = 0;
    if (!ja_conferiu && getenv("FT3_DIAG_CONSERVA") != NULL) {
        ja_conferiu = 1;
        real h = 0.0;
        for (int k = 0; k < n && h <= 0.0; k++) h = _h_em(sfd[0], pos[k]);
        real hd = h * h * h;
        for (int dim = 0; dim < DIM; dim++) {
            real lado_marc = 0.0, escala = 0.0;
            for (int k = 0; k < n; k++) {
                lado_marc += F[k][dim];
                escala    += fabs(F[k][dim]);
            }
            real lado_grade = 0.0;
            const hig_facet_snapshot *hfs = sfd_get_snapshot(sfd[dim]);
            for (int flid = 0; flid < hfs->n; flid++)
                lado_grade += dp_get_value(dpF[dim], flid) * hd;
            fprintf(stderr, "FT3 conserva dim=%d: grade=%.6e  vertices=%.6e  "
                    "erro=%.3e  erro/escala=%.3e\n", dim, lado_grade, lado_marc,
                    fabs(lado_grade - lado_marc),
                    escala > 0.0 ? fabs(lado_grade - lado_marc) / escala : 0.0);
        }
    }

    free(lids); free(pesos); free(pos); free(F); free(peso);
}

// Gancho chamado a cada passo, no mesmo ponto em que o corpo rigido injeta.
static void _aplica_tensao(higflow_solver *ns, void *vctx)
{
    ft3_ctx *ctx = (ft3_ctx *) vctx;

    if (ctx->advecta) {
        // Nao ha' campo analitico aqui: a velocidade vem da malha.  Fica para a
        // fase seguinte, junto com o teste que a exercita -- o B2 e' a gota
        // ESTATICA, e mover a superficie sem oraculo de movimento seria
        // exatamente o que os degraus anteriores existem para evitar.
        fprintf(stderr, "FT3: adveccao acoplada ainda nao implementada "
                        "(B2 e' a gota estatica)\n");
        abort();
    }

    // Diagnostico periodico: a superficie nao deve se mexer no B2, entao area,
    // volume e topologia sao constantes -- qualquer deriva aqui e' defeito.
    if ((ctx->passo % 50) == 0) {
        print0f("=+= FT3 passo %d: nv=%d nt=%d area=%.8f volume=%.8f chi=%d =+=\n",
                ctx->passo, ft3_num_vertices(ctx->sup), ft3_num_triangulos(ctx->sup),
                (double) ft3_area(ctx->sup), (double) ft3_volume(ctx->sup),
                ft3_euler(ctx->sup));
    }
    ctx->passo++;

    ft3_espalha_tensao_solver(ctx->sup, ctx->sigma, ns->sfdF, ns->dpFU);
}

//! Instala o front-tracking 3D (tensao superficial) no solver.
extern "C" void front_tracking_3d_instala(higflow_solver *ns, ft3_superficie *sup,
                                          real sigma)
{
    ft3_ctx *ctx = (ft3_ctx *) malloc(sizeof *ctx);
    ctx->sup     = sup;
    ctx->sigma   = sigma;
    ctx->passo   = 0;
    { const char *s = getenv("FT3_ADVECTA"); ctx->advecta = (s != NULL) ? atoi(s) : 0; }
    higflow_set_fronteira_imersa(ns, _aplica_tensao, ctx);
}
