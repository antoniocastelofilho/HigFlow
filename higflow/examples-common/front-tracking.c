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
    int        advecta;   // 0 = frente fixa (Laplace puro); 1 = advecta
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

// Gancho: roda antes do preditor.  Opcionalmente advecta a frente no campo de
// velocidade atual (para medir correntes parasitas), depois espalha a forca.
static void _aplica_tensao(higflow_solver *ns, void *vctx)
{
    ft_ctx *ctx = (ft_ctx *) vctx;
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
    higflow_set_fronteira_imersa(ns, _aplica_tensao, ctx);
}
