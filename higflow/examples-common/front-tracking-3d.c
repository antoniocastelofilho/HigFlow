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
#include "criterio-esteira-pura.h"
#include "hig-flow-fronteira-imersa.h"   // fi_suporte_facetas, fi_suporte_capacidade
#include "hig-flow-front-tracking-3d.h"
#include "hig-flow-amr-criterio.h"   // o nucleo do criterio, compartilhado com o VOF
#include "hig-mesh-snapshot.h"           // hig_facet_snapshot (para o oraculo)
#include "hig-flow-eval.h"               // compute_facet_value_at_point (vorticidade)
#include "hig-flow-discret.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct {
    ft3_superficie *sup;
    real            sigma;
    int             advecta;   // 0 = superficie fixa (Laplace puro)
    int             passo;
    double          t_fracvol; // cronometro da fracao (FT3_CRONO_FRACVOL)
    real            vol_ini;   // volume no instante em que a adveccao comecou
    int             cirurgias; // operacoes de cirurgia acumuladas
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

// ---------------------------------------------------------------------------
// MULTIFASICO: a superficie manda na fracao, e as propriedades saem dela.
//
// Mesma forma do 2D: em vez de advectar a fracao com PLIC, ela e' RECALCULADA
// da geometria a cada passo, e dai' para frente rho(x), mu(x) e o momento de
// coeficiente variavel sao o MESMO codigo que o VOF usa.  A comparacao entre os
// dois passa a isolar exatamente a representacao da interface.
//
// CONVENCAO DO SOLVER, conferida no codigo e nao suposta:
//     dens = (1 - fracvol)*dens0 + fracvol*dens1
// ou seja fracvol=1 -> FASE 1, e fase 1 e' DENTRO.  `ft3_fracao_na_caixa`
// devolve justamente a fracao de dentro, entao as duas convencoes batem sem
// inversao.  No 2D esta suposicao foi feita ao contrario uma vez e produziu uma
// bolha PESADA num meio leve, que afundou ate' a base.
// ---------------------------------------------------------------------------
static void _preenche_fracvol(higflow_solver *ns, ft3_superficie *sup,
                              real *vol_no_campo, real faixa[2])
{
    real soma_vol = 0.0, fmin = 1e300, fmax = -1e300;
    // GRADE por chamada.  A superficie pode ter mudado desde o passo anterior
    // (adveccao, cirurgia), e uma grade sobrevivente apontaria para triangulos
    // que nao existem mais.  Construi-la custa O(nt); a alternativa -- invalidar
    // por versao -- exige que TODO caminho que mexe na superficie lembre de
    // marcar, e esquecer um nao da' erro, da' resultado errado.
    ft3_grade *grade = ft3_grade_cria(sup);
    sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
    mp_mapper  *mp  = sd_get_domain_mapper(sdm);
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(sdm); !higcit_isfinished(it);
         higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        int clid = mp_lookup(mp, hig_get_cid(c));
        if (clid < 0) continue;
        Point lo, hi;
        hig_get_lowpoint(c, lo);
        hig_get_highpoint(c, hi);
        real f = ft3_fracao_na_caixa_g(sup, grade, lo, hi);
        if (f < 0.0) f = 0.0;
        if (f > 1.0) f = 1.0;
        dp_set_value(ns->ed.mult.dpfracvol, clid, f);
        real v = 1.0;
        for (int d = 0; d < DIM; d++) v *= (hi[d] - lo[d]);
        soma_vol += f * v;
        if (f > fmax) fmax = f;
        if (f < fmin) fmin = f;
    }
    higcit_destroy(it);
    ft3_grade_destroi(grade);
    dp_sync(ns->ed.mult.dpfracvol);
    if (vol_no_campo) *vol_no_campo = soma_vol;
    if (faixa) { faixa[0] = fmin; faixa[1] = fmax; }
}

// ---------------------------------------------------------------------------
// ADVECCAO ACOPLADA: a superficie anda com a velocidade INTERPOLADA da malha.
//
// Mesma estrutura do 2D, inclusive o FILTRO DE POSSE, que la' foi o que separou
// interpolacao correta de interpolacao incompleta.  Em np=1 nao ha' franja e o
// filtro nunca recusa; ele esta' aqui porque o caminho paralelo e' o mesmo e
// escrever depois custaria o mesmo trabalho duas vezes.
//
// O ORACULO E' O VOLUME.  Campo discretamente livre de divergencia preserva o
// volume fechado, entao deriva de volume mede diretamente se a interpolacao esta'
// certa -- e' o mesmo oraculo do B1 3D, agora com a velocidade vindo da malha em
// vez de uma formula.
// ---------------------------------------------------------------------------
typedef struct {
    sim_facet_domain     *sfd[DIM];
    distributed_property *dpu[DIM];
} ft3_interp;

static real _pior_desvio_unidade_3d = 0.0;

static void _campo_da_malha_3d(const Point *x, int n, real t, void *vctx, Point *u)
{
    ft3_interp *ic = (ft3_interp *) vctx;
    (void) t;
    const int capac = fi_suporte_capacidade();
    int  *lids  = (int  *) malloc((size_t) capac * sizeof(int));
    real *pesos = (real *) malloc((size_t) capac * sizeof(real));

    real *val_loc = (real *) calloc((size_t) n * DIM, sizeof(real));
    real *cnt_loc = (real *) calloc((size_t) n * DIM, sizeof(real));

    for (int k = 0; k < n; k++) {
        const real h = _h_em(ic->sfd[0], x[k]);
        if (h <= 0.0) continue;              // este rank nao enxerga: contribui 0
        for (int dim = 0; dim < DIM; dim++) {
            const int m = fi_suporte_facetas(ic->sfd[dim], dim, x[k], h,
                                             lids, pesos, capac);
            real val = 0.0, soma = 0.0;
            for (int i = 0; i < m; i++) {
                val  += dp_get_value(ic->dpu[dim], lids[i]) * pesos[i];
                soma += pesos[i];
            }
            // PARTICAO DA UNIDADE.  Os pesos de Roma somam 1 em aritmetica
            // exata; soma diferente de 1 quer dizer suporte INCOMPLETO, e a
            // interpolacao nao vale.  A tolerancia e' folgada de proposito: o
            // desvio observado em serie e' da ordem de 1e-14.
            if (fabs(soma - 1.0) < 1e-9) {
                val_loc[k*DIM + dim] = val;
                cnt_loc[k*DIM + dim] = 1.0;
            }
            const real desvio = fabs(soma - 1.0);
            if (desvio > _pior_desvio_unidade_3d) _pior_desvio_unidade_3d = desvio;
        }
    }

    int ntasks = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    real *val_g = val_loc, *cnt_g = cnt_loc;
    if (ntasks > 1) {
        val_g = (real *) malloc((size_t) n * DIM * sizeof(real));
        cnt_g = (real *) malloc((size_t) n * DIM * sizeof(real));
        MPI_Allreduce(val_loc, val_g, n * DIM, MPI_HIGREAL, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(cnt_loc, cnt_g, n * DIM, MPI_HIGREAL, MPI_SUM, MPI_COMM_WORLD);
    }

    for (int k = 0; k < n; k++)
        for (int dim = 0; dim < DIM; dim++) {
            const real c = cnt_g[k*DIM + dim];
            if (c < 0.5) {
                // NENHUM processo tem o suporte inteiro deste vertice nesta
                // direcao.  Nao ha' velocidade correta a devolver, e devolver
                // zero seria pregar o vertice no lugar em silencio -- que e'
                // exatamente o defeito que o 2D levou meio dia para achar.
                fprintf(stderr, "FT3: vertice %d sem suporte completo em dim=%d "
                        "(x = %g %g %g).  Abortando em vez de devolver zero.\n",
                        k, dim, (double) x[k][0], (double) x[k][1], (double) x[k][2]);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            u[k][dim] = val_g[k*DIM + dim] / c;
        }

    if (ntasks > 1) { free(val_g); free(cnt_g); }
    free(val_loc); free(cnt_loc); free(lids); free(pesos);
}

// Gancho chamado a cada passo, no mesmo ponto em que o corpo rigido injeta.
static void _aplica_tensao(higflow_solver *ns, void *vctx)
{
    ft3_ctx *ctx = (ft3_ctx *) vctx;

    if (ctx->advecta) {
        // ORDEM.  Neste ponto `ns->dpu` e' u^n -- a velocidade final do passo
        // anterior, JA' PROJETADA e portanto discretamente livre de divergencia.
        // Entao: move com u^n, faz a cirurgia, e espalha a forca nas posicoes
        // NOVAS, que sao as que o preditor vai ver.
        ft3_interp ic;
        for (int d = 0; d < DIM; d++) {
            ic.sfd[d] = psfd_get_local_domain(ns->psfdu[d]);
            ic.dpu[d] = ns->dpu[d];
        }
        if (ctx->vol_ini <= 0.0) ctx->vol_ini = ft3_volume(ctx->sup);
        ft3_advecta_lote(ctx->sup, _campo_da_malha_3d, &ic, ns->par.t, ns->par.dt);
        ctx->cirurgias += ft3_cirurgia(ctx->sup);
    }

    // Diagnostico periodico: a superficie nao deve se mexer no B2, entao area,
    // volume e topologia sao constantes -- qualquer deriva aqui e' defeito.
    if ((ctx->passo % 50) == 0) {
        real v = ft3_volume(ctx->sup);
        // O ORACULO DA ADVECCAO: campo livre de divergencia preserva o volume
        // fechado.  A deriva mede a interpolacao, e nada mais.
        real deriva = (ctx->vol_ini > 0.0) ? (v - ctx->vol_ini) / ctx->vol_ini : 0.0;
        print0f("=+= FT3 passo %d: nv=%d nt=%d area=%.8f volume=%.8f chi=%d "
                "deriva=%+.3e cirurgias=%d unidade=%.2e =+=\n",
                ctx->passo, ft3_num_vertices(ctx->sup), ft3_num_triangulos(ctx->sup),
                (double) ft3_area(ctx->sup), (double) v, ft3_euler(ctx->sup),
                (double) deriva, ctx->cirurgias, (double) _pior_desvio_unidade_3d);
    }
    ctx->passo++;

    // A fracao tem de ser preenchida ANTES de higflow_compute_viscosity /
    // density_multiphase, e e' exatamente aqui que o gancho foi posto no passo
    // multifasico.
    if (ns->contr.flowtype == MULTIPHASE) {
        real t0 = 0.0, t1 = 0.0;
        const int cronometra = (getenv("FT3_CRONO_FRACVOL") != NULL);
        if (cronometra) t0 = MPI_Wtime();
        real vol_campo = 0.0, faixa[2] = {0.0, 0.0};
        _preenche_fracvol(ns, ctx->sup, &vol_campo, faixa);
        // ORACULO DO ACOPLAMENTO: o volume que o CAMPO carrega tem de ser o
        // volume da SUPERFICIE.  Sem isto, "a fracao foi calculada" e "a fracao
        // chegou ao solver" sao duas afirmacoes diferentes, e so' a primeira
        // estaria medida.
        if ((ctx->passo % 50) == 1) {
            real vs = ft3_volume(ctx->sup);
            print0f("=+= FT3 fracvol: volume no campo=%.8f  da superficie=%.8f  "
                    "erro=%.3e  faixa=[%.4f , %.4f] =+=\n",
                    (double) vol_campo, (double) vs,
                    (double) (vs != 0.0 ? fabs(vol_campo - vs)/vs : 0.0),
                    (double) faixa[0], (double) faixa[1]);
        }
        if (cronometra) {
            t1 = MPI_Wtime();
            ctx->t_fracvol += (t1 - t0);
            if ((ctx->passo % 50) == 0)
                print0f("=+= FT3 fracvol: %.4f s acumulados em %d passos "
                        "(%.4f s/passo) =+=\n",
                        ctx->t_fracvol, ctx->passo, ctx->t_fracvol / ctx->passo);
        }
    }

    ft3_espalha_tensao_solver(ctx->sup, ctx->sigma, ns->sfdF, ns->dpFU);
}

//! Preenche `dpfracvol` a partir da GEOMETRIA da superficie, uma vez, para a
//! CONDICAO INICIAL.  Sem isto o campo nasce com o degrau suavizado que a
//! funcao `fracvol` do exemplo devolve, e o primeiro quadro -- aquele que fixa
//! a massa de referencia -- carrega o erro dele.
//!
//! Importa para a comparacao FT contra VOF: o VOF inicializa a fracao por
//! subdivisao subcelular, com erro de volume 6.1e-05 em N=16, enquanto o
//! degrau erra 1.2e-02.  Se um dos dois nascer com 1.2% de massa a mais, toda
//! comparacao de conservacao ou de velocidade de subida mede a diferenca de
//! condicao inicial, nao a diferenca de metodo.
//!
//! Aqui a fracao vem da forma fechada do corte plano-caixa aplicada a' propria
//! superficie triangulada (ft3_fracao_na_caixa), que e' a maquinaria do FT --
//! nao a do VOF.  O erro residual passa a ser o da TRIANGULACAO, que e' a
//! representacao que o metodo de fato usa.
extern "C" void front_tracking_3d_fracvol_inicial(higflow_solver *ns,
                                                  ft3_superficie *sup)
{
    if (ns->contr.flowtype != MULTIPHASE) return;
    real vol_campo = 0.0, faixa[2] = {0.0, 0.0};
    _preenche_fracvol(ns, sup, &vol_campo, faixa);
    real vol_sup = ft3_volume(sup);
    print0f("=+=+=+= FT3 fracvol inicial da SUPERFICIE: vol campo=%.8f  "
            "vol superficie=%.8f  dif=%.3e  fracao em [%.4f, %.4f] =+=+=+=\n",
            (double) vol_campo, (double) vol_sup,
            (double) (vol_sup != 0.0 ? fabs(vol_campo-vol_sup)/vol_sup : 0.0),
            (double) faixa[0], (double) faixa[1]);
}

//! Instala o front-tracking 3D (tensao superficial) no solver.
extern "C" void front_tracking_3d_instala(higflow_solver *ns, ft3_superficie *sup,
                                          real sigma)
{
    ft3_ctx *ctx = (ft3_ctx *) malloc(sizeof *ctx);
    ctx->sup     = sup;
    ctx->sigma   = sigma;
    ctx->passo   = 0;
    ctx->t_fracvol = 0.0;
    ctx->vol_ini   = 0.0;
    ctx->cirurgias = 0;
    { const char *s = getenv("FT3_ADVECTA"); ctx->advecta = (s != NULL) ? atoi(s) : 0; }
    higflow_set_fronteira_imersa(ns, _aplica_tensao, ctx);
}

// ---------------------------------------------------------------------------
// CRITERIO DE REFINO EM 3D: escreve um .amr multinivel a partir da geometria
// da superficie.
//
// Porte do criterio 2D (examples-common/front-tracking.c), com a MESMA regra --
// celula com 0<fracvol<1 e' SEMENTE, e uma celula recebe o nivel L se estiver a
// menos de `limiar[L-1]` da semente mais proxima, pela REGRA DAS CELULAS
// MINIMAS, com os limiares DERIVADOS e nao digitados -- e tres diferencas que a
// terceira dimensao obriga:
//
// 1. A TABELA DE NIVEIS SAI POR DILATACAO, nao por busca exaustiva.  Em 2D a
//    forca bruta e' O(ncel x nsem) e passa; em 3D, 250 mil celulas contra
//    milhares de sementes dao ~1e9 operacoes POR REMALHAMENTO.  Aqui cada
//    semente carimba a sua vizinhanca de raio ceil(thr[0]/h) -- medido, 5
//    celulas com dois niveis, caixa de 1331 --, o que da' O(nsem x R^3).
//
//    A dilatacao e' EXATA, nao uma aproximacao: o nivel e' monotono decrescente
//    na distancia, logo o maximo sobre as sementes e' o nivel da semente mais
//    proxima.  `amr_tabela_exaustiva` existe ao lado para o arreio exigir
//    igualdade bit a bit entre as duas -- mesma disciplina da grade espacial.
//
// 2. A VORTICIDADE DEIXA DE SER ESCALAR.  Em 2D basta omega_z; em 3D e'
//    |rot u| com as tres componentes.  Sao 12 avaliacoes de faceta por celula
//    contra 4 em 2D.
//
// 3. O FORMATO DO .amr e' generico em DIM -- conferido no parser
//    higio_read_amr_info --, mas o numero de campos por retalho dobra:
//    initcell[3] e patchsize[3].  INDICE BASE 1, como em 2D: com base 0 a
//    contagem de celulas finas sai CERTA e as posicoes ERRADAS, e esse e' o modo
//    de falha que nao aparece em nenhum total.
// ---------------------------------------------------------------------------

//! Escreve o .amr multinivel.  `ns` nulo no ARRANQUE: sem campo de velocidade
//! nao ha' esteira a medir, e a malha inicial sai so' da geometria.
extern "C" void ft3_escreve_amr_criterio(higflow_solver *ns, ft3_superficie *sup,
                                         real lx, real ly, real lz,
                                         int nx, int ny, int nz,
                                         int niveis, int cel_min,
                                         const char *caminho)
{
    const real hx = lx / nx, hy = ly / ny, hz = lz / nz;
    real h_base = hx;
    if (hy < h_base) h_base = hy;
    if (hz < h_base) h_base = hz;

    real *thr = (real *) malloc((size_t) niveis * sizeof(real));
    amr_limiares(h_base, niveis, cel_min, thr);

    // SEMENTES: centros das celulas base cortadas pela superficie.  A grade
    // espacial e' obrigatoria aqui -- sem ela cada consulta e' O(nt), e sao
    // nx*ny*nz consultas.  A grade e' criada uma vez e morre no fim; a
    // superficie nao muda sob ela dentro desta funcao.
    ft3_grade *g = ft3_grade_cria(sup);
    const long ncel = (long) nx * ny * nz;
    real *sem = (real *) malloc((size_t) ncel * 3 * sizeof(real));
    long nsem = 0;
    // So' a envoltoria da superficie pode ter celula cortada: varrer o dominio
    // inteiro desperdica o grosso das consultas.  A folga de duas celulas cobre
    // o arredondamento dos indices.
    real elo[DIM], ehi[DIM];
    {
        const int nv = ft3_num_vertices(sup);
        const Point *pv = ft3_posicoes(sup);
        for (int d = 0; d < DIM; d++) { elo[d] = 1e300; ehi[d] = -1e300; }
        for (int q = 0; q < nv; q++)
            for (int d = 0; d < DIM; d++) {
                if (pv[q][d] < elo[d]) elo[d] = pv[q][d];
                if (pv[q][d] > ehi[d]) ehi[d] = pv[q][d];
            }
    }
    int i0 = (int) (elo[0]/hx) - 2, i1 = (int) (ehi[0]/hx) + 2;
    int j0 = (int) (elo[1]/hy) - 2, j1 = (int) (ehi[1]/hy) + 2;
    int k0 = (int) (elo[2]/hz) - 2, k1 = (int) (ehi[2]/hz) + 2;
    if (i0 < 0) i0 = 0;
    if (i1 >= nx) i1 = nx - 1;
    if (j0 < 0) j0 = 0;
    if (j1 >= ny) j1 = ny - 1;
    if (k0 < 0) k0 = 0;
    if (k1 >= nz) k1 = nz - 1;
    for (int kk = k0; kk <= k1; kk++)
      for (int jj = j0; jj <= j1; jj++)
        for (int ii = i0; ii <= i1; ii++) {
            real lo[DIM], hi[DIM];
            lo[0] = ii*hx; hi[0] = (ii+1)*hx;
            lo[1] = jj*hy; hi[1] = (jj+1)*hy;
            lo[2] = kk*hz; hi[2] = (kk+1)*hz;
            real f = ft3_fracao_na_caixa_g(sup, g, lo, hi);
            if (f > 0.001 && f < 0.999) {
                sem[3*nsem+0] = lo[0] + 0.5*hx;
                sem[3*nsem+1] = lo[1] + 0.5*hy;
                sem[3*nsem+2] = lo[2] + 0.5*hz;
                nsem++;
            }
        }
    ft3_grade_destroi(g);

    signed char *tab = (signed char *) malloc((size_t) ncel);
    amr_tabela(sem, nsem, hx, hy, hz, nx, ny, nz, thr, niveis, tab);

    // -----------------------------------------------------------------
    // A ESTEIRA, pelo MAXIMO com a banda geometrica.  Em 3D a vorticidade e'
    // |rot u| com as tres componentes, nao o escalar omega_z do 2D.
    //
    // HISTERESE (FT3_VORT_HISTERESE): quem JA' estava refinado so' engrossa se
    // |omega| cair abaixo da METADE do limiar.  Em 2D, sem isso, mediu-se
    // FLAPPING -- celulas cruzando o limiar seco a cada ciclo e a malha
    // derivando.  Sair da malha fina tem de ser mais dificil que entrar.
    // -----------------------------------------------------------------
    static signed char *tab_ant = NULL;
    static long ncel_ant = 0;
    {
        const char *e1 = getenv("FT3_VORT1"), *e2 = getenv("FT3_VORT2");
        const real v1 = (e1 != NULL) ? atof(e1) : 0.0;
        const real v2 = (e2 != NULL) ? atof(e2) : 0.0;
        const int  hist = (getenv("FT3_VORT_HISTERESE") != NULL);
        if (v1 > 0.0 && ns != NULL) {
            sim_domain *sdp = psd_get_local_domain(ns->psdp);
            sim_facet_domain *sfdu[DIM];
            for (int d = 0; d < DIM; d++)
                sfdu[d] = psfd_get_local_domain(ns->psfdu[d]);
            const hig_mesh_snapshot *hms = sd_get_snapshot(sdp);
            for (int clid = 0; clid < hms->n; clid++) {
                Point cc, cd;
                hms_center(hms, clid, cc);
                hms_delta(hms, clid, cd);
                // d u_a / d x_b por diferenca centrada nas facetas
                real du[DIM][DIM];
                for (int a = 0; a < DIM; a++) {
                    for (int b = 0; b < DIM; b++) {
                        Point pp;
                        POINT_ASSIGN(pp, cc); pp[b] = cc[b] + cd[b];
                        real up = compute_facet_value_at_point(sfdu[a], cc, pp, 1.0,
                                                               ns->dpu[a], ns->stn);
                        POINT_ASSIGN(pp, cc); pp[b] = cc[b] - cd[b];
                        real um = compute_facet_value_at_point(sfdu[a], cc, pp, 1.0,
                                                               ns->dpu[a], ns->stn);
                        du[a][b] = (up - um) / (2.0 * cd[b]);
                    }
                }
                // |rot u| pela funcao COMPARTILHADA com o arreio, para o oraculo
                // testar o codigo que RODA: e' aqui que sinal trocado ou eixo
                // permutado passa sem sintoma -- a malha sai refinada no lugar
                // quase certo e nada denuncia.
                const real om = esteira_omega_grad(du);
                int ii = (int)(cc[0]/hx), jj = (int)(cc[1]/hy), kk = (int)(cc[2]/hz);
                if (ii < 0) ii = 0;  if (ii >= nx) ii = nx - 1;
                if (jj < 0) jj = 0;  if (jj >= ny) jj = ny - 1;
                if (kk < 0) kk = 0;  if (kk >= nz) kk = nz - 1;
                const long q = ii + (long) nx * (jj + (long) ny * kk);
                signed char alvo = 0;
                if (v2 > 0.0 && om > v2) alvo = (signed char) niveis;
                else if (om > v1)        alvo = 1;
                if (hist && tab_ant != NULL && ncel_ant == ncel && tab_ant[q] > alvo) {
                    if (om > 0.5 * v1) alvo = tab_ant[q];
                }
                if (alvo > tab[q]) tab[q] = alvo;   // MAXIMO com a banda
            }
        }
        if (ncel_ant != ncel) { free(tab_ant); tab_ant = NULL; }
        if (tab_ant == NULL) tab_ant = (signed char *) malloc((size_t) ncel);
        memcpy(tab_ant, tab, (size_t) ncel);
        ncel_ant = ncel;
    }

    // -----------------------------------------------------------------
    // Escrita.  SO' O RANK 0 ESCREVE, E TODOS ESPERAM: as sementes saem da
    // SUPERFICIE, que e' replicada e portanto completa em todo rank, logo nao
    // ha' reuniao a fazer -- mas sem a guarda todos truncam o mesmo arquivo, e
    // sem a barreira alguem o rele' antes de estar pronto.
    // -----------------------------------------------------------------
    int meurank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &meurank);
    if (meurank != 0) {
        MPI_Barrier(MPI_COMM_WORLD);
        free(tab); free(sem); free(thr);
        return;
    }
    if (amr_escreve_arquivo(caminho, lx, ly, lz, nx, ny, nz,
                                niveis, tab) != 0) {
        MPI_Barrier(MPI_COMM_WORLD);
        free(tab); free(sem); free(thr);
        return;
    }
    fprintf(stderr, "FT3 amr: %ld sementes de %ld celulas base, niveis=%d, "
                    "limiares=", nsem, ncel, niveis);
    for (int l = 0; l < niveis; l++) fprintf(stderr, " %.5f", (double) thr[l]);
    fprintf(stderr, "  (>=%d celulas finas por lado)\n", cel_min);
    free(tab); free(sem); free(thr);
    MPI_Barrier(MPI_COMM_WORLD);
}
