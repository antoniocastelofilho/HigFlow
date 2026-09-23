// Oraculo para a transferencia de solucao entre malhas.
//
// Duas partes, e a primeira sozinha nao bastaria:
//
//   PARTE 1, malha IGUAL.  Colher de A e plantar em A' identica.  Todo valor
//   tem de voltar bit a bit igual -- a transferencia e' copia, nao faz conta.
//   Mas: MEDIDO, o particionador e' determinista no conjunto de celulas, entao
//   A' sai com a MESMA particao e nenhum valor troca de rank.  Uma
//   transferencia completamente quebrada passa nesta parte, porque cada rank
//   acha em casa tudo o que pede.  Tentei duas alavancas para forcar particao
//   diferente com a mesma malha -- peso de grupo (1 contra 9) e decomposicao em
//   8 contra 4 arvores -- e as DUAS deram particao identica, conferido pela
//   contagem de celulas proprias por rank.
//
//   PARTE 2, malha MUDADA.  E' a unica que move dado, e e' o caso real.  B
//   refina uma regiao que A tinha grossa.  Entao:
//     - as posicoes que existem nas duas tem de voltar exatas, e muitas mudam
//       de rank;
//     - as posicoes NOVAS tem de ser reportadas como nao achadas, e tem de
//       estar TODAS dentro da regiao que mudou.  Buraco fora dela seria perda
//       silenciosa de campo.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <petsc.h>

#include "higtree.h"
#include "higtree-iterator.h"
#include "pdomain.h"
#include "domain.h"
#include "lbal.h"
#include "hig-flow-remalha.h"

#define NC       16
#define NTIRAS    4
#define FRANJA    5
#define XCORTE  0.75      // a regiao que B refina e A nao: x > XCORTE e y >= YCORTE
#define YCORTE  0.50
#define FOLGA   1.0e-9

static int g_constante = 0;      // 1: campo uniforme, para o oraculo de facetas

static real campo_centro(const Point x)
{
    if (g_constante) return 7.25;
    real v = 1.0 + x[0] + 10.0 * x[1];
#if DIM == 3
    v += 100.0 * x[2];
#endif
    return v;
}
static real campo_faceta(int dim, const Point x)
{
    if (g_constante) return 7.25 + 100.0 * dim;   // constante, distinta por direcao
    return 1000.0 * (dim + 1) + campo_centro(x);
}

// Fora da regiao que mudou, nao achar e' buraco.
static int na_regiao_nova(const Point x)
{
    return (x[0] > XCORTE - FOLGA) && (x[1] > YCORTE - FOLGA);
}

typedef struct {
    partition_graph      *pg;
    sim_domain           *sd;
    psim_domain          *psd;
    sim_facet_domain     *sfd[DIM];
    psim_facet_domain    *psfd[DIM];
    distributed_property *dpc;
    distributed_property *dpu[DIM];
} Malha;

static Malha monta(int rank, int extra)
{
    Malha M;
    M.pg = pg_create(MPI_COMM_WORLD);
    pg_set_fringe_size(M.pg, FRANJA);
    load_balancer *lb = lb_create(MPI_COMM_WORLD, 1);

    if (rank == 0) {
        for (int t = 0; t < NTIRAS; t++) {
            Point lo, hi;
            POINT_ASSIGN_SCALAR(lo, 0.0);
            POINT_ASSIGN_SCALAR(hi, 1.0);
            lo[0] = (real) t       / NTIRAS;
            hi[0] = (real) (t + 1) / NTIRAS;
            hig_cell *raiz = hig_create_root(lo, hi);
            int nc[DIM];
            for (int d = 0; d < DIM; d++) nc[d] = NC;
            nc[0] = NC / NTIRAS;
            hig_refine_uniform(raiz, nc);

            int dois[DIM];
            for (int d = 0; d < DIM; d++) dois[d] = 2;
            higcit_celliterator *it;
            for (it = higcit_create_all_leaves(raiz); !higcit_isfinished(it);
                 higcit_nextcell(it)) {
                hig_cell *c = higcit_getcell(it);
                Point cc;
                hig_get_center(c, cc);
                const int base  = (cc[1] < YCORTE);
                const int nova  = extra && (cc[0] > XCORTE) && (cc[1] > YCORTE);
                if (base || nova) hig_refine_uniform(c, dois);
            }
            higcit_destroy(it);
            lb_add_input_tree(lb, raiz, true, 0);
        }
    }
    lb_calc_partition(lb, M.pg);

    M.sd = sd_create(NULL);
    sd_set_interpolator_order(M.sd, 2);
    for (unsigned i = 0; i < lb_get_num_local_trees(lb); ++i)
        sd_add_higtree(M.sd, lb_get_local_tree(lb, i, NULL));
    lb_destroy(lb);
    M.psd = psd_create(M.sd, M.pg);
    psd_synced_mapper(M.psd);
    M.dpc = psd_create_property(M.psd);

    for (int dim = 0; dim < DIM; dim++) {
        M.sfd[dim] = sfd_create(NULL, dim);
        sfd_set_interpolator_order(M.sfd[dim], 2);
        sfd_copy_higtrees_from_center_domain(M.sfd[dim], M.sd);
        sfd_adjust_facet_ids(M.sfd[dim]);
        M.psfd[dim] = psfd_create(M.sfd[dim], M.psd);
        psfd_compute_sfbi(M.psfd[dim]);
        psfd_synced_mapper(M.psfd[dim]);
        M.dpu[dim] = psfd_create_property(M.psfd[dim]);
    }
    return M;
}

// Percorre as entidades PROPRIAS chamando `f(lid, posicao, ctx)`.
typedef void (*Visita)(int lid, const Point x, void *ctx);

static void percorre_centro(Malha *M, Visita f, void *ctx)
{
    mp_mapper *m = sd_get_domain_mapper(M->sd);
    const int nloc = M->dpc->pdata->local_count;
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(M->sd); !higcit_isfinished(it);
         higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        const int lid = mp_lookup(m, hig_get_cid(c));
        if (lid < 0 || lid >= nloc) continue;
        Point cc;
        hig_get_center(c, cc);
        f(lid, cc, ctx);
    }
    higcit_destroy(it);
}

static void percorre_faceta(Malha *M, int dim, Visita f, void *ctx)
{
    mp_mapper *m = M->sfd[dim]->fm;
    const int nloc = M->dpu[dim]->pdata->local_count;
    for (int k = 0; k < sd_get_num_higtrees(M->sd); k++) {
        hig_cell *root = sd_get_higtree(M->sd, k);
        Point blo, bhi;
        POINT_ASSIGN_SCALAR(blo, -1.0); POINT_ASSIGN_SCALAR(bhi, 2.0);
        higfit_facetiterator *fit;
        for (fit = higfit_create_bounding_box_facets(root,
                    M->sfd[dim]->dimofinterest, blo, bhi);
             !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            hig_facet *ff = higfit_getfacet(fit);
            const int lid = mp_lookup(m, hig_get_fid(ff));
            if (lid < 0 || lid >= nloc) continue;
            Point fc;
            hig_get_facet_center(ff, fc);
            f(lid, fc, ctx);
        }
        higfit_destroy(fit);
    }
}

typedef struct { distributed_property *dp; int dim; } Semeia;
static void _semeia_c(int lid, const Point x, void *v)
{ dp_set_value(((Semeia *) v)->dp, lid, campo_centro(x)); }
static void _semeia_f(int lid, const Point x, void *v)
{ Semeia *s = (Semeia *) v; dp_set_value(s->dp, lid, campo_faceta(s->dim, x)); }

typedef struct {
    distributed_property *dp;
    const char *achado;
    int   dim;            // -1 = centro
    long  erradas, buracos_fora, novas;
} Confere;

static void _confere(int lid, const Point x, void *v)
{
    Confere *C = (Confere *) v;
    const real esperado = (C->dim < 0) ? campo_centro(x) : campo_faceta(C->dim, x);
    if (C->achado != NULL && !C->achado[lid]) {
        C->novas++;
        if (!na_regiao_nova(x)) C->buracos_fora++;   // buraco fora da regiao que mudou
        return;
    }
    if (dp_get_value(C->dp, lid) != esperado) C->erradas++;
}

static void checa(int rank, const char *o_que, int ok)
{
    int todos = 0;
    MPI_Allreduce(&ok, &todos, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (rank == 0) printf("  %-54s %s\n", o_que, todos ? "ok" : "FALHOU");
}

// A integral do campo: soma valor * volume nas celulas PROPRIAS.  E' a
// quantidade que a interpolacao conservativa tem de preservar, qualquer que
// seja o campo -- e' isso que "conservativa" quer dizer.
typedef struct { distributed_property *dp; double soma; } Integral;

static void _integra_c(int lid, const Point x, void *v)
{
    (void) x;
    Integral *I = (Integral *) v;
    I->soma += dp_get_value(I->dp, lid);     // peso entra fora, ver abaixo
}

// Para celula e faceta o peso e' o volume (ou a area) da entidade.  Como a
// malha e' dyadica e o valor e' constante por entidade, e' preciso o tamanho --
// entao a soma e' feita percorrendo a malha diretamente, e nao pelo visitante.
static double integral_centro(Malha *M)
{
    mp_mapper *m = sd_get_domain_mapper(M->sd);
    const int nloc = M->dpc->pdata->local_count;
    double s = 0.0;
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(M->sd); !higcit_isfinished(it);
         higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        const int lid = mp_lookup(m, hig_get_cid(c));
        if (lid < 0 || lid >= nloc) continue;
        Point cl, ch;
        hig_get_lowpoint(c, cl);
        hig_get_highpoint(c, ch);
        double vol = 1.0;
        for (int d = 0; d < DIM; d++) vol *= (ch[d] - cl[d]);
        s += dp_get_value(M->dpc, lid) * vol;
    }
    higcit_destroy(it);
    return s;
}

static double integral_faceta(Malha *M, int dim)
{
    mp_mapper *m = M->sfd[dim]->fm;
    const int nloc = M->dpu[dim]->pdata->local_count;
    double s = 0.0;
    for (int k = 0; k < sd_get_num_higtrees(M->sd); k++) {
        hig_cell *root = sd_get_higtree(M->sd, k);
        Point blo, bhi;
        POINT_ASSIGN_SCALAR(blo, -1.0); POINT_ASSIGN_SCALAR(bhi, 2.0);
        higfit_facetiterator *fit;
        for (fit = higfit_create_bounding_box_facets(root,
                    M->sfd[dim]->dimofinterest, blo, bhi);
             !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            hig_facet *f = higfit_getfacet(fit);
            const int lid = mp_lookup(m, hig_get_fid(f));
            if (lid < 0 || lid >= nloc) continue;
            hig_cell *c = hig_get_facet_cell(f);
            Point cl, ch;
            hig_get_lowpoint(c, cl);
            hig_get_highpoint(c, ch);
            double area = 1.0;               // area da faceta: o produto das
            for (int d = 0; d < DIM; d++)    // extensoes MENOS a normal
                if (d != dim) area *= (ch[d] - cl[d]);
            s += dp_get_value(M->dpu[dim], lid) * area;
        }
        higfit_destroy(fit);
    }
    return s;
}

static double global(double x)
{
    double g = 0.0;
    MPI_Allreduce(&x, &g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return g;
}

// Uma transferencia completa: colhe de `de`, planta em `para`, interpola o que
// faltou, e afirma as tres coisas que importam.
static int transfere(int rank, int ntasks, const char *nome, Malha *de, Malha *para)
{
    // Semeia a origem e mede a integral dela.
    { Semeia s = { de->dpc, -1 }; percorre_centro(de, _semeia_c, &s); }
    for (int d = 0; d < DIM; d++) { Semeia s = { de->dpu[d], d }; percorre_faceta(de, d, _semeia_f, &s); }

    double Ic_antes = global(integral_centro(de));
    double If_antes[DIM];
    for (int d = 0; d < DIM; d++) If_antes[d] = global(integral_faceta(de, d));

    // SONDA: quantas facetas PROPRIAS existem no plano x = 1,0 (a fronteira do
    // dominio) na origem?  Se forem zero, a colheita nao as publica, e o pai do
    // plano do meio adjacente a' borda nunca sera' achado.
    if (getenv("REMALHA_ONDE") != NULL) {
        mp_mapper *mm = de->sfd[0]->fm;
        const int nl = de->dpu[0]->pdata->local_count;
        long na_borda = 0, total = 0;
        for (int k = 0; k < sd_get_num_higtrees(de->sd); k++) {
            hig_cell *root = sd_get_higtree(de->sd, k);
            Point blo, bhi;
            POINT_ASSIGN_SCALAR(blo, -1.0); POINT_ASSIGN_SCALAR(bhi, 2.0);
            higfit_facetiterator *fit;
            for (fit = higfit_create_bounding_box_facets(root,
                        de->sfd[0]->dimofinterest, blo, bhi);
                 !higfit_isfinished(fit); higfit_nextfacet(fit)) {
                hig_facet *ff = higfit_getfacet(fit);
                const int lid = mp_lookup(mm, hig_get_fid(ff));
                if (lid < 0 || lid >= nl) continue;
                Point fc;
                hig_get_facet_center(ff, fc);
                total++;
                {
                    const char *xs = getenv("REMALHA_PLANO");
                    const double xp = xs ? atof(xs) : 1.0;
                    if (fabs(fc[0] - xp) < 1e-9) na_borda++;
                }
            }
            higfit_destroy(fit);
        }
        printf("     [sonda] origem dim 0: %ld facetas proprias, %ld no plano pedido\n",
               total, na_borda);
        fflush(stdout);
    }

    rem_colheita *hc = rem_colhe_centro(de->sd, de->dpc);
    rem_colheita *hf[DIM];
    for (int d = 0; d < DIM; d++) hf[d] = rem_colhe_faceta(de->sfd[d], de->dpu[d]);

    // Zera o destino, para que "ficou com valor" signifique "foi preenchido".
    for (int i = 0; i < para->dpc->pdata->total_count; i++) dp_set_value(para->dpc, i, 0.0);
    for (int d = 0; d < DIM; d++)
        for (int i = 0; i < para->dpu[d]->pdata->total_count; i++) dp_set_value(para->dpu[d], i, 0.0);

    char *ac = (char *) calloc((size_t) para->dpc->pdata->total_count, 1);
    rem_planta_centro(hc, para->sd, para->dpc, ac);
    long exatos_c = 0;
    { Confere C = { para->dpc, ac, -1, 0, 0, 0 }; percorre_centro(para, _confere, &C);
      exatos_c = C.erradas; }

    long faltam = rem_interpola_centro(hc, para->sd, para->dpc, ac);
    for (int d = 0; d < DIM; d++) {
        char *af = (char *) calloc((size_t) para->dpu[d]->pdata->total_count, 1);
        rem_planta_faceta(hf[d], para->sfd[d], para->dpu[d], af);
        faltam += rem_interpola_faceta(hf[d], para->sfd[d], para->dpu[d], af);
        free(af);
    }
    free(ac);

    double Ic_dep = global(integral_centro(para));
    double If_dep[DIM];
    for (int d = 0; d < DIM; d++) If_dep[d] = global(integral_faceta(para, d));

    long g[2] = { faltam, exatos_c }, gs[2];
    MPI_Allreduce(g, gs, 2, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    long refino = 0, engross = 0, gr[2], gg[2];
    rem_interpolou(&refino, &engross);
    gr[0] = refino; gr[1] = engross;
    MPI_Allreduce(gr, gg, 2, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0) printf("\n%s  (np = %d)\n", nome, ntasks);

    const double rel_c = fabs(Ic_dep - Ic_antes) / (fabs(Ic_antes) + 1e-300);
    double rel_f = 0.0;
    for (int d = 0; d < DIM; d++) {
        const double r = fabs(If_dep[d] - If_antes[d]) / (fabs(If_antes[d]) + 1e-300);
        if (r > rel_f) rel_f = r;
    }

    // Conferir o campo em TODA entidade do destino -- inclusive as que a
    // interpolacao preencheu.  Com campo constante, as tres rotas tem de
    // devolver a constante exata: filha recebe o da mae, mae recebe a media das
    // filhas, e a faceta do meio recebe a media das duas paralelas.
    long fora_c = 0, fora_f = 0;
    if (g_constante) {
        Confere C = { para->dpc, NULL, -1, 0, 0, 0 };
        percorre_centro(para, _confere, &C);
        fora_c = C.erradas;
        for (int d = 0; d < DIM; d++) {
            Confere F = { para->dpu[d], NULL, d, 0, 0, 0 };
            percorre_faceta(para, d, _confere, &F);
            fora_f += F.erradas;
        }
        // ONDE falham.  "Falha" sem sitio nao e' diagnostico.
        if (getenv("REMALHA_ONDE") != NULL) {
            for (int d = 0; d < DIM; d++) {
                mp_mapper *mm = para->sfd[d]->fm;
                const int nl = para->dpu[d]->pdata->local_count;
                int mostrados = 0;
                for (int k = 0; k < sd_get_num_higtrees(para->sd) && mostrados < 6; k++) {
                    hig_cell *root = sd_get_higtree(para->sd, k);
                    Point blo, bhi;
                    POINT_ASSIGN_SCALAR(blo, -1.0); POINT_ASSIGN_SCALAR(bhi, 2.0);
                    higfit_facetiterator *fit;
                    for (fit = higfit_create_bounding_box_facets(root,
                                para->sfd[d]->dimofinterest, blo, bhi);
                         !higfit_isfinished(fit) && mostrados < 6; higfit_nextfacet(fit)) {
                        hig_facet *ff = higfit_getfacet(fit);
                        const int lid = mp_lookup(mm, hig_get_fid(ff));
                        if (lid < 0 || lid >= nl) continue;
                        Point fc;
                        hig_get_facet_center(ff, fc);
                        const real esp = campo_faceta(d, fc);
                        const real got = dp_get_value(para->dpu[d], lid);
                        if (got == esp) continue;
                        hig_cell *cc = hig_get_facet_cell(ff);
                        Point cl, ch;
                        hig_get_lowpoint(cc, cl); hig_get_highpoint(cc, ch);
                        printf("     FALHA dim=%d faceta (%.4f,%.4f) celula %.4fx%.4f "
                               "valor %.3f esperado %.3f\n",
                               d, fc[0], fc[1], ch[0]-cl[0], ch[1]-cl[1], got, esp);
                        mostrados++;
                    }
                    higfit_destroy(fit);
                }
            }
            fflush(stdout);
        }
    }
    long gf[2] = { fora_c, fora_f }, gfs[2];
    MPI_Allreduce(gf, gfs, 2, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    checa(rank, "1. nada ficou sem valor", gs[0] == 0);
    checa(rank, "2. o que veio exato esta' bit a bit igual", gs[1] == 0);
    checa(rank, "3. integral de centro preservada", rel_c < 1e-12);
    // A INTEGRAL DE FACETA NAO E' INVARIANTE, e exigi-la era erro meu.  Refinar
    // cria faceta no plano do MEIO, que no nivel grosso era interior: ela soma
    // fluxo que antes nao existia.  Somar u*A sobre TODAS as facetas portanto
    // muda, e deve mudar.  O que a transferencia tem de garantir e' que o campo
    // constante saia constante -- e isso vale nas tres rotas.
    if (g_constante) {
        checa(rank, "4. campo constante reproduzido em toda celula", gfs[0] == 0);
        checa(rank, "5. campo constante reproduzido em toda faceta", gfs[1] == 0);
    }
    // Sem isto, uma transferencia que nunca interpolasse passaria em 1-4.
    checa(rank, "6. a interpolacao foi exercitada", gg[0] + gg[1] > 0);
    if (rank == 0)
        printf("     (%ld por refino, %ld por engrossamento; erro relativo na\n"
               "      integral: centro %.2e, faceta %.2e)\n",
               gg[0], gg[1], rel_c, rel_f);

    for (int d = 0; d < DIM; d++) rem_destroi(hf[d]);
    rem_destroi(hc);
    (void) rel_f;
    return (gs[0] != 0) || (gs[1] != 0) || !(rel_c < 1e-12) ||
           (gg[0] + gg[1] == 0) || (gfs[0] != 0) || (gfs[1] != 0);
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    PetscInitialize(&argc, &argv, NULL, NULL);
    int rank, ntasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    Malha A = monta(rank, 0);      // sem a regiao extra
    Malha B = monta(rank, 1);      // com a regiao extra refinada

    int falhou = 0;
    // As DUAS direcoes.  So' refinar deixaria o caminho de engrossamento sem
    // teste nenhum, e ele e' metade do codigo.
    for (int passo = 0; passo < 2; passo++) {
        g_constante = passo;
        if (rank == 0)
            printf("\n========== campo %s ==========\n",
                   g_constante ? "CONSTANTE (afirma toda entidade, inclusive as interpoladas)"
                               : "analitico (afirma exatidao do que casou, e a integral de centro)");
        falhou |= transfere(rank, ntasks, "A -> B  (refina)", &A, &B);
        falhou |= transfere(rank, ntasks, "B -> A  (engrossa)", &B, &A);
    }

    if (rank == 0)
        printf("\n%s\n", falhou ? "FALHOU" : "todas as afirmacoes passaram");

    PetscFinalize();
    MPI_Finalize();
    return falhou;
}
