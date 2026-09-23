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

static real campo_centro(const Point x)
{
    real v = 1.0 + x[0] + 10.0 * x[1];
#if DIM == 3
    v += 100.0 * x[2];
#endif
    return v;
}
static real campo_faceta(int dim, const Point x)
{
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

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    PetscInitialize(&argc, &argv, NULL, NULL);
    int rank, ntasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    Malha A = monta(rank, 0);
    { Semeia s = { A.dpc, -1 }; percorre_centro(&A, _semeia_c, &s); }
    for (int d = 0; d < DIM; d++) { Semeia s = { A.dpu[d], d }; percorre_faceta(&A, d, _semeia_f, &s); }

    rem_colheita *hc = rem_colhe_centro(A.sd, A.dpc);
    long colisoes = rem_colisoes();
    rem_colheita *hf[DIM];
    for (int d = 0; d < DIM; d++) { hf[d] = rem_colhe_faceta(A.sfd[d], A.dpu[d]); colisoes += rem_colisoes(); }

    int falhou = 0;
    for (int parte = 0; parte < 2; parte++) {
        const int extra = parte;      // 0: mesma malha.  1: malha mudada.
        Malha D = monta(rank, extra);

        char *ac = (char *) calloc((size_t) D.dpc->pdata->total_count, 1);
        long perdidas = rem_planta_centro(hc, D.sd, D.dpc, ac);
        long fora = rem_vieram_de_outro_rank();
        Confere Cc = { D.dpc, ac, -1, 0, 0, 0 };
        percorre_centro(&D, _confere, &Cc);

        long erradas = Cc.erradas, buracos = Cc.buracos_fora, novas = Cc.novas;
        for (int d = 0; d < DIM; d++) {
            char *af = (char *) calloc((size_t) D.dpu[d]->pdata->total_count, 1);
            perdidas += rem_planta_faceta(hf[d], D.sfd[d], D.dpu[d], af);
            fora += rem_vieram_de_outro_rank();
            Confere Cf = { D.dpu[d], af, d, 0, 0, 0 };
            percorre_faceta(&D, d, _confere, &Cf);
            erradas += Cf.erradas; buracos += Cf.buracos_fora; novas += Cf.novas;
            free(af);
        }
        free(ac);

        long g[6] = { colisoes, perdidas, fora, erradas, buracos, novas }, gs[6];
        MPI_Allreduce(g, gs, 6, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

        if (rank == 0)
            printf("\nPARTE %d -- malha %s, np = %d\n", parte + 1,
                   extra ? "MUDADA (B refina x>0,75 e y>0,5)" : "IGUAL", ntasks);

        checa(rank, "1. sem colisao de chave na colheita", gs[0] == 0);
        checa(rank, "2. todo valor achado voltou bit a bit igual", gs[3] == 0);
        checa(rank, "3. nenhum buraco fora da regiao que mudou", gs[4] == 0);
        if (!extra) {
            checa(rank, "4. malha igual: nada ficou sem valor", gs[1] == 0);
            if (rank == 0)
                printf("     (particao igual: %ld valores mudaram de rank --\n"
                       "      por isso esta parte NAO afirma a redistribuicao)\n", gs[2]);
            falhou |= (gs[0] != 0) || (gs[3] != 0) || (gs[4] != 0) || (gs[1] != 0);
        } else {
            checa(rank, "4. malha mudada: houve celula nova, e foi reportada", gs[5] > 0);
            if (ntasks > 1)
                checa(rank, "5. houve valor vindo de OUTRO rank (nao vazio)", gs[2] > 0);
            if (rank == 0)
                printf("     (%ld valores mudaram de rank, %ld posicoes novas)\n",
                       gs[2], gs[5]);
            falhou |= (gs[0] != 0) || (gs[3] != 0) || (gs[4] != 0) || (gs[5] == 0) ||
                      (ntasks > 1 && gs[2] == 0);
        }
    }

    rem_destroi(hc);
    for (int d = 0; d < DIM; d++) rem_destroi(hf[d]);
    PetscFinalize();
    MPI_Finalize();
    return falhou;
}
