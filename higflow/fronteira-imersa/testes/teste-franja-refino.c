// Afirma o VALOR na franja, nao a forma dela.
//
// A guarda da HiGTree (pdomain.c) confere que a troca de franja e' SIMETRICA:
// que o remetente manda tantos elementos quantos o destinatario dimensionou.
// Isso pega o defeito do example2d_DynamicMeshAdapt, onde as contagens diferem.
// Nao pega nada sobre o CONTEUDO: uma troca pode ser perfeitamente simetrica e
// ainda assim depositar o valor errado em cada faceta.
//
// Este teste poe um campo ANALITICO nas facetas proprias, envenena as de franja
// com um valor impossivel, sincroniza, e confere faceta por faceta contra a
// funcao avaliada no centro dela.  Erro de valor aparece como diferenca; franja
// nunca preenchida aparece como veneno sobrevivente.  Os dois sao reportados
// com a POSICAO, porque "diverge" nao e' diagnostico e "diverge aqui" e'.
//
// A malha tem a interface de refino no plano x = 0,5 -- uma linha inteira
// atravessando o dominio, que qualquer particao nao trivial cruza.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "higtree.h"
#include "higtree-io.h"
#include "higtree-iterator.h"
#include "pdomain.h"
#include "domain.h"
#include "lbal.h"
#include <petsc.h>

#define NC      16
#define FRANJA  5
#define VENENO  (-987654.0)

// Campo analitico: distinto em cada faceta, e nao simetrico em nenhum eixo, para
// que uma troca de indices errada nao possa passar por coincidencia.
static real campo(const Point x)
{
    real v = 1.0 + x[0] + 10.0 * x[1];
#if DIM == 3
    v += 100.0 * x[2];
#endif
    return v;
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    PetscInitialize(&argc, &argv, NULL, NULL);
    int rank, ntasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    partition_graph *pg = pg_create(MPI_COMM_WORLD);
    pg_set_fringe_size(pg, FRANJA);
    load_balancer *lb = lb_create(MPI_COMM_WORLD, 1);
    if (rank == 0) {
        Point lo, hi;
        POINT_ASSIGN_SCALAR(lo, 0.0);
        POINT_ASSIGN_SCALAR(hi, 1.0);
        hig_cell *raiz = hig_create_root(lo, hi);
        int nc[DIM];
        for (int d = 0; d < DIM; d++) nc[d] = NC;
        hig_refine_uniform(raiz, nc);

        // Meio dominio refinado: a interface fica no plano x = 0,5, inteira.
        int dois[DIM];
        for (int d = 0; d < DIM; d++) dois[d] = 2;
        higcit_celliterator *it;
        for (it = higcit_create_all_leaves(raiz); !higcit_isfinished(it);
             higcit_nextcell(it)) {
            hig_cell *c = higcit_getcell(it);
            Point cc;
            hig_get_center(c, cc);
            if (cc[0] < 0.5) hig_refine_uniform(c, dois);
        }
        higcit_destroy(it);
        lb_add_input_tree(lb, raiz, true, 0);
    }
    lb_calc_partition(lb, pg);

    sim_domain *sd = sd_create(NULL);
    sd_set_interpolator_order(sd, 2);
    for (unsigned i = 0; i < lb_get_num_local_trees(lb); ++i)
        sd_add_higtree(sd, lb_get_local_tree(lb, i, NULL));
    lb_destroy(lb);
    psim_domain *psd = psd_create(sd, pg);
    psd_synced_mapper(psd);

    int falhou = 0;
    for (int dim = 0; dim < DIM; dim++) {
        sim_facet_domain *sfd = sfd_create(NULL, dim);
        sfd_set_interpolator_order(sfd, 2);
        sfd_copy_higtrees_from_center_domain(sfd, sd);
        sfd_adjust_facet_ids(sfd);
        psim_facet_domain *psfd = psfd_create(sfd, psd);
        psfd_compute_sfbi(psfd);
        psfd_synced_mapper(psfd);
        distributed_property *dp = psfd_create_property(psfd);

        const int nloc = dp->pdata->local_count;
        const int ntot = dp->pdata->total_count;

        // Proprias recebem o campo; franja recebe veneno.  Quem preenche a
        // franja e' o dp_sync, e so' ele.
        for (int lid = 0; lid < ntot; lid++) dp_set_value(dp, lid, VENENO);
        {
            sim_domain *cd = sfd->cdom;
            for (int k = 0; k < sd_get_num_higtrees(cd); k++) {
                hig_cell *root = sd_get_higtree(cd, k);
                Point blo, bhi;
                POINT_ASSIGN_SCALAR(blo, -1.0);
                POINT_ASSIGN_SCALAR(bhi,  2.0);
                higfit_facetiterator *fit;
                for (fit = higfit_create_bounding_box_facets(root,
                            sfd->dimofinterest, blo, bhi);
                     !higfit_isfinished(fit); higfit_nextfacet(fit)) {
                    hig_facet *f = higfit_getfacet(fit);
                    int lid = mp_lookup(sfd->fm, hig_get_fid(f));
                    if (lid >= 0 && lid < nloc) {
                        Point fc;
                        hig_get_facet_center(f, fc);
                        dp_set_value(dp, lid, campo(fc));
                    }
                }
                higfit_destroy(fit);
            }
        }

        dp_sync(dp);

        // Agora a afirmacao: TODA faceta visivel vale o campo no centro dela.
        // Um teste que passa por nao exercitar o caso e' pior que nenhum.
        // Estes dois contadores afirmam que ha' franja, e que ha' franja NA
        // interface de refino -- sem eles, "passou" nao quer dizer nada.
        long n_franja = 0, n_franja_int = 0;
        long n_veneno = 0, n_errada = 0, n_conf = 0;
        real pior = 0.0;
        Point onde;
        POINT_ASSIGN_SCALAR(onde, 0.0);
        {
            sim_domain *cd = sfd->cdom;
            for (int k = 0; k < sd_get_num_higtrees(cd); k++) {
                hig_cell *root = sd_get_higtree(cd, k);
                Point blo, bhi;
                POINT_ASSIGN_SCALAR(blo, -1.0);
                POINT_ASSIGN_SCALAR(bhi,  2.0);
                higfit_facetiterator *fit;
                for (fit = higfit_create_bounding_box_facets(root,
                            sfd->dimofinterest, blo, bhi);
                     !higfit_isfinished(fit); higfit_nextfacet(fit)) {
                    hig_facet *f = higfit_getfacet(fit);
                    int lid = mp_lookup(sfd->fm, hig_get_fid(f));
                    if (lid < 0 || lid >= ntot) continue;
                    Point fc;
                    hig_get_facet_center(f, fc);
                    const real v = dp_get_value(dp, lid);
                    n_conf++;
                    if (lid >= nloc) {
                        n_franja++;
                        if (fabs(fc[0] - 0.5) < 1.5 / NC) n_franja_int++;
                    }
                    if (v == VENENO) {
                        n_veneno++;
                        if (n_veneno == 1) POINT_ASSIGN(onde, fc);
                        continue;
                    }
                    const real e = fabs(v - campo(fc));
                    // NaN nao e' maior que nada: testar por !(<=) e nao por (>).
                    if (!(e <= pior)) { pior = e; POINT_ASSIGN(onde, fc); }
                    if (!(e < 1e-12)) n_errada++;
                }
                higfit_destroy(fit);
            }
        }

        long g[5] = {0, 0, 0, 0, 0}, gs[5];
        g[0] = n_veneno; g[1] = n_errada; g[2] = n_conf;
        g[3] = n_franja; g[4] = n_franja_int;
        MPI_Allreduce(g, gs, 5, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
        // Vacuidade e' falha, nao sucesso.
        if (ntasks > 1 && gs[4] == 0) {
            if (rank == 0)
                printf("dim %d: VACUO -- nenhuma faceta de franja na interface "
                       "de refino; o teste nao exercita o caso\n", dim);
            falhou = 1;
        }

        if (gs[0] || gs[1]) {
            falhou = 1;
            // Cada rank que viu defeito diz ONDE, porque o sitio e' o dado.
            if (n_veneno || n_errada)
                printf("  dim %d rank %d: %ld veneno, %ld valor errado "
                       "(pior %.3e) em (%.4f,%.4f)\n",
                       dim, rank, n_veneno, n_errada, pior, onde[0], onde[1]);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0)
            printf("dim %d: %ld conferidas (%ld de franja, %ld delas na interface), "
                   "%ld veneno, %ld erradas -> %s\n",
                   dim, gs[2], gs[3], gs[4], gs[0], gs[1],
                   (gs[0] || gs[1]) ? "FALHOU" : "ok");
        // --------------------------------------------------------------
        // Segunda afirmacao: o CONDICIONAMENTO do estencil.
        //
        // Um WLS de ordem 2 reproduz campo linear exatamente POR CONSTRUCAO --
        // e' um ajuste de minimos quadrados com base linear sobre dado linear.
        // Logo, "reproduziu o campo linear" nao prova que o operador e' sao:
        // ele pode reproduzir a base e ainda assim ter coeficientes enormes,
        // que amplificam tudo que nao esta' na base.  A medida que importa e'
        // soma|w|, o fator de amplificacao.  Numa interpolacao sa ele fica na
        // ordem de 1.  Se num sitio ele explode, aquele sitio e' um detonador:
        // qualquer ruido ali vira um pico em um passo.
        //
        // E' exatamente o que o solver faz: sfd_get_stencil no centro da
        // faceta deslocado de +-h em cada direcao (hig-flow-step.c:890).
        {
            sim_stencil *stn = stn_create();
            real pior_amp = 0.0;
            Point onde_amp;
            POINT_ASSIGN_SCALAR(onde_amp, 0.0);
            real pior_lin = 0.0;
            sim_domain *cd = sfd->cdom;
            for (int k = 0; k < sd_get_num_higtrees(cd); k++) {
                hig_cell *root = sd_get_higtree(cd, k);
                Point blo, bhi;
                POINT_ASSIGN_SCALAR(blo, -1.0);
                POINT_ASSIGN_SCALAR(bhi,  2.0);
                higfit_facetiterator *fit;
                for (fit = higfit_create_bounding_box_facets(root,
                            sfd->dimofinterest, blo, bhi);
                     !higfit_isfinished(fit); higfit_nextfacet(fit)) {
                    hig_facet *f = higfit_getfacet(fit);
                    int lid = mp_lookup(sfd->fm, hig_get_fid(f));
                    if (lid < 0 || lid >= nloc) continue;   // so' as proprias
                    Point fc;
                    hig_get_facet_center(f, fc);
                    hig_cell *cel = hig_get_facet_cell(f);
                    Point cl, ch;
                    hig_get_lowpoint(cel, cl);
                    hig_get_highpoint(cel, ch);
                    for (int d = 0; d < DIM; d++) {
                        const real h = ch[d] - cl[d];
                        for (int s = -1; s <= 1; s += 2) {
                            Point p;
                            POINT_ASSIGN(p, fc);
                            p[d] += s * h;
                            stn_reset(stn);
                            sfd_get_stencil(sfd, fc, p, 1.0, stn);
                            const int ne = stn_get_numelems(stn);
                            int  *ids = stn_get_ids(stn);
                            real *w   = stn_get_vals(stn);
                            real amp = 0.0, val = 0.0;
                            for (int j = 0; j < ne; j++) {
                                amp += fabs(w[j]);
                                if (ids[j] >= 0 && ids[j] < ntot)
                                    val += w[j] * dp_get_value(dp, ids[j]);
                            }
                            if (!(amp <= pior_amp)) {
                                pior_amp = amp;
                                POINT_ASSIGN(onde_amp, p);
                            }
                            // Campo linear: o estencil (mais o termo de
                            // contorno, que aqui nao ha') deve dar campo(p).
                            if (ne > 0) {
                                const real e = fabs(val - campo(p));
                                if (!(e <= pior_lin)) pior_lin = e;
                            }
                        }
                    }
                }
                higfit_destroy(fit);
            }
            stn_destroy(stn);

            struct { double v; int r; } loc, glo;
            loc.v = pior_amp; loc.r = rank;
            MPI_Allreduce(&loc, &glo, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
            double lin_g = 0.0;
            MPI_Allreduce(&pior_lin, &lin_g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            if (rank == glo.r)
                printf("dim %d: AMPLIFICACAO maxima soma|w| = %.4e em "
                       "(%.4f,%.4f) no rank %d ; erro no campo linear = %.3e\n",
                       dim, glo.v, onde_amp[0], onde_amp[1], rank, lin_g);
            MPI_Barrier(MPI_COMM_WORLD);
        }

        dp_destroy(dp);
    }

    if (rank == 0)
        printf("\n%s (np = %d)\n", falhou ? "FALHOU" : "todas as afirmacoes passaram",
               ntasks);
    PetscFinalize();
    MPI_Finalize();
    return falhou;
}
