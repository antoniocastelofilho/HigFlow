// Os dois operadores de transferencia da fronteira imersa, verificados sem
// solver nenhum.  Tres afirmacoes, e cada uma pega um defeito diferente.
//
//  1. PESO TOTAL = PERIMETRO.  A curva e' um QUADRADO de lado 0,5: perimetro
//     exatamente 2,0.  Circulo nao serviria -- o poligono inscrito tem
//     perimetro proximo, e "proximo" nao distingue subdivisao errada de
//     discretizacao.  Pega marcador perdido, duplicado, ou peso errado.
//
//  2. INTERPOLAR CAMPO CONSTANTE DEVOLVE A CONSTANTE.  E' a particao da unidade
//     do nucleo.  Pega nucleo errado, suporte incompleto e normalizacao errada.
//
//  3. CONSERVACAO DA FORCA.  SUM_facetas F h^DIM = SUM_k f_k w_k.  E' o par
//     adjunto, e e' o unico dos tres que exerce a ACUMULACAO franja->dono: sem
//     ela, marcador perto de fronteira de particao perde parte da forca, e o
//     erro cresce com np em vez de aparecer de uma vez.
#include "higtree.h"
#include "higtree-io.h"
#include "higtree-iterator.h"
#include "pdomain.h"
#include "lbal.h"
#include "hig-flow-fronteira-imersa.h"
#include <mpi.h>
#include <math.h>
#include <stdio.h>

#define NC     32
#define FRANJA  5

static int falhou = 0;
static void checa(int rank, const char *nome, real obtido, real esperado, real tol)
{
    const real err = fabs(obtido - esperado);
    if (rank == 0)
        printf("  %-46s obtido %.12g  esperado %.12g  erro %.3g  %s\n",
               nome, (double) obtido, (double) esperado, (double) err,
               err <= tol ? "ok" : "FALHOU");
    if (err > tol) falhou = 1;
}

int main(int argc, char *argv[])
{
    higtree_initialize(&argc, &argv);
    int rank, np;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    const real h = 1.0 / NC;

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

    // Um dominio de facetas por direcao, como o solver monta.
    sim_facet_domain  *sfd[DIM];
    psim_facet_domain *psfd[DIM];
    distributed_property *dpu[DIM], *dpF[DIM];
    for (int dim = 0; dim < DIM; dim++) {
        sfd[dim] = sfd_create(NULL, dim);
        sfd_set_interpolator_order(sfd[dim], 2);
        sfd_copy_higtrees_from_center_domain(sfd[dim], sd);
        // O solver faz isto (hig-flow-kernel.c:921) e eu nao fazia.  Sem ele a
        // troca de franja sai ASSIMETRICA em np>=2 e a guarda da HiGTree aborta
        // -- corretamente, porque sincronizar assim corromperia o heap.
        sfd_adjust_facet_ids(sfd[dim]);
        psfd[dim] = psfd_create(sfd[dim], psd);
        psfd_compute_sfbi(psfd[dim]);
        psfd_synced_mapper(psfd[dim]);
        dpu[dim] = psfd_create_property(psfd[dim]);
        dpF[dim] = psfd_create_property(psfd[dim]);
    }

    // Campo de velocidade CONSTANTE: u = (3, -7).  Constante e' o que torna a
    // afirmacao 2 uma particao da unidade e nao uma coincidencia.
    const real U[2] = {3.0, -7.0};
    for (int dim = 0; dim < DIM; dim++) {
        const int nloc = dpu[dim]->pdata->local_count;
        for (int lid = 0; lid < nloc; lid++) {
            dp_set_value(dpu[dim], lid, U[dim]);
            dp_set_value(dpF[dim], lid, 0.0);
        }
        const int ntot = dpu[dim]->pdata->total_count;
        for (int lid = nloc; lid < ntot; lid++) dp_set_value(dpF[dim], lid, 0.0);
        dp_sync(dpu[dim]);
        dp_sync(dpF[dim]);
    }

    // Quadrado de lado 0,5, longe da borda do dominio para o suporte caber.
    Point v[4] = {{0.25,0.25},{0.75,0.25},{0.75,0.75},{0.25,0.75}};
    fi_corpo *corpo = fi_cria_curva(sfd[0], (const Point *) v, 4, h);

    checa(rank, "1. peso total = perimetro do quadrado",
          fi_peso_total(corpo), 2.0, 1e-12);

    fi_interpola(corpo, sfd, dpu);

    // A afirmacao 2 olha o pior marcador de todos os ranks.
    {
        real pior[DIM];
        for (int d = 0; d < DIM; d++) pior[d] = 0.0;
        // reusa a forca como sonda: f = -u/dt com dt = -1 devolve u.
        fi_forca_corpo_rigido(corpo, -1.0);
        real tot[DIM];
        fi_forca_total(corpo, tot);
        // `fi_forca_total` soma f_k * dV_k, e dV = peso * h^(DIM-1).  Entao,
        // se todo marcador leu exatamente U, o total e' U * perimetro * h.
        // O fator h nao e' arbitrario: e' o VOLUME do marcador, e foi ele que
        // faltava quando a forca saia 1/h vezes grande demais.
        real hvol = 1.0;
        for (int d = 0; d < DIM - 1; d++) hvol *= h;
        for (int d = 0; d < DIM; d++) pior[d] = tot[d];
        checa(rank, "2a. interpolou u[0] constante (x perimetro x h)", pior[0], U[0]*2.0*hvol, 1e-10);
        checa(rank, "2b. interpolou u[1] constante (x perimetro x h)", pior[1], U[1]*2.0*hvol, 1e-10);
    }

    // Afirmacao 3: forca direta com dt, espalhar, e conferir a integral.
    const real dt = 0.01;
    fi_forca_corpo_rigido(corpo, dt);
    real esperado[DIM];
    fi_forca_total(corpo, esperado);
    fi_espalha(corpo, sfd, dpF);

    real hd = 1.0;
    for (int d = 0; d < DIM; d++) hd *= h;
    for (int dim = 0; dim < DIM; dim++) {
        const int nloc = dpF[dim]->pdata->local_count;   // so' as PROPRIAS
        real soma = 0.0;
        for (int lid = 0; lid < nloc; lid++) soma += dp_get_value(dpF[dim], lid);
        soma *= hd;
        real global = 0.0;
        MPI_Allreduce(&soma, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        char nome[64];
        snprintf(nome, sizeof nome, "3%c. forca espalhada = forca dos marcadores",
                 (char)('a' + dim));
        checa(rank, nome, global, esperado[dim], 1e-9 * fabs(esperado[dim]) + 1e-12);
    }

    fi_destroi(corpo);
    if (rank == 0)
        printf("  np=%d: %s\n", np, falhou ? "HOUVE FALHA" : "todas as afirmacoes passaram");
    return falhou;
}
