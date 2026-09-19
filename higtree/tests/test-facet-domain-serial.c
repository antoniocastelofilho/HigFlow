// Montagem SERIAL de um `sim_facet_domain`, do zero ate' uma consulta de estencil.
//
// POR QUE ESTE TESTE EXISTE.  Ate' 2026-09-19 esta sequencia nao era possivel: o
// `sfd->sfbi[]` -- de onde o iterador de facetas tira quais facetas sao do bloco --
// so' era preenchido por `psfd_compute_sfbi`, no caminho PARTICIONADO.  `sfd_create`
// deixava o vetor nulo, e a primeira consulta estourava.  Nao era sequencia
// faltando: varrendo higtree e higflow, nao havia um unico chamador serial de
// `sfd_set_sfbi`.
//
// A consequencia era arquitetural, nao um incomodo de teste: metade do contrato de
// dominio de facetas nao podia ser exercitada sem MPI, e interface que so' existe
// acoplada ao particionamento e' exatamente a dificil de substituir.  A metade local
// passou a viver no domain.c como `sfd_compute_sfbi`; este teste e' o que impede que
// ela volte a depender do particionamento sem ninguem perceber.
//
// O que ele afirma, em ordem:
//   1. antes de `sfd_compute_sfbi`, o bloco e' nulo -- o estado de partida
//   2. depois, nao e'
//   3. o numero de facetas mapeadas e' o das facetas INTERIORES, calculado da malha
//   4. uma consulta de estencil responde, com pesos que somam 1 (particao da unidade)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "higtree.h"
#include "higtree-iterator.h"
#include "domain.h"
#include "testing.h"

#define NC 4                     // celulas por direcao

int main(void) {
    Point lo, hi;
    POINT_ASSIGN_SCALAR(lo, 0.0);
    POINT_ASSIGN_SCALAR(hi, 1.0);
    hig_cell *raiz = hig_create_root(lo, hi);
    int nc[DIM];
    for (int d = 0; d < DIM; d++) nc[d] = NC;
    hig_refine_uniform(raiz, nc);

    sim_domain *sd = sd_create(NULL);
    sd_set_interpolator_order(sd, 2);
    sd_add_higtree(sd, raiz);
    mp_mapper *mc = sd_get_domain_mapper(sd);
    higcit_celliterator *cit = sd_get_domain_celliterator(sd);
    mp_assign_from_celliterator(mc, cit, 0);
    higcit_destroy(cit);

    sim_facet_domain *sfd = sfd_create(NULL, 0);   // facetas normais a x
    sfd_set_interpolator_order(sfd, 2);
    sfd_copy_higtrees_from_center_domain(sfd, sd);
    sfd_adjust_facet_ids(sfd);

    t_case("sfbi_nulo_antes_de_computar");
    T_CHECK_MSG(sfd_get_sfbi(sfd, 0) == NULL,
        "sfd_create deveria deixar sfbi[0] nulo, e nao deixou -- o estado de "
        "partida deste teste nao e' mais o que ele supoe");

    t_case("montagem_serial_preenche_o_bloco_de_faceta");
    sfd_compute_sfbi(sfd);
    T_CHECK_MSG(sfd_get_sfbi(sfd, 0) != NULL,
        "sfd_compute_sfbi nao preencheu sfbi[0]: a montagem serial voltou a "
        "depender do caminho particionado");

    t_case("facetas_interiores_mapeadas");
    mp_mapper *mf = sfd_get_domain_mapper(sfd);
    higfit_facetiterator *fit = sfd_get_domain_facetiterator(sfd);
    int nf = mp_assign_from_facetiterator(mf, fit, 0);
    higfit_destroy(fit);

    // Sem condicao de contorno registrada, so' as facetas INTERIORES normais a x
    // pertencem ao bloco: (NC-1) por linha, vezes NC nas demais direcoes.
    int esperado = NC - 1;
    for (int d = 1; d < DIM; d++) esperado *= NC;
    T_CHECK_MSG(nf == esperado,
        "facetas mapeadas: obtido %d, esperado %d (as interiores normais a x "
        "numa malha %d^%d)", nf, esperado, NC, DIM);

    t_case("consulta_de_estencil_responde_e_soma_um");
    sim_stencil *stn = stn_create();
    Point centro, x;
    POINT_ASSIGN_SCALAR(centro, 0.5);
    POINT_ASSIGN_SCALAR(x, 0.5);
    sfd_get_stencil(sfd, centro, x, 1.0, stn);

    int n = stn_get_numelems(stn);
    T_CHECK_MSG(n > 0, "o estencil voltou vazio (numelems = %d)", n);

    // Particao da unidade: interpolar a funcao constante 1 tem de dar 1, entao os
    // pesos somam 1.  Nao consulta referencia gravada.
    real *w = stn_get_vals(stn);
    real soma = 0.0;
    for (int i = 0; i < n; i++) soma += w[i];
    T_NEAR(soma, 1.0, 1.0e-12, "soma dos pesos do estencil (particao da unidade)");

    stn_destroy(stn);
    return t_end();
}
