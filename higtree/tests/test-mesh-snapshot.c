// A FRONTEIRA DAS CONSULTAS: o instantaneo, e os dois backends que o preenchem.
//
// POR QUE ESTE TESTE EXISTE.  Das 907 chamadas de consulta dentro de laco quente
// em higflow/src, 777 -- 86% -- sao leitura de centro, tamanho e indice.  Nenhuma
// precisa da arvore.  O `hig_mesh_snapshot` as atende com arranjo plano preenchido
// uma vez por producao, e com isso o backend passa a ENTREGAR o que as consultas
// precisam, em vez de entregar uma estrutura que elas saibam navegar.
//
// O que se afirma aqui, em ordem:
//
//   instantaneo_do_mtree_bate_com_a_arvore   o arranjo diz o mesmo que as
//       consultas antigas diriam, celula a celula.  Sem isto, mover a fronteira
//       seria trocar uma resposta por outra, nao pela mesma.
//
//   indice_e_o_id_local                      center[i*DIM] e' a celula de id i.
//       E' o que apaga as 190 chamadas de identificador do laco.
//
//   instantaneo_ignora_a_franja              a franja da' suporte a' interpolacao
//       e nao pertence ao dominio sobre o qual se resolve (clausula C6).
//
//   t8code_preenche_o_mesmo_conjunto         [so' com -DHIGTREE_COM_T8CODE]
//       o t8code preenche o instantaneo DIRETO da floresta, sem materializar
//       octree nenhum, e o conjunto de celulas e' o mesmo que o MTree produz.
//       E' este caso que mostra a fronteira movida: os dois backends entregam o
//       mesmo dado por caminhos que nao se parecem.
//
// A comparacao e' por CONJUNTO, nao posicao a posicao: os dois backends numeram
// as celulas em ordens diferentes, porque a ordem e' consequencia de como cada um
// percorre.  Exigir a mesma ordem seria exigir a mesma implementacao.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "higtree.h"
#include "higtree-iterator.h"
#include "domain.h"
#include "hig-mesh-snapshot.h"
#include "utils.h"
#include "testing.h"

#ifdef HIGTREE_COM_T8CODE
#include "t8code/t8-mesh-producer.h"
#endif

#define NC   4

static hig_cell *bloco(real x0, real x1) {
    Point lo, hi;
    POINT_ASSIGN_SCALAR(lo, 0.0);
    POINT_ASSIGN_SCALAR(hi, 1.0);
    lo[0] = x0; hi[0] = x1;
    hig_cell *r = hig_create_root(lo, hi);
    int nc[DIM];
    for (int d = 0; d < DIM; d++) nc[d] = NC;
    nc[0] = 2;
    hig_refine_uniform(r, nc);
    return r;
}

// Malha equivalente a' do produtor do t8code: raiz 4 por direcao, a celula que
// contem 0,375 refinada, e uma neta dela refinada de novo.
static sim_domain *dominio_nao_graduado(void) {
    Point lo, hi;
    POINT_ASSIGN_SCALAR(lo, 0.0);
    POINT_ASSIGN_SCALAR(hi, 1.0);
    hig_cell *raiz = hig_create_root(lo, hi);
    int nc[DIM];
    for (int d = 0; d < DIM; d++) nc[d] = NC;
    hig_refine_uniform(raiz, nc);
    for (int passo = 0; passo < 2; passo++) {
        Point p;
        POINT_ASSIGN_SCALAR(p, passo == 0 ? 0.375 : 0.3125);
        hig_cell *c = hig_get_cell_with_point(raiz, p);
        if (c == NULL) return NULL;
        int n2[DIM];
        for (int d = 0; d < DIM; d++) n2[d] = 2;
        hig_refine_uniform(c, n2);
    }
    sim_domain *sd = sd_create(NULL);
    sd_add_higtree(sd, raiz);
    mp_mapper *m = sd_get_domain_mapper(sd);
    higcit_celliterator *cit = sd_get_domain_celliterator(sd);
    mp_assign_from_celliterator(m, cit, 0);
    higcit_destroy(cit);
    return sd;
}

int main(int argc, char *argv[]) {
    higtree_initialize(&argc, &argv);   // o produtor do t8code exige MPI

    sim_domain *sd = dominio_nao_graduado();
    T_CHECK_MSG(sd != NULL, "nao consegui montar o dominio nao graduado");
    mp_mapper *m = sd_get_domain_mapper(sd);

    hig_mesh_snapshot *s = hms_from_domain(sd);
    T_CHECK_MSG(s != NULL, "hms_from_domain devolveu NULL");

    // ------------------------------------------------------------------
    t_case("instantaneo_do_mtree_bate_com_a_arvore");
    {
        int visitadas = 0, divergentes = 0;
        char primeira[256]; primeira[0] = '\0';
        higcit_celliterator *cit;
        for (cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit);
             higcit_nextcell(cit)) {
            hig_cell *c = higcit_getcell(cit);
            const int i = mp_lookup(m, hig_get_cid(c));
            Point ce, de;
            hig_get_center(c, ce);
            hig_get_delta(c, de);
            visitadas++;
            for (int d = 0; d < DIM; d++) {
                if (fabs(s->center[i * DIM + d] - ce[d]) > 1e-15 ||
                    fabs(s->delta[i * DIM + d]  - de[d]) > 1e-15) {
                    if (!divergentes) {
                        snprintf(primeira, sizeof primeira,
                                 "celula %d, direcao %d: instantaneo (%.17g, %.17g) "
                                 "contra arvore (%.17g, %.17g)", i, d,
                                 (double) s->center[i * DIM + d],
                                 (double) s->delta[i * DIM + d],
                                 (double) ce[d], (double) de[d]);
                    }
                    divergentes++;
                    break;
                }
            }
        }
        higcit_destroy(cit);
        T_CHECK_MSG(visitadas == s->n,
            "o iterador visitou %d celulas, o instantaneo tem %d", visitadas, s->n);
        T_CHECK_MSG(divergentes == 0,
            "%d celula(s) divergem entre o instantaneo e a arvore.  %s",
            divergentes, primeira);
    }

    // ------------------------------------------------------------------
    t_case("indice_e_o_id_local");
    {
        // Para toda celula, o id que o mapeador devolve indexa a propria celula
        // no instantaneo.  E' o que permite tirar hig_get_cid do laco.
        int fora = 0;
        higcit_celliterator *cit;
        for (cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit);
             higcit_nextcell(cit)) {
            hig_cell *c = higcit_getcell(cit);
            const int i = mp_lookup(m, hig_get_cid(c));
            if (i < 0 || i >= s->n) { fora++; continue; }
            Point ce;
            hig_get_center(c, ce);
            for (int d = 0; d < DIM; d++) {
                if (fabs(s->center[i * DIM + d] - ce[d]) > 1e-15) { fora++; break; }
            }
        }
        higcit_destroy(cit);
        T_CHECK_MSG(fora == 0,
            "%d celula(s) cujo id nao indexa a propria posicao no instantaneo", fora);
    }

    // ------------------------------------------------------------------
    t_case("instantaneo_ignora_a_franja");
    {
        // Mesma malha, mas com um bloco acrescentado como FRANJA.  O instantaneo
        // nao pode crescer: a franja da' suporte, nao pertence ao dominio (C6).
        sim_domain *sf = sd_create(NULL);
        sd_add_higtree(sf, bloco(0.0, 0.5));
        sd_add_fringe_higtree(sf, bloco(0.5, 1.0));
        mp_mapper *mf = sd_get_domain_mapper(sf);
        higcit_celliterator *cit = sd_get_domain_celliterator(sf);
        unsigned livre = mp_assign_from_celliterator(mf, cit, 0);
        higcit_destroy(cit);
        for (unsigned k = 0; k < sd_get_num_fringe_higtrees(sf); ++k) {
            cit = higcit_create_all_leaves(sd_get_fringe_higtree(sf, k));
            livre = mp_assign_from_celliterator(mf, cit, livre);
            higcit_destroy(cit);
        }
        hig_mesh_snapshot *sfs = hms_from_domain(sf);
        T_CHECK_MSG(sfs != NULL, "hms_from_domain devolveu NULL com franja");
        int esperado = 2;
        for (int d = 1; d < DIM; d++) esperado *= NC;
        T_CHECK_MSG(sfs != NULL && sfs->n == esperado,
            "o instantaneo tem %d celulas e o bloco local tem %d: a franja entrou",
            sfs ? sfs->n : -1, esperado);
        hms_destroy(sfs);
    }

    // ------------------------------------------------------------------
#ifdef HIGTREE_COM_T8CODE
    t_case("t8code_preenche_o_mesmo_conjunto");
    {
        hig_mesh_snapshot *t8 = t8_preenche_instantaneo();
        T_CHECK_MSG(t8 != NULL,
            "t8_preenche_instantaneo devolveu NULL -- a floresta nao tem o salto");
        if (t8 != NULL) {
            char detalhe[256];
            const int sem_par = hms_mesmo_conjunto(s, t8, 1e-12, detalhe);
            T_CHECK_MSG(sem_par == 0,
                "o t8code preencheu conjunto diferente do MTree: %d celula(s) "
                "sem par.  %s", sem_par, detalhe);
            hms_destroy(t8);
        }
    }
#endif

    hms_destroy(s);
    return t_end();
}
