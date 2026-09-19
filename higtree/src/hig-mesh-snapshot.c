// A fronteira das consultas.  Ver hig-mesh-snapshot.h para o porque.

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "hig-mesh-snapshot.h"
#include "higtree.h"
#include "higtree-iterator.h"
#include "mapper.h"

hig_mesh_snapshot *
hms_create(int n)
{
    hig_mesh_snapshot *s = (hig_mesh_snapshot *) malloc(sizeof *s);
    if (s == NULL) return NULL;
    s->n = n;
    s->dim = DIM;
    s->center = (real *) calloc((size_t) n * DIM, sizeof *s->center);
    s->delta  = (real *) calloc((size_t) n * DIM, sizeof *s->delta);
    if (s->center == NULL || s->delta == NULL) { hms_destroy(s); return NULL; }
    return s;
}

void
hms_destroy(hig_mesh_snapshot *s)
{
    if (s == NULL) return;
    free(s->center);
    free(s->delta);
    free(s);
}

hig_mesh_snapshot *
hms_from_domain(sim_domain *sd)
{
    mp_mapper *m = sd_get_domain_mapper(sd);

    // Conta primeiro: o iterador de dominio ve' so' as celulas LOCAIS (C6), e e'
    // essa contagem que define o tamanho do instantaneo.
    int n = 0;
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(sd); !higcit_isfinished(it);
         higcit_nextcell(it)) {
        n++;
    }
    higcit_destroy(it);

    hig_mesh_snapshot *s = hms_create(n);
    if (s == NULL) return NULL;

    for (it = sd_get_domain_celliterator(sd); !higcit_isfinished(it);
         higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        const int i = mp_lookup(m, hig_get_cid(c));
        if (i < 0 || i >= n) {          // mapeador incoerente com o iterador
            higcit_destroy(it);
            hms_destroy(s);
            return NULL;
        }
        Point ce, de;
        hig_get_center(c, ce);
        hig_get_delta(c, de);
        for (int d = 0; d < DIM; d++) {
            s->center[i * DIM + d] = ce[d];
            s->delta[i * DIM + d]  = de[d];
        }
    }
    higcit_destroy(it);
    return s;
}

int
hms_mesmo_conjunto(const hig_mesh_snapshot *a, const hig_mesh_snapshot *b,
                   real tol, char *detalhe)
{
    if (detalhe) detalhe[0] = '\0';
    if (a == NULL || b == NULL) {
        if (detalhe) snprintf(detalhe, 256, "instantaneo nulo");
        return -1;
    }
    if (a->n != b->n) {
        if (detalhe)
            snprintf(detalhe, 256, "contagens diferentes: %d contra %d", a->n, b->n);
        return a->n > b->n ? a->n - b->n : b->n - a->n;
    }

    // O(n^2) de proposito: isto e' verificacao, nao caminho quente, e uma busca
    // linear nao tem estrutura de indice para sair errada.
    char *usado = (char *) calloc((size_t) b->n, 1);
    if (usado == NULL) return -1;

    int sem_par = 0;
    for (int i = 0; i < a->n; i++) {
        int achou = -1;
        for (int j = 0; j < b->n && achou < 0; j++) {
            if (usado[j]) continue;
            int bate = 1;
            for (int d = 0; d < DIM && bate; d++) {
                if (fabs(a->center[i * DIM + d] - b->center[j * DIM + d]) > tol)
                    bate = 0;
            }
            if (bate) achou = j;
        }
        if (achou < 0) {
            if (sem_par == 0 && detalhe) {
                snprintf(detalhe, 256,
                         "celula %d de a, centro (%.6f, %.6f), sem par em b", i,
                         (double) a->center[i * DIM],
                         (double) a->center[i * DIM + (DIM > 1 ? 1 : 0)]);
            }
            sem_par++;
            continue;
        }
        usado[achou] = 1;
        for (int d = 0; d < DIM; d++) {
            if (fabs(a->delta[i * DIM + d] - b->delta[achou * DIM + d]) > tol) {
                if (sem_par == 0 && detalhe) {
                    snprintf(detalhe, 256,
                             "celula em (%.6f, ...) com tamanhos diferentes: "
                             "%.9f contra %.9f na direcao %d",
                             (double) a->center[i * DIM],
                             (double) a->delta[i * DIM + d],
                             (double) b->delta[achou * DIM + d], d);
                }
                sem_par++;
                break;
            }
        }
    }
    free(usado);
    return sem_par;
}
