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
    s->low  = (real *) calloc((size_t) n * DIM, sizeof *s->low);
    s->high = (real *) calloc((size_t) n * DIM, sizeof *s->high);
    if (s->low == NULL || s->high == NULL) { hms_destroy(s); return NULL; }
    return s;
}

void
hms_destroy(hig_mesh_snapshot *s)
{
    if (s == NULL) return;
    free(s->low);
    free(s->high);
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
        // A CAIXA, copiada como a celula a guarda -- centro e delta saem dela
        // por derivacao, com as mesmas contas do `hig_get_center`.
        hig_get_lowpoint(c, ce);
        hig_get_highpoint(c, de);
        for (int d = 0; d < DIM; d++) {
            s->low[i * DIM + d]  = ce[d];
            s->high[i * DIM + d] = de[d];
        }
    }
    higcit_destroy(it);
    return s;
}

hig_facet_snapshot *
hfs_create(int n)
{
    hig_facet_snapshot *s = (hig_facet_snapshot *) calloc(1, sizeof *s);
    if (s == NULL) return NULL;
    s->n = n;
    s->low  = (real *) calloc((size_t) n * DIM, sizeof *s->low);
    s->high = (real *) calloc((size_t) n * DIM, sizeof *s->high);
    s->dim  = (signed char *) calloc((size_t) n, sizeof *s->dim);
    s->dir  = (signed char *) calloc((size_t) n, sizeof *s->dir);
    if (s->low == NULL || s->high == NULL || s->dim == NULL || s->dir == NULL) {
        hfs_destroy(s);
        return NULL;
    }
    return s;
}

void
hfs_destroy(hig_facet_snapshot *s)
{
    if (s == NULL) return;
    free(s->low);
    free(s->high);
    free(s->dim);
    free(s->dir);
    free(s);
}

hig_facet_snapshot *
hfs_from_facet_domain(sim_facet_domain *sfd)
{
    mp_mapper *m = sfd_get_domain_mapper(sfd);

    int n = 0;
    higfit_facetiterator *fit;
    for (fit = sfd_get_domain_facetiterator(sfd); !higfit_isfinished(fit);
         higfit_nextfacet(fit)) {
        n++;
    }
    higfit_destroy(fit);

    hig_facet_snapshot *s = hfs_create(n);
    if (s == NULL) return NULL;

    for (fit = sfd_get_domain_facetiterator(sfd); !higfit_isfinished(fit);
         higfit_nextfacet(fit)) {
        hig_facet *f = higfit_getfacet(fit);
        const int i = mp_lookup(m, hig_get_fid(f));
        if (i < 0 || i >= n) {          // mapeador incoerente com o iterador
            higfit_destroy(fit);
            hfs_destroy(s);
            return NULL;
        }
        // A caixa da CELULA da faceta, copiada.  Centro e tamanho da faceta saem
        // dela pelas mesmas contas do `hig_get_facet_center`/`_delta`.
        hig_cell *c = hig_get_facet_cell(f);
        Point lo, hi;
        hig_get_lowpoint(c, lo);
        hig_get_highpoint(c, hi);
        for (int d = 0; d < DIM; d++) {
            s->low[i * DIM + d]  = lo[d];
            s->high[i * DIM + d] = hi[d];
        }
        s->dim[i] = (signed char) f->dim;
        s->dir[i] = (signed char) f->dir;
    }
    higfit_destroy(fit);
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
        Point ca, da;
        hms_center(a, i, ca);
        hms_delta(a, i, da);
        int achou = -1;
        for (int j = 0; j < b->n && achou < 0; j++) {
            if (usado[j]) continue;
            Point cb;
            hms_center(b, j, cb);
            int bate = 1;
            for (int d = 0; d < DIM && bate; d++) {
                if (fabs(ca[d] - cb[d]) > tol)
                    bate = 0;
            }
            if (bate) achou = j;
        }
        if (achou < 0) {
            if (sem_par == 0 && detalhe) {
                snprintf(detalhe, 256,
                         "celula %d de a, centro (%.6f, %.6f), sem par em b", i,
                         (double) ca[0],
                         (double) ca[DIM > 1 ? 1 : 0]);
            }
            sem_par++;
            continue;
        }
        usado[achou] = 1;
        Point db;
        hms_delta(b, achou, db);
        for (int d = 0; d < DIM; d++) {
            if (fabs(da[d] - db[d]) > tol) {
                if (sem_par == 0 && detalhe) {
                    snprintf(detalhe, 256,
                             "celula em (%.6f, ...) com tamanhos diferentes: "
                             "%.9f contra %.9f na direcao %d",
                             (double) ca[0],
                             (double) da[d],
                             (double) db[d], d);
                }
                sem_par++;
                break;
            }
        }
    }
    free(usado);
    return sem_par;
}
