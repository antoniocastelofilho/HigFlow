// Refino adaptativo compartilhado pelos exemplos bifasicos (VOF e front-tracking).
//
// Veio de example2d_DynamicMeshAdapt/mesh_adapt_function.c, que e' o unico lugar
// do sistema onde o dominio ja' era reconstruido no meio da corrida com
// multifase -- o higflow_reconstroi_dominio da fronteira imersa aborta em
// qualquer coisa que nao seja NEWTONIANO (hig-flow-kernel.c:91).
//
// O CRITERIO: celulas com 0,001 < fracvol < 0,999 sao SEMENTES de interface; uma
// celula qualquer recebe o nivel L se estiver a menos de `limiar[L-1]` da semente
// mais proxima.  Os limiares sao DISTANCIAS e sao decrescentes -- o ultimo e' o
// do nivel mais fino.
//
// A REGRA DAS CELULAS MINIMAS, e por que ela e' derivada e nao digitada:
// o VOF precisa de pelo menos N celulas do nivel MAIS FINO de cada lado da
// interface para reconstruir a normal e a curvatura (o estencil de altura
// precisa da coluna inteira dentro do nivel fino).  Um limiar digitado a olho
// satisfaz ou nao satisfaz essa condicao sem avisar, e o sintoma seria curvatura
// ruim -- nao um erro.  Aqui o limiar SAI da regra (malha_adapt_limiares) e a
// regra e' CONFERIDA POR MEDIDA na malha de verdade (malha_adapt_mede_banda).

// Era INCLUIDO no exemplo (depois do cabecalho do kernel); como unidade de
// compilacao propria precisa declarar o solver por conta.
#include "hig-flow-kernel.h"
#include "hig-flow-io.h"
#include "hig-flow-bc.h"
#include "malha-adaptativa.h"
#include "higtree.h"
#include "higtree-io.h"
#include "higtree-iterator.h"
#include "pdomain.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glib.h>
#include <float.h>
#include <mpi.h>

typedef struct { Point center; } InterfaceSeedAdapt;

static real dist_sq_c(Point p1, Point p2) {
    real d = 0.0;
    for (int i = 0; i < DIM; i++) d += (p1[i]-p2[i])*(p1[i]-p2[i]);
    return d;
}

// --- Cell level cache (O(1) lookup, avoids walking to the root every time) ---
static int get_cached_level(GHashTable *cache, hig_cell *c) {
    if (c == NULL) return -1;
    gpointer val = g_hash_table_lookup(cache, c);
    if (val != NULL) return GPOINTER_TO_INT(val);
    int level = 0;
    hig_cell *p = hig_get_parent(c);
    while (p != NULL) { level++; p = hig_get_parent(p); }
    level--;
    g_hash_table_insert(cache, c, GINT_TO_POINTER(level));
    return level;
}

// --- 2D spatial hash for fast neighbour-seed queries ---
typedef struct {
    int nx, ny;
    real bin_w, bin_h, ox, oy;
    int *start;   // size nx*ny+1 (offset into seeds array)
    int *seeds;   // flat array of seed indices per bin
} SeedHash;

static void build_seed_hash(SeedHash *sh, InterfaceSeedAdapt *seeds, int n,
                             real bin_w, real bin_h, Point lo, Point hi)
{
    real sx = hi[0] - lo[0], sy = hi[1] - lo[1];
    sh->nx = (int)ceil(sx / bin_w); if (sh->nx < 1) sh->nx = 1;
    sh->ny = (int)ceil(sy / bin_h); if (sh->ny < 1) sh->ny = 1;
    sh->bin_w = bin_w; sh->bin_h = bin_h;
    sh->ox = lo[0]; sh->oy = lo[1];

    int nb = sh->nx * sh->ny;
    int *cnt = (int *) calloc(nb, sizeof(int));
    for (int s = 0; s < n; s++) {
        int bx = (int)((seeds[s].center[0] - lo[0]) / bin_w);
        int by = (int)((seeds[s].center[1] - lo[1]) / bin_h);
        if (bx < 0) bx = 0; if (bx >= sh->nx) bx = sh->nx - 1;
        if (by < 0) by = 0; if (by >= sh->ny) by = sh->ny - 1;
        cnt[by * sh->nx + bx]++;
    }
    sh->start = (int *) malloc((nb + 1) * sizeof(int));
    int offset = 0;
    for (int i = 0; i < nb; i++) { sh->start[i] = offset; offset += cnt[i]; }
    sh->start[nb] = offset;
    sh->seeds = (int *) malloc(n * sizeof(int));
    int *pos = (int *) malloc(nb * sizeof(int));
    memcpy(pos, sh->start, nb * sizeof(int));
    for (int s = 0; s < n; s++) {
        int bx = (int)((seeds[s].center[0] - lo[0]) / bin_w);
        int by = (int)((seeds[s].center[1] - lo[1]) / bin_h);
        if (bx < 0) bx = 0; if (bx >= sh->nx) bx = sh->nx - 1;
        if (by < 0) by = 0; if (by >= sh->ny) by = sh->ny - 1;
        sh->seeds[pos[by * sh->nx + bx]++] = s;
    }
    free(pos); free(cnt);
}

static void free_seed_hash(SeedHash *sh) {
    free(sh->start); free(sh->seeds); memset(sh, 0, sizeof(*sh));
}

// Smallest squared distance from a point p to any seed, via the spatial hash.
static real min_dist2_to_seeds(Point p, InterfaceSeedAdapt *seeds,
                                SeedHash *sh)
{
    int bx = (int)((p[0] - sh->ox) / sh->bin_w);
    int by = (int)((p[1] - sh->oy) / sh->bin_h);
    real d2_min = DBL_MAX;
    for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {
            int cx = bx + dx, cy = by + dy;
            if (cx < 0 || cx >= sh->nx || cy < 0 || cy >= sh->ny) continue;
            int idx = cy * sh->nx + cx;
            for (int i = sh->start[idx]; i < sh->start[idx + 1]; i++) {
                real d2 = dist_sq_c(p, seeds[sh->seeds[i]].center);
                if (d2 < d2_min) d2_min = d2;
            }
        }
    }
    return d2_min;
}

// =====================================================================
// Collect interface seeds from the local multiphase domain.
// Returns a malloc'd array and writes its length to *out_count.
// At step 0 the analytical fracvol is used; afterwards the distributed
// property is read.
// =====================================================================
static InterfaceSeedAdapt *collect_interface_seeds_local(
    higflow_solver *ns, int *out_count)
{
    sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
    mp_mapper *mp = sd_get_domain_mapper(sdm);

    int seed_cap = 1000;
    int seed_count = 0;
    InterfaceSeedAdapt *seeds = (InterfaceSeedAdapt *) malloc(seed_cap * sizeof(*seeds));

    higcit_celliterator *it = sd_get_domain_celliterator(sdm);
    while (!higcit_isfinished(it)) {
        hig_cell *c = higcit_getcell(it);
        int clid = mp_lookup(mp, hig_get_cid(c));
        if (clid < 0) {
            higcit_nextcell(it);
            continue;
        }

        real val;
        if (ns->par.step == 0) {
            Point xc, delta;
            hig_get_center(c, xc);
            hig_get_delta(c, delta);
            // DUAS APIs convivem: a antiga (ponteiros soltos, que o
            // example2d_DynamicMeshAdapt usa) e a de OBJETO DE PROBLEMA (que o VOF e
            // o front-tracking usam).  Quem registra pelo objeto deixa get_fracvol
            // NULO, e chama-lo da' SEGV -- foi o que aconteceu ao trazer o modulo.
            if (ns->ed.mult.problem != NULL)
                val = ns->ed.mult.problem->fracvol(xc, delta, ns->par.t);
            else if (ns->ed.mult.get_fracvol != NULL)
                val = ns->ed.mult.get_fracvol(xc, delta, ns->par.t);
            else
                val = 0.0;
        } else {
            val = dp_get_value(ns->ed.mult.dpfracvol, clid);
        }

        if (val > 0.001 && val < 0.999) {
            if (seed_count >= seed_cap) {
                seed_cap *= 2;
                seeds = (InterfaceSeedAdapt *) realloc(seeds, seed_cap * sizeof(*seeds));
            }
            hig_get_center(c, seeds[seed_count].center);
            seed_count++;
        }
        higcit_nextcell(it);
    }
    higcit_destroy(it);

    *out_count = seed_count;
    return seeds;
}

// =====================================================================
// Gather interface seeds from all MPI ranks into a single global array.
// Every rank ends up with the same seed set.  Caller frees the result.
// =====================================================================
static InterfaceSeedAdapt *gather_seeds_mpi(InterfaceSeedAdapt *local,
                                            int local_count,
                                            int *out_total)
{
    int ntasks;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    int *counts = (int *) malloc(ntasks * sizeof(int));
    int *disps  = (int *) malloc(ntasks * sizeof(int));
    MPI_Allgather(&local_count, 1, MPI_INT, counts, 1, MPI_INT,
                  MPI_COMM_WORLD);

    int total = 0;
    for (int r = 0; r < ntasks; r++) total += counts[r];
    for (int r = 0, d = 0; r < ntasks; r++) {
        disps[r] = d;
        d += counts[r];
    }

    InterfaceSeedAdapt *global = NULL;
    if (total > 0) {
        global = (InterfaceSeedAdapt *) malloc(total * sizeof(*global));
        int sz = (int)sizeof(InterfaceSeedAdapt);
        int *byte_counts = (int *) malloc(ntasks * sizeof(int));
        int *byte_disps  = (int *) malloc(ntasks * sizeof(int));
        for (int r = 0; r < ntasks; r++) {
            byte_counts[r] = counts[r] * sz;
            byte_disps[r]  = disps[r] * sz;
        }
        MPI_Allgatherv(local, local_count * sz, MPI_BYTE,
                       global, byte_counts, byte_disps, MPI_BYTE,
                       MPI_COMM_WORLD);
        free(byte_counts);
        free(byte_disps);
    }

    free(counts);
    free(disps);
    *out_total = total;
    return global;
}

// =====================================================================
// Adapt a tree given an explicit set of interface seeds.
// Refines cells near any seed and coarsens those far away.
// =====================================================================
static void adapt_tree_with_seeds(real *thresholds, hig_cell *root,
                                  InterfaceSeedAdapt *seeds, int seed_count)
{
    int num_levels = 0;
    real thr_sq[10], thr_coarse_sq[10];
    while (num_levels < 10 && thresholds[num_levels] >= 0) num_levels++;
    for (int i = 0; i < num_levels; i++) {
        thr_sq[i] = thresholds[i] * thresholds[i];
        thr_coarse_sq[i] = (thresholds[i] + 0.005) *
                           (thresholds[i] + 0.005);
    }
    real max_search = (num_levels > 0) ? thresholds[0] : 0.0;

    GHashTable *level_cache = g_hash_table_new(g_direct_hash,
                                                g_direct_equal);
    GHashTable *ref_set     = g_hash_table_new(g_direct_hash,
                                                g_direct_equal);

    Point bbox_lo = {DBL_MAX, DBL_MAX};
    Point bbox_hi = {-DBL_MAX, -DBL_MAX};
    for (int s = 0; s < seed_count; s++) {
        for (int d = 0; d < DIM; d++) {
            if (seeds[s].center[d] - max_search < bbox_lo[d])
                bbox_lo[d] = seeds[s].center[d] - max_search;
            if (seeds[s].center[d] + max_search > bbox_hi[d])
                bbox_hi[d] = seeds[s].center[d] + max_search;
        }
    }

    SeedHash seed_hash;
    if (seed_count > 0) {
        real bh = max_search + 0.005;
        build_seed_hash(&seed_hash, seeds, seed_count, bh, bh,
                        bbox_lo, bbox_hi);
    }

    int refine_cap = (seed_count > 0) ? seed_count * 4 : 1000;
    int merge_cap  = (seed_count > 0) ? seed_count * 2 : 1000;

    // STAGE A: REFINEMENT
    if (seed_count > 0 && num_levels > 0) {
        for (int pass = 0; pass < num_levels + 1; pass++) {
            int refine_count = 0;
            hig_cell **to_refine = (hig_cell * *) malloc(refine_cap * sizeof(hig_cell*));
            g_hash_table_remove_all(ref_set);

            higcit_celliterator *it_bb =
                higcit_create_bounding_box(root, bbox_lo, bbox_hi);
            while (!higcit_isfinished(it_bb)) {
                hig_cell *neigh = higcit_getcell(it_bb);
                if (hig_get_number_of_children(neigh) == 0) {
                    Point cn;
                    hig_get_center(neigh, cn);
                    real d2 = min_dist2_to_seeds(cn, seeds, &seed_hash);
                    int target_level = 0;
                    for (int l = num_levels - 1; l >= 0; l--) {
                        if (d2 <= thr_sq[l]) {
                            target_level = l + 1;
                            break;
                        }
                    }
                    if (target_level > 0) {
                        int cl = get_cached_level(level_cache, neigh);
                        if (cl < target_level &&
                            !g_hash_table_contains(ref_set, neigh)) {
                            if (refine_count >= refine_cap) {
                                refine_cap *= 2;
                                to_refine = (hig_cell * *) realloc(to_refine,
                                    refine_cap * sizeof(hig_cell*));
                            }
                            to_refine[refine_count++] = neigh;
                            g_hash_table_add(ref_set, neigh);
                        }
                    }
                }
                higcit_nextcell(it_bb);
            }
            higcit_destroy(it_bb);

            if (refine_count == 0) {
                free(to_refine);
                break;
            }
            for (int i = 0; i < refine_count; i++) {
                int nc[DIM] = {2, 2};
                hig_refine_uniform(to_refine[i], nc);
            }
            free(to_refine);
        }
    }

    // STAGE B: COARSENING
    if (seed_count > 0 && num_levels > 0) {
        for (int pass = 0; pass < 4; pass++) {
            int merge_count = 0;
            hig_cell **parents_to_merge =
                (hig_cell **) malloc(merge_cap * sizeof(hig_cell*));

            higcit_celliterator *it_all =
                higcit_create_all_higtree(root);
            while (!higcit_isfinished(it_all)) {
                hig_cell *c = higcit_getcell(it_all);
                int nc = hig_get_number_of_children(c);
                if (nc > 0) {
                    int all_leaves = 1;
                    for (int k = 0; k < nc; k++) {
                        if (hig_get_number_of_children(
                                hig_get_child(c, k)) > 0) {
                            all_leaves = 0;
                            break;
                        }
                    }
                    if (all_leaves) {
                        int plvl = get_cached_level(level_cache, c);
                        int safe = 1;
                        for (int k = 0; k < nc && safe; k++) {
                            hig_cell *ch = hig_get_child(c, k);
                            Point cch;
                            hig_get_center(ch, cch);
                            real d2 = min_dist2_to_seeds(cch, seeds,
                                                         &seed_hash);
                            real thr = (plvl < num_levels)
                                           ? thr_coarse_sq[plvl]
                                           : -1.0;
                            if (d2 <= thr) safe = 0;
                        }
                        if (safe) {
                            if (merge_count >= merge_cap) {
                                merge_cap *= 2;
                                parents_to_merge = (hig_cell * *) realloc(
                                    parents_to_merge,
                                    merge_cap * sizeof(hig_cell*));
                            }
                            parents_to_merge[merge_count++] = c;
                        }
                    }
                }
                higcit_nextcell(it_all);
            }
            higcit_destroy(it_all);

            if (merge_count == 0) {
                free(parents_to_merge);
                break;
            }
            for (int i = 0; i < merge_count; i++) {
                int nc = hig_get_number_of_children(parents_to_merge[i]);
                for (int k = 0; k < nc; k++)
                    g_hash_table_remove(level_cache,
                        hig_get_child(parents_to_merge[i], k));
                hig_merge_children(parents_to_merge[i]);
            }
            free(parents_to_merge);
        }
    }

    g_hash_table_destroy(level_cache);
    g_hash_table_destroy(ref_set);
    if (seed_count > 0) free_seed_hash(&seed_hash);
}

// =====================================================================
// Backwards-compatible wrapper: collect seeds from the local domain and
// adapt the given tree.  Used by the serial and in-place paths.
// =====================================================================
static void adapt_tree_core(higflow_solver *ns, real *thresholds,
                            hig_cell *root)
{
    int seed_count = 0;
    InterfaceSeedAdapt *seeds = collect_interface_seeds_local(ns,
                                                              &seed_count);
    adapt_tree_with_seeds(thresholds, root, seeds, seed_count);
    free(seeds);
}

// Build a fresh adapted tree (clone of the live mesh). Caller owns the result.
hig_cell *higflow_make_adapted_tree_params(higflow_solver *ns, real *thresholds)
{
    sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
    hig_cell *root_o = sd_get_higtree(sdm, 0);
    if (!root_o) return NULL;
    hig_cell *root_c = hig_clone(root_o);
    if (root_c) adapt_tree_core(ns, thresholds, root_c);
    return root_c;
}

// Adapt the live mesh in place.
void higflow_refine_tree_inplace(higflow_solver *ns, real *thresholds)
{
    sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
    hig_cell *root = sd_get_higtree(sdm, 0);
    if (root) adapt_tree_core(ns, thresholds, root);
}

// =====================================================================
// Build a globally adapted tree identical on all ranks.
//
// Reads the original AMR file (identical on every rank), gathers the
// interface seeds from all ranks, and adapts the tree with the same
// seed set everywhere.  The returned tree can be passed to the load
// balancer for parallel partitioning.
// =====================================================================
hig_cell *higflow_make_global_adapted_tree(higflow_solver *ns,
                                           real *thresholds,
                                           const char *amr_filename)
{
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    FILE *fd = fopen(amr_filename, "r");
    if (!fd) {
        if (myrank == 0)
            fprintf(stderr, "Cannot open AMR file %s\n", amr_filename);
        return NULL;
    }
    higio_amr_info *mi = higio_read_amr_info(fd);
    fclose(fd);
    if (!mi) return NULL;

    hig_cell *root = higio_read_from_amr_info(mi);
    higio_amr_info_destroy(mi);
    if (!root) return NULL;

    int local_count = 0;
    InterfaceSeedAdapt *local_seeds =
        collect_interface_seeds_local(ns, &local_count);

    int global_count = 0;
    InterfaceSeedAdapt *global_seeds =
        gather_seeds_mpi(local_seeds, local_count, &global_count);
    free(local_seeds);

    adapt_tree_with_seeds(thresholds, root, global_seeds, global_count);

    int leaves_after = 0;
    higcit_celliterator *it = higcit_create_all_leaves(root);
    for (; !higcit_isfinished(it); higcit_nextcell(it)) leaves_after++;
    higcit_destroy(it);

    print0f("===> AMR: seeds=%d leaves=%d\n",
            global_count, leaves_after);

    free(global_seeds);
    return root;
}

// =====================================================================
// BC higtree refinement for AMR cases.
// Registered via higflow_set_bc_refine_hook() so each boundary higtree
// is refined in place to match the adjacent internal mesh before the
// sim_boundary is created.  _bc_domain_root is set just before calling
// higflow_initialize_boundaries_yaml() and cleared immediately after.
// =====================================================================
static hig_cell *_bc_domain_root = NULL;

static int _bc_level(hig_cell *c) {
    int l = 0;
    while (c) { c = hig_get_parent(c); l++; }
    return l;
}

static void _refine_bc_tree(hig_cell *bc_root, int bc_id) {
    if (!_bc_domain_root) return;
    const real eps = 1e-7;
    bool changed = true;
    while (changed) {
        changed = false;
        higcit_celliterator *it = higcit_create_all_leaves(bc_root);
        for (; !higcit_isfinished(it); higcit_nextcell(it)) {
            hig_cell *bc_leaf = higcit_getcell(it);
            Point center; hig_get_center(bc_leaf, center);
            Point q; q[0] = center[0]; q[1] = center[1];
            if      (bc_id == 0) q[0] += eps;
            else if (bc_id == 1) q[1] -= eps;
            else if (bc_id == 2) q[0] -= eps;
            else                 q[1] += eps;
            hig_cell *dom = hig_get_cell_with_point(_bc_domain_root, q);
            if (!dom) continue;
            int bc_lev  = _bc_level(bc_leaf);
            int dom_lev = _bc_level(dom);
            if (dom_lev > bc_lev) {
                int nc[DIM];
                if (bc_id == 0 || bc_id == 2) { nc[0]=1; nc[1]=2; }
                else                            { nc[0]=2; nc[1]=1; }
                hig_refine_uniform(bc_leaf, nc);
                changed = true;
                higcit_destroy(it);
                break;
            }
        }
        if (!changed) higcit_destroy(it);
    }
}

// =====================================================================
// A REGRA DAS CELULAS MINIMAS, derivada -- e o oraculo que a confere.
// =====================================================================

//! Preenche `thr` com os limiares de distancia para `niveis` niveis de refino,
//! de modo que o nivel MAIS FINO cubra ao menos `cel_min` celulas finas de cada
//! lado da interface.  `thr` deve caber `niveis+1` reais (o ultimo e' o -1 que
//! encerra a lista).
//!
//! A margem nao e' enfeite: a semente e' o CENTRO de uma celula de interface, e
//! pode estar a ate' h_base/2 da interface de verdade; e a celula candidata e'
//! julgada pelo centro DELA, a ate' h_fino/2 da borda.  Sem somar as duas, a
//! banda medida sai menor que a pedida, e justamente nas celulas que importam.
void malha_adapt_limiares(real h_base, int niveis, int cel_min, real *thr)
{
    real h_fino = h_base;
    for (int l = 0; l < niveis; l++) h_fino *= 0.5;
    // MARGEM DE UM h_base INTEIRO, e a medida e' que decidiu isto.  Com
    // 0,5*h_base a banda medida dava 4,02 celulas finas onde a regra pedia 5.
    // A razao e' que o criterio decide por CELULA BASE: o teste e' no centro
    // dela (incerteza h_base/2) e a celula ainda se estende outro h_base/2 alem.
    const real margem = h_base + h_fino;
    const real banda_fina = cel_min * h_fino + margem;
    // Limiares DECRESCENTES: thr[niveis-1] e' o do nivel mais fino.  Os niveis
    // intermediarios recebem bandas progressivamente maiores, para a malha
    // crescer em degraus e nao num salto so'.
    for (int l = niveis - 1; l >= 0; l--) {
        int passos = niveis - 1 - l;                 // 0 no mais fino
        real h_nivel = h_base;
        for (int k = 0; k < l + 1; k++) h_nivel *= 0.5;
        // cada nivel mais grosso ganha a banda do anterior mais cel_min celulas
        // DELE, que e' o que mantem o degrau com largura suficiente.
        thr[l] = banda_fina;
        for (int k = 0; k < passos; k++) thr[l] += cel_min * h_nivel;
    }
    thr[niveis] = -1.0;
}

//! MEDE, na malha de verdade, a menor distancia entre uma semente de interface e
//! uma celula que NAO esteja no nivel mais fino -- e a devolve em unidades de
//! h_fino.  E' o oraculo da regra: se a banda fina tem ao menos `cel_min`
//! celulas de cada lado, esse numero e' >= cel_min.
//!
//! Devolve -1 se nao houver semente (sem interface no rank) ou se toda a malha
//! estiver no nivel fino (banda infinita, regra trivialmente satisfeita).
real malha_adapt_mede_banda(higflow_solver *ns, int niveis, real h_base)
{
    int n_sem = 0;
    InterfaceSeedAdapt *sem = collect_interface_seeds_local(ns, &n_sem);
    if (getenv("FT_DIAG_MALHA") != NULL)
        fprintf(stderr, "  [banda] sementes=%d passo=%d\n", n_sem, ns->par.step);
    if (n_sem == 0) { free(sem); return -1.0; }

    real h_fino = h_base;
    for (int l = 0; l < niveis; l++) h_fino *= 0.5;

    // O DOMINIO CERTO E' O DA PRESSAO, e essa distincao custou uma medida vazia.
    // O sdmult do HiGFlow e' UNIFORME no nivel mais fino POR DESENHO
    // (higflow_create_amr_info_mult: "this uniform mesh has the same size as the
    // last level of the original mesh"), entao la' nao existe celula grossa e a
    // busca devolvia -1 sempre -- um oraculo que nunca podia falhar, que e' o
    // mesmo que nao ter oraculo.  As SEMENTES continuam vindo do sdmult, onde o
    // fracvol vive; o tamanho de celula e' julgado no dominio ADAPTADO.
    sim_domain *sdm = psd_get_local_domain(ns->psdp);
    real d2min = 1e300, dmin_v = 1e300, dmax_v = 0.0;
    long ncel_v = 0, ngrossa_v = 0;
    Point pior_c = {0,0}; real pior_d = 0.0;
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(sdm); !higcit_isfinished(it);
         higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        Point d; hig_get_delta(c, d);
        if (d[0] < dmin_v) dmin_v = d[0];
        if (d[0] > dmax_v) dmax_v = d[0];
        ncel_v++;
        // "nao esta' no nivel mais fino" = celula maior que h_fino com folga
        if (d[0] < 1.5 * h_fino) continue;
        ngrossa_v++;
        Point cc; hig_get_center(c, cc);
        real dd = 0.0;
        for (int k = 0; k < n_sem; k++) {
            real s = 0.0;
            for (int i = 0; i < DIM; i++) {
                real t = cc[i] - sem[k].center[i];
                s += t * t;
            }
            if (k == 0 || s < dd) dd = s;
        }
        if (dd < d2min) { d2min = dd; pior_c[0]=cc[0]; pior_c[1]=cc[1]; pior_d=d[0]; }
    }
    higcit_destroy(it);
    free(sem);
    if (getenv("FT_DIAG_MALHA") != NULL)
        fprintf(stderr, "  [banda] celulas=%ld delta[%.5f,%.5f] h_fino=%.5f "
                "grossas=%ld d2min=%.4e  pior em (%.4f,%.4f) delta=%.5f\n", ncel_v,
                (double) dmin_v, (double) dmax_v, (double) h_fino, ngrossa_v,
                (double) d2min, (double) pior_c[0], (double) pior_c[1],
                (double) pior_d);
    if (d2min > 1e299) return -1.0;
    // O QUE DESCONTAR, e eu errei isto na primeira versao.  A distancia medida
    // e' de CENTRO a CENTRO.  Para saber ate' onde a regiao fina se estende sao
    // duas correcoes, e nenhuma e' h_base:
    //   - a celula grossa mais proxima: sua BORDA esta' a meia-largura DELA do
    //     centro (pior_d/2), nao h_base/2 -- ela pode ser de nivel intermediario;
    //   - a semente: e' o centro de uma celula do dominio da fracao, que e'
    //     uniforme no nivel fino, entao a interface real esta' a ate' h_fino/2.
    // Subtrair h_base/2 (o que eu fazia) penalizava a medida em quase o dobro:
    // relatava 3,52 celulas onde a banda de verdade tinha ~4,9.
    real dist = sqrt(d2min) - 0.5 * pior_d - 0.5 * h_fino;
    if (dist < 0.0) dist = 0.0;
    return dist / h_fino;
}


// =====================================================================
// Reconstrucao do dominio com a malha adaptada (compartilhada).
// =====================================================================
//! Reconstroi o dominio com a malha adaptada.  Parametrizado pelo caminho do
//! .amr base e pelos limiares, para servir aos DOIS exemplos bifasicos -- VOF e
//! front-tracking -- com o MESMO criterio.  Como o criterio le' `fracvol`, e o
//! front-tracking agora ESCREVE o fracvol da geometria da frente, os dois
//! adaptam identicamente e a comparacao isola a representacao da interface.
void malha_adapt_reconstroi(higflow_solver *ns, higflow_solver *ns2,
                            int myrank, int ntasks, int cache, int order_center,
                            real *limiares, const char *amr_base)
{
    (void)cache;
    (void)ntasks;

    // Global adapted tree identical on every rank.
    hig_cell *root = higflow_make_global_adapted_tree(
        ns, limiares, amr_base);

    // Sanity check: every rank must see the same adapted tree.
    long local_leaves = root ? hig_get_number_of_leaves(root) : 0;
    long global_leaves;
    MPI_Allreduce(&local_leaves, &global_leaves, 1, MPI_LONG, MPI_SUM,
                  MPI_COMM_WORLD);
    print0f("===> AMR: rank %d leaves=%ld total_leaves=%ld\n",
            myrank, local_leaves, global_leaves);

    // Uniform multiphase tree built from the same AMR file.
    FILE *fd = fopen(amr_base, "r");
    higio_amr_info *mi = fd ? higio_read_amr_info(fd) : NULL;
    if (fd) fclose(fd);
    higio_amr_info *mi_mult = mi ? higflow_create_amr_info_mult(mi) : NULL;
    hig_cell *mult_root = mi_mult ? higio_read_from_amr_info(mi_mult) : NULL;

    // Partition graph and load balancer (group 0 = flow, group 1 = multiphase).
    print0f("===> AMR: rank %d creating partition graph\n", myrank);
    partition_graph *pg = pg_create(MPI_COMM_WORLD);
    pg_set_fringe_size(pg, 5);
    load_balancer *lb = lb_create(MPI_COMM_WORLD, 2);
    if (myrank == 0 && root)
        lb_add_input_tree(lb, root, true, 0);
    if (myrank == 0 && mult_root)
        lb_add_input_tree(lb, mult_root, true, 1);
    print0f("===> AMR: rank %d calling lb_calc_partition\n", myrank);
    lb_calc_partition(lb, pg);
    print0f("===> AMR: rank %d lb_calc_partition done\n", myrank);

    // Add the local trees returned by the LB to ns2's domains.
    print0f("===> AMR: rank %d adding local trees to domains\n", myrank);
    int numhigs0 = lb_get_num_local_trees_in_group(lb, 0);
    for (int h = 0; h < numhigs0; h++) {
        hig_cell *t = lb_get_local_tree_in_group(lb, h, 0);
        sd_add_higtree(ns2->sdp, t);
        sd_add_higtree(ns2->sdF, t);
        if (ns2->ed.mult.contr.viscoelastic_either == true)
            sd_add_higtree(ns2->ed.sdED, t);
    }
    int numhigs1 = lb_get_num_local_trees_in_group(lb, 1);
    for (int h = 0; h < numhigs1; h++) {
        hig_cell *t = lb_get_local_tree_in_group(lb, h, 1);
        sd_add_higtree(ns2->ed.mult.sdmult, t);
    }
    lb_destroy(lb);

    // Partitioned domains and stencils.
    print0f("===> AMR: rank %d creating partitioned domains\n", myrank);
    higflow_create_partitioned_domain(ns2, pg, order_center);
    print0f("===> AMR: rank %d creating multiphase domain\n", myrank);
    higflow_create_partitioned_domain_multiphase(ns2, pg, order_center);
    print0f("===> AMR: rank %d creating stencils\n", myrank);
    higflow_create_stencil(ns2);
    higflow_create_stencil_multiphase(ns2);
    if (ns2->ed.mult.contr.viscoelastic_either == true)
        higflow_create_stencil_for_extra_domain(ns2);

    // Distributed properties.
    print0f("===> AMR: rank %d creating distributed properties\n", myrank);
    higflow_create_distributed_properties(ns2);

    // Boundaries refined to match the adapted internal mesh.
    print0f("===> AMR: rank %d initializing boundaries\n", myrank);
    _bc_domain_root = sd_get_higtree(psd_get_local_domain(ns2->psdp), 0);
    higflow_set_bc_refine_hook(_refine_bc_tree);
    higflow_initialize_boundaries_yaml(ns2);
    higflow_set_bc_refine_hook(NULL);
    _bc_domain_root = NULL;

    // Solver creation is intentionally left to the caller.  Creating it
    // here and again in main overwrites the solver handle and leaves the
    // distributed properties mapped to the discarded solver.
    print0f("===> AMR: rank %d creating solver\n", myrank);
    higflow_create_solver(ns2);
    print0f("===> AMR: rank %d solver created\n", myrank);
}
