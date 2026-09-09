
#include "higtree.h"
#include "higtree-io.h"
#include "higtree-iterator.h"
#include "pdomain.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <glib.h>
#include <float.h>

/*! Seed for interface refinement — stores cell centre */
typedef struct { Point center; } InterfaceSeedAdapt;

static real dist_sq_c(Point p1, Point p2) {
    real d = 0.0;
    for (int i = 0; i < DIM; i++) d += (p1[i]-p2[i])*(p1[i]-p2[i]);
    return d;
}

/*! Cell-level cache (GHashTable) for O(1) lookup. */
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

#if DIM == 2
/* ------------------------------------------------------------------ */
/*  2-D spatial hash                                                  */
/* ------------------------------------------------------------------ */
typedef struct {
    int nx, ny;   real bin_w, bin_h, ox, oy;
    int *start;   int *seeds;
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
    int *cnt = calloc(nb, sizeof(int));
    for (int s = 0; s < n; s++) {
        int bx = (int)((seeds[s].center[0] - lo[0]) / bin_w);
        int by = (int)((seeds[s].center[1] - lo[1]) / bin_h);
        if (bx < 0) bx = 0; if (bx >= sh->nx) bx = sh->nx - 1;
        if (by < 0) by = 0; if (by >= sh->ny) by = sh->ny - 1;
        cnt[by * sh->nx + bx]++;
    }
    sh->start = malloc((nb + 1) * sizeof(int));
    int offset = 0;
    for (int i = 0; i < nb; i++) { sh->start[i] = offset; offset += cnt[i]; }
    sh->start[nb] = offset;
    sh->seeds = malloc(n * sizeof(int));
    int *pos = malloc(nb * sizeof(int));
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

static real min_dist2_to_seeds(Point p, InterfaceSeedAdapt *seeds,
                                SeedHash *sh)
{
    int bx = (int)((p[0] - sh->ox) / sh->bin_w);
    int by = (int)((p[1] - sh->oy) / sh->bin_h);
    real d2_min = DBL_MAX;
    for (int dy = -1; dy <= 1; dy++)
        for (int dx = -1; dx <= 1; dx++) {
            int cx = bx + dx, cy = by + dy;
            if (cx < 0 || cx >= sh->nx || cy < 0 || cy >= sh->ny) continue;
            int idx = cy * sh->nx + cx;
            for (int i = sh->start[idx]; i < sh->start[idx + 1]; i++) {
                real d2 = dist_sq_c(p, seeds[sh->seeds[i]].center);
                if (d2 < d2_min) d2_min = d2;
            }
        }
    return d2_min;
}

#elif DIM == 3
/* ------------------------------------------------------------------ */
/*  3-D spatial hash                                                  */
/* ------------------------------------------------------------------ */
typedef struct {
    int nx, ny, nz;   real bin_w, bin_h, bin_d, ox, oy, oz;
    int *start;       int *seeds;
} SeedHash;

static void build_seed_hash(SeedHash *sh, InterfaceSeedAdapt *seeds, int n,
                             real bin_w, real bin_h, real bin_d,
                             Point lo, Point hi)
{
    real sx = hi[0] - lo[0], sy = hi[1] - lo[1], sz = hi[2] - lo[2];
    sh->nx = (int)ceil(sx / bin_w); if (sh->nx < 1) sh->nx = 1;
    sh->ny = (int)ceil(sy / bin_h); if (sh->ny < 1) sh->ny = 1;
    sh->nz = (int)ceil(sz / bin_d); if (sh->nz < 1) sh->nz = 1;
    sh->bin_w = bin_w; sh->bin_h = bin_h; sh->bin_d = bin_d;
    sh->ox = lo[0]; sh->oy = lo[1]; sh->oz = lo[2];

    int nb = sh->nx * sh->ny * sh->nz;
    int *cnt = calloc(nb, sizeof(int));
    for (int s = 0; s < n; s++) {
        int bx = (int)((seeds[s].center[0] - lo[0]) / bin_w);
        int by = (int)((seeds[s].center[1] - lo[1]) / bin_h);
        int bz = (int)((seeds[s].center[2] - lo[2]) / bin_d);
        if (bx < 0) bx = 0; if (bx >= sh->nx) bx = sh->nx - 1;
        if (by < 0) by = 0; if (by >= sh->ny) by = sh->ny - 1;
        if (bz < 0) bz = 0; if (bz >= sh->nz) bz = sh->nz - 1;
        cnt[(bz * sh->ny + by) * sh->nx + bx]++;
    }
    sh->start = malloc((nb + 1) * sizeof(int));
    int offset = 0;
    for (int i = 0; i < nb; i++) { sh->start[i] = offset; offset += cnt[i]; }
    sh->start[nb] = offset;
    sh->seeds = malloc(n * sizeof(int));
    int *pos = malloc(nb * sizeof(int));
    memcpy(pos, sh->start, nb * sizeof(int));
    for (int s = 0; s < n; s++) {
        int bx = (int)((seeds[s].center[0] - lo[0]) / bin_w);
        int by = (int)((seeds[s].center[1] - lo[1]) / bin_h);
        int bz = (int)((seeds[s].center[2] - lo[2]) / bin_d);
        if (bx < 0) bx = 0; if (bx >= sh->nx) bx = sh->nx - 1;
        if (by < 0) by = 0; if (by >= sh->ny) by = sh->ny - 1;
        if (bz < 0) bz = 0; if (bz >= sh->nz) bz = sh->nz - 1;
        sh->seeds[pos[(bz * sh->ny + by) * sh->nx + bx]++] = s;
    }
    free(pos); free(cnt);
}

static void free_seed_hash(SeedHash *sh) {
    free(sh->start); free(sh->seeds); memset(sh, 0, sizeof(*sh));
}

static real min_dist2_to_seeds(Point p, InterfaceSeedAdapt *seeds,
                                SeedHash *sh)
{
    int bx = (int)((p[0] - sh->ox) / sh->bin_w);
    int by = (int)((p[1] - sh->oy) / sh->bin_h);
    int bz = (int)((p[2] - sh->oz) / sh->bin_d);
    real d2_min = DBL_MAX;
    for (int dz = -1; dz <= 1; dz++)
        for (int dy = -1; dy <= 1; dy++)
            for (int dx = -1; dx <= 1; dx++) {
                int cx = bx + dx, cy = by + dy, cz = bz + dz;
                if (cx < 0 || cx >= sh->nx) continue;
                if (cy < 0 || cy >= sh->ny) continue;
                if (cz < 0 || cz >= sh->nz) continue;
                int idx = ((cz * sh->ny + cy) * sh->nx) + cx;
                for (int i = sh->start[idx]; i < sh->start[idx + 1]; i++) {
                    real d2 = dist_sq_c(p, seeds[sh->seeds[i]].center);
                    if (d2 < d2_min) d2_min = d2;
                }
            }
    return d2_min;
}
#endif  // DIM == 2 / 3

/* =====================================================================
 *  Mesh-adaptation driver
 *
 *  1. Collect interface seeds (cells with 0.001 < fracvol < 0.999)
 *  2. Build a spatial hash for O(1) neighbour queries
 *  3. Refine pass: cells within thresholds distance → subdivide
 *  4. Coarsening pass: cells far from interface → merge children
 * ===================================================================== */
hig_cell *higflow_make_adapted_tree_params(higflow_solver *ns, real *thresholds)
{
    int num_levels = 0;
    real thr_sq[10], thr_coarse_sq[10];
    while (num_levels < 10 && thresholds[num_levels] >= 0) num_levels++;
    for (int i = 0; i < num_levels; i++) {
        thr_sq[i] = thresholds[i] * thresholds[i];
        thr_coarse_sq[i] = (thresholds[i] + 0.005) * (thresholds[i] + 0.005);
    }
    real max_search = (num_levels > 0) ? thresholds[0] : 0.0;

    sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
    mp_mapper *mp = sd_get_domain_mapper(sdm);
    int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    hig_cell *root_o = sd_get_higtree(sdm, 0);
    if (!root_o) return NULL;
    hig_cell *root_c = hig_clone(root_o);
    if (!root_c) return NULL;

    /* ---- 1. Collect seeds ---- */
    int seed_cap = 1000, seed_count = 0;
    InterfaceSeedAdapt *seeds = malloc(seed_cap * sizeof(*seeds));
    higcit_celliterator *it = sd_get_domain_celliterator(sdm);
    while (!higcit_isfinished(it)) {
        hig_cell *c = higcit_getcell(it);
        int clid = mp_lookup(mp, hig_get_cid(c));
        if (clid >= 0) {
            real val;
            if (ns->par.step == 0) {
                Point xc, delta;
                hig_get_center(c, xc); hig_get_delta(c, delta);
                val = ns->ed.mult.get_fracvol(xc, delta, ns->par.t);
            } else {
                val = dp_get_value(ns->ed.mult.dpfracvol, clid);
            }
            if (val > 0.001 && val < 0.999) {
                if (seed_count >= seed_cap) {
                    seed_cap *= 2;
                    seeds = realloc(seeds, seed_cap * sizeof(*seeds));
                }
                hig_get_center(c, seeds[seed_count].center);
                seed_count++;
            }
        }
        higcit_nextcell(it);
    }
    higcit_destroy(it);

    /* ---- 2. Build level cache and refinement set ---- */
    GHashTable *level_cache = g_hash_table_new(g_direct_hash, g_direct_equal);
    GHashTable *ref_set     = g_hash_table_new(g_direct_hash, g_direct_equal);

    /* ---- Bounding box around all seeds ---- */
    Point bbox_lo, bbox_hi;
    for (int d = 0; d < DIM; d++) { bbox_lo[d] = DBL_MAX; bbox_hi[d] = -DBL_MAX; }
    for (int s = 0; s < seed_count; s++) {
        for (int d = 0; d < DIM; d++) {
            if (seeds[s].center[d] - max_search < bbox_lo[d])
                bbox_lo[d] = seeds[s].center[d] - max_search;
            if (seeds[s].center[d] + max_search > bbox_hi[d])
                bbox_hi[d] = seeds[s].center[d] + max_search;
        }
    }

    /* ---- Spatial hash ---- */
    SeedHash seed_hash;
    if (seed_count > 0) {
        real bin_size = max_search + 0.005;
#if DIM == 2
        build_seed_hash(&seed_hash, seeds, seed_count,
                        bin_size, bin_size, bbox_lo, bbox_hi);
#elif DIM == 3
        build_seed_hash(&seed_hash, seeds, seed_count,
                        bin_size, bin_size, bin_size, bbox_lo, bbox_hi);
#endif
    }

    int refine_cap = (seed_count > 0) ? seed_count * 4 : 1000;
    int merge_cap  = (seed_count > 0) ? seed_count * 2 : 1000;

    /* ==================================================================
     *  STEP A: REFINEMENT
     * ================================================================== */
    if (seed_count > 0 && num_levels > 0) {
        for (int pass = 0; pass < num_levels + 1; pass++) {
            int refine_count = 0;
            hig_cell **to_refine = malloc(refine_cap * sizeof(hig_cell*));
            g_hash_table_remove_all(ref_set);

            higcit_celliterator *it_bb =
                higcit_create_bounding_box(root_c, bbox_lo, bbox_hi);
            while (!higcit_isfinished(it_bb)) {
                hig_cell *neigh = higcit_getcell(it_bb);
                if (hig_get_number_of_children(neigh) == 0) {
                    Point cn; hig_get_center(neigh, cn);
                    real d2 = min_dist2_to_seeds(cn, seeds, &seed_hash);
                    int target_level = 0;
                    for (int l = num_levels - 1; l >= 0; l--) {
                        if (d2 <= thr_sq[l] + 1e-14) { target_level = l + 1; break; }
                    }
                    int cl = get_cached_level(level_cache, neigh);
                    if (target_level > cl) {
                        if (g_hash_table_contains(ref_set, neigh)) continue;
                        if (refine_count >= refine_cap) {
                            refine_cap *= 2;
                            to_refine = realloc(to_refine, refine_cap * sizeof(hig_cell*));
                        }
                        to_refine[refine_count++] = neigh;
                        g_hash_table_insert(ref_set, neigh, GINT_TO_POINTER(1));
                    }
                }
                higcit_nextcell(it_bb);
            }
            higcit_destroy(it_bb);

            for (int i = 0; i < refine_count; i++) {
                int nc[DIM];
                for (int d = 0; d < DIM; d++) nc[d] = 2;
                hig_refine_uniform(to_refine[i], nc);
            }
            free(to_refine);
        }
    }

    /* ==================================================================
     *  STEP B: COARSENING
     * ================================================================== */
    for (int pass = 0; pass < 4; pass++) {
        int merge_count = 0;
        hig_cell **to_merge = malloc(merge_cap * sizeof(hig_cell*));
        g_hash_table_remove_all(ref_set);

        higcit_celliterator *it_bb2 =
            higcit_create_bounding_box(root_c, bbox_lo, bbox_hi);
        while (!higcit_isfinished(it_bb2)) {
            hig_cell *c = higcit_getcell(it_bb2);
            if (hig_get_number_of_children(c) > 0) {
                int all_leaves = 1;
                hig_cell *ch = hig_get_child(c, 0);
                if (ch == NULL) all_leaves = 0;
                for (int ci = 1; ci < 4; ci++) {
                    if (hig_get_child(c, ci) == NULL) { all_leaves = 0; break; }
                }
                if (!all_leaves) { higcit_nextcell(it_bb2); continue; }

                Point cc; hig_get_center(c, cc);
                real d2 = min_dist2_to_seeds(cc, seeds, &seed_hash);
                int coarse_level = 0;
                for (int l = 0; l < num_levels; l++) {
                    if (d2 <= thr_coarse_sq[l] + 1e-14) { coarse_level = l; break; }
                }
                int cl = get_cached_level(level_cache, c);
                if (cl > 0 && cl > coarse_level) {
                    if (g_hash_table_contains(ref_set, c)) continue;
                    if (merge_count >= merge_cap) {
                        merge_cap *= 2;
                        to_merge = realloc(to_merge, merge_cap * sizeof(hig_cell*));
                    }
                    to_merge[merge_count++] = c;
                    g_hash_table_insert(ref_set, c, GINT_TO_POINTER(1));
                }
            }
            higcit_nextcell(it_bb2);
        }
        higcit_destroy(it_bb2);

        for (int i = 0; i < merge_count; i++) {
            hig_merge_children(to_merge[i]);
        }
        free(to_merge);
    }

    /* ---- Cleanup ---- */
    if (seed_count > 0) free_seed_hash(&seed_hash);
    g_hash_table_destroy(level_cache);
    g_hash_table_destroy(ref_set);
    free(seeds);

    return root_c;
}
