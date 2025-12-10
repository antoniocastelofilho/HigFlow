
#include "higtree.h"
#include "higtree-io.h"
#include "higtree-iterator.h"
#include "pdomain.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Estrutura para sementes da interface
typedef struct {
    Point center;
} InterfaceSeedAdapt;

// Função auxiliar para nível
int get_cell_level_c(hig_cell *c) {
    int level = 0;
    hig_cell *p = hig_get_parent(c);
    while (p != NULL) {
        level++;
        p = hig_get_parent(p);
    }
    return level - 1;
}

// Distância Euclidiana
real dist_sq_c(Point p1, Point p2) {
    real d = 0.0;
    for(int i=0; i<DIM; i++) d += (p1[i]-p2[i])*(p1[i]-p2[i]);
    return d;
}

// =========================================================================================
// FUNÇÃO DE ADAPTAÇÃO PARAMETRIZADA
// Recebe lista de distancias: thresholds[0] -> Level 1, thresholds[1] -> Level 2, ...
// =========================================================================================
hig_cell * higflow_make_adapted_tree_params(higflow_solver *ns, int num_levels, real *thresholds) {
    sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
    mp_mapper *mp = sd_get_domain_mapper(sdm);
    
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // 1. Clonagem (Adaptação na cópia)
    hig_cell *root_original = sd_get_higtree(sdm, 0);
    if (!root_original) return NULL;
    hig_cell *root_copy = hig_clone(root_original);
    if (!root_copy) return NULL;

    // Coletar Sementes (Interface)
    int seed_cap = 1000;
    int seed_count = 0;
    InterfaceSeedAdapt *seeds = (InterfaceSeedAdapt *)malloc(seed_cap * sizeof(InterfaceSeedAdapt));
    
    higcit_celliterator *it = sd_get_domain_celliterator(sdm);
    while (!higcit_isfinished(it)) {
        hig_cell *c = higcit_getcell(it);
        int clid = mp_lookup(mp, hig_get_cid(c));
        if (clid >= 0) {
            real val = dp_get_value(ns->ed.mult.dpfracvol, clid);
            if (val > 0.001 && val < 0.999) {
                if (seed_count >= seed_cap) {
                    seed_cap *= 2;
                    seeds = (InterfaceSeedAdapt *)realloc(seeds, seed_cap * sizeof(InterfaceSeedAdapt));
                }
                hig_get_center(c, seeds[seed_count].center);
                seed_count++;
            }
        }
        higcit_nextcell(it);
    }
    higcit_destroy(it);

    // =====================================================================
    // ETAPA A: REFINAMENTO
    // =====================================================================
    
    if (seed_count > 0) {
        // Iterate up to num_levels + 1 passes to ensure propagation
        for (int pass = 0; pass < num_levels + 1; pass++) {
            int refine_cap = 1000;
            int refine_count = 0;
            hig_cell **to_refine = (hig_cell **)malloc(refine_cap * sizeof(hig_cell*));

            // Max search distance is the largest threshold (level 1)
            real max_search = (num_levels > 0) ? thresholds[0] : 0.0;
            
            for(int s=0; s<seed_count; s++) {
                Point center_s;
                POINT_ASSIGN(center_s, seeds[s].center);
                Point lo, hi;
                for(int d=0; d<DIM; d++) {
                    lo[d] = center_s[d] - max_search;
                    hi[d] = center_s[d] + max_search;
                }
                
                higcit_celliterator *it_bb = higcit_create_bounding_box(root_copy, lo, hi);
                while(!higcit_isfinished(it_bb)) {
                    hig_cell *neigh = higcit_getcell(it_bb);
                    
                    if(hig_get_number_of_children(neigh) == 0) { // Folha
                        Point center_n;
                        hig_get_center(neigh, center_n);
                        real dist = sqrt(dist_sq_c(center_s, center_n));
                        
                        int target_level = 0;
                        // Check thresholds from highest level down to 1
                        for (int l = num_levels - 1; l >= 0; l--) {
                            if (dist <= thresholds[l]) {
                                target_level = l + 1;
                                break;
                            }
                        }
                        
                        int current_lvl = get_cell_level_c(neigh);
                        
                        if (current_lvl < target_level) {
                            if (refine_count >= refine_cap) {
                                refine_cap *= 2;
                                to_refine = (hig_cell **)realloc(to_refine, refine_cap * sizeof(hig_cell*));
                            }
                            int added = 0;
                            for(int k=0; k<refine_count; k++) if(to_refine[k] == neigh) { added=1; break; }
                            
                            if(!added) {
                                to_refine[refine_count++] = neigh;
                            }
                        }
                    }
                    higcit_nextcell(it_bb);
                }
                higcit_destroy(it_bb);
            }
            
            if(refine_count == 0 && pass > 0) { 
                 if (pass > num_levels) break; 
            }
            
               for(int i=0; i<refine_count; i++) {
                if(hig_get_number_of_children(to_refine[i]) == 0) {
                     int nc[DIM];
                     for(int d=0; d<DIM; d++) nc[d] = 2;
                     hig_refine_uniform(to_refine[i], nc);
                }
            }
            free(to_refine);
        }
    }

    // =====================================================================
    // ETAPA B: ENGROSSAMENTO (4 Passes)
    // =====================================================================
    
    float coarsen_hys = 0.005; // Margem de segurança
    
    for (int pass = 0; pass < 4; pass++) { 
        int merge_cap = 1000;
        int merge_count = 0;
        hig_cell **parents_to_merge = (hig_cell **)malloc(merge_cap * sizeof(hig_cell*));
        
        higcit_celliterator *it_all = higcit_create_all_higtree(root_copy);
        while (!higcit_isfinished(it_all)) {
            hig_cell *c = higcit_getcell(it_all);
            
            if (hig_get_number_of_children(c) > 0) {
                int all_leaves = 1;
                for(int k=0; k<hig_get_number_of_children(c); k++) {
                    hig_cell *ch = hig_get_child(c, k);
                    if(hig_get_number_of_children(ch) > 0) {
                        all_leaves = 0;
                        break;
                    }
                }
                
                if (all_leaves) {
                    int parent_lvl = get_cell_level_c(c); 
                    
                    int safe_to_merge = 1;
                    
                    if (seed_count > 0) {
                        for(int k=0; k<hig_get_number_of_children(c); k++) {
                            hig_cell *ch = hig_get_child(c, k);
                            Point center;
                            hig_get_center(ch, center);
                            
                            real min_d2 = 1.0e20;
                            for(int s=0; s<seed_count; s++) {
                                real d2 = dist_sq_c(center, seeds[s].center);
                                if (d2 < min_d2) min_d2 = d2;
                            }
                            real dist = sqrt(min_d2);
                            
                            real threshold = 0.0;
                            // logic: if parent_lvl is K, children are at K+1.
                            // To keep children at K+1, dist must be <= thresholds[K].
                            // To MERGE (go back to K), dist must be > thresholds[K] + hys.
                            // BUT wait, thresholds array is 0-indexed: index 0 is Level 1.
                            // If parent is Level 0, it has children at Level 1. threshold is thresholds[0].
                            // If parent is Level 1, it has children at Level 2. threshold is thresholds[1].
                            
                            if (parent_lvl < num_levels) {
                                threshold = thresholds[parent_lvl] + coarsen_hys;
                            } else {
                                // parent level >= num_levels? Should not happen if we only refine up to num_levels.
                                // But if it does, we should probably merge.
                                threshold = -1.0; // Always merge? Or huge?
                                // If threshold is small, we merge if dist > threshold.
                                // If we are WAY deep, we want to merge. 
                                // Set threshold to something that allows merge.
                                threshold = 1.0e20; // Wait, if dist <= threshold we DON'T merge.
                                // If we want to merge, make condition fail. 
                                // dist <= huge is TRUE -> safe_to_merge = 0 (Don't merge).
                                // Wait, logic:
                                // if (dist <= threshold) safe_to_merge = 0;
                                // We want to merge if we are too deep.
                                // Actually, if parent_lvl >= num_levels, we definitely want to merge to reduce level.
                                // So we want safe_to_merge = 1.
                                // So ensure (dist <= threshold) is FALSE.
                                threshold = -1.0; 
                            }
                            
                            if (dist <= threshold) {
                                safe_to_merge = 0; 
                                break;
                            }
                        }
                    }
                    
                    if (safe_to_merge) {
                        if (merge_count >= merge_cap) {
                            merge_cap *= 2;
                            parents_to_merge = (hig_cell **)realloc(parents_to_merge, merge_cap * sizeof(hig_cell*));
                        }
                        parents_to_merge[merge_count++] = c;
                    }
                }
            }
            higcit_nextcell(it_all);
        }
        higcit_destroy(it_all);
        
        if(merge_count > 0) {
            for(int i=0; i<merge_count; i++) {
                hig_merge_children(parents_to_merge[i]);
            }
        } else {
            break;
        }
        free(parents_to_merge);
    }
    
    free(seeds);
    return root_copy;
}
