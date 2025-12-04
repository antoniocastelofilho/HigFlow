// Função auxiliar para pegar o nível da célula
int get_cell_level(hig_cell *c) {
    int level = 0;
    hig_cell *p = hig_get_parent(c);
    while (p != NULL) {
        level++;
        p = hig_get_parent(p);
    }
    return level;
}

#define MAX_REF_LEVEL 3

#include "higtree.h"
#include "higtree-io.h"
#include "higtree-iterator.h"
#include "pdomain.h"
#include <math.h>

// Estrutura para armazenar posições da interface
typedef struct {
    Point center;
} InterfaceSeed;

// Função auxiliar para calcular o tamanho da menor célula (Nível Máximo)
real get_h_min(hig_cell *root) {
    Point delta;
    hig_get_delta(root, delta);
    real h = delta[0];
    return h;
    // for (int d = 1; d < DIM; d++) if (delta[d] < h) h = delta[d];
    // Divide pelo nível máximo (2^MAX_LEVEL)
    // for(int i = 0; i < max_level; i++) h *= 0.5;
}

// Função auxiliar para calcular distância ao quadrado entre pontos
real dist_sq(Point p1, Point p2) {
    real d = 0.0;
    for(int i=0; i<DIM; i++) d += (p1[i]-p2[i])*(p1[i]-p2[i]);
    return d;
}

void higflow_save_refined_mesh_preview(higflow_solver *ns, int frame_id) {
    sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
    mp_mapper *mp = sd_get_domain_mapper(sdm);
    
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // 1. Clonagem
    hig_cell *root_original = sd_get_higtree(sdm, 0);
    if (root_original == NULL) return;
    
    hig_cell *root_copy = hig_clone(root_original);
    if (root_copy == NULL) return;

    if(myrank == 0) printf("Preview: Refinando malha (DIM=%d) com Buffer e Grading (Frame %d)...\n", DIM, frame_id);

    
    float search_dist[MAX_REF_LEVEL]; // Tamanho fixo, mas usando a constante em cálculos
    search_dist[0] = 0.1;
    search_dist[1] = 0.15;
    search_dist[2] = 0.2;

    float refiment_levels[MAX_REF_LEVEL];
    for (int i = 0; i < MAX_REF_LEVEL; i++) refiment_levels[i] = 1;

    // Loop Iterativo de Propagação
    for (int pass = 0; pass < MAX_REF_LEVEL; pass++) {
        int capacity = 1000;
        int count = 0;
        hig_cell **cells_to_refine = (hig_cell **)malloc(capacity * sizeof(hig_cell*));


        higcit_celliterator *it = higcit_create_all_higtree(root_copy);
        while (!higcit_isfinished(it)) {
            hig_cell *c_copy = higcit_getcell(it);
            int childrens = hig_get_number_of_children(c_copy);
            
            // if (get_cell_level(c_copy) < MAX_REF_LEVEL) {
            if (childrens < pow(2, refiment_levels[pass])) {
                Point center;
                Point size;
                hig_get_center(c_copy, center);
                hig_get_delta(c_copy, size);

                int should_refine = 0;
                real val_center = -1.0;

                // 1. Valor no Centro
                hig_cell *c_orig = sd_get_cell_with_point(sdm, center);
                if (c_orig) {
                    int clid = mp_lookup(mp, hig_get_cid(c_orig));
                    if (clid >= 0) val_center = dp_get_value(ns->ed.mult.dpfracvol, clid);
                }

                if (val_center > 0.001 && val_center < 0.999) {
                    should_refine = 1;
                } 
                // Achar as céluas vizinhas à interface
                else if (val_center >= 0.0) { 
                    // Distância de busca (Buffer + Grading)
                    // --- SONDAGEM GENERALIZADA (DIMENSIONAL-AGNOSTIC) ---
                    // Itera por cada dimensão (x, y, [z])
                    for (int dim = 0; dim < DIM; dim++) {
                        // Itera por direções (-1 e +1)
                        for (int dir = -1; dir <= 1; dir += 2) {
                            
                            // Cria ponto de sonda copiando o centro manualmente
                            Point probe;
                            for(int k=0; k<DIM; k++) probe[k] = center[k];
                            
                            // Aplica o deslocamento na dimensão atual
                            //  Pega a bounding box
                            probe[dim] += dir * search_dist[pass];

                            // Verifica a sonda
                            hig_cell *c_probe = sd_get_cell_with_point(sdm, probe);
                            if (c_probe) {
                                int pid = mp_lookup(mp, hig_get_cid(c_probe));
                                if (pid >= 0) {
                                    real val_probe = dp_get_value(ns->ed.mult.dpfracvol, pid);
                                    
                                    // Detecta mudança brusca (Fluido A <-> B) ou toque na Interface
                                    if (fabs(val_probe - val_center) > 0.001 || 
                                       (val_probe > 0.001 && val_probe < 0.999)) {
                                        should_refine = 1;
                                        // Break duplo para sair dos loops de direção/dimensão
                                        goto found_refinement; 
                                    }
                                }
                            }
                        }
                    }
                    found_refinement:; // Label para sair do loop aninhado
                }
                if (should_refine) {
                    if (count >= capacity) {
                        capacity *= 2;
                        cells_to_refine = (hig_cell **)realloc(cells_to_refine, capacity * sizeof(hig_cell*));
                    }
                    cells_to_refine[count++] = c_copy;
                }
            }
            higcit_nextcell(it);
        }
        higcit_destroy(it);

        if (count == 0) {
            free(cells_to_refine);
            break; 
        }

        // Aplica o refino
        for (int i = 0; i < count; i++) {
            int numcells[DIM];
            hig_get_cells_per_dim(cells_to_refine[i], numcells);
            // Generaliza a divisão para DIM (ex: 2x2 em 2D, 2x2x2 em 3D)
            // Assume que a macro aceita args extras ou ignore se for 2D
            POINT_ASSIGN_INTS(numcells, 1, 1, 1)
            // Para ser 100% seguro em POINT_ASSIGN_INTS com DIM variável:
            for(int d=0; d<DIM; d++) numcells[d] = 2; 
            
            hig_refine_uniform(cells_to_refine[i], numcells);
        }
        free(cells_to_refine);
    } 

    // Salva VTK
    char filename[256];
    sprintf(filename, "preview_mesh_r%d_f%d.vtk", myrank, frame_id);
    FILE *fd = fopen(filename, "w");
    if (fd) {
        // Tenta usar a função disponível na sua versão
        higio_print_in_vtk2d(fd, root_copy); 
        // Nota: Se estiver em 3D, precisará trocar para higio_print_in_vtk3d(fd, root_copy) 
        // ou a função genérica hig_print_vtk(root_copy, fd) se disponível.
        fclose(fd);
    }

    hig_destroy(root_copy);
    if(myrank == 0) printf("Preview: Concluido.\n");
    // exit(0);
}

// void higflow_refine_transition_zone(higflow_solver *ns) {
//     // Obtém o domínio real da simulação
//     sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
//     mp_mapper *mp = sd_get_domain_mapper(sdm);
//
//     int myrank;
//     MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//
//     // --- PASSO 1: Coleta de Sementes (Interface) ---
//     // Necessário para calcular a distância
//     int seed_cap = 1000;
//     int seed_count = 0;
//     InterfaceSeed *seeds = (InterfaceSeed *)malloc(seed_cap * sizeof(InterfaceSeed));
//
//     higcit_celliterator *it = sd_get_domain_celliterator(sdm);
//     while (!higcit_isfinished(it)) {
//         hig_cell *c = higcit_getcell(it);
//         int clid = mp_lookup(mp, hig_get_cid(c));
//         if (clid >= 0) {
//             real val = dp_get_value(ns->ed.mult.dpfracvol, clid);
//             // Apenas mistura real (exclui paredes)
//             if (val > 0.001 && val < 0.999) {
//                 if (seed_count >= seed_cap) {
//                     seed_cap *= 2;
//                     seeds = (InterfaceSeed *)realloc(seeds, seed_cap * sizeof(InterfaceSeed));
//                 }
//                 hig_get_center(c, seeds[seed_count].center);
//                 seed_count++;
//             }
//         }
//         higcit_nextcell(it);
//     }
//     higcit_destroy(it);
//
//     if (seed_count == 0) {
//         free(seeds);
//         return; // Nada a fazer se não há interface
//     }
//
//     // --- PASSO 2: Identifica Células de Transição para Refinar ---
//     int refine_cap = 1000;
//     int refine_count = 0;
//     hig_cell **cells_to_refine = (hig_cell **)malloc(refine_cap * sizeof(hig_cell*));
//
//     it = sd_get_domain_celliterator(sdm);
//     while (!higcit_isfinished(it)) {
//         hig_cell *c = higcit_getcell(it);
//
//         // Verifica o nível atual
//         int current_lvl = get_cell_level(c);
//
//         // Só nos interessa refinar se estiver ABAIXO do nível de transição
//         if (current_lvl < TRANSITION_LEVEL) {
//             Point center;
//             hig_get_center(c, center);
//
//             // Calcula distância mínima à interface
//             real min_dist_sq = 1.0e20;
//             for (int k = 0; k < seed_count; k++) {
//                 real d2 = dist_sq(center, seeds[k].center);
//                 if (d2 < min_dist_sq) min_dist_sq = d2;
//             }
//             real dist = sqrt(min_dist_sq);
//
//             // Verifica se está na BANDA DE TRANSIÇÃO
//             if (dist > SEARCH_DIST && dist <= SEARCH_DIST_TRANSITION) {
//                 if (refine_count >= refine_cap) {
//                     refine_cap *= 2;
//                     cells_to_refine = (hig_cell **)realloc(cells_to_refine, refine_cap * sizeof(hig_cell*));
//                 }
//                 cells_to_refine[refine_count++] = c;
//             }
//         }
//         higcit_nextcell(it);
//     }
//     higcit_destroy(it);
//     free(seeds);
//
//     // --- PASSO 3: Aplica Refinamento e Atualiza Solver ---
//     if (refine_count > 0) {
//         if(myrank == 0) printf("Transition Refine: Refinando %d celulas na zona de transicao...\n", refine_count);
//
//         // 3.1 Refina
//         for(int i=0; i<refine_count; i++) {
//             int nc[DIM];
//             // Refino padrão (ex: dividir em 2x2)
//             for(int d=0; d<DIM; d++) nc[d] = 2; 
//             hig_refine_uniform(cells_to_refine[i], nc);
//         }
//     }
//
//     free(cells_to_refine);
// }
