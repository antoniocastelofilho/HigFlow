/*
 * ===========================================================================
 * LÓGICA DE PRÉ-VISUALIZAÇÃO DE REFINAMENTO DE MALHA (PREVIEW MESH ADAPT)
 * ===========================================================================
 *
 * Esta função gera uma cópia virtual da malha e simula o refinamento
 * adaptativo baseado na distância da interface, sem alterar a simulação real.
 *
 * 1. PREPARAÇÃO (SETUP)
 * - Clona a estrutura da malha original (HigTree) para 'root_copy'.
 * - Define os raios de influência física para cada nível de refinamento:
 * * Nível 3 (Fino):   Aplicado se dist. da interface <= 0.03
 * * Nível 2 (Médio):  Aplicado se dist. da interface <= 0.05
 * * Nível 1 (Grosso): Aplicado se dist. da interface <= 0.08
 *
 * 2. IDENTIFICAÇÃO DE SEMENTES (INTERFACE SEEDS)
 * - Varre o domínio original buscando células onde a física da interface
 * acontece. Critério: Fração de Volume entre 0.001 e 0.999.
 * - Armazena as coordenadas (centros) dessas células numa lista.
 *
 * 3. PROPAGAÇÃO DE REFINAMENTO (LOOP DE PASSADAS)
 * - Executa (MAX_LEVEL + 1) iterações para garantir que o refinamento
 * cresça gradualmente (ex: uma célula Nível 0 vira 1, depois 2...).
 *
 * PARA CADA PASSADA (PASS):
 * a. PARA CADA SEMENTE (Ponto de Interface):
 * - Cria uma Bounding Box (Caixa) ao redor da semente com o
 * tamanho do maior raio de influência (0.08).
 * - Itera apenas sobre as células vizinhas DENTRO desta caixa.
 *
 * PARA CADA VIZINHO (Cell Neighbor):
 * - Calcula a Distância Euclidiana (d) até a Semente.
 * - Determina o Nível Alvo (Target) hierarquicamente:
 * Se (d <= 0.03) ENTÃO Alvo = 3 (Zona Crítica)
 * Senão Se (d <= 0.05) ENTÃO Alvo = 2 (Zona Média)
 * Senão Se (d <= 0.08) ENTÃO Alvo = 1 (Zona Transição)
 *
 * - Se (Nível Atual < Nível Alvo):
 * Marca o vizinho para ser refinado.
 *
 * b. APLICAÇÃO:
 * - Refina todas as células marcadas (divisão uniforme).
 * - Se nenhuma célula precisar de refino, encerra o processo.
 *
 * 4. EXPORTAÇÃO
 * - Salva a malha clonada em VTK para visualização (Paraview).
 * - Destrói a cópia e libera memória.
 * ===========================================================================
 */#include "higtree.h"
#include "higtree-io.h"
#include "higtree-iterator.h"
#include "pdomain.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Definição do nível máximo de refinamento na interface
#define MAX_REF_LEVEL 3

// Estrutura para armazenar o centro das células que contêm a interface
typedef struct {
    Point center;
} InterfaceSeed;

// Função auxiliar para obter o nível atual de uma célula na árvore
// Conta o número de pais que uma célula tem
int get_cell_level(hig_cell *c) {
    int level = 0;
    hig_cell *p = hig_get_parent(c);
    while (p != NULL) {
        level++;
        p = hig_get_parent(p);
    }
    return level;
}

void higflow_save_refined_mesh_preview(higflow_solver *ns, int frame_id) {
    // Acessa o domínio da simulação multifásica para ler a propriedade real
    sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
    mp_mapper *mp = sd_get_domain_mapper(sdm);

    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // 1. CLONAGEM DA MALHA (PREVIEW)
    // Trabalhamos numa cópia para não alterar a malha da simulação
    hig_cell *root_original = sd_get_higtree(sdm, 0);
    if (root_original == NULL) return;

    hig_cell *root_copy = hig_clone(root_original);
    if (root_copy == NULL) return;

    if (myrank == 0) {
        printf("Preview: Refinando malha (BBox) - Frame %d...\n", frame_id);
    }

    // 2. DEFINIÇÃO DE DISTÂNCIAS (CRITÉRIO FÍSICO)
    // Define o raio de influência de cada nível a partir da interface.
    float search_dist[MAX_REF_LEVEL];

    // Zona Grossa: Células a até 0.08 unidades devem ser Nível 1
    search_dist[0] = 0.08;
    // Zona Média: Células a até 0.05 unidades devem ser Nível 2
    search_dist[1] = 0.005;
    // Zona Fina: Células a até 0.03 unidades devem ser Nível 3
    search_dist[2] = 0.003;

    // A maior distância define o tamanho da Bounding Box de busca
    real max_search_dist = search_dist[0];

    // 3. IDENTIFICAÇÃO DAS CÉLULAS DE INTERFACE ("SEMENTES")
    // Varremos toda a malha original para encontrar a física da interface.
    int seed_cap = 1000;
    int seed_count = 0;
    InterfaceSeed *seeds = (InterfaceSeed *)malloc(
        seed_cap * sizeof(InterfaceSeed));

    higcit_celliterator *it = sd_get_domain_celliterator(sdm);
    while (!higcit_isfinished(it)) {
        hig_cell *c = higcit_getcell(it);
        int clid = mp_lookup(mp, hig_get_cid(c));

        if (clid >= 0) {
            // dpfracvol (Fração de Volume):
            // 0.0 -> Apenas Fluido A; 1.0 -> Apenas Fluido B
            // 0.0 < val < 1.0 -> Célula contém a INTERFACE
            real val = dp_get_value(ns->ed.mult.dpfracvol, clid);

            // Filtro numérico: considera interface se estritamente entre 0 e 1
            if (val > 0.001 && val < 0.999) {
                if (seed_count >= seed_cap) {
                    seed_cap *= 2;
                    seeds = (InterfaceSeed *)realloc(seeds,
                        seed_cap * sizeof(InterfaceSeed));
                }
                hig_get_center(c, seeds[seed_count].center);
                seed_count++;
            }
        }
        higcit_nextcell(it);
    }
    higcit_destroy(it);

    if (seed_count == 0) {
        if (myrank == 0) printf("Preview: Nenhuma interface encontrada.\n");
        free(seeds);
        hig_destroy(root_copy);
        return;
    }

    // 4. LOOP DE REFINAMENTO (MULTIBLOCO / PASSADAS)
    // O refinamento em árvore precisa ser gradual (0 -> 1 -> 2 -> 3).
    // Executamos várias passadas para propagar corretamente.
    int max_passes = MAX_REF_LEVEL + 1;

    for (int pass = 0; pass < max_passes; pass++) {
        int refine_cap = 1000;
        int refine_count = 0;
        hig_cell **to_refine = (hig_cell **)malloc(
            refine_cap * sizeof(hig_cell*));

        // --- Iteração sobre as Sementes (Interfaces) ---
        for (int i = 0; i < seed_count; i++) {
            Point center_seed;
            POINT_ASSIGN(center_seed, seeds[i].center);

            // Criação da Bounding Box (BBox) em torno da semente
            // Otimização: Verificamos apenas células no raio máximo.
            Point lo, hi;
            for (int d = 0; d < DIM; d++) {
                lo[d] = center_seed[d] - max_search_dist;
                hi[d] = center_seed[d] + max_search_dist;
            }

            // --- Iterador Geométrico (BBox) ---
            // Itera apenas sobre células da CÓPIA dentro da caixa
            higcit_celliterator *it_bb = higcit_create_bounding_box(
                root_copy, lo, hi);

            while (!higcit_isfinished(it_bb)) {
                hig_cell *c_neighbor = higcit_getcell(it_bb);

                // Só refinamos folhas (células sem filhos)
                if (hig_get_number_of_children(c_neighbor) == 0) {
                    Point center_neigh;
                    hig_get_center(c_neighbor, center_neigh);

                    // Distância euclidiana entre interface e vizinho
                    real dist = co_distance(center_seed, center_neigh);

                    // Lógica de Gradiente (Grading)
                    // Define o nível ideal baseado na proximidade.
                    int target_level = 0;

                    // Verifica do mais restritivo (fino) ao permissivo
                    if (dist <= search_dist[2]) {
                        target_level = 3; // Muito perto -> Fino
                    }
                    else if (dist <= search_dist[1]) {
                        target_level = 2; // Perto -> Médio
                    }
                    else if (dist <= search_dist[0]) {
                        target_level = 1; // Longe -> Grosso
                    }

                    // Verifica se célula está aquém do nível desejado
                    int current_lvl = get_cell_level(c_neighbor);
                    if (current_lvl < target_level) {
                        if (refine_count >= refine_cap) {
                            refine_cap *= 2;
                            to_refine = (hig_cell **)realloc(to_refine,
                                refine_cap * sizeof(hig_cell*));
                        }
                        // --- Aplicação do Refinamento ---
                        int nc[DIM];
                        // Divide em 2 em cada direção (sobe +1 nível)
                        for (int d = 0; d < DIM; d++) nc[d] = 2;
                        hig_refine_uniform(c_neighbor, nc);
                    }
                }
                higcit_nextcell(it_bb);
            }
            higcit_destroy(it_bb);
        }
    }

    free(seeds);

    // 5. EXPORTAÇÃO (VTK)
    char filename[256];
    sprintf(filename, "preview_mesh_r%d_f%d.vtk", myrank, frame_id);
    FILE *fd = fopen(filename, "w");
    if (fd) {
        higio_print_in_vtk2d(fd, root_copy);
        fclose(fd);
    }

    hig_destroy(root_copy);
    if (myrank == 0) printf("Preview: Concluido.\n");
    // exit(0);
}
