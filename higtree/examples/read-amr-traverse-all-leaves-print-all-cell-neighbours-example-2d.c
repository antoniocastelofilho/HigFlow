/**
 * @file amr_neighbor_analyzer.c
 * @brief AMR Grid Neighbor Analysis Tool
 * 
 * This program reads an AMR (Adaptive Mesh Refinement) 2D grid file and
 * analyzes the neighbor relationships between all leaf cells in the grid
 * hierarchy. It demonstrates HIG tree traversal and neighbor detection
 * capabilities.
 */

#include<stdio.h>
#include<stdlib.h>

#include "higtree.h"
#include "higtree-io.h"

/**
 * @brief Main function for AMR neighbor analysis
 * 
 * This program:
 * 1. Reads an AMR 2D grid from file
 * 2. Iterates through all leaf cells in the grid
 * 3. For each leaf cell, finds and prints all its neighbors
 * 4. Outputs cell ID and neighbor IDs for analysis
 * 
 * @param argc Number of command line arguments
 * @param argv Array of command line arguments
 * @return int Exit status (0 on success)
 * 
 * @note Usage: ./program amr2d.file
 * @note The program expects an AMR 2D format file as input
 * 
 * @example
 * Input file: grid.amr
 * Output:
 * neighbours of cell 0: 1 2 3
 * neighbours of cell 1: 0 4 5
 * ...
 */
int main(int argc, char *argv[]) {
    /* Validate command line arguments */
    /* Example mesh mesh2d-2.amr */
    if (argc <= 1) {
        printf("usage: %s amr2d.file\n", argv[0]);
        exit(0);
    }
    
    /* Open and read AMR 2D grid file */
    FILE * fdin = fopen(argv[1], "r");
    hig_cell * root = higio_read_from_amr2d(fdin);
    fclose(fdin);

    /* Create iterator for all leaf cells in the AMR grid */
    higcit_celliterator *it;
    
    /* Iterate through all leaf cells in the grid */
    for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
        /* Get current leaf cell */
        hig_cell *c = higcit_getcell(it);
        
        /* Print cell identifier */
        printf("neighbours of cell %d:", hig_get_cid(c));
        
        /* Create iterator for neighbors of current cell */
        higcit_celliterator *itn;
        
        /* Iterate through all neighbors of current cell */
        for(itn = higcit_create_neighbours(c);
              !higcit_isfinished(itn); higcit_nextcell(itn)) {
            /* Get neighbor cell */
            hig_cell *n = higcit_getcell(itn);
            
            /* Print neighbor cell identifier */
            printf(" %d", hig_get_cid(n));
        }
        
        /* Complete the neighbor list for current cell */
        printf("\n");
        
        /* Clean up neighbor iterator */
        higcit_destroy(itn);
    }
    
    /* Clean up leaf cell iterator */
    higcit_destroy(it);

    /* Free the AMR grid hierarchy */
    hig_destroy(root);
    
    return 0;
}
