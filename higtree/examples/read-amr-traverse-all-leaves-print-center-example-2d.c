/**
 * @file read-amr-traverse-all-leaves-print-center-example-2d.c
 * @brief AMR Grid Cell Center Coordinate Extractor
 * 
 * This program reads an AMR (Adaptive Mesh Refinement) 2D grid file and
 * extracts the center coordinates of all leaf cells in the grid hierarchy.
 * It demonstrates HIG tree traversal and coordinate extraction capabilities.
 */

#include<stdio.h>
#include<stdlib.h>

#include "higtree.h"
#include "higtree-io.h"

/**
 * @brief Main function for AMR cell center analysis
 * 
 * This program performs the following operations:
 * 1. Reads an AMR 2D grid from a specified file
 * 2. Iterates through all leaf cells in the grid hierarchy
 * 3. Calculates and prints the center coordinates of each leaf cell
 * 4. Outputs cell ID with corresponding (x,y) center coordinates
 * 
 * @param argc Number of command line arguments
 * @param argv Array of command line argument strings
 *        - argv[0]: Program name
 *        - argv[1]: Path to AMR 2D input file
 * 
 * @return int Program exit status:
 *         - 0: Success
 *         - -1: Error (invalid arguments or file access)
 * 
 * @note The program requires exactly one command line argument (input file path)
 * @note Input file must be in AMR 2D format readable by higio_read_from_amr2d()
 * @note Output format: "center of cell <id>: (<x>, <y>)"
 * 
 * @warning Program exits with error code -1 if no input file is provided
 * @warning Ensure the input file exists and is in correct AMR 2D format
 * 
 * @example
 * Command line: ./read-amr-traverse-all-leaves-print-center-example-2d grid.amr
 * Output:
 * center of cell 0: (0.25, 0.25)
 * center of cell 1: (0.75, 0.25)
 * center of cell 2: (0.25, 0.75)
 * center of cell 3: (0.75, 0.75)
 * 
 * @see higio_read_from_amr2d()
 * @see higcit_create_all_leaves()
 * @see hig_get_center()
 */
int main(int argc, char *argv[]) {
    /* Validate command line arguments */
    if (argc <= 1) {
        printf("usage: %s amr2d.file\n", argv[0]);
        exit(-1);
    }
    
    /* Open input file for reading */
    FILE * fdin = fopen(argv[1], "r");
    
    /* Read AMR 2D grid structure from file */
    hig_cell * root = higio_read_from_amr2d(fdin);
    fclose(fdin);

    /* Create iterator for traversing all leaf cells in the AMR grid */
    higcit_celliterator *it;
    
    /* Iterate through all leaf cells in the grid hierarchy */
    for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
        /* Get current leaf cell from iterator */
        hig_cell *c = higcit_getcell(it);
        
        /* Calculate center coordinates of the cell */
        Point center;
        hig_get_center(c, center);
        
        /* Print cell ID and center coordinates */
        printf("center of cell %d: (%lg, %lg)\n", hig_get_cid(c), center[0], center[1]);
    }
    
    /* Clean up cell iterator to free resources */
    higcit_destroy(it);

    /* Destroy the AMR grid hierarchy and free allocated memory */
    hig_destroy(root);
    
    return 0;
}
