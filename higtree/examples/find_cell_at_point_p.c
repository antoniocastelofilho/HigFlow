/**
 * @file amr_boundary_extractor.c
 * @brief AMR Grid Boundary Face Extractor and Analyzer
 * 
 * This program demonstrates the creation and analysis of an Adaptive Mesh 
 * Refinement (AMR) grid with specific refinement patterns. It creates a root
 * grid, refines cells containing specified points, and then extracts and
 * analyzes boundary faces between grid cells.
 * 
 * The program:
 * 1. Creates a base grid and refines cells around specific points
 * 2. Outputs the grid to VTK format for visualization
 * 3. Analyzes and extracts boundary faces between grid cells
 * 4. Prints boundary face coordinates for further analysis
 */

#include<stdio.h>

#include "utils.h"
#include "higtree.h"
#include "higtree-io.h"

/**
 * @brief Main function for AMR boundary analysis
 * 
 * This program performs the following operations:
 * 1. Creates a root grid with specified dimensions
 * 2. Refines cells containing predefined points uniformly
 * 3. Exports the refined grid to VTK format
 * 4. Iterates through all leaf cells to identify boundary faces
 * 5. Prints coordinates of boundary faces between cells
 * 
 * @param argc Number of command line arguments (unused)
 * @param argv Array of command line arguments (unused)
 * @return int Always returns 0 indicating successful execution
 * 
 * @note The grid spans from (0,0) to (8,8) with initial 2x2 configuration
 * @note Refinement points are hardcoded within the program
 * @note VTK output is written to /tmp/t.vtk for visualization
 * @note Boundary faces are printed to stdout in format: id = ((x1,y1), (x2,y2))
 * 
 * @warning The VTK file path is hardcoded to /tmp/t.vtk
 * @warning Refinement points and grid parameters are hardcoded
 * 
 * @example
 * Output example:
 * 0 = ((0.000000, 0.000000), (0.000000, 8.000000))
 * 1 = ((0.000000, 0.000000), (8.000000, 0.000000))
 * 2 = ((4.000000, 0.000000), (4.000000, 8.000000))
 * 
 * @see hig_create_root()
 * @see hig_get_cell_with_point()
 * @see hig_refine_uniform()
 * @see higio_print_in_vtk2d()
 * @see higcit_create_all_leaves()
 */
int main(int argc, char *argv[]) {
    /* Define grid boundaries and dimensions */
    Point l, h;
    int nc[DIM];
    POINT_ASSIGN_INTS(nc, 2, 2);           /**< 2x2 refinement pattern */
    POINT_ASSIGN_REALS(l, 0.0, 0.0);       /**< Lower bounds of grid */
    POINT_ASSIGN_REALS(h, 8.0, 8.0);       /**< Upper bounds of grid */
    
    /* Predefined points around which cells will be refined */
    Point ps[] = {
        {4.0, 4.0},  /**< Center point */
        {2.0, 2.0},  /**< Bottom-left quadrant */
        {2.0, 6.0},  /**< Top-left quadrant */  
        {1.0, 5.0},  /**< Additional refinement point */
        {3.0, 5.0}   /**< Additional refinement point */
    };
    int numps = 5;  /**< Number of refinement points */

    /* Create root grid and refine cells containing specified points */
    hig_cell *root = hig_create_root(l, h);
    hig_cell *c;
    for (int i = 0; i < numps; i++) {
        /* Find cell containing the point and refine it uniformly */
        c = hig_get_cell_with_point(root, ps[i]);
        hig_refine_uniform(c, nc);
    }

    /* Export the refined grid to VTK format for visualization */
    FILE *fd = fopen("/tmp/t.vtk", "w");
    higio_print_in_vtk2d(fd, root);
    fclose(fd);

    /* Analyze boundary faces between grid cells */
    higcit_celliterator *it;
    int cnt = 0;  /**< Counter for boundary faces */
    
    /* Iterate through all leaf cells in the refined grid */
    for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        
        /* Get cell boundaries and center */
        Point cl, ch, ccenter;
        hig_get_lowpoint(c, cl);
        hig_get_highpoint(c, ch);
        hig_get_center(c, ccenter);

        /* Check each dimension for boundary faces */
        for(int i = 0; i < DIM; i++) {
            Point fcenter;
            Point ofcenter;
            POINT_ASSIGN(fcenter, ccenter);
            POINT_ASSIGN(ofcenter, ccenter);
            
            /* Check both low and high faces in current dimension */
            for(int j = 0; j < 2; j++) {
                /* Set face center and offset for neighbor detection */
                fcenter[i] = ((j==0)?cl[i]:ch[i]);
                ofcenter[i] = ((j==0)?cl[i]-EPSDELTA:ch[i]+EPSDELTA);
                
                int print = 0;  /**< Flag indicating if this is a boundary face */
                
                if (j == 1) {
                    /* Always consider high faces as potential boundaries */
                    print = 1;
                } else {
                    /* Check if low face is a boundary by looking for neighbors */
                    hig_cell *ofc = hig_get_cell_with_point(root, ofcenter);
                    if (ofc == NULL) {
                        /* No cell found at offset - this is a boundary */
                        print = 1;
                    } else {
                        /* Check if neighbor cell aligns properly */
                        Point ofl;
                        hig_get_lowpoint(ofc, ofl);
                        for(int k = 0; k < DIM; k++) {
                            if (k != i && FLT_NE(cl[k], ofl[k])) {
                                /* Cells don't align in other dimensions */
                                //print = 1;  // Currently commented out
                                break;
                            }
                        }
                        if (!print) {
                            /* Additional check with high points */
                            Point ofh;
                            hig_get_highpoint(ofc, ofh);
                            for(int k = 0; k < DIM; k++) {
                                if (k != i && FLT_NE(ch[k], ofh[k])) {
                                    /* Cells don't align in other dimensions */
                                    //print = 1;  // Currently commented out
                                    break;
                                }
                            }
                        }
                    }
                }
                
                /* If this is a boundary face, print its coordinates */
                if (print) {
                    Point fl, fh;
                    POINT_ASSIGN(fl, cl);
                    POINT_ASSIGN(fh, ch);
                    
                    /* Set face boundaries based on direction */
                    if (j == 0) {
                        fh[i] = cl[i];  /**< Low face */
                    } else {
                        fl[i] = ch[i];  /**< High face */
                    }
                    
                    /* Print boundary face coordinates */
                    printf("%d = ((%f, %f), (%f, %f))\n", cnt++, fl[0], fl[1], fh[0], fh[1]);
                }
            }
        }
    }
    
    /* Clean up resources */
    higcit_destroy(it);
    return 0;
}
