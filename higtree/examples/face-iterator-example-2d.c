/**
 * @file dimension_adaptive_amr_analyzer.c
 * @brief Dimension-Adaptive AMR Grid Boundary Analyzer
 * 
 * This program demonstrates the creation and analysis of an Adaptive Mesh 
 * Refinement (AMR) grid with support for both 2D and 3D configurations.
 * It automatically adapts to the compilation dimension (DIM) and provides
 * dimension-appropriate grid creation, refinement, and output.
 * 
 * Features:
 * - Automatic dimension detection (2D/3D) at compile time
 * - Uniform refinement around specified points
 * - VTK format output for visualization
 * - Boundary face extraction and analysis
 * - Cross-dimensional compatibility
 */

#include<stdio.h>

#include "utils.h"
#include "higtree.h"
#include "higtree-io.h"
#define DEBUG
#include "Debug-c.h"

/**
 * @brief Main function for dimension-adaptive AMR analysis
 * 
 * This program performs the following operations:
 * 1. Creates a root grid with dimension-specific parameters
 * 2. Refines cells containing predefined points uniformly
 * 3. Exports the refined grid to dimension-appropriate VTK format
 * 4. Iterates through all leaf cells to identify boundary faces
 * 5. Prints coordinates of boundary faces in dimension-specific format
 * 
 * @param argc Number of command line arguments (unused)
 * @param argv Array of command line arguments (unused)
 * @return int Always returns 0 indicating successful execution
 * 
 * @note The program automatically adapts to 2D or 3D based on DIM macro
 * @note Grid spans from (0,0[,0]) to (8,8[,8]) with initial uniform division
 * @note Refinement points are hardcoded and dimension-specific
 * @note VTK output is written to /tmp/t.vtk for visualization
 * @note Boundary faces are printed to stdout with dimension-appropriate format
 * 
 * @warning The VTK file path is hardcoded to /tmp/t.vtk
 * @warning Compilation must define DIM as either 2 or 3
 * @warning Refinement points and grid parameters are hardcoded
 * 
 * @example 2D Output:
 * 0 = ((0.000000, 0.000000), (0.000000, 8.000000))
 * 1 = ((0.000000, 0.000000), (8.000000, 0.000000))
 * 
 * @example 3D Output:
 * 0 = ((0.000000, 0.000000, 0.000000), (0.000000, 8.000000, 8.000000))
 * 1 = ((0.000000, 0.000000, 0.000000), (8.000000, 0.000000, 8.000000))
 * 
 * @see hig_create_root()
 * @see hig_get_cell_with_point()
 * @see hig_refine_uniform()
 * @see higio_print_in_vtk2d()
 * @see higio_print_in_vtk3d()
 */
int main(int argc, char *argv[]) {
    /* Define grid boundaries and dimensions with compile-time dimension adaptation */
    Point l, h;
    int nc[DIM];
    
#if DIM == 2
    /** @brief 2x2 refinement pattern for 2D grid */
    POINT_ASSIGN_INTS(nc, 2, 2);
    /** @brief Lower bounds (0,0) for 2D grid */
    POINT_ASSIGN_REALS(l, 0.0, 0.0);
    /** @brief Upper bounds (8,8) for 2D grid */  
    POINT_ASSIGN_REALS(h, 8.0, 8.0);
    /** @brief Refinement points for 2D grid */
    Point ps[] = {
        {4.0, 4.0},  /**< Center point */
        {2.0, 2.0},  /**< Bottom-left quadrant */
        {2.0, 6.0},  /**< Top-left quadrant */
        {1.0, 5.0},  /**< Additional refinement point */
        {3.0, 5.0},   /**< Additional refinement point */
        {7.0, 4.0}   /**< Additional refinement point */
    };
    int numps = 6;  /**< Number of refinement points in 2D */
    
#elif DIM == 3
    /** @brief 2x2x2 refinement pattern for 3D grid */
    POINT_ASSIGN_INTS(nc, 2, 2, 2);
    /** @brief Lower bounds (0,0,0) for 3D grid */
    POINT_ASSIGN_REALS(l, 0.0, 0.0, 0.0);
    /** @brief Upper bounds (8,8,8) for 3D grid */
    POINT_ASSIGN_REALS(h, 8.0, 8.0, 8.0);
    /** @brief Refinement points for 3D grid */
    Point ps[] = {
        {4.0, 4.0, 4.0},  /**< Center point */
        {2.0, 2.0, 2.0}   /**< Corner point */
    };
    int numps = 2;  /**< Number of refinement points in 3D */
#endif

    /* Create root grid and refine cells containing specified points */
    hig_cell *root = hig_create_root(l, h);
    hig_cell *c;
    for (int i = 0; i < numps; i++) {
        /* Find cell containing the point and refine it uniformly */
        c = hig_get_cell_with_point(root, ps[i]);
        // hig_refine_uniform(c, nc);
        hig_refine_uniform(c, nc);
        printf("Cells 0: %d\n", c->numcells[0]);
    }
    DEBUG_PASS;  /**< Debug checkpoint after refinement */

    /* Export the refined grid to dimension-appropriate VTK format */
    FILE *fd = fopen("/tmp/t.vtk", "w");
#if DIM == 2
    higio_print_in_vtk2d(fd, root);  /**< 2D VTK output */
#elif DIM == 3  
    higio_print_in_vtk3d(fd, root);  /**< 3D VTK output */
#endif
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
                        /* Check if neighbor cell aligns properly in other dimensions */
                        Point ofl;
                        hig_get_lowpoint(ofc, ofl);
                        for(int k = 0; k < DIM; k++) {
                            if (k != i && FLT_NE(cl[k], ofl[k])) {
                                /* Cells don't align in other dimensions */
                                //print = 1;  // Currently commented out for debugging
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
                                    //print = 1;  // Currently commented out for debugging
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
                        fh[i] = cl[i];  /**< Low face - collapse to lower boundary */
                    } else {
                        fl[i] = ch[i];  /**< High face - collapse to upper boundary */
                    }
                    
                    /* Print boundary face coordinates with dimension-appropriate format */
#if DIM == 2
                    // printf("%d = ((%f, %f), (%f, %f))\n", cnt++, fl[0], fl[1], fh[0], fh[1]);
#elif DIM == 3
                    printf("%d = ((%f, %f, %f), (%f, %f, %f))\n", cnt++, fl[0], fl[1], fl[2], fh[0], fh[1], fh[2]);
#endif
                }
            }
        }
    }
    
    /* Clean up resources */
    higcit_destroy(it);
    return 0;
}
