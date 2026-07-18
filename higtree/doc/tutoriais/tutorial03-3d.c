#include<stdio.h>

#include "coord.h"
#include "higtree-io.h"
#include "higtree.h"

int main (int argc, char *argv[]) {
    higtree_initialize(&argc, &argv);
    Point lo, hi;
    POINT_ASSIGN_SCALAR(lo, -1.0); //defines the point (-1,-1,-1)
    POINT_ASSIGN_SCALAR(hi, 1.0);  //defines the point (1, ,1 ,1)

    //creates the root of the cell. it's a cube centered at the origin with width=2
    hig_cell *root = hig_create_root(lo, hi); 

    //refines into 2 cells in the x direction, 3 in the y direction and 4 in the z direction
    int nc[3] = {2,3,4}; 
    hig_refine_uniform(root, nc);

    //refines even further one cell in the mesh
    Point a;
    POINT_ASSIGN_SCALAR(a, 0.5);
    hig_cell *c = hig_get_cell_with_point(root, a);
    nc[0] = 2; nc[1] = 2; nc[2] = 2;
    hig_refine_uniform(c, nc);

    //Saves the mesh created
    FILE *f = fopen("VTKS/tut03.vtk", "w");
    higio_print_in_vtk3d(f, root);

    fclose(f);
    hig_destroy(c);
}
