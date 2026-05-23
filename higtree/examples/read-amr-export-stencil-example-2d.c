/*
*Reads ONE amr file, sets Dirichlet boundary conditions and export the stencil of each cell
*file_coords_out exports a .csv with the relation between id (global), coordinate and delta.
*file_weights_out exports a .txt with a python dictionary. The keys are the coordinates of the that is not in the mesh
 and the value associated is a list with tuples containing the id and the weight of the point.
*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "utils.h"
#include "higtree.h"
#include "higtree-io.h"
#include "mapper.h"
#include "wls.h"
#include "glib.h"
#include "domain.h"
#include "solver.h"

int boundary_func = 0;
real f(Point p){
    switch (boundary_func){
    case 0:
        return 1.0;
        break;
    }
}

int rhs_func = 0;
real g(Point p){
    switch(rhs_func){
    case 0:
        return 0.0;
        break;
    }
}


void set_dirichlet_boundary_conditions(sim_domain *sd) {
	higcit_celliterator *it;
	int numbcs = sd_get_num_bcs(sd, DIRICHLET);
	for (int i = 0; i < numbcs; i++) {
		sim_boundary *bc = sd_get_bc(sd, DIRICHLET, i);
		mp_mapper *bm = sb_get_mapper(bc);
		for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
			hig_cell *bcell = higcit_getcell(it);
			Point bccenter;
			hig_get_center(bcell, bccenter);
			int bcgid = mp_lookup(bm, hig_get_cid(bcell));
			sb_set_value(bc, bcgid, f(bccenter));
		}
		higcit_destroy(it);
	}
}

int main(int argc,char *argv[]) {
    int write_vtk = 1;

    assert(DIM == 2 || DIM == 3);
	higtree_initialize(&argc, &argv);
	FILE *file_in = fopen(argv[1], "r");

	char *p = strtok(argv[1], ".");

	char file_coords_out[100]  = "/media/kainas/ubuntu/documents/data_stencils/file_coords_";
    char file_weights_out[100] = "/media/kainas/ubuntu/documents/data_stencils/file_weights_";
	char vtk_file_name[100] = "/media/kainas/ubuntu/documents/vtk_files/";
	
	strcat(file_coords_out, p);
	strcat(file_coords_out, ".csv");


	strcat(file_weights_out,p);
	strcat(file_weights_out, ".txt");
	

	strcat(vtk_file_name, p);
	strcat(vtk_file_name, ".vtk");

    FILE *vtk_file = fopen(vtk_file_name, "w");

    higio_amr_info *amr_info = higio_read_amr_info(file_in);
    hig_cell *root = higio_read_from_amr_info(amr_info);
	higio_print_in_vtk2d(vtk_file, root);

	fclose(vtk_file);
    fclose(file_in);

    higcit_celliterator *it;
    mp_mapper *m = mp_create();
    
    sim_domain *sd = sd_create(m);
    sd_add_higtree(sd, root);

	it = higcit_create_all_leaves(root);
	int size = mp_assign_from_celliterator(m, it, 0);
    higcit_destroy(it);
	m = sd_get_domain_mapper(sd);


    sd_set_interpolator_order(sd, 3);

	

    for(int dim = 0; dim < DIM; dim++) {
		for(int dir = 0; dir < 2; dir++) {
			hig_cell *bcg = higio_read_bc_from_amr(amr_info, dim, dir);
			mp_mapper *bm = mp_create();
			it = higcit_create_all_leaves(bcg);
			mp_assign_from_celliterator(bm, it, 0);
			higcit_destroy(it);
			sim_boundary *bc = sb_create(bcg, DIRICHLET, bm);
			sd_add_boundary(sd, bc);
		}
	}
    set_dirichlet_boundary_conditions(sd);

    higio_amr_info_destroy(amr_info);


    sim_stencil *stn = stn_create();

	printf("writing files...\n");

	FILE *file_coords = fopen(file_coords_out, "w");
	FILE *file_weights; 
	file_weights = fopen(file_weights_out, "w");
	fprintf(file_coords, "gid,x,y,dx,dy\n");
	
	fprintf(file_weights,"{");

	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		Point ccenter;
		Point cdelta;
		hig_get_center(c, ccenter);
		real cx = ccenter[0];
		real cy = ccenter[1];
		hig_get_delta(c, cdelta);
		real dx = cdelta[0];
		real dy = cdelta[1];
		dy *= 1.0;

		Point center = {cx     , cy     };
		Point left   = {cx - dx, cy     };
		Point right  = {cx + dx, cy     };
		Point top    = {cx,      cy + dy};
		Point bottom = {cx,      cy - dy};

		int cgid = mp_lookup(m, hig_get_cid(c));
		fprintf(file_coords, "%d,%f,%f,%f,%f\n", cgid, cx, cy, dx, dy);
		
		stn_set_rhs(stn, g(center));
		for(int i = 0; i<4; i++) {
			stn_reset(stn);

			Point current;
			switch (i)
			{
			case 0:
				POINT_ASSIGN_REALS(current, right[0], right[1]); 
				break;
			case 1:
				POINT_ASSIGN_REALS(current, left[0], left[1]); 
				break;
			case 2:
				POINT_ASSIGN_REALS(current, top[0], top[1]); 
				break;
			case 3:
				POINT_ASSIGN_REALS(current, bottom[0], bottom[1]); 
				break;
			}

			bool cache_found = false;
			_cache_wls_choice *cached = NULL;
			cache_found = ptm_lookup(sd->cwls, current, (void **)&cached);

			if(cache_found == false) {
				sd_get_stencil(sd, center, current, 1.0, stn);

				int *ids = stn_get_ids(stn); //global ids
				real *vals = stn_get_vals(stn);
				int numelems = stn_get_numelems(stn);

				if(numelems > 1){
					fprintf(file_weights, "(%.12f, %.12f): [", current[0], current[1]);
					for(int j = 0; j < numelems-1; j++) {
						fprintf(file_weights, "[%d, %.12f],", ids[j], vals[j]);
					}
					fprintf(file_weights, "[%d, %.12f]],\n\n", ids[numelems-1], vals[numelems-1]);
				}
			}
		} //for that iterates files
	} //for that iterates tree


	fseek(file_weights, -3, SEEK_CUR);
	fprintf(file_weights, "}");

	higcit_destroy(it);
	stn_destroy(stn);

	fclose(file_coords);
	fclose(file_weights);

	printf("done\n");
    return 0;
}	