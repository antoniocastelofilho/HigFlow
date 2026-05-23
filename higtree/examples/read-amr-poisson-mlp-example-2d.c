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
#include "mlp-test/mlp.h"


#define DEBUG
#include"Debug-c.h"

/** 
 * @brief Must executed in the beginning ofthe code.
 * Loads the mlp model to do the inference of the interpolation weights
 * @param model_path path to the .pt containing the model
*/
extern int init_mlp_model(const char* model_path);

/**
 * @brief Makes the inference
 * @param max_n Max size of an interpolation scheme
 * @param features_por_no number of feature per point (using 3 currently)
 * @param input_features a vector containing
 * @param output_weights vector containing the returned weights 
 * @return 0 if there's no error, -1 if there's error
 */
extern int run_mlp_inference(const float* input_features, int max_n, int features_por_no, float* output_weights);

/**
 * @brief Function that is applied to all the points in the boundary
 *@param p Point that the function is applied 
*/
int boundary_func = 0;
real f(Point p){
    switch (boundary_func){
    case 0:
        return 0.0;
        break;
    }
}

int func = 0;
real u(real x, real y) {
    switch(func) {
		case 0:	return sin(x)*cos(y);                    // Solução suave
		case 1: return 4*x*x*y+2*x*x-4*y*x-5*x;        // Polinômio
		case 2: return cos(x*y)+sin(x+y);              // Combinação trigonométrica
		case 3: return cos(x-y)*sin(x+y);              // Produto trigonométrico
		case 4: return cos(x+sin(y));                  // Composição trigonométrica
		case 5: return x*x+y*y-x*y-1.1234;            // Forma quadrática
	}
}

real g(real x, real y){
    switch(func) {
		case 0: return -2.0*sin(x)*cos(y);            // Laplaciano de sin(x)cos(y)
		case 1: return 4 + 8*y;                       // Laplaciano do polinômio
		case 2: return (-x*x-y*y)*cos(x*y)-2*sin(x+y); // Laplaciano complexo
		case 3: return -4*cos(x-y)*sin(x+y);          // Laplaciano do produto
		case 4: return -(1+cos(y)*cos(y))*cos(x+sin(y))+sin(y)*(sin(x+sin(y))); // Laplaciano composto
		case 5: return 4.0;                           // Laplaciano constante
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

void mlp_calculate_stencil(sim_domain *sd, ) {
	static _wls_item_list items;
	bool use_neumann = true;
	bool use_dirichlet = true;

	// TODO: Make cache work at this level.

	// Find if point is inside the domain and if it is at a cell
	// center:
	// if it is at a cell center, already add the element to the stencil.
	point_location in_domain;
	if(funcs->find_in_center(specific, x, alpha, &in_domain, stn)) {
		return;
	}

	if(in_domain == ON_BOUNDARY) {

		// If on domain boundary, there is a chance it is at the center of a
		// Dirichlet boundary condition.
		if(dirichlet_find_in_center(d, x, alpha, stn)) {
			return;
		}

		if(use_neumann) {
			// If on a Neumann boundary, the point value will be determined by a
			// Neumann boundary condition.
			if(get_stencil_neumann_boundary_any_order(d, x, alpha, stn, funcs, &items, specific, from_sfd)) {
				return;
			}
		}

		if(use_dirichlet) {
			// If on a Dirichlet boundary, the point value will be determined by a
			// Dirichlet boundary condition.
			if(get_stencil_dirichlet_boundary(d, x, alpha, stn, funcs, &items, specific, from_sfd)) {
				return;
			}
		}
	
	} else if(in_domain == OUTSIDE_DOMAIN) {

		if(use_neumann) {
			// If outside the domain, the point value may be determined by a
			// Neumann boundary condition.
			if(get_stencil_neumann_any_order(d, x, alpha, stn, funcs, &items, specific, from_sfd)) {
				return;
			}
		}

		if(use_dirichlet) { // modificação daniel
			// If outside the domain, the point value may be determined by a
			// Dirichlet boundary condition.
			if(get_stencil_dirichlet_any_order(d, x, alpha, stn, funcs, &items, specific, from_sfd)) {
				return;
			}
		}
	}

	// If not at center, nor given by a Neumann BC, use the general case:
}




int main(int argc,char *argv[]) {
    assert(DIM == 2 || DIM == 3);
	higtree_initialize(&argc, &argv);
    //First argument is the amr file
	FILE *file_in = fopen(argv[1], "r");


    higio_amr_info *amr_info = higio_read_amr_info(file_in);
    hig_cell *root = higio_read_from_amr_info(amr_info);

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
	printf("Testing with WLS interpolator")
	printf("creating solver")
	DEBUG_DIFF_TIME_START
    solver *slv = slv_create(SOLVER_ANY, 0, size);
    slv_set_maxnonzeros(slv, 300);

	printf("calculating stencils")
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

		/* Define pontos para estêncil de 5 pontos */
		Point center = {cx     , cy     };
		Point left   = {cx - dx, cy     };
		Point right  = {cx + dx, cy     };
		Point top    = {cx,      cy + dy};
		Point bottom = {cx,      cy - dy};

		/* Monta estêncil do Laplaciano discreto */
		stn_reset(stn);
		stn_set_rhs(stn, g(cx, cy));  // Lado direito
		
		/* Coeficientes do esquema de diferenças finitas:
		   ∇²u ≈ (u_{i-1,j} + u_{i+1,j} - 2u_{i,j})/dx² + (u_{i,j-1} + u_{i,j+1} - 2u_{i,j})/dy² */
		sd_get_stencil(sd, center, center, (-2.0/(dx*dx)-2.0/(dy*dy)), stn);  // Ponto central
		sd_get_stencil(sd, center, left,    (1.0/(dx*dx)),              stn);  // Vizinho esquerdo
		sd_get_stencil(sd, center, right,   (1.0/(dx*dx)),              stn);  // Vizinho direito
		sd_get_stencil(sd, center, top,     (1.0/(dy*dy)),              stn);  // Vizinho superior
		sd_get_stencil(sd, center, bottom,  (1.0/(dy*dy)),              stn);  // Vizinho inferior

		/* Adiciona linha na matriz do sistema */
		int *ids = stn_get_ids(stn);
		real *vals = stn_get_vals(stn);
		int numelems = stn_get_numelems(stn);

		int cgid = mp_lookup(m, hig_get_cid(c));
		slv_set_Ai(slv, cgid, numelems, ids, vals);
		slv_set_bi(slv, cgid, stn_get_rhs(stn));
	}
    higcit_destroy(it);
	DEBUG_DIFF_TIME;
	
	/* Montagem e solução do sistema linear */
	printf("assembling...\n");
	slv_assemble(slv);
	DEBUG_DIFF_TIME;

	printf("solving...\n");
	slv_solve(slv);
	DEBUG_DIFF_TIME;
    /* Cálculo de erros */
	real norm_inf = -1.0;
	real norm_2 = 0.0;

	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		int cgid = mp_lookup(m, hig_get_cid(c));
		Point ccenter;
		hig_get_center(c, ccenter);
		real cx = ccenter[0];
		real cy = ccenter[1];

		/* Erro absoluto no ponto */
		real error = fabs(u(cx, cy) - slv_get_xi(slv, cgid));
		norm_2 += error*error;
		norm_inf = ((norm_inf<error) ? error : norm_inf);
	}
	norm_2 /= (real) size;
	norm_2 = sqrt(norm_2);

	/* Exibe resultados */
	printf("WLS: norm_inf = %f\n", norm_inf);
	printf("WLS: norm_2 = %f\n", norm_2);
	printf("{%e, %e}\n\n", norm_inf, norm_2);
	higcit_destroy(it);





	printf("Testing with MLP interpolator")
	printf("creating solver")
	DEBUG_DIFF_TIME_START
    solver *slv_mlp = slv_create(SOLVER_ANY, 0, size);
    slv_set_maxnonzeros(slv_mlp, 300);

	printf("calculating stencils")
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

		/* Define pontos para estêncil de 5 pontos */
		Point center = {cx     , cy     };
		Point left   = {cx - dx, cy     };
		Point right  = {cx + dx, cy     };
		Point top    = {cx,      cy + dy};
		Point bottom = {cx,      cy - dy};

		/* Monta estêncil do Laplaciano discreto */
		stn_reset(stn);
		stn_set_rhs(stn, g(cx, cy));  // Lado direito
		
		/* Coeficientes do esquema de diferenças finitas:
		   ∇²u ≈ (u_{i-1,j} + u_{i+1,j} - 2u_{i,j})/dx² + (u_{i,j-1} + u_{i,j+1} - 2u_{i,j})/dy² */
		sd_get_stencil(sd, center, center, (-2.0/(dx*dx)-2.0/(dy*dy)), stn);  // Ponto central
		sd_get_stencil(sd, center, left,    (1.0/(dx*dx)),              stn);  // Vizinho esquerdo
		sd_get_stencil(sd, center, right,   (1.0/(dx*dx)),              stn);  // Vizinho direito
		sd_get_stencil(sd, center, top,     (1.0/(dy*dy)),              stn);  // Vizinho superior
		sd_get_stencil(sd, center, bottom,  (1.0/(dy*dy)),              stn);  // Vizinho inferior

		/* Adiciona linha na matriz do sistema */
		int *ids = stn_get_ids(stn);
		real *vals = stn_get_vals(stn);
		int numelems = stn_get_numelems(stn);

		int cgid = mp_lookup(m, hig_get_cid(c));
		slv_set_Ai(slv_mlp, cgid, numelems, ids, vals);
		slv_set_bi(slv_mlp, cgid, stn_get_rhs(stn));
	}
    higcit_destroy(it);
	DEBUG_DIFF_TIME;
	
	/* Montagem e solução do sistema linear */
	printf("assembling...\n");
	slv_assemble(slv_mlp);
	DEBUG_DIFF_TIME;

	printf("solving...\n");
	slv_solve(slv_mlp);
	DEBUG_DIFF_TIME;
    /* Cálculo de erros */
	real norm_inf = -1.0;
	real norm_2 = 0.0;

	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		int cgid = mp_lookup(m, hig_get_cid(c));
		Point ccenter;
		hig_get_center(c, ccenter);
		real cx = ccenter[0];
		real cy = ccenter[1];

		/* Erro absoluto no ponto */
		real error = fabs(u(cx, cy) - slv_get_xi(slv_mlp, cgid));
		norm_2 += error*error;
		norm_inf = ((norm_inf<error) ? error : norm_inf);
	}
	norm_2 /= (real) size;
	norm_2 = sqrt(norm_2);

	/* Exibe resultados */
	printf("MLP: norm_inf = %f\n", norm_inf);
	printf("MLP: norm_2 = %f\n", norm_2);
	printf("{%e, %e}\n\n", norm_inf, norm_2);
	higcit_destroy(it);




    return 0;
}	