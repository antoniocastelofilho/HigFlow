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
#include "nn-weights.h"


#define DEBUG
#include"Debug-c.h"


extern int init_mlp_model(const char* model_path);


extern void nn_inference(int numpts, wls_item items[], int max_n, real w[]);

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
			sb_set_value(bc, bcgid, u(bccenter[0], bccenter[1]));
		}
		higcit_destroy(it);
	}
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
	double norm_inf,norm_2;

	printf("Testing with WLS interpolator\n");
	printf("creating solver\n");
	DEBUG_DIFF_TIME_START;
    solver *slv = slv_create(SOLVER_ANY, 0, size);
    slv_set_maxnonzeros(slv, 300);


	printf("calculating stencils\n");
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
	
	/* Montagem e solução do sistema linear */
	printf("assembling...\n");
	slv_assemble(slv);


	printf("solving...\n");
	slv_solve(slv);
	DEBUG_DIFF_TIME_FINISH("Finished wls");
    /* Cálculo de erros */
	norm_inf = -1.0;
	norm_2 = 0.0;

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





	
	sd_use_cache(sd, 0);
	init_mlp_model("../../nn_models/mlp_relativo.pt");
	wls_set_type(sd->inter.wls, NEURALNETWORK);
	printf("Testing with MLP interpolator\n");
	printf("creating solver\n");
	DEBUG_DIFF_TIME_START;
    solver *slv_mlp = slv_create(SOLVER_ANY, 0, size);
    slv_set_maxnonzeros(slv_mlp, 300);

	printf("calculating stencils\n");
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

		// Define pontos para estêncil de 5 pontos
		Point center = {cx     , cy     };
		Point left   = {cx - dx, cy     };
		Point right  = {cx + dx, cy     };
		Point top    = {cx,      cy + dy};
		Point bottom = {cx,      cy - dy};

		// Monta estêncil do Laplaciano discreto
		stn_reset(stn);
		stn_set_rhs(stn, g(cx, cy));  // Lado direito
		
		// Coeficientes do esquema de diferenças finitas:
		//∇²u ≈ (u_{i-1,j} + u_{i+1,j} - 2u_{i,j})/dx² + (u_{i,j-1} + u_{i,j+1} - 2u_{i,j})/dy²
		sd_get_stencil(sd, center, center, (-2.0/(dx*dx)-2.0/(dy*dy)), stn);  // Ponto central
		sd_get_stencil(sd, center, left,    (1.0/(dx*dx)),              stn);  // Vizinho esquerdo
		sd_get_stencil(sd, center, right,   (1.0/(dx*dx)),              stn);  // Vizinho direito
		sd_get_stencil(sd, center, top,     (1.0/(dy*dy)),              stn);  // Vizinho superior
		sd_get_stencil(sd, center, bottom,  (1.0/(dy*dy)),              stn);  // Vizinho inferior

		// Adiciona linha na matriz do sistema
		int *ids = stn_get_ids(stn);
		real *vals = stn_get_vals(stn);
		int numelems = stn_get_numelems(stn);

		int cgid = mp_lookup(m, hig_get_cid(c));
		slv_set_Ai(slv_mlp, cgid, numelems, ids, vals);
		slv_set_bi(slv_mlp, cgid, stn_get_rhs(stn));
	}
    higcit_destroy(it);
	
	// Montagem e solução do sistema linear
	printf("assembling...\n");
	slv_assemble(slv_mlp);

	printf("solving...\n");
	slv_solve(slv_mlp);
	DEBUG_DIFF_TIME_FINISH("Finished mlp");
    // Cálculo de erros
	norm_inf = -1.0;
	norm_2 = 0.0;

	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		int cgid = mp_lookup(m, hig_get_cid(c));
		Point ccenter;
		hig_get_center(c, ccenter);
		real cx = ccenter[0];
		real cy = ccenter[1];

		// Erro absoluto no ponto
		real error = fabs(u(cx, cy) - slv_get_xi(slv_mlp, cgid));
		norm_2 += error*error;
		norm_inf = ((norm_inf<error) ? error : norm_inf);
	}
	norm_2 /= (real) size;
	norm_2 = sqrt(norm_2);

	// Exibe resultados
	printf("MLP: norm_inf = %f\n", norm_inf);
	printf("MLP: norm_2 = %f\n", norm_2);
	printf("{%e, %e}\n\n", norm_inf, norm_2);
	higcit_destroy(it);

    return 0;
}	