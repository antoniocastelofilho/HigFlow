#include "hig-flow-vof-finite-difference-normal-curvature.h"
#include "hig-flow-vof-hf-padrao.h"
#include "hig-flow-vof-9-cells.h"
#include "hig-flow-vof-adap-hf.h"


real direction_hf_horizontal(sim_domain *sdp, higflow_solver *ns, Point center, Point delta, int *aux_dir){
	real frac, fracvol;
	Point p;
	real hf_horizontal=0.0;
	
	p[0] = center[0]+delta[0];
	//p[1] = center[1];
	real hf_hor_dir = 0.0;
	int i = 0;
	do {
		p[1] = center[1] + (i*delta[1]-delta[1]);
		real frac = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
		hf_hor_dir += frac;
		i++;
	} while (i < 3);
	p[0] = center[0]-delta[0];
	real hf_hor_esq = 0.0;
	i = 0;
	do {
		p[1] = center[1] + (i*delta[1]-delta[1]);
		frac = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
		hf_hor_esq += frac;
		i++;
	} while (i < 3);
	
	hf_horizontal = fabs(hf_hor_dir - hf_hor_esq);

	// aproveitando esse bloco de código para definir um valor auxiliar que indicarah o sentido do vetor normal
	p[1] = center[1];
	p[0] = center[0]+delta[0];
	real frac_right = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
	p[0] = center[0]-delta[0];
	real frac_left = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
	if (frac_left == 1.0 && frac_right == 0.0) {
		*aux_dir = 0; // normal aponta para esquerda
	} else {
		*aux_dir = 1; // normal aponta para direita
	}
	
	return hf_horizontal;
}


real direction_hf_vertical(sim_domain *sdp, higflow_solver *ns, Point center, Point delta, int *aux_dir){
	real frac, fracvol;
	Point p;
	real hf_vertical=0.0;
	
	p[1] = center[1]+delta[1];
	//p[1] = center[1];
	real hf_ver_sup = 0.0;
	int i = 0;
	do {
		p[0] = center[0] + (i*delta[0]-delta[0]);
		real frac = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
		hf_ver_sup += frac;
		i++;
	} while (i < 3);
	p[1] = center[1]-delta[1];
	real hf_ver_inf = 0.0;
	i = 0;
	do {
		p[0] = center[0] + (i*delta[0]-delta[0]);
		frac = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
		hf_ver_inf += frac;
		i++;
	} while (i < 3);
	
	hf_vertical = fabs(hf_ver_sup - hf_ver_inf);
	
	// aproveitando esse bloco de código para definir um valor auxiliar que indicarah o sentido do vetor normal
	p[0] = center[0];
	p[1] = center[1]+delta[1];
	real frac_top = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
	p[1] = center[1]-delta[1];
	real frac_bottom = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
	if (frac_bottom == 1.0 && frac_top == 0.0) {
		*aux_dir = 2; // normal aponta para baixo
	} else {
		*aux_dir = 3; // normal aponta para cima
	}
	
	return hf_vertical;
}



void vertical_collumn_padrao(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *vertical, int *aux){
	int status;
	real fracvol, fracvol_aux;
	Point pp, ppt, ppb;
	*aux=0;
	// Vertical Up
	status = get_frac_vol(sdp, ns, 1, center, p, delta, &fracvol);
	if (status != 1) {
		if (status == -1) {
			printf("Up - Vertical cells with different sizes \n");
			return;
		} else {
			//printf("Up - Vertical cell out of domain \n");
			fracvol = 0.0;
		}
	}
	pp[0] = p[0];
	*vertical = fracvol;
	int i = 0;
	do {
		i++;
		pp[1] = center[1] + i * delta[1];
		status = get_frac_vol(sdp, ns, 1, center, pp, delta, &fracvol);
		if (status == 1) {
			*vertical += fracvol;
		} else if (status == 0) {
			//printf("Up - Vertical cell out of domain \n");
			fracvol = 0.0;
		} else {
			printf("Up - Vertical cells with different sizes \n");
			return;
		}
	} while (i <3 && status == 1);
	fracvol_aux = fracvol;
	ppt[0]=pp[0];ppt[1]=pp[1]+0.5*delta[1]; 
	real fracvol_top = fracvol;             
	// Vertical Down
	i = 0;
	do {
		i++;
		pp[1] = p[1] - i * delta[1];
		status = get_frac_vol(sdp, ns, 1, center, pp, delta, &fracvol);
		if (status == 1) {
			*vertical += fracvol;
		} else if (status == 0) {
			//printf("Down - Vertical cell out of domain \n");
			fracvol = 0.0;return;
		} else {
			printf("Down - Vertical cells with different sizes \n");
			return;
		}
	} while (i <3 && status == 1);

	if (fracvol == fracvol_aux) {
		//printf("The phases are not different \n");
		return;
	}

	if (fracvol == 1.0 && fracvol_top == 0.0) {
		*aux = -1;
	} else {
		*aux = 1;
	}
}

void horizontal_row_padrao(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *horizontal, int *aux){
	int status;
	real fracvol, fracvol_aux;
	Point pp, ppr, ppl;
	*aux=0;
	// Horizontal Right
	status = get_frac_vol(sdp, ns, 0, center, p, delta, &fracvol);
	if (status != 1) {
		if (status == -1) {
			printf("Right - Horizontal cells with different sizes \n");
			return;
		} else {
			//printf("Right - Horizontal cell out of domain \n");
			fracvol = 0.0;
		}
	}
	pp[1] = p[1];
	*horizontal=fracvol;
	int i = 0;
	do {
		i++;
		pp[0] = center[0] + i*delta[0];
		status=get_frac_vol(sdp, ns, 0, center, pp, delta, &fracvol);
		if (status == 1) {
			*horizontal += fracvol;
		} else if (status == 0) {
			//printf("Right - Horizontal cell out of domain \n");
			fracvol = 0.0;
		} else {
			printf("Right - Horizontal cells with different sizes \n");
			return;
		}
	} while (i <3 && status == 1);
	fracvol_aux = fracvol;
	ppr[0]=pp[0]+0.5*delta[0];
	ppr[1]=pp[1];
	real fracvol_right = fracvol;
	// Horizontal Left
	i = 0;
	//fracvol_aux=0;
	do {
		i++;
		pp[0] = p[0] - i*delta[0];
		status = get_frac_vol(sdp, ns, 0, center, pp, delta, &fracvol);
		if (status == 1) {
			*horizontal += fracvol;
		} else if (status == 0) {
			//printf("Left - Horizontal cell out of domain \n");
			fracvol = 0.0;
		} else {
			printf("Left - Horizontal cells with different sizes \n");
			return;
		}
	}while (i <3 && status == 1);

	if (fracvol == fracvol_aux){
		//printf("The phases are not different:\n");
		return;
	}

	if (fracvol == 1.0 && fracvol_right == 0.0) {
		*aux = -1;
	} else {
		*aux = 1;
	}
}


//void calculate_interfacial_force(sim_domain *sdp, higflow_solver *ns, int clid, Point center, Point IF){
	/////////////////////////////////////////////////////////////
	//// Set the curvature in the distributed curvature property
	//real curvature = compute_value_at_point(sdp, center, center, 1.0, ns->ed.mult.dpcurvature, ns->ed.stn);
	//real Normalx = compute_value_at_point(sdp, center, center, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
	//real Normaly = compute_value_at_point(sdp, center, center, 1.0,	ns->ed.mult.dpnormal[1], ns->ed.stn);
	//IF[0] = curvature * Normalx;
	//IF[1] = curvature * Normaly;

////	arquivoCurv(NULL,center[0],center[1],curvature);
	//for (int i = 0; i < DIM; i++) {
		//dp_set_value(ns->ed.mult.dpIF[i], clid, IF[i]);
	//}

//}

void higflow_compute_curvature_and_normal_with_HF_2D(higflow_solver *ns) {
	if (ns->contr.flowtype == 2) {
		real IF[DIM];
		// Get the local sub-domain for the cells
		sim_domain *sdp = psd_get_local_domain(ns->ed.mult.psdmult);

		// Get the map for the domain properties
		mp_mapper *mp = sd_get_domain_mapper(sdp);

		// Loop for each cell
		higcit_celliterator *it;

		for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
			// Get the cell
			hig_cell *c = higcit_getcell(it);

			// Get the cell identifier
			int clid = mp_lookup(mp, hig_get_cid(c));

			// Get the center of the cell
			Point center;
			hig_get_center(c, center);

			// Get the delta of the cell
			Point delta;
			hig_get_delta(c, delta);

			// Case bi-dimensional
			Point p;

			p[0] = center[0];
			p[1] = center[1];
			real frac = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			

//			arquivoFrac(NULL,p[0],p[1],frac);
//			continue;

			dp_set_value(ns->ed.mult.dpcurvature, clid, 0.0);/* Set the curvature in the distributed curvature property*/
			//for (int i = 0; i < DIM; i++) {
				//dp_set_value(ns->ed.mult.dpnormal[i], clid, 0.0);
			//}

			if (frac == 0.0 || frac == 1.0){
				continue;
			}

			//normal calculation ===========================
			//shirani_9_cells(sdp, ns, clid, center, p, delta);
			//Point Normal;
			//Normal[0] = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
			//Normal[1] = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
			
			//if(fabs(Normal[0])>fabs(Normal[1])){
				//calculate_exact_normal_x_dominant_2D(ns, clid, p);
			//}else{
				//calculate_exact_normal_y_dominant_2D(ns, clid, p);
			//}
			
			int aux_dir;
			real hf_horizontal = direction_hf_horizontal(sdp, ns, center, delta, &aux_dir); //height deriction verification
			real hf_vertical = direction_hf_vertical(sdp, ns, center, delta, &aux_dir); //height deriction verification
			
			real hm, hb, ht, vm, vl, vr;
			int aux_mh, aux_b, aux_t, aux_mv, aux_l, aux_r;
			if(hf_horizontal > hf_vertical){
				// Middle
	//			printf("Horizontal: going middle\n");
				horizontal_row_padrao(sdp, ns, center, p, delta, &hm, &aux_mh);
				// Top
				p[1] = center[1] + delta[1];
	//			printf("Horizontal: going top\n");
				horizontal_row_padrao(sdp, ns, center, p, delta, &ht, &aux_t);
				// Bottom
				p[1] = center[1] - delta[1];
	//			printf("Horizontal: going BOTTOM\n");
				horizontal_row_padrao(sdp, ns, center, p, delta, &hb, &aux_b);
				
				real H[3];
				H[0] = hb;
				H[1] = hm;
				H[2] = ht;
				//calculate_normal_cell_central_2nd_order_finite_difference_Horizontal(ns,clid,H[0],H[1],H[2],delta[0],delta[1],aux_mh);
				calculate_curvature_cell_central_2nd_order_finite_difference_Horizontal(ns,clid,H[0],H[1],H[2],delta[0],delta[1],aux_mh);
				//p[0] = center[0];
				//p[1] = center[1];
				//calculate_exact_normal_x_dominant_2D(ns, clid, p, aux_mh);
				//shirani_9_cells(sdp, ns, clid, center, p, delta);
			}else{
				// Middle
	//			printf("Vertical: going middle\n");
				vertical_collumn_padrao(sdp, ns, center, p, delta, &vm, &aux_mv);
				// Right
				p[0] = center[0] + delta[0];
				vertical_collumn_padrao(sdp, ns, center, p, delta, &vr, &aux_r);
	//			printf("Vertical: RIGTH: vr=%lf auxvr=%d origvr=%lf\n",vr,aux_r,orig_r);
				// Left
				p[0] = center[0] - delta[0];
				vertical_collumn_padrao(sdp, ns, center, p, delta, &vl, &aux_l);
				
				real V[3];
				V[0] = vl;
				V[1] = vm;
				V[2] = vr;
				//calculate_normal_cell_central_2nd_order_finite_difference_Vertical(ns, clid, V[0], V[1], V[2], delta[0], delta[1], aux_mv);
				calculate_curvature_cell_central_2nd_order_finite_difference_Vertical(ns, clid, V[0], V[1], V[2], delta[0], delta[1], aux_mv);
			}
			
		}
		// Destroy the iterator
		higcit_destroy(it);
		// Sync the distributed pressure property
		dp_sync(ns->ed.mult.dpcurvature);
		//for (int i = 0; i < DIM; i++) {
			//dp_sync(ns->ed.mult.dpnormal[i]);
		//}
	}
}
