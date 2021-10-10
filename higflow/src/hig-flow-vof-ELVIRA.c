#include "hig-flow-vof-finite-difference-normal-curvature.h"
#include "hig-flow-vof-adap-hf.h"
#include "hig-flow-vof-ELVIRA.h"
#include "hig-flow-vof-plic.h"

void ELVIRA_vertical_collumn(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *vertical, int *aux){
	int status;
	real fracvol, fracvol_aux;
	Point pp, ppt, ppb;
	*aux=0;
	// Vertical Up
	status = get_frac_vol(sdp, ns, 1, center, p, delta, &fracvol);
	if (status != 1) {
		if (status == -1) {
			printf("Up - Vertical cells with different sizes \n");
		} else {
			//printf("Up - Vertical cell out of domain \n");
		}
		return;
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
			return;
		} else {
			printf("Up - Vertical cells with different sizes \n");
			return;
		}
	} while (i <2 && status == 1);
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
			return;
		} else {
			printf("Down - Vertical cells with different sizes \n");
			return;
		}
	} while (i <2 && status == 1);
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

void ELVIRA_horizontal_row(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *horizontal, int *aux){
	int status;
	real fracvol, fracvol_aux;
	Point pp, ppr, ppl;
	*aux=0;
	// Horizontal Right
	status = get_frac_vol(sdp, ns, 0, center, p, delta, &fracvol);
	if (status != 1) {
		if (status == -1) {
			printf("Right - Horizontal cells with different sizes \n");
		} else {
			//printf("Right - Horizontal cell out of domain \n");
		}
		return;
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
			return;
		} else {
			printf("Right - Horizontal cells with different sizes \n");
			return;
		}
	} while (i <2 && status == 1);
	fracvol_aux = fracvol;
	ppr[0]=pp[0]+0.5*delta[0];
	ppr[1]=pp[1];
	real fracvol_right = fracvol;
	// Horizontal Left
	i = 0;
	// fracvol_aux = 0;
	do {
		i++;
		pp[0] = p[0] - i*delta[0];
		status = get_frac_vol(sdp, ns, 0, center, pp, delta, &fracvol);
		if (status == 1) {
			*horizontal += fracvol;
		} else if (status == 0) {
			//printf("Left - Horizontal cell out of domain \n");
			return;
		} else {
			printf("Left - Horizontal cells with different sizes \n");
			return;
		}
	} while (i <2 && status == 1);
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

void higflow_compute_curvature_interfacial_force_normal_multiphase_2D_EL(higflow_solver *ns) {
	if (ns->contr.flowtype == 2) {
		real IF[DIM];
		// Get the local sub-domain for the cells
		sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
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
			for (int i = 0; i < DIM; i++) {
				dp_set_value(ns->ed.mult.dpIF[i], clid, 0.0);
				dp_set_value(ns->ed.mult.dpnormal[i], clid, 0.0);
			}
			if (frac == 0.0 || frac == 1.0){
				continue;
			}
			
			real area = frac*delta[0]*delta[1];
			
			//normal calculation ======================================
			real hm, hb, ht, vm, vl, vr;
			int aux_mh, aux_b, aux_t, aux_mv, aux_l, aux_r;
			
			// Middle
//			printf("Horizontal: going middle\n");
			ELVIRA_horizontal_row(sdp, ns, center, p, delta, &hm, &aux_mh);
			// Top
			p[1] = center[1] + delta[1];
//			printf("Horizontal: going top\n");
			ELVIRA_horizontal_row(sdp, ns, center, p, delta, &ht, &aux_t);
			// Bottom
			p[1] = center[1] - delta[1];
//			printf("Horizontal: going BOTTOM\n");
			ELVIRA_horizontal_row(sdp, ns, center, p, delta, &hb, &aux_b);
			
			real H[3];
			H[0] = hb;
			H[1] = hm;
			H[2] = ht;
			
			//==========================================================
			////horizontal-regressive finite difference
			ELVIRA_calculate_normal_cell_regressive_1st_order_finite_difference_Horizontal(ns, clid, hm, hb, delta[0], delta[1], aux_mh);
			
			Point Normal;
			Normal[0] = compute_value_at_point(sdp, center, center, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
			Normal[1] = compute_value_at_point(sdp, center, center, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
			
			real distance1 = distance_from_center(Normal,delta,area);
			
			//real area1 = area_left_line_origin_center(Normal, delta, 0.010217);
			real area1 = area_left_line_origin_center(Normal, delta, distance1);
			
			real erro1 = sqrt(pow(area1-area, 2));
			
			printf("distance1 = %.16lf\t", distance1);
			printf("area1       = %.16lf\n", area1);
			printf("erro1       = %.16lf\n", erro1);
			//==========================================================
			////horizontal-central finite difference
			ELVIRA_calculate_normal_cell_central_2nd_order_finite_difference_Horizontal
			(ns, clid, ht, hb, delta[0], delta[1], aux_mh);
			
			Normal[0] = compute_value_at_point(sdp, center, center, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
			Normal[1] = compute_value_at_point(sdp, center, center, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
			
			real distance2 = distance_from_center(Normal,delta,area);
			
			//real area2 = area_left_line_origin_center(Normal, delta, 0.010040);
			real area2 = area_left_line_origin_center(Normal, delta, distance2);
			
			real erro2 = sqrt(pow(area2-area, 2));
			
			printf("distance2 = %.16lf\t", distance2);
			printf("area2       = %.16lf\n", area2);
			printf("erro2       = %.16lf\n", erro2);
			//==========================================================
			
			////horizontal-progressive finite difference
			ELVIRA_calculate_normal_cell_progressive_1st_order_finite_difference_Horizontal
			(ns, clid, ht, hm, delta[0], delta[1], aux_mh);
			
			Normal[0] = compute_value_at_point(sdp, center, center, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
			Normal[1] = compute_value_at_point(sdp, center, center, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
			
			real distance3 = distance_from_center(Normal,delta,area);
			
			//real area3 = area_left_line_origin_center(Normal, delta, 0.009830);
			real area3 = area_left_line_origin_center(Normal, delta, distance3);
			
			real erro3 = sqrt(pow(area3-area, 2));
			
			printf("distance3 = %.16lf\t", distance3);
			printf("area3       = %.16lf\n", area3);
			printf("erro3       = %.16lf\n", erro3);
			getchar();
			
			//==========================================================
			p[0] = center[0];
			p[1] = center[1];
			// Middle
//			printf("Vertical: going middle\n");
			ELVIRA_vertical_collumn(sdp, ns, center, p, delta, &vm, &aux_mv);
			// Right
			p[0] = center[0] + delta[0];
			ELVIRA_vertical_collumn(sdp, ns, center, p, delta, &vr, &aux_r);
//			printf("Vertical: RIGTH: vr=%lf auxvr=%d origvr=%lf\n",vr,aux_r,orig_r);
			// Left
			p[0] = center[0] - delta[0];
			ELVIRA_vertical_collumn(sdp, ns, center, p, delta, &vl, &aux_l);
			
			real V[3];
			V[0] = vl;
			V[1] = vm;
			V[2] = vr;
			
			////vertical-progressive finite difference
			ELVIRA_calculate_normal_cell_progressive_1st_order_finite_difference_Vertical
			(ns, clid, vr, vm, delta[0], delta[1], aux_mv);
			
			Normal[0] = compute_value_at_point(sdp, center, center, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
			Normal[1] = compute_value_at_point(sdp, center, center, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
			
			real distance4 = distance_from_center(Normal,delta,area);
			//printf("%lf\t", distance4);
			
			real area4 = area_left_line_origin_center(Normal, delta, distance4);
			//printf("%.16lf\n", area4);
			
			real erro4 = sqrt(pow(area4-area, 2));
			//==========================================================
			
			////vertical-central finite difference
			ELVIRA_calculate_normal_cell_central_2nd_order_finite_difference_Vertical
			(ns, clid, vr, vl, delta[0], delta[1], aux_mv);
			
			Normal[0] = compute_value_at_point(sdp, center, center, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
			Normal[1] = compute_value_at_point(sdp, center, center, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
			
			real distance5 = distance_from_center(Normal,delta,area);
			//printf("%lf\t", distance5);
			
			//real area1 = area_left_line_origin_center(Normal, delta, 0.010217);
			real area5 = area_left_line_origin_center(Normal, delta, distance5);
			//printf("%.16lf\n", area5);
			
			real erro5 = sqrt(pow(area5-area, 2));
			//==========================================================
			
			////vertical-regressive finite difference
			ELVIRA_calculate_normal_cell_regressive_1st_order_finite_difference_Vertical
			(ns, clid, vm, vl, delta[0], delta[1], aux_mv);
			
			Normal[0] = compute_value_at_point(sdp, center, center, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
			Normal[1] = compute_value_at_point(sdp, center, center, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
			
			real distance6 = distance_from_center(Normal,delta,area);
			//printf("%lf\t", distance6);
			
			//real area1 = area_left_line_origin_center(Normal, delta, 0.010217);
			real area6 = area_left_line_origin_center(Normal, delta, distance6);
			//printf("%.16lf\n", area6);
			//printf("\n");
			
			real erro6 = sqrt(pow(area6-area, 2));
			//==========================================================
			
			p[0] = center[0];
			p[1] = center[1];
			
			if(fabs(erro1)<fabs(erro2) && fabs(erro1)<fabs(erro3) && fabs(erro1)<fabs(erro4)
			&& fabs(erro1)<fabs(erro5) && fabs(erro1)<fabs(erro6)){
				ELVIRA_calculate_normal_cell_regressive_1st_order_finite_difference_Horizontal
				(ns, clid, hm, hb, delta[0], delta[1], aux_mh);
			} else if(fabs(erro2)<fabs(erro1) && fabs(erro2)<fabs(erro3) && fabs(erro2)<fabs(erro4)
			&& fabs(erro2)<fabs(erro5) && fabs(erro2)<fabs(erro6)){
				ELVIRA_calculate_normal_cell_central_2nd_order_finite_difference_Horizontal
				(ns, clid, ht, hb, delta[0], delta[1], aux_mh);
			} else if(fabs(erro3)<fabs(erro1) && fabs(erro3)<fabs(erro2) && fabs(erro3)<fabs(erro4)
			&& fabs(erro3)<fabs(erro5) && fabs(erro3)<fabs(erro6)){
				ELVIRA_calculate_normal_cell_progressive_1st_order_finite_difference_Horizontal
				(ns, clid, ht, hm, delta[0], delta[1], aux_mh);
			} else if(fabs(erro4)<fabs(erro1) && fabs(erro4)<fabs(erro2) && fabs(erro4)<fabs(erro3)
			&& fabs(erro4)<fabs(erro5) && fabs(erro4)<fabs(erro6)){
				ELVIRA_calculate_normal_cell_progressive_1st_order_finite_difference_Vertical
			(ns, clid, vr, vm, delta[0], delta[1], aux_mv);
			} else if(fabs(erro5)<fabs(erro1) && fabs(erro5)<fabs(erro2) && fabs(erro5)<fabs(erro3)
			&& fabs(erro5)<fabs(erro4) && fabs(erro5)<fabs(erro6)){
				ELVIRA_calculate_normal_cell_central_2nd_order_finite_difference_Vertical
			(ns, clid, vr, vl, delta[0], delta[1], aux_mv);
			} else{
				ELVIRA_calculate_normal_cell_regressive_1st_order_finite_difference_Vertical
			(ns, clid, vm, vl, delta[0], delta[1], aux_mv);
			}
			
			
		}
		// Destroy the iterator
		higcit_destroy(it);
		// Sync the distributed pressure property
		dp_sync(ns->ed.mult.dpcurvature);
		for (int i = 0; i < DIM; i++) {
//			dp_sync(ns->ed.mult.dpIF[i]);
			dp_sync(ns->ed.mult.dpnormal[i]);
		}
	}
}
