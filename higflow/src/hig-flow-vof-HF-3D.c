#if DIM == 3
#include "hig-flow-vof-finite-difference-normal-curvature_3D.h"
#include "hig-flow-vof-HF-3D.h"

void fraction_correction_at_get_3D(real *fracvol){
	if(fabs(*fracvol - 1.0)<1.0e-8){
		*fracvol = 1.0;
	}else if(fabs(*fracvol)<1.0e-8){
		*fracvol=0.0;
	}
	
}
//======================================================================
int get_frac_vol_3D(sim_domain *sd, higflow_solver *ns, int dim, 
	Point Center, Point P,Point Delta,real *fracvol){
	hig_cell *c = sd_get_cell_with_point(sd,P);
	if (c == NULL){
		return 0;
	} else {
		Point Delta2;
		hig_get_delta(c, Delta2);
		
		if (fabs(Delta[dim] - Delta2[dim])>1.0e-12){
			printf("Different sizes dc=[%.18lf %.18lf %.18lf] dp= [%.18lf %.18lf %.18lf] \t ",
			Delta[0],Delta[1],Delta[2],Delta2[0],Delta2[1],Delta2[2]);
			return -1;
		} else {
			*fracvol=compute_value_at_point(sd, Center, P, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			fraction_correction_at_get_3D(fracvol);
			return 1;
		}
	}
}

//======================================================================
//verifying which is the dominant deriction=============================
//======================================================================
real direction_hf_x(sim_domain *sdp, higflow_solver *ns, Point center, Point delta){
	real frac, fracvol;
	Point p;
	real hf_x=0.0;
	
	p[0] = center[0]+delta[0];
	real hf_hor_dir = 0.0;
	int k = 0;
	do {
		p[2] = center[2] + (k*delta[2]-delta[2]);
		int i = 0;
		do {
			p[1] = center[1] + (i*delta[1]-delta[1]);
			real frac = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			hf_hor_dir += frac;
			i++;
		} while (i < 3);
		k++;
	} while (k < 3);
	
	p[0] = center[0]-delta[0];
	real hf_hor_esq = 0.0;
	k = 0;
	do{
		p[2] = center[2] + (k*delta[2]-delta[2]);
		int i = 0;
		do {
			p[1] = center[1] + (i*delta[1]-delta[1]);
			frac = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			hf_hor_esq += frac;
			i++;
		} while (i < 3);
		k++;
	}while (k < 3);
	
	hf_x = fabs(hf_hor_dir - hf_hor_esq);
	return hf_x;
}
//======================================================================
real direction_hf_y(sim_domain *sdp, higflow_solver *ns, Point center, Point delta){
	real frac, fracvol;
	Point p;
	real hf_y=0.0;
	
	p[1] = center[1]+delta[1];
	real hf_ver_sup = 0.0;
	int k = 0;
	do{
		p[2] = center[2] + (k*delta[2]-delta[2]);
		int i = 0;
		do {
			p[0] = center[0] + (i*delta[0]-delta[0]);
			real frac = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			hf_ver_sup += frac;
			i++;
		} while (i < 3);
		k++;
	} while (k < 3);
		
	p[1] = center[1]-delta[1];
	real hf_ver_inf = 0.0;
	k = 0;
	do{
		p[2] = center[2] + (k*delta[2]-delta[2]);
		int i = 0;
		do {
			p[0] = center[0] + (i*delta[0]-delta[0]);
			frac = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			hf_ver_inf += frac;
			i++;
		} while (i < 3);
		k++;
	} while (k < 3);
	
	hf_y = fabs(hf_ver_sup - hf_ver_inf);
	return hf_y;
}
//======================================================================
real direction_hf_z(sim_domain *sdp, higflow_solver *ns, Point center, Point delta){
	real frac, fracvol;
	Point p;
	real hf_z=0.0;
	
	p[2] = center[2]+delta[2];
	real hf_z_sup = 0.0;
	int k = 0;
	do{
		p[1] = center[1] + (k*delta[1]-delta[1]);
		int i = 0;
		do {
			p[0] = center[0] + (i*delta[0]-delta[0]);
			real frac = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			hf_z_sup += frac;
			i++;
		} while (i < 3);
		k++;
	} while (k < 3);
		
	p[2] = center[2]-delta[2];
	real hf_z_inf = 0.0;
	k = 0;
	do{
		p[1] = center[1] + (k*delta[1]-delta[1]);
		int i = 0;
		do {
			p[0] = center[0] + (i*delta[0]-delta[0]);
			frac = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			hf_z_inf += frac;
			i++;
		} while (i < 3);
		k++;
	} while (k < 3);
	
	hf_z = fabs(hf_z_sup - hf_z_inf);
	//printf("%lf %lf %lf\n",hf_ver_dir, hf_ver_esq, hf_z);
	return hf_z;
}
//======================================================================
//arranging the HF stencil==============================================
//======================================================================

void x_block(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *x_row, int *aux){
	int status;
	real fracvol, fracvol_aux;
	Point pp, ppr, ppl;
	*aux=0;
	// x_row Right
	status = get_frac_vol_3D(sdp, ns, 0, center, p, delta, &fracvol);
	if (status != 1) {
		if (status == -1) {
			printf("Right - Horizontal cells with different sizes \n");
		} else {
			//printf("Right - Horizontal cell out of domain \n");
			fracvol = 0.0;
		}
		return;
	}
	pp[1] = p[1];
	pp[2] = p[2];
	*x_row=fracvol;
	int i = 0;
	do {
		i++;
		pp[0] = center[0] + i*delta[0];
		status=get_frac_vol_3D(sdp, ns, 0, center, pp, delta, &fracvol);
		if (status == 1) {
			*x_row += fracvol;
		} else if (status == 0) {
			//printf("Right - Horizontal cell out of domain \n");
			fracvol = 0.0;
			//return;
		} else {
			printf("Right - Horizontal cells with different sizes \n");
			return;
		}
	} while (i <4 && status == 1);
	fracvol_aux = fracvol;
	//ppr[0]=pp[0]+0.5*delta[0];
	//ppr[1]=pp[1];
	real fracvol_right = fracvol;
	// x_row Left
	i = 0;
	//fracvol_aux=0;
	do {
		i++;
		pp[0] = p[0] - i*delta[0];
		status = get_frac_vol_3D(sdp, ns, 0, center, pp, delta, &fracvol);
		if (status == 1) {
			*x_row += fracvol;
		} else if (status == 0) {
			//printf("Left - Horizontal cell out of domain \n");
			fracvol = 0.0;
			//return;
		} else {
			printf("Left - Horizontal cells with different sizes \n");
			return;
		}
	}while (i <4 && status == 1);

	//if (fracvol == fracvol_aux){
		////printf("The phases are not different:\n");
		//return;
	//}

	if (fracvol == 1.0 && fracvol_right == 0.0) {
		*aux = -1;
	} else {
		*aux = 1;
	}
}
//======================================================================

void y_block(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *y_row, int *aux){
	int status;
	real fracvol, fracvol_aux;
	Point pp, ppr, ppl;
	*aux=0;
	// x_row Right
	status = get_frac_vol_3D(sdp, ns, 0, center, p, delta, &fracvol);
	if (status != 1) {
		if (status == -1) {
			printf("Right - Horizontal cells with different sizes \n");
		} else {
			//printf("Right - Horizontal cell out of domain \n");
			fracvol = 0.0;
		}
		return;
	}
	pp[0] = p[0];
	pp[2] = p[2];
	*y_row=fracvol;
	int i = 0;
	do {
		i++;
		pp[1] = center[1] + i*delta[1];
		status=get_frac_vol_3D(sdp, ns, 0, center, pp, delta, &fracvol);
		if (status == 1) {
			*y_row += fracvol;
		} else if (status == 0) {
			//printf("Right - Horizontal cell out of domain \n");
			fracvol = 0.0;
			//return;
		} else {
			printf("Right - Horizontal cells with different sizes \n");
			return;
		}
	} while (i <4 && status == 1);
	fracvol_aux = fracvol;
	//ppr[0]=pp[0]+0.5*delta[0];
	//ppr[1]=pp[1];
	real fracvol_right = fracvol;
	// x_row Left
	i = 0;
	//fracvol_aux=0;
	do {
		i++;
		pp[1] = p[1] - i*delta[1];
		status = get_frac_vol_3D(sdp, ns, 0, center, pp, delta, &fracvol);
		if (status == 1) {
			*y_row += fracvol;
		} else if (status == 0) {
			//printf("Left - Horizontal cell out of domain \n");
			fracvol = 0.0;
			//return;
		} else {
			printf("Left - Horizontal cells with different sizes \n");
			return;
		}
	}while (i <4 && status == 1);

	//if (fracvol == fracvol_aux){
		////printf("The phases are not different:\n");
		//return;
	//}

	if (fracvol == 1.0 && fracvol_right == 0.0) {
		*aux = -1;
	} else {
		*aux = 1;
	}
}
//======================================================================

void z_block(sim_domain *sdp, higflow_solver *ns, Point center, Point p, Point delta, real *z_row, int *aux){
	int status;
	real fracvol, fracvol_aux;
	Point pp, ppr, ppl;
	*aux=0;
	// z_row upper
	status = get_frac_vol_3D(sdp, ns, 0, center, p, delta, &fracvol);
	if (status != 1) {
		if (status == -1) {
			printf("Right - Horizontal cells with different sizes \n");
		} else {
			//printf("Right - Horizontal cell out of domain \n");
			fracvol = 0.0;
		}
		return;
	}
	pp[1] = p[1];
	pp[0] = p[0];
	*z_row=fracvol;
	int i = 0;
	do {
		i++;
		pp[2] = center[2] + i*delta[2];
		status=get_frac_vol_3D(sdp, ns, 0, center, pp, delta, &fracvol);
		if (status == 1) {
			*z_row += fracvol;
		} else if (status == 0) {
			//printf("Right - Horizontal cell out of domain \n");
			fracvol = 0.0;
			//return;
		} else {
			printf("Right - Horizontal cells with different sizes \n");
			return;
		}
	} while (i <4 && status == 1);
	fracvol_aux = fracvol;
	//ppr[0]=pp[0]+0.5*delta[0];
	//ppr[1]=pp[1];
	real fracvol_right = fracvol;
	// z_row bottom
	i = 0;
	//fracvol_aux=0;
	do {
		i++;
		pp[2] = p[2] - i*delta[2];
		status = get_frac_vol_3D(sdp, ns, 0, center, pp, delta, &fracvol);
		if (status == 1) {
			*z_row += fracvol;
		} else if (status == 0) {
			//printf("Left - Horizontal cell out of domain \n");
			fracvol = 0.0;
			//return;
		} else {
			printf("Left - Horizontal cells with different sizes \n");
			return;
		}
	}while (i <4 && status == 1);

	//if (fracvol == fracvol_aux){
		////printf("The phases are not different:\n");
		//return;
	//}

	if (fracvol == 1.0 && fracvol_right == 0.0) {
		*aux = -1;
	} else {
		*aux = 1;
	}
}
//======================================================================
//======================================================================

void higflow_compute_curvature_interfacial_force_normal_multiphase_3D_HF_padrao(higflow_solver *ns) {
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
			p[2] = center[2];
			real frac = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			

//			arquivoFrac(NULL,p[0],p[1],frac);
//			continue;

			dp_set_value(ns->ed.mult.dpcurvature, clid, 0.0);/* Set the curvature in the distributed curvature property*/
			//for (int i = 0; i < DIM; i++) {
				////dp_set_value(ns->ed.mult.dpIF[i], clid, 0.0);
				//dp_set_value(ns->ed.mult.dpnormal[i], clid, 0.0);
			//}

			//if (frac < 1.0e-12 || fabs(frac - 1.0) < 1.0e-12){
			if (frac == 0.0 || frac == 1.0){
				continue;
			}


			//printf("%lf %lf %lf %lf\n",p[0],p[1],p[2], frac);
			
			//shirani_125_cells(sdp, ns, clid, center, p, delta);
			//shirani_125_cells(ns);
			//Point Normal;
			//Normal[0] = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
			//Normal[1] = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
			//Normal[2] = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpnormal[2], ns->ed.stn);
			//printf("%lf %lf %lf\n", Normal[0], Normal[1], Normal[2]);
			//shirani_125_cells_curvature(sdp, ns, clid, center, p, delta);
			
			

			real hf_x = direction_hf_x(sdp, ns, center, delta); //height deriction verification
			real hf_y = direction_hf_y(sdp, ns, center, delta); //height deriction verification
			real hf_z = direction_hf_z(sdp, ns, center, delta); //height deriction verification
			//printf("%lf %lf %lf\n", hf_x, hf_y, hf_z)
			
			real hxm, hxu, hxb, hxr, hxl, hxbl, hxbr, hxul,hxur;
			int  aux_xm, aux_xu, aux_xb, aux_xr, aux_xl, aux_xbl, 
			aux_xbr, aux_xur, aux_xul;
			
			real hym, hyu, hyb, hyr, hyl, hybl, hybr, hyul,hyur;
			int  aux_ym, aux_yu, aux_yb, aux_yr, aux_yl, aux_ybl, 
			aux_ybr, aux_yur, aux_yul;
			
			real hzm, hzu, hzb, hzr, hzl, hzbl, hzbr, hzul,hzur;
			int  aux_zm, aux_zu, aux_zb, aux_zr, aux_zl, aux_zbl, 
			aux_zbr, aux_zur, aux_zul;
			
			if((hf_x > hf_y) && (hf_x > hf_z)){
			//if(Normal[0]>Normal[1] && Normal[0]>Normal[2]){
				// Middle
				//=========================
				x_block(sdp, ns, center, p, delta, &hxm, &aux_xm);
				// Up
				//=========================
				p[2] = center[2] + delta[2];
				x_block(sdp, ns, center, p, delta, &hxu, &aux_xu);
				// Bottom
				//=========================
				p[2] = center[2] - delta[2];
				x_block(sdp, ns, center, p, delta, &hxb, &aux_xb);
				// Right
				//=========================
				p[1] = center[1] + delta[1];
				p[2] = center[2];
				x_block(sdp, ns, center, p, delta, &hxr, &aux_xr);
				// Left
				//=========================
				p[1] = center[1] - delta[1];
				x_block(sdp, ns, center, p, delta, &hxl, &aux_xl);
				//Mix stencil rows
				// Up-Right
				//=========================
				p[1] = center[1] + delta[1];
				p[2] = center[2] + delta[2];
				x_block(sdp, ns, center, p, delta, &hxur, &aux_xur);
				// Up-Left
				//=========================
				p[1] = center[1] - delta[1];
				p[2] = center[2] + delta[2];
				x_block(sdp, ns, center, p, delta, &hxul, &aux_xul);
				//=========================
				// Bottom-Right
				p[1] = center[1] + delta[1];
				p[2] = center[2] - delta[2];
				x_block(sdp, ns, center, p, delta, &hxbr, &aux_xbr);
				//=========================
				// Bottom-Left
				p[1] = center[1] - delta[1];
				p[2] = center[2] - delta[2];
				x_block(sdp, ns, center, p, delta, &hxbl, &aux_xbl);
	
				real H[9];
				H[0] = hxl;
				H[1] = hxr;
				H[2] = hxb;
				H[3] = hxu;
				H[4] = hxm;
				H[5] = hxul;
				H[6] = hxur;
				H[7] = hxbl;
				H[8] = hxbr;
				//Curvature
				calculate_HF_curvature_x_dominant(ns,clid,H[0],H[1],H[2],H[3],H[4],
				H[5],H[6],H[7],H[8],delta[0],delta[1],delta[2],aux_xm);
				//Normal vector
				//calculate_HF_normal_x_dominant(ns,clid,H[0],H[1],H[2],H[3],delta[0],delta[1],delta[2],aux_xm);
				////Theoretical normal
				//p[0] = center[0];
				//p[1] = center[1];
				//p[2] = center[2];
				//calculate_exact_normal_x_dominant(ns, clid,  p, aux_xm);
				//printf("%lf %lf %lf \n",p[0],p[1],p[2]);
			}
			else if((hf_y > hf_x) && (hf_y > hf_z)){
			//else if(Normal[1]>Normal[0] && Normal[1]>Normal[2]){
				// Middle
				//=========================
				y_block(sdp, ns, center, p, delta, &hym, &aux_ym);
				// Up
				//=========================
				p[2] = center[2] + delta[2];
				y_block(sdp, ns, center, p, delta, &hyu, &aux_yu);
				// Bottom
				//=========================
				p[2] = center[2] - delta[2];
				y_block(sdp, ns, center, p, delta, &hyb, &aux_yb);
				// Right
				//=========================
				p[0] = center[0] + delta[0];
				p[2] = center[2];
				y_block(sdp, ns, center, p, delta, &hyr, &aux_yr);
				// Left
				//=========================
				p[0] = center[0] - delta[0];
				y_block(sdp, ns, center, p, delta, &hyl, &aux_yl);
				//Mix stencil rows
				// Up-Right
				//=========================
				p[0] = center[0] + delta[0];
				p[2] = center[2] + delta[2];
				y_block(sdp, ns, center, p, delta, &hyur, &aux_yur);
				// Up-Left
				//=========================
				p[0] = center[0] - delta[0];
				p[2] = center[2] + delta[2];
				y_block(sdp, ns, center, p, delta, &hyul, &aux_yul);
				//=========================
				// Bottom-Right
				p[0] = center[0] + delta[0];
				p[2] = center[2] - delta[2];
				y_block(sdp, ns, center, p, delta, &hybr, &aux_ybr);
				//=========================
				// Bottom-Left
				p[0] = center[0] - delta[0];
				p[2] = center[2] - delta[2];
				y_block(sdp, ns, center, p, delta, &hybl, &aux_ybl);
	
				real H[9];
				H[0] = hyl;
				H[1] = hyr;
				H[2] = hyb;
				H[3] = hyu;
				H[4] = hym;
				H[5] = hyul;
				H[6] = hyur;
				H[7] = hybl;
				H[8] = hybr;
				//Curvature
				calculate_HF_curvature_y_dominant(ns,clid,H[0],H[1],H[2],H[3],H[4],
				H[5],H[6],H[7],H[8],delta[0],delta[1],delta[2],aux_ym);
				//Normal vector
				//calculate_HF_normal_y_dominant(ns,clid,H[0],H[1],H[2],H[3],delta[0],delta[1],delta[2],aux_ym);
				////Theoretical normal
				//p[0] = center[0];
				//p[1] = center[1];
				//p[2] = center[2];
				//calculate_exact_normal_y_dominant(ns, clid,  p, aux_ym);
				//printf("%lf %lf %lf \n",p[0],p[1],p[2]);
			}
			else{
			//else if((hf_z > (hf_x)) && (hf_z > (hf_y))){
			//else if(Normal[2]>Normal[0] && Normal[2]>Normal[1]){
				// Middle
				//=========================
				z_block(sdp, ns, center, p, delta, &hzm, &aux_zm);
				// Up
				//=========================
				p[0] = center[0] + delta[0];
				z_block(sdp, ns, center, p, delta, &hzu, &aux_zu);
				// Bottom
				//=========================
				p[0] = center[0] - delta[0];
				z_block(sdp, ns, center, p, delta, &hzb, &aux_zb);
				// Right
				//=========================
				p[1] = center[1] + delta[1];
				p[0] = center[0];
				z_block(sdp, ns, center, p, delta, &hzr, &aux_zr);
				// Left
				//=========================
				p[1] = center[1] - delta[1];
				z_block(sdp, ns, center, p, delta, &hzl, &aux_zl);
				//Mix stencil rows
				// Up-Right
				//=========================
				p[1] = center[1] + delta[1];
				p[0] = center[0] + delta[0];
				z_block(sdp, ns, center, p, delta, &hzur, &aux_zur);
				// Up-Left
				//=========================
				p[1] = center[1] - delta[1];
				p[0] = center[0] + delta[0];
				z_block(sdp, ns, center, p, delta, &hzul, &aux_zul);
				//=========================
				// Bottom-Right
				p[1] = center[1] + delta[1];
				p[0] = center[0] - delta[0];
				z_block(sdp, ns, center, p, delta, &hzbr, &aux_zbr);
				//=========================
				// Bottom-Left
				p[1] = center[1] - delta[1];
				p[0] = center[0] - delta[0];
				z_block(sdp, ns, center, p, delta, &hzbl, &aux_zbl);
	
				real H[9];
				H[0] = hzl;
				H[1] = hzr;
				H[2] = hzb;
				H[3] = hzu;
				H[4] = hzm;
				H[5] = hzul;
				H[6] = hzur;
				H[7] = hzbl;
				H[8] = hzbr;
				//Curvature
				calculate_HF_curvature_z_dominant(ns,clid,H[0],H[1],H[2],H[3],H[4],
				H[5],H[6],H[7],H[8],delta[0],delta[1],delta[2],aux_zm);
				//Normal vector
				//calculate_HF_normal_z_dominant(ns,clid,H[0],H[1],H[2],H[3],delta[0],delta[1],delta[2],aux_zm);
				////Theoretical normal
				//p[0] = center[0];
				//p[1] = center[1];
				//p[2] = center[2];
				//calculate_exact_normal_z_dominant(ns, clid,  p, aux_zm);
				//printf("%lf %lf %lf \n",p[0],p[1],p[2]);
			}
		}
		// Destroy the iterator
		higcit_destroy(it);
		// Sync the distributed pressure property
		dp_sync(ns->ed.mult.dpcurvature);
	}
}
#endif
