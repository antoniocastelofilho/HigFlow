#if DIM == 3
//#include "hig-flow-vof-finite-difference-normal-curvature.h"
#include "hig-flow-vof-finite-difference-normal-curvature_3D.h"
//#include "hig-flow-vof-HF-3D_adap.h"
//#include "hig-flow-vof-HF-3D.h"

int cell_location_2(sim_domain *sd, higflow_solver *ns, int dim, Point Center, Point P){
	hig_cell *c = sd_get_cell_with_point(sd,P);
	if (c == NULL){
		return 0;
	} else {
		return 1;
	}
}

void shirani_125_cells(higflow_solver *ns){
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

			Point p;
			p[0] = center[0];
			p[1] = center[1];
			p[2] = center[2];
			real frac = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			

//			arquivoFrac(NULL,p[0],p[1],frac);
//			continue;

			for (int i = 0; i < DIM; i++) {
				dp_set_value(ns->ed.mult.dpnormal[i], clid, 0.0);
			}

			if (frac < 1.0e-10 || fabs(frac - 1.0) < 1.0e-10){
			//if (frac == 0.0 || frac == 1.0){
				continue;
			}
		
			Point pp, Normal;
			
			real fracvol, gradFij1, gradFij2, gradFij3, gradFij, gradFik1, gradFik2,
				gradFik3, gradFik, gradFji1, gradFji2, gradFji3, gradFji, gradFjk1,
				gradFjk2, gradFjk3, gradFjk, gradFki1, gradFki2, gradFki3, gradFki,
				gradFkj1, gradFkj2, gradFkj3, gradFkj, gradFi, gradFj, gradFk;
				
			int status;
		
		
		
			pp[0]=p[0];
			pp[1]=p[1];
			pp[2]=p[2];
			real f_c = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			
			pp[0]=p[0]+delta[0];
			pp[1]=p[1]-delta[1];
			pp[2]=p[2]-delta[2];
			real f_fld = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_fld = 0.0; }
			
			pp[0]=p[0]-delta[0];
			pp[1]=p[1]-delta[1];
			pp[2]=p[2]-delta[2];
			real f_bld = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_bld = 0.0; }
			
			pp[0]=p[0]+delta[0];
			pp[1]=p[1];
			pp[2]=p[2]-delta[2];
			real f_fd = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_fd = 0.0; }
			
			pp[0]=p[0]-delta[0];
			pp[1]=p[1];
			pp[2]=p[2]-delta[2];
			real f_bd = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_bd = 0.0; }
			
			pp[0]=p[0]+delta[0];
			pp[1]=p[1]+delta[1];
			pp[2]=p[2]-delta[2];
			real f_frd = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_frd = 0.0; }
			
			pp[0]=p[0]-delta[0];
			pp[1]=p[1]+delta[1];
			pp[2]=p[2]-delta[2];
			real f_brd = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_brd = 0.0; }
			
			pp[0]=p[0]+delta[0];
			pp[1]=p[1]-delta[1];
			pp[2]=p[2];
			real f_fl = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_fl = 0.0; }
			
			pp[0]=p[0]-delta[0];
			pp[1]=p[1]-delta[1];
			pp[2]=p[2];
			real f_bl = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_bl = 0.0; }
			
			pp[0]=p[0]+delta[0];
			pp[1]=p[1];
			pp[2]=p[2];
			real f_f = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_f = 0.0; }
			
			pp[0]=p[0]-delta[0];
			pp[1]=p[1];
			pp[2]=p[2];
			real f_b = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_b = 0.0; }
			
			pp[0]=p[0]+delta[0];
			pp[1]=p[1]+delta[1];
			pp[2]=p[2];
			real f_fr = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_fr = 0.0; }
			
			pp[0]=p[0]-delta[0];
			pp[1]=p[1]+delta[1];
			pp[2]=p[2];
			real f_br = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_br = 0.0; }
			
			pp[0]=p[0]+delta[0];
			pp[1]=p[1]-delta[1];
			pp[2]=p[2]+delta[2];
			real f_flu = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_flu = 0.0; }
			
			pp[0]=p[0]-delta[0];
			pp[1]=p[1]-delta[1];
			pp[2]=p[2]+delta[2];
			real f_blu = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_blu = 0.0; }
			
			pp[0]=p[0]+delta[0];
			pp[1]=p[1];
			pp[2]=p[2]+delta[2];
			real f_fu = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_fu = 0.0; }
			
			pp[0]=p[0]-delta[0];
			pp[1]=p[1];
			pp[2]=p[2]+delta[2];
			real f_bu = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_bu = 0.0; }
			
			pp[0]=p[0]+delta[0];
			pp[1]=p[1]+delta[1];
			pp[2]=p[2]+delta[2];
			real f_fru = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_fru = 0.0; }
			
			pp[0]=p[0]-delta[0];
			pp[1]=p[1]+delta[1];
			pp[2]=p[2]+delta[2];
			real f_bru = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_bru = 0.0; }
			
			pp[0]=p[0];
			pp[1]=p[1]+delta[1];
			pp[2]=p[2]-delta[2];
			real f_rd = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_rd = 0.0; }
			
			pp[0]=p[0];
			pp[1]=p[1]-delta[1];
			pp[2]=p[2]-delta[2];
			real f_ld = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_ld = 0.0; }
			
			pp[0]=p[0];
			pp[1]=p[1]+delta[1];
			pp[2]=p[2];
			real f_r = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_r = 0.0; }
			
			pp[0]=p[0];
			pp[1]=p[1]-delta[1];
			pp[2]=p[2];
			real f_l = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_l = 0.0; }
			
			pp[0]=p[0];
			pp[1]=p[1]+delta[1];
			pp[2]=p[2]+delta[2];
			real f_ru = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_ru = 0.0; }
			
			pp[0]=p[0];
			pp[1]=p[1]-delta[1];
			pp[2]=p[2]+delta[2];
			real f_lu = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_lu = 0.0; }
			
			pp[0]=p[0];
			pp[1]=p[1];
			pp[2]=p[2]+delta[2];
			real f_u = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_u = 0.0; }
			
			pp[0]=p[0];
			pp[1]=p[1];
			pp[2]=p[2]-delta[2];
			real f_d = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			status = cell_location_2(sdp, ns, 0, center, pp); 
			if(status == 0){f_d = 0.0; }
			
			
			
			//at x direction
			//==================================================================
			//plane i-j
			gradFij1 = (f_fld-f_bld + 2.0*(f_fd-f_bd) + f_frd-f_brd)/(delta[0]);
			
			gradFij2 = (f_fl-f_bl + 2.0*(f_f-f_b) + f_fr-f_br)/(delta[0]);
			
			gradFij3 = (f_flu-f_blu + 2.0*(f_fu-f_bu) + f_fru-f_bru)/(delta[0]);
			
			gradFij = (gradFij1 + gradFij2 + gradFij3)/24.0;
			
			//plane i-k
			gradFik1 = (f_fld-f_bld + 2.0*(f_fl-f_bl) + f_flu-f_blu)/(delta[0]);
			
			gradFik2 = (f_fd-f_bd + 2.0*(f_f-f_b) + f_fd-f_bd)/(delta[0]);
			
			gradFik3 = (f_frd-f_brd + 2.0*(f_fr-f_br) + f_fru-f_bru)/(delta[0]);
			
			gradFik = (gradFik1 + gradFik2 + gradFik3)/24.0;
			
			//----------------------------------------------
					gradFi = (gradFij + gradFik)/2.0;
			//----------------------------------------------
			
			
			//at y direction
			//==================================================================
			//plane j-i
			gradFji1 = (f_frd-f_fld + 2.0*(f_rd-f_ld) + f_brd-f_bld)/delta[1];
			
			gradFji2 = (f_fr-f_fl + 2.0*(f_r-f_l) + f_br-f_bl)/delta[1];
			
			gradFji3 = (f_fru-f_flu + 2.0*(f_ru-f_lu) + f_bru-f_blu)/delta[1];
			
			gradFji = (gradFji1 + gradFji2 + gradFji3)/24.0;
			
			//plane j-k
			gradFjk1 = (f_brd - f_bld + 2.0*(f_br - f_bl) + f_bru - f_blu)/delta[1];
			
			gradFjk2 = (f_rd - f_ld + 2.0*(f_r - f_l) + f_ru - f_lu)/delta[1];
			
			gradFjk3 = (f_frd - f_fld + 2.0*(f_fr - f_fl) + f_fru - f_flu)/delta[1];
			
			gradFjk = (gradFjk1 + gradFjk2 + gradFjk3)/24.0;
			
			//----------------------------------------------
			 gradFj = (gradFji + gradFjk)/2.0;
			//----------------------------------------------
			
			//at z direction
			//==================================================================
			//plane k-i
			gradFki1 = (f_blu-f_bld + 2.0*(f_lu - f_ld) + f_flu-f_fld)/delta[2];
			
			gradFki2 = (f_lu-f_ld + 2.0*(f_u - f_d) + f_fu-f_fd)/delta[2];
			
			gradFki3 = (f_bru-f_brd + 2.0*(f_ru - f_rd) + f_fru-f_frd)/delta[2];
			
			gradFki = (gradFki1 + gradFki2 + gradFki3)/24.0;
			
			//plane k-j
			gradFkj1 = (f_blu-f_bld + 2.0*(f_lu-f_ld) + f_bru-f_brd)/delta[2];
			
			gradFkj2 = (f_lu-f_ld + 2.0*(f_u-f_d) + f_ru-f_rd)/delta[2];
			
			gradFkj3 = (f_flu-f_fld + 2.0*(f_fu-f_fd) + f_fru-f_frd)/delta[2];
			
			gradFkj = (gradFkj1 + gradFkj2 + gradFkj3)/24.0;
			
			//----------------------------------------------
			 gradFk = (gradFki + gradFkj)/2.0;
			//----------------------------------------------
			
			
			//Normal
			//=================================================================
			real norm_grad = sqrt(pow(gradFi,2)+pow(gradFj,2)+pow(gradFk,2));
			Normal[0] = gradFi/(norm_grad+1.0e-14);
			Normal[1] = gradFj/(norm_grad+1.0e-14);
			Normal[2] = gradFk/(norm_grad+1.0e-14);
			
			//printf("%lf %lf %lf %lf\n",p[0],p[1],p[2], f_c);
			
			//printf("%lf %lf %lf\n", Normal[0], Normal[1], Normal[2]);
			
			//////========================================================
			////// Calculating exact normalo vector=======================
			//////========================================================
			//p[0] = center[0];
			//p[1] = center[1];
			//p[2] = center[2];
			//if(Normal[0]>Normal[1] && Normal[0]>Normal[2]){
				//calculate_exact_normal_x_dominant(ns, clid,  p, 1);
			//}else if(Normal[1]>Normal[0] && Normal[1]>Normal[2]){
				//calculate_exact_normal_y_dominant(ns, clid,  p, 1);
			//}else{
				//calculate_exact_normal_z_dominant(ns, clid,  p, 1);
			//}
			//////========================================================
			
			
			for (int i=0; i<DIM; i++){
				dp_set_value(ns->ed.mult.dpnormal[i], clid, Normal[i]);
			}
		}
		// Destroy the iterator
		higcit_destroy(it);
		// Sync the distributed pressure property
		for (int i = 0; i < DIM; i++) {
			dp_sync(ns->ed.mult.dpnormal[i]);
		}
	}
}
#endif
