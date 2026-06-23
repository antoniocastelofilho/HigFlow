#if DIM == 3
#include "hig-flow-vof-finite-difference-normal-curvature_3D.h"

int cell_location(sim_domain *sd, higflow_solver *ns, int dim, Point Center, Point P){
	hig_cell *c = sd_get_cell_with_point(sd,P);
	if (c == NULL){
		return 0;
	} else {
		return 1;
	}
}

void mehta(higflow_solver *ns){
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
	
		Point p;
		p[0] = center[0];
		p[1] = center[1];
		p[2] = center[2];
		real frac = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		
	
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
		int status;
		
		real fracvol, gradFij1, gradFij2, gradFik1, gradFik2, gradFjk1, gradFjk2, gradFi, gradFj, gradFk;
		
		pp[0]=p[0];
		pp[1]=p[1];
		pp[2]=p[2];
		real f_c = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		
		pp[0]=p[0]+delta[0];
		pp[1]=p[1]-delta[1];
		pp[2]=p[2]-delta[2];
		real f_fld = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_fld = 0.0; }
		
		pp[0]=p[0]-delta[0];
		pp[1]=p[1]-delta[1];
		pp[2]=p[2]-delta[2];
		real f_bld = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_bld = 0.0; }
		
		pp[0]=p[0]+delta[0];
		pp[1]=p[1];
		pp[2]=p[2]-delta[2];
		real f_fd = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_fd = 0.0; }
		
		pp[0]=p[0]-delta[0];
		pp[1]=p[1];
		pp[2]=p[2]-delta[2];
		real f_bd = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_bd = 0.0; }
		
		pp[0]=p[0]+delta[0];
		pp[1]=p[1]+delta[1];
		pp[2]=p[2]-delta[2];
		real f_frd = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_frd = 0.0; }
		
		pp[0]=p[0]-delta[0];
		pp[1]=p[1]+delta[1];
		pp[2]=p[2]-delta[2];
		real f_brd = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_brd = 0.0; }
		
		pp[0]=p[0]+delta[0];
		pp[1]=p[1]-delta[1];
		pp[2]=p[2];
		real f_fl = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_fl = 0.0; }
		
		pp[0]=p[0]-delta[0];
		pp[1]=p[1]-delta[1];
		pp[2]=p[2];
		real f_bl = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_bl = 0.0; }
		
		pp[0]=p[0]+delta[0];
		pp[1]=p[1];
		pp[2]=p[2];
		real f_f = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_f = 0.0; }
		
		pp[0]=p[0]-delta[0];
		pp[1]=p[1];
		pp[2]=p[2];
		real f_b = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_b = 0.0; }
		
		pp[0]=p[0]+delta[0];
		pp[1]=p[1]+delta[1];
		pp[2]=p[2];
		real f_fr = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_fr = 0.0; }
		
		pp[0]=p[0]-delta[0];
		pp[1]=p[1]+delta[1];
		pp[2]=p[2];
		real f_br = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_br = 0.0; }
		
		pp[0]=p[0]+delta[0];
		pp[1]=p[1]-delta[1];
		pp[2]=p[2]+delta[2];
		real f_flu = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_flu = 0.0; }
		
		pp[0]=p[0]-delta[0];
		pp[1]=p[1]-delta[1];
		pp[2]=p[2]+delta[2];
		real f_blu = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_blu = 0.0; }
		
		pp[0]=p[0]+delta[0];
		pp[1]=p[1];
		pp[2]=p[2]+delta[2];
		real f_fu = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_fu = 0.0; }
		
		pp[0]=p[0]-delta[0];
		pp[1]=p[1];
		pp[2]=p[2]+delta[2];
		real f_bu = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_bu = 0.0; }
		
		pp[0]=p[0]+delta[0];
		pp[1]=p[1]+delta[1];
		pp[2]=p[2]+delta[2];
		real f_fru = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_fru = 0.0; }
		
		pp[0]=p[0]-delta[0];
		pp[1]=p[1]+delta[1];
		pp[2]=p[2]+delta[2];
		real f_bru = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_bru = 0.0; }
		
		pp[0]=p[0];
		pp[1]=p[1]+delta[1];
		pp[2]=p[2]-delta[2];
		real f_rd = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_rd = 0.0; }
		
		pp[0]=p[0];
		pp[1]=p[1]-delta[1];
		pp[2]=p[2]-delta[2];
		real f_ld = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_ld = 0.0; }
		
		pp[0]=p[0];
		pp[1]=p[1]+delta[1];
		pp[2]=p[2];
		real f_r = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_r = 0.0; }
		
		pp[0]=p[0];
		pp[1]=p[1]-delta[1];
		pp[2]=p[2];
		real f_l = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_l = 0.0; }
		
		pp[0]=p[0];
		pp[1]=p[1]+delta[1];
		pp[2]=p[2]+delta[2];
		real f_ru = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_ru = 0.0; }
		
		pp[0]=p[0];
		pp[1]=p[1]-delta[1];
		pp[2]=p[2]+delta[2];
		real f_lu = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_lu = 0.0; }
		
		pp[0]=p[0];
		pp[1]=p[1];
		pp[2]=p[2]+delta[2];
		real f_u = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_u = 0.0; }
		
		pp[0]=p[0];
		pp[1]=p[1];
		pp[2]=p[2]-delta[2];
		real f_d = compute_value_at_point(sdp, center, pp, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
		status = cell_location(sdp, ns, 0, center, pp); 
		if(status == 0){f_d = 0.0; }
		
		
		
		//at x direction
		//==================================================================
		//plane j-k
		gradFjk1 = (f_fru+2.0*f_fr+f_frd + 2.0*f_fu+4.0*f_f+2.0*f_fd + f_flu+2.0*f_fl+f_fld);
		gradFjk2 = (f_bru+2.0*f_br+f_brd + 2.0*f_bu+4.0*f_b+2.0*f_bd + f_blu+2.0*f_bl+f_bld);
		//----------------------------------------------
				gradFi = (gradFjk1 - gradFjk2)/(32.0*delta[0]);
		//----------------------------------------------
		
		//at y direction
		//==================================================================
		//plane i-k
		gradFik1 = (f_fru+2.0*f_fr+f_frd + 2.0*f_ru+4.0*f_r+2.0*f_rd + f_bru+2.0*f_br+f_brd);
		gradFik2 = (f_flu+2.0*f_fl+f_fld +2.0*f_lu+4.0*f_l+2.0*f_ld +  f_blu+2.0*f_bl+f_bld);
		//----------------------------------------------
				gradFj = (gradFik1 - gradFik2)/(32.0*delta[1]);
		//----------------------------------------------
		
		//at z direction
		//==================================================================
		//plane i-j
		gradFij1 = (f_fru+2.0*f_fu+f_flu + 2.0*f_ru+4.0*f_u+2.0*f_lu + f_bru+2.0*f_bu+f_blu);
		gradFij2 = (f_frd+2.0*f_fd+f_fld + 2.0*f_r+4.0*f_d+2.0*f_ld +  f_brd+2.0*f_bd+f_bld);
		//----------------------------------------------
				gradFk = (gradFij1 - gradFij2)/(32.0*delta[2]);
		//----------------------------------------------
		
		
		//Normal
		//=================================================================
		real norm_grad = sqrt(pow(gradFi,2)+pow(gradFj,2)+pow(gradFk,2));
		Normal[0] = gradFi/(norm_grad+1.0e-14);
		Normal[1] = gradFj/(norm_grad+1.0e-14);
		Normal[2] = gradFk/(norm_grad+1.0e-14);
		
		
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
#endif
