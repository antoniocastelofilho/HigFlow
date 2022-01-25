#include "hig-flow-vof-finite-difference-normal-curvature.h"

real norm_vec(real *vec, int n){
	real norm=0.0;
	for(int i=0;i<n;i++){
		norm += vec[i]*vec[i];
	}
	return norm = sqrt(norm);
}

void calculate_normal_cell_central_2nd_order_finite_difference_Horizontal(higflow_solver *ns, int clid, real Hk_,real H,real Hk,real dx,real dy,int auxh){
	real Hy=-auxh*(-Hk_+Hk)*dx/(2*dy);
	Point Normal;
	Normal[0]=auxh;
	Normal[1]=-auxh*Hy;
	real norm=norm_vec(Normal,2);
	Normal[0]=Normal[0]/norm;
	Normal[1]=Normal[1]/norm;
	for (int i=0; i<DIM; i++){
		dp_set_value(ns->ed.mult.dpnormal[i], clid, Normal[i]);
	}
}

void ELVIRA_calculate_normal_cell_progressive_1st_order_finite_difference_Horizontal(higflow_solver *ns, int clid, Point Normal, real Ht, real Hm, real dx, real dy, int auxh){
	real Hy=-auxh*(Ht-Hm)*dx/(dy);
	Normal[0]=auxh;
	Normal[1]=-auxh*Hy;
	real norm=norm_vec(Normal,2);
	Normal[0]=Normal[0]/norm;
	Normal[1]=Normal[1]/norm;
}

void ELVIRA_calculate_normal_cell_regressive_1st_order_finite_difference_Horizontal(higflow_solver *ns, int clid, Point Normal, real Hm,real Hb,real dx,real dy,int auxh){
	real Hy=-auxh*(Hm-Hb)*dx/(dy);
	Normal[0]=auxh;
	Normal[1]=-auxh*Hy;
	real norm=norm_vec(Normal,2);
	Normal[0]=Normal[0]/norm;
	Normal[1]=Normal[1]/norm;
}

void ELVIRA_calculate_normal_cell_central_2nd_order_finite_difference_Horizontal(higflow_solver *ns, int clid, Point Normal, real Ht,real Hb,real dx,real dy,int auxh){
	real Hy=-auxh*(Ht-Hb)*dx/(2.0*dy);
	Normal[0]=auxh;
	Normal[1]=-auxh*Hy;
	real norm=norm_vec(Normal,2);
	Normal[0]=Normal[0]/norm;
	Normal[1]=Normal[1]/norm;
}

void ELVIRA_calculate_normal_cell_central_2nd_order_finite_difference_Vertical(higflow_solver *ns, int clid, Point Normal, real Vr,real Vl,real dx,real dy,int auxv){
	real Vx=-auxv*(Vr-Vl)*dy/(2*dx);
	Normal[0]=-auxv*Vx;
	Normal[1]=auxv;
	real norm=norm_vec(Normal,2);
	Normal[0]=Normal[0]/norm;
	Normal[1]=Normal[1]/norm;
}

void ELVIRA_calculate_normal_cell_progressive_1st_order_finite_difference_Vertical(higflow_solver *ns, int clid, Point Normal, real Vr,real Vm,real dx,real dy,int auxv){
	real Vx=-auxv*(Vr-Vm)*dy/(dx);
	Normal[0]=-auxv*Vx;
	Normal[1]=auxv;
	real norm=norm_vec(Normal,2);
	Normal[0]=Normal[0]/norm;
	Normal[1]=Normal[1]/norm;
}

void ELVIRA_calculate_normal_cell_regressive_1st_order_finite_difference_Vertical(higflow_solver *ns, int clid, Point Normal,real Vm,real Vl,real dx,real dy,int auxv){
	real Vx=-auxv*(Vm-Vl)*dy/(dx);
	Normal[0]=-auxv*Vx;
	Normal[1]=auxv;
	real norm=norm_vec(Normal,2);
	Normal[0]=Normal[0]/norm;
	Normal[1]=Normal[1]/norm;
}

void calculate_curvature_cell_central_2nd_order_finite_difference_Horizontal(higflow_solver *ns, int clid, real Hk_,real H,real Hk,real dx,real dy,int auxh){
	real Hy=-auxh*(-Hk_+Hk)*dx/(2*dy);
	real Hyy=-auxh*(Hk_-2*H+Hk)*dx/(dy*dy);
	real curv = auxh*Hyy/pow((1+Hy*Hy),1.5);
	dp_set_value(ns->ed.mult.dpcurvature, clid, curv);// Set the curvature in the distributed curvature property
}

void calculate_normal_cell_central_2nd_order_finite_difference_Vertical(higflow_solver *ns, int clid, real Vk_,real V,real Vk,real dx,real dy,int auxv){
	real Vx=-auxv*(-Vk_+Vk)*dy/(2*dx);
	Point Normal;
	Normal[0]=-auxv*Vx;
	Normal[1]=auxv;
	real norm=norm_vec(Normal,2);
	Normal[0]=Normal[0]/norm;
	Normal[1]=Normal[1]/norm;
	for (int i=0; i<DIM; i++){
		dp_set_value(ns->ed.mult.dpnormal[i], clid, Normal[i]);
	}
}

void calculate_curvature_cell_central_2nd_order_finite_difference_Vertical(higflow_solver *ns, int clid, real Vk_,real V,real Vk,real dx,real dy,int auxv){

	real Vx=-auxv*(-Vk_+Vk)*dy/(2*dx);
	real Vxx=-auxv*(Vk_-2*V+Vk)*dy/(dx*dx);
	real curv = auxv*Vxx/pow((1+Vx*Vx),1.5);

	dp_set_value(ns->ed.mult.dpcurvature, clid, curv);// Set the curvature in the distributed curvature property
}

void calculate_normal_cell_progressive_2nd_order_finite_difference_Horizontal(higflow_solver *ns, int clid, real H,real Hk,real Hkk,real dx,real dy,int auxh){
	real Hy=-auxh*(-1.5*H+2*Hk-0.5*Hkk)*dx/(dy);
	Point Normal;
	Normal[0]=auxh;
	Normal[1]=-auxh*Hy;
	real norm=norm_vec(Normal,2);
	Normal[0]=Normal[0]/norm;
	Normal[1]=Normal[1]/norm;

	for (int i=0; i<DIM; i++){
		dp_set_value(ns->ed.mult.dpnormal[i], clid, Normal[i]);
	}
}

void calculate_curvature_cell_progressive_1st_order_finite_difference_Horizontal(higflow_solver *ns, int clid, real H,real Hk,real Hkk,real dx,real dy,int auxh){
	real Hy=-auxh*(-1.5*H+2*Hk-0.5*Hkk)*dx/(dy);
	real Hyy=-auxh*(H-2*Hk+Hkk)*dx/(dy*dy);
	
	real curv = auxh*Hyy/pow((1+Hy*Hy),1.5);
	
	dp_set_value(ns->ed.mult.dpcurvature, clid, curv);// Set the curvature in the distributed curvature property
}

void calculate_normal_cell_progressive_2nd_order_finite_difference_Vertical(higflow_solver *ns, int clid, real V,real Vk,real Vkk,real dx,real dy,int auxv){
	real Vx=-auxv*(-1.5*V+2*Vk-0.5*Vkk)*dy/(dx);
	Point Normal;
	Normal[0]=-auxv*Vx;
	Normal[1]=auxv;
	real norm=norm_vec(Normal,2);
	Normal[0]=Normal[0]/norm;
	Normal[1]=Normal[1]/norm;

	for (int i=0; i<DIM; i++){
		dp_set_value(ns->ed.mult.dpnormal[i], clid, Normal[i]);
	}
}

void calculate_curvature_cell_progressive_1st_order_finite_difference_Vertical(higflow_solver *ns, int clid, real V,real Vk,real Vkk,real dx,real dy,int auxv){
	real Vx=-auxv*(-1.5*V+2*Vk-0.5*Vkk)*dy/(dx);
	real Vxx=-auxv*(V-2*Vk+Vkk)*dy/(dx*dx);
	
	real curv = auxv*Vxx/pow((1+Vx*Vx),1.5);
	
	dp_set_value(ns->ed.mult.dpcurvature, clid, curv);// Set the curvature in the distributed curvature property
}

void calculate_normal_cell_regressive_2nd_order_finite_difference_Horizontal(higflow_solver *ns, int clid, real Hk_k_,real Hk_,real H,real dx,real dy,int auxh){
	real Hy=-auxh*(0.5*Hk_k_-2*Hk_+1.5*H)*dx/(dy);
	Point Normal;
	Normal[0]=auxh;
	Normal[1]=-auxh*Hy;
	real norm=norm_vec(Normal,2);
	Normal[0]=Normal[0]/norm;
	Normal[1]=Normal[1]/norm;

	for (int i=0; i<DIM; i++){
		dp_set_value(ns->ed.mult.dpnormal[i], clid, Normal[i]);
	}
}

void calculate_curvature_cell_regressive_1st_order_finite_difference_Horizontal(higflow_solver *ns, int clid, real Hk_k_,real Hk_,real H,real dx,real dy,int auxh){
	real Hy=-auxh*(0.5*Hk_k_-2*Hk_+1.5*H)*dx/(dy);
	real Hyy=-auxh*(Hk_k_- 2*Hk_ + H)*dx/(dy*dy);
	
	real curv = auxh*Hyy/pow((1+Hy*Hy),1.5);
	
	dp_set_value(ns->ed.mult.dpcurvature, clid, curv);// Set the curvature in the distributed curvature property
}

void calculate_normal_cell_regressive_2nd_order_finite_difference_Vertical(higflow_solver *ns, int clid, real Vk_k_,real Vk_,real V,real dx,real dy,int auxv){
	real Vx=-auxv*(0.5*Vk_k_-2*Vk_+1.5*V)*dy/(dx);
	Point Normal;
	Normal[0]=-auxv*Vx;
	Normal[1]=auxv;
	real norm=norm_vec(Normal,2);
	Normal[0]=Normal[0]/norm;
	Normal[1]=Normal[1]/norm;

	for (int i=0; i<DIM; i++){
		dp_set_value(ns->ed.mult.dpnormal[i], clid, Normal[i]);
	}
}

void calculate_curvature_cell_regressive_1st_order_finite_difference_Vertical(higflow_solver *ns, int clid, real Vk_k_,real Vk_,real V,real dx,real dy,int auxv){
	real Vx=-auxv*(0.5*Vk_k_-2*Vk_+1.5*V)*dy/(dx);
	real Vxx=-auxv*(Vk_k_- 2*Vk_ + V)*dy/(dx*dx);
	real curv = auxv*Vxx/pow((1.0+Vx*Vx),1.5);
	dp_set_value(ns->ed.mult.dpcurvature, clid, curv);// Set the curvature in the distributed curvature property
}

void calculate_exact_normal_x_dominant_2D(higflow_solver *ns, int clid, Point p){
   Point Nt;
   if(p[0] > 0.5){
      Nt[0] = sqrt(pow(0.25,2)-(pow(p[1]-0.5,2)));
      Nt[1] = p[1]-0.5;
      real norm=norm_vec(Nt,2);
      Nt[0] = -Nt[0]/norm;
      Nt[1] = -Nt[1]/norm;
   }else{
      Nt[0] = (-sqrt(pow(0.25,2)-(pow(p[1]-0.5,2))));
      Nt[1] = p[1]-0.5;
      real norm=norm_vec(Nt,2);
      Nt[0] = -Nt[0]/norm;
      Nt[1] = -Nt[1]/norm;
   }	
}

void calculate_exact_normal_y_dominant_2D(higflow_solver *ns, int clid, Point p){
		 Point Nt;
		 if(p[1] > 0.5){
			 Nt[0] = p[0]-0.5;
			 Nt[1] = sqrt(pow(0.25,2)-(pow(p[0]-0.5,2)));
			 real norm=norm_vec(Nt,2);
			 Nt[0] = -Nt[0]/norm;
			 Nt[1] = -Nt[1]/norm;
		}else{
			 Nt[0] = p[0]-0.5;
			 Nt[1] = -sqrt(pow(0.25,2)-(pow(p[0]-0.5,2)));
			 real norm=norm_vec(Nt,3);
			 Nt[0] = -Nt[0]/norm;
			 Nt[1] = -Nt[1]/norm;
		}	
}
