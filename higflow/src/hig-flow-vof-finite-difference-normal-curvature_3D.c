#if DIM == 3
#include "hig-flow-vof-finite-difference-normal-curvature.h"


void calculate_HF_curvature_x_dominant(higflow_solver *ns, int clid, real Xl,
	 real Xr, real Yb, real Yu, real XY, real Xul, real Xur, real Xbl, 
	 real Xbr, real dx, real dy, real dz, int aux){
		real hy = -aux*(Xr-Xl)*dx/(2*dy);
		real hz = -aux*(Yu-Yb)*dx/(2*dz);
		real hyy = -aux*(Xr-2*XY+Xl)*dx/(dy*dy);
		real hzz = -aux*(Yu-2*XY+Yb)*dx/(dz*dz);
		real hyz = -aux*(Xur-Xul-Xbr+Xbl)*dx/(2*dy*2*dz);
		real curv = aux*(hyy+hzz+hyy*pow(hz,2)+hzz*pow(hy,2)-2*hyz*hy*hz)/
					pow((1+pow(hy,2)+pow(hz,2)), 1.5);
		
		dp_set_value(ns->ed.mult.dpcurvature, clid, curv);// Set the curvature in the distributed curvature property
		
		//if(fabs(curv-4.0)>0.5){
			//printf("%lf\n", curv);
		//}
}

void calculate_HF_curvature_y_dominant(higflow_solver *ns, int clid, real Xl,
	 real Xr, real Yb, real Yu, real XY, real Xul, real Xur, real Xbl, 
	 real Xbr, real dx, real dy, real dz, int aux){
		real hx = -aux*(Xr-Xl)*dy/(2*dx);
		real hz = -aux*(Yu-Yb)*dy/(2*dz);
		real hxx = -aux*(Xr-2*XY+Xl)*dy/(dx*dx);
		real hzz = -aux*(Yu-2*XY+Yb)*dy/(dz*dz);
		real hxz = -aux*(Xur-Xul-Xbr+Xbl)*dy/(2*dy*2*dz);
		real curv = aux*(hxx+hzz+hxx*pow(hz,2)+hzz*pow(hx,2)-2*hxz*hx*hz)/
					pow((1+pow(hx,2)+pow(hz,2)), 1.5);
		//printf("%lf\n", curv);
		dp_set_value(ns->ed.mult.dpcurvature, clid, curv);// Set the curvature in the distributed curvature property
		
		//if(fabs(curv-4.0)>0.5){
			//printf("%lf\n", curv);
		//}
}

void calculate_HF_curvature_z_dominant(higflow_solver *ns, int clid, real Xl,
	 real Xr, real Yb, real Yu, real XY, real Xul, real Xur, real Xbl, 
	 real Xbr, real dx, real dy, real dz, int aux){
		real hx = -aux*(Xr-Xl)*dz/(2*dx);
		real hy = -aux*(Yu-Yb)*dz/(2*dy);
		real hxx = -aux*(Xr-2*XY+Xl)*dz/(dx*dx);
		real hyy = -aux*(Yu-2*XY+Yb)*dz/(dy*dy);
		real hxy = -aux*(Xur-Xul-Xbr+Xbl)*dz/(2*dy*2*dx);
		real curv = aux*(hxx+hyy+hxx*pow(hy,2)+hyy*pow(hx,2)-2*hxy*hx*hy)/
					pow((1+pow(hx,2)+pow(hy,2)), 1.5);
		//printf("%lf\n", curv);
		dp_set_value(ns->ed.mult.dpcurvature, clid, curv);// Set the curvature in the distributed curvature property
		
		//if(fabs(curv-4.0)>0.5){
			//printf("%lf\n", curv);
		//}
}

void calculate_HF_normal_x_dominant(higflow_solver *ns, int clid, real Xl,
	 real Xr, real Yb, real Yu,real dx,real dy, real dz,int auxh){
		real hy = -auxh*(Xr-Xl)*dx/(2*dy);
		real hz = -auxh*(Yu-Yb)*dx/(2*dz);
		Point Normal;
		Normal[0]=auxh;
		Normal[1]=-auxh*hy;
		Normal[2]=-auxh*hz;
		real norm=norm_vec(Normal,3)+1e-14;
		Normal[0]=Normal[0]/norm;
		Normal[1]=Normal[1]/norm;
		Normal[2]=Normal[2]/norm;
		
		//printf("%lf %lf %lf\n", Normal[0], Normal[1], Normal[2]);
		
		for (int i=0; i<DIM; i++){
			dp_set_value(ns->ed.mult.dpnormal[i], clid, Normal[i]);
		}
}

void calculate_HF_normal_y_dominant(higflow_solver *ns, int clid, real Xl,
	 real Xr, real Yb, real Yu,real dx,real dy, real dz,int auxh){
		real hx = -auxh*(Xr-Xl)*dy/(2*dx);
		real hz = -auxh*(Yu-Yb)*dy/(2*dz);
		Point Normal;
		Normal[0]=-auxh*hx;
		Normal[1]=auxh;
		Normal[2]=-auxh*hz;
		real norm=norm_vec(Normal,3)+1e-14;
		Normal[0]=Normal[0]/norm;
		Normal[1]=Normal[1]/norm;
		Normal[2]=Normal[2]/norm;
		
		//printf("%lf %lf %lf\n", Normal[0], Normal[1], Normal[2]);
		
		for (int i=0; i<DIM; i++){
			dp_set_value(ns->ed.mult.dpnormal[i], clid, Normal[i]);
		}
}

void calculate_HF_normal_z_dominant(higflow_solver *ns, int clid, real Xl,
	 real Xr, real Yb, real Yu,real dx,real dy, real dz,int auxh){
		real hy = -auxh*(Xr-Xl)*dz/(2*dy);
		real hx = -auxh*(Yu-Yb)*dz/(2*dx);
		Point Normal;
		Normal[0]=-auxh*hx;
		Normal[1]=-auxh*hy;
		Normal[2]=auxh;
		real norm=norm_vec(Normal,3)+1e-14;
		Normal[0]=Normal[0]/norm;
		Normal[1]=Normal[1]/norm;
		Normal[2]=Normal[2]/norm;
		
		//printf("%lf %lf %lf\n", Normal[0], Normal[1], Normal[2]);
		
		for (int i=0; i<DIM; i++){
			dp_set_value(ns->ed.mult.dpnormal[i], clid, Normal[i]);
		}
}


void calculate_exact_normal_x_dominant(higflow_solver *ns, int clid, Point p, int auxh){
		 
		 Point Nt;
		 if(p[0] > 0.5){
			 Nt[0] = sqrt(pow(0.25,2)-(pow(p[1]-0.5,2)+pow(p[2]-0.5,2)));
			 //printf(">=== %lf %lf %lf\n", pow(0.25,2),(pow(p[0]-0.5,2)+pow(p[1]-0.5,2)), Nt[0]);
			 Nt[1] = p[1]-0.5;
			 Nt[2] = p[2]-0.5;
			 real norm=norm_vec(Nt,3)+1e-14;;
			 Nt[0] = -Nt[0]/norm;
			 Nt[1] = -Nt[1]/norm;
			 Nt[2] = -Nt[2]/norm;
		}else{
			 Nt[0] = (-sqrt(pow(0.25,2)-(pow(p[1]-0.5,2)+pow(p[2]-0.5,2))));
			 //printf("%lf %lf %lf\n", pow(0.25,2) , (pow(p[0]-0.5,2)+pow(p[1]-0.5,2)), Nt[0]);
			 Nt[1] = p[1]-0.5;
			 Nt[2] = p[2]-0.5;
			 real norm=norm_vec(Nt,3)+1e-14;;
			 Nt[0] = -Nt[0]/norm;
			 Nt[1] = -Nt[1]/norm;
			 Nt[2] = -Nt[2]/norm;
		}	
			printf("%lf %lf %lf\n", Nt[0], Nt[1], Nt[2]);
}

void calculate_exact_normal_y_dominant(higflow_solver *ns, int clid, Point p, int auxh){
		 
		 Point Nt;
		 if(p[1] > 0.5){
			 Nt[0] = p[0]-0.5;
			 Nt[1] = sqrt(pow(0.25,2)-(pow(p[0]-0.5,2)+pow(p[2]-0.5,2)));
			 //printf("%lf %lf %lf\n", pow(0.25,2),(pow(p[0]-0.5,2)+pow(p[1]-0.5,2)), Nt[1]);
			 Nt[2] = p[2]-0.5;
			 real norm=norm_vec(Nt,3)+1e-14;;
			 Nt[0] = -Nt[0]/norm;
			 Nt[1] = -Nt[1]/norm;
			 Nt[2] = -Nt[2]/norm;
		}else{
			 Nt[0] = p[0]-0.5;
			 Nt[1] = -sqrt(pow(0.25,2)-(pow(p[0]-0.5,2)+pow(p[2]-0.5,2)));
			 //printf("%lf %lf %lf\n", pow(0.25,2) , (pow(p[0]-0.5,2)+pow(p[1]-0.5,2)), Nt[1]);
			 Nt[2] = p[2]-0.5;
			 real norm=norm_vec(Nt,3)+1e-14;;
			 Nt[0] = -Nt[0]/norm;
			 Nt[1] = -Nt[1]/norm;
			 Nt[2] = -Nt[2]/norm;
		}	
		
		printf("%lf %lf %lf\n", Nt[0], Nt[1], Nt[2]);
}

void calculate_exact_normal_z_dominant(higflow_solver *ns, int clid, Point p, int auxh){
		 
		 Point Nt;
		 if(p[2] > 0.5){
			 Nt[0] = p[0]-0.5;
			 Nt[1] = p[1]-0.5;
			 Nt[2] = sqrt(pow(0.25,2)-(pow(p[0]-0.5,2)+pow(p[1]-0.5,2)));
			 //printf("%lf %lf %lf\n", pow(0.25,2) , (pow(p[0]-0.5,2)+pow(p[1]-0.5,2)), Nt[2]);
			 real norm=norm_vec(Nt,3)+1e-14;;
			 Nt[0] = -Nt[0]/norm;
			 Nt[1] = -Nt[1]/norm;
			 Nt[2] = -Nt[2]/norm;
		}else{
			 Nt[0] = p[0]-0.5;
			 Nt[1] = p[1]-0.5;
			 Nt[2] = -sqrt(pow(0.25,2)-(pow(p[0]-0.5,2)+pow(p[1]-0.5,2)));
			 //printf("%lf %lf %lf\n", pow(0.25,2) , (pow(p[0]-0.5,2)+pow(p[1]-0.5,2)), Nt[2]);
			 real norm=norm_vec(Nt,3)+1e-14;;
			 Nt[0] = -Nt[0]/norm;
			 Nt[1] = -Nt[1]/norm;
			 Nt[2] = -Nt[2]/norm;
		}	
		
		printf("%lf %lf %lf\n", Nt[0], Nt[1], Nt[2]);
}
		 
//Progressive HF method===================================================
void calculate_progressive_HF_curvature_x_dominant(higflow_solver *ns, int clid, real Xl,
	 real Xr, real XY, real Xul, real Xur, real Xu, real Xuul, real Xuur, 
	 real Xuu, real dx, real dy, real dz, int aux){
		real hy = -aux*(Xr-Xl)*dx/(2*dy);
		real hz = -aux*(-1.5*XY+2*Xu-0.5*Xuu)*dx/(dz);
		real hyy = -aux*(Xr-2*XY+Xl)*dx/(dy*dy);
		real hzz = -aux*(XY-2*Xu+Xuu)*dx/(dz*dz);
		real hyz = -aux*(-1.5*(Xr-Xl)+2*(Xur-Xul)-0.5*(Xuur-Xuul))*dx/(2*dy*dz);
		real curv = aux*(hyy+hzz+hyy*pow(hz,2)+hzz*pow(hy,2)-2*hyz*hy*hz)/
					pow((1+pow(hy,2)+pow(hz,2)), 1.5);
		//printf("%lf\n", curv);
		dp_set_value(ns->ed.mult.dpcurvature, clid, curv);// Set the curvature in the distributed curvature property
}

//horizintal is dominant with progressivite in z
void calculate_progressive_HF_normal_x_dominant(higflow_solver *ns, int clid, real Xl,
	 real Xr, real XY, real Yu, real Yuu, real dx,real dy, real dz,int auxh){
		real hy = -auxh*(Xr-Xl)*dx/(2*dy);
		real hz = -auxh*(-1.5*XY+2*Yu-0.5*Yuu)*dx/(dz);
		Point Normal;
		Normal[0]=auxh;
		Normal[1]=-auxh*hy;
		Normal[2]=-auxh*hz;
		real norm=norm_vec(Normal,3)+1e-14;
		Normal[0]=Normal[0]/norm;
		Normal[1]=Normal[1]/norm;
		Normal[2]=Normal[2]/norm;
		
		//printf("%lf %lf %lf\t", Normal[0], Normal[1], Normal[2]);
		
		for (int i=0; i<DIM; i++){
			dp_set_value(ns->ed.mult.dpnormal[i], clid, Normal[i]);
		}
}

//Regressive HF method===================================================
void calculate_regressive_HF_curvature_x_dominant(higflow_solver *ns, int clid, real Xl,
	 real Xr, real XY, real Xbl, real Xbr, real Xb, real Xbbl, real Xbbr, 
	 real Xbb, real dx, real dy, real dz, int aux){
		real hy = -aux*(Xr-Xl)*dx/(2*dy);
		real hz = -aux*(1.5*XY-2*Xb+0.5*Xbb)*dx/(dz);
		real hyy = -aux*(Xr-2*XY+Xl)*dx/(dy*dy);
		real hzz = -aux*(XY-2*Xb+Xbb)*dx/(dz*dz);
		real hyz = -aux*(1.5*(Xr-Xl)-2*(Xbr-Xbl)+0.5*(Xbbr-Xbbl))*dx/(2*dy*dz);
		real curv = aux*(hyy+hzz+hyy*pow(hz,2)+hzz*pow(hy,2)-2*hyz*hy*hz)/
					pow((1+pow(hy,2)+pow(hz,2)), 1.5);
		//printf("%lf\n", curv);
		dp_set_value(ns->ed.mult.dpcurvature, clid, curv);// Set the curvature in the distributed curvature property
}

void calculate_regressive_HF_normal_x_dominant(higflow_solver *ns, int clid, real Xl,
	 real Xr, real XY, real Yb, real Ybb, real dx,real dy, real dz,int auxh){
		real hy = -auxh*(Xr-Xl)*dx/(2*dy);
		real hz = -auxh*(1.5*XY-2*Yb+0.5*Ybb)*dx/(dz);
		Point Normal;
		Normal[0]=auxh;
		Normal[1]=-auxh*hy;
		Normal[2]=-auxh*hz;
		real norm=norm_vec(Normal,3)+1e-14;
		Normal[0]=Normal[0]/norm;
		Normal[1]=Normal[1]/norm;
		Normal[2]=Normal[2]/norm;
		
		//printf("%lf %lf %lf\t", Normal[0], Normal[1], Normal[2]);
		
		for (int i=0; i<DIM; i++){
			dp_set_value(ns->ed.mult.dpnormal[i], clid, Normal[i]);
		}
}
#endif
