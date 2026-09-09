#if DIM == 3
#include "hig-flow-vof-plic-3D.h"
#include "hig-flow-vof-plic.h"
#include "hig-flow-vof-finite-difference-normal-curvature.h"

real newton_raphson(real x0, real A, real B, real C, real D){
	
	real erro, fx, dfx, x1;
	erro = 0.0;
	int i=0;
	do{
		fx = A*pow(x0,3)+B*pow(x0,2)+C*x0+D;
		dfx = 3*A*pow(x0,2)+2*B*x0+C;
		x1 = x0-(fx/dfx);
		erro = fabs(x0-x1);
		x0 = x1;
		i++;
	}while(erro > 1e-10 && i < 1000);
			
	return x0;
}
//======================================================================
real solver_equation_second_order(real volume, real n_x, real n_y, real n_z,
 real dx, real dy, real dz){
	
	real A, B, C, Value;

	A = 3.0*dx;
	B = -3.0*n_x*pow(dx,2);
	C = pow(n_x,2)*pow(dx,3)-6*n_y*n_z*volume;
		
	if(A==0){
		Value=-C/B;
		return Value;
	}
	Value=(-B+sqrt(B*B-4.0*A*C))/(2.0*A);
	
	return Value;
}
//======================================================================
real solver_equation_third_order(real volume, real n_x, real n_y, real n_z,
 real dx, real dy, real dz){
		 
	real A, B, C, D, erro, fx, dfx, x1, frac, Value;
	real ig1,ig2,ig3;
	Point p0;

	A = -1.0;
	B = 3.0*(n_x*dx+n_y*dy);
	C = -3.0*(pow(n_x,2)*pow(dx,2)+pow(n_y,2)*pow(dy,2));
	D = pow(n_x,3)*pow(dx,3)+pow(n_y,3)*pow(dy,3)-6.0*n_x*n_y*n_z*volume;
	
	frac = volume/(dx*dy*dz);
	ig1 = dx*sqrt(3)*frac; //diagonal do cubo
	Value = newton_raphson(ig1, A, B, C, D);
		   
	return Value;
}
//======================================================================
real solver_equation_third_order_2(real volume, real n_x, real n_y, real n_z,
 real dx, real dy, real dz){
		 
	real A, B, C, D, erro, fx, dfx, x1, frac, Value;
	real ig1,ig2,ig3;
	Point p0;

	A = -2.0;
	B = 3.0*(n_x*dx+n_y*dy+n_z*dz);
	C = -3.0*(pow(n_x,2)*pow(dx,2)+pow(n_y,2)*pow(dy,2)+pow(n_z,2)*pow(dz,2));
	D = pow(n_x,3)*pow(dx,3)+pow(n_y,3)*pow(dy,3)+pow(n_z,3)*pow(dz,3)-6.0*n_x*n_y*n_z*volume;
		   
	 
	frac = volume/(dx*dy*dz);
	ig1 = dx*sqrt(3)*frac; //diagonal do cubo
	Value = newton_raphson(ig1, A, B, C, D);
		   	
	return Value;
}
//======================================================================
real solver_equation_third_order_3(real volume, real n_x, real n_y, real n_z,
 real dx, real dy, real dz){
		 
	real Value = 0.5*(n_x*dx+n_y*dy)+(n_z*volume)/(dx*dy);
	
	return Value;
}
//======================================================================

//if only one argument of the normal vector is non-zero ================
real solver_equation_n_nonzero(real volume, real n1, real dy, real dz){

	real Value=n1*volume/(dy*dz);
	
	return Value;
}
//======================================================================

//calculating the volume when all arguments of the normal vector is non-zero
real volume_3D(real n_x, real n_y, real n_z, real dx, real dy, real dz, real d, real aux_d){
	
	real aux_x, aux_y, aux_z, aux_x_y, aux_x_z, aux_y_z;
	n_x = fabs(n_x);
	n_y = fabs(n_y);
	n_z = fabs(n_z);
	real n1 = n_x+1e-14;
	real n2 = n_y+1e-14;
	real n3 = n_z+1e-14;	
	d = fabs(d);
	real volume;
	volume = 0.0;
	real tol_n = 1e-8;
	
	aux_x = d - n1*dx;
	aux_y = d - n2*dy;
	aux_z = d - n3*dz;
	aux_x_y = d - n1*dx - n2*dy;
	aux_x_z = d - n1*dx - n3*dz;
	aux_y_z = d - n2*dy - n3*dz;
			
	volume = (1/(6.0*n1*n2*n3))*(pow(d,3) - H(aux_x)*pow(d-n1*dx,3) - H(aux_y)*pow(d-n2*dy,3) 
	- H(aux_z)*pow(d-n3*dz,3) + H(aux_x_y)*pow(d-n1*dx-n2*dy,3)
	+ H(aux_x_z)*pow(d-n1*dx-n3*dz,3) + H(aux_y_z)*pow(d-n2*dy-n3*dz,3)); 
	
	if(aux_d<0){
		volume=dx*dy*dz-volume;
	}
	
	if(n_x<=tol_n && (n_y>tol_n && n_z>tol_n) || n_y<=tol_n  && 
	(n_x>tol_n && n_z>tol_n) || n_z<=tol_n  && (n_y>tol_n && n_x>tol_n)){
		volume = (dz/(2.0*n1*n2))*(pow(d,2) - H(aux_x)*pow(d-n1*dx,2) 
		- H(aux_y)*pow(d-n2*dy,2));
		if(aux_d<0){
			volume = dx*dy*dz-volume;
		}
	}
	return volume;
}
//======================================================================
//checking if there is an extra triangle.
real H(real x){
	
	real value = (x <= 0.0) ? 0.0 : 1.0;
	
	return value;
}
//======================================================================
real trans_center_to_p0(Point Delta,Point Normal,real d_from_center){
	
		real d_from_p0 =  -Normal[0] * 0.5 * Delta[0] - Normal[1] * 0.5 * Delta[1]
						- Normal[2] * 0.5 * Delta[2] - d_from_center;
		
	return d_from_p0;
}
//======================================================================
real trans_p0_to_center(Point Delta, Point Normal, real d_from_p0) {
	
	real nx = Normal[0];
	real ny = Normal[1];
	real nz = Normal[2];
	real n1 = fabs(nx);
	real n2 = fabs(ny);
	real n3 = fabs(nz);
	real tol_n = 1e-8;
	real n;
	
	real d_from_center =  - Normal[0] * 0.5 * Delta[0] - Normal[1] * 0.5 * Delta[1]
						- Normal[2] * 0.5 * Delta[2] - d_from_p0;
	
	return d_from_center;
}
//======================================================================
// O vetor normal deve apontar para fora da interface.
real distance_from_center_3D(Point Normal,Point Delta,real VOLUME){
	
	real tol_volume = 1e-8, tol_n = 1e-8;
	real nx = Normal[0], ny = Normal[1], nz = Normal[2];
	real n_x = fabs(nx), n_y = fabs(ny), n_z = fabs(nz);
	real dx = Delta[0], dy = Delta[1], dz = Delta[2];
	Point Value;
	real volume, d;
		
	real n;
	n = maximo(n_x, n_y, n_z);
	n = (n==n_x) ? nx : ((n==n_y) ? ny : nz);
	//==================================================================
	//first case========================================================
	//==================================================================
	if ((n_x < tol_n && n_y < tol_n) || (n_x < tol_n && 
	n_z < tol_n) || (n_y < tol_n && n_z < tol_n)){
		
		real n1, n_1;	
		n1 = (n_x < tol_n && n_y < tol_n) ? nz : ((n_x < tol_n && n_z < tol_n) ? ny : nx);
		n_1 = (n_x < tol_n && n_y < tol_n) ? -n_z : ((n_x < tol_n && n_z < tol_n) ? -n_y : -n_x);
		
		d  = solver_equation_n_nonzero(VOLUME, n1, dy, dz);
		d = 0.5*Delta[0] - sign(n)*d;
		return d;
		
	}else if(n_x < tol_n || n_y < tol_n || n_z < tol_n){
	//==================================================================
	//Second case=======================================================
	//==================================================================
		real n1, n2, n3, n_1, n_2, n_3, f;
		n1 = (n_x < tol_n) ? ny : ((n_y < tol_n) ? nx : nx);
		n2 = (n_x < tol_n) ? nz : ((n_y < tol_n) ? nz : ny);
		n3 = 0.0;
		n_1 = (n_x < tol_n) ? -n_y : ((n_y < tol_n) ? -n_x : -n_x);
		n_2 = (n_x < tol_n) ? -n_z : ((n_y < tol_n) ? -n_z : -n_y);
		n_3 = 0.0;
		Normal[0] = n_1;
		Normal[1] = n_2;
		Normal[2] = n_3;
		n = n2;
		
		if(n2 < 0.0 && VOLUME < 0.5*(dx*dy*dz) || n2 > 0.0 && VOLUME > 0.5*(dx*dy*dz)){
		//Parte inferior================================================
		
			volume = (VOLUME < 0.5*(dx*dy*dz)) ? VOLUME : dx*dy*dz - VOLUME;
			
			f = volume/(dx*dy*dz); 
			d = sqrt((2*volume*n_1*n_2)/dx);
			d = -sign(n) * trans_p0_to_center(Delta, Normal, d);
			return d;
		}
		else if(n2 < 0.0 && VOLUME > 0.5*(dx*dy*dz) || n2 > 0.0 && VOLUME < 0.5*(dx*dy*dz)){
		//Parte superior================================================
				
			volume = (VOLUME < 0.5*(dx*dy*dz)) ? VOLUME : dx*dy*dz - VOLUME;
			
			f = volume/(dx*dy*dz);
			d = sqrt((2*volume*n_1*n_2)/dx);
			d = sign(n) * trans_p0_to_center(Delta, Normal, d);
			return d;
		}
	}else if(n_x>tol_n && n_y>tol_n && n_z>tol_n){
	//==================================================================
	//Third case =======================================================
	//==================================================================
		real n1, n2, n3, n_1, n_2, n_3, f;
		
		n_1 = -n_x;
		n_2 = -n_y;
		n_3 = -n_z;
		Normal[0] = n_1;
		Normal[1] = n_2;
		Normal[2] = n_3;
		n = nz;
		
		real pow_nx = pow(n_x,3)*pow(dx,3);
		real pow_ny = pow(n_y,3)*pow(dy,3);
		real pow_nz = pow(n_z,3)*pow(dz,3);
		//Parte inferior================================================
		//==============================================================
		if(nz < 0.0 && VOLUME < 0.5*(dx*dy*dz) || nz > 0.0 && VOLUME > 0.5*(dx*dy*dz)){
			
			volume = (nz < 0.0 && VOLUME < 0.5*(dx*dy*dz)) ? VOLUME : dx*dy*dz - VOLUME;
			
			f = volume/(dx*dy*dz);		
			real c = 6.0*n_x*n_y*n_z*volume;
			
			if((c < pow_nx) && (c < pow_ny) && (c < pow_nz)){
			////Este caso funciona quando a piramide esta completamente dentro da celula	
				d = cbrt(6*volume*n_x*n_y*n_z);
				d = -sign(n) * trans_p0_to_center(Delta, Normal, d);
				return d;
			}else if((c > pow_nx) && (c < pow_ny) && (c < pow_nz)){
			////Este caso funciona quando tem-se apenas uma piramide fora da celula
				d = solver_equation_second_order(volume, n_x, n_y, n_z, dx, dy, dz);
				d = -sign(n) * trans_p0_to_center(Delta, Normal, d);
				return d;
			}else if((c < pow_nx) && (c > pow_ny) && (c < pow_nz)){
			////Este caso funciona quando tem-se apenas uma piramide fora da celula
				d = solver_equation_second_order(volume, n_y, n_x, n_z, dy, dx, dz);
				d = -sign(n) * trans_p0_to_center(Delta, Normal, d);
				return d;
			}else if((c < pow_nx) && (c < pow_ny) && (c > pow_nz)){
			////Este caso funciona quando tem-se apenas uma piramide fora da celula
				d = solver_equation_second_order(volume, n_z, n_y, n_x, dz, dy, dx);
				d = -sign(n) * trans_p0_to_center(Delta, Normal, d);
				return d;
			}else if((c > pow_nx) && (c > pow_ny) && (c < pow_nz)){
			////Este caso funciona quando tem-se duas piramides fora da celula
				d = solver_equation_third_order(volume, n_x, n_y, n_z, dx, dy, dz);
				d = -sign(n) * trans_p0_to_center(Delta, Normal, d);
				real d2 = solver_equation_third_order_3(volume, n_x, n_y, n_z, dx, dy, dz);
				if(d2>n_x*dx && d2>n_y*dy && d2<n_z*dz){
					d = -sign(n) * trans_p0_to_center(Delta, Normal, d2);
				}
				return d;
			}else if((c > pow_nx) && (c < pow_ny) && (c > pow_nz)){
			////Este caso funciona quando tem-se duas piramides fora da celula
				d = solver_equation_third_order(volume, n_x, n_z, n_y, dx, dz, dy);
				d = -sign(n) * trans_p0_to_center(Delta, Normal, d);
				real d2 = solver_equation_third_order_3(volume, n_x, n_z, n_y, dx, dz, dy);
				if(d2>n_x*dx && d2<n_y*dy && d2>n_z*dz){
						d = -sign(n) * trans_p0_to_center(Delta, Normal, d2);
				}
				return d;
			}else if((c < pow_nx) && (c > pow_ny) && (c > pow_nz)){
			////Este caso funciona quando tem-se duas piramides fora da celula
				d = solver_equation_third_order(volume, n_z, n_y, n_x, dz, dy, dx);
				d = -sign(n) * trans_p0_to_center(Delta, Normal, d);
				real d2 = solver_equation_third_order_3(volume, n_z, n_y, n_x, dz, dy, dx);
				if(d2<n_x*dx && d2>n_y*dy && d2>n_z*dz){
						d = -sign(n) * trans_p0_to_center(Delta, Normal, d2);
				}
				return d;
			}else if((c > pow_nx) && (c > pow_ny) && (c > pow_nz)){
			////Este caso funciona quando tem-se tres piramides fora da celula
				d = solver_equation_third_order_2(volume, n_x, n_y, n_z, dx, dy, dz);
				d = -sign(n) * trans_p0_to_center(Delta, Normal, d);
				return d;
			}
			
		//Parte superior================================================
		//==============================================================
		}else if(nz < 0.0 && VOLUME > 0.5*(dx*dy*dz) || nz > 0.0 && VOLUME < 0.5*(dx*dy*dz)){
			
			volume = (nz > 0.0 && VOLUME < 0.5*(dx*dy*dz)) ? VOLUME : dx*dy*dz - VOLUME;
			
			f = volume/(dx*dy*dz);
			real c = 6.0*n_x*n_y*n_z*volume;
			
			if((c < pow_nx) && (c < pow_ny) && (c < pow_nz)){
			////Este caso funciona quando a piramide esta completamente dentro da celula	
				d = cbrt(6*volume*n_x*n_y*n_z);
				d = sign(n) * trans_p0_to_center(Delta, Normal, d);
				return d;
			}else if((c > pow_nx) && (c < pow_ny) && (c < pow_nz)){
			////Este caso funciona quando tem-se apenas uma piramide fora da celula
				d = solver_equation_second_order(volume, n_x, n_y, n_z, dx, dy, dz);
				d = sign(n) * trans_p0_to_center(Delta, Normal, d);
				return d;
			}else if((c < pow_nx) && (c > pow_ny) && (c < pow_nz)){
			////Este caso funciona quando tem-se apenas uma piramide fora da celula
				d = solver_equation_second_order(volume, n_y, n_x, n_z, dy, dx, dz);
				d = sign(n) * trans_p0_to_center(Delta, Normal, d);
				return d;
			}else if((c < pow_nx) && (c < pow_ny) && (c > pow_nz)){
			////Este caso funciona quando tem-se apenas uma piramide fora da celula
				d = solver_equation_second_order(volume, n_z, n_y, n_x, dz, dy, dx);
				d = sign(n) * trans_p0_to_center(Delta, Normal, d);
				return d;
			}else if((c > pow_nx) && (c > pow_ny) && (c < pow_nz)){
			////Este caso funciona quando tem-se duas piramides fora da celula
				d = solver_equation_third_order(volume, n_x, n_y, n_z, dx, dy, dz);
				d = sign(n) * trans_p0_to_center(Delta, Normal, d);
				
				real d2 = solver_equation_third_order_3(volume, n_x, n_y, n_z, dx, dy, dz);
				if(d2>n_x*dx && d2>n_y*dy && d2<n_z*dz){
						d = sign(n) * trans_p0_to_center(Delta, Normal, d2);
				}
				return d;
			}else if((c > pow_nx) && (c < pow_ny) && (c > pow_nz)){
			////Este caso funciona quando tem-se duas piramides fora da celula
				d = solver_equation_third_order(volume, n_x, n_z, n_y, dx, dz, dy);
				d = sign(n) * trans_p0_to_center(Delta, Normal, d);
				
				real d2 = solver_equation_third_order_3(volume, n_x, n_z, n_y, dx, dz, dy);
				if(d2>n_x*dx && d2<n_y*dy && d2>n_z*dz){
						d = sign(n) * trans_p0_to_center(Delta, Normal, d2);
				}
				return d;
			}else if((c < pow_nx) && (c > pow_ny) &&
			(c > pow_nz)){
			////Este caso funciona quando tem-se duas piramides fora da celula
				d = solver_equation_third_order(volume, n_z, n_y, n_x, dz, dy, dx);
				d = sign(n) * trans_p0_to_center(Delta, Normal, d);
				
				real d2 = solver_equation_third_order_3(volume, n_z, n_y, n_x, dz, dy, dx);
				if(d2<n_x*dx && d2>n_y*dy && d2>n_z*dz){
						d = sign(n) * trans_p0_to_center(Delta, Normal, d2);
				}
				return d;
			}else if((c > pow_nx) && (c > pow_ny) &&
			(c > pow_nz)){
			////Este caso funciona quando tem-se tres piramides fora da celula
				d = solver_equation_third_order_2(volume, n_x, n_y, n_z, dx, dy, dz);
				d = sign(n) * trans_p0_to_center(Delta, Normal, d);
				return d;
			}
			
		}
	}
}
//======================================================================

real parallel_case_volume(Point Normal,Point Delta,real d_from_center,real tol_n){
	real nx = Normal[0], ny = Normal[1], nz = Normal[2];
	real dx = Delta[0], dy = Delta[1], dz = Delta[2];
	int n_x, n_y, n_z, n_1, n_2, n_3;
	real volume, d;
	
	n_1 = fabs(Normal[0]);
	n_2 = fabs(Normal[1]);
	n_3 = fabs(Normal[2]);
	Normal[0] = -n_1;
	Normal[1] = -n_2;
	Normal[2] = -n_3;
	
	d = fabs(d_from_center);
	d = trans_center_to_p0(Delta,Normal,d);
	
	if(fabs(nx)<tol_n && fabs(ny)<tol_n){
		n_x = 0;
		n_y = 0;
		n_z = 1;
	}else if(fabs(nx)<tol_n && fabs(nz)<tol_n){
		n_x = 0;
		n_y = 1;
		n_z = 0;
	}else if(fabs(ny)<tol_n && fabs(nz)<tol_n){
		n_x = 1;
		n_y = 0;
		n_z = 0;
	}else{
		printf("Nenhuma das condicoes foi atendida");
	}
		
	
	volume = (n_x==1) ? d*(dy*dz) : ((n_y==1) ? d*(dx*dz) : d*(dy*dx));

	return volume;
}
//======================================================================

real volume_left_line_origin_center(Point Normal,Point Delta,real d_from_center){
	
	real aux_d = d_from_center;
	d_from_center = fabs(d_from_center);
	real d, volume;
	real dx  = Delta[0], dy  = Delta[1], dz  = Delta[2];
	real nx = Normal[0], ny = Normal[1], nz = Normal[2];
	real n_x = fabs(nx), n_y = fabs(ny), n_z = fabs(nz);
	real tol_n = 1e-8;
	
	real dmax   = n_x * 0.5 * Delta[0] + n_y * 0.5 * Delta[1] + n_z * 0.5 * Delta[2];
	real volmax = dx*dy*dz;
	
	if((d_from_center >= dmax) && aux_d<0) {
		 volume = dx*dy*dz;
	} else if((d_from_center >= dmax) && aux_d>0) {
		 volume = 0.0;
	}else{
		if((n_x<=tol_n && n_y<=tol_n ) || (n_x<=tol_n && n_z<=tol_n ) || (n_y<=tol_n && n_z<=tol_n )) {
			volume = parallel_case_volume(Normal, Delta, d_from_center, tol_n);
			if(aux_d<0){
				volume=dx*dy*dz-volume;
			}
		}else if(n_x<=tol_n && (n_y>tol_n && n_z>tol_n) || n_y<=tol_n  && 
		(n_x>tol_n && n_z>tol_n) || n_z<=tol_n  && (n_y>tol_n && n_x>tol_n)){
			
			real n1, n2, n3, n_1, n_2, n_3;
			n1 = (n_x < tol_n) ? ny : nx;
			n2 = (n_x < tol_n || n_y < tol_n) ? nz : ny;
			n3 = 0.0;
			n_1 = (n_x < tol_n) ? -n_y : -n_x;
			n_2 = (n_x < tol_n || n_y < tol_n) ? -n_z : -n_y;
			n_3 = 0.0;
			Normal[0] = n_1;
			Normal[1] = n_2;
			Normal[2] = n_3;
			
			d = trans_center_to_p0(Delta, Normal, d_from_center);
			volume = volume_3D(n_1, n_2, n_3, dx, dy, dz, d, aux_d);
		}else if(n_x>tol_n && n_y>tol_n && n_z>tol_n){
			
			real n1, n2, n3, n_1, n_2, n_3;
			n_1 = -n_x;
			n_2 = -n_y;
			n_3 = -n_z;
			Normal[0] = n_1;
			Normal[1] = n_2;
			Normal[2] = n_3;
			
			d = trans_center_to_p0(Delta, Normal, d_from_center);
			volume = volume_3D(nx, ny, nz, dx, dy, dz, d, aux_d);
		}
	}
	
	volume = minimo(volume,volmax);
	
	return volume;
}
//======================================================================
real maximo(real nx, real ny, real nz){
		real max;		
		max = (nx > ny && nx > nz) ? nx : ((ny > nx && ny > nz) ? ny : nz);
		
		return max;
}
//======================================================================
real minimo(real vol1, real vol2){
		real min;
		min = (vol1 < vol2) ? vol1 : vol2;
		
		return min;
}
//======================================================================
void higflow_compute_distance_multiphase_3D(higflow_solver *ns) {
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
			
			Point Normal;
			Normal[0] = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
			Normal[1] = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
			Normal[2] = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpnormal[2], ns->ed.stn);
						
			real tol_n = 1e-8;
			if (fabs(Normal[0]) < tol_n && fabs(Normal[1]) < tol_n && fabs(Normal[2]) < tol_n){
				continue;
			}
			
			//real frac = compute_value_at_point(sdp, center, center, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			real volume = frac*delta[0]*delta[1]*delta[2];
			
			real distance3D = distance_from_center_3D(Normal,delta,volume);
			
			dp_set_value(ns->ed.mult.dpdistance, clid, distance3D);

			//arquivoDist(NULL,center[0],center[1],distance);
		}
		// Destroy the iterator
		higcit_destroy(it);
		// Sync the distributed pressure property
		dp_sync(ns->ed.mult.dpdistance);
	}
}
#endif
