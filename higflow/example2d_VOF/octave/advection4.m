function[]=advection4()
	source('plic.m');
	source('reconstruction4.m');
	source('draws.m');
	
	
	global Frac Frac_aux Normalx Normaly 
	global NormalAux Kappa;
	
	[X,Y,Xh,Yh]=create_domain();
	
	[Nx Ny]=size(Frac);
	Frac_aux=Frac;
	
	delta = [X(2)-X(1) Y(2)-Y(1)];dx=delta(1),dy=delta(2);
	
	calculation_normal_curvature(Xh,Yh,delta,X,Y);
	build_distance(Nx,Ny,dx,dy);
	%draw_interface(X,Y,Xh,Yh,dx,dy,Nx,Ny);
	%drawFrac(Xh,Yh);
	
	dt = 0.2*dx;
	Nt = 1280;
	
	disp('ADVECTION 4');
	for K=1:Nt
		
		%figure(1);clf;hold on;
		%title('*******OLD*********');
		
		%draw_interface(X,Y,Xh,Yh,dx,dy,Nx,Ny);
		%drawFrac(Xh,Yh);
		%draw_x_flux_rectangle(Xh,Yh,(K-1)*dt,dt);
		
		
		%figure(2);shg;clf;hold on;
		%title('******* ADV-4*********');
		
		disp('=================');
		disp('advecting fraction...');fflush(stdout);
		for I = 2:Nx-1
			for J = 2:Ny-1
				x = Xh(I);y = Yh(J);
				ur = u(x + 0.5*dx,y,(K - 1)*dt); 
				ul = u(x - 0.5*dx,y,(K - 1)*dt);
				vr = v(x,y + 0.5*dy,(K - 1)*dt); 
				vl = v(x,y - 0.5*dy,(K - 1)*dt);
				
				advecx(I,J,ur,ul,dx,dy,dt,x,y);
			end
		end
		disp('DONE advection...');fflush(stdout);
		
		disp('calculating interface...');fflush(stdout);
		Normalx=zeros(Nx,Ny);Normaly=zeros(Nx,Ny);
		Kappa=zeros(Nx,Ny);NormalAux=zeros(Nx,Ny);
		Frac=Frac_aux;
		calculation_normal_curvature(Xh,Yh,delta,X,Y);
		build_distance(Nx,Ny,dx,dy);
		disp('DONE interface');fflush(stdout);
		%drawCircle(X,Y);
		%drawDomain(X,Y);
		%draw_interface(X,Y,Xh,Yh,dx,dy,Nx,Ny);
		%drawFrac(Xh,Yh);
		%drawKappa(X,Y,Xh,Yh);
		disp('saving...');fflush(stdout);
		save FRAC4_512x256.txt Frac;
		
		
		printf('K = %d\n\n',K);
		disp('=================');
		
		fflush(stdout);
		%pause;
		
	end
	
	%save Frac2.txt Frac;
	figure(2);shg;clf;hold on;
	drawCircle(X,Y);
	drawDomain(X,Y);
	draw_interface(X,Y,Xh,Yh,dx,dy,Nx,Ny);
	drawFrac(Xh,Yh);
	drawKappa(X,Y,Xh,Yh);
	title('256x256');
	
end

%%%%%%%%%%%%%%%%%%%%%%%
function[]=advecx(I,J,ur,ul,dx,dy,dt,x,y)
	
	global Frac Frac_aux distance Normalx Normaly;
	
	tol_u = 1e-8;
	%%% advection right face
	if (abs(ur) > tol_u)
			
		if(Frac(I,J)==1||Frac(I,J)==0||(Normalx(I,J)==0&&Normaly(I,J)==0))
			
			Ar=Frac(I,J)*dy*ur*dt;
		else
			d = distance(I,J);
			n_x_r = Normalx(I,J);
			n_y_r = Normaly(I,J);
			d_center_new_r = d - 0.5*(dx - ur*dt)*n_x_r;
			Ar = area_left_line_origin_center...
			(n_x_r,n_y_r,ur*dt,dy,d_center_new_r);
		end
	else
		
	end
	
	%%% advection left face
	if (abs(ul) > tol_u)
			
		if(Frac(I-1,J)==0||Frac(I-1,J)==1||(Normalx(I-1,J)==0&&Normaly(I-1,J)==0))
			
			Al=Frac(I-1,J)*dy*ul*dt;
			
		else
			d = distance(I-1,J);
			n_x_l = Normalx(I-1,J);
			n_y_l = Normaly(I-1,J);
			d_center_new_l = d - 0.5*(dx - ul*dt)*n_x_l;
			Al = area_left_line_origin_center...
			(n_x_l,n_y_l,ul*dt,dy,d_center_new_l);
		end
	else
		
	end
	
	
	Fracr=Ar/(ur*dt*dy);
	Fracl=Al/(ul*dt*dy);
	Frac_aux(I,J)=Frac(I,J) - dt*(Fracr*ur - Fracl*ul)/dx;
	Frac_aux(I,J)=fraction_correction(Frac_aux(I,J),I,J);
	
	
end

%%%%%%%%%%%%%%%%%%%%%%%
function[ ]=advecy(I,J,vr,vl,dx,dy,dt,x,y)
	global Frac distance Normalx Normaly;
	
	return;
	
	%%% advection right face
	
	%%% advection left face
	
end

%%%%%%%%%%%%%%%%%%%%%%%%
function[Frac_new] = fraction_correction(frac,I,J)
	
	tol_frac = 1e-6;
	Frac_new = frac;
	
	if (Frac_new < tol_frac)
		Frac_new = 0.0;
	elseif (Frac_new > 1.0 - tol_frac)
		Frac_new = 1.0;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%
function [value]=u(x,y,t)
	value = 1.0;
	value = -4.0*y*(y - 1.0);
	
end

%%%%%%%%%%%%%%%%%%%%%%%%%
function [value]=v(x,y,t)
	value = 0.0;
end