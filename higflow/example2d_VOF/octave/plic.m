1;
function [value] = solve_equation(Rx,Ry,nx,ny,area,dx,dy)

	A = (1-Rx-Ry)/abs(nx*ny);
	B = 2.0*(Rx*dx/abs(ny)+Ry*dy/abs(nx));
	C = -2.0*area-Rx*dx^2*abs(nx/ny)-Ry*dy^2*abs(ny/nx);
	
	if(A==0)
		value=-C/B;
		return;
	end
	
	x1=(-B+sqrt(B^2-4.0*A*C))/(2.0*A);
	x2=(-B-sqrt(B^2-4.0*A*C))/(2.0*A);
  
	value=[x1 x2];
end

function[d]=parallel_left_line_origin_center_DISTANCE(nx,ny,dx,dy,area,tol_n)
	
	if(abs(nx)<tol_n)
		n_x = 0;
		if(ny>0)
			n_y = 1;
		else
			n_y = -1;
		end
	else
		n_y=0;
		if(nx>0)
			n_x=1;
		else
			n_x=-1;
		end
	end
	
	if(n_x==-1)
		d=-area/dy;
	elseif(n_x==1)
		d=-(area/dy-dx);
	elseif(n_y==-1)
		d=-area/dx;
	elseif(n_y==1)
		d=-(area/dx-dy);
	end
	
end

function[d]=distance_from_center(nx,ny,dx,dy,AREA)
	
	tol_area=1e-8;
	tol_n=1e-8;
	
	if(abs(nx)<tol_n||abs(ny)<tol_n)
		
		d=parallel_left_line_origin_center_DISTANCE...
		(nx,ny,dx,dy,AREA,tol_n);
		d=trans_bl_2_center(dx,...
		dy,nx,ny,d);
		
	elseif(nx*ny>0)
		
		dLT=nx*(-0.5*dx)+ny*(0.5*dy);
		aLT=area_left_line_origin_center...
		(nx,ny,dx,dy,dLT);
		
		dRB=nx*(0.5*dx)+ny*(-0.5*dy);
		aRB=area_left_line_origin_center...
		(nx,ny,dx,dy,dRB);
		
		area=AREA;
		
		if(abs(area-aLT)<=tol_area)
			d=dLT;
			return;
		elseif(abs(area-aRB)<=tol_area)
			d=dRB;
			return;
		end
		
		if(nx>0)
			area=dx*dy-AREA;
			aLT=dx*dy-aLT;
			aRB=dx*dy-aRB;
		end
		
		amax=max(aLT,aRB);
		amin=min(aLT,aRB);
		
		if(area>amax)
			Rx=1;Ry=1;
			d=solve_equation(Rx,Ry,nx,ny,...
			area,dx,dy);
			d=min(d);
		elseif(area<amax&&area>amin)
			if(aRB>aLT)
				Rx=0;Ry=1;
			else
				Rx=1;Ry=0;
			end
			d=solve_equation(Rx,Ry,nx,ny,...
			area,dx,dy);
		elseif(area<amin)
			Rx=0;Ry=0;
			d=solve_equation(Rx,Ry,nx,ny,...
			area,dx,dy);
			d=max(d);
		end
		d=sign(nx)*d;
		d=trans_bl_2_center(dx,...
		dy,nx,ny,d);
	elseif(nx*ny<0)
		
		dLB=nx*(0.5*dx)+ny*(0.5*dy);
		aLB=area_left_line_origin_center...
		(nx,ny,dx,dy,dLB);
		
		dRT=nx*(-0.5*dx)+ny*(-0.5*dy);
		aRT=area_left_line_origin_center...
		(nx,ny,dx,dy,dRT);
		printf('aRTold=%f aLBold=%f\n',aRT,aLB);
		area=AREA;
		
		if(abs(area-aRT)<=tol_area)
			d=dRT;
			return;
		elseif(abs(area-aLB)<=tol_area)
			d=dLB;
			return;
		end
		
		if(nx<0)
			area=dx*dy-AREA;
			aRT=dx*dy-aRT;
			aLB=dx*dy-aLB;
		end
		
		amax=max(aRT,aLB);
		amin=min(aRT,aLB);
		printf('aRTold2=%f aLBold2=%f\n',aRT,aLB);
    
		if(area>amax)
			Rx=1;Ry=1;
			d=solve_equation(Rx,Ry,nx,ny,...
			area,dx,dy);
			d=min(d); 
		elseif(area>amin&&area<amax)
			if(aRT>aLB)
				Rx=0;Ry=1;
			else
				Rx=1;Ry=0;
			end
			d=solve_equation(Rx,Ry,nx,ny,...
			area,dx,dy);
      printf('dRT=%f dLB=%f\n',dRT,dLB);
      printf('aRT=%f aLB=%f\n',aRT,aLB);
      printf('Rx=%d Ry=%d nx=%f ny=%f\n',Rx,Ry,nx,ny);
      printf("Passei Aquisson: dvec[0]=%f dvec[1]=%f area=%f"
      ,d,d,area);
		elseif(area<amin)
			Rx=0;Ry=0;
			d=solve_equation(Rx,Ry,nx,ny,...
			area,dx,dy);
			d=max(d);
		end
		d=-sign(nx)*d;
		d=trans_br_2_center(dx,...
		dy,nx,ny,d);
		printf("d=%f\n",d);
	end
	
	
	
	
end


function[area]=parallel_left_line_origin_center_AREA(nx,ny,dx,dy,d_from_center,tol_n)
	
	if(abs(nx)<tol_n)
		n_x=0;
		if(ny>0)
			n_y=1;
		else
			n_y=-1;
		end
	else
		n_y=0;
		if(nx>0)
			n_x=1;
		else
			n_x=-1;
		end
	end
	
	d = trans_center_2_bl(dx,dy,nx,ny,d_from_center);
	
	if(n_x==-1)
		area=-d*dy;
	elseif(n_x==1)
		area=dy*(dx-d);
	elseif(n_y==-1)
		area=-d*dx;
	elseif(n_y==1)
		area=dx*(dy-d);
	end
	
end

function [area] = area_left_line_origin_center(nx,ny,dx,dy,d_from_center)
	
	n_x   = abs(nx);
	n_y   = abs(ny);
	
	tol_n = 1e-8;
	
	srt = left_right([0.5*dx 0.5*dy],nx,ny,d_from_center);
	srb = left_right([0.5*dx -0.5*dy],nx,ny,d_from_center);
	slt = left_right([-0.5*dx 0.5*dy],nx,ny,d_from_center);
	slb = left_right([-0.5*dx -0.5*dy],nx,ny,d_from_center);
	
	if(srt<=0&&srb<=0&&slt<=0&&slb<=0)
		area = dx*dy;
		%disp('$$$$$$$$$$  FULL $$$$$$$$$$');
		return
	elseif(srt>=0&&srb>=0&&slt>=0&&slb>=0)
		area = 0.0;
		%disp('$$$$$$$$  EMPTY  $$$$$$$$$$' );
		return;
	end
	
	if(n_x <= tol_n||n_y <= tol_n)
		[area] = parallel_left_line_origin_center_AREA(nx,ny,dx,dy,d_from_center,tol_n);
		return;
	end
	
	if(nx*ny>0)
		d = trans_center_2_bl(dx,dy,nx,ny,d_from_center);
	else
		d = trans_center_2_br(dx,dy,nx,ny,d_from_center);
	end
	
	d      = abs(d);
	
	aux1   = 1/(2.0*n_x*n_y);
	auxx   = d - n_x*dx;
	auxx_1 = (auxx)^2.0;
	auxy   = d - n_y*dy;
	auxy_1 = (auxy)^2.0;
	
	area   = aux1*(d^2 - R(auxx)*auxx_1 - R(auxy)*auxy_1);
	
	if(ny>0)
		area=dx*dy-area;
	end
	
end

function [side] = left_right(P,n_x,n_y,d)
	side = sign(-n_x*P(1) - n_y*P(2) + d);
end

function [d_from_bl] = trans_center_2_bl(dx,dy,n_x,n_y,d_from_center)
	d_from_bl = d_from_center +	n_x*0.5*dx + n_y*0.5*dy;
end

function [d_from_center] = trans_bl_2_center(dx,dy,nx,ny,d_from_bl)
	d_from_center = d_from_bl - nx*0.5*dx - ny*0.5*dy;
end

function [d_from_br] = trans_center_2_br(dx,dy,n_x,n_y,d_from_center)
	d_from_br = d_from_center - n_x*0.5*dx + n_y*0.5*dy;
end

function[d_from_center]=trans_br_2_center(dx,dy,n_x,n_y,d_from_br)
	d_from_center = d_from_br	+ n_x*0.5*dx - n_y*0.5*dy;
end

function [value] = R(x)
	if (x<= 0.0)
		value = 0.0;
	else
		value = 1.0;
	end	
end