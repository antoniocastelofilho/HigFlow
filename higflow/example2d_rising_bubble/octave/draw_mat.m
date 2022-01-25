function[]=draw_mat()

	global nx ny XXh YYh hplot
	
	xf=1; yf=1;
	nx=120;ny=120;dx=xf/nx;dy=yf/ny;

	x = 0:dx:xf; y = 0:dy:yf;
	
	xh=dx/2:dx:xf-dx/2; yh=dy/2:dy:yf-dy/2;
	
	[XX YY] = meshgrid(x,y);
	[XXh YYh] = meshgrid(xh,yh);
	
	dt=0.0001;
	xc=1.0;yc=1.0;r=1/6;
	
	x=0.5*x;y=0.5*y;xh=0.5*xh;yh=0.5*yh;dx=0.5*dx;dy=0.5*dy;
	XX=0.5*XX;YY=0.5*YY;XXh=0.5*XXh;YYh=0.5*YYh;
	xc=0.5*xc;yc=0.5*yc;r=0.5*r;
	xf=0.5*xf;yf=0.5*yf;
	
	i=0:100-1;
	i=0;
	
	for step=i
		
	figure(2);shg;clf;hold on;
%shg;clf;hold on;
		
		%plot(XX,YY,'k');
		%plot(XX',YY','k');
		x01=0:0.1:xf;
		y01=0:0.1:yf;
		
		[xx01 yy01]=meshgrid(x01,y01);
		hplot=plot(xx01,yy01,'k','linestyle','--');
		hplot=[hplot;plot(xx01',yy01','k','linestyle','--')];

		hplot=[hplot;plot([0 xf],[0.5*yf 0.5*yf],'color','y','linewidth',2)];
		
		%drawCircle();
		
		str_title=strcat('step:   ',num2str(step));
		%title(str_title);
		
		%[frac vx vy d normal kappa iforce gp]=loadFiles(step);
		[frac d normal]=loadFiles(step);
		frac(:,1:2)=0.5*frac(:,1:2);
		T=step*dt;
		draw_interface(x,y,xh,yh,dx,dy,max(size(frac)),frac,normal,0.5*d);
		%draw_exact_interface(xc,yc,r,T);
		%draw_press(frac,iforce,max(size(frac)),dx,dy, gp)
		
		%rotate(hplot,[0 0 1],90);
		
		%draw_curv(frac,kappa,max(size(frac)),dx,dy);
		
		
		
	end
	%draw_frac(frac,max(size(frac)),dx,dy,0);
	%draw_v(frac,vx,max(size(frac)),dx,dy,1)
	%draw_v(frac,vy,max(size(frac)),dx,dy,2)
	
	%axis([0.2 0.8 0.2 0.8])
	%axis equal;
	%print (figure(1),'circulo.pdf', '-dpdflatexstandalone');
end

%function[frac vx vy d normal kappa press iforce]=loadFiles(step)
%function[frac vx vy d normal kappa iforce gp]=loadFiles(step)
function[frac d normal]=loadFiles(step)
	
	
	str=strcat('../DATA/',num2str(step),'_frac_visc.txt');
	frac=load(str,'-ascii');
	
%	str=strcat('../DATA/',num2str(step),'_vx.txt');
%	vx=1;%vx=load(str,'-ascii');
	
%	str=strcat('../DATA/',num2str(step),'_vy.txt');
%	vy=1;%load(str,'-ascii');
	
	str=strcat('../DATA/',num2str(step),'_d.txt');
	d=load(str,'-ascii');
	
	str=strcat('../DATA/',num2str(step),'_norm.txt');
	normal=load(str,'-ascii');
	
%	str=strcat('../DATA/',num2str(step),'_curv.txt');
%	kappa=1;%load(str,'-ascii');
	
%	str=strcat('../DATA/',num2str(step),'_press.txt');
%	press=load(str,'-ascii');
	
%	str=strcat('../DATA/',num2str(step),'_if.txt');
%	iforce=load(str,'-ascii');
	
%	str=strcat('../DATA/',num2str(step),'_gp.txt');
%	gp = 1;
%	gp=load(str,'-ascii');
	
	
end

function[]=drawCircle()
  
  xc=0.5;yc=0.5;r=0.25;
  t=[0:0.01:2*pi];
  x=xc+r*cos(t);
  y=yc+r*sin(t);
  plot(x,y,'color','c','linewidth',2);
  
end

function[]=draw_frac(frac,n,dx,dy,ind)
	for i = 1:n
		x=frac(i,1);
		y=frac(i,2);
		f=frac(i,3);
		if ((f == 0 || f == 1))
			continue;
		end
		text(x,y+0.0*dy,num2str(f),'color','b');
		
		if(ind==1)
			text(x,y+0.25*dy,num2str(i-1),'color','r');
		end
		
	end
end

function[]=draw_v(frac,v,n,dx,dy,dim)
	
	color=0;
	d=[0 0];
	
	if(dim==1)
		d(dim)=dx;color='c';
	else
		d(dim)=dy;color='m';
	end
	
	for i = 1:n
		x=frac(i,1);
		y=frac(i,2);
		f=frac(i,3);
		ul=v(i,1);
		ur=v(i,2);
		if ((f == 0 || f == 1))
			continue;
		end
		pr=[x y] +0.5*d;
		pl=[x y] -0.5*d;
		
		text(pr(1),pr(2),num2str(ur),'color',color);
		text(pl(1),pl(2),num2str(ul),'color',color);
		
		if(dim==1)
			text(x+0.25*dx,y,num2str((ur-ul)/dx),'color','k');
		else
			text(x-0.25*dx,y,num2str((ur-ul)/dy),'color','k');
		end
		
	end
	
end

function[]=draw_curv(frac,kappa,m,dx,dy)
	
	for i = 1:m
		x=frac(i,1);
		y=frac(i,2);
		f=frac(i,3);
		
		curv=kappa(i,1);
		
		if(f==1||f==0)
			continue;
		end
		
		text(x,y-0.25*dy,num2str(curv),'color','g');
		
		
	end
end


function[]=draw_press(frac,iforce,m,dx,dy, gp)
	
	global nx ny XXh YYh
	d = abs(iforce - gp);
	d1=d(:,1).*d(:,1) + d(:,2).*d(:,2);
	
	
	sify=reshape(iforce(:,1),ny,nx);
	sifx=reshape(iforce(:,2),ny,nx);
	
	sgp = size(gp);
	
	gy=reshape(gp(:,1),ny,nx);
	gx=reshape(gp(:,2),ny,nx);
	
	normif=sqrt(sify.^2  + sifx.^2);
	
	normgp=sqrt(gy.^2  + gx.^2);
	
%	pressnew=reshape(press,ny,nx);
	
	sxh=size(XXh);
	syh=size(YYh);
%	sp=size(pressnew);
	
	
	
%	normifinf = (max(max(abs(normif))));
%	normgpinf = max((max(abs(normgp))));

	
	
	surfc(XXh,YYh,normgp);
	title(sprintf('Pressure gradient'));
	figure(2);clf;
	surfc(XXh,YYh,normif);
	title(sprintf('Interfacial force'));
%	figure(3);clf;
%	surfc(XXh,YYh, normifgp);
%	title(sprintf('Capillary pressure'));
	%contourf(XXh,YYh,reshape(press,ny,nx));
	%colorbar();
	%	disp('pausing...');
%	pause;
	return;
	
	press_aux=reshape(press,nx,ny);
	for i=1:0
		for j=1:ny
		plot(XX(i),YY(j));
		pause;
		text(XX(i),YY(j)-0.25*dy,num2str(press_aux(i,j)));
		pause;
		end
	end
	
	
	
	
	for i = 1:m
		x=frac(i,1);
		y=frac(i,2);
		f=frac(i,3);
		
		p=press(i,1);
		
		if(f==1||f==0)
			%continue;
		end
		
		text(x,y-0.25*dy,num2str(p),'color','g');
		plot(x,y,'color','g');
		pause;
		
	end
end

function[ ]=draw_interface(X,Y,Xh,Yh,dx,dy,n,Frac,Normal,distance)
	
	for i=1:n
		
		nx=Normal(i,1);ny=Normal(i,2);
		center=[Frac(i,1),Frac(i,2)];
		dc=distance(i);
		drawCell(nx,ny,dx,dy,center,dc,Frac(i,3),'k','k','k');
		
	end
	
end

function[]=draw_exact_interface(xc,yc,r,T)
	
	
	
	s=0:0.2:2*pi;
	x0=xc+r*cos(s);
	y0=yc+r*sin(s);
	
	u=@(t)1;u=@(t)sin(pi*t);
	v=@(t)0;v=@(t)0.4*cos(pi*t);
	
	x= x0 + quad(u,0,T);
	y= y0  + quad(v,0,T);
	
	%plot(x,y,'ob');
	
	y=r*sin(s)+xc;
	x=r*cos(s) +yc -4.0*y.*(y - 1.0).*T;
	
	plot(x,y,'ob');
	
end


function[]=drawCell(nx,ny,dx,dy,center,dc,frac,colorp,colorn,colora)
	global hplot;
	
	if(frac==0||frac==1)
		return;
	end
		
	if(nx*ny>0)
		
		[d] = trans_center_2_bl(dx,dy,nx,ny,dc);
		
		Rx=R(abs(d)-abs(nx)*dx);
		Ry=R(abs(d)-abs(ny)*dy);
		
		if(Rx==1&&Ry==1)
			P0=center+[(dc-0.5*ny*dy)/nx,0.5*dy]; 
			P1=center+[0.5*dx,(dc-0.5*nx*dx)/ny];
		elseif(Rx==0&&Ry==1)
			P0=center+[(dc-0.5*ny*dy)/nx,0.5*dy];
			P1=center+[(dc+0.5*ny*dy)/nx,-0.5*dy];
		elseif(Rx==1&&Ry==0)
			P0=center+[-0.5*dx,(dc+0.5*nx*dx)/ny];
			P1=center+[0.5*dx,(dc-0.5*nx*dx)/ny];
		elseif(Rx==0&&Ry==0)
			P0=center+[-0.5*dx,(dc+0.5*nx*dx)/ny];
			P1=center+[(dc+0.5*ny*dy)/nx,-0.5*dy];
		end
		vx=[P0(1),P1(1)];vy=[P0(2),P1(2)];
		PM=0.5*(P0+P1);
		hplot=[hplot;plot(vx,vy,'color',colorp,'linewidth',1.5)];
		axis equal;
		PH=PM+min(dx,dy)*[nx ny];
		%plot([PM(1),PH(1)],[PM(2),PH(2)],colora);
		%plot(PH(1),PH(2),'marker','.','color',colora);
		
	elseif(nx*ny<0)
		
		[d] = trans_center_2_br(dx,dy,nx,ny,dc);
		
		Rx=R(abs(d)-abs(nx)*dx);
		Ry=R(abs(d)-abs(ny)*dy);
		if(Rx==1&&Ry==1)
			P0=center+[(dc-0.5*ny*dy)/nx,0.5*dy]; 
			P1=center+[-0.5*dx,(dc+0.5*nx*dx)/ny];
		elseif(Rx==0&&Ry==1)
			P0=center+[(dc-0.5*ny*dy)/nx,0.5*dy]; 
			P1=center+[(dc+0.5*ny*dy)/nx,-0.5*dy];
		elseif(Rx==1&&Ry==0)
			P0=center+[0.5*dx,(dc-0.5*nx*dx)/ny];
			P1=center+[-0.5*dx,(dc+0.5*nx*dx)/ny];
		elseif(Rx==0&&Ry==0)
			P0=center+[0.5*dx,(dc-0.5*nx*dx)/ny];
			P1=center+[(dc+0.5*ny*dy)/nx,-0.5*dy];
		end
		vx=[P0(1),P1(1)];vy=[P0(2),P1(2)];
		PM=0.5*(P0+P1);
		hplot=[hplot;plot(vx,vy,'color',colorn,'linewidth',1.5)];
		axis equal;
		PH=PM+min(dx,dy)*[nx ny];
		%plot([PM(1),PH(1)],[PM(2),PH(2)],colora);
		%plot(PH(1),PH(2),'marker','.','color',colora);
	end
	
end

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
		
		d=parallel_left_line_origin_center_DISTANCE(nx,ny,dx,dy,AREA,tol_n);
		d=trans_bl_2_center(dx,dy,nx,ny,d);
		
	elseif(nx*ny>0)
		
		dLT=nx*(-0.5*dx)+ny*(0.5*dy);
		aLT=area_left_line_origin_center(nx,ny,dx,dy,dLT);
		
		dRB=nx*(0.5*dx)+ny*(-0.5*dy);
		aRB=area_left_line_origin_center(nx,ny,dx,dy,dRB);
		
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
			d=solve_equation(Rx,Ry,nx,ny,area,dx,dy);
			d=min(d);
		elseif(area<amax&&area>amin)
			if(aRB>aLT)
				Rx=0;Ry=1;
			else
				Rx=1;Ry=0;
			end
			d=solve_equation(Rx,Ry,nx,ny,area,dx,dy);
		elseif(area<amin)
			Rx=0;Ry=0;
			d=solve_equation(Rx,Ry,nx,ny,area,dx,dy);
			d=max(d);
		end
		d=sign(nx)*d;
		d=trans_bl_2_center(dx,dy,nx,ny,d);
	elseif(nx*ny<0)
		
		dLB=nx*(0.5*dx)+ny*(0.5*dy);
		aLB=area_left_line_origin_center(nx,ny,dx,dy,dLB);
		
		dRT=nx*(-0.5*dx)+ny*(-0.5*dy);
		aRT=area_left_line_origin_center(nx,ny,dx,dy,dRT);
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
			d=solve_equation(Rx,Ry,nx,ny,area,dx,dy);
			d=min(d); 
		elseif(area>amin&&area<amax)
			if(aRT>aLB)
				Rx=0;Ry=1;
			else
				Rx=1;Ry=0;
			end
			d=solve_equation(Rx,Ry,nx,ny,area,dx,dy);
      printf('dRT=%f dLB=%f\n',dRT,dLB);
      printf('aRT=%f aLB=%f\n',aRT,aLB);
      printf('Rx=%d Ry=%d nx=%f ny=%f\n',Rx,Ry,nx,ny);
      printf('Passei Aquisson: dvec[0]=%f dvec[1]=%f area=%f',d,d,area);
		elseif(area<amin)
			Rx=0;Ry=0;
			d=solve_equation(Rx,Ry,nx,ny,area,dx,dy);
			d=max(d);
		end
		d=-sign(nx)*d;
		d=trans_br_2_center(dx,dy,nx,ny,d);
		printf('d=%f\n',d);
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
