function[]=draw_compare()
	source('plic.m');
	global nx ny XXh YYh hplot
	xf=1; yf=1;
	%close all;
	nx=16;ny=16;dx=xf/nx;dy=yf/ny;

	x = 0:dx:xf;
	y = 0:dy:yf;
	
	xh=dx/2:dx:xf-dx/2;
	yh=dy/2:dy:yf-dy/2;
	
	[XX YY] = meshgrid(x,y);
	[XXh YYh] = meshgrid(xh,yh);
	dt=0.0001;
	xc=0.5;yc=0.75;r=0.15;
	
	i=0;
	figure(1);
	shg;clf;hold on;
	x01=0:xf:xf;
	y01=0:yf:yf;
	
	[xx01 yy01]=meshgrid(x01,y01);
	hplot=plot(xx01,yy01,'k','linestyle','-');
	hplot=[hplot;plot(xx01',yy01','k','linestyle','-')];
	
	hplot=[hplot;plot([0 xf],[0.5*yf 0.5*yf],'color','b','linewidth',2)];
	
	plot(XX,YY,'--k')
	plot(XX',YY','--k')
	%plot(XXh',YYh','--g')
	%plot(XXh,YYh,'--g')
	
	str_title=strcat("step:   ",num2str(i));
	[frac vx vy d normal kappa iforce gp]=loadFiles(i);
	draw_interface(x,y,xh,yh,dx,dy,max(size(frac)),frac,normal,d,'m','m','m');
	%axis([0.3 0.8 0.5 1.0])
	%draw_frac(frac,nx*ny,dx,dy,0)
	return
	i=20000;
	str_title=strcat("step:   ",num2str(i));
	[frac vx vy d normal kappa iforce gp]=loadFiles(i);
	draw_interface(x,y,xh,yh,dx,dy,max(size(frac)),frac,normal,d,'r','r','r');
	
end

function[ ]=draw_interface(X,Y,Xh,Yh,dx,dy,n,Frac,Normal,distance,colorp,colorn,colora)
	for i=1:n
		nx=Normal(i,1);ny=Normal(i,2);
		center=[Frac(i,1),Frac(i,2)];
		dc=distance(i);
		drawCell(nx,ny,dx,dy,center,dc,Frac(i,3),colorp,colorn,colora);
	end
end

function[]=drawCell(nx,ny,dx,dy,center,dc,frac,colorp,colorn,colora)
	global hplot;
	if(frac==0||frac==1)
		return;
	end
	plot(center(1),center(2),'ob')
	if(nx*ny>=0)
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
	end
end

function[frac vx vy d normal kappa iforce gp]=loadFiles(step)
	str=strcat("../DATA/",num2str(step),"_frac_visc.txt");
	frac=load(str,"-ascii");
	
	str=strcat("../DATA/",num2str(step),"_vx.txt");
	vx=1;%vx=load(str,"-ascii");
	
	str=strcat("../DATA/",num2str(step),"_vy.txt");
	vy=1;%load(str,"-ascii");
	
	str=strcat("../DATA/",num2str(step),"_d.txt");
	d=load(str,"-ascii");
	
	str=strcat("../DATA/",num2str(step),"_norm.txt");
	normal=load(str,"-ascii");
	
	str=strcat("../DATA/",num2str(step),"_curv.txt");
	kappa=1;%load(str,"-ascii");

	str=strcat("../DATA/",num2str(step),"_if.txt");
	iforce=1;%iforce=load(str,"-ascii");
	
	str=strcat("../DATA/",num2str(step),"_gp.txt");
	gp=1;%gp=load(str,"-ascii");
end

function[]=drawCirclezalesak(xc,yc,r)
  t=[0:0.01:2*pi];
  x=xc+r*cos(t);
  y=yc+r*sin(t);
  plot(x,y,'color','c','linewidth',2);
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
	
	normif=sqrt(sify.**2  + sifx.**2);
	
	normgp=sqrt(gy.**2  + gx.**2);
	
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

function[] = draw_ref_inter()
ref = load('a20.txt');
x = ref(:,1);
y = ref(:,2);
plot(x,y,'b--', 'linewidth', 2)
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
