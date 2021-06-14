source('plic.m');

function[]=drawDomain(X,Y)
	global Frac;
	[XX YY]=meshgrid(X,Y);
	plot(XX,YY,'color','b','linewidth',0.5);
	plot(XX',YY','color','b','linewidth',0.5);
	drawCircle(X,Y);
	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=drawFrac(Xh,Yh)
	global Frac;
	Nxh=max(size(Xh));Nyh=max(size(Yh));
	dy=Yh(2)-Yh(1);
	for i=1:Nxh
		for j=1:Nyh
			if(Frac(i,j)>0&&Frac(i,j)<1)
				text(Xh(i),Yh(j)-0.25*dy,num2str(Frac(i,j)),'color','m');
			end
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=drawKappa(X,Y,Xh,Yh)
	global Kappa;
	global Normalx Normaly NormalAux;
	global Frac;
	
	Nx=size(Kappa,1);
	Ny=size(Kappa,2);
	dx=Xh(2)-Xh(1);
	dy=Yh(2)-Yh(1);
	
	PP1=[0 0];
	PP2=[0 0];
	
	for i=1:Nx
		for j=1:Ny
			
			if(abs(Kappa(i,j))>0)
				kappa=num2str(Kappa(i,j));
				if(NormalAux(i,j)==1)
					color='b';
				else
					color='k';
				end
				text(Xh(i),Yh(j),kappa,'color',color);
			end
			
			if(Frac(i,j)>0&&Frac(i,j)<1&&Kappa(i,j)==0)
				v1=[Xh(i)-0.5*dx,Yh(j)-0.5*dy];
				v2=[v1(1)+dx,v1(2)];
				v3=[v2(1),v2(2)+dy];
				v4=[v3(1)-dx,v3(2)];
				v=[v1;v2;v3;v4];
				fill(v(:,1),v(:,2),'g');
			end
			
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[ ]=draw_interface(X,Y,Xh,Yh,dx,dy,Nx,Ny)
	
	global Frac;
	global Normalx;
	global Normaly;
	global distance;
	
	for I=1:Nx
		for J=1:Ny
			nx=Normalx(I,J);ny=Normaly(I,J);
			center=[Xh(I),Yh(J)];
			dc=distance(I,J);
			drawCell(nx,ny,dx,dy,center,dc,Frac(I,J),...
			'k','r','b');
			
		end
	end
	
	
end

function[]=draw_x_flux_rectangle(Xh,Yh,t,dt)
	
	global Frac;
	
	dx=Xh(2)-Xh(1);
	dy=Yh(2)-Yh(1);
	
	[Nx Ny]=size(Frac);
	
	for I=1:Nx
		for J=1:Ny
			
			if(Frac(I,J)==1||Frac(I,J)==0)
				continue
			end
			
			x=Xh(I);y=Yh(J);
			
			xr=x;yr=y;ur=u(x+dx/2,y,t);
			v1=[xr+dx/2,yr-dy/2];v2=[v1(1),yr+dy/2 ];
			v3=[v1(1)-ur*dt,yr+dy/2];v4=[v1(1)-ur*dt,yr-dy/2];
			v=[v1;v2;v3;v4];fill(v(:,1),v(:,2),'y');
			
			xl=x-dx;yl=y;ul=u(x-dx/2,y,t);
			v1=[xl+dx/2,yl-dy/2];v2=[v1(1),yl+dy/2 ];
			v3=[v1(1)-ur*dt,yl+dy/2];v4=[v1(1)-ul*dt,yl-dy/2];
			v=[v1;v2;v3;v4];fill(v(:,1),v(:,2),'g');
			
			
		end
	end
	
end


function[]=drawCell(nx,ny,dx,dy,center,dc,frac,...
colorp,colorn,colora)
	
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
		plot(vx,vy,'color',colorp,'linewidth',1.5);
		PH=PM+min(dx,dy)*[nx ny];
		plot([PM(1),PH(1)],[PM(2),PH(2)],colora);
		plot(PH(1),PH(2),'marker','.','color',colora);
		
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
		plot(vx,vy,'color',colorn,'linewidth',1.5);
		PH=PM+min(dx,dy)*[nx ny];
		plot([PM(1),PH(1)],[PM(2),PH(2)],colora);
		plot(PH(1),PH(2),'marker','.','color',colora);
	end
	
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=drawCircle(X,Y)
	cCircle=[0.5 0.5];radius=0.25;
	c=cCircle;r=radius;
	t=[0:0.01:2*pi];
	x=c(1)+r*cos(t);
	y=c(2)+r*sin(t);
	plot(x,y,'marker','*','color','c');
	plot([0 X(end)],[0.5*Y(end) 0.5*Y(end)],'color','m','linewidth',2);
end