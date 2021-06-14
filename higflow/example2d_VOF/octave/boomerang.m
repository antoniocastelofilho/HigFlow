function[]=boomerang()
	
	source('draws.m');
	source('reconstruction4.m');
  
	
	global Frac Normalx Normaly Kappa distance;
	
	[X,Y,Xh,Yh,t]=create_domain_boomerang();
	delta=[X(2)-X(1),Y(2)-Y(1)];
	[Nx Ny]=size(Frac);
	
	figure(1);clf;hold on;
	drawDomain(X,Y);
	t=0;
	%drawCurve(0.5,0.5,0.25,t,min(delta));
	
	calculation_normal_curvature(Xh,Yh,delta,X,Y);
	build_distance(Nx,Ny,delta(1),delta(2),Xh,Yh);
	
	draw_interface(X,Y,Xh,Yh,delta(1),delta(2),Nx,Ny);
	%drawKappa(X,Y,Xh,Yh);
	drawdistance(Xh,Yh)
	return;
	FracAUX=Frac;
	
	
	figure(3);clf;hold on;
	Frac=load('FRAC4_512x256.txt','-ascii');
	
	calculation_normal_curvature(Xh,Yh,delta,X,Y);
	build_distance(Nx,Ny,delta(1),delta(2),Xh,Yh);
	
	draw_interface(X,Y,Xh,Yh,delta(1),delta(2),Nx,Ny);
	drawKappa(X,Y,Xh,Yh);
	drawDomain(X,Y);
	
	
	Frac=Frac-FracAUX;
	drawFrac(Xh,Yh);
	drawCurve(0.5,0.5,0.25,t,min(delta));
	
end

function[]=drawCurve(c1,c2,r,t,dx)
	
	s=[0:4.0*dx:2*pi];
	
	y=r*sin(s)+c2;
	x=r*cos(s) +c1 -4.0*y.*(y - 1.0).*t;
	
	plot(x,y,'color','b','marker','.');
	
	ns=max(size(s));
	
	
	for i=1:ns
		value=kappa(r,s(i),t,c1,c2);
		str=num2str(value);
		text(x(i)+0.05*dx,y(i),str,'color','g');
	end
	
	
end


function[value]=kappa(r,s,t,c1,c2)
	
	value=(-r*(-r*t*(4*c2 + 4*r*sin(s))*cos(s) - ...
	4*r*t*(c2 + r*sin(s) - 1)*cos(s) - ...
	r*sin(s))*sin(s) - r*(-8*r**2*t*cos(s)**2 + ...
	r*t*(4*c2 + 4*r*sin(s))*sin(s) + ...
	4*r*t*(c2 + r*sin(s) - 1)*sin(s) - ...
	r*cos(s))*cos(s))/(r**2*cos(s)**2 + ...
	(-r*t*(4*c2 + 4*r*sin(s))*cos(s) - ...
	4*r*t*(c2 + r*sin(s) - 1)*cos(s) -...
	r*sin(s))**2)**(3/2);
	
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=drawdistance(Xh,Yh)
	global distance Frac;
	
	Nx=size(distance,1);
	Ny=size(distance,2);
	dx=Xh(2)-Xh(1);
	dy=Yh(2)-Yh(1);
	
	for i=1:Nx
		for j=1:Ny
			if(Frac(i,j)>0.0&&Frac(i,j)<1.0)
				text(Xh(i),Yh(j),num2str(distance(i,j)));
			end
			
		end
	end
end