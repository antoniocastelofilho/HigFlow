function[]=create_domain_ju()
	
	global delta
	global Nx Ny
	
	Nx=128;Ny=128;
	xf=75e-3;yf=75e-3;
	
	delta=[xf/Nx yf/Ny];
	dt = 0.00005;
	
	X=[0:delta(1):xf];Y=[0:delta(2):yf];
	Xh=[0.5*delta(1):delta(1):xf-0.5*delta(1)];
	Yh=[0.5*delta(2):delta(2):yf-0.5*delta(2)];
	
	Nxh=max(size(Xh));
	Nyh=max(size(Yh));
	
	cont = 0;
	for step = 0:100:8700
		cont = cont + 1;
		%disp('fraction initialization...');
		[Frac]=loadFiles(step);
		
		[rxnum] = findRadius(Xh,Yh,Frac);
		
		rxnumvec(cont) = rxnum;
		T(cont)        = step*dt;
		
	end
	figure(1)
	figure(1, 'DefaultAxesFontSize',14)
	plot(T,rxnumvec,'b');
	xlabel('t[s]', 'Color', 'k')
	ylabel('y[m]', 'Color', 'k')
	legend('64x64','location', 'northwest')
end

function[frac]=loadFiles(step)
	str=strcat('data/',num2str(step),'_frac.txt');
	frac=load(str,'-ascii');
end

function[r]=findRadius(Xh,Yh,Frac)
	
	global Nx Ny
	
	r=0;
	dx=Xh(2)-Xh(1);dy=Yh(2)-Yh(1);
	
	j = ceil(0.5*length(Yh));
	
	for i=1:Nx
		k = (i-1)*Ny + j;
		if(Frac(k,3) > 0.0)
			r=r+Frac(k,3);
		end
	end
	rx=0.5*r*dx;
	
	i= ceil(0.5*length(Xh));
	
	r=0;
	for j=1:Ny
		k = (i-1)*Ny + j;
		if(Frac(k,3) > 0.0)
			r=r+Frac(k,3);
		end
	end
	
	ry=0.5*r*dy;
	r=37.5e-3+ry;
	
	
	
end