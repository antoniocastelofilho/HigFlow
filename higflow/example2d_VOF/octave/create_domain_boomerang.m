function[X,Y,Xh,Yh,t]=create_domain_boomerang()
	
	global Frac Kappa Normalx Normaly NormalAux;
	global distance;
	global TIME
	
	
	Nx=80;Ny=40;
	xf=2;yf=1;
	
	delta=[xf/Nx yf/Ny];
	
	X=[0:delta(1):xf];Y=[0:delta(2):yf];
	Xh=[0.5*delta(1):delta(1):xf-0.5*delta(1)];
	Yh=[0.5*delta(2):delta(2):yf-0.5*delta(2)];
	
	dt=0.2*delta(1);
	t=1280*dt;
	TIME=t;
	
	TIME=0.0;
	
	Nxh=max(size(Xh))
	Nyh=max(size(Yh));
	
	Frac=zeros(Nxh,Nyh);
	Kappa=Frac;Normalx=Frac;
	Normaly=Frac;NormalAux=Frac;
	distance = Frac;
	disp('fraction initialization...');fflush(stdout);tic
	initilization_fraction(Xh,Yh);
	%Frac=load('boomerang512x256.txt','-ascii');
	disp('done');fflush(stdout);toc
	fflush(stdout);
	%save boomerang512x256.txt Frac;
	
	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[ ]=initilization_fraction(Xh,Yh)
	global Frac;
	radius=0.25;r=radius;
	
	Nx=max(size(Xh));Ny=max(size(Yh));
	delta=[Xh(2)-Xh(1),Yh(2)-Yh(1)];
	
	for i=1:Nx
		for j=1:Ny
			center=[Xh(i) Yh(j)];
			value=getFracVol(center,delta);
			Frac(i,j)=value;
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[ value]=square_case(f0,f1,f2,f3)
	
	value=-1;
	
	if((f0 <= 0.0) && (f1 <= 0.0) && (f2 <= 0.0) && (f3 <= 0.0))
		value=0;
	elseif((f0 > 0.0) && (f1 <= 0.0) && (f2 <= 0.0) && (f3 <= 0.0))
		value=1;
	elseif((f0 <= 0.0) && (f1 > 0.0) && (f2 <= 0.0) && (f3 <= 0.0))
		value=2;
	elseif((f0 <= 0.0) && (f1 <= 0.0) && (f2 > 0.0) && (f3 <= 0.0))
		value=3;
	elseif((f0 <= 0.0) && (f1 <= 0.0) && (f2 <= 0.0) && (f3 > 0.0))
		value=4;
	elseif((f0 > 0.0) && (f1 > 0.0) && (f2 <= 0.0) && (f3 <= 0.0))
		value=5;
	elseif((f0 > 0.0) && (f1 <= 0.0) && (f2 > 0.0) && (f3 <= 0.0))
		value=6;
	elseif((f0 <= 0.0) && (f1 > 0.0) && (f2 <= 0.0) && (f3 > 0.0))
		value=7;
	elseif((f0 <= 0.0) && (f1 <= 0.0) && (f2 > 0.0) && (f3 > 0.0))
		value=8;
	elseif((f0 <= 0.0) && (f1 > 0.0) && (f2 > 0.0) && (f3 <= 0.0))
		value=9;
	elseif((f0 > 0.0) && (f1 <= 0.0) && (f2 <= 0.0) && (f3 > 0.0))
		value=10;
	elseif((f0 <= 0.0) && (f1 > 0.0) && (f2 > 0.0) && (f3 > 0.0))
		value=11;
	elseif((f0 > 0.0) && (f1 <= 0.0) && (f2 > 0.0) && (f3 > 0.0))
		value=12;
	elseif((f0 > 0.0) && (f1 > 0.0) && (f2 <= 0.0) && (f3 > 0.0))
		value=13;
	elseif((f0 > 0.0) && (f1 > 0.0) && (f2 > 0.0) && (f3 <= 0.0))
		value=14;
	elseif((f0 > 0.0) && (f1 > 0.0) && (f2 > 0.0) && (f3 > 0.0))
		value=15;
	end
	
end

function[ p]=Intersec(p0,p1,f0,f1)
	
	lambda=-f0/(f1-f0);
	p=p0+lambda*(p1-p0);
	
end

function[value]=func(p)
	
	global TIME;
	
	cCircle=[0.5 0.5];radius=0.25;
	
	c=cCircle;r=radius;t=TIME;
	
	%value=-(p(1)-c(1))^2/(r^2)-(p(2)-c(2))^2/(r^2)+1;
	%y=r*sin(s)+c2;
	%x=r*cos(s) +c1 -4.0*y.*(y - 1.0).*t;
	
	value=-(p(1)-c(1)+4.0*p(2)*(p(2) - 1.0)*t)^2/(r^2)-(p(2)-c(2))^2/(r^2)+1;
	
end
function[value]=getFracVol(center,delta)
	
	value=-1;
	
	%bottom left corner
	p0(1)=center(1)-0.5*delta(1);
	p0(2)=center(2)-0.5*delta(2);
	f0=func(p0);
	
	%bottom rigth corner
	p1(1)=center(1)+0.5*delta(1);
	p1(2)=center(2)-0.5*delta(2);
	f1=func(p1);
	
	%top left corner
	p2(1)=center(1)-0.5*delta(1);
	p2(2)=center(2)+0.5*delta(2);
	f2=func(p2);
	
	%top right corner
	p3(1)=center(1)+0.5*delta(1);
	p3(2)=center(2)+0.5*delta(2);
	f3=func(p3);
	
	case_=square_case(f0,f1,f2,f3);
	
	if(case_==0)
		value=0;
	elseif(case_==15)
		value=delta(1)*delta(2);
	else
		value=0;
		N=8;
		delta_new=delta/N;
		x0=center(1)-delta(1)/2;
		y0=center(2)-delta(2)/2;
		xh=[x0+delta_new(1)/2:delta_new(1):x0+delta(1)-delta_new(1)/2];
		yh=[y0+delta_new(2)/2:delta_new(2):y0+delta(2)-delta_new(2)/2];
		for i=1:N
			for j=1:N
				center_new=[xh(i),yh(j)];
				value=value+i_vN(center_new,delta_new);
			end
		end
	end
	
	value=value/(delta(1)*delta(2));
	
end

function[value]=i_vN(center,delta)
	value=-1;
	
	%bottom left corner
	p0(1)=center(1)-0.5*delta(1);
	p0(2)=center(2)-0.5*delta(2);
	f0=func(p0);
	
	%bottom rigth corner
	p1(1)=center(1)+0.5*delta(1);
	p1(2)=center(2)-0.5*delta(2);
	f1=func(p1);
	
	%top left corner
	p2(1)=center(1)-0.5*delta(1);
	p2(2)=center(2)+0.5*delta(2);
	f2=func(p2);
	
	%top right corner
	p3(1)=center(1)+0.5*delta(1);
	p3(2)=center(2)+0.5*delta(2);
	f3=func(p3);
	
	case_=square_case(f0,f1,f2,f3);
	
	if(case_==0)
		value=0;
	elseif(case_==1)
		p01=Intersec(p0,p1,f0,f1);
		p02=Intersec(p0,p2,f0,f2);
		value = 0.5*(p01(1)-p0(1))*(p02(2)-p0(2));
	elseif(case_==2)
		p01=Intersec(p0,p1,f0,f1);
		p13=Intersec(p1,p3,f1,f3);
		value = 0.5*(p1(1)-p01(1))*(p13(2)-p1(2));
	elseif(case_==3)
		p02=Intersec(p0,p2,f0,f2);
		p23=Intersec(p2,p3,f2,f3);
		value = 0.5*(p23(1)-p2(1))*(p2(2)-p02(2));
	elseif(case_==4)
		p13=Intersec(p1,p3,f1,f3);
		p23=Intersec(p2,p3,f2,f3);
		value = 0.5*(p3(1)-p23(1))*(p3(2)-p13(2));
	elseif(case_==5)
		p02=Intersec(p0,p2,f0,f2);
		p13=Intersec(p1,p3,f1,f3);
		value = 0.5*((p02(2)-p0(2))+(p13(2)-p1(2)))*(p1(1)-p0(1));
	elseif(case_==6)
		p01=Intersec(p0,p1,f0,f1);
		p23=Intersec(p2,p3,f2,f3);
		value = 0.5*((p01(1)-p0(1))+(p23(1)-p2(1)))*(p2(2)-p0(2));
	elseif(case_==7)
		p01=Intersec(p0,p1,f0,f1);
		p23=Intersec(p2,p3,f2,f3);
		value = 0.5*((p1(1)-p01(1))+(p3(1)-p23(1)))*(p3(2)-p1(2));
	elseif(case_==8)
		p02=Intersec(p0,p2,f0,f2);
		p13=Intersec(p1,p3,f1,f3);
		value = 0.5*((p2(2)-p02(2))+(p3(2)-p13(2)))*(p3(1)-p2(1));
	elseif(case_==9)
		disp('----malha pouco redinada------');%exit;
	elseif(case_==10)
		disp('----malha pouco refinada------');%exit;
	elseif(case_==11)
		p01=Intersec(p0,p1,f0,f1);
		p02=Intersec(p0,p2,f0,f2);
		value = delta(1)*delta(2) - 0.5*(p01(1)-p0(1))*(p02(2)-p0(2));
	elseif(case_==12)
		p01=Intersec(p0,p1,f0,f1);
		p13=Intersec(p1,p3,f1,f3);
		value = delta(1)*delta(2) - 0.5*(p1(1)-p01(1))*(p13(2)-p1(2));
	elseif(case_==13)
		p02=Intersec(p0,p2,f0,f2);
		p23=Intersec(p2,p3,f2,f3);
		value = delta(1)*delta(2) - 0.5*(p23(1)-p2(1))*(p2(2)-p02(2));
	elseif(case_==14)
		p13=Intersec(p1,p3,f1,f3);
		p23=Intersec(p2,p3,f2,f3);
		value = delta(1)*delta(2) - 0.5*(p3(1)-p23(1))*(p3(2)-p13(2));
	elseif(case_==15)
		value=delta(1)*delta(2);
	end
	
end
