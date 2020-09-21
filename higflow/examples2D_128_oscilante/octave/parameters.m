function[]=parameters()
	
	rho1=787.88;
	rho2=1.1768;
	R=22.56e-3;
	n=2;
	sigma=0.02361
	omega=sqrt((n^3-n)*sigma/(rho1+rho2)/R^3)
	T=2*pi/omega
	printf('T=%.18f\n',T);
	
	
end