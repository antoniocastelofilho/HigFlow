/*UNIFORM RECTANGULAR MESH*/

#include <stdio.h>


void create_domain_amr(double x0,double xf,
double y0,double yf,int Nx,int Ny)
{
	double dx=(xf-x0)/Nx;
	double dy=(yf-y0)/Ny;
	
	FILE*fp=fopen("domain/ch-d-0.amr","w");
	
	fprintf(fp,"%.12lf %.12lf %.12lf %.12lf\n",x0,xf,y0
	,yf);
	
	fprintf(fp,"1\n");
	
	fprintf(fp,"%.12lf %.12lf 1\n",dx,dy);
	
	fprintf(fp,"1 1 %d %d",Nx,Ny);
	
	fclose(fp);
}

void create_boundary_amr(double x0,double xf,
double y0,double yf,int Nx,int Ny)
{
	
	double dx=(xf-x0)/Nx;
	double dy=(yf-y0)/Ny;
	
	/*LEFT*/
	
	FILE*fp=fopen("bc/ch-bc-0.amr"
	,"w");
	
	fprintf(fp,"%.12lf %.12lf %.12lf %.12lf\n",x0,x0,y0,yf);
	
	fprintf(fp,"1\n");
	
	fprintf(fp,"0 %.12lf 1\n",dy);
	
	fprintf(fp,"1 1 1 %d",Ny);
	
	fclose(fp);
	
	/*TOP*/
	
	fp=fopen("bc/ch-bc-1.amr"
	,"w");
	
	fprintf(fp,"%.12lf %.12lf %.12lf %.12lf\n",x0,xf,yf,yf);
	
	fprintf(fp,"1\n");
	
	fprintf(fp,"%.12lf 0 1\n",dx);
	
	fprintf(fp,"1 1 %d 1",Nx);
	
	fclose(fp);
	
	/*RIGHT*/
	
	fp=fopen("bc/ch-bc-2.amr"
	,"w");
	
	fprintf(fp,"%.12lf %.12lf %.12lf %.12lf\n",xf,xf,y0,yf);
	
	fprintf(fp,"1\n");
	
	fprintf(fp,"0 %.12lf 1\n",dy);
	
	fprintf(fp,"1 1 1 %d",Ny);
	
	fclose(fp);
	
	/*BOTTOM*/
	
	fp=fopen("bc/ch-bc-3.amr"
	,"w");
	
	fprintf(fp,"%.12lf %.12lf %.12lf %.12lf\n",x0,xf,y0,y0);
	
	fprintf(fp,"1\n");
	
	fprintf(fp,"%.12lf 0 1\n",dx);
	
	fprintf(fp,"1 1 %d 1",Nx);
	
	fclose(fp);
}


int main()
{
	/* DOMAIN*/
	
	double x0,xf,y0,yf;
	
	x0=0.0;xf=16.0;
	y0=0.0;yf=8.0;
	
	/*MESH RESOLUTION*/
	
	int Nx,Ny;
	Nx=64; Ny=32;
	printf("x0 = %.12lf , y0 = %.12lf , xf = %.12lf , yf = %.12lf,  nx= %d , ny = %d \n", x0,y0,xf,yf,Nx,Ny);
	create_domain_amr(x0,xf,y0,yf,Nx,Ny);
	create_boundary_amr(x0,xf,y0,yf,Nx,Ny);
	
	return 0;
}
