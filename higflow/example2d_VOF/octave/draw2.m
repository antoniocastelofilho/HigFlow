function[]=draw2()
	
	figure(1);clf;hold on;
	
	frac=load('data/0.000000_frac.txt','-ascii');
	fracaux=load('data/0.000000_fracaux.txt','-ascii');
	
	Nx=80;Ny=40;
	
	X=[0:2/Nx:2];dx=2/Nx;
	Xh=[0.5*dx:dx:2-0.5*dx];
	
	Y=[0:1/Ny:1];dy=1/Ny;
	Yh=[0.5*dy:dy:1-0.5*dy];
	
	draw_domain(X,Y);
	
	drawFrac(frac,dx,dy)
	%drawFracAux(fracaux,dx,dy)
	
	return;
	figure(2);clf;hold on;
	
	draw_domain(X,Y);
	
	frac_resh=reshape(frac(:,1),80,40);
	
	drawFrac_resh(frac_resh,dx,dy,Xh,Yh)
	
end

function[]=drawFrac_resh(frac,dx,dy,Xh,Yh)
	
	Nx=size(frac,1);
	Ny=size(frac,2);
	
	for I=1:Nx
		for J=1:Ny
			tc=num2str(frac(I,J));
			text(Xh(I),Yh(J),tc);
			pause();
			
		end
	end
	
end



function[]=drawFrac(frac,dx,dy)
	
	
	N=max(size(frac));
	
	for(i=1:N)
		if(frac(i,3)==0||frac(i,3)==1)
			%continue;
		end
		tc=num2str(frac(i,3));
		text(frac(i,1),frac(i,2),tc);
		pause();
		%plot(frac(i,1),frac(i,2),'.r');
		
	end
	
	
	
end


function[]=drawFracAux(frac,dx,dy)
	
	
	N=max(size(frac));
	
	for(i=1:N)
		if(frac(i,3)==0)%||frac(i,3)==1)
			continue;
		end
		tc=num2str(frac(i,3));
		text(frac(i,1),frac(i,2)-0.25*dy,tc,'color','r');
		%plot(frac(i,1),frac(i,2),'.r');
		
	end
	
	
end






function[]=draw_domain(X,Y)
	
	
	[XX YY]=meshgrid(X,Y);
	
	plot(XX,YY,'k');
	plot(XX',YY','k');
	
	return;
	xc=0.5;yc=0.5;t=[0:0.01:2*pi];r=0.25;
	x=xc+r*cos(t);
	y=yc+r*sin(t);
	plot(x,y,'b');
	
end


function []=draw_norm(vector)
  
  %%figure();
  %%hold on;
  
  draw_domain();
  n = max(size(vector));
  vector_new = vector;
  h=2/80;
  
  for i = 1:n
    if (norm(vector_new(i,3:4)) ~= 0.0)
      vector_new(i,3:4) = h*vector_new(i,3:4)/norm(vector_new(i,3:4)); 
      x1 = vector_new(i,1);
      y1 = vector_new(i,2);
      x2 = x1 + vector_new(i,3);
      y2 = y1 + vector_new(i,4);
    else
      x1 = vector_new(i,1);
      y1 = vector_new(i,2);
      x2 = x1 + vector_new(i,3);
      y2 = y1 + vector_new(i,4);
    end
    plot([x1 x2],[y1 y2]);
    plot(x2,y2,'marker','*'); 
  end
  
  
  
  
end





