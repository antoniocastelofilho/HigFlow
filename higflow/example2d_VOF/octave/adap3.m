1.0;
function[value]=is_inside(I,J)
	
	 global Frac;
	
	 [Nx Ny]=size(Frac);
	
	 value=true;
	
	 if(I>Nx||I<1||J>Ny||J<1)	
		 value=false;
     end
end

function[vertical,aux,orig]=vertical_collumn(I,...
J,Xh,Yh,info)
	
	vertical=0;
	aux=0;
	orig=0;
	
	if(is_inside(I,J)==false)
		return;
	end
	
	 dx=Xh(2)-Xh(1);
	 dy=Yh(2)-Yh(1);
	
	 global Frac;
	 j=1;
	
	 vertical=Frac(I,J);
	 while(is_inside(I,J+j)==true&&Frac(I,J+j)>0&&Frac(I,J+j)<1)
		 vertical=vertical+Frac(I,J+j); 
		 j=j+1;
	 end
	 
	 if(is_inside(I,J+j)==false)
		 disp('*********top boundary achieved***********');
		 vertical=0;
		 return;
	 end
	
	 vertical=vertical+Frac(I,J+j);
	 vt=Frac(I,J+j);
	 jt=J+j;
	
	 j=-1;
	 
	 %bottom
	 while(is_inside(I,J+j)==true&&Frac(I,J+j)>0&&Frac(I,J+j)<1&&Frac(I,J+j)~=vt)
		vertical=vertical+Frac(I,J+j);
		j=j-1;
	end
	
	
	if(is_inside(I,J+j)==false||Frac(I,J+j)==vt)
		vertical=0;
		
		
		%disp('**********phases are not distinct or boundary achieved***********');
		%printf('is_inside: %d vb=vt=%d\n',is_inside(I,J+j),vt);
		%printf('info: %s Frac= %f x=%f y=%f\n',info,Frac(I,J),Xh(I),Yh(J));
		
		if(strcmp(info,'MIDDLE')==1)
			fill_cell(Xh(I),Yh(J),dx,dy,Frac(I,J),'c','M');
		elseif(strcmp(info,'RIGHT')==1)
			fill_cell(Xh(I-1),Yh(J),dx,dy,Frac(I-1,J),'c','R');
		elseif(strcmp(info,'LEFT')==1)
			fill_cell(Xh(I+1),Yh(J),dx,dy,Frac(I+1,J),'c','L');
		end
		
		return;
		
	end
	
	jb=J+j;
	vertical=vertical+Frac(I,J+j);
	aux=1;
	orig=Yh(jt)+dy/2;
	
	 if(vt==0)
		 aux=-1;
		 orig=Yh(jb)-dy/2;
	 end
end

function[]=fill_cell(x,y,dx,dy,frac,color,info);
	
	v1=[x-0.5*dx,y-0.5*dy];
	v2=[x+0.5*dx,y-0.5*dy];
	v3=[x+0.5*dx,y+0.5*dy];
	v4=[x-0.5*dx,y+0.5*dy];
	v=[v1;v2;v3;v4];
	fill(v(:,1),v(:,2),color);
	%text(x-0.5*dx,y-0.5*dy,num2str(frac));
	text(x,y+0.25*dy,info);
	
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[horizontal,aux,orig]=horizontal_row(I,J,Xh,Yh,color)
	
	aux=0;
	horizontal=0;
	orig=0;
	
	if(is_inside(I,J)==false)
		return;
	end
	
	dx=Xh(2)-Xh(1);
	dy=Yh(2)-Yh(1);
	
	global Frac;
	
	i=1;
	
	%%%right
	horizontal=Frac(I,J);
	while(is_inside(I+i,J)==true&&Frac(I+i,J)>0&&Frac(I+i,J)<1)
		horizontal=horizontal+Frac(I+i,J);
		i=i+1;
	end
	
	%left
	if(is_inside(I+i,J)==false)
		%disp('horizontal boundary');
		horizontal=0;
		return;
	end
	
	horizontal=horizontal+Frac(I+i,J);
	
	hr=Frac(I+i,J);
	ir=I+i;
	i=-1;
	
	while(is_inside(I+i,J)==true&&Frac(I+i,J)>0&&Frac(I+i,J)<1&&Frac(I+i,J)!=hr)
		horizontal=horizontal+Frac(I+i,J); 
		i=i-1;
	end
	
	if(is_inside(I+i,J)==false||Frac(I+i,J)==hr)
		horizontal=0;
		%disp('phases are not distinct or boundary achieved');%
		%printf('is_inside: %d hl=hr=%d\n',is_inside(I+i,J),hr);
		return;
	end
	
	horizontal=horizontal+Frac(I+i,J);
	il=I+i;
	
	aux=1;
	orig=Xh(ir)+dx/2;
	
	if(hr==0)
		aux=-1;
		orig=Xh(il)-dx/2;
	end
end

function[Vnew]=...
set_common_orig_vertical(V,...
auxv,origv,I,Xh,dy)
	
	Vnew=V;
	n=max(size(V));
	
	if(auxv>0)
		
		orig=max(origv);
		delta=orig*ones(1,n)-origv;
		delta=delta/dy;
		
		Vnew=Vnew+delta;
		
	else
		
		orig=min(origv);
		delta=-orig*ones(1,n)+origv;
		delta=delta/dy;
		
		Vnew=Vnew+delta;
		
	end
	
end

function[Hnew]=...
set_common_orig_horizontal(H,...
auxh,origh,J,Yh,dx)
	
	Hnew=H;
	n=max(size(H));
	
	if(auxh>0)
		
		orig=max(origh);
		delta=orig*ones(1,n)-origh;
		delta=delta/dx;
		Hnew=Hnew+delta;
		
	else
		
		orig=min(origh);
		delta=-orig*ones(1,n)+origh;
		delta=delta/dx;
		Hnew=Hnew+delta;
		
	end
	
end










