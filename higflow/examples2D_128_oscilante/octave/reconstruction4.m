source('plic.m');
source('adap3.m');
source('finite_difference_normal_curvature.m');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = build_distance(Nx,Ny,dx,dy,Xh,Yh)
	global Frac Normalx Normaly distance;
	for I=1:Nx
		for J=1:Ny
			if(Frac(I,J)==1||Frac(I,J)==0||(Normalx(I,J)==0&&Normaly(I,J)==0))
				continue;
			end
			nx=Normalx(I,J);ny=Normaly(I,J);
			AREA=Frac(I,J)*dx*dy;
      if(nx*ny<0)
        printf('x=%f y=%f\n',Xh(I),Yh(J));
      end
			distance(I,J)=distance_from_center(nx,ny,dx,dy,AREA);
		end
	end
end


function[ ]=calculation_normal_curvature(Xh,Yh,delta,X,Y)
	global Frac Kappa Normalx Normaly NormalAux;
	
	H=zeros(1,3);V=zeros(1,3);
	Nx=size(Frac,1);
	Ny=size(Frac,2);
	dx=delta(1);dy=delta(2);
	for I=1:Nx
		for J=1:Ny
			if(Frac(I,J)==0||Frac(I,J)==1)
				continue;
			end
			[H_middle auxh origh]=horizontal_row...
			(I,J,Xh,Yh,'MIDDLE');
			[H_top auxh_t origh_t]=horizontal_row...
			(I,J+1,Xh,Yh,'TOP');
			[H_bottom auxh_b origh_b]=horizontal_row...
			(I,J-1,Xh,Yh,'BOTTOM');
			if(abs(auxh+auxh_t+auxh_b)==3)
				H=set_common_orig_horizontal([H_middle,...
				H_top,H_bottom],...
				auxh,[origh,origh_t,origh_b],J,Yh,dx);
				H_middle=H(1);H_top=H(2);H_bottom=H(3);
				calculation_normal_cell_central_2nd_order_finite_difference_Horizontal...
				(H_bottom,H_middle,H_top,dx,dy,I,J,auxh)
				calculation_curvature_cell_central_2nd_order_finite_difference_Horizontal...
				(H_bottom,H_middle,H_top,dx,dy,I,J,auxh)
				NormalAux(I,J)=2;%Naux
				continue;
			else
				[V_middle auxv origv]=vertical_collumn...
				(I,J,Xh,Yh,'MIDDLE');
				[V_right auxv_r origv_r]=vertical_collumn...
				(I+1,J,Xh,Yh,'RIGHT');
				[V_left auxv_l origv_l]=vertical_collumn...
				(I-1,J,Xh,Yh,'LEFT');
				if(abs(auxv+auxv_l+auxv_r)==3)
					V = set_common_orig_vertical...
					([V_middle,V_right,V_left],...
					auxv,[origv,origv_r,origv_l],I,Xh,dy);
					V_middle=V(1);V_right=V(2);V_left=V(3);
					calculation_normal_cell_central_2nd_order_finite_difference_Vertical...
					(V_left,V_middle,V_right,dx,dy,I,J,auxv)
					calculation_curvature_cell_central_2nd_order_finite_difference_Vertical...
					(V_left,V_middle,V_right,dx,dy,I,J,auxv)
					NormalAux(I,J)=1;%Naux
					continue;
				end
			end
			if(abs(auxh+auxh_t)==2)
				[H_top_top auxh_t_t origh_t_t]=...
				horizontal_row(I,J+2,Xh,Yh,'g');
				if(abs(auxh+auxh_t+auxh_t_t)==3)
					H=set_common_orig_horizontal([H_middle,...
					H_top,H_top_top],...
					auxh,[origh,origh_t,origh_t_t],J,Yh,dx);
					H_middle=H(1);H_top=H(2);H_top_top=H(3);
					calculation_normal_cell_progressive_2nd_order_finite_difference_Horizontal...
					(H_middle,H_top,H_top_top,dx,dy,I,J,auxh)
					calculation_curvature_cell_progressive_1st_order_finite_difference_Horizontal...
					(H_middle,H_top,H_top_top,dx,dy,I,J,auxh)
					NormalAux(I,J)=2;%Naux
					continue;
				end
			end
			if(abs(auxh+auxh_b)==2)
				[H_bottom_bottom auxh_b_b origh_b_b]=...
				horizontal_row(I,J-2,Xh,Yh,'g');
				if(abs(auxh+auxh_b+auxh_b_b)==3)
					H=set_common_orig_horizontal([H_middle,...
					H_bottom,H_bottom_bottom],...
					auxh,[origh,origh_b,origh_b_b],J,Yh,dx);
					H_middle=H(1);H_bottom=H(2);H_bottom_bottom=H(3);
					calculation_normal_cell_regressive_2nd_order_finite_difference_Horizontal...
					(H_bottom_bottom,H_bottom,H_middle,dx,dy,I,J,auxh)
					calculation_curvature_cell_regressive_1st_order_finite_difference_Horizontal...
					(H_bottom_bottom,H_bottom,H_middle,dx,dy,I,J,auxh)
					NormalAux(I,J)=2;%Naux
					continue;
				end
			end
			if(abs(auxv+auxv_r)==2)
				[V_right_right auxv_r_r origv_r_r]=...
				vertical_collumn(I+2,J,Xh,Yh,'g');
				if(abs(auxv+auxv_r+auxv_r_r)==3)
					V = set_common_orig_vertical([V_middle,...
					V_right,V_right_right],...
					auxv,[origv,origv_r,origv_r_r],I,Xh,dy);
					V_middle=V(1);V_right=V(2);V_right_right=V(3);
					calculation_normal_cell_progressive_2nd_order_finite_difference_Vertical...
					(V_middle,V_right,V_right_right,dx,dy,I,J,auxv)
					calculation_curvature_cell_progressive_1st_order_finite_difference_Vertical...
					(V_middle,V_right,V_right_right,dx,dy,I,J,auxv)
					NormalAux(I,J)=1;%Naux
					continue;
				end
			end
			if(abs(auxv+auxv_l)==2)
				[V_left_left auxv_l_l origv_l_l]=...
				vertical_collumn(I-2,J,Xh,Yh,'g');
				if(abs(auxv+auxv_l+auxv_l_l)==3)
					V = set_common_orig_vertical([V_middle,...
					V_left,V_left_left],...
					auxv,[origv,origv_l,origv_l_l],I,Xh,dy);
					V_middle=V(1);V_left=V(2);V_left_left=V(3);
					calculation_normal_cell_regressive_2nd_order_finite_difference_Vertical...
					(V_left_left,V_left,V_middle,dx,dy,I,J,auxv)
					calculation_curvature_cell_regressive_1st_order_finite_difference_Vertical...
					(V_left_left,V_left,V_middle,dx,dy,I,J,auxv)
					NormalAux(I,J)=1;%Naux
					continue;
				end
			end
			
			
			
			if(abs(auxh+auxh_t)==2)
			H=set_common_orig_horizontal([H_middle,H_top],...
			auxh,[origh,origh_t],J,Yh,dx);
			H_middle=H(1);H_top=H(2);
			Hy=-auxh*(H_top-H_middle)*dx/(dy);
			normal=-auxh*[-1,Hy];normal=normal/norm(normal);
			Normalx(I,J)=normal(1);Normaly(I,J)=normal(2);
			NormalAux(I,J)=2;%Naux
			disp('A1');fflush(stdout);
		elseif(abs(auxh+auxh_b)==2)
			H=set_common_orig_horizontal([H_middle,H_bottom],...
			auxh,[origh,origh_b],J,Yh,dx);
			H_middle=H(1);H_bottom=H(2);
			Hy=-auxh*(H_middle-H_bottom)*dx/(dy);
			normal=-auxh*[-1,Hy];normal=normal/norm(normal);
			Normalx(I,J)=normal(1);Normaly(I,J)=normal(2);
			NormalAux(I,J)=2;%Naux
			disp('A2');fflush(stdout);
		elseif(abs(auxv+auxv_r)==2)
			V=set_common_orig_vertical([V_middle,V_right],...
			auxv,[origv,origv_r],I,Xh,dy);
			V_middle=V(1);V_right=V(2);
			Vx=-auxv*(V_right-V_middle)*dy/(dx);
			normal=auxv*[-Vx,1];normal=normal/norm(normal);
			Normalx(I,J)=normal(1);Normaly(I,J)=normal(2);
			NormalAux(I,J)=1;%Naux
			disp('A3');fflush(stdout);
		elseif(abs(auxv+auxv_l)==2)
			V=set_common_orig_vertical([V_middle,V_left],...
			auxv,[origv,origv_l],I,Xh,dy);
			V_middle=V(1);V_left=V(2);
			Vx=-auxv*(V_middle-V_left)*dy/(dx);
			normal=auxv*[-Vx,1];normal=normal/norm(normal);
			Normalx(I,J)=normal(1);Normaly(I,J)=normal(2);
			NormalAux(I,J)=1;%Naux
			disp('A4');fflush(stdout);
		end
			
		end
	end
	
end


