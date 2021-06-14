1;

%%%%%%%%%%%% 3-points CENTRAL HORIZONTAL %%%%%%%%%%%%

function[ ]=calculation_curvature_cell_central_2nd_order_finite_difference_Horizontal(Hk_,H,Hk,dx,dy,I,J,auxh)
	
	global Kappa
	
	Hy=-auxh*(-Hk_+Hk)*dx/(2*dy);
	Hyy=-auxh*(Hk_-2*H+Hk)*dx/(dy^2);
	Kappa(I,J)=auxh*Hyy/((1+Hy^2)^1.5);
	
end

function[ ]=calculation_normal_cell_central_2nd_order_finite_difference_Horizontal(Hk_,H,Hk,dx,dy,I,J,auxh)
	
	global Normalx Normaly
	
	Hy=-auxh*(-Hk_+Hk)*dx/(2*dy);
	normal=-auxh*[-1,Hy];normal=normal/norm(normal);
	Normalx(I,J)=normal(1);Normaly(I,J)=normal(2);
	
end

%%%%%%%%%%%% 3-points CENTRAL VERTICAL   %%%%%%%%%%%%

function[ ]=calculation_curvature_cell_central_2nd_order_finite_difference_Vertical(Vk_,V,Vk,dx,dy,I,J,auxv)
	
	global Kappa
	
	Vx=-auxv*(-Vk_+Vk)*dy/(2*dx);
	Vxx=-auxv*(Vk_-2*V+Vk)*dy/(dx^2);
	Kappa(I,J)=auxv*Vxx/((1+Vx^2)^1.5);
	
end


function[ ]=calculation_normal_cell_central_2nd_order_finite_difference_Vertical(Vk_,V,Vk,dx,dy,I,J,auxv)
	
	global Normalx Normaly
	
	Vx=-auxv*(-Vk_+Vk)*dy/(2*dx);
	normal=auxv*[-Vx,1];normal=normal/norm(normal);
	Normalx(I,J)=normal(1);Normaly(I,J)=normal(2);
	
end

%%%%%%%%%% 3-points PROGRESIVE HORIZONTAL %%%%%%%%%

function[ ]=calculation_curvature_cell_progressive_1st_order_finite_difference_Horizontal(H,Hk,Hkk,dx,dy,I,J,auxh)
	
	global Kappa
	
	Hy=-auxh*(-1.5*H+2*Hk-0.5*Hkk)*dx/(dy);
	Hyy=-auxh*(H-2*Hk+Hkk)*dx/(dy^2);
	Kappa(I,J)=auxh*Hyy/((1+Hy^2)^1.5);
	
end

function[ ]=calculation_normal_cell_progressive_2nd_order_finite_difference_Horizontal(H,Hk,Hkk,dx,dy,I,J,auxh)
	
	global Normalx Normaly
	
	Hy=-auxh*(-1.5*H+2*Hk-0.5*Hkk)*dx/(dy);
	normal=-auxh*[-1,Hy];normal=normal/norm(normal);
	Normalx(I,J)=normal(1);Normaly(I,J)=normal(2);
	
end

%%%%%%%%%%% 3-points PROGRESSIVE VERTICAL %%%%%%%%%%

function[ ]=calculation_curvature_cell_progressive_1st_order_finite_difference_Vertical(V,Vk,Vkk,dx,dy,I,J,auxv)
	
	global Kappa
	
	Vx=-auxv*(-1.5*V+2*Vk-0.5*Vkk)*dy/(dx);
	Vxx=-auxv*(V-2*Vk+Vkk)*dy/(dx^2);
	Kappa(I,J)=auxv*Vxx/((1+Vx^2)^1.5);
	
end

function[ ]=calculation_normal_cell_progressive_2nd_order_finite_difference_Vertical(V,Vk,Vkk,dx,dy,I,J,auxv)
	
	global Normalx Normaly
	
	Vx=-auxv*(-1.5*V+2*Vk-0.5*Vkk)*dy/(dx);
	normal=auxv*[-Vx,1];normal=normal/norm(normal);
	Normalx(I,J)=normal(1);Normaly(I,J)=normal(2);
	
end

%%%%%%%%%% 3-points REGRESSIVE HORIZONTAL %%%%%%%%%

function[ ]=calculation_curvature_cell_regressive_1st_order_finite_difference_Horizontal(Hk_k_,Hk_,H,dx,dy,I,J,auxh)
	
	global Kappa
	
	Hy=-auxh*(0.5*Hk_k_-2*Hk_+1.5*H)*dx/(dy);
	Hyy=-auxh*(Hk_k_- 2*Hk_ + H)*dx/(dy^2);
	Kappa(I,J)=auxh*Hyy/((1+Hy^2)^1.5);
	
end

function[ ]=calculation_normal_cell_regressive_2nd_order_finite_difference_Horizontal(Hk_k_,Hk_,H,dx,dy,I,J,auxh)
	
	global Normalx Normaly
	
	Hy=-auxh*(0.5*Hk_k_-2*Hk_+1.5*H)*dx/(dy);
	normal=-auxh*[-1,Hy];normal=normal/norm(normal);
	Normalx(I,J)=normal(1);Normaly(I,J)=normal(2);
	
end


%%%%%%%%%% 3-points REGRESSIVE VERTICAL %%%%%%%%%

function[ ]=calculation_curvature_cell_regressive_1st_order_finite_difference_Vertical(Vk_k_,Vk_,V,dx,dy,I,J,auxv)
	
	global Kappa
	
	Vx=-auxv*(0.5*Vk_k_-2*Vk_+1.5*V)*dy/(dx);
	Vxx=-auxv*(Vk_k_- 2*Vk_ + V)*dy/(dx^2);
	Kappa(I,J)=auxv*Vxx/((1+Vx^2)^1.5);
	
end

function[ ]=calculation_normal_cell_regressive_2nd_order_finite_difference_Vertical(Vk_k_,Vk_,V,dx,dy,I,J,auxv)
	
	global Normalx Normaly
	
	Vx=-auxv*(0.5*Vk_k_-2*Vk_+1.5*V)*dy/(dx);
	normal=auxv*[-Vx,1];normal=normal/norm(normal);
	Normalx(I,J)=normal(1);Normaly(I,J)=normal(2);
	
end

