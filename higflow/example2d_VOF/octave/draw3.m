function[]=draw3()
source('draws.m');
global Frac;
figure(1);clf;hold on;
[X,Y,Xh,Yh,t]=create_domain_boomerang();
drawDomain(X,Y);
drawFrac(Xh,Yh);
end
