function[]=kinetic_energy()
global  Kin
load_initialize_arquichives();
plot_figure_kinetic_energy();
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=plot_figure_kinetic_energy()
    global Kin
    h = figure(1);
    plot(Kin.t,Kin.Ene,'r','LineWidth',1.5)
    axis([0 20 0 0.02])
    legend('Kinetic Energy','Location','southeast');
    title('Parameters - Re = 1.0, \rho_{0}=\rho_{1}=1 and \mu_{0} = \mu_{1}=1')
    xlabel('Time');
    ylabel('Kinetic energy');
    grid minor
    hold off
%    savefig(h,'Kini_Ener.fig')
    %saveas(gcf,'Kini_Ener','epsc')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=load_initialize_arquichives()
global Kin
initialize_kinetic_energy();
load_arquichives_kinetic_energy();
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=initialize_kinetic_energy()
global Kin
Kin = struct('t',{0},'Ene',{0});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=load_arquichives_kinetic_energy()
global Kin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str            = strcat('../DATA/Kinetic.txt');
aux          = load(str,'-ascii');
Kin.t        = aux(:,1);
Kin.Ene = aux(:,2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%