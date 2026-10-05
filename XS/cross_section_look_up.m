% Develop a look-up table for 3D PIC code (Deca et al. 2018)
% of the electron-impact ionisation on H2O
% -----------------------------------------------------------

close all

%read_comet
% Energy grid: Einc (eV)
% Total ionisation cross section of H2O due to e- impact in cm^2
sigma = sige_ioni(1,:);

figure
f=plot(Einc,sige_ioni(1,:));
hold on
grid on
axis([0 300 0 2.5e-16])

% Set an energy grid

Eth = 12.6;

% E1 = [13:10:60];
% E2 = [50:10:100];
% E3 = [90:10:300];
% 
% 
% [m ind_min] = min(abs(Einc-min(E1)));
% [m ind_max] = min(abs(Einc-max(E1)));
% p1 = polyfit(Einc(ind_min:ind_max),sigma(ind_min:ind_max),1);
% 
% sigma_fit = polyval(p1,E1);
% f1=plot(E1,sigma_fit);
% 
% [m ind_min] = min(abs(Einc-min(E2)));
% [m ind_max] = min(abs(Einc-max(E2)));
% p2 = polyfit(Einc(ind_min:ind_max),sigma(ind_min:ind_max),1);
% 
% sigma_fit = polyval(p2,E2);
% f2=plot(E2,sigma_fit);
% 
% [m ind_min] = min(abs(Einc-min(E3)));
% [m ind_max] = min(abs(Einc-max(E3)));
% p3 = polyfit(Einc(ind_min:ind_max),sigma(ind_min:ind_max),1);
% 
% sigma_fit = polyval(p3,E3);
% f3=plot(E3,sigma_fit);

E_look = [13:1:50, 52:2:100, 120:20:200];

for ie=1:length(E_look)-1
    [m,ind_min]=min(abs(Einc-E_look(ie)));
    [m,ind_max]=min(abs(Einc-E_look(ie+1)));
    sigma_look(ie) = integrale([Einc(ind_min:ind_max)],...
        [sigma(ind_min:ind_max)])/(Einc(ind_max)-Einc(ind_min));
    f0=plot([E_look(ie),E_look(ie+1)],[sigma_look(ie) sigma_look(ie)],'r');
    a_var(ie,:) = [E_look(ie),E_look(ie+1),sigma_look(ie)];
end

save -ascii e_H2O_ioni_cross_section.txt  a_var

