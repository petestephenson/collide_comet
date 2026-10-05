%% Importing Cross Sections
fid=fopen('/Users/pps18/Documents/Cross_Sections/Ionization_XS/comet_H2O_CO_CO2_ioni.txt');

%H2O XSections
formatSpec='%f%f%f%f%f%f%f';
startRow=15;
H2O_XS=textscan(fid, formatSpec, 'Delimiter', ',', 'HeaderLines' ,startRow);
H2O_XS=cell2mat(H2O_XS);

% CO XSections
startRow=16;
formatSpec='%f%f%f%f%f%f%f%f';
CO_XS=textscan(fid, formatSpec, 'Delimiter', ',', 'HeaderLines' ,startRow);
CO_XS=cell2mat(CO_XS);

% CO2 Xsections
startRow=18;
formatSpec='%f%f%f%f%f%f%f%f%f%f%f%f';
CO2_XS=textscan(fid, formatSpec, 'Delimiter', ',', 'HeaderLines' ,startRow);
fclose(fid);
CO2_XS=cell2mat(CO2_XS);
CO2_XS=CO2_XS(:,[1,3:end]);
% O2 XSections
O2_XS=cell(1,5);
fid=fopen('/Users/pps18/Documents/Cross_Sections/Ionization_XS/O2_ionisation.txt');
startRow=9;
formatSpec='%f%f%f%f%f%f%f%f%f%f%f%f';

for i_proc=1:5
    O2_XS{i_proc}=textscan(fid,formatSpec, 'Delimiter', '\t', 'HeaderLines' ,startRow);
    
    O2_XS{i_proc}=cell2mat(O2_XS{i_proc});
    O2_XS{i_proc}=reshape(O2_XS{i_proc}',numel(O2_XS{i_proc})/2,2);

    O2_XS{i_proc}=O2_XS{i_proc}(~isnan(O2_XS{i_proc}(:,1)),:);
end
fclose(fid);

O2_XS_2=[O2_XS{1},zeros(size(O2_XS{1},1),4)];
for i_product=2:5
    ix_xs=arrayfun(@(x) find(x(:,1)==O2_XS_2(:,1)),O2_XS{i_product}(:,1));
    O2_XS_2(ix_xs,i_product+1)=O2_XS{i_product}(:,2);
end

O2_XS=O2_XS_2;
clear O2_XS_2
%% Plotting XS
figure('Renderer', 'painters', 'Position', [10 10 900 400])
set(0,'DefaultLineLineWidth',2);
ax1=plot(H2O_XS(:,1),sum(H2O_XS(:,2:end),2),'DisplayName','H2O');
hold on
ax2=plot(CO_XS(:,1),sum(CO_XS(:,2:end),2),'DisplayName','CO');
hold on
ax3=plot(CO2_XS(:,1),sum(CO2_XS(:,2:end),2),'DisplayName','CO2');
hold on
ax4=plot(O2_XS(:,1),sum(O2_XS(:,2:end),2),'DisplayName','O2');
hold on
set(gca,{'xscale','yscale'},{'log','log'});
set(gca,'Xlim',[10,8e2]);
legend('location','southeast');
set(gca,'Fontsize',15);
xlabel('Energy/eV');
ylabel('Ionization Cross Section/cm^2','Interpreter','tex');

