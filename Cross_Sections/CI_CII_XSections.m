%% Calculating the CII 1335 Cross sections.

% Eqns used are from Shirai et al. (2001)
nEng=300;
E=logspace(0.9,4,nEng);
sig_COCII1335=COCII1335(E);
sig_CO2CII1335=CO2CII1335(E);


%% Importing Cross sections for CI 1561 and 1657 from Mumma et al. (1972)

fid=fopen('/Users/pps18/Documents/Cross_Sections/xs_CO2_CI1657.txt');
formatSpec = ['%f %f'];
startRow=1;
sig_CO2CI1657=textscan(fid,formatSpec,'Delimiter', ' ','HeaderLines' ,startRow);
fclose(fid);
fid=fopen('/Users/pps18/Documents/Cross_Sections/XS_CO2_CI1561.txt');
formatSpec = ['%f %f'];
startRow=1;
sig_CO2CI1561=textscan(fid,formatSpec,'Delimiter', ' ','HeaderLines' ,startRow);
fclose(fid);

%NB factor of 1e18 had been omitted in the text file. Reinserted here
%Column 1 is energy/eV. Column 2 is XS/cm^2.
sig_CO2CI1561=[sig_CO2CI1561{1},sig_CO2CI1561{2}*1e-18];
sig_CO2CI1657=[sig_CO2CI1657{1},sig_CO2CI1657{2}*1e-18];%cm^2
%% CO 4PG total emission XS of the v'=0 level.
%From table 3 in Beegle et al. (1999)
CO_4PG_energy=[8,9,10,11,12,13,14,15,16,18,20,22,24,26,28,30,35,40,45,50,55,60,...
    65,70,80,90,100,110,120,140,160,180,200,220,240,260,280,300,350,400,...
    450,500,550,600,650,700,750];
sig_CO_4PG=[0.341,1.770,3.073,4.000,4.862,5.422,5.755,6.119,6.368,...
    6.713,6.930,6.985,6.967,6.873,6.703,6.496,5.995,5.632,5.308,5.023,...
    4.773,4.553,4.359,4.187,3.895,3.658,3.460,3.291,3.145,2.898,2.692,...
    2.513,2.354,2.210,2.080,1.961,1.852,1.754,1.548,1.389,1.269,1.178,...
    1.110,1.058,1.017,0.985,0.958]*1e-18;

sig_CO_4PG_02=sig_CO_4PG*0.2786; %Band at 1653.2A
%% Setting all Cross Sections to the same Energy grid
% sig_CO2CI1561_new=zeros([nEng,1]);
% sig_CO2CI1657=zeros([nEng,1]);

%interpolating to the regular grid in logspace. XS set to zero for energies
%lower than the threshold. 
sig_CO_4PG_02=interp1(CO_4PG_energy,sig_CO_4PG_02,E,'linear');
sig_CO2CI1561=interp1(sig_CO2CI1561(:,1),sig_CO2CI1561(:,2),E,'linear');
sig_CO2CI1657=interp1(sig_CO2CI1657(:,1),sig_CO2CI1657(:,2),E,'linear');
sig_CO_4PG_02(E<sig_CO_4PG_02(1)|isnan(sig_CO_4PG_02))=0;
sig_CO2CI1561(E<sig_CO2CI1561(1,1)|isnan(sig_CO2CI1561))=0;
sig_CO2CI1657(E<sig_CO2CI1657(1,1)|isnan(sig_CO2CI1657))=0;

%% Writing Data File containing the Cross sections for CO and CO2.
sig_CO=[E',sig_COCII1335',sig_CO_4PG_02'];
sig_CO2=[E',sig_CO2CI1561',sig_CO2CI1657',sig_CO2CII1335'];
headers_CO={'Energy/eV','XS CII1335/cm^2','XS 4PG (0,2) 1653/cm^2'};
headers_CO2={'Energy/eV','XS CI11561/cm^2','XS CI1657/cm^2','XS CII1335/cm^2'};
csvwrite_with_headers('/Users/pps18/Documents/Cross_Sections/sig_CICII_fromCO.txt',sig_CO,headers_CO);
csvwrite_with_headers('/Users/pps18/Documents/Cross_Sections/sig_CICII_fromCO2.txt',sig_CO2,headers_CO2);
%%
% fid=fopen('/Users/pps18/Documents/Cross_Sections/sig_OI_fromCO2.txt','rt');
% xs_CO2=fscanf(fid,'%f%f%f',[3,inf]);
% fclose(fid);
% lines_CO2=[1304,1356];
% sig_CO2_1356=zeros(300,1);
% % for i=1:300
% sig_CO2_1356=interp1(xs_CO2(1,:), xs_CO2(3,:),E,'linear');
% sig_CO2_1356(isnan(sig_CO2_1356))=0;
%     
% % end

%% Plot of Atomic XSections
% 
% % ax1=loglog(E,sig_CO2CI1561,':','DisplayName','CO2 CI 1561');
% % hold on
% % ax2=loglog(E,sig_CO2CI1657,'--','DisplayName','CO2 CI 1657');
% % hold on
% % ax3=loglog(E,sig_COCII1335,'-.','DisplayName','CO CII 1335');
% % hold on
% % ax4=loglog(E,sig_CO2CII1335,'-','DisplayName','CO2 CII 1335');
% % set([ax1,ax2,ax3,ax4],'LineWidth', 2);
% ratio1=sig_CO2CI1561./sig_CO2_1356;
% ratio1(isnan(ratio1))=0;
% ratio2=sig_CO2CI1657./sig_CO2_1356;
% ratio2(isnan(ratio2))=0;
% ax1=semilogx(E,ratio1,'DisplayName','CO2 1561/1356');
% hold on
% ax2=semilogx(E,ratio2,'DisplayName','CO2 1657/1356');
% set([ax1,ax2],'LineWidth', 2);
% % set(gca,'xlim',[10,500]);
% set(gca,'ylim',[0,2.4]);
% xlab=xlabel('Energy/eV');
% ylab=ylabel('Cross Section/$cm^{2}$');
% set(xlab,'Interpreter','Latex');
% set(ylab,'Interpreter','Latex');
% set(gca, 'fontsize', 15);
% % set(ax1,'Xscale', log);
% % set(gca,'yscale',log);
% 
% legend;


%% Shirai et al. Functions
% Line 69 from the table in Shirai et al. (2001)
function sig_CO2CII1335=CO2CII1335(E)
c1=4.310e-3;
c2=7.070e-1;
c3=8.960e-2;
c4=6.720e-1;
E_th=4.3e-2;
sig_CO2CII1335=zeros(size(E));
E=E/1e3;
for i=1:length(E)
    if E(i)>E_th
        sig_CO2CII1335(i)=f2(E(i)-E_th,c1,c2,c3,c4);
    end
end

end


%Line 40 from the table in Shirai et al. (2001)
function sig_COCII1335=COCII1335(E)
c1=1.900e-2;
c2=9.190e-1;
c3=7.900e-2;
c4=8.540e-1;
E_th=3.16e-2;
sig_COCII1335=zeros(size(E));
E=E/1e3;
for i=1:length(E)
    if E(i)>E_th
        sig_COCII1335(i)=f2(E(i)-E_th,c1,c2,c3,c4);
    end
end

end
%Functions below are as defined in Shirai et al. (2001)
function f_3=f3(x,c1,c2,c3,c4,c5,c6)
f_1=f1(x,c1,c2);
factor=1+(x/c3)^(c2+c4)+(x/c5)^(c2+c4);
f_3=f_1/factor;
end
function f_2=f2(x,c1,c2,c3,c4)
f_1=f1(x,c1,c2);
factor=1+(x/c3)^(c2+c4);
f_2=f_1/factor;
end
function f_1=f1(x,c1,c2)
sig0=1e-16;%cm^2
E_R=13.61e-3;%eV Rydberg Constant
f_1=sig0*c1*(x/E_R)^c2;
end
