%% -------- Checking for OS ----------------------
OS_name = system_dependent('getos');
i_windows = strcmp('Microsoft',OS_name(1:9));
if i_windows
    one_drive = 'C:/Users/pete4/OneDrive - Imperial College London/';
else
    one_drive = '/Users/pps18/OneDrive - Imperial College London/';
end

%% Open E and B fields

EB_file = 'EB_fields_210F5.nc';
% EB_file = 'EB_fields_210F1_trim.nc';
ncid = netcdf.open(EB_file);
finfo = ncinfo(EB_file);
E_ave = netcdf.getVar(ncid,0);
B_ave = netcdf.getVar(ncid,1);
% load('E_ave.mat'); 
% load('B_ave.mat');
%% Analyse Trajectories
% i_pl = 1;
% pl_no = [1:10];
% nPls = numel(pl_no);
% X = cell(nPls, 1);
% V = cell(nPls, 1);
% for i_pl = 1:nPls
%     pl_traj = csvread([one_drive,'Fortran/Energy_Degradation/'...
%             'testing123/ioni/pl_traj_ioni',...
%             num2str(pl_no(i_pl))], 0,0);
% 
%     plot_range = pl_traj(2:end,1)~=0;
%     plot_range = [plot_range;false];
% 
%     X{i_pl} = pl_traj(plot_range,1:3);
%     V{i_pl} = pl_traj(plot_range,4:6);
% end
%%
run_no = 1738693;
% traj_file = ['C:\Users\pete4\trajectories1597812.nc'];
traj_file = ['C:\Users\pete4\trajectories',num2str(run_no), '.nc'];

ncid = netcdf.open(traj_file);
finfo = ncinfo(traj_file);
nPls = numel(finfo.Variables)/2;
nPls = 1;
X = cell(nPls, 1);
V = cell(nPls, 1);
for i_pl = 1:nPls
    i_pl = 2;
%     i_pl =6;
    x_id = netcdf.inqVarID(ncid, ['X',num2str(i_pl)]);
    v_id = netcdf.inqVarID(ncid, ['V',num2str(i_pl)]);
    X_temp = netcdf.getVar(ncid, x_id)';
    V_temp = netcdf.getVar(ncid, v_id)';
    plot_range = X_temp(2:end,1)~=0;
%     plot_range(1000:end) = 0;
    plot_range = [plot_range;false];
    nT = find(plot_range, 1, 'last');
    
    X{i_pl} = X_temp(plot_range,:);
    V{i_pl} = V_temp(plot_range,:);
%     if sum(plot_range > 1e5)


%     for i_part = 1:ceil(nT/6e4)
%         i_cell = i_pl+i_part - 1;
%         cell_range = find(and([1:nT]'>6e4*(i_part-1), [1:nT]'<=6e4*i_part));
% %         
%         X{i_cell} = X_temp(cell_range,:);
%         V{i_cell} = V_temp(cell_range,:);
% 
%     end
%     end
end
X = X(i_pl); V = V(i_pl);

%%
% for i_pl = 1:4%5%numel(X)
pls = analyse_traj(X,V, E_ave, B_ave, 1,1);
% end
%%
% save(['pls_Cometary_210F1.mat'], 'pls')


%%
% load('pls_Cometary_210F1.mat');
% load('pls_SW_210F1.mat');
%%
i_pl = 1;
%%

% i_pl = 1;
%% ------ Plot 3D trajectory with Energy Variation ---------------

i_tot_eng = 1;
i_eng_para = 0;
i_eng_perp = 0;
i_pot = 0;
i_time = 0;

clims = [-Inf, Inf];
% E_max = max(color_array(1,:));
% E_min = min(color_array(1,:));

% x_comet = [134,97,97];
x_comet = [134,113,113];

r_com = 17.;
figure('units', 'centimeters', 'position', [2, 2, 15, 15]);
colormap(jet);
% for i_pl = 1:1%9:9
if i_tot_eng
    color_array = [pls(i_pl).Eng'; pls(i_pl).Eng'];
    clabel_name = 'Total Energy [eV]';
elseif i_eng_para
    color_array = [pls(i_pl).Eng_para';pls(i_pl).Eng_para'];
    clabel_name = 'Parallel Energy [eV]';
elseif i_eng_perp 
    color_array = [pls(i_pl).Eng_perp';pls(i_pl).Eng_perp'];
    clabel_name = 'Perp Energy [eV]';
elseif i_pot
    color_array = [pls(i_pl).V_pot_para';pls(i_pl).V_pot_para'];
    clabel_name = 'Potential Energy [eV]';
elseif i_time
    color_array = [pls(i_pl).T';pls(i_pl).T'];
    clabel_name = 'Time [s]';
end

x_plot = (pls(i_pl).X - x_comet) * 10;
ax0 = surface( [x_plot(:,1)';x_plot(:,1)'],...
    [x_plot(:,2)';x_plot(:,2)'],...
    [x_plot(:,3)'; x_plot(:,3)'],...
    color_array, 'facecolor','none','edgecolor','interp','linewidth',0.9);
% end
c = colorbar;
c.Label.String= clabel_name;
c.Label.FontSize = 15;
caxis(clims);
[x_67P,y_67P,z_67P] = sphere;
x_67P = (x_67P *r_com);% Should be *1.7/10 inside brackets
y_67P = (y_67P *r_com);
z_67P = (z_67P *r_com);
hold on

set(gca, 'fontsize', 15);

view([-1,-1,0.2]);

grid on;
daspect([1,1,1]);
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
%% ------ Plot E and B components as function of time
set(0, 'DefaultLineLineWidth', 1.3);
dim_name = {'x','y','z'};
figure;

subplot(2,1,1)
for i_dim = 1:3
    plot(pls(i_pl).T,pls(i_pl).E_pl(:,i_dim), 'DisplayName', ['E_', dim_name{i_dim}]);
    ylabel('E [Vm^{-1}]');
    hold on
end
set(gca, 'Fontsize', 15);
legend;
subplot(2,1,2)
for i_dim = 1:3
    plot(pls(i_pl).T,pls(i_pl).B_pl(:,i_dim), 'DisplayName', ['B_', dim_name{i_dim}]);
    ylabel('B [nT]');
    hold on
end
set(gca, 'Fontsize', 15);

xlabel('T [s]');
legend;
%% ------ Plotting lab velocity and drifts ------------------------
set(0, 'DefaultLineLineWidth', 1.3);
dim_name = {'x','y','z'};
figure;
% set(gca,'DefaultLineLineWidth',2)
for i_dim = 1:3
    subplot(3,1,i_dim);
    plot(pls(i_pl).T,pls(i_pl).V(:,i_dim), 'DisplayName', 'V_{tot}');
    hold on
    plot(pls(i_pl).T,pls(i_pl).V_perp_gyro(:,i_dim), 'DisplayName', 'v_{perp}');
    
    plot(pls(i_pl).T,pls(i_pl).V_para(:,i_dim), 'DisplayName', 'v_{para}');
    plot(pls(i_pl).T,pls(i_pl).v_ExB(:,i_dim), 'DisplayName', 'ExB');
%     hold on
    plot(pls(i_pl).T,pls(i_pl).v_curv(:,i_dim), 'DisplayName', 'Curv');
%     plot(pls(i_pl).T,pls(i_pl).v_gradB(:,i_dim), 'DisplayName', 'gradB');
    ylabel(['V_', dim_name{i_dim}, ' [ms^{-1}]']);
end
xlabel('T [s]');
legend;

%% Plotting rotated velocity and drifts
set(0, 'DefaultLineLineWidth', 1.3);
dim_name = {'B','E','{ExB}'};
figure;

    % set(gca,'DefaultLineLineWidth',2)
plot_line = pls(i_pl).V_gyro_rot;% - pls(i_pl).v_curv_rot;
line_label = 'V_{tot} - v_{ExB}';% -v_{curv}';
for i_dim = 1:3
    subplot(3,1,i_dim);
%     set(gca,'DefaultLineLineWidth',2)

     plot(pls(i_pl).T,plot_line(:,i_dim),...
     'DisplayName', line_label );
    hold on
%     plot(pls(i_pl).T,V_gyro2(:,i_dim),...
%      'DisplayName', 'V_{tot} 2');
%     plot(pls(i_pl).T,pls(i_pl).v_perp_rot(:,i_dim), 'DisplayName', 'v_{perp}');
%     hold on
%     plot(pls(i_pl).T,pls(i_pl).v_para_rot(:,i_dim), 'DisplayName', 'v_{para}');
    plot(pls(i_pl).T,pls(i_pl).v_ExB_rot(:,i_dim), 'DisplayName', 'ExB drift');
%     plot(pls(i_pl).T,v_ExB_gyro(:,i_dim), 'DisplayName', 'ExB');
%     hold on
    plot(pls(i_pl).T,pls(i_pl).v_ExB_rot(:,i_dim) + pls(i_pl).v_curv_rot(:,i_dim), 'DisplayName', 'ExB Curv');
%     plot(pls(i_pl).T,pls(i_pl).v_gradB_rot(:,i_dim), 'DisplayName', 'gradB');
    ylabel(['V_', dim_name{i_dim}, ' [ms^{-1}]']);
    xlim([3,4.3]);

set(gca, 'Fontsize',13);
end
xlabel('T [s]');
legend;
%% 
[EHigh, ELow] = envelope(pls(i_pl).E_rot(:,2),10,'peak'); 
[BHigh, BLow] = envelope(sqrt(sum(pls(i_pl).B_pl.^2,2)),10,'peak'); 
figure;
subplot(2,1,1)
plot(pls(i_pl).T, pls(i_pl).E_rot(:,2));
hold on 
plot(pls(i_pl).T, EHigh);
plot(pls(i_pl).T, ELow);
% xlim([5,5.2]);
subplot(2,1,2)
plot(pls(i_pl).T, sqrt(sum(pls(i_pl).B_pl.^2,2)));
hold on
plot(pls(i_pl).T, BHigh);
plot(pls(i_pl).T, BLow);
% xlim([5,5.2]);


%% Plot Perpendicular Energy
% ix_fix = 1;
ix_fix = find(pls(i_pl).T > 0.6, 1, 'first');

delt_para = pls(i_pl).Eng_para - (pls(i_pl).Eng_para(ix_fix)...
    + pls(i_pl).V_pot_para(ix_fix) - pls(i_pl).V_pot_para);
delt_perp = pls(i_pl).Eng_perp_gyro - (pls(i_pl).Eng_perp_gyro(ix_fix)...
    + pls(i_pl).V_pot_curv(ix_fix) - pls(i_pl).V_pot_curv);
figure;
% yyaxis left
% plot(pls(i_pl).T, delt_para);
hold on
% ylim([-20,5]);
% yyaxis right
plot(pls(i_pl).T, pls(i_pl).Eng_perp, 'DisplayName', 'EngPerp');
plot(pls(i_pl).T, pls(i_pl).Eng_perp_gyro, 'DisplayName', 'EngPerp');
hold on

% ylim([0,25]);
plot(pls(i_pl).T,pls(i_pl).Eng_perp_gyro(ix_fix)...
    + pls(i_pl).V_pot_curv(ix_fix) - pls(i_pl).V_pot_curv, 'DisplayName', 'Pot_{curv}');
xlabel('T [s]');
ylabel('Energy [eV]');
legend;

%% Parallel Energy
Data = pls(i_pl).V_pot_para;
[Maxima,MaxIdx] = findpeaks(Data);
DataInv = 1.01*max(Data) - Data;
[Minima,MinIdx] = findpeaks(DataInv);
%%
% ix_fix = 40000;
% ix_fix = 1;
ix_fix = find(pls(i_pl).T > 0.6, 1, 'first');
[PotHigh, PotLow] = envelope(pls(i_pl).Eng_para(ix_fix)...
    + pls(i_pl).V_pot_para(ix_fix) - pls(i_pl).V_pot_para,50,'peak');
[EngHigh, EngLow] = envelope(pls(i_pl).Eng_para,100,'peak');
figure;
hold on
subplot(3,1,1)
% hax=axes; 
hold on
plot(pls(i_pl).T,pls(i_pl).Eng_para(ix_fix)...
    + pls(i_pl).V_pot_para(ix_fix) - pls(i_pl).V_pot_para, 'DisplayName', 'Pot_{para}');
plot(pls(i_pl).T, pls(i_pl).Eng_para,...
    'DisplayName', 'Eng_{para}');
plot(pls(i_pl).T, EngHigh,'DisplayName', 'EngMax');
% plot(pls(i_pl).T, EngLow,'DisplayName', 'EngMin');
plot(pls(i_pl).T, PotHigh-PotLow,'DisplayName', 'DPot');
% plot(pls(i_pl).T, PotLow,'DisplayName', 'Pot Min');
xlabel('T [s]');
ylabel('Energy [eV]');
legend;
yheight = [-40,100];
% for i_line = 1:numel(Maxima)
% plot([pls(i_pl).T(MaxIdx),pls(i_pl).T(MaxIdx)],yheight, 'color', 'k',...
%     'HandleVisibility', 'off','Linewidth', 0.5);
% plot([pls(i_pl).T(MinIdx),pls(i_pl).T(MinIdx)],yheight, 'color', 'r',...
%     'HandleVisibility', 'off','Linewidth', 0.5);
% end

% xlim([0.5,2]);

Deng_para = pls(i_pl).Eng_para-( pls(i_pl).Eng_para(ix_fix)...
    + pls(i_pl).V_pot_para(ix_fix) - pls(i_pl).V_pot_para);

subplot(3,1,2);
hold on
% xlim([0.5,2]);
plot(pls(i_pl).T, pls(i_pl).Eng_perp, 'DisplayName', 'EngPerp');
plot(pls(i_pl).T, pls(i_pl).Eng_perp_gyro, 'DisplayName', 'EngPerp');
hold on

% ylim([0,25]);
% plot(pls(i_pl).T,pls(i_pl).Eng_perp_gyro(ix_fix)...
%     + pls(i_pl).V_pot_curv(ix_fix) - pls(i_pl).V_pot_curv, 'DisplayName', 'Pot_{curv}');
plot(pls(i_pl).T,- PotLow, 'DisplayName', '- V_{min}');
xlabel('T [s]');
ylabel('Energy [eV]');
yheight = [0,25];
% plot([pls(i_pl).T(MaxIdx),pls(i_pl).T(MaxIdx)],yheight, 'color', 'k',...
%     'HandleVisibility', 'off','Linewidth', 0.5);
% plot([pls(i_pl).T(MinIdx),pls(i_pl).T(MinIdx)],yheight, 'color', 'r',...
%     'HandleVisibility', 'off','Linewidth', 0.5);
legend;

subplot(3,1,3)
magB = sqrt(sum(pls(i_pl).B_pl.^2,2));
plot(pls(i_pl).T, magB);
hold on;
% plot([pls(i_pl).T(MaxIdx),pls(i_pl).T(MaxIdx)],yheight, 'color', 'k',...
%     'HandleVisibility', 'off','Linewidth', 0.5);
% plot([pls(i_pl).T(MinIdx),pls(i_pl).T(MinIdx)],yheight, 'color', 'r',...
%     'HandleVisibility', 'off','Linewidth', 0.5);
ylim([0,1.5e-8]);
% xlim([0.5,2]);
% figure; 
% plot(pls(i_pl).T,Deng_para);
% xlabel('T [s]');
% ylabel('Energy [eV]');

%%
Deng_perp = pls(i_pl).Eng_perp_gyro-(pls(i_pl).Eng_perp_gyro(ix_fix)...
    +  (pls(i_pl).V_pot_curv(ix_fix) - pls(i_pl).V_pot_curv));

figure;
plot(pls(i_pl).T, Deng_para+Deng_perp);
hold on;
% plot(pls(i_pl).T, -Deng_perp);
% plot(pls(i_pl).T, 


%%

% mag_pressure = 
plot_array = pls(i_pl).Eng_gyro;
figure;
subplot(2,1,1)
plot(pls(i_pl).T,pls(i_pl).Eng_para(ix_fix)...
    + pls(i_pl).V_pot_para(ix_fix) - pls(i_pl).V_pot_para, 'DisplayName', 'Pot_{para}');
hold on
% plot(pls(i_pl).T, plot_array, 'DisplayName', 'Eng_{gyro}');
% plot(pls(i_pl).T,pls(i_pl).Eng, 'DisplayName', 'Energy');
% plot(pls(i_pl).T,pls(i_pl).Eng_gyro(ix_fix)...
%      + pls(i_pl).V_pot_curv(ix_fix) - pls(i_pl).V_pot_curv...
%      + pls(i_pl).V_pot_para(ix_fix) - pls(i_pl).V_pot_para,...
%      'DisplayName', 'Pot_{tot}');
plot(pls(i_pl).T, pls(i_pl).Eng_para,...
    'DisplayName', 'Eng_{para}');
% plot(pls(i_pl).T, pls(i_pl).Eng_perp_gyro,...
%     'DisplayName', 'Eng_{perp}');
legend;
subplot(2,1,2);
% yyaxis left
% plot(pls(i_pl).T,sum(pls(i_pl).B_pl.^2,2));
% yyaxis right
% plot(pls(i_pl).T,delt_para, 'DisplayName','D Para');
hold on
% plot(pls(i_pl).T,delt_perp, 'DisplayName', 'D Perp');
plot(pls(i_pl).T,delt_perp + delt_para, 'DisplayName', 'D Perp');
legend;



%% Variation of ExB energy
figure;
yyaxis left
plot(pls(i_pl).T, sqrt(sum(pls(i_pl).v_circ_rot.^2,2)),...
     'DisplayName', 'Circular');
 yyaxis right
 plot(pls(i_pl).T,sqrt(sum(pls(i_pl).B_pl.^2,2)), 'DisplayName', 'Bfield');
ylabel('B Field [T]');
% plot(pls(i_pl).T, pls(i_pl).v_circ_rot(:,i_dim)
%%
pls(i_pl).Eng_ExB = me * sum(pls(i_pl).v_ExB.^2,2)/(2*qe);

figure;
hold on
plot(pls(i_pl).T,pls(i_pl).Eng_ExB, 'DisplayName', 'ExB Energy');
plot(pls(i_pl).T,pls(i_pl).Eng_perp, 'DisplayName', 'Perp Energy');

legend;

%% 
figure;
hold on 
subplot(2,1,1)
plot(pls(i_pl).T, pls(i_pl).MagMom);
ylabel('MagMoment [eV/T]');
subplot(2,1,2)
plot(pls(i_pl).T, pls(i_pl).PitchAngle);
xlabel('T [s]');
ylabel('Pitch Angle');



%%
figure;
yyaxis left
plot(pls(i_pl).T,pls(i_pl).GyroFreq, 'DisplayName', 'GyroFreq');
hold on
ylabel('Gyrofreq [s^{-1}]');
yyaxis right
plot(pls(i_pl).T,sqrt(sum(pls(i_pl).B_pl.^2,2)), 'DisplayName', 'Bfield');
ylabel('B Field [T]');
legend;


