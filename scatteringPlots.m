
% --------------- Load Scattering XS -----------------------
% From Cho et al. (2004)
addpath('C:\Users\pete4\OneDrive - Imperial College London\Energy_Degradation')
scat_XS = csvread('XS\Scat_XS\e_H2O_elas_scattering.csv',1,0);
scat_theta = scat_XS(:,1);
scat_err = scat_XS(:,3:2:end);
scat_XS = scat_XS(:,2:2:end);
scat_Eng = [4,6,8,10,15,20,30,40,50];

scat_XS(1,1:3) = [2.5, 2.5, 2.5] .* scat_XS(2,1:3);

%
Faure_1eV = csvread('XS\Scat_XS\Faure_2004_1eV.csv',1,0);
Faure_025eV = csvread('XS\Scat_XS\Faure_2004_0_25eV.csv',1,0);

scat_XS_1eV = interp1(Faure_1eV(:,1),Faure_1eV(:,2), scat_theta, 'linear', 'extrap'); 
scat_XS_025eV = interp1(Faure_025eV(:,1),Faure_025eV(:,2), scat_theta, 'linear', 'extrap');

scat_Eng = [0.25, 1, scat_Eng];
scat_XS = [scat_XS_025eV, scat_XS_1eV, scat_XS] * 1e-16;

%% Integrated Cross Section

int_XS = trapz(scat_theta * pi / 180, 2 * pi * sind(scat_theta) .* scat_XS);

nEng = size(scat_Eng,2);
nTheta = size(scat_theta,1);
XS_cdf = zeros(nTheta,nEng);

for i_theta = 1:nTheta
    XS_cdf(i_theta,:) = trapz(scat_theta(1:i_theta) * pi / 180, 2 * pi * sind(scat_theta(1:i_theta)) .* scat_XS(1:i_theta,:),1);
end

XS_cdf = XS_cdf./int_XS;


%% Add Scattering cross sections at high and low energy
load('xsections.mat');
scat_Eng = [0, scat_Eng, XS_Eng(end)];
scat_XS_isot = ones(size(scat_theta));
% scat_XS_highEng = ones(size(scat_theta));

isot_cdf = zeros(size(scat_theta));
for i_theta = 1:nTheta
    isot_cdf(i_theta) = trapz(scat_theta(1:i_theta) * pi / 180, 2 * pi * sind(scat_theta(1:i_theta)) .* scat_XS_isot(1:i_theta,:),1);
end
isot_cdf = isot_cdf ./ isot_cdf(end);

XS_cdf = [isot_cdf, XS_cdf, XS_cdf(:,end)];
scat_XS = [ones(nTheta, 1), scat_XS, scat_XS(:,end)];
%% -------------- Plot Differential cross section -------------------------
f = figure('units', 'centimeters','Position', [1,1,26,12]); 
ax = axes(f);
set(0, 'DefaultLineLineWidth',1.2);
% colorOrder = #00876c
% #379469
% #58a066
% #78ab63
% #98b561
% #b8bf62
% #dac767
% #deb256
% #e09d4b
% #e18745
% #e06f45
% #dc574a
% #d43d51
colorOrder = [166,206,227;...
    31,120,180;...
    178,223,138;...
    51,160,44;...
    251,154,153;...
    227,26,28;...
    253,191,111;...
    255,127,0;...
    202,178,214;...
    106,61,154;...
    177,89,40]/255;
% colorOrder = [64,0,75;...
%     118,42,131;...
%     153,112,171;...
%     194,165,207;...
%     231,212,232;...
%     247,247,247;...
%     217,240,211;...
%     166,219,160;...
%     90,174,97;...
%     27,120,55;...
%     0,68,27]/255;
% colorOrder = [ '#00876c';'#379469';'#58a066';'#78ab63';'#98b561';...
% '#b8bf62';'#dac767';'#deb256';'#e09d4b';'#e18745';'#e06f45';'#dc574a';'#d43d51'];
R_scat = 0.75;
engEx = 30;%eV
ixExample = find(scat_Eng == engEx, 1, 'last');
exScat_theta = interp1(XS_cdf(:,ixExample),scat_theta, R_scat);


% ---------- plot options ----------------------------
lef = [0.09, 0.58 0.7];
bot = [0.15, 0.15, 0.1];
wid = 0.4 * ones(1,2);%, 0.3, 0.3];
heig = [0.8, 0.8, 0.8];
cb_wid = 0.01;
cb_heig = 0.8;
yLims = {[1e-18,1e-12], [1e-20, 3e-16]};
xLims = {[1e-2,3e3], [10, 3e2]};

plot_range = {1:6, 7:12};
nPlots = 2;
legCols = [1:2];
ax = cell(nPlots, 1);

%--------------------------------------------------------------------------
%---------------- Plot Cumalative dist Function ---------------------------
%--------------------------------------------------------------------------

iPlot = 1;
ax{iPlot} = subplot('Position', [lef(iPlot), bot(iPlot), wid(iPlot), heig(iPlot)]);

for i_eng = 2:nEng
    plot(scat_theta, scat_XS(:,i_eng), 'DisplayName', [num2str(scat_Eng(i_eng')),'eV'], 'Color', colorOrder(i_eng,:));
hold on
end
legend('NumColumns',3, 'Location', 'northeast');

ylabel('d\sigma / d\Omega  [cm^{2} \cdot sr^{-1}]');
xlabel('Scattering Angle [{\circ}]');
set(gca, 'FontSize', 13);
set(gca, 'yscale', 'log');
xlim([0,180]);
set(ax{iPlot}, 'XTick', [0,30,60,90,120,150,180]);
ylim(ax{iPlot}, [5e-18,1e-13]);
%--------------------------------------------------------------------------
%---------------- Plot Cumalative dist Function ---------------------------
%--------------------------------------------------------------------------
iPlot = 2;
ax{iPlot} = subplot('Position', [lef(iPlot), bot(iPlot), wid(iPlot), heig(iPlot)]);
hold(ax{iPlot}, 'on');
for i_eng = [2:2:nEng, nEng] % Plotting subset of energies for clarity
    plot(ax{iPlot}, scat_theta, XS_cdf(:,i_eng), 'DisplayName', [num2str(scat_Eng(i_eng)), 'eV'],...
        'LineWidth', 1.5, 'Color', colorOrder(i_eng,:));
end
plot(ax{iPlot}, [0, exScat_theta], [R_scat, R_scat],...
    'LineWidth', 1.5, 'Color', colorOrder(ixExample,:), 'LineStyle', '--',...
    'HandleVisibility', 'off');
plot(ax{iPlot}, [exScat_theta, exScat_theta],[0, R_scat], ...
    'LineWidth', 1.5, 'Color', colorOrder(ixExample,:), 'LineStyle', '--',...
    'HandleVisibility', 'off');
t1 = text(5,R_scat + 0.04, 'R_{scatter}', 'Color', colorOrder(ixExample,:),...
    'FontSize', 13, 'FontWeight', 'bold');
t2 = text(exScat_theta + 3, 0.05, '\theta_{scatter}', 'Color', colorOrder(ixExample,:),...
    'FontSize', 13, 'FontWeight', 'bold');
ylabel(ax{iPlot}, 'Cumulative Distribution Function');
xlabel(ax{iPlot},'Scattering Angle [{\circ}]');
set(ax{iPlot}, 'FontSize', 13);
% set(ax{iPlot}, 'yscale', 'log');
xlim(ax{iPlot}, [0,180]);
set(ax{iPlot}, 'XTick', [0,30,60,90,120,150,180]);
% ylim(ax{iPlot}, [5e-18,1e-13]);
legend(ax{iPlot}, 'Location', 'southeast');
%%


pl_eng = 0.1;%eV

XS_pl = interp1(scat_Eng', scat_XS', pl_eng,'linear','extrap');

ix_Eng = find(scat_Eng< pl_eng, 1, 'last');



figure;
hold on
yyaxis left
plot(scat_theta, scat_XS(:,ix_Eng), 'DisplayName', num2str(scat_Eng(ix_Eng)), 'LineWidth', 1.5);
plot(scat_theta, scat_XS(:,ix_Eng+1),'DisplayName', num2str(scat_Eng(ix_Eng+1)), 'LineWidth', 1.5);
ylabel('Diff XS [cm^2sr^-1]');
% plot(Faure_025eV(:,1),Faure_025eV(:,2), 'DisplayName', 'Faure 0.25eV')
% plot(Faure_1eV(:,1),Faure_1eV(:,2), 'DisplayName', 'Faure 1eV')
% plot(scat_theta, XS_pl, 'DisplayName', 'Pl XS');
set(gca, 'yscale', 'log');
yyaxis right
plot(scat_theta, XS_cdf(:,ix_Eng), 'DisplayName', ['cdf',num2str(scat_Eng(ix_Eng))], 'LineWidth', 1.5);
plot(scat_theta, XS_cdf(:,ix_Eng+1),'DisplayName', ['cdf ', num2str(scat_Eng(ix_Eng+1))], 'LineWidth', 1.5);
ylabel('Cumulative Distribution Function');
legend;
xlabel('Scatter Angle [deg]');


%%
load('xsections.mat');


ix_ext_low = XS_Eng <0.03;
ix_ext_high = XS_Eng > 100;
ext_int_XS =  [XS(ix_ext_low,1); int_XS'; XS(ix_ext_high,1)];
ext_Eng = [XS_Eng(ix_ext_low); scat_Eng(2:end-1)'; XS_Eng(ix_ext_high)];
lwidth = 1.3;
f = figure;
ax1 = axes(f);
hold on
l = cell(1,1);

l{1} = plot(ax1,XS_Eng, sum(XS(:,1:2),2),'DisplayName', 'Elastic + J=0\rightarrow 1', 'LineWidth', lwidth);
l{2} = plot(XS_Eng, sum(XS(:,1),2),'DisplayName', 'Elastic', 'LineWidth', lwidth);
l{3} = plot(XS_Eng, sum(XS(:,2),2),'DisplayName', 'J=0\rightarrow 1',  'LineWidth', lwidth);
l{4} = plot(ext_Eng, ext_int_XS, 'DisplayName', '\int d\sigma/d\Omega', 'LineWidth', lwidth); % Integrated Differential XS
% l{5} = plot(scat_Eng, int_XS, 'DisplayName', 'Integrated');
% xlim([0,50]);
for i_l = 1:numel(l)
    l{i_l}.LineWidth = 2;
end

legend;
set(gca, 'yscale', 'log');
set(gca, 'xscale', 'log');
ylabel('Cross Section [cm^2]');
xlabel('Energy [eV]');
set(ax1, 'FontSize', 14);
xlim(ax1, [1e-2,3e2]);


