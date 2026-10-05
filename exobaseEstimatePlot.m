% exobase estimate
sigma = [6.68e-21, 5e-20, 6.68e-19]'; %m^2 - inel Xsecs at 3, 10 and 60~eV
u = 6e2;%ms^-1
Q = logspace(25,29);
% Q = [2e26, 9e27];
Q_sims = [1e26, 3e26, 1.5e27]' .* ones(size(Q));
r_exo = Q .* sigma .* 1e-3 ./ (4 *pi * u); %km
% r_sims = Q_sims .* sigma .* 1e-3 ./ (4 *pi * u); %km

r_surf = 1.7 * ones(size(Q));%km
rmin = min(r_exo, [], 'all');
rmax = max(r_exo, [], 'all');
rRange = linspace(rmin, rmax, numel(Q));
lWidth = 1.5;
f = figure;
ax = axes(f);
hold(ax, 'on');
lStyle = {'-','-','-',':'};
LineName = {'Inel XS at 10 eV', 'Momentum Transfer', 'Inel XS at 10 eV boosted'};
for iEng = [1,3] %1:numel(sigma)
    plot(ax, Q, r_exo(iEng,:), 'LineWidth', lWidth, 'LineStyle',lStyle{iEng},...
        'DisplayName', LineName{iEng});
%     scatter(ax, Q_sims, r_sims(iEng,:), 50,'x', 'k', 'LineWidth', 2);
end
plot(ax, Q, r_surf, 'k', 'LineWidth', lWidth, 'HandleVisibility','off');
plot(ax, Q_sims, rRange, 'LineStyle','--', 'Color','k', 'LineWidth', lWidth,...
    'HandleVisibility','off');
scatter(ax, Q_sims(:,1), [15, 55, 195],160,'k','Marker','x',...
    'LineWidth', lWidth, 'DisplayName','Test Pl Sims');
% ,[5,5,5], 'Marker','none',...
%     'LineStyle', 'none', 'LineWidth', lWidth);
set(ax, 'XScale', 'log');
set(ax, 'YScale', 'log');
xlim(ax, [Q(1), Q(end)]);
ylim(ax, [rmin, rmax]);
set(ax, 'FontSize', 14);
xlabel(ax, 'Outgassing Rate [s^{-1}]');
legend(ax, 'Location', 'southeast');
ylabel(ax, 'Exobase Radius [km]');

%% ------------ Adapted location from Gilet et al 2020 --------------------
% Read Gilet et al 
engDegDir = 'C:\Users\pps18\OneDrive - Imperial College London\Energy_Degradation\';
coldElPoints = csvread([engDegDir, 'Gilet_cold_Detections_v2.csv'],1,0);
acquiredPoints = csvread([engDegDir, 'Gilet_Acquired_Spectra_v2.csv'],1,0);
acquiredPoints = [acquiredPoints; coldElPoints];
% ------- clean data points from vertical line ----------------
% find(coldElPoints(:,1))

%%

rSims = [15, 55, 195];
rSimFit = Q_sims * 6.68e-15 / ( 4 * pi * 1e10 ); %Based on exobase position
rExtrap = logspace(0,3);
r_Range = linspace(rSimFit(1), rSimFit(3));
boostFactor = 6.68e-15 / 5e-16;
Qsims = [1e26, 3e26, 1.5e27]' .* ones(size(rExtrap));
rSimFit = Qsims * 6.68e-15 / ( 4 * pi * 1e10 ); %Based on exobase position
r_surf = 1.7 * ones(size(rExtrap));%km
lef = 0.15;
bot = 0.15;
wid = 0.8;
heig = 0.8;

lWidth = 2;

f = figure('Position',[100,100, 600,500]);
ax = axes(f, "Position",[lef, bot, wid, heig]);



hold(ax, 'on');
scatter(ax, acquiredPoints(:,1), acquiredPoints(:,2), 5,...
    [125,125,145]/255, 'filled', 'DisplayName', 'Acquired MIP Spectra');
% scatter(ax, acquiredPoints(:,1) * boostFactor, acquiredPoints(:,2), 5,...
%     [125,125,145]/255, 'filled', 'DisplayName', 'Acquired MIP Spectra');
scatter(ax, coldElPoints(:,1) , coldElPoints(:,2), 5,...
    'b', 'filled', 'DisplayName', 'cold e^{-}');
plot(ax, rExtrap, rExtrap, 'r','LineStyle', '-', 'LineWidth', lWidth, 'HandleVisibility','off');
% plot(ax, r_Range, r_Range, 'r', 'LineWidth', lWidth, 'HandleVisibility','off');
% plot(ax, rExtrap, rExtrap, 'r','LineStyle', ':', 'LineWidth', lWidth, 'HandleVisibility','off');
% for iQ = 1:3
% plot(ax, r_surf, rExtrap, 'k', 'LineWidth', lWidth, 'HandleVisibility','off');
% plot(ax, rSimFit, rExtrap, 'LineStyle','--', 'Color','k', 'LineWidth', lWidth,...
%         'HandleVisibility','off');
% end
% xlim(ax, [1.99,2.01]);
xlabel(ax, 'L_{en} [km]');
ylabel(ax, 'r_{S/C} [km]');
grid(ax, 'on');

set(ax, {'XScale','YScale'}, {'log', 'log'});

% ax2 = axes(f, "Position",[lef, bot, wid, heig]);
% ax2.XAxisLocation = 'top';
% ax2.Color = 'none';
% xlabel(ax2, 'Q [s^{-1}]');
% set(ax2, {'XScale','YScale'}, {'log', 'log'});
% ax2.YTick = {};
% xlim(ax2, QLims);
% set(ax2, 'Fontsize', 15);

% linkaxes([ax, ax2], 'x')
% xlim(ax, [1e0, 1e4])
QLims = ax.XLim * 4 * pi * 1e10 / 6.68e-15;


set(ax, 'XScale', 'log');
set(ax, 'YScale', 'log');
legend(ax, 'Location', 'southeast');
set(ax, 'Fontsize', 15)

% set(ax, 'Fontsize', 15)

