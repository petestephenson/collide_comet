
function plot_photoelectron_EngDist(E_grid, E_th, diffIoniRate)
%PLOT_PHOTOELECTRON_ENGDIST Plots photoelectron energy distribution.
%   PLOT_PHOTOELECTRON_ENGDIST(E_GRID, E_TH, DIFFIONIRATE) plots the
%   photoelectron energy distribution.
%
%   Inputs:
%   - E_grid: Energy grid.
%   - E_th: Threshold energy for each cross-section.
%   - diffIoniRate: Differential ionization rate.
LineName = {'OH+', 'O+', 'H+', 'H2O+'};
f = figure;
ax = axes(f);
hold(ax, 'on');
for i_xs = 1:4
    plot(ax, E_grid - E_th(i_xs), diffIoniRate(:,i_xs), 'DisplayName', LineName{i_xs});
end
ylabel(ax, 'Differential Ioni Rate [s^{-1} eV^{-1}]')
xlabel(ax, 'Photoelectron Energy [eV]');
set(ax, 'XScale', 'log');
set(ax, 'YScale', 'log');
set(ax, 'FontSize', 14);
legend(ax);
end
