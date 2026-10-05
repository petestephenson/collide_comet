
function plot_PE_eng_dist(E_grid, maxwellDist, elEngDist)
%PLOT_PE_ENG_DIST Plots photoelectron energy distribution.
%   PLOT_PE_ENG_DIST(E_GRID, MAXWELLDIST, ELENGDIST) plots the photoelectron
%   energy distribution.
%
%   Inputs:
%   - E_grid: Energy grid.
%   - maxwellDist: Maxwellian distribution.
%   - elEngDist: Electron energy distribution.

LineName = {'OH+', 'O+', 'H+', 'H2O+'};

lwidth = 1.3;

f = figure;
ax = axes(f);
hold(ax, 'on');
for i_xs = 1:4
    plot(ax, E_grid, elEngDist(:,i_xs), 'DisplayName', LineName{i_xs}, ...
        'LineWidth', lwidth);
end
plot(ax, E_grid, sum(elEngDist,2), 'DisplayName', 'Total',...
    'LineWidth', lwidth);
plot(ax, E_grid, maxwellDist, 'DisplayName', 'Maxwellian Dist',...
    'LineWidth', lwidth);
legend(ax);
xlabel(ax, 'Energy [eV]');
ylabel(ax, 'Photoelectron Energy Distribution [electrons s^{-1} eV^{-1}]');
set(ax, 'xscale', 'log');
set(ax, 'yscale', 'log');
ylim(ax, [1e-13, 1e-6]);
xlim(ax, [1e-1, 1e3]);
set(ax, 'FontSize', 13);

end
