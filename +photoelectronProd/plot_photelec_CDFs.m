
function plot_photelec_CDFs(E_grid, pe_CDF, maxwell_CDF)
%PLOT_PHOTELEC_CDFS Plots photoelectron and Maxwellian cumulative distribution functions.
%   PLOT_PHOTELEC_CDFS(E_GRID, PE_CDF, MAXWELL_CDF) plots the cumulative 
%   distribution functions for photoelectron and Maxwellian distributions.
%
%   Inputs:
%   - E_grid: Energy grid.
%   - pe_CDF: Cumulative distribution function for photoelectron.
%   - maxwell_CDF: Cumulative distribution function for Maxwellian.


f = figure;
ax = axes(f);
hold(ax, 'on');
plot(ax, E_grid, pe_CDF, 'DisplayName', 'Photoelectron', 'LineWidth', 1.3);
plot(ax, E_grid, maxwell_CDF, 'DisplayName', 'Maxwellian', 'LineWidth', 1.3);
ylabel(ax, 'CDF');
xlabel(ax, 'Energy [eV]');
set(ax, 'FontSize', 13);
xlim(ax, [0,70]);
ylim(ax, [0, 1.1]);
legend(ax);


end
