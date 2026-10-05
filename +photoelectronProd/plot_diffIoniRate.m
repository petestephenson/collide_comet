
function  plot_diffIoniRate(E_grid, XS_int, diffIoniRate)
%PLOT_DIFFIONIRATE Plots differential ionization rate.
%   PLOT_DIFFIONIRATE(E_GRID, XS_INT, DIFFIONIRATE) plots the differential
%   ionization rate against energy.
%
%   Inputs:
%   - E_grid: Energy grid.
%   - XS_int: Interpolated photoionization cross-sections.
%   - diffIoniRate: Differential ionization rate.

LineName = {'OH+', 'O+', 'H+', 'H2O+'};
f = figure;
ax = axes(f);
hold(ax, 'on');
for i_xs = 1:4
    plot(ax, E_grid, XS_int(:,i_xs), 'DisplayName', LineName{i_xs}, 'LineWidth', 1.2);
end
ylabel(ax, 'Differential Ioni Rate [s^{-1} eV^{-1}]')
xlabel(ax, 'Energy [eV]');
set(ax, 'XScale', 'log');
set(ax, 'YScale', 'log');
set(ax, 'FontSize', 14);
legend(ax);
end