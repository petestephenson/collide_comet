
function plot_photoelectron_fraction(E_grid, x_grid, elEngDist, Q, b)
totElEngDist = permute(squeeze(sum( elEngDist,2)), [2,1]);
[X_mesh, Y_mesh] = meshgrid(E_grid, x_grid);
f = figure;
ax = axes(f);
surf(ax, Y_mesh, X_mesh, totElEngDist./totElEngDist(1,:), 'EdgeColor', 'none');

view(ax, 2);
cb = colorbar(ax);
ylabel(ax, 'Energy [eV]');
xlabel(ax, 'X LOS [km]');
cb.Label.String = 'Transmissivity';
set(gca, 'FontSize', 14);
title(['Transmissivity for Q =', num2str(Q),...
    ' s^{-1} at b =',num2str(b) ,' km']); 
end