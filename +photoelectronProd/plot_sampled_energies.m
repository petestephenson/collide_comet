
function plot_sampled_energies(E_grid, elEngDist, maxwellDist, E_PE, E_maxwell)
% PLOT_SAMPLED_ENERGIES Plots sampled energies.
%   PLOT_SAMPLED_ENERGIES(E_GRID, ELENGDIST, MAXWELLDIST, E_PE, E_MAXWELL)
%   plots sampled energies along with theoretical and Maxwellian distributions.
%
%   Inputs:
%   - E_grid: Energy grid.
%   - elEngDist: Electron energy distribution.
%   - maxwellDist: Maxwellian distribution.
%   - E_PE: Random photoelectron energies.
%   - E_maxwell: Random Maxwellian energies.maxEng = 200; % eV

nBins = 2001;
BinEdges = linspace(0,maxEng, nBins);
BinCounts = histcounts(E_PE, BinEdges);
BinCountsMax = histcounts(E_maxwell, BinEdges);

f = figure;
ax = axes(f);
hold(ax, 'on');

h1 = histogram(ax, 'BinEdges', BinEdges, 'BinCounts', BinCounts, 'DisplayName',...
    'PE Random Dist', 'EdgeColor', 'none');
h2 = histogram(ax, 'BinEdges', BinEdges, 'BinCounts', BinCountsMax, 'DisplayName',...
    'Maxwell Dist', 'EdgeColor', 'none');
totelEngDist = sum(elEngDist,2);
scaledEngDist = totelEngDist * max(BinCounts) / max(totelEngDist);

totMaxwell = sum(maxwellDist, 2);
scaledMaxwellian = totMaxwell * max(BinCountsMax) / max(totMaxwell);
plot(ax, E_grid,  scaledEngDist, 'k', 'DisplayName', 'Theoretical Photo el',...
    'LineWidth', 1.5);
plot(ax, E_grid,  scaledMaxwellian, 'k', 'DisplayName', 'Maxwellian',...
    'LineWidth', 1.5);
set(gca, 'YScale', 'log');
ylabel(ax, 'Bin Counts');
xlabel(ax, 'Energy [eV]');

set(ax, 'FontSize',14);
legend(ax, 'Location', 'SouthWest');
set(gca, 'XScale', 'log');

end 