
function [maxwellDist, elEngDist, totIoniRate] = calcEngDist(E_grid, E_th, diffIoniRate)
%CALENGDIST Calculates energy distribution.
%   [MAXWELLDIST, ELENGDIST, TOTIONIRATE] = CALENGDIST(E_GRID, E_TH, DIFFIONIRATE)
%   calculates the energy distribution.
%
%   Inputs:
%   - E_grid: Energy grid.
%   - E_th: Threshold energy for each cross-section.
%   - diffIoniRate: Differential ionization rate.
%
%   Outputs:
%   - maxwellDist: Maxwellian distribution.
%   - elEngDist: Electron energy distribution.
%   - totIoniRate: Total ionization rate.


if ndims(diffIoniRate)
   diffIoniRate = reshape(diffIoniRate, [size(diffIoniRate),1]);
end
kBT = 10 * 1.6e-19; % Electron temperature
lwidth = 1.3;
nEng = numel(E_grid);
nLos = size(diffIoniRate,3);
elEngDist = zeros(nEng, 4, nLos);
for i_xs = 1:4
    elEngDist(:,i_xs,:) = interp1(E_grid-E_th(i_xs), diffIoniRate(:,i_xs,:), E_grid);
end
totIoniRate = trapz(E_grid,  sum(elEngDist,2, 'omitnan'));

maxwellDist = 2 .* sqrt(E_grid .* 1.6e-19 / pi) .* ((1 ./ kBT).^1.5 ) .* ...
    exp(- E_grid .* 1.6e-19 ./ kBT) .* 1.6e-19 .* totIoniRate; 
end
