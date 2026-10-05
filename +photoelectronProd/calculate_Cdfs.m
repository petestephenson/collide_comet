
function [pe_CDF, maxwell_CDF] = calculate_Cdfs(E_grid, elEngDist, maxwellDist)
%CALCULATE_CDFS Calculates cumulative distribution functions.
%   [PE_CDF, MAXWELL_CDF] = CALCULATE_CDFS(E_GRID, ELENGDIST, MAXWELLDIST)
%   calculates the cumulative distribution functions.
%
%   Inputs:
%   - E_grid: Energy grid.
%   - elEngDist: Electron energy distribution.
%   - maxwellDist: Maxwellian distribution.
%
%   Outputs:
%   - pe_CDF: Cumulative distribution function for photoelectron.
%   - maxwell_CDF: Cumulative distribution function for Maxwellian.

nEng = numel(E_grid);

pe_CDF = zeros(nEng,1);
maxwell_CDF = zeros(nEng,1);

for iEng = 2:nEng
    pe_CDF(iEng) = trapz(E_grid(1:iEng), sum(elEngDist(1:iEng,:),2, 'omitnan'));
    maxwell_CDF(iEng) = trapz(E_grid(1:iEng), sum(maxwellDist(1:iEng,:),2));
end
pe_CDF = pe_CDF / pe_CDF(end);
maxwell_CDF = maxwell_CDF/ maxwell_CDF(end);
end
