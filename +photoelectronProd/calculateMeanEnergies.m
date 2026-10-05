
function [meanPE_eng,meanMax_eng] =  calculateMeanEnergies(E_grid, maxwellDist, elEngDist)
%CALCULATEMEANENERGIES Calculates mean energies.
%   [MEANPE_ENG, MEANMAX_ENG] = CALCULATEMEANENERGIES(E_GRID, MAXWELLDIST, ELENGDIST)
%   calculates the mean energies.
%
%   Inputs:
%   - E_grid: Energy grid.
%   - maxwellDist: Maxwellian distribution.
%   - elEngDist: Electron energy distribution.
%
%   Outputs:
%   - meanPE_eng: Mean photoelectron energy.
%   - meanMax_eng: Mean Maxwellian energy.



meanPE_eng = trapz(E_grid, sum(elEngDist,2, 'omitnan') .* E_grid) ./...
    trapz(E_grid, sum(elEngDist,2, 'omitnan'));
meanMax_eng = trapz(E_grid, maxwellDist .* E_grid) ./ trapz(E_grid, maxwellDist);
end