function [E_PE, E_maxwell] = generate_random_energies(nRand,E_grid, pe_CDF, maxwell_CDF)
% GENERATE_RANDOM_ENERGIES Generates random energies based on cumulative distribution functions.
%   [E_PE, E_MAXWELL] = GENERATE_RANDOM_ENERGIES(NRAND, E_GRID, PE_CDF, MAXWELL_CDF)
%   generates random energies based on the provided cumulative distribution functions.
%
%   Inputs:
%   - nRand: Number of random energies to generate.
%   - E_grid: Energy grid.
%   - pe_CDF: Cumulative distribution function for photoelectron.
%   - maxwell_CDF: Cumulative distribution function for Maxwellian.
%
%   Outputs:
%   - E_PE: Random photoelectron energies.
%   - E_maxwell: Random Maxwellian energies.
rEng = rand(nRand, 1);
ixend = find(diff(pe_CDF) < 1e-6, 1, 'first');
ixMaxend = find(diff(maxwell_CDF) < 1e-6, 1, 'first');
E_PE= interp1(pe_CDF(1:ixend), E_grid(1:ixend), rEng);
E_maxwell = interp1(maxwell_CDF(1:ixMaxend), E_grid(1:ixMaxend), rEng);

end
