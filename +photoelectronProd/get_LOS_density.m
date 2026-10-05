
function [x_grid, dxgrid, n_grid, losColDensity] = get_LOS_density(Q, u_outgas, b, nX, maxR)
% GET_LOS_DENSITY Calculates line of sight (LOS) density.
%   [X_GRID, DXGRID, N_GRID, LOSCOLDENSITY] = GET_LOS_DENSITY(Q, U_OUTGAS, B, NX, MAXR)
%   calculates the line of sight (LOS) density.
%
%   Inputs:
%   - Q: Outgassing rate.
%   - u_outgas: Outgassing velocity.
%   - b: Impact parameter.
%   - nX: Number of points in the grid.
%   - maxR: Maximum distance.
%
%   Outputs:
%   - x_grid: Grid along the LOS.
%   - dxgrid: Grid spacing.
%   - n_grid: Density grid.
%   - losColDensity: Cumulative column density along LOS.

disp('NOTE: LOS calc should be updated using Beth et al. (2019) analytical solns.')
x_grid = linspace(-maxR, maxR, nX);

r_pos = sqrt(x_grid.^2 + b.^2);
n_grid = Q * 1e-12 ./ (4 * pi * u_outgas .* r_pos.^2); %cm^-3
dxgrid = [0,diff(x_grid)];

losColDensity = cumtrapz(x_grid, n_grid,2)* 1e5;

end