function solarFlux = interp_sol_flux(sp_Eng, I_Eng, E_grid)
    %INTERP_SOL_FLUX Interpolates solar flux.
%   SOLARFLUX = INTERP_SOL_FLUX(SP_ENG, I_ENG, E_GRID) interpolates solar
%   flux based on energy grid.
%
%   Inputs:
%   - sp_Eng: Energy.
%   - I_Eng: Photon flux.
%   - E_grid: Energy grid.
%
%   Output:
%   - solarFlux: Interpolated solar flux.

    solarFlux = interp1(sp_Eng / 1.6e-19, I_Eng .* 1.6e-19, E_grid);
end