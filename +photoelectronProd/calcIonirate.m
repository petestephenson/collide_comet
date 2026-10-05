function [E_grid, XS_int, diffIoniRate] = calcIonirate(H2O_XS, XS_Eng, XS,E_th, E_grid, solarFlux)
%CALCIONIRATE Calculates ionization rate.
%   [E_GRID, XS_INT, DIFFIONIRATE] = CALCIONIRATE(H2O_XS, XS_ENG, XS, E_TH, E_GRID, SOLARFLUX) 
%   calculates the ionization rate.
%
%   Inputs:
%   - H2O_XS: H2O photoionization cross-section data.
%   - XS_Eng: Energy grid.
%   - XS: Photoionization cross-sections.
%   - E_th: Threshold energy for each cross-section.
%   - E_grid: Energy grid.
%   - solarFlux: Solar flux.
%
%   Outputs:
%   - E_grid: Energy grid.
%   - XS_int: Interpolated photoionization cross-sections.
%   - diffIoniRate: Differential ionization rate.

if ndims(solarFlux) ==2
    solarFlux = permute(reshape(solarFlux, [size(solarFlux), 1]), [1,3,2]);
    end
h = 6.63e-34;
c = 3e8;
nEng = 1001;
XS_lam = H2O_XS(:,1);
XS = H2O_XS(:, 6:9); %[OH+, O+, H+, H2O+] - cm^2
LineName = {'OH+', 'O+', 'H+', 'H2O+'};
E_th = [18.11, 18.65, 18.72, 12.616];
XS_Eng = h * c * 1e10./(flip(XS_lam) * 1.6e-19); %J
XS = flip(XS);

% SolarFlux = interp1(sp_Eng / 1.6e-19, I_Eng .* 1.6e-19, E_grid);
XS_int = interp1(XS_Eng, XS, E_grid);
solarFlux(isnan(solarFlux)) = 0;
XS_int(isnan(XS_int)) = 0;
diffIoniRate = solarFlux .* XS_int;
end
