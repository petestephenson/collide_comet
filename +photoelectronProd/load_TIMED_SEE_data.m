
function [sp_wav, sp_Eng, sp_flux, I_lam, I_Eng] = load_TIMED_SEE_data(srcDir)
%LOAD_TIMED_SEE_DATA Imports Timed See solar flux at 1AU.
%   [SP_WAV, SP_ENG, SP_FLUX, I_LAM, I_ENG] = LOAD_TIMED_SEE_DATA(SRCDIR) 
%   imports Timed See solar flux data from the specified directory and
%   calculates related parameters.
%
%   Inputs:
%   - srcDir: Directory containing the Timed See solar flux data.
%
%   Outputs:
%   - sp_wav: Spectral wavelength.
%   - sp_Eng: Energy corresponding to each wavelength.
%   - sp_flux: Spectral flux.
%   - I_lam: Photon flux per unit wavelength.
%   - I_Eng: Photon flux per unit energy.

% -------------- Import Timed See solar flux at 1AU -----------------------
h = 6.63e-34;
c = 3e8;

dofyear_str = '2015030';
timedseeID = netcdf.open([srcDir,'see__L3_',dofyear_str,'_012_01.ncdf']);
% CDF with Variables 4 and 5 being SP.WAV - Spectral wavelength and SP.FLUX
% - Spectral Flux

sp_wav = double(netcdf.getVar(timedseeID, 4)); 
sp_flux = double(netcdf.getVar(timedseeID, 5)); % Units: Wm^-2nm^-1

sp_Eng = double(h*c*1e9./flip(sp_wav)); %J
I_lam = double(flip(sp_flux)*1e-4./sp_Eng); %photons cm^-2 nm^-1
I_Eng = double(flip(sp_wav).*I_lam./sp_Eng); %photons s^-1 cm^-2 J^-1

end