
function plot_TIMED_SEE_data(sp_wav, sp_flux)
%PLOT_TIMED_SEE_DATA Plots Timed See solar flux data.
%   PLOT_TIMED_SEE_DATA(SP_WAV, SP_FLUX) plots the Timed See solar flux 
%   data with spectral wavelength on the x-axis and spectral flux on the 
%   y-axis.
%
%   Inputs:
%   - sp_wav: Spectral wavelength.
%   - sp_flux: Spectral flux.

f = figure('Position', [100, 100, 800, 350]);
ax = axes(f);
plot(ax, sp_wav * 10, sp_flux / 10, 'LineWidth', 1.2);
set(ax, 'yscale', 'log');
Ang = char(197);
xlabel(ax, ['Wavelength [', Ang, ']']);
ylabel(ax, ['Spectral Flux [W m^{-2}' , Ang,'^{-1}]']);
set(ax, 'FontSize', 14);
end
function plot_photon_flux(sp_Eng, I_Eng)
%PLOT_PHOTON_FLUX Plots photon flux.
%   PLOT_PHOTON_FLUX(SP_ENG, I_ENG) plots the photon flux as a function of
%   energy.
%
%   Inputs:
%   - sp_Eng: Energy. [eV]
%   - I_Eng: Photon flux.
f = figure;
ax = axes(f);
hold(ax, 'on');
plot(sp_Eng/1.6e-19, I_Eng, 'LineWidth', 1.2);

xlabel(ax, 'Energy [eV]');
ylabel(ax, 'Photon Flux [photons cm^{-2} s^{-1} J^{-1}]');

set(ax, 'FontSize',14);
set(ax, 'YScale', 'log');
set(ax, 'XScale', 'log');
ylim(ax, [1e26, Inf]);

end