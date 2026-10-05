
function plot_XSections(XS_Eng, E_th, XS)
%PLOT_XSECTIONS Plots photoionization cross-sections.
%   PLOT_XSECTIONS(XS_ENG, E_TH, XS) plots the photoionization 
%   cross-sections against the produced photoelectron energy.
%
%   Inputs:
%   - XS_Eng: Energy grid.
%   - E_th: Threshold energy for each cross-section.
%   - XS: Photoionization cross-sections.

f = figure;
ax = axes(f);
hold(ax, 'on');
for i_xs = 1:4
    plot(ax, XS_Eng-E_th(i_xs), XS(:,i_xs), 'DisplayName', LineName{i_xs});
end
plot(ax, XS_Eng, sum(XS,2, 'omitnan'), 'DisplayName', 'Tot');
ylabel(ax, 'Cross Section [cm^2]')
xlabel(ax, 'Produced Photoelectron energy [eV]')
legend(ax);
end
function E_grid = get_energy_grid(minEng, maxEng, nEng)
    E_grid = linspace(minEng, maxEng, nEng)';
end