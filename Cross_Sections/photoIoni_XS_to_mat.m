hvXS = importdata('C:\Users\pps18\OneDrive - Imperial College London\Cross_Sections/photoionization_H2O.csv');
wav_grid = [0:100]; %nm
XS = zeros(numel(wav_grid), 4);

for iLine = 1:4
    [lineWav, ixSort] = sort(hvXS.data(:,2*iLine-1));
    lineXS = hvXS.data(ixSort,2*iLine);
    ixNan = isnan(lineWav);
    XS(:,iLine) = interp1(lineWav(~ixNan), lineXS(~ixNan), wav_grid);
end

%%
processName = {'H_2O^+', 'OH^{+}', 'H^{+}', 'O^{+}'};
E_th = [12.62, 18.12, 16.95, 19.00];%eV 'H_2O^+', 'OH^{+}', 'H^{+}', 'O^{+}'
f = figure;
ax = axes(f);
hold(ax, 'on');
for iLine = 1:4 
    scatter(ax, wav_grid, XS(:,iLine), 20, 'filled', 'DisplayName',processName{iLine});
end
set(ax, 'YScale','log');
xlabel(ax, 'Wavelength [nm]');
ylabel(ax, 'Cross Section [cm^{-2}]');
legend(ax, 'Location', 'southeast');
%%
save('photoioni_XS_H2O.mat', "E_th", "XS", "wav_grid", "processName")