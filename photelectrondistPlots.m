OS_name = system_dependent('getos');
i_windows = strcmp('Microsoft',OS_name(1:9));
if i_windows
    one_drive = 'C:/Users/pete4/OneDrive - Imperial College London/';
    lat_folder = 'C:\Users\pete4\Documents\PhD\latPSA\';
else
    one_drive = '/Users/pps18/OneDrive - Imperial College London/';
    lat_folder = '/Volumes/Seagate Expansion Drive/latPSA/latPSA/';
end

% -------------- Import Timed See solar flux at 1AU -----------------------
h = 6.63e-34;
c = 3e8;

dofyear_str = '2015030';
timedseeID = netcdf.open([one_drive,'TIMEDSEE_data/see__L3_',dofyear_str,...
    '_012_01.ncdf']);
% CDF with Variables 4 and 5 being SP.WAV - Spectral wavelength and SP.FLUX
% - Spectral Flux
%% Plot 1 - Solar Flux
sp_wav = netcdf.getVar(timedseeID, 4); 
sp_flux = netcdf.getVar(timedseeID, 5); % Units: Wm^-2nm^-1
f = figure('units', 'centimeters', 'Position', [1,1,15,9]);
ax = axes(f);
plot(ax, sp_wav, sp_flux, 'LineWidth', 1.2);
set(ax, 'yscale', 'log');
xlabel(ax, 'Wavelength [nm]');
ylabel(ax, 'Spectral Flux [Wm^{-2}nm^{-1}]');
set(ax, 'FontSize', 14);

%%
sp_Eng = h*c*1e9./flip(sp_wav); %J

I_lam = flip(sp_flux)*1e-4./sp_Eng; %photons cm^-2 nm^-1
I_E = flip(sp_wav).*I_lam./sp_Eng; %photons s^-1 cm^-2 J^-1
% I_lam_old = sp_flux*1e-4./Eng; 
% I_E_old = sp_wav.*I_lam./Eng; %photons s^-1 cm^-2 J^-1

% XS_lam = h*c*1e9./(flip(excite_XS(:,1)'));
f = figure;
ax = axes(f);
hold(ax, 'on');
plot(sp_Eng/1.6e-19, I_E, 'LineWidth', 1.2);
% plot(Eng/1.6e-19, I_E_old, 'DisplayName', 'old');
xlabel(ax, 'Energy [eV]');
ylabel(ax, 'Photon Flux [photons cm^{-2} s^{-1} J^{-1}]');
% legend(ax);
set(ax, 'FontSize',14);
set(ax, 'YScale', 'log');
set(ax, 'XScale', 'log');
ylim(ax, [1e26, Inf]);
%% -------------- Import photoionization cross section ---------------------
% 
% XS_path = [one_drive, '\Energy_Degradation\XS\H2O_Photoionization\'];
% XS_1b1 = csvread([XS_path, '1b1.csv']);
% XS_1b2 = csvread([XS_path, '1b2.csv']);
% XS_2a1 = csvread([XS_path, '2a1.csv']);
% XS_3a1 = csvread([XS_path, '3a1.csv']);
% 
% nEng = 10001;
% XS_Eng = linspace(0,80, nEng);
% XS = zeros(nEng, 4);
% LineName = {'1b1', '1b2', '2a1', '3a1'};
% XS(:,1) = interp1(XS_1b1(:,1), XS_1b1(:,2), XS_Eng,'linear', 'extrap');
% XS(XS_Eng < XS_1b1(1,1),1) = 0;
% XS(:,2) = interp1(XS_1b2(:,1), XS_1b2(:,2), XS_Eng,'linear', 'extrap');
% XS(XS_Eng < XS_1b2(1,1),2) = 0;
% XS(:,3) = interp1(XS_2a1(:,1), XS_2a1(:,2), XS_Eng,'linear', 'extrap');
% XS(XS_Eng < XS_2a1(1,1),3) = 0;
% XS(:,4) = interp1(XS_3a1(:,1), XS_3a1(:,2), XS_Eng,'linear', 'extrap');
% XS(XS_Eng < XS_3a1(1,1),4) = 0;
% 
% XS = XS * 1e-18; %cm^2 - conversion for Mb
% lam_th = [984, 720.7, 411.4, 907];%\AA
% E_th = h*c*1e10./(lam_th * 1.6e-19);
% 
% XS(XS_Eng<E_th(1),1) = 0;
% lwidth = 1.2;
% % E_th = [12.6, 18.6, 32.6, 14.8];%eV Ionization potential
% f = figure;
% ax = axes(f);
% hold(ax, 'on');
% plot(ax, XS_Eng, sum(XS,2, 'omitnan'), 'DisplayName', 'Tot','LineWidth', lwidth);
% for i_xs = 1:4
%     plot(ax, XS_Eng, XS(:,i_xs), 'DisplayName', LineName{i_xs},...
%         'LineWidth', lwidth);
% end
% 
% xlabel(ax, 'Energy [eV]');
% ylabel(ax, 'Absorption Cross Section [cm^2]');
% legend(ax);
% set(ax, 'FontSize', 14);
% xlim(ax, [0,80]);
%%

H2O_XS = readmatrix([XS_path, 'H2O_ionization_XS.txt'], 'NumHeaderLines', 56);
%% Plot 2 ------ Photoionization XS------------------------
nEng = 1001;
XS_lam = H2O_XS(:,1);
XS = H2O_XS(:, 6:9); %[OH+, O+, H+, H2O+] - cm^2
LineName = {'OH^+', 'O^+', 'H^+', 'H_2O^+'};
E_th = [18.11, 18.65, 18.72, 12.616];
XS_Eng = h * c * 1e10./(flip(XS_lam) * 1.6e-19); %J
XS = flip(XS);
E_grid = linspace(0, 200, nEng)';
SolarFlux = interp1(sp_Eng / 1.6e-19, I_E .* 1.6e-19, E_grid);
XS_int = interp1(log10(XS_Eng), log10(XS), log10(E_grid));
XS_int = 10.^(XS_int);
SolarFlux(isnan(SolarFlux)) = 0;
XS_int(isnan(XS_int)) = 0;
diffIoniRate = SolarFlux .* XS_int;

f = figure('units', 'centimeters', 'Position', [1,1,15,10]);
ax = axes(f);
hold(ax, 'on');
for i_xs = 1:4
    plot(ax, E_grid, XS_int(:,i_xs), 'DisplayName', LineName{i_xs}, 'LineWidth', 1.2);
end
ylabel(ax, 'Ionization Cross Section [cm^2]')
xlabel(ax, 'Energy [eV]');
ylim(ax, [1e-20, Inf]);
set(ax, 'XScale', 'log');
set(ax, 'YScale', 'log');
set(ax, 'FontSize', 14);
legend(ax);%, 'Location', 'southwest');
% I_lam = sp_flux*1e-4./Eng; % Units: Wm^-2nm^-1
%% Differential Photoionization rate @1au
% SolarFlux = interp1(Eng / 1.6e-19, I_E .* 1.6e-19, XS_Eng)'; %Converts flux to eV^-1
% diffIoniRate = SolarFlux .* XS;
f = figure;
ax = axes(f);
hold(ax, 'on');
for i_xs = 1:4
    plot(ax, E_grid, diffIoniRate(:,i_xs), 'DisplayName', LineName{i_xs});
end
%% Photoelectron energy Dist
f = figure;
ax = axes(f);
hold(ax, 'on');
for i_xs = 1:4
    plot(ax, E_grid - E_th(i_xs), diffIoniRate(:,i_xs), 'DisplayName', LineName{i_xs});
end
legend(ax);

%% ---- interpolate onto a fixed energy grid. ----------
% NB all at 1au. 
% E_grid = linspace(0, 1000, nEng);
kBT = 10 * 1.6e-19;
lwidth = 1.3;
elEngDist = zeros(nEng, 4);
for i_xs = 1:4
    elEngDist(:,i_xs) = interp1(E_grid-E_th(i_xs), diffIoniRate(:,i_xs), E_grid);
end
totIoniRate = trapz(E_grid,  sum(elEngDist,2, 'omitnan'));

MaxwellDist = 2 .* sqrt(E_grid .* 1.6e-19 / pi) .* ((1 ./ kBT).^1.5 ) .* ...
    exp(- E_grid .* 1.6e-19 ./ kBT) * 1.6e-19 * totIoniRate; 
f = figure;
ax = axes(f);
hold(ax, 'on');
for i_xs = 1:4
    plot(ax, E_grid, elEngDist(:,i_xs), 'DisplayName', LineName{i_xs}, ...
        'LineWidth', lwidth);
end
plot(ax, E_grid, sum(elEngDist,2), 'DisplayName', 'Total',...
    'LineWidth', lwidth);
plot(ax, E_grid, MaxwellDist, 'DisplayName', 'Maxwell',...
    'LineWidth', lwidth);
legend(ax);
xlabel(ax, 'Energy [eV]');
ylabel(ax, 'Diff. Production Rate [electrons s^{-1} eV^{-1}]');
set(ax, 'xscale', 'log');
set(ax, 'yscale', 'log');
ylim(ax, [1e-13, 1e-6]);
xlim(ax, [1e-1, 1e3]);
set(ax, 'FontSize', 13);
%% ---------- Calculate Mean electron energy ------------------------------


meanPE_eng = trapz(E_grid, sum(elEngDist,2, 'omitnan') .* E_grid) ./...
    trapz(E_grid, sum(elEngDist,2, 'omitnan'));
meanMax_eng = trapz(E_grid, MaxwellDist .* E_grid) ./ trapz(E_grid, MaxwellDist);
%% ---------- Convert Total to CDF ----------------------------------------
pe_CDF = zeros(nEng,1);
maxwell_cdf = zeros(nEng,1);
for iEng = 2:nEng
    pe_CDF(iEng) = trapz(E_grid(1:iEng), sum(elEngDist(1:iEng,:),2, 'omitnan'));
    maxwell_cdf(iEng) = trapz(E_grid(1:iEng), sum(MaxwellDist(1:iEng,:),2));
end


% pe_CDF = pe_CDF / pe_CDF(end);
%%

% pe_CDF(isnan(pe_CDF)) = max(pe_CDF);
pe_CDF = pe_CDF / pe_CDF(end);
maxwell_cdf = maxwell_cdf / maxwell_cdf(end);
f = figure;
ax = axes(f);
hold(ax, 'on');
plot(ax, E_grid, pe_CDF, 'DisplayName', 'Photoelectron', 'LineWidth', 1.3);
plot(ax, E_grid, maxwell_cdf, 'DisplayName', 'Maxwellian', 'LineWidth', 1.3);
ylabel(ax, 'CDF');
xlabel(ax, 'Energy [eV]');
set(ax, 'FontSize', 13);
xlim(ax, [0,70]);
ylim(ax, [0, 1.1]);
legend(ax);

%%
nRand = 1000000;
rEng = rand(nRand, 1);
% ixend = find(pe_CDF == 1, 1, 'first');
ixend = find(diff(pe_CDF) < 1e-6, 1, 'first');
ixMaxend = find(diff(maxwell_cdf) < 1e-6, 1, 'first');
E_rand = interp1(pe_CDF(1:ixend), E_grid(1:ixend), rEng);
E_Maxwell = interp1(maxwell_cdf(1:ixMaxend), E_grid(1:ixMaxend), rEng);
% %%
% figure;
% plot(diff(pe_CDF));
%% Testing the random generation of photoelectron energies
BinEdges = linspace(0,200, 2001);
BinCounts = histcounts(E_rand, BinEdges);
BinCountsMax = histcounts(E_Maxwell, BinEdges);
% BinCounts = ceil(log10(BinCounts));
% BinCounts(BinCounts<0) = 0;


f = figure;
ax = axes(f);
hold(ax, 'on');
% discretize(E_rand,
% hist(E_rand, 50);

h1 = histogram(ax, 'BinEdges', BinEdges, 'BinCounts', BinCounts, 'DisplayName',...
    'PE Random Dist', 'EdgeColor', 'none');
h2 = histogram(ax, 'BinEdges', BinEdges, 'BinCounts', BinCountsMax, 'DisplayName',...
    'Maxwell Dist', 'EdgeColor', 'none');
totelEngDist = sum(elEngDist,2);
scaledEngDist = totelEngDist * max(BinCounts) / max(totelEngDist);
plot(ax, E_grid,  scaledEngDist, 'k', 'DisplayName', 'Theoretical Photo el',...
    'LineWidth', 1.5);
set(gca, 'YScale', 'log');
ylabel(ax, 'Bin Counts');
xlabel(ax, 'Energy [eV]');
% xlim(ax, [0,80]);
% ylim(ax, [1e2, 2e5]);
set(ax, 'FontSize',14);
legend(ax, 'Location', 'SouthWest');
set(gca, 'XScale', 'log');

%% Write Photoelectron Energy Distribution

% nEngSave = find(E_grid == 200);
% pe_CDF_save = pe_CDF(1:nEngSave)/pe_CDF(nEngSave); %Scale to 1
% E_grid_save = E_grid(1:nEngSave);
% % % nScatEng = numel(scat_Eng);
% cdf_name = 'PE_dist.nc';
% nccreate(cdf_name, 'PE_CDF',...
%     'Dimensions',{'nEng', nEngSave}, 'Format','classic');
% nccreate(cdf_name, 'PE_Eng',...
%     'Dimensions', {'nEng', nEngSave}, 'Format','classic');
% ncwrite(cdf_name, 'PE_CDF', pe_CDF_save);
% ncwrite(cdf_name, 'PE_Eng', E_grid_save);

%%
% % nccreate(cdf_name, 'E_th',...
% %     'Dimensions', {'nXS', nXS}, 'Format','classic');
% % nccreate(cdf_name, 'Scat_CDF',...
% %     'Dimensions', {'nScatEng', nScatEng, 'nTheta', nTheta}, 'Format', 'classic');
% % nccreate(cdf_name, 'Scat_Theta',...
% %     'Dimensions', {'nTheta', nTheta}, 'Format', 'classic');
% % nccreate(cdf_name, 'Scat_Eng',...
% %     'Dimensions', {'nScatEng', nScatEng}, 'Format', 'classic');

% ncwrite(cdf_name, 'XS_Eng', E_grid_save);
% ncwrite(cdf_name, 'E_th', E_th);
% ncwrite(cdf_name, 'Scat_CDF', XS_cdf');
% ncwrite(cdf_name, 'Scat_Theta', scat_theta);
% ncwrite(cdf_name, 'Scat_Eng', scat_Eng);


%% ------------- Photoabsorption -----------------------------------------
% ------------- Import XS -----------------------------------------------
H2O_abs_XS = readmatrix([XS_path, 'H2O_absorption_XS.txt'], 'NumHeaderLines', 2);
%% Pull out required data - wavlength and abs XSection

absLam = H2O_abs_XS(:,1);
absEng = h * c * 1e10 ./ (flip(absLam) * 1.6e-19);
absXS = flip(H2O_abs_XS(:,3)) * 1e-18;
figure;
plot(absEng, absXS);
%% Lambert beer Law
% integral of dFlux(lambda, z)/ dz = - sum( sig_i(lambda,z), n_i(z)) * flux(z)  
% --------- Parameters of neutral Dens ------------------------------------
Q = 1.5e27; %s^-1
u_outgas = 1e3; %ms^-1

% ---------- construct an example number density profile ------------------
nX = 250;
x_grid = linspace(-200, 200, nX);
y_pos =  10; %km
z_pos = 10; %km
r_pos = sqrt(x_grid.^2 + y_pos.^2 + z_pos.^2);
n_grid = Q * 1e-12 ./ (4 * pi * u_outgas .* r_pos.^2); %cm^-3
dxgrid = [0,diff(x_grid)];
nEng = 200;
E_grid = linspace(12.64, 124,nEng);

% ------------- 
SolarFlux = interp1(sp_Eng / 1.6e-19, I_E .* 1.6e-19, E_grid)';
intAbsXS = interp1(absEng, absXS, E_grid)';

% ------------- 
optDepth = intAbsXS .* cumsum(n_grid .* dxgrid * 1e5,2);
RedSolarFlux = SolarFlux .* exp(-optDepth);

% ------------- 
XS_int = interp1(XS_Eng, XS, E_grid);
SolarFlux(isnan(SolarFlux)) = 0;
XS_int(isnan(XS_int)) = 0;

% -------------
diffIoniRate = permute(XS_int,[1,3,2]) .* RedSolarFlux;

elEngDist = zeros(nEng, nX, 4);
elEngGrid = linspace(0,120, nEng)';
for i_xs = 1:4
    elEngDist(:,:,i_xs) = interp1(E_grid-E_th(i_xs), diffIoniRate(:,:,i_xs), elEngGrid);
end
elEngDist(isnan(elEngDist)) = 0;

% Resultant photoelectron energy distribution at different positions. 
totElEngDist = permute(sum(elEngDist,3), [2,1]);
[X_mesh, Y_mesh] = meshgrid(elEngGrid, x_grid);
f = figure;
ax = axes(f);
surf(ax, Y_mesh, X_mesh, totElEngDist./totElEngDist(1,:), 'EdgeColor', 'none');

view(ax, 2);
cb = colorbar(ax);
ylabel(ax, 'Energy [eV]');
xlabel(ax, 'X [km]');
cb.Label.String = 'Transmissivity';
set(gca, 'FontSize', 14);
title(['Transmissivity for Q =', num2str(Q),...
    ' s^{-1} at (y,z) = (',num2str(y_pos),',', num2str(z_pos),')  km']); 

%%
nc_file = 'C:/Users/pete4/plas_dens1.nc';
ncid = netcdf.open(nc_file);

pe_eng_id = netcdf.inqVarID(ncid, 'Photoel Eng Dist');
PE_Eng_Dist = netcdf.getVar(ncid, pe_eng_id);

PE_grid_id = netcdf.inqVarID(ncid, 'Energy Boundaries');
PE_grid = netcdf.getVar(ncid, PE_grid_id);
dPE_grid = diff(PE_grid);
BinCounts = PE_Eng_Dist./dPE_grid;
f = figure;
ax0 = axes(f);
scaledDist = totelEngDist * max(BinCounts) / max(totelEngDist);
histogram(ax0, 'BinEdges', PE_grid','BinCounts', BinCounts ,...
                      'DisplayName','BinCounts');
hold(ax0, 'on');
plot(ax0, E_grid, scaledDist, 'LineWidth', 1.3)
