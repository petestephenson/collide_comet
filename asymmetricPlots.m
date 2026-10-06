

EngDegDir = 'C:/';
plasDensDir = 'E:/TestPlSimData/';


% peRefNo = '8289790';
% % swRefNo = '8292503';
% swRefNo = '9655074';
% % peCasesNo = {'8289793', '8289804', '8289809','8289809', '8289790'};
% peCasesNo = {'9351241', '9351250','9351253','9351253','9351238'};
% % swCasesNo = { '8292504', '8292506', '8292507', '8292508','8292503'};
% % swCasesNo = {'9361336', '9361342','9361346','9361350','9361329'};
% swCasesNo = {'9655077','9655080','9691056','9707749', '9655074'};


% peRefNo = '681347';
% swRefNo = '699652';
% peRefNo = '777453';
peRefNo = '2857936';
% peRefNo = '1944692';
% swRefNo = '760331';
swRefNo = '2857233';
% peCasesNo = {'681360', '681371', '681400', '681420', '681435', '681347'};
% peCasesNo = {'681360', '681389', '681420', '681435', '681347'};
% peCasesNo = {'777455', '777458', '777462', '777463', '777453'};
% peCasesNo = {'1944473', '1944688', '1944689', '1944690', '1944691'};
peCasesNo = {'1944473',... % -X
            '1944689', ...% +Y
            '1944691',...% +Z
            '1944692',...% -Z
            '2857936'};% Ref
            %'777453'};%Ref

% {'-X', '+Y', '+Z, '-Z','Ref'}
%             '1944688', ...% +X
%             '1944690', ...% -Y

% swCasesNo = {'699654','699664', '699666', '699673', '699681'};
% swCasesNo = {'760332', '760361', '760365', '760367', '760331'};
% swCasesNo = {'2767385', '760361', '760365', '760367', '760331'};
% swCasesNo = {'2767385', '2777452', '2777459', '2777465', '2777466'};

swCasesNo = {'2767385',... % Dayside
            '2777459', ....% +Y
            '2777466', ...% +Z
            '2777467', ...% -Z
            '2857233'}; ...% Ref
% '2777452', ... % +X (nightside)
% '2777465', ...% -Y





%% Load  reference case
% figDir = 'D:/TestPlFigs/Total/';
% figDir = 'D:/TestPlFigs/Total_Jan25/';
figDir = 'E:/TestPlFigs/Total_Jan26/';

plot_lim = 300;
% figDir = 'E:/TestPlFigs/Secondaries/';
% figDir = 'E:/TestPlFigs/Photoelectrons/';
% figDir = 'E:/TestPlFigs/SolarWind/';
% figDir = 'E:/TestPlFigs/SW2/';
% vange
vangleList = {[1,0,0],[0,-1,0],[0,0,1]};
viewNamelist = { 'YZ','XZ','XY'};
slices = {[0],[0],[0]};
[XS, XS_Eng, E_th] = asymPlots.loadXSections();
% load data from
inc_gens = [1,2];
E_cold = 1;
peRefName = [plasDensDir, 'plas_dens', peRefNo, '.nc'];
swRefName = [plasDensDir, 'plas_dens', swRefNo, '.nc'];
[eDepthRadGridref, eDepthGridref, eDepthDistref]  = asymPlots.load_electron_depth_dist(peRefName);


% tGrid = read_timestep_grid(peRefName);
[dtGrid, dtCount] = asymPlots.read_timestep_grid(peRefName);

%%
f = figure();
ax = axes(f);
histogram(ax, 'BinEdges', dtGrid', 'BinCounts', dtCount');
set(ax, 'xScale', 'log')
%%
[gridsRef, plas_dens_E_PE_ref] = asymPlots.loadRunData(peRefName);
[dum, plas_dens_E_SW_ref] = asymPlots.loadRunData(swRefName);


plas_dens_E_ref = plas_dens_E_PE_ref + plas_dens_E_SW_ref;

%%

[pe_bulk_vel, pe_square_vel] = asymPlots.load_velocity_arrays(peRefName);
[sw_bulk_vel, sw_square_vel] = asymPlots.load_velocity_arrays(swRefName);
%%

dE_grid = diff(gridsRef{4});
[plas_dens_tot,bulk_vel_tot, square_vel_tot] = asymPlots.combine_run_velocities(pe_bulk_vel, pe_square_vel,plas_dens_E_PE_ref, ...
    sw_bulk_vel, sw_square_vel, plas_dens_E_SW_ref, dE_grid,[1,2]);
me = asymPlots.get_electron_mass(peRefName);
temperatureRef = asymPlots.calculateTemperature(square_vel_tot, bulk_vel_tot, 0, me);
square_vel_sw_only = permute(sw_square_vel(:,:,:,:,1),[2,3,4,1]);
bulk_vel_sw_only = permute(sw_bulk_vel(:,:,:,:,1),[2,3,4,1]);
temperatureSW_Ref = asymPlots.calculateTemperature(square_vel_sw_only, bulk_vel_sw_only, 0, me);


%%
% plas_dens_E_ref = plas_dens_E_SW;
[qSpherical, qDayside, qUniform, qVec] = asymPlots.loadOutgassingProperties(swRefName);


qDayside = 0;
qUniform = qSpherical;
coldEProd_Ref =  asymPlots.calculateColdeProduction(plas_dens_E_ref, gridsRef, E_cold, qDayside, qUniform,qVec,XS_Eng, XS,E_th, inc_gens);

eImpactIoniref = asymPlots.calcEIoniFreq(plas_dens_E_ref, gridsRef{4}, XS_Eng, XS, inc_gens);
eImpactProdref = asymPlots.calculateEIIProduction(plas_dens_E_ref, gridsRef, qDayside, qUniform,qVec,XS_Eng, XS,inc_gens);


dE_grid = gridsRef{4}(2:end) -gridsRef{4}(1:end-1);
minEng = 0; maxEng = 1;
% calculatePlasmaDensities(grids, plas_dens_E,inc_gens, minEng, maxEng)
[plas_dens_ref, plas_dens_cold_ref] = asymPlots.calculatePlasmaDensities(gridsRef, plas_dens_E_ref, inc_gens, minEng, maxEng);
[X_mesh, Y_mesh, Z_mesh] = meshgrid(gridsRef{1}(1:end-1), gridsRef{2}(1:end-1), gridsRef{3}(1:end-1));

%% Plot Temperature
for iAngle = 1:2
    vangle = vangleList{iAngle};
    [f,ax0, cb] = asymPlots.makeTemperaturePlot(temperatureRef,X_mesh,Y_mesh,Z_mesh, slices, 500,  vangle   , [0, 40]);
%     [f,ax0, cb] = asymPlots.makeTemperaturePlot(temperatureSW_Ref,X_mesh,Y_mesh,Z_mesh, slices, plot_lim,  vangle, [0, 40]);
end

%% Plot temperature vs Radius
[f, ax0, cb, rMids, temp_vs_r, fitresult] = asymPlots.makeTemperature_vs_r_Plot(temperatureRef, plas_dens_ref,  X_mesh, Y_mesh, Z_mesh);
% [f, ax0, cb] = asymPlots.makeTemperature_vs_r_Plot(temperatureSW_Ref, X_mesh, Y_mesh, Z_mesh);
% print(X_mesh)


%%
iCase = 3;

peAsymNo = peCasesNo{iCase};
swAsymNo = swCasesNo{iCase};
peAsymName = [plasDensDir, 'plas_dens', peAsymNo, '.nc'];
swAsymName = [plasDensDir, 'plas_dens', swAsymNo, '.nc'];

[gridsAsym, plas_dens_E_PE_Asym] = asymPlots.loadRunData(swAsymName);


% module load intel-suite/2019.4
% module load netcdf/4.5.2-fortran
%
% module list %%
% for iCase = 5:5
for iCase = [1,2,3,4,5]

    peAsymNo = peCasesNo{iCase};
    swAsymNo = swCasesNo{iCase};
    peAsymName = [plasDensDir, 'plas_dens', peAsymNo, '.nc'];
    swAsymName = [plasDensDir, 'plas_dens', swAsymNo, '.nc'];

    [gridsAsym, plas_dens_E_PE_Asym] = asymPlots.loadRunData(peAsymName);
    [dum, plas_dens_E_SW_Asym] = asymPlots.loadRunData(swAsymName);
    me = asymPlots.get_electron_mass(peAsymName);
    plas_dens_E_Asym = plas_dens_E_PE_Asym + plas_dens_E_SW_Asym;
    [qSpherical, qDayside, qUniform, qVec] = asymPlots.loadOutgassingProperties(swAsymName);
%     plas_dens_E_Asym = plas_dens_E_SW_Asym;
    [XS, XS_Eng, E_th] = asymPlots.loadXSections();
    eImpactIoniAsym = asymPlots.calcEIoniFreq(plas_dens_E_Asym, gridsAsym{4}, XS_Eng, XS, inc_gens);

%     [plas_dens_Asym, plas_dens_cold_Asym] = asymPlots.calculatePlasmaDensities(gridsAsym, plas_dens_E_Asym, dE_grid, inc_gens, minEng, maxEng);
    [plas_dens_Asym, plas_dens_cold_Asym] = asymPlots.calculatePlasmaDensities(gridsAsym, plas_dens_E_Asym, inc_gens, minEng, maxEng);
    if iCase == 5
        qDayside = 0.0;
        qUniform = qSpherical;
    end
    coldEProd_Asym =  asymPlots.calculateColdeProduction(plas_dens_E_Asym, gridsAsym, E_cold, qDayside, qUniform,qVec,XS_Eng, XS,E_th, inc_gens);
    eImpactProdAsym = asymPlots.calculateEIIProduction(plas_dens_E_Asym, gridsAsym, qDayside, qUniform,qVec,XS_Eng, XS,inc_gens);

    % ---- Calculate Temperatures ------
    [pe_bulk_vel, pe_square_vel] = asymPlots.load_velocity_arrays(peAsymName);
    [sw_bulk_vel, sw_square_vel] = asymPlots.load_velocity_arrays(swAsymName);
    bulk_vel_sw_only = permute(sw_bulk_vel(:,:,:,:,1),[2,3,4,1]);
    square_vel_sw_only = permute(sw_square_vel(:,:,:,:,1),[2,3,4,1]);
    dE_grid = diff(gridsAsym{4});
    [plas_dens_tot,bulk_vel_tot, square_vel_tot] = asymPlots.combine_run_velocities(pe_bulk_vel, pe_square_vel,plas_dens_E_PE_Asym, ...
        sw_bulk_vel, sw_square_vel, plas_dens_E_SW_Asym, dE_grid,[1,2]);
    me = asymPlots.get_electron_mass(peAsymName);
    temperatureAsym = asymPlots.calculateTemperature(square_vel_tot, bulk_vel_tot, 0,me);
    temperatureSW_Asym = asymPlots.calculateTemperature(square_vel_sw_only, bulk_vel_sw_only, 0,me);

    disp(['Case ', num2str(iCase),'. qvec: ']);
%     qVec
    % Plasma Density

    plot_array =  plas_dens_Asym;
    clims = [-1,3];
    for iAngle = 1:3
        vangle = vangleList{iAngle};
        viewName = viewNamelist{iAngle};
        iLogArray = 1;
        asymPlots.makePlasDensPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle, clims);
        print(gcf, [figDir, 'plas_dens_PE_SW_Case',num2str(iCase),'_',viewName], '-djpeg');
    end
%
%     % Plasma Density Ratio
%      plot_array =  plas_dens_Asym./plas_dens_ref;
%     for iAngle = 1:3
%         vangle = vangleList{iAngle};
%
%         viewName = viewNamelist{iAngle};
%         iLogArray = 1;
%         asymPlots.makePlasDensRatioPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle);
%         print(gcf, [figDir, 'plas_dens_ratio_PE_SW_Case',num2str(iCase),'_',viewName], '-djpeg');
%     end
%
%     % Cold Plasma Density
%     plot_array = plas_dens_cold_Asym;
%     for iAngle = 1:3
%         vangle = vangleList{iAngle};
%         viewName = viewNamelist{iAngle};
%         iLogArray = 1;
%         asymPlots.makeColdPlasDensPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle);
%         print(gcf, [figDir, 'cold_dens_PE_SW_Case',num2str(iCase),'_',viewName], '-djpeg');
%     end
%
%     % Cold Density Ratio
%     plot_array = plas_dens_cold_Asym./plas_dens_cold_ref;
%     for iAngle = 1:3
%         vangle = vangleList{iAngle};
%         viewName = viewNamelist{iAngle};
%         iLogArray = 1;
%          asymPlots.makeColdDensRatioPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle);
%         print(gcf, [figDir, 'cold_Ratio_PE_SW_Case',num2str(iCase),'_',viewName], '-djpeg');
%     end
%
%     % EII freq
%     plot_array = eImpactIoniAsym;
%     for iAngle = 1:3
%         vangle = vangleList{iAngle};
%         viewName = viewNamelist{iAngle};
%         iLogArray = 1;
%         [f,ax0,cb] = asymPlots.makeEIIFreqPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle);
%         print(gcf, [figDir, 'Eii_freq_PE_SW_Case',num2str(iCase),'_',viewName], '-djpeg');
%     end
%     % EII freq Ratio
%     plot_array = eImpactIoniAsym./eImpactIoniref;
%     for iAngle = 1:3
%         vangle = vangleList{iAngle};
%         viewName = viewNamelist{iAngle};
%         iLogArray = 0;
%         [f,ax0,cb] = asymPlots.makeEIIFreqRatioPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle);
%         print(gcf, [figDir, 'Eii_Ratio_PE_SW_Case',num2str(iCase),'_',viewName], '-djpeg');
%     end
%
%     plot_array = coldEProd_Asym;
%     for iAngle = 1:3
%         vangle = vangleList{iAngle};
%         viewName = viewNamelist{iAngle};
%         iLogArray = 1;
%         [f,ax0,cb] = asymPlots.makecoldEProdPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle);
%         print(gcf, [figDir, 'coldEProd_Case',num2str(iCase),'_',viewName], '-djpeg');
%     end
%     plot_array = coldEProd_Asym./coldEProd_Ref;
%     for iAngle = 1:3
%         vangle = vangleList{iAngle};
%         viewName = viewNamelist{iAngle};
%         iLogArray = 1;
%         [f,ax0,cb] = asymPlots.makecoldEProdRatioPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle);
%         print(gcf, [figDir, 'coldEProdRatio_Case',num2str(iCase),'_',viewName], '-djpeg');
%     end

%     plot_array = eImpactProdAsym;
%     for iAngle = 1:3
%         vangle = vangleList{iAngle};
%         viewName = viewNamelist{iAngle};
%         iLogArray = 1;
%         [f,ax0,cb] = asymPlots.makeEIIproductionPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle);
%         print(gcf, [figDir, 'eiiProd_Case',num2str(iCase),'_',viewName], '-djpeg');
%     end


%     plot_array = temperatureSW_Asym ./temperatureSW_Ref;
%     for iAngle = 1:3
%         vangle = vangleList{iAngle};
%         viewName = viewNamelist{iAngle};
%         iLogArray = 1;
%         [f,ax0,cb] = asymPlots.makeTemperatureRatioPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, vangle);
%
%
%         print(gcf, [figDir, 'SW_tempRatio_Case',num2str(iCase),'_',viewName], '-djpeg');
%     end


%     plot_array = temperatureAsym ./temperatureRef;
%     for iAngle = 1:3
%         vangle = vangleList{iAngle};
%         viewName = viewNamelist{iAngle};
%         iLogArray = 1;
%         [f,ax0,cb] = asymPlots.makeTemperatureRatioPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, vangle);
%
%
%         print(gcf, [figDir, 'totTempRatio_Case',num2str(iCase),'_',viewName], '-djpeg');
%     end



end


%%
% plot_array = plas_dens_Asym./plas_dens_ref;
% plot_lim = 300;
% slices = {[0],[0],[0]};
% iLogArray = 1;
% % vangle = [1,0,0];
% % vangle = [0,-1,0];
% % vangle = [0,0,1];
% [f,ax0] = makePlot(plot_array, X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle);
%%



figDir = 'E:/TestPlFigs/GroupedFigs_Jan26/';
vangleList = {[1,0,0],[0,-1,0],[0,0,1]};
viewNamelist = { 'YZ','XZ','XY'};
%%
for iAngle = 1:3
    [f,ax] = asymPlots.make_SW_temp_plot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    viewName = viewNamelist{iAngle};
    print(f, [figDir, 'swTemp_',viewName], '-djpeg');
end
%%
% make_SW_temp_plot(peCasesNo, swCasesNo, plasDensDir);
for iAngle = 1:3
    [f,ax] = asymPlots.make_SW_tempRatioplot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    viewName = viewNamelist{iAngle};
    print(f, [figDir, 'swTempRatio_',viewName], '-djpeg');
end
%%



dimNames = {'X', 'Y', 'Z'};
for iDim = 1:3
    [f, ax] = asymPlots.make_SW_edf_plot(peCasesNo, swCasesNo, plasDensDir, iDim);
    print(f, [figDir, 'sw_edf', dimNames{iDim}], '-djpeg');
%
end


%%
dimNames = {'X', 'Y', 'Z'};
for iDim = 1:3
    [f,ax0] = asymPlots.make_SW_edfRatioplot(peCasesNo, swCasesNo, plasDensDir, iDim);
     print(f, [figDir, 'sw_edfRatio_', dimNames{iDim}], '-djpeg');
end

%%

for iAngle = 1:3
    [f,ax] = asymPlots.make_grouped_EIIFreqPlot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    viewName = viewNamelist{iAngle};
    print(f, [figDir, 'eiiFreq_Total', viewName], '-djpeg');
end

%%

for iAngle = 1:3
    [f,ax] = asymPlots.make_grouped_EIIFreqRatio_Plot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    viewName = viewNamelist{iAngle};
    print(f, [figDir, 'eiiFreqRatio_Total', viewName], '-djpeg');
end
%%
for iAngle = 1:3
    [f,ax] = asymPlots.make_grouped_nDensityPlot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    viewName = viewNamelist{iAngle};
    print(f, [figDir, 'nDensity', viewName], '-djpeg');
end

%%
for iAngle = 1:3
%     make_grouped_EIIProductionPlot
    [f,ax] = asymPlots.make_grouped_EIIProductionPlot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    viewName = viewNamelist{iAngle};
    print(f, [figDir, 'eiiProduction_Total', viewName], '-djpeg');
end

%%
for iAngle = 1:3
%     make_grouped_EIIProductionPlot
    [f,ax] = asymPlots.make_grouped_EIIProductionRatio_Plot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    viewName = viewNamelist{iAngle};
    print(f, [figDir, 'eiiProductionRatio_Total', viewName], '-djpeg');
end
%%
% for iAngle = 1:3
%     make_grouped_EIIProductionPlot
[f,ax] = asymPlots.make_electrondens_radial_Plot(peCasesNo, swCasesNo, plasDensDir);
% viewName = viewNamelist{iAngle};
%     print(f, [figDir, 'totalElecDensities', viewName], '-djpeg');
% end
%%

[f, ax] = asymPlots.make_electrondens_axisCone_Plot(peCasesNo, swCasesNo, plasDensDir, 3);

%%

[f, ax] = asymPlots.make_electrondens_axisCone_MC_Plot(peCasesNo, swCasesNo, plasDensDir, 1);

%%
[f, ax, prof] = asymPlots.make_electrondens_axisCone_MC_scaledExobase_Plot( ...
    peCasesNo, swCasesNo, plasDensDir, 3, 2e5, true);

%%

[f, ax, prof] = asymPlots.make_electrondens_angleToQVec_shell_MC_AllCases_Plot( ...
    peCasesNo, swCasesNo, plasDensDir, 20, 10, 2e5, 36);

%%
opts = struct();
opts.plotQuantity = "total";
opts.statistic = "median";
opts.nMC = 3e5;
opts.nRadBins = 30;
opts.nAngleBins = 20;
opts.rMin = 0;
opts.rMax = 300;
opts.plotLog10 = true;
opts.clims = [];

[f, ax, out] = asymPlots.make_electrondens_rAngle_QVec_2d_MC_AllCases_Plot_v2( ...
    peCasesNo, swCasesNo, plasDensDir, opts);

% print(f, [figDir, 'cold_electron_dens_vs_qvec'], '-djpeg');
%%
opts = struct();
opts.plotQuantity = "cold";
opts.statistic = "median";
opts.nMC = 3e5;
opts.nRadBins = 30;
opts.nAngleBins = 20;
opts.rMin = 0;
opts.rMax = 300;
opts.plotLog10 = true;
opts.clims = [];

[f, ax, out] = asymPlots.make_electrondens_rAngle_QVec_2d_MC_AllCases_Plot( ...
    peCasesNo, swCasesNo, plasDensDir, opts);

% print(f, [figDir, 'cold_electron_dens_vs_qvec'], '-djpeg');
%%

opts = struct();
opts.plotQuantity = "fraction";
opts.statistic = "median";
opts.nMC = 3e7;
opts.nRadBins = 30;
opts.nAngleBins = 20;
opts.rMin = 0;
opts.rMax = 300;
opts.plotLog10 = false;
opts.clims = [];

[f, ax, out] = asymPlots.make_electrondens_rAngle_QVec_2d_MC_AllCases_Plot( ...
    peCasesNo, swCasesNo, plasDensDir, opts);

% print(f, [figDir, 'cold_fraction_vs_qvec'], '-djpeg');

%%

[f, ax, prof] = asymPlots.make_electrondens_angleToQVec_shell_Plot( ...
    peCasesNo, swCasesNo, plasDensDir, 2, 50, 20);

%%
[f, ax, prof] = asymPlots.make_electrondens_angleToQVec_shell_MC_Plot( ...
    peCasesNo, swCasesNo, plasDensDir, 3, 50, 20, 5e5, 36);

%%
[f, ax, prof] = asymPlots.make_electrondens_angleToQVec_shell_MC_AllCases_Plot( ...
    peCasesNo, swCasesNo, plasDensDir, 50, 20, 5e5, 20);
% print(f, [figDir, 'densities_vs_qvec_50km'], '-djpeg');
%%

[f, ax, prof] = asymPlots.make_electrondens_angleToQVec_shell_MC_AllCases_Plot_v2( ...
    peCasesNo, swCasesNo, plasDensDir, 50, 20, 5e5, 20);
% print(f, [figDir, 'densities_vs_qvec_50km'], '-djpeg');

%%
rEdges = 0:5:300;

opts = struct();
opts.caseNameList = {'-X', '+Y', '+Z', '-Z', 'Ref'};
opts.rMax = 300;
opts.yScaleDensity = "log";
opts.yScaleFraction = "log";
opts.fractionYLim = [1e-3, 1];

[f, ax, prof] = asymPlots.make_electrondens_radial_AllCases_Plot( ...
    peCasesNo, swCasesNo, plasDensDir, rEdges, opts);
print(f, [figDir, 'densities_vs_r'], '-djpeg');
%%
for iAngle = 1:3
%     make_grouped_EIIProductionPlot
    [f,ax] = asymPlots.make_totalelectrons_Plot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    viewName = viewNamelist{iAngle};
    print(f, [figDir, 'totalElecDensities', viewName], '-djpeg');
end
%%
for iAngle = 1:3
%     make_grouped_EIIProductionPlot
    [f,ax] = asymPlots.make_pe_electrons_Plot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    viewName = viewNamelist{iAngle};
    print(f, [figDir, 'pe_ElecDensities', viewName], '-djpeg');
end
%%
for iAngle = 1:3
%     make_grouped_EIIProductionPlot
    [f,ax] = asymPlots.make_pe_electrons_Ratio_Plot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    viewName = viewNamelist{iAngle};
    print(f, [figDir, 'pe_ElecDensities_ratio', viewName], '-djpeg');
end
%%
for iAngle = 1:3
%     make_grouped_EIIProductionPlot
    [f,ax] = asymPlots.make_swelectrons_Plot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    viewName = viewNamelist{iAngle};
    print(f, [figDir, 'swElecDensities', viewName], '-djpeg');
end


%%
for iAngle = 1:3
%     make_grouped_EIIProductionPlot
    [f,ax] = asymPlots.make_sw_warm_electrons_Plot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    viewName = viewNamelist{iAngle};
    print(f, [figDir, 'sw_warmElecDensities_3evmin', viewName], '-djpeg');
end


%%
for iAngle = 1:3
%     make_grouped_EIIProductionPlot
    [f,ax] = asymPlots.make_swelectrons_Ratio_Plot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    viewName = viewNamelist{iAngle};
    print(f, [figDir, 'swElecDensityRatio', viewName], '-djpeg');
end
%%
for iAngle = 1:3
%     make_grouped_EIIProductionPlot
    [f,ax] = asymPlots.make_sw_seccoldelectrons_Plot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    viewName = viewNamelist{iAngle};
    print(f, [figDir, 'sw_sec_cold_ElecDensities', viewName], '-djpeg');
end
%%

for iAngle = 1:3
%     make_grouped_EIIProductionPlot
    [f,ax] = asymPlots.make_sw_seccoldelectrons_Ratio_Plot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    viewName = viewNamelist{iAngle};
    print(f, [figDir, 'sw_sec_cold_ElecDens_ratios', viewName], '-djpeg');
end
%%
for iAngle = 1:3
%     make_grouped_EIIProductionPlot
    [f,ax] = asymPlots.make_totalelectrons_Ratio_Plot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    viewName = viewNamelist{iAngle};
    print(f, [figDir, 'totalElecDensitiy_ratio', viewName], '-djpeg');
end
%%
for iAngle = 1:3
%     make_grouped_EIIProductionPlot
    [f,ax] = asymPlots.make_secondaryelectrons_Plot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    viewName = viewNamelist{iAngle};
    print(f, [figDir, 'secElecDensities_Total', viewName], '-djpeg');
end

%%
for iAngle = 1:3
%     make_grouped_EIIProductionPlot
    [f,ax] = asymPlots.make_secondaryelectrons_Ratio_Plot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    viewName = viewNamelist{iAngle};
    print(f, [figDir, 'secElecDensitiyRatio_Total', viewName], '-djpeg');
end

%%
for iAngle = 1:3
%     make_grouped_EIIProductionPlot
    [f,ax] = asymPlots.make_coldelectrons_Plot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    viewName = viewNamelist{iAngle};
    print(f, [figDir, 'coldElecDensities_Total', viewName,'V2'], '-djpeg');
end
%%
for iAngle = 1:3
%     make_grouped_EIIProductionPlot
    [f,ax] = asymPlots.make_coldelectrons_Ratio_Plot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    viewName = viewNamelist{iAngle};
    print(f, [figDir, 'coldElecDensitiyRatio_Total', viewName], '-djpeg');
end

%%

for iAngle = 1:3
    [f,ax] = asymPlots.make_grouped_cold_e_prod_Plot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    viewName = viewNamelist{iAngle};
%     print(f, [figDir, 'coldElecProduction_Total', viewName], '-djpeg');
end
%%
[f, ax, summary] = asymPlots.make_grouped_cold_e_prod_Plot_v3( ...
    peCasesNo, swCasesNo, plasDensDir, 2, 50, 10, 36);

% print(f, [figDir, 'coldElecProduction_Total_50km', viewName], '-djpeg');
%%

opts = struct();
opts.caseNameList = {'-x', '+y', '+z', '-z', 'sym'};
opts.nBins = 60;
opts.normalization = 'probability';
% opts = struct();
opts.rRange = [0, 55];
opts.nBins = 50;
% opts.normalization = 'probability';
% opts = struct();
% opts.normalization = 'probability';
opts.yScale = 'log';
opts.yMinLog = 1e-5;

[f, ax, out] = asymPlots.make_cold_e_prod_volumeWeighted_Hist_AllCases( ...
    peCasesNo, swCasesNo, plasDensDir, opts);

%%
relAxis = {'X', 'Y', 'Z'};
for iAngle = 1:3
    % Change fontsize to 14

    [f,ax] = asymPlots.make_grouped_ElDepth2d_Plot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    viewName = relAxis{iAngle};
%     print(f, [figDir, 'elDepth2D_', viewName,'_v2'], '-djpeg');
end

%%
opts = struct();
opts.iTau = 3;
opts.medianLevel = 0.5;
% opts.clims = [-1.5, 1.5];
opts.clims = [-2, 1];
opts.rMax = 300;
% opts.cmap = parula;
opts.depthThresh = 1;
opts.showDepthBoundary = true;
opts.boundaryLineSpec = 'k-';
opts.boundaryLineWidth = 2.0;

[f, ax, out] = asymPlots.make_eDepth_QVec_2d_AllCases_Plot( ...
    peCasesNo, swCasesNo, plasDensDir, ...
    [30], 20, opts);
print(f, [figDir, 'elDepth2Ds_cartesian_test'], '-djpeg');

%%
% opts = struct();
% opts.iTau = 3;
% opts.medianLevel = 0.5;
% opts.clims = [-2, 0.5];
% opts.rMax = 300;
%
% [f, ax, out] = asymPlots.make_eDepth_QVec_2d_AllCases_Plot( ...
%     peCasesNo, swCasesNo, plasDensDir, opts);

%%
% opts = struct();
% opts.iTau = 3;
% opts.medianLevel = 0.5;
% opts.exoThreshLog10 = 0;
% opts.smoothHull = true;
% opts.smoothMethod = "sgolay";
% opts.smoothSpan = 14;
% opts.sgolayOrder = 2;
% opts.hullResampleN = 50;
f = figure();
ax = axes(f);
hold(ax, 'on')
for iAngle = 2:2
    exoAll = asymPlots.compute_grouped_ElDepthMedianExobase(peCasesNo, swCasesNo, plasDensDir, iAngle);
    exo3 = exoAll(3).exo;
    plot(ax, exo3.x_smooth, exo3.y_smooth);
    plot(ax, exo3.envelope.x, exo3.envelope.y);



    axis equal
end
% Example: smoothed hull for case 3
% exo3 = outAll(3).exo;
% plot(exo3.x_smooth, exo3.y_smooth);
% axis equal


%%
for iAngle = 1:3
%     [f,ax] = asymPlots.make_grouped_Temperature_Plot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    % Default: combined temperature
%         make_grouped_Temperature_Plot(peCasesNo, swCasesNo, plasDensDir, 1);

    % SW-only temperature (generation 1)
%     [f,ax] = asymPlots.make_grouped_Temperature_Plot(peCasesNo, swCasesNo, plasDensDir, iAngle);
    [f,ax] = asymPlots.make_grouped_Temperature_Plot(peCasesNo, swCasesNo, plasDensDir, iAngle, 'PlotMode', 'sw');

    viewName = viewNamelist{iAngle};
%     print(f, [figDir, 'totTemps_', viewName,'_v2'], '-djpeg');
    print(f, [figDir, 'swTemps_', viewName,'_v2'], '-djpeg');

end

%%
ixXs = 3;
[f,ax] = asymPlots.plot_eDepth_dist(peCasesNo, swCasesNo, plasDensDir, ixXs);
% print(f, [figDir, 'eDepthDist_rbelow200km'], '-djpeg');

%%


ixXs = 3;
[f,ax] = asymPlots.plot_eDepth_dist_asym(peCasesNo, swCasesNo, plasDensDir, ixXs);
% print(f, [figDir, 'eDepthDist_rbelow200km'], '-djpeg');

%%
% exoOutline =
% plot(ax{1},

%%
% [qSpherical, qDayside, qUniform, qvec] = loadOutgassingProperties(ncid);

% Write function to Calculate neutral density at a point
%
% function [f,ax] = make_SW_edf_plot(peCasesNo, swCasesNo, plasDensDir, iDim)
%     f = figure(Position=[0,0,1000, 600]);
%     left = [0.33, 0.61, 0.33, 0.61, 0.05];
%     bottom = [0.57, 0.57, 0.1, 0.1,0.57];
%     cbpos = [0.9,0.1,0.01,0.8];
%     wid = 0.23;
%     heig = 0.4;
%     ax = cell(5,1);
%     cb = cell(5,1);
%     vangleList = {[1,0,0],[0,-1,0],[0,0,1]};
%     slices = {[0],[0],[0]};
%
%     plot_lim=300;
%
%     swRefNo = swCasesNo{5};
%
%     swRefName = [plasDensDir, 'plas_dens', swRefNo, '.nc'];
%     [qSpherical, qDayside, qUniform, qVec] = loadOutgassingProperties(swRefName);
%     [gridsRef, ~] = loadRunData(swRefName);
%     [dum, plas_dens_E_SW_ref] = loadRunData(swRefName);
% %     ixGrid = zeros(3,1);
% %     ixGrid(1) = find(slices{1} <= gridsRef{1}, 1, "first");
% %     ixGrid(2) = find(slices{2} <= gridsRef{2}, 1, "first");
% %     ixGrid(3) = find(slices{3} <= gridsRef{3}, 1, "first");
%     clims = [-4,1];
%     for iCase = [1,2,3,4,5]
%         ax{iCase} = axes(f, 'Position',[left(iCase), bottom(iCase), wid, heig]);
%
%
%
%         peAsymNo = peCasesNo{iCase};
%         swAsymNo = swCasesNo{iCase};
%         peAsymName = [plasDensDir, 'plas_dens', peAsymNo, '.nc'];
%         swAsymName = [plasDensDir, 'plas_dens', swAsymNo, '.nc'];
%         [gridsAsym, ~] = loadRunData(peAsymName);
% %         [spatial_mesh, energy_mesh] = meshgrid(gridsAsym{iDim}(1:end-1), gridsAsym{2}(1:end-1), gridsAsym{4}(1:end-1));
%
%         [dum, plas_dens_E_SW_Asym] = loadRunData(swAsymName);
%
%         [f,ax{iCase}, cb{iCase}] = make_edf_plot(plas_dens_E_SW_Asym(:,:,:,:,1), gridsAsym,iDim, slices,clims, f, ax{iCase});
% %         me = get_electron_mass(peAsymName);
%         if iCase == 5
%             qDayside = 0.0;
%             qUniform = qSpherical;
%         end
%         % ---- Calculate Temperatures ------
%
%
%
% %         me = get_electron_mass(peAsymName);
% %
% %         [f,ax{iCase}, cb{iCase}] = makeTemperaturePlot(temperatureSW_Asym,X_mesh,Y_mesh,Z_mesh, slices, plot_lim,  vangle, f, ax{iCase});
%         if or((iCase == 1), (iCase ==2))
%             xlabel(ax{iCase}, '');
%         end
% %         if mod(iCase,2) == 1
%         delete(cb{iCase});
% %         end
%
%     end
%     cbPlot = colorbar(ax{1}, 'Position', cbpos);
%     plot_str = 'f(E) [cm^{-3}s^{-1}eV^{-1}';
%     cbPlot.Label.String = plot_str;
%     cbPlot.FontSize = 14;
% %     cbPlot.Location = 'eastoutside';
% %     colormap(ax0, jet);
% end
% function [f,ax0, cb] = make_edf_plot(plas_dens_E, grids, iDim, slices, clims,  varargin)
%     if length(varargin) >= 1 && ~isempty(varargin{1})
%         f = varargin{1};
%     else
%         f = figure();
%     end
%
%     if length(varargin) >= 2 && ~isempty(varargin{2})
%         ax0 = varargin{2};
%     else
%         ax0 = axes(f);
%     end
%     dimName = {'X', 'Y', 'Z'};
%     ixGrid = cell(3,1);
%     ixGrid{1} = find(slices{1} <= grids{1}, 1, "first");
%     ixGrid{2} = find(slices{2} <= grids{2}, 1, "first");
%     ixGrid{3} = find(slices{3} <= grids{3}, 1, "first");
%     ixGrid{iDim} = 1:numel(grids{iDim})-1;
% %     nDims = cellfun(@numel, grids) ;
%     smoothedPlasDens_E = smooth_plas_dens_E(plas_dens_E);
%
%     [spatialMesh, energyMesh] = meshgrid(grids{iDim}(1:end-1),grids{4}(1:end-1));
%     plot_array = squeeze(log10(smoothedPlasDens_E(ixGrid{1},ixGrid{2}, ixGrid{3},:)));
%
%     surf(ax0, spatialMesh, energyMesh,plot_array', 'EdgeColor',"none", 'FaceColor',"flat");
%     ylabel(ax0,'Energy [eV]');
%     xlabel(ax0, [dimName{iDim}, ' [km]']);
%     view(ax0, 2);
%     clim(ax0, clims);
%     cb = colorbar(ax0);
%     plot_str = 'f(E) [cm^{-3}s^{-1}eV^{-1}';
%     cb.Label.String = plot_str;
%     cb.FontSize = 14;
%     cb.Location = 'eastoutside';
%     colormap(ax0, jet);
%
% end
% function [f,ax] = make_SW_temp_plot(peCasesNo, swCasesNo, plasDensDir, iAngle)
%
%     f = figure(Position=[0,0,1000, 700]);
%     left = [0.4, 0.7, 0.4, 0.7, 0.1];
%     bottom = [0.6, 0.6, 0.15, 0.15,0.6];
%     wid = 0.25;
%     heig = 0.3;
%     ax = cell(5,1);
%     cb = cell(5,1);
%     vangleList = {[1,0,0],[0,-1,0],[0,0,1]};
%     viewNamelist = { 'YZ','XZ','XY'};
%     slices = {[0],[0],[0]};
%     plot_lim=300;
% %     iAngle = 2;
%     vangle = vangleList{iAngle};
%     swRefNo = swCasesNo{5};
%     swRefName = [plasDensDir, 'plas_dens', swRefNo, '.nc'];
%
%     [qSpherical, qDayside, qUniform, qVec] = loadOutgassingProperties(swRefName);
%     for iCase = [1,2,3,4,5]
%         ax{iCase} = axes(f, 'Position',[left(iCase), bottom(iCase), wid, heig]);
%
%
%
%         peAsymNo = peCasesNo{iCase};
%         swAsymNo = swCasesNo{iCase};
%         peAsymName = [plasDensDir, 'plas_dens', peAsymNo, '.nc'];
%         swAsymName = [plasDensDir, 'plas_dens', swAsymNo, '.nc'];
%         [gridsAsym, ~] = loadRunData(peAsymName);
%         [X_mesh, Y_mesh, Z_mesh] = meshgrid(gridsAsym{1}(1:end-1), gridsAsym{2}(1:end-1), gridsAsym{3}(1:end-1));
%
%         me = get_electron_mass(peAsymName);
%         if iCase == 5
%             qDayside = 0.0;
%             qUniform = qSpherical;
%         end
%         % ---- Calculate Temperatures ------
%
%         [sw_bulk_vel, sw_square_vel] = load_velocity_arrays(swAsymName);
%         bulk_vel_sw_only = permute(sw_bulk_vel(:,:,:,:,1),[2,3,4,1]);
%         square_vel_sw_only = permute(sw_square_vel(:,:,:,:,1),[2,3,4,1]);
%
%         me = get_electron_mass(peAsymName);
%     %     temperatureAsym = calculateTemperature(square_vel_tot, bulk_vel_tot, 0,me);
%         temperatureSW_Asym = calculateTemperature(square_vel_sw_only, bulk_vel_sw_only, 0,me);
%
%         [f,ax{iCase}, cb{iCase}] = makeTemperaturePlot(temperatureSW_Asym,X_mesh,Y_mesh,Z_mesh, slices, plot_lim,  vangle, f, ax{iCase});
%         if or((iCase == 1), (iCase ==2))
%             xlabel(ax{iCase}, '');
%         end
%         if mod(iCase,2) == 1
%             delete(cb{iCase});
%         end
%
%     end
% end
%
% function [f,ax] = make_SW_tempRatioplot(peCasesNo, swCasesNo, plasDensDir, iAngle)
%
%     f = figure();
%     left = [0.1,0.55, 0.1, 0.55];
%     bottom = [0.6, 0.6, 0.15, 0.15];
%     wid = 0.4;
%     heig = 0.3;
%     ax = cell(5,1);
%     cb = cell(5,1);
%     vangleList = {[1,0,0],[0,-1,0],[0,0,1]};
%     viewNamelist = { 'YZ','XZ','XY'};
%     slices = {[0],[0],[0]};
%     plot_lim=300;
% %     iAngle = 2;
%     vangle = vangleList{iAngle};
%     swRefNo = swCasesNo{5};
%     peRefNo = peCasesNo{5};
%     swRefName = [plasDensDir, 'plas_dens', swRefNo, '.nc'];
%     peRefName = [plasDensDir, 'plas_dens', peRefNo, '.nc'];
%
%     [qSpherical, qDayside, qUniform, qVec] = loadOutgassingProperties(swRefName);
%     [sw_bulk_vel, sw_square_vel] = load_velocity_arrays(swRefName);
%     bulk_vel_sw_only = permute(sw_bulk_vel(:,:,:,:,1),[2,3,4,1]);
%     square_vel_sw_only = permute(sw_square_vel(:,:,:,:,1),[2,3,4,1]);
%     me = get_electron_mass(peRefName);
%     temperatureSW_Ref = calculateTemperature(square_vel_sw_only, bulk_vel_sw_only, 0,me);
%     for iCase = [1,2,3,4]
%         ax{iCase} = axes(f, 'Position',[left(iCase), bottom(iCase), wid, heig]);
%
%
%
%         peAsymNo = peCasesNo{iCase};
%         swAsymNo = swCasesNo{iCase};
%         peAsymName = [plasDensDir, 'plas_dens', peAsymNo, '.nc'];
%         swAsymName = [plasDensDir, 'plas_dens', swAsymNo, '.nc'];
%         [gridsAsym, ~] = loadRunData(peAsymName);
%         [X_mesh, Y_mesh, Z_mesh] = meshgrid(gridsAsym{1}(1:end-1), gridsAsym{2}(1:end-1), gridsAsym{3}(1:end-1));
%
%         me = get_electron_mass(peAsymName);
%         if iCase == 5
%             qDayside = 0.0;
%             qUniform = qSpherical;
%         end
%         % ---- Calculate Temperatures ------
%
%         [sw_bulk_vel, sw_square_vel] = load_velocity_arrays(swAsymName);
%         bulk_vel_sw_only = permute(sw_bulk_vel(:,:,:,:,1),[2,3,4,1]);
%         square_vel_sw_only = permute(sw_square_vel(:,:,:,:,1),[2,3,4,1]);
%
%         me = get_electron_mass(peAsymName);
%     %     temperatureAsym = calculateTemperature(square_vel_tot, bulk_vel_tot, 0,me);
%         temperatureSW_Asym = calculateTemperature(square_vel_sw_only, bulk_vel_sw_only, 0,me);
%         plot_array = temperatureSW_Asym ./ temperatureSW_Ref;
%         [f,ax{iCase}, cb{iCase}] = makeTemperatureRatioPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim,  vangle, f, ax{iCase});
%         if mod(iCase,2) == 1
%             delete(cb{iCase})
%         end
%         if iCase <3
%             xlabel(ax{iCase}, '')
%         end
%     end
% end

% function nDens = calcNeutralDens(x,y,z, qUniform, qDayside, qVec)
%     u_gas = 1.0; %kms^-1
% %     Qpt = qUniform;
%     rPts = sqrt(x.^2+y.^2+z.^2);
%
%     xNorm = x./rPts;
%     yNorm = y./rPts;
%     zNorm = z./rPts;
%
%     cos2Max = xNorm.*qVec(1) +  yNorm.*qVec(2) +  yNorm.*qVec(3);
%     cos2Max(cos2Max < 0) = 0.0;
%     Qpts = qUniform + 4.* qDayside.* cos2Max;
%
%
%     nDens = Qpts ./ ( 4* pi* u_gas * (rPts.^2)*1e15); % s^-1/(km^3s^-1)
% end
%
%
% function [qSpherical, qDayside, qUniform, qvec] = loadOutgassingProperties(runFname)
%     disp(runFname);
%     ncid = netcdf.open(runFname);
%     qVecID = netcdf.inqVarID(ncid, 'Q max vector');
%     qUniID = netcdf.inqVarID(ncid, 'OGR Uniform');
%     qDaysideID = netcdf.inqVarID(ncid, 'OGR Dayside');
%     qID = netcdf.inqVarID(ncid, 'OGR');
%     qvec = netcdf.getVar(ncid, qVecID);
%     qUniform = netcdf.getVar(ncid, qUniID);
%     qDayside = netcdf.getVar(ncid, qDaysideID);
%     qSpherical = netcdf.getVar(ncid, qID);
%     netcdf.close(ncid);
% end
% function [f,ax0,cb] = makecoldEProdPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle)
%
%     plot_str = 'P_{e,cold} [cm^{-3}s^{-1}]';
%     if iLogArray
%         plot_array = log10(abs(plot_array));
%         plot_str = ['log_{10}(',plot_str,')'];
%         clims = [-4,1];
%     else
%         clims = [0,10];
%     end
%
%     [f,ax0,cb] = makePlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle);
%     cb.Label.String = plot_str;
%     clim(ax0, clims);
%
% end
%
% function [f,ax0,cb] = makecoldEProdRatioPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle)
%
%     plot_str = 'P_{e,cold,Asym}/P_{e,c,Sym}';
%     if iLogArray
%         plot_array = log10(abs(plot_array));
%         plot_str = ['log_{10}(',plot_str,')'];
%         clims = [-1,1];
%     else
%         clims = [0,10];
%     end
%
%     [f,ax0,cb] = makePlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle);
%     cb.Label.String = plot_str;
%     clim(ax0, clims);
%
% end
% function [f,ax0,cb] = makeEIIFreqRatioPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle)
%
%     plot_str = '\nu_{i,e,A} / \nu_{i,e,S}';
%     if iLogArray
%         plot_array = log10(abs(plot_array));
%         plot_str = ['log_{10}(',plot_str,')'];
%         clims = [-1,1];
%     else
%         clims = [0,10];
%     end
%
%     [f,ax0,cb] = makePlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle);
%     cb.Label.String = plot_str;
%     clim(ax0, clims);
%
% end

% function [f,ax0,cb] = makeEIIproductionPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle)
%
%     plot_str = 'P_{e,EII} [cm^{-3}s^{-1}]';
%     if iLogArray
%         plot_array = log10(abs(plot_array));
%         plot_str = ['log_{10}(',plot_str,')'];
%     end
%
%     [f,ax0,cb] = makePlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle);
%     cb.Label.String = plot_str;
%     clim(ax0, [-2,2]);
% %     clim(ax0, [-9,-6]);
%
% end
%
% function [f,ax0,cb] = makeEIIFreqPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle)
%
%     plot_str = '\nu_{ioni,e} [s^{-1}]';
%     if iLogArray
%         plot_array = log10(abs(plot_array));
%         plot_str = ['log_{10}(',plot_str,')'];
%     end
%
%     [f,ax0,cb] = makePlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle);
%     cb.Label.String = plot_str;
%     clim(ax0, [-8,-5]);
% %     clim(ax0, [-9,-6]);
%
% end

% function [f,ax0,cb] = makeColdDensRatioPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle)
%
%     plot_str = 'N_{e,c,A}/N_{e,c,sym}';
%     if iLogArray
%         plot_array = log10(abs(plot_array));
%         plot_str = ['log_{10}(',plot_str,')'];
%         clims = [-2,2];
%     else
%         clims = [0,10];
%     end
%
%     [f,ax0,cb] = makePlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle);
%     cb.Label.String = plot_str;
%     clim(ax0, clims);
% end
%
% function [f,ax0,cb] = makePlasDensRatioPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle)
%
%     plot_str = 'N_{e,Asym}/N_{e,sym}';
%     if iLogArray
%         plot_array = log10(abs(plot_array));
%         plot_str = ['log_{10}(',plot_str,')'];
%         clims = [-2,2];
%     else
%         clims = [0,10];
%     end
%
%     [f,ax0,cb] = makePlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle);
%     cb.Label.String = plot_str;
%     clim(ax0,clims)
%
% end
% function [f,ax0,cb] = makeColdPlasDensPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle)
%
%     plot_str = 'N_{e,cold} [cm^{-3}]';
%     if iLogArray
%         plot_array = log10(abs(plot_array));
%         plot_str = ['log_{10}(',plot_str,')'];
%         clims = [-2,1];
%
%     else
%         clims = [-inf, inf];
%     end
%
%     [f,ax0,cb] = makePlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle);
%     cb.Label.String = plot_str;
%     clim(ax0, clims);
%
% end
% function [f,ax0,cb] = makePlasDensPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle)
%
%     plot_str = 'N_e [cm^{-3}]';
%     if iLogArray
%         plot_array = log10(abs(plot_array));
%         plot_str = ['log_{10}(',plot_str,')'];
%         clims = [-1,2];
%     else
%         clims = [0,40];
%     end
%
%     [f,ax0,cb] = makePlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle);
%     cb.Label.String = plot_str;
%     clim(ax0, clims);
%
% end
% function [f,ax0, cb] = makeTemperaturePlot(temperature,X_mesh,Y_mesh,Z_mesh, slices, plot_lim,  vangle, varargin)
%     if length(varargin) >= 1 && ~isempty(varargin{1})
%         f = varargin{1};
%     else
%         f = figure();
%     end
%
%     if length(varargin) >= 2 && ~isempty(varargin{2})
%         ax0 = varargin{2};
%     else
%         ax0 = axes(f);
%     end
%
%     hold(ax0, 'on');
%
%     plot_str = 'T_e [eV]';
%
%     clims = [0,40];
%     iLogArray = 0;
%     [f,ax0,cb] = makePlot(temperature,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle, f, ax0);
%     cb.Label.String = plot_str;
%     clim(ax0, clims);
% end
% function [f,ax0, cb] = makeTemperatureRatioPlot(tempRatio,X_mesh,Y_mesh,Z_mesh, slices, plot_lim,  vangle, varargin)
%     if length(varargin) >= 1 && ~isempty(varargin{1})
%         f = varargin{1};
%     else
%         f = figure();
%     end
%
%     if length(varargin) >= 2 && ~isempty(varargin{2})
%         ax0 = varargin{2};
%     else
%         ax0 = axes(f);
%     end
%     plot_str = 'T_{e,asym} / T_{e,sym}';
%
%     clims = [0,2];
%     iLogArray = 0;
%     [f,ax0,cb] = makePlot(tempRatio,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle,f, ax0);
%     cb.Label.String = plot_str;
%     clim(ax0, clims);
% end
%
% function [f, ax0, cb] = makePlot(plot_array, X_mesh, Y_mesh, Z_mesh, slices, plot_lim, iLogArray, vangle, varargin)
%     % Parse optional inputs
%     if length(varargin) >= 1 && ~isempty(varargin{1})
%         f = varargin{1};
%     else
%         f = figure();
%     end
%
%     if length(varargin) >= 2 && ~isempty(varargin{2})
%         ax0 = varargin{2};
%     else
%         ax0 = axes(f);
%     end
%
%     hold(ax0, 'on');
%
%
%     % --- Create Plot Array ------------------------
%
%     plot_array = permute(plot_array,[2,1,3]);
%
%     % set viewing angle
%     view(ax0, vangle);
%
%     % Plot slices of plasma density
%     nslice = sum(cellfun(@numel, slices));
%     ax_slice = slice(ax0, X_mesh,Y_mesh,Z_mesh,...
%         plot_array, slices{1},slices{2},slices{3});
%
%     for i_slice = 1:nslice
%         ax_slice(i_slice).EdgeColor = 'none';
%     end
%
%     % Axes labels
%     xlabel(ax0, 'X [km]');
%     ylabel(ax0, 'Y [km]');
%     zlabel(ax0, 'Z [km]');
%     set(ax0, 'FontSize', 14);
%
%     % Axis limits
%     xlim(ax0, [-plot_lim,plot_lim]);
%     ylim(ax0, [-plot_lim,plot_lim]);
%     zlim(ax0, [-plot_lim,plot_lim]);
%
%     daspect(ax0, [1,1,1])
%     % Set colormap
%     rb_map = zeros(256, 3);
%     rb_map = redblue(256);
%
%     customMap = [255,255,255;...
%         2,5,95;...
%         2,241,229;...
%         23,142,5;...
%         251,246,13;...
%         197,62,6;...
%         105,11,13;...
%         212, 72, 73;...
%         0,0,0];
%     customMap = customMap / 255;
%     cMap = zeros(256,3);
%
%     nStep = round(1+ ([0:size(customMap,1)-1] * (size(cMap,1)-1) / (size(customMap,1)-1)));
%
%     cMap(nStep, :) = customMap;
%
%     for iStep = 1:numel(nStep)-1
%         cMap(nStep(iStep):nStep(iStep+1),:) = interp1(nStep(iStep:iStep+1), customMap(iStep:iStep + 1,:), nStep(iStep):nStep(iStep+1));
%
%     end
%     cb = colorbar(ax0);
%     colormap(ax0, rb_map);
%
%     cb.Location = 'eastoutside';
% end
%%
%
% function [XS, XS_Eng, E_th] = loadXSections()
% % app.oneDrive = [basePath, 'OneDrive - Imperial College London\'];
% %             app.EngDegDir = [app.oneDrive, 'Energy_Degradation\'];
%     XSFolder ='C:/Users/pete4/OneDrive - Imperial College London/Energy_Degradation/XS/';
%     XSFile = 'Coll_XS.nc';
%     ncid = netcdf.open([XSFolder, XSFile]);
%
%     xs_id = netcdf.inqVarID(ncid, 'XS');
%     xsEng_id = netcdf.inqVarID(ncid, 'XS_Eng');
%     Eth_id = netcdf.inqVarID(ncid, 'E_th');
%     XS = netcdf.getVar(ncid, xs_id);
%     XS_Eng = netcdf.getVar(ncid, xsEng_id);
%     E_th = netcdf.getVar(ncid, Eth_id);
%
%     netcdf.close(ncid);
% end
% function [grids, plas_dens_E] = loadRunData(runFname)
%     ncid = netcdf.open(runFname);
%     x_id = netcdf.inqVarID(ncid, 'X Boundaries');
%     x_grid = netcdf.getVar(ncid,x_id);
%     y_id = netcdf.inqVarID(ncid, 'Y Boundaries');
%     y_grid = netcdf.getVar(ncid,y_id);
%     z_id = netcdf.inqVarID(ncid, 'Z Boundaries');
%     z_grid = netcdf.getVar(ncid,z_id);
%     E_id = netcdf.inqVarID(ncid, 'Energy Boundaries');
%     E_grid = netcdf.getVar(ncid,E_id);
%
%     grids = {x_grid, y_grid, z_grid, E_grid};
%
%     plas_id = netcdf.inqVarID(ncid, 'Plasma Density');
%     plas_dens_E = netcdf.getVar(ncid,plas_id);
%
%     netcdf.close(ncid);
% end


%
% function eImpactIoni = calcEIoniFreq(plas_dens_E, E_grid, XS_Eng, XS, inc_gens)
% % Calculate e-impact ioni rate
%     qe = 1.6e-19;
%     me = 9.11e-31;
%     % Mid points of energy grid
%     E_mid = 0.5 * (E_grid(1:end-1) + E_grid(2:end));
%     E_mid = reshape(E_mid, 1, 1, 1, []);
%     ix_ioni = 7;
%
%     % Get coll XS.
%     XS_interped = reshape(interp1(XS_Eng, sum(XS(:,ix_ioni:end),2), E_mid), ...
%         1,1,1,[]);
%     % Convert from edf to pl flux
%     elFlux = sum(plas_dens_E(:,:,:,:,inc_gens), 5) .* ...
%         sqrt(E_mid * 2 * qe ./ me) * 1e2;
%     % Integrate over energy
%     eImpactIoni = trapz(reshape(E_mid,[],1), XS_interped .* elFlux, 4);
% end
%
% function coldElecProd = calculateColdeProduction(plas_dens_E, grids, E_cold, qDayside, qUniform,qVec,XS_Eng, XS,E_th, inc_gens)
%
%     % Load cross-section data
% %     load([EngDegDir, '\xsections.mat']);
%     ix_ioni = 7;
%
%     nEng = numel(grids{4}) - 1;
%     nXS = numel(E_th);
%
%         % Calculating sum(Int f(E) * XS(E) n_H2O(x) dE)
%     xMid = reshape(0.5 * (grids{1}(2:end) + grids{1}(1:end-1)), [], 1);
%     yMid = reshape(0.5 * (grids{2}(2:end) + grids{2}(1:end-1)), 1, []);
%     zMid = reshape(0.5 * (grids{3}(2:end) + grids{3}(1:end-1)), 1, 1, []);
%     rMid = sqrt(xMid.^2 + yMid.^2 + zMid.^2);
% %     nDens = Q * 1e-15 / (4 * pi * u_gas .* rMid.^2);
%     nDens = calcNeutralDens(xMid, yMid, zMid, qUniform, qDayside, qVec);
%
%     E_grid_mid = 0.5 * (grids{4}(2:end) + grids{4}(1:end-1));
%     dE_grid = (grids{4}(2:end) - grids{4}(1:end-1));
%
%     % Initialize XS_total
%     XS_total = zeros(1, nEng);
%
%     for i_xs = 1:nXS
%         if E_th(i_xs) > 0
%             % Determine Limits for each process
%             E_lims = [max(E_th(i_xs), E_cold), E_cold + E_th(i_xs)];
%
%             if i_xs >= ix_ioni
%                 E_lims(2) = E_lims(2) + E_cold;
%             end
%
%             ix_E_th = find(grids{4} > E_lims(1), 1, 'first') - 1;
%             ix_E_cold = find(grids{4} > E_lims(2), 1, 'first') - 1;
%
%             eng_low = 0.5 * (E_lims(1) + grids{4}(ix_E_th + 1));
%             eng_upp = 0.5 * (E_lims(2) + grids{4}(ix_E_cold));
%
%             low_factor = (grids{4}(ix_E_th + 1) - E_lims(1)) / ...
%                 (grids{4}(ix_E_th + 1) - grids{4}(ix_E_th));
%             upp_factor = (E_lims(2) - grids{4}(ix_E_cold)) / ...
%                 (grids{4}(ix_E_cold + 1) - grids{4}(ix_E_cold));
%
%             XS_Eng2 = E_grid_mid';
%             XS_Eng2([ix_E_th, ix_E_cold]) = [eng_low, eng_upp];
%             XS_interp = interp1(XS_Eng, XS(:, i_xs), XS_Eng2);
%
%             XS_interp([ix_E_th, ix_E_cold]) = XS_interp([ix_E_th, ix_E_cold]) .* ...
%                 [low_factor, upp_factor];
%             XS_interp([1:ix_E_th-1, ix_E_cold + 1:nEng]) = 0;
%
%             XS_total = XS_total + XS_interp;
%         end
%     end
%     qe = 1.6e-19;
%     me = 9.11e-31;
%     v_eng = sqrt(2 * XS_Eng2 * abs(qe) / me);
%     temp = reshape(XS_total .* dE_grid' .* v_eng * 1e2, 1, 1, 1, []);
%
% %     for iRun = 1:2
%     coldEProdPerMol = sum(plas_dens_E(:,:,:,:,inc_gens) .* temp, [4, 5]);
%     coldElecProd = coldEProdPerMol .* nDens;
% %     end
% end
% end

%
% function eiiProd = calculateEIIProduction(plas_dens_E, grids, qDayside, qUniform,qVec,XS_Eng, XS,inc_gens)
% %
% %     % Load cross-section data
% % %     load([EngDegDir, '\xsections.mat']);
% %     ix_ioni = 7;
% %
% %     nEng = numel(grids{4}) - 1;
% %     nXS = numel(E_th);
% %
%         % Calculating sum(Int f(E) * XS(E) n_H2O(x) dE)
%     E_grid  = grids{4};
%     xMid = reshape(0.5 * (grids{1}(2:end) + grids{1}(1:end-1)), [], 1);
%     yMid = reshape(0.5 * (grids{2}(2:end) + grids{2}(1:end-1)), 1, []);
%     zMid = reshape(0.5 * (grids{3}(2:end) + grids{3}(1:end-1)), 1, 1, []);
% %     rMid = sqrt(xMid.^2 + yMid.^2 + zMid.^2);
% %     nDens = Q * 1e-15 / (4 * pi * u_gas .* rMid.^2);
%     nDens = calcNeutralDens(xMid, yMid, zMid, qUniform, qDayside, qVec);
%     eImpactIoni = calcEIoniFreq(plas_dens_E, E_grid, XS_Eng, XS, inc_gens);
%     eiiProd = nDens .* eImpactIoni; % cm^-3s^-1
%
% end
%
%
% function [plas_dens, plas_dens_cold] = calculatePlasmaDensities(grids, plas_dens_E, dE_grid, inc_gens, minEng, maxEng)
%     % Initialize cell arrays
% %     inc_eng = cell(2, 3);
% %     plas_dens_all = cell(2, 3);
% %     plas_dens_cold_all = cell(2, 3);
%
% %     for iOGR = 1:3
% %         for iRun = 1:2
%     ix_min_eng = find(grids{4} > minEng, 1, 'first') - 1;
%     ix_max_eng = find(grids{4} <= maxEng, 1, 'last') - 1;
%     inc_eng = [ix_min_eng:ix_max_eng]; % Energy Range to plot
% %         end
%
% %         for iRun = 1:2
% %     print(inc_eng)
%     plas_dens = sum(plas_dens_E(:,:,:,:,inc_gens) ...
%         .* reshape(dE_grid, 1, 1, 1, []), [4, 5]);
%     plas_dens_cold = sum(plas_dens_E(:,:,:,inc_eng, inc_gens) ...
%         .* reshape(dE_grid(inc_eng), 1, 1, 1, []), [4, 5]);
%
% %     nDims = size(plas_dens);
%     plot_array_smooth = zeros(nDims);
%     plot_array_cold_smooth = zeros(nDims);
%     temp_array = plas_dens([1, 1:end, end], [1, 1:end, end], [1, 1:end, end]);
%     temp_array_cold = plas_dens_cold([1, 1:end, end], [1, 1:end, end], [1, 1:end, end]);
%     for iX = 1:3
%         ix_X = [1:nDims(1)] + iX - 1;
%         for iY = 1:3
%             ix_Y = [1:nDims(2)] + iY - 1;
%             for iZ = 1:3
%                 ix_Z = [1:nDims(3)] + iZ - 1;
%                 plot_array_smooth = plot_array_smooth + ...
%                     temp_array(ix_X, ix_Y, ix_Z);
%                 plot_array_cold_smooth = plot_array_cold_smooth + ...
%                     temp_array_cold(ix_X, ix_Y, ix_Z);
%             end
%         end
%     end
%     plas_dens = plot_array_smooth / 27;
%     plas_dens_cold = plot_array_cold_smooth / 27;
% %         end
% %     end
% end

% function smoothedPlasDens = smooth_plas_dens_E(plas_dens_E)
%     nDims = size(plas_dens_E);
%     plot_array_smooth = zeros(size(plas_dens_E));
%
%
%     temp_array = plas_dens_E([1, 1:end, end], [1, 1:end, end], [1, 1:end, end],:);
%
%     for iX = 1:3
%         ix_X = [1:nDims(1)] + iX - 1;
%         for iY = 1:3
%             ix_Y = [1:nDims(2)] + iY - 1;
%             for iZ = 1:3
%                 ix_Z = [1:nDims(3)] + iZ - 1;
%                 plot_array_smooth = plot_array_smooth + ...
%                     temp_array(ix_X, ix_Y, ix_Z,:);
%
%             end
%         end
%     end
%     smoothedPlasDens = plot_array_smooth / 27;
%
% end
%
% function Temp = calculateTemperature(square_vel, bulk_vel,i_dim, me)
%     qe = 1.6e-19;%C
%     if i_dim == 0 % Magnitude of Velocity
%         sq_vel_temp = sum(square_vel,4, 'omitnan'); % Finds magnitude of bulk vel
%         bulk_vel_temp = sqrt(sum(bulk_vel.^2,4, 'omitnan'));
%
%     else % Velocity Component
%         sq_vel_temp = square_vel(:,:,:,i_dim);
%         bulk_vel_temp = bulk_vel(:,:,:,i_dim);
%
%     end
%
%     sq_vel_temp = sq_vel_temp - (bulk_vel_temp.^2);
%
%     Temp = me * (sq_vel_temp) / qe;
%
%     if i_dim == 0
%         Temp= Temp / 3; %Taking 1/3 * Tr(Temp)
%     end
%
% end
%
% function [eDepthRadGrid, eDepthGrid, eDepthDist] = load_electron_depth_dist(runFName)
%     ncid = netcdf.open(runFName);
%     eDepth_id = netcdf.inqVarID(ncid, 'Electron Depth Dist');
%     eDepthDist = netcdf.getVar(ncid, eDepth_id);
%
%
%     eDepthRadGrid_id = netcdf.inqVarID(ncid, "elDepthRadGrid");
%     eDepthRadGrid = netcdf.getVar(ncid, eDepthRadGrid_id);
%
%     eDepthGrid_id = netcdf.inqVarID(ncid, 'elDepthGrid');
%     eDepthGrid = netcdf.getVar(ncid, eDepthGrid_id);
%
%     netcdf.close(ncid)
%
% end
%
%
% function [bulk_vel, square_vel] = load_velocity_arrays(runFName)
%     ncid = netcdf.open(runFName);
%     bulk_v_id = netcdf.inqVarID(ncid, 'Bulk Velocity');
%     bulk_vel = netcdf.getVar(ncid,bulk_v_id);
%
%     % Read square velocities for temperature
%     sq_v_id = netcdf.inqVarID(ncid, 'Square Velocity');
%     square_vel = netcdf.getVar(ncid,sq_v_id);
%     square_vel(isnan(square_vel)) = 0;
%     netcdf.close(ncid)
% end
%
%
% function [tot_plas_dens, ave_bulk_vel, ave_square_vel] = combine_run_velocities(pe_bulk_vel, pe_square_vel,plas_dens_E_PE_ref, ...
%     sw_bulk_vel, sw_square_vel, plas_dens_E_SW_ref, dE_grid, inc_gens)
%
%     % --- Calculate plas_dens for each generation -----
%     plas_dens_PE_gen = sum(plas_dens_E_PE_ref(:,:,:,:,:)...
%                     .* reshape(dE_grid(:),1,1,1,[]),4);
%     plas_dens_SW_gen = sum(plas_dens_E_SW_ref(:,:,:,:,:)...
%                     .* reshape(dE_grid(:),1,1,1,[]),4);
%
%     % --- Calculate total sw_square and bulk velocities ----
%     sq_vel_temp = permute(pe_square_vel,[2,3,4,1,5]);
%     sq_vel_PE = sum(sq_vel_temp(:,:,:,:,inc_gens) .* ...
%             plas_dens_PE_gen(:,:,:,:,inc_gens),5, 'omitnan');
%
%     bulk_vel_temp = permute(pe_bulk_vel,[2,3,4,1,5]);
%     bulk_vel_PE = sum(bulk_vel_temp(:,:,:,:,inc_gens) .*...
%             plas_dens_PE_gen(:,:,:,:,inc_gens),5, 'omitnan');
%
%     sq_vel_temp = permute(sw_square_vel,[2,3,4,1,5]);
%     sq_vel_SW = sum(sq_vel_temp(:,:,:,:,inc_gens) .* ...
%             plas_dens_SW_gen(:,:,:,:,inc_gens),5, 'omitnan');
%
%     bulk_vel_temp = permute(sw_bulk_vel,[2,3,4,1,5]);
%     bulk_vel_SW = sum(bulk_vel_temp(:,:,:,:,inc_gens) .*...
%             plas_dens_SW_gen(:,:,:,:,inc_gens),5, 'omitnan');
%
%     plas_dens_PE  = sum(plas_dens_PE_gen(:,:,:,:,inc_gens), 5, 'omitnan');
%     plas_dens_SW  = sum(plas_dens_SW_gen(:,:,:,:,inc_gens), 5, 'omitnan');
%     %---- Calculate total -------
%
%     tot_square_vel = sq_vel_PE   + sq_vel_SW;
%     tot_bulk_vel   = bulk_vel_PE + bulk_vel_SW;
%
%     tot_plas_dens  = plas_dens_PE + plas_dens_SW;
%
%     ave_square_vel = tot_square_vel ./ tot_plas_dens;
%     ave_bulk_vel   = tot_bulk_vel   ./ tot_plas_dens;
%
% end
%
% function me = get_electron_mass(runFName)
%     ncid = netcdf.open(runFName);
%     me_id = netcdf.inqVarID(ncid, 'Electron Mass');
%     me = netcdf.getVar(ncid, me_id);
%     netcdf.close(ncid)
% end

function [dtGrid, dtCount] = read_timestep_grid(runFName)
    ncid = netcdf.open(runFName);
    dtGrid_id = netcdf.inqVarID(ncid, "dT Boundaries");
    dtCount_id = netcdf.inqVarID(ncid, "dT Count");
    dtGrid = netcdf.getVar(ncid, dtGrid_id);
    dtCount = netcdf.getVar(ncid, dtCount_id);
    netcdf.close(ncid);
end

% function
