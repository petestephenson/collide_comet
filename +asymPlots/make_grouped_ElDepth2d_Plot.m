function [f,ax] = make_grouped_ElDepth2d_Plot(peCasesNo, swCasesNo, plasDensDir, iAngle)

    f = figure(Position=[0,0,1000, 600]);
    left = [0.35, 0.65, 0.35, 0.65, 0.07];
    bottom = [0.57, 0.57, 0.1, 0.1,0.57];
    cbpos = [0.9,0.1,0.01,0.8];
    wid = 0.23;
    heig = 0.35;
    ax = cell(5,1);
    cb = cell(5,1);
    vangleList = {[1,0,0],[0,-1,0],[0,0,1]};
    viewNamelist = { 'YZ','XZ','XY'};
%     slices = {[0],[0],[0]};
    plot_lim = 300;
    clims = [-2,0.5];
    inc_gens = [1,2];
%     vangle = vangleList{iAngle};
    swRefNo = swCasesNo{5};
    swRefName = [plasDensDir, 'plas_dens', swRefNo, '.nc'];


%     [XS, XS_Eng, E_th] = asymPlots.loadXSections();
%     [qSpherical, qDayside, qUniform, qVec] = asymPlots.loadOutgassingProperties(swRefName);
    for iCase = [1,2,3,4,5]
        ax{iCase} = axes(f, 'Position',[left(iCase), bottom(iCase), wid, heig]);
        


        peAsymNo = peCasesNo{iCase};
        swAsymNo = swCasesNo{iCase};
        peAsymName = [plasDensDir, 'plas_dens', peAsymNo, '.nc'];
        swAsymName = [plasDensDir, 'plas_dens', swAsymNo, '.nc'];
%         [gridsAsym, plas_dens_E_Asym_pe] = asymPlots.loadRunData(peAsymName);
%         [gridsAsym, plas_dens_E_Asym_sw] = asymPlots.loadRunData(swAsymName);
%         plas_dens_E_Asym = plas_dens_E_Asym_sw + plas_dens_E_Asym_pe;
%         [X_mesh, Y_mesh, Z_mesh] = meshgrid(gridsAsym{1}(1:end-1), gridsAsym{2}(1:end-1), gridsAsym{3}(1:end-1));
%         [qSpherical, qDayside, qUniform, qVec] = asymPlots.loadOutgassingProperties(swAsymName);
%         me = asymPlots.get_electron_mass(peAsymName);
%         asymPlots.load_electron_depth_dist
        [eDepthRadGrid_sw, eDepthGrid_sw, eDepthDist_sw] = asymPlots.load_electron_depth_dist(swAsymName);
        [eDepthRadGrid_pe, eDepthGrid_pe, eDepthDist_pe] = asymPlots.load_electron_depth_dist(peAsymName);
        eDepthDist_all = eDepthDist_sw + eDepthDist_pe;
        [qSpherical, qDayside, qUniform, qVec] = asymPlots.loadOutgassingProperties(swAsymName);
        [f,ax{iCase}, cb{iCase}] = asymPlots.make_2dElDepthPlot(eDepthRadGrid_pe, eDepthGrid_pe, eDepthDist_all, iAngle, f, ax{iCase}, clims); 
        % ---- Calculate EII Freq ------
        if or((iCase == 1), (iCase ==2))
            xlabel(ax{iCase}, '');
        end

        delete(cb{iCase});

    end
    cbPlot = colorbar(ax{1}, 'Position', cbpos);
    plot_str = 'log_{10}(\tau_e [cm^{-3}s^{-1}])';
    cbPlot.Label.String = plot_str;
end