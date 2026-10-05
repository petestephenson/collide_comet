function [f,ax] = make_grouped_EIIProductionRatio_Plot(peCasesNo, swCasesNo, plasDensDir, iAngle)

    f = figure(Position=[0,0,1000, 600]);
    left = [0.1,0.5, 0.1, 0.5];
    bottom = [0.6, 0.6, 0.15, 0.15];
    cbpos = [0.9,0.1,0.01,0.8];
    wid = 0.36;
    heig = 0.35;
    ax = cell(5,1);
    cb = cell(5,1);
    vangleList = {[1,0,0],[0,-1,0],[0,0,1]};
    viewNamelist = { 'YZ','XZ','XY'};
    slices = {[0],[0],[0]};
    plot_lim = 300;
    clims = [-1,1];
    inc_gens = [1,2];
    vangle = vangleList{iAngle};
    swRefNo = swCasesNo{5};
    swRefName = [plasDensDir, 'plas_dens', swRefNo, '.nc'];
    [XS, XS_Eng, E_th] = asymPlots.loadXSections();
    [qSpherical, qDayside, qUniform, qVec] = asymPlots.loadOutgassingProperties(swRefName);
    peRefNo = peCasesNo{5};
    swRefNo = swCasesNo{5};
    peRefName = [plasDensDir, 'plas_dens', peRefNo, '.nc'];
    swRefName = [plasDensDir, 'plas_dens', swRefNo, '.nc'];
    [gridsRef, plas_dens_E_Ref_pe] = asymPlots.loadRunData(peRefName);
    [gridsRef, plas_dens_E_Ref_sw] = asymPlots.loadRunData(swRefName);
    plas_dens_E_Ref = plas_dens_E_Ref_sw + plas_dens_E_Ref_pe;
%     eImpactIoniRef = asymPlots.calcEIoniFreq(plas_dens_E_Ref, gridsRef{4}, XS_Eng, XS, inc_gens);
    eiiProdRef = asymPlots.calculateEIIProduction(plas_dens_E_Ref, gridsRef, qDayside, qUniform,qVec,XS_Eng, XS,inc_gens);
    for iCase = [1,2,3,4]
        ax{iCase} = axes(f, 'Position',[left(iCase), bottom(iCase), wid, heig]);

        peAsymNo = peCasesNo{iCase};
        swAsymNo = swCasesNo{iCase};
        peAsymName = [plasDensDir, 'plas_dens', peAsymNo, '.nc'];
        swAsymName = [plasDensDir, 'plas_dens', swAsymNo, '.nc'];
        [gridsAsym, plas_dens_E_Asym_pe] = asymPlots.loadRunData(peAsymName);
        [gridsAsym, plas_dens_E_Asym_sw] = asymPlots.loadRunData(swAsymName);
        plas_dens_E_Asym = plas_dens_E_Asym_sw + plas_dens_E_Asym_pe;
        [X_mesh, Y_mesh, Z_mesh] = meshgrid(gridsAsym{1}(1:end-1), gridsAsym{2}(1:end-1), gridsAsym{3}(1:end-1));

        me = asymPlots.get_electron_mass(peAsymName);
        if iCase == 5
            qDayside = 0.0;
            qUniform = qSpherical;
        end
        % ---- Calculate EII Freq ------
%         eImpactIoniAsym = asymPlots.calcEIoniFreq(plas_dens_E_Asym, gridsAsym{4}, XS_Eng, XS, inc_gens);
        eiiProdAsym = asymPlots.calculateEIIProduction(plas_dens_E_Asym, gridsAsym, qDayside, qUniform,qVec,XS_Eng, XS,inc_gens);
        plot_array = eiiProdAsym ./ eiiProdRef;
%         [f,ax{iCase}, cb{iCase}] = asymPlots.makeTemperaturePlot(temperatureSW_Asym,X_mesh,Y_mesh,Z_mesh, slices, plot_lim,  vangle, f, ax{iCase});
        [f,ax{iCase}, cb{iCase}] = asymPlots.makeEIIproductionPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices,...
                plot_lim, true,  vangle, clims, f, ax{iCase});
        if or((iCase == 1), (iCase ==2))
            xlabel(ax{iCase}, '');
        end

        delete(cb{iCase});
    end
    cbPlot = colorbar(ax{1}, 'Position', cbpos);
    plot_str = 'log_{10}(P_{e, Asym}/P_{e, Sym})';
    cbPlot.Label.String = plot_str;
end