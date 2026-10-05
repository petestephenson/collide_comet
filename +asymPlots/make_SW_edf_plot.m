
function [f,ax] = make_SW_edf_plot(peCasesNo, swCasesNo, plasDensDir, iDim)
    f = figure(Position=[0,0,1000, 600]);
    left = [0.33, 0.61, 0.33, 0.61, 0.05];
    bottom = [0.57, 0.57, 0.1, 0.1,0.57];
    cbpos = [0.9,0.1,0.01,0.8];
    wid = 0.23;
    heig = 0.4;
    ax = cell(5,1);
    cb = cell(5,1);
    vangleList = {[1,0,0],[0,-1,0],[0,0,1]};
    slices = {[0],[0],[0]};
    yplotlim = 350;
    plot_lim=300;

    swRefNo = swCasesNo{5};
    
    swRefName = [plasDensDir, 'plas_dens', swRefNo, '.nc'];
    [qSpherical, qDayside, qUniform, qVec] = asymPlots.loadOutgassingProperties(swRefName);
    [gridsRef, ~] = asymPlots.loadRunData(swRefName);
    [dum, plas_dens_E_SW_ref] = asymPlots.loadRunData(swRefName);
%     ixGrid = zeros(3,1);
%     ixGrid(1) = find(slices{1} <= gridsRef{1}, 1, "first");
%     ixGrid(2) = find(slices{2} <= gridsRef{2}, 1, "first");
%     ixGrid(3) = find(slices{3} <= gridsRef{3}, 1, "first");
    clims = [-4,1];
    for iCase = [1,2,3,4,5]
        ax{iCase} = axes(f, 'Position',[left(iCase), bottom(iCase), wid, heig]);
        


        peAsymNo = peCasesNo{iCase};
        swAsymNo = swCasesNo{iCase};
        peAsymName = [plasDensDir, 'plas_dens', peAsymNo, '.nc'];
        swAsymName = [plasDensDir, 'plas_dens', swAsymNo, '.nc'];
        [gridsAsym, ~] = asymPlots.loadRunData(peAsymName);
%         [spatial_mesh, energy_mesh] = meshgrid(gridsAsym{iDim}(1:end-1), gridsAsym{2}(1:end-1), gridsAsym{4}(1:end-1));

        [dum, plas_dens_E_SW_Asym] = asymPlots.loadRunData(swAsymName);

        [f,ax{iCase}, cb{iCase}] = asymPlots.make_edf_plot(plas_dens_E_SW_Asym(:,:,:,:,1), gridsAsym,iDim, slices,clims, f, ax{iCase});
%         me = get_electron_mass(peAsymName);
        if iCase == 5
            qDayside = 0.0;
            qUniform = qSpherical;
        end
        % ---- Calculate Temperatures ------
        % Set xlim to -plot_lim, plot_lim
        xlim(ax{iCase}, [-plot_lim, plot_lim]);
        ylim(ax{iCase}, [0, yplotlim]);


%         me = get_electron_mass(peAsymName);
% 
%         [f,ax{iCase}, cb{iCase}] = makeTemperaturePlot(temperatureSW_Asym,X_mesh,Y_mesh,Z_mesh, slices, plot_lim,  vangle, f, ax{iCase});
        if or((iCase == 1), (iCase ==2))
            xlabel(ax{iCase}, '');
        end
%         if mod(iCase,2) == 1
        delete(cb{iCase});
%         end

    end
    cbPlot = colorbar(ax{1}, 'Position', cbpos);
    plot_str = 'log_{10}(f(E) [cm^{-3}s^{-1}eV^{-1}])';
    cbPlot.Label.String = plot_str;
    cbPlot.FontSize = 14;           
%     cbPlot.Location = 'eastoutside';
%     colormap(ax0, jet);
end