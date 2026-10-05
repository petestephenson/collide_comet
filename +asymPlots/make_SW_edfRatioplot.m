
function [f,ax] = make_SW_edfRatioplot(peCasesNo, swCasesNo, plasDensDir, iDim)
    f = figure(Position=[0,0,1000, 600]);
%     left = [0.33, 0.61, 0.33, 0.61, 0.05];
%     bottom = [0.57, 0.57, 0.1, 0.1,0.57];
    left = [0.07,0.5, 0.07, 0.5];
    bottom = [0.6, 0.6, 0.15, 0.15];
    cbpos = [0.9,0.1,0.01,0.8];
    wid = 0.36;
    heig = 0.35;
    ax = cell(5,1);
    cb = cell(5,1);
    vangleList = {[1,0,0],[0,-1,0],[0,0,1]};
    slices = {[0],[0],[0]};
    
    plot_lim=300;

    swRefNo = swCasesNo{5};
    
    swRefName = [plasDensDir, 'plas_dens', swRefNo, '.nc'];
    [qSpherical, qDayside, qUniform, qVec] = asymPlots.loadOutgassingProperties(swRefName);
    [gridsRef, ~] = asymPlots.loadRunData(swRefName);
    [dum, plas_dens_E_SW_ref] = asymPlots.loadRunData(swRefName);
%     plas_dens_E_SW_ref_smooth = asymPlots.smooth_plas_dens_E(plas_dens_E_SW_ref(:,:,:,:,1));

    clims = [-1,1];
    for iCase = [1,2,3,4]
        ax{iCase} = axes(f, 'Position',[left(iCase), bottom(iCase), wid, heig]);
        


        peAsymNo = peCasesNo{iCase};
        swAsymNo = swCasesNo{iCase};
        peAsymName = [plasDensDir, 'plas_dens', peAsymNo, '.nc'];
        swAsymName = [plasDensDir, 'plas_dens', swAsymNo, '.nc'];
        [gridsAsym, ~] = asymPlots.loadRunData(peAsymName);
%         [spatial_mesh, energy_mesh] = meshgrid(gridsAsym{iDim}(1:end-1), gridsAsym{2}(1:end-1), gridsAsym{4}(1:end-1));

        [dum, plas_dens_E_SW_Asym] = asymPlots.loadRunData(swAsymName);

        [f,ax{iCase}, cb{iCase}] = asymPlots.make_edfRatioplot(plas_dens_E_SW_ref(:,:,:,:,1),plas_dens_E_SW_Asym(:,:,:,:,1),...
            gridsAsym,iDim, slices,clims, f, ax{iCase});
%         me = get_electron_mass(peAsymName);
        if iCase == 5
            qDayside = 0.0;
            qUniform = qSpherical;
        end
       

        if or((iCase == 1), (iCase ==2))
            xlabel(ax{iCase}, '');
        end
%         if mod(iCase,2) == 1
        delete(cb{iCase});
%         end

    end
    cbPlot = colorbar(ax{1}, 'Position', cbpos);
    plot_str = 'log_{10}(f_{asym}(E)/f_{sym}(E))';
    cbPlot.Label.String = plot_str;
    cbPlot.FontSize = 14;           

end