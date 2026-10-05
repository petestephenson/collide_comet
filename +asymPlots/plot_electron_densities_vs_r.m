function [f,ax] = plot_electron_densities_vs_r(plas_dens_all,X_mesh_all, Y_mesh_all, Z_mesh_all)
    
    f = figure('Position',[10,10,1000,500]);
    caseNameList = {'Ref','-X', '+Y', '+Z', '-Z'};
    nRuns = numel(plas_dens_all);
    
    ax = axes(f);
    
    for iOGR = 1:nRuns
        plot_array = plas_dens_all{iOGR};%+ plas_dens_cold_all{2, iOGR};
        rMesh = sqrt(X_mesh_all{iOGR}.^2 + Y_mesh_all{iOGR}.^2 + Z_mesh_all{iOGR}.^2);

        scatter(ax, rMesh(:), plot_array(:), 10, 'filled', 'DisplayName', caseNameList{iOGR});
        hold(ax, 'on');
    end