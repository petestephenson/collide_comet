function [f,ax0, cb] = make_edfRatioplot(plas_dens_E_ref,plas_dens_E_asym, grids, iDim, slices, clims,  varargin)
    if length(varargin) >= 1 && ~isempty(varargin{1})
        f = varargin{1};
    else
        f = figure();
    end
    
    if length(varargin) >= 2 && ~isempty(varargin{2})
        ax0 = varargin{2};
    else
        ax0 = axes(f);
    end
    dimName = {'X', 'Y', 'Z'};
    ixGrid = cell(3,1);
    ixGrid{1} = find(slices{1} <= grids{1}, 1, "first");
    ixGrid{2} = find(slices{2} <= grids{2}, 1, "first");
    ixGrid{3} = find(slices{3} <= grids{3}, 1, "first");
    ixGrid{iDim} = 1:numel(grids{iDim})-1;
%     nDims = cellfun(@numel, grids) ;
    smoothedPlasDens_E_ref = asymPlots.smooth_plas_dens_E(plas_dens_E_ref);
    smoothedPlasDens_E_asym = asymPlots.smooth_plas_dens_E(plas_dens_E_asym);
    plot_array = smoothedPlasDens_E_asym./smoothedPlasDens_E_ref;
    [spatialMesh, energyMesh] = meshgrid(grids{iDim}(1:end-1),grids{4}(1:end-1));
    plot_array = squeeze(log10(plot_array (ixGrid{1},ixGrid{2}, ixGrid{3},:)));

    surf(ax0, spatialMesh, energyMesh,plot_array', 'EdgeColor',"none", 'FaceColor',"flat");
    ylabel(ax0,'Energy [eV]');
    xlabel(ax0, [dimName{iDim}, ' [km]']);
    view(ax0, 2);
    clim(ax0, clims);
    cb = colorbar(ax0);
    plot_str = 'f_{asym}(E)/f_{sym}(E)';
    cb.Label.String = plot_str;
    cb.FontSize = 14;           
    cb.Location = 'eastoutside';
    colormap(ax0, jet);

end