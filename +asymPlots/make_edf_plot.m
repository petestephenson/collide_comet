function [f,ax0, cb] = make_edf_plot(plas_dens_E, grids, iDim, slices, clims,  varargin)
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
    smoothedPlasDens_E = asymPlots.smooth_plas_dens_E(plas_dens_E);
    
    [spatialMesh, energyMesh] = meshgrid(grids{iDim}(1:end-1),grids{4}(1:end-1));
    plot_array = squeeze(log10(smoothedPlasDens_E(ixGrid{1},ixGrid{2}, ixGrid{3},:)));

    surf(ax0, spatialMesh, energyMesh,plot_array', 'EdgeColor',"none", 'FaceColor',"flat");
    ylabel(ax0,'Energy [eV]');
    xlabel(ax0, [dimName{iDim}, ' [km]']);
    view(ax0, 2);
    clim(ax0, clims);
    cb = colorbar(ax0);
    plot_str = 'log_{10}(f(E) [cm^{-3}s^{-1}eV^{-1})';
    cb.Label.String = plot_str;
    cb.FontSize = 14;           
    cb.Location = 'eastoutside';
    colormap(ax0, jet);

end