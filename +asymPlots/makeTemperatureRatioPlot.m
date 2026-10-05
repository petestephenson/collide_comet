function [f,ax0, cb] = makeTemperatureRatioPlot(tempRatio,X_mesh,Y_mesh,Z_mesh, slices, plot_lim,  vangle, varargin)
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
    plot_str = 'T_{e,asym} / T_{e,sym}';
    
    clims = [0,2];
    iLogArray = 0;
    [f,ax0,cb] = asymPlots.makePlot(tempRatio,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle,f, ax0);
    cb.Label.String = plot_str;
    clim(ax0, clims);
end