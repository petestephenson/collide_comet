

function [f,ax0,cb] = makeEIIFreqPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices,...
        plot_lim, iLogArray,  vangle, clims, varargin)
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
    
    plot_str = '\nu_{ioni,e} [s^{-1}]';
    if iLogArray
        plot_array = log10(abs(plot_array));
        plot_str = ['log_{10}(',plot_str,')'];
    end

    [f,ax0,cb] = asymPlots.makePlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle, f, ax0);
    cb.Label.String = plot_str;
    clim(ax0, clims);
%     clim(ax0, [-9,-6]);

end