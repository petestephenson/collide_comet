
function [f,ax0,cb] = makeColdDensRatioPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle)
    
    plot_str = 'N_{e,c,A}/N_{e,c,sym}';
    if iLogArray
        plot_array = log10(abs(plot_array));
        plot_str = ['log_{10}(',plot_str,')'];
        clims = [-2,2];
    else
        clims = [0,10];
    end

    [f,ax0,cb] = makePlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle);
    cb.Label.String = plot_str;
    clim(ax0, clims);
end