function [f,ax0,cb] = makeColdPlasDensPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle)
    
    plot_str = 'N_{e,cold} [cm^{-3}]';
    if iLogArray
        plot_array = log10(abs(plot_array));
        plot_str = ['log_{10}(',plot_str,')'];
        clims = [-2,1];

    else
        clims = [-inf, inf];
    end

    [f,ax0,cb] = makePlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle);
    cb.Label.String = plot_str;
    clim(ax0, clims);

end