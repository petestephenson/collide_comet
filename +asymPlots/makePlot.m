
function [f, ax0, cb] = makePlot(plot_array, X_mesh, Y_mesh, Z_mesh, slices, plot_lim, iLogArray, vangle, varargin)
    % Parse optional inputs
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

    hold(ax0, 'on');

            
    % --- Create Plot Array ------------------------

    plot_array = permute(plot_array,[2,1,3]);

    % set viewing angle
    view(ax0, vangle);
   
    % Plot slices of plasma density
    nslice = sum(cellfun(@numel, slices));
    ax_slice = slice(ax0, X_mesh,Y_mesh,Z_mesh,...
        plot_array, slices{1},slices{2},slices{3});
    
    for i_slice = 1:nslice
        ax_slice(i_slice).EdgeColor = 'none';
    end

    % Axes labels
    xlabel(ax0, 'X [km]');
    ylabel(ax0, 'Y [km]');
    zlabel(ax0, 'Z [km]');
    set(ax0, 'FontSize', 14);
    
    % Axis limits
    xlim(ax0, [-plot_lim,plot_lim]);
    ylim(ax0, [-plot_lim,plot_lim]);
    zlim(ax0, [-plot_lim,plot_lim]);
    
    daspect(ax0, [1,1,1])
    % Set colormap
    rb_map = zeros(256, 3);
    rb_map = redblue(256);

    customMap = [255,255,255;...
        2,5,95;...
        2,241,229;...
        23,142,5;...
        251,246,13;...
        197,62,6;...
        105,11,13;...
        212, 72, 73;...
        0,0,0];
    customMap = customMap / 255;
    cMap = zeros(256,3);
                
    nStep = round(1+ ([0:size(customMap,1)-1] * (size(cMap,1)-1) / (size(customMap,1)-1)));
    
    cMap(nStep, :) = customMap;
    
    for iStep = 1:numel(nStep)-1
        cMap(nStep(iStep):nStep(iStep+1),:) = interp1(nStep(iStep:iStep+1), customMap(iStep:iStep + 1,:), nStep(iStep):nStep(iStep+1));
    
    end
    cb = colorbar(ax0);
    colormap(ax0, rb_map);

    cb.Location = 'eastoutside'; 
end