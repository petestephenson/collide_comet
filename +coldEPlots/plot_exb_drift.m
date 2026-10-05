function plot_exb_drift(vExB, velHat, Qmesh, XPicMesh_all, YPicMesh_all, ZPicMesh_all)
f = figure('Position',[10,10,1200,400]);
lef = [0.08,0.375,0.67];
bot = 0.15;
wid = 0.25;
heig = 0.8;
ax = cell(1,3);

plot_lim = 350;

clims = [0,0,0;8e5,8e5,8e5];
iLog = 0;
i_dim = 0;
slices = {[],0,0};
PicVangle = [0,-1,0];
nslice = sum(cellfun(@numel, slices));

for iOGR = 1:3

    ax{iOGR} = axes(f, 'Position',[lef(iOGR), bot, wid, heig]);
    hold(ax{iOGR}, 'on');
    
    dim_name = {'X','Y','Z'};

    
    if i_dim == 0 % Magnitude of Velocity
        vExBPlot = sqrt(sum(vExB{iOGR}.^2,4)); % Finds magnitude of bulk vel
        varStr = ['v_{ExB}'];
    else % Velocity Component
        vExBPlot = vExB{iOGR}(:,:,:,i_dim);
        varStr = ['v_{ExB,', dim_name{i_dim},'}'];
    end
    unitStr = '[ms^{-1}]';
    % --- Create Plot Array ------------------------
    plot_array = vExBPlot;
    plot_str = [varStr,' ', unitStr];
    
    if iLog
        plot_array = log10(abs(plot_array));
        plot_str = ['log_{10}(', plot_str, ')'];
    end

    plot_array = permute(plot_array, [2,1,3]);
    plotVelHat = velHat{iOGR};

    
    ax_slice = slice(ax{iOGR}, XPicMesh_all{iOGR},YPicMesh_all{iOGR},ZPicMesh_all{iOGR},...
        plot_array, slices{1},slices{2},slices{3});
    
    for i_slice = 1:nslice
        ax_slice(i_slice).EdgeColor = 'none';
    end
         
%     try
    ArrowLength = plot_lim / 15;
    plotVelHat = reshape(plotVelHat, [], 3);            
    QmeshPlot = reshape(Qmesh{iOGR}, [], 3);
    
    ix_nan = isnan(plotVelHat(:,1));
    QmeshPlot = QmeshPlot(~ix_nan,:);
    plotVelHat = plotVelHat(~ix_nan,:);
    slice_dim = find(cellfun(@(x) ~isempty(x), slices),1);
    QmeshPlot(:, slice_dim) = slices{slice_dim};
    if size(QmeshPlot, 1) > 0

        quivplot = quiver3D(ax{iOGR},QmeshPlot, plotVelHat * ArrowLength, 'arrowColors','k',...
            'arrowRadius', 0.06);
    end
%     end
    ax{iOGR}.LineWidth = 2;
    ax{iOGR}.TickLength = [0.02,0.02];
    ax{iOGR}.ZAxis.TickDirection = 'out';
    ax{iOGR}.XAxis.TickDirection = 'out';
    

    xlim(ax{iOGR}, [-150,350]);
    ylim(ax{iOGR}, [-plot_lim,plot_lim]);
    zlim(ax{iOGR}, [-150,350]);

    xlabel(ax{iOGR}, 'X [km]');
    if iOGR == 1
        ylabel(ax{iOGR}, 'Y [km]');
        zlabel(ax{iOGR}, 'Z [km]');
    end

    rb_map = zeros(256, 3);
    rb_map = redblue(256);
    colormap(ax{iOGR}, rb_map);
    if iOGR == 3
        cb = colorbar(ax{iOGR});
        cb.Label.String = plot_str;
        cb.Position = [0.93, bot - 0.02, 0.02, heig];
    end

    caxis(ax{iOGR}, clims(:,iOGR));      

    
    view(ax{iOGR}, PicVangle);
    daspect(ax{iOGR}, [1,1,1])

    set(ax{iOGR}, 'Fontsize', 14);
end

end