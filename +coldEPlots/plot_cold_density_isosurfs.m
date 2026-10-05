function [f,ax] = plot_cold_density_isosurfs(plas_dens_cold_all, X_mesh_all, Y_mesh_all, Z_mesh_all)
    f = figure('Position',[10,10,1000,500]);
    lef = [0.08,0.40,0.72];
    bot = 0.1;
    wid = 0.27;
    heig = 0.8;
    ax = cell(1,3);
    
    plot_lim = 350;
   
    plotRuns = [1,2];
    r_67P = 1.7;%km
    densLevel = [1e-2, 1e-1,1e0,1e1];

    az = -20; el = 10;
    
    ixCol = [3,2,1,4];
    nDens = numel(densLevel);
    p = cell(nDens, 3);
    for iOGR = 1:3    
        ax{iOGR} = axes(f, 'Position',[lef(iOGR), bot, wid, heig]);
        hold(ax{iOGR}, 'on');
        color_order = get(ax{iOGR},'colororder');
        plot_array = plas_dens_cold_all{1, iOGR};%+ plas_dens_cold_all{2, iOGR};
  
 
        plot_array = permute(plot_array, [2,1,3]);

        for iDens = 1:nDens
            plotLevelStr = num2str(densLevel(iDens));
            plotStr = ['N_{e} = ', plotLevelStr, 'cm^{-3}'];
    %         plotStr = ['P_{e,cold}/N_{e,cold}  = ', plotLevelStr, ' s^{-1}'];
            p{iDens, iOGR} = patch(ax{iOGR}, isosurface(X_mesh_all{plotRuns(1), iOGR},Y_mesh_all{plotRuns(1),iOGR},...
                Z_mesh_all{plotRuns(1), iOGR},...
                    plot_array,densLevel(iDens)),...
                    'DisplayName', plotStr);
                p{iDens, iOGR}.FaceVertexCData = densLevel(iDens);
            %     p{iDens, iOGR}.FaceColor = 'flat';
                iCol = ixCol(iDens);
                p{iDens, iOGR}.FaceColor = color_order(iCol,:);
                p{iDens, iOGR}.EdgeColor = color_order(iCol,:);
                p{iDens, iOGR}.FaceAlpha = 0.2;
                p{iDens, iOGR}.EdgeAlpha = 0;
        end
    %             grid on
        % ------ Plot comet nucleus for context --------
        [x_com, y_com, z_com] = sphere;
        x_com = x_com * r_67P;
        y_com = y_com * r_67P;
        z_com = z_com * r_67P;
        
        surf(ax{iOGR}, x_com, y_com, z_com,'FaceColor', 'k',...
            'Edgecolor','none', 'HandleVisibility','off');
        xlabel(ax{iOGR}, 'X [km]');
    
        if iOGR == 1
            zlabel(ax{iOGR}, 'Z [km]');
            ylabel(ax{iOGR}, 'Y [km]');
        
        end
        set(ax{iOGR}, 'FontSize', 14);
        
        xlim(ax{iOGR}, [-150,350]);
    
        ylim(ax{iOGR}, [-plot_lim,plot_lim]);
    
        zlim(ax{iOGR}, [-150,350]);
        grid(ax{iOGR}, 'on');
        
        daspect(ax{iOGR}, [1,1,1]);
        if iOGR == 1
            Leg = legend(ax{iOGR}, 'Orientation', 'horizontal', 'Position', [0.3, 0.8, 0.45, 0.08]);
        end
        view(ax{iOGR}, az, el);
    
    end


end