function [f,ax] = plot_cold_e_fraction(plas_dens_cold_all, plas_dens_all, X_mesh_all, Y_mesh_all, Z_mesh_all)
    %% ----------------- Cold electron Fraction -------------------------------
    f = figure('Position',[10,10,1000,500]);
    lef = [0.08,0.4,0.705];
    bot = 0.1;
    wid = 0.28;
    heig = 0.8;
    ax = cell(1,3);
    
    plot_lim = 350;
    clims = [0,0,0;0.05,0.1,0.3];
    % clims = [-inf,-inf,-inf;inf,inf,inf];
    iLog = 0;
    slices = {[],0,0};
    vangle = [0,-1,0];
    plotRuns = [1,2];
    for iOGR = 1:3    
        ax{iOGR} = axes(f, 'Position',[lef(iOGR), bot, wid, heig]);
    
        % plot_array = (plas_dens_cold_all{1, iOGR} + plas_dens_cold_all{2, iOGR}) ./ ...
        %     (plas_dens_all{1, iOGR} + plas_dens_all{2, iOGR});
        plot_array = (plas_dens_cold_all{1, iOGR} + plas_dens_cold_all{2, iOGR}) ./ ...
            (plas_dens_all{1, iOGR} + plas_dens_all{2, iOGR});
        % plot_array = (plas_dens_cold_all{1, iOGR}) ./ ...
        %     (plas_dens_all{1, iOGR});
        plot_array = permute(plot_array, [2,1,3]);
        plot_str = 'N_{e,cold} / N_{e}';
    %     plot_str = 'Fraction Cold';
    
        if iLog == 1
            plot_array = log10(abs(plot_array));
            plot_str = ['log_{10}(',plot_str,')'];
        end
               
                
        nslice = sum(cellfun(@numel, slices));
        
               
    %            plot_array = permute( 2 * plas_dens{1} + 1 * plas_dens{2},[2,1,3]);
        ax_slice = slice(ax{iOGR}, X_mesh_all{plotRuns(1), iOGR},Y_mesh_all{plotRuns(1), iOGR},...
            Z_mesh_all{plotRuns(1), iOGR},...
            plot_array, slices{1},slices{2},slices{3});
        view(ax{iOGR}, vangle);        
        for i_slice = 1:nslice
            ax_slice(i_slice).EdgeColor = 'none';
        end
        xlabel(ax{iOGR}, 'X [km]');
        if iOGR ==1
            ylabel(ax{iOGR}, 'Y [km]');
            zlabel(ax{iOGR}, 'Z [km]');
        else
            ax{iOGR}.ZTickLabel = {};
        end
        ax{iOGR}.LineWidth = 2;
        ax{iOGR}.TickLength = [0.02,0.02];
        ax{iOGR}.ZAxis.TickDirection = 'out';
        ax{iOGR}.XAxis.TickDirection = 'out';
        set(ax{iOGR}, 'FontSize', 14);
        
        xlim(ax{iOGR}, [-150,350]);
        ylim(ax{iOGR}, [-plot_lim,plot_lim]);
        zlim(ax{iOGR}, [-150,350]);
                
        daspect(ax{iOGR}, [1,1,1])
        
        rb_map = zeros(256, 3);
        rb_map = redblue(256);
        
        rb_map = colorblindmap(256, 'heat');
        colormap(ax{iOGR}, rb_map);
    %     colormap(ax{iOGR}, jet);
        caxis(ax{iOGR}, clims(:,iOGR));
        
    %     if iOGR == 3
            cb = colorbar(ax{iOGR});
            cb.Label.String = plot_str;
            cb.FontSize = 14;
            cb.Location = 'northoutside';
            cb.Position(2) = cb.Position(2) - 0.01;
    %         cb.Position = [0.89, bot, 0.02, heig];
    %     end
    end

    % Copy figure to clipboard
    print(f, '-clipboard', '-dbitmap');
end