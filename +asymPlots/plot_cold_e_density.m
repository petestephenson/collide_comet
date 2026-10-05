function [f,ax] = plot_cold_e_density(plas_dens_cold_all,X_mesh_all, Y_mesh_all, Z_mesh_all)
    f = figure('Position',[10,10,1200,400]);
    nRuns = numel(plas_dens_cold_all);
%     lef = [0.08,0.35,0.62];
    lef = linspace(0.08, 1-0.5*((2+0.5)/nRuns ), nRuns);
    caseNameList = {'Ref','-X', '+Y', '+Z', '-Z'};
    bot = 0.15;
    wid = 1/(nRuns+2);
    heig = 0.8;
    ax = cell(1,3);
    
    plot_lim = 350;
    clims = zeros(2, nRuns);
    clims(1,:) = -3;
    clims(2,:) = 1;

    % clims = [0,0,0;0.05,10,35];
    iLog = 1;
    slices = {[],0,0};
    vangle = [0,-1,0];
    plotRuns = [1];
    
    for iOGR = 1:nRuns    
        ax{iOGR} = axes(f, 'Position',[lef(iOGR), bot, wid, heig]);
    
        plot_array = plas_dens_cold_all{iOGR};%+ plas_dens_cold_all{2, iOGR};
        plot_array = permute(plot_array, [2,1,3]);
        plot_str = 'N_e [cm^{-3}]';
    
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
                
        daspect(ax{iOGR}, [1,1,1]);
        
        rb_map = colorblindmap(256, 'heat');
        
        
        colormap(ax{iOGR}, rb_map);
        caxis(ax{iOGR}, clims(:,iOGR));
        title(ax{iOGR}, caseNameList{iOGR});
        if iOGR == 3
            cb = colorbar(ax{iOGR});
            cb.Label.String = plot_str;
            cb.FontSize = 14;
    %         cb.Location = 'northoutside';
            cb.Position = [0.89, bot, 0.02, heig];
        end
    end
    % Copy figure to clipboard
    sgtitle('Cold e Density')
    print(f, '-clipboard', '-dbitmap');
end