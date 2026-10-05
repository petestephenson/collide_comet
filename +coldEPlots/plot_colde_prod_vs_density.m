function [f,ax] = plot_colde_prod_vs_density(plas_dens_cold_all, coldElecProd, X_mesh_all, Y_mesh_all, Z_mesh_all)
% ---------- Plot cold electron Production vs. Density -------------------

f = figure('Position',[10,10,1000,500]);
lef = [0.08,0.40,0.72];
bot = 0.1;
wid = 0.26;
heig = 0.8;
ax = cell(1,3);

plot_lim = 350;

clims = [-3,-3,-3;0,0,0];
iLog = 1;
slices = {[],0,0};
vangle = [0,-1,0];
plotRuns = [1,2];
r_67P = 1.7;%km
minVal = [1e-2,1e-2,1e-1];

iRun = 1;
for iOGR = 1:3    
    ax{iOGR} = axes(f, 'Position',[lef(iOGR), bot, wid, heig]);
    hold(ax{iOGR}, 'on');
                        
    plot_str = 'Cold Elec Production / Cold Elec Density [s^{-1}]';
%     plot_array = coldElecProd{1, iOGR} + coldElecProd{2,iOGR};
    plot_array = (coldElecProd{1, iOGR} + coldElecProd{2,iOGR})./...
        (plas_dens_cold_all{1,iOGR} + plas_dens_cold_all{2,iOGR});
    ixNan = (plas_dens_cold_all{1,iOGR} + plas_dens_cold_all{2,iOGR}) < minVal(iOGR);
    plot_array(ixNan) = NaN;
    plot_array = permute(plot_array ,[2,1,3]);
    if iLog
        plot_array = log10(abs(plot_array)); 
        plot_str = ['log_{10}(',plot_str,')'];
    end
    
    nslice = sum(cellfun(@numel, slices));
    
    ax_slice = slice(ax{iOGR}, X_mesh_all{iRun,iOGR}, Y_mesh_all{iRun,iOGR},Z_mesh_all{iRun,iOGR},...
        plot_array, slices{1},slices{2},slices{3});
    
    for i_slice = 1:nslice
        ax_slice(i_slice).EdgeColor = 'none';
    end

    xlim(ax{iOGR}, [-150,350]);
    ylim(ax{iOGR}, [-plot_lim,plot_lim]);
    zlim(ax{iOGR}, [-150,350]);
    
    rb_map = zeros(256, 3);
    rb_map = redblue(256);
%             colormap(ax{iOGR}, rb_map);
    colormap(ax{iOGR}, jet);
    
    cb = colorbar(ax{iOGR});
    if iOGR == 2
        cb.Label.String = plot_str;
    end
    cb.FontSize = 14;
    cb.Location = 'northoutside';
    cb.Position(2) = cb.Position(2) - 0.01;
    cb.Color = 'k';
         
    
    daspect(ax{iOGR}, [1,1,1])
    caxis(ax{iOGR}, clims(:,iOGR));

    view(ax{iOGR}, vangle);

    xlabel(ax{iOGR}, 'X [km]');
    ylabel(ax{iOGR}, 'Y [km]');
    if iOGR ==1
        zlabel(ax{iOGR}, 'Z [km]');
    end
    set(ax{iOGR}, 'FontSize', 14);
    ax{iOGR}.LineWidth = 2;
    ax{iOGR}.TickLength = [0.02,0.02];
    ax{iOGR}.ZAxis.TickDirection = 'out';
    ax{iOGR}.XAxis.TickDirection = 'out';
end



end