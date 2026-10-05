function [f,ax] = plot_cold_e_production(plas_dens_cold_all, coldElecProd, vExB, velHat,Qmesh, X_mesh_all, Y_mesh_all, Z_mesh_all,...
    XPicMesh_all, YPicMesh_all, ZPicMesh_all)
f = figure('Position',[10,10,1000,500]);
lef = [0.08,0.40,0.72];
bot = 0.1;
wid = 0.26;
heig = 0.8;
ax = cell(1,3);
iRun = 1;

plot_lim = 350;
clims = [-3,-3,0;0,0,8e5];
iLog = 1;
slices = {[],0,0};
vangle = [0,-1,0];
plotRuns = [1,2];
r_67P = 1.7;%km
minVal = [1e-2,1e-2,1e-1];

  
ax{1} = axes(f, 'Position',[lef(1), bot, wid, heig]);
hold(ax{1}, 'on');
                    
plot_str = 'P_{e, cold} [cm^{-3} s^{-1}]';
plot_array = coldElecProd{1, 2} + coldElecProd{2, 2};

ixNan = (plas_dens_cold_all{1,2} + plas_dens_cold_all{2,2}) < minVal(2);
plot_array(ixNan) = NaN;
plot_array = permute(plot_array ,[2,1,3]);
iLog = 1;
if iLog
    plot_array = log10(abs(plot_array)); 
    plot_str = ['log_{10}(',plot_str,')'];
end
    
nslice = sum(cellfun(@numel, slices));

ax_slice = slice(ax{1}, X_mesh_all{iRun,2}, Y_mesh_all{iRun,2},Z_mesh_all{iRun,2},...
    plot_array, slices{1},slices{2},slices{3});

for i_slice = 1:nslice
    ax_slice(i_slice).EdgeColor = 'none';
end

xlim(ax{1}, [-150,350]);
ylim(ax{1}, [-plot_lim,plot_lim]);
zlim(ax{1}, [-150,350]);
    

rb_map = colorblindmap(256, 'heat');
colormap(ax{1}, rb_map);

cb1 = colorbar(ax{1});
% if iOGR == 2
cb1.Label.String = plot_str;
% end
cb1.FontSize = 14;
cb1.Location = 'northoutside';
cb1.Position(2) = cb1.Position(2) - 0.01;
cb1.Color = 'k';
         
 
ax{2} = axes(f, 'Position',[lef(2), bot, wid, heig]);
hold(ax{2}, 'on');

                        
plot_str = 'P_{e,cold} / n_{e,cold} [s^{-1}]';

plot_array = (coldElecProd{1, 2} + coldElecProd{2,2})./...
    (plas_dens_cold_all{1,2} + plas_dens_cold_all{2,2});
ixNan = (plas_dens_cold_all{1,2} + plas_dens_cold_all{2,2}) < minVal(2);
plot_array(ixNan) = NaN;
plot_array = permute(plot_array ,[2,1,3]);

plot_array = log10(abs(plot_array)); 
plot_str = ['log_{10}(',plot_str,')'];


nslice = sum(cellfun(@numel, slices));

ax_slice = slice(ax{2}, X_mesh_all{iRun,2}, Y_mesh_all{iRun,2},Z_mesh_all{iRun,2},...
    plot_array, slices{1},slices{2},slices{3});

for i_slice = 1:nslice
    ax_slice(i_slice).EdgeColor = 'none';
end

xlim(ax{2}, [-150,350]);
ylim(ax{2}, [-plot_lim,plot_lim]);
zlim(ax{2}, [-150,350]);

rb_map = zeros(256, 3);
rb_map = redblue(256);

rb_map = colorblindmap(256, 'heat');
colormap(ax{2}, rb_map);

cb2 = colorbar(ax{2});

cb2.Label.String = plot_str;
cb2.FontSize = 14;
cb2.Location = 'northoutside';
cb2.Position(2) = cb2.Position(2) - 0.01;
cb2.Color = 'k';
     
ax{3} = axes(f, 'Position',[lef(3), bot, wid, heig]);
hold(ax{3}, 'on');

vExBPlot = sqrt(sum(vExB{2}.^2,4)); % Finds magnitude of bulk vel
varStr = ['v_{ExB}'];

unitStr = '[ms^{-1}]';
% --- Create Plot Array ------------------------
plot_array = vExBPlot;
plot_str = [varStr,' ', unitStr];


plot_array = permute(plot_array, [2,1,3]);
plotVelHat = velHat{2};



ax_slice = slice(ax{3}, XPicMesh_all{2},YPicMesh_all{2},ZPicMesh_all{2},...
    plot_array, slices{1},slices{2},slices{3});

for i_slice = 1:nslice
    ax_slice(i_slice).EdgeColor = 'none';
end
     
%     try
ArrowLength = plot_lim / 15;
plotVelHat = reshape(plotVelHat, [], 3);            
QmeshPlot = reshape(Qmesh{2}, [], 3);

ix_nan = isnan(plotVelHat(:,1));
QmeshPlot = QmeshPlot(~ix_nan,:);
plotVelHat = plotVelHat(~ix_nan,:);
slice_dim = find(cellfun(@(x) ~isempty(x), slices),1);
QmeshPlot(:, slice_dim) = slices{slice_dim};
if size(QmeshPlot, 1) > 0

    quivplot = quiver3D(ax{3},QmeshPlot, plotVelHat * ArrowLength, 'arrowColors','k',...
        'arrowRadius', 0.06);
end
%     end
ax{3}.LineWidth = 2;


rb_map = zeros(256, 3);
rb_map = redblue(256);
rb_map = colorblindmap(256, 'heat');
colormap(ax{3}, rb_map);
% if iOGR == 3
cb3 = colorbar(ax{3});
cb3.Location = 'northoutside';
cb3.Position(2) = cb3.Position(2) - 0.01;
cb3.Color = 'k';
cb3.Label.String = 'v_{ExB} [ms^{-1}]';



for iPlot = 1:3
    xlim(ax{iPlot}, [-150,350]);
    ylim(ax{iPlot}, [-plot_lim,plot_lim]);
    zlim(ax{iPlot}, [-150,350]);
    set(ax{iPlot}, 'Fontsize', 14);
    daspect(ax{iPlot}, [1,1,1])
    caxis(ax{iPlot}, clims(:,iPlot));
    
    view(ax{iPlot}, vangle);
    
    xlabel(ax{iPlot}, 'X [km]');
    ylabel(ax{iPlot}, 'Y [km]');
    if iPlot ==1
        zlabel(ax{iPlot}, 'Z [km]');
    end
    set(ax{iPlot}, 'FontSize', 14);
    ax{iPlot}.LineWidth = 2;
    ax{iPlot}.TickLength = [0.02,0.02];
    ax{iPlot}.ZAxis.TickDirection = 'out';
    ax{iPlot}.XAxis.TickDirection = 'out';
end

end