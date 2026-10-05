
%% Create and save Pic Field datasets at mat files
% global PicData XPicMesh YPicMesh ZPicMesh PicVangle
% global PicGrids


%  Q = 1e26 sims: PE: 4178608 SW: 4123267
% save('Pic_Q1e26.mat', 'PicData', 'XPicMesh', 'YPicMesh', 'ZPicMesh',...
%     'PicVangle', 'PicGrids');



%  Q = 3e26 sims: PE: 5222768 SW: 5222826
% save('Pic_Q3e26.mat', 'PicData', 'XPicMesh', 'YPicMesh', 'ZPicMesh',...
%     'PicVangle', 'PicGrids');

%  Q = 1.5e27 sims: PE: 4417108 SW: 4563470
% save('Pic_Q15e26.mat', 'PicData', 'XPicMesh', 'YPicMesh', 'ZPicMesh',...
%     'PicVangle', 'PicGrids');
%%
% [PicData_all, XPicMesh_all, YPicMesh_all, ZPicMesh_all] = deal(cell(3,1));
% PicGrids_all = cell(3,3);
% OGRFile = ["Pic_Q1e26.mat", "Pic_Q3e26.mat", "Pic_Q15e26.mat"];
% for iOGR = 1:3
%     load(OGRFile(iOGR))
%     PicData_all{iOGR} = PicData; 
%     XPicMesh_all{iOGR} = XPicMesh; 
%     YPicMesh_all{iOGR} = YPicMesh; 
%     ZPicMesh_all{iOGR} = ZPicMesh; 
%     PicGrids_all(iOGR,:) = PicGrids; 
% end
% 
% save('All_PiCFields.mat', 'PicData_all', 'XPicMesh_all', 'YPicMesh_all', 'ZPicMesh_all', 'PicGrids_all', '-v7.3');
load("All_PiCFields.mat")

%% Calculate quiver plots and ExB vels
nPicArrows = 40;
plot_lim = 400;
plot_dim  = 2;
slices = {[],0,0};
vExB = cell(3,1);
velHat = cell(3,1);
Qmesh = cell(3,1);
for iOGR = 1:3
    vExB{iOGR} = cross(PicData_all{iOGR}.E_Pic, PicData_all{iOGR}.B_Pic, 4) ...
        ./ sum(PicData_all{iOGR}.B_Pic.^2, 4);

    Qmesh{iOGR} = zeros([size(XPicMesh_all{iOGR}),3]);
    Qmesh{iOGR}(:,:,:,1) = XPicMesh_all{iOGR};
    Qmesh{iOGR}(:,:,:,2) = YPicMesh_all{iOGR};
    Qmesh{iOGR}(:,:,:,3) = ZPicMesh_all{iOGR};
    Qmesh{iOGR} = permute(Qmesh{iOGR}, [2,1,3,4]);
    
    % --- need to reduce no. of arrows for visibilty -------
    ix_slice = find(slices{plot_dim} < PicGrids_all{iOGR, plot_dim}, 1, 'first') - 1;

    arrow_pos = -plot_lim + (2*plot_lim / (nPicArrows+1)) * [1:nPicArrows];
    arrow_select = {1:numel(PicGrids_all{iOGR, 1},1)-1, 1:numel(PicGrids_all{iOGR, 2}, 1)-1, 1:numel(PicGrids_all{iOGR, 3})-1};
    arrow_select{plot_dim} = ix_slice;

    
    for dim = 1:3
        if dim ~= plot_dim
            ix_arrows = arrayfun(@(x) find(x < PicGrids_all{iOGR, dim},1, 'first'), arrow_pos);
            ix_arrows = unique(ix_arrows);
            arrow_select{dim} = ix_arrows;
            
        end
    end
    Qmesh{iOGR} = Qmesh{iOGR}(arrow_select{1}, arrow_select{2}, arrow_select{3},:);
    
    velHat{iOGR} = vExB{iOGR};
    velHat{iOGR}(:,:,:,plot_dim) = 0;
    velHat{iOGR} = velHat{iOGR} ./ sqrt(sum(velHat{iOGR}.^2, 4));
    velHat{iOGR} = velHat{iOGR}(arrow_select{1}, arrow_select{2}, arrow_select{3},:);
    % ---------------- Remove arrwos from low density region -----------

end

%% Create plot
f = figure('Position',[10,10,1200,400]);
lef = [0.08,0.375,0.67];
bot = 0.15;
wid = 0.25;
heig = 0.8;
ax = cell(1,3);

plot_lim = 350;
% clims = [-inf,-inf,-inf;inf,inf,inf];
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
%     try
%         plotvelHat = permute(velHat{iOGR}, [1,2,3,4]);
%         Qmesh{iOGR} = permute(Qmesh{iOGR}, [1,2,3,4]);
%     end   
    
    
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
%     xlim(ax, [-plot_lim,plot_lim]);
%     ylim(ax, [-plot_lim,plot_lim]);
%     zlim(ax, [-plot_lim,plot_lim]);
    
    %             colormap(ax, jet);
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