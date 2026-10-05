%%
xPicComet = {[100.45, 145, 145], [134,113,113]};
dxEBPic = {2200/287, 10.0};
fName = {'EB_fields_Q1e26.nc', 'EB_fields_210F1_trim.nc'};
oneDrive =  'C:\Users\pete4\OneDrive - Imperial College London';% Description
EngDegDir = '\Energy_Degradation';
PiCFieldDir = 'EB_fields';
%  Load Test particle model outputs and B field at location.
% ls([oneDrive, '\',EngDegDir, '\', PiCFieldDir, '\','EB_fields*.nc'])
% FilePath = [oneDrive, '\',EngDegDir, '\', PiCFieldDir, '\','EB_fields_Q1e26.nc'];
% xPicComet = ;
% 
PicData = struct();
PicData.E_Pic = cell(2,1);
PicData.B_Pic = cell(2,1);
PicData.XPicMesh = cell(2,1);
PicData.YPicMesh = cell(2,1);
PicData.ZPicMesh = cell(2,1);
PicData.vExB = cell(2,1);
for iOGR = 1:2
    FilePath = [oneDrive, '\',EngDegDir, '\', PiCFieldDir, '\',fName{iOGR}];
    ncid = netcdf.open(FilePath);
    finfo = ncinfo(FilePath);
    PicData.E_Pic{iOGR} = netcdf.getVar(ncid,0);
    PicData.B_Pic{iOGR} = netcdf.getVar(ncid,1);
    [~, nPicX] = netcdf.inqDim(ncid,0);
    [~, nPicY] = netcdf.inqDim(ncid,1);
    [~, nPicZ] = netcdf.inqDim(ncid,2);
    if size(PicData.E_Pic,1) == 3
        PicData.E_Pic{iOGR} = permute(PicData.E_Pic{iOGR}, [3,4,2,1]);
        PicData.B_Pic{iOGR} = permute(PicData.B_Pic{iOGR}, [3,4,2,1]);
    end
    x_pic = ([1:nPicX] - xPicComet{iOGR}(1)) * dxEBPic{iOGR};
    y_pic = ([1:nPicY] - xPicComet{iOGR}(2)) * dxEBPic{iOGR};
    z_pic = ([1:nPicZ] - xPicComet{iOGR}(3)) * dxEBPic{iOGR};
    PicGrids = {x_pic, y_pic, z_pic};
    [PicData.XPicMesh{iOGR}, PicData.YPicMesh{iOGR}, PicData.ZPicMesh{iOGR}] = meshgrid(x_pic, y_pic, z_pic);
    PicData.vExB{iOGR} = cross(PicData.E_Pic{iOGR}, PicData.B_Pic{iOGR}, 4)...
        ./ sum(PicData.B_Pic{iOGR}.^2, 4);
    
end
%%

% dxEBPic = 10;
% %%
% ncid = netcdf.open(FilePath);
% finfo = ncinfo(FilePath);
% % xPicComet = [149,160,160];
% % dxEBPic = 10.7;
% 
% PicData.E_Pic = netcdf.getVar(ncid,0);
% PicData.B_Pic = netcdf.getVar(ncid,1);
% [~, nPicX] = netcdf.inqDim(ncid,0);
% [~, nPicY] = netcdf.inqDim(ncid,1);
% [~, nPicZ] = netcdf.inqDim(ncid,2);
% if size(PicData.E_Pic,1) == 3
%     PicData.E_Pic = permute(PicData.E_Pic, [3,4,2,1]);
%     PicData.B_Pic = permute(PicData.B_Pic, [3,4,2,1]);
% end
% x_pic = ([1:nPicX] - xPicComet(1)) * dxEBPic;
% y_pic = ([1:nPicY] - xPicComet(2)) * dxEBPic;
% z_pic = ([1:nPicZ] - xPicComet(3)) * dxEBPic;
% PicGrids = {x_pic, y_pic, z_pic};
% [XPicMesh, YPicMesh, ZPicMesh] = meshgrid(x_pic, y_pic, z_pic);
% 



%% --------------- Save necessary Parameters from n_e v T_E plots --------

global plas_dens_E X_mesh Y_mesh Z_mesh
global bulk_vel grids pl_count pl_count_E plas_dens_gen
global dE_grid
global square_vel 



% -------------- No Collisions -------------------------------------------
% PE: 4178608.   SW: 4123267
% save('Q_1e26_no_coll', 'plas_dens_E', "plas_dens_gen", "square_vel", "dE_grid", ...
% "pl_count_E", "pl_count", "bulk_vel", "grids", "X_mesh", "Y_mesh", "Z_mesh");


% -------------- With Collisions -----------------------------------------
% PE: 4417108.   SW: 4563470
% save('Q_15e26_w_coll', 'plas_dens_E', "plas_dens_gen", "square_vel", "dE_grid", ...
% "pl_count_E", "pl_count", "bulk_vel", "grids", "X_mesh", "Y_mesh", "Z_mesh");

%%
load("Q_1e26_w_coll.mat");
runData = struct();
runData(1).plas_dens_E = plas_dens_E;
runData(1).plas_dens_E = plas_dens_E;
runData(1).plas_dens_gen = plas_dens_gen;
runData(1). square_vel = square_vel;
runData(1).dE_grid = dE_grid;
runData(1).pl_count_E = pl_count_E;
runData(1).pl_count = pl_count;
runData(1).bulk_vel = bulk_vel;
runData(1).grids = grids;
runData(1).X_mesh = X_mesh;
runData(1).Y_mesh = Y_mesh;
runData(1).Z_mesh = Z_mesh;

load("Q_15e26_w_coll.mat");

runData(2).plas_dens_E = plas_dens_E;
runData(2).plas_dens_gen = plas_dens_gen;
runData(2). square_vel = square_vel;
runData(2).dE_grid = dE_grid;
runData(2).pl_count_E = pl_count_E;
runData(2).pl_count = pl_count;
runData(2).bulk_vel = bulk_vel;
runData(2).grids = grids;
runData(2).X_mesh = X_mesh;
runData(2).Y_mesh = Y_mesh;
runData(2).Z_mesh = Z_mesh;

clear plas_dens_E plas_dens_gen square_vel dE_grid
clear pl_count_E pl_count bulk_vel grids X_mesh Y_mesh Z_mesh
%% Compute plot plas Density
% plot_plas_dens = 
inc_gens = [1,2];
for iColl = 1:2
    runData(iColl).plasDens = cell(2,1);
    for iRun = 1:2
        inc_eng = true(size(runData(iColl).dE_grid{iRun}));
        runData(iColl).plasDens{iRun} = sum(runData(iColl).plas_dens_E{iRun}(:,:,:,inc_eng, inc_gens)...
                .* reshape(runData(iColl).dE_grid{iRun}(inc_eng),1,1,1,[]),[4,5]);                    
    end
    runData(iColl).plotPlasDens = runData(iColl).plasDens{1} + runData(iColl).plasDens{2};
end

%% Compute plot Bulk velocity
inc_gens = [1,2];
for iColl = 1:2
    arrayDims = size(runData(iColl).bulk_vel{1});
    tot_bulk_vel = zeros(arrayDims([2,3,4,1]));
    tot_plas_dens = zeros(arrayDims([2,3,4]));
    bulk_vel_temp = cell(2,1);

    for iRun = 1:2
        bulk_vel_temp{iRun} = permute(runData(iColl).bulk_vel{iRun},[2,3,4,1,5]);
        bulk_vel_temp{iRun} = sum(bulk_vel_temp{iRun}(:,:,:,:,inc_gens) .*...
            runData(iColl).plas_dens_gen{iRun}(:,:,:,:,inc_gens),5, 'omitnan');
        
        plas_dens_temp = sum(runData(iColl).plas_dens_gen{iRun}(:,:,:,:,inc_gens), 5, 'omitnan');
            
%         if iSumRuns
        % Compute total for PE + SW runs 
        tot_bulk_vel  = tot_bulk_vel + bulk_vel_temp{iRun};
        tot_plas_dens = tot_plas_dens + plas_dens_temp;
%         end
        bulk_vel_temp{iRun} = bulk_vel_temp{iRun} ./ plas_dens_temp;
    end

    % Compute bulk velocity for sum of PE and SW pops
    tot_bulk_vel = tot_bulk_vel./ tot_plas_dens;
    runData(iColl).bulkVPlot = tot_bulk_vel;

end
clear plas_dens_temp bulk_vel_temp arrayDims tot_bulk_vel tot_plas_dens

%% Compute Temperatrue of the populations
inc_gens = [1,2];
iTempdim = 0;
Temp = cell(2,1);
sq_vel_temp = cell(2,1);
bulk_vel_temp = cell(2,1);
for iColl = 1:2
    arrayDims = size(runData(iColl).square_vel{1});
    tot_square_vel = zeros(arrayDims([2,3,4,1]));
    tot_bulk_vel = zeros(arrayDims([2,3,4,1]));
    tot_plas_dens = zeros(arrayDims([2,3,4]));

    for iRun = 1:2
        sq_vel_temp{iRun} = permute(runData(iColl).square_vel{iRun},[2,3,4,1,5]);
        sq_vel_temp{iRun} = sum(sq_vel_temp{iRun}(:,:,:,:,inc_gens) .* ...
            runData(iColl).plas_dens_gen{iRun}(:,:,:,:,inc_gens),5, 'omitnan');
        bulk_vel_temp{iRun} = permute(runData(iColl).bulk_vel{iRun},[2,3,4,1,5]);
        bulk_vel_temp{iRun} = sum(bulk_vel_temp{iRun}(:,:,:,:,inc_gens) .*...
            runData(iColl).plas_dens_gen{iRun}(:,:,:,:,inc_gens),5, 'omitnan');
        plas_dens_temp = sum(runData(iColl).plas_dens_gen{iRun}(:,:,:,:,inc_gens), 5, 'omitnan');

        tot_square_vel = tot_square_vel + sq_vel_temp{iRun};
        tot_bulk_vel  = tot_bulk_vel + bulk_vel_temp{iRun};
        tot_plas_dens = tot_plas_dens + plas_dens_temp;
    end

    tot_square_vel = tot_square_vel ./ tot_plas_dens;
    tot_bulk_vel = tot_bulk_vel./ tot_plas_dens;
    runData(iColl).totBulkVel = tot_bulk_vel;
    runData(iColl).totSquareVel = tot_square_vel;
end

clear tot_square_vel tot_bulk_vel tot_plas_dens bulk_vel_temp sq_vel_temp
clear arrayDims

%% Create no Coll vs Coll Plots 
% cmap = zeros(256, 3);
% cmap = redblue(256);
cmap = colorblindmap(256, 'heat');
iVeldim = 0;
iVelLog = 0;
nVarrows = 7; % no. of velocity vvectors in each dim
arrowColor = 'k';
dim_name = {'X', 'Y', 'Z'};

iTempdim = 0;

ax = cell(6,1);
cb = cell(6,1);

xLimits = [-400,400];
yLimits = [-400,400];
zLimits = [-400,400];
plasLim = [-1, 2];
% bulkVLim = [-inf, inf];
climsBulkV = [0, 2e6];
tempLim = [0, 40];
plotLim = 300;

f = figure('units', 'centimeters','Position', [2,2,30,16]); 

lef = [0.1, 0.1, 0.39, 0.39, 0.72, 0.72] - 0.04;
bot = [0.55, 0.1, 0.55, 0.1, 0.55, 0.1];
wid = zeros(6,1) + 0.28;
heig = [0.4, 0.4, 0.4, 0.4, 0.4, 0.4];
vangle = [0,0,1];
slices = {[],[],[0]};
nslice = sum(cellfun(@numel, slices));
% ------------------------------------------------------------------------
%--------------- Plasma Density -------------------------------------------
% ------------------------------------------------------------------------
iLog = 1;

for iColl = 1:2
    iPlot = iColl;
    plot_array = runData(iColl).plotPlasDens;
    if iLog 
        plot_array = log10(plot_array);
    end
    plot_array = permute(plot_array,[2,1,3]);
    ax{iPlot} = subplot('Position', [lef(iPlot), bot(iPlot), wid(iPlot), heig(iPlot)]);
    ax_slice = slice(ax{iPlot}, runData(iColl).X_mesh{1},runData(iColl).Y_mesh{1},runData(iColl).Z_mesh{1},...
                plot_array, slices{1},slices{2},slices{3});
    for i_slice = 1:nslice
        ax_slice(i_slice).EdgeColor = 'none';
    end

    xlabel(ax{iPlot}, 'X [km]');
    ylabel(ax{iPlot}, 'Y [km]');
    zlabel(ax{iPlot}, 'Z [km]');
    set(ax{iPlot}, 'FontSize', 14);
    view(ax{iPlot}, vangle);

    xlim(ax{iPlot}, xLimits);
    ylim(ax{iPlot}, yLimits);
    zlim(ax{iPlot}, zLimits);
    
    daspect(ax{iPlot}, [1,1,1])
    

    cb{iPlot} = colorbar(ax{iPlot});
    colormap(ax{iPlot}, cmap);
    caxis(ax{iPlot}, plasLim);
end

% ------------------------------------------------------------------------
% ---------------- Plot Bulk velocity ------------------------------------
% ------------------------------------------------------------------------
me = 9.11e-31; %kg
qe = 1.6e-19;%C 

for iColl = 1:2
    iPlot = iColl + 2;
    ax{iPlot} = subplot('Position', [lef(iPlot), bot(iPlot), wid(iPlot), heig(iPlot)]);
    hold(ax{iPlot}, 'on');

    %-------- Create Quiver Plot -------------------
    plot_dim = find(~cellfun(@isempty, slices),1,'first');
    % grids = {grids{iRun, 1}, grids{iRun, 2}, grids{iRun, 3}};
    ix_slice = find(slices{plot_dim} < runData(iColl).grids{iRun, plot_dim},1, 'first') -1;
    % try
    Qmesh = zeros([size(runData(iColl).X_mesh{2}),3]);
    Qmesh(:,:,:,1) = runData(iColl).X_mesh{2};
    Qmesh(:,:,:,2) = runData(iColl).Y_mesh{2};
    Qmesh(:,:,:,3) = runData(iColl).Z_mesh{2};
    Qmesh = permute(Qmesh, [2,1,3,4]);
    
    % --- need to reduce no. of arrows for visibilty ------- 
    arrow_pos = -plotLim + (2*plotLim/(nVarrows+1))*[1:nVarrows];
    arrow_select = {1:numel(runData(iColl).grids{iRun, 1})-1,...
        1:numel(runData(iColl).grids{iRun, 2})-1, 1:numel(runData(iColl).grids{iRun,3})-1};
    arrow_select{plot_dim} = ix_slice;
    for dim = 1:3
        if dim ~= plot_dim
            ix_arrows = arrayfun(@(x) find(x< runData(iColl).grids{iRun, dim},1, 'first'), arrow_pos);
            ix_arrows = unique(ix_arrows);
            arrow_select{dim} = ix_arrows;
    
        end
    end
    Qmesh = Qmesh(arrow_select{1}, arrow_select{2}, arrow_select{3},:);
%     ix_exclude_arrows = plot_cold_dens(arrow_select{1}, arrow_select{2}, arrow_select{3})<-2;
    velHat = runData(iColl).bulkVPlot;
    velHat(:,:,:,plot_dim) = 0;
    velHat = velHat ./ sqrt(sum(velHat.^2, 4));
    velHat = velHat(arrow_select{1}, arrow_select{2}, arrow_select{3},:);
%     for i_dim = 1:3 
%         temp = velHat(:,:,:, i_dim);
%         temp(ix_exclude_arrows) = NaN;
%         velHat(:,:,:,i_dim) = temp;
%     end
% end
% ---- Choose Velocity Dimension to Plot ------
% for iRun = plotRuns
if iVeldim == 0 % Magnitude of Velocity
    bulkVPlot = sqrt(sum(runData(iColl).bulkVPlot.^2,4)); % Finds magnitude of bulk vel
    varStr = ['U'];
else % Velocity Component
    bulkVPlot = runData(iColl).bulkVPlot(:,:,:,iVeldim);
    varStr = ['U_{', dim_name{iVeldim},'}'];
end
% end
unitStr = '[ms^{-1}]';
plot_array = bulkVPlot;
% --- Create Plot Array ------------------------

plot_array = permute(plot_array,[2,1,3]);
plot_str = [varStr,' ', unitStr];

nslice = sum(cellfun(@numel, slices));

if iVelLog
    plot_array = log10(abs(plot_array));
    plot_str = ['log_{10}(', plot_str, ')'];
end

% try
    plotVelHat = permute(velHat, [1,2,3,4]);
    Qmesh = permute(Qmesh, [1,2,3,4]);
% end     
% permute(cold_plas_dens,[2,1,3])
% i_exclude = permute(cold_plas_dens{iRun},[2,1,3])<10^-2;
% plot_array(i_exclude) = NaN;
ax_slice = slice(ax{iPlot}, runData(iColl).X_mesh{2},runData(iColl).Y_mesh{2},...
    runData(iColl).Z_mesh{2}, plot_array, slices{1},slices{2},slices{3});

for i_slice = 1:nslice
    ax_slice(i_slice).EdgeColor = 'none';
end

% try
    plotVelHat = reshape(plotVelHat, [], 3);            
    Qmesh = reshape(Qmesh, [], 3);

   ix_nan = isnan(plotVelHat(:,1));
    Qmesh = Qmesh(~ix_nan,:);
    plotVelHat = plotVelHat(~ix_nan,:);
    slice_dim = find(cellfun(@(x) ~isempty(x), slices),1);
    Qmesh(:, slice_dim) = slices{slice_dim};
    if size(Qmesh, 1) > 0
        quivplot = quiver3D(ax{iPlot},Qmesh, plotVelHat*plotLim/10, 'arrowColors',arrowColor,...
            'arrowRadius', 0.1);
    end
% end

xlim(ax{iPlot}, xLimits);
ylim(ax{iPlot}, yLimits);
zlim(ax{iPlot}, zLimits);
% cb2 = colorbar(ax{iPlot}, 'Position', [lef(2) + wid(2) + 0.005, bot(2), cb_wid, cb_heig]);
cb{iPlot} = colorbar(ax{iPlot});%, 'Position', [cb_lef(4), cb_bot(4), cb_wid, cb_heig]);
% cb4.Location = 'eastoutside';
% cb = colorbar(ax{iPlot});
% colormap(ax{iPlot}, jet);



colormap(ax{iPlot}, cmap);
cb{iPlot}.Label.String = plot_str;
% cb2.Location = 'eastoutside';

caxis(ax{iPlot}, climsBulkV);

view(ax{iPlot}, vangle);
daspect(ax{iPlot}, [1,1,1])

xlabel(ax{iPlot}, 'X [km]');
ylabel(ax{iPlot}, 'Y [km]');
zlabel(ax{iPlot}, 'Z [km]');
set(gca, 'Fontsize', 14);


%     xlabel(ax{iPlot}, 'X [km]');
%     ylabel(ax{iPlot}, 'Y [km]');
%     zlabel(ax{iPlot}, 'Z [km]');
%     set(ax{iPlot}, 'FontSize', 14);
%     view(ax{iPlot}, vangle);
% 
%     xlim(ax{iPlot}, xLimits);
%     ylim(ax{iPlot}, yLimits);
%     zlim(ax{iPlot}, zLimits);
%     
%     daspect(ax{iPlot}, [1,1,1])
%     
%     cmap = zeros(256, 3);
%     cmap = redblue(256);
% 
%     cb{iPlot} = colorbar(ax{iPlot});
%     colormap(ax{iPlot}, cmap);
%     caxis(ax{iPlot}, plasLim);
end

% ------------------------------------------------------------------------
% ---------------- Plot Temperature --------------------------------------
% ------------------------------------------------------------------------
for iColl = 1:2
    iPlot = iColl + 4;
    ax{iPlot} = subplot('Position', [lef(iPlot), bot(iPlot), wid(iPlot), heig(iPlot)]);
    hold(ax{iPlot}, 'on');

    
if iTempdim == 0 % Magnitude of Velocity
    sqVPlot = sqrt(sum(runData(iColl).totSquareVel.^2,4)); % Finds magnitude of Temperature
    bulkVPlot = sum(runData(iColl).totBulkVel.^2,4); % Finds magnitude of Temperature
    Temp = sqVPlot - bulkVPlot;
    varStr = ['T_e'];
else % Velocity Component
    sqVPlot = runData(iColl).totSquareVel(:,:,:,iTempdim);
    bulkVPlot = runData(iColl).totBulkVel(:,:,:,iTempdim);
    Temp = sqVPlot - bulkVPlot;
    varStr = ['T_{e,', dim_name{iTempdim},'}'];
end
Temp = me * Temp / abs(qe);
                
if iTempdim == 0 
    Temp = Temp / 3; %Taking 1/3 * Tr(Temp)
end
% end
unitStr = '[eV]';
plot_array = Temp;
% --- Create Plot Array ------------------------

plot_array = permute(plot_array,[2,1,3]);
plot_str = [varStr,' ', unitStr];

nslice = sum(cellfun(@numel, slices));

% if iVelLog
%     plot_array = log10(abs(plot_array));
%     plot_str = ['log_{10}(', plot_str, ')'];
% end

% permute(cold_plas_dens,[2,1,3])
% i_exclude = permute(cold_plas_dens{iRun},[2,1,3])<10^-2;
% plot_array(i_exclude) = NaN;
ax_slice = slice(ax{iPlot}, runData(iColl).X_mesh{2},runData(iColl).Y_mesh{2},...
    runData(iColl).Z_mesh{2}, plot_array, slices{1},slices{2},slices{3});

for i_slice = 1:nslice
    ax_slice(i_slice).EdgeColor = 'none';
end


xlim(ax{iPlot}, xLimits);
ylim(ax{iPlot}, yLimits);
zlim(ax{iPlot}, zLimits);
% cb2 = colorbar(ax{iPlot}, 'Position', [lef(2) + wid(2) + 0.005, bot(2), cb_wid, cb_heig]);
cb{iPlot} = colorbar(ax{iPlot});%, 'Position', [cb_lef(4), cb_bot(4), cb_wid, cb_heig]);
% cb4.Location = 'eastoutside';
% cb = colorbar(ax{iPlot});
% colormap(ax{iPlot}, jet);

colormap(ax{iPlot}, cmap);
cb{iPlot}.Label.String = plot_str;
% cb2.Location = 'eastoutside';

caxis(ax{iPlot}, tempLim);

view(ax{iPlot}, vangle);
daspect(ax{iPlot}, [1,1,1])

xlabel(ax{iPlot}, 'X [km]');
ylabel(ax{iPlot}, 'Y [km]');
zlabel(ax{iPlot}, 'Z [km]');
set(gca, 'Fontsize', 14);

end

%% Plot Electorn Density for populations
%% Compute plot plas Density
% plot_plas_dens = 
inc_gens = [1,2];
nPops = 3;

for iOGR = 1:2
    runData(iOGR).plotPlasDensPop = cell(nPops,1);
    runData(iOGR).plasDens = cell(2,1);
    for iRun = 1:2
        inc_eng = true(size(runData(iOGR).dE_grid{iRun}));
        runData(iOGR).plasDens{iRun} = sum(runData(iOGR).plas_dens_E{iRun}(:,:,:,inc_eng, inc_gens)...
                .* reshape(runData(iOGR).dE_grid{iRun}(inc_eng),1,1,1,[]),4);                    
    end
    runData(iOGR).plotPlasDensPop{1} = runData(iOGR).plasDens{1}(:,:,:,:,1);
    runData(iOGR).plotPlasDensPop{2} = runData(iOGR).plasDens{2}(:,:,:,:,1);
    runData(iOGR).plotPlasDensPop{3} = runData(iOGR).plasDens{1}(:,:,:,:,2) + runData(iOGR).plasDens{2}(:,:,:,:,2);
end

%% calc eXb streamlines
% iOGR = 1;
nStreamlines  = 10;
sPosX = zeros(nStreamlines,1)-200;
sPosY = zeros(nStreamlines,1);
sPosZ = linspace(-200,200, nStreamlines);
vExBstreamlines = cell(2,1);
for iOGR = 1:2
    vExBstreamlines{iOGR} = streamline(PicData.XPicMesh{iOGR},PicData.YPicMesh{iOGR},PicData.ZPicMesh{iOGR}, permute(PicData.vExB{iOGR}(:,:,:,1),[2,1,3]),...
        zeros(size(permute(PicData.vExB{iOGR}(:,:,:,1),[2,1,3]))),permute(PicData.vExB{iOGR}(:,:,:,3),[2,1,3]),...
        sPosX, sPosY, sPosZ);
end
%% Plot electron Populations
% cmap = zeros(256, 3);
% cmap = redblue(256);
cmap = colorblindmap(256, 'heat');
arrowColor = 'k';
dim_name = {'X', 'Y', 'Z'};

iTempdim = 0;

ax = cell(6,1);
cb = cell(6,1);

xLimits = [-380,380];
yLimits = [-380,380];
zLimits = [-380,380];
plasLim = [-1, 2; 0, 3.5];
plasTicks = {[0.1,0.3 1,3 10,30, 100],...
    [1,3,10,30,100,300,1000,3000]};
% bulkVLim = [-inf, inf];
climsBulkV = [0, 2e6];
tempLim = [0, 40];
plotLim = 300;

f = figure('units', 'centimeters','Position', [2,2,30,16]); 

lef = [0.1, 0.36,  0.62,0.1, 0.36, 0.62] - 0.04;
bot = [0.55, 0.55, 0.55 , 0.1, 0.1, 0.1];
wid = zeros(6,1) + 0.25;
heig = 0.4 * ones(6,1);
cblef = 0.82;
cbbot = [0.55,0.1];
cbwid = 0.02;
cbheig = heig(1);

vangle = [0,-1,0];
slices = {[],[0],[]};
nslice = sum(cellfun(@numel, slices));
% ------------------------------------------------------------------------
%--------------- Plasma Density -------------------------------------------
% ------------------------------------------------------------------------
iLog = 1;

for iPlot = 1:6
    iOGR = ceil(iPlot/3);
    iPop = mod(iPlot-1,3)+1;
    plot_array = runData(iOGR).plotPlasDensPop{iPop};
    if iLog 
        plot_array = log10(plot_array);
    end
    plot_array = permute(plot_array,[2,1,3]);
    ax{iPlot} = subplot('Position', [lef(iPlot), bot(iPlot), wid(iPlot), heig(iPlot)], 'LineWidth', 1.2);
    ax_slice = slice(ax{iPlot}, runData(iOGR).X_mesh{1},runData(iOGR).Y_mesh{1},runData(iOGR).Z_mesh{1},...
                plot_array, slices{1},slices{2},slices{3});
    for i_slice = 1:nslice
        ax_slice(i_slice).EdgeColor = 'none';
    end
%    for iOGR = 1:2
    vExBstreamlines{iOGR} = streamline(ax{iPlot},PicData.XPicMesh{iOGR},PicData.YPicMesh{iOGR},PicData.ZPicMesh{iOGR}, permute(PicData.vExB{iOGR}(:,:,:,1),[2,1,3]),...
        zeros(size(permute(PicData.vExB{iOGR}(:,:,:,1),[2,1,3]))),permute(PicData.vExB{iOGR}(:,:,:,3),[2,1,3]),...
        sPosX, sPosY, sPosZ);
    vExBstreamlines{iOGR} = streamline(ax{iPlot},PicData.XPicMesh{iOGR},PicData.YPicMesh{iOGR},PicData.ZPicMesh{iOGR}, -permute(PicData.vExB{iOGR}(:,:,:,1),[2,1,3]),...
        zeros(size(permute(PicData.vExB{iOGR}(:,:,:,1),[2,1,3]))),-permute(PicData.vExB{iOGR}(:,:,:,3),[2,1,3]),...
        sPosX, sPosY, sPosZ);
%     end
    if iOGR == 2
        xlabel(ax{iPlot}, 'X [km]');
    else
        ax{iPlot}.XTickLabel = {''};
    end
    ylabel(ax{iPlot}, 'Y [km]');
    if iPop == 1
        zlabel(ax{iPlot}, 'Z [km]');
%     else
%         ax{iPlot}.ZTickLabel = {''};
    end
    set(ax{iPlot}, 'FontSize', 14);
    view(ax{iPlot}, vangle);

    xlim(ax{iPlot}, xLimits);   
    ylim(ax{iPlot}, yLimits);
    zlim(ax{iPlot}, zLimits);
    
    daspect(ax{iPlot}, [1,1,1])
    ax{iPlot}.LineWidth = 1;
    ax{iPlot}.TickDir = 'out';
    
    if iPop ==3
        cb{iPlot} = colorbar(ax{iPlot}, 'Position', [cblef, cbbot(iOGR), cbwid, cbheig]);
%         cb{iPlot}.Label.String = 'log_{10}(n_e [cm^{-3}])';
        cb{iPlot}.Label.String = 'n_e [cm^{-3}]';
%         if iOGR == 1
        cb{iPlot}.Ticks = log10(plasTicks{iOGR});
        cb{iPlot}.TickLabels = plasTicks{iOGR};
%         end

    end
    colormap(ax{iPlot}, cmap);
    caxis(ax{iPlot}, plasLim(iOGR,:));
end
%% Plot total densities
% cmap = zeros(256, 3);
% cmap = redblue(256);
cmap = colorblindmap(256, 'heat');
arrowColor = 'k';
dim_name = {'X', 'Y', 'Z'};

iTempdim = 0;

ax = cell(6,1);
cb = cell(6,1);

xLimits = [-380,380];
yLimits = [-380,380];
zLimits = [-380,380];
plasLim = [-1, 2; 0, 3.5];
plasTicks = {[0.1,0.3 1,3 10,30, 100],...
    [1,3,10,30,100,300,1000,3000]};
% bulkVLim = [-inf, inf];
climsBulkV = [0, 2e6];
tempLim = [0, 40];
plotLim = 300;

f = figure('units', 'centimeters','Position', [2,2,20,12]); 

lef = [0.12, 0.6];%,  0.62,0.1, 0.36, 0.62] - 0.04;
bot = [0.12, 0.12];%, 0.55 , 0.1, 0.1, 0.1];
wid = zeros(6,1) + 0.35;
heig = 0.7 * ones(6,1);
cblef = [0.44,0.87];
cbbot = [0.78];
cbwid = 0.01;
cbheig = heig(1);

vangle = [0,-1,0];
slices = {[],[0],[]};
nslice = sum(cellfun(@numel, slices));
% ------------------------------------------------------------------------
%--------------- Plasma Density -------------------------------------------
% ------------------------------------------------------------------------
iLog = 1;

for iPlot = 1:2
    iOGR = iPlot;
%     iPop = mod(iPlot-1,3)+1;
    plot_array = runData(iOGR).plotPlasDensPop{1}+runData(iOGR).plotPlasDensPop{2}+runData(iOGR).plotPlasDensPop{3};
    if iLog 
        plot_array = log10(plot_array);
    end
    plot_array = permute(plot_array,[2,1,3]);
    ax{iPlot} = subplot('Position', [lef(iPlot), bot(iPlot), wid(iPlot), heig(iPlot)], 'LineWidth', 1.2);
    ax_slice = slice(ax{iPlot}, runData(iOGR).X_mesh{1},runData(iOGR).Y_mesh{1},runData(iOGR).Z_mesh{1},...
                plot_array, slices{1},slices{2},slices{3});
    for i_slice = 1:nslice
        ax_slice(i_slice).EdgeColor = 'none';
    end
%    for iOGR = 1:2
    vExBstreamlines{iOGR} = streamline(ax{iPlot},PicData.XPicMesh{iOGR},PicData.YPicMesh{iOGR},PicData.ZPicMesh{iOGR}, permute(PicData.vExB{iOGR}(:,:,:,1),[2,1,3]),...
        zeros(size(permute(PicData.vExB{iOGR}(:,:,:,1),[2,1,3]))),permute(PicData.vExB{iOGR}(:,:,:,3),[2,1,3]),...
        sPosX, sPosY, sPosZ);
    vExBstreamlines{iOGR} = streamline(ax{iPlot},PicData.XPicMesh{iOGR},PicData.YPicMesh{iOGR},PicData.ZPicMesh{iOGR}, -permute(PicData.vExB{iOGR}(:,:,:,1),[2,1,3]),...
        zeros(size(permute(PicData.vExB{iOGR}(:,:,:,1),[2,1,3]))),-permute(PicData.vExB{iOGR}(:,:,:,3),[2,1,3]),...
        sPosX, sPosY, sPosZ);
%     end
%     if iOGR == 2
    xlabel(ax{iPlot}, 'X [km]');
%     else
%         ax{iPlot}.XTickLabel = {''};
%     end
    ylabel(ax{iPlot}, 'Y [km]');
%     if iPlot == 1
    zlabel(ax{iPlot}, 'Z [km]');
%     else
%         ax{iPlot}.ZTickLabel = {''};
%     end
    set(ax{iPlot}, 'FontSize', 14);
    view(ax{iPlot}, vangle);

    xlim(ax{iPlot}, xLimits);   
    ylim(ax{iPlot}, yLimits);
    zlim(ax{iPlot}, zLimits);
    
    daspect(ax{iPlot}, [1,1,1])
    ax{iPlot}.LineWidth = 1;
    ax{iPlot}.TickDir = 'out';
    
%     if iPop ==3
%     cb{iPlot} = colorbar(ax{iPlot}, 'Position', [cblef(iOGR), cbbot(iOGR), cbwid, cbheig]);
    cb{iPlot} = colorbar(ax{iPlot}, 'northoutside');
    cb{iPlot}.Position(4)= 0.01;
    cb{iPlot}.Position(2)= cbbot(1);
%         cb{iPlot}.Label.String = 'log_{10}(n_e [cm^{-3}])';
    cb{iPlot}.Label.String = 'n_e [cm^{-3}]';
%         if iOGR == 1
    cb{iPlot}.Ticks = log10(plasTicks{iOGR});
    cb{iPlot}.TickLabels = plasTicks{iOGR};
%         end

%     end
    colormap(ax{iPlot}, cmap);
    caxis(ax{iPlot}, plasLim(iOGR,:));
end

%% plot Plasma Fraction

% cmap = zeros(256, 3);
% cmap = redblue(256);
cmap = colorblindmap(256, 'heat');
arrowColor = 'k';
dim_name = {'X', 'Y', 'Z'};

iTempdim = 0;

ax = cell(6,1);
cb = cell(6,1);

xLimits = [-380,380];
yLimits = [-380,380];
zLimits = [-380,380];
% plasLim = [-1, 2; 0, 3.5];
plasLim = [0,1;0,1];
% plasTicks = {[0.1,0.3 1,3 10,30, 100],...
%     [1,3,10,30,100,300,1000,3000]};
% bulkVLim = [-inf, inf];
climsBulkV = [0, 2e6];
tempLim = [0, 40];
plotLim = 300;

f = figure('units', 'centimeters','Position', [2,2,30,16]); 

lef = [0.1, 0.36,  0.62,0.1, 0.36, 0.62] - 0.04;
bot = [0.55, 0.55, 0.55 , 0.1, 0.1, 0.1];
wid = zeros(6,1) + 0.25;
heig = 0.4 * ones(6,1);
cblef = 0.82;
cbbot = [0.55,0.1];
cbwid = 0.02;
cbheig = heig(1);

vangle = [0,-1,0];
% vangle = [0,0,1];
slices = {[0],[0],[0]};
nslice = sum(cellfun(@numel, slices));
% ------------------------------------------------------------------------
%--------------- Plasma Density -------------------------------------------
% ------------------------------------------------------------------------
iLog = 0;

totPlotPlasDens = cell(2, 1);
for iOGR = 1:2
    totPlotPlasDens{iOGR} = runData(iOGR).plotPlasDensPop{1}+runData(iOGR).plotPlasDensPop{2}+runData(iOGR).plotPlasDensPop{3};
end
for iPlot = 1:6
    iOGR = ceil(iPlot/3);
    iPop = mod(iPlot-1,3)+1;
    plot_array = runData(iOGR).plotPlasDensPop{iPop}./totPlotPlasDens{iOGR};
    if iLog 
        plot_array = log10(plot_array);
    end
    plot_array = permute(plot_array,[2,1,3]);
    ax{iPlot} = subplot('Position', [lef(iPlot), bot(iPlot), wid(iPlot), heig(iPlot)], 'LineWidth', 1.2);
    ax_slice = slice(ax{iPlot}, runData(iOGR).X_mesh{1},runData(iOGR).Y_mesh{1},runData(iOGR).Z_mesh{1},...
                plot_array, slices{1},slices{2},slices{3});
    for i_slice = 1:nslice
        ax_slice(i_slice).EdgeColor = 'none';
    end
%    for iOGR = 1:2
    vExBstreamlines{iOGR} = streamline(ax{iPlot},PicData.XPicMesh{iOGR},PicData.YPicMesh{iOGR},PicData.ZPicMesh{iOGR}, permute(PicData.vExB{iOGR}(:,:,:,1),[2,1,3]),...
        zeros(size(permute(PicData.vExB{iOGR}(:,:,:,1),[2,1,3]))),permute(PicData.vExB{iOGR}(:,:,:,3),[2,1,3]),...
        sPosX, sPosY, sPosZ);
    vExBstreamlines{iOGR} = streamline(ax{iPlot},PicData.XPicMesh{iOGR},PicData.YPicMesh{iOGR},PicData.ZPicMesh{iOGR}, -permute(PicData.vExB{iOGR}(:,:,:,1),[2,1,3]),...
        zeros(size(permute(PicData.vExB{iOGR}(:,:,:,1),[2,1,3]))),-permute(PicData.vExB{iOGR}(:,:,:,3),[2,1,3]),...
        sPosX, sPosY, sPosZ);
%     end
    if iOGR == 2
        xlabel(ax{iPlot}, 'X [km]');
    else
        ax{iPlot}.XTickLabel = {''};
    end
    ylabel(ax{iPlot}, 'Y [km]');
    if iPop == 1
        zlabel(ax{iPlot}, 'Z [km]');
%     else
%         ax{iPlot}.ZTickLabel = {''};
    end
    set(ax{iPlot}, 'FontSize', 14);
    view(ax{iPlot}, vangle);

    xlim(ax{iPlot}, xLimits);   
    ylim(ax{iPlot}, yLimits);
    zlim(ax{iPlot}, zLimits);
    
    daspect(ax{iPlot}, [1,1,1])
    ax{iPlot}.LineWidth = 1;
    ax{iPlot}.TickDir = 'out';
    
    if iPop ==3
        cb{iPlot} = colorbar(ax{iPlot}, 'Position', [cblef, cbbot(iOGR), cbwid, cbheig]);
%         cb{iPlot}.Label.String = 'log_{10}(n_e [cm^{-3}])';
        cb{iPlot}.Label.String = 'n_e/n_{Tot}';
%         if iOGR == 1
%         cb{iPlot}.Ticks = log10(plasTicks{iOGR});
%         cb{iPlot}.TickLabels = plasTicks{iOGR};
%         end

    end
%     colormap(ax{iPlot}, cmap);
    colormap(ax{iPlot}, jet);
    caxis(ax{iPlot}, plasLim(iOGR,:));
end
%% Calculate ionization Frequency From each population
load('xsections.mat');
inc_gens = [1,2];
nPops = 3;

qe = 1.6e-19;
me = 9.11e-31;
for iOGR = 1:2
    runData(iOGR).EIoni = cell(2,1);
    runData(iOGR).plotEIoni = cell(nPops,1);

    
%     runData(iOGR).plasDens = cell(2,1);
    for iRun = 1:2
        E_grid = runData(iOGR).grids{iRun,4};
        E_mid = 0.5 * (E_grid(1:end-1) + E_grid(2:end));
        E_mid = reshape(E_mid, 1, 1, 1, []);
        ix_ioni = 7;

        XS_interped = reshape(interp1(XS_Eng, sum(XS(:,7:end),2), E_mid), ...
        1,1,1,[]);               
        elFlux = runData(iOGR).plas_dens_E{iRun}(:,:,:,:,inc_gens) .* ...
            sqrt(E_mid * 2 * qe ./ me) * 1e2; 
        runData(iOGR).EIoni{iRun} = trapz(reshape(E_mid,[],1), XS_interped .* elFlux, 4);
%         inc_eng = true(size(runData(iOGR).dE_grid{iRun}));
%         runData(iOGR).plasDens{iRun} = sum(runData(iOGR).plas_dens_E{iRun}(:,:,:,inc_eng, inc_gens)...
%                 .* reshape(runData(iOGR).dE_grid{iRun}(inc_eng),1,1,1,[]),4);                    
    end
    runData(iOGR).plotEIoni{1} = runData(iOGR).EIoni{1}(:,:,:,:,1);
    runData(iOGR).plotEIoni{2} = runData(iOGR).EIoni{2}(:,:,:,:,1);
    runData(iOGR).plotEIoni{3} = runData(iOGR).EIoni{1}(:,:,:,:,2) + runData(iOGR).EIoni{2}(:,:,:,:,2);
end

% for iRun = plotRuns
%     E_mid = 0.5 * (E_grid(1:end-1) + E_grid(2:end));
%     E_mid = reshape(E_mid, 1, 1, 1, []);
%     ix_ioni = 7;
%     
%     XS_interped = reshape(interp1(XS_Eng{iRun}, sum(XS{iRun}(:,7:end),2), E_mid), ...
%         1,1,1,[]);               
%     elFlux = sum(plas_dens_E{iRun}(:,:,:,:,inc_gens), 5) .* ...
%         sqrt(E_mid * 2 * app.qe ./ app.me(iRun)) * 1e2; 
%     eImpactIoni{iRun} = trapz(reshape(E_mid,[],1), XS_interped .* elFlux, 4);
% end
%% Plot e-impact ionization freq
% cmap = zeros(256, 3);
% cmap = redblue(256);
cmap = colorblindmap(256, 'heat');
arrowColor = 'k';
dim_name = {'X', 'Y', 'Z'};

iTempdim = 0;

ax = cell(6,1);
cb = cell(6,1);

xLimits = [-380,380];
yLimits = [-380,380];
zLimits = [-380,380];
plasLim = [-9, -6; -8, -5.5];
plasTicks = {[1e-9, 3e-9, 1e-8,3e-8, 1e-7, 3e-7, 1e-6],...
    [1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6]};
% bulkVLim = [-inf, inf];
climsBulkV = [0, 2e6];
tempLim = [0, 40];
plotLim = 300;

f = figure('units', 'centimeters','Position', [2,2,30,16]); 

lef = [0.1, 0.36,  0.62,0.1, 0.36, 0.62] - 0.04;
bot = [0.55, 0.55, 0.55 , 0.1, 0.1, 0.1];
wid = zeros(6,1) + 0.25;
heig = 0.4 * ones(6,1);
cblef = 0.82;
cbbot = [0.55,0.1];
cbwid = 0.02;
cbheig = heig(1);
vangle = [0,-1,0];
slices = {[],[0],[]};
nslice = sum(cellfun(@numel, slices));
% ------------------------------------------------------------------------
%--------------- Plasma Density -------------------------------------------
% ------------------------------------------------------------------------
iLog = 1;

for iPlot = 1:6
    iOGR = ceil(iPlot/3);
    iPop = mod(iPlot-1,3)+1;
    plot_array = runData(iOGR).plotEIoni{iPop};
    if iLog 
        plot_array = log10(plot_array);
    end
    plot_array = permute(plot_array,[2,1,3]);
    ax{iPlot} = subplot('Position', [lef(iPlot), bot(iPlot), wid(iPlot), heig(iPlot)], 'LineWidth', 1.2);
    ax_slice = slice(ax{iPlot}, runData(iOGR).X_mesh{1},runData(iOGR).Y_mesh{1},runData(iOGR).Z_mesh{1},...
                plot_array, slices{1},slices{2},slices{3});
    for i_slice = 1:nslice
        ax_slice(i_slice).EdgeColor = 'none';
    end

    if iOGR == 2
        xlabel(ax{iPlot}, 'X [km]');
    else
        ax{iPlot}.XTickLabel = {''};
    end
    ylabel(ax{iPlot}, 'Y [km]');
    if iPop == 1
        zlabel(ax{iPlot}, 'Z [km]');
%     else
%         ax{iPlot}.ZTickLabel = {''};
    end
    set(ax{iPlot}, 'FontSize', 14);
    view(ax{iPlot}, vangle);

    xlim(ax{iPlot}, xLimits);   
    ylim(ax{iPlot}, yLimits);
    zlim(ax{iPlot}, zLimits);
    
    daspect(ax{iPlot}, [1,1,1])
    ax{iPlot}.LineWidth = 1;
    ax{iPlot}.TickDir = 'out';
   vExBstreamlines{iOGR} = streamline(ax{iPlot},PicData.XPicMesh{iOGR},PicData.YPicMesh{iOGR},PicData.ZPicMesh{iOGR}, permute(PicData.vExB{iOGR}(:,:,:,1),[2,1,3]),...
        zeros(size(permute(PicData.vExB{iOGR}(:,:,:,1),[2,1,3]))),permute(PicData.vExB{iOGR}(:,:,:,3),[2,1,3]),...
        sPosX, sPosY, sPosZ);
    vExBstreamlines{iOGR} = streamline(ax{iPlot},PicData.XPicMesh{iOGR},PicData.YPicMesh{iOGR},PicData.ZPicMesh{iOGR}, -permute(PicData.vExB{iOGR}(:,:,:,1),[2,1,3]),...
        zeros(size(permute(PicData.vExB{iOGR}(:,:,:,1),[2,1,3]))),-permute(PicData.vExB{iOGR}(:,:,:,3),[2,1,3]),...
        sPosX, sPosY, sPosZ);
    if iPop ==3
        cb{iPlot} = colorbar(ax{iPlot}, 'Position', [cblef, cbbot(iOGR), cbwid, cbheig], 'FontSize',15);
%         cb{iPlot}.Label.String = 'log_{10}(\nu^{ioni}_{e} [s^{-1}])';
        cb{iPlot}.Label.String = '\nu^{ioni}_{e} [s^{-1}]';
        cb{iPlot}.Ticks = log10(plasTicks{iOGR});
        cb{iPlot}.TickLabels = plasTicks{iOGR};
    end
    colormap(ax{iPlot}, cmap);
    caxis(ax{iPlot}, plasLim(iOGR,:));
end

%% Plot e-impact ionization freq
% cmap = zeros(256, 3);
% cmap = redblue(256);
cmap = colorblindmap(256, 'heat');
arrowColor = 'k';
dim_name = {'X', 'Y', 'Z'};

iTempdim = 0;

ax = cell(6,1);
cb = cell(6,1);

xLimits = [-380,380];
yLimits = [-380,380];
zLimits = [-380,380];
plasLim = [-9, -6; -8, -5.5];
plasTicks = {[1e-9, 3e-9, 1e-8,3e-8, 1e-7, 3e-7, 1e-6],...
    [1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6]};
% bulkVLim = [-inf, inf];
climsBulkV = [0, 2e6];
tempLim = [0, 40];
plotLim = 300;

f = figure('units', 'centimeters','Position', [2,2,30,16]); 

lef = [0.1, 0.36,  0.62,0.1, 0.36, 0.62] - 0.04;
bot = [0.55, 0.55, 0.55 , 0.1, 0.1, 0.1];
wid = zeros(6,1) + 0.25;
heig = 0.4 * ones(6,1);
cblef = 0.82;
cbbot = [0.55,0.1];
cbwid = 0.02;
cbheig = heig(1);
vangle = [0,-1,0];
slices = {[],[0],[]};
nslice = sum(cellfun(@numel, slices));
% ------------------------------------------------------------------------
%--------------- Plasma Density -------------------------------------------
% ------------------------------------------------------------------------
iLog = 1;

for iPlot = 1:6
    iOGR = ceil(iPlot/3);
    iPop = mod(iPlot-1,3)+1;
    plot_array = runData(iOGR).plotEIoni{iPop};
    if iLog 
        plot_array = log10(plot_array);
    end
    plot_array = permute(plot_array,[2,1,3]);
    ax{iPlot} = subplot('Position', [lef(iPlot), bot(iPlot), wid(iPlot), heig(iPlot)], 'LineWidth', 1.2);
    ax_slice = slice(ax{iPlot}, runData(iOGR).X_mesh{1},runData(iOGR).Y_mesh{1},runData(iOGR).Z_mesh{1},...
                plot_array, slices{1},slices{2},slices{3});
    for i_slice = 1:nslice
        ax_slice(i_slice).EdgeColor = 'none';
    end

    if iOGR == 2
        xlabel(ax{iPlot}, 'X [km]');
    else
        ax{iPlot}.XTickLabel = {''};
    end
    ylabel(ax{iPlot}, 'Y [km]');
    if iPop == 1
        zlabel(ax{iPlot}, 'Z [km]');
%     else
%         ax{iPlot}.ZTickLabel = {''};
    end
    set(ax{iPlot}, 'FontSize', 14);
    view(ax{iPlot}, vangle);

    xlim(ax{iPlot}, xLimits);   
    ylim(ax{iPlot}, yLimits);
    zlim(ax{iPlot}, zLimits);
    
    daspect(ax{iPlot}, [1,1,1])
    ax{iPlot}.LineWidth = 1;
    ax{iPlot}.TickDir = 'out';
   vExBstreamlines{iOGR} = streamline(ax{iPlot},PicData.XPicMesh{iOGR},PicData.YPicMesh{iOGR},PicData.ZPicMesh{iOGR}, permute(PicData.vExB{iOGR}(:,:,:,1),[2,1,3]),...
        zeros(size(permute(PicData.vExB{iOGR}(:,:,:,1),[2,1,3]))),permute(PicData.vExB{iOGR}(:,:,:,3),[2,1,3]),...
        sPosX, sPosY, sPosZ);
    vExBstreamlines{iOGR} = streamline(ax{iPlot},PicData.XPicMesh{iOGR},PicData.YPicMesh{iOGR},PicData.ZPicMesh{iOGR}, -permute(PicData.vExB{iOGR}(:,:,:,1),[2,1,3]),...
        zeros(size(permute(PicData.vExB{iOGR}(:,:,:,1),[2,1,3]))),-permute(PicData.vExB{iOGR}(:,:,:,3),[2,1,3]),...
        sPosX, sPosY, sPosZ);
    if iPop ==3
        cb{iPlot} = colorbar(ax{iPlot}, 'Position', [cblef, cbbot(iOGR), cbwid, cbheig], 'FontSize',15);
%         cb{iPlot}.Label.String = 'log_{10}(\nu^{ioni}_{e} [s^{-1}])';
        cb{iPlot}.Label.String = '\nu^{ioni}_{e} [s^{-1}]';
        cb{iPlot}.Ticks = log10(plasTicks{iOGR});
        cb{iPlot}.TickLabels = plasTicks{iOGR};
    end
    colormap(ax{iPlot}, cmap);
    caxis(ax{iPlot}, plasLim(iOGR,:));
end
%% Save Counts Data for comp to r and B

global Xedges Yedges nCounts
global xMid yMed lQ uQ

% save('IoniFreqVsR_Q1e26', 'Xedges', 'Yedges', 'nCounts', 'xMid',...
%     'yMed', 'lQ', 'uQ')
% save('IoniFreqVsB_Q1e26', 'Xedges', 'Yedges', 'nCounts', 'xMid',...
%     'yMed', 'lQ', 'uQ')
% save('IoniFreqVsB_Q15e26', 'Xedges', 'Yedges', 'nCounts', 'xMid',...
%     'yMed', 'lQ', 'uQ')
% save('IoniFreqVsR_Q15e26', 'Xedges', 'Yedges', 'nCounts', 'xMid',...
%     'yMed', 'lQ', 'uQ')
%%
Xedges_all = cell(2,2);
Yedges_all = cell(2,2);
nCounts_all = cell(2,2);
xMid_all = cell(2,2);
yMed_all = cell(2,2);
lQ_all = cell(2,2);
uQ_all = cell(2,2);

load('IoniFreqVsR_Q1e26')
Xedges_all{1,1} = Xedges;
Yedges_all{1,1} = Yedges;
nCounts_all{1,1} = nCounts;
xMid_all{1,1} = xMid;
yMed_all{1,1} = yMed;
lQ_all{1,1} = lQ;
uQ_all{1,1} = uQ;
load('IoniFreqVsR_Q15e26')
Xedges_all{2,1} = Xedges;
Yedges_all{2,1} = Yedges;
nCounts_all{2,1} = nCounts;
xMid_all{2,1} = xMid;
yMed_all{2,1} = yMed;
lQ_all{2,1} = lQ;
uQ_all{2,1} = uQ;
load('IoniFreqVsB_Q1e26')
Xedges_all{1,2} = Xedges;
Yedges_all{1,2} = Yedges;
nCounts_all{1,2} = nCounts;
xMid_all{1,2} = xMid;
yMed_all{1,2} = yMed;
lQ_all{1,2} = lQ;
uQ_all{1,2} = uQ;
load('IoniFreqVsB_Q15e26')
Xedges_all{2,2} = Xedges;
Yedges_all{2,2} = Yedges;
nCounts_all{2,2} = nCounts;
xMid_all{2,2} = xMid;
yMed_all{2,2} = yMed;
lQ_all{2,2} = lQ;
uQ_all{2,2} = uQ;
%%
cmap = colorblindmap(256, 'heat');
iProp = 1; % vs COmetodist
f = figure('Position', [100,100,1100,450]);

lef = [0.1, 0.52];
bot = 0.15;
wid = 0.35;
heig = 0.8 ;
cblef = 0.89;
cbbot = bot;
cbwid = 0.02;
cbheig = heig(1);
vangle = [0,-1,0];

ElIoniLims = [-9, -6;-9, -5.2];
yTickPlot = {[1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6],...
    [1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6]};
slices = {[],[0],[]};
nslice = sum(cellfun(@numel, slices));
xStr = 'r [km]';
yStr = 'log_{10}(\nu^{ioni}_{e} [s^{-1}])';
cbStr = 'Fraction of Counts';
for iOGR = 1:2
    ax{iOGR,1} = axes(f, 'Position', [lef(iOGR), bot, wid(1), heig]);
    hold(ax{iOGR,1},'on');
    grid(ax{iOGR,1}, 'on')
    histogram2(ax{iOGR,1}, 'XBinEdges',Xedges_all{iOGR,1},'YBinEdges',Yedges_all{iOGR,1},...
        'BinCounts',nCounts_all{iOGR,1}, 'DisplayStyle','tile', 'HandleVisibility','off');

    ebars = errorbar(ax{iOGR,1}, xMid_all{iOGR,1},yMed_all{iOGR,1},...
        lQ_all{iOGR,1} - yMed_all{iOGR,1},uQ_all{iOGR,1} -yMed_all{iOGR,1}, 'LineStyle', 'None', ...
        'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k',...
        'LineWidth', 1.3, 'Color', 'k'); 
    
    colormap(ax{iOGR,1}, cmap)
    if iOGR == 2
        cb = colorbar(ax{iOGR,1}, 'Position', [cblef, cbbot, cbwid, cbheig]);
        cb.Label.String = 'Fraction of points';
        cb.Label.FontSize = 14;
    else
        ylabel(ax{iOGR,1}, yStr);
    end
    xlabel(ax{iOGR,1}, xStr);
    
    caxis(ax{iOGR,1}, [0,0.2]);
    set(ax{iOGR,1}, 'FontSize', 14);
    ylim(ax{iOGR,1}, ElIoniLims(iOGR,:));
    yticks(ax{iOGR,1}, log10(yTickPlot{iOGR}))
    yticklabels(ax{iOGR,1}, yTickPlot{iOGR});


end
%%
cmap = colorblindmap(256, 'heat');
iProp = 2; % vs Mag Field
f = figure('Position', [100,100,1100,450]);

lef = [0.1, 0.5];
bot = 0.15;
wid = 0.35;
heig = 0.8 ;
cblef = 0.87;
cbbot = bot;
cbwid = 0.02;
cbheig = heig(1);
vangle = [0,-1,0];

ElIoniLims = [-9, -6;-9, -5.2];
yTickPlot = {[1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6],...
    [1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6]};
slices = {[],[0],[]};
nslice = sum(cellfun(@numel, slices));
xStr = 'B [nT]';
yStr = 'log_{10}(\nu^{ioni}_{e} [s^{-1}])';
cbStr = 'Fraction of Counts';
for iOGR = 1:2
    ax{iOGR,1} = axes(f, 'Position', [lef(iOGR), bot, wid(1), heig]);
    hold(ax{iOGR,1},'on');
    grid(ax{iOGR,1}, 'on')
    histogram2(ax{iOGR,1}, 'XBinEdges',Xedges_all{iOGR,iProp},'YBinEdges',Yedges_all{iOGR,iProp},...
        'BinCounts',nCounts_all{iOGR,iProp}, 'DisplayStyle','tile', 'HandleVisibility','off');

    ebars = errorbar(ax{iOGR,1}, xMid_all{iOGR,iProp},yMed_all{iOGR,iProp},...
        lQ_all{iOGR,iProp} - yMed_all{iOGR,iProp},uQ_all{iOGR,iProp} -yMed_all{iOGR,iProp}, 'LineStyle', 'None', ...
        'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k',...
        'LineWidth', 1.3, 'Color', 'k'); 

    colormap(ax{iOGR,1}, cmap);
    if iOGR == 2
        cb = colorbar(ax{iOGR,1}, 'Position', [cblef, cbbot, cbwid, cbheig]);
        cb.Label.String = 'Fraction of points';
        cb.Label.FontSize = 14;
    else
        ylabel(ax{iOGR,1}, yStr);
    end
    xlabel(ax{iOGR,1}, xStr);
    
    caxis(ax{iOGR,1}, [0,0.2]);
    set(ax{iOGR,1}, 'FontSize', 14);
    ylim(ax{iOGR,1}, ElIoniLims(iOGR,:));
    yticks(ax{iOGR,1}, log10(yTickPlot{iOGR}))
    yticklabels(ax{iOGR,1}, yTickPlot{iOGR});
end

% %%
% oneDrive =  'C:\Users\pete4\OneDrive - Imperial College London';% Description
% EngDegDir = '\Energy_Degradation';
% PiCFieldDir = 'EB_fields';
% %  Load Test particle model outputs and B field at location.
% % ls([oneDrive, '\',EngDegDir, '\', PiCFieldDir, '\','EB_fields*.nc'])
% FilePath = [oneDrive, '\',EngDegDir, '\', PiCFieldDir, '\','EB_fields_Q1e26.nc'];
% xPicComet = [100.45, 145, 145];
% dxEBPic = 2200/287;
% % FilePath = [oneDrive, '\',EngDegDir, '\', PiCFieldDir, '\','EB_fields_210F1_trim.nc'];
% % xPicComet = [134,113,113];
% % dxEBPic = 10;
% ncid = netcdf.open(FilePath);
% finfo = ncinfo(FilePath);
% % xPicComet = [149,160,160];
% % dxEBPic = 10.7;
% global PicData
% PicData = struct();
% PicData.E_Pic = netcdf.getVar(ncid,0);
% PicData.B_Pic = netcdf.getVar(ncid,1);
% [~, nPicX] = netcdf.inqDim(ncid,0);
% [~, nPicY] = netcdf.inqDim(ncid,1);
% [~, nPicZ] = netcdf.inqDim(ncid,2);
% if size(PicData.E_Pic,1) == 3
%     PicData.E_Pic = permute(PicData.E_Pic, [3,4,2,1]);
%     PicData.B_Pic = permute(PicData.B_Pic, [3,4,2,1]);
% end
% x_pic = ([1:nPicX] - xPicComet(1)) * dxEBPic;
% y_pic = ([1:nPicY] - xPicComet(2)) * dxEBPic;
% z_pic = ([1:nPicZ] - xPicComet(3)) * dxEBPic;
% PicGrids = {x_pic, y_pic, z_pic};
% [XPicMesh, YPicMesh, ZPicMesh] = meshgrid(x_pic, y_pic, z_pic);
% PicData.XPicMesh = XPicMesh;
% PicData.YPicMesh = YPicMesh;
% PicData.ZPicMesh = ZPicMesh;
% %%
% 
% PicData.xPicComet = xPicComet;
% PicData.dxEBPic = dxEBPic;
% 
% %%
% % plotStr = 'Bmag';
% plotStr = 'r';
% % function plotElImpactvsProp(app, tile_no, plotRuns, plotStr)
% % plot histogram of e-impact ioni vs. a property of the coma
% iOGR = 1;
% x_grid = runData(iOGR).grids{1,1};
% y_grid = runData(iOGR).grids{1,2};
% z_grid = runData(iOGR).grids{1,3};
% E_grid = runData(iOGR).grids{1,4};
% % runData
% eIoni = cell(2,1);
% % calculate e-impact ioni freq.
% 
% % -------- Get midpoints of cells in grid --------------------
% xMid = 0.5 * (x_grid(2:end) + x_grid(1:end-1));
% yMid = 0.5 * (y_grid(2:end) + y_grid(1:end-1));
% zMid = 0.5 * (z_grid(2:end) + z_grid(1:end-1));
% plot_array = permute(runData(iOGR).plotEIoni{1} + runData(iOGR).plotEIoni{2} + runData(iOGR).plotEIoni{3}, [2,1,3]);
% [XMidmesh, YMidmesh, ZMidmesh] = meshgrid(xMid, yMid, zMid);
% 
% % -------- Generate Sobol Set of Random Points ---------------
% % ------- Interpolate From Grid points -----------------------
% % Histogram is populated by points randomly selected from the
% % domain using a sobol set
% P = sobolset(3);
% XSob = net(P, 1e6);
% xRange = diff(x_grid([1,end]));
% yRange = diff(y_grid([1,end]));
% zRange = diff(z_grid([1,end]));
% XSob = [x_grid(1), y_grid(1), z_grid(1)] +...
%     [xRange, yRange, zRange] .* XSob;
% eIoniSob = interp3(XMidmesh, YMidmesh, ZMidmesh,...
%     plot_array, XSob(:,1), XSob(:,2), XSob(:,3));
% tic;
% if strcmp(plotStr, 'r')
%     rPt = sqrt(sum(XSob.^2, 2));
%     xPlot = rPt;
%     xStr = 'r [km]';
% elseif strcmp(plotStr, 'Bmag')
%     % Interpolate magnetic field to position of cell centre
%     global BSob
%     interpFieldstoSobolSet(XSob)
%     % Calculate magnitude of the magnetic field
%     Bmag = sqrt(sum(BSob .^ 2, 2));
%     % Set as x Value for plots.
%     xPlot = Bmag * 1e9;
%     xStr = 'B [nT]';
% end
% toc;
% % -------- Plot points across domain -------------------------
% % t = findobj(app.RunPlotsTab, 'Type', 'tiledlayout');
% f = figure();
% ax0 = axes(f)
% hold(ax0, 'on');
% plot_lim = 500;
% 
% yStr = '\nu_{ioni}^{e} [s^{-1}]';
% 
% yPlot = eIoniSob;       
% % ----- Ensure correct shape ----------------------
% xPlot = reshape(xPlot, [], 1);
% yPlot = reshape(yPlot, [], 1);
% 
% 
% % ---- Median and Quartiles ----------------------
% nBins = 50;
% yMed = zeros(nBins, 1);
% lQ = zeros(nBins, 1);
% uQ = zeros(nBins, 1);
% 
% 
% % if app.ElImpactXlogButton.Value
% %     xPlot = log10(xPlot);
% %     xStr = ['log_{10}(', xStr, ')'];
% %     xMid = log10(xMid);
% 
% % end
% 
% % if app.ElImpactYlogButton.Value
%     yPlot = log10(yPlot);
%     yStr = ['log_{10}(', yStr, ')'];
% % end
% % Clean infinite Points 
% ixInf = or(isinf(xPlot), isinf(yPlot));
% xPlot = xPlot(~ixInf); yPlot = yPlot(~ixInf);
% % ---------- Calculate Median, Upper, Lower Quartiles -------
% 
% xEdges = linspace(0, max(xPlot, [], 'all'), nBins + 1);
% xMid  = 0.5 * (xEdges(1:end-1) + xEdges(2:end));
% ixBin = discretize(xPlot, xEdges);
% 
% for iBin = 1:nBins
%     yMed(iBin) = prctile(yPlot(ixBin == iBin), 50,  'all');
%     lQ(iBin) = prctile(yPlot(ixBin == iBin), 25, 'all');
%     uQ(iBin) = prctile(yPlot(ixBin == iBin), 75, 'all');
% end
% 
% % ------ Create Limits of Histogram and bin edges -----------
% xHistRange  = [min(xPlot, [], 'all'), max(xPlot, [], 'all')];
% yHistRange  =  [min(yPlot, [], 'all'), max(yPlot, [], 'all')];
% yLims = [-10, -5.2];
% yHistRange(~isinf(yLims)) = yLims(~isinf(yLims));
% 
% xInc = and(xPlot >= xHistRange(1), xPlot <= xHistRange(2));
% yInc =  and(yPlot >= yHistRange(1), yPlot <= yHistRange(2));
% xInc = and(xInc, yInc);
% 
% xPlot = xPlot(xInc);
% yPlot = yPlot(xInc);
% 
% 
% 
% % ---------- Create 2D Histogram ----------------------------
% 
% 
%% ------ Discretize plot ------------------------------------
% [nCounts, Xedges, Yedges] = histcounts2(xPlot,yPlot,nBins); 
% nCountsX = sum(nCounts, 2);



evsB1e26 = load('IoniFreq_testPlHistogramVsB_1e26.mat');
evsB15e26 = load('IoniFreq_testPlHistogramVsB_15e26.mat');
evsB = {evsB1e26, evsB15e26};
evsr1e26 = load('IoniFreq_testPlHistogramVsr_1e26.mat');
evsr15e26 = load('IoniFreq_testPlHistogramVsr_15e26.mat');
evsr = {evsr1e26, evsr15e26};
%%
cmap = colorblindmap(256, 'heat');
lef = [0.1, 0.52];
% bot = 0.15;
bot  = [0.16, 0.78];
wid = 0.36;
% heig = 0.8 ;
heig = [0.55, 0.15];
% wid = 0.35;
% heig = 0.8 ;
cblef = 0.89;
cbbot = bot(1);
cbwid = 0.02;
cbheig = heig(1);
vangle = [0,-1,0];

ElIoniLims = [-9, -6;-9, -5.2];
BLims = [4,14;0,30];
% yTickPlot = {[1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6],...
%     [1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6]};
slices = {[],[0],[]};
nslice = sum(cellfun(@numel, slices));
xStr = 'B [nT]';
yStr = '\nu^{ioni}_{e} [s^{-1}]';
cbStr = 'Fraction of Counts';
f = figure('Position', [100,100,1100,400]);
ax =cell(2,2);
yTickPlot = {[1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6],...
    [1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6]};
rCoeff = nan(2,1);
pFits = cell(2,1);
for iOGR = 1:2
    Xedges = evsB{iOGR}.Xedges;
    Yedges = evsB{iOGR}.Yedges;
    nCounts = evsB{iOGR}.nCounts;
    xMid= evsB{iOGR}.xMid;
    yMed = evsB{iOGR}.yMed;
    lQ = evsB{iOGR}.lQ;
    uQ = evsB{iOGR}.uQ;
    ax{iOGR,1} = axes(f, 'Position', [lef(iOGR), bot(1), wid(1), heig(1)]);
    hold(ax{iOGR,1},'on');
    grid(ax{iOGR,1}, 'on')
    nCountsX = sum(nCounts,2);
    nCounts = nCounts ./sum(nCounts,2);
    histogram2(ax{iOGR,1}, 'XBinEdges',Xedges,'YBinEdges',Yedges,...
            'BinCounts',nCounts, 'DisplayStyle','tile', 'HandleVisibility','off');
    ebars = errorbar(ax{iOGR,1}, xMid,yMed,lQ - yMed,uQ -yMed, 'LineStyle', 'None', ...
        'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k',...
        'LineWidth', 1.3, 'Color', 'k'); 
    colormap(ax{iOGR,1}, cmap);
    if iOGR == 2
        cb = colorbar(ax{iOGR,1}, 'Position', [cblef, cbbot, cbwid, cbheig]);
        cb.Label.String = 'Fraction of points';
        cb.Label.FontSize = 14;
    end
    xlabel(ax{iOGR,1}, xStr);
    if iOGR == 1
        ylabel(ax{iOGR,1}, yStr);
    end
    
    ylim(ax{iOGR,1},ElIoniLims(iOGR,:));
    xlim(ax{iOGR,1},BLims(iOGR,:))
%     caxis(ax{iOGR,1},[0, 6000]);
    caxis(ax{iOGR,1},[0, 0.2]);
    set(ax{iOGR,1}, 'FontSize', 14);

    % Fitting curve 1/r
%     if iOGR ==1
    y = yMed(xMid>6.5);
    x = xMid(xMid>6.5)';
%     else
%         y=yMed(and(xMid<400, xMid>100));
%         x = xMid(and(xMid<400, xMid>100))';
%     end
    yfit = @(b,x) b(1).*(x)+b(2) ;  % Objective Function
    CF = @(b) sum((y-yfit(b,x)).^2);        % Cost Function
    b0 = rand(2,1)*10;                      % Initial Parameter Estimates
    [B, fv] = fminsearch(CF, b0);           % Estimate Parameters7
    xFit = linspace(6.5, max(xMid));
    yFit = yfit(B,xFit);
%     rCoefTemp = corrcoef(x, y);
    rCoefTemp = corrcoef(y,yfit(B,x))
    rCoeff(iOGR) = rCoefTemp(2,1);
    pFits{iOGR} = B;
    h = plot(ax{iOGR,1}, xFit, yFit, 'x--', 'LineWidth',1.9, 'Color',gca().ColorOrder(6,:));
    h.MarkerIndices = h.MarkerIndices(1:4:end);

    ylim(ax{iOGR,1}, ElIoniLims(iOGR,:));
    yticks(ax{iOGR,1}, log10(yTickPlot{iOGR}))
    yticklabels(ax{iOGR,1}, yTickPlot{iOGR});
    ax{iOGR,2} = axes(f, "Position",[lef(iOGR), bot(2), wid, heig(2)]);
    histogram(ax{iOGR,2}, 'BinEdges', Xedges, 'BinCounts',...
            nCountsX, 'HandleVisibility', 'off');
%     xlim(ax{iOGR,2}, xLims(iOGR,:));
    %                 ylim(ax{iOGR,2}, [0, 1.2 * max(nCountsXPlot{iOGR}, [], 'all')])
    xlim(ax{iOGR,2},BLims(iOGR,:))
    if iOGR == 1 
        ylabel(ax{iOGR,2}, 'Bin Counts');
    end
    ax{iOGR,2}.XTickLabel = {''};
    set(ax{iOGR,2}, 'FontSize', fSize);
    ax{iOGR, 2}.TickDir = 'out';
    ax{iOGR, 2}.TickLength = [0.01,0.0250];
    ax{iOGR, 2}.LineWidth = 1.2;
    ax{iOGR, 2}.Box = 'off';
    ax{iOGR, 1}.Box = 'off';
    hold(ax{iOGR, 2}, 'on');
    
    grid(ax{iOGR,2}, 'on');
end
%%
intStarts = [9.75, 6.095];
intEnds = [42.1775 30.8028];
% intH = [-8.5, -7.5];
intH = [0.5e4, 3e4];
intStr = {'Early', 'Late'};
fSize =14;
cmap = colorblindmap(256, 'heat');
lef = [0.1, 0.52];
% bot = 0.15;
bot  = [0.16, 0.78];
wid = 0.35;
% heig = 0.8 ;
heig = [0.55, 0.15];
cblef = 0.89;
cbbot = bot(1);
cbwid = 0.02;
cbheig = heig(1);
vangle = [0,-1,0];
rLims = [0,600;0,800];
ElIoniLims = [-9, -6;-9, -5.2];
% yTickPlot = {[1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6],...
%     [1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6]};
slices = {[],[0],[]};
nslice = sum(cellfun(@numel, slices));
xStr = 'r [km]';
yStr = '\nu^{ioni}_{e} [s^{-1}]';
cbStr = 'Fraction of Counts';
f = figure('Position', [100,100,1100,400]);
ax =cell(2,1);
yTickPlot = {[1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6],...
    [1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6]};
rCoeff = nan(2,1);
pFits = cell(2,1);
for iOGR = 1:2
    Xedges = evsr{iOGR}.Xedges;
    Yedges = evsr{iOGR}.Yedges;
    nCounts = evsr{iOGR}.nCounts;
    xMid= evsr{iOGR}.xMid;
    yMed = evsr{iOGR}.yMed;
    lQ = evsr{iOGR}.lQ;
    uQ = evsr{iOGR}.uQ;
    ax{iOGR,1} = axes(f, 'Position', [lef(iOGR), bot(1), wid(1), heig(1)]);
    hold(ax{iOGR,1},'on');
    grid(ax{iOGR,1}, 'on')
    nCountsX = sum(nCounts,2);
    nCounts = nCounts ./sum(nCounts,2);
    histogram2(ax{iOGR,1}, 'XBinEdges',Xedges,'YBinEdges',Yedges,...
            'BinCounts',nCounts, 'DisplayStyle','tile', 'HandleVisibility','off');
    ebars = errorbar(ax{iOGR,1}, xMid,yMed,lQ - yMed,uQ -yMed, 'LineStyle', 'None', ...
        'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k',...
        'LineWidth', 1.3, 'Color', 'k'); 
    colormap(ax{iOGR,1}, cmap);
     if iOGR == 2
        cb = colorbar(ax{iOGR,1}, 'Position', [cblef, cbbot, cbwid, cbheig]);
        cb.Label.String = 'Fraction of points';
        cb.Label.FontSize = 14;
    end
    xlabel(ax{iOGR,1}, xStr);
    if iOGR == 1
        ylabel(ax{iOGR,1}, yStr);
    end
    
    ylim(ax{iOGR,1},ElIoniLims(iOGR,:));
    xlim(ax{iOGR,1},rLims(iOGR,:))
%     caxis(ax{iOGR,1},[0, 6000]);
    caxis(ax{iOGR,1},[0, 0.2]);
    set(ax{iOGR,1}, 'FontSize', 14);

    % Fitting curve 1/r
%     if iOGR ==1
%         y=yMed;
%         x = xMid';
%     else
    y=yMed(and(xMid<400, xMid>100));
    x = xMid(and(xMid<400, xMid>100))';
%     end
    yfit = @(b,x) log10(b./(x)) ;  % Objective Function
    CF = @(b) sum((y-yfit(b,x)).^2);        % Cost Function
    b0 = rand(1,1)*10;                      % Initial Parameter Estimates
    [B, fv] = fminsearch(CF, b0);           % Estimate Parameters7
%     xFit = linspace(min(xMid), max(xMid));
    xFit = linspace(100, max(xMid));
    yFit = yfit(B,xFit);
    h = plot(ax{iOGR,1}, xFit, yFit, 'x--', 'LineWidth',1.9, 'Color',gca().ColorOrder(6,:));
    h.MarkerIndices = h.MarkerIndices(1:3:end);
    ylim(ax{iOGR,1}, ElIoniLims(iOGR,:));
    yticks(ax{iOGR,1}, log10(yTickPlot{iOGR}))
    yticklabels(ax{iOGR,1}, yTickPlot{iOGR});
%     rCoefTemp = corrcoef(xMid, yMed)
    rCoefTemp = corrcoef(yMed,yfit(B,xMid));
    rCoeff(iOGR) = rCoefTemp(2,1);
    pFits{iOGR} = B;
    
    ax{iOGR,2} = axes(f, "Position",[lef(iOGR), bot(2), wid, heig(2)]);
    histogram(ax{iOGR,2}, 'BinEdges', Xedges, 'BinCounts',...
            nCountsX, 'HandleVisibility', 'off');
%     xlim(ax{iOGR,2}, xLims(iOGR,:));
    %                 ylim(ax{iOGR,2}, [0, 1.2 * max(nCountsXPlot{iOGR}, [], 'all')])
    xlim(ax{iOGR,2},rLims(iOGR,:))
    if iOGR == 1 
        ylabel(ax{iOGR,2}, 'Bin Counts');
    end
    xlim(ax{iOGR,2},rLims(iOGR,:))
    ax{iOGR,2}.XTickLabel = {''};
    set(ax{iOGR,2}, 'FontSize', fSize);
    ax{iOGR, 2}.TickDir = 'out';
    ax{iOGR, 2}.TickLength = [0.01,0.0250];
    ax{iOGR, 2}.LineWidth = 1.2;
    ax{iOGR, 2}.Box = 'off';
    ax{iOGR, 1}.Box = 'off';
    hold(ax{iOGR, 2}, 'on')
    for iInt = 1:2
        
        plot(ax{iOGR,2}, [intStarts(iInt), intEnds(iInt)], repmat(intH(iInt),2),...
                'Color', 'k', 'Marker', '|', 'MarkerSize',10, 'LineStyle', '-', 'LineWidth',1.2)
        xStrPos = mean([intStarts(iInt), intEnds(iInt)]);
%         text(ax{iOGR,2},xStrPos+15*iOGR, intH(iInt)+0.2+0.05*iOGR, intStr{iInt}, 'HorizontalAlignment','center', 'FontSize',15)
        text(ax{iOGR,2},xStrPos+15*iOGR, intH(iInt)+1.25e4, intStr{iInt}, 'HorizontalAlignment','center', 'FontSize',15)
    end
    grid(ax{iOGR,2}, 'on');
end

%%
intStarts = [9.75, 6.095];
intEnds = [42.1775 30.8028];
intH = [-8.5, -7.5];
intStr = {'Early', 'Late'};

cmap = colorblindmap(256, 'heat');
lef = [0.1, 0.52];
% bot = 0.15;
bot  = [0.16, 0.78];
wid = 0.35;
% heig = 0.8 ;
heig = [0.55, 0.15];
cblef = 0.89;
cbbot = bot(1);
cbwid = 0.02;
cbheig = heig(1);
vangle = [0,-1,0];
rLims = [0,60;0,100];
ElIoniLims = [-9, -6;-9, -5.2];
% yTickPlot = {[1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6],...
%     [1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6]};
slices = {[],[0],[]};
nslice = sum(cellfun(@numel, slices));
xStr = 'r [km]';
yStr = '\nu^{ioni}_{e} [s^{-1}]';
cbStr = 'Fraction of Counts';
f = figure('Position', [100,100,1100,400]);
ax =cell(2,1);
yTickPlot = {[1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6],...
    [1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6]};
rCoeff = nan(2,1);
pFits = cell(2,1);
for iOGR = 1:2
    Xedges = evsr{iOGR}.Xedges;
    Yedges = evsr{iOGR}.Yedges;
    nCounts = evsr{iOGR}.nCounts;
    xMid= evsr{iOGR}.xMid;
    yMed = evsr{iOGR}.yMed;
    lQ = evsr{iOGR}.lQ;
    uQ = evsr{iOGR}.uQ;
    ax{iOGR,1} = axes(f, 'Position', [lef(iOGR), bot(1), wid(1), heig(1)]);
    hold(ax{iOGR,1},'on');
    grid(ax{iOGR,1}, 'on')
    nCountsX = sum(nCounts,2);
    nCounts = nCounts ./sum(nCounts,2);
    histogram2(ax{iOGR,1}, 'XBinEdges',Xedges,'YBinEdges',Yedges,...
            'BinCounts',nCounts, 'DisplayStyle','tile', 'HandleVisibility','off');
    ebars = errorbar(ax{iOGR,1}, xMid,yMed,lQ - yMed,uQ -yMed, 'LineStyle', 'None', ...
        'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k',...
        'LineWidth', 1.3, 'Color', 'k'); 
    colormap(ax{iOGR,1}, cmap);
     if iOGR == 2
        cb = colorbar(ax{iOGR,1}, 'Position', [cblef, cbbot, cbwid, cbheig]);
        cb.Label.String = 'Fraction of points';
        cb.Label.FontSize = 14;
    end
    xlabel(ax{iOGR,1}, xStr);
    if iOGR == 1
        ylabel(ax{iOGR,1}, yStr);
    end
    
    ylim(ax{iOGR,1},ElIoniLims(iOGR,:));
    xlim(ax{iOGR,1},rLims(iOGR,:))
%     caxis(ax{iOGR,1},[0, 6000]);
    caxis(ax{iOGR,1},[0, 0.2]);
    set(ax{iOGR,1}, 'FontSize', 14);

    % Fitting curve 1/r
    if iOGR ==1
        y=yMed;
        x = xMid';
    else
        y=yMed(and(xMid<400, xMid>100));
        x = xMid(and(xMid<400, xMid>100))';
    end
    yfit = @(b,x) log10(b./(x)) ;  % Objective Function
    CF = @(b) sum((y-yfit(b,x)).^2);        % Cost Function
    b0 = rand(1,1)*10;                      % Initial Parameter Estimates
    [B, fv] = fminsearch(CF, b0);           % Estimate Parameters7
    xFit = linspace(min(xMid), max(xMid));
    yFit = yfit(B,xFit);
    h = plot(ax{iOGR,1}, xFit, yFit, 'x--', 'LineWidth',1.9, 'Color',gca().ColorOrder(6,:));
    h.MarkerIndices = h.MarkerIndices(1:3:end);

    yFit2 = polyfit(xMid(xMid<60), yMed(xMid<60),1);
    yFitTest = polyval(yFit2, xMid);
    plot(ax{iOGR,1}, xMid, yFitTest)
    ylim(ax{iOGR,1}, ElIoniLims(iOGR,:));
    yticks(ax{iOGR,1}, log10(yTickPlot{iOGR}))
    yticklabels(ax{iOGR,1}, yTickPlot{iOGR});
    rCoefTemp = corrcoef(xMid, yMed);
    rCoeff(iOGR) = rCoefTemp(2,1);
    pFits{iOGR} = B;
    for iInt = 1:2
        
        plot(ax{iOGR,1}, [intStarts(iInt), intEnds(iInt)], repmat(intH(iInt),2),...
                'Color', 'k', 'Marker', '|', 'MarkerSize',10, 'LineStyle', '-', 'LineWidth',1.2)
        xStrPos = mean([intStarts(iInt), intEnds(iInt)]);
        text(xStrPos+15*iOGR, intH(iInt)+0.2+0.05*iOGR, intStr{iInt}, 'HorizontalAlignment','center', 'FontSize',15)
    end
    ax{iOGR,2} = axes(f, "Position",[lef(iOGR), bot(2), wid, heig(2)]);
    histogram(ax{iOGR,2}, 'BinEdges', Xedges, 'BinCounts',...
            nCountsX, 'HandleVisibility', 'off');
%     xlim(ax{iOGR,2}, xLims(iOGR,:));
    %                 ylim(ax{iOGR,2}, [0, 1.2 * max(nCountsXPlot{iOGR}, [], 'all')])
    xlim(ax{iOGR,2},rLims(iOGR,:))
    if iOGR == 1 
        ylabel(ax{iOGR,2}, 'Bin Counts');
    end
    xlim(ax{iOGR,2},rLims(iOGR,:))
    ax{iOGR,2}.XTickLabel = {''};
    set(ax{iOGR,2}, 'FontSize', fSize);
    ax{iOGR, 2}.TickDir = 'out';
    ax{iOGR, 2}.TickLength = [0.01,0.0250];
    ax{iOGR, 2}.LineWidth = 1.2;
    ax{iOGR, 2}.Box = 'off';
    ax{iOGR, 1}.Box = 'off';
    grid(ax{iOGR,2}, 'on');
end
% ylim(ax0, [app.ElIoniRadialCBMin.Value, app.ElIoniRadialCBMax.Value]);
%%

% function interpFieldstoSobolSet(XSob)
%     global PicData BSob
% %     iRun = plotRuns(1);
% %     if numel(plotRuns) == 1
% %         fileName = app.EBFieldSet{iRun};
% %     else
% %         if strcmp(app.EBFieldSet{1},app.EBFieldSet{2})
% %             fileName = app.EBFieldSet{iRun};
% %         else 
% %             warning('Test Pl. Runs use different PiC fields');
% %             fileName = app.EBFieldSet{iRun};
% %         end
% %     end
%     
% 
% %     app.SelectPiCFileDropDown.Value = fileName;
% %     tic;
% %     loadPicFile(app);
% %     toc;
%     % ------ Convert Sobol set to PiC location -------------------
%     XSobPic = (XSob ./ PicData.dxEBPic) + PicData.xPicComet;
%     BSob = zeros(size(XSobPic, 1), 3);
%     % ------- Interpolate to Position of sobol set ---------------
%     for i_dim = 1:3
%         BSob(:,i_dim) = interp3(permute(PicData.B_Pic(:,:,:,i_dim),[2,1,3]), XSobPic(:,1), XSobPic(:,2), XSobPic(:,3));
%     end
% end