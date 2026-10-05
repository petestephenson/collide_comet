%% --------------- Save necessary Parameters from n_e v T_E plots --------

% global plas_dens_E X_mesh Y_mesh Z_mesh
% global bulk_vel grids pl_count pl_count_E plas_dens_gen
% global dE_grid
% global square_vel 



% -------------- No Collisions -------------------------------------------
% PE: 3821679.   SW: 4715285
% save('Q_1e26_no_coll', 'plas_dens_E', "plas_dens_gen", "square_vel", "dE_grid", ...
% "pl_count_E", "pl_count", "bulk_vel", "grids", "X_mesh", "Y_mesh", "Z_mesh");


% -------------- With Collisions -----------------------------------------
% PE: 4178608.   SW: 4123267
% save('Q_1e26_w_coll', 'plas_dens_E', "plas_dens_gen", "square_vel", "dE_grid", ...
% "pl_count_E", "pl_count", "bulk_vel", "grids", "X_mesh", "Y_mesh", "Z_mesh");

load("Q_1e26_no_coll.mat");
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

load("Q_1e26_w_coll.mat");

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
    
    rb_map = zeros(256, 3);
    rb_map = redblue(256);

    cb{iPlot} = colorbar(ax{iPlot});
    colormap(ax{iPlot}, rb_map);
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

rb_map = zeros(256, 3);
rb_map = redblue(256);

colormap(ax{iPlot}, rb_map);
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
%     rb_map = zeros(256, 3);
%     rb_map = redblue(256);
% 
%     cb{iPlot} = colorbar(ax{iPlot});
%     colormap(ax{iPlot}, rb_map);
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

rb_map = zeros(256, 3);
rb_map = redblue(256);

colormap(ax{iPlot}, rb_map);
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


