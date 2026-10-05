% global X_mesh Y_mesh Z_mesh
% global EHall EAmbi vPot
% 
% 
% 
% save('fieldsQ15e26_noColl', 'X_mesh', 'Y_mesh', ...
%     'Z_mesh', 'EHall', 'EAmbi', 'vPot');
% save('fieldsQ15e26_wColl', 'X_mesh', 'Y_mesh', ...
%     'Z_mesh', 'EHall', 'EAmbi', 'vPot');

%% --------- combine separate properties into 1 set of cells --------------
X_mesh_all = cell(2,1);
Y_mesh_all = cell(2,1);
Z_mesh_all = cell(2,1);
EHall_all = cell(2,1);
EAmbi_all = cell(2,1);
vPot_all = cell(2,1);

load('fieldsQ15e26_noColl.mat');
X_mesh_all{1} = X_mesh{1};
Y_mesh_all{1} = Y_mesh{1};
Z_mesh_all{1} = Z_mesh{1};
EHall_all{1} = EHall;
EAmbi_all{1} = EAmbi;
vPot_all{1} = vPot;

load('fieldsQ15e26_wColl.mat');
X_mesh_all{2} = X_mesh{1};
Y_mesh_all{2} = Y_mesh{1};
Z_mesh_all{2} = Z_mesh{1};
EHall_all{2} = EHall;
EAmbi_all{2} = EAmbi;
vPot_all{2} = vPot;

%% --------- Plot Ambipolar Fields ----------------------------------------
% Plot coll / non-coll
cmap = colorblindmap(256, 'heat');
i_dim = 2;
iLog = 0;
smallPot = 2e-4;
% vangle = [1,0,0]; %YZ
% vangle = [0,-1,0]; %XZ
vangle = [0,0,1]; %XY
clims = [0, 1];
% clims = [0, 0.5];
slices = {0,0,0};
plot_array = cell(2,1);
plotRuns = [1,2];
plot_lim = 400;%km

dimName = {'X', 'Y', 'Z'};

if i_dim == 0
    plot_str = 'E_{Ambi,coll} / E_{Ambi, no coll}';
else
    plot_str = ['E_{Ambi,', dimName{i_dim},',coll} / E_{Ambi,',...
        dimName{i_dim},', no coll}'];
end


for iRun = 1:2
    if i_dim == 0
        plot_array{iRun} = sqrt(sum(EAmbi_all{iRun}.^2,4));
%         dimName = '';
        
    else
        plot_array{iRun} = EAmbi_all{iRun}(:,:,:,i_dim);
%         dimName = dimName{i_dim};
    end
end
iSmallPot = plot_array{1} < smallPot;
plot_array = plot_array{2} ./ plot_array{1};
if iLog 
    plot_array = log10(abs(plot_array));
    plot_str = ['log_{10}(',plot_str,')'];
end

% plot_array(iSmallPot) = NaN;
nslice = sum(cellfun(@numel, slices));
f = figure;
ax1 = axes(f);

ax_slice = slice(ax1, X_mesh_all{plotRuns(1)},Y_mesh_all{plotRuns(1)},...
    Z_mesh_all{plotRuns(1)},...
    plot_array, slices{1},slices{2},slices{3});
            
for i_slice = 1:nslice
    ax_slice(i_slice).EdgeColor = 'none';
end
xlabel(ax1, 'X [km]');
ylabel(ax1, 'Y [km]');
zlabel(ax1, 'Z [km]');
set(ax1, 'FontSize', 14);

xlim(ax1, [-plot_lim,plot_lim]);
ylim(ax1, [-plot_lim,plot_lim]);
zlim(ax1, [-plot_lim,plot_lim]);

daspect(ax1, [1,1,1])
cb = colorbar(ax1);
colormap(ax1, cmap);
%             colormap(ax0, rb_map);
%             colormap(ax0, cMap);
cb.Label.String = plot_str;
cb.FontSize = 14;

caxis(ax1, clims);
cb.Location = 'eastoutside';
view(ax1, vangle);

r_pt = sqrt(X_mesh_all{1}.^2 + Y_mesh_all{1}.^2 + Z_mesh_all{1}.^2);

%% --------- Plot Hall Fields ----------------------------------------
% Plot coll / non-coll
i_dim = 0;
iLog = 0;
smallPot = 2e-4;
% vangle = [1,0,0]; %YZ
% vangle = [0,-1,0]; %XZ
vangle = [0,0,1]; %XY
clims = [1, 1.5];
% clims = [-2, inf];
slices = {0,0,0};
plot_array = cell(2,1);
plotRuns = [1,2];
plot_lim = 400;%km

dimName = {'X', 'Y', 'Z'};

if i_dim == 0
    plot_str = 'E_{Hall,coll} / E_{Hall, no coll}';
else
    plot_str = ['E_{Hall,', dimName{i_dim},',coll} / E_{Hall,',...
        dimName{i_dim},', no coll}'];
end


for iRun = 1:2
    if i_dim == 0
        plot_array{iRun} = sqrt(sum(EHall_all{iRun}.^2,4));
        dimName = '';
        
    else
        plot_array{iRun} = EHall_all{iRun}(:,:,:,i_dim);
        dimName = dimname{i_dim};
    end
end
iSmallField = plot_array{1} < smallPot;
plot_array = plot_array{2} ./ plot_array{1};
if iLog 
    plot_array = log10(plot_array);
    plot_str = ['log_{10}(',plot_str,')'];
end

plot_array(iSmallField) = NaN;
nslice = sum(cellfun(@numel, slices));
f = figure;
ax1 = axes(f);

ax_slice = slice(ax1, X_mesh_all{plotRuns(1)},Y_mesh_all{plotRuns(1)},...
    Z_mesh_all{plotRuns(1)},...
    plot_array, slices{1},slices{2},slices{3});
            
for i_slice = 1:nslice
    ax_slice(i_slice).EdgeColor = 'none';
end
xlabel(ax1, 'X [km]');
ylabel(ax1, 'Y [km]');
zlabel(ax1, 'Z [km]');
set(ax1, 'FontSize', 14);

xlim(ax1, [-plot_lim,plot_lim]);
ylim(ax1, [-plot_lim,plot_lim]);
zlim(ax1, [-plot_lim,plot_lim]);

daspect(ax1, [1,1,1])
cb = colorbar(ax1);
colormap(ax1, cmap);
%             colormap(ax0, rb_map);
%             colormap(ax0, cMap);
cb.Label.String = plot_str;
cb.FontSize = 14;

caxis(ax1, clims);
cb.Location = 'eastoutside';
view(ax1, vangle);

r_pt = sqrt(X_mesh_all{1}.^2 + Y_mesh_all{1}.^2 + Z_mesh_all{1}.^2);
%% --------- Plot Pot Well ----------------------------------------
% Plot coll / non-coll
i_dim = 0;
iLog = 0;
smallPot = 5;
% vangle = [1,0,0]; %YZ
% vangle = [0,-1,0]; %XZ
vangle = [0,0,1]; %XY

clims = [0, 1.5];
% clims = [-1,0];
slices = {0,0,0};
plot_array = cell(2,1);
plotRuns = [1,2];
plot_lim = 400;%km

dimName = {'X', 'Y', 'Z'};

% if i_dim == 0
plot_str = 'V_{coll} / V_{no coll}';
% else
%     plot_str = ['V_{,', dimName{i_dim},',coll} / E_{,',...
%         dimName{i_dim},', no coll}'];
% end


for iRun = 1:2
    plot_array{iRun} = vPot_all{iRun}(:,:,:);
end
% ismallPot = plot_array{1} < smallPot;
% ismallPot = plot_array{1} < smallPot;
plot_array = plot_array{2} ./ plot_array{1};
plot_array(plot_array < 0) = nan;
% plot_array(ismallPot) = NaN;
if iLog 
    plot_array = log10(plot_array);
    plot_str = ['log_{10}(',plot_str,')'];
end


nslice = sum(cellfun(@numel, slices));
f = figure;
ax1 = axes(f);
% plot
plot_array = permute(plot_array, [2,1,3]);
ax_slice = slice(ax1, X_mesh_all{plotRuns(1)},Y_mesh_all{plotRuns(1)},...
    Z_mesh_all{plotRuns(1)},...
    plot_array, slices{1},slices{2},slices{3});
            
for i_slice = 1:nslice
    ax_slice(i_slice).EdgeColor = 'none';
end
xlabel(ax1, 'X [km]');
ylabel(ax1, 'Y [km]');
zlabel(ax1, 'Z [km]');
set(ax1, 'FontSize', 14);

xlim(ax1, [-plot_lim,plot_lim]);
ylim(ax1, [-plot_lim,plot_lim]);
zlim(ax1, [-plot_lim,plot_lim]);

daspect(ax1, [1,1,1])
cb = colorbar(ax1);
colormap(ax1, cmap);
%             colormap(ax0, rb_map);
%             colormap(ax0, cMap);
cb.Label.String = plot_str;
cb.FontSize = 14;

caxis(ax1, clims);
cb.Location = 'eastoutside';
view(ax1, vangle);

r_pt = sqrt(X_mesh_all{1}.^2 + Y_mesh_all{1}.^2 + Z_mesh_all{1}.^2);

%% All Field Plot 


f = figure('units', 'centimeters','Position', [2,2,30,16]); 
fSize = 14;
lef = [0.1, 0.1, 0.41, 0.41,0.71, 0.71] - 0.04;
bot = [0.51, 0.1, 0.51, 0.1, 0.51, 0.1];
wid = zeros(6,1) + 0.28;
heig = ones(6,1) * 0.35;
cb_wid = 0.01;
cb_heig = 0.85;
% ax = axes(f);
vangle = [0,0,1];
slices = {0,0,0};
% clims = [0,8];

% plas_dens = cell(2,1);
iRun = 2;
plotRuns = [1,2];
inc_gens = [1,2];
ax_plot = cell(6,1);
cb = cell(6,1);
cb_heigh = 0.01;
Xlimit = 380;
Ylimit = 300;
set(gca, 'FontSize', 14)
% -------------- Plot Hall fields -----------------------------------
clims = [0,8e-3];
i_dim = 0;
for i_plot = 1:2
    iRun = i_plot;
    ax_plot{i_plot} = subplot('Position', [lef(i_plot), bot(i_plot), wid(i_plot), heig(i_plot)]);
    plot_lim = 380;
    



    % --- Create Plot Array ------------------------

    % plot_array = permute(smooth_plas_dens{3-mod(i_plot,3)},[2,1,3]);
    % plot_str = 'log_{10}(N_e [cm^{-3}])';

    % ----- If comparing Runs -------
    slices = {0,0,0};
    nslice = sum(cellfun(@numel, slices));

    if i_dim == 0
        plot_array = sqrt(sum(EHall_all{iRun}.^2,4));
        dimName = '';
    else
        plot_array = EHall_all{iRun}(:,:,:,i_dim);
%         dimName = dimname{i_dim};
    end
    %            plot_array = permute( 2 * plas_dens{1} + 1 * plas_dens{2},[2,1,3]);
    ax_slice = slice(ax_plot{i_plot}, X_mesh_all{1},Y_mesh_all{1},Z_mesh_all{1},...
        plot_array, slices{1},slices{2},slices{3});

    for i_slice = 1:nslice
        ax_slice(i_slice).EdgeColor = 'none';
    end
    %     xlabel(ax1, 'X [km]');
    %     ylabel(ax1, 'Y [km]');
    %     zlabel(ax1, 'Z [km]');
    %     set(ax1, 'FontSize', 14);

    xlim(ax_plot{i_plot}, [-Xlimit,Xlimit]);
    ylim(ax_plot{i_plot}, [-Ylimit,Ylimit]);
    zlim(ax_plot{i_plot}, [-Ylimit,Ylimit]);

    daspect(ax_plot{i_plot}, [1,1,1])

    % rb_map = zeros(256, 3);
    % rb_map = redblue(256);
    % 
    % cb = colorbar(ax);
    
    
    colormap(ax_plot{i_plot}, cmap);
    % cb.Label.String = plot_str;
    % cb.FontSize = 14;

    view(ax_plot{i_plot}, vangle);
    caxis(ax_plot{i_plot}, clims);
    if mod(i_plot,2) == 1
        cb{i_plot} = colorbar(ax_plot{i_plot});
        cb{i_plot}.Location = 'northoutside'; 
        ax_plot{i_plot}.Position = [lef(i_plot), bot(i_plot), wid(i_plot), heig(i_plot)];
        cb{i_plot}.Position(4) = 0.01;
        cb{i_plot}.Position(2) = bot(1) + heig(1) + 0.02;
        cb{i_plot}.Label.String = 'Hall Field [V m^{-1}]';
        cb{i_plot}.Label.FontWeight = 'bold';
        cb{i_plot}.Label.FontSize = fSize;
        cb{i_plot}.FontSize = fSize;
    else 
        xlabel(ax_plot{i_plot}, 'X [km]');
    end
    
    ylabel(ax_plot{i_plot}, 'Y [km]');
    
    %     % cb.Location = 'eastoutside';
end
% -------------- Plot Ambipolar fields -----------------------------------
clims = [-1e-3,1e-3];
i_dim = 2;
for i_plot = 3:4
    iRun = i_plot - 2;
    ax_plot{i_plot} = subplot('Position', [lef(i_plot), bot(i_plot), wid(i_plot), heig(i_plot)]);
    plot_lim = 380;
    



    % --- Create Plot Array ------------------------

    % plot_array = permute(smooth_plas_dens{3-mod(i_plot,3)},[2,1,3]);
    % plot_str = 'log_{10}(N_e [cm^{-3}])';

    % ----- If comparing Runs -------
    slices = {0,0,0};
    nslice = sum(cellfun(@numel, slices));

    if i_dim == 0
        plot_array = sqrt(sum(EAmbi_all{iRun}.^2,4));
%         dimName = '';
    else
        plot_array = EAmbi_all{iRun}(:,:,:,i_dim);
%         dimName = dimname{i_dim};
    end
    %            plot_array = permute( 2 * plas_dens{1} + 1 * plas_dens{2},[2,1,3]);
    ax_slice = slice(ax_plot{i_plot}, X_mesh_all{1},Y_mesh_all{1},Z_mesh_all{1},...
        plot_array, slices{1},slices{2},slices{3});

    for i_slice = 1:nslice
        ax_slice(i_slice).EdgeColor = 'none';
    end
    %     xlabel(ax1, 'X [km]');
    %     ylabel(ax1, 'Y [km]');
    %     zlabel(ax1, 'Z [km]');
    %     set(ax1, 'FontSize', 14);

    xlim(ax_plot{i_plot}, [-Xlimit,Xlimit]);
    ylim(ax_plot{i_plot}, [-Ylimit,Ylimit]);
    zlim(ax_plot{i_plot}, [-Ylimit,Ylimit]);


    daspect(ax_plot{i_plot}, [1,1,1])

    % rb_map = zeros(256, 3);
    % rb_map = redblue(256);
    % 
    % cb = colorbar(ax1);
    colormap(ax_plot{i_plot}, cmap);
    
    if mod(i_plot,2) == 1
        cb{i_plot} = colorbar(ax_plot{i_plot});
        cb{i_plot}.Location = 'northoutside'; 
        ax_plot{i_plot}.Position = [lef(i_plot), bot(i_plot), wid(i_plot), heig(i_plot)];
        cb{i_plot}.Position(4) = 0.01;
        cb{i_plot}.Position(2) = bot(1) + heig(1) + 0.02;
        cb{i_plot}.Label.String = 'E_{Ambi, Y} [V m^{-1}]';
        cb{i_plot}.Label.FontWeight = 'bold';
        cb{i_plot}.Label.FontSize = fSize;
        cb{i_plot}.FontSize = fSize;
    else 
        xlabel(ax_plot{i_plot}, 'X [km]');
    end
    
    ylabel(ax_plot{i_plot}, 'Y [km]');
    % cb.Label.String = plot_str;
    % cb.FontSize = 14;

    view(ax_plot{i_plot}, vangle);
    caxis(ax_plot{i_plot}, clims);
    %     % cb.Location = 'eastoutside';
end


% -------------- Plot Potential Well -----------------------------------
clims = [-10,120];
clims = [-20, 200];
% i_dim = 2;
for i_plot = 5:6
    iRun = i_plot - 4;
    ax_plot{i_plot} = subplot('Position', [lef(i_plot), bot(i_plot), wid(i_plot), heig(i_plot)]);
    plot_lim = 380;



    % --- Create Plot Array ------------------------

    % plot_array = permute(smooth_plas_dens{3-mod(i_plot,3)},[2,1,3]);
    % plot_str = 'log_{10}(N_e [cm^{-3}])';

    % ----- If comparing Runs -------
    slices = {0,0,0};
    nslice = sum(cellfun(@numel, slices));
    plot_array = permute(vPot_all{iRun}, [2,1,3]);
%     if i_dim == 0
%         plot_array = sqrt(sum(EAmbi_all{iRun}.^2,4));
% %         dimName = '';
%     else
%         plot_array = EAmbi_all{iRun}(:,:,:,i_dim);
% %         dimName = dimname{i_dim};
%     end
    %            plot_array = permute( 2 * plas_dens{1} + 1 * plas_dens{2},[2,1,3]);
    ax_slice = slice(ax_plot{i_plot}, X_mesh_all{1},Y_mesh_all{1},Z_mesh_all{1},...
        plot_array, slices{1},slices{2},slices{3});

    for i_slice = 1:nslice
        ax_slice(i_slice).EdgeColor = 'none';
    end
    %     xlabel(ax1, 'X [km]');
    %     ylabel(ax1, 'Y [km]');
    %     zlabel(ax1, 'Z [km]');
    %     set(ax1, 'FontSize', 14);
    xlim(ax_plot{i_plot}, [-Xlimit,Xlimit]);
    ylim(ax_plot{i_plot}, [-Ylimit,Ylimit]);
    zlim(ax_plot{i_plot}, [-Ylimit,Ylimit]);


    daspect(ax_plot{i_plot}, [1,1,1])

    % rb_map = zeros(256, 3);
    % rb_map = redblue(256);
    % 
    % cb = colorbar(ax1);
    colormap(ax_plot{i_plot}, cmap);
    if mod(i_plot,2) == 1
        cb{i_plot} = colorbar(ax_plot{i_plot});
        cb{i_plot}.Location = 'northoutside'; 
        ax_plot{i_plot}.Position = [lef(i_plot), bot(i_plot), wid(i_plot), heig(i_plot)];
        cb{i_plot}.Position(4) = 0.01;
        cb{i_plot}.Position(2) = bot(1) + heig(1) + 0.02;
        cb{i_plot}.Label.String = 'Electric Potential [V]';
        cb{i_plot}.Label.FontWeight = 'bold';
        cb{i_plot}.Label.FontSize = fSize;
        cb{i_plot}.FontSize = fSize;
    else 
        xlabel(ax_plot{i_plot}, 'X [km]');
    end
    
    ylabel(ax_plot{i_plot}, 'Y [km]');
    % cb.Label.String = plot_str;
    % cb.FontSize = 14;

    view(ax_plot{i_plot}, vangle);
    caxis(ax_plot{i_plot}, clims);
    %     % cb.Location = 'eastoutside';
end

for i_plot = 1:6
    ax_plot{i_plot}.FontSize = 13;
end

%%
f2 = figure();
ax2 = axes(f2);
hold(ax2, 'on');
xpos = -40;
zpos = 0;
lStr = {'No Colls', 'Colls'};
for iRun = 1:2
    plot_array = permute(vPot_all{iRun}, [2,1,3]);
    xgrid = X_mesh_all{iRun}(1,:,1);
    ygrid = Y_mesh_all{iRun}(:,1,1);
    zgrid = Z_mesh_all{iRun}(1,1,:);
    ix_x = find(xgrid>=xpos,1,'first');
    ix_z = find(zgrid>=zpos,1,'first');
    dy = diff(ygrid(1:2));
    plot(ax2, ygrid+dy/2, -plot_array(:, ix_x, ix_z), 'LineWidth', 1.2,...
        'DisplayName',lStr{iRun});
end
legend(ax2, 'Location', 'best')
ylabel(ax2, 'Ambipolar Potential [V]')
xlabel('Y [km]');
xlim(ax2, [-300, 300]);
set(ax2, 'FontSize',15)
