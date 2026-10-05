%% ----------------------- Load Field data points -------------------------
EB_Folder = 'C:\Users\pete4\OneDrive - Imperial College London\Energy_Degradation\EB_fields';
% EB_file = 'EB_fields_210F5.nc';
% EB_file = 'EB_fields_210F1_trim.nc';
EB_file = 'EB_fields_Q1e26.nc';
EB_file = [EB_Folder,'\', EB_file];
ncid = netcdf.open(EB_file);
finfo = ncinfo(EB_file);
E_ave = netcdf.getVar(ncid,0);
B_ave = netcdf.getVar(ncid,1);
[~, nX] = netcdf.inqDim(ncid,0);
[~, nY] = netcdf.inqDim(ncid,1);
[~, nZ] = netcdf.inqDim(ncid,2);

netcdf.close(ncid);

%% ------------------------------------------------------------------------
% ----------------------- Create field Plots ------------------------------
% -------------------------------------------------------------------------

f = figure('units', 'centimeters','Position', [2,2,30,16]); 

lef = [0.08, 0.39, 0.72, 0.08, 0.39, 0.72] - 0.04;
bot = [0.58, 0.58, 0.58, 0.1, 0.1, 0.1];
wid = zeros(6,1) + 0.29;
heig = 0.39 * ones(6,1);
cb_wid = 0.01;
cb_heig = 0.85;
% ax = axes(f);
ax_plot = cell(6,1);
x_comet = [100.45,145,145];
dx_EB = 2200 / 287;
grid = cell(3,1);
grid{1} = ([1:nX] - x_comet(1)) * dx_EB; grid{2} = ([1:nY]- x_comet(2)) * dx_EB;
grid{3} = ([1:nZ] - x_comet(3)) * dx_EB;
[X,Y,Z] = meshgrid(grid{1}, grid{2},grid{3});
slices = {[],[],0};
plotB = permute(B_ave * 1e9, [2,1,3,4]);
plotE = permute(E_ave, [2,1,3,4]);

BLims = [-4, 4;...
    6,7;...
    -2.5,2.5];
ELims = [-5, 5;...
    -1.5, 1.5;...
    -3, 3] * 1e-3;
    
xlims = [-770,770];
ylims = [-770,770];
zlims = [-770,770];
dimName = {'X', 'Y', 'Z'};
for i_plot = 1:3
    ax_plot{i_plot} = subplot('Position',  [lef(i_plot), bot(i_plot), wid(i_plot), heig(i_plot)]);
    
    ax_slice = slice(X, Y, Z, plotB(:,:,:,i_plot), slices{1},slices{2},slices{3});
    ax_slice.EdgeColor = 'none';
%     xlabel(ax_plot{i_plot},'X [km]');
    
    zlabel(ax_plot{i_plot},'Z [km]');
    set(ax_plot{i_plot}, 'Fontsize', 15);
    cb = colorbar(ax_plot{i_plot}, 'Location', 'eastoutside');
    cb.Position(3) = 0.01;
    cb.Label.String = ['B_{',dimName{i_plot},'} [nT]'];
    colormap(ax_plot{i_plot}, jet);
    view(ax_plot{i_plot}, [0,0,1]);
    daspect(ax_plot{i_plot}, [1,1,1]);
    caxis(ax_plot{i_plot}, BLims(i_plot,:));
    xlim(ax_plot{i_plot}, xlims);
    ylim(ax_plot{i_plot}, ylims);
    zlim(ax_plot{i_plot}, zlims);
    if i_plot == 1
        ylabel(ax_plot{i_plot},'Y [km]');
    end
    
end

for i_plot = 4:6
    ax_plot{i_plot} = subplot('Position',  [lef(i_plot), bot(i_plot), wid(i_plot), heig(i_plot)]);
    
    ax_slice = slice(X, Y, Z, plotE(:,:,:,i_plot-3), slices{1},slices{2},slices{3});
    ax_slice.EdgeColor = 'none';
    xlabel(ax_plot{i_plot},'X [km]');
    if i_plot == 4
        ylabel(ax_plot{i_plot},'Y [km]');
    end
    zlabel(ax_plot{i_plot},'Z [km]');
    set(ax_plot{i_plot}, 'Fontsize', 15);
    cb = colorbar(ax_plot{i_plot}, 'Location', 'eastoutside');
    cb.Position(3) = 0.01;
    cb.Label.String = ['E_{',dimName{i_plot - 3},'} [V m^{-1}]'];
    colormap(ax_plot{i_plot}, jet);
    view(ax_plot{i_plot}, [0,0,1]);
    daspect(ax_plot{i_plot}, [1,1,1]);
    caxis(ax_plot{i_plot}, ELims(i_plot - 3,:));
    xlim(ax_plot{i_plot}, xlims);
    ylim(ax_plot{i_plot}, ylims);
    zlim(ax_plot{i_plot}, zlims);
end