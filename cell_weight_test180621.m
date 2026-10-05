global cell_weight
global X_mesh Y_mesh Z_mesh
global x_grid y_grid z_grid
global dx_grid dy_grid dz_grid

iRun = 1;

vol_grid = dx_grid{iRun} .* dy_grid{iRun}' .* ...
    reshape(dz_grid{iRun},1,1,[]);
x_mid = 0.5 * (x_grid{iRun}(1:end-1) + x_grid{iRun}(2:end));
y_mid = 0.5 * (y_grid{iRun}(1:end-1) + y_grid{iRun}(2:end));
z_mid = 0.5 * (z_grid{iRun}(1:end-1) + z_grid{iRun}(2:end));

r_cell = sqrt(x_mid.^2 .* y_mid'.^2 .* ...
    reshape(z_mid,1,1,[]).^2);

%% Check weighting

nDensTest = cell_weight{iRun} .* r_cell.^2 ./ (vol_grid);


plot_array = permute(nDensTest,[2,1,3]);
slices = {0,0,0};
nslice = 3;
plot_lim = 400;
f = figure;
ax0 = axes(f);
ax_slice = slice(ax0, X_mesh{1},Y_mesh{1},Z_mesh{1},...
    plot_array, slices{1},slices{2},slices{3});

for i_slice = 1:nslice
    ax_slice(i_slice).EdgeColor = 'none';
end
xlabel(ax0, 'X [km]');
ylabel(ax0, 'Y [km]');
zlabel(ax0, 'Z [km]');
set(ax0, 'FontSize', 14);

xlim(ax0, [-plot_lim,plot_lim]);
ylim(ax0, [-plot_lim,plot_lim]);
zlim(ax0, [-plot_lim,plot_lim]);

daspect(ax0, [1,1,1])

rb_map = zeros(256, 3);
rb_map = redblue(256);
%             colormap(ax0, rb_map);

cb = colorbar(ax0);
colormap(ax0, jet);
cb.Label.String = plot_str;
cb.FontSize = 14;