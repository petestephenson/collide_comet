function [f,ax] = make_isocontour_plot(plas_dens_cold_all,PicData_all,vExB, fieldProps, X_mesh_all,Y_mesh_all,Z_mesh_all, XPicMesh_all, YPicMesh_all, ZPicMesh_all)

iRun = 2;
minVal = [1e-2,1e-2,1e-1];
plot_lim = 350;

plot_array = plas_dens_cold_all{1,iRun};% plas_dens_cold_all{2,iRun};
ixNan = plot_array < minVal(2);
plot_array(ixNan) = NaN;
plot_array = permute(plot_array ,[2,1,3]);
iLog = 0;
if iLog
    plot_array = log10(abs(plot_array)); 
    plot_str = ['log_{10}(',plot_str,')'];
end
% 

isoLevels = [1e-2, 1e-1,1e0];
isoAlpha = [0.5, 0.5, 0.5];

colorIndex = [1,2,3,4];

colors = colorblindmap(numel(isoLevels)+2, 'blue2');


disp(['Size of colors: ', num2str(size(colors))]);

colors = colors(end-2:end,:);

set(0,'DefaultAxesFontSize',15)
f = figure('Position',[10,10,800,500]);
ax = axes(f);
hold(ax, 'on');
% Add inset plot to the figure 
ax_inset = axes(f, 'Position',[0.57,0.65,0.3,0.3]);
hold(ax_inset, 'on');

for iLevel = 1:numel(isoLevels)
    % # print isoLevels(iLevel) and iLevel
    disp(['Plotting isocontour ', num2str(iLevel), ' of ', num2str(numel(isoLevels))]);
    [faces, verts] = isosurface(X_mesh_all{iRun,2}, Y_mesh_all{iRun,2},Z_mesh_all{iRun,2},...
        plot_array, isoLevels(iLevel));
    
    % Check if verts is not empty before plotting
    if ~isempty(verts)
        
        color_data = repmat(colors(iLevel, :), size(verts, 1), 1); % Color from colormap
        
        
        p = patch(ax, 'Vertices', verts, 'Faces', faces, 'FaceColor', 'interp', 'FaceVertexCData', color_data, 'EdgeColor', 'none','DisplayName', ['n_{e,c}=', num2str(isoLevels(iLevel)), 'cm^{-3}']);
        p.FaceAlpha = isoAlpha(iLevel);

        p = patch(ax_inset, 'Vertices', verts, 'Faces', faces, 'FaceColor', 'interp', 'FaceVertexCData', color_data, 'EdgeColor', 'none','DisplayName', ['n_{e,c}=', num2str(isoLevels(iLevel)), 'cm^{-3}']);
        p.FaceAlpha = isoAlpha(iLevel);
    end
end

% Calculate streamlines using vExB drift velocities
% stream3(X,Y,Z,U,V,W,startx,starty,startz)
% With U,V,W the components of the velocity field and startx, starty, startz
% the starting points of the streamlines.
% Create meshgrid of starting points in y and z from -plot_lim to plot_lim
% with 10 points in each direction
nLines = 3;

streamLims = 100;
startx = linspace(-streamLims, streamLims, nLines);
starty = linspace(-streamLims, streamLims, nLines);
startz = zeros(size(starty));
[startXmesh, startYmesh] = meshgrid(startx, starty);
startZmesh = zeros(size(startXmesh));




l1 = stream3(XPicMesh_all{iRun}, YPicMesh_all{iRun},ZPicMesh_all{iRun},...
    vExB{2}(:,:,:,1), vExB{2}(:,:,:,2), vExB{2}(:,:,:,3),...
    startXmesh, startYmesh, startZmesh);
% Calculate streamlines in -ve direction
l2 = stream3(XPicMesh_all{iRun}, YPicMesh_all{iRun},ZPicMesh_all{iRun},...
    -vExB{2}(:,:,:,1), -vExB{2}(:,:,:,2), -vExB{2}(:,:,:,3),...
    startXmesh, startYmesh, startZmesh);

% plot each of the streamlines
for iLine = 1:(nLines*nLines)
    plot3(ax, l1{iLine}(:,1), l1{iLine}(:,2), l1{iLine}(:,3), 'k', 'HandleVisibility', 'off', 'LineWidth', 1.5);
    plot3(ax, l2{iLine}(:,1), l2{iLine}(:,2), l2{iLine}(:,3), 'k', 'HandleVisibility', 'off', 'LineWidth', 1.5);

    % plot the streamlines in the inset plot
    plot3(ax_inset, l1{iLine}(:,1), l1{iLine}(:,2), l1{iLine}(:,3), 'k', 'HandleVisibility', 'off', 'LineWidth', 1.5);
    plot3(ax_inset, l2{iLine}(:,1), l2{iLine}(:,2), l2{iLine}(:,3), 'k', 'HandleVisibility', 'off', 'LineWidth', 1.5);
end



% Plot streamlines for the B field from PicData_all{2}.B_Pic. Using xPicMesh_all{2} etc. 
%  Use start points from a meshgrid 
% in x and z from -bLims to bLims with 10 points in each direction
nBLines = 3;
bLims = [0, 200];
startx = linspace(bLims(1), bLims(2), nBLines);
startz = linspace(bLims(1), bLims(2), nBLines);
[startXmesh, startZmesh] = meshgrid(startx, startz);
startYmesh = zeros(size(startXmesh));
bDataMesh = permute(PicData_all{iRun}.B_Pic, [2,1,3,4,]);
% Calculate streamlines in +ve direction
lB1 = stream3(XPicMesh_all{iRun}, YPicMesh_all{iRun},ZPicMesh_all{iRun},...
    bDataMesh(:,:,:,1), bDataMesh(:,:,:,2), bDataMesh(:,:,:,3),...
    startXmesh, startYmesh, startZmesh);   
% Calculate streamlines in -ve direction
lB2 = stream3(XPicMesh_all{iRun}, YPicMesh_all{iRun},ZPicMesh_all{iRun},...
    -bDataMesh(:,:,:,1), -bDataMesh(:,:,:,2), -bDataMesh(:,:,:,3),...
    startXmesh, startYmesh, startZmesh);


% plot each of the streamlines with color 'm' with line width 2
for iLine = 1:(nBLines*nBLines)
    plot3(ax, lB1{iLine}(:,1), lB1{iLine}(:,2), lB1{iLine}(:,3), 'm', 'HandleVisibility', 'off', 'LineWidth', 2);
    plot3(ax, lB2{iLine}(:,1), lB2{iLine}(:,2), lB2{iLine}(:,3), 'm', 'HandleVisibility', 'off', 'LineWidth', 2);
    % plot the streamlines in the inset plot
    plot3(ax_inset, lB1{iLine}(:,1), lB1{iLine}(:,2), lB1{iLine}(:,3), 'm', 'HandleVisibility', 'off', 'LineWidth', 2);
    plot3(ax_inset, lB2{iLine}(:,1), lB2{iLine}(:,2), lB2{iLine}(:,3), 'm', 'HandleVisibility', 'off', 'LineWidth', 2);
end


slices = {[],[],[0,200]};
sliceAlphas = [0.8, 0.8];

% Plot slices of E ambipolar y component
% Use colormap from colorblindmap
rb_map = colorblindmap(256, 'heat');
colormap(ax, rb_map);


sPot = slice(ax, fieldProps.X_mesh_field_all{2}, fieldProps.Y_mesh_field_all{2},fieldProps.Z_mesh_field_all{2},...
    permute(fieldProps.vPot_all{2}, [2,1,3]), slices{1},slices{2},slices{3});
% Set edgecolor of sPot to none
for i_slice = 1:numel(sPot)
    sPot(i_slice).EdgeColor = 'none';
    sPot(i_slice).FaceAlpha = sliceAlphas(i_slice);
    % Set handle visibilty to off so that it doesn't appear in legend
    sPot(i_slice).HandleVisibility = 'off';
end
% Add colorbar
cb = colorbar(ax);
caxis(ax, [0,50]);
cb.Label.String = 'Ambipolar Potential [V]';
cb.FontSize = 14;

% Set xyz limits and aspect ratio
xlim(ax, [-150,350]);
ylim(ax, [-plot_lim,plot_lim]);
zlim(ax, [-150,350]);
daspect(ax, [1,1,1])

% add legend for isocontours
leg = legend(ax, 'show', 'Location', 'northoutside', 'Orientation', 'horizontal');

leg.Position(1) = leg.Position(1) - 0.05;

% add xyz lables
xlabel(ax, 'X [km]');
ylabel(ax, 'Y [km]');
zlabel(ax, 'Z [km]');

view(ax, -60, 11)
view(ax_inset, -22, -7)

% Remove labels from inset plot
% and the ticks from the axes
xlabel(ax_inset, '');
ylabel(ax_inset, '');
zlabel(ax_inset, '');
ax_inset.XTick = [];
ax_inset.YTick = [];
ax_inset.ZTick = [];
% Set xyz limits and aspect ratio
xlim(ax_inset, [-150,350]);
ylim(ax_inset, [-plot_lim,plot_lim]);
zlim(ax_inset, [-150,350]);
daspect(ax_inset, [1,1,1])





end