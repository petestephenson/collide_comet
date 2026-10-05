
nGrid = 500;
plotLim = 60;
[X_sph, Y_sph] = sphere(50);
X_sph = X_sph * 1.7;
Y_sph = Y_sph * 1.7;
Z_Sph = Z_Sph * 1.7;
xGrid = linspace(-plotLim, plotLim, nGrid);
yGrid = xGrid; 
% zGrid = xGrid;

[xmesh, ymesh] = meshgrid(xGrid, yGrid);
rmesh = sqrt(xmesh.^2 + ymesh.^2);

cData = (1 ./ rmesh) .^ 0.1;
f = figure;
ax = axes(f);
hold(ax, 'on');
surf(ax, X_sph, Y_sph, Z_Sph, 'FaceColor', 'k');
s1 = surf(ax, xmesh, ymesh, zeros(size(xmesh)), cData);
% s1 = slice(ax, xmesh, ymesh, zmesh, cData, 0, 0, 0);
for iSlice = 1:numel(s1)
    s1(iSlice).EdgeColor = 'none';
end
cmap = colorblindmap(256, 'purple');
colormap(ax, cmap);
% colorbar(ax);
caxis(ax,[min(cData*0.9,[],'all'),max(cData*0.75,[],'all')])
daspect(ax, [1,1,1]);
view(ax, 2);