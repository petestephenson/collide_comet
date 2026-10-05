function plot_elDepth_XZ(elDepthMedians, radialGrid, angularGrid)
% Define the radial and angular axes
rMax = 300;
thetaTicks = linspace(0,pi,6);
xThetaTicks = rMax * sin(thetaTicks);
yThetaTicks = rMax * cos(thetaTicks);


% Create a meshgrid for the polar coordinates
[theta, rho] = meshgrid(angularGrid, radialGrid);
% Convert polar coordinates to Cartesian coordinates
x = rho .* sin(theta);
y = rho .* cos(theta);


% Generate some data for the plot (replace with your own data)
data = elDepthMedians(:,:,1);


% Create a 2D plot
f = figure;
ax  = axes(f);
hold(ax, "on")
pcolor(ax, x, y, data);

shading interp; % Interpolate the colors
colormap jet; % Set the colormap

% Set the aspect ratio to equal to avoid distortion
axis equal;
ax.XTick = ax.XTick(ax.XTick >=0);
% Set the x-axis limits

xlim([0,300]);
ylim([-300, 300]);

% Set the y-axis limits

%     ylim([0,300]);

% Add labels and title
ylabel('r [km]');
xlabel('y');


hold on;
angleTicks = linspace(0, pi, 7); % Adjust the number of lines as needed
angleLabels = {'0°', '30°', '60°', '90°', '120°', '150°', '180°'}; % Adjust the labels as needed
for i = 1:length(angleTicks)
    angle = angleTicks(i);
    xline = 300 * sin(angle);
    yline = 300 * cos(angle);
    plot([-xline, xline], [-yline, yline], 'Color', [0.2, 0.2, 0.2], 'LineStyle', '--');
    text(1.1 * xline - abs(0.05*xline), 1.1* yline, angleLabels{i}, 'Color', [0,0,0]);
end
% Add a colorbar
    colorbar;
    ax.Box = 'off';
    ax.Color = 'none';
% ax.YLabel = ''
    ax.XAxis.Visible = 'off';
    cb = colorbar();
    cb.Label.String = 'log_{10}(Median \tau)';
    caxis(ax,[-2,0.5])
    set(gca, 'FontSize', 12);


end