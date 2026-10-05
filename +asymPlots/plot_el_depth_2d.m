function [f,ax, xInExo] = plot_el_depth_2d(ElDepthDist_all, iDim)
f = figure('Position',[10,10,1200,400]);
nRuns = numel(ElDepthDist_all);
%     lef = [0.08,0.35,0.62];s
lef = linspace(0.08, 1-0.5*((2+0.5)/nRuns ), nRuns);
% lef = [0.08];
caseNameList = {'Ref','-X', '+Y', '+Z', '-Z'};
axisName = {'X', 'Y', 'Z'};
bot = 0.15;
wid = 1/(nRuns+2);
heig = 0.8;
ax = cell(1,nRuns);

nRad = size(ElDepthDist_all{1}, 3);
nElDepthThetaBounds = size(ElDepthDist_all{1},4)+1;
nRuns = numel(ElDepthDist_all);
maxrProd = 300;
elDepthRadGrid = linspace(0, maxrProd, nRad + 1);
elDepthRadMids = 0.5* (elDepthRadGrid(1:end-1)+elDepthRadGrid(2:end));

iTau = 3;
ixExobase = cell(nRuns,1);
for iOGR = 1:nRuns
    elDepthDist2 = ElDepthDist_all{iOGR};
%     elDepthDist2 = sum(elDepthDist2, [4]);
    medianTau = cumsum(elDepthDist2(2:end,iTau,:,:,iDim), 1);
    medianTau = squeeze(medianTau ./medianTau(end,:,:,:));
    elDepthGrid = linspace(-5,4,51);
    elDepthGrid = elDepthGrid(2:end);
    [r,c,v] = find(medianTau >= 0.5);
    rowIndices = nan(size(elDepthDist2, [3,4]));
    % elDe
%     for iDim = 1:3
    for iPhi = 1:size(medianTau, 3)
        for col = 1:size(medianTau, 2)
            ixMedian = find(squeeze(medianTau(:, col, iPhi)) >= 0.5, 1, 'first');
            if ~isempty(ixMedian)
                rowIndices(col, iPhi) = ixMedian;
        %             valIndices(col) = v;
            end
        end
    end

%     end
    elDepthMedians = nan(size(rowIndices));
    elDepthMedians(rowIndices>0) = elDepthGrid(rowIndices(rowIndices>0));

    ixInExobase = elDepthMedians>0;
    ixExobase{iOGR} = ixInExobase;
    ax{iOGR} = axes(f, 'Position',[lef(iOGR), bot, wid, heig]);

    rMax = 300;

    radialGrid = elDepthRadGrid; % Linear grid from 0 to 300 km
    angularAxis = acos(linspace(-1, 1, nElDepthThetaBounds)); % Linear grid from -1 to 1
    radialGrid = 0.5 * (radialGrid(2:end)+ radialGrid(1:end-1));
    angularAxis = 0.5*(angularAxis(2:end)+angularAxis(1:end-1));
    
    
    % Create a meshgrid for the polar coordinates
    [theta, rho] = meshgrid(angularAxis, radialGrid);
    % Convert polar coordinates to Cartesian coordinates
    x = rho .* cos(theta);
    y = rho .* sin(theta);


% Create a 2D plot
    data = elDepthMedians(:,:);
    hold(ax{iOGR}, "on")

    pcolor(ax{iOGR}, x, y, data);
    xsInExo = x(ixInExobase);
    ysInExo = y(ixInExobase);
    k = convhull(xsInExo, ysInExo);
    plot(ax{iOGR}, xsInExo(k), ysInExo(k), 'k-')
    shading interp; % Interpolate the colors
    colormap jet; % Set the colormap
    
    % Set the aspect ratio to equal to avoid distortion
    axis equal;
    ax{iOGR}.XTick = ax{iOGR}.XTick(ax{iOGR}.XTick >=0);
    % Set the x-axis limits
    xlim([-300, 300]);
    ylim([0, 300]);
    
    % ax.XTicks = xticks(xticks >= 0);
    % Add labels and title
    xlabel('r [km]');
    ylabel('y');
    
    % Add faint grey lines at regular angles with labels
    hold on;
    angleTicks = linspace(0, pi, 7); % Adjust the number of lines as needed
    angleLabels = {'0°', '30°', '60°', '90°', '120°', '150°', '180°'}; % Adjust the labels as needed
    for i = 1:length(angleTicks)
        angle = angleTicks(i);
        xline = 300 * cos(angle);
        yline = 300 * sin(angle);
        plot([-xline, xline], [-yline, yline], 'Color', [0.2, 0.2, 0.2], 'LineStyle', '--');
        text(1.1 * xline - abs(0.05*xline), 1.1* yline, angleLabels{i}, 'Color', [0,0,0]);
    end
    % Add a colorbar
%     colorbar;
    ax{iOGR}.Box = 'off';
    ax{iOGR}.Color = 'none';
    % ax.YLabel = ''
    ax{iOGR}.YAxis.Visible = 'off';
    if iOGR == nRuns    
        cb = colorbar('Position', [0.9,bot, 0.02, heig]);
        cb.Label.String = 'log_{10}(Median \tau)';
    end
    caxis(ax{iOGR},[-2,0.5]);
    % set(gcf, 'Color', 'none');
    % Optional: Customize the plot appearance
    % For example, you can set the font size of the labels and title:
    set(gca, 'FontSize', 12);
    t1 = title(ax{iOGR}, caseNameList{iOGR});
    t1.Position(2) = t1.Position(2)+50;
    
end
sgtitle(['Median e Depths Relative to ',axisName{iDim},' axis']);

% 
% % f = figure();
% % ax = axes(f);
% hold(ax, 'on')
% % for iX= 1:nAsym
% %     for iY = 1:nAsym
% %         for iZ = 1:nAsym
% plot(ax, elDepthRadMids, elDepthMedians(:));
% %         end
% %     end
% % end
% %%
% iTau = 3;
% medianTau = cumsum(elDepthDist(2:end,iTau,:,:,:), 1);
% medianTau = squeeze(medianTau ./medianTau(end,:,:,:,:));
% nAsym=1;
% % medianIndex = find(medianTau >= 0.5, 1, "first");
% elDepthGrid = linspace(-5,4,51);
% elDepthGrid = elDepthGrid(2:end);
% [r,c,v] = find(medianTau >= 0.5);
% % binEdges = 1:size(medianTau, 1) + 1; % Assuming the bins are equally spaced
% % medianValues = mean(binEdges([medianIndex, medianIndex+1]));
% 
% 
% % Initialize the arrays to store the row, column, and value of the first true element
% rowIndices = zeros(size(medianTau, 2),size(medianTau, 3),nAsym);
% colIndices = 1:size(medianTau, 2);
% % valIndices = zeros( size(medianTau, 2), size(medianTau, 3));
% 
% % Iterate over the columns and find the first true value greater than or equal to 0.5
% for iAsym = 1:nAsym
%     for col = 1:size(medianTau, 2)
%         for theta = 1:size(medianTau, 3)
%             ixMedian = find(squeeze(medianTau(:, col, theta, iAsym)) >= 0.5, 1, 'first');
%             if ~isempty(ixMedian)
%                 rowIndices(col, theta, iAsym) = ixMedian;
%     %             valIndices(col) = v;
%             end
%         end
%     
%     end
% end
% % Display the results
% disp('First true values greater than or equal to 0.5:');
% disp('Column | Row | Value');
% disp([colIndices, rowIndices, valIndices]);
%%
% elDepthMedians = nan(size(rowIndices));
% elDepthMedians(rowIndices>0) = elDepthGrid(rowIndices(rowIndices>0));
%%
% % Define the radial and angular axes
% rMax = 300;
% % nElDepthThetaBounds 1= size()
% radialGrid = elDepthRadGrid; % Linear grid from 0 to 300 km
% angularAxis = acos(linspace(-1, 1, nElDepthThetaBounds)); % Linear grid from -1 to 1
% radialGrid = 0.5 * (radialGrid(2:end)+ radialGrid(1:end-1));
% angularAxis = 0.5*(angularAxis(2:end)+angularAxis(1:end-1));
% 
% thetaTicks = linspace(0,pi,6);
% xThetaTicks = rMax * cos(thetaTicks);
% yThetaTicks = rMax * sin(thetaTicks);
% 
% 
% % Create a meshgrid for the polar coordinates
% [theta, rho] = meshgrid(angularAxis, radialGrid);
% % Convert polar coordinates to Cartesian coordinates
% x = rho .* cos(theta);
% y = rho .* sin(theta);
% 
% 
% % Generate some data for the plot (replace with your own data)
% % data = cos(theta) .* rho; % Example data, replace with your own variable
% % data = rowIndices;
% data = elDepthMedians(:,:,1);
% 
% 
% % Create a 2D plot
% f = figure;
% ax  = axes(f);
% hold(ax, "on")
% pcolor(ax, x, y, data);
% if nAsym ==2
%     y2 = - rho .* sin(theta);
%     data2 = elDepthMedians(:,:,2);
%     pcolor(ax, x, y2, data2);
% end
% 
% shading interp; % Interpolate the colors
% colormap jet; % Set the colormap
% 
% % Set the aspect ratio to equal to avoid distortion
% axis equal;
% ax.XTick = ax.XTick(ax.XTick >=0);
% % Set the x-axis limits
% xlim([-300, 300]);
% 
% % Set the y-axis limits
% if nAsym ==2
%     ylim([-300, 300]);
% else
%     ylim([0,300]);
% end
% % ax.XTicks = xticks(xticks >= 0);
% % Add labels and title
% xlabel('r [km]');
% ylabel('y');
% 
% % Add faint grey lines at regular angles with labels
% hold on;
% angleTicks = linspace(0, pi, 7); % Adjust the number of lines as needed
% angleLabels = {'0°', '30°', '60°', '90°', '120°', '150°', '180°'}; % Adjust the labels as needed
% for i = 1:length(angleTicks)
%     angle = angleTicks(i);
%     xline = 300 * cos(angle);
%     yline = 300 * sin(angle);
%     plot([-xline, xline], [-yline, yline], 'Color', [0.2, 0.2, 0.2], 'LineStyle', '--');
%     text(1.1 * xline - abs(0.05*xline), 1.1* yline, angleLabels{i}, 'Color', [0,0,0]);
% end
% % Add a colorbar
% colorbar;
% ax.Box = 'off';
% ax.Color = 'none';
% % ax.YLabel = ''
% ax.YAxis.Visible = 'off';
% cb = colorbar();
% cb.Label.String = 'log_{10}(Median \tau)';
% caxis(ax,[-2,0.5])
% % set(gcf, 'Color', 'none');
% % Optional: Customize the plot appearance
% % For example, you can set the font size of the labels and title:
% set(gca, 'FontSize', 12);
% % title('2D Plot with Radial and Angular Axes');
% 
% %% 2nd orientation
% 
% [angularAxis, radialGrid] = makePlotGrids(nElDepthThetaBounds, elDepthRadGrid);
% plot_elDepth_XZ(elDepthMedians, radialGrid, angularAxis)
