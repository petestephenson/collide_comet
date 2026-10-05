function [f,ax] = plot_el_depth_2d2ref(ElDepthDist_all, iDim)
f = figure('Position',[10,10,1200,400]);
nRuns = numel(ElDepthDist_all)-1;
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
% nRuns = numel(ElDepthDist_all);
maxrProd = 300;
elDepthRadGrid = linspace(0, maxrProd, nRad + 1);
elDepthRadMids = 0.5* (elDepthRadGrid(1:end-1)+elDepthRadGrid(2:end));

iTau = 3;

for iOGR = 0:nRuns
    elDepthDist2 = ElDepthDist_all{iOGR+1};
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
    if iOGR == 0
        elDepthMedianRef = elDepthMedians;
        continue
    end
    ixInExobase = elDepthMedians>0;
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
    data = elDepthMedians-elDepthMedianRef;
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
        cb.Label.String = 'log_{10}(Median \tau/ Reference)';
    end
    caxis(ax{iOGR},[-1,1]);
    % set(gcf, 'Color', 'none');
    % Optional: Customize the plot appearance
    % For example, you can set the font size of the labels and title:
    set(gca, 'FontSize', 12);
    t1 = title(ax{iOGR}, caseNameList{iOGR+1});
    t1.Position(2) = t1.Position(2)+50;
    
end
sgtitle(['Median e Depths/ Ref Relative to ',axisName{iDim},' axis']);
print(f, '-clipboard', '-dbitmap');
