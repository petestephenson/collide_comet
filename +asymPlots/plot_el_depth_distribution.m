function [f,ax] = plot_el_depth_distribution(ElDepthDist_all)
    f = figure('Position',[10,10,1000,500]);
    caseNameList = {'Ref','-X', '+Y', '+Z', '-Z'};
    nRuns = numel(ElDepthDist_all);
    nRad = size(ElDepthDist_all{1}, 3);
    ElDepthGrid_all = cell(1, nRuns);
    lef = linspace(0.08, 1-0.5*((2+0.5)/nRuns ), nRuns);
    bot = 0.15;
    wid = 1/(nRuns+2);
    heig = 0.8;
    ax = cell(1,nRuns);
    % iRun = 1;

    maxrProd = zeros(nRuns)+300;
    MinDepth = zeros(1,nRuns);
    MaxDepth = zeros(1,nRuns);
    MinDepth(:) = -5;
    MaxDepth(:) = 4;
    
    i_xs = 3;
    xlims = zeros(2, nRuns);
    ylims = zeros(2, nRuns);
    xlims(1,:) = -4;
    xlims(2,:) = 2;
    ylims(1,:) = 0;
    ylims(2,:) = 6e20;
    
    % ylims = [0,0,0;8e18,14e18, 1e20];

%     CDF_Depth = cell(3,2); % containts CDF of electron depth (summed over cometodist)
%     cdf_depth_wr = cell(3,2);
    fSize = 14;
    for iOGR = 1:nRuns

        ax{iOGR} = axes(f, 'Position', [lef(iOGR),bot, wid, heig]);
        

        
        plotElDist = sum(ElDepthDist_all{1, iOGR}(:,:,:,:,1),4) ;

        ElDepthGrid_all{1,iOGR} = MinDepth(iOGR) + (MaxDepth(iOGR) - MinDepth(iOGR)) * [0:size(plotElDist,1)]'...
                    / size(plotElDist,1);

        plotMid = 0.5 * (ElDepthGrid_all{1,iOGR}(1:end-1) + ElDepthGrid_all{1,iOGR}(2:end));
        barWidth = ElDepthGrid_all{1,iOGR}(2:end) - ElDepthGrid_all{1,iOGR}(1:end-1);
        plotEdges = ElDepthGrid_all{1,iOGR};
        
        hold(ax{iOGR}, 'on');
        
        plotRadGrid = linspace(0, maxrProd(iOGR), nRad + 1);
        
        
        plotDist = permute(plotElDist(:,i_xs,:),[1,3,2]);
        cmap = colorblindmap(256, 'heat');
        colormap(ax{iOGR}, cmap);
        barPlot = bar(ax{iOGR}, plotMid, plotDist(:,1:end-1),1,...
            'stacked', 'FaceColor', 'flat');
        for iRad = 1:nRad-1
            barPlot(iRad).CData = repmat(plotRadGrid(iRad),size(plotElDist,1) ,1);
            barPlot(iRad).HandleVisibility = 'off';

        end
        barPlot(end).HandleVisibility = 'on';
        barPlot(end).DisplayName = 'PiC Fields';
%         disp('testing')
%         size(plotMid)
%         size(sum(plotRadDist(:,i_xs,1:end-1),3))

%         plot(ax{iOGR}, plotMid, sum(plotRadDist(:,i_xs,1:end-1),3), 'k', 'LineWidth',1.5,...
%             'DisplayName', 'Radial Outflow');

        grid(ax{iOGR}, 'on');
        xlabel(ax{iOGR}, 'log_{10}(Electron Depth)');

        if iOGR == 1
            ylabel(ax{iOGR}, 'Weight of Electrons [electrons s^{-1}]');
        end
        
        ax{iOGR}.CLim = plotRadGrid([1,end]);
        xlim(ax{iOGR}, xlims(:,iOGR));
        ylim(ax{iOGR}, ylims(:,iOGR))
        set(ax{iOGR}, 'FontSize', fSize);
        if iOGR == 3
            cb = colorbar(ax{iOGR}, 'Position', [0.9,bot, 0.02, heig]);
            cb.Label.String = 'Cometocentric Distance at Production [km]';
            set(cb, 'FontSize', fSize);
            
        end
        % Set caxis for colorbar between 0 and maxrProd
%         caxis(ax{iOGR}, [0, maxrProd(iOGR)]);
        if iOGR == 1
            legend(ax{iOGR})
        end
        title(caseNameList{iOGR})
    end



        % Save figure to file E:\TestPlFigs\coldEReplots\elDepthDist.png
        
end
