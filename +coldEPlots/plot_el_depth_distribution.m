function [f,ax] = plot_el_depth_distribution(innerData, vInnerData, ElDepthDist_all, RadialElDepthDist)
    f = figure('Position',[10,10,1000,500]);
    lef = [0.08,0.36,0.64];
    bot = 0.15;
    wid = 0.25;
    heig = 0.8;
    ax = cell(1,3);
    % iRun = 1;

    maxrProd = [300,300,300];
    MinDepth = [-5,-5,-5];
    MaxDepth = [4,4,4];

    i_xs = 3;
    xlims = [-4,-4,-4;4,4,4];
    % ylims = [0,0,0;8e18,14e18, 1e20];
    ylims = [0,0,0;4e18,1e19, 5e19];
    % i_xs = 3;
    xlims = [-4,-4,-4;2,2,2];
    % ylims = [0,0,0;6e18,12e18, 1e20];
    CDF_Depth = cell(3,2); % containts CDF of electron depth (summed over cometodist)
    cdf_depth_wr = cell(3,2);
    fSize = 14;
    for iOGR = 1:3

        ax{iOGR} = axes(f, 'Position', [lef(iOGR),bot, wid, heig]);
        

        % plotElDist = innerData.ElDepthDist_all{iRun,iOGR};

        % plotElDist = ElDepthDist_all{1, iOGR};
        plotElDist = ElDepthDist_all{1, iOGR}+ElDepthDist_all{2, iOGR};
        plotElDist(:,:,1:5) = innerData.ElDepthDist_all{1,iOGR}; %Distribution using fields
        plotElDist(:,:,1) = vInnerData.ElDepthDist_all{1, iOGR};
        % plotElDist = ElDepthDist_all{2, iOGR};

    
        plotRadDist = RadialElDepthDist{iOGR};

        disp(size(plotRadDist))
        if iOGR == 1
            
            plotRadDist = reshape(plotRadDist,2,50,12,30);
            plotRadDist = reshape(sum(plotRadDist, 1), 50,12,30);
            % plotElDist(:,:,2:end) = plotElDist(:,:,1:end-1);
        end
        plotRadDist = plotRadDist /2.5;
        nRad = size(plotElDist,3);
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
        disp('testing')
        size(plotMid)
        size(sum(plotRadDist(:,i_xs,1:end-1),3))

        plot(ax{iOGR}, plotMid, sum(plotRadDist(:,i_xs,1:end-1),3), 'k', 'LineWidth',1.5,...
            'DisplayName', 'Radial Outflow');

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
        if iOGR == 1
            legend(ax{iOGR})
        end
        % Calculate CDF of e depths
        depthVals = sum(plotDist, 2);
        depthValsRad = reshape(sum(plotRadDist(:, i_xs,:), 3), [],1);
    %     depthValsRad = reshape(sum(plotRadDist(:, i_xs,:), 3), [],1);
        depthVals([1,end]) = 0;
        depthValsRad([1,end]) = 0;
        CDF_Depth{iOGR,1} = cumsum(depthVals./sum(depthVals)) ;
        CDF_Depth{iOGR,2} = cumsum(depthValsRad ./sum(depthValsRad)) ;
        plotDist(1,:) = 0;
        plotRadDist(1,:,:) = 0;

        cdf_depth_wr{iOGR, 1} = cumsum(plotDist./sum(plotDist, 1), 1);

        cdf_depth_wr{iOGR, 2} = squeeze(cumsum(plotRadDist(:,i_xs,:)./sum(plotRadDist(:,i_xs,:),1)));
    end



        % Save figure to file E:\TestPlFigs\coldEReplots\elDepthDist.png
        
end
