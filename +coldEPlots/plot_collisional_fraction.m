function [f,ax]  = plot_collisional_fraction(innerData, vInnerData, ElDepthDist_all, RadialElDepthDist)


    exponentialModel = fittype('a * exp(-b * x) + c', 'independent', 'x', 'dependent', 'y');

    options = fitoptions(exponentialModel);
    options.StartPoint = [1.0, 0.1,0.];

    f = figure('Position',[10,10,1000,500]);
    lef = [0.07,0.38,0.69];
    bot = 0.15;
    wid = 0.27;
    heig = 0.8;
    ax = cell(1,3);
    iRun = 1;
    maxrProd = [300,300,300];
    fSize = 17;
    MinDepth = [-5,-5,-5];
    MaxDepth = [4,4,4];
    rexo = zeros(3,1);
    i_xs = 3;

    ylims = [0,0,0;1,1,1];
    fSize = 14;
    for iOGR = 1:3
        ax{iOGR} = axes(f, 'Position', [lef(iOGR),bot, wid, heig]);
        hold(ax{iOGR}, 'on');
        plotElDist = ElDepthDist_all{iRun, iOGR};
        plotElDist(:,:,1:5) = innerData.ElDepthDist_all{iRun,iOGR}; %Distribution using fields
        plotElDist(:,:,1) = vInnerData.ElDepthDist_all{iRun, iOGR};
        plotRadDist = RadialElDepthDist{iOGR}; % Distribution from Radial outflow
        

        if iOGR == 1
            plotElDist(2:end,:,:) = plotElDist(1:end-1,:,:);
            plotElDist(2:end,:,1) = plotElDist(1:end-1,:,1);
        end
        nRad = size(plotElDist,3);
        plotRadGrid = linspace(0, maxrProd(iOGR), nRad + 1);
        ElDepthGrid{iRun,iOGR} = MinDepth(iOGR) + (MaxDepth(iOGR) - MinDepth(iOGR)) * [0:size(plotElDist,1)]'...
                    / size(plotElDist,1);

        plotMid = 0.5 * (ElDepthGrid{iRun,iOGR}(1:end-1) + ElDepthGrid{iRun,iOGR}(2:end));
        
        ix_0 = find(ElDepthGrid{iRun,iOGR} < 0, 1, 'last');
    %     print(ix_0)
        WeightCollisional = sum(plotElDist(ix_0:end,i_xs,:),1);
        totWeight = sum(plotElDist(2:end,3,:),1);
        

        fracCollisional = reshape(WeightCollisional./ totWeight,[],1);
      
        plotRadMid =  0.5 * (plotRadGrid(2:end) + plotRadGrid(1:end-1));

        bar(ax{iOGR}, plotRadMid(1:end-1), fracCollisional(1:end-1),1);
        if iOGR == 2
            xlabel(ax{iOGR}, 'Cometocentric Distance at e^{-} creation [km]');
        end
        if iOGR == 1
            ylabel(ax{iOGR}, 'Collisional Fraction of Electrons');
        end
        grid(ax{iOGR}, 'on');
        ylim(ax{iOGR}, ylims(:,iOGR));
        set(ax{iOGR}, 'FontSize', fSize);

        % p
        plot(ax{iOGR}, [0,300], 1./exp([1.0, 1.0]))
        fitResult = fit(plotRadMid', fracCollisional, exponentialModel, options);
        yFit = feval(fitResult, plotRadMid);
        plot(ax{iOGR}, plotRadMid, yFit, 'r', 'DisplayName', 'Fitted Curve');
        rexo(iOGR) = interp1(yFit,plotRadMid,1/exp(1), 'cubic');

    end

end