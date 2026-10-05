function [f,ax] = plot_cold_e_frac_vs_r(plas_dens_cold_all, plas_dens_all, X_mesh_all, Y_mesh_all, Z_mesh_all)

    f = figure();
    numBinsX = 400;
    numBinsY = 300;
    cminValue = 0;
    cmaxValue = 0.2;
    lef = [0.08,0.389,0.7];
    bot = 0.13;
    wid = 0.28;
    heig = 0.8;
    ax = cell(1,3);
    plotRuns = [1,2];
    rb_map = colorblindmap(256, 'heat');
    ax = cell(1,3);
    for iOGR = 1:3
        ax{iOGR} = axes(f, 'Position',[lef(iOGR), bot, wid, heig]);
        plot_array = (plas_dens_cold_all{1, iOGR}) ./ ...
            (plas_dens_all{1, iOGR});
        plot_array = permute(plot_array, [2,1,3]);
        plot_str = 'N_{e,cold} / N_{e}';


        rMesh = sqrt(X_mesh_all{plotRuns(1), iOGR}.^2+Y_mesh_all{plotRuns(1), iOGR}.^2+...
            Z_mesh_all{plotRuns(1), iOGR}.^2);
        rMeshFlat = reshape(rMesh, 1,[]);
        plot_array = reshape(plot_array, 1,[]);
    
        % Create a 2D histogram
        histogram2 = histcounts2(rMeshFlat, plot_array, [numBinsX, numBinsY]);
        histogram2 = histogram2 ./ sum(histogram2,2);
        colormap(rb_map)
    
        
        [xGrid, yGrid] = meshgrid(linspace(min(rMeshFlat), max(rMeshFlat), numBinsX), linspace(min(plot_array), max(plot_array), numBinsY));
    
        ax{iOGR}.XDir = "normal";
        ax{iOGR}.YDir = "normal";
        xlabel(ax{iOGR}, 'CometoDist [km]');
        if iOGR == 1
            ylabel(ax{iOGR}, 'Cold e frac');
        end
        cb = colorbar(ax{iOGR}, "north");
        clim(ax{iOGR}, [cminValue, cmaxValue]);
        ylim(ax{iOGR}, [0,0.1])
    %     hold(ax{iOGR}, 'on')
        rMesh_min = min(rMeshFlat);
        rMesh_max = max(rMeshFlat);
        binEdges = linspace(rMesh_min, rMesh_max, numBinsX + 1);
        %   Find the median values in each bin
    %     [~, binIndices] = histc(rMeshFlat, binEdges);
    %     medianValues = accumarray(binIndices, plot_array, [], @median);
        % Find the median values in each bin
        medianValues = zeros(1, numBinsX);
        uQvals = zeros(1, numBinsX);
        for i = 1:numBinsX
            % Find indices of elements in the current bin
            indices = rMeshFlat >= binEdges(i) & rMeshFlat < binEdges(i+1);
            
            % Calculate the median for elements in the current bin
            medianValues(i) = median(plot_array(indices));
            uQvals(i) = quantile(plot_array(indices), 0.75);
        end
        hold(ax{iOGR}, 'on');
        % Plot the median values against the bin centers
        binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;
        scatter(ax{iOGR},binCenters, medianValues, 'o', 'filled');
        scatter(ax{iOGR},binCenters, uQvals, 'ro', 'filled');
        xlim(ax{iOGR}, [0,150]);
        grid(ax{iOGR}, 'on')
        ylim(ax{iOGR}, [0,0.2]);

        %
        % plot(ax{iOGR}, [rexo(iOGR), rexo(iOGR)], [0, 1]);
    end



end