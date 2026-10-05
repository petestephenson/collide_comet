function [f, ax, out] = make_cold_e_prod_volumeWeighted_Hist_AllCases( ...
    peCasesNo, swCasesNo, plasDensDir, opts)
%MAKE_COLD_E_PROD_VOLUMEWEIGHTED_HIST_ALLCASES
%
% Plot a single volume-weighted histogram of cold-electron production rates
% for each outgassing case.
%
% Each grid cell contributes to the histogram with weight dV = dx*dy*dz.
%
% By default this plots:
%
%   x = log10(P_e [cm^-3 s^-1])
%   y = volume fraction
%
% Usage:
%
%   opts = struct();
%   opts.caseNameList = {'-x', '+y', '+z', '-z', 'sym'};
%   opts.nBins = 60;
%   opts.normalization = 'probability';
%   opts.rRange = [];          % all radii
%   opts.E_cold = 1;           % eV
%
%   [f, ax, out] = asymPlots.make_cold_e_prod_volumeWeighted_Hist_AllCases( ...
%       peCasesNo, swCasesNo, plasDensDir, opts);

    % ------------------------------------------------------------
    % Defaults
    % ------------------------------------------------------------
    if nargin < 4 || isempty(opts)
        opts = struct();
    end

    if ~isfield(opts, 'caseNameList')
        opts.caseNameList = {'-x', '+y', '+z', '-z', 'sym'};
    end

    if ~isfield(opts, 'E_cold')
        opts.E_cold = 1; % eV
    end

    if ~isfield(opts, 'inc_gens')
        opts.inc_gens = [1, 2];
    end

    if ~isfield(opts, 'nBins')
        opts.nBins = 60;
    end

    if ~isfield(opts, 'yScale')
        opts.yScale = 'linear';   % options: 'linear' or 'log'
    end

    if ~isfield(opts, 'yMinLog')
        opts.yMinLog = 1e-6;      % minimum plotted y value for log scale
    end
    if ~isfield(opts, 'binEdges')
        opts.binEdges = [];
    end

    if ~isfield(opts, 'showPercentiles')
        opts.showPercentiles = true;
    end

    if ~isfield(opts, 'percentiles')
        opts.percentiles = [50, 75, 90];
    end

    if ~isfield(opts, 'percentileLineStyles')
        opts.percentileLineStyles = {'-', '--', ':'};
    end
    if ~isfield(opts, 'normalization')
        % Options include:
        %   'probability'  -> volume fraction
        %   'count'        -> volume in km^3 per bin
        %   'pdf'          -> volume-weighted PDF
        opts.normalization = 'probability';
    end

    if ~isfield(opts, 'rRange')
        % [] means use all grid cells.
        % Otherwise use [rMin, rMax] in km.
        opts.rRange = [];
    end

    if ~isfield(opts, 'cometaryOnly')
        opts.cometaryOnly = false;
    end

    if ~isfield(opts, 'plotLog10')
        opts.plotLog10 = true;
    end

    if ~isfield(opts, 'lineWidth')
        opts.lineWidth = 1.8;
    end

    if ~isfield(opts, 'xlim')
        opts.xlim = [];
    end

    if ~isfield(opts, 'ylim')
        opts.ylim = [];
    end

    if ~isfield(opts, 'titleStr')
        opts.titleStr = 'Volume-weighted cold electron production';
    end

    nCases = numel(peCasesNo);

    caseColors = lines(nCases);

    % ------------------------------------------------------------
    % Load cross sections once
    % ------------------------------------------------------------
    [XS, XS_Eng, E_th] = asymPlots.loadXSections();

    % ------------------------------------------------------------
    % First pass: calculate production values and collect global range
    % ------------------------------------------------------------
    out = struct();
    out.opts = opts;
    out.cases = struct([]);

    allPlotVals = [];

    for iCase = 1:nCases

        caseStr = getCaseName(opts.caseNameList, iCase);

        fprintf('Calculating cold electron production histogram data for %s (%d/%d)...\n', ...
            caseStr, iCase, nCases);

        peAsymNo = peCasesNo{iCase};
        swAsymNo = swCasesNo{iCase};

        peAsymName = fullfile(plasDensDir, ['plas_dens', peAsymNo, '.nc']);
        swAsymName = fullfile(plasDensDir, ['plas_dens', swAsymNo, '.nc']);

        [gridsAsym, plas_dens_E_Asym_pe] = asymPlots.loadRunData(peAsymName);
        [~,         plas_dens_E_Asym_sw] = asymPlots.loadRunData(swAsymName);

        if opts.cometaryOnly
            plas_dens_E_Asym_sw(:,:,:,:,1) = 0.0;
        end

        plas_dens_E_Asym = plas_dens_E_Asym_sw + plas_dens_E_Asym_pe;

        [qSpherical, qDayside, qUniform, qVec] = ...
            asymPlots.loadOutgassingProperties(swAsymName);

        isRefCase = strcmpi(string(caseStr), "sym") || ...
                    strcmpi(string(caseStr), "ref") || ...
                    contains(lower(string(caseStr)), "sym") || ...
                    contains(lower(string(caseStr)), "ref");

        if isRefCase
            qDayside = 0.0;
            qUniform = qSpherical;
        end

        coldElecProd = asymPlots.calculateColdeProduction( ...
            plas_dens_E_Asym, ...
            gridsAsym, ...
            opts.E_cold, ...
            qDayside, ...
            qUniform, ...
            qVec, ...
            XS_Eng, ...
            XS, ...
            E_th, ...
            opts.inc_gens);

        [cellVol, rCell] = makeCellVolumeGrid(gridsAsym, size(coldElecProd));

        if ~isequal(size(cellVol), size(coldElecProd))
            if isequal(size(cellVol), size(permute(coldElecProd, [2, 1, 3])))
                coldElecProd = permute(coldElecProd, [2, 1, 3]);
            else
                error( ...
                    ['Size mismatch for %s: coldElecProd is %s but cellVol is %s. ', ...
                     'Check grid orientation.'], ...
                    caseStr, mat2str(size(coldElecProd)), mat2str(size(cellVol)));
            end
        end

        valid = isfinite(coldElecProd) & ...
                isfinite(cellVol) & ...
                coldElecProd > 0 & ...
                cellVol > 0;

        if ~isempty(opts.rRange)
            rMin = opts.rRange(1);
            rMax = opts.rRange(2);

            valid = valid & ...
                    isfinite(rCell) & ...
                    rCell >= rMin & ...
                    rCell < rMax;
        end

        prodVals = coldElecProd(valid);
        volVals = cellVol(valid);

        if opts.plotLog10
            plotVals = log10(prodVals);
        else
            plotVals = prodVals;
        end

        allPlotVals = [allPlotVals; plotVals(:)]; %#ok<AGROW>

        out.cases(iCase).caseName = caseStr;
        out.cases(iCase).peAsymNo = peAsymNo;
        out.cases(iCase).swAsymNo = swAsymNo;
        out.cases(iCase).prodVals = prodVals;
        out.cases(iCase).plotVals = plotVals;
        out.cases(iCase).volVals = volVals;
        out.cases(iCase).totalVolume = sum(volVals, 'omitnan');
        out.cases(iCase).nCellsUsed = numel(prodVals);
        out.cases(iCase).percentiles = opts.percentiles;
        out.cases(iCase).weightedPercentileVals = weightedPercentiles( ...
            plotVals, ...
            volVals, ...
            opts.percentiles);
        

    end
    
    % ------------------------------------------------------------
    % Histogram bin edges
    % ------------------------------------------------------------
    if isempty(opts.binEdges)

        allPlotVals = allPlotVals(isfinite(allPlotVals));

        if isempty(allPlotVals)
            error('No valid cold electron production values found.');
        end

        xMin = floor(min(allPlotVals));
        xMax = ceil(max(allPlotVals));

        if xMin == xMax
            xMin = xMin - 0.5;
            xMax = xMax + 0.5;
        end

        binEdges = linspace(xMin, xMax, opts.nBins + 1);

    else

        binEdges = opts.binEdges;

    end

    out.binEdges = binEdges;
    out.binCentres = 0.5 .* (binEdges(1:end-1) + binEdges(2:end));

    % ------------------------------------------------------------
    % Plot
    % ------------------------------------------------------------
    f = figure( ...
        'Position', [100, 100, 850, 550], ...
        'Color', 'w');

    ax = axes(f);
    hold(ax, 'on');
    box(ax, 'on');

    set(ax, ...
        'FontSize', 13, ...
        'LineWidth', 1.1, ...
        'TickDir', 'out');

    for iCase = 1:nCases

        % histogram(ax, ...
        %     out.cases(iCase).plotVals, ...
        %     binEdges, ...
        %     'Weights', out.cases(iCase).volVals, ...
        %     'Normalization', opts.normalization, ...
        %     'DisplayStyle', 'stairs', ...
        %     'LineWidth', opts.lineWidth, ...
        %     'EdgeColor', caseColors(iCase, :), ...
        %     'DisplayName', out.cases(iCase).caseName);


        % for iCase = 1:nCases

        [binCountsWeighted, binCentres] = weightedHistCountsManual( ...
            out.cases(iCase).plotVals, ...
            out.cases(iCase).volVals, ...
            binEdges, ...
            opts.normalization);
        if strcmpi(opts.yScale, 'log')
            binCountsWeighted(binCountsWeighted <= 0 | ~isfinite(binCountsWeighted)) = NaN;
        end
        stairs(ax, ...
            binCentres, ...
            binCountsWeighted, ...
            'LineWidth', opts.lineWidth, ...
            'Color', caseColors(iCase, :), ...
            'DisplayName', out.cases(iCase).caseName);

        out.cases(iCase).histCounts = binCountsWeighted;
        out.cases(iCase).histBinCentres = binCentres;

        % end
        if opts.showPercentiles

            pctVals = out.cases(iCase).weightedPercentileVals;

            yLimNow = ylim(ax);

            for iPct = 1:numel(opts.percentiles)

                pct = opts.percentiles(iPct);
                pctVal = pctVals(iPct);

                if ~isfinite(pctVal)
                    continue
                end

                lineStyle = opts.percentileLineStyles{ ...
                    mod(iPct - 1, numel(opts.percentileLineStyles)) + 1};

                xline(ax, ...
                    pctVal, ...
                    'LineStyle', lineStyle, ...
                    'Color', caseColors(iCase, :), ...
                    'LineWidth', 1.0, ...
                    'HandleVisibility', 'off');

            end

            ylim(ax, yLimNow);

        end
    end


    
    set(ax, 'YScale', opts.yScale);

    if strcmpi(opts.yScale, 'log')
        ylimVals = ylim(ax);

        if ylimVals(1) <= 0 || ~isfinite(ylimVals(1))
            ylimVals(1) = opts.yMinLog;
        end

        ylim(ax, ylimVals);
    end

    if opts.plotLog10
        xlabel(ax, 'log_{10}(P_e [cm^{-3} s^{-1}])', ...
            'Interpreter', 'tex');
    else
        xlabel(ax, 'P_e [cm^{-3} s^{-1}]', ...
            'Interpreter', 'tex');
        set(ax, 'XScale', 'log');
    end

    switch lower(opts.normalization)
        case 'probability'
            ylabel(ax, 'Volume fraction');
        case 'count'
            ylabel(ax, 'Cell volume per bin [km^3]');
        case 'pdf'
            ylabel(ax, 'Volume-weighted PDF');
        otherwise
            ylabel(ax, opts.normalization, 'Interpreter', 'none');
    end

    title(ax, opts.titleStr, 'Interpreter', 'none');

    if ~isempty(opts.rRange)
        subtitle(ax, sprintf('%.0f \\leq r < %.0f km', ...
            opts.rRange(1), opts.rRange(2)), ...
            'Interpreter', 'tex');
    end

    if ~isempty(opts.xlim)
        xlim(ax, opts.xlim);
    end

    if ~isempty(opts.ylim)
        ylim(ax, opts.ylim);
    end

    legend(ax, ...
        'Location', 'best', ...
        'Interpreter', 'none');

    if opts.showPercentiles

        hPct = gobjects(numel(opts.percentiles), 1);

        for iPct = 1:numel(opts.percentiles)

            lineStyle = opts.percentileLineStyles{ ...
                mod(iPct - 1, numel(opts.percentileLineStyles)) + 1};

            hPct(iPct) = plot(ax, NaN, NaN, ...
                'LineStyle', lineStyle, ...
                'Color', [0, 0, 0], ...
                'LineWidth', 1.0, ...
                'DisplayName', sprintf('%dth percentile', opts.percentiles(iPct)));

        end

    end
    out.figure = f;
    out.ax = ax;

end


% ========================================================================
function caseStr = getCaseName(caseNameList, iCase)

    if iCase <= numel(caseNameList)
        caseStr = caseNameList{iCase};
    else
        caseStr = sprintf('Case %d', iCase);
    end

end


% ========================================================================
function [cellVol, rCell] = makeCellVolumeGrid(gridsAsym, targetSize)
%MAKECELLVOLUMEGRID
%
% Build a cell-volume grid and radius grid from grid edges.
%
% Uses ndgrid so output size is [nx, ny, nz], which commonly matches
% model arrays indexed as (x, y, z).

    x_edges = gridsAsym{1};
    y_edges = gridsAsym{2};
    z_edges = gridsAsym{3};

    dx = diff(x_edges);
    dy = diff(y_edges);
    dz = diff(z_edges);

    x_centres = x_edges(1:end-1) + 0.5 .* dx;
    y_centres = y_edges(1:end-1) + 0.5 .* dy;
    z_centres = z_edges(1:end-1) + 0.5 .* dz;

    [dX, dY, dZ] = ndgrid(dx, dy, dz);
    cellVol = dX .* dY .* dZ;

    [Xc, Yc, Zc] = ndgrid(x_centres, y_centres, z_centres);
    rCell = sqrt(Xc.^2 + Yc.^2 + Zc.^2);

    % If the target uses meshgrid-style [ny, nx, nz], transpose x/y.
    if nargin > 1 && ~isempty(targetSize)

        if ~isequal(size(cellVol), targetSize)

            if isequal(size(permute(cellVol, [2, 1, 3])), targetSize)
                cellVol = permute(cellVol, [2, 1, 3]);
                rCell = permute(rCell, [2, 1, 3]);
            end

        end

    end

end

function [histVals, binCentres] = weightedHistCountsManual(vals, weights, binEdges, normalization)
%WEIGHTEDHISTCOUNTSMANUAL
% Manual weighted histogram for MATLAB versions without histogram weights.
%
% vals    : data values, e.g. log10(P_e)
% weights : weights for each value, e.g. cell volume dV
%
% normalization:
%   'count'       -> weighted volume per bin
%   'probability' -> volume fraction per bin
%   'pdf'         -> volume-weighted probability density

    vals = vals(:);
    weights = weights(:);

    good = isfinite(vals) & isfinite(weights) & weights > 0;

    vals = vals(good);
    weights = weights(good);

    nBins = numel(binEdges) - 1;
    binCentres = 0.5 .* (binEdges(1:end-1) + binEdges(2:end));

    histVals = zeros(nBins, 1);

    if isempty(vals)
        histVals(:) = NaN;
        return
    end

    % Use discretize if available.
    binIdx = discretize(vals, binEdges);

    goodBin = isfinite(binIdx) & binIdx >= 1 & binIdx <= nBins;

    binIdx = binIdx(goodBin);
    weights = weights(goodBin);

    if isempty(binIdx)
        histVals(:) = NaN;
        return
    end

    histVals = accumarray( ...
        binIdx, ...
        weights, ...
        [nBins, 1], ...
        @sum, ...
        0);

    switch lower(normalization)

        case 'count'
            % Leave as total weighted volume per bin.
            % Units are km^3 if weights are cell volumes in km^3.

        case 'probability'
            totalWeight = sum(histVals, 'omitnan');

            if totalWeight > 0
                histVals = histVals ./ totalWeight;
            else
                histVals(:) = NaN;
            end

        case 'pdf'
            totalWeight = sum(histVals, 'omitnan');
            binWidths = diff(binEdges(:));

            if totalWeight > 0
                histVals = histVals ./ totalWeight ./ binWidths;
            else
                histVals(:) = NaN;
            end

        otherwise
            warning('Unknown normalization "%s". Using count.', normalization);

    end

end

function qVals = weightedPercentiles(vals, weights, percentiles)
%WEIGHTEDPERCENTILES Weighted percentiles of vals.
%
% vals        : data values, e.g. log10(P_e)
% weights     : weights, e.g. cell volume dV
% percentiles : values in percent, e.g. [50, 75, 90]

    vals = vals(:);
    weights = weights(:);
    percentiles = percentiles(:);

    good = isfinite(vals) & ...
           isfinite(weights) & ...
           weights > 0;

    vals = vals(good);
    weights = weights(good);

    qVals = NaN(size(percentiles));

    if isempty(vals)
        return
    end

    [valsSorted, sortIdx] = sort(vals);
    weightsSorted = weights(sortIdx);

    cumWeights = cumsum(weightsSorted);
    totalWeight = cumWeights(end);

    if totalWeight <= 0 || ~isfinite(totalWeight)
        return
    end

    cumProb = cumWeights ./ totalWeight;

    for iPct = 1:numel(percentiles)

        targetProb = percentiles(iPct) ./ 100;

        ix = find(cumProb >= targetProb, 1, 'first');

        if ~isempty(ix)
            qVals(iPct) = valsSorted(ix);
        end

    end

end