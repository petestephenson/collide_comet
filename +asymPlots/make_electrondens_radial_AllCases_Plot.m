function [f, ax, prof] = make_electrondens_radial_AllCases_Plot( ...
    peCasesNo, swCasesNo, plasDensDir, rEdges, opts)
%MAKE_ELECTRONDENS_RADIAL_ALLCASES_PLOT
%
% Plot volume-weighted radial profiles for all outgassing cases.
%
% Panels:
%   1. Total electron density vs radius
%   2. Cold electron density vs radius
%   3. Cold / total electron density vs radius
%
% This is a radial-only replacement for the angle-to-QVec MC function.
%
% For each simulation:
%   1. Load total and cold electron density.
%   2. Build cell-volume and radius grids.
%   3. Bin grid cells by radius.
%   4. Calculate volume-weighted mean total density, cold density,
%      and cold electron fraction in each radial bin.
%
% Usage:
%
%   rEdges = 0:5:300;
%
%   [f, ax, prof] = asymPlots.make_electrondens_radial_AllCases_Plot( ...
%       peCasesNo, swCasesNo, plasDensDir, rEdges);
%
% Optional:
%
%   opts = struct();
%   opts.caseNameList = {'-X', '+Y', '+Z', '-Z', 'Ref'};
%   opts.rMax = 300;
%   opts.cometaryOnly = false;
%   opts.plotMarkers = true;
%
%   [f, ax, prof] = asymPlots.make_electrondens_radial_AllCases_Plot( ...
%       peCasesNo, swCasesNo, plasDensDir, 0:5:300, opts);

    % ------------------------------------------------------------
    % Defaults
    % ------------------------------------------------------------
    if nargin < 4 || isempty(rEdges)
        rEdges = 0:5:300; % km
    end

    if nargin < 5 || isempty(opts)
        opts = struct();
    end

    if ~isfield(opts, "caseNameList")
        opts.caseNameList = {'-X', '+Y', '+Z', '-Z', 'Ref'};
    end

    if ~isfield(opts, "cometaryOnly")
        opts.cometaryOnly = false;
    end

    if ~isfield(opts, "inc_gens")
        opts.inc_gens = [1, 2];
    end

    if ~isfield(opts, "minEng")
        opts.minEng = 0;
    end

    if ~isfield(opts, "maxEng")
        opts.maxEng = 1;
    end

    if ~isfield(opts, "rMax")
        opts.rMax = [];
    end

    if ~isfield(opts, "plotMarkers")
        opts.plotMarkers = true;
    end

    if ~isfield(opts, "markerStep")
        opts.markerStep = 3;
    end

    if ~isfield(opts, "lineWidth")
        opts.lineWidth = 1.5;
    end

    if ~isfield(opts, "yScaleDensity")
        opts.yScaleDensity = "log";
    end

    if ~isfield(opts, "yScaleFraction")
        opts.yScaleFraction = "log";
    end

    if ~isfield(opts, "fractionYLim")
        opts.fractionYLim = [1e-3, 1];
    end

    if ~isfield(opts, "showLegend")
        opts.showLegend = true;
    end

    rEdges = rEdges(:).';
    rCentres = 0.5 .* (rEdges(1:end-1) + rEdges(2:end));
    nRBins = numel(rCentres);

    nCases = numel(peCasesNo);

    caseColors = lines(nCases);
    caseMarkers = {'o', 'pentagram', '^', 'v', 's', 'd', 'p', 'h'};
    caseLineStyles = {'-', '-', '-', '-', '--', ':', '-.'};

    % ------------------------------------------------------------
    % Create figure
    % ------------------------------------------------------------
    f = figure( ...
        'Position', [100, 100, 1450, 500], ...
        'Color', 'w');

    t = tiledlayout(f, 1, 3, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    ax = gobjects(1, 3);

    ax(1) = nexttile(t, 1);
    ax(2) = nexttile(t, 2);
    ax(3) = nexttile(t, 3);

    for iAx = 1:3
        hold(ax(iAx), 'on');
        box(ax(iAx), 'on');

        set(ax(iAx), ...
            'FontSize', 12, ...
            'LineWidth', 1.1, ...
            'TickDir', 'out');
    end

    set(ax(1), 'YScale', opts.yScaleDensity);
    set(ax(2), 'YScale', opts.yScaleDensity);
    set(ax(3), 'YScale', opts.yScaleFraction);

    % ------------------------------------------------------------
    % Output storage
    % ------------------------------------------------------------
    prof = struct();

    prof.opts = opts;
    prof.rEdges = rEdges;
    prof.rCentres = rCentres;
    prof.nRBins = nRBins;
    prof.profiles = struct([]);

    % ------------------------------------------------------------
    % Main case loop
    % ------------------------------------------------------------
    for iCase = 1:nCases

        if iCase <= numel(opts.caseNameList)
            caseStr = opts.caseNameList{iCase};
        else
            caseStr = sprintf('Case %d', iCase);
        end

        fprintf('Processing radial electron-density profile for %s (%d/%d)...\n', ...
            caseStr, iCase, nCases);

        peAsymNo = peCasesNo{iCase};
        swAsymNo = swCasesNo{iCase};

        peAsymName = fullfile(plasDensDir, ['plas_dens', peAsymNo, '.nc']);
        swAsymName = fullfile(plasDensDir, ['plas_dens', swAsymNo, '.nc']);

        % --------------------------------------------------------
        % Load and combine PE + SW electron distribution functions
        % --------------------------------------------------------
        [gridsAsym, plas_dens_E_Asym_pe] = ...
            asymPlots.loadRunData(peAsymName);

        [~, plas_dens_E_Asym_sw] = ...
            asymPlots.loadRunData(swAsymName);

        if opts.cometaryOnly
            plas_dens_E_Asym_sw(:,:,:,:,1) = 0.0;
        end

        plas_dens_E_Asym = plas_dens_E_Asym_sw + plas_dens_E_Asym_pe;

        % --------------------------------------------------------
        % Calculate total and cold electron densities
        % --------------------------------------------------------
        [plas_dens_total, plas_dens_cold] = ...
            asymPlots.calculatePlasmaDensities( ...
                gridsAsym, ...
                plas_dens_E_Asym, ...
                opts.inc_gens, ...
                opts.minEng, ...
                opts.maxEng);

        % --------------------------------------------------------
        % Build radius and volume grids
        % --------------------------------------------------------
        [rGrid, cellVol] = makeRadiusVolumeGrid_local( ...
            gridsAsym, ...
            size(plas_dens_total));

        % Make sure density arrays match the radius/volume grids.
        [plas_dens_total, plas_dens_cold] = matchDensityOrientation_local( ...
            plas_dens_total, ...
            plas_dens_cold, ...
            size(rGrid), ...
            caseStr);

        % --------------------------------------------------------
        % Radial volume-weighted profiles
        % --------------------------------------------------------
        [mean_total, mean_cold, frac_cold, volPerBin, nCellsPerBin] = ...
            radialVolumeWeightedProfiles_local( ...
                rGrid, ...
                cellVol, ...
                plas_dens_total, ...
                plas_dens_cold, ...
                rEdges);

        % --------------------------------------------------------
        % Plot styles
        % --------------------------------------------------------
        plotColor = caseColors(iCase, :);

        lineStyle = caseLineStyles{ ...
            mod(iCase - 1, numel(caseLineStyles)) + 1};

        markerStyle = caseMarkers{ ...
            mod(iCase - 1, numel(caseMarkers)) + 1};

        if ~opts.plotMarkers
            markerStyle = 'none';
        end

        markerIdx = 1:opts.markerStep:numel(rCentres);

        % --------------------------------------------------------
        % Plot radial profiles
        % --------------------------------------------------------
        goodTotal = isfinite(mean_total) & mean_total > 0;
        goodCold  = isfinite(mean_cold)  & mean_cold  > 0;
        goodFrac  = isfinite(frac_cold)  & frac_cold  > 0;

        plot(ax(1), ...
            rCentres(goodTotal), ...
            mean_total(goodTotal), ...
            'LineStyle', lineStyle, ...
            'Marker', markerStyle, ...
            'MarkerIndices', markerIdx(markerIdx <= nnz(goodTotal)), ...
            'Color', plotColor, ...
            'MarkerFaceColor', plotColor, ...
            'LineWidth', opts.lineWidth, ...
            'MarkerSize', 5, ...
            'DisplayName', caseStr);

        plot(ax(2), ...
            rCentres(goodCold), ...
            mean_cold(goodCold), ...
            'LineStyle', lineStyle, ...
            'Marker', markerStyle, ...
            'MarkerIndices', markerIdx(markerIdx <= nnz(goodCold)), ...
            'Color', plotColor, ...
            'MarkerFaceColor', plotColor, ...
            'LineWidth', opts.lineWidth, ...
            'MarkerSize', 5, ...
            'DisplayName', caseStr);

        plot(ax(3), ...
            rCentres(goodFrac), ...
            frac_cold(goodFrac), ...
            'LineStyle', lineStyle, ...
            'Marker', markerStyle, ...
            'MarkerIndices', markerIdx(markerIdx <= nnz(goodFrac)), ...
            'Color', plotColor, ...
            'MarkerFaceColor', plotColor, ...
            'LineWidth', opts.lineWidth, ...
            'MarkerSize', 5, ...
            'DisplayName', caseStr);

        % --------------------------------------------------------
        % Store output
        % --------------------------------------------------------
        prof.profiles(iCase).caseName = caseStr;
        prof.profiles(iCase).peAsymNo = peAsymNo;
        prof.profiles(iCase).swAsymNo = swAsymNo;

        prof.profiles(iCase).rCentres = rCentres;
        prof.profiles(iCase).rEdges = rEdges;

        prof.profiles(iCase).mean_total = mean_total;
        prof.profiles(iCase).mean_cold = mean_cold;
        prof.profiles(iCase).frac_cold = frac_cold;

        prof.profiles(iCase).volPerBin = volPerBin;
        prof.profiles(iCase).nCellsPerBin = nCellsPerBin;

    end

    % ------------------------------------------------------------
    % Axis formatting
    % ------------------------------------------------------------
    xlabel(ax(1), 'Radius [km]', 'Interpreter', 'tex');
    xlabel(ax(2), 'Radius [km]', 'Interpreter', 'tex');
    xlabel(ax(3), 'Radius [km]', 'Interpreter', 'tex');

    ylabel(ax(1), 'Total n_e [cm^{-3}]', 'Interpreter', 'tex');
    ylabel(ax(2), 'Cold n_e [cm^{-3}]', 'Interpreter', 'tex');
    ylabel(ax(3), 'Cold / total', 'Interpreter', 'tex');

    title(ax(1), 'Total electron density');
    title(ax(2), 'Cold electron density');
    title(ax(3), 'Cold electron fraction');

    if ~isempty(opts.rMax)
        for iAx = 1:3
            xlim(ax(iAx), [0, opts.rMax]);
        end
    else
        for iAx = 1:3
            xlim(ax(iAx), [rEdges(1), rEdges(end)]);
        end
    end

    if ~isempty(opts.fractionYLim)
        ylim(ax(3), opts.fractionYLim);
    end

    if opts.showLegend
        legend(ax(3), ...
            'Location', 'bestoutside', ...
            'Interpreter', 'none');
    end

end


% ========================================================================
function [rGrid, cellVol] = makeRadiusVolumeGrid_local(gridsAsym, targetSize)
%MAKERADIUSVOLUMEGRID_LOCAL
%
% Build radius and cell-volume grids from the simulation grid edges.
%
% Uses ndgrid, giving arrays with size [nx, ny, nz].

    x_edges = gridsAsym{1};
    y_edges = gridsAsym{2};
    z_edges = gridsAsym{3};

    dx = diff(x_edges);
    dy = diff(y_edges);
    dz = diff(z_edges);

    x_centres = x_edges(1:end-1) + 0.5 .* dx;
    y_centres = y_edges(1:end-1) + 0.5 .* dy;
    z_centres = z_edges(1:end-1) + 0.5 .* dz;

    [Xc, Yc, Zc] = ndgrid(x_centres, y_centres, z_centres);
    [dX, dY, dZ] = ndgrid(dx, dy, dz);

    rGrid = sqrt(Xc.^2 + Yc.^2 + Zc.^2);
    cellVol = dX .* dY .* dZ;

    % If the model uses meshgrid-style [ny, nx, nz], transpose x/y.
    if nargin > 1 && ~isempty(targetSize)

        if ~isequal(size(rGrid), targetSize)

            if isequal(size(permute(rGrid, [2, 1, 3])), targetSize)
                rGrid = permute(rGrid, [2, 1, 3]);
                cellVol = permute(cellVol, [2, 1, 3]);
            end

        end

    end

end


% ========================================================================
function [dens_total, dens_cold] = matchDensityOrientation_local( ...
    dens_total, dens_cold, targetSize, caseStr)
%MATCHDENSITYORIENTATION_LOCAL
%
% Ensure density arrays match the radius/cell-volume grid orientation.

    if isequal(size(dens_total), targetSize) && ...
       isequal(size(dens_cold), targetSize)
        return
    end

    if isequal(size(permute(dens_total, [2, 1, 3])), targetSize) && ...
       isequal(size(permute(dens_cold,  [2, 1, 3])), targetSize)

        dens_total = permute(dens_total, [2, 1, 3]);
        dens_cold  = permute(dens_cold,  [2, 1, 3]);

        return
    end

    error( ...
        ['Size mismatch for %s. total density is %s, cold density is %s, ', ...
         'but target grid size is %s.'], ...
        caseStr, ...
        mat2str(size(dens_total)), ...
        mat2str(size(dens_cold)), ...
        mat2str(targetSize));

end


% ========================================================================
function [mean_total, mean_cold, frac_cold, volPerBin, nCellsPerBin] = ...
    radialVolumeWeightedProfiles_local( ...
        rGrid, cellVol, dens_total, dens_cold, rEdges)
%RADIALVOLUMEWEIGHTEDPROFILES_LOCAL
%
% Calculate volume-weighted radial profiles:
%
%   mean_total(r) = sum(n_total * dV) / sum(dV)
%   mean_cold(r)  = sum(n_cold  * dV) / sum(dV)
%   frac_cold(r)  = mean_cold / mean_total

    rVals = rGrid(:);
    volVals = cellVol(:);
    totalVals = dens_total(:);
    coldVals = dens_cold(:);

    good = isfinite(rVals) & ...
           isfinite(volVals) & ...
           isfinite(totalVals) & ...
           isfinite(coldVals) & ...
           volVals > 0 & ...
           totalVals > 0 & ...
           coldVals >= 0;

    rVals = rVals(good);
    volVals = volVals(good);
    totalVals = totalVals(good);
    coldVals = coldVals(good);

    nBins = numel(rEdges) - 1;

    mean_total = NaN(nBins, 1);
    mean_cold = NaN(nBins, 1);
    frac_cold = NaN(nBins, 1);
    volPerBin = NaN(nBins, 1);
    nCellsPerBin = zeros(nBins, 1);

    if isempty(rVals)
        warning('No valid density cells found for radial binning.');
        return
    end

    [~, ~, binIdx] = histcounts(rVals, rEdges);

    inBin = binIdx > 0;

    binIdx = binIdx(inBin);
    volVals = volVals(inBin);
    totalVals = totalVals(inBin);
    coldVals = coldVals(inBin);

    if isempty(binIdx)
        warning('No cells fell inside the requested radial bins.');
        return
    end

    volPerBin = accumarray( ...
        binIdx, ...
        volVals, ...
        [nBins, 1], ...
        @sum, ...
        NaN);

    totalWeighted = accumarray( ...
        binIdx, ...
        totalVals .* volVals, ...
        [nBins, 1], ...
        @sum, ...
        NaN);

    coldWeighted = accumarray( ...
        binIdx, ...
        coldVals .* volVals, ...
        [nBins, 1], ...
        @sum, ...
        NaN);

    nCellsPerBin = accumarray( ...
        binIdx, ...
        ones(size(binIdx)), ...
        [nBins, 1], ...
        @sum, ...
        0);

    mean_total = totalWeighted ./ volPerBin;
    mean_cold = coldWeighted ./ volPerBin;
    frac_cold = mean_cold ./ mean_total;

    badTotal = ~isfinite(mean_total) | mean_total <= 0 | ...
               ~isfinite(volPerBin) | volPerBin <= 0 | ...
               nCellsPerBin <= 0;

    badCold = ~isfinite(mean_cold) | mean_cold <= 0 | ...
              ~isfinite(volPerBin) | volPerBin <= 0 | ...
              nCellsPerBin <= 0;

    badFrac = ~isfinite(frac_cold) | frac_cold <= 0 | ...
              mean_total <= 0 | nCellsPerBin <= 0;

    mean_total(badTotal) = NaN;
    mean_cold(badCold) = NaN;
    frac_cold(badFrac) = NaN;

end