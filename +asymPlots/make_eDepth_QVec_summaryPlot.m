function [f, ax, out] = make_eDepth_QVec_summaryPlot( ...
    peCasesNo, swCasesNo, plasDensDir, r0, dr, nAngleBins, iCase2D, opts)

%MAKE_EDEPTH_QVEC_SUMMARYPLOT
% Summary plots for electron depth:
%
%   Panel 1:
%       median electron depth at r0 +/- dr/2 vs angle to max outgassing
%       for all cases.
%
%   Panel 2:
%       rectangular radius-angle map for one case:
%           x = angle to max outgassing
%           y = radius
%           color = log10(median electron depth)
%
% This function uses elDepthDist in the same format as make_2dElDepthPlot:
%       elDepthDist(depth, tauIndex, radius, theta, iDim)

    if nargin < 4 || isempty(r0)
        r0 = 50;
    end

    if nargin < 5 || isempty(dr)
        dr = 10;
    end

    if nargin < 6 || isempty(nAngleBins)
        nAngleBins = [];
    end

    if nargin < 7 || isempty(iCase2D)
        iCase2D = 5;
    end

    % if nargin < 8 || isempty(iDim)
    %     iDim = 3;
    % end

    if nargin < 8 || isempty(opts)
        opts = struct();
    end

    if ~isfield(opts, "iTau")
        opts.iTau = 3;
    end

    if ~isfield(opts, "medianLevel")
        opts.medianLevel = 0.5;
    end

    if ~isfield(opts, "clims")
        opts.clims = [-2, 0.5];
    end

    if ~isfield(opts, "caseNameList")
        opts.caseNameList = {'-X', '+Y', '+Z', '-Z', 'Ref'};
    end

    if ~isfield(opts, "plotLog10")
        opts.plotLog10 = true;
    end

    nCases = numel(peCasesNo);

    % ----------------------------
    % Figure setup
    % ----------------------------
    f = figure( ...
        'Position', [100, 100, 1250, 520], ...
        'Color', 'w');

    t = tiledlayout(f, 1, 2, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    ax = gobjects(1, 2);

    ax(1) = nexttile(t, 1);
    ax(2) = nexttile(t, 2);

    hold(ax(1), 'on');
    box(ax(1), 'on');
    grid(ax(1), 'on');

    set(ax(1), ...
        'FontSize', 12, ...
        'LineWidth', 1.1, ...
        'TickDir', 'out', ...
        'YScale', 'log');

    set(ax(2), ...
        'FontSize', 12, ...
        'LineWidth', 1.1, ...
        'TickDir', 'out');

    caseColors = lines(nCases);

    out = struct();
    out.r0 = r0;
    out.dr = dr;
    % out.iDim = iDim;
    out.iDimMode = 'case_asymmetry_axis';
    out.iCase2D = iCase2D;
    out.opts = opts;

    % ============================================================
    % Panel 1: median electron depth at fixed radius vs angle
    % ============================================================
    for iCase = 1:nCases

        if iCase <= numel(opts.caseNameList)
            caseStr = opts.caseNameList{iCase};
        else
            caseStr = sprintf('Case %d', iCase);
        end

        [iDimCase, axisSignCase] = get_case_asymmetry_dim(iCase, caseStr);
        fprintf('Processing shell profile for %s (%d/%d)...\n', ...
            caseStr, iCase, nCases);

        peAsymNo = peCasesNo{iCase};
        swAsymNo = swCasesNo{iCase};

        peAsymName = fullfile(plasDensDir, ['plas_dens', peAsymNo, '.nc']);
        swAsymName = fullfile(plasDensDir, ['plas_dens', swAsymNo, '.nc']);

        [eDepthRadEdges, eDepthGrid, eDepthDist_all] = ...
            load_total_eDepth_dist(peAsymName, swAsymName);
        [medianDepth, rMids, thetaMids] = compute_medianDepth_grid_local( ...
            eDepthRadEdges, ...
            eDepthGrid, ...
            eDepthDist_all, ...
            iDimCase, ...
            opts.iTau, ...
            opts.medianLevel);
        out.cases(iCase).iDim = iDimCase;
        out.cases(iCase).axisSign = axisSignCase;
        thetaDeg = rad2deg(thetaMids);
        if axisSignCase < 0
            thetaDeg = 180 - thetaDeg;
        end

        % Optionally re-bin theta if requested.
        if isempty(nAngleBins)
            angle_centres = thetaDeg(:);
            shellProfile = median_shell_profile_from_native_theta( ...
                medianDepth, rMids, r0, dr);
        else
            angle_edges = linspace(0, 180, nAngleBins + 1);
            angle_centres = angle_edges(1:end-1) + 0.5 .* diff(angle_edges);

            nativeProfile = median_shell_profile_from_native_theta( ...
                medianDepth, rMids, r0, dr);

            shellProfile = rebin_angle_profile_median( ...
                thetaDeg(:), nativeProfile(:), angle_edges);
        end

        good = isfinite(shellProfile) & shellProfile > 0;

        plot(ax(1), angle_centres(good), shellProfile(good), 'o-', ...
            'Color', caseColors(iCase, :), ...
            'MarkerFaceColor', caseColors(iCase, :), ...
            'LineWidth', 1.5, ...
            'MarkerSize', 5, ...
            'DisplayName', caseStr);

        out.cases(iCase).caseName = caseStr;
        out.cases(iCase).peAsymNo = peAsymNo;
        out.cases(iCase).swAsymNo = swAsymNo;
        out.cases(iCase).rMids = rMids;
        out.cases(iCase).thetaDeg = thetaDeg;
        out.cases(iCase).medianDepth = medianDepth;
        out.cases(iCase).shell_angle_centres = angle_centres;
        out.cases(iCase).shellProfile = shellProfile;

    end

    xlabel(ax(1), 'Angle to maximum outgassing / QVec [deg]');
    ylabel(ax(1), 'Median electron depth');
    title(ax(1), sprintf('Median electron depth at r = %.1f \\pm %.1f km', ...
        r0, 0.5 * dr));

    xlim(ax(1), [0, 180]);
    xticks(ax(1), 0:30:180);

    legend(ax(1), ...
        'Location', 'best', ...
        'Interpreter', 'none');

    % ============================================================
    % Panel 2: rectangular radius-angle plot for one case
    % ============================================================
    if iCase2D <= numel(opts.caseNameList)
        caseStr2D = opts.caseNameList{iCase2D};
    else
        caseStr2D = sprintf('Case %d', iCase2D);
    end
    [iDim2D, axisSign2D] = get_case_asymmetry_dim(iCase2D, caseStr2D);

    peAsymNo = peCasesNo{iCase2D};
    swAsymNo = swCasesNo{iCase2D};

    peAsymName = fullfile(plasDensDir, ['plas_dens', peAsymNo, '.nc']);
    swAsymName = fullfile(plasDensDir, ['plas_dens', swAsymNo, '.nc']);

    [eDepthRadEdges, eDepthGrid, eDepthDist_all] = ...
        load_total_eDepth_dist(peAsymName, swAsymName);

    [medianDepth2D, rMids2D, thetaMids2D] = compute_medianDepth_grid_local( ...
        eDepthRadEdges, ...
        eDepthGrid, ...
        eDepthDist_all, ...
        iDim2D, ...
        opts.iTau, ...
        opts.medianLevel);

    thetaDeg2D = rad2deg(thetaMids2D);
    if axisSign2D < 0
        thetaDeg2D = 180 - thetaDeg2D;
    end

    out.case2D.iDim = iDim2D;
    out.case2D.axisSign = axisSign2D;

    plotVals = medianDepth2D;

    if opts.plotLog10
        plotVals(plotVals <= 0) = NaN;
        plotVals = log10(plotVals);
    end

    imagesc(ax(2), thetaDeg2D, rMids2D, plotVals);
    set(ax(2), 'YDir', 'normal');

    xlabel(ax(2), 'Angle to maximum outgassing / QVec [deg]');
    ylabel(ax(2), 'Radius [km]');
    % title(ax(2), sprintf('Radius-angle electron depth map, %s', caseStr2D));
    title(ax(2), sprintf('Radius-angle electron depth map, %s, iDim = %d, sign = %+d', ...
        caseStr2D, iDim2D, axisSign2D));    

    xlim(ax(2), [0, 180]);
    xticks(ax(2), 0:30:180);

    cb = colorbar(ax(2));

    if opts.plotLog10
        cb.Label.String = 'log_{10}(median electron depth)';
        clim(ax(2), opts.clims);
    else
        cb.Label.String = 'Median electron depth';
    end

    out.case2D.caseName = caseStr2D;
    out.case2D.rMids = rMids2D;
    out.case2D.thetaDeg = thetaDeg2D;
    out.case2D.medianDepth = medianDepth2D;
    out.case2D.plotVals = plotVals;

    title(t, 'Electron depth summary vs radius and angle');

end


% ========================================================================
function [eDepthRadEdges, eDepthGrid, eDepthDist_all] = ...
    load_total_eDepth_dist(peAsymName, swAsymName)

    [eDepthRadEdges_sw, eDepthGrid_sw, eDepthDist_sw] = ...
        asymPlots.load_electron_depth_dist(swAsymName);

    [eDepthRadEdges_pe, eDepthGrid_pe, eDepthDist_pe] = ...
        asymPlots.load_electron_depth_dist(peAsymName);

    % Basic consistency checks.
    if ~isequal(eDepthRadEdges_sw, eDepthRadEdges_pe)
        warning('SW and PE electron-depth radial grids differ. Using SW radial grid.');
    end

    if ~isequal(eDepthGrid_sw, eDepthGrid_pe)
        warning('SW and PE electron-depth grids differ. Using SW electron-depth grid.');
    end

    eDepthRadEdges = eDepthRadEdges_sw;
    eDepthGrid = eDepthGrid_sw;

    eDepthDist_all = eDepthDist_sw + eDepthDist_pe;

end


% ========================================================================
function [medianDepth, rMids, thetaMids] = compute_medianDepth_grid_local( ...
    elDepthRadEdges, elDepthGrid, elDepthDist, iDim, iTau, medianLevel)

% This is the same basic logic used in make_2dElDepthPlot.
%
% elDepthDist is assumed to have dimensions:
%   depth x tauIndex x radius x theta x dimension
%
% The output medianDepth has dimensions:
%   radius x theta

    rMids = 0.5 .* (elDepthRadEdges(1:end-1) + elDepthRadEdges(2:end));

    nTheta = size(elDepthDist, 4);

    thetaBounds = acos(linspace(-1, 1, nTheta + 1));
    thetaMids = 0.5 .* (thetaBounds(1:end-1) + thetaBounds(2:end));

    slab = squeeze(elDepthDist(2:end, iTau, :, :, iDim));

    % After squeeze, expected size is:
    %   nDepth x nRadius x nTheta
    if ndims(slab) ~= 3
        error('Expected slab to be 3D after squeeze, but size is %s.', ...
            mat2str(size(slab)));
    end

    cumSlab = cumsum(slab, 1);

    denom = cumSlab(end, :, :);
    denom(denom == 0) = NaN;

    frac = cumSlab ./ denom;

    nR = size(frac, 2);
    nTheta = size(frac, 3);

    medianDepth = NaN(nR, nTheta);

    for it = 1:nTheta
        for ir = 1:nR
            ix = find(frac(:, ir, it) >= medianLevel, 1, 'first');

            if ~isempty(ix)
                medianDepth(ir, it) = elDepthGrid(ix + 1);
            end
        end
    end

end


% ========================================================================
function shellProfile = median_shell_profile_from_native_theta( ...
    medianDepth, rMids, r0, dr)

% medianDepth is radius x theta.
% This takes the median across radial bins within r0 +/- dr/2.

    rMask = rMids >= (r0 - 0.5 * dr) & ...
            rMids <  (r0 + 0.5 * dr);

    if ~any(rMask)
        warning('No radial bins found in shell r = %.3f +/- %.3f.', ...
            r0, 0.5 * dr);

        shellProfile = NaN(1, size(medianDepth, 2));
        return
    end

    vals = medianDepth(rMask, :);
    vals(vals <= 0) = NaN;

    shellProfile = median(vals, 1, 'omitnan');

end


% ========================================================================
function yRebinned = rebin_angle_profile_median(thetaDeg, yVals, angle_edges)

    nBins = numel(angle_edges) - 1;

    binIdx = discretize(thetaDeg, angle_edges);

    good = ~isnan(binIdx) & isfinite(yVals) & yVals > 0;

    if ~any(good)
        yRebinned = NaN(nBins, 1);
        return
    end

    yRebinned = accumarray( ...
        binIdx(good), ...
        yVals(good), ...
        [nBins, 1], ...
        @median, ...
        NaN);

end

function [iDim, axisSign] = get_case_asymmetry_dim(iCase, caseStr)
%GET_CASE_ASYMMETRY_DIM
% Return the electron-depth dimension corresponding to the asymmetry axis.
%
% Outputs:
%   iDim     = 1, 2, or 3 for X, Y, or Z
%   axisSign = +1 or -1 for positive/negative asymmetry direction
%
% Assumes:
%   iDim = 1 -> X asymmetry axis
%   iDim = 2 -> Y asymmetry axis
%   iDim = 3 -> Z asymmetry axis

    caseStr = string(caseStr);

    switch caseStr
        case "+X"
            iDim = 1;
            axisSign = +1;

        case "-X"
            iDim = 1;
            axisSign = -1;

        case "+Y"
            iDim = 2;
            axisSign = +1;

        case "-Y"
            iDim = 2;
            axisSign = -1;

        case "+Z"
            iDim = 3;
            axisSign = +1;

        case "-Z"
            iDim = 3;
            axisSign = -1;

        case {"Ref", "REF", "ref"}
            iDim = 3;
            axisSign = +1;

        otherwise
            switch iCase
                case 1
                    iDim = 1;  % -X
                    axisSign = -1;

                case 2
                    iDim = 2;  % +Y
                    axisSign = +1;

                case 3
                    iDim = 3;  % +Z
                    axisSign = +1;

                case 4
                    iDim = 3;  % -Z
                    axisSign = -1;

                case 5
                    iDim = 3;  % Ref
                    axisSign = +1;

                otherwise
                    iDim = 3;
                    axisSign = +1;
            end
    end

end