function [f, ax, out] = make_eDepth_QVec_2d_AllCases_Plot( ...
    peCasesNo, swCasesNo, plasDensDir, r0, dr, opts)

%MAKE_EDEPTH_QVEC_2D_ALLCASES_PLOT
% Plot:
%   Tile 1: 1D median electron depth at r0 +/- dr/2 vs angle to max outgassing
%   Tiles 2-6: 2D radius-angle median electron-depth maps for each case.
%
% Each 2D panel shows:
%   x = angle to maximum outgassing direction [deg]
%   y = radius [km]
%   color = log10(median electron depth)
%
% Example:
%   opts = struct();
%   opts.iTau = 3;
%   opts.medianLevel = 0.5;
%   opts.clims = [-2, 0.5];
%   opts.rMax = 300;
%
%   [f, ax, out] = asymPlots.make_eDepth_QVec_2d_AllCases_Plot( ...
%       peCasesNo, swCasesNo, plasDensDir, 50, 10, opts);

    if nargin < 4 || isempty(r0)
        r0 = 50;
    end

    if nargin < 5 || isempty(dr)
        dr = 10;
    end

    r0 = r0(:)';
    dr = dr(:)';

    if isscalar(dr) && numel(r0) > 1
        dr = repmat(dr, size(r0));
    end

    if isscalar(r0) && numel(dr) > 1
        r0 = repmat(r0, size(dr));
    end

    if numel(r0) ~= numel(dr)
        error('r0 and dr must either have the same length or one must be scalar.');
    end

    nShells = numel(r0);


    if nargin < 6 || isempty(opts)
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

    if ~isfield(opts, "cmap")
        opts.cmap = [];
    end

    if ~isfield(opts, "rMax")
        opts.rMax = [];
    end
    if ~isfield(opts, "depthThresh")
        opts.depthThresh = 1;
    end

    if ~isfield(opts, "showDepthBoundary")
        opts.showDepthBoundary = true;
    end

    if ~isfield(opts, "boundaryLineSpec")
        opts.boundaryLineSpec = 'k-';
    end

    if ~isfield(opts, "boundaryLineWidth")
        opts.boundaryLineWidth = 1.0;
    end

    if ~isfield(opts, "showExobase")
        opts.showExobase = true;
    end

    if ~isfield(opts, "exobaseLineSpec")
        opts.exobaseLineSpec = 'k--';
    end

    if ~isfield(opts, "exobaseLineWidth")
        opts.exobaseLineWidth = 1.0;
    end

    nCases = numel(peCasesNo);

    % ----------------------------
    % Figure layout
    % ----------------------------
    f = figure( ...
        'Position', [100, 100, 1500, 800], ...
        'Color', 'w');

    t = tiledlayout(f, 2, 3, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    % Desired tile arrangement:
    %
    %   Tile 1: Reference case
    %   Tile 2: First non-reference case
    %   Tile 3: Second non-reference case
    %   Tile 4: 1D summary panel
    %   Tile 5: Third non-reference case
    %   Tile 6: Fourth non-reference case

    ax = gobjects(nCases + 1, 1);

    summaryTile = 4;
    refTile = 1;
    nonRefTiles = [2, 3, 5, 6];

    caseNameListStr = string(opts.caseNameList);
    refCaseIdx = find(strcmpi(caseNameListStr, "Ref"), 1, "first");

    if isempty(refCaseIdx)
        warning('No reference case named "Ref" found. Using final case as reference.');
        refCaseIdx = nCases;
    end

    nonRefCaseIdx = setdiff(1:nCases, refCaseIdx, 'stable');

    if numel(nonRefCaseIdx) > numel(nonRefTiles)
        error('Too many non-reference cases for the fixed 2 x 3 layout.');
    end

    caseTile = NaN(1, nCases);
    caseTile(refCaseIdx) = refTile;

    for ii = 1:numel(nonRefCaseIdx)
        caseTile(nonRefCaseIdx(ii)) = nonRefTiles(ii);
    end

    ax2DList = gobjects(nCases, 1);
    out = struct();
    out.opts = opts;
    out.r0 = r0;
    out.dr = dr;
    out.nShells = nShells;

    caseColors = lines(nCases);
    shellLineStyles = {'-', '--', ':', '-.'};
    shellMarkers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};

    % Reference-case summary plotting:
    % show electron depth relative to X, Y, and Z using different symbols.
    refVecDims = [1, 2, 3];
    refVecNames = {'X', 'Y', 'Z'};
    refVecMarkers = {'o', 'pentagram', '^'};
    refStyle  = {'-', '--', ':'};
    refColor = [0, 0, 0];
    refAlpha = 0.5;

    if nShells > numel(shellLineStyles)
        warning('More shells than predefined line styles; styles will repeat.');
    end



    % ----------------------------
    % Bottom-left tile: 1D summary panel
    % ----------------------------
    ax(1) = nexttile(t, summaryTile);
    hold(ax(1), 'on');
    box(ax(1), 'on');
    grid(ax(1), 'off');

    set(ax(1), ...
        'FontSize', 12, ...
        'LineWidth', 1.1, ...
        'TickDir', 'out', ...
        'YScale', 'log');

    xlabel(ax(1), 'Angle to max outgassing [deg]');
    ylabel(ax(1), 'Median electron depth');
    % title(ax(1), 'Median electron depth at selected radii');

    xlim(ax(1), [0, 180]);
    xticks(ax(1), 0:30:180);

    % ----------------------------
    % Main case loop
    % ----------------------------
    for iCase = 1:nCases

        if iCase <= numel(opts.caseNameList)
            caseStr = opts.caseNameList{iCase};
        else
            caseStr = sprintf('Case %d', iCase);
        end

        fprintf('Processing electron-depth maps for %s (%d/%d)...\n', ...
            caseStr, iCase, nCases);

        peAsymNo = peCasesNo{iCase};
        swAsymNo = swCasesNo{iCase};

        peAsymName = fullfile(plasDensDir, ['plas_dens', peAsymNo, '.nc']);
        swAsymName = fullfile(plasDensDir, ['plas_dens', swAsymNo, '.nc']);
        [qSpherical, qDayside, qUniform, qVec] = ...
            asymPlots.loadOutgassingProperties(peAsymName);
        [eDepthRadEdges, eDepthGrid, eDepthDist_all] = ...
            load_total_eDepth_dist(peAsymName, swAsymName);

        [iDimCase, axisSignCase] = get_case_asymmetry_dim(iCase, caseStr);

        [medianDepth, rMids, thetaMids] = compute_medianDepth_grid_local( ...
            eDepthRadEdges, ...
            eDepthGrid, ...
            eDepthDist_all, ...
            iDimCase, ...
            opts.iTau, ...
            opts.medianLevel);

        thetaDeg = rad2deg(thetaMids);

        % For negative outgassing axes, make 0 deg point toward max outgassing.
        if axisSignCase < 0
            thetaDeg = 180 - thetaDeg;

            % Keep x-axis increasing from 0 to 180 deg.
            [thetaDeg, sortIdx] = sort(thetaDeg);
            medianDepth = medianDepth(:, sortIdx);
        end
        % ========================================================
        % Add this case to the 1D summary panel
        % ========================================================
        isRefCase = strcmpi(string(caseStr), "Ref");

        shellProfiles = NaN(nShells, numel(thetaDeg));

        if isRefCase

            % Do not plot the reference case here by itself.
            % The reference profiles are instead plotted alongside each
            % asymmetric case using the corresponding asymmetry direction.
            % continue
            tempNo = 0;

        else

            % ----------------------------------------------------
            % First plot the asymmetric case, as before.
            % ----------------------------------------------------
            for iShell = 1:nShells

                shellProfile = median_shell_profile_from_native_theta( ...
                    medianDepth, rMids, r0(iShell), dr(iShell));

                shellProfiles(iShell, :) = shellProfile;

                goodShell = isfinite(shellProfile) & shellProfile > 0;

                lineStyle = shellLineStyles{ ...
                    mod(iShell - 1, numel(shellLineStyles)) + 1};

                markerStyle = shellMarkers{ ...
                    mod(iShell - 1, numel(shellMarkers)) + 1};

                if nShells == 1
                    displayName = sprintf('%s', caseStr);
                else
                    displayName = sprintf('%s, r = %.0f km', ...
                        caseStr, r0(iShell));
                end

                plot(ax(1), ...
                    thetaDeg(goodShell), ...
                    shellProfile(goodShell), ...
                    'LineStyle', lineStyle, ...
                    'Marker', shellMarkers{iCase}, ...
                    'Color', caseColors(iCase, :), ...
                    'MarkerFaceColor', caseColors(iCase, :), ...
                    'LineWidth', 1, ...
                    'MarkerSize', 3, ...
                    'DisplayName', displayName);

            end

            % ----------------------------------------------------
            % Now plot the symmetric/reference case in the same
            % axis direction as this asymmetric case.
            % ----------------------------------------------------
            if ~isempty(refCaseIdx)

                peRefNo = peCasesNo{refCaseIdx};
                swRefNo = swCasesNo{refCaseIdx};

                peRefName = fullfile(plasDensDir, ['plas_dens', peRefNo, '.nc']);
                swRefName = fullfile(plasDensDir, ['plas_dens', swRefNo, '.nc']);

                [eDepthRadEdgesRef, eDepthGridRef, eDepthDistRef_all] = ...
                    load_total_eDepth_dist(peRefName, swRefName);

                [medianDepthRef, rMidsRef, thetaMidsRef] = ...
                    compute_medianDepth_grid_local( ...
                        eDepthRadEdgesRef, ...
                        eDepthGridRef, ...
                        eDepthDistRef_all, ...
                        iDimCase, ...
                        opts.iTau, ...
                        opts.medianLevel);

                thetaDegRef = rad2deg(thetaMidsRef);

                % Use the same angular convention as the asymmetric case.
                % This means 0 deg points along the labelled asymmetric direction,
                % e.g. -X, +Y, +Z, or -Z.
                if axisSignCase < 0
                    thetaDegRef = 180 - thetaDegRef;

                    [thetaDegRef, sortIdxRef] = sort(thetaDegRef);
                    medianDepthRef = medianDepthRef(:, sortIdxRef);
                end

                for iShell = 1:nShells

                    shellProfileRef = median_shell_profile_from_native_theta( ...
                        medianDepthRef, ...
                        rMidsRef, ...
                        r0(iShell), ...
                        dr(iShell));

                    goodShellRef = isfinite(shellProfileRef) & shellProfileRef > 0;

                    lineStyle = shellLineStyles{ ...
                        mod(iShell - 1, numel(shellLineStyles)) + 1};

                    if nShells == 1
                        displayNameRef = sprintf('Ref %s', caseStr);
                    else
                        displayNameRef = sprintf('Ref %s, r = %.0f km', ...
                            caseStr, r0(iShell));
                    end

                    plot(ax(1), ...
                        thetaDegRef(goodShellRef), ...
                        shellProfileRef(goodShellRef), ...
                        'LineStyle', '--', ...
                        'Marker', shellMarkers{iCase}, ...
                        'Color', caseColors(iCase, :), ...
                        'LineWidth', 1.0, ...
                        'markerSize', 3, ...
                        'DisplayName', displayNameRef);

                end
            end
        end
        % % ========================================================
        % % Add this case to the 1D summary panel
        % % ========================================================
        % isRefCase = strcmpi(string(caseStr), "Ref");

        % shellProfiles = NaN(nShells, numel(thetaDeg));

        % if isRefCase

        %     % For the reference case, plot profiles relative to X, Y, Z.
        %     % These are shown in black, with different symbols for each vector.
        %     refShellProfiles = cell(numel(refVecDims), nShells);

        %     for iVec = 1:numel(refVecDims)

        %         iDimRef = refVecDims(iVec);

        %         [medianDepthRef, rMidsRef, thetaMidsRef] = ...
        %             compute_medianDepth_grid_local( ...
        %                 eDepthRadEdges, ...
        %                 eDepthGrid, ...
        %                 eDepthDist_all, ...
        %                 iDimRef, ...
        %                 opts.iTau, ...
        %                 opts.medianLevel);

        %         thetaDegRef = rad2deg(thetaMidsRef);

        %         for iShell = 1:nShells

        %             shellProfileRef = median_shell_profile_from_native_theta( ...
        %                 medianDepthRef, ...
        %                 rMidsRef, ...
        %                 r0(iShell), ...
        %                 dr(iShell));

        %             refShellProfiles{iVec, iShell} = shellProfileRef;

        %             goodShell = isfinite(shellProfileRef) & shellProfileRef > 0;

        %             lineStyle = shellLineStyles{ ...
        %                 mod(iShell - 1, numel(shellLineStyles)) + 1};

        %             if nShells == 1
        %                 displayName = sprintf('Ref: angle to %s', ...
        %                     refVecNames{iVec});
        %             else
        %                 displayName = sprintf('Ref: angle to %s, r = %.0f km', ...
        %                     refVecNames{iVec}, r0(iShell));
        %             end
        %             refAlphaColor = [refColor, refAlpha];
        %             plot(ax(1), ...
        %                 thetaDegRef(goodShell), ...
        %                 shellProfileRef(goodShell), ...
        %                 'LineStyle', refStyle{iVec}, ...
        %                 'Marker', refVecMarkers{iVec}, ...
        %                 'Color', refAlphaColor, ...
        %                 'MarkerFaceColor', refColor, ...
        %                 'LineWidth', 1, ...
        %                 'MarkerSize', 4, ...
        %                 'DisplayName', displayName);

        %         end
        %     end

        % else

        %     % For asymmetric cases, keep the existing behaviour:
        %     % angle is relative to the maximum outgassing direction.
        %     for iShell = 1:nShells

        %         shellProfile = median_shell_profile_from_native_theta( ...
        %             medianDepth, rMids, r0(iShell), dr(iShell));

        %         shellProfiles(iShell, :) = shellProfile;

        %         goodShell = isfinite(shellProfile) & shellProfile > 0;

        %         lineStyle = shellLineStyles{ ...
        %             mod(iShell - 1, numel(shellLineStyles)) + 1};

        %         if nShells == 1
        %             displayName = sprintf('%s', caseStr);
        %         else
        %             displayName = sprintf('%s, r = %.0f km', ...
        %                 caseStr, r0(iShell));
        %         end

        %         plot(ax(1), ...
        %             thetaDeg(goodShell), ...
        %             shellProfile(goodShell), ...
        %             'LineStyle', lineStyle, ...
        %             'Marker', shellMarkers{ ...
        %                 mod(iShell - 1, numel(shellMarkers)) + 1}, ...
        %             'Color', caseColors(iCase, :), ...
        %             'MarkerFaceColor', caseColors(iCase, :), ...
        %             'LineWidth', 1, ...
        %             'MarkerSize', 3, ...
        %             'DisplayName', displayName);

        %     end
        % end

        % ========================================================
        % 2D radius-angle panel for this case
        % ========================================================
        % ax2D = nexttile(t, iCase + 1);
        % ax(iCase + 1) = ax2D;
        ax2D = nexttile(t, caseTile(iCase));
        ax(iCase + 1) = ax2D;
        ax2DList(iCase) = ax2D;
        hold(ax2D, 'on');
        box(ax2D, 'on');

        set(ax2D, ...
            'FontSize', 12, ...
            'LineWidth', 1.1, ...
            'TickDir', 'out');

        plotVals = medianDepth;

        if opts.plotLog10
            plotVals(plotVals <= 0) = NaN;
            plotVals = log10(plotVals);
        end

        imagesc(ax2D, thetaDeg, rMids, plotVals);
        set(ax2D, 'YDir', 'normal');

        % --------------------------------------------------------
        % Overlay estimated exobase distance
        %
        % Asymmetric cases:
        %   qLocal(theta) = qUniform + 4*qDayside*max(cos(theta), 0)
        %
        % Reference case:
        %   use qSpherical as a constant exobase distance.
        % --------------------------------------------------------
        if opts.showExobase

            if isRefCase

                qPlot = qSpherical;

                if numel(qPlot) > 1
                    qPlot = qPlot(1);
                end

                exobaseTheta = [0, 180];

                sigma_val = 6.7e-17; % cm^2
                u_gas = 1.0;         % km/s

                qRadius = 100 .* qPlot .* sigma_val .* 1e-10 ./ (4 .* pi .* u_gas); % km

                exobaseRadius = qRadius .* ones(size(exobaseTheta));

            else

                % thetaDeg is already measured relative to the max outgassing
                % direction, with negative-axis cases flipped earlier.
                cosMax = cosd(thetaDeg);

                qPlot = calc_local_q_from_angle( ...
                    cosMax, ...
                    qUniform, ...
                    qDayside, ...
                    qVec);

                exobaseTheta = thetaDeg;
                % exobaseRadius = qPlot;
                
                sigma_val = 6.7e-17; % cm^2
                u_gas = 1.0; % km/s
                exobaseRadius =  100 * qPlot * sigma_val *1e-10 / (4*pi* u_gas); % km

                exobaseLabel = 'q_{local}';

            end

            goodQ = isfinite(exobaseTheta) & isfinite(exobaseRadius);

            plot(ax2D, ...
                exobaseTheta(goodQ), ...
                exobaseRadius(goodQ), ...
                opts.exobaseLineSpec, ...
                'LineWidth', opts.exobaseLineWidth, ...
                'DisplayName', exobaseLabel);

        end
        % --------------------------------------------------------
        % Overlay largest-radius boundary where medianDepth > threshold
        % --------------------------------------------------------
        if opts.showDepthBoundary

            rBoundary = largest_radius_above_threshold( ...
                medianDepth, ...
                rMids, ...
                opts.depthThresh);

            goodBoundary = isfinite(rBoundary);
            rBoundarySmooth = smoothdata(rBoundary, 'movmedian', 6, 'omitnan');

            % plot(ax2D, thetaDeg(goodBoundary), rBoundarySmooth(goodBoundary), ...
            %     opts.boundaryLineSpec, ...
            %     'LineWidth', opts.boundaryLineWidth);
        end

        xlim(ax2D, [0, 180]);
        xticks(ax2D, 0:30:180);

        if ~isempty(opts.rMax)
            ylim(ax2D, [0, opts.rMax]);
        end

        if opts.plotLog10
            clim(ax2D, opts.clims);
        end

        if ~isempty(opts.cmap)
            colormap(ax2D, opts.cmap);
        else
            try
                colormap(ax2D, redblue);
            catch
                % Use MATLAB default if redblue is unavailable.
            end
        end

        title(ax2D, sprintf('%s', caseStr), ...
            'Interpreter', 'none');

        xlabel(ax2D, 'Angle to max outgassing [deg]');
        ylabel(ax2D, 'Radius [km]');

        out.cases(iCase).caseName = caseStr;
        out.cases(iCase).peAsymNo = peAsymNo;
        out.cases(iCase).swAsymNo = swAsymNo;
        out.cases(iCase).iDim = iDimCase;
        out.cases(iCase).axisSign = axisSignCase;
        out.cases(iCase).rMids = rMids;
        out.cases(iCase).thetaDeg = thetaDeg;
        out.cases(iCase).medianDepth = medianDepth;
        out.cases(iCase).plotVals = plotVals;
        % out.cases(iCase).shellProfile = shellProfile;
        out.cases(iCase).shellProfiles = shellProfiles;
        out.cases(iCase).shell_r0 = r0;
        out.cases(iCase).shell_dr = dr;
        if opts.showExobase
            out.cases(iCase).qSpherical = qSpherical;
            out.cases(iCase).qDayside = qDayside;
            out.cases(iCase).qUniform = qUniform;
            out.cases(iCase).qVec = qVec;
            out.cases(iCase).exobaseTheta = exobaseTheta;
            out.cases(iCase).exobaseRadius = exobaseRadius;
        else
            out.cases(iCase).qSpherical = [];
            out.cases(iCase).qDayside = [];
            out.cases(iCase).qUniform = [];
            out.cases(iCase).qVec = [];
            out.cases(iCase).exobaseTheta = [];
            out.cases(iCase).exobaseRadius = [];
        end
        % if isRefCase
        %     out.cases(iCase).refVecDims = refVecDims;
        %     out.cases(iCase).refVecNames = refVecNames;
        %     out.cases(iCase).refShellProfiles = refShellProfiles;
        % end
        if isRefCase
            out.cases(iCase).refVecDims = refVecDims;
            out.cases(iCase).refVecNames = refVecNames;
            out.cases(iCase).refShellProfiles = [];
        end
        if opts.showDepthBoundary
            out.cases(iCase).rBoundary_depthThresh = rBoundary;
        else
            out.cases(iCase).rBoundary_depthThresh = [];
        end

        out.cases(iCase).depthThresh = opts.depthThresh;

    end

    % # two column legend on the summary panel
    legend(ax(1), ...
        'Location', 'best', ...
        'Interpreter', 'none', ...
        'NumColumns', 2);

    % ----------------------------
    % Force all 2D colour limits to match
    % ----------------------------
    valid2DAxes = ax2DList(isgraphics(ax2DList));

    if opts.plotLog10
        sharedClims = opts.clims;
    else
        % Optional automatic range if plotting linear values.
        allVals = [];

        for iCase = 1:nCases
            if isfield(out.cases, 'plotVals') && ...
                    iCase <= numel(out.cases) && ...
                    isfield(out.cases(iCase), 'plotVals')

                vals = out.cases(iCase).plotVals;
                allVals = [allVals; vals(:)]; %#ok<AGROW>
            end
        end

        allVals = allVals(isfinite(allVals));

        if isempty(allVals)
            sharedClims = [0, 1];
        else
            sharedClims = [min(allVals), max(allVals)];
        end
    end

    for iAx = 1:numel(valid2DAxes)
        clim(valid2DAxes(iAx), sharedClims);
    end
    % ----------------------------
    % Shared colorbar for 2D panels
    % ----------------------------
    % cb = colorbar(ax(end));
    % cb.Layout.Tile = 'east';
    % cb = colorbar(ax2DList(refCaseIdx));
    % cb.Layout.Tile = 'east';
    cb = colorbar(valid2DAxes(1));
    cb.Layout.Tile = 'east';

    if opts.plotLog10
        cb.Label.String = 'log_{10}(median electron depth)';
    else
        cb.Label.String = 'Median electron depth';
    end

    % title(t, 'Median electron depth vs radius and angle to maximum outgassing');

    out.colorbar = cb;

end
% ========================================================================
function [eDepthRadEdges, eDepthGrid, eDepthDist_all] = ...
    load_total_eDepth_dist(peAsymName, swAsymName)

    [eDepthRadEdges_sw, eDepthGrid_sw, eDepthDist_sw] = ...
        asymPlots.load_electron_depth_dist(swAsymName);

    [eDepthRadEdges_pe, eDepthGrid_pe, eDepthDist_pe] = ...
        asymPlots.load_electron_depth_dist(peAsymName);

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

% elDepthDist dimensions:
%   depth x tauIndex x radius x theta x dimension
%
% Output:
%   medianDepth(radius, theta)

    rMids = 0.5 .* (elDepthRadEdges(1:end-1) + elDepthRadEdges(2:end));

    nTheta = size(elDepthDist, 4);

    thetaBounds = acos(linspace(-1, 1, nTheta + 1));
    thetaMids = 0.5 .* (thetaBounds(1:end-1) + thetaBounds(2:end));

    slab = squeeze(elDepthDist(2:end, iTau, :, :, iDim));

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
function [iDim, axisSign] = get_case_asymmetry_dim(iCase, caseStr)


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
            iDim = 1;
            axisSign = -1;

        otherwise
            switch iCase
                case 1
                    iDim = 1;
                    axisSign = -1;

                case 2
                    iDim = 2;
                    axisSign = +1;

                case 3
                    iDim = 3;
                    axisSign = +1;

                case 4
                    iDim = 3;
                    axisSign = -1;

                case 5
                    iDim = 1;
                    axisSign = -1;

                otherwise
                    iDim = 1;
                    axisSign = -1;
            end
    end
end


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


function rBoundary = largest_radius_above_threshold(medianDepth, rMids, depthThresh)
%LARGEST_RADIUS_ABOVE_THRESHOLD
% For each angle bin, find the largest radius where medianDepth > depthThresh.
%
% Inputs:
%   medianDepth : nRadius x nTheta array
%   rMids       : nRadius x 1 or 1 x nRadius radial midpoint array
%   depthThresh : threshold value, e.g. 1
%
% Output:
%   rBoundary   : 1 x nTheta array of outermost threshold-crossing radii

    rMids = rMids(:);

    nTheta = size(medianDepth, 2);

    rBoundary = NaN(1, nTheta);

    for it = 1:nTheta

        vals = medianDepth(:, it);

        good = isfinite(vals) & vals > depthThresh;

        if any(good)
            rBoundary(it) = max(rMids(good));
        end

    end

end

function [Qpts] = calc_local_q_from_angle(cos2Max, qUniform, qDayside, qVec)
%CALC_LOCAL_Q_FROM_ANGLE Calculate local exobase distance.
%
% qVec is included for compatibility with existing function signatures,
% but is not used here because cos2Max has already been supplied.

    %#ok<INUSD>

    cos2Max_flat = cos2Max;
    cos2Max_flat(cos2Max_flat < 0) = 0.0;

    Qpts = qUniform + 4 .* qDayside .* cos2Max_flat;

end