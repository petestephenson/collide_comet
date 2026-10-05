function [f, ax, out] = make_electrondens_rAngle_QVec_2d_AllCases_Plot( ...
    peCasesNo, swCasesNo, plasDensDir, opts)

%MAKE_ELECTRONDENS_RANGLE_QVEC_2D_ALLCASES_PLOT
% Plot 2D radius-angle maps of electron density using Cartesian grid cells.
%
% Each panel shows one outgassing case:
%
%   x-axis = angle to QVec / highest outgassing direction [deg]
%   y-axis = cometocentric radius [km]
%   color  = total n_e, cold n_e, or cold / total
%
% The 2D map is produced by binning Cartesian grid-cell centres into
% radius-angle bins and taking the median or mean value in each bin.
%
% Usage:
%
%   opts = struct();
%   opts.plotQuantity = "cold";       % "total", "cold", or "fraction"
%   opts.statistic = "median";        % "median" or "mean"
%   opts.nRadBins = 80;
%   opts.nAngleBins = 72;
%   opts.rMax = 300;
%   opts.clims = [];                  % auto, or e.g. [-2 4]
%
%   [f, ax, out] = asymPlots.make_electrondens_rAngle_QVec_2d_AllCases_Plot( ...
%       peCasesNo, swCasesNo, plasDensDir, opts);

    if nargin < 4 || isempty(opts)
        opts = struct();
    end

    % ----------------------------
    % Options
    % ----------------------------
    if ~isfield(opts, "plotQuantity")
        opts.plotQuantity = "cold";   % "total", "cold", "fraction"
    end

    if ~isfield(opts, "statistic")
        opts.statistic = "median";    % "median" or "mean"
    end

    if ~isfield(opts, "nRadBins")
        opts.nRadBins = 80;
    end

    if ~isfield(opts, "nAngleBins")
        opts.nAngleBins = 72;
    end

    if ~isfield(opts, "rMin")
        opts.rMin = 0;
    end

    if ~isfield(opts, "rMax")
        opts.rMax = [];
    end

    if ~isfield(opts, "clims")
        opts.clims = [];
    end

    if ~isfield(opts, "plotLog10")
        opts.plotLog10 = [];
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

    if ~isfield(opts, "cmap")
        opts.cmap = [];
    end

    if ~isfield(opts, "showBoundary")
        opts.showBoundary = false;
    end

    if ~isfield(opts, "boundaryThresh")
        opts.boundaryThresh = [];
    end

    if ~isfield(opts, "boundaryLineSpec")
        opts.boundaryLineSpec = 'k-';
    end

    if ~isfield(opts, "boundaryLineWidth")
        opts.boundaryLineWidth = 2.0;
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

    if ~isfield(opts, "exobaseSigma")
        opts.exobaseSigma = 6.7e-17; % cm^2
    end

    if ~isfield(opts, "exobaseUGas")
        opts.exobaseUGas = 1.0; % km/s
    end

    opts.plotQuantity = lower(string(opts.plotQuantity));
    opts.statistic = lower(string(opts.statistic));

    if ~ismember(opts.plotQuantity, ["total", "cold", "fraction"])
        error('opts.plotQuantity must be "total", "cold", or "fraction".');
    end

    if ~ismember(opts.statistic, ["median", "mean"])
        error('opts.statistic must be "median" or "mean".');
    end

    % Default log plotting:
    %   total/cold density -> log10
    %   fraction           -> linear
    if isempty(opts.plotLog10)
        opts.plotLog10 = opts.plotQuantity ~= "fraction";
    end

    nCases = numel(peCasesNo);

    % ----------------------------
    % Figure layout
    % ----------------------------
    f = figure( ...
        'Position', [100, 100, 1450, 650], ...
        'Color', 'w');

    t = tiledlayout(f, 2, 3, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    ax = gobjects(nCases, 1);

    % ------------------------------------------------------------
    % Tile ordering:
    %   Put Ref in the top-left tile.
    %   Remaining cases follow the natural order.
    % ------------------------------------------------------------
    caseNameListStr = string(opts.caseNameList);
    refCaseIdx = find(strcmpi(caseNameListStr, "Ref"), 1, "first");

    if isempty(refCaseIdx)
        warning('No reference case named "Ref" found. Using final case as reference.');
        refCaseIdx = nCases;
    end

    nonRefCaseIdx = setdiff(1:nCases, refCaseIdx, 'stable');

    plotOrder = [refCaseIdx, nonRefCaseIdx];

    caseTile = NaN(1, nCases);

    for iTile = 1:numel(plotOrder)
        caseTile(plotOrder(iTile)) = iTile;
    end

    out = struct();
    out.opts = opts;

    % ----------------------------
    % Main case loop
    % ----------------------------
    for iCase = 1:nCases

        if iCase <= numel(opts.caseNameList)
            caseStr = opts.caseNameList{iCase};
        else
            caseStr = sprintf('Case %d', iCase);
        end

        fprintf('Processing %s (%d/%d)...\n', caseStr, iCase, nCases);

        % ax(iCase) = nexttile(t, iCase);
        ax(iCase) = nexttile(t, caseTile(iCase));
        hold(ax(iCase), 'on');
        box(ax(iCase), 'on');

        set(ax(iCase), ...
            'FontSize', 12, ...
            'LineWidth', 1.1, ...
            'TickDir', 'out');

        peAsymNo = peCasesNo{iCase};
        swAsymNo = swCasesNo{iCase};

        peAsymName = fullfile(plasDensDir, ['plas_dens', peAsymNo, '.nc']);
        swAsymName = fullfile(plasDensDir, ['plas_dens', swAsymNo, '.nc']);

        % --------------------------------------------------------
        % Load and calculate densities
        % --------------------------------------------------------
        [gridsAsym, plas_dens_E_Asym_pe] = asymPlots.loadRunData(peAsymName);
        [~,         plas_dens_E_Asym_sw] = asymPlots.loadRunData(swAsymName);

        if opts.cometaryOnly
            plas_dens_E_Asym_sw(:,:,:,:,1) = 0.0;
        end

        plas_dens_E_Asym = plas_dens_E_Asym_sw + plas_dens_E_Asym_pe;

        [plas_dens_asym, plas_dens_cold_asym] = asymPlots.calculatePlasmaDensities( ...
            gridsAsym, ...
            plas_dens_E_Asym, ...
            opts.inc_gens, ...
            opts.minEng, ...
            opts.maxEng);

        % --------------------------------------------------------
        % Load QVec: highest outgassing direction
        % --------------------------------------------------------
        % [~, ~, ~, QVec] = asymPlots.loadOutgassingProperties(peAsymName);
        [qSpherical, qDayside, qUniform, QVec] = ...
            asymPlots.loadOutgassingProperties(peAsymName);

        QVec = QVec(:);

        if numel(QVec) ~= 3 || norm(QVec) == 0
            warning('Skipping %s: invalid QVec.', caseStr);
            continue
        end

        Qhat = QVec ./ norm(QVec);

        % --------------------------------------------------------
        % Build Cartesian cell-centre grid
        % --------------------------------------------------------
        x_edges = gridsAsym{1};
        y_edges = gridsAsym{2};
        z_edges = gridsAsym{3};

        x_centres = x_edges(1:end-1) + 0.5 .* diff(x_edges);
        y_centres = y_edges(1:end-1) + 0.5 .* diff(y_edges);
        z_centres = z_edges(1:end-1) + 0.5 .* diff(z_edges);

        [X, Y, Z] = meshgrid(x_centres, y_centres, z_centres);

        rGrid = sqrt(X.^2 + Y.^2 + Z.^2);

        rSafe = rGrid;
        rSafe(rSafe == 0) = NaN;

        cosAngle = ...
            (X .* Qhat(1) + Y .* Qhat(2) + Z .* Qhat(3)) ./ rSafe;

        cosAngle(cosAngle > 1) = 1;
        cosAngle(cosAngle < -1) = -1;

        angleGrid = acosd(cosAngle);

        % --------------------------------------------------------
        % Match density orientation to meshgrid output.
        % meshgrid gives size [ny, nx, nz].
        % --------------------------------------------------------
        dens_total = permute(plas_dens_asym, [2, 1, 3]);
        dens_cold  = permute(plas_dens_cold_asym, [2, 1, 3]);

        if ~isequal(size(rGrid), size(dens_total))
            error( ...
                ['Size mismatch for %s: rGrid is %s but total density is %s. ', ...
                 'Check whether permute(plas_dens_asym,[2,1,3]) is correct.'], ...
                caseStr, mat2str(size(rGrid)), mat2str(size(dens_total)));
        end

        if ~isequal(size(rGrid), size(dens_cold))
            error( ...
                ['Size mismatch for %s: rGrid is %s but cold density is %s. ', ...
                 'Check whether permute(plas_dens_cold_asym,[2,1,3]) is correct.'], ...
                caseStr, mat2str(size(rGrid)), mat2str(size(dens_cold)));
        end

        % --------------------------------------------------------
        % Choose plotted quantity
        % --------------------------------------------------------
        switch opts.plotQuantity
            case "total"
                valGrid = dens_total;
                quantityLabel = 'Total n_e [cm^{-3}]';

            case "cold"
                valGrid = dens_cold;
                quantityLabel = 'Cold n_e [cm^{-3}]';

            case "fraction"
                valGrid = dens_cold ./ dens_total;
                quantityLabel = 'Cold / total n_e';
        end

        % --------------------------------------------------------
        % Radius bin limits
        % --------------------------------------------------------
        if isempty(opts.rMax)
            rMaxUse = max(rGrid(:), [], 'omitnan');
        else
            rMaxUse = opts.rMax;
        end

        r_edges = linspace(opts.rMin, rMaxUse, opts.nRadBins + 1);
        angle_edges = linspace(0, 180, opts.nAngleBins + 1);

        r_centres = r_edges(1:end-1) + 0.5 .* diff(r_edges);
        angle_centres = angle_edges(1:end-1) + 0.5 .* diff(angle_edges);

        % --------------------------------------------------------
        % Bin Cartesian cells into radius-angle histogram
        % --------------------------------------------------------
        [mapVals, nPerBin] = bin_cartesian_r_angle_values( ...
            rGrid, ...
            angleGrid, ...
            valGrid, ...
            r_edges, ...
            angle_edges, ...
            opts.statistic, ...
            opts.plotQuantity);

        plotVals = mapVals;

        if opts.plotLog10
            plotVals(plotVals <= 0) = NaN;
            plotVals = log10(plotVals);
        end

        imagesc(ax(iCase), angle_centres, r_centres, plotVals);
        set(ax(iCase), 'YDir', 'normal');

        % --------------------------------------------------------
        % Overlay estimated exobase distance
        %
        % Reference case:
        %   use qSpherical as a constant horizontal line.
        %
        % Asymmetric cases:
        %   qLocal(theta) = qUniform + 4*qDayside*max(cos(theta), 0)
        %   then convert q to radius using the same scaling.
        % --------------------------------------------------------
        isRefCase = strcmpi(string(caseStr), "Ref");

        if opts.showExobase

            sigma_val = opts.exobaseSigma; % cm^2
            u_gas = opts.exobaseUGas;      % km/s

            if isRefCase

                qPlot = qSpherical;

                if numel(qPlot) > 1
                    qPlot = qPlot(1);
                end

                exobaseTheta = [0, 180];

                qRadius = 100 .* qPlot .* sigma_val .* 1e-10 ./ ...
                    (4 .* pi .* u_gas); % km

                exobaseRadius = qRadius .* ones(size(exobaseTheta));

            else

                % angle_centres are already measured relative to QVec / max
                % outgassing direction.
                cosMax = cosd(angle_centres);

                qLocal = calc_local_q_from_angle( ...
                    cosMax, ...
                    qUniform, ...
                    qDayside, ...
                    QVec);

                exobaseTheta = angle_centres;

                exobaseRadius = 100 .* qLocal .* sigma_val .* 1e-10 ./ ...
                    (4 .* pi .* u_gas); % km

            end

            goodExobase = isfinite(exobaseTheta) & isfinite(exobaseRadius);

            plot(ax(iCase), ...
                exobaseTheta(goodExobase), ...
                exobaseRadius(goodExobase), ...
                opts.exobaseLineSpec, ...
                'LineWidth', opts.exobaseLineWidth);

        end

        xlim(ax(iCase), [0, 180]);
        xticks(ax(iCase), 0:30:180);

        ylim(ax(iCase), [opts.rMin, rMaxUse]);

        xlabel(ax(iCase), 'Angle to highest outgassing / QVec [deg]');
        ylabel(ax(iCase), 'Radius [km]');

        title(ax(iCase), sprintf('%s', caseStr), ...
            'Interpreter', 'none');

        if ~isempty(opts.clims)
            clim(ax(iCase), opts.clims);
        end

        if ~isempty(opts.cmap)
            colormap(ax(iCase), opts.cmap);
        else
            try
                colormap(ax(iCase), redblue);
            catch
                % Use default colormap if redblue is unavailable.
            end
        end

        % --------------------------------------------------------
        % Optional boundary line
        % --------------------------------------------------------
        if opts.showBoundary

            if isempty(opts.boundaryThresh)
                if opts.plotQuantity == "fraction"
                    boundaryThresh = 0.1;
                else
                    boundaryThresh = 1;
                end
            else
                boundaryThresh = opts.boundaryThresh;
            end

            rBoundary = largest_radius_above_threshold( ...
                mapVals, ...
                r_centres, ...
                boundaryThresh);

            goodBoundary = isfinite(rBoundary);

            plot(ax(iCase), angle_centres(goodBoundary), rBoundary(goodBoundary), ...
                opts.boundaryLineSpec, ...
                'LineWidth', opts.boundaryLineWidth);

        else
            rBoundary = [];
            boundaryThresh = [];
        end

        % --------------------------------------------------------
        % Store outputs
        % --------------------------------------------------------
        out.cases(iCase).caseName = caseStr;
        out.cases(iCase).peAsymNo = peAsymNo;
        out.cases(iCase).swAsymNo = swAsymNo;
        out.cases(iCase).QVec = QVec;
        out.cases(iCase).Qhat = Qhat;
        out.cases(iCase).r_edges = r_edges;
        out.cases(iCase).r_centres = r_centres;
        out.cases(iCase).angle_edges = angle_edges;
        out.cases(iCase).angle_centres = angle_centres;
        out.cases(iCase).mapVals = mapVals;
        out.cases(iCase).plotVals = plotVals;
        out.cases(iCase).nPerBin = nPerBin;
        out.cases(iCase).rBoundary = rBoundary;
        out.cases(iCase).boundaryThresh = boundaryThresh;
        out.cases(iCase).qSpherical = qSpherical;
        out.cases(iCase).qDayside = qDayside;
        out.cases(iCase).qUniform = qUniform;

        if opts.showExobase
            out.cases(iCase).exobaseTheta = exobaseTheta;
            out.cases(iCase).exobaseRadius = exobaseRadius;
        else
            out.cases(iCase).exobaseTheta = [];
            out.cases(iCase).exobaseRadius = [];
        end
    end

    % ----------------------------
    % Shared colorbar
    % ----------------------------
    % cb = colorbar(ax(end));
    cb = colorbar(ax(refCaseIdx));
    cb.Layout.Tile = 'east';

    if opts.plotLog10
        cb.Label.String = ['log_{10}(', quantityLabel, ')'];
    else
        cb.Label.String = quantityLabel;
    end

    out.colorbar = cb;

    title(t, sprintf('%s electron-density map vs radius and angle to QVec', ...
        upperFirst(char(opts.plotQuantity))), ...
        'Interpreter', 'none');

end


% ========================================================================
function [mapVals, nPerBin] = bin_cartesian_r_angle_values( ...
    rGrid, angleGrid, valGrid, r_edges, angle_edges, statistic, plotQuantity)

%BIN_CARTESIAN_R_ANGLE_VALUES
% Bin Cartesian cell-centre values into radius-angle bins.

    rVals = rGrid(:);
    aVals = angleGrid(:);
    vVals = valGrid(:);

    valid = isfinite(rVals) & ...
            isfinite(aVals) & ...
            isfinite(vVals);

    switch string(plotQuantity)
        case {"total", "cold"}
            valid = valid & vVals > 0;

        case "fraction"
            valid = valid & vVals >= 0 & vVals <= 1;
    end

    rVals = rVals(valid);
    aVals = aVals(valid);
    vVals = vVals(valid);

    rBin = discretize(rVals, r_edges);
    aBin = discretize(aVals, angle_edges);

    good = ~isnan(rBin) & ~isnan(aBin);

    rBin = rBin(good);
    aBin = aBin(good);
    vVals = vVals(good);

    nRadBins = numel(r_edges) - 1;
    nAngleBins = numel(angle_edges) - 1;

    nPerBin = accumarray( ...
        [rBin, aBin], ...
        ones(size(vVals)), ...
        [nRadBins, nAngleBins], ...
        @sum, ...
        0);

    switch string(statistic)
        case "median"
            mapVals = accumarray( ...
                [rBin, aBin], ...
                vVals, ...
                [nRadBins, nAngleBins], ...
                @median, ...
                NaN);

        case "mean"
            mapVals = accumarray( ...
                [rBin, aBin], ...
                vVals, ...
                [nRadBins, nAngleBins], ...
                @mean, ...
                NaN);
    end

end


% ========================================================================
function rBoundary = largest_radius_above_threshold(mapVals, r_centres, thresh)

%LARGEST_RADIUS_ABOVE_THRESHOLD
% For each angle bin, find largest radius where mapVals > thresh.

    r_centres = r_centres(:);

    nAngle = size(mapVals, 2);
    rBoundary = NaN(1, nAngle);

    for ia = 1:nAngle

        vals = mapVals(:, ia);

        good = isfinite(vals) & vals > thresh;

        if any(good)
            rBoundary(ia) = max(r_centres(good));
        end

    end

end


% ========================================================================
function strOut = upperFirst(strIn)

    if isempty(strIn)
        strOut = strIn;
        return
    end

    strOut = strIn;
    strOut(1) = upper(strOut(1));

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