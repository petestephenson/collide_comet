function [f, ax, out] = make_electrondens_rAngle_QVec_2d_MC_AllCases_Plot_v2( ...
    peCasesNo, swCasesNo, plasDensDir, opts)

%MAKE_ELECTRONDENS_RANGLE_QVEC_2D_MC_ALLCASES_PLOT_V2
% MC-sampled 2D radius-angle maps of electron density.
%
% Asymmetric cases:
%   x-axis = angle to QVec / highest outgassing direction [deg]
%
% Reference case:
%   x-axis = angle to +X, +Y, and +Z axes [deg], plotted as 3 panels
%
% y-axis = cometocentric radius [km]
% color  = total n_e, cold n_e, or cold / total

    if nargin < 4 || isempty(opts)
        opts = struct();
    end

    if ~isfield(opts, "plotQuantity")
        opts.plotQuantity = "cold";
    end

    if ~isfield(opts, "statistic")
        opts.statistic = "median";
    end

    if ~isfield(opts, "nMC")
        opts.nMC = 5e5;
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
        opts.rMax = 300;
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

    if ~isfield(opts, "rngSeed")
        opts.rngSeed = 1;
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

    opts.plotQuantity = lower(string(opts.plotQuantity));
    opts.statistic = lower(string(opts.statistic));

    if ~ismember(opts.plotQuantity, ["total", "cold", "fraction"])
        error('opts.plotQuantity must be "total", "cold", or "fraction".');
    end

    if ~ismember(opts.statistic, ["median", "mean", "geomean"])
        error('opts.statistic must be "median", "mean", or "geomean".');
    end

    if isempty(opts.plotLog10)
        opts.plotLog10 = opts.plotQuantity ~= "fraction";
    end

    if opts.rMin < 0
        error('opts.rMin must be >= 0.');
    end

    if opts.rMax <= opts.rMin
        error('opts.rMax must be greater than opts.rMin.');
    end

    rng(opts.rngSeed);

    nCases = numel(peCasesNo);

    r_edges = linspace(opts.rMin, opts.rMax, opts.nRadBins + 1);
    angle_edges = linspace(0, 180, opts.nAngleBins + 1);

    r_centres = r_edges(1:end-1) + 0.5 .* diff(r_edges);
    angle_centres = angle_edges(1:end-1) + 0.5 .* diff(angle_edges);

    % --------------------------------------------------------
    % Count panels.
    % Reference cases get 3 panels: +X, +Y, +Z.
    % Asymmetric cases get 1 panel: QVec.
    % --------------------------------------------------------
    nPanels = 0;

    for iCase = 1:nCases
        if iCase <= numel(opts.caseNameList)
            caseStrTmp = opts.caseNameList{iCase};
        else
            caseStrTmp = sprintf('Case %d', iCase);
        end

        isRefCaseTmp = strcmpi(strtrim(string(caseStrTmp)), "ref") || ...
                       contains(lower(string(caseStrTmp)), "ref");

        if isRefCaseTmp
            nPanels = nPanels + 3;
        else
            nPanels = nPanels + 1;
        end
    end

    nRows = 2;
    nCols = ceil(nPanels / nRows);

    f = figure( ...
        'Position', [100, 100, 320 * nCols, 650], ...
        'Color', 'w');

    t = tiledlayout(f, nRows, nCols, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    ax = gobjects(nPanels, 1);
    iPanel = 0;

    out = struct();
    out.opts = opts;
    out.r_edges = r_edges;
    out.r_centres = r_centres;
    out.angle_edges = angle_edges;
    out.angle_centres = angle_centres;

    % ====================================================================
    % Main case loop
    % ====================================================================
    for iCase = 1:nCases

        if iCase <= numel(opts.caseNameList)
            caseStr = opts.caseNameList{iCase};
        else
            caseStr = sprintf('Case %d', iCase);
        end

        isRefCase = strcmpi(strtrim(string(caseStr)), "ref") || ...
                    contains(lower(string(caseStr)), "ref");

        fprintf('Processing MC 2D map for %s (%d/%d)...\n', ...
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

        [plas_dens_asym, plas_dens_cold_asym] = asymPlots.calculatePlasmaDensities( ...
            gridsAsym, ...
            plas_dens_E_Asym, ...
            opts.inc_gens, ...
            opts.minEng, ...
            opts.maxEng);

        [~, ~, ~, QVec] = asymPlots.loadOutgassingProperties(peAsymName);
        QVec = QVec(:);

        if numel(QVec) ~= 3 || norm(QVec) == 0
            warning('Skipping %s: invalid QVec.', caseStr);
            continue
        end

        Qhat = QVec ./ norm(QVec);

        x_edges_grid = gridsAsym{1};
        y_edges_grid = gridsAsym{2};
        z_edges_grid = gridsAsym{3};

        x_centres = x_edges_grid(1:end-1) + 0.5 .* diff(x_edges_grid);
        y_centres = y_edges_grid(1:end-1) + 0.5 .* diff(y_edges_grid);
        z_centres = z_edges_grid(1:end-1) + 0.5 .* diff(z_edges_grid);

        % Match density orientation to meshgrid output:
        % meshgrid(x,y,z) has size [ny, nx, nz].
        dens_total = permute(plas_dens_asym, [2, 1, 3]);
        dens_cold  = permute(plas_dens_cold_asym, [2, 1, 3]);

        expectedSize = [numel(y_centres), numel(x_centres), numel(z_centres)];

        if ~isequal(size(dens_total), expectedSize)
            error( ...
                ['Size mismatch for %s: expected density size %s, ', ...
                 'but total density is %s.'], ...
                caseStr, mat2str(expectedSize), mat2str(size(dens_total)));
        end

        if ~isequal(size(dens_cold), expectedSize)
            error( ...
                ['Size mismatch for %s: expected density size %s, ', ...
                 'but cold density is %s.'], ...
                caseStr, mat2str(expectedSize), mat2str(size(dens_cold)));
        end

        % --------------------------------------------------------
        % MC sample full radial range uniformly in volume
        % --------------------------------------------------------
        [xMC, yMC, zMC, rMC] = generateShellMCPoints( ...
            opts.nMC, opts.rMin, opts.rMax);

        insideGrid = ...
            xMC >= x_edges_grid(1) & xMC <= x_edges_grid(end) & ...
            yMC >= y_edges_grid(1) & yMC <= y_edges_grid(end) & ...
            zMC >= z_edges_grid(1) & zMC <= z_edges_grid(end);

        xMC = xMC(insideGrid);
        yMC = yMC(insideGrid);
        zMC = zMC(insideGrid);
        rMC = rMC(insideGrid);

        nValidDomain = numel(rMC);

        if nValidDomain == 0
            warning('Skipping %s: no MC points lie inside the grid domain.', caseStr);
            continue
        end

        % --------------------------------------------------------
        % Unit position vectors of MC points
        % --------------------------------------------------------
        rSafe = rMC;
        rSafe(rSafe == 0) = NaN;

        xHatMC = xMC ./ rSafe;
        yHatMC = yMC ./ rSafe;
        zHatMC = zMC ./ rSafe;

        % --------------------------------------------------------
        % Nearest grid-cell lookup
        % --------------------------------------------------------
        ix = nearestIndex1D(x_centres, xMC);
        iy = nearestIndex1D(y_centres, yMC);
        iz = nearestIndex1D(z_centres, zMC);

        linIdx = sub2ind(size(dens_total), iy, ix, iz);

        totalSample = dens_total(linIdx);
        coldSample  = dens_cold(linIdx);

        switch opts.plotQuantity
            case "total"
                sampleVals = totalSample;
                quantityLabel = 'Total n_e [cm^{-3}]';

            case "cold"
                sampleVals = coldSample;
                quantityLabel = 'Cold n_e [cm^{-3}]';

            case "fraction"
                sampleVals = coldSample ./ totalSample;
                quantityLabel = 'Cold / total n_e';
        end

        % --------------------------------------------------------
        % Define angle axes.
        %
        % Reference case:
        %   plot angle to +X, +Y, +Z.
        %
        % Asymmetric cases:
        %   plot angle to QVec only.
        % --------------------------------------------------------
        if isRefCase

            angleDirs = [
                1, 0, 0
                0, 1, 0
                0, 0, 1
            ];

            angleNames = {
                '+X axis'
                '+Y axis'
                '+Z axis'
            };

            xLabels = {
                'Angle to +X axis [deg]'
                'Angle to +Y axis [deg]'
                'Angle to +Z axis [deg]'
            };

        else

            angleDirs = Qhat(:).';
            angleNames = {'QVec'};
            xLabels = {'Angle to highest outgassing / QVec [deg]'};

        end

        % --------------------------------------------------------
        % Plot one panel per requested angle direction
        % --------------------------------------------------------
        for iAngle = 1:size(angleDirs, 1)

            angleDir = angleDirs(iAngle, :).';
            angleDir = angleDir ./ norm(angleDir);

            cosAngle = ...
                xHatMC .* angleDir(1) + ...
                yHatMC .* angleDir(2) + ...
                zHatMC .* angleDir(3);

            cosAngle(cosAngle > 1) = 1;
            cosAngle(cosAngle < -1) = -1;

            angleMC = acosd(cosAngle);

            [mapVals, nPerBin] = bin_MC_r_angle_values( ...
                rMC, ...
                angleMC, ...
                sampleVals, ...
                r_edges, ...
                angle_edges, ...
                opts.statistic, ...
                opts.plotQuantity);

            plotVals = mapVals;

            if opts.plotLog10
                plotVals(plotVals <= 0) = NaN;
                plotVals = log10(plotVals);
            end

            iPanel = iPanel + 1;

            ax(iPanel) = nexttile(t, iPanel);
            hold(ax(iPanel), 'on');
            box(ax(iPanel), 'on');

            set(ax(iPanel), ...
                'FontSize', 12, ...
                'LineWidth', 1.1, ...
                'TickDir', 'out');

            imagesc(ax(iPanel), angle_centres, r_centres, plotVals);
            set(ax(iPanel), 'YDir', 'normal');

            xlim(ax(iPanel), [0, 180]);
            xticks(ax(iPanel), 0:30:180);
            ylim(ax(iPanel), [opts.rMin, opts.rMax]);

            xlabel(ax(iPanel), xLabels{iAngle});
            ylabel(ax(iPanel), 'Radius [km]');

            if isRefCase
                title(ax(iPanel), sprintf('%s: %s', caseStr, angleNames{iAngle}), ...
                    'Interpreter', 'none');
            else
                title(ax(iPanel), caseStr, 'Interpreter', 'none');
            end

            if ~isempty(opts.clims)
                clim(ax(iPanel), opts.clims);
            end

            if ~isempty(opts.cmap)
                colormap(ax(iPanel), opts.cmap);
            else
                try
                    colormap(ax(iPanel), redblue);
                catch
                    % Use MATLAB default if redblue is unavailable.
                end
            end

            % --------------------------------------------------------
            % Optional largest-radius boundary
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

                plot(ax(iPanel), angle_centres(goodBoundary), rBoundary(goodBoundary), ...
                    opts.boundaryLineSpec, ...
                    'LineWidth', opts.boundaryLineWidth);

            else
                boundaryThresh = [];
                rBoundary = [];
            end

            out.panels(iPanel).caseName = caseStr;
            out.panels(iPanel).angleName = angleNames{iAngle};
            out.panels(iPanel).angleDir = angleDir;
            out.panels(iPanel).peAsymNo = peAsymNo;
            out.panels(iPanel).swAsymNo = swAsymNo;
            out.panels(iPanel).QVec = QVec;
            out.panels(iPanel).Qhat = Qhat;
            out.panels(iPanel).nMC_requested = opts.nMC;
            out.panels(iPanel).nMC_insideGrid = nValidDomain;
            out.panels(iPanel).mapVals = mapVals;
            out.panels(iPanel).plotVals = plotVals;
            out.panels(iPanel).nPerBin = nPerBin;
            out.panels(iPanel).rBoundary = rBoundary;
            out.panels(iPanel).boundaryThresh = boundaryThresh;

        end

    end

    % Trim unused axes if any case was skipped.
    ax = ax(1:iPanel);

    cb = colorbar(ax(end));
    cb.Layout.Tile = 'east';

    if opts.plotLog10
        cb.Label.String = ['log_{10}(', quantityLabel, ')'];
    else
        cb.Label.String = quantityLabel;
    end

    out.colorbar = cb;

    title(t, sprintf( ...
        ['MC %s electron-density map vs radius and angle. ', ...
         'Reference case uses +X/+Y/+Z axes; asymmetric cases use QVec. N = %d per case'], ...
        char(opts.plotQuantity), opts.nMC), ...
        'Interpreter', 'none');

end


% ========================================================================
function [xMC, yMC, zMC, rMC] = generateShellMCPoints(nMC, rMin, rMax)
%GENERATESHELLMCPOINTS Generate points uniformly in spherical shell volume.
%
% Uniform in volume:
%   r = (rMin^3 + U * (rMax^3 - rMin^3))^(1/3)

    u_r = rand(nMC, 1);
    u_mu = rand(nMC, 1);
    u_phi = rand(nMC, 1);

    rMC = (rMin.^3 + u_r .* (rMax.^3 - rMin.^3)).^(1/3);

    mu = -1 + 2 .* u_mu;
    phi = 2 .* pi .* u_phi;

    sinTheta = sqrt(1 - mu.^2);

    xMC = rMC .* sinTheta .* cos(phi);
    yMC = rMC .* sinTheta .* sin(phi);
    zMC = rMC .* mu;

end


% ========================================================================
function idx = nearestIndex1D(gridCentres, queryVals)
%NEARESTINDEX1D Return nearest index in a monotonic 1-D grid.

    gridCentres = gridCentres(:);
    queryVals = queryVals(:);

    nGrid = numel(gridCentres);

    if nGrid < 2
        idx = ones(size(queryVals));
        return
    end

    idx = interp1( ...
        gridCentres, ...
        1:nGrid, ...
        queryVals, ...
        'nearest', ...
        'extrap');

    idx = round(idx);

    idx(idx < 1) = 1;
    idx(idx > nGrid) = nGrid;

end


% ========================================================================
function [mapVals, nPerBin] = bin_MC_r_angle_values( ...
    rVals, angleVals, sampleVals, r_edges, angle_edges, statistic, plotQuantity)

%BIN_MC_R_ANGLE_VALUES
% Bin MC-sampled values into radius-angle bins.

    rVals = rVals(:);
    angleVals = angleVals(:);
    sampleVals = sampleVals(:);

    valid = isfinite(rVals) & ...
            isfinite(angleVals) & ...
            isfinite(sampleVals);

    switch string(plotQuantity)
        case {"total", "cold"}
            valid = valid & sampleVals > 0;

        case "fraction"
            valid = valid & sampleVals >= 0 & sampleVals <= 1;
    end

    rVals = rVals(valid);
    angleVals = angleVals(valid);
    sampleVals = sampleVals(valid);

    rBin = discretize(rVals, r_edges);
    aBin = discretize(angleVals, angle_edges);

    good = ~isnan(rBin) & ~isnan(aBin);

    rBin = rBin(good);
    aBin = aBin(good);
    sampleVals = sampleVals(good);

    nRadBins = numel(r_edges) - 1;
    nAngleBins = numel(angle_edges) - 1;

    nPerBin = accumarray( ...
        [rBin, aBin], ...
        ones(size(sampleVals)), ...
        [nRadBins, nAngleBins], ...
        @sum, ...
        0);

    switch string(statistic)

        case "median"
            mapVals = accumarray( ...
                [rBin, aBin], ...
                sampleVals, ...
                [nRadBins, nAngleBins], ...
                @median, ...
                NaN);

        case "mean"
            mapVals = accumarray( ...
                [rBin, aBin], ...
                sampleVals, ...
                [nRadBins, nAngleBins], ...
                @mean, ...
                NaN);

        case "geomean"
            logVals = log10(sampleVals);

            logMean = accumarray( ...
                [rBin, aBin], ...
                logVals, ...
                [nRadBins, nAngleBins], ...
                @mean, ...
                NaN);

            mapVals = 10.^logMean;
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