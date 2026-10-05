function [f, ax, prof] = make_electrondens_axisCone_MC_Plot(peCasesNo, swCasesNo, plasDensDir, plotCase, nMC)
%MAKE_ELECTRONDENS_AXISCONE_MC_PLOT
% Plot Monte-Carlo sampled electron density profiles along +/-X, +/-Y,
% and +/-Z semi-axes for one simulation case.
%
% For each semi-axis, this function:
%
%   1. Generates nMC points uniformly in volume inside a cone.
%   2. Finds the nearest grid-cell centre to each MC point.
%   3. Samples total and cold electron densities at those nearest cells.
%   4. Calculates radial profiles from the sampled values.
%
% Outputs:
%   f    - figure handle
%   ax   - axes handles
%   prof - structure containing sampled radial profiles
%
% Usage:
%   [f, ax, prof] = make_electrondens_axisCone_MC_Plot( ...
%       peCasesNo, swCasesNo, plasDensDir, 1, 2e5);
%
% Notes:
%   - nMC is the number of Monte Carlo points per cone.
%   - The MC sampling is uniform in cone volume, so the resulting radial
%     bin means approximate volume-weighted means.
%   - The nearest-grid lookup assumes gridsAsym{1:3} define monotonic
%     structured grid edges.

    if nargin < 4 || isempty(plotCase)
        plotCase = 1;
    end

    if nargin < 5 || isempty(nMC)
        nMC = 2e5;
    end

    % ----------------------------
    % Plot/settings
    % ----------------------------
    cometaryOnly = true;

    plot_lim = 300;          % km
    nBins = 50;
    coneAngleDeg = 10;       % cone half-angle

    caseNameList = {'-X', '+Y', '+Z', '-Z','Ref'};

    inc_gens = [1, 2];
    minEng = 0;
    maxEng = 1;

    % Optional: set a fixed seed for reproducibility.
    % Comment this out if you want different MC realisations each time.
    rng(1);

    % ----------------------------
    % Semi-axis definitions
    % ----------------------------
    axisLabels = {'+X', '-X', '+Y', '-Y', '+Z', '-Z'};

    axisVecs = [
         1,  0,  0;
        -1,  0,  0;
         0,  1,  0;
         0, -1,  0;
         0,  0,  1;
         0,  0, -1
    ];

    baseColors = lines(3);

    lineColors = [
        baseColors(1, :);
        baseColors(1, :);
        baseColors(2, :);
        baseColors(2, :);
        baseColors(3, :);
        baseColors(3, :)
    ];

    lineStyles = {'-', '--', '-', '--', '-', '--'};

    % ----------------------------
    % Create figure/subplots
    % ----------------------------
    f = figure( ...
        'Position', [100, 100, 1350, 500], ...
        'Color', 'w' ...
    );

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
        grid(ax(iAx), 'on');

        set(ax(iAx), ...
            'FontSize', 12, ...
            'LineWidth', 1.1, ...
            'TickDir', 'out');
    end

    set(ax(1), 'YScale', 'log');
    set(ax(2), 'YScale', 'log');
    set(ax(3), 'YScale', 'linear');

    % ----------------------------
    % Radial bins
    % ----------------------------
    r_edges = linspace(0, plot_lim, nBins + 1);
    r_centres = r_edges(1:end-1) + 0.5 .* diff(r_edges);

    % ----------------------------
    % Load selected case
    % ----------------------------
    peAsymNo = peCasesNo{plotCase};
    swAsymNo = swCasesNo{plotCase};

    peAsymName = fullfile(plasDensDir, ['plas_dens', peAsymNo, '.nc']);
    swAsymName = fullfile(plasDensDir, ['plas_dens', swAsymNo, '.nc']);

    [gridsAsym, plas_dens_E_Asym_pe] = asymPlots.loadRunData(peAsymName);
    [~,         plas_dens_E_Asym_sw] = asymPlots.loadRunData(swAsymName);

    if cometaryOnly
        plas_dens_E_Asym_sw(:,:,:,:,1) = 0.0;
    end

    plas_dens_E_Asym = plas_dens_E_Asym_sw + plas_dens_E_Asym_pe;

    [plas_dens_asym, plas_dens_cold_asym] = asymPlots.calculatePlasmaDensities( ...
        gridsAsym, ...
        plas_dens_E_Asym, ...
        inc_gens, ...
        minEng, ...
        maxEng ...
    );

    % ----------------------------
    % Build cell-centre grids
    % ----------------------------
    x_edges = gridsAsym{1};
    y_edges = gridsAsym{2};
    z_edges = gridsAsym{3};

    dx = diff(x_edges);
    dy = diff(y_edges);
    dz = diff(z_edges);

    x_centres = x_edges(1:end-1) + 0.5 .* dx;
    y_centres = y_edges(1:end-1) + 0.5 .* dy;
    z_centres = z_edges(1:end-1) + 0.5 .* dz;

    % meshgrid ordering gives arrays of size:
    %   [numel(y_centres), numel(x_centres), numel(z_centres)]
    [X_mesh, Y_mesh, Z_mesh] = meshgrid(x_centres, y_centres, z_centres);
    rMesh = sqrt(X_mesh.^2 + Y_mesh.^2 + Z_mesh.^2);

    % Match density orientation to meshgrid output
    dens_total = permute(plas_dens_asym, [2, 1, 3]);
    dens_cold  = permute(plas_dens_cold_asym, [2, 1, 3]);

    if ~isequal(size(rMesh), size(dens_total))
        error( ...
            ['Size mismatch: rMesh is %s but total density is %s. ', ...
             'Check whether permute(plas_dens_asym,[2,1,3]) is correct.'], ...
            mat2str(size(rMesh)), ...
            mat2str(size(dens_total)) ...
        );
    end

    if ~isequal(size(rMesh), size(dens_cold))
        error( ...
            ['Size mismatch: rMesh is %s but cold density is %s. ', ...
             'Check whether permute(plas_dens_cold_asym,[2,1,3]) is correct.'], ...
            mat2str(size(rMesh)), ...
            mat2str(size(dens_cold)) ...
        );
    end

    % ----------------------------
    % Storage for outputs
    % ----------------------------
    prof = struct();

    prof.plotCase = plotCase;
    prof.nMC = nMC;
    prof.coneAngleDeg = coneAngleDeg;
    prof.r_edges = r_edges;
    prof.r_centres = r_centres;
    prof.axisLabels = axisLabels;

    % ----------------------------
    % Main MC cone loop
    % ----------------------------
    for iDir = 1:numel(axisLabels)

        thisLabel = axisLabels{iDir};
        thisAxisVec = axisVecs(iDir, :);

        % Generate MC points uniformly in cone volume.
        [xMC, yMC, zMC, rMC] = generateConeMCPoints( ...
            nMC, ...
            thisAxisVec, ...
            coneAngleDeg, ...
            plot_lim ...
        );

        % Find nearest grid-cell centre for each MC point.
        ix = nearestIndex1D(x_centres, xMC);
        iy = nearestIndex1D(y_centres, yMC);
        iz = nearestIndex1D(z_centres, zMC);

        % Remove MC points whose nearest point is outside the usable grid.
        validIndex = ...
            ix >= 1 & ix <= numel(x_centres) & ...
            iy >= 1 & iy <= numel(y_centres) & ...
            iz >= 1 & iz <= numel(z_centres);

        ix = ix(validIndex);
        iy = iy(validIndex);
        iz = iz(validIndex);
        rSample = rMC(validIndex);

        % Convert nearest subscripts into linear indices.
        %
        % dens_total has size [ny, nx, nz], so subscript order is:
        %   row    -> iy
        %   column -> ix
        %   page   -> iz
        linIdx = sub2ind(size(dens_total), iy, ix, iz);

        totalSample = dens_total(linIdx);
        coldSample  = dens_cold(linIdx);

        % Calculate radial profiles from the sampled points.
        [mean_total, mean_cold, frac_cold, nPerBin] = calcMCProfiles( ...
            rSample, ...
            totalSample, ...
            coldSample, ...
            r_edges, ...
            nBins, ...
            plot_lim ...
        );

        % Store profiles
        safeField = matlab.lang.makeValidName(thisLabel);

        prof.(safeField).label = thisLabel;
        prof.(safeField).axisVec = thisAxisVec;
        prof.(safeField).mean_total = mean_total;
        prof.(safeField).mean_cold = mean_cold;
        prof.(safeField).frac_cold = frac_cold;
        prof.(safeField).nPerBin = nPerBin;
        prof.(safeField).nValidMC = numel(rSample);

        % Plot profiles
        plot(ax(1), r_centres, mean_total, ...
            'LineWidth', 2.2, ...
            'Color', lineColors(iDir, :), ...
            'LineStyle', lineStyles{iDir}, ...
            'DisplayName', sprintf('%s, N=%d', thisLabel, numel(rSample)));

        plot(ax(2), r_centres, mean_cold, ...
            'LineWidth', 2.2, ...
            'Color', lineColors(iDir, :), ...
            'LineStyle', lineStyles{iDir}, ...
            'DisplayName', sprintf('%s, N=%d', thisLabel, numel(rSample)));

        plot(ax(3), r_centres, frac_cold, ...
            'LineWidth', 2.2, ...
            'Color', lineColors(iDir, :), ...
            'LineStyle', lineStyles{iDir}, ...
            'DisplayName', sprintf('%s, N=%d', thisLabel, numel(rSample)));

    end

    % ----------------------------
    % Axis labels/titles
    % ----------------------------
    xlabel(ax(1), 'r [km]', 'Interpreter', 'tex');
    xlabel(ax(2), 'r [km]', 'Interpreter', 'tex');
    xlabel(ax(3), 'r [km]', 'Interpreter', 'tex');

    ylabel(ax(1), 'Total n_e [cm^{-3}]', 'Interpreter', 'tex');
    ylabel(ax(2), 'Cold n_e [cm^{-3}]', 'Interpreter', 'tex');
    ylabel(ax(3), 'Cold / total', 'Interpreter', 'tex');

    title(ax(1), sprintf('Total n_e, MC cones, %d^\\circ', coneAngleDeg), ...
        'Interpreter', 'tex');

    title(ax(2), sprintf('Cold n_e, MC cones, %d^\\circ', coneAngleDeg), ...
        'Interpreter', 'tex');

    title(ax(3), sprintf('Cold fraction, MC cones, %d^\\circ', coneAngleDeg), ...
        'Interpreter', 'tex');

    for iAx = 1:3
        xlim(ax(iAx), [0, plot_lim]);
    end

    ylim(ax(3), [0, 1]);

    legend(ax(3), ...
        'Location', 'bestoutside', ...
        'Interpreter', 'none');

    if plotCase <= numel(caseNameList)
        caseStr = caseNameList{plotCase};
    else
        caseStr = sprintf('Case %d', plotCase);
    end

    title(t, sprintf('%s: MC-sampled electron density along semi-axes', caseStr), ...
        'Interpreter', 'none');

end


function [xMC, yMC, zMC, rMC] = generateConeMCPoints(nMC, axisVec, coneAngleDeg, rMax)
%GENERATECONEMCPOINTS Generate points uniformly in volume inside a cone.
%
% The cone has apex at the origin, extends to radius rMax, and has
% half-angle coneAngleDeg around axisVec.
%
% Uniform volume sampling:
%   r = rMax * U^(1/3)
%
% Uniform solid-angle sampling inside cone:
%   cos(theta) is uniform between cos(coneAngle) and 1
%   phi is uniform between 0 and 2*pi

    axisVec = axisVec(:);
    axisVec = axisVec ./ norm(axisVec);

    cosCone = cosd(coneAngleDeg);

    u_r = rand(nMC, 1);
    u_mu = rand(nMC, 1);
    u_phi = rand(nMC, 1);

    rMC = rMax .* u_r.^(1/3);

    mu = cosCone + (1 - cosCone) .* u_mu;
    theta = acos(mu);
    phi = 2 .* pi .* u_phi;

    % Local cone coordinates, where cone axis is initially +Z.
    xLocal = rMC .* sin(theta) .* cos(phi);
    yLocal = rMC .* sin(theta) .* sin(phi);
    zLocal = rMC .* cos(theta);

    localPts = [xLocal, yLocal, zLocal];

    % Rotate local +Z axis onto requested axisVec.
    R = rotationFromZ(axisVec);

    pts = localPts * R.';

    xMC = pts(:, 1);
    yMC = pts(:, 2);
    zMC = pts(:, 3);

end


function R = rotationFromZ(targetVec)
%ROTATIONFROMZ Rotation matrix that maps [0;0;1] onto targetVec.

    zVec = [0; 0; 1];

    targetVec = targetVec(:);
    targetVec = targetVec ./ norm(targetVec);

    c = dot(zVec, targetVec);

    if c > 1 - 1e-12
        R = eye(3);
        return
    end

    if c < -1 + 1e-12
        % 180 degree rotation about X maps +Z to -Z.
        R = [
            1,  0,  0;
            0, -1,  0;
            0,  0, -1
        ];
        return
    end

    v = cross(zVec, targetVec);
    s = norm(v);

    vx = [
          0,   -v(3),  v(2);
          v(3),  0,   -v(1);
         -v(2), v(1),  0
    ];

    % Rodrigues' rotation formula.
    R = eye(3) + vx + vx * vx .* ((1 - c) ./ s.^2);

end


function idx = nearestIndex1D(gridCentres, queryVals)
%NEARESTINDEX1D Return nearest index in a monotonic 1-D grid.
%
% This avoids needing knnsearch and is much faster for structured grids.
%
% gridCentres must be monotonic increasing.

    gridCentres = gridCentres(:);
    queryVals = queryVals(:);

    nGrid = numel(gridCentres);

    if nGrid < 2
        idx = ones(size(queryVals));
        return
    end

    % Use interp1 to estimate the nearest fractional index.
    fracIdx = interp1( ...
        gridCentres, ...
        1:nGrid, ...
        queryVals, ...
        'nearest', ...
        'extrap' ...
    );

    idx = round(fracIdx);

    % Clamp to grid range.
    idx(idx < 1) = 1;
    idx(idx > nGrid) = nGrid;

end


function [mean_total, mean_cold, frac_cold, nPerBin] = calcMCProfiles( ...
    rSample, totalSample, coldSample, r_edges, nBins, plot_lim)
%CALCMCPROFILES Calculate radial profiles from MC sampled density values.
%
% Since the MC points are uniformly sampled in cone volume, the mean in each
% radial bin approximates the volume-weighted mean density in that cone.

    rSample = rSample(:);
    totalSample = totalSample(:);
    coldSample = coldSample(:);

    valid = isfinite(rSample) & ...
            isfinite(totalSample) & ...
            isfinite(coldSample) & ...
            totalSample > 0 & ...
            coldSample >= 0 & ...
            rSample <= plot_lim;

    rSample = rSample(valid);
    totalSample = totalSample(valid);
    coldSample = coldSample(valid);

    [~, ~, bin_idx] = histcounts(rSample, r_edges);

    inBin = bin_idx > 0;

    bin_idx = bin_idx(inBin);
    totalSample = totalSample(inBin);
    coldSample = coldSample(inBin);

    nPerBin = accumarray( ...
        bin_idx, ...
        ones(size(bin_idx)), ...
        [nBins, 1], ...
        @sum, ...
        0 ...
    );

    total_sum = accumarray( ...
        bin_idx, ...
        totalSample, ...
        [nBins, 1], ...
        @sum, ...
        NaN ...
    );

    cold_sum = accumarray( ...
        bin_idx, ...
        coldSample, ...
        [nBins, 1], ...
        @sum, ...
        NaN ...
    );

    mean_total = total_sum ./ nPerBin;
    mean_cold  = cold_sum  ./ nPerBin;

    frac_cold = mean_cold ./ mean_total;

    badTotal = ~isfinite(mean_total) | mean_total <= 0 | nPerBin <= 0;
    badCold  = ~isfinite(mean_cold)  | mean_cold < 0  | nPerBin <= 0;
    badFrac  = ~isfinite(frac_cold)  | mean_total <= 0 | nPerBin <= 0;

    mean_total(badTotal) = NaN;
    mean_cold(badCold) = NaN;
    frac_cold(badFrac) = NaN;

end