function [f, ax, prof] = make_electrondens_axisCone_MC_scaledExobase_Plot( ...
    peCasesNo, swCasesNo, plasDensDir, plotCase, nMC, scaleByLocalExobase)
%MAKE_ELECTRONDENS_AXISCONE_MC_SCALEDEXOBASE_PLOT
% Plot Monte-Carlo sampled electron density profiles along +/-X, +/-Y,
% and +/-Z semi-axes for one simulation case.
%
% Optional distance scaling:
%
%   xProfile = r
%
% or
%
%   xProfile = r / qLocal
%
% where qLocal is calculated from the local direction relative to qVec:
%
%   qLocal = qUniform + 4*qDayside*max(cos(theta_Q), 0)
%
% This follows the form:
%
%   Qpts = qUniform + 4*qDayside*cos2Max_flat
%
% with negative values flattened to zero.
%
% Usage:
%
%   [f, ax, prof] = make_electrondens_axisCone_MC_scaledExobase_Plot( ...
%       peCasesNo, swCasesNo, plasDensDir, 1, 2e5, true);

    if nargin < 4 || isempty(plotCase)
        plotCase = 1;
    end

    if nargin < 5 || isempty(nMC)
        nMC = 2e5;
    end

    if nargin < 6 || isempty(scaleByLocalExobase)
        scaleByLocalExobase = false;
    end

    % ----------------------------
    % Plot/settings
    % ----------------------------
    cometaryOnly = true;

    plot_lim = 300;          % physical max radius sampled [km]
    nBins = 50;
    coneAngleDeg = 10;       % cone half-angle

    caseNameList = {'-X', '+Y', '+Z', '-Z','Ref'};

    inc_gens = [1, 2];
    minEng = 0;
    maxEng = 1;

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
    % Load outgassing properties for local exobase scaling
    % ----------------------------
    [qSpherical, qDayside, qUniform, qVec] = asymPlots.loadOutgassingProperties(peAsymName);

    qVec = qVec(:);

    if numel(qVec) ~= 3
        error('qVec must have 3 components, but has size %s.', mat2str(size(qVec)));
    end

    if norm(qVec) == 0
        error('qVec has zero magnitude.');
    end

    qHat = qVec ./ norm(qVec);

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

    [X_mesh, Y_mesh, Z_mesh] = meshgrid(x_centres, y_centres, z_centres);
    rMesh = sqrt(X_mesh.^2 + Y_mesh.^2 + Z_mesh.^2);

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
    % Pre-generate MC points for all cones, so we can determine bin range
    % before plotting if using scaled distances.
    % ----------------------------
    mc = struct();

    allProfileDist = [];

    for iDir = 1:numel(axisLabels)

        thisAxisVec = axisVecs(iDir, :);

        [xMC, yMC, zMC, rMC] = generateConeMCPoints( ...
            nMC, ...
            thisAxisVec, ...
            coneAngleDeg, ...
            plot_lim ...
        );

        insideGrid = ...
            xMC >= x_edges(1) & xMC <= x_edges(end) & ...
            yMC >= y_edges(1) & yMC <= y_edges(end) & ...
            zMC >= z_edges(1) & zMC <= z_edges(end);

        xMC = xMC(insideGrid);
        yMC = yMC(insideGrid);
        zMC = zMC(insideGrid);
        rMC = rMC(insideGrid);

        qLocal = calcLocalQFromPosition(xMC, yMC, zMC, qUniform, qDayside, qHat);
        

        if scaleByLocalExobase
            sigma_val = 6.7e-17; % cm^2
            u_gas = 1.0; % km/s
            rexo = 100 * qLocal * sigma_val *1e-10 / (4*pi* u_gas); % km
        
            profileDist = rMC ./ rexo;
            % print the average exobase distance for this cone
            fprintf('Average exobase distance for cone %s: %.2f km\n', axisLabels{iDir}, mean(rexo));
        else
            profileDist = rMC;
        end

        mc(iDir).xMC = xMC;
        mc(iDir).yMC = yMC;
        mc(iDir).zMC = zMC;
        mc(iDir).rMC = rMC;
        mc(iDir).qLocal = qLocal;
        mc(iDir).profileDist = profileDist;

        allProfileDist = [allProfileDist; profileDist(:)];

    end

    validProfileDist = isfinite(allProfileDist) & allProfileDist >= 0;
    allProfileDist = allProfileDist(validProfileDist);

    if isempty(allProfileDist)
        error('No valid MC profile distances were generated.');
    end

    if scaleByLocalExobase
        profile_lim = max(allProfileDist);
        xLabelStr = 'r / q_{local}';
    else
        profile_lim = plot_lim;
        xLabelStr = 'r [km]';
    end

    profile_edges = linspace(0, profile_lim, nBins + 1);
    profile_centres = profile_edges(1:end-1) + 0.5 .* diff(profile_edges);

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
    % Storage for outputs
    % ----------------------------
    prof = struct();

    prof.plotCase = plotCase;
    prof.nMC = nMC;
    prof.coneAngleDeg = coneAngleDeg;
    prof.scaleByLocalExobase = scaleByLocalExobase;
    prof.qSpherical = qSpherical;
    prof.qUniform = qUniform;
    prof.qDayside = qDayside;
    prof.qVec = qVec;
    prof.qHat = qHat;
    prof.profile_edges = profile_edges;
    prof.profile_centres = profile_centres;
    prof.axisLabels = axisLabels;

    % ----------------------------
    % Main MC cone loop
    % ----------------------------
    for iDir = 1:numel(axisLabels)

        thisLabel = axisLabels{iDir};
        thisAxisVec = axisVecs(iDir, :);

        xMC = mc(iDir).xMC;
        yMC = mc(iDir).yMC;
        zMC = mc(iDir).zMC;
        rMC = mc(iDir).rMC;
        qLocal = mc(iDir).qLocal;
        profileDist = mc(iDir).profileDist;

        ix = nearestIndex1D(x_centres, xMC);
        iy = nearestIndex1D(y_centres, yMC);
        iz = nearestIndex1D(z_centres, zMC);

        validIndex = ...
            ix >= 1 & ix <= numel(x_centres) & ...
            iy >= 1 & iy <= numel(y_centres) & ...
            iz >= 1 & iz <= numel(z_centres) & ...
            isfinite(profileDist) & ...
            profileDist >= 0 & ...
            qLocal > 0;

        ix = ix(validIndex);
        iy = iy(validIndex);
        iz = iz(validIndex);
        rSample = rMC(validIndex);
        qSample = qLocal(validIndex);
        profileSample = profileDist(validIndex);

        linIdx = sub2ind(size(dens_total), iy, ix, iz);

        totalSample = dens_total(linIdx);
        coldSample  = dens_cold(linIdx);

        [mean_total, mean_cold, frac_cold, nPerBin] = calcMCProfiles( ...
            profileSample, ...
            totalSample, ...
            coldSample, ...
            profile_edges, ...
            nBins, ...
            profile_lim ...
        );

        safeField = matlab.lang.makeValidName(thisLabel);

        prof.(safeField).label = thisLabel;
        prof.(safeField).axisVec = thisAxisVec;
        prof.(safeField).mean_total = mean_total;
        prof.(safeField).mean_cold = mean_cold;
        prof.(safeField).frac_cold = frac_cold;
        prof.(safeField).nPerBin = nPerBin;
        prof.(safeField).nValidMC = numel(profileSample);
        prof.(safeField).rSample = rSample;
        prof.(safeField).qSample = qSample;
        prof.(safeField).profileSample = profileSample;

        plot(ax(1), profile_centres, mean_total, ...
            'LineWidth', 2.2, ...
            'Color', lineColors(iDir, :), ...
            'LineStyle', lineStyles{iDir}, ...
            'DisplayName', sprintf('%s, N=%d', thisLabel, numel(profileSample)));

        plot(ax(2), profile_centres, mean_cold, ...
            'LineWidth', 2.2, ...
            'Color', lineColors(iDir, :), ...
            'LineStyle', lineStyles{iDir}, ...
            'DisplayName', sprintf('%s, N=%d', thisLabel, numel(profileSample)));

        plot(ax(3), profile_centres, frac_cold, ...
            'LineWidth', 2.2, ...
            'Color', lineColors(iDir, :), ...
            'LineStyle', lineStyles{iDir}, ...
            'DisplayName', sprintf('%s, N=%d', thisLabel, numel(profileSample)));

    end

    % ----------------------------
    % Axis labels/titles
    % ----------------------------
    xlabel(ax(1), xLabelStr, 'Interpreter', 'tex');
    xlabel(ax(2), xLabelStr, 'Interpreter', 'tex');
    xlabel(ax(3), xLabelStr, 'Interpreter', 'tex');

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
        xlim(ax(iAx), [0, profile_lim]);
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

    if scaleByLocalExobase
        titleStr = sprintf( ...
            '%s: MC electron density along semi-axes, distance scaled by q_{local}', ...
            caseStr);
    else
        titleStr = sprintf( ...
            '%s: MC electron density along semi-axes', ...
            caseStr);
    end

    title(t, titleStr, 'Interpreter', 'tex');

end


function qLocal = calcLocalQFromPosition(x, y, z, qUniform, qDayside, qHat)
%CALCLOCALQFROMPOSITION Calculate local exobase distance from position.
%
% Uses:
%
%   qLocal = qUniform + 4*qDayside*max(cos(theta_Q), 0)
%
% where theta_Q is the angle between the local radial direction and qVec.
%
% This mirrors the supplied form:
%
%   Qpts = qUniform + 4*qDayside*cos2Max_flat
%
% with negative values set to zero.

    r = sqrt(x.^2 + y.^2 + z.^2);

    rSafe = r;
    rSafe(rSafe == 0) = NaN;

    xHat = x ./ rSafe;
    yHat = y ./ rSafe;
    zHat = z ./ rSafe;

    cosToQ = ...
        xHat .* qHat(1) + ...
        yHat .* qHat(2) + ...
        zHat .* qHat(3);

    cosToQ(cosToQ < 0) = 0.0;

    qLocal = qUniform + 4 .* qDayside .* cosToQ;

end


function [Qpts] = calc_local_q_from_angle(cos2Max, qUniform, qDayside, qVec)
%CALC_LOCAL_Q_FROM_ANGLE Calculate local exobase distance.
%
% qVec is included for compatibility with your existing function signature,
% but is not used here because cos2Max has already been supplied.

    %#ok<INUSD>

    cos2Max_flat = cos2Max;
    cos2Max_flat(cos2Max_flat < 0) = 0.0;

    Qpts = qUniform + 4 .* qDayside .* cos2Max_flat;

end


function [xMC, yMC, zMC, rMC] = generateConeMCPoints(nMC, axisVec, coneAngleDeg, rMax)
%GENERATECONEMCPOINTS Generate points uniformly in volume inside a cone.

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

    xLocal = rMC .* sin(theta) .* cos(phi);
    yLocal = rMC .* sin(theta) .* sin(phi);
    zLocal = rMC .* cos(theta);

    localPts = [xLocal, yLocal, zLocal];

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

    R = eye(3) + vx + vx * vx .* ((1 - c) ./ s.^2);

end


function idx = nearestIndex1D(gridCentres, queryVals)
%NEARESTINDEX1D Return nearest index in a monotonic 1-D grid.

    gridCentres = gridCentres(:);
    queryVals = queryVals(:);

    nGrid = numel(gridCentres);

    if nGrid < 2
        idx = ones(size(queryVals));
        return
    end

    fracIdx = interp1( ...
        gridCentres, ...
        1:nGrid, ...
        queryVals, ...
        'nearest', ...
        'extrap' ...
    );

    idx = round(fracIdx);

    idx(idx < 1) = 1;
    idx(idx > nGrid) = nGrid;

end


function [mean_total, mean_cold, frac_cold, nPerBin] = calcMCProfiles( ...
    profileSample, totalSample, coldSample, profile_edges, nBins, profile_lim)
%CALCMCPROFILES Calculate profiles from MC sampled density values.

    profileSample = profileSample(:);
    totalSample = totalSample(:);
    coldSample = coldSample(:);

    valid = isfinite(profileSample) & ...
            isfinite(totalSample) & ...
            isfinite(coldSample) & ...
            totalSample > 0 & ...
            coldSample >= 0 & ...
            profileSample >= 0 & ...
            profileSample <= profile_lim;

    profileSample = profileSample(valid);
    totalSample = totalSample(valid);
    coldSample = coldSample(valid);

    [~, ~, bin_idx] = histcounts(profileSample, profile_edges);

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
    mean_cold  = cold_sum ./ nPerBin;

    frac_cold = mean_cold ./ mean_total;

    badTotal = ~isfinite(mean_total) | mean_total <= 0 | nPerBin <= 0;
    badCold  = ~isfinite(mean_cold)  | mean_cold < 0  | nPerBin <= 0;
    badFrac  = ~isfinite(frac_cold)  | mean_total <= 0 | nPerBin <= 0;

    mean_total(badTotal) = NaN;
    mean_cold(badCold) = NaN;
    frac_cold(badFrac) = NaN;

end