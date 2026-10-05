function [f, ax, prof] = make_electrondens_angleToQVec_shell_MC_Plot( ...
    peCasesNo, swCasesNo, plasDensDir, plotCase, r0, dr, nMC, nAngleBins)
%MAKE_ELECTRONDENS_ANGLETOQVEC_SHELL_MC_PLOT
% Plot electron density as a function of angle to QVec for one simulation,
% restricted to a radial shell, using Monte Carlo sampling of space.
%
% Panels:
%   1. Total electron density vs angle to QVec
%   2. Cold electron density vs angle to QVec
%   3. Cold / total electron density vs angle to QVec
%
% The radial shell is:
%
%   r0 - dr/2 <= r <= r0 + dr/2
%
% The MC points are sampled uniformly in shell volume. Each MC point is then
% assigned the total/cold density of the nearest grid-cell centre.
%
% Usage:
%
%   [f, ax, prof] = make_electrondens_angleToQVec_shell_MC_Plot( ...
%       peCasesNo, swCasesNo, plasDensDir, 1, 50, 10, 2e5, 36);
%
% Inputs:
%   plotCase    - simulation case index
%   r0          - central shell radius [km]
%   dr          - shell width [km]
%   nMC         - number of MC points in shell
%   nAngleBins  - number of angular bins between 0 and 180 deg

    if nargin < 4 || isempty(plotCase)
        plotCase = 1;
    end

    if nargin < 5 || isempty(r0)
        r0 = 50;   % km
    end

    if nargin < 6 || isempty(dr)
        dr = 10;   % km
    end

    if nargin < 7 || isempty(nMC)
        nMC = 2e5;
    end

    if nargin < 8 || isempty(nAngleBins)
        nAngleBins = 36;
    end

    % ----------------------------
    % Settings
    % ----------------------------
    cometaryOnly = false;

    caseNameList = {'-X', '+Y', '+Z', '-Z','Ref'};%{'Ref', '-X', '+Y', '+Z', '-Z'};

    inc_gens = [1, 2];
    minEng = 0;
    maxEng = 1;

    rMin = r0 - 0.5 * dr;
    rMax = r0 + 0.5 * dr;

    if rMin < 0
        error('r0 - dr/2 must be >= 0. Current rMin = %.3f km.', rMin);
    end

    % Optional reproducibility.
    rng(1);

    % ----------------------------
    % Angular bins
    % ----------------------------
    angle_edges = linspace(0, 180, nAngleBins + 1);
    angle_centres = angle_edges(1:end-1) + 0.5 .* diff(angle_edges);

    % ----------------------------
    % Create figure
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
    % Load QVec
    % ----------------------------
    [~, ~, ~, QVec] = asymPlots.loadOutgassingProperties(peAsymName);

    QVec = QVec(:);

    if numel(QVec) ~= 3
        error('QVec must have 3 components, but has size %s.', mat2str(size(QVec)));
    end

    if norm(QVec) == 0
        error('QVec has zero magnitude.');
    end

    Qhat = QVec ./ norm(QVec);

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

    % Match density orientation to meshgrid output.
    % meshgrid gives size [ny, nx, nz].
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
    % Monte Carlo sample shell volume
    % ----------------------------
    [xMC, yMC, zMC, rMC] = generateShellMCPoints(nMC, rMin, rMax);

    % Keep only MC points inside the actual grid domain.
    % This avoids clamping out-of-domain points onto boundary cells.
    insideGrid = ...
        xMC >= x_edges(1) & xMC <= x_edges(end) & ...
        yMC >= y_edges(1) & yMC <= y_edges(end) & ...
        zMC >= z_edges(1) & zMC <= z_edges(end);

    xMC = xMC(insideGrid);
    yMC = yMC(insideGrid);
    zMC = zMC(insideGrid);
    rMC = rMC(insideGrid);

    nValidDomain = numel(rMC);

    if nValidDomain == 0
        error('No MC points lie inside the grid domain for this shell.');
    end

    % ----------------------------
    % Angle to QVec for MC points
    % ----------------------------
    rSafe = rMC;
    rSafe(rSafe == 0) = NaN;

    xHatMC = xMC ./ rSafe;
    yHatMC = yMC ./ rSafe;
    zHatMC = zMC ./ rSafe;

    cosAngle = ...
        xHatMC .* Qhat(1) + ...
        yHatMC .* Qhat(2) + ...
        zHatMC .* Qhat(3);

    cosAngle(cosAngle > 1) = 1;
    cosAngle(cosAngle < -1) = -1;

    angleToQ_MC = acosd(cosAngle);

    % ----------------------------
    % Nearest grid-cell lookup
    % ----------------------------
    ix = nearestIndex1D(x_centres, xMC);
    iy = nearestIndex1D(y_centres, yMC);
    iz = nearestIndex1D(z_centres, zMC);

    linIdx = sub2ind(size(dens_total), iy, ix, iz);

    totalSample = dens_total(linIdx);
    coldSample  = dens_cold(linIdx);

    % ----------------------------
    % Angular profiles from MC samples
    % ----------------------------
    [mean_total, mean_cold, frac_cold, nPerBin] = calcMCVolumeAngleProfiles( ...
        angleToQ_MC, ...
        totalSample, ...
        coldSample, ...
        angle_edges, ...
        nAngleBins ...
    );

    % ----------------------------
    % Plot
    % ----------------------------
    plot(ax(1), angle_centres, mean_total, ...
        'LineWidth', 2.2, ...
        'Color', 'k', ...
        'DisplayName', 'Total');

    plot(ax(2), angle_centres, mean_cold, ...
        'LineWidth', 2.2, ...
        'Color', 'k', ...
        'DisplayName', 'Cold');

    plot(ax(3), angle_centres, frac_cold, ...
        'LineWidth', 2.2, ...
        'Color', 'k', ...
        'DisplayName', 'Cold / total');

    % Optional: show bin counts on a right axis for diagnostics.
    % Uncomment if useful.
    %
    % yyaxis(ax(3), 'right');
    % plot(ax(3), angle_centres, nPerBin, 'Color', [0.5, 0.5, 0.5], ...
    %     'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'N per bin');
    % ylabel(ax(3), 'MC samples per bin');

    % ----------------------------
    % Axis formatting
    % ----------------------------
    xlabel(ax(1), 'Angle to QVec [deg]', 'Interpreter', 'tex');
    xlabel(ax(2), 'Angle to QVec [deg]', 'Interpreter', 'tex');
    xlabel(ax(3), 'Angle to QVec [deg]', 'Interpreter', 'tex');

    ylabel(ax(1), 'Total n_e [cm^{-3}]', 'Interpreter', 'tex');
    ylabel(ax(2), 'Cold n_e [cm^{-3}]', 'Interpreter', 'tex');
    ylabel(ax(3), 'Cold / total', 'Interpreter', 'tex');

    title(ax(1), 'Total electron density');
    title(ax(2), 'Cold electron density');
    title(ax(3), 'Cold electron fraction');

    for iAx = 1:3
        xlim(ax(iAx), [0, 180]);
        xticks(ax(iAx), 0:30:180);
    end

    ylim(ax(3), [0, 1]);
    % set(ax(3), 'YScale', 'log');

    legend(ax(3), ...
        'Location', 'best', ...
        'Interpreter', 'none');

    if plotCase <= numel(caseNameList)
        caseStr = caseNameList{plotCase};
    else
        caseStr = sprintf('Case %d', plotCase);
    end

    title(t, sprintf( ...
        '%s: MC electron density vs angle to QVec, r = %.1f \\pm %.1f km, N = %d', ...
        caseStr, r0, 0.5 * dr, nValidDomain), ...
        'Interpreter', 'tex');

    % ----------------------------
    % Store outputs
    % ----------------------------
    prof = struct();

    prof.plotCase = plotCase;
    prof.QVec = QVec;
    prof.Qhat = Qhat;

    prof.r0 = r0;
    prof.dr = dr;
    prof.rMin = rMin;
    prof.rMax = rMax;

    prof.nMC_requested = nMC;
    prof.nMC_insideGrid = nValidDomain;

    prof.angle_edges = angle_edges;
    prof.angle_centres = angle_centres;

    prof.mean_total = mean_total;
    prof.mean_cold = mean_cold;
    prof.frac_cold = frac_cold;
    prof.nPerBin = nPerBin;

    prof.angleToQ_MC = angleToQ_MC;
    prof.totalSample = totalSample;
    prof.coldSample = coldSample;
    prof.rMC = rMC;

end


function [xMC, yMC, zMC, rMC] = generateShellMCPoints(nMC, rMin, rMax)
%GENERATESHELLMCPOINTS Generate points uniformly in spherical shell volume.
%
% Uniform in volume requires:
%
%   r = (rMin^3 + U * (rMax^3 - rMin^3))^(1/3)
%
% Directions are sampled uniformly over the sphere by taking cos(theta)
% uniformly distributed between -1 and 1.

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


function idx = nearestIndex1D(gridCentres, queryVals)
%NEARESTINDEX1D Return nearest index in a monotonic 1-D grid.
%
% This avoids knnsearch and is fast for structured grids.

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
        'extrap' ...
    );

    idx = round(idx);

    idx(idx < 1) = 1;
    idx(idx > nGrid) = nGrid;

end


function [mean_total, mean_cold, frac_cold, nPerBin] = calcMCVolumeAngleProfiles( ...
    angleVals, totalVals, coldVals, angle_edges, nAngleBins)
%CALCMCVOLUMEANGLEPROFILES
% Calculate angular profiles from MC-sampled density values.
%
% Because the MC positions are sampled uniformly in shell volume, a simple
% mean in each angular bin approximates the volume-weighted mean density.

    angleVals = angleVals(:);
    totalVals = totalVals(:);
    coldVals  = coldVals(:);

    valid = isfinite(angleVals) & ...
            isfinite(totalVals) & ...
            isfinite(coldVals) & ...
            totalVals > 0 & ...
            coldVals >= 0;

    angleVals = angleVals(valid);
    totalVals = totalVals(valid);
    coldVals  = coldVals(valid);

    [~, ~, bin_idx] = histcounts(angleVals, angle_edges);

    inBin = bin_idx > 0;

    bin_idx = bin_idx(inBin);
    totalVals = totalVals(inBin);
    coldVals  = coldVals(inBin);

    nPerBin = accumarray( ...
        bin_idx, ...
        ones(size(bin_idx)), ...
        [nAngleBins, 1], ...
        @sum, ...
        0 ...
    );

    total_sum = accumarray( ...
        bin_idx, ...
        totalVals, ...
        [nAngleBins, 1], ...
        @sum, ...
        NaN ...
    );

    cold_sum = accumarray( ...
        bin_idx, ...
        coldVals, ...
        [nAngleBins, 1], ...
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