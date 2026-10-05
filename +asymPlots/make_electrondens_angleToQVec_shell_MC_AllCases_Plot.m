function [f, ax, prof] = make_electrondens_angleToQVec_shell_MC_AllCases_Plot( ...
    peCasesNo, swCasesNo, plasDensDir, r0, dr, nMC, nAngleBins)
%MAKE_ELECTRONDENS_ANGLETOQVEC_SHELL_MC_ALLCASES_PLOT
% Plot MC-sampled electron density vs angle to QVec for all simulations.
%
% For each simulation:
%   1. Load total and cold electron density.
%   2. Load QVec.
%   3. Generate MC points uniformly in a radial shell.
%   4. Map each MC point to the nearest grid-cell centre.
%   5. Bin sampled densities by angle to QVec.
%
% Panels:
%   1. Total electron density vs angle to QVec
%   2. Cold electron density vs angle to QVec
%   3. Cold / total electron density vs angle to QVec
%
% Usage:
%
%   [f, ax, prof] = make_electrondens_angleToQVec_shell_MC_AllCases_Plot( ...
%       peCasesNo, swCasesNo, plasDensDir, 50, 10, 2e5, 36);
%
% Inputs:
%   r0          - central shell radius [km]
%   dr          - shell width [km]
%   nMC         - number of MC points per simulation
%   nAngleBins  - number of angular bins between 0 and 180 deg

    if nargin < 4 || isempty(r0)
        r0 = 50;   % km
    end

    if nargin < 5 || isempty(dr)
        dr = 10;   % km
    end

    if nargin < 6 || isempty(nMC)
        nMC = 2e5;
    end

    if nargin < 7 || isempty(nAngleBins)
        nAngleBins = 36;
    end

    % ----------------------------
    % Settings
    % ----------------------------
    cometaryOnly = false;

    caseNameList ={'-X', '+Y', '+Z', '-Z','Ref'};
    
    nCases = numel(peCasesNo);

    inc_gens = [1, 2];
    minEng = 0;
    maxEng = 1;

    rMin = r0 - 0.5 * dr;
    rMax = r0 + 0.5 * dr;

    if rMin < 0
        error('r0 - dr/2 must be >= 0. Current rMin = %.3f km.', rMin);
    end

    % Fixed seed for reproducibility.
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
        'Position', [100, 100, 1450, 500], ...
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

    caseColors = lines(nCases);

    % ----------------------------
    % Output storage
    % ----------------------------
    prof = struct();

    prof.r0 = r0;
    prof.dr = dr;
    prof.rMin = rMin;
    prof.rMax = rMax;
    prof.nMC_requested = nMC;
    prof.nAngleBins = nAngleBins;
    prof.angle_edges = angle_edges;
    prof.angle_centres = angle_centres;

    % ----------------------------
    % Main case loop
    % ----------------------------
    for iCase = 1:nCases

        if iCase <= numel(caseNameList)
            caseStr = caseNameList{iCase};
        else
            caseStr = sprintf('Case %d', iCase);
        end

        fprintf('Processing %s (%d/%d)...\n', caseStr, iCase, nCases);

        peAsymNo = peCasesNo{iCase};
        swAsymNo = swCasesNo{iCase};

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
        % Load QVec for this case
        % ----------------------------
        [~, ~, ~, QVec] = asymPlots.loadOutgassingProperties(peAsymName);

        QVec = QVec(:);

        if numel(QVec) ~= 3
            warning('Skipping %s: QVec must have 3 components.', caseStr);
            continue
        end

        if norm(QVec) == 0
            warning('Skipping %s: QVec has zero magnitude.', caseStr);
            continue
        end

        Qhat = QVec ./ norm(QVec);

        % ----------------------------
        % Build cell-centre grid
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
                ['Size mismatch for %s: rMesh is %s but total density is %s. ', ...
                 'Check whether permute(plas_dens_asym,[2,1,3]) is correct.'], ...
                caseStr, ...
                mat2str(size(rMesh)), ...
                mat2str(size(dens_total)) ...
            );
        end

        if ~isequal(size(rMesh), size(dens_cold))
            error( ...
                ['Size mismatch for %s: rMesh is %s but cold density is %s. ', ...
                 'Check whether permute(plas_dens_cold_asym,[2,1,3]) is correct.'], ...
                caseStr, ...
                mat2str(size(rMesh)), ...
                mat2str(size(dens_cold)) ...
            );
        end

        % ----------------------------
        % MC sample shell volume
        % ----------------------------
        [xMC, yMC, zMC, rMC] = generateShellMCPoints(nMC, rMin, rMax);

        % Keep only MC points inside the actual grid domain.
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
            warning('Skipping %s: no MC points lie inside the grid domain.', caseStr);
            continue
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
        % [mean_total, mean_cold, frac_cold, nPerBin] = calcMCVolumeAngleProfiles( ...
        %     angleToQ_MC, ...
        %     totalSample, ...
        %     coldSample, ...
        %     angle_edges, ...
        %     nAngleBins ...
        % );
        % [mean_total, std_total, mean_cold, std_cold, frac_cold, std_frac, nPerBin] = ...
        %     calcMCVolumeAngleProfiles( ...
        %         angleToQ_MC, ...
        %         totalSample, ...
        %         coldSample, ...
        %         angle_edges, ...
        %         nAngleBins ...
        %     );
        [mean_total, errLow_total, errHigh_total, ...
            mean_cold,  errLow_cold,  errHigh_cold, ...
            frac_cold,  errLow_frac,  errHigh_frac, ...
            nPerBin] = calcMCVolumeAngleProfiles( ...
                angleToQ_MC, ...
                totalSample, ...
                coldSample, ...
                angle_edges, ...
                nAngleBins ...
            );

        % ----------------------------
        % Store profile
        % ----------------------------
        fieldName = matlab.lang.makeValidName(caseStr);

        prof.cases(iCase).caseName = caseStr;
        prof.cases(iCase).peAsymNo = peAsymNo;
        prof.cases(iCase).swAsymNo = swAsymNo;
        prof.cases(iCase).QVec = QVec;
        prof.cases(iCase).Qhat = Qhat;
        prof.cases(iCase).nMC_insideGrid = nValidDomain;
        prof.cases(iCase).mean_total = mean_total;
        % prof.cases(iCase).std_total = std_total;

        prof.cases(iCase).mean_cold = mean_cold;
        % prof.cases(iCase).std_cold = std_cold;

        prof.cases(iCase).frac_cold = frac_cold;
        % prof.cases(iCase).std_frac = std_frac;


        prof.cases(iCase).errLow_total = errLow_total;
        prof.cases(iCase).errHigh_total = errHigh_total;

        prof.cases(iCase).errLow_cold = errLow_cold;
        prof.cases(iCase).errHigh_cold = errHigh_cold;

        prof.cases(iCase).errLow_frac = errLow_frac;
        prof.cases(iCase).errHigh_frac = errHigh_frac;
        prof.cases(iCase).nPerBin = nPerBin;

        prof.byName.(fieldName) = prof.cases(iCase);

        % ----------------------------
        % Plot this case
        % ----------------------------
        errorbar(ax(1), angle_centres, mean_total, errLow_total, errHigh_total, ...
            'o-', ...
            'LineWidth', 1.5, ...
            'MarkerSize', 5, ...
            'CapSize', 4, ...
            'Color', caseColors(iCase, :), ...
            'MarkerFaceColor', caseColors(iCase, :), ...
            'DisplayName', caseStr);

        errorbar(ax(2), angle_centres, mean_cold, errLow_cold, errHigh_cold, ...
            'o-', ...
            'LineWidth', 1.5, ...
            'MarkerSize', 5, ...
            'CapSize', 4, ...
            'Color', caseColors(iCase, :), ...
            'MarkerFaceColor', caseColors(iCase, :), ...
            'DisplayName', caseStr);

        errorbar(ax(3), angle_centres, frac_cold, errLow_frac, errHigh_frac, ...
            'o-', ...
            'LineWidth', 1.5, ...
            'MarkerSize', 5, ...
            'CapSize', 4, ...
            'Color', caseColors(iCase, :), ...
            'MarkerFaceColor', caseColors(iCase, :), ...
            'DisplayName', caseStr);

    end

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

    % ylim(ax(3), [0, 0.2]);
    ylim(ax(3), [1e-3, 0.2]);
    set(ax(3), 'YScale', 'log');


    legend(ax(3), ...
        'Location', 'bestoutside', ...
        'Interpreter', 'none');

    title(t, sprintf( ...
        'MC electron density vs angle to QVec, r = %.1f \\pm %.1f km, N = %d per case', ...
        r0, 0.5 * dr, nMC), ...
        'Interpreter', 'tex');

end


function [xMC, yMC, zMC, rMC] = generateShellMCPoints(nMC, rMin, rMax)
%GENERATESHELLMCPOINTS Generate points uniformly in spherical shell volume.
%
% Uniform in volume requires:
%
%   r = (rMin^3 + U * (rMax^3 - rMin^3))^(1/3)
%
% Directions are sampled uniformly over the sphere.

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


% function [mean_total, mean_cold, frac_cold, nPerBin] = calcMCVolumeAngleProfiles( ...
%     angleVals, totalVals, coldVals, angle_edges, nAngleBins)
% %CALCMCVOLUMEANGLEPROFILES
% % Calculate angular profiles from MC-sampled density values.
% %
% % Because the MC positions are sampled uniformly in shell volume, a simple
% % mean in each angular bin approximates the volume-weighted mean density.

%     angleVals = angleVals(:);
%     totalVals = totalVals(:);
%     coldVals  = coldVals(:);

%     valid = isfinite(angleVals) & ...
%             isfinite(totalVals) & ...
%             isfinite(coldVals) & ...
%             totalVals > 0 & ...
%             coldVals >= 0;

%     angleVals = angleVals(valid);
%     totalVals = totalVals(valid);
%     coldVals  = coldVals(valid);

%     [~, ~, bin_idx] = histcounts(angleVals, angle_edges);

%     inBin = bin_idx > 0;

%     bin_idx = bin_idx(inBin);
%     totalVals = totalVals(inBin);
%     coldVals  = coldVals(inBin);

%     nPerBin = accumarray( ...
%         bin_idx, ...
%         ones(size(bin_idx)), ...
%         [nAngleBins, 1], ...
%         @sum, ...
%         0 ...
%     );

%     total_sum = accumarray( ...
%         bin_idx, ...
%         totalVals, ...
%         [nAngleBins, 1], ...
%         @sum, ...
%         NaN ...
%     );

%     cold_sum = accumarray( ...
%         bin_idx, ...
%         coldVals, ...
%         [nAngleBins, 1], ...
%         @sum, ...
%         NaN ...
%     );

%     mean_total = total_sum ./ nPerBin;
%     mean_cold  = cold_sum ./ nPerBin;

%     frac_cold = mean_cold ./ mean_total;

%     badTotal = ~isfinite(mean_total) | mean_total <= 0 | nPerBin <= 0;
%     badCold  = ~isfinite(mean_cold)  | mean_cold < 0  | nPerBin <= 0;
%     badFrac  = ~isfinite(frac_cold)  | mean_total <= 0 | nPerBin <= 0;

%     mean_total(badTotal) = NaN;
%     mean_cold(badCold) = NaN;
%     frac_cold(badFrac) = NaN;

% end
function [mean_total, errLow_total, errHigh_total, ...
          mean_cold,  errLow_cold,  errHigh_cold, ...
          frac_cold,  errLow_frac,  errHigh_frac, ...
          nPerBin] = calcMCVolumeAngleProfiles( ...
              angleVals, totalVals, coldVals, angle_edges, nAngleBins)
%CALCMCVOLUMEANGLEPROFILES
% Calculate angular profiles from MC-sampled density values.
%
% For total and cold density, uncertainties are calculated in log10-space
% and converted back to asymmetric linear error bars. This is better for
% plotting on log-scaled axes because it avoids non-positive lower bounds.
%
% For cold/total fraction, uncertainties are calculated from the point-wise
% fraction distribution in linear space, with clipping to [0, 1].

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

    fracVals = coldVals ./ totalVals;

    [~, ~, bin_idx] = histcounts(angleVals, angle_edges);

    inBin = bin_idx > 0;

    bin_idx   = bin_idx(inBin);
    totalVals = totalVals(inBin);
    coldVals  = coldVals(inBin);
    fracVals  = fracVals(inBin);

    nPerBin = accumarray( ...
        bin_idx, ...
        ones(size(bin_idx)), ...
        [nAngleBins, 1], ...
        @sum, ...
        0 ...
    );

    % ------------------------------------------------------------
    % Log-space density statistics
    % ------------------------------------------------------------
    [mean_total, errLow_total, errHigh_total] = binnedLogStats( ...
        bin_idx, totalVals, nAngleBins);

    [mean_cold, errLow_cold, errHigh_cold] = binnedLogStats( ...
        bin_idx, coldVals, nAngleBins);

    % ------------------------------------------------------------
    % Fraction statistics
    % ------------------------------------------------------------
    % Keep the central fraction as ratio of binned density means.
    frac_cold = mean_cold ./ mean_total;

    frac_mean_pointwise = accumarray( ...
        bin_idx, ...
        fracVals, ...
        [nAngleBins, 1], ...
        @mean, ...
        NaN ...
    );

    frac_std_pointwise = accumarray( ...
        bin_idx, ...
        fracVals, ...
        [nAngleBins, 1], ...
        @std, ...
        NaN ...
    );

    frac_low  = max(frac_mean_pointwise - frac_std_pointwise, 0);
    frac_high = min(frac_mean_pointwise + frac_std_pointwise, 1);

    errLow_frac  = frac_cold - frac_low;
    errHigh_frac = frac_high - frac_cold;

    % Guard against cases where ratio-of-means and pointwise mean differ
    % enough that the lower/upper ranges cross the plotted central value.
    errLow_frac(errLow_frac < 0 | ~isfinite(errLow_frac)) = NaN;
    errHigh_frac(errHigh_frac < 0 | ~isfinite(errHigh_frac)) = NaN;

    % ------------------------------------------------------------
    % Mask poorly populated or invalid bins
    % ------------------------------------------------------------
    badTotal = ~isfinite(mean_total) | mean_total <= 0 | nPerBin <= 1;
    badCold  = ~isfinite(mean_cold)  | mean_cold <= 0  | nPerBin <= 1;
    badFrac  = ~isfinite(frac_cold)  | frac_cold < 0 | frac_cold > 1 | nPerBin <= 1;

    mean_total(badTotal) = NaN;
    errLow_total(badTotal) = NaN;
    errHigh_total(badTotal) = NaN;

    mean_cold(badCold) = NaN;
    errLow_cold(badCold) = NaN;
    errHigh_cold(badCold) = NaN;

    frac_cold(badFrac) = NaN;
    errLow_frac(badFrac) = NaN;
    errHigh_frac(badFrac) = NaN;

end


function [yMid, errLow, errHigh] = binnedLogStats(bin_idx, vals, nBins)
%BINNEDLOGSTATS
% Return geometric mean and asymmetric 1-sigma bounds in linear space.
%
% yMid   = 10 ^ mean(log10(vals))
% yLow   = 10 ^ (mean(log10(vals)) - std(log10(vals)))
% yHigh  = 10 ^ (mean(log10(vals)) + std(log10(vals)))
%
% This gives multiplicative error bars, which are appropriate for log axes.

    vals = vals(:);
    bin_idx = bin_idx(:);

    valid = isfinite(vals) & vals > 0 & isfinite(bin_idx) & bin_idx > 0;

    vals = vals(valid);
    bin_idx = bin_idx(valid);

    logVals = log10(vals);

    logMean = accumarray( ...
        bin_idx, ...
        logVals, ...
        [nBins, 1], ...
        @mean, ...
        NaN ...
    );

    logStd = accumarray( ...
        bin_idx, ...
        logVals, ...
        [nBins, 1], ...
        @std, ...
        NaN ...
    );

    yMid  = 10.^logMean;
    yLow  = 10.^(logMean - logStd);
    yHigh = 10.^(logMean + logStd);

    errLow  = yMid - yLow;
    errHigh = yHigh - yMid;

    bad = ~isfinite(yMid) | ~isfinite(errLow) | ~isfinite(errHigh) | ...
          yMid <= 0 | errLow < 0 | errHigh < 0;

    yMid(bad) = NaN;
    errLow(bad) = NaN;
    errHigh(bad) = NaN;

end