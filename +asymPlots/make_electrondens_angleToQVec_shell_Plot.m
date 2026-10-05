function [f, ax, prof] = make_electrondens_angleToQVec_shell_Plot( ...
    peCasesNo, swCasesNo, plasDensDir, plotCase, r0, dr)
%MAKE_ELECTRONDENS_ANGLETOQVEC_SHELL_PLOT
% Plot electron density as a function of angle to QVec for one simulation,
% restricted to a radial shell.
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
% Uses volume-weighted angular-bin means:
%
%   <n>_V = sum(n_i dV_i) / sum(dV_i)

    if nargin < 4 || isempty(plotCase)
        plotCase = 1;
    end

    if nargin < 5 || isempty(r0)
        r0 = 50;   % km
    end

    if nargin < 6 || isempty(dr)
        dr = 10;   % km
    end

    % ----------------------------
    % Settings
    % ----------------------------
    cometaryOnly = false;

    nAngleBins = 36;     % 5 degree bins
    caseNameList = {'-X', '+Y', '+Z', '-Z','Ref'};

    inc_gens = [1, 2];
    minEng = 0;
    maxEng = 1;

    rMin = r0 - 0.5 * dr;
    rMax = r0 + 0.5 * dr;

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
    % Build spatial mesh using cell centres
    % ----------------------------
    x_edges = gridsAsym{1};
    y_edges = gridsAsym{2};
    z_edges = gridsAsym{3};

    dx = diff(x_edges);
    dy = diff(y_edges);
    dz = diff(z_edges);

    x = x_edges(1:end-1) + 0.5 .* dx;
    y = y_edges(1:end-1) + 0.5 .* dy;
    z = z_edges(1:end-1) + 0.5 .* dz;

    [X_mesh, Y_mesh, Z_mesh] = meshgrid(x, y, z);
    rMesh = sqrt(X_mesh.^2 + Y_mesh.^2 + Z_mesh.^2);

    [dX, dY, dZ] = meshgrid(dx, dy, dz);
    dV = dX .* dY .* dZ;

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
    % Calculate angle to QVec
    % ----------------------------
    rSafe = rMesh;
    rSafe(rSafe == 0) = NaN;

    Xhat = X_mesh ./ rSafe;
    Yhat = Y_mesh ./ rSafe;
    Zhat = Z_mesh ./ rSafe;

    cosAngle = ...
        Xhat .* Qhat(1) + ...
        Yhat .* Qhat(2) + ...
        Zhat .* Qhat(3);

    cosAngle(cosAngle > 1) = 1;
    cosAngle(cosAngle < -1) = -1;

    angleToQ = acosd(cosAngle);

    % ----------------------------
    % Radial shell mask
    % ----------------------------
    shellMask = rMesh >= rMin & rMesh <= rMax;

    nShell = nnz(shellMask);

    if nShell == 0
        warning('No grid cells found in shell %.2f <= r <= %.2f km.', rMin, rMax);
    end

    % ----------------------------
    % Angular profiles within shell
    % ----------------------------
    [mean_total, mean_cold, frac_cold, nPerBin] = calcVolumeWeightedAngleProfiles( ...
        angleToQ, ...
        dens_total, ...
        dens_cold, ...
        dV, ...
        shellMask, ...
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

    legend(ax(3), ...
        'Location', 'best', ...
        'Interpreter', 'none');

    if plotCase <= numel(caseNameList)
        caseStr = caseNameList{plotCase};
    else
        caseStr = sprintf('Case %d', plotCase);
    end

    title(t, sprintf( ...
        '%s: Electron density vs angle to QVec, r = %.1f \\pm %.1f km', ...
        caseStr, r0, 0.5 * dr), ...
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
    prof.nShell = nShell;

    prof.angle_edges = angle_edges;
    prof.angle_centres = angle_centres;

    prof.mean_total = mean_total;
    prof.mean_cold = mean_cold;
    prof.frac_cold = frac_cold;
    prof.nPerBin = nPerBin;

end


function [mean_total, mean_cold, frac_cold, nPerBin] = calcVolumeWeightedAngleProfiles( ...
    angleToQ, dens_total, dens_cold, dV, mask, angle_edges, nAngleBins)
%CALCVOLUMEWEIGHTEDANGLEPROFILES
% Calculate volume-weighted electron-density profiles as a function of
% angle to QVec.

    angleVals = angleToQ(mask);
    totalVals = dens_total(mask);
    coldVals  = dens_cold(mask);
    volVals   = dV(mask);

    angleVals = angleVals(:);
    totalVals = totalVals(:);
    coldVals  = coldVals(:);
    volVals   = volVals(:);

    valid = isfinite(angleVals) & ...
            isfinite(totalVals) & ...
            isfinite(coldVals) & ...
            isfinite(volVals) & ...
            totalVals > 0 & ...
            coldVals >= 0 & ...
            volVals > 0;

    angleVals = angleVals(valid);
    totalVals = totalVals(valid);
    coldVals  = coldVals(valid);
    volVals   = volVals(valid);

    [~, ~, bin_idx] = histcounts(angleVals, angle_edges);

    inBin = bin_idx > 0;

    bin_idx = bin_idx(inBin);
    totalVals = totalVals(inBin);
    coldVals  = coldVals(inBin);
    volVals   = volVals(inBin);

    nPerBin = accumarray( ...
        bin_idx, ...
        ones(size(bin_idx)), ...
        [nAngleBins, 1], ...
        @sum, ...
        0 ...
    );

    vol_sum = accumarray( ...
        bin_idx, ...
        volVals, ...
        [nAngleBins, 1], ...
        @sum, ...
        NaN ...
    );

    total_weighted_sum = accumarray( ...
        bin_idx, ...
        totalVals .* volVals, ...
        [nAngleBins, 1], ...
        @sum, ...
        NaN ...
    );

    cold_weighted_sum = accumarray( ...
        bin_idx, ...
        coldVals .* volVals, ...
        [nAngleBins, 1], ...
        @sum, ...
        NaN ...
    );

    mean_total = total_weighted_sum ./ vol_sum;
    mean_cold  = cold_weighted_sum ./ vol_sum;

    frac_cold = mean_cold ./ mean_total;

    badTotal = ~isfinite(mean_total) | mean_total <= 0 | nPerBin <= 0;
    badCold  = ~isfinite(mean_cold)  | mean_cold < 0  | nPerBin <= 0;
    badFrac  = ~isfinite(frac_cold)  | mean_total <= 0 | nPerBin <= 0;

    mean_total(badTotal) = NaN;
    mean_cold(badCold) = NaN;
    frac_cold(badFrac) = NaN;

end