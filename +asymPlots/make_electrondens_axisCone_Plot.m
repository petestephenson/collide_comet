function [f, ax] = make_electrondens_axisCone_Plot(peCasesNo, swCasesNo, plasDensDir, plotCase)
%MAKE_ELECTRONDENS_AXISCONE_PLOT
% Plot volume-weighted mean total and cold electron density profiles
% along +/-X, +/-Y, and +/-Z semi-axes for one simulation case.
%
% Points are selected if they lie within coneAngleDeg of each semi-axis.

    if nargin < 4
        plotCase = 1;
    end

    % ----------------------------
    % Plot settings
    % ----------------------------
    cometaryOnly = false;
    plot_lim = 300;          % km
    nBins = 30;
    coneAngleDeg = 30;       % half-angle around each semi-axis

    caseNameList = {'-X', '+Y', '+Z', '-Z','Ref'};

    inc_gens = [1, 2];
    minEng = 0;
    maxEng = 1;

    % ----------------------------
    % Axis directions to sample
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

    axisColors = lines(3);

    lineColors = [
        axisColors(1, :);
        axisColors(1, :);
        axisColors(2, :);
        axisColors(2, :);
        axisColors(3, :);
        axisColors(3, :)
    ];

    lineStyles = {'-', '--', '-', '--', '-', '--'};

    cosCone = cosd(coneAngleDeg);

    % ----------------------------
    % Create figure/subplots
    % ----------------------------
    f = figure( ...
        'Position', [100, 100, 1100, 500], ...
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
            'TickDir', 'out', ...
            'YScale', 'log');
    end

    % ----------------------------
    % Radial bins
    % ----------------------------
    r_edges = linspace(0, plot_lim, nBins + 1);
    r_centres = r_edges(1:end-1) + 0.5 * diff(r_edges);

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
    % Build spatial mesh using cell centres
    % ----------------------------
    x_edges = gridsAsym{1};
    y_edges = gridsAsym{2};
    z_edges = gridsAsym{3};

    dx = diff(x_edges);
    dy = diff(y_edges);
    dz = diff(z_edges);

    x = x_edges(1:end-1) + 0.5 * dx;
    y = y_edges(1:end-1) + 0.5 * dy;
    z = z_edges(1:end-1) + 0.5 * dz;

    [X_mesh, Y_mesh, Z_mesh] = meshgrid(x, y, z);
    rMesh = sqrt(X_mesh.^2 + Y_mesh.^2 + Z_mesh.^2);

    [dX, dY, dZ] = meshgrid(dx, dy, dz);
    dV = dX .* dY .* dZ;

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
    % Unit position vector at each grid cell
    % ----------------------------
    rSafe = rMesh;
    rSafe(rSafe == 0) = NaN;

    Xhat = X_mesh ./ rSafe;
    Yhat = Y_mesh ./ rSafe;
    Zhat = Z_mesh ./ rSafe;

    % ----------------------------
    % Plot profiles along each semi-axis
    % ----------------------------
    for iDir = 1:numel(axisLabels)

        thisVec = axisVecs(iDir, :);

        cosTheta = ...
            thisVec(1) .* Xhat + ...
            thisVec(2) .* Yhat + ...
            thisVec(3) .* Zhat;

        mask = cosTheta >= cosCone;

        % [mean_total, mean_cold] = calcVolumeWeightedProfiles( ...
        %     rMesh, ...
        %     dens_total, ...
        %     dens_cold, ...
        %     dV, ...
        %     mask, ...
        %     r_edges, ...
        %     nBins, ...
        %     plot_lim ...
        % );

        [mean_total, mean_cold, frac_cold] = calcVolumeWeightedProfiles( ...
            rMesh, ...
            dens_total, ...
            dens_cold, ...
            dV, ...
            mask, ...
            r_edges, ...
            nBins, ...
            plot_lim ...
        );
        plot(ax(1), r_centres, mean_total, ...
            'LineWidth', 2.2, ...
            'Color', lineColors(iDir, :), ...
            'LineStyle', lineStyles{iDir}, ...
            'DisplayName', axisLabels{iDir});

        plot(ax(2), r_centres, mean_cold, ...
            'LineWidth', 2.2, ...
            'Color', lineColors(iDir, :), ...
            'LineStyle', lineStyles{iDir}, ...
            'DisplayName', axisLabels{iDir});

        plot(ax(3), r_centres, frac_cold, ...
            'LineWidth', 2.2, ...
            'Color', lineColors(iDir, :), ...
            'LineStyle', lineStyles{iDir}, ...
            'DisplayName', axisLabels{iDir});
    end

    % ----------------------------
    % Axis formatting
    % ----------------------------
    xlabel(ax(1), 'r [km]', 'Interpreter', 'tex');
    xlabel(ax(2), 'r [km]', 'Interpreter', 'tex');
    xlabel(ax(3), 'r [km]', 'Interpreter', 'tex');

    ylabel(ax(1), 'Total n_e [cm^{-3}]', 'Interpreter', 'tex');
    ylabel(ax(2), 'Cold n_e [cm^{-3}]', 'Interpreter', 'tex');
    ylabel(ax(3), 'Cold / total', 'Interpreter', 'tex');

    title(ax(1), sprintf('Total electron density, \\pm axes, %d^\\circ cone', coneAngleDeg));
    title(ax(2), sprintf('Cold electron density, \\pm axes, %d^\\circ cone', coneAngleDeg));
    title(ax(3), sprintf('Cold electron fraction, \\pm axes, %d^\\circ cone', coneAngleDeg));

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

    title(t, sprintf('%s: Electron density along semi-axes', caseStr), ...
        'Interpreter', 'none');

end

function [mean_total, mean_cold, frac_cold] = calcVolumeWeightedProfiles( ...
    rMesh, dens_total, dens_cold, dV, mask, r_edges, nBins, plot_lim)
%CALCVOLUMEWEIGHTEDPROFILES Calculate volume-weighted mean density profiles.
%
% For each radial bin:
%
%   <n>_V = sum(n_i dV_i) / sum(dV_i)
%
% The cold electron fraction is:
%
%   <n_cold>_V / <n_total>_V

    rVals = rMesh(mask);
    totalVals = dens_total(mask);
    coldVals = dens_cold(mask);
    volVals = dV(mask);

    rVals = rVals(:);
    totalVals = totalVals(:);
    coldVals = coldVals(:);
    volVals = volVals(:);

    valid = isfinite(rVals) & ...
            isfinite(totalVals) & ...
            isfinite(coldVals) & ...
            isfinite(volVals) & ...
            totalVals > 0 & ...
            coldVals >= 0 & ...
            volVals > 0 & ...
            rVals <= plot_lim;

    rVals = rVals(valid);
    totalVals = totalVals(valid);
    coldVals = coldVals(valid);
    volVals = volVals(valid);

    [~, ~, bin_idx] = histcounts(rVals, r_edges);

    inBin = bin_idx > 0;

    bin_idx = bin_idx(inBin);
    totalVals = totalVals(inBin);
    coldVals = coldVals(inBin);
    volVals = volVals(inBin);

    vol_sum = accumarray( ...
        bin_idx, ...
        volVals, ...
        [nBins, 1], ...
        @sum, ...
        NaN ...
    );

    total_weighted_sum = accumarray( ...
        bin_idx, ...
        totalVals .* volVals, ...
        [nBins, 1], ...
        @sum, ...
        NaN ...
    );

    cold_weighted_sum = accumarray( ...
        bin_idx, ...
        coldVals .* volVals, ...
        [nBins, 1], ...
        @sum, ...
        NaN ...
    );

    mean_total = total_weighted_sum ./ vol_sum;
    mean_cold  = cold_weighted_sum  ./ vol_sum;

    frac_cold = mean_cold ./ mean_total;

    mean_total(~isfinite(mean_total) | mean_total <= 0) = NaN;
    mean_cold(~isfinite(mean_cold) | mean_cold < 0) = NaN;
    frac_cold(~isfinite(frac_cold) | mean_total <= 0) = NaN;

end