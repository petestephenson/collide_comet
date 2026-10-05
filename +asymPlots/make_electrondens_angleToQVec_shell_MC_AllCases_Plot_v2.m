function [f, ax, prof] = make_electrondens_angleToQVec_shell_MC_AllCases_Plot_v2( ...
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
        % grid(ax(iAx), 'on');

        set(ax(iAx), ...
            'FontSize', 12, ...
            'LineWidth', 1.1, ...
            'TickDir', 'out');
    end

    set(ax(1), 'YScale', 'log');
    set(ax(2), 'YScale', 'log');
    set(ax(3), 'YScale', 'log');

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
    % prof.profiles = struct([]);
    prof.profiles = {};

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
        % Unit position vectors for MC points
        % ----------------------------
        rSafe = rMC;
        rSafe(rSafe == 0) = NaN;

        xHatMC = xMC ./ rSafe;
        yHatMC = yMC ./ rSafe;
        zHatMC = zMC ./ rSafe;

        % ----------------------------
        % Define angle axes
        %
        % Reference case:
        %   plot relative to +X, +Y, and +Z.
        %
        % Asymmetric cases:
        %   plot relative to QVec only.
        % ----------------------------
        isRefCase = strcmpi(strtrim(string(caseStr)), "ref") || ...
                    contains(lower(string(caseStr)), "ref");

        % if isRefCase

        %     angleDirs = [
        %         1, 0, 0
        %         0, 1, 0
        %         0, 0, 1
        %     ];

        %     angleLabels = {
        %         'Ref: +X'
        %         'Ref: +Y'
        %         'Ref: +Z'
        %     };

        %     angleLineStyles = {
        %         '-'
        %         '--'
        %         ':'
        %     };

        % else

        %     angleDirs = Qhat(:).';
        %     angleLabels = {caseStr};
        %     angleLineStyles = {'-'};

        % end
        if isRefCase

            % ------------------------------------------------------------
            % For the reference case, plot profiles relative to each
            % asymmetric case axis/direction, e.g. -X, +Y, +Z, -Z.
            % ------------------------------------------------------------
            asymCaseIdx = find(~strcmpi(string(caseNameList), "Ref"));

            nAsymDirs = numel(asymCaseIdx);

            angleDirs = NaN(nAsymDirs, 3);
            angleLabels = cell(nAsymDirs, 1);
            angleLineStyles = cell(nAsymDirs, 1);
            angleMarkers = cell(nAsymDirs, 1);
            angleColors = NaN(nAsymDirs, 3);

            refLineStyles = {'--', '--', '--', '--', ':', '-.'};

            for iDir = 1:nAsymDirs

                iAsymCase = asymCaseIdx(iDir);
                asymCaseStr = caseNameList{iAsymCase};

                angleDirs(iDir, :) = caseAxisDirectionFromName(asymCaseStr);

                angleLabels{iDir} = sprintf('Ref relative to %s', asymCaseStr);

                angleLineStyles{iDir} = refLineStyles{ ...
                    mod(iDir - 1, numel(refLineStyles)) + 1};

                % Match marker style to the axis.
                if contains(asymCaseStr, 'X')
                    angleMarkers{iDir} = 'o';
                elseif contains(asymCaseStr, 'Y')
                    angleMarkers{iDir} = 'pentagram';
                elseif contains(asymCaseStr, 'Z')
                    angleMarkers{iDir} = '^';
                else
                    angleMarkers{iDir} = 's';
                end

                % Use the same colour as the corresponding asymmetric case,
                % but with reduced alpha later.
                angleColors(iDir, :) = caseColors(iAsymCase, :);

            end

        else

            angleDirs = Qhat(:).';
            angleLabels = {caseStr};
            angleLineStyles = {'-'};
            angleMarkers = {'o'};
            angleColors = caseColors(iCase, :);

        end
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
        for iAngleAxis = 1:size(angleDirs, 1)

            angleDir = angleDirs(iAngleAxis, :).';
            angleDir = angleDir ./ norm(angleDir);

            cosAngle = ...
                xHatMC .* angleDir(1) + ...
                yHatMC .* angleDir(2) + ...
                zHatMC .* angleDir(3);

            cosAngle(cosAngle > 1) = 1;
            cosAngle(cosAngle < -1) = -1;

            angleMC = acosd(cosAngle);

            [mean_total, errLow_total, errHigh_total, ...
                mean_cold,  errLow_cold,  errHigh_cold, ...
                frac_cold,  errLow_frac,  errHigh_frac, ...
                nPerBin] = calcMCVolumeAngleProfiles( ...
                    angleMC, ...
                    totalSample, ...
                    coldSample, ...
                    angle_edges, ...
                    nAngleBins ...
                );

            % ----------------------------
            % Store profile
            % ----------------------------
            if isRefCase
                fieldName = matlab.lang.makeValidName(angleLabels{iAngleAxis});
                profLabel = angleLabels{iAngleAxis};
            else
                fieldName = matlab.lang.makeValidName(caseStr);
                profLabel = caseStr;
            end

            profEntry = struct();

            profEntry.caseName = caseStr;
            profEntry.profileLabel = profLabel;
            profEntry.angleDir = angleDir;
            profEntry.peAsymNo = peAsymNo;
            profEntry.swAsymNo = swAsymNo;
            profEntry.QVec = QVec;
            profEntry.Qhat = Qhat;
            profEntry.nMC_insideGrid = nValidDomain;

            profEntry.mean_total = mean_total;
            profEntry.mean_cold = mean_cold;
            profEntry.frac_cold = frac_cold;

            profEntry.errLow_total = errLow_total;
            profEntry.errHigh_total = errHigh_total;

            profEntry.errLow_cold = errLow_cold;
            profEntry.errHigh_cold = errHigh_cold;

            profEntry.errLow_frac = errLow_frac;
            profEntry.errHigh_frac = errHigh_frac;

            profEntry.nPerBin = nPerBin;

            prof.profiles{end + 1} = profEntry; %#ok<AGROW>
            prof.byName.(fieldName) = profEntry;

            % ----------------------------
            % Plot this profile
            % ----------------------------
            
            % if isRefCase
            %     plotColor = [0,0,0];
            %     plotStyle = angleLineStyles{iAngleAxis};
            %     displayName = angleLabels{iAngleAxis};
            %     lalpha = 0.7;
            %     switch iAngleAxis
            %         case 1
            %             plotMarker = 'o';  % +X
            %         case 2
            %             plotMarker = "pentagram";  % +Y
            %         case 3
            %             plotMarker = '^';  % +Z
            %     end

            % else
            %     plotColor = caseColors(iCase, :);
            %     plotStyle = '-';
            %     displayName = caseStr;
            %     plotMarker = 'o';
            %     lalpha = 1.0;
            %     % add a switch depending on the caseStr to assign different markers if desired
            %     if contains(caseStr, 'X')    
            %         plotMarker = 'o';
            %     elseif contains(caseStr, 'Y')     %     
            %         plotMarker = 'pentagram';
            %     elseif contains(caseStr, 'Z')     %     plotMarker = '^';
            %         plotMarker = '^';
            %     end

            % end
            if isRefCase

                plotColor = angleColors(iAngleAxis, :);
                plotStyle = angleLineStyles{iAngleAxis};
                plotMarker = angleMarkers{iAngleAxis};
                displayName = angleLabels{iAngleAxis};

                % Keep reference curves visually secondary.
                lalpha = 0.55;
                patchAlpha = 0.06;

            else

                plotColor = caseColors(iCase, :);
                plotStyle = '-';
                displayName = caseStr;
                lalpha = 1.0;
                patchAlpha = 0.13;

                if contains(caseStr, 'X')
                    plotMarker = 'o';
                elseif contains(caseStr, 'Y')
                    plotMarker = 'pentagram';
                elseif contains(caseStr, 'Z')
                    plotMarker = '^';
                else
                    plotMarker = 's';
                end

            end
            plotColorAlpha = [plotColor, lalpha];
            patchAlpha = 0.1;
            % errorbar(ax(1), angle_centres, mean_total, errLow_total, errHigh_total, ...
            %     'LineStyle', plotStyle, ...
            %     'Marker', plotMarker, ...
            %     'LineWidth', 1.5, ...
            %     'MarkerSize', 5, ...
            %     'CapSize', 4, ...
            %     'Color', plotColorAlpha, ...
            %     'MarkerFaceColor', plotColor, ...
            %     'DisplayName', displayName);

            % errorbar(ax(2), angle_centres, mean_cold, errLow_cold, errHigh_cold, ...
            %     'LineStyle', plotStyle, ...
            %     'Marker', plotMarker, ...
            %     'LineWidth', 1.5, ...
            %     'MarkerSize', 5, ...
            %     'CapSize', 4, ...
            %     'Color', plotColorAlpha, ...
            %     'MarkerFaceColor', plotColor, ...
            %     'DisplayName', displayName);

            % errorbar(ax(3), angle_centres, frac_cold, errLow_frac, errHigh_frac, ...
            %     'LineStyle', plotStyle, ...
            %     'Marker', plotMarker, ...
            %     'LineWidth', 1.5, ...
            %     'MarkerSize', 5, ...
            %     'CapSize', 4, ...
            %     'Color', plotColorAlpha, ...
            %     'MarkerFaceColor', plotColor, ...
            %     'DisplayName', displayName);
            % if isRefCase
            %     patchAlpha = 0.07;
            % else
            %     patchAlpha = 0.13;
            % end

            plotProfileWithPatch( ...
                ax(1), ...
                angle_centres, ...
                mean_total, ...
                errLow_total, ...
                errHigh_total, ...
                plotColor, ...
                plotStyle, ...
                plotMarker, ...
                displayName, ...
                patchAlpha, ...
                lalpha);

            plotProfileWithPatch( ...
                ax(2), ...
                angle_centres, ...
                mean_cold, ...
                errLow_cold, ...
                errHigh_cold, ...
                plotColor, ...
                plotStyle, ...
                plotMarker, ...
                displayName, ...
                patchAlpha, ...
                lalpha);

            plotProfileWithPatch( ...
                ax(3), ...
                angle_centres, ...
                frac_cold, ...
                errLow_frac, ...
                errHigh_frac, ...
                plotColor, ...
                plotStyle, ...
                plotMarker, ...
                displayName, ...
                patchAlpha, ...
                lalpha);
        end
    end

    % ----------------------------
    % Axis formatting
    % ----------------------------
    % xlabel(ax(1), 'Angle to QVec [deg]', 'Interpreter', 'tex');
    % xlabel(ax(2), 'Angle to QVec [deg]', 'Interpreter', 'tex');
    % xlabel(ax(3), 'Angle to QVec [deg]', 'Interpreter', 'tex');
    xlabel(ax(1), 'Angle to reference direction [deg]', 'Interpreter', 'tex');
    xlabel(ax(2), 'Angle to reference direction [deg]', 'Interpreter', 'tex');
    xlabel(ax(3), 'Angle to reference direction [deg]', 'Interpreter', 'tex');

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

    % title(t, sprintf( ...
    %     'MC electron density vs angle to QVec, r = %.1f \\pm %.1f km, N = %d per case', ...
    %     r0, 0.5 * dr, nMC), ...
    %     'Interpreter', 'tex');

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

function h = plotProfileWithPatch(ax, x, y, errLow, errHigh, plotColor, ...
    plotStyle, plotMarker, displayName, patchAlpha, lineAlpha, yMinPatch)
%PLOTPROFILEWITHPATCH Plot central profile with semi-transparent uncertainty patch.

    if nargin < 12 || isempty(yMinPatch)
        yMinPatch = 1e-30;
    end

    x = x(:);
    y = y(:);
    errLow = errLow(:);
    errHigh = errHigh(:);

    yLow = y - errLow;
    yHigh = y + errHigh;

    good = isfinite(x) & isfinite(y) & ...
           isfinite(yLow) & isfinite(yHigh);

    
    % For log axes, keep valid central values, but clip the lower patch edge.
    if strcmpi(ax.YScale, 'log')
        good = good & y > 0 & yHigh > 0;

        yLow(yLow <= yMinPatch | ~isfinite(yLow)) = yMinPatch;
        yLow = max(yLow, yMinPatch);
    end
    x = x(good);
    y = y(good);
    yLow = yLow(good);
    yHigh = yHigh(good);

    if isempty(x)
        h = gobjects(1, 1);
        return
    end

    % Sort by x so the patch is well behaved.
    [x, sortIdx] = sort(x);
    y = y(sortIdx);
    yLow = yLow(sortIdx);
    yHigh = yHigh(sortIdx);

    % Semi-transparent uncertainty envelope.
    patch( ...
        ax, ...
        [x; flipud(x)], ...
        [yLow; flipud(yHigh)], ...
        plotColor, ...
        'FaceAlpha', patchAlpha, ...
        'EdgeColor', 'none', ...
        'HandleVisibility', 'off');

    % Central profile.
    h = plot( ...
        ax, ...
        x, ...
        y, ...
        'LineStyle', plotStyle, ...
        'Marker', plotMarker, ...
        'LineWidth', 1.5, ...
        'MarkerSize', 5, ...
        'Color', [plotColor, lineAlpha], ...
        'MarkerFaceColor', plotColor, ...
        'DisplayName', displayName, ...
        'MarkerIndices', 1:3:numel(x));

end

function axisDir = caseAxisDirectionFromName(caseStr)
%CASEAXISDIRECTIONFROMNAME Return unit vector for case label.
%
% Examples:
%   '-X' -> [-1,  0,  0]
%   '+Y' -> [ 0, +1,  0]
%   '+Z' -> [ 0,  0, +1]
%   '-Z' -> [ 0,  0, -1]

    caseStr = strtrim(string(caseStr));

    switch upper(caseStr)
        case "+X"
            axisDir = [1, 0, 0];

        case "-X"
            axisDir = [-1, 0, 0];

        case "+Y"
            axisDir = [0, 1, 0];

        case "-Y"
            axisDir = [0, -1, 0];

        case "+Z"
            axisDir = [0, 0, 1];

        case "-Z"
            axisDir = [0, 0, -1];

        otherwise
            warning('Unknown case axis "%s". Using +Z.', caseStr);
            axisDir = [0, 0, 1];
    end

end