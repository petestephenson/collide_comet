function [f, ax] = plot_eDepth_dist(peCasesNo, swCasesNo, plasDensDir, ixXs)
%PLOT_EDEPTH_DIST Plot electron-depth distribution vs tau_e for case pairs.
%Adds a smaller top axis showing weighted median tau_e and horizontal IQR bars.
%
% Inputs
%   peCasesNo   cell array of strings, length >= 5
%   swCasesNo   cell array of strings, length >= 5
%   plasDensDir directory string (with or without trailing file separator)
%   ixXs        scalar index or vector of indices into X dimension
%
% Outputs
%   f   figure handle
%   ax  main (bottom) axes handle

    % ------------------- checks -------------------
    if nargin < 4 || isempty(ixXs)
        error('ixXs must be provided (scalar or vector index into X dimension).');
    end
    if numel(peCasesNo) < 5 || numel(swCasesNo) < 5
        error('peCasesNo and swCasesNo must have at least 5 entries.');
    end

    % Ensure trailing file separator
    if ~isempty(plasDensDir) && ~endsWith(plasDensDir, filesep)
        plasDensDir = [plasDensDir, filesep];
    end

    % ------------------- settings -------------------
    maxRThreshold = 200;
    asymName = {'-x', '+y', '+z', '-z', 'sym'};
    nCases = 5;

    % ------------------- figure/axes -------------------
    f = figure('Position', [0, 0, 500, 400]);

    % Manually-sized axes so top is smaller than bottom
    axTop = axes('Parent', f, 'Position', [0.12 0.73 0.83 0.20]); hold(axTop,'on');
    ax    = axes('Parent', f, 'Position', [0.12 0.12 0.83 0.55]); hold(ax,'on');

    set([ax axTop], 'XScale', 'log');
    grid(ax, 'off'); grid(axTop, 'off');

    ylabel(ax,    'Electron distribution');
    xlabel(ax,    '$\tau_e$', 'Interpreter', 'latex');
    ylabel(axTop, 'Median $\tau_e$ (IQR)', 'Interpreter', 'latex');

    % Top axis is categorical in Y, but shares log X with bottom
    ylim(axTop, [0.5 nCases+0.5]);
    set(axTop, 'YTick', 1:nCases, 'YTickLabel', asymName);
    set(axTop, 'TickDir', 'out');
    set(ax,    'TickDir', 'out');

    linkaxes([axTop ax], 'x');

    % Use consistent colors for both axes
    cols = lines(nCases);

    % ------------------- loop cases -------------------
    for iCase = 1:nCases
        peAsymName = [plasDensDir, 'plas_dens', peCasesNo{iCase}, '.nc'];
        swAsymName = [plasDensDir, 'plas_dens', swCasesNo{iCase}, '.nc'];

        % Load SW + PE and combine
        [~, eDepthGrid_sw, eDepthDist_sw] = asymPlots.load_electron_depth_dist(swAsymName);
        [eDepthRadGrid_pe, eDepthGrid_pe, eDepthDist_pe] = asymPlots.load_electron_depth_dist(peAsymName); %#ok<ASGLU>

%         eDepthDist_all = eDepthDist_sw + eDepthDist_pe;
        eDepthDist_all = eDepthDist_pe;

        % Radial selection (use PE radial grid; assume same as SW)
        eDepthRadMids = 0.5 * (eDepthRadGrid_pe(2:end) + eDepthRadGrid_pe(1:end-1));
        ixR = (eDepthRadMids < maxRThreshold);

        % Expect: (nTau, nX, nR, nZ, nSpecies)
        d = eDepthDist_all(:, ixXs, ixR, :, 1);

        % Integrate over X,R,Z
        w = squeeze(sum(d, [2 3 4]));
        w = w(:);

        % Tau midpoints (use PE tau grid consistently)
        tau = sqrt(eDepthGrid_pe(2:end) .* eDepthGrid_pe(1:end-1));
        tau = tau(:);

        % Align lengths safely
        n = min(numel(tau), numel(w));
        tau = tau(1:n);
        w   = w(1:n);

        % Drop first bin (to match your original plotting choice)
        if n >= 2
            tau = tau(2:end);
            w   = w(2:end);
        end

        % Normalise weights robustly
        w(~isfinite(w)) = 0;
        w(w < 0) = 0;
        sw = sum(w);
        if sw > 0
            w = w ./ sw;
        end

        % Weighted median + quartiles
        [medTau, q1Tau, q3Tau] = local_weighted_quantiles(tau, w, [0.50 0.25 0.75]);

        % -------- bottom: distribution (consistent color) --------
        plot(ax, tau, w, ...
            'Color', cols(iCase,:), ...
            'LineWidth', 1.5, ...
            'DisplayName', asymName{iCase});

        % -------- top: horizontal IQR bar + median (same color) --------
        y0 = iCase;
        plot(axTop, [q1Tau q3Tau], [y0 y0], ...
            'Color', cols(iCase,:), ...
            'LineWidth', 2);

        plot(axTop, medTau, y0, 'o', ...
            'Color', cols(iCase,:), ...
            'MarkerFaceColor', cols(iCase,:), ...
            'MarkerEdgeColor', cols(iCase,:), ...
            'MarkerSize', 6, ...
            'LineWidth', 1.0);
    end

    % Limits + legend
    xlim(ax, [1e-4 1e3]);
    legend(ax, 'Location', 'best');


    xline(ax,    1,'k--', 'LineWidth', 1.5);
    xline(axTop, 1, 'k--', 'LineWidth', 1.5);
    
    % Make top cleaner
    set(axTop, 'XTickLabel', []); % hide x tick labels on top axis
    box(axTop, 'off'); box(ax, 'off');
end

% ======================= local helper(s) ==========================
function [q50, q25, q75] = local_weighted_quantiles(x, w, probs)
% Compute weighted quantiles for x with nonnegative weights w.
% probs e.g. [0.50 0.25 0.75] (returned in that same order)

    x = x(:); w = w(:);
    mask = isfinite(x) & isfinite(w) & (w >= 0);
    x = x(mask); w = w(mask);

    if isempty(x) || sum(w) <= 0
        q50 = NaN; q25 = NaN; q75 = NaN;
        return;
    end

    [xs, idx] = sort(x);
    ws = w(idx);
    cdf = cumsum(ws) ./ sum(ws);

    qs = nan(size(probs));
    for i = 1:numel(probs)
        p = probs(i);
        j = find(cdf >= p, 1, 'first');
        if isempty(j)
            qs(i) = xs(end);
        else
            qs(i) = xs(j);
        end
    end

    q50 = qs(1);
    q25 = qs(2);
    q75 = qs(3);
end


% function [f, ax] = plot_eDepth_dist(peCasesNo, swCasesNo, plasDensDir, ixXs)
% %PLOT_EDEPTH_DIST Plot electron-depth distribution vs tau_e for case pairs.
% %
% % Inputs
% %   peCasesNo   cell array of strings, length >= 5
% %   swCasesNo   cell array of strings, length >= 5
% %   plasDensDir directory string (with or without trailing file separator)
% %   ixXs        scalar index or vector of indices into X dimension
% %
% % Outputs
% %   f   figure handle
% %   ax  axes handle
% 
%     if nargin < 4 || isempty(ixXs)
%         error('ixXs must be provided (scalar or vector index into X dimension).');
%     end
%     if numel(peCasesNo) < 5 || numel(swCasesNo) < 5
%         error('peCasesNo and swCasesNo must have at least 5 entries.');
%     end
% 
%     % Ensure trailing file separator
%     if ~isempty(plasDensDir) && ~endsWith(plasDensDir, filesep)
%         plasDensDir = [plasDensDir, filesep];
%     end
% 
% 
%     maxRThreshold = 200;
%     f = figure('Position', [0, 0, 600, 400]);
%     ax = axes('Parent', f);
%     hold(ax, 'on');
% %     Turn grid off
%     grid(ax, 'off');
%     xlabel(ax, '$\tau_e$', 'Interpreter', 'latex');
%     ylabel(ax, 'Electron distribution');
%     set(ax, 'XScale', 'log');
% %     set(ax, 'YScale', 'log');
%     asymName = {'-x', '+y', '+z', '-z', 'sym'};
%     for iCase = 1:5
%         peAsymNo = peCasesNo{iCase};
%         swAsymNo = swCasesNo{iCase};
% 
%         peAsymName = [plasDensDir, 'plas_dens', peAsymNo, '.nc'];
%         swAsymName = [plasDensDir, 'plas_dens', swAsymNo, '.nc'];
% 
%         % Correct: load SW from swAsymName, PE from peAsymName
%         [eDepthRadGrid_sw, eDepthGrid_sw, eDepthDist_sw] = asymPlots.load_electron_depth_dist(swAsymName);
%         [eDepthRadGrid_pe, eDepthGrid_pe, eDepthDist_pe] = asymPlots.load_electron_depth_dist(peAsymName);
% 
% %         eDepthDist_all = eDepthDist_sw;
%         eDepthDist_all = eDepthDist_sw + eDepthDist_pe;
%         eDepthRadMids = 0.5* ( eDepthRadGrid_pe(2:end) +  eDepthRadGrid_pe(1:end-1) );
%         ixR = eDepthRadMids< maxRThreshold;
%         % Expect something like: (nTau, nX, nY, nZ, nSpecies)
%         % Take species index 1, integrate over selected X and all Y,Z.
%         d = eDepthDist_all(:, ixXs, ixR, :, 1);
% 
%         % Sum over X,Y,Z dims robustly (dims 2,3,4 of d)
%         eDepth_plot = squeeze(sum(d, [2 3 4]));  % yields nTau x 1 (or 1 x nTau)
%         eDepth_plot = eDepth_plot / sum(eDepth_plot(2:end));
% 
%         % Midpoints of tau bins (use one grid consistently)
%         eDepthGrid = eDepthGrid_pe;
%         eDepthMids = sqrt(eDepthGrid(2:end) .* eDepthGrid(1:end-1));
% 
%         % Make sure vectors align
%         eDepth_plot = eDepth_plot(:);
%         eDepthMids = eDepthMids(:);
% 
%         n = min(numel(eDepthMids), numel(eDepth_plot));
%         plot(ax, eDepthMids(2:n), eDepth_plot(2:n), 'DisplayName', asymName{iCase}, 'LineWidth',1.5);
% 
%         
%     end
%     xlim(ax, [1e-4, 1e3]);
%     legend(ax, 'Location', 'best');
% %     grid(ax, 'on');
% end
