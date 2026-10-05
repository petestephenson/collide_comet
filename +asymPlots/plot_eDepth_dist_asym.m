function [f, axTop, axBot] = plot_eDepth_dist_asym(peCasesNo, swCasesNo, plasDensDir, ixXs)
%PLOT_EDEPTH_DIST_ASYM 3x2 panel: top=median+IQR, bottom=distributions.
% Splits into asymmetries along X, Y, Z and includes symmetric reference case in each.
% For each case, plots TOTAL plus the +/- split, with markers to distinguish them.

    % ------------------- checks -------------------
    if nargin < 4 || isempty(ixXs)
        error('ixXs must be provided (scalar or vector index into X dimension).');
    end
    if numel(peCasesNo) < 5 || numel(swCasesNo) < 5
        error('peCasesNo and swCasesNo must have at least 5 entries.');
    end
    if ~isempty(plasDensDir) && ~endsWith(plasDensDir, filesep)
        plasDensDir = [plasDensDir, filesep];
    end

    % ------------------- settings -------------------
    maxRThreshold = 200;

    asymName = {'-x', '+y', '+z', '-z', 'sym'};
    dimNo    = [ 1,    2,    3,    3,    1];

    groups(1).title   = 'X asymmetry';
    groups(1).caseIdx = [1 5];          % -x, sym
    groups(2).title   = 'Y asymmetry';
    groups(2).caseIdx = [2 5];          % +y, sym
    groups(3).title   = 'Z asymmetry';
    groups(3).caseIdx = [3 4 5];        % +z, -z, sym

    nCols = 3;

    % Marker + style scheme for the 3 curves per case
    mk.total = 'o';
    mk.plus  = 's';
    mk.minus = 'd';

    ls.total = '-';
    ls.plus  = '--';
    ls.minus = '-.';

    % ------------------- figure + manual layout -------------------
    f = figure('Position', [0, 0, 1200, 600]);

    left   = 0.07;
    right  = 0.02;
    top    = 0.06;
    bottom = 0.10;
    hgap   = 0.04;

    topH   = 0.22;
    botH   = 0.56;

    colW = (1 - left - right - (nCols-1)*hgap) / nCols;

    axTop = gobjects(1,nCols);
    axBot = gobjects(1,nCols);

    for c = 1:nCols
        x0 = left + (c-1)*(colW + hgap);

        axTop(c) = axes('Parent', f, 'Position', [x0, 1-top-topH, colW, topH]); hold(axTop(c),'on');
        axBot(c) = axes('Parent', f, 'Position', [x0, bottom,    colW, botH]); hold(axBot(c),'on');

        set([axTop(c) axBot(c)], 'XScale', 'log');
        box(axTop(c),'off'); box(axBot(c),'off');
        set([axTop(c) axBot(c)], 'TickDir', 'out');

        title(axTop(c), groups(c).title);

        ylabel(axTop(c), 'Median \tau_e (IQR)', 'Interpreter','tex');
        ylabel(axBot(c), 'Electron distribution');
        xlabel(axBot(c), '$\tau_e$', 'Interpreter','latex');

        set(axTop(c), 'XTickLabel', []); % cleaner top axis
        linkaxes([axTop(c) axBot(c)], 'x');

        xline(axTop(c), 1, 'k--', 'LineWidth', 1.2);
        xline(axBot(c), 1, 'k--', 'LineWidth', 1.2);
    end

    % consistent colors across all cases (1..5)
    colsAll = lines(5);

    % ------------------- compute + plot -------------------
    for c = 1:nCols
        caseIdx  = groups(c).caseIdx;
        nInGroup = numel(caseIdx);

        ylim(axTop(c), [0.5 nInGroup+0.5]);
        set(axTop(c), 'YTick', 1:nInGroup, 'YTickLabel', asymName(caseIdx));

        % vertical offsets so TOTAL / + / - are separated on the top axis
        dyT = 0.20;   % total
        dyP = 0.00;   % plus
        dyM = -0.20;  % minus

        for k = 1:nInGroup
            iCase = caseIdx(k);

            peAsymName = [plasDensDir, 'plas_dens', peCasesNo{iCase}, '.nc'];
            swAsymName = [plasDensDir, 'plas_dens', swCasesNo{iCase}, '.nc'];

            % Load SW + PE and combine
            [~, eDepthGrid_sw, eDepthDist_sw] = asymPlots.load_electron_depth_dist(swAsymName); %#ok<ASGLU>
            [eDepthRadGrid_pe, eDepthGrid_pe, eDepthDist_pe] = asymPlots.load_electron_depth_dist(peAsymName); %#ok<ASGLU>

%             eDepthDist_all = eDepthDist_sw + eDepthDist_pe;
            eDepthDist_all = eDepthDist_pe;

            % Radial selection
            eDepthRadMids = 0.5 * (eDepthRadGrid_pe(2:end) + eDepthRadGrid_pe(1:end-1));
            ixR = (eDepthRadMids < maxRThreshold);

            nTheta = size(eDepthDist_all, 4);
            nTHalf = floor(nTheta / 2);

            % Expect: (nTau, nX, nR, nTheta, nSpecies)
            d1 = eDepthDist_all(:, ixXs, ixR, 1:nTHalf,      dimNo(iCase));
            d2 = eDepthDist_all(:, ixXs, ixR, nTHalf+1:end,  dimNo(iCase));
%             d1 = eDepthDist_all(:, ixXs, ixR, 1:nTHalf,      1);
%             d2 = eDepthDist_all(:, ixXs, ixR, nTHalf+1:end,  1);

            % Integrate over X,R,Theta  -> RAW weights
            w1raw = squeeze(sum(d1, [2 3 4])); w1raw = w1raw(:);
            w2raw = squeeze(sum(d2, [2 3 4])); w2raw = w2raw(:);
            wtraw = w1raw + w2raw;

            % Tau midpoints
            tau = sqrt(eDepthGrid_pe(2:end) .* eDepthGrid_pe(1:end-1));
            tau = tau(:);

            % Align lengths safely
            n = min([numel(tau), numel(w1raw), numel(w2raw), numel(wtraw)]);
            tau   = tau(1:n);
            w1raw = w1raw(1:n);
            w2raw = w2raw(1:n);
            wtraw = wtraw(1:n);

            % Drop first bin
            if n >= 2
                tau   = tau(2:end);
                w1raw = w1raw(2:end);
                w2raw = w2raw(2:end);
                wtraw = wtraw(2:end);
            end

            % Clean / clamp
            w1raw(~isfinite(w1raw)) = 0; w1raw(w1raw < 0) = 0;
            w2raw(~isfinite(w2raw)) = 0; w2raw(w2raw < 0) = 0;
            wtraw(~isfinite(wtraw)) = 0; wtraw(wtraw < 0) = 0;

            % Normalise EACH distribution (so each integrates to 1)
            w1 = w1raw; s1 = sum(w1); if s1 > 0, w1 = w1 ./ s1; end
            w2 = w2raw; s2 = sum(w2); if s2 > 0, w2 = w2 ./ s2; end
            wt = wtraw; st = sum(wt); if st > 0, wt = wt ./ st; end

            % Weighted median + quartiles (TOTAL + +/-)
            [med_t, q1_t, q3_t] = local_weighted_quantiles(tau, wt, [0.50 0.25 0.75]);
            [med_p, q1_p, q3_p] = local_weighted_quantiles(tau, w1, [0.50 0.25 0.75]);
            [med_m, q1_m, q3_m] = local_weighted_quantiles(tau, w2, [0.50 0.25 0.75]);

            cc = colsAll(iCase,:);

            % ---------- bottom: 3 distributions per case ----------
            plot(axBot(c), tau, wt, ...
                'Color', cc, 'LineWidth', 1.6, 'LineStyle', ls.total, ...
                'Marker', mk.total, 'MarkerSize', 4, ...
                'DisplayName', sprintf('%s, total', asymName{iCase}));

            plot(axBot(c), tau, w1, ...
                'Color', cc, 'LineWidth', 1.4, 'LineStyle', ls.plus, ...
                'Marker', mk.plus, 'MarkerSize', 4, ...
                'DisplayName', sprintf('%s, >0', asymName{iCase}));

            plot(axBot(c), tau, w2, ...
                'Color', cc, 'LineWidth', 1.4, 'LineStyle', ls.minus, ...
                'Marker', mk.minus, 'MarkerSize', 4, ...
                'DisplayName', sprintf('%s, <0', asymName{iCase}));

            % ---------- top: 3 IQR bars + 3 medians ----------
            y0 = k;

            % TOTAL (solid IQR)
            plot(axTop(c), [q1_t q3_t], [y0+dyT y0+dyT], 'Color', cc, 'LineWidth', 2.2, 'LineStyle', '-');
            plot(axTop(c), med_t, y0+dyT, mk.total, 'Color', cc, ...
                'MarkerFaceColor', cc, 'MarkerEdgeColor', cc, 'MarkerSize', 6);

            % + (dashed IQR)
            plot(axTop(c), [q1_p q3_p], [y0+dyP y0+dyP], 'Color', cc, 'LineWidth', 2.0, 'LineStyle', '--');
            plot(axTop(c), med_p, y0+dyP, mk.plus, 'Color', cc, ...
                'MarkerFaceColor', cc, 'MarkerEdgeColor', cc, 'MarkerSize', 6);

            % - (dotted IQR)
            plot(axTop(c), [q1_m q3_m], [y0+dyM y0+dyM], 'Color', cc, 'LineWidth', 2.0, 'LineStyle', ':');
            plot(axTop(c), med_m, y0+dyM, mk.minus, 'Color', cc, ...
                'MarkerFaceColor', 'w', 'MarkerEdgeColor', cc, 'MarkerSize', 6);
        end

        % Limits + legend
        xlim(axBot(c), [1e-4 1e3]);
        xlim(axTop(c), [1e-4 1e3]);

        legend(axBot(c), 'Location', 'best');
        grid(axBot(c), 'off');
        grid(axTop(c), 'off');
    end
end

% ======================= local helper(s) ==========================
function [q50, q25, q75] = local_weighted_quantiles(x, w, probs)
% Compute weighted quantiles for x with nonnegative weights w.

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
