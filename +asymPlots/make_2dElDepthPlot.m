function [f, ax, cb, out] = make_2dElDepthPlot( ...
    elDepthRadEdges, elDepthGrid, elDepthDist, iDim, f, ax, clims, opts)
%MAKE_2DELDEPTHPLOT Plot log10(medianDepth) in polar (x,y) and overlay a SMOOTHED convex-hull exobase.
%
% Exobase mask: log10(medianDepth) > opts.exoThreshLog10 (default 0)
% Exobase outline: convex hull of masked points, then periodic smoothing of the hull polygon.

    if nargin < 8 || isempty(opts), opts = struct(); end
    if ~isfield(opts, "iTau"),           opts.iTau = 3; end
    if ~isfield(opts, "medianLevel"),    opts.medianLevel = 0.5; end
    if ~isfield(opts, "exoThreshLog10"), opts.exoThreshLog10 = 0; end
    if ~isfield(opts, "rMax"),           opts.rMax = max(elDepthRadEdges(:)); end
    if ~isfield(opts, "cmap"),           opts.cmap = []; end
    if ~isfield(opts, "showAngleRays"),  opts.showAngleRays = true; end

    % Smoothed hull options
    if ~isfield(opts, "smoothHull"),     opts.smoothHull = true; end
    if ~isfield(opts, "smoothMethod"),   opts.smoothMethod = "sgolay"; end % "movmean" | "sgolay"
    if ~isfield(opts, "smoothSpan"),     opts.smoothSpan = 14; end          % odd recommended
    if ~isfield(opts, "sgolayOrder"),    opts.sgolayOrder = 2; end
    if ~isfield(opts, "hullResampleN"),  opts.hullResampleN = 50; end     % number of points on smoothed curve

    if nargin < 5 || isempty(f) || ~ishandle(f)
        f = figure('Position', [10, 10, 1000, 650]);
    end
    if nargin < 6 || isempty(ax) || ~ishandle(ax)
        ax = axes(f);
    end
    if nargin < 7 || isempty(clims)
        clims = [-2, 2];
    end

    % ---- Compute medianDepth(r,theta) ----
    [medianDepth, rMids, thetaMids] = compute_medianDepth_grid( ...
        elDepthRadEdges, elDepthGrid, elDepthDist, iDim, opts.iTau, opts.medianLevel);

    log10MedianDepth = log10(medianDepth);
    inExo = log10MedianDepth > opts.exoThreshLog10;

    % ---- Mesh (x,y) ----
    [theta, rho] = meshgrid(thetaMids, rMids);
    x = rho .* cos(theta);
    y = rho .* sin(theta);

    % ---- Plot field ----
    cla(ax); hold(ax, "on");
    h = pcolor(ax, x, y, log10MedianDepth);
    set(h, 'EdgeColor', 'none');
    axis(ax, "equal");
    xlim(ax, [-opts.rMax, opts.rMax]);
    ylim(ax, [0, opts.rMax]);

    xt = ax.XTick;
    ax.XTick = xt(xt >= 0);
    xlabel(ax, 'x [km]');
    ax.YAxis.Visible = 'off';
    ax.Box = 'off';
    ax.Color = 'none';

    if ~isempty(opts.cmap)
        colormap(ax, opts.cmap);
    else
        try, colormap(ax, redblue); catch, end
    end
    clim(ax, clims);

    % ---- Smoothed convex-hull exobase ----
    exo = compute_smoothed_exobase_hull_xy(x, y, inExo, opts);
    if ~isempty(exo)
        plot(ax, exo.x_smooth, exo.y_smooth, 'k-', 'LineWidth', 1.7);
    end

    % ---- Optional: angle rays ----
    if opts.showAngleRays
        angleTicks  = linspace(0, pi, 7);
        angleLabels = {'0°','30°','60°','90°','120°','150°','180°'};
        for i = 1:numel(angleTicks)
            ang = angleTicks(i);
            xline = opts.rMax * cos(ang);
            yline = opts.rMax * sin(ang);
            plot(ax, [-xline, xline], [-yline, yline], '--', 'Color', [0.2 0.2 0.2]);
            text(1.08*xline - abs(0.05*xline), 1.08*yline, angleLabels{i}, ...
                 'Color', [0 0 0], 'Parent', ax);
        end
    end

    % ---- Colorbar ----
    cb = colorbar(ax, 'Position', [0.90, 0.15, 0.02, 0.80]);
    cb.Label.String = 'log_{10}(Median depth)';

    % ---- Outputs ----
    out = struct();
    out.medianDepth      = medianDepth;
    out.log10MedianDepth = log10MedianDepth;
    out.inExo            = inExo;
    out.rMids            = rMids;
    out.thetaMids        = thetaMids;
    out.x                = x;
    out.y                = y;
    out.exo              = exo; % includes raw + smoothed hull curves (may be [])
end


function [medianDepth, rMids, thetaMids] = compute_medianDepth_grid( ...
    elDepthRadEdges, elDepthGrid, elDepthDist, iDim, iTau, medianLevel)

    rMids = 0.5 * (elDepthRadEdges(1:end-1) + elDepthRadEdges(2:end));

    nTheta = size(elDepthDist, 4);
    thetaBounds = acos(linspace(-1, 1, nTheta + 1));
    thetaMids = 0.5 * (thetaBounds(1:end-1) + thetaBounds(2:end));

    slab = squeeze(elDepthDist(2:end, iTau, :, :, iDim));  % (nDepth x nR x nTheta)
    cumSlab = cumsum(slab, 1);
    denom = cumSlab(end, :, :);
    denom(denom == 0) = NaN;
    frac = cumSlab ./ denom;

    nR = size(frac, 2);
    medianDepth = NaN(nR, nTheta);

    for it = 1:nTheta
        for ir = 1:nR
            ix = find(frac(:, ir, it) >= medianLevel, 1, 'first');
            if ~isempty(ix)
                medianDepth(ir, it) = elDepthGrid(ix + 1); % +1 because slab started at 2:end
            end
        end
    end
end


function exo = compute_smoothed_exobase_hull_xy(x, y, inExo, opts)
%COMPUTE_SMOOTHED_EXOBASE_HULL_XY
% 1) Take convex hull of inExo points.
% 2) Convert hull polygon to arc-length parameter.
% 3) Resample to uniform arc-length spacing.
% 4) Apply periodic smoothing (movmean or Savitzky-Golay).
%
% Returns [] if insufficient points.

    xs = x(inExo);
    ys = y(inExo);

    if numel(xs) < 3
        exo = [];
        return;
    end

    good = isfinite(xs) & isfinite(ys);
    xs = xs(good); ys = ys(good);
    if numel(xs) < 3
        exo = [];
        return;
    end

    P = [xs(:), ys(:)];
    if rank(P - mean(P,1)) < 2
        exo = [];
        return;
    end

    k = convhull(xs, ys); % closed polygon (last repeats first)
    xr = xs(k);
    yr = ys(k);

    % --- Arc-length parameter along the closed hull polygon ---
    dx = diff(xr);
    dy = diff(yr);
    ds = hypot(dx, dy);
    s = [0; cumsum(ds)];          % s length = numel(xr)
    L = s(end);
    if ~isfinite(L) || L == 0
        exo = [];
        return;
    end

    % --- Resample to uniform arc length ---
    N = max(50, opts.hullResampleN);
    sU = linspace(0, L, N+1)';    % include closure point
    xU = interp1(s, xr, sU, 'linear');
    yU = interp1(s, yr, sU, 'linear');

    % Drop the duplicated last point for smoothing, then re-close later
    xU = xU(1:end-1);
    yU = yU(1:end-1);

    % --- Periodic smoothing: wrap the data so endpoints blend ---
    span = max(3, round(opts.smoothSpan));
    if mod(span,2) == 0, span = span + 1; end % keep odd
    halfw = floor(span/2);

    xW = [xU(end-halfw+1:end); xU; xU(1:halfw)];
    yW = [yU(end-halfw+1:end); yU; yU(1:halfw)];

    method = string(opts.smoothMethod);
    switch lower(method)
        case "sgolay"
            % Requires Signal Processing Toolbox for sgolayfilt.
            try
                ord = min(opts.sgolayOrder, span-2);
                xWs = sgolayfilt(xW, ord, span);
                yWs = sgolayfilt(yW, ord, span);
            catch
                % Fallback to movmean if sgolayfilt unavailable
                xWs = movmean(xW, span);
                yWs = movmean(yW, span);
            end
        otherwise % "movmean"
            xWs = movmean(xW, span);
            yWs = movmean(yW, span);
    end

    % Extract the centre (original-length) portion and re-close
    xS = xWs(halfw+1 : halfw+numel(xU));
    yS = yWs(halfw+1 : halfw+numel(yU));

    xS(end+1) = xS(1);
    yS(end+1) = yS(1);

    exo = struct();
    exo.x_raw = xr;
    exo.y_raw = yr;
    exo.x_smooth = xS;
    exo.y_smooth = yS;
    exo.k = k;
end


% function [f, ax, cb, out] = make_2dElDepthPlot( ...
%     elDepthRadEdges, elDepthGrid, elDepthDist, iDim, f, ax, clims, opts)
% %MAKE_2DELDEPTHPLOT Plot log10(medianDepth) in polar (x,y) and overlay a convex-hull exobase.
% %
% % Exobase overlay is computed as the convex hull of points where
% % log10(medianDepth) > opts.exoThreshLog10 (default 0).
% 
%     if nargin < 8 || isempty(opts), opts = struct(); end
%     if ~isfield(opts, "iTau"),           opts.iTau = 3; end
%     if ~isfield(opts, "medianLevel"),    opts.medianLevel = 0.5; end
%     if ~isfield(opts, "exoThreshLog10"), opts.exoThreshLog10 = 0; end
%     if ~isfield(opts, "rMax"),           opts.rMax = max(elDepthRadEdges(:)); end
%     if ~isfield(opts, "cmap"),           opts.cmap = []; end
%     if ~isfield(opts, "showAngleRays"),  opts.showAngleRays = true; end
% 
%     if nargin < 5 || isempty(f) || ~ishandle(f)
%         f = figure('Position', [10, 10, 1000, 650]);
%     end
%     if nargin < 6 || isempty(ax) || ~ishandle(ax)
%         ax = axes(f);
%     end
%     if nargin < 7 || isempty(clims)
%         clims = [-2, 2];
%     end
% 
%     % ---- Compute medianDepth(r,theta) ----
%     [medianDepth, rMids, thetaMids] = compute_medianDepth_grid( ...
%         elDepthRadEdges, elDepthGrid, elDepthDist, iDim, opts.iTau, opts.medianLevel);
% 
%     log10MedianDepth = log10(medianDepth);
%     inExo = log10MedianDepth > opts.exoThreshLog10;
% 
%     % ---- Mesh (x,y) ----
%     [theta, rho] = meshgrid(thetaMids, rMids);
%     x = rho .* cos(theta);
%     y = rho .* sin(theta);
% 
%     % ---- Plot field ----
%     cla(ax); hold(ax, "on");
%     % --- Edges (not mids) ---
%     rEdges = elDepthRadEdges(:);          % (nR+1)x1
%     nTheta = size(elDepthDist, 4);
%     thetaEdges = acos(linspace(-1, 1, nTheta+1));  % (1)x(nTheta+1), spans 0..pi
%     
%     % --- Make edge mesh in polar then to Cartesian ---
%     [thE, rE] = meshgrid(thetaEdges, rEdges);      % (nR+1)x(nTheta+1)
%     xE = rE .* cos(thE);
%     yE = rE .* sin(thE);
%     
%     % --- Bin values: one per (r,theta) bin ---
%     C = log10MedianDepth;                 % (nR)x(nTheta) from your computation
%     
%     % --- Plot as flat bins that fill the domain ---
%     Cpad = C;
%     Cpad(end+1, :) = C(end, :);      % duplicate last radial bin
%     Cpad(:, end+1) = Cpad(:, end);   % duplicate last theta bin
% 
%     h = pcolor(ax, xE, yE, Cpad);
%     set(h, 'EdgeColor', 'none');
%     shading(ax, 'flat');                  % <- constant per bin (no interp)
%     
%     axis(ax, 'equal');
%     xlim(ax, [-opts.rMax, opts.rMax]);
%     ylim(ax, [0, opts.rMax]);
% 
%     xt = ax.XTick;
%     ax.XTick = xt(xt >= 0);
%     xlabel(ax, 'x [km]');
%     ax.YAxis.Visible = 'off';
%     ax.Box = 'off';
%     ax.Color = 'none';
% 
%     if ~isempty(opts.cmap)
%         colormap(ax, opts.cmap);
%     else
%         try
%             colormap(ax, redblue);
%         catch
%         end
%     end
%     clim(ax, clims);
% 
%     % ---- Convex-hull exobase (smooth outline) ----
%     [exoHull, xs, ys] = compute_exobase_hull_xy(x, y, inExo);
%     if ~isempty(exoHull)
%         plot(ax, exoHull.x, exoHull.y, 'k-', 'LineWidth', 0.5);
%     end
%     scatter(ax, xs, ys, 5, 'k', 'filled', 'MarkerFaceAlpha', 1);
% 
%     % ---- Optional: angle rays ----
%     if opts.showAngleRays
%         angleTicks  = linspace(0, pi, 7);
%         angleLabels = {'0°','30°','60°','90°','120°','150°','180°'};
%         for i = 1:numel(angleTicks)
%             ang = angleTicks(i);
%             xline = opts.rMax * cos(ang);
%             yline = opts.rMax * sin(ang);
%             plot(ax, [-xline, xline], [-yline, yline], '--', 'Color', [0.2 0.2 0.2]);
%             text(1.08*xline - abs(0.05*xline), 1.08*yline, angleLabels{i}, ...
%                  'Color', [0 0 0], 'Parent', ax);
%         end
%     end
% 
%     % ---- Colorbar ----
%     cb = colorbar(ax, 'Position', [0.90, 0.15, 0.02, 0.80]);
%     cb.Label.String = 'log_{10}(Median \tau)';
% 
%     % ---- Outputs ----
%     out = struct();
%     out.medianDepth      = medianDepth;
%     out.log10MedianDepth = log10MedianDepth;
%     out.inExo            = inExo;
%     out.rMids            = rMids;
%     out.thetaMids        = thetaMids;
%     out.x                = x;
%     out.y                = y;
%     out.exoHull          = exoHull; % empty if insufficient points
% end
% 
% 
% function [medianDepth, rMids, thetaMids] = compute_medianDepth_grid( ...
%     elDepthRadEdges, elDepthGrid, elDepthDist, iDim, iTau, medianLevel)
% 
%     rMids = 0.5 * (elDepthRadEdges(1:end-1) + elDepthRadEdges(2:end));
% 
%     nTheta = size(elDepthDist, 4);
%     thetaBounds = acos(linspace(-1, 1, nTheta + 1));
%     thetaMids = 0.5 * (thetaBounds(1:end-1) + thetaBounds(2:end));
% 
%     slab = squeeze(elDepthDist(2:end, iTau, :, :, iDim));  % (nDepth x nR x nTheta)
%     cumSlab = cumsum(slab, 1);
%     denom = cumSlab(end, :, :);
%     denom(denom == 0) = NaN;
%     frac = cumSlab ./ denom;
% 
%     nR = size(frac, 2);
%     medianDepth = NaN(nR, nTheta);
% 
%     for it = 1:nTheta
%         for ir = 1:nR
%             ix = find(frac(:, ir, it) >= medianLevel, 1, 'first');
%             if ~isempty(ix)
%                 medianDepth(ir, it) = elDepthGrid(ix + 1); % +1 because slab started at 2:end
%             end
%         end
%     end
% end
% 
% 
% function [exoHull, xs, ys] = compute_exobase_hull_xy(x, y, inExo)
% %COMPUTE_EXOBASE_HULL_XY Convex hull of in-exobase points in (x,y).
% %
% % Returns struct with fields:
% %   .x, .y  : hull polygon coordinates (closed)
% %   .k      : hull indices into the in-exo point list
% %
% % Returns [] if not enough points.
% 
%     xs = x(inExo);
%     ys = y(inExo);
% 
%     % Need at least 3 non-collinear points for a hull
%     if numel(xs) < 3
%         exoHull = [];
%         return;
%     end
% 
%     % Remove NaNs/Infs
%     good = isfinite(xs) & isfinite(ys);
%     xs = xs(good);
%     ys = ys(good);
%     if numel(xs) < 3
%         exoHull = [];
%         return;
%     end
% 
%     % If points are collinear (rank < 2), convhull can error; guard it.
%     P = [xs(:), ys(:)];
%     if rank(P - mean(P,1)) < 2
%         exoHull = [];
%         return;
%     end
% 
%     k = convhull(xs, ys); % k closes the polygon (last == first)
%     exoHull = struct();
%     exoHull.k = k;
%     exoHull.x = xs(k);
%     exoHull.y = ys(k);
% end
% 
% 
% 
% % function [f, ax, cb, out] = make_2dElDepthPlot( ...
% %     elDepthRadEdges, elDepthGrid, elDepthDist, iDim, f, ax, clims, opts)
% % %MAKE_2DELDEPTHPLOT Plot median electron-depth metric in a polar (x,y) view + exobase contour.
% % %
% % % Inputs
% % %   elDepthRadEdges : (nR+1) radial bin edges [km]
% % %   elDepthGrid     : (nE)   grid of electron-depth values (same indexing as elDepthDist's depth axis)
% % %   elDepthDist     : array with dims at least: (depthIdx, tauIdx, rIdx, thetaIdx, dimIdx)
% % %   iDim            : which dimIdx to plot (scalar)
% % %   f, ax           : optional existing figure/axes handles (pass [] to create)
% % %   clims           : [cmin cmax] color limits for log10(medianDepth)
% % %   opts            : optional struct, fields:
% % %       .iTau          (default 3)   index along tauIdx
% % %       .medianLevel   (default 0.5) target cumulative fraction
% % %       .exoThreshLog10(default 0)   exobase threshold: log10(medianDepth) > this
% % %       .rMax          (default max(elDepthRadEdges)) for plot limits / angle rays
% % %       .cmap          (default [])  colormap name/array; [] keeps current or uses redblue if available
% % %       .showAngleRays (default true)
% % %
% % % Outputs
% % %   f, ax, cb : handles
% % %   out       : struct with computed products:
% % %       .medianDepth       (nR x nTheta)
% % %       .log10MedianDepth  (nR x nTheta)
% % %       .inExo             (nR x nTheta) logical
% % %       .rMids             (nR x 1)
% % %       .thetaMids         (1 x nTheta)
% % %       .x, .y             (nR x nTheta) mesh
% % %       .exo               exobase curve in polar and cart:
% % %           .rho (1 x nTheta), .theta (1 x nTheta), .x (1 x nTheta), .y (1 x nTheta)
% % 
% %     if nargin < 8 || isempty(opts), opts = struct(); end
% %     if ~isfield(opts, "iTau"),           opts.iTau = 3; end
% %     if ~isfield(opts, "medianLevel"),    opts.medianLevel = 0.5; end
% %     if ~isfield(opts, "exoThreshLog10"), opts.exoThreshLog10 = 0; end
% %     if ~isfield(opts, "rMax"),           opts.rMax = max(elDepthRadEdges(:)); end
% %     if ~isfield(opts, "cmap"),           opts.cmap = []; end
% %     if ~isfield(opts, "showAngleRays"),  opts.showAngleRays = true; end
% % 
% %     if nargin < 5 || isempty(f) || ~ishandle(f)
% %         f = figure('Position', [10, 10, 1000, 650]);
% %     end
% %     if nargin < 6 || isempty(ax) || ~ishandle(ax)
% %         ax = axes(f);
% %     end
% %     if nargin < 7 || isempty(clims)
% %         clims = [-2, 2];
% %     end
% % 
% %     % --------- Compute medianDepth (nR x nTheta) ----------
% %     [medianDepth, rMids, thetaMids] = compute_medianDepth_grid( ...
% %         elDepthRadEdges, elDepthGrid, elDepthDist, iDim, opts.iTau, opts.medianLevel);
% % 
% %     log10MedianDepth = log10(medianDepth);
% %     inExo = log10MedianDepth > opts.exoThreshLog10;
% % 
% %     % --------- Build polar mesh (x,y) ----------
% %     [theta, rho] = meshgrid(thetaMids, rMids);
% %     x = rho .* cos(theta);
% %     y = rho .* sin(theta);
% % 
% %     % --------- Plot ----------
% %     cla(ax);
% %     hold(ax, "on");
% % 
% %     h = pcolor(ax, x, y, log10MedianDepth);
% %     set(h, 'EdgeColor', 'none'); % cleaner than shading interp for pcolor grids
% %     axis(ax, "equal");
% %     xlim(ax, [-opts.rMax, opts.rMax]);
% %     ylim(ax, [0, opts.rMax]);
% % 
% %     % Use only non-negative x-ticks (matching your original)
% %     xt = ax.XTick;
% %     ax.XTick = xt(xt >= 0);
% % 
% %     xlabel(ax, 'x [km]');
% %     % y-axis is radial up in your plot; you hid it originally:
% %     ax.YAxis.Visible = 'off';
% %     ax.Box = 'off';
% %     ax.Color = 'none';
% % 
% %     % Colormap handling (keep your redblue preference if available)
% %     if ~isempty(opts.cmap)
% %         colormap(ax, opts.cmap);
% %     else
% %         try
% %             colormap(ax, redblue);
% %         catch
% %             % fallback: leave as-is
% %         end
% %     end
% %     clim(ax, clims);
% % 
% %     % --------- Exobase curve as distance vs (r, theta) ----------
% %     exo = compute_exobase_distance_r_theta(inExo, rMids, thetaMids);
% % 
% %     % Plot exobase curve (only where defined)
% %     valid = ~isnan(exo.rho);
% %     plot(ax, exo.x(valid), exo.y(valid), 'k-', 'LineWidth', 1.5);
% % 
% %     % --------- Optional: angle rays + labels ----------
% %     if opts.showAngleRays
% %         angleTicks  = linspace(0, pi, 7);
% %         angleLabels = {'0°','30°','60°','90°','120°','150°','180°'};
% %         for i = 1:numel(angleTicks)
% %             ang = angleTicks(i);
% %             xline = opts.rMax * cos(ang);
% %             yline = opts.rMax * sin(ang);
% %             plot(ax, [-xline, xline], [-yline, yline], '--', 'Color', [0.2 0.2 0.2]);
% %             text(1.08*xline - abs(0.05*xline), 1.08*yline, angleLabels{i}, 'Color', [0 0 0], ...
% %                  'Parent', ax);
% %         end
% %     end
% % 
% %     % --------- Colorbar ----------
% %     cb = colorbar(ax, 'Position', [0.90, 0.15, 0.02, 0.80]);
% %     cb.Label.String = 'log_{10}(Median \tau)';
% % 
% %     % --------- Pack outputs ----------
% %     out = struct();
% %     out.medianDepth      = medianDepth;
% %     out.log10MedianDepth = log10MedianDepth;
% %     out.inExo            = inExo;
% %     out.rMids            = rMids;
% %     out.thetaMids        = thetaMids;
% %     out.x                = x;
% %     out.y                = y;
% %     out.exo              = exo;
% % end
% % 
% % 
% % function [medianDepth, rMids, thetaMids] = compute_medianDepth_grid( ...
% %     elDepthRadEdges, elDepthGrid, elDepthDist, iDim, iTau, medianLevel)
% % %COMPUTE_MEDIANDEPTH_GRID Convert elDepthDist -> medianDepth(r,theta) via cumulative fraction threshold.
% % 
% %     % Radial midpoints
% %     rMids = 0.5 * (elDepthRadEdges(1:end-1) + elDepthRadEdges(2:end));
% % 
% %     % Theta midpoints: assume theta bounds are sampled uniformly in cos(theta) from -1..1
% %     nTheta = size(elDepthDist, 4);
% %     thetaBounds = acos(linspace(-1, 1, nTheta + 1));
% %     thetaMids = 0.5 * (thetaBounds(1:end-1) + thetaBounds(2:end)); % 1 x nTheta
% % 
% %     % Pull the slice we need: (depthIdx, rIdx, thetaIdx)
% %     % Expecting elDepthDist(depthIdx, tauIdx, rIdx, thetaIdx, dimIdx)
% %     slab = squeeze(elDepthDist(2:end, iTau, :, :, iDim)); % drop first depth bin to match your original
% %     % slab: (nDepth x nR x nTheta)
% % 
% %     % Cumulative fraction along depthIdx
% %     cumSlab = cumsum(slab, 1);
% %     denom = cumSlab(end, :, :);
% %     % Avoid divide-by-zero
% %     denom(denom == 0) = NaN;
% %     frac = cumSlab ./ denom;
% % 
% %     % Find first depth index where frac >= medianLevel, for each (r,theta)
% %     nR = size(frac, 2);
% %     medianDepth = NaN(nR, nTheta);
% % 
% %     for it = 1:nTheta
% %         for ir = 1:nR
% %             ix = find(frac(:, ir, it) >= medianLevel, 1, 'first');
% %             if ~isempty(ix)
% %                 % ix indexes elDepthGrid(2:end) because slab started at 2:end
% %                 medianDepth(ir, it) = elDepthGrid(ix + 1);
% %             end
% %         end
% %     end
% % end
% % 
% % 
% % function exo = compute_exobase_distance_r_theta(inExo, rMids, thetaMids)
% % %COMPUTE_EXOBASE_DISTANCE_R_THETA Compute exobase radius (distance) vs theta.
% % %
% % % Defines the "exobase distance" at each theta as the OUTERMOST radius (largest r)
% % % for which inExo(r,theta) is true. If none are true, returns NaN for that theta.
% % 
% %     nTheta = numel(thetaMids);
% %     rho = NaN(1, nTheta);
% % 
% %     for it = 1:nTheta
% %         idx = find(inExo(:, it), 1, 'last'); % outermost in-exobase point
% %         if ~isempty(idx)
% %             rho(it) = rMids(idx);
% %         end
% %     end
% % 
% %     exo = struct();
% %     exo.theta = reshape(thetaMids, 1, []);
% %     exo.rho   = rho;
% %     exo.x     = rho .* cos(exo.theta);
% %     exo.y     = rho .* sin(exo.theta);
% % end
% % 
% % 
% % 
% % % function [f,ax,cb] = make_2dElDepthPlot(elDepthRadGrid, elDepthGrid, elDepthDist,iDim, f, ax, clims)
% % % % function [f,ax, xInExo] = plot_el_depth_2d(ElDepthDist_all, iDim)
% % % % f = figure('Position',[10,10,1200,400]);
% % % % nRuns = numel(ElDepthDist_all);
% % % %     lef = [0.08,0.35,0.62];s
% % % % lef = linspace(0.08, 1-0.5*((2+0.5)/nRuns ), nRuns);
% % % % lef = [0.08];
% % % % caseNameList = {'Ref','-X', '+Y', '+Z', '-Z'};
% % % % axisName = {'X', 'Y', 'Z'};
% % % bot = 0.15;
% % % % wid = 1/(nRuns+2);
% % % heig = 0.8;
% % % % ax = cell(1,nRuns);
% % % % nRad = size(ElDepthDist_all{1}, 3);
% % % 
% % % nElDepthThetaBounds = size(elDepthDist,4)+1;
% % % % nRuns = numel(ElDepthDist_all);
% % % maxrProd = 300;
% % % % elDepthRadGrid = linspace(0, maxrProd, nRad + 1);
% % % elDepthRadMids = 0.5* (elDepthRadGrid(1:end-1)+elDepthRadGrid(2:end));
% % % 
% % % iTau = 3;
% % % % ixExobase = cell(nRuns,1);
% % % % for iOGR = 1:nRuns
% % % elDepthDist2 = elDepthDist;
% % % %     elDepthDist2 = sum(elDepthDist2, [4]);
% % % medianTau = cumsum(elDepthDist2(2:end,iTau,:,:,iDim), 1);
% % % medianTau = squeeze(medianTau ./medianTau(end,:,:,:));
% % % % elDepthGrid = linspace(-5,4,51);
% % % % elDepthGrid = elDepthGrid(2:end);
% % % [r,c,v] = find(medianTau >= 0.5);
% % % rowIndices = nan(size(elDepthDist2, [3,4]));
% % %     % elDe
% % % %     for iDim = 1:3
% % % for iPhi = 1:size(medianTau, 3)
% % %     for col = 1:size(medianTau, 2)
% % %         ixMedian = find(squeeze(medianTau(:, col, iPhi)) >= 0.5, 1, 'first');
% % %         if ~isempty(ixMedian)
% % %             rowIndices(col, iPhi) = ixMedian;
% % %     %             valIndices(col) = v;
% % %         end
% % %     end
% % % end
% % % 
% % % %     end
% % % elDepthMedians = nan(size(rowIndices));
% % % elDepthMedians(rowIndices>0) = elDepthGrid(rowIndices(rowIndices>0));
% % % 
% % % ixInExobase = log10(elDepthMedians)>0;
% % % ixExobase = ixInExobase;
% % % % ax = axes(f, 'Position',[lef(iOGR), bot, wid, heig]);
% % % 
% % % rMax = 300;
% % % 
% % % radialGrid = elDepthRadGrid; % Linear grid from 0 to 300 km
% % % angularAxis = acos(linspace(-1, 1, nElDepthThetaBounds)); % Linear grid from -1 to 1
% % % radialGrid = 0.5 * (radialGrid(2:end)+ radialGrid(1:end-1));
% % % angularAxis = 0.5*(angularAxis(2:end)+angularAxis(1:end-1));
% % % 
% % % 
% % % % Create a meshgrid for the polar coordinates
% % % [theta, rho] = meshgrid(angularAxis, radialGrid);
% % % % Convert polar coordinates to Cartesian coordinates
% % % x = rho .* cos(theta);
% % % y = rho .* sin(theta);
% % % 
% % % 
% % % % Create a 2D plot
% % % data = log10(elDepthMedians(:,:));
% % % hold(ax, "on")
% % % 
% % % pcolor(ax, x, y, data);
% % % xsInExo = x(ixInExobase);
% % % ysInExo = y(ixInExobase);
% % % k = convhull(xsInExo, ysInExo);
% % % plot(ax, xsInExo(k), ysInExo(k), 'k-')
% % % shading interp; % Interpolate the colors
% % % colormap redblue; % Set the colormap
% % % 
% % % % Set the aspect ratio to equal to avoid distortion
% % % axis equal;
% % % ax.XTick = ax.XTick(ax.XTick >=0);
% % % % Set the x-axis limits
% % % xlim([-300, 300]);
% % % ylim([0, 300]);
% % % 
% % % % ax.XTicks = xticks(xticks >= 0);
% % % % Add labels and title
% % % xlabel('r [km]');
% % % ylabel('y');
% % % 
% % % % Add faint grey lines at regular angles with labels
% % % hold on;
% % % angleTicks = linspace(0, pi, 7); % Adjust the number of lines as needed
% % % angleLabels = {'0°', '30°', '60°', '90°', '120°', '150°', '180°'}; % Adjust the labels as needed
% % % for i = 1:length(angleTicks)
% % %     angle = angleTicks(i);
% % %     xline = 300 * cos(angle);
% % %     yline = 300 * sin(angle);
% % %     plot(ax, [-xline, xline], [-yline, yline], 'Color', [0.2, 0.2, 0.2], 'LineStyle', '--');
% % %     text(1.1 * xline - abs(0.05*xline), 1.1* yline, angleLabels{i}, 'Color', [0,0,0]);
% % % end
% % % % Add a colorbar
% % % %     colorbar;
% % % ax.Box = 'off';
% % % ax.Color = 'none';
% % % % ax.YLabel = ''
% % % ax.YAxis.Visible = 'off';
% % % % if iOGR == nRuns    
% % % cb = colorbar('Position', [0.9,bot, 0.02, heig]);
% % % cb.Label.String = 'log_{10}(Median \tau)';
% % % clim(ax, clims)
% % % % end
% % % % caxis(ax{iOGR},[-2,0.5]);
% % % % set(gcf, 'Color', 'none');
% % % % Optional: Customize the plot appearance
% % % % For example, you can set the font size of the labels and title:
% % % % set(gca, 'FontSize', 12);
% % % % t1 = title(ax{iOGR}, caseNameList{iOGR});
% % % % t1.Position(2) = t1.Position(2)+50;
% % % 
% % % % end
% % % % sgtitle(['Median e Depths Relative to ',axisName{iDim},' axis']);
% % % 
% % % 
% % % end