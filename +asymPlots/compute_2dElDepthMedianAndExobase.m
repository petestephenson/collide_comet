function out = compute_2dElDepthMedianAndExobase( ...
    elDepthRadEdges, elDepthGrid, elDepthDist, iDim, opts)
%COMPUTE_2DELDEPTHMEDIANANDEXOBASE Compute median depths and (smoothed) exobase hull with no plotting.
%
% Adds theta-padding so you have explicit 0 and pi columns (useful for later pcolor/surface with no sliver).
% Padding strategy: duplicate the first/last theta-bin values at the boundaries.

    if nargin < 5 || isempty(opts), opts = struct(); end
    if ~isfield(opts, "iTau"),           opts.iTau = 3; end
    if ~isfield(opts, "medianLevel"),    opts.medianLevel = 0.5; end
    if ~isfield(opts, "exoThreshLog10"), opts.exoThreshLog10 = 0; end

    % Smoothed hull options
    if ~isfield(opts, "smoothHull"),     opts.smoothHull = true; end
    if ~isfield(opts, "smoothMethod"),   opts.smoothMethod = "sgolay"; end % "movmean" | "sgolay"
    if ~isfield(opts, "smoothSpan"),     opts.smoothSpan = 14; end         % odd recommended
    if ~isfield(opts, "sgolayOrder"),    opts.sgolayOrder = 2; end
    if ~isfield(opts, "hullResampleN"),  opts.hullResampleN = 50; end

    % Theta padding options
    if ~isfield(opts, "padTheta"),       opts.padTheta = true; end
    if ~isfield(opts, "thetaPadMode"),   opts.thetaPadMode = "duplicate"; end % "duplicate" (recommended)

    % Envelope options (NEW)
    if ~isfield(opts, "smoothEnvelope"), opts.smoothEnvelope = true; end
    if ~isfield(opts, "envelopeMethod"), opts.envelopeMethod = "movmean"; end % "movmean" | "sgolay"
    if ~isfield(opts, "envelopeSpan"),   opts.envelopeSpan = 9; end           % odd recommended
    if ~isfield(opts, "envelopeOrder"),  opts.envelopeOrder = 2; end          % sgolay order if used
    if ~isfield(opts, "envelopeFill"),   opts.envelopeFill = "interp"; end    % "interp" | "none"


    % ---- Median depth grid on theta MIDPOINTS ----
    [medianDepth, rMids, thetaMids] = compute_medianDepth_grid( ...
        elDepthRadEdges, elDepthGrid, elDepthDist, iDim, opts.iTau, opts.medianLevel);

    log10MedianDepth = log10(medianDepth);
    inExo = log10MedianDepth > opts.exoThreshLog10;

    % ---- Midpoint mesh in (x,y) ----
    [theta, rho] = meshgrid(thetaMids, rMids);
    x = rho .* cos(theta);
    y = rho .* sin(theta);

    % ---- Pad in theta with explicit 0 and pi columns (and duplicate values) ----
    if opts.padTheta
        [thetaMids_pad, medianDepth_pad, log10MedianDepth_pad, inExo_pad, x_pad, y_pad] = ...
            pad_theta_midpoint_fields(rMids, thetaMids, medianDepth, log10MedianDepth, inExo);
    else
        thetaMids_pad        = thetaMids;
        medianDepth_pad      = medianDepth;
        log10MedianDepth_pad = log10MedianDepth;
        inExo_pad            = inExo;
        x_pad                = x;
        y_pad                = y;
    end

    % ---- Hull uses padded fields (so it can reach the ends if needed) ----
    if opts.smoothHull
        exo = compute_smoothed_exobase_hull_xy(x_pad, y_pad, inExo_pad, opts);
    else
        exo = compute_raw_exobase_hull_xy(x_pad, y_pad, inExo_pad);
    end

    % (2) Envelope r(theta): outermost r where inExo is true per theta bin
    exo.envelope = compute_exobase_envelope_rt(rMids, thetaMids_pad, inExo_pad, opts);

    % ---- Pack outputs ----
    out = struct();
    out.medianDepth         = medianDepth_pad;
    out.log10MedianDepth    = log10MedianDepth_pad;
    out.inExo               = inExo_pad;
    out.rMids               = rMids;
    out.thetaMids           = thetaMids_pad;
    out.x                   = x_pad;
    out.y                   = y_pad;
    out.exo                 = exo;

    % Keep unpadded versions too
    out.medianDepth_unpad      = medianDepth;
    out.log10MedianDepth_unpad = log10MedianDepth;
    out.inExo_unpad            = inExo;
    out.thetaMids_unpad        = thetaMids;
end


function [thetaMids_pad, median_pad, log10median_pad, inExo_pad, x_pad, y_pad] = ...
    pad_theta_midpoint_fields(rMids, thetaMids, median, log10median, inExo)
%PAD_THETA_MIDPOINT_FIELDS
% Create explicit theta columns at 0 and pi by duplicating the first/last theta-bin values.
%
% Inputs:
%   rMids       (nR x 1)
%   thetaMids   (1 x nT) in radians
%   median      (nR x nT)
%   log10median (nR x nT)
%   inExo       (nR x nT)
%
% Outputs:
%   thetaMids_pad (1 x (nT+2)) = [0, thetaMids, pi]
%   *_pad         (nR x (nT+2)) with duplicated end columns
%   x_pad, y_pad  (nR x (nT+2))

    nT = numel(thetaMids);

    thetaMids_pad = [pi, reshape(thetaMids,1,[]), 0];

    % Duplicate first/last columns so boundary columns have sane values
    median_pad = [median(:,1), median, median(:,nT)];
    log10median_pad = [log10median(:,1), log10median, log10median(:,nT)];
    inExo_pad = [inExo(:,1), inExo, inExo(:,nT)];

    [thPad, rPad] = meshgrid(thetaMids_pad, rMids);
    x_pad = rPad .* cos(thPad);
    y_pad = rPad .* sin(thPad);
end


function [medianDepth, rMids, thetaMids] = compute_medianDepth_grid( ...
    elDepthRadEdges, elDepthGrid, elDepthDist, iDim, iTau, medianLevel)

    rMids = 0.5 * (elDepthRadEdges(1:end-1) + elDepthRadEdges(2:end));

    nTheta = size(elDepthDist, 4);
    thetaBounds = acos(linspace(-1, 1, nTheta + 1)); % 0..pi
    thetaMids = 0.5 * (thetaBounds(1:end-1) + thetaBounds(2:end)); % midpoints (radians)

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


function exo = compute_raw_exobase_hull_xy(x, y, inExo)
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

    k = convhull(xs, ys);
    exo = struct();
    exo.x_raw = xs(k);
    exo.y_raw = ys(k);
    exo.k = k;
end


function exo = compute_smoothed_exobase_hull_xy(x, y, inExo, opts)
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

    k = convhull(xs, ys);
    xr = xs(k);
    yr = ys(k);

    ds = hypot(diff(xr), diff(yr));
    s = [0; cumsum(ds)];
    L = s(end);
    if ~isfinite(L) || L == 0
        exo = [];
        return;
    end

    N = max(20, round(opts.hullResampleN));
    sU = linspace(0, L, N+1)';
    xU = interp1(s, xr, sU, 'linear');
    yU = interp1(s, yr, sU, 'linear');

    xU = xU(1:end-1);
    yU = yU(1:end-1);

    span = max(3, round(opts.smoothSpan));
    if mod(span,2) == 0, span = span + 1; end
    halfw = floor(span/2);

    xW = [xU(end-halfw+1:end); xU; xU(1:halfw)];
    yW = [yU(end-halfw+1:end); yU; yU(1:halfw)];

    method = string(opts.smoothMethod);
    switch lower(method)
        case "sgolay"
            try
                ord = min(round(opts.sgolayOrder), span-2);
                xWs = sgolayfilt(xW, ord, span);
                yWs = sgolayfilt(yW, ord, span);
            catch
                xWs = movmean(xW, span);
                yWs = movmean(yW, span);
            end
        otherwise
            xWs = movmean(xW, span);
            yWs = movmean(yW, span);
    end

    xS = xWs(halfw+1 : halfw+numel(xU));
    yS = yWs(halfw+1 : halfw+numel(yU));

    xS(end+1) = xS(1);
    yS(end+1) = yS(1);

    exo = struct();
    exo.k = k;
    exo.x_raw = xr;
    exo.y_raw = yr;
    exo.x_resamp = [xU; xU(1)];
    exo.y_resamp = [yU; yU(1)];
    exo.x_smooth = xS;
    exo.y_smooth = yS;
end


function exoRT = compute_exobase_envelope_rt(rMids, thetaMids, inExo, opts)
% Compute r_exo(theta) = max r where inExo is true, then smooth vs theta.

    if nargin < 4 || isempty(opts), opts = struct(); end
    if ~isfield(opts,"smoothTheta"), opts.smoothTheta = true; end
    if ~isfield(opts,"smoothSpan"),  opts.smoothSpan  = 1; end

    nT = numel(thetaMids);
    r_exo = NaN(1,nT);

    for it = 1:nT
        idx = find(inExo(:,it), 1, 'last');   % outermost r
        if ~isempty(idx)
            r_exo(it) = rMids(idx);
        end
    end

    % Optional smoothing in theta (not periodic for 0..pi)
    if opts.smoothTheta
        span = max(3, round(opts.smoothSpan));
        if mod(span,2)==0, span = span+1; end
        r_exo_s = movmean(r_exo, span, 'omitnan');  % simple + robust
    else
        r_exo_s = r_exo;
    end

    % Convert to xy for plotting later
    x = r_exo_s .* cos(thetaMids);
    y = r_exo_s .* sin(thetaMids);

    exoRT = struct('theta',thetaMids,'r',r_exo,'r_smooth',r_exo_s,'x',x,'y',y);
end

