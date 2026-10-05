function [f, ax0, cb, rMids, temp_vs_r, fitresult, ...
          temp_var_vs_r, temp_var_unweighted, temp_std_vs_r] = ...
    makeTemperature_vs_r_Plot_v2(temperature, plas_dens, X_mesh, Y_mesh, Z_mesh, varargin)
% makeTemperature_vs_r_Plot
% Plots T_e vs radius with a density-weighted mean curve, an unweighted mean,
% and a rational (rat22) fit. Also returns binned variance (weighted & unweighted).
%
% Outputs (first 6 unchanged for backward-compat):
%   1 f, 2 ax0, 3 cb, 4 rMids, 5 temp_vs_r (weighted mean), 6 fitresult,
%   7 temp_var_vs_r (weighted population variance),
%   8 temp_var_unweighted (unweighted population variance),
%   9 temp_std_vs_r (sqrt of weighted variance)
%
% Notes:
% - Weighted stats use plas_dens as weights (population variance: sum(w*(x-mu)^2)/sum(w)).
% - Unweighted variance is population variance (var(x,1)).

    % ---------- Figure/Axes ----------
    if numel(varargin) >= 1 && ~isempty(varargin{1}), f = varargin{1}; else, f = figure(); end
    if numel(varargin) >= 2 && ~isempty(varargin{2}), ax0 = varargin{2}; else, ax0 = axes(f); end
    hold(ax0, 'on');

    % ---------- Data prep ----------
    % Match your original orientation change
    temperature = permute(temperature, [2,1,3]);
    plas_dens   = permute(plas_dens,   [2,1,3]);

    % Radius grid
    r_Mesh = sqrt(X_mesh.^2 + Y_mesh.^2 + Z_mesh.^2);

    % Flatten everything for fast vectorized binning
    r   = r_Mesh(:);
    T   = temperature(:);
    W   = plas_dens(:);
    Cx  = X_mesh(:);                 % for scatter coloring

    % Keep only finite, non-negative weights
    good = isfinite(r) & isfinite(T) & isfinite(W) & (W >= 0);
    r = r(good); T = T(good); W = W(good); Cx = Cx(good);

    % ---------- Scatter ----------
    sc = scatter(ax0, r, T, 5, Cx, 'filled'); %#ok<NASGU>
    xlabel(ax0, 'r [km]');
    ylabel(ax0, 'T_e [eV]');
    cb = colorbar(ax0); cb.Label.String = 'X [km]';

    % ---------- Bin setup ----------
    nBins = 100;
    rMin = min(r);
    rMax = max(r);
    edges = linspace(rMin, rMax, nBins+1);
    rMids = 0.5*(edges(1:end-1) + edges(2:end));

    % Bin index for each point
    b = discretize(r, edges);

    % ---------- Weighted mean ----------
    % mean_w = sum(w*x)/sum(w)
    num_w = accumarray(b, W.*T, [nBins,1], @sum, NaN);
    den_w = accumarray(b, W,    [nBins,1], @sum, NaN);
    temp_vs_r = num_w ./ den_w;                 % weighted mean

    % ---------- Unweighted mean ----------
    temp_vs_r_unweighted = accumarray(b, T, [nBins,1], @(x) mean(x,'omitnan'), NaN);

    % ---------- Weighted variance (population) ----------
    % var_w = sum(w*(x-mu)^2)/sum(w)
    mu_per_point = temp_vs_r(b);                        % map bin mean to each point
    sq_term = W.*(T - mu_per_point).^2;
    num_var_w = accumarray(b, sq_term, [nBins,1], @sum, NaN);
    temp_var_vs_r = num_var_w ./ den_w;                 % weighted population variance
    temp_std_vs_r = sqrt(temp_var_vs_r);

    % ---------- Unweighted variance (population) ----------
    % var(x,1) => normalize by N (population)
    temp_var_unweighted = accumarray(b, T, [nBins,1], @(x) var(x,1,'omitnan'), NaN);
    % ---------- Plot binned curves ----------
    % Make all vectors columns to avoid implicit expansion creating NxN masks
    rMids                 = rMids(:);
    temp_vs_r             = temp_vs_r(:);
    temp_vs_r_unweighted  = temp_vs_r_unweighted(:);
    temp_var_vs_r         = temp_var_vs_r(:);
    temp_std_vs_r         = temp_std_vs_r(:);
    
    plot(ax0, rMids, temp_vs_r, 'k-',  'LineWidth', 1.6, 'DisplayName', '⟨T_e⟩_w (density-weighted)');
    plot(ax0, rMids, temp_vs_r_unweighted, 'k--', 'LineWidth', 1.4, 'DisplayName', '⟨T_e⟩ (unweighted)');
    
    % Build a 1×N logical mask safely (column vector)
    ok = isfinite(rMids) & isfinite(temp_vs_r) & isfinite(temp_std_vs_r);
    
    if any(ok)
        errorbar(ax0, rMids(ok), temp_vs_r(ok), temp_std_vs_r(ok), ...
                 'LineStyle','none','Color',[0 0 0 0.25],'CapSize',0, ...
                 'DisplayName','1\sigma (weighted)');
    end

    % ---------- Plot binned curves ----------
%     plot(ax0, rMids, temp_vs_r, 'k-',  'LineWidth', 1.6, 'DisplayName', '⟨T_e⟩_w (density-weighted)');
%     plot(ax0, rMids, temp_vs_r_unweighted, 'k--', 'LineWidth', 1.4, 'DisplayName', '⟨T_e⟩ (unweighted)');

    % Optional: 1σ band for weighted mean (light error bars)
    % Only draw where finite:
%     ok = isfinite(rMids) & isfinite(temp_vs_r) & isfinite(temp_std_vs_r);
%     if any(ok)
%         % Thin errorbars to avoid clutter
%         errorbar(ax0, rMids(ok), temp_vs_r(ok), temp_std_vs_r(ok), ...
%                  'LineStyle','none','Color',[0 0 0 0.25],'CapSize',0, ...
%                  'DisplayName','1\sigma (weighted)');
%     end

    % ---------- Fit (rat22) on finite weighted means ----------
    fit_ok = isfinite(rMids) & isfinite(temp_vs_r);
    xData = rMids(fit_ok);
    yData = temp_vs_r(fit_ok);

    if numel(xData) >= 5 && numel(yData) >= 5
        ft2 = fittype('rat22');
        opts2 = fitoptions('Method','NonlinearLeastSquares');
        opts2.Display = 'Off';
        opts2.StartPoint = [0.117417650855806, 0.296675873218327, 0.318778301925882, 0.424166759713807, 0.507858284661118];
        % you had a lower bound on the 3rd coeff:
        opts2.Lower = [-Inf, -Inf, 1e3, -Inf, -Inf];

        [fitresult, ~] = fit(xData(:), yData(:), ft2, opts2);

        testRGrid = linspace(rMin, rMax, 200);
        plot(ax0, testRGrid, fitresult(testRGrid), 'b-', 'LineWidth', 1.6, 'DisplayName', 'Rational fit (rat22)');
    else
        fitresult = [];
    end

    legend(ax0, 'Location','northwest');
    box(ax0,'on');

end
