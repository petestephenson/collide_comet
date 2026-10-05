function [f,ax0, cb, rMids, temp_vs_r, fitresult] = makeTemperature_vs_r_Plot(temperature,plas_dens, X_mesh,Y_mesh,Z_mesh, varargin)
    if length(varargin) >= 1 && ~isempty(varargin{1})
        f = varargin{1};
    else
        f = figure();
    end
    
    if length(varargin) >= 2 && ~isempty(varargin{2})
        ax0 = varargin{2};
    else
        ax0 = axes(f);
    end

    hold(ax0, 'on');

    plot_str = 'X [km]';
    
%     clims = [0,40];
    iLogArray = 0;
    temperature = permute(temperature,[2,1,3]);
    plas_dens = permute(plas_dens,[2,1,3]);
    r_Mesh = sqrt(X_mesh.^2 + Y_mesh.^2 + Z_mesh.^2);
    % # Print shapes of r_Mesh, temperature, X_mesh
    disp(size(r_Mesh));
    disp(size(temperature));
    disp(size(X_mesh));
    scatter(ax0, reshape(r_Mesh, [], 1), reshape(temperature, [], 1), 5, reshape(X_mesh, [], 1));
    % scatter(ax0, r_Mesh,temperature, 5,X_mesh);
    % Label x axis as r [km] and y axis as T_e [eV]
    xlabel('r [km]');
    ylabel('T_e [eV]');
    cb = colorbar;
%     [f,ax0,cb] = asymPlots.makePlot(temperature,X_mesh,Y_mesh,Z_mesh, slices, plot_lim, iLogArray, vangle, f, ax0);
    cb.Label.String = plot_str;

    % Bin the data using r_Mesh and calculate average temperature in each bin
    temp_vs_r = zeros(1, 100);
    temp_std_vs_r = zeros(1,100);
    temp_vs_r_unweighted = zeros(1, 100);
    rBins = linspace(0, max(r_Mesh, [], "all"), 101);
    rMids = rBins(1:100) + diff(rBins)/2;
    for i = 1:100
        rBin = rBins(i);
        rBinNext = rBins(i+1);
        rBinIndices = find(r_Mesh >= rBin & r_Mesh < rBinNext);
        
        temp_vs_r(i) = sum(temperature(rBinIndices).*plas_dens(rBinIndices))/sum(plas_dens(rBinIndices));
        temp_vs_r_unweighted(i) = mean(temperature(rBinIndices));
    end
    plot(ax0, rMids, temp_vs_r, 'LineWidth', 1.5, 'Color', 'k', 'DisplayName', 'Mean Temperature vs r');
    plot(ax0, rMids, temp_vs_r_unweighted, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '--');
    % PLot a horizontal line at y=20 eV

    [xData, yData] = prepareCurveData( rMids, temp_vs_r );

% Set up fittype and options.
    ft = fittype( 'exp2' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Normalize = 'on';
    opts.StartPoint = [23.6066960946882 -0.0717843407161476 -0.00442084004287311 -5.23126230352778];
    % Set lower bounds for the parameters
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    % plot(ax0, fitresult, xData, yData, 'predobs', 0.99);


    testRGrid = linspace(0, 500, 100);
    testTemp = fitresult(testRGrid);
    % plot(ax0, testRGrid, testTemp, 'LineWidth', 1.5, 'Color', 'r');

    ft2 = fittype( 'rat22' );
    opts2 = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts2.Display = 'Off';
    opts2.StartPoint = [0.117417650855806 0.296675873218327 0.318778301925882 0.424166759713807 0.507858284661118];
    opts2.Lower = [-Inf, -Inf, 1e3, -Inf, -Inf];
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft2, opts2 );
    testTemp2 = fitresult(testRGrid);
    plot(ax0, testRGrid, testTemp2, 'LineWidth', 1.5, 'Color', 'b', 'DisplayName', 'Rational Fit');
%     yline(ax0, 20, ':', 'Color', 'k', 'LineWidth', 1.5);
    % Add legend to the plot
    legend(ax0, 'Location', 'northwest');
%     %  Make a custom fit for rMids vs temp_vs_r with a constant term and two exponential terms
%     %  f(x) = a + b*exp(c*x) + d*exp(e*x)
%     %  where a, b, c, d, e are the coefficients
%     %  x is rMids
%     %  y is temp_vs_r
%     %  Coefficients are initialized to 20, 20, -0.01, 20, -0.01

%     %  Fit the data
%     opts = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [ 23.6, -0.07, -0.0044, -3.5],...
%       'Lower', [-Inf, -Inf, -Inf, -Inf], 'Upper', [Inf, Inf, Inf, Inf], 'Algorithm', 'Trust-Region',...
%        'DiffMinChange', 1e-8, 'DiffMaxChange', 0.1, 'TolFun', 1e-8, 'TolX', 1e-8);
% % fitresult = fit(x', y', ft, opts);

%     ft = fittype('b*exp(c*x) + d*exp(e*x)', 'independent', 'x', 'coefficients', {'b', 'c','d','e'});
%     fitresult = fit(rMids', temp_vs_r', ft, opts);
%     disp(fitresult);
    % fitCoeffs = fit(rMids', temp_vs_r', 'b*exp(c*x) + d*exp(e*x)', 'StartPoint', [ 20, -0.1, -1, -3.566], ...
    %     'Method', 'NonlinearLeastSquares', 'Algorithm', 'Trust-Region', 'DiffMinChange', 1e-8, 'MaxIter', 1000);
    % Print the coefficients
    % disp(fitCoeffs);
    % plot(fitresult, rMids, temp_vs_r);
    %  Calculate the fit values
    % fitValues = fitCoeffs(rMids);

    %  Plot the fit
    % plot(ax0, rMids, fitValues, 'LineWidth', 1.5, 'Color', 'r');



%     clim(ax0, clims);
end