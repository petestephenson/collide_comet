function [f, ax, summary] = make_grouped_cold_e_prod_Plot_v3( ...
    peCasesNo, swCasesNo, plasDensDir, iAngle, r0, dr, nAngleBins)

%MAKE_GROUPED_COLD_E_PROD_PLOT
%
% Plot grouped cold-electron production maps for each outgassing case.
%
% Panels:
%   1-5: 2D slices of cold electron production.
%   6  : Summary plot of shell-averaged cold electron production
%        at r = r0 +/- dr/2, binned by angle to outgassing direction.
%
% Summary plot:
%   - Asymmetric cases are plotted relative to their own qVec.
%   - Symmetric/reference case is plotted relative to each asymmetric
%     outgassing direction: -X, +Y, +Z, -Z.
%
% Usage:
%   [f, ax, summary] = asymPlots.make_grouped_cold_e_prod_Plot( ...
%       peCasesNo, swCasesNo, plasDensDir, 1, 50, 10, 36);

    if nargin < 4 || isempty(iAngle)
        iAngle = 1;
    end

    if nargin < 5 || isempty(r0)
        r0 = 50; % km
    end

    if nargin < 6 || isempty(dr)
        dr = 10; % km
    end

    if nargin < 7 || isempty(nAngleBins)
        nAngleBins = 36;
    end

    rMin = r0 - 0.5 * dr;
    rMax = r0 + 0.5 * dr;

    if rMin < 0
        error('r0 - dr/2 must be >= 0. Current rMin = %.3f km.', rMin);
    end

    % ------------------------------------------------------------
    % Settings
    % ------------------------------------------------------------
    E_cold = 1; % eV

    asymName = {'-x', '+y', '+z', '-z', 'sym'};

    % Directions used for plotting the symmetric case in the summary panel.
    refAxisNames = { ...
        'sym, -x', ...
        'sym, +y', ...
        'sym, +z', ...
        'sym, -z'};

    refAxisDirs = { ...
        [-1, 0, 0], ...
        [ 0, 1, 0], ...
        [ 0, 0, 1], ...
        [ 0, 0,-1]};

    caseColors = lines(5);
    caseMarkers = {'o', 'pentagram', '^', 'v', 's'};

    left = [0.35, 0.65, 0.35, 0.65, 0.07];
    bottom = [0.57, 0.57, 0.10, 0.10, 0.57];

    cbpos = [0.90, 0.10, 0.01, 0.80];

    wid = 0.23;
    heig = 0.35;

    vangleList = { ...
        [1, 0, 0], ...
        [0,-1, 0], ...
        [0, 0, 1]};

    viewNamelist = {'YZ', 'XZ', 'XY'}; %#ok<NASGU>

    slices = {[0], [0], [0]};
    plot_lim = 300;
    clims = [-3, 3];

    inc_gens = [1, 2];

    vangle = vangleList{iAngle};

    % ------------------------------------------------------------
    % Figure and axes
    % ------------------------------------------------------------
    f = figure( ...
        'Position', [0, 0, 1200, 700], ...
        'Color', 'w');

    ax = cell(6, 1);
    cb = cell(5, 1);

    % 2D map axes.
    for iCase = 1:5
        ax{iCase} = axes(f, ...
            'Position', [left(iCase), bottom(iCase), wid, heig]);
    end

    % Bottom-left summary axis.
    summaryAx = 6;

    ax{summaryAx} = axes(f, ...
        'Position', [0.07, 0.10, wid, heig]);

    hold(ax{summaryAx}, 'on');
    box(ax{summaryAx}, 'on');

    set(ax{summaryAx}, ...
        'FontSize', 12, ...
        'LineWidth', 1.1, ...
        'TickDir', 'out', ...
        'YScale', 'log');

    xlabel(ax{summaryAx}, 'Angle to outgassing direction [deg]');
    ylabel(ax{summaryAx}, 'Mean P_e [cm^{-3} s^{-1}]', ...
        'Interpreter', 'tex');

    title(ax{summaryAx}, sprintf('r = %.0f \\pm %.0f km', ...
        r0, 0.5 * dr), ...
        'Interpreter', 'tex');

    xlim(ax{summaryAx}, [0, 180]);
    xticks(ax{summaryAx}, 0:30:180);

    % ------------------------------------------------------------
    % Summary output storage
    % ------------------------------------------------------------
    summary = struct();

    summary.r0 = r0;
    summary.dr = dr;
    summary.rMin = rMin;
    summary.rMax = rMax;
    summary.nAngleBins = nAngleBins;

    summary.angleEdges = linspace(0, 180, nAngleBins + 1);
    summary.angleCentres = 0.5 .* ...
        (summary.angleEdges(1:end-1) + summary.angleEdges(2:end));

    summary.cases = struct([]);

    % ------------------------------------------------------------
    % Load cross sections
    % ------------------------------------------------------------
    swRefNo = swCasesNo{5};
    swRefName = [plasDensDir, 'plas_dens', swRefNo, '.nc'];

    [XS, XS_Eng, E_th] = asymPlots.loadXSections();

    % This reference load is kept in case you need it later.
    [qSphericalRef, qDaysideRef, qUniformRef, qVecRef] = ...
        asymPlots.loadOutgassingProperties(swRefName); %#ok<ASGLU>

    % ------------------------------------------------------------
    % Main case loop
    % ------------------------------------------------------------
    for iCase = 1:5

        fprintf('Processing cold electron production for %s (%d/5)...\n', ...
            asymName{iCase}, iCase);

        peAsymNo = peCasesNo{iCase};
        swAsymNo = swCasesNo{iCase};

        peAsymName = [plasDensDir, 'plas_dens', peAsymNo, '.nc'];
        swAsymName = [plasDensDir, 'plas_dens', swAsymNo, '.nc'];

        [gridsAsym, plas_dens_E_Asym_pe] = ...
            asymPlots.loadRunData(peAsymName);

        [gridsAsym, plas_dens_E_Asym_sw] = ...
            asymPlots.loadRunData(swAsymName);

        plas_dens_E_Asym = plas_dens_E_Asym_sw + plas_dens_E_Asym_pe;

        % --------------------------------------------------------
        % Grid for 2D plotting
        % --------------------------------------------------------
        [X_mesh, Y_mesh, Z_mesh] = meshgrid( ...
            gridsAsym{1}(1:end-1), ...
            gridsAsym{2}(1:end-1), ...
            gridsAsym{3}(1:end-1));

        % --------------------------------------------------------
        % Cell-centre grid for summary calculation
        % --------------------------------------------------------
        x_edges = gridsAsym{1};
        y_edges = gridsAsym{2};
        z_edges = gridsAsym{3};

        dx = diff(x_edges);
        dy = diff(y_edges);
        dz = diff(z_edges);

        x_centres = x_edges(1:end-1) ;%+ 0.5 .* dx;
        y_centres = y_edges(1:end-1) ;%+ 0.5 .* dy;
        z_centres = z_edges(1:end-1) ;%+ 0.5 .* dz;

        [Xc, Yc, Zc] = meshgrid(x_centres, y_centres, z_centres);
        [dX, dY, dZ] = meshgrid(dx, dy, dz);

        cellVol = dX .* dY .* dZ; % km^3

        % --------------------------------------------------------
        % Load outgassing properties
        % --------------------------------------------------------
        [qSpherical, qDayside, qUniform, qVec] = ...
            asymPlots.loadOutgassingProperties(swAsymName);

        me = asymPlots.get_electron_mass(peAsymName); %#ok<NASGU>

        if iCase == 5
            qDayside = 0.0;
            qUniform = qSpherical;
        end

        % --------------------------------------------------------
        % Calculate cold-electron production
        % --------------------------------------------------------
        coldElecProd = asymPlots.calculateColdeProduction( ...
            plas_dens_E_Asym, ...
            gridsAsym, ...
            E_cold, ...
            qDayside, ...
            qUniform, ...
            qVec, ...
            XS_Eng, ...
            XS, ...
            E_th, ...
            inc_gens);

        % --------------------------------------------------------
        % 2D cold-electron production panel
        % --------------------------------------------------------
        axes(ax{iCase}); %#ok<LAXES>

        [f, ax{iCase}, cb{iCase}] = asymPlots.makecoldEProdPlot( ...
            coldElecProd, ...
            X_mesh, Y_mesh, Z_mesh, ...
            slices, ...
            plot_lim, ...
            1, ...
            vangle, ...
            clims, ...
            f, ...
            ax{iCase});

        title(ax{iCase}, asymName{iCase}, ...
            'Interpreter', 'none');

        if iCase == 1 || iCase == 2
            xlabel(ax{iCase}, '');
        end

        delete(cb{iCase});

        % --------------------------------------------------------
        % Summary profile
        % --------------------------------------------------------
        if iCase < 5

            % ----------------------------------------------------
            % Asymmetric cases:
            % plot relative to this case's own qVec.
            % ----------------------------------------------------
            qHat = qVec(:);

            if norm(qHat) == 0
                warning('Skipping summary profile for %s: qVec has zero magnitude.', ...
                    asymName{iCase});

                summary.cases(iCase).caseName = asymName{iCase};
                summary.cases(iCase).qHat = [];
                summary.cases(iCase).prodProfile = [];
                summary.cases(iCase).nCellsProfile = [];

            else

                qHat = qHat ./ norm(qHat);

                [prodProfile, nCellsProfile] = coldProdShellAngleProfile( ...
                    coldElecProd, ...
                    Xc, Yc, Zc, ...
                    cellVol, ...
                    qHat, ...
                    rMin, rMax, ...
                    summary.angleEdges);

                goodProfile = isfinite(prodProfile) & prodProfile > 0;

                xPlot = summary.angleCentres(goodProfile);
                yPlot = prodProfile(goodProfile);

                plot(ax{summaryAx}, ...
                    xPlot, ...
                    yPlot, ...
                    'LineStyle', '-', ...
                    'Marker', caseMarkers{iCase}, ...
                    'MarkerIndices', 1:3:numel(xPlot), ...
                    'Color', caseColors(iCase, :), ...
                    'MarkerFaceColor', caseColors(iCase, :), ...
                    'LineWidth', 1.2, ...
                    'MarkerSize', 4, ...
                    'DisplayName', asymName{iCase});

                summary.cases(iCase).caseName = asymName{iCase};
                summary.cases(iCase).qHat = qHat;
                summary.cases(iCase).prodProfile = prodProfile;
                summary.cases(iCase).nCellsProfile = nCellsProfile;

            end

        else

            % ----------------------------------------------------
            % Symmetric case:
            % plot relative to each asymmetric outgassing direction.
            % ----------------------------------------------------
            summary.cases(iCase).caseName = asymName{iCase};
            summary.cases(iCase).refProfiles = struct([]);

            for iRef = 1:numel(refAxisDirs)

                qHatRef = refAxisDirs{iRef}(:);
                qHatRef = qHatRef ./ norm(qHatRef);

                [prodProfileRef, nCellsProfileRef] = coldProdShellAngleProfile( ...
                    coldElecProd, ...
                    Xc, Yc, Zc, ...
                    cellVol, ...
                    qHatRef, ...
                    rMin, rMax, ...
                    summary.angleEdges);

                goodProfileRef = isfinite(prodProfileRef) & ...
                                 prodProfileRef > 0;

                xPlot = summary.angleCentres(goodProfileRef);
                yPlot = prodProfileRef(goodProfileRef);

                plot(ax{summaryAx}, ...
                    xPlot, ...
                    yPlot, ...
                    'LineStyle', '--', ...
                    'Marker', caseMarkers{iRef}, ...
                    'MarkerIndices', 1:3:numel(xPlot), ...
                    'Color', caseColors(iRef, :), ...
                    'MarkerFaceColor', 'none', ...
                    'LineWidth', 1.2, ...
                    'MarkerSize', 4, ...
                    'DisplayName', refAxisNames{iRef});

                summary.cases(iCase).refProfiles(iRef).refName = ...
                    refAxisNames{iRef};

                summary.cases(iCase).refProfiles(iRef).qHat = qHatRef;
                summary.cases(iCase).refProfiles(iRef).prodProfile = ...
                    prodProfileRef;

                summary.cases(iCase).refProfiles(iRef).nCellsProfile = ...
                    nCellsProfileRef;

            end

        end

    end

    % ------------------------------------------------------------
    % Shared colourbar for 2D panels
    % ------------------------------------------------------------
    cbPlot = colorbar(ax{1}, 'Position', cbpos);
    cbPlot.Label.String = 'log_{10}(P_e [cm^{-3}s^{-1}])';
    cbPlot.Label.Interpreter = 'tex';

    % ------------------------------------------------------------
    % Summary legend
    % ------------------------------------------------------------
    legend(ax{summaryAx}, ...
        'Location', 'best', ...
        'Interpreter', 'none', ...
        'NumColumns', 2);

end


% ========================================================================
function [prodProfile, nCellsProfile] = coldProdShellAngleProfile( ...
    coldElecProd, Xc, Yc, Zc, cellVol, qHat, rMin, rMax, angleEdges)
%COLDPRODSHELLANGLEPROFILE
%
% Volume-weighted mean cold-electron production in a radial shell,
% binned by angle to qHat.
%
% coldElecProd units are assumed to be cm^-3 s^-1.
% prodProfile has the same units.

    qHat = qHat(:);

    if norm(qHat) == 0
        error('qHat must have non-zero magnitude.');
    end

    qHat = qHat ./ norm(qHat);

    % ------------------------------------------------------------
    % Match coldElecProd orientation to meshgrid output if needed.
    % meshgrid(x_centres, y_centres, z_centres) gives size [ny, nx, nz].
    % ------------------------------------------------------------
    if ~isequal(size(coldElecProd), size(Xc))

        if isequal(size(permute(coldElecProd, [2, 1, 3])), size(Xc))

            coldElecProd = permute(coldElecProd, [2, 1, 3]);

        else

            error( ...
                ['coldElecProd size %s does not match grid size %s, ', ...
                 'and permute(coldElecProd,[2,1,3]) also does not match.'], ...
                mat2str(size(coldElecProd)), ...
                mat2str(size(Xc)));

        end
    end

    if ~isequal(size(cellVol), size(Xc))
        error('cellVol size %s does not match grid size %s.', ...
            mat2str(size(cellVol)), mat2str(size(Xc)));
    end

    % ------------------------------------------------------------
    % Radius and angular coordinate
    % ------------------------------------------------------------
    r = sqrt(Xc.^2 + Yc.^2 + Zc.^2);

    rSafe = r;
    rSafe(rSafe == 0) = NaN;

    xHat = Xc ./ rSafe;
    yHat = Yc ./ rSafe;
    zHat = Zc ./ rSafe;

    cosAng = ...
        xHat .* qHat(1) + ...
        yHat .* qHat(2) + ...
        zHat .* qHat(3);

    cosAng(cosAng > 1) = 1;
    cosAng(cosAng < -1) = -1;

    angleDeg = acosd(cosAng);

    % ------------------------------------------------------------
    % Select shell and valid values
    % ------------------------------------------------------------
    shellMask = r >= rMin & r < rMax;

    good = shellMask & ...
           isfinite(angleDeg) & ...
           isfinite(coldElecProd) & ...
           isfinite(cellVol) & ...
           coldElecProd > 0 & ...
           cellVol > 0;

    angleVals = angleDeg(good);
    prodVals = coldElecProd(good);
    volVals = cellVol(good);

    nAngleBins = numel(angleEdges) - 1;

    prodProfile = NaN(nAngleBins, 1);
    nCellsProfile = zeros(nAngleBins, 1);

    if isempty(angleVals)
        warning('No valid cells found in requested radial shell.');
        return
    end

    [~, ~, binIdx] = histcounts(angleVals, angleEdges);

    inBin = binIdx > 0;

    binIdx = binIdx(inBin);
    prodVals = prodVals(inBin);
    volVals = volVals(inBin);

    if isempty(binIdx)
        warning('No cells fell inside the angular bins.');
        return
    end

    % ------------------------------------------------------------
    % Volume-weighted mean production per angle bin
    % ------------------------------------------------------------
    weightedProdSum = accumarray( ...
        binIdx, ...
        prodVals .* volVals, ...
        [nAngleBins, 1], ...
        @sum, ...
        NaN);

    volSum = accumarray( ...
        binIdx, ...
        volVals, ...
        [nAngleBins, 1], ...
        @sum, ...
        NaN);

    nCellsProfile = accumarray( ...
        binIdx, ...
        ones(size(binIdx)), ...
        [nAngleBins, 1], ...
        @sum, ...
        0);

    prodProfile = weightedProdSum ./ volSum;

    bad = ~isfinite(prodProfile) | prodProfile <= 0 | ...
          ~isfinite(volSum) | volSum <= 0 | ...
          nCellsProfile <= 0;

    prodProfile(bad) = NaN;

end