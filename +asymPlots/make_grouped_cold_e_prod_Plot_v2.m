function [f, ax, summary] = make_grouped_cold_e_prod_Plot_v2( ...
        peCasesNo, swCasesNo, plasDensDir, iAngle, r0, dr, nAngleBins)

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
                error('r0 - dr/2 must be >= 0.');
        end
        E_cold = 1;% eV - Make default value and allow as an input variable

        f = figure(Position=[0,0,1000, 600]);
        left = [0.35, 0.65, 0.35, 0.65, 0.07];
        bottom = [0.57, 0.57, 0.1, 0.1,0.57];
        asymName = {'-x', '+y', '+z', '-z', 'sym'};
        cbpos = [0.9,0.1,0.01,0.8];
        wid = 0.23;
        heig = 0.35;
        ax = cell(6, 1);
        cb = cell(5, 1);

        summaryAx = 6;
        ax{summaryAx} = axes(f, 'Position', [0.07, 0.1, wid, heig]);
        hold(ax{summaryAx}, 'on');
        box(ax{summaryAx}, 'on');

        set(ax{summaryAx}, ...
        'FontSize', 12, ...
        'LineWidth', 1.1, ...
        'TickDir', 'out', ...
        'YScale', 'log');

        xlabel(ax{summaryAx}, 'Angle to max outgassing [deg]');
        ylabel(ax{summaryAx}, 'Mean P_e [cm^{-3} s^{-1}]');
        title(ax{summaryAx}, sprintf('r = %.0f \\pm %.0f km', r0, 0.5 * dr));

        xlim(ax{summaryAx}, [0, 180]);
        xticks(ax{summaryAx}, 0:30:180);

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

        vangleList = {[1,0,0],[0,-1,0],[0,0,1]};
        viewNamelist = { 'YZ','XZ','XY'};
        slices = {[0],[0],[0]};
        plot_lim = 300;
        clims = [-3,3];
        inc_gens = [1,2];
        vangle = vangleList{iAngle};
        swRefNo = swCasesNo{5};
        swRefName = [plasDensDir, 'plas_dens', swRefNo, '.nc'];
        [XS, XS_Eng, E_th] = asymPlots.loadXSections();
        [qSpherical, qDayside, qUniform, qVec] = asymPlots.loadOutgassingProperties(swRefName);

        for iCase = [1,2,3,4,5]
                ax{iCase} = axes(f, 'Position',[left(iCase), bottom(iCase), wid, heig]);



                peAsymNo = peCasesNo{iCase};
                swAsymNo = swCasesNo{iCase};
                peAsymName = [plasDensDir, 'plas_dens', peAsymNo, '.nc'];
                swAsymName = [plasDensDir, 'plas_dens', swAsymNo, '.nc'];
                [gridsAsym, plas_dens_E_Asym_pe] = asymPlots.loadRunData(peAsymName);
                [gridsAsym, plas_dens_E_Asym_sw] = asymPlots.loadRunData(swAsymName);
                plas_dens_E_Asym = plas_dens_E_Asym_sw + plas_dens_E_Asym_pe;
                [X_mesh, Y_mesh, Z_mesh] = meshgrid(gridsAsym{1}(1:end-1), gridsAsym{2}(1:end-1), gridsAsym{3}(1:end-1));
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


                [qSpherical, qDayside, qUniform, qVec] = asymPlots.loadOutgassingProperties(swAsymName);
                me = asymPlots.get_electron_mass(peAsymName);
                
                if iCase == 5
                        qDayside = 0.0;
                        qUniform = qSpherical;
                end
                % ---- Calculate EII Freq ------
                %         eImpactIoniAsym = asymPlots.calcEIoniFreq(plas_dens_E_Asym, gridsAsym{4}, XS_Eng, XS, inc_gens);
                coldElecProd = asymPlots.calculateColdeProduction(plas_dens_E_Asym, gridsAsym, E_cold, qDayside, qUniform,qVec,XS_Eng, XS,E_th, inc_gens);
                
                % ------------------------------------------------------------
                % 1D summary: cold electron production in shell r0 +/- dr/2
                % ------------------------------------------------------------
                qHat = qVec(:);

                if norm(qHat) == 0 || iCase == 5
                % For the symmetric case, use the current viewing direction
                % or choose a fixed reference direction.
                        qHat = vangle(:);
                end

                qHat = qHat ./ norm(qHat);

                [prodProfile, nCellsProfile] = coldProdShellAngleProfile( ...
                coldElecProd, ...
                Xc, Yc, Zc, ...
                cellVol, ...
                qHat, ...
                rMin, rMax, ...
                summary.angleEdges);

                goodProfile = isfinite(prodProfile) & prodProfile > 0;

                caseColor = lines(5);
                plot(ax{summaryAx}, ...
                summary.angleCentres(goodProfile), ...
                prodProfile(goodProfile), ...
                'LineStyle', '-', ...
                'Marker', 'o', ...
                'MarkerIndices', 1:3:nnz(goodProfile), ...
                'Color', caseColor(iCase, :), ...
                'MarkerFaceColor', caseColor(iCase, :), ...
                'LineWidth', 1.2, ...
                'MarkerSize', 4, ...
                'DisplayName', asymName{iCase});

                summary.cases(iCase).caseName = asymName{iCase};
                summary.cases(iCase).qHat = qHat;
                summary.cases(iCase).prodProfile = prodProfile;
                summary.cases(iCase).nCellsProfile = nCellsProfile;
                
                %         eiiProdAsym = asymPlots.calculateEIIProduction(plas_dens_E_Asym, gridsAsym, qDayside, qUniform,qVec,XS_Eng, XS,inc_gens);
                %         [f,ax{iCase}, cb{iCase}] = asymPlots.makeTemperaturePlot(temperatureSW_Asym,X_mesh,Y_mesh,Z_mesh, slices, plot_lim,  vangle, f, ax{iCase});
                [f,ax{iCase}, cb{iCase}] = asymPlots.makecoldEProdPlot(coldElecProd,X_mesh,Y_mesh,Z_mesh, slices,...
                        plot_lim, 1,  vangle, clims, f, ax{iCase});
                if or((iCase == 1), (iCase ==2))
                        xlabel(ax{iCase}, '');
                end
                %         if mod(iCase,2) == 1
                delete(cb{iCase});
                %         end

        end
        legend(ax{summaryAx}, ...
                'Location', 'best', ...
                'Interpreter', 'none');
        cbPlot = colorbar(ax{1}, 'Position', cbpos);
        plot_str = 'log_{10}(P_e [cm^{-3}s^{-1}])';
        cbPlot.Label.String = plot_str;
end


    
function [prodProfile, nCellsProfile] = coldProdShellAngleProfile( ...
    coldElecProd, Xc, Yc, Zc, cellVol, qHat, rMin, rMax, angleEdges)
%COLDPRODSHELLANGLEPROFILE
% Volume-weighted mean cold-electron production in a radial shell,
% binned by angle to qHat.
%
% coldElecProd units are assumed to be cm^-3 s^-1.
% The returned prodProfile is also cm^-3 s^-1.

    qHat = qHat(:) ./ norm(qHat);

    % Make sure coldElecProd has the same orientation as meshgrid output.
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

    r = sqrt(Xc.^2 + Yc.^2 + Zc.^2);

    rSafe = r;
    rSafe(rSafe == 0) = NaN;

    xHat = Xc ./ rSafe;
    yHat = Yc ./ rSafe;
    zHat = Zc ./ rSafe;

    cosAng = xHat .* qHat(1) + yHat .* qHat(2) + zHat .* qHat(3);
    cosAng = max(min(cosAng, 1), -1);

    angleDeg = acosd(cosAng);

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

    [~, ~, binIdx] = histcounts(angleVals, angleEdges);

    inBin = binIdx > 0;

    binIdx = binIdx(inBin);
    prodVals = prodVals(inBin);
    volVals = volVals(inBin);

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

    bad = ~isfinite(prodProfile) | prodProfile <= 0 | nCellsProfile <= 0;
    prodProfile(bad) = NaN;

end