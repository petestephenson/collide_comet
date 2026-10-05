
function coldElecProd = calculateColdeProduction(plas_dens_E, grids, E_cold, qDayside, qUniform,qVec,XS_Eng, XS,E_th, inc_gens)

    % Load cross-section data
%     load([EngDegDir, '\xsections.mat']);
    ix_ioni = 7;

    nEng = numel(grids{4}) - 1;
    nXS = numel(E_th);

        % Calculating sum(Int f(E) * XS(E) n_H2O(x) dE)
    xMid = reshape(0.5 * (grids{1}(2:end) + grids{1}(1:end-1)), [], 1);
    yMid = reshape(0.5 * (grids{2}(2:end) + grids{2}(1:end-1)), 1, []);
    zMid = reshape(0.5 * (grids{3}(2:end) + grids{3}(1:end-1)), 1, 1, []);
    rMid = sqrt(xMid.^2 + yMid.^2 + zMid.^2);
%     nDens = Q * 1e-15 / (4 * pi * u_gas .* rMid.^2);
    nDens = asymPlots.calcNeutralDens(xMid, yMid, zMid, qUniform, qDayside, qVec);

    E_grid_mid = 0.5 * (grids{4}(2:end) + grids{4}(1:end-1));
    dE_grid = (grids{4}(2:end) - grids{4}(1:end-1));

    % Initialize XS_total
    XS_total = zeros(1, nEng);

    for i_xs = 1:nXS
        if E_th(i_xs) > 0
            % Determine Limits for each process
            E_lims = [max(E_th(i_xs), E_cold), E_cold + E_th(i_xs)];

            if i_xs >= ix_ioni
                E_lims(2) = E_lims(2) + E_cold;
            end

            ix_E_th = find(grids{4} > E_lims(1), 1, 'first') - 1;
            ix_E_cold = find(grids{4} > E_lims(2), 1, 'first') - 1;

            eng_low = 0.5 * (E_lims(1) + grids{4}(ix_E_th + 1));
            eng_upp = 0.5 * (E_lims(2) + grids{4}(ix_E_cold));

            low_factor = (grids{4}(ix_E_th + 1) - E_lims(1)) / ...
                (grids{4}(ix_E_th + 1) - grids{4}(ix_E_th));
            upp_factor = (E_lims(2) - grids{4}(ix_E_cold)) / ...
                (grids{4}(ix_E_cold + 1) - grids{4}(ix_E_cold));

            XS_Eng2 = E_grid_mid';
            XS_Eng2([ix_E_th, ix_E_cold]) = [eng_low, eng_upp];
            XS_interp = interp1(XS_Eng, XS(:, i_xs), XS_Eng2);

            XS_interp([ix_E_th, ix_E_cold]) = XS_interp([ix_E_th, ix_E_cold]) .* ...
                [low_factor, upp_factor];
            XS_interp([1:ix_E_th-1, ix_E_cold + 1:nEng]) = 0;

            XS_total = XS_total + XS_interp;
        end
    end
    qe = 1.6e-19;
    me = 9.11e-31;
    v_eng = sqrt(2 * XS_Eng2 * abs(qe) / me);
    temp = reshape(XS_total .* dE_grid' .* v_eng * 1e2, 1, 1, 1, []);
    
%     for iRun = 1:2
    coldEProdPerMol = sum(plas_dens_E(:,:,:,:,inc_gens) .* temp, [4, 5]);
    coldElecProd = coldEProdPerMol .* nDens;
%     end
end