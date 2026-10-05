
function eiiProd = calculateEIIProduction(plas_dens_E, grids, qDayside, qUniform,qVec,XS_Eng, XS,inc_gens)

    E_grid  = grids{4};
    xMid = reshape(0.5 * (grids{1}(2:end) + grids{1}(1:end-1)), [], 1);
    yMid = reshape(0.5 * (grids{2}(2:end) + grids{2}(1:end-1)), 1, []);
    zMid = reshape(0.5 * (grids{3}(2:end) + grids{3}(1:end-1)), 1, 1, []);
%     rMid = sqrt(xMid.^2 + yMid.^2 + zMid.^2);
%     nDens = Q * 1e-15 / (4 * pi * u_gas .* rMid.^2);
    nDens = asymPlots.calcNeutralDens(xMid, yMid, zMid, qUniform, qDayside, qVec);
    eImpactIoni = asymPlots.calcEIoniFreq(plas_dens_E, E_grid, XS_Eng, XS, inc_gens);
    eiiProd = nDens .* eImpactIoni; % cm^-3s^-1

end

