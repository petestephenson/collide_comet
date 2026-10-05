
function eImpactIoni = calcEIoniFreq(plas_dens_E, E_grid, XS_Eng, XS, inc_gens)
% Calculate e-impact ioni rate
    qe = 1.6e-19;
    me = 9.11e-31;
    % Mid points of energy grid
    E_mid = 0.5 * (E_grid(1:end-1) + E_grid(2:end));
    E_mid = reshape(E_mid, 1, 1, 1, []);
    ix_ioni = 7;
    
    % Get coll XS.
    XS_interped = reshape(interp1(XS_Eng, sum(XS(:,ix_ioni:end),2), E_mid), ...
        1,1,1,[]);  
    % Set XS to 0 for nans
    XS_interped(isnan(XS_interped)) = 0;
    % Convert from edf to pl flux
    elFlux = sum(plas_dens_E(:,:,:,:,inc_gens), 5) .* ...
        sqrt(E_mid * 2 * qe ./ me) * 1e2; 
    % Integrate over energy
    eImpactIoni = trapz(reshape(E_mid,[],1), XS_interped .* elFlux, 4);
end
