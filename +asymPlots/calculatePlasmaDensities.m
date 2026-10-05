
function [plas_dens, plas_dens_cold] = calculatePlasmaDensities(grids, plas_dens_E,inc_gens, minEng, maxEng)
    dE_grid = grids{4}(2:end) - grids{4}(1:end-1);
    ix_min_eng = find(grids{4} > minEng, 1, 'first') - 1;
    ix_max_eng = find(grids{4} <= maxEng, 1, 'last') - 1;
    inc_eng = [ix_min_eng:ix_max_eng]; % Energy Range to plot

    plas_dens = sum(plas_dens_E(:,:,:,:,inc_gens) ...
        .* reshape(dE_grid, 1, 1, 1, []), [4, 5]);
    plas_dens_cold = sum(plas_dens_E(:,:,:,inc_eng, inc_gens) ...
        .* reshape(dE_grid(inc_eng), 1, 1, 1, []), [4, 5]);
    nDims = size(plas_dens);
    plot_array_smooth = zeros(nDims);
    plot_array_cold_smooth = zeros(nDims);
    temp_array = plas_dens([1, 1:end, end], [1, 1:end, end], [1, 1:end, end]);
    temp_array_cold = plas_dens_cold([1, 1:end, end], [1, 1:end, end], [1, 1:end, end]);
    for iX = 1:3
        ix_X = [1:nDims(1)] + iX - 1;
        for iY = 1:3
            ix_Y = [1:nDims(2)] + iY - 1;
            for iZ = 1:3
                ix_Z = [1:nDims(3)] + iZ - 1;
                plot_array_smooth = plot_array_smooth + ...
                    temp_array(ix_X, ix_Y, ix_Z);
                plot_array_cold_smooth = plot_array_cold_smooth + ...
                    temp_array_cold(ix_X, ix_Y, ix_Z);
            end
        end
    end
    plas_dens = plot_array_smooth / 27;
    plas_dens_cold = plot_array_cold_smooth / 27;

end