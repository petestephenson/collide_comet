
function smoothedPlasDens = smooth_plas_dens_E(plas_dens_E)
    nDims = size(plas_dens_E);
    plot_array_smooth = zeros(size(plas_dens_E));

   
    temp_array = plas_dens_E([1, 1:end, end], [1, 1:end, end], [1, 1:end, end],:);
    
    for iX = 1:3
        ix_X = [1:nDims(1)] + iX - 1;
        for iY = 1:3
            ix_Y = [1:nDims(2)] + iY - 1;
            for iZ = 1:3
                ix_Z = [1:nDims(3)] + iZ - 1;
                plot_array_smooth = plot_array_smooth + ...
                    temp_array(ix_X, ix_Y, ix_Z,:);
                
            end
        end
    end
    smoothedPlasDens = plot_array_smooth / 27;

end