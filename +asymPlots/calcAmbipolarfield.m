% Calculate ambipolar field from the electron pressure and the x,y, z grids
function eAmbi = calcAmbipolarfield(pressure, plas_dens, x_grid, y_grid, z_grid)

    % Calculate ambipolar field
    % input pressure in Pa
    % input plasma density in m^-3
    % P = n*kB*T
    % n = plas_dens
    % T = temperature
    kb = 1.38e-23;
    e = 1.602e-19;
%     plas_dens = plas_dens; % convert to m^-3
    % pressure = pressure; % convert to J
    % Calculate gradient of pressure in x, y, z
    [dPdx, dPdy, dPdz] = gradient(pressure/kb, x_grid, y_grid, z_grid);

    % Calculate ambipolar field
    % Make array of zeros with shape of (3,)+ shape of pressure
%     eAmbi = zeros(3, size(pressure));
    eAmbi = zeros(3, size(pressure, 1), size(pressure, 2), size(pressure, 3));

    eAmbi(1, :, :, :) = -dPdx./(1e6*plas_dens*e);
    eAmbi(2, :, :, :) = -dPdy./(1e6*plas_dens*e);
    eAmbi(3, :, :, :) = -dPdz./(1e6*plas_dens*e);
    % Divide by 1e3 to convert X mesh from km to m
    eAmbi = eAmbi/1e3;
end


