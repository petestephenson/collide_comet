% Write function to calculate electron pressure from input denisty and temperature arrays
function [pressure] = calculatePressure(plas_dens, temperature)
    % Calculate electron pressure
    % input plasma density in cm^-3
    % input temperature in eV
    % P = n*kB*T
    % n = plas_dens
    % T = temperature
    kb = 1.38e-23;
    plas_dens = plas_dens*1e6; % convert to m^-3
    temperature = temperature*1.602e-19; % convert to J
    
    pressure = plas_dens.*temperature*kb;
end