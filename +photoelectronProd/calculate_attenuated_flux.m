
function attenSolarFlux = calculate_attenuated_flux(losColDensity,solarFlux, intAbsXS)
    optDepth = intAbsXS .* losColDensity;
    attenSolarFlux = solarFlux .* exp(-optDepth);


end