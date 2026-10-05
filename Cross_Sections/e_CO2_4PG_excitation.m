band_0 = csvread('/Users/pps18/Documents/Cross_Sections/CO2_4PG_excitation_bands/CO2_4PG_excitation_band_0.txt',1,0);
band_1 = csvread('/Users/pps18/Documents/Cross_Sections/CO2_4PG_excitation_bands/CO2_4PG_excitation_band_1.txt',1,0);
band_2 = csvread('/Users/pps18/Documents/Cross_Sections/CO2_4PG_excitation_bands/CO2_4PG_excitation_band_2.txt',1,0);
band_3 = csvread('/Users/pps18/Documents/Cross_Sections/CO2_4PG_excitation_bands/CO2_4PG_excitation_band_3.txt',1,0);
band_4 = csvread('/Users/pps18/Documents/Cross_Sections/CO2_4PG_excitation_bands/CO2_4PG_excitation_band_4.txt',1,0);
energy_grid = logspace(log10(4),log10(400),100);
bands = {band_0,band_1, band_2, band_3, band_4};
excitation_bands = cellfun(@(x) interp1(x(:,1),x(:,2),energy_grid'), bands,'UniformOutput',false);
excitation_bands = [excitation_bands{:}];
% 
% test = band_0(2:end,1)-band_0(1:end-1,1);
% find(test<=0)
% 
figure;
% plot(band_2(:,1),band_2(:,2))
for i = 1:5
    plot(energy_grid,excitation_bands(:,i));
    hold on
end

XS_at_300 = interp1(energy_grid',excitation_bands,300);

writematrix([energy_grid',excitation_bands],'/Users/pps18/Documents/Cross_Sections/CO2_4PG_excitation_bands/CO2_4PG_band_shapes.txt')