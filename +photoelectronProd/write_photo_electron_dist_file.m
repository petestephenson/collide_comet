
function write_photo_electron_dist_file(E_grid, pe_CDF, saveDir)
% WRITE_PHOTO_ELECTRON_DIST_FILE Writes photoelectron distribution to a file.
%   WRITE_PHOTO_ELECTRON_DIST_FILE(E_GRID, PE_CDF, SAVE_DIR) writes the photoelectron
%   distribution to a netCDF file at the specified save directory.
%
%   Inputs:
%   - E_grid: Energy grid.
%   - pe_CDF: Cumulative distribution function for photoelectron.
%   - saveDir: Directory to save the file.
maxEng = 200; % eV
nEngSave = find(E_grid >= maxEng);
pe_CDF_save = pe_CDF(1:nEngSave)/pe_CDF(nEngSave); %Scale to 1
E_grid_save = E_grid(1:nEngSave);

cdf_name = [saveDir, 'PE_dist.nc'];

nccreate(cdf_name, 'PE_CDF',...
    'Dimensions',{'nEng', nEngSave}, 'Format','classic');
nccreate(cdf_name, 'PE_Eng',...
    'Dimensions', {'nEng', nEngSave}, 'Format','classic');
ncwrite(cdf_name, 'PE_CDF', pe_CDF_save);
ncwrite(cdf_name, 'PE_Eng', E_grid_save);
end
