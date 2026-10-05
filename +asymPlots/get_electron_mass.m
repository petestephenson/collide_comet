
function me = get_electron_mass(runFName)
    ncid = netcdf.open(runFName);
    me_id = netcdf.inqVarID(ncid, 'Electron Mass');
    me = netcdf.getVar(ncid, me_id);
    netcdf.close(ncid)
end