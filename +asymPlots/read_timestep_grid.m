
function [dtGrid, dtCount] = read_timestep_grid(runFName)
    ncid = netcdf.open(runFName);
    dtGrid_id = netcdf.inqVarID(ncid, "dT Boundaries");
    dtCount_id = netcdf.inqVarID(ncid, "dT Count");
    dtGrid = netcdf.getVar(ncid, dtGrid_id);
    dtCount = netcdf.getVar(ncid, dtCount_id);
    netcdf.close(ncid);
end