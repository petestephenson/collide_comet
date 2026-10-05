function [grids, plas_dens_E] = loadRunData(runFname)
    ncid = netcdf.open(runFname);
    x_id = netcdf.inqVarID(ncid, 'X Boundaries'); 
    x_grid = netcdf.getVar(ncid,x_id);
    y_id = netcdf.inqVarID(ncid, 'Y Boundaries'); 
    y_grid = netcdf.getVar(ncid,y_id);
    z_id = netcdf.inqVarID(ncid, 'Z Boundaries');
    z_grid = netcdf.getVar(ncid,z_id);
    E_id = netcdf.inqVarID(ncid, 'Energy Boundaries');
    E_grid = netcdf.getVar(ncid,E_id);

    grids = {x_grid, y_grid, z_grid, E_grid};
    
    plas_id = netcdf.inqVarID(ncid, 'Plasma Density');
    plas_dens_E = netcdf.getVar(ncid,plas_id);

    netcdf.close(ncid);
end