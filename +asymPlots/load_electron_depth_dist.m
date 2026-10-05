
function [eDepthRadGrid, eDepthGrid, eDepthDist] = load_electron_depth_dist(runFName)
    ncid = netcdf.open(runFName);
    eDepth_id = netcdf.inqVarID(ncid, 'Electron Depth Dist');
    eDepthDist = netcdf.getVar(ncid, eDepth_id);


    eDepthRadGrid_id = netcdf.inqVarID(ncid, "elDepthRadGrid");
    eDepthRadGrid = netcdf.getVar(ncid, eDepthRadGrid_id);
    
    eDepthGrid_id = netcdf.inqVarID(ncid, 'elDepthGrid');
    eDepthGrid = netcdf.getVar(ncid, eDepthGrid_id);
    eDepthGrid = logspace(-5,4,size(eDepthGrid, 1));
    netcdf.close(ncid)

end