
function [plThetaGrid, plThetaDist] = load_electron_theta_dist(runFName)
    ncid = netcdf.open(runFName);
    plThetaDist_id = netcdf.inqVarID(ncid, 'Particle Theta Dist');
    plThetaDist = netcdf.getVar(ncid, plThetaDist_id);

    plThetaGrid_id = netcdf.inqVarID(ncid, 'elDepthThetaGrid');
    plThetaGrid = netcdf.getVar(ncid, plThetaGrid_id);
%     eDepthRadGrid_id = netcdf.inqVarID(ncid, "elDepthThetaGrid");
%     eDepthRadGrid = netcdf.getVar(ncid, eDepthRadGrid_id);
%     
%     eDepthGrid_id = netcdf.inqVarID(ncid, 'elDepthGrid');
%     eDepthGrid = netcdf.getVar(ncid, eDepthGrid_id);
%     eDepthGrid = logspace(-5,4,size(eDepthGrid, 1));
    netcdf.close(ncid)

end