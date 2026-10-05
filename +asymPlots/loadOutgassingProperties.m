function [qSpherical, qDayside, qUniform, qvec] = loadOutgassingProperties(runFname)
%     disp(runFname);
    ncid = netcdf.open(runFname);
    qVecID = netcdf.inqVarID(ncid, 'Q max vector');
    qUniID = netcdf.inqVarID(ncid, 'OGR Uniform');
    qDaysideID = netcdf.inqVarID(ncid, 'OGR Dayside');
    qID = netcdf.inqVarID(ncid, 'OGR');
    qvec = netcdf.getVar(ncid, qVecID);
    qUniform = netcdf.getVar(ncid, qUniID);
    qDayside = netcdf.getVar(ncid, qDaysideID);
    qSpherical = netcdf.getVar(ncid, qID);
    netcdf.close(ncid);
end


% function [qSpherical, qDayside, qUniform, qvec] = loadOutgassingProperties(runFname)
%     disp(runFname);
%     ncid = netcdf.open(runFname);
%     qVecID = netcdf.inqVarID(ncid, 'Q max vector');
%     qUniID = netcdf.inqVarID(ncid, 'OGR Uniform');
%     qDaysideID = netcdf.inqVarID(ncid, 'OGR Dayside');
%     qID = netcdf.inqVarID(ncid, 'OGR');
%     qvec = netcdf.getVar(ncid, qVecID);
%     qUniform = netcdf.getVar(ncid, qUniID);
%     qDayside = netcdf.getVar(ncid, qDaysideID);
%     qSpherical = netcdf.getVar(ncid, qID);
%     netcdf.close(ncid);
% end