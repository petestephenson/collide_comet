function [bulk_vel, square_vel] = load_velocity_arrays(runFName)
    ncid = netcdf.open(runFName);
    bulk_v_id = netcdf.inqVarID(ncid, 'Bulk Velocity');
    bulk_vel = netcdf.getVar(ncid,bulk_v_id);
  
    % Read square velocities for temperature
    sq_v_id = netcdf.inqVarID(ncid, 'Square Velocity');
    square_vel = netcdf.getVar(ncid,sq_v_id);
    square_vel(isnan(square_vel)) = 0;
    netcdf.close(ncid)
end
