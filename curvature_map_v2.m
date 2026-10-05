%% Generate curvature map

%% Import E and B fields
EB_file = 'EB_fields_210F1_trim.nc';
ncid = netcdf.open(EB_file);
finfo = ncinfo(EB_file);
E_ave = netcdf.getVar(ncid,0);
B_ave = netcdf.getVar(ncid,1);
nX = size(E_ave,1); nY = size(E_ave,2); nZ = size(E_ave,3);
%% Define start positions of lines
% Lines should be predominantly in the y direction with some draping
% need array in x and z with set y values. 
% May require some adjustment for points in the cavity behind the comet. 
% Possibly second array here. 

x_start = linspace(2, numel(
