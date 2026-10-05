%% -------- Checking for OS ----------------------
OS_name = system_dependent('getos');
i_windows = strcmp('Microsoft',OS_name(1:9));
if i_windows
%     one_drive = 'C:/Users/pps18/OneDrive - Imperial College London/';
    one_drive = 'C:/Users/pete4/OneDrive - Imperial College London/';
else
    one_drive = '/Users/pete4/OneDrive - Imperial College London/';
end
fieldDir = [one_drive, 'Energy_Degradation/EB_fields/Asym/'];
%%
iFields = 1;
iMoments = 0;
iSym = 0;
iAsym = 1;
caseName = 'Ref';
% if iSym 
%     E_in = dlmread([fieldDir,'0118/E_Jan_cycle40000.vtk'],' ',9,0); %288 x 224 x 224
%     B_in = dlmread([fieldDir,'0118/B_Jan_cycle40000.vtk'],' ',9,0); %288 x 224 x 224
if iFields
    E_in = readmatrix([fieldDir,'Case', caseName, '/E_Jan_cycle40000.vtk'],'NumHeaderLines', 9, 'ExpectedNumVariables',3, 'Delimiter', ' ', 'FileType', 'text'); %288 x 224 x 224
    B_in = readmatrix([fieldDir,'Case', caseName, '/B_Jan_cycle40000.vtk'],'NumHeaderLines', 9, 'ExpectedNumVariables',3, 'Delimiter', ' ', 'FileType', 'text'); %288 x 224 x 224
elseif iMoments
    rho_in = readmatrix([fieldDir,'0118/rho4S_Jan_cycle40000.vtk'],'NumHeaderLines', 10, 'ExpectedNumVariables',1, 'Delimiter', ' ', 'FileType', 'text'); %288 x 224 x 224
end
% %     rho_in = dlmread([one_drive,'Energy_Degradation/PICFields/Q1e26new/rho4S_Jan_cycle40000.vtk'],' ',[10,0,23887881,0]);
% elseif iAsym
% %     E_in = dlmread([fieldDir,'0203/E_Jan_cycle40000.vtk'],' ',9,0); %288 x 224 x 224
%     if iFields
%         E_in = readmatrix([fieldDir,'0203/E_Jan_cycle40000.vtk'],'NumHeaderLines', 9, 'ExpectedNumVariables',3, 'Delimiter', ' ', 'FileType', 'text'); %288 x 224 x 224
%         B_in = readmatrix([fieldDir,'0203/B_Jan_cycle40000.vtk'],'NumHeaderLines', 9, 'ExpectedNumVariables',3, 'Delimiter', ' ', 'FileType', 'text'); %288 x 224 x 224
%     elseif iMoments   
%         rho_in = readmatrix([fieldDir,'0203/rho4S_Jan_cycle40000.vtk'],'NumHeaderLines', 10, 'ExpectedNumVariables',1, 'Delimiter', ' ', 'FileType', 'text'); %288 x 224 x 224
%     end
% %     B_in = dlmread([fieldDir,'0203/B_Jan_cycle40000.vtk'],' ',9,0); %288 x 224 x 224
% end

nX = 288; nY = 288; nZ = 288;
% i_Q1e25 = 0;
% i_Q1e26 = 0;
% i_Q1e26new = 1;
% 
% if i_Q1e25
%     E_in = dlmread([one_drive,'Energy_Degradation/PICFields/Q1e25/E_Jan_run210F5_average.vtk'],' ',9,0);
%     B_in = dlmread([one_drive,'Energy_Degradation/PICFields/Q1e25/B_Jan_run210F5_average.vtk'],' ',9,0);
%     nX = 288; nY = 192; nZ = 192;
% 
% elseif i_Q1e26
%     E_in = dlmread([one_drive,'Energy_Degradation/PICFields/Q1e26/E_Jan_run210F1_average.vtk'],' ',9,0);
%     B_in = dlmread([one_drive,'Energy_Degradation/PICFields/Q1e26/B_Jan_run210F1_average.vtk'],' ',9,0);   
%     nX = 288; nY = 224; nZ = 224;
% elseif i_Q1e26new
%     E_in = dlmread([one_drive,'Energy_Degradation/PICFields/Q1e26new/E_Jan_cycle40000.vtk'],' ',9,0);
%     B_in = dlmread([one_drive,'Energy_Degradation/PICFields/Q1e26new/B_Jan_cycle40000.vtk'],' ',9,0);   
%     rho_in = dlmread([one_drive,'Energy_Degradation/PICFields/Q1e26new/rho4S_Jan_cycle40000.vtk'],' ',[10,0,23887881,0]);
%     nX = 288; nY = 288; nZ = 288;
% end


%% Density and pressure matrices
rho2 = reshape([nan;nan;rho_in], [], 4);
% rho2 = rho2,[],4)
rhoComEl = reshape(rho2(3:end,3), nX, nY, nZ);
rhoComIons = reshape(rho2(3:end,4), nX, nY, nZ);
rhoSWEl = reshape(rho2(3:end,1), nX, nY, nZ);
rhoSWIons = reshape(rho2(3:end,2), nX, nY, nZ);
if iSym 
    save('Rho_Q_3e26_sym.mat', "rhoComEl", "rhoComIons", "rhoSWEl", "rhoSWIons");
elseif iAsym
    save('Rho_Q_3e26_Asym.mat', "rhoComEl", "rhoComIons", "rhoSWEl", "rhoSWIons");
end
%% EB Fields 
B2 = reshape(B_in, nX, nY, nZ, 3);
E2 = reshape(E_in, nX, nY, nZ, 3);
% if i_Q1e26new
%     rho_ave = reshape(rho_in, nX, nY, nZ);
% end
u_sw = 300e3; %ms^-1 - SW velocity
By_up = 6e-9; %T - By upstream in SW
Ez_up = - u_sw * By_up;

By_PIC_max = max(B2(1,:,:,2),[],'all');
By_PIC_min = min(B2(1,:,:,2),[],'all');
Ez_PIC_max = max(E2(1,:,:,3),[],'all');
Ez_PIC_min = min(E2(1,:,:,3),[],'all');

if By_PIC_max ~= By_PIC_min
    error('B_y varies on upstream boundary')
end
if Ez_PIC_max ~= Ez_PIC_min
     error('E_z varies on upstream boundary')
end

E_factor = Ez_up / Ez_PIC_max;
B_factor = By_up / By_PIC_max;

E_ave = E2 * E_factor;
B_ave = B2 * B_factor;
% clear E2 B2 By_PIC_max By_PIC_min Ez_PIC_max Ez_PIC_min E_factor B_factor

%%
% rho_ave = rhoComEl;

dim_name = {'','X','Y','Z'};
% nX = size(B_ave, 1); nY = size(B_ave, 2); nZ = size(B_ave, 3);
bounds = 1200;
i_plot_E = 0;
i_plot_B = 1;
i_plotRho = 0;
i_dim = 2;

i_slice = 3;
%---- viewing angle---
i_xy_plane = 0;
i_xz_plane = 1;

if i_xy_plane
    v_angle = [0,0,1];
elseif i_xz_plane
    v_angle = [0,-1,0];
else 
    v_angle = [1,0,0];
end


if i_plot_E
    title_str = ['E_{', dim_name{i_dim+1},'}'];
    cstr = ['E_{', dim_name{i_dim+1},'} [Vm^{-1}]'];
elseif i_plot_B
    title_str = ['B_{', dim_name{i_dim+1},'}'];
    cstr = ['B_{', dim_name{i_dim+1},'} [nT]'];
elseif i_plotRho
    title_str = ['\rho_e'];
    cstr = ['\rho_e [cm^{-3}]'];
end

xslice = [0];
yslice = [0];%(nY/2 + 1)];
zslice = [0];%[(nZ/2 +1), (nZ/2 +1)+ 10];
% x_comet = [134,113,113];
% x_comet = [134,113,113];
% x_comet = [134,97,97];
x_comet = [134, 145, 145];
% x_comet = []
% dx_EB = 10 * (10/15);
% dx_EB = 10;
dx_EB = 15 * 227 /288;
nslice = numel(xslice) + numel(yslice) + numel(zslice);
figure;
if i_plot_E
    plot_array = permute(E_ave(:,:,:,:),[2,1,3,4]);
elseif i_plot_B
    plot_array = permute(B_ave(:,:,:,:),[2,1,3,4])*1e9;
elseif i_plotRho 
    plot_array = log10(abs(permute(rho_ave(:,:,:),[2,1,3])));
end
if or(i_plot_B, i_plot_E)
    if i_dim == 0
        plot_array = sqrt(sum(plot_array.^2, 4));
    else 
        plot_array = plot_array(:,:,:,i_dim);
    end
end
x_grid2 = ([1:nX] - x_comet(1)) * dx_EB;
y_grid2 = ([1:nY] - x_comet(2)) * dx_EB;
z_grid2 = ([1:nZ] - x_comet(3)) * dx_EB;
[X2,Y2,Z2] = meshgrid(x_grid2, y_grid2, z_grid2);

ax_slice = slice(X2,Y2,Z2,plot_array, xslice,yslice,zslice);
for i_slice = 1:nslice
    ax_slice(i_slice).EdgeColor = 'none';
end
title(title_str);
colormap('jet');
% rb_map = zeros(256, 3);
% rb_map = redblue(256);
% colormap(rb_map);
c = colorbar;
c.Location = 'west';
xlim([-bounds, bounds]);
ylim([-bounds, bounds]);
zlim([-bounds, bounds]);
% caxis([14,18]);
c.Label.String = cstr;
xlabel('X [km]');ylabel('Y [km]');zlabel('Z [km]');
daspect([1,1,1]);
view(v_angle);
set(gca, 'Fontsize', 15);

%% Trim downstream edges for 201F1 run
% E_ave2 = E_ave(1:end-5,:,:,:);
% B_ave2 = B_ave(1:end-5,:,:,:);

%% Create netcdf file with E and B field data

% if i_Q1e25
%     cdf_name = 'EB_fields_210F5.nc';
% elseif i_Q1e26
%     cdf_name = 'EB_fields_210F1_trim.nc';
% elseif i_Q1e26new
%     cdf_name = 'EB_fields_Q1e26.nc';
% end
% if iSym
%     cdf_name = 'EB_fields/EB_fields_Q3e26.nc';
% 
% % elseif iAsym
% cdf_name = ['EB_fields/EB_fields_Q5e26_Asym_Case', caseName, '.nc'];
% % end
% % 
% nccreate(cdf_name, 'E_field',...
%     'Dimensions',{'nX',nX, 'nY', nY, 'nZ', nZ, 'nDims',3}, 'Format','classic');
% nccreate(cdf_name, 'B_field',...
%     'Dimensions',{'nX',nX, 'nY', nY, 'nZ', nZ, 'nDims',3}, 'Format','classic');
% 
% ncwrite(cdf_name, 'E_field', E_ave);
% ncwrite(cdf_name, 'B_field', B_ave);


%%
% clear;
EB_Folder = 'EB_fields';
% EB_file = 'EB_fields_210F5.nc';
% EB_file = 'EB_fields_210F1_trim.nc';
% EB_file = 'EB_fields_Q3e26_Asym.nc';
caseName = 'Ref';
EB_file = ['EB_fields_Q5e26_Asym_Case', caseName, '.nc'];
EB_file = [EB_Folder,'\', EB_file];
ncid = netcdf.open(EB_file);
finfo = ncinfo(EB_file);
E_ave = netcdf.getVar(ncid,0);
B_ave = netcdf.getVar(ncid,1);
[~, nX] = netcdf.inqDim(ncid,0);
[~, nY] = netcdf.inqDim(ncid,1);
[~, nZ] = netcdf.inqDim(ncid,2);
% z_grid = netcdf.getVar(ncid,2);
% E_grid = netcdf.getVar(ncid,3);

netcdf.close(ncid);

