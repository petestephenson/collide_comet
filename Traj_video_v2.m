%% Read Obj shape file for Rosetta

shape_67P = readObj('1567260122412-ESA_Rosetta_OSIRIS_67P_SHAP2P.obj');


%% Load PIC E and B fields?
EB_folder = 'C:\Users\pete4\OneDrive - Imperial College London\Energy_Degradation\EB_fields\';
EB_file = [EB_folder,'EB_fields_210F5.nc'];
% EB_file = [EB_folder, 'EB_fields_210F1_trim.nc'];
ncid = netcdf.open(EB_file);
finfo = ncinfo(EB_file);
E_ave = netcdf.getVar(ncid,0);
B_ave = netcdf.getVar(ncid,1);
nX = size(B_ave,1);
nY = size(B_ave,2);
nZ = size(B_ave,3);
%% Create an arrow

introvideo = VideoWriter('cometintrov2');
introvideo.FrameRate = 15;
f = figure('Position', [100,100,1500,1000]);
ax0 = axes(f);
hold(ax0, 'on');
r_67P = 100;
% figure;
x_comet = [134, 97, 97];
% x_comet = [134,113,113];
patch('Faces', shape_67P.f.v, 'Vertices', shape_67P.v *r_67P, ...
    'EdgeColor', 'none','FaceColor',[0.5,0.5,0.5]);
daspect([1,1,1]);

% [hline, hhead] = arrow3d([0,0,0], [1,1,1]);
% hold on
%--------- Set limits of plot ---------------------------
xlims = [-x_comet(1)*10,1000];
ylims = [-x_comet(2)*10,x_comet(2)*10];
xlim(xlims);
ylim(ylims);
zlim([-x_comet(3)*10, x_comet(3)*10]);


set(ax0, 'color', 'k');
% c_photon = [1,0.6,0.3];
c_photon = [238,130,238]/255;
% h2 = mArrow3([0,0,0],[20,100,100], 'color', c_photon, 'tipWidth', 20);
light('Position', [xlims(1),0,0], 'Style', 'infinite');

%------------ Set up photons ------------------------
nPhotons = 5;
nPhotPts = 50;
photons = struct();

% Create photon
l_photon = 200;
h_photon = 15;
lArrow = 40;
photon_spacing = 300;
photon_shift = 1300; %How far to shift photons in - x direction initially
grid_int = 2;


y_in = [-grid_int:grid_int]*photon_spacing;
z_in = zeros(size(y_in));
nPhotons = numel(y_in);

%------------- initial photon positions -----------------------------------
for i_photon = 1:nPhotons
%     photons(i_photon).X_init = [500*(2* rand() -1),y_in(i_photon),z_in(i_photon)]; %+ [xlims(1),0,0]
    photons(i_photon).X_init = [- l_photon - photon_shift,y_in(i_photon),z_in(i_photon)]; %+ [xlims(1),0,0]
    photons(i_photon).X = zeros(nPhotPts, 3);
    photons(i_photon).X(:,1) = linspace(xlims(1), xlims(1) + l_photon, nPhotPts);
    photons(i_photon).X(:,2) = h_photon * cos(4 * pi *...
        (photons(i_photon).X(:,1)-xlims(1)) / l_photon);
    photons(i_photon).X(:,3) = h_photon * sin(4 * pi *...
        (photons(i_photon).X(:,1)-xlims(1)) / l_photon);
    photons(i_photon).X = photons(i_photon).X + photons(i_photon).X_init; 

    photons(i_photon).Line =  plot3(ax0, photons(i_photon).X(:,1),...
        photons(i_photon).X(:,2),photons(i_photon).X(:,3), 'color', c_photon,...
        'Linewidth', 2);
%     photon(i_photon).End = [x_p(end), y_p(end), z_p(end)];
%     p1_arrow = 
    photons(i_photon).Arrow = mArrow3(photons(i_photon).X(end,:),photons(i_photon).X(end,:)...
        + [lArrow,0,0], 'color', c_photon, 'tipWidth', 20);
    set(photons(i_photon).Arrow,'facelighting','none');
%     photons(i_photon).line =  plot3(ax0, x_p, y_p,z_p, 'color', c_photon,...
%     'Linewidth', 2);

end
% ----------------------------------------------------------------
% ------------------ Set up field lines --------------------------
% ----------------------------------------------------------------
% if i_plot_B_lines
nBLines = 5;
B_line_range = 200;
x_B_start = linspace(xlims(1) + 30, xlims(1) + 30 + B_line_range, nBLines);
y_B_start = ylims(1) + 200;

x_s = (x_B_start/10) + x_comet(1);
y_s = (y_B_start/10) + x_comet(2);

zslice = 0;
[~,ix_z] = min(abs(([1:nZ]-x_comet(3))*10-zslice));
[str_X,str_Y] = meshgrid([1:nY], [1:nX]);
B_temp = permute(B_ave, [2,1,3,4]);
hs = stream2(str_X, str_Y, B_ave(:,:,ix_z,2), B_ave(:,:,ix_z,1), repmat(y_s, size(x_s)),x_s);
hn = stream2(str_X, str_Y, -B_ave(:,:,ix_z,2), -B_ave(:,:,ix_z,1), repmat(y_s, size(x_s)),x_s);
plot_lines = find(and(~cellfun(@isempty,hs), ~cellfun(@isempty,hn)));

B_lines = struct();
for i_line = plot_lines
    lin_s = (hs{i_line} - x_comet([2,1]))*10;
    lin_n = (hn{i_line} - x_comet([2,1]))*10;
    lin_s = [lin_s, repmat(zslice, size(lin_s,1), 1)];
    lin_n = [lin_n, repmat(zslice, size(lin_n,1), 1)];
    B_lines(i_line).Data = [flipud(lin_n); lin_s(2:end,:)];
    B_lines(i_line).Line = line(ax0,B_lines(i_line).Data(:,2),B_lines(i_line).Data(:,1),...
        B_lines(i_line).Data(:,3),'Linewidth',0.8,'color','w');
%     B_lines(i_line,1).Line = line(ax0,lin_s(:,2),lin_s(:,1),lin_s(:,3),'Linewidth',0.8,'color','w');
%     B_lines(i_line,2).Line = line(ax0,lin_n(:,2),lin_n(:,1),lin_n(:,3),'Linewidth',0.8,'color','w');
end

% ----------------------------------------------------------------
% -------------- Create a heatwave -------------------------------
% ----------------------------------------------------------------
l_heat = 300;
[x_heat, y_heat] = meshgrid(linspace(xlims(1), xlims(1) + l_heat), linspace(ylims(1), ylims(2)));
z_heat = 40 +zeros(size(x_heat));
% [255, 160,19]/255
c_heat = repmat(permute([255, 94,19]/255,[1,3,2]), size(x_heat));
% c_heat = permute
% Opacity of heat wave 
x_max = x_heat(1,end);
x_min = x_heat(1, 1);

alpha_heat = sigmoid(x_heat, (x_min + x_max) /2,1/40);

% figure;
% plot(x_heat(1,:), alpha_heat);
heat_wave = surf(ax0, x_heat, y_heat, z_heat, c_heat, 'EdgeColor', 'none');
% cb = colorbar('Location', 'west');
heat_wave.AlphaData = alpha_heat;
heat_wave.FaceColor = 'interp';
heat_wave.FaceAlpha = 'interp';
heat_wave.EdgeColor = 'none';
heat_wave.FaceLighting = 'none';

% ----------------------------------------------------------------
% --------------- Set up SW pls. ---------------------------------
% ----------------------------------------------------------------
% Create photon
SWPls = struct();
size_SW_pls = 20;
col_SW_pls = 'r';
% lArrow = 50;
SW_spacing = 200;
% grid_int = floor(ceil(sqrt(nPhotons))/2);
grid_int = 3;
% initial photon positions grid
% [y_SW, z_SW] = meshgrid([-grid_int:grid_int]*photon_spacing,...
%     [-grid_int:grid_int]*photon_spacing);
y_SW = [-grid_int:grid_int]*SW_spacing;
z_SW = zeros(size(y_SW));
nSWPls = numel(y_SW);
[x_pl, y_pl, z_pl] = sphere;
x_pl = x_pl * size_SW_pls; y_pl = y_pl * size_SW_pls; z_pl = z_pl * size_SW_pls;

ix_SW_Bline = 5; % Choose which B field line to add pls to. 
for i_sw = 1:nSWPls
    
    SWPls(i_sw).X_centre(1) = interp1(B_lines(ix_SW_Bline).Data(:,1), B_lines(ix_SW_Bline).Data(:,2),...
        y_SW(i_sw));
    SWPls(i_sw).X_centre(2) = y_SW(i_sw);
    SWPls(i_sw).X_centre(3) = 0;
    
    SWPls(i_sw).X = zeros([size(x_pl),3]);
    
    SWPls(i_sw).X(:,:,1) = SWPls(i_sw).X_centre(1) + x_pl ;%+ 200*(2* rand() -1);
    SWPls(i_sw).X(:,:,2) = SWPls(i_sw).X_centre(2) + y_pl;
    SWPls(i_sw).X(:,:,3) = SWPls(i_sw).X_centre(3) + z_pl;
    
    SWPls(i_sw).Surface = surf(SWPls(i_sw).X(:,:,1), SWPls(i_sw).X(:,:,2),...
        SWPls(i_sw).X(:,:,3), 'FaceColor', col_SW_pls, 'Edgecolor','none');
    set(SWPls(i_sw).Surface, 'facelighting', 'gouraud');
%     photons(i_photon).Line =  plot3(ax0, photons(i_photon).X(:,1),...
%         photons(i_photon).X(:,2),photons(i_photon).X(:,3), 'color', c_photon,...
%         'Linewidth', 2);
    
%     photons(i_photon).X_init = [500*(2* rand() -1),y_in(i_photon),z_in(i_photon)]; %+ [xlims(1),0,0]
%     photons(i_photon).X = zeros(nPhotPts, 3);

end

%--------- Set up neutral gas -----------------------------------
size_cloud = 500;
x_gas = linspace(-size_cloud, size_cloud);
% y_gas = linspace(-100, 100);
% z_gas = linspace(-100, 100);
[Y_gas,X_gas,Z_gas] = meshgrid(x_gas, x_gas, x_gas);
r_gas = sqrt(Y_gas.^2 + X_gas.^2 + Z_gas.^2);%km
Q_gas = 1e25; %s^-1
u_gas = 1; %kms^-1

n_gas = Q_gas *1e-15 ./ (4 * pi * u_gas * r_gas.^2); 

r_gas_lin = linspace(0,size_cloud);
theta_gas = linspace(0,360);
[r_gas,theta_gas] = meshgrid(r_gas_lin, theta_gas);
x_gas = r_gas .* cosd(theta_gas);
y_gas = zeros(size(x_gas));
z_gas = r_gas .* sind(theta_gas);
% r_gas = sqrt(y_gas.^2 + x_gas.^2 + z_gas.^2);%km
pos_gas = [permute(x_gas,[3,4,1,2]),permute(y_gas,[3,4,1,2]), permute(z_gas,[3,4,1,2])];
c_gas_slice = [114,0,88]/255;
% c_gas_slice = [148,0,211]/255;
gas_step_init = 150; % Step after which neutral gas appears
nGasSteps = 100; % number of steps for gas to grow

% ix_major
% gas_slice = slice(ax0, Y_gas,X_gas, Z_gas, log10(n_gas), [],[],[0]);
% set(gas_slice, 'facelighting','none');
% gas_slice.EdgeColor = 'none';
% n_gas_slice = gas_slice.CData;
% gas_slice.CData(:,:) = 1;
% alpha_mat = 1 * (log10(n_gas_slice) - min(log10(n_gas_slice), [], 'all')) ./ ...
%                 (max(log10(n_gas_slice), [], 'all') - min(log10(n_gas_slice), [], 'all'));
% gas_slice.AlphaData = alpha_mat;
% gas_slice.FaceColor = 'interp';
% gas_slice.FaceAlpha = 'interp';
gas_map = [138,43,226] .* [1:6]'/(6*256);
% gas_map = [0 0.3 0.3
%     0 0.4 0.4
%     0 0 0.5
%     0 0 0.6
%     0 0 0.8
%     0 0 1.0];
% -------------------------------------------------------------------------
% ---------------------- Arrows to indicate gas outflow -------------------
% -------------------------------------------------------------------------
nGasArrow = 4;
gasArrowRad = 1.2 * r_67P;
gasArrowTheta = linspace(0,2*pi,nGasArrow + 1)';
gasArrowTheta = gasArrowTheta(1:end-1);
gasArrTipLength = 50; % Length of Arrow Tip
gasArrTipAngle = 75;
gasArrowMaxLength = 200;
nGasArrowGrowSteps = 40;
nGasArrowFadeSteps = 20;
dGasArrowLength = (gasArrowMaxLength - gasArrTipLength) / nGasArrowGrowSteps;
GasArrowFadeStart = gas_step_init + nGasArrowGrowSteps + 10;
gasArrows = struct();
gasArrowColor = [148,0,211]/255;

for i_arrow = 1:nGasArrow
    gasArrows(i_arrow).X_min = gasArrowRad * [cos(gasArrowTheta(i_arrow)),sin(gasArrowTheta(i_arrow)),...
        0];
    gasArrows(i_arrow).V_hat = [cos(gasArrowTheta(i_arrow)),sin(gasArrowTheta(i_arrow)),...
        0];
    gasArrows(i_arrow).length = min(gasArrTipLength, gasArrowMaxLength);
    gasArrows(i_arrow).X_max = (gasArrowRad + gasArrTipLength) * [cos(gasArrowTheta(i_arrow)),...
        sin(gasArrowTheta(i_arrow)), 0];
end
% -------------------------------------------------------------------------
% ---------------------- Create Photoelectrons ----------------------------
% -------------------------------------------------------------------------

photo_els = struct();
r_photo_el = r_67P * 1.1;
nPhotoEls = 10;
thetaPhEl = linspace(0, 2*pi, nPhotoEls);
StepInitPhEl = 320;
for i_photo_el = 1:nPhotoEls
    photo_els(i_photo_el).X_centre = r_photo_el * [cos(thetaPhEl(i_photo_el)),...
        sin(thetaPhEl(i_photo_el)),0];
    photo_els(i_photo_el).V_hat = photo_els(i_photo_el).X_centre ./ ...
                sqrt(sum(photo_els(i_photo_el).X_centre.^2));
    photo_els(i_photo_el).X = zeros([size(x_pl),3]);
    photo_els(i_photo_el).X(:,:,1) = photo_els(i_photo_el).X_centre(1) + x_pl;
    photo_els(i_photo_el).X(:,:,2) = photo_els(i_photo_el).X_centre(2) + y_pl;
    photo_els(i_photo_el).X(:,:,3) = photo_els(i_photo_el).X_centre(3) + z_pl;
end
% camtarget(photons(floor(nPhotons/2)).X(end,:));

campos([0,0,500]);
view(ax0,0,-90);
camva(5);

% Move photons
dx_ph = 8; % Distance change per frame for photon
dx_sw = 3; % Distance moved by SW electron
dx_ph_el = 0.5; % Distance change per frame for photoelectron
dx_ph_ion = 1; % Distance change per frame for photoelectron
dx_heatwave = dx_ph;
col_ph_ion = 'b';
col_ph_el = 'b';
nSteps = 600;
% nAbsPhotons = 0;
% nPhotoEls = 0;
open(introvideo);
for i_step = 1:nSteps
    % -------------------------------------------------------------------
    % ------------------- Plot Photons from upstream --------------------
    % -------------------------------------------------------------------
    for i_photon = 1:nPhotons
        delete(photons(i_photon).Line);
        delete(photons(i_photon).Arrow);
        photons(i_photon).X(:,1) = photons(i_photon).X(:,1) + dx_ph;
        r_photon = sqrt(sum(photons(i_photon).X(end,:).^2,2));
%         if r_photon > 1.5 * r_67P 
            photons(i_photon).Line =  plot3(ax0, photons(i_photon).X(:,1),...
            photons(i_photon).X(:,2),photons(i_photon).X(:,3), 'color', c_photon,...
            'Linewidth', 2);
        
            
            photons(i_photon).Arrow = mArrow3(photons(i_photon).X(end,:),photons(i_photon).X(end,:)...
                + [lArrow,0,0], 'color', c_photon, 'tipLength', lArrow);
            set(photons(i_photon).Arrow,'facelighting','none');
%             nAbsPhotons = nAbsPhotons + 1;
%         end
    end
    
    % ---------------------------------------------------------------------
    %------------------- Plot heatwave ------------------------------------
    delete(heat_wave);
    x_heat = x_heat + dx_heatwave;
    heat_wave = surf(x_heat, y_heat, z_heat, c_heat, 'EdgeColor', 'none');
%     cb = colorbar('Location', 'west');
    heat_wave.AlphaDataMapping = 'none';
    heat_wave.AlphaData = alpha_heat;
%     heat_wave.FaceColor = [1,0,0];
    heat_wave.FaceAlpha = 'flat';
    heat_wave.EdgeColor = 'none';
    heat_wave.FaceLighting = 'none';
    heat_wave.AmbientStrength = 1;


    
    % ---------------------------------------------------------------------
    % ------------------ Plot Magnetic Field lines ------------------------
    % ---------------------------------------------------------------------
    
   
    x_s = x_s + (dx_sw/10);
%     y_s = (y_B_start/10) + x_comet(2);

    zslice = 0;
%     [~,ix_z] = min(abs(([1:nZ]-x_comet(3))*10-zslice));
%     [str_X,str_Y] = meshgrid([1:nY], [1:nX]);
%     B_temp = permute(B_ave, [2,1,3,4]);
    hs = stream2(str_X, str_Y, B_ave(:,:,ix_z,2), B_ave(:,:,ix_z,1), repmat(y_s, size(x_s)),x_s);
    hn = stream2(str_X, str_Y, -B_ave(:,:,ix_z,2), -B_ave(:,:,ix_z,1), repmat(y_s, size(x_s)),x_s);
    plot_lines = find(and(~cellfun(@isempty,hs), ~cellfun(@isempty,hn)));
    for i_line = plot_lines
        delete(B_lines(i_line).Line);
%         delete(B_lines(i_line).Line);
        
        lin_s = (hs{i_line} - x_comet([2,1]))*10;
        lin_n = (hn{i_line} - x_comet([2,1]))*10;
        lin_s = [lin_s, repmat(zslice, size(lin_s,1), 1)];
        lin_n = [lin_n, repmat(zslice, size(lin_n,1), 1)];
        B_lines(i_line).Data = [flipud(lin_n); lin_s(2:end,:)];
        B_lines(i_line).Line = line(ax0,B_lines(i_line).Data(:,2),B_lines(i_line).Data(:,1),...
            B_lines(i_line).Data(:,3),'Linewidth',0.8,'color','w');
%         B_lines(i_line,1).Line = line(ax0,lin_s(:,2),lin_s(:,1),lin_s(:,3),'Linewidth',0.8,'color','w');
%         B_lines(i_line,2).Line = line(ax0,lin_n(:,2),lin_n(:,1),lin_n(:,3),'Linewidth',0.8,'color','w');
    end
    
    
    % -------------------------------------------------------------------
    % ------------------- Plot neutral Gas Cloud ------------------------
    % -------------------------------------------------------------------
    
    if exist('gas_slice')
        delete(gas_slice)
    end
    
    if i_step > gas_step_init
        [az, el] = view();
        cam_vec = [- cosd(el) * cosd(az + 90), - cosd(el) * sind(az+ 90),...
            sind(el)];
    %     cam_vec = ([0,0,0] - campos())/norm([0,0,0] - campos());
        % -----------Calculate rotation used for neutral gas slice ------
        
        rotvec = vrrotvec([0,1,0],cam_vec);
        rot_mat = vrrotvec2mat(rotvec);
        pos_gas_plot = permute(sum(rot_mat .* pos_gas,2),[3,4,1,2]);
        % -------- Determine size of neutral gas slice --------------------
        % For growth of cloud 
        r_gas_max = r_67P + (size_cloud - r_67P) * (i_step - gas_step_init)/nGasSteps;
        pos_gas_plot = pos_gas_plot(:,r_gas(1,:) <r_gas_max,:);
        
        % ------------ Plot the slice for the neutral gas ----------------
        gas_slice = slice(ax0, Y_gas,X_gas, Z_gas, log10(n_gas), pos_gas_plot(:,:,1),...
            pos_gas_plot(:,:,2),pos_gas_plot(:,:,3));
        set(gas_slice,'facelighting','none');
        
        % ---------- Calculate transparency for neutral gas slice ------
        n_gas_slice = gas_slice.CData;
        
        alpha_mat = 0.2 + 0.8 * (log10(n_gas_slice) - min(log10(n_gas_slice), [], 'all')) ./ ...
                        (max(log10(n_gas_slice), [], 'all') - min(log10(n_gas_slice), [], 'all'));
%         gas_slice.CData(:,:) = 1;
        gas_slice.CData = repmat(permute([114,0,88]/255,[1,3,2]), size(n_gas_slice));
%         colormap( gas_map);
        gas_slice.AlphaDataMapping = 'none';
        gas_slice.AlphaData = alpha_mat;
        gas_slice.FaceColor = 'interp';
        gas_slice.FaceAlpha = 'interp';
        gas_slice.EdgeColor = 'none';
        gas_slice.AmbientStrength = 1;
        
        % Create arrows to indicate outflow of gas
        for i_arrow = 1:nGasArrow
            % Delete any old arrows
            if isfield(gasArrows(i_arrow), 'Surface')
                delete(gasArrows(i_arrow).Surface);
            end
            
            %Find new start and end position of arrows
            if i_step > gas_step_init + nGasArrowGrowSteps
                gasArrows(i_arrow).X_min = gasArrows(i_arrow).X_min +...
                    dGasArrowLength * gasArrows(i_arrow).V_hat; 
            end
            gasArrows(i_arrow).X_max = gasArrows(i_arrow).X_max +...
                    dGasArrowLength * gasArrows(i_arrow).V_hat; 
            
            % Plot New arrows
            gasArrows(i_arrow).Surface = mArrow3(gasArrows(i_arrow).X_min,...
                gasArrows(i_arrow).X_max, 'tipLength', gasArrTipLength, 'tipAngle',...
                gasArrTipAngle, 'FaceColor', gasArrowColor);
            set(gasArrows(i_arrow).Surface, 'FaceLighting', 'none');
            if i_step > GasArrowFadeStart
                gasArrowAlpha = max(0, 1-(i_step - GasArrowFadeStart)/nGasArrowFadeSteps);
                gasArrows(i_arrow).Surface.FaceAlpha = gasArrowAlpha;
            end
        end
        
        
    end
    
    % -------------------------------------------------------------------
    % ------------------- Create pls. in the solar wind -----------------
    % -------------------------------------------------------------------
    
    for i_sw = 1:nSWPls
        delete(SWPls(i_sw).Surface);
        SWPls(i_sw).X_centre(1) = interp1(B_lines(ix_SW_Bline).Data(:,1), B_lines(ix_SW_Bline).Data(:,2),...
            y_SW(i_sw));
        SWPls(i_sw).X_centre(2) = y_SW(i_sw);
        SWPls(i_sw).X_centre(3) = 0;

        SWPls(i_sw).X = zeros([size(x_pl),3]);

        SWPls(i_sw).X(:,:,1) = SWPls(i_sw).X_centre(1) + x_pl ;%+ 200*(2* rand() -1);
        SWPls(i_sw).X(:,:,2) = SWPls(i_sw).X_centre(2) + y_pl;
        SWPls(i_sw).X(:,:,3) = SWPls(i_sw).X_centre(3) + z_pl;
        
        SWPls(i_sw).X(:,:,1) = SWPls(i_sw).X(:,:,1) + dx_sw;
       SWPls(i_sw).Surface = surf(SWPls(i_sw).X(:,:,1), SWPls(i_sw).X(:,:,2),...
        SWPls(i_sw).X(:,:,3), 'FaceColor', col_SW_pls, 'Edgecolor','none');
        set(SWPls(i_sw).Surface, 'facelighting', 'gouraud'); 
        
    end
    
    
    %--------------------------------------------------------------------
    %------------- Rotate to foucs on Comet -----------------------------
    %--------------------------------------------------------------------
%     if i_step == 40
%         % Initial view angle before pan
%         [az_pan_init, el_pan_init] = view();
%     end
%     if and(i_step > 40, i_step<50)
%         
%         va = camva();
%         campan(az_pan_init/10, el_pan_init/10); 
%         camva(va + 1);
% %     campan
%     
%     end
    % -------------------------------------------------------------------
    % ------------------- Generate Electrons in the coma ----------------
    % -------------------------------------------------------------------
    
    
%     if i_step > 50
%         camtarget([0,0,0]);
        % focus camera on comet nucleus 
%         camva(50);
%     for i_photon = 1:nPhotons
%         photons(i_photon).X(:,1) = photons(i_photon).X(:,1) + dx_ph;
%         r_photon = sqrt(sum(photons(i_photon).X(end,:).^2,2));
%         if and(r_photon < 3 * r_67P, rand()<0.03)
%             % Delete photon and create electron                
%             delete(photons(i_photon).Line);
%             delete(photons(i_photon).Arrow);
%             nPhotoEls = nPhotoEls + 1;
%             i_photo_el = nPhotoEls;
%             photo_els(i_photo_el).X_centre = photons(i_photon).X(end,:); 
%             photo_ions(i_photo_el).X_centre = photons(i_photon).X(end,:);
%             photo_els(i_photo_el).V_hat = photo_els(i_photo_el).X_centre ./ ...
%                 sqrt(sum(photo_els(i_photo_el).X_centre.^2));
%             photo_ions(i_photo_el).V_hat = photo_ions(i_photo_el).X_centre ./ ...
%                 sqrt(sum(photo_ions(i_photo_el).X_centre.^2));
%             
%             photo_els(i_photo_el).X = zeros([size(x_pl),3]);
%             photo_els(i_photo_el).X(:,:,1) = photo_els(i_photo_el).X_centre(1) + x_pl;
%             photo_els(i_photo_el).X(:,:,2) = photo_els(i_photo_el).X_centre(2) + y_pl;
%             photo_els(i_photo_el).X(:,:,3) = photo_els(i_photo_el).X_centre(3) + z_pl;
% 
%             photo_ions(i_photo_el).X = zeros([size(x_pl),3]);
%             photo_ions(i_photo_el).X(:,:,1) = photo_ions(i_photo_el).X_centre(1) + 2 * x_pl;
%             photo_ions(i_photo_el).X(:,:,2) = photo_ions(i_photo_el).X_centre(2) + 2 * y_pl;
%             photo_ions(i_photo_el).X(:,:,3) = photo_ions(i_photo_el).X_centre(3) + 2 * z_pl;
% %                 if rand() >0.5
%             photo_els(i_photo_el).Surface = surf(photo_els(i_photo_el).X(:,:,1),...
%                 photo_els(i_photo_el).X(:,:,2), photo_els(i_photo_el).X(:,:,3),...
%                 'FaceColor', col_ph_el, 'Edgecolor','none');
%                 set(photo_els(i_photo_el).Surface, 'facelighting', 'gouraud'); 
%             photo_ions(i_photo_el).Surface = surf(photo_ions(i_photo_el).X(:,:,1),...
%                 photo_ions(i_photo_el).X(:,:,2), photo_ions(i_photo_el).X(:,:,3),...
%                 'FaceColor', col_ph_el, 'Edgecolor','none');
%                 set(photo_ions(i_photo_el).Surface, 'facelighting', 'gouraud');
% 
%             photons(i_photon) = photons(end);
% 
%             
%             nPhotons = nPhotons - 1;
% %             i_photon = i_photon - 1;
%             
% 
%         end
%     end
%     end
    photons = photons(1:nPhotons);
%         end
    if i_step > StepInitPhEl 
        
    for i_photo_el = 1:nPhotoEls
        if i_step > StepInitPhEl + 1
            delete(photo_els(i_photo_el).Surface);
        end
%         delete(photo_ions(i_photo_el).Surface);
        photo_els(i_photo_el).Surface = surf(photo_els(i_photo_el).X(:,:,1),...
            photo_els(i_photo_el).X(:,:,2), photo_els(i_photo_el).X(:,:,3),...
            'FaceColor', col_ph_el, 'Edgecolor','none');
            set(photo_els(i_photo_el).Surface, 'facelighting', 'gouraud');
%         photo_ions(i_photo_el).Surface = surf(photo_ions(i_photo_el).X(:,:,1),...
%             photo_ions(i_photo_el).X(:,:,2), photo_ions(i_photo_el).X(:,:,3),...
%             'FaceColor', col_ph_ion, 'Edgecolor','none');
%             set(photo_ions(i_photo_el).Surface, 'facelighting', 'gouraud');

%         photo_ions(i_photo_el).X = photo_ions(i_photo_el).X +...
%             reshape(photo_ions(i_photo_el).V_hat * dx_ph_ion,1,1,3);
        photo_els(i_photo_el).X = photo_els(i_photo_el).X +...
            reshape(photo_els(i_photo_el).V_hat * dx_ph_el,1,1,3);


    end
    end
    
    drawnow;
    
    frame = getframe(f);
    writeVideo(introvideo, frame);
    % ------------- Program Camera Movement ----------------------------
%     if i_step <30
%         camtarget(photons(floor(nPhotons/2)).X(end,:));
% %         if i_step == 1
% %             camva(70);
% %         end
%         campos([0,-500,0]);
% %         camzoom(min(3, 1 + 0.1* i_step))
%         if i_step > 5
% %             camzoom(1.05);
%             va = camva();
%             camva(va*0.96);
%         end
% %     else
% %         camlookat(photons(floor(nPhotons/2)).Line);
%     end
%     drawnow;
    
%     Plot neutral gas
end

close(introvideo);


%%

% ----------------------------------------------------------------
% -------------- Create a heatwave -------------------------------
% ----------------------------------------------------------------
l_heat = 100;
[x_heat, y_heat] = meshgrid(linspace(xlims(1), xlims(1) + l_heat), linspace(ylims(1), ylims(2)));
z_heat = zeros(size(x_heat));

c_heat = repmat(permute([255, 94,19]/255,[1,3,2]), size(x_heat));
% c_heat = permute
% Opacity of heat wave 
x_max = x_heat(1,end);
x_min = x_heat(1, 1);

alpha_heat = sigmoid(x_heat, (x_min + x_max) /2,1/20);

figure;
% plot(x_heat(1,:), alpha_heat);
heat_wave = surf(x_heat, y_heat, z_heat, c_heat, 'EdgeColor', 'none');
cb = colorbar('Location', 'west');
heat_wave.AlphaData = alpha_heat;
heat_wave.FaceColor = 'interp';
heat_wave.FaceAlpha = 'interp';
heat_wave.EdgeColor = 'none';
view(2)