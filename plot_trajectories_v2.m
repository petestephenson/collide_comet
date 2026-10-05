%% -------- Checking for OS ----------------------
OS_name = system_dependent('getos');
i_windows = strcmp('Microsoft',OS_name(1:9));
if i_windows
    one_drive = 'C:/Users/pete4/OneDrive - Imperial College London/';
else
    one_drive = '/Users/pps18/OneDrive - Imperial College London/';
end



%% Open E and B fields
load('E_ave.mat');
load('B_ave.mat');

nFieldX = size(E_ave, 1);
nFieldY = size(E_ave, 2);
nFieldZ = size(E_ave, 3);
%%
i_pl = 1;
pl_no = [1];
pl_traj = csvread([one_drive,'Fortran/Energy_Degradation/'...
        'testing123/ioni/pl_traj_ioni',...
        num2str(pl_no(i_pl))], 0,0);

plot_range = pl_traj(2:end,1)~=0;
plot_range = [plot_range;false];

X = pl_traj(plot_range,1:3);
V = pl_traj(plot_range,4:6);

pls = analyse_traj(X,V, E_ave, B_ave);
%%




gradB = zeros(size(B_ave));
% for i_dim = 1:3
[gradB(:,:,:,1),gradB(:,:,:,2),gradB(:,:,:,3)] = ...
    gradient(sqrt(sum(B_ave.^2,4)));
% end
gradB = gradB / 1e4; %Converts from 10km to 1m
%%
qe = 1.6e-19; me = 9.11e-31;

B_temp = permute(B_ave, [2,1,3,4]);
E_temp = permute(E_ave, [2,1,3,4]);
gradB_temp = permute(gradB, [2,1,3,4]);
x_comet = [134, 97, 97];
pl_no = [1];
nPls = numel(pl_no);
% pl_traj = cell(nPls,1);
pls = struct();
%%
for i_pl = 1:nPls
    pl_traj = csvread([one_drive,'Fortran/Energy_Degradation/'...
        'testing123/ioni/pl_traj_ioni',...
        num2str(pl_no(i_pl))], 0,0);
%     pl_traj = csvread([one_drive,'Fortran/Energy_Degradation/'...
%         'testing123/cold/pl_traj_cold',...
%         num2str(pl_no(i_pl))], 0,0);
    plot_range = pl_traj(2:end,1)~=0;
    plot_range = [plot_range;false];
% pls = cell(2,1);
   X = pl_traj(plot_range,1:3);
   V = pl_traj(plot_range,4:6);
   
    pls(i_pl) = analyse_traj(X,V, E_ave, B_ave);

   pls(i_pl).X = pl_traj(plot_range,1:3);
   pls(i_pl).V = pl_traj(plot_range,4:6);
   
   
   pls(i_pl).dt = zeros(size(pls(i_pl).X,1),1);
   pls(i_pl).dt(1:end-1) = 1e4 * (pls(i_pl).X(2:end,1) - pls(i_pl).X(1:end-1,1))./pls(i_pl).V(1:end-1,1);
   pls(i_pl).T = [0;cumsum(pls(i_pl).dt(1:end-1))]; 
   pls(i_pl).Eng = 0.5 * 9.11e-31 * sum(pls(i_pl).V.^2,2)/1.6e-19;
   pls(i_pl).E_pl = zeros(size(pls(i_pl).X));
   pls(i_pl).B_pl = zeros(size(pls(i_pl).X));
   pls(i_pl).gradB = zeros(size(pls(i_pl).X));
   pls(i_pl).X_gyro = zeros(size(pls(i_pl).X));
   pls(i_pl).V_gyro = zeros(size(pls(i_pl).X));
   for i_dim = 1:3
        pls(i_pl).E_pl(:, i_dim) = interp3(...
            E_temp(:,:,:,i_dim),pls(i_pl).X(:,1), pls(i_pl).X(:,2),pls(i_pl).X(:,3));
        pls(i_pl).B_pl(:, i_dim) = interp3(...
            B_temp(:,:,:,i_dim),pls(i_pl).X(:,1), pls(i_pl).X(:,2),pls(i_pl).X(:,3));
        pls(i_pl).gradB(:,i_dim) = interp3(...
            gradB_temp(:,:,:,i_dim), pls(i_pl).X(:,1), pls(i_pl).X(:,2),pls(i_pl).X(:,3));
        [envHigh, envLow] = envelope(pls(i_pl).V(:,i_dim),10,'peak');
        pls(i_pl).V_gyro(:,i_dim)  = (envHigh+envLow)/2;
        [envHigh, envLow] = envelope(pls(i_pl).X(:,i_dim),10,'peak');
        pls(i_pl).X_gyro(:,i_dim)  = (envHigh+envLow)/2;
        
   end

   pls(i_pl).b_hat = pls(i_pl).B_pl...
    ./ sqrt(sum(pls(i_pl).B_pl.^2,2));
   pls(i_pl).e_hat = pls(i_pl).E_pl...
    ./ sqrt(sum(pls(i_pl).E_pl.^2,2));
%    pls(i_pl).X = (pls(i_pl).X - x_comet) * 10;
%    pls(i_pl).gyroX = movmean(pls(i_pl).X,20);
    
   
   pls(i_pl).V_para_noisy = dot(pls(i_pl).V,  pls(i_pl).b_hat,2) .* pls(i_pl).b_hat ;
   pls(i_pl).V_perp_noisy = pls(i_pl).V - pls(i_pl).V_para_noisy;

   pls(i_pl).Eng_para = me * sum(pls(i_pl).V_para_noisy.^2,2)/(2*qe);
   pls(i_pl).Eng_perp = me * sum(pls(i_pl).V_perp_noisy.^2,2)/(2*qe);
   
   [envHigh, envLow] = envelope(pls(i_pl).Eng_perp,10,'peak');
    pls(i_pl).Eng_perp_gyro =  (envHigh+envLow)/2;
   [envHigh, envLow] = envelope(pls(i_pl).Eng,10,'peak');
    pls(i_pl).Eng_gyro =  (envHigh+envLow)/2;
   
   pls(i_pl).V_para = dot(pls(i_pl).V_gyro,  pls(i_pl).b_hat,2) .* pls(i_pl).b_hat ;
   pls(i_pl).V_perp = pls(i_pl).V_gyro - pls(i_pl).V_para;
end
%% Retrieve B field lines for each point in the trajectory

x_grid = [1:nFieldX];
y_grid = [1:nFieldY];
z_grid = [1:nFieldZ];

[meshX, meshY, meshZ] = meshgrid(y_grid, x_grid, z_grid);
stepsize = 0.2;


i_pl = 1;
%%
hs = stream3(meshX,meshY,meshZ,B_ave(:,:,:,2),B_ave(:,:,:,1), B_ave(:,:,:,3),...
   pls(i_pl).X(:,2), pls(i_pl).X(:,1), pls(i_pl).X(:,3), stepsize);
hn = stream3(meshX,meshY,meshZ,-B_ave(:,:,:,2),-B_ave(:,:,:,1), -B_ave(:,:,:,3),...
   pls(i_pl).X(:,2), pls(i_pl).X(:,1), pls(i_pl).X(:,3), stepsize);

%%
% pls(i_pl).b_hat = pls(i_pl).B_pl...
%     ./ sqrt(sum(pls(i_pl).B_pl.^2,2));
field_line = cell(size(hn));
nT = numel(hn);
%%
% L = cell(nT,1);
% R = cell(nT,1);
% k = cell(nT,1);
%%

pls(i_pl).k_Rc = zeros(nT,3);
pls(i_pl).Rc = zeros(nT,1);
AvePts = 4;
trange = find(~cellfun(@isempty, hn));

for i_t = trange
    
    field_line{i_t} = [flipud(hn{i_t});hs{i_t}(2:end,:)];
    pl_loc = size(hn{i_t},1);
    low_bound = max([1,pl_loc-AvePts]);
    up_bound = min([size(field_line{i_t},1),pl_loc+AvePts]);
    
    [L,R,k] = curvature(field_line{i_t}(low_bound:up_bound,[2,1,3]));
    % central points form a straight line as driven by same vector B
    % Average curvature from points nearby
    k = k ./sqrt(sum(k.^2,2));
    pls(i_pl).k_Rc(i_t,:) = mean(k(~isnan(R),:),1);

    pls(i_pl).Rc(i_t) = mean(R(~isnan(R)),1) * 1e4;
end
pls(i_pl).k_Rc = pls(i_pl).k_Rc ./ ...
    sqrt(sum(pls(i_pl).k_Rc.^2,2));
%% Compute the Curvature drift

pls(i_pl).GyroFreq = sqrt(sum(pls(i_pl).B_pl.^2,2)) * qe /me;
pls(i_pl).v_curv = cross(pls(i_pl).k_Rc, pls(i_pl).b_hat,2);
pls(i_pl).v_curv = pls(i_pl).v_curv .* sum(pls(i_pl).V_para.^2,2);
pls(i_pl).v_curv = pls(i_pl).v_curv ./ ...
    (pls(i_pl).GyroFreq .* pls(i_pl).Rc);

pls(i_pl).v_gradB = cross(pls(i_pl).b_hat, pls(i_pl).gradB, 2);
pls(i_pl).v_gradB = pls(i_pl).v_gradB .* ...
    sum(pls(i_pl).V_perp.^2,2)./ ...
    (2 * pls(i_pl).GyroFreq .* sqrt(sum(pls(i_pl).B_pl.^2,2)));
pls(i_pl).v_ExB = cross(pls(i_pl).E_pl, pls(i_pl).B_pl,2) ./...
    sum(pls(i_pl).B_pl.^2,2);
%% Rotate B field along the x-axis
pls(i_pl).rot_theta = - atan2(pls(i_pl).b_hat(:,3),...
        pls(i_pl).b_hat(:,2)) * 180 / pi;
pls(i_pl).rot_a = - atan2(...
    cos(pls(i_pl).rot_theta * pi / 180).* pls(i_pl).b_hat(:,2) -...
    sin(pls(i_pl).rot_theta * pi /180).* pls(i_pl).b_hat(:,3),...
    pls(i_pl).b_hat(:,1)) * 180 / pi;
%%
pls(i_pl).Rot_mat = zeros(3,3,nT);
pls(i_pl).v_para_rot = zeros(nT,3);
pls(i_pl).v_perp_rot = zeros(nT,3);
pls(i_pl).v_circ_rot = zeros(nT,3);
pls(i_pl).E_rot = zeros(nT,3);

for i_t = 1:nT
    pls(i_pl).Rot_mat(:,:,i_t) = ...
        (rotz(pls(i_pl).rot_a(i_t)) * rotx(pls(i_pl).rot_theta(i_t)));
    pls(i_pl).E_rot(i_t,:) = pls(i_pl).Rot_mat(:,:,i_t) *...
        pls(i_pl).E_pl(i_t,:)';
end
pls(i_pl).rot_E_theta = - atan2(pls(i_pl).E_rot(:,3),...
    pls(i_pl).E_rot(:,2)) * 180 / pi;

pls(i_pl).v_circ = pls(i_pl).V - pls(i_pl).V_gyro;

for i_t = 1:nT
    pls(i_pl).Rot_mat(:,:,i_t) = rotx(pls(i_pl).rot_E_theta(i_t)) *...
        (rotz(pls(i_pl).rot_a(i_t)) * rotx(pls(i_pl).rot_theta(i_t)));
    pls(i_pl).V_gyro_rot(i_t,:) = pls(i_pl).Rot_mat(:,:,i_t) *...
        pls(i_pl).V_gyro(i_t,:)';
    pls(i_pl).v_para_rot(i_t,:) = pls(i_pl).Rot_mat(:,:,i_t) *...
        pls(i_pl).V_para(i_t,:)';
    pls(i_pl).v_perp_rot(i_t,:) = pls(i_pl).Rot_mat(:,:,i_t) *...
        pls(i_pl).V_perp(i_t,:)'; 
    pls(i_pl).v_ExB_rot(i_t,:) = pls(i_pl).Rot_mat(:,:,i_t) *...
        pls(i_pl).v_ExB(i_t,:)';
    pls(i_pl).v_curv_rot(i_t,:) = pls(i_pl).Rot_mat(:,:,i_t) *...
        pls(i_pl).v_curv(i_t,:)';
    pls(i_pl).v_gradB_rot(i_t,:) = pls(i_pl).Rot_mat(:,:,i_t) *...
        pls(i_pl).v_gradB(i_t,:)';
    pls(i_pl).v_circ_rot(i_t,:) = pls(i_pl).Rot_mat(:,:,i_t) *...
        pls(i_pl).v_circ(i_t,:)';
end
pls(i_pl).MagMom = pls(i_pl).Eng_perp_gyro ./sqrt(sum(pls(i_pl).B_pl.^2,2));
pls(i_pl).PitchAngle = atan2(sqrt(sum(pls(i_pl).V_perp_noisy.^2,2)),...
        sqrt(sum(pls(i_pl).V_para_noisy.^2,2))) * 180 / pi;
    
%% Calculate potential change
pls(i_pl).dX = (pls(i_pl).X(2:end,:) - pls(i_pl).X(1:end-1,:)) * 1e4;% in m
pls(i_pl).V_pot_para = zeros(nT,1);
pls(i_pl).V_pot_curv = zeros(nT,1);
for i_t = 1:nT-1
    pls(i_pl).V_pot_para(i_t + 1) = pls(i_pl).V_pot_para(i_t) + ...
        dot(pls(i_pl).b_hat(i_t,:), pls(i_pl).dX(i_t,:)) .*...
        pls(i_pl).E_rot(i_t,1);
    pls(i_pl).V_pot_curv(i_t + 1) = pls(i_pl).V_pot_curv(i_t) + ...
        dot(pls(i_pl).E_pl(i_t,:), pls(i_pl).v_curv(i_t,:)) .* ...
        pls(i_pl).dt(i_t);
end


%%
ix_fix = 2906;

delt_para = pls(i_pl).Eng_para - (pls(i_pl).Eng_para(ix_fix)...
    + pls(i_pl).V_pot_para(ix_fix) - pls(i_pl).V_pot_para);
delt_perp = pls(i_pl).Eng_perp_gyro - (pls(i_pl).Eng_perp_gyro(ix_fix)...
    + pls(i_pl).V_pot_curv(ix_fix) - pls(i_pl).V_pot_curv);
figure;
yyaxis left
plot(pls(i_pl).T, delt_para);
hold on
ylim([-20,5]);
yyaxis right 
plot(pls(i_pl).T, pls(i_pl).Eng_perp_gyro);
hold on
ylim([0,25]);
plot(pls(i_pl).T,pls(i_pl).Eng_perp_gyro(ix_fix)...
    + pls(i_pl).V_pot_curv(ix_fix) - pls(i_pl).V_pot_curv, 'DisplayName', 'Pot_{curv}');
%%

% mag_pressure = 
plot_array = pls(i_pl).Eng_gyro;
figure;
subplot(2,1,1)
% plot(pls(i_pl).T,pls(i_pl).Eng_para(ix_fix)...
%     + pls(i_pl).V_pot_para(ix_fix) - pls(i_pl).V_pot_para, 'DisplayName', 'Pot_{para}');
hold on
% plot(pls(i_pl).T, plot_array, 'DisplayName', 'Eng_{gyro}');
% plot(pls(i_pl).T,pls(i_pl).Eng, 'DisplayName', 'Energy');
plot(pls(i_pl).T,pls(i_pl).Eng_gyro(ix_fix)...
     + pls(i_pl).V_pot_curv(ix_fix) - pls(i_pl).V_pot_curv...
     + pls(i_pl).V_pot_para(ix_fix) - pls(i_pl).V_pot_para,...
     'DisplayName', 'Pot_{tot}');
plot(pls(i_pl).T, pls(i_pl).Eng_para,...
    'DisplayName', 'Eng_{para}');
% plot(pls(i_pl).T, pls(i_pl).Eng_perp_gyro,...
%     'DisplayName', 'Eng_{perp}');
legend;
subplot(2,1,2);
% yyaxis left
% plot(pls(i_pl).T,sum(pls(i_pl).B_pl.^2,2));
% yyaxis right
% plot(pls(i_pl).T,delt_para, 'DisplayName','D Para');
hold on
% plot(pls(i_pl).T,delt_perp, 'DisplayName', 'D Perp');
plot(pls(i_pl).T,delt_perp + delt_para, 'DisplayName', 'D Perp');
legend;


%% Plotting lab velocity and drifts
dim_name = {'x','y','z'};
figure;
for i_dim = 1:3
    subplot(3,1,i_dim);
    plot(pls(i_pl).T,pls(i_pl).V_gyro(:,i_dim), 'DisplayName', 'V_{tot}');
    hold on
    plot(pls(i_pl).T,pls(i_pl).V_perp(:,i_dim), 'DisplayName', 'v_{perp}');
    
    plot(pls(i_pl).T,pls(i_pl).V_para(:,i_dim), 'DisplayName', 'v_{para}');
    plot(pls(i_pl).T,pls(i_pl).v_ExB(:,i_dim), 'DisplayName', 'ExB');
%     hold on
    plot(pls(i_pl).T,pls(i_pl).v_curv(:,i_dim), 'DisplayName', 'Curv');
    plot(pls(i_pl).T,pls(i_pl).v_gradB(:,i_dim), 'DisplayName', 'gradB');
    ylabel(['V_', dim_name{i_dim}, ' [ms^{-1}]']);
end
xlabel('T [s]');
legend;

%% Plotting rotated velocity and drifts
dim_name = {'x','y','z'};
figure;
for i_dim = 1:3
    subplot(3,1,i_dim);
    plot(pls(i_pl).T,pls(i_pl).V_gyro_rot(:,i_dim) - ...
        pls(i_pl).v_ExB_rot(:,i_dim), 'DisplayName', 'V_{tot} - v_{ExB}');
    hold on
%     plot(pls(i_pl).T,pls(i_pl).v_perp_rot(:,i_dim), 'DisplayName', 'v_{perp}');
%     hold on
%     plot(pls(i_pl).T,pls(i_pl).v_para_rot(:,i_dim), 'DisplayName', 'v_{para}');
%     plot(pls(i_pl).T,pls(i_pl).v_ExB_rot(:,i_dim), 'DisplayName', 'ExB');
%     hold on
    plot(pls(i_pl).T,pls(i_pl).v_curv_rot(:,i_dim), 'DisplayName', 'Curv');
    plot(pls(i_pl).T,pls(i_pl).v_gradB_rot(:,i_dim), 'DisplayName', 'gradB');
   
    ylabel(['V_', dim_name{i_dim}, ' [ms^{-1}]']);
end
xlabel('T [s]');
legend;
%% Variation of ExB energy
figure;
yyaxis left
plot(pls(i_pl).T, sqrt(sum(pls(i_pl).v_circ_rot.^2,2)),...
     'DisplayName', 'Circular');
 yyaxis right
 plot(pls(i_pl).T,sqrt(sum(pls(i_pl).B_pl.^2,2)), 'DisplayName', 'Bfield');
ylabel('B Field [T]');
% plot(pls(i_pl).T, pls(i_pl).v_circ_rot(:,i_dim)
%%
pls(i_pl).Eng_ExB = me * sum(pls(i_pl).v_ExB.^2,2)/(2*qe);

figure;
hold on
plot(pls(i_pl).T,pls(i_pl).Eng_ExB, 'DisplayName', 'ExB Energy');
plot(pls(i_pl).T,pls(i_pl).Eng_perp, 'DisplayName', 'Perp Energy');

legend;

%% 
figure;
hold on 
subplot(2,1,1)
plot(pls(i_pl).T, pls(i_pl).MagMom);
ylabel('MagMoment [eV/T]');
subplot(2,1,2)
plot(pls(i_pl).T, pls(i_pl).PitchAngle);
xlabel('T [s]');
ylabel('Pitch Angle');



%%
figure;
yyaxis left
plot(pls(i_pl).T,pls(i_pl).GyroFreq, 'DisplayName', 'GyroFreq');
hold on
ylabel('Gyrofreq [s^{-1}]');
yyaxis right
plot(pls(i_pl).T,sqrt(sum(pls(i_pl).B_pl.^2,2)), 'DisplayName', 'Bfield');
ylabel('B Field [T]');
legend;







%% Gyroaverage of perpendicular velocity
fs = 1e4;
pls(i_pl).tvperp = timeseries(pls(i_pl).V_perp, pls(i_pl).T); 
timevec = [0:1/fs: floor(pls(i_pl).T(end)*1e4)/1e4];
L = numel(timevec);
if mod(L,2) == 1
    timevec = timevec(1:end-1);
end
pls(i_pl).tvperp = resample(pls(i_pl).tvperp, timevec);

pls(i_pl).V_gyro = zeros(numel(timevec),3);
for i_dim = 1:3
pls(i_pl).V_gyro(:,i_dim) = lowpass(pls(i_pl).tvperp.Data(:,i_dim),...
    1e-2,1e4);
end

%% Trying to fft
L = numel(timevec);
f = fs * (0:(L/2))/L;
Y = fft(pls(i_pl).tvperp.Data(:,2));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
figure;
plot(f,P1);
%%

figure;
plot(pls(i_pl).T, pls(i_pl).V(:,1));
hold on
plot(pls(i_pl).T, v_mean);


%%
figure; 
plot(timevec, pls(i_pl).V_gyro);
figure;
pspectrum(pls(i_pl).V_gyro(:,i_dim),1e4);
hold on 
pspectrum(pls(i_pl).tvperp.Data(:, i_dim), 1e4);




















%% Interpolate position for E and B fields


i_pl = 1; 

field_pl_x = x_comet + (pls(i_pl).X/10);
nT = size(pls(i_pl).X,1);

        % for i_t = 1:nT
%     for i_dim = 1:3
%         pls(i_pl).E_pl(i_t,i_dim) = lininterp3([1:nFieldX],[1:nFieldY], [1:nFieldZ],...
%             E_ave(:,:,:,i_dim), field_pl_x(i_t,i_dim), field_pl_x(i_t, i_dim),field_pl_x(i_t,i_dim));
%         pls(i_pl).B_pl(i_t,i_dim) = lininterp3([1:nFieldX],[1:nFieldY], [1:nFieldZ],...
%             B_ave(:,:,:,i_dim), field_pl_x(i_t,i_dim), field_pl_x(i_t, i_dim),field_pl_x(i_t,i_dim));
%     end
% end

%% ---------- Calculate Trajectory Properties -----------------------------
figure; 
plot(pls(1).T, pls(1).X);
hold on
plot(pls(1).T, movmean(pls(1).X,20));

%% Plot Pl energy
figure; 
for i_pl = 1:nPls
   plot(pls(i_pl).Eng, 'displayname', num2str(i_pl));
    hold on
    
end
legend;
%% ---------- Trajectory Plots --------------------------------------------
x_comet = [134, 97, 97];
r_com = 1.7;


i_tot_eng = 1;
i_time = 0;
if any([i_tot_eng, i_time])
    i_colbar = 1;
else
    i_colbar = 0;
    
end

if i_tot_eng
    col_str = 'Eng';
    clabel_name = 'Total Energy [eV]';
elseif i_time
    col_str = 'T';
    clabel_name = 'Time [s]';
end


ax = cell(nPls,1);
figure('units', 'centimeters', 'position', [2, 2, 15, 15]);
colororder = get(gca, 'ColorOrder');
for i_pl = 14:14 %1:nPls
    if i_colbar
        coltemp = getfield(pls(i_pl), col_str);
        color_array = [coltemp';coltemp'];
        ax{i_pl} = surface( [pls(i_pl).X(:,1)';pls(i_pl).X(:,1)'],...
        [pls(i_pl).X(:,2)';pls(i_pl).X(:,2)'],...
        [pls(i_pl).X(:,3)'; pls(i_pl).X(:,3)'],...
        color_array, 'facecolor','none','edgecolor','interp','linewidth',0.9);
%         ax{i_pl} = surface( [pls(i_pl).gyroX(:,1)';pls(i_pl).gyroX(:,1)'],...
%         [pls(i_pl).gyroX(:,2)';pls(i_pl).gyroX(:,2)'],...
%         [pls(i_pl).gyroX(:,3)'; pls(i_pl).gyroX(:,3)'],...
%         color_array, 'facecolor','none','edgecolor','interp','linewidth',0.9,...
%         'DisplayName', 'Pl');
    else
        plot3(pls(i_pl).gyroX(:,1),pls(i_pl).gyroX(:,2),pls(i_pl).gyroX(:,3),...
            'LineWidth',2, 'DisplayName',['Generation',num2str(i_pl)]);
    end
    
    hold on
end
if i_colbar
    c = colorbar;
    colormap(jet);
    c.Label.String= clabel_name;
    c.Label.FontSize = 15;
%     caxis([0,15]);
end
 legend;
% caxis(col_lims);
[x_67P,y_67P,z_67P] = sphere;
x_67P = (x_67P *r_com);% Should be *1.7/10 inside brackets
y_67P = (y_67P *r_com);
z_67P = (z_67P *r_com);
surf(x_67P, y_67P, z_67P, 'FaceColor','k','FaceAlpha', 0.5, 'DisplayName', 'Nucleus')
daspect([1,1,1]);
set(gca, 'Fontsize', 15);
hold on

xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');

grid on 

xlim([-100,100]);
ylim([-100,100]);
zlim([-100,100]);
view([1,0,0]);
% plot_range = pl1(2:end,1)~=0;
% plot_range = [plot_range;false];
% % pls = cell(2,1);
% 
% 
% pls(1).x_pl = pl_traj(plot_range,1:3);
% pls(1).v_pl = pl_traj(plot_range,4:6);

% pls{1}


%%
x_comet = [134, 97, 97];
r_com = 1.7;

plot_range = t_pl(:,1) ~=0;

%% Plotting time as a function of step
% t_plot = t_pl(plot_range);
% figure('units', 'centimeters', 'position', [2, 2, 15, 15]);
% yyaxis left
% plot(t_plot(2:end)-t_plot(1:end-1),'DisplayName','Time Interval [s]');
% ylabel('Simulation Time [s]')
% yyaxis right
% plot(1./sqrt(sum(B_pl.^2,2)), 'DisplayName','Magnitude of B');
% ylabel('1/Magnitude of B [nT^{-1}]')
% set(gca, 'Fontsize', 15);
% xlabel('Timestep')
% 
%%

figure; 
for i_pl = 1:1
    avedt = movmean(pls(i_pl).dt,20);
    plot(movmean(pls(i_pl).T,20), movmean(pls(i_pl).X(:,1),20));
end

%%
fs = 1;
ts = timeseries(pls(i_pl).V(:,1),pls(i_pl).T);
T_max = pls(i_pl).T(end);
tvec = linspace(0,T_max, 100001);
ty = resample(ts, tvec);

figure;
% plot(tvec, test)
hold on
plot(tvec, movmean(ty.Data,10));
%% Calc ave velocity
i_pl = 15;
ave_vx = (pls(i_pl).X(end,1) - pls(i_pl).X(1,1))/(pls(i_pl).T(end)-pls(i_pl).T(1))
%%
%%
v_sw = 400e3;
figure;

max_steps = 0;

for i_pl = 1:1%:nPls
    
    pls(i_pl).gyroV = zeros(size(pls(i_pl).V));
    for i_dim = 1:3
        pls(i_pl).gyroV(:, i_dim) = movmean(pls(i_pl).V(:, i_dim), 20);
    end
    hold on;
    max_steps = max(size(pls(i_pl).gyroV,1), max_steps);
    plot(pls(i_pl).gyroV(:,1));
end
plot(v_sw*ones(max_steps,1),'k', 'LineWidth', 2, 'linestyle', ':');
%% ---------------- Plot trajectory -----------------------------

i_tot_eng = 1;
i_eng_para = 0;
i_eng_perp = 0;
i_pot = 0;
i_time = 0;
if i_tot_eng
    color_array = [Eng_pl(plot_range)'; Eng_pl(plot_range)'];
    clabel_name = 'Total Energy [eV]';
elseif i_eng_para
    color_array = [Eng_para(plot_range)';Eng_para(plot_range)'];
    clabel_name = 'Parallel Energy [eV]';
elseif i_eng_perp 
    color_array = [Eng_perp(plot_range)';Eng_perp(plot_range)'];
    clabel_name = 'Perp Energy [eV]';
elseif i_pot
    color_array = [V_pot(plot_range)';V_pot(plot_range)'];
    clabel_name = 'Potential Energy [eV]';
elseif i_time
    color_array = [t_pl(plot_range)';t_pl(plot_range)'];
    clabel_name = 'Time [s]';
end
E_max = max(color_array(1,:));
E_min = min(color_array(1,:));

x_plot = (x_pl(plot_range,:) - x_comet) * 10;

figure('units', 'centimeters', 'position', [2, 2, 15, 15]);
ax0 = surface( [x_plot(:,1)';x_plot(:,1)'],...
    [x_plot(:,2)';x_plot(:,2)'],...
    [x_plot(:,3)'; x_plot(:,3)'],...
    color_array, 'facecolor','none','edgecolor','interp','linewidth',0.9);
c = colorbar;
c.Label.String= clabel_name;
c.Label.FontSize = 15;
caxis([E_min, E_max]);
[x_67P,y_67P,z_67P] = sphere;
x_67P = (x_67P *r_com);% Should be *1.7/10 inside brackets
y_67P = (y_67P *r_com);
z_67P = (z_67P *r_com);
hold on
% surf(x_67P, y_67P, z_67P, 'FaceColor','k','FaceAlpha', 0.5)
set(gca, 'fontsize', 15);
% xlim([min(x_pl(:,1),[],1), max(x_pl(:,1),[],1)]);
% ylim([min(x_pl(:,2),[],1), max(x_pl(:,2),[],1)]);
% zlim([180, 181]);
view([1,0,0]);
% caxis([0,25]);
% hold on
% [x_67P,y_67P,z_67P] = sphere;
% x_67P = (x_67P *1.7/10) + x_comet(1);% Should be *1.7/10 inside brackets
% y_67P = (y_67P *1.7/10) + x_comet(2);
% z_67P = (z_67P *1.7/10) + x_comet(3);
% 
% % surf(x_67P, y_67P, z_67P, 'FaceColor','k','FaceAlpha', 0.5);
% hold on
% i_layer = 90;
% [y,x,z] = meshgrid(i_layer,1:288,1:192);
% E_mag = sqrt(sum(E_ave.^2,4));
% ax1 = quiver3(x,y,z, E_ave(:,i_layer,:,1),E_ave(:,i_layer,:,2),E_ave(:,i_layer,:,3));
% ax2 = surf(reshape(x,288,192), reshape(y,288,192), reshape(z,288,192), reshape( E_mag(:,i_layer,:),288,192) );
% colorbar;
% xlim([-30,30]);
% ylim([-40,30]);
% zlim([-10,50]);

grid on;
% daspect([1,1,1]);
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
% view(3);
%% Drifts
%-------------------- E x B ---------------------------------------------
exb = cross(E_pl, B_pl, 2);
vexb = exb./sum(B_pl.^2,2);
v_n =  dot(v_pl, exb,2)./sqrt(sum(exb.^2,2));
v_n2 = movmean(v_n,40);
figure('units', 'centimeters', 'position', [2, 2, 15, 15]);
plot(t_pl(plot_range), v_n2(plot_range), 'displayname', 'Electron speed in ExB');
hold on
plot(t_pl(plot_range),sqrt(sum(vexb(plot_range,:).^2,2)), 'DisplayName', 'Theoretical ExB drift');
legend;
ylabel('Electron Velocity [ms^{-1}]')
xlabel('Time [s]')
set(gca, 'Fontsize', 15);
max(v_n2);
max(norm(vexb,2))

%---------------------- Grad B ------------------------------------------

%%  -------------------- Plot Energies ----------------------------------
figure; 
% yyaxis left
plot(t_pl(plot_range),Eng_perp(plot_range), 'DisplayName', 'E_{perp}',...
    'linewidth', 2,'linestyle','-','color','r');
hold on
plot(t_pl(plot_range), Eng_para(plot_range), 'DisplayName', 'E_{para}',...
    'linewidth', 2,'linestyle','-','color','g');
plot(t_pl(plot_range), Eng_pl(plot_range), 'DisplayName', 'Total Eng',...
    'linewidth', 2,'linestyle','-','color','k');
% plot(t_pl(plot_range), V_pot(plot_range),...
%     'DisplayName','Field Aligned Potential', 'linewidth', 2,'linestyle','-','color','b');
ylabel('Energy [eV]');
% yyaxis right 
% plot(t_pl(plot_range), coll(plot_range ~= 0),'Displayname', 'Collision','color', 'k', 'linestyle','-');
% plot(sqrt(sum(B_pl(plot_range).^2,2)),'DisplayName', 'Bfield');
legend;
xlabel('Time [s]');
set(gca, 'Fontsize',15);

%%
dE = Eng_pl - Eng_pl(1);
figure; 
plot(movmean(dE(plot_range),20)/Eng_pl(1));

%%  -------------- Plasma Density -----------------
% /Users/pps18/OneDrive - Imperial College London/Energy_Degradation/Runs/1096396
pl_count_in = csvread([one_drive,'Energy_Degradation/Runs/1096396/pl_count.txt']);
plasma_dens_in = csvread([one_drive, 'Energy_Degradation/Runs/1096396/plas_dens.txt']);
bulk_vel_in = csvread([one_drive, 'Energy_Degradation/Runs/1096396/bulk_velocity.txt']);
square_vel_in = csvread([one_drive, 'Energy_Degradation/Runs/1096396/square_velocity.txt']);
%%
nGridE = 40;
nGridX = 70; nGridY = 60; nGridZ = 60;
dx_lin = 1; dy_lin = 1; dz_lin = 1; %km
nLinX = 10; nLinY = 10; nLinZ =10;
nLinE = 10; dE_lin = 0.5; %eV
gridx_max = 1540; gridy_max = 1540; gridz_max = 1540;%km
gridE_max = 400;%eV

%--- plas_dens - actually f(x,y,z,E) for primary and secondary elecs ------
plas_dens_E = reshape(plasma_dens_in, nGridZ, nGridE,2,nGridX,nGridY);
plas_dens_E = permute(plas_dens_E, [4,5,1,2,3]);
%--- pl_count -------------------------------------------------------------
pl_count = reshape(pl_count_in, nGridZ, nGridE,2,nGridX,nGridY);
pl_count = permute(pl_count, [4,5,1,2,3]);
%--- bulk_vel ---------------------------------------------------------
bulk_vel = reshape(bulk_vel_in, nGridY, nGridZ, 2, 3, nGridX);
bulk_vel = permute(bulk_vel, [4,5,1,2,3]);
%--- square_vel ---------------------------------------------------------
square_vel = reshape(square_vel_in, nGridY, nGridZ, 2, 3, nGridX);
square_vel = permute(square_vel, [4,5,1,2,3]);


%% -------------- Define the x, y, z grid spacing -------------------------
x_grid = zeros(nGridX+1,1);
nHalfX = nGridX/2;
y_grid = zeros(nGridY+1,1);
nHalfY = nGridY/2;
z_grid = zeros(nGridZ+1,1);
nHalfZ = nGridZ/2;

gridx_min = dx_lin * nLinX;
for i_x = 1:nLinX
    x_grid(nHalfX +1 + i_x) = dx_lin * i_x;
end
    
for i_x = nLinX + 1:nHalfX    
    x_grid(nHalfX +1 + i_x) = exp(log(gridx_min) + ...
        (log(gridx_max) - log(gridx_min))*(i_x - nLinX)/(nHalfX - nLinX));
end

x_grid(1:nHalfX) = -flipud(x_grid(nHalfX+2:end));

gridy_min = dy_lin * nLinY;
for i_y = 1:nLinY
    y_grid(nHalfY +1 + i_y) = dy_lin * i_y;
end
    
for i_y = nLinY + 1:nHalfY    
    y_grid(nHalfY +1 + i_y) = exp(log(gridy_min) + ...
        (log(gridy_max) - log(gridy_min))*(i_y - nLinY)/(nHalfY - nLinY));
end
y_grid(1:nHalfY) = -flipud(y_grid(nHalfY+2:end));
gridz_min = dz_lin * nLinZ;
for i_x = 1:nLinZ
    z_grid(nHalfZ +1 + i_x) = dz_lin * i_x;
end
    
for i_x = nLinZ + 1:nHalfZ    
    z_grid(nHalfZ +1 + i_x) = exp(log(gridz_min) + ...
        (log(gridz_max) - log(gridz_min))*(i_x - nLinZ)/(nHalfZ - nLinZ));
end
z_grid(1:nHalfZ) = -flipud(z_grid(nHalfZ+2:end));

%% ---------------- Define the Energy grid spacing ------------------------
E_grid = zeros(nGridE+1,1);

for i_E = 1:nLinE
    E_grid(i_E + 1) = i_E * dE_lin;
    
end
gridE_min = dE_lin * nLinE;

for i_E = nLinE + 1:nGridE
    E_grid(i_E + 1) = exp(log(gridE_min) + (log(gridE_max) - log(gridE_min))...
        * ( i_E - nLinE ) / (nGridE - nLinE));
end

%% ----- Volume of grid cells --------------------------------------------

dE_grid = E_grid(2:end) - E_grid(1:end-1);
dx_grid = x_grid(2:end) - x_grid(1:end-1);
dy_grid = y_grid(2:end) - y_grid(1:end-1);
dz_grid = z_grid(2:end) - z_grid(1:end-1);
grid_vol = dx_grid .* reshape(dy_grid, 1, []) .* reshape(dz_grid, 1,1, []);


%%
plas_dens = sum(plas_dens_E .* reshape(dE_grid,1,1,1,[]),[4,5]);
pl_count_wo_E = sum(pl_count,[4,5]);

%% 
E_cold = 1;
ix_cold = find(0.5*(E_grid(1:end-1)+E_grid(2:end))<E_cold);
plas_dens_cold = sum(plas_dens_E(:,:,:,ix_cold,:) ...
    .* reshape(dE_grid(ix_cold),1,1,1,[]),[4,5]);

%%
plas_dens = pl_count_new;%./ grid_vol;
plas_dens = permute(plas_dens,[3,2,1]);
plas_dens = sum(plas_dens,3);


[X,Y] = meshgrid(x_grid);
plas_dens2 = [plas_dens; plas_dens(end,:)];
plas_dens2 = [plas_dens2, plas_dens2(:,end)];
% histogram2('XBinEdges',x_grid,'YBinEdges',x_grid,'BinCounts',plas_dens,'DisplayStyle','tile')
% xlabel('X [km]');
% ylabel('Y [km]');
% % set(gca, 'zscale', 'log');
% set(gca, 'Fontsize', 15);
% colorbar;

f = figure; 
ax1 = surf(X,Y, log10(plas_dens2), 'facecolor', 'flat');
ax1.EdgeColor = 'none';
view(2)
colorbar;
xlabel('X [km]');
ylabel('Y [km]');

% set(gca, 'Xlim', [-200,200]);
% set(gca, 'Ylim', [-100,100]);
% ('BinCounts',test, 'Edges', {x_grid, x_grid});
% f = figure; 
% ax1 = surf(X,Y, log10(plas_dens2), 'facecolor', 'flat');
% ax1.EdgeColor = 'none';
% view(2)
% colorbar;
% xlabel('X [km]');
% ylabel('Y [km]');

%% Slice data
[X,Y,Z] = meshgrid(x_grid(1:end-1), y_grid(1:end-1), z_grid(1:end-1));

i_plasma_dens = 0;
i_plas_dens_cold = 1;
i_pl_count = 0;

if i_plasma_dens
    plot_array = permute(plas_dens,[2,1,3]);
elseif i_plas_dens_cold
    plot_array = permute(plas_dens_cold,[2,1,3]);
elseif i_pl_count
    plot_array = permute(pl_count_wo_E,[2,1,3]);
end

f2 = figure('units','centimeters','position', [1,1,20,15]);

xslice = [0];
yslice = [0];
zslice = [0];
nslice = numel(xslice) + numel(yslice) + numel(zslice);

ax_slice = slice(X,Y,Z,log10(plot_array), xslice,yslice,zslice);
% ax_slice = slice(X,Y,Z,log(plas_dens), xslice,yslice,zslice);
for i_slice = 1:nslice
    ax_slice(i_slice).EdgeColor = 'none';
end
% f2.CurrentAxes.Children(1).Edgecolor = 'none'
% ax_slice = EdgeColor = 'none';
c = colorbar;
colormap(jet);
% caxis([0,Inf]);
c.Label.String = 'log_{10}(N_e [cm^{-3}])';
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
% xlim([-20,20]);
% ylim([-20,20]);
% zlim([-20,20]);
% view([0,1,0]);
set(gca, 'Fontsize', 15);
%%
bulk_v_mag = sqrt(sum(bulk_v_new.^2,4));
%% Checking the potential change
ix_coll = coll ~= 0;
E_coll = E_pl(ix_coll,:);
dx_coll = x_pl(ix_coll,:) - x_pl(find(ix_coll)-1,:);
dV = dot(E_coll, dx_coll, 2);

%% 
nEng = 50;
Eng_grid = zeros(nEng,1);
nLinE = 10;
dELin_grid = 0.5;
E_gridmax = 200;
Eng_grid(1:nLinE) = dELin_grid * [1:nLinE];
Eng_grid(nLinE + 1:end) = exp(log(nLinE * dELin_grid) + ...
    ( log(E_gridmax)-log(nLinE * dELin_grid) ) * [1:nEng - nLinE]/(nEng - 1 - nLinE));
figure;
plot(Eng_grid)%, 'interp', 'flat')


%%
nLinx = 10;
dxlin = 1;
nX = 100;
nHalfX = nX/2;
x_grid = zeros(nX +1 ,1);
gridx_min = dxlin * nLinx;
gridx_max = 1540;
for i_x = 1:nLinx
    x_grid(nHalfX +1 + i_x) = dxlin * i_x;
    x_grid(nHalfX +1 - i_x) = - dxlin * i_x;
end
    
for i_x = nLinx + 1:nHalfX    
    x_grid(nHalfX +1 + i_x) = exp(log(gridx_min) + ...
        (log(gridx_max) - log(gridx_min))*(i_x - nLinx)/(nHalfX - nLinx));
        
%     x_grid(nHalfX +1 - i_x) = -exp(log(gridx_min) + ...
%         (log(gridx_max) - log(gridx_min))*(i_x - nLinx -1)/(nHalfX - 1 - nLinx));
end
figure;
plot(x_grid)


