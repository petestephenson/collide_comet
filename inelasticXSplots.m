%% ------ Simplified Collision cross sections for PIC simulation ----------
OS_name = system_dependent('getos');
i_windows = strcmp('Microsoft',OS_name(1:9));
if i_windows
    one_drive = 'C:/Users/pete4/OneDrive - Imperial College London/';
else
    one_drive = '/Users/pps18/OneDrive - Imperial College London/';
end
%%
% load('xsections.mat');
load('xsections_w_scatter.mat');
% XS - 103 x 12
% XS_eng - 103 x 1
% E_th - 1 x 12
nXS = numel(E_th);
nEng = numel(XS_Eng);
%% -------- Extend the Dissociation XS with straight line in loglog -------
ix_diss = 6;
ix_elec = 5;
last_point = find(XS(:, ix_diss)>0, 1, 'last');
ix_100 = find(XS_Eng == 100, 1 ,'first');
grad_XS = log(XS(end,ix_elec)./XS(ix_100,ix_elec)) ./ log(XS_Eng(end)./XS_Eng(ix_100));
grad_log_XS =  log(XS(2:end,ix_diss)./XS(1:end-1,ix_diss)) ./ log(XS_Eng(2:end)./XS_Eng(1:end-1));
XS(last_point + 1 :end , ix_diss) = exp(log(XS(last_point, ix_diss)) + grad_XS ...
    * (log(XS_Eng(last_point + 1 :end)/XS_Eng(last_point)))); 
%% -------- Loss Rates ----------------------------------------------------
% Excitation: Loss Rate = XS * E_th
% Ionization: Loss Rate = XS * (E_th + 0.5*E)
% Double Ionization: LR = XS * (E_th + [E/3])
ix_exc = 1:6; %Excitaiton
ix_ioni = 7:11; %Single Ionization indices
ix_ioni2 = 12; % Double Ionization indices

LR = zeros(nEng, nXS);
% ------------- Calculate Loss Rate --------------------------------------
LR(:, ix_exc) = XS(:, ix_exc) .* E_th(ix_exc);
LR(:, ix_ioni) = 0.5 * XS(:, ix_ioni) .* (E_th(ix_ioni) + XS_Eng);
LR(:, ix_ioni2) = XS(:, ix_ioni2) .* (E_th(ix_ioni2) + 2 * XS_Eng) / 3;

d_Eng_ioni = zeros(nEng,6);
d_Eng_ioni(:,ix_ioni-6) = 0.5 .* (E_th(ix_ioni) + XS_Eng);
d_Eng_ioni(:,ix_ioni2-6) =  (E_th(ix_ioni2) + 2 * XS_Eng)/3;




%% ------------------ Plot Cross Sections ---------------------------------
process_names = {'Elastic', 'Rotation', 'Vibration 010', 'Vibration Stretch',...
    'Electronic', 'Dissociation', 'H_{2}O^{+}', 'OH^{+}','O^{+}','H_{2}^{+}','H^{+}', 'O^{2+}'};
f = figure('units', 'centimeters','Position', [1,1,26,12]); 

% ---------- plot options ----------------------------
lef = [0.09, 0.58 0.7];
bot = [0.15, 0.15, 0.1];
wid = 0.4 * ones(1,2);%, 0.3, 0.3];
heig = [0.8, 0.8, 0.8];
cb_wid = 0.01;
cb_heig = 0.8;
yLims = {[1e-18,1e-12], [1e-20, 3e-16]};
xLims = {[1e-2,3e3], [10, 3e2]};

plot_range = {1:6, 7:12};
nPlots = 2;
legCols = [1:2];
ax = cell(nPlots, 1);

for iPlot = 1:nPlots
    ax{iPlot} = subplot('Position', [lef(iPlot), bot(iPlot), wid(iPlot), heig(iPlot)]);

%     figure('units', 'centimeters', 'position', [1,1,25,15]);
    hold on

    for i_xs = plot_range{iPlot}
        plot(XS_Eng, XS(:,i_xs), 'DisplayName',process_names{i_xs}, 'linewidth', 3);
    end
    
    xlabel('Energy [eV]');
    ylabel('Cross Section [cm^2]') 
    set(gca, 'xscale','log');
    set(gca, 'yscale','log');
%     ylim([1e-20,3e-16]);
    xlim(xLims{iPlot});
    set(gca, 'fontsize', 14);
    grid on;

    legend('NumColumns', legCols(iPlot));
    ylim(yLims{iPlot});
end
ax{2}.Legend.Position(1:2) = [0.7965    0.3979];

%% -------------------- Plot Loss Rates -----------------------------------
i_plot_exc = 1; %Plots excitation Loss Rates
i_plot_ioni = 1; % Plots Ionization Loss Rates

process_names = {'Elastic', 'Rotation', 'Vib 010', 'Vib Stretch',...
    'Electronic', 'Dissociation', 'Ionization'};
plot_range = 1:7;

LRplot(:,7) = sum(LR(:,7:12),2);
if or(i_plot_exc, i_plot_ioni)
    figure('units', 'centimeters', 'position', [1,1,22,12]);
    hold on
    for i_xs = plot_range
        plot(XS_Eng, LRplot(:,i_xs), 'DisplayName',process_names{i_xs}, 'linewidth', 3);
    end
    xlabel('Energy [eV]');
    ylabel('Loss Rate [eV \cdot cm^2]') 
    set(gca, 'xscale','log');
    set(gca, 'yscale','log');
    ylim([5e-19,3e-14]);
%     ylim([1e-18,3e-15]);
    xlim([1e-2,3e2]);
    set(gca, 'fontsize', 15);
    grid on;

    legend('Location', 'southwest');
end
