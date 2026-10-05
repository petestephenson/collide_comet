
function [f,ax] = make_SW_tempRatioplot(peCasesNo, swCasesNo, plasDensDir, iAngle)

    f = figure();
    left = [0.1,0.55, 0.1, 0.55];
    bottom = [0.6, 0.6, 0.15, 0.15];
    wid = 0.4;
    heig = 0.3;
    ax = cell(5,1);
    cb = cell(5,1);
    vangleList = {[1,0,0],[0,-1,0],[0,0,1]};
    viewNamelist = { 'YZ','XZ','XY'};
    slices = {[0],[0],[0]};
    plot_lim=300;
%     iAngle = 2;
    vangle = vangleList{iAngle};
    swRefNo = swCasesNo{5};
    peRefNo = peCasesNo{5};
    swRefName = [plasDensDir, 'plas_dens', swRefNo, '.nc'];
    peRefName = [plasDensDir, 'plas_dens', peRefNo, '.nc'];

    [qSpherical, qDayside, qUniform, qVec] = asymPlots.loadOutgassingProperties(swRefName);
    [sw_bulk_vel, sw_square_vel] = asymPlots.load_velocity_arrays(swRefName);
    bulk_vel_sw_only = permute(sw_bulk_vel(:,:,:,:,1),[2,3,4,1]);
    square_vel_sw_only = permute(sw_square_vel(:,:,:,:,1),[2,3,4,1]);
    me = asymPlots.get_electron_mass(peRefName);
    temperatureSW_Ref = asymPlots.calculateTemperature(square_vel_sw_only, bulk_vel_sw_only, 0,me);
    for iCase = [1,2,3,4]
        ax{iCase} = axes(f, 'Position',[left(iCase), bottom(iCase), wid, heig]);
        


        peAsymNo = peCasesNo{iCase};
        swAsymNo = swCasesNo{iCase};
        peAsymName = [plasDensDir, 'plas_dens', peAsymNo, '.nc'];
        swAsymName = [plasDensDir, 'plas_dens', swAsymNo, '.nc'];
        [gridsAsym, ~] = asymPlots.loadRunData(peAsymName);
        [X_mesh, Y_mesh, Z_mesh] = meshgrid(gridsAsym{1}(1:end-1), gridsAsym{2}(1:end-1), gridsAsym{3}(1:end-1));

        me = asymPlots.get_electron_mass(peAsymName);
        if iCase == 5
            qDayside = 0.0;
            qUniform = qSpherical;
        end
        % ---- Calculate Temperatures ------

        [sw_bulk_vel, sw_square_vel] = asymPlots.load_velocity_arrays(swAsymName);
        bulk_vel_sw_only = permute(sw_bulk_vel(:,:,:,:,1),[2,3,4,1]);
        square_vel_sw_only = permute(sw_square_vel(:,:,:,:,1),[2,3,4,1]);

        me = asymPlots.get_electron_mass(peAsymName);
    %     temperatureAsym = calculateTemperature(square_vel_tot, bulk_vel_tot, 0,me);
        temperatureSW_Asym = asymPlots.calculateTemperature(square_vel_sw_only, bulk_vel_sw_only, 0,me);
        plot_array = temperatureSW_Asym ./ temperatureSW_Ref;
        [f,ax{iCase}, cb{iCase}] = asymPlots.makeTemperatureRatioPlot(plot_array,X_mesh,Y_mesh,Z_mesh, slices, plot_lim,  vangle, f, ax{iCase});
        if mod(iCase,2) == 1 
            delete(cb{iCase})
        end
        if iCase <3
            xlabel(ax{iCase}, '')
        end
    end
end