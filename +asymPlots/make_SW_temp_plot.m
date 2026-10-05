function [f,ax] = make_SW_temp_plot(peCasesNo, swCasesNo, plasDensDir, iAngle)

    f = figure(Position=[0,0,1000, 600]);
    left = [0.35, 0.65, 0.35, 0.65, 0.07];
    bottom = [0.57, 0.57, 0.1, 0.1,0.57];
    cbpos = [0.9,0.1,0.01,0.8];
    wid = 0.23;
    heig = 0.35;
    ax = cell(5,1);
    cb = cell(5,1);
    
    clims = [0,60];

    vangleList = {[1,0,0],[0,-1,0],[0,0,1]};
    viewNamelist = { 'YZ','XZ','XY'};
    slices = {[0],[0],[0]};
    plot_lim=300;
%     iAngle = 2;
    vangle = vangleList{iAngle};
    swRefNo = swCasesNo{5};
    swRefName = [plasDensDir, 'plas_dens', swRefNo, '.nc'];

    [qSpherical, qDayside, qUniform, qVec] = asymPlots.loadOutgassingProperties(swRefName);
    for iCase = [1,2,3,4,5]
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

        [f,ax{iCase}, cb{iCase}] = asymPlots.makeTemperaturePlot(temperatureSW_Asym,X_mesh,Y_mesh,Z_mesh, slices, plot_lim,  vangle,clims, f, ax{iCase});
        if or((iCase == 1), (iCase ==2))
            xlabel(ax{iCase}, '');
        end
%         if mod(iCase,2) == 1
        delete(cb{iCase});
%         end

    end
    cbPlot = colorbar(ax{1}, 'Position', cbpos);
    plot_str = 'T_{e} [eV])';
    cbPlot.Label.String = plot_str;
end