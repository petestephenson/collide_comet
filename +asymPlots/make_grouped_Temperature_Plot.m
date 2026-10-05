function [f, ax] = make_grouped_Temperature_Plot(peCasesNo, swCasesNo, plasDensDir, iAngle, varargin)
    %MAKE_GROUPED_TEMPERATURE_PLOT Make a 2x2(+1) grouped temperature plot.
    %
    % Inputs
    %   peCasesNo    cellstr / string array, length >= 5 (PE run IDs)
    %   swCasesNo    cellstr / string array, length >= 5 (SW run IDs)
    %   plasDensDir  char/string path to run files (ending with / or \)
    %   iAngle       index into view list (1..3)
    %
    % Name-value pairs
    %   'PlotMode'   'combined' (default) or 'sw'   -> plot combined Te or SW-only Te
    %   'SWGenIdx'   (default 1)                    -> which SW generation index to use for SW-only
    %
    % Outputs
    %   f   figure handle
    %   ax  cell array of axes handles (1x5)
    
        % ---------------- Options ----------------
        p = inputParser;
        p.addParameter('PlotMode', 'combined', @(s) ischar(s) || isstring(s));
        p.addParameter('SWGenIdx', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
        p.parse(varargin{:});
        plotMode = lower(string(p.Results.PlotMode));
        swGenIdx = p.Results.SWGenIdx;
    
        % ---------------- Figure / layout ----------------
        f = figure('Position', [0, 0, 1000, 600]);
    
        left   = [0.35, 0.65, 0.35, 0.65, 0.07];
        bottom = [0.57, 0.57, 0.10, 0.10, 0.57];
        cbpos  = [0.90, 0.10, 0.01, 0.80];
        wid    = 0.23;
        heig   = 0.35;
    
        plot_lim = 300;
        clims    = [0, 40];
        slices   = {[0], [0], [0]};  % keep as-is for asymPlots.makeTemperaturePlot
    
        vangleList = {[1,0,0], [0,-1,0], [0,0,1]};
        vangle     = vangleList{iAngle};
    
        ax = cell(5, 1);
        cb = cell(5, 1);
        exoAll = cell(3,1);
        for iAsymDim = 1:3
            exoAll{iAsymDim} = asymPlots.compute_grouped_ElDepthMedianExobase(peCasesNo, swCasesNo, plasDensDir, iAsymDim);
        end
        % ---------------- Case loop ----------------
        for iCase = 1:5
            ax{iCase} = axes(f, 'Position', [left(iCase), bottom(iCase), wid, heig]);
            
            switch plotMode
                case "combined"
                    [Te, gridsAsym] = calcCombinedElectronTemperature( ...
                        peCasesNo{iCase}, swCasesNo{iCase}, plasDensDir);
    
                case "sw"
                    % Use helper to compute SW-only temperature
                    [~, gridsAsym, extras] = calcCombinedElectronTemperature( ...
                        peCasesNo{iCase}, swCasesNo{iCase}, plasDensDir, ...
                        'ComputeSWOnly', true, 'SWGenIdx', swGenIdx);
                    Te = extras.Te_sw;
    
                otherwise
                    error("PlotMode must be 'combined' or 'sw'. Got: %s", plotMode);
            end
    
            [X_mesh, Y_mesh, Z_mesh] = meshgrid( ...
                gridsAsym{1}(1:end-1), gridsAsym{2}(1:end-1), gridsAsym{3}(1:end-1));
    
            [f, ax{iCase}, cb{iCase}] = asymPlots.makeTemperaturePlot( ...
                Te, X_mesh, Y_mesh, Z_mesh, slices, plot_lim, vangle, clims, f, ax{iCase});
            for iAsymDim = 1:3
                if iAngle ~= iAsymDim
                    rVals = exoAll{iAsymDim}(iCase).exo.envelope.r;
                    rValsFull = [rVals,rVals(1:end-1:-1)];
                    thetaVals =  exoAll{iAsymDim}(iCase).exo.envelope.theta;
                    thetaValsFull = [thetaVals, -thetaVals(1:end-1:-1)];
                end
            end

            % Clean up labels for top row
            if iCase == 1 || iCase == 2
                xlabel(ax{iCase}, '');
            end
    
            % Remove per-panel colorbars (use one shared)
            delete(cb{iCase});
        end
    
        % ---------------- Shared colorbar ----------------
        cbPlot = colorbar(ax{1}, 'Position', cbpos);
        cbPlot.Label.String = 'T_{e} [eV]';
    end
    
    
    function [Te_tot, grids, extras] = calcCombinedElectronTemperature(peInput, swInput, plasDensDir, varargin)
    %CALCCOMBINEDELECTRONTEMPERATURE Compute electron temperature from combined PE+SW runs.
    %
    % Name-value pairs
    %   'IncGens'         (default [1 2])    generation indices to include in combine_run_velocities
    %   'ComputeSWOnly'   (default false)    additionally compute SW-only temperature
    %   'SWGenIdx'        (default 1)        which SW generation index to use for SW-only temperature
    %
    % Outputs
    %   Te_tot   Combined temperature array (eV)
    %   grids    Grids cell array from asymPlots.loadRunData
    %   extras   struct with optional diagnostics:
    %              .Te_sw   SW-only temperature (if requested)
    %              .me      electron mass used
    %              .bulk_vel_tot, .square_vel_tot
    %              .plas_dens_E_PE, .plas_dens_E_SW
    
        % -------- Parse options --------
        p = inputParser;
        p.addParameter('IncGens', [1 2], @(x) isnumeric(x) && isvector(x) && ~isempty(x));
        p.addParameter('ComputeSWOnly', false, @(x) islogical(x) && isscalar(x));
        p.addParameter('SWGenIdx', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
        p.parse(varargin{:});
    
        inc_gens      = p.Results.IncGens;
        computeSWOnly = p.Results.ComputeSWOnly;
        swGenIdx      = p.Results.SWGenIdx;
    
        % -------- Resolve file paths --------
        peFile = resolveRunFile(peInput, plasDensDir);
        swFile = resolveRunFile(swInput, plasDensDir);
    
        % -------- Load velocities --------
        [pe_bulk_vel, pe_square_vel] = asymPlots.load_velocity_arrays(peFile);
        [sw_bulk_vel, sw_square_vel] = asymPlots.load_velocity_arrays(swFile);
    
        % -------- Load density-vs-energy arrays + grids --------
        [grids, plas_dens_E_PE] = asymPlots.loadRunData(peFile);
        [~,     plas_dens_E_SW] = asymPlots.loadRunData(swFile);
    
        dE_grid = diff(grids{4});
    
        % -------- Combine and compute temperature --------
        [~, bulk_vel_tot, square_vel_tot] = asymPlots.combine_run_velocities( ...
            pe_bulk_vel, pe_square_vel, plas_dens_E_PE, ...
            sw_bulk_vel, sw_square_vel, plas_dens_E_SW, ...
            dE_grid, inc_gens);
    
        me     = asymPlots.get_electron_mass(peFile);
        Te_tot = asymPlots.calculateTemperature(square_vel_tot, bulk_vel_tot, 0, me);
    
        % -------- Diagnostics / optional SW-only --------
        extras = struct();
        extras.me = me;
        extras.bulk_vel_tot = bulk_vel_tot;
        extras.square_vel_tot = square_vel_tot;
        extras.plas_dens_E_PE = plas_dens_E_PE;
        extras.plas_dens_E_SW = plas_dens_E_SW;
    
        if computeSWOnly
            % Match your original convention: select generation index and permute dims
            square_vel_sw = permute(sw_square_vel(:,:,:,:,swGenIdx), [2,3,4,1]);
            bulk_vel_sw   = permute(sw_bulk_vel(:,:,:,:,swGenIdx),   [2,3,4,1]);
            extras.Te_sw  = asymPlots.calculateTemperature(square_vel_sw, bulk_vel_sw, 0, me);
        end
    end
    
    
    function runFile = resolveRunFile(inputVal, plasDensDir)
    % Resolve either:
    %   - a full path (char/string) that ends in .nc
    %   - or a case id that needs building as <plasDensDir>/plas_dens<id>.nc
    
        if nargin < 2
            plasDensDir = "";
        end
    
        s = string(inputVal);
    
        % If it looks like a file path already, use it
        if endsWith(lower(s), ".nc") || contains(s, filesep) || contains(s, "/") || contains(s, "\")
            runFile = char(s);
            return;
        end
    
        % Otherwise build from directory + case id
        runFile = char(fullfile(plasDensDir, "plas_dens" + s + ".nc"));
    end
    