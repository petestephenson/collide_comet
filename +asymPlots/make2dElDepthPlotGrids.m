function [angularAxis, radialGrid] = make2dElDepthPlotGrids(nElDepthThetaBounds, elDepthRadGrid)
%     rMax = 300;
    radialGrid = elDepthRadGrid; % Linear grid from 0 to 300 km
    angularAxis = acos(linspace(-1, 1, nElDepthThetaBounds)); % Linear grid from -1 to 1
    radialGrid = 0.5 * (radialGrid(2:end)+ radialGrid(1:end-1));
    angularAxis = 0.5*(angularAxis(2:end)+angularAxis(1:end-1));
    
end
