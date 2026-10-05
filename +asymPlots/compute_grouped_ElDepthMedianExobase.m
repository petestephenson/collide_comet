function outAll = compute_grouped_ElDepthMedianExobase( ...
    peCasesNo, swCasesNo, plasDensDir, iAngle, opts)
%COMPUTE_GROUPED_ELDEPTHMEDIANEXOBASE Compute median depths + exobase hull for multiple PE/SW cases.
%
% Mirrors the data-loading logic in make_grouped_ElDepth2d_Plot, but does NO plotting.
%
% Inputs
%   peCasesNo, swCasesNo : cell arrays / string arrays of case IDs (length >= 5)
%   plasDensDir          : directory containing plas_dens<case>.nc (include trailing slash or not)
%   iAngle               : dim index passed through to median/exobase code (your "iAngle")
%   opts                 : (optional) options struct passed to compute_2dElDepthMedianAndExobase
%
% Output
%   outAll : 1xN struct array (N = numel(casesUsed)) with fields:
%       .iCase
%       .peCaseNo, .swCaseNo
%       .peFile, .swFile
%       .elDepthRadEdges, .elDepthGrid, .elDepthDist_all
%       .medianDepth, .log10MedianDepth, .inExo, .rMids, .thetaMids, .x, .y, .exo
%
% Notes
% - Uses asymPlots.load_electron_depth_dist(...) exactly as in your plotting function.
% - Assumes both PE and SW runs share compatible grids; uses the PE grid for edges/grid.

    if nargin < 5 || isempty(opts), opts = struct(); end

    % Normalise directory
    plasDensDir = char(plasDensDir);
    if ~isempty(plasDensDir) && ~endsWith(plasDensDir, filesep)
        plasDensDir = [plasDensDir, filesep];
    end

    casesUsed = 1:5;  % matches your loop; change if needed
    nCases = numel(casesUsed);

    outAll = repmat(struct(), 1, nCases);

    for ii = 1:nCases
        iCase = casesUsed(ii);

        peCaseNo = peCasesNo{iCase};
        swCaseNo = swCasesNo{iCase};

        peFile = [plasDensDir, 'plas_dens', char(peCaseNo), '.nc'];
        swFile = [plasDensDir, 'plas_dens', char(swCaseNo), '.nc'];

        % Load electron-depth distributions
        % (Your original code had pe/sw variable names swapped; here we keep the INTENT:
        %  load from peFile and swFile respectively.)
        [eRadEdges_pe, eGrid_pe, eDist_pe] = asymPlots.load_electron_depth_dist(peFile);
        [eRadEdges_sw, eGrid_sw, eDist_sw] = asymPlots.load_electron_depth_dist(swFile);

        % Basic compatibility checks (lightweight)
        if numel(eRadEdges_pe) ~= numel(eRadEdges_sw) || any(eRadEdges_pe(:) ~= eRadEdges_sw(:))
            warning('Case %d: radial edges differ between PE and SW; proceeding with PE edges.', iCase);
        end
        if numel(eGrid_pe) ~= numel(eGrid_sw) || any(eGrid_pe(:) ~= eGrid_sw(:))
            warning('Case %d: depth grid differs between PE and SW; proceeding with PE grid.', iCase);
        end

        % Combine PE + SW contributions
        eDist_all = eDist_pe + eDist_sw;

        % Compute median depths + exobase hull (no plotting)
        out = asymPlots.compute_2dElDepthMedianAndExobase(eRadEdges_pe, eGrid_pe, eDist_all, iAngle, opts);

        % Package per-case output
        outAll(ii).iCase           = iCase;
        outAll(ii).peCaseNo        = peCaseNo;
        outAll(ii).swCaseNo        = swCaseNo;
        outAll(ii).peFile          = peFile;
        outAll(ii).swFile          = swFile;

        outAll(ii).elDepthRadEdges = eRadEdges_pe;
        outAll(ii).elDepthGrid     = eGrid_pe;
        outAll(ii).elDepthDist_all = eDist_all;  %#ok<*NASGU> (keep if you want downstream reuse)

        % Computed products
        outAll(ii).medianDepth      = out.medianDepth;
        outAll(ii).log10MedianDepth = out.log10MedianDepth;
        outAll(ii).inExo            = out.inExo;
        outAll(ii).rMids            = out.rMids;
        outAll(ii).thetaMids        = out.thetaMids;
        outAll(ii).x                = out.x;
        outAll(ii).y                = out.y;
        outAll(ii).exo              = out.exo;
    end
end
