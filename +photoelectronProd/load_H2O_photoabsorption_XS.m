
function [absLam, absEng, absXS] = load_H2O_photoabsorption_XS(srcDir)
% LOAD_H2O_PHOTOABSORPTION_XS Loads H2O photoabsorption cross-section data.
%   [ABS_LAM, ABS_ENG, ABS_XS] = LOAD_H2O_PHOTOABSORPTION_XS(SRC_DIR) loads H2O
%   photoabsorption cross-section data from the specified source directory.
%
%   Inputs:
%   - srcDir: Source directory path.
%
%   Outputs:
%   - absLam: Wavelength grid.
%   - absEng: Energy grid.
%   - absXS: Photoabsorption cross-sections.
h = 6.63e-34;
c = 3e8;
XS_path = [srcDir, 'Energy_Degradation/XS/H2O_Photoionization/'];
H2O_abs_XS = readmatrix([XS_path, 'H2O_absorption_XS.txt'], 'NumHeaderLines', 2);
absLam = H2O_abs_XS(:,1);
absEng = h * c * 1e10 ./ (flip(absLam) * 1.6e-19);
absXS = flip(H2O_abs_XS(:,3)) * 1e-18;
end
