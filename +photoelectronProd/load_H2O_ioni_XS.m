
function [XS_Eng, XS, E_th] = load_H2O_ioni_XS(srcDir)
% LOAD_H2O_IONI_XS Loads H2O photoionization cross-section data.
%   [XS_ENG, XS, E_TH] = LOAD_H2O_IONI_XS(SRC_DIR) loads H2O
%   photoionization cross-section data from the specified source directory.
%
%   Inputs:
%   - srcDir: Source directory path.
%
%   Outputs:
%   - XS_Eng: Energy grid.
%   - XS: Photoionization cross-sections.
%   - E_th: Threshold energy for each cross-section.

h = 6.63e-34;
c = 3e8;
H2O_XS = readmatrix([srcDir, 'Energy_Degradation/XS/H2O_Photoionization/H2O_ionization_XS.txt'], 'NumHeaderLines', 56);
XS_lam = H2O_XS(:,1);
XS = H2O_XS(:, 6:9); %[OH+, O+, H+, H2O+] - cm^2
LineName = {'OH+', 'O+', 'H+', 'H2O+'};
E_th = [18.11, 18.65, 18.72, 12.616];
XS_Eng = h * c * 1e10./(flip(XS_lam) * 1.6e-19); %J
XS = flip(XS);
end