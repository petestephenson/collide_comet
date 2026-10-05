
function [XS_Eng, XS, E_th] = load_PE_XSections(srcDir)
%LOAD_PE_XSECTIONS Loads photoionization cross-section data.
%   [XS_ENG, XS, E_TH] = LOAD_PE_XSECTIONS(SRCDIR) loads photoionization
%   cross-section data from the specified directory and calculates related
%   parameters.
%
%   Inputs:
%   - srcDir: Directory containing the photoionization cross-section data.
%
%   Outputs:
%   - XS_Eng: Energy grid.[eV]
%   - XS: Photoionization cross-sections. [cm^2]
%   - E_th: Threshold energy for each cross-section.[eV]

h = 6.63e-34;
c = 3e8;

XS_path = [srcDir, '\Energy_Degradation\XS\H2O_Photoionization\'];
XS_1b1 = csvread([XS_path, '1b1.csv']);
XS_1b2 = csvread([XS_path, '1b2.csv']);
XS_2a1 = csvread([XS_path, '2a1.csv']);
XS_3a1 = csvread([XS_path, '3a1.csv']);

nEng = 10001;
XS_Eng = linspace(0,80, nEng);
XS = zeros(nEng, 4);
LineName = {'1b1', '1b2', '2a1', '3a1'};
XS(:,1) = interp1(XS_1b1(:,1), XS_1b1(:,2), XS_Eng,'linear', 'extrap');
XS(XS_Eng < XS_1b1(1,1),1) = 0;
XS(:,2) = interp1(XS_1b2(:,1), XS_1b2(:,2), XS_Eng,'linear', 'extrap');
XS(XS_Eng < XS_1b2(1,1),2) = 0;
XS(:,3) = interp1(XS_2a1(:,1), XS_2a1(:,2), XS_Eng,'linear', 'extrap');
XS(XS_Eng < XS_2a1(1,1),3) = 0;
XS(:,4) = interp1(XS_3a1(:,1), XS_3a1(:,2), XS_Eng,'linear', 'extrap');
XS(XS_Eng < XS_3a1(1,1),4) = 0;

XS = XS * 1e-18; %cm^2 - conversion for Mb
lam_th = [984, 720.7, 411.4, 907];%\AA
E_th = h*c*1e10./(lam_th * 1.6e-19);

XS(XS_Eng<E_th(1),1) = 0;

end