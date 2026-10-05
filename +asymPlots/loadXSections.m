function [XS, XS_Eng, E_th] = loadXSections()
% app.oneDrive = [basePath, 'OneDrive - Imperial College London\'];
%             app.EngDegDir = [app.oneDrive, 'Energy_Degradation\'];
    XSFolder ='C:/Users/pete4/OneDrive - Imperial College London/Energy_Degradation/XS/';
    XSFile = 'Coll_XS.nc';
    ncid = netcdf.open([XSFolder, XSFile]);
    
    xs_id = netcdf.inqVarID(ncid, 'XS');
    xsEng_id = netcdf.inqVarID(ncid, 'XS_Eng');
    Eth_id = netcdf.inqVarID(ncid, 'E_th');
    XS = netcdf.getVar(ncid, xs_id);
    XS_Eng = netcdf.getVar(ncid, xsEng_id);
    E_th = netcdf.getVar(ncid, Eth_id);
    
    netcdf.close(ncid);
end