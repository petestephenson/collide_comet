function nDens = calcNeutralDens(x,y,z, qUniform, qDayside, qVec)
    u_gas = 1.0; %kms^-1
%     Qpt = qUniform;
    rPts = sqrt(x.^2+y.^2+z.^2);

    xNorm = x./rPts;
    yNorm = y./rPts;
    zNorm = z./rPts;

    cos2Max = xNorm.*qVec(1) +  yNorm.*qVec(2) +  zNorm.*qVec(3);
    cos2Max_unedited = cos2Max;
    cos2Max(cos2Max < 0) = 0.0;
%     ftemp = figure();
%     axtemp = axes(ftemp);
    
    Qpts = qUniform + 4.* qDayside.* cos2Max;
    cos2MaxFlat = reshape(cos2Max_unedited, [], 1);
    QptsFlat = reshape(Qpts, [], 1);
    angleFlat = acosd(cos2MaxFlat);
%     scatter(axtemp, angleFlat, QptsFlat)
    % # set xlabel as cos2Max and ylabel as Qpts
%     xlabel('angle2max');
%     ylabel('Qpts');
%     xMid = reshape(0.5 * (grids{1}(2:end) + grids{1}(1:end-1)), [], 1);
%     yMid = reshape(0.5 * (grids{2}(2:end) + grids{2}(1:end-1)), 1, []);
%     zMid = reshape(0.5 * (grids{3}(2:end) + grids{3}(1:end-1)), 1, 1, []);
%     rMid = sqrt(x.^2 + y.^2 + y.^2);

    
    nDens = Qpts ./ ( 4* pi* u_gas * (rPts.^2)*1e15); % s^-1/(km^3s^-1)
end