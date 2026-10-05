function [Qpts, cos2Max] = calc_local_q_from_pos(x,y,z, qUniform, qDayside, qVec)
    % u_gas = 1.0; %kms^-1
%     Qpt = qUniform;
    rPts = sqrt(x.^2+y.^2+z.^2);

    xNorm = x./rPts;
    yNorm = y./rPts;
    zNorm = z./rPts;

    cos2Max = xNorm.*qVec(1) +  yNorm.*qVec(2) +  zNorm.*qVec(3);
    % cos2max_unedited = cos2Max;
    cos2Max_flat(cos2Max < 0) = 0.0;
%     ftemp = figure();
%     axtemp = axes(ftemp);
    
    Qpts = qUniform + 4.* qDayside.* cos2Max_flat;
    % cos2MaxFlat = reshape(cos2Max_unedited, [], 1);
    