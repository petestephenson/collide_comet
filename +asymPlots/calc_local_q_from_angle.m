function [Qpts] = calc_local_q_from_angle(cos2Max, qUniform, qDayside, qVec)
    % u_gas = 1.0; %kms^-1

    % cos2Max 
    % cos2max_unedited = cos2Max;
    cos2Max_flat(cos2Max < 0) = 0.0;
%     ftemp = figure();
%     axtemp = axes(ftemp);
    
    Qpts = qUniform + 4.* qDayside.* cos2Max_flat;
    % cos2MaxFlat = reshape(cos2Max_unedited, [], 1);
    