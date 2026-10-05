
function Temp = calculateTemperature(square_vel, bulk_vel,i_dim, me)
    qe = 1.6e-19;%C
    if i_dim == 0 % Magnitude of Velocity
        sq_vel_temp = sum(square_vel,4, 'omitnan'); % Finds magnitude of bulk vel
        bulk_vel_temp = sqrt(sum(bulk_vel.^2,4, 'omitnan'));

    else % Velocity Component
        sq_vel_temp = square_vel(:,:,:,i_dim);
        bulk_vel_temp = bulk_vel(:,:,:,i_dim);

    end
   
    sq_vel_temp = sq_vel_temp - (bulk_vel_temp.^2);
  
    Temp = me * (sq_vel_temp) / qe;
    
    if i_dim == 0 
        Temp= Temp / 3; %Taking 1/3 * Tr(Temp)
    end

end
