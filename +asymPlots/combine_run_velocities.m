

function [tot_plas_dens, ave_bulk_vel, ave_square_vel] = combine_run_velocities(pe_bulk_vel, pe_square_vel,plas_dens_E_PE_ref, ...
    sw_bulk_vel, sw_square_vel, plas_dens_E_SW_ref, dE_grid, inc_gens)

    % --- Calculate plas_dens for each generation -----
    plas_dens_PE_gen = sum(plas_dens_E_PE_ref(:,:,:,:,:)...
                    .* reshape(dE_grid(:),1,1,1,[]),4);
    plas_dens_SW_gen = sum(plas_dens_E_SW_ref(:,:,:,:,:)...
                    .* reshape(dE_grid(:),1,1,1,[]),4);

    % --- Calculate total sw_square and bulk velocities ----
    sq_vel_temp = permute(pe_square_vel,[2,3,4,1,5]);
    sq_vel_PE = sum(sq_vel_temp(:,:,:,:,inc_gens) .* ...
            plas_dens_PE_gen(:,:,:,:,inc_gens),5, 'omitnan');
    
    bulk_vel_temp = permute(pe_bulk_vel,[2,3,4,1,5]);
    bulk_vel_PE = sum(bulk_vel_temp(:,:,:,:,inc_gens) .*...
            plas_dens_PE_gen(:,:,:,:,inc_gens),5, 'omitnan');

    sq_vel_temp = permute(sw_square_vel,[2,3,4,1,5]);
    sq_vel_SW = sum(sq_vel_temp(:,:,:,:,inc_gens) .* ...
            plas_dens_SW_gen(:,:,:,:,inc_gens),5, 'omitnan');

    bulk_vel_temp = permute(sw_bulk_vel,[2,3,4,1,5]);
    bulk_vel_SW = sum(bulk_vel_temp(:,:,:,:,inc_gens) .*...
            plas_dens_SW_gen(:,:,:,:,inc_gens),5, 'omitnan');

    plas_dens_PE  = sum(plas_dens_PE_gen(:,:,:,:,inc_gens), 5, 'omitnan');
    plas_dens_SW  = sum(plas_dens_SW_gen(:,:,:,:,inc_gens), 5, 'omitnan');
    %---- Calculate total -------

    tot_square_vel = sq_vel_PE   + sq_vel_SW;
    tot_bulk_vel   = bulk_vel_PE + bulk_vel_SW;

    tot_plas_dens  = plas_dens_PE + plas_dens_SW;

    ave_square_vel = tot_square_vel ./ tot_plas_dens;
    ave_bulk_vel   = tot_bulk_vel   ./ tot_plas_dens;
    
end