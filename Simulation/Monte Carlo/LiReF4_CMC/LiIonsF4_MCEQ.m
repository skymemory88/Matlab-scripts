function [energies,lattice,E_0,lat_mom]=LiIonsF4_MCEQ(params,ion,inter,lattice,T,E_0,lat_mom,field_change)

N = params.N_Er; % Nb of spins
iterMax = params.NiterEQ * N;
estm = 0; % initiate an estimator
energies=zeros(1,iterMax);% array for the energy per spin

for j=1:N
    lattice{j}.mom=lat_mom(1:3,j)';
    lattice{j}.beta=lat_mom(4,j);
    lattice{j}.alpha=lat_mom(5,j);
    lattice{j}.energy=lat_mom(6,j);
end

for iterations=1:iterMax
    % Rotate one moment
    chosen=0;
    while(chosen==0)
        new_ion=randi(size(lattice,2),1,1);
        [alpha_new, beta_new]=random_angles;
        [new_mom, new_E_Zee,field_change]=cl_eff_mod(alpha_new, beta_new, params.field, lattice{new_ion},field_change);
        chosen=1;
    end
    
    % Computes the energy variation
    dE=energy_update(ion,params.field,lattice,params.L,inter,new_ion,new_mom,new_E_Zee,lattice{new_ion}.energy);

    % Metropolis-Hastings
    acc=metropolis(dE,T);
%     % Glauber rate
%     acc=glauber(dE,T);
    if acc==1
        lattice{new_ion}.mom=new_mom;
        lattice{new_ion}.alpha=alpha_new;
        lattice{new_ion}.beta=beta_new;
        lattice{new_ion}.energy=new_E_Zee;
        lat_mom(1:3,new_ion)=new_mom;
        lat_mom(4,new_ion)=alpha_new;
        lat_mom(5,new_ion)=beta_new;
        lat_mom(6,new_ion)=new_E_Zee;
    end
    E_0=E_0+acc*dE;
    energies(iterations) = E_0/N;
    % Convergence test
    if mod(iterations,N) == 0
        old_estm = estm;
        estm = var(energies(iterations-N+1:iterations))/N;
        if (abs(old_estm-estm) <= estm*0.0015) % When the estimator stops changing, then interrupt the calculation
            energies(iterations:end) = energies(iterations);
            disp('Convergence reached, iteration terminated.') % Stop the iteration if the energy converges
            break
        elseif iterations == iterMax
            disp('Maximum limit reached, iteration terminated.');
        end
    end
end
energies = energies(1:N:end);
clearvars -except energies lattice E_0 lat_mom
return