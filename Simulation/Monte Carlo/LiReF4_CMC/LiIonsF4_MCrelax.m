function [energies,mag,maltx,malty,lattice,E_0,lat_mom,params]=LiIonsF4_MCrelax(ion,params,inter,lattice,field,T,E_0,lat_mom)

L=params.L;
N=4*(L^3); % number of sites
Niter=params.Nitermeas;

energies=zeros(1,Niter);
mag=zeros(3,4,Niter);
maltx=zeros(Niter,1);
malty=zeros(Niter,1);


for j=1:N
    lattice{j}.mom=lat_mom(1:3,j)';
    lattice{j}.beta=lat_mom(4,j);
    lattice{j}.alpha=lat_mom(5,j);
    lattice{j}.energy=lat_mom(6,j);
end

for iterations=1:Niter
    
    % Rotate one moment
    chosen=0;
    while(chosen==0)
        new_ion=randi(size(lattice,2),1,1);
        [alpha_new, beta_new]=random_angles;
        [new_mom, new_E_Zee,params]=cl_eff_mod(alpha_new, beta_new, field, lattice{new_ion},params);
        chosen=1;
    end
    
    % Computes the energy variation
    dE=energy_update(ion,field,lattice,L,inter,new_ion,new_mom,new_E_Zee,lattice{new_ion}.energy);
    
    % Metropolis-Hastings
    acc=metropolis(dE,T);
    % Glauber rate
    %acc=glauber(dE,T);
    
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
    E=E_0+acc*dE;
    
    energies(iterations)=E/N;
    
    % magnetizations and magnetizations squared per type of site
    mag(:,1,iterations)=sum(lat_mom(1:3,1:4:end)/N,2);
    mag(:,2,iterations)=sum(lat_mom(1:3,2:4:end)/N,2);
    mag(:,3,iterations)=sum(lat_mom(1:3,3:4:end)/N,2);
    mag(:,4,iterations)=sum(lat_mom(1:3,4:4:end)/N,2);
    
    for m=1:4
        maltx(iterations)=maltx(iterations)+params.C(m,1)*mag(1,m,iterations);
        malty(iterations)=malty(iterations)+params.C(m,2)*mag(2,m,iterations);
    end
    
    E_0=E;
    
end

return
