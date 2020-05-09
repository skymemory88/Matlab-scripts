function [energies,Egs,Egs2,mag,magsq,sq0,sqx,sqy,lattice,E_0,lat_mom,params]=LiIonsF4_MC(ion,params,inter,lattice,field,T,E_0,lat_mom)
persistent phasex phasey phase0 % local variables

L=params.L;
N=4*(L^3); % number of sites
Niter=params.Nitermeas;
N_meas=floor(Niter/N);
%N_meas=1;

if isempty(phasex)
    phase0=ones(N);
    qx=2*pi/params.L;
    qy=qx;
    pos=zeros(3,N);
    for i=1:N
        pos(:,i)=lattice{i}.position;
    end
    phasex=exp(1i*qx*pos(1,:))'*exp(1i*qx*pos(1,:));
    phasey=exp(1i*qy*pos(2,:))'*exp(1i*qy*pos(2,:));
   
end

energies=zeros(1,Niter);

% arrays to write the magnetizations
mag=zeros([3,4,N_meas]);
magsq=zeros([3,3,4,N_meas]);
% form factors
sq0=zeros([3,3,4,4,N_meas]);
sqx=zeros([3,3,4,4,N_meas]);
sqy=zeros([3,3,4,4,N_meas]);
% energies per spin and energies per spin squared
Egs=zeros([1,N_meas]);
Egs2=zeros([1,N_meas]);
% counter
p=1;

for j=1:N
    % disp(j)
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
    %fprintf('dE = %f meV, acc = %d \n',dE,acc)
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
         
    if rem(iterations,N)==0
        Egs(p)=energies(iterations);
        Egs2(p)=(energies(iterations))^2;
        
        for i=1:3
            for j=1:3
                for m=1:4
                    for mp=1:4
                        sq0(i,j,m,mp,p)=lat_mom(i,m:4:end)*phase0(m:4:end,mp:4:end)*lat_mom(j,mp:4:end)';
                        sqx(i,j,m,mp,p)=lat_mom(i,m:4:end)*phasex(m:4:end,mp:4:end)*lat_mom(j,mp:4:end)';
                        sqy(i,j,m,mp,p)=lat_mom(i,m:4:end)*phasey(m:4:end,mp:4:end)*lat_mom(j,mp:4:end)';
                    end
                end
            end
        end
        
        % magnetizations and magnetizations squared per type of site
        mag(:,1,p)=sum(lat_mom(1:3,1:4:end)/N,2);
        mag(:,2,p)=sum(lat_mom(1:3,2:4:end)/N,2);
        mag(:,3,p)=sum(lat_mom(1:3,3:4:end)/N,2);
        mag(:,4,p)=sum(lat_mom(1:3,4:4:end)/N,2);
        
        for a=1:3
            for b=1:3
                for m=1:4                   
                    magsq(a,b,m,p)=sum((lat_mom(a,m:4:end)/N).*(lat_mom(b,1:4:end)/N),2);
                end
            end
        end       
        p=p+1;
    end

    E_0=E;
    
end

return
