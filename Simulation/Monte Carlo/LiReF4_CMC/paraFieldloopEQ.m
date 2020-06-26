function [relaxE,EQlat_mom,lattice,params]=paraFieldloopEQ(ion,params,inter,lattice,lat_mom)

N=4*(params.L)^3;
wforce = gcp(); % Start a new parallel pool if it doesn't exist

% Initiate containers for equilibrium configurations as well as equilibrium energies
relaxE=zeros([size(params.field,1),params.NiterEQ],wforce.NumWorkers);
EQlat_mom=cell(wforce.NumWorkers,size(params.field,1));

% Calculate the energy of the current configuration as the initial state for all field points
E_0 = energy0(ion,params,inter,lattice,lat_mom);
disp(['E_0=',num2str(E_0/N)]);

% Seed the first (T = min(T)) thermalization process with the initial configuration
relaxE(1,:,:) = E_0;
EQlat_mom{:,1} = lat_mom; 

% Carry out thermalization in parallel across all cores and save the data in the end
spmd
    for jj = 1:size(params.field,1)
        field_change = true;
        if labindex == 1
            disp(['T = ',num2str(params.temp),' Field = [',num2str(params.field(jj,1)),',',num2str(params.field(jj,2)),',',num2str(params.field(jj,3)),']']);
        end
        [relaxE(jj,:,labindex),~,~,EQlat_mom{labindex,jj}]=LiIonsF4_MCEQ(ion,params.L,params.NiterEQ,inter,lattice,params.field(jj,:),params.temp,relaxE(1,end,labindex),EQlat_mom{labindex,jj},field_change);
        if ii < size(params.field,1)
            relaxE(jj+1,:,labindex) = relaxE(jj,:,labindex); % Seed the thermalization process of the next temperature point
            EQlat_mom{labindex,jj+1} = EQlat_mom{labindex,jj}; % Seed the tehermalization process of the next temperature point
        end
    end
end
clearvars -except relaxE EQlat_mom lattice params
end