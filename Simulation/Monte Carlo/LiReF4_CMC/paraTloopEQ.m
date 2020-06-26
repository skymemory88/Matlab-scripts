function [frelaxE,fEQlat_mom,lattice,params]=paraTloopEQ(ion,params,inter,lattice,lat_mom)

N = params.N_Er;
wforce = gcp(); % Start a new parallel pool if it doesn't exist
threads = wforce.NumWorkers;

% Initiate containers for equilibrium configurations as well as equilibrium energies
relaxE = zeros(size(params.temp,2),params.NiterEQ/N);
EQlat_mom = cell(1,size(params.temp,2));
EQlat_mom{1} = lat_mom;

% Initate containers for the final configurations
frelaxE = zeros(size(params.temp,2),params.NiterEQ/N);
fEQlat_mom = cell(1,size(params.temp,2));

% Calculate the energy of the current configuration as the initial state for all field points
E_0 = energy0(ion,params,inter,lattice,lat_mom);
disp(['E_0 = ',num2str(E_0/N)]);

% Seed the first (T = min(T)) thermalization process with the initial configuration
relaxE(1,:) = E_0;

% Carry out thermalization in parallel across all cores and save the data in the end
spmd
    for ii = 1:size(params.temp,2)
        field_change = false; % trigger for Ising_basis()
        if labindex == 1
            disp(['T = ',num2str(params.temp(ii)),' Field = [',num2str(params.field(1)),',',num2str(params.field(2)),',',num2str(params.field(3)),']']);
        end
        [relaxE(ii,:),~,~,EQlat_mom{ii}]=LiIonsF4_MCEQ(ion,params.L,params.NiterEQ,inter,lattice,params.field,params.temp(ii),relaxE(ii,end)*N,EQlat_mom{ii},field_change);
        if ii < size(params.temp,2)
            relaxE(ii+1,:) = relaxE(ii,:); % Seed the thermalization process of the next temperature point
            EQlat_mom{ii+1} = EQlat_mom{ii}; % Seed the tehermalization process of the next temperature point
        end
    end
end

% collect the results from all parallel threads
relaxE = {relaxE{:}};
EQlat_mom = {EQlat_mom{:}};

% pick the configuration of each temperature from the parallel pool that has the lowest final total energy
for ii = 1:size(params.temp,2)
    finals = double.empty(0,threads);
    for kk = 1:threads
        finals(1,kk) = relaxE{kk}(ii,end);
    end
    [~,idx] = min(finals);
    frelaxE(ii,:) = relaxE{idx}(ii,:);
    fEQlat_mom{ii} = EQlat_mom{idx}(ii);
end
    clearvars -except params lattice fEQlat_mom frelaxE N
end