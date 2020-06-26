function [relaxE,EQlat_mom,lattice,params]=TloopEQ(ion,params,inter,lattice,lat_mom)

N=4*(params.L)^3;

% Initiate containers for equilibrium configurations as well as equilibrium energies
relaxE=zeros(size(params.temp,2),params.NiterEQ/params.N_Er);
EQlat_mom=cell(length(params.temp),1);

% Calculate the energy of the current configuration as the initial state for all field points
E_0 = energy0(ion,params,inter,lattice,lat_mom);
disp(['E_0/N = ',num2str(E_0/N)]);

% Seed the first (T = min(T)) thermalization process with the initial configuration
relaxE(1,:) = E_0/N;
EQlat_mom{1,1} = lat_mom;

% Carry out thermalization
for ii = 1:size(params.temp,2)
    field_change = false; % trigger for Ising_basis()
    disp(['T = ',num2str(params.temp(ii)),' Field = [',num2str(params.field(1)),',',num2str(params.field(2)),',',num2str(params.field(3)),']']);
    [relaxE(ii,:),lattice,~,EQlat_mom{ii,1}]=LiIonsF4_MCEQ(ion,params.L,params.NiterEQ,inter,lattice,params.field,params.temp(ii),relaxE(ii,end)*N,EQlat_mom{ii,1},field_change);
    if ii < size(params.temp,2)
        relaxE(ii+1,:) = relaxE(ii,:); % Seed the thermalization process of the next temperature point
        EQlat_mom{ii+1,1} = EQlat_mom{ii,1}; % Seed the tehermalization process of the next temperature point
    end
end
clearvars -except relaxE EQlat_mom lattice params
end