function [relaxE,EQlat_mom,lattice,params]=paraTloopEQ(ion,params,inter,lattice,lat_mom)

temp = params.temp;
field = params.field;
relaxE=zeros([size(temp,2),params.NiterEQ]);
L = params.L;
N=4*(L)^3;
EQiterLim = params.NiterEQ;

E_0=energy0(ion,params,inter,lattice,lat_mom);
disp(['E_0 = ',num2str(E_0/N)]);
EQlat_mom=cell(length(temp),1);

wforce = gcp(); % Start a new parallel pool if it doesn't exist
para_opts = parforOptions(wforce, 'RangePartitionMethod','auto');

parfor (ii = 1:size(temp,2), para_opts)
    field_change = true; % trigger for Ising_basis()
    worker = getCurrentTask();
    disp(['T = ',num2str(temp(ii)),' Field = [',num2str(field(1)),',',num2str(field(2)),',',num2str(field(3)),'] ', num2str(worker.ID,'Worker ID: %u')]);
    [energies,~,~,new_lat_mom]=LiIonsF4_MCEQ(ion,L,EQiterLim,inter,lattice,field,temp(ii),E_0,lat_mom,field_change);
    EQlat_mom{ii}=new_lat_mom;
    relaxE(ii,:)=energies;
end

end