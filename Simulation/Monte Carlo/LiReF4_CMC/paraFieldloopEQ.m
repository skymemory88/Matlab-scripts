function [relaxE,EQlat_mom,lattice,params]=paraFieldloopEQ(ion,params,inter,lattice,lat_mom)

temp = params.temp;
field = params.field;
relaxE=zeros([size(field,1),params.NiterEQ]);
L = params.L;
N=4*(L)^3;
EQiterLim = params.NiterEQ;

% Calculate the energy of the current configuration as the initial state for all field points
E_0=energy0(ion,params,inter,lattice,lat_mom);
disp(['E_0=',num2str(E_0/N)]);
EQlat_mom=cell(1,size(field,1));

wforce = gcp(); % Start a new parallel pool if it doesn't exist
para_opts = parforOptions(wforce, 'RangePartitionMethod','auto');

parfor (jj = 1:size(field,1),para_opts)
    field_change = true;
    worker = getCurrentTask();
    disp(['T = ',num2str(temp),' Field = [',num2str(field(jj,1)),',',num2str(field(jj,2)),',',num2str(field(jj,3)),'], ', num2str(worker.ID,'Worker ID: %u')]);
    [energies,~,~,EQlat_mom{1,jj}]=LiIonsF4_MCEQ(ion,L,EQiterLim,inter,lattice,field(jj,:),temp,E_0,lat_mom,field_change);
%     EQlat_mom{1,jj}=new_lat_mom;
    relaxE(jj,:)=energies;
end

end