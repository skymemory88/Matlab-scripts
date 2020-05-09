function [relaxE,EQlat_mom,lattice,params]=FieldloopEQ(ion,params,inter,lattice,lat_mom)

temp=params.temp;
field=params.field;
relaxE=zeros([size(field,1),params.NiterEQ]);
N=4*(params.L)^3;

E_0=energy0(ion,params,inter,lattice,lat_mom);
disp(['E_0=',num2str(E_0/N)]);

EQlat_mom=cell(1,size(field,1));

for j=1:size(field,1)
    t=tic;
    
    for k=1:N
        [lattice{k}.VV,~]=Ising_basis(field(j,:),lattice{k});
    end
    disp(['T = ',num2str(temp),' Field = [',num2str(field(j,1)),',',num2str(field(j,2)),',',num2str(field(j,3)),']']);
    [energies,lattice,E_0,lat_mom,params]=LiIonsF4_MCEQ(ion,params,inter,lattice,field(j,:),temp,E_0,lat_mom);
    EQlat_mom{1,j}=lat_mom;
    relaxE(j,:)=energies;
    
    toc(t);
end


end