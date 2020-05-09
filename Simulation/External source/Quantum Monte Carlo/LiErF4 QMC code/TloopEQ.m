function [relaxE,EQlat_mom,lattice,params]=TloopEQ(ion,params,inter,lattice,lat_mom)

temp=params.temp;
field=params.field;
relaxE=zeros([size(temp,2),params.NiterEQ]);
N=4*(params.L)^3;

E_0=energy0(ion,params,inter,lattice,lat_mom);
disp(['E_0 = ',num2str(E_0/N)]);

EQlat_mom=cell(length(temp),1);

for i=1:size(temp,2)
    t=tic;
    
    for k=1:N %loop added in reference to FieldloopEQ.m --Yikai
        [lattice{k}.VV,~]=Ising_basis(field,lattice{k});
    end
    disp(['T = ',num2str(temp(i)),' Field = [',num2str(field(1)),',',num2str(field(2)),',',num2str(field(3)),']']);
    [energies,lattice,E_0,lat_mom,params]=LiIonsF4_MCEQ(ion,params,inter,lattice,field,temp(i),E_0,lat_mom);
    EQlat_mom{i}=lat_mom;
    relaxE(i,:)=energies;       

    toc(t);

end

end