function [bestE,bestE2,C_fdt,angl,mmagsq,msq0,msqx,msqy,mmag,malt,lat_mom,params]=Tloop(ion,params,inter,lattice,EQlat_mom)

temp = params.temp;
field = params.field;

% angles in the XY plane
angl=zeros([size(temp,2),size(lattice,2)]);

% energy per spin, squared energy per spin and fluctuations
bestE=zeros(1,size(temp,2));
bestE2=zeros(1,size(temp,2));
C_fdt=zeros(1,size(temp,2));

% mean magnetization squared, magnetization and alternated magnetization
mmagsq=zeros(3,3,4,size(temp,2));
mmag=zeros(3,4,size(temp,2));
malt=zeros(1,size(temp,2));

% form factors 
msq0=zeros(3,3,4,4,size(temp,2));
msqx=zeros(3,3,4,4,size(temp,2));
msqy=zeros(3,3,4,4,size(temp,2));

% E_0=energy0(ion,params,inter,lattice,EQlat_mom); % Original script
% --Yikai
% disp(['E_0 = ',num2str(E_0/(4*(params.L)^3))]);

for i=1:size(temp,2)
    t=tic;
    
    E_0=energy0(ion,params,inter,lattice,EQlat_mom{i,1});
    disp(['E_0 = ',num2str(E_0/(4*(params.L)^3))]);

    for k=1:length(lattice) % Loop added in reference to Fieldloop.m --Yikai
        [lattice{k}.VV,~]=Ising_basis(field,lattice{k});
    end
    
    disp(['T = ',num2str(temp(i)), ', Field = [', num2str(field(1)), ', ', num2str(field(2)), ', ', num2str(field(3)),']']);
%     [~,Egs,Egs2,mag,magsq,sq0,sqx,sqy,lattice,E_0,lat_mom,params]=LiIonsF4_MC(ion,params,inter,lattice,field,temp(i),E_0,EQlat_mom);
    [~,Egs,Egs2,mag,magsq,sq0,sqx,sqy,lattice,~,lat_mom,params]=LiIonsF4_MC(ion,params,inter,lattice,field,temp(i),E_0,EQlat_mom{i,1});
    
    % mean values over the MC steps
    bestE(i)=mean(Egs);
    bestE2(i)=mean(Egs2);
    mmagsq(:,:,:,i)=mean(magsq,4);
    mmag(:,:,i)=mean(mag,3);
    msq0(:,:,:,:,i)=mean(sq0,5);
    msqx(:,:,:,:,i)=mean(sqx,5);
    msqy(:,:,:,:,i)=mean(sqy,5);
    
    % alternated magnetizations
    altx=zeros(size(mag,3),1);
    alty=zeros(size(mag,3),1);
    for m=1:4
        altx=altx+params.C(m,1)*squeeze(mag(1,m,:));
        alty=alty+params.C(m,2)*squeeze(mag(2,m,:));
    end
    malt(i)=mean(sqrt(altx.^2+alty.^2));
    
    
    if(temp(i)~=0)
        C_fdt(i)=(bestE2(i)-(bestE(i)^2))/((1/11.6)*(temp(i)^2));
    end
    
    for j=1:size(lattice,2)
        angl(i,j)=atan2(lattice{j}.mom(2),lattice{j}.mom(1));
    end
    
    disp(['magnetX_alt = ',num2str(squeeze(mag(1,:,i))*params.C(:,1))]);
    disp(['magnetY_alt = ',num2str(squeeze(mag(2,:,i))*params.C(:,2))]);
    disp(['E = ',num2str(bestE(i))]);
    
    toc(t);
  
end

end