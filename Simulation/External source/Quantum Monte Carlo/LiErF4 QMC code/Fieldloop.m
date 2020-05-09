function [bestE,bestE2,C_fdt,angl,mmagsq,msq0,msqx,msqy,mmag,malt,lat_mom,params]=Fieldloop(ion,params,inter,lattice,EQlat_mom)

temp = params.temp;
field = params.field;

% angles in the XY plane
angl=zeros([size(field,1),size(lattice,2)]);

% energy per spin, squared energy per spin and fluctuations
bestE=zeros(1,size(field,1));
bestE2=zeros(1,size(field,1));
C_fdt=zeros(1,size(field,1));

% mean magnetization squared, magnetization and alternated magnetization
mmagsq=zeros(3,3,4,size(field,1));
mmag=zeros(3,4,size(field,1));
malt=zeros(1,size(field,1));

% form factors 
msq0=zeros(3,3,4,4,size(field,1));
msqx=zeros(3,3,4,4,size(field,1));
msqy=zeros(3,3,4,4,size(field,1));

% E_0=energy0(ion,params,inter,lattice,EQlat_mom);
% disp(['E_0 = ',num2str(E_0/(4*(params.L)^3))]);
    
for i=1:size(field,1)
    t=tic;
    
    E_0=energy0(ion,params,inter,lattice,EQlat_mom{1,i});
    disp(['E_0 = ',num2str(E_0/(4*(params.L)^3))]);
    
    for k=1:length(lattice)
        [lattice{k}.VV,~]=Ising_basis(field(i,:),lattice{k});
    end
    
    disp(['T = ',num2str(temp),' Field = [', num2str(field(i,1)), ', ', num2str(field(i,2)), ', ', num2str(field(i,3)), ']']);
    [~,Egs,Egs2,mag,magsq,sq0,sqx,sqy,lattice,~,lat_mom,params]=LiIonsF4_MC(ion,params,inter,lattice,field(i,:)',temp,E_0,EQlat_mom{1,i});
    
    % mean values over the MC steps
    bestE(i)=mean(Egs);
    bestE2(i)=mean(Egs2);
    mmagsq(:,:,:,i)=mean(magsq,4);
    mmag(:,:,i)=mean(mag,3);
    msq0(:,:,:,:,i)=mean(sq0,5);
    msqx(:,:,:,:,i)=mean(sqx,5); %Originally commented out
    msqy(:,:,:,:,i)=mean(sqy,5); %Originally commented out
    
    % alternated magnetizations
    altx=zeros(size(mag,3),1);
    alty=zeros(size(mag,3),1);
    for m=1:4
        altx=altx+params.C(m,1)*squeeze(mag(1,m,:));
        alty=alty+params.C(m,2)*squeeze(mag(2,m,:));
    end
    malt(i)=mean(sqrt(altx.^2+alty.^2));
    
    
    if(temp~=0)
        C_fdt(i)=(bestE2(i)-(bestE(i)^2))/((1/11.6)*(temp^2));
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