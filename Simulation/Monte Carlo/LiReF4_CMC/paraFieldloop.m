function [bestE,bestE2,C_fdt,angl,mmagsq,msq0,msqx,msqy,mmag,malt,lat_mom,params]=paraFieldloop(ion,params,inter,lattice,EQlat_mom)

temp = params.temp;
field = params.field;
L = params.L;
N=4*L^3;
iterLim = params.Nitermeas;

% angles in the XY plane
angl=zeros([size(field,1), size(lattice,2)]);

% spin moments at each temperature
lat_mom = cell(1,size(field,1));

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

% wforce = gcp(); % Start a new parallel pool if it doesn't exist
% para_opts = parforOptions(wforce, 'RangePartitionMethod','auto');

% parfor (ii=1:size(field,1), para_opts)
profile on
for ii = 1:size(field,1) % for debugging
    t=tic;
    field_change = true;
%     worker = getCurrentTask();
%     disp(['T = ',num2str(temp),' Field = [', num2str(field(ii,1)), ', ', num2str(field(ii,2)), ', ', num2str(field(ii,3)), '] ', num2str(worker.ID,'Worker ID: %u')]);
    disp(['T = ',num2str(temp),' Field = [', num2str(field(ii,1)), ', ', num2str(field(ii,2)), ', ', num2str(field(ii,3)), ']']);
    E_0=energy0(ion,params,inter,lattice,EQlat_mom{1,ii});
    
    [~,Egs,Egs2,mag,magsq,sq0,sqx,sqy,templattice,~,lat_mom{1,ii}]=LiIonsF4_MC(ion,L,iterLim,inter,lattice,field(ii,:)',temp,E_0,EQlat_mom{1,ii},field_change);
    
    % mean values over the MC steps
    bestE(ii)=mean(Egs);
    bestE2(ii)=mean(Egs2);
    mmagsq(:,:,:,ii)=mean(magsq,4);
    mmag(:,:,ii)=mean(mag,3);
    msq0(:,:,:,:,ii)=mean(sq0,5);
    msqx(:,:,:,:,ii)=mean(sqx,5);
    msqy(:,:,:,:,ii)=mean(sqy,5);
    
    % alternated magnetizations (vectorized from the original loop --Yikai)
    altx = sum(params.C(1:4,1).*squeeze(mag(1,1:4,:)));
    alty = sum(params.C(1:4,2).*squeeze(mag(2,1:4,:)));
    malt(ii)=mean(sqrt(altx.^2+alty.^2));
    
    % Specific heat from fluctuation-dissipation theorem
    if(temp~=0)
        C_fdt(ii)=(bestE2(ii)-(bestE(ii)^2))/((1/11.6)*(temp^2)); % kB~1/11.6 [meV/K]
    end
    
    for jj = 1:N
        angl(jj,ii)=atan2(templattice{jj}.mom(2),templattice{jj}.mom(1));
    end
%     angl(:,ii) = arrayfun(@(lat) atan2(lat.mom(2),lat.mom(1)),templattice); % Vectorized from original
    
    disp(['magnetX_alt = ',num2str(squeeze(mag(1,:,ii))*params.C(:,1))]);
    disp(['magnetY_alt = ',num2str(squeeze(mag(2,:,ii))*params.C(:,2))]);
    
    toc(t);
end
profile viewer
end