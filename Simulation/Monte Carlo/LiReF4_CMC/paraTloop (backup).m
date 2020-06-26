% function [relaxE,bestE,bestE2,C_v,angl,mmagsq,msq0,msqx,msqy,mmag,malt,lat_mom,params]=paraTloop(ion,params,inter,lattice,EQlat_mom)
function [relaxE,bestE,bestE2,C_v,malt,lat_mom,params]=paraTloop(ion,params,inter,lattice,EQlat_mom)

temp = params.temp;

% % angles in the XY plane
% angl=zeros([size(temp,2),size(lattice,2)]);

% energy per spin, squared energy per spin and fluctuations
relaxE=zeros(size(params.temp,2),params.Nitermeas/params.N_Er);
bestE=zeros(1,size(temp,2));
bestE2=zeros(1,size(temp,2));
C_v=zeros(1,size(temp,2));

% % mean magnetization squared, magnetization and alternated magnetization
% mmagsq=zeros(3,3,4,size(temp,2));
% mmag=zeros(3,4,size(temp,2));
% malt=zeros(1,size(temp,2));
% 
% % form factors 
% msq0=zeros(3,3,4,4,size(temp,2));
% msqx=zeros(3,3,4,4,size(temp,2));
% msqy=zeros(3,3,4,4,size(temp,2));

wforce = gcp(); % Start a new parallel pool if it doesn't exist
numWkrs = wforce.NumWorkers;
para_opts = parforOptions(wforce, 'RangePartitionMethod','auto');

field_change = params.field_changed;  % trigger for Ising_basis(), for tempscan, use only once
parfor (ii=1:size(temp,2),para_opts)
% parfor (ii=5:size(temp,2)-6,para_opts)
% for ii=1:size(temp,2) % single thread for-loop for debugging
    t=tic;  
    worker = getCurrentTask();
    disp(['T = ',num2str(temp(ii)), ', Field = [', num2str(params.field(1)), ', ', num2str(params.field(2)), ', ', num2str(params.field(3)), '] ', num2str(worker.ID,'Worker ID: %u')]);
  
    E_0=energy0(ion,params,inter,lattice,EQlat_mom{ii});
%     [relaxE(ii,:),Egs,Egs2,mag,~,~,~,~,templattice,~,lat_mom{ii,1}]=LiIonsF4_MC(ion,params.L,iterLim,inter,lattice,params.field,temp(ii),E_0,EQlat_mom{ii},field_change);
%     [relaxE(ii,:),Egs,Egs2,mag,~,lat_mom{ii,1}]=LiIonsF4_MC(ion,params,inter,lattice,temp(ii),E_0,EQlat_mom{ii},numWkrs,worker.ID,field_change);
    [relaxE(ii,:),Egs,Egs2,mag,~,lat_mom{ii,1}]=LiIonsF4_MC(ion,params,inter,lattice,temp(ii),E_0,EQlat_mom{ii},numWkrs,worker.ID,field_change);    
    % mean values over the MC steps
    bestE(ii)=mean(Egs);
    bestE2(ii)=mean(Egs2);
%     mmagsq(:,:,:,ii)=mean(magsq,4);
%     mmag(:,:,ii)=mean(mag,3);
%     msq0(:,:,:,:,ii)=mean(sq0,5);
%     msqx(:,:,:,:,ii)=mean(sqx,5);
%     msqy(:,:,:,:,ii)=mean(sqy,5);
    
    % alternating magnetizations
    altx = sum(params.C(1:4,1).*squeeze(mag(1,1:4,:)));
    alty = sum(params.C(1:4,2).*squeeze(mag(2,1:4,:)));
    malt(ii)=mean(sqrt(altx.^2+alty.^2));
    
    if(temp(ii)~=0)
        C_v(ii)=(bestE2(ii)-(bestE(ii)^2))/((1/11.6)*(temp(ii)^2)); % kB~1/11.6 [meV/K]
    end
    
%     for jj=1:N
%         angl(ii,jj)=atan2(templattice{jj}.mom(2),templattice{jj}.mom(1));
%     end
%     angl(:,ii) = arrayfun(@(lat) atan2(lat.mom(2),lat.mom(1)),templattice); % Vectorized from original
    
    disp(['magnetX_alt = ',num2str(squeeze(mag(1,:,end))*params.C(:,1))]);
    disp(['magnetY_alt = ',num2str(squeeze(mag(2,:,end))*params.C(:,2))]);
    
    toc(t);
  
end

end