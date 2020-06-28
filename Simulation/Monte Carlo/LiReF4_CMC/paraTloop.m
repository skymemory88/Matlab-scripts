% function [relaxE,bestE,bestE2,C_v,angl,mmagsq,msq0,msqx,msqy,mmag,malt,lat_mom,params]=paraTloop(ion,params,inter,lattice,EQlat_mom)
function [relaxE,bestE,bestE2,C_v,malt,lat_mom,params,tempMark]=paraTloop(ion,params,inter,lattice,EQlat_mom)

temp = params.temp;

wforce = gcp(); % Start a new parallel pool if there isn't one
numWkrs = wforce.NumWorkers;
% fprintf('Total number of workers: %d.\n',numWkrs); %Checkpoint
% numWkrs = 1; % For debugging

if mod(size(temp,2),numWkrs) == 0
    block = int8(size(temp,2)/numWkrs);
%     block = size(temp,2); % for debugging
else
    error('Number of parallel threads mismatch the temperature dimension, parallel tempering will fail!');
end

% container for parallel tempering
index = 1;
half_interval = 50*params.N_Er;
tempMark = ones(ceil(params.Nitermeas/half_interval),1);
prob = rand(size(tempMark));

% % angles in the XY plane
% angl=zeros([size(temp,2),size(lattice,2)]);

% energy per spin, squared energy per spin and fluctuations
relaxE = double.empty(block,params.Nitermeas/params.N_Er,0);
bestE = double.empty(block,0);
bestE2 = double.empty(block,0);
C_v = double.empty(block,0);

% mean magnetization squared, magnetization and alternated magnetization
% mmagsq=zeros(3,3,4,size(temp,2));
% mmag=zeros(3,4,size(temp,2));
malt = double.empty(block,0);

% % form factors 
% msq0=zeros(3,3,4,4,size(temp,2));
% msqx=zeros(3,3,4,4,size(temp,2));
% msqy=zeros(3,3,4,4,size(temp,2));

field_change = params.field_changed;  % trigger for Ising_basis(), for tempscan, use only once
spmd
    tempMark = tempMark.*labindex;
    for ii=labindex:numWkrs:size(temp,2)
%     for ii=1:size(temp,2) % single thread for-loop for debugging
        t=tic;  
        fprintf('T = %1$.3f, Field = [%2$.2f, %3$.2f, %4$.2f]\n',temp(ii),params.field(1),params.field(2),params.field(3));
        
        E_0 = energy0(ion,params,inter,lattice,EQlat_mom{ii});
    %     [relaxE(ii,:),Egs,Egs2,mag,~,~,~,~,templattice,~,lat_mom{ii,1}]=LiIonsF4_MC(ion,params.L,iterLim,inter,lattice,params.field,temp(ii),E_0,EQlat_mom{ii},field_change);
    %     [relaxE(ii,:),Egs,Egs2,mag,~,lat_mom{ii,1}]=LiIonsF4_MC(ion,params,inter,lattice,temp(ii),E_0,EQlat_mom{ii},numWkrs,worker.ID,field_change);
        [relaxE(index,:,1),Egs,Egs2,mag,~,lat_mom,tempMark]=LiIonsF4_MC(ion,params,inter,lattice,temp(ii),E_0,EQlat_mom{ii},field_change,numWkrs,half_interval,tempMark,prob);    
        % mean values over the MC steps
        bestE(index,1) = mean(Egs);
        bestE2(index,1) = mean(Egs2);
        
%         bestE(ii)=mean(Egs);
%         bestE2(ii)=mean(Egs2);
    %     mmagsq(:,:,:,ii)=mean(magsq,4);
    %     mmag(:,:,ii)=mean(mag,3);
    %     msq0(:,:,:,:,ii)=mean(sq0,5);
    %     msqx(:,:,:,:,ii)=mean(sqx,5);
    %     msqy(:,:,:,:,ii)=mean(sqy,5);

        % alternating magnetizations
        altx = sum(params.C(1:4,1).*squeeze(mag(1,1:4,:)));
        alty = sum(params.C(1:4,2).*squeeze(mag(2,1:4,:)));
%         malt(ii)=mean(sqrt(altx.^2+alty.^2));
        malt(index,1) = mean(sqrt(altx.^2+alty.^2));

        if(temp(ii)~=0)
%             C_v(ii)=(bestE2(ii)-(bestE(ii)^2))/((1/11.6)*(temp(ii)^2)); % kB~1/11.6 [meV/K]
            C_v(index,1) = (bestE2(index,1)-(bestE(index,1)^2))/((1/11.6)*(temp(ii)^2)); % kB~1/11.6 [meV/K]
        end

    %     for jj=1:N
    %         angl(ii,jj)=atan2(templattice{jj}.mom(2),templattice{jj}.mom(1));
    %     end
    %     angl(:,ii) = arrayfun(@(lat) atan2(lat.mom(2),lat.mom(1)),templattice); % Vectorized from original

        disp(['magnetX_alt = ',num2str(squeeze(mag(1,:,end))*params.C(:,1))]);
        disp(['magnetY_alt = ',num2str(squeeze(mag(2,:,end))*params.C(:,2))]);

        index = index + 1;
        toc(t);

    end
end
clearvars -except relaxE bestE bestE2 C_v malt lat_mom params tempMark
end