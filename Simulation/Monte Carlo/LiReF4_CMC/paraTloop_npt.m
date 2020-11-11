% function [relaxE,bestE,bestE2,C_v,angl,mmagsq,msq0,msqx,msqy,mmag,malt,lat_mom,params]=paraTloop(ion,params,inter,lattice,EQlat_mom)
function [Egs,Egs2,bestE,bestE2,acc_rate,C_v,malt,lat_mom,params,tempf,meas]=paraTloop(ion,params,inter,lattice,EQlat_mom)

temp = params.temp;

wforce = gcp(); % Start a new parallel pool if there isn't one
numWkrs = wforce.NumWorkers;
fprintf('Total number of workers: %d.\n',numWkrs); %Checkpoint
% numWkrs = 1; % For debugging

if mod(size(temp,2),numWkrs) == 0
    block = int8(size(temp,2)/numWkrs);
%     block = size(temp,2); % for debugging
else
    error('Number of parallel threads mismatch the temperature dimension, parallel tempering will fail!');
end

% container for parallel tempering
index = 1;
acc_rate = double.empty(block,params.Nitermeas,0); % container for acceptance rate
tempMark = ones(block,params.Nitermeas+1);
% marker = double.empty(block,int8(params.Nitermeas/params.pt_intv)+1,numWkrs,0);
prob = rand(size(tempMark));

% % angles in the XY plane
% angl=zeros([size(temp,2),size(lattice,2)]);

% energy per spin, squared energy per spin and fluctuations
bestE = double.empty(block,0);
bestE2 = double.empty(block,0);
C_v = double.empty(block,0);
% mean magnetization squared, magnetization and alternated magnetization
% mmagsq=zeros(3,3,4,size(temp,2));
% mmag=zeros(3,4,size(temp,2));
malt = double.empty(block,0);
tempf = double.empty(block,0);

% % form factors 
% msq0=zeros(3,3,4,4,size(temp,2));
% msqx=zeros(3,3,4,4,size(temp,2));
% msqy=zeros(3,3,4,4,size(temp,2));

field_change = params.field_changed;  % trigger for Ising_basis(), for tempscan, use only once
spmd
    for ii=1:block
%     for ii=1:size(temp,2) % single thread for-loop for debugging
        t=tic;  
        % the temperature assignment is such that each core moves forward by one point but separated from each other by a stride of the size of the total number of cores.
        fprintf('T = %1$.3f, Field = [%2$.2f, %3$.2f, %4$.2f]\n',temp(ii+(labindex-1)*block),params.field(1),params.field(2),params.field(3));
        tempMark(index,1) = labindex;
        E_0 = energy0(ion,params,inter,lattice,EQlat_mom{ii+(labindex-1)*block});
        disp('Initial energy calculated.')
    %     [relaxE(ii+(labindex-1)*block,:),Egs,Egs2,mag,~,~,~,~,templattice,~,lat_mom{ii+(labindex-1)*block,1}]=LiIonsF4_MC(ion,params.L,iterLim,inter,lattice,params.field,temp(ii+(labindex-1)*block),E_0,EQlat_mom{ii+(labindex-1)*block},field_change);
    %     [relaxE(ii+(labindex-1)*block,:),Egs,Egs2,mag,~,lat_mom{ii+(labindex-1)*block,1}]=LiIonsF4_MC(ion,params,inter,lattice,temp(ii+(labindex-1)*block),E_0,EQlat_mom{ii+(labindex-1)*block},numWkrs,worker.ID,field_change);
        [Egs,Egs2,acc_rate(index,:,1),mag,~,lat_mom,tempMark(index,:),tempf(index,1),mark,meas]=LiIonsF4_MC(ion,params,inter,lattice,temp(ii+(labindex-1)*block),E_0,EQlat_mom{ii+(labindex-1)*block},field_change,numWkrs,tempMark(index,:),prob(index,:));    
        % mean values over the MC steps
        bestE(index,1) = mean(Egs);
        bestE2(index,1) = mean(Egs2);
        Egs = Egs';
        Egs2 = Egs2';
        
%         bestE(ii+(labindex-1)*block)=mean(Egs);
%         bestE2(ii+(labindex-1)*block)=mean(Egs2);
%         mmagsq(:,:,:,ii+(labindex-1)*block)=mean(magsq,4);
%         mmag(:,:,ii+(labindex-1)*block)=mean(mag,3);
%         msq0(:,:,:,:,ii+(labindex-1)*block)=mean(sq0,5);
%         msqx(:,:,:,:,ii+(labindex-1)*block)=mean(sqx,5);
%         msqy(:,:,:,:,ii+(labindex-1)*block)=mean(sqy,5);

        % alternating magnetizations
        altx = sum(params.C(1:4,1).*squeeze(mag(1,1:4,:)));
        alty = sum(params.C(1:4,2).*squeeze(mag(2,1:4,:)));
%         malt(ii+(labindex-1)*block)=mean(sqrt(altx.^2+alty.^2));
        malt(index,1) = mean(sqrt(altx.^2+alty.^2));

        if(temp(ii+(labindex-1)*block)~=0)
%             C_v(ii+(labindex-1)*block)=(bestE2(ii+(labindex-1)*block)-(bestE(ii+(labindex-1)*block)^2))/((1/11.6)*(temp(ii+(labindex-1)*block)^2)); % kB~1/11.6 [meV/K]
            C_v(index,1) = params.N_Er*(bestE2(index,1)-(bestE(index,1)^2))/((1/11.6)*(temp(ii+(labindex-1)*block)^2)); % kB~1/11.6 [meV/K]
        end

    %     for jj=1:N
    %         angl(ii,jj)=atan2(templattice{jj}.mom(2),templattice{jj}.mom(1));
    %     end
    %     angl(:,ii) = arrayfun(@(lat) atan2(lat.mom(2),lat.mom(1)),templattice); % Vectorized from original

        disp(['magnetX_alt = ',num2str(squeeze(mag(1,:,end))*params.C(:,1))]);
        disp(['magnetY_alt = ',num2str(squeeze(mag(2,:,end))*params.C(:,2))]);

        index = index + 1;
        toc(t);
%         fprintf('Time of parallel tempering attempt: %d.\n', mark)
    end
end

% for ii = 1:size(marker,3)
%     marker(:,:,ii,1) = tempMark{ii};
% end

clearvars -except Egs Egs2 acc_rate bestE bestE2 C_v malt lat_mom params tempf mark meas
end