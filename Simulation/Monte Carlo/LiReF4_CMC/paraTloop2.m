% function [relaxE,bestE,bestE2,C_v,angl,mmagsq,msq0,msqx,msqy,mmag,malt,lat_mom,params]=paraTloop(ion,params,inter,lattice,EQlat_mom)
function [Egs,Egs2,bestE,bestE2,acc_rate,C_v,malt,lat_mom,params,tempf]=paraTloop2(ion,params,inter,lattice,EQlat_mom)

temp = params.temp;
N_meas = params.Nitermeas;
N_stats = params.Niterstat;
wforce = gcp(); % Start a new parallel pool if there isn't one
numWkrs = wforce.NumWorkers;
% numWkrs = 1; % For debugging
fprintf('Total number of workers: %d.\n',numWkrs); %Checkpoint

% Check point to make sure each worker gets equal amount of tasks
if mod(size(temp,2),numWkrs) == 0
    block = int8(size(temp,2)/numWkrs);
%     block = size(temp,2); % for debugging
else
    error('Number of parallel threads mismatch the temperature dimension, parallel tempering will fail!');
end

% container for parallel tempering
index = 1;
acc_rate = double.empty(block,params.meas_intv,N_meas,0); % container for acceptance rate
tempMark = ones(block,ceil(N_meas * params.meas_intv / params.pt_intv)+1);
% marker = double.empty(block,int8(params.Nitermeas/params.pt_intv)+1,numWkrs,0);
prob = rand(size(tempMark));

% % angles in the XY plane
% angl=zeros([size(temp,2),size(lattice,2)]);

% energy per spin, squared energy per spin and fluctuations
Egs = double.empty(block,N_meas,0);
Egs2 = double.empty(block,N_meas,0);
bestE = double.empty(block,0); % Ground state energy
bestE2 = double.empty(block,0); 
C_v = double.empty(block,0); % Specific heat from F-D theorem

% % mean magnetization squared, magnetization and alternated magnetization
% mmagsq = zeros(3,3,4,size(temp,2));
% mmag = zeros(3,4,size(temp,2));
malt = double.empty(block,0); % Alternating moments as order parameter
tempf = double.empty(block,0);

% % form factors 
% msq0 = zeros(3,3,4,4,size(temp,2));
% msqx = zeros(3,3,4,4,size(temp,2));
% msqy = zeros(3,3,4,4,size(temp,2));

field_change = params.field_changed;  % trigger for Ising_basis(), for tempscan, use only once
spmd
    for ii = 1:block
%     for ii = 1:size(temp,2) % single thread for-loop for debugging
%         t=tic;  
        % the temperature assignment is such that each core moves forward by one point but separated from each other by a stride of the size of the total number of cores.
        fprintf('T = %1$.3f, Field = [%2$.2f, %3$.2f, %4$.2f]\n',temp(ii+(labindex-1)*block),params.field(1),params.field(2),params.field(3));
        tempMark(index,1) = labindex;
        E_0 = energy0(ion,params,inter,lattice,EQlat_mom{ii+(labindex-1)*block});
%         disp('Initial energy calculated.')

        % containers for statistical points at each measurement
        altx_temp = zeros(1,N_stats,N_meas);
        alty_temp = zeros(1,N_stats,N_meas);
        lat_mom = EQlat_mom{ii+(labindex-1)*block};
        for jj = 1:N_meas
            [E_0,~,~,acc_rate(index,:,jj,1),~,lattice,lat_mom,~,temp(ii+(labindex-1)*block),~] = rand_walk(params.meas_intv,params.pt_intv,ion,params,inter,lattice,temp(ii+(labindex-1)*block),E_0,lat_mom,field_change,numWkrs,tempMark(index,:),prob(index,:));
            [E_0,Egs(index,jj,1),Egs2(index,jj,1),~,mag,lattice,lat_mom,~,temp(ii+(labindex-1)*block),~] = rand_walk(N_stats,inf,ion,params,inter,lattice,temp(ii+(labindex-1)*block),E_0,lat_mom,field_change,numWkrs,tempMark(index,:),prob(index,:));     
            altx_temp(1,:,jj) = sum(params.C(1:4,1).*permute(mag(1,1:4,:),[2 3 1])); % alternating magnetizations in x
            alty_temp(1,:,jj) = sum(params.C(1:4,2).*permute(mag(2,1:4,:),[2 3 1])); % alternating magnetizations in y
        end
        tempf(index,1) = temp(ii+(labindex-1)*block);
        bestE(index,1) = mean(Egs(index,:,1),2); % Ground state energy (mean values over the MC steps)
        bestE2(index,1) = mean(Egs2(index,:,1),2);
%         Egs = squeeze(Egs)';
%         Egs2 = squeeze(Egs2)';
%         Egs = Egs';
%         Egs2 = Egs2';
        acc_rate = mean(squeeze(acc_rate),3);
        altx = squeeze(sum(altx_temp,3)); % alternating magnetizations in x
        alty = squeeze(sum(alty_temp,3)); % alternating magnetizations in y
%         malt(ii+(labindex-1)*block) = mean(sqrt(altx.^2+alty.^2));
        malt(index,1) = mean(sqrt(altx.^2+alty.^2));

        if temp(ii+(labindex-1)*block) ~= 0
%             C_v(ii+(labindex-1)*block) = (bestE2(ii+(labindex-1)*block)-(bestE(ii+(labindex-1)*block)^2))/((1/11.6)*(temp(ii+(labindex-1)*block)^2)); % kB~1/11.6 [meV/K]
            C_v(index,1) = params.N_Er*(bestE2(index,1)-(bestE(index,1)^2))/((1/11.6)*(temp(ii+(labindex-1)*block)^2)); % kB~1/11.6 [meV/K]
        end

    %     for jj=1:N
    %         angl(ii,jj) = atan2(templattice{jj}.mom(2),templattice{jj}.mom(1));
    %     end
    %     angl(:,ii) = arrayfun(@(lat) atan2(lat.mom(2),lat.mom(1)),templattice); % Vectorized from original

        disp(['magnetX_alt = ',num2str(squeeze(mag(1,:,end))*params.C(:,1))]);
        disp(['magnetY_alt = ',num2str(squeeze(mag(2,:,end))*params.C(:,2))]);

        index = index + 1;
%         toc(t);
%         fprintf('Time of parallel tempering attempt: %d.\n', mark)
    end
end

% for ii = 1:size(marker,3)
%     marker(:,:,ii,1) = tempMark{ii};
% end

clearvars -except Egs Egs2 acc_rate bestE bestE2 C_v malt lat_mom params tempf
end