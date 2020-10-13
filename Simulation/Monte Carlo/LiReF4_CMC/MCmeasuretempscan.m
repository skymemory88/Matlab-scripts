function MCmeasuretempscan(filein_id)
store1 = ('G:\My Drive\File sharing\PhD program\Research projects\LiErF4 project\Quantum Monte Carlo\Test\');
% store2 = ('G:\My Drive\File sharing\PhD program\Research projects\LiErF4 project\Quantum Monte Carlo\Iteration limit');
EQname = sprintf(['resultsEQ_',filein_id,'.mat']);
EQs = fullfile(store1,EQname);
load(EQs);

% [relaxE,bestE,bestE2,C_v,~,mmagsq,msq0,~,~,mmag,malt,lat_mom,params] = paraTloop(ion,params,inter,lattice,EQlat_mom);
[relaxE,relaxE2,bestE,bestE2,acc_rate,C_v,malt,lat_mom,params,TT] = paraTloop(ion,params,inter,lattice,EQlat_mom);

% Gather all the data from different workers into arrays
TT = [TT{:}]';
bestE = [bestE{:}]';
bestE2 = [bestE2{:}]'; 
relaxE = {relaxE{:}};
relaxE2 = {relaxE2{:}};
malt = [malt{:}]';
C_v = [C_v{:}]';
lat_mom = [lat_mom{:}];
acc_rate = {acc_rate{:}};

% % order all the arrays' elements by parallel worker indices (core 1, core 2, ... core n, core 1, core 2 ... core n)
TT = TT(TT~=0);
bestE = bestE(:);
bestE2 = bestE2(:);
% relaxE = relaxE(:)';
% relaxE2 = relaxE2(:)';
malt = malt(:);
C_v = C_v(:);

% order all arrays' elements by asending temperature
% [~,sIdx] = sort(TT,'ascend');
% bestE = bestE(sIdx);
% bestE2 = bestE2(sIdx);
% % relaxE = relaxE(sIdx)';
% malt = malt(sIdx);
% C_v = C_v(sIdx);

% specific heat (derivative of the energy with respect to temperature) --Yikai
C_fdt = diff(bestE)./diff(TT);

% % susceptibility 
% chi_v = zeros(3,3,length(params.temp));
% chi = zeros(3,3,length(params.temp));
% for a=1:3
%     for b=1:3
%         for t=1:length(params.temp)
%             if(params.temp(t)~=0)
%                 chi_v(a,b,t)=(mean(mmagsq(a,b,:,t))-mean(mmag(:,a,t))*mean(mmag(:,b,t)))/params.temp(t);
%                 chi(a,b,t)=mean(mean(msq0(a,b,:,:,t)))/params.temp(t);
%             end
%         end
%     end
% end

cd('G:\My Drive\File sharing\PhD program\Research projects\LiErF4 project\Quantum Monte Carlo\Test');
filename=sprintf(['results_tempscan_',params.jobid,'.mat']);
% save(filename,'relaxE','bestE','bestE2','C_fdt','C_v','chi_v','chi','mmagsq','angl','msq0','msqx','msqy','mmag','malt','params','lat_mom','-v7.3');
save(filename,'relaxE','relaxE2','bestE','bestE2','acc_rate','C_v','C_fdt','malt','lat_mom','params','TT','-v7.3');
cd('G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\Monte Carlo\LiReF4_CMC');
clearvars
end