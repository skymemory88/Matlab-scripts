function ContMCtempscan(filein_id,start,finish)
store1 = ('G:\My Drive\File sharing\PhD program\Research projects\LiErF4 project\Quantum Monte Carlo\Test\');
% store2 = ('G:\My Drive\File sharing\PhD program\Research projects\LiErF4 project\Quantum Monte Carlo\Iteration limit');
EQname1 = sprintf(['results_tempscan_',filein_id, num2str(start,'_%d.mat')]);
EQname2 = sprintf(['resultsEQ_',filein_id,'.mat']);
EQs1 = fullfile(store1,EQname1);
EQs2 = fullfile(store1,EQname2);
load(EQs1,'-mat','lat_mom','params','bestE','bestE2');
load(EQs2,'-mat','ion','inter','lattice');

[relaxE,bbestE,bbestE2,C_v,malt,lat_mom,params] = paraTloop(ion,params,inter,lattice,lat_mom);

bestE = (bbestE + bestE)./2;
bestE2 = (bbestE2 + bestE2)./2;

% specific heat (derivative of the energy with respect to temperature) --Yikai
C_fdt = diff(bestE)./diff(params.temp);

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
filename=sprintf(['results_tempscan_',params.jobid, num2str(finish,'_%d.mat')]);
% save(filename,'relaxE','bestE','bestE2','C_fdt','C_v','chi_v','chi','mmagsq','angl','msq0','msqx','msqy','mmag','malt','params','lat_mom','-v7.3');
save(filename,'relaxE','bestE','bestE2','C_fdt','C_v','params','malt','lat_mom','-v7.3');
cd('G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\Monte Carlo\LiReF4_CMC');
clearvars
end