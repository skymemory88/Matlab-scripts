function MCmeasuretempscan2(filein_id)
store1 = ('G:\My Drive\File sharing\PhD program\Research projects\LiErF4 project\Quantum Monte Carlo\Tc scaling\');
% store2 = ('G:\My Drive\File sharing\PhD program\Research projects\LiErF4 project\Quantum Monte Carlo\Iteration limit');
EQname = sprintf(['resultsEQ_',filein_id,'.mat']);
EQs = fullfile(store1,EQname);
load(EQs);
% cd('G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\External source\Monte Carlo\Parallel-CMC code');

% for t=1:length(params.temp) % Appears to be duplicated loop
%     [bestE(t),bestE2(t),C_fdt(t),angl(t,:),mmagsq(t,:,:,:),msq0(t,:,:,:,:),msqx(t,:,:,:,:),msqy(t,:,:,:,:),mmag(t,:,:),malt(t),lat_mom,params]=Tloop(ion,params,inter,lattice,EQlat_mom{t,1},temp(t));
% end

[bestE,bestE2,C_fdt,angl,mmagsq,msq0,msqx,msqy,mmag,malt,lat_mom,params]=paraTloop(ion,params,inter,lattice,EQlat_mom);

% specific heat (derivative of the energy with respect to temperature) --Yikai
C_v = diff(bestE)./diff(params.temp);

% susceptibility 
chi_fdt=zeros(3,3,length(params.temp));
chi=zeros(3,3,length(params.temp));
for a=1:3
    for b=1:3
        for t=1:length(params.temp)
            if(params.temp(t)~=0)
%                 chi_fdt(a,b,t)=(mean(mmagsq(t,a,b,:))-mean(mmag(t,:,a))*mean(mmag(t,:,b)))/params.temp(t);
%                 chi(a,b,t)=mean(mean(msq0(t,a,b,:,:)))/params.temp(t);
                chi_fdt(a,b,t)=(mean(mmagsq(a,b,:,t))-mean(mmag(:,a,t))*mean(mmag(:,b,t)))/params.temp(t);
                chi(a,b,t)=mean(mean(msq0(a,b,:,:,t)))/params.temp(t);
            end
        end
    end
end

cd('G:\My Drive\File sharing\PhD program\Research projects\LiErF4 project\Quantum Monte Carlo\Tc scaling (Ising Basis)');
filename=sprintf(['results_tempscan_',params.jobid,'.mat']);
save(filename,'bestE','bestE2','C_fdt','C_v','chi_fdt','chi','mmagsq','angl','msq0','msqx','msqy','mmag','malt','params','lat_mom');
cd('G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\External source\Monte Carlo\Parallel-CMC code');
end