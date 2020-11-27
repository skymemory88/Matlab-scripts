function MCmeasuretempscan(filein_id,N_meas,N_stats,meas_intv,pt_intv)
store1 = ('G:\My Drive\File sharing\PhD program\Research projects\LiErF4 project\Quantum Monte Carlo\Test\');
% store2 = ('G:\My Drive\File sharing\PhD program\Research projects\LiErF4 project\Quantum Monte Carlo\Iteration limit');
EQname = sprintf(['resultsEQ_',filein_id,'.mat']);
EQs = fullfile(store1,EQname);
load(EQs);

t0 = tic;
params.Nitermeas = N_meas; % Number of measurements
params.meas_intv = meas_intv; % interval size in between measurements
params.Niterstat = N_stats; % Number of statistical points for each measurement
params.pt_intv = pt_intv; % Interval size in between parallel tempering attempts

% [relaxE,bestE,bestE2,C_v,~,mmagsq,msq0,~,~,mmag,malt,lat_mom,params] = paraTloop(ion,params,inter,lattice,EQlat_mom);
[relaxE,relaxE2,bestE,bestE2,acc_rate,Marker,C_v,C_fdt,malt,lat_mom,params] = paraTloop3(ion,params,inter,lattice,EQlat_mom);

% % susceptibility 
% chi_v = zeros(3,3,length(params.temp));
% chi = zeros(3,3,length(params.temp));
% for a=1:3
%     for b=1:3
%         for t=1:length(params.temp)
%             if(params.temp(t)~=0)
%                 chi_v(a,b,t) = (mean(mmagsq(a,b,:,t))-mean(mmag(:,a,t))*mean(mmag(:,b,t)))/params.temp(t);
%                 chi(a,b,t) = mean(mean(msq0(a,b,:,:,t)))/params.temp(t);
%             end
%         end
%     end
% end

filepath = 'G:\My Drive\File sharing\PhD program\Research projects\LiErF4 project\Quantum Monte Carlo\Test';
filename = sprintf(['results_tempscan_',params.jobid,'_%1$u_%2$u_%3$u_%4$u.mat'],params.Nitermeas, params.Niterstat, params.meas_intv, params.pt_intv);
fileobj = fullfile(filepath,filename);
% save(fileobj,'relaxE','bestE','bestE2','C_fdt','C_v','chi_v','chi','mmagsq','angl','msq0','msqx','msqy','mmag','malt','params','lat_mom','-v7.3');
% save(fileobj,'relaxE','relaxE2','bestE','bestE2','acc_rate','C_v','C_fdt','malt','lat_mom','params','TT','temp','-v7.3'); % For paraTloop()
save(fileobj,'relaxE','relaxE2','bestE','bestE2','acc_rate','Marker','C_v','C_fdt','malt','lat_mom','params','-v7.3');

t_intv = datevec(toc(t0)/(60*60*24));
fprintf('Total time cost: %1$u days %2$u hours %3$u minutes and %4$.2f seconds\n', t_intv(3), t_intv(4), t_intv(5), t_intv(6));

clearvars
end