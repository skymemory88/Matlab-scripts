function MC_process(file_code,measN,meas_intv)
file_name=sprintf(['results_tempscan_',file_code,'_%u_%u.mat'],measN,meas_intv);
file_path = 'G:\My Drive\File sharing\PhD program\Research projects\LiErF4 project\Quantum Monte Carlo\Test';
file_obj = fullfile(file_path, file_name);
load(file_obj,'acc_rate','bestE','bestE2','C_fdt','C_v','lat_mom','malt','meas','params','relaxE','relaxE2','TT');

[temp, idx] = sort(TT,'ascend');
C_v = C_v(idx);
% C_fdt = C_fdt(idx(2:end)-1);
malt = malt(idx);
bestE = bestE(idx);
% bestE2 = bestE2(idx);
acc_rate = cell2mat(acc_rate');

figure
plot(temp, C_v, '-o');
ylabel('Specific heat')
hold on
yyaxis right
plot(temp(2:end),C_fdt(idx(2:end)-1),'-s');
xlabel('Temperature (K)')
ylabel('Specific heat')
legend('C_v','C_{fdt}')
title('Specific heat from F-D theorem and dF/dT')

figure
plot(temp, malt, '-o');
xlabel('Temperature (K)')
ylabel('Order parameter')
title('Alternating moments as order parameter')

figure
plot(acc_rate(C_v == max(C_v),:));
xlabel('Iteration (LSW)')
ylabel('Acceptance rate')

figure
plot(bestE,'-o')
xlabel('Iteration (LSW)')
ylabel('Global energy')

return