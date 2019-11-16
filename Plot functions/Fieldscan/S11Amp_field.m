%First two lines for Windows and last two for MacOS
addpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\Fieldscan')
cd('G:\My Drive\File sharing\PhD projects\LiHoF4 project\Data\Experiment\LiHoF4\SC127\S11_Amp_Field')
% addpath('/Users/yikaiyang/Google Drive/File sharing/Programming scripts/Matlab/Plot functions/Fieldscan')
% cd('/Users/yikaiyang/Google Drive/File sharing/PhD projects/LiHoF4 project/Data/Experiment/LiHoF4/SC127/S11_Amp_Field');
order = 4;
figure
hold on
lgd = {''};
Hmin = 1.5;
Hmax = 3;
dBmax = [];
dBmin = [];

% load('DPPH_Down_4.000_-20dBm_0dB.mat');
% dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
% H0 = H0(H0 >= Hmin & H0 <= Hmax);
% dB0 = medfilt1(dB0,order);
% plot(H0,dB0,'o-');
% lgd = [lgd,'Down 4K -20dBm 0dB'];
% dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), -20];
% dBmin = [dBmin; unique(H0(dB0 == min(dB0))), min(dB0), -20];

% load('Up_0.080_-30dBm_0dB.mat');
% dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
% H0 = H0(H0 >= Hmin & H0 <= Hmax);
% % dB0 = medfilt1(dB0,order);
% dB0 = smoothdata(dB0,'loess');
% plot(H0,dB0,'o-');
% lgd = [lgd,'Up 84mK -30dBm 0dB'];
% dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), -30];
% dBmin = [dBmin; unique(H0(dB0 == min(dB0))), min(dB0), -30];

load('Up_0.400_-30dBm_0dB.mat');
dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
H0 = H0(H0 >= Hmin & H0 <= Hmax);
% dB0 = medfilt1(dB0,order);
dB0 = smoothdata(dB0,'loess');
plot(H0(1:20:end),dB0(1:20:end),'>-');
lgd = [lgd, 'Up 400mK -30dBm 0dB'];
dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), -30];
dBmin = [dBmin; unique(H0(dB0 == min(dB0))), min(dB0), -30];

% load('Down_0.080_-20dBm.mat');
% dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
% H0 = H0(H0 >= Hmin & H0 <= Hmax);
% % dB0 = medfilt1(dB0,order);
% dB0 = smoothdata(dB0,'loess');
% plot(H0,dB0,'s-');
% lgd = [lgd,'Down 80mK -20dBm 0dB']; 
% dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), -20];
% dBmin = [dBmin; unique(H0(dB0 == min(dB0))), min(dB0), -20];

% load('Down_0.250_-20dBm_30dB.mat');
% dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
% H0 = H0(H0 >= Hmin & H0 <= Hmax);
% % dB0 = medfilt1(dB0,order);
% dB0 = smoothdata(dB0,'loess');
% plot(H0,dB0,'>-');
% lgd = [lgd, 'Down 250mK -20dBm 30dB'];
% dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), -20];
% dBmin = [dBmin; unique(H0(dB0 == min(dB0))), min(dB0), -20];

% load('Up_0.400_-20dBm_0dB.mat');
% dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
% H0 = H0(H0 >= Hmin & H0 <= Hmax);
% % dB0 = medfilt1(dB0,order);
% dB0 = smoothdata(dB0,'loess');
% plot(H0,dB0,'>-');
% lgd = [lgd, 'Up 400mK -20dBm 0dB'];
% dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), -20];
% dBmin = [dBmin; unique(H0(dB0 == min(dB0))), min(dB0), -20];

load('Up_0.400_-20dBm_30dB.mat');
dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
H0 = H0(H0 >= Hmin & H0 <= Hmax);
% dB0 = medfilt1(dB0,order);
dB0 = smoothdata(dB0,'loess');
plot(H0(1:20:end),dB0(1:20:end),'o-');
% plot(H0(1:20:end),smoothdata(dB0(1:20:end)),'o-');
lgd = [lgd, 'Up 400mK -20dBm 30dB'];
dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), -20];
dBmin = [dBmin; unique(H0(dB0 == min(dB0))), min(dB0), -20];

% load('Up_0.080_-10dBm.mat');
% dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
% H0 = H0(H0 > Hmin & H0 <= Hmax);
% % dB0 = medfilt1(dB0,order);
% dB0 = smoothdata(dB0,'loess');
% plot(H0,dB0,'*-');
% lgd = [lgd,'Up 80mK -10dBm 0dB'];
% dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), -10];
% dBmin = [dBmin; unique(H0(dB0 == min(dB0))), min(dB0), -10];

% load('Up_0.250_-10dBm_30dB.mat');
% dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
% H0 = H0(H0 >= Hmin & H0 <= Hmax);
% % dB0 = medfilt1(dB0,order);
% dB0 = smoothdata(dB0,'loess');
% plot(H0(1:20:end),smoothdata(dB0(1:20:end)),'s-');
% lgd = [lgd, 'Down 250mK -10dBm 30dB'];
% dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), -10];
% dBmin = [dBmin; unique(H0(dB0 == min(dB0))), min(dB0), -10];

% load('Down_0.400_-10dBm_0dB.mat');
% dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
% H0 = H0(H0 >= Hmin & H0 <= Hmax);
% % dB0 = medfilt1(dB0,order);
% dB0 = smoothdata(dB0,'loess');
% plot(H0,dB0,'<-');
% lgd = [lgd, 'Down 400mK -10dBm 0dB'];
% dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), -10];
% dBmin = [dBmin; unique(H0(dB0 == min(dB0))), min(dB0), -10];

load('Down_0.400_-10dBm_30dB.mat');
dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
H0 = H0(H0 >= Hmin & H0 <= Hmax);
% dB0 = medfilt1(dB0,order);
dB0 = smoothdata(dB0,'loess');
plot(H0(1:20:end),dB0(1:20:end),'s-');
% plot(H0(1:20:end),smoothdata(dB0(1:20:end)),'s-');
lgd = [lgd, 'Down 400mK -10dBm 30dB'];
dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), -10];
dBmin = [dBmin; unique(H0(dB0 == min(dB0))), min(dB0), -10];

% load('Up_0.150_0dBm.mat');
% dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
% H0 = H0(H0 >= Hmin & H0 <= Hmax);
% % dB0 = medfilt1(dB0,order);
% dB0 = smoothdata(dB0,'loess');
% plot(H0,dB0,'d-');
% lgd = [lgd, 'Up 150mK 0dBm 30dB'];
% dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), 0];
% dBmin = [dBmin; unique(H0(dB0 == min(dB0))), min(dB0), 0];

% load('Down_0.290_0dBm.mat');
% dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
% H0 = H0(H0 >= Hmin & H0 <= Hmax);
% % dB0 = medfilt1(dB0,order);
% dB0 = smoothdata(dB0,'loess');
% plot(H0,dB0,'d-');
% lgd = [lgd, 'Down 290mK 0dBm 30dB'];
% dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), 0];
% dBmin = [dBmin; unique(H0(dB0 == min(dB0))), min(dB0), 0];
 
load('Down_0.400_0dBm.mat');
dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
H0 = H0(H0 >= Hmin & H0 <= Hmax);
% dB0 = medfilt1(dB0,order);
dB0 = smoothdata(dB0,'loess');
plot(H0(1:20:end),dB0(1:20:end),'d-');
% plot(H0(1:20:end),smoothdata(dB0(1:20:end)),'d-');
lgd = [lgd, 'Down 400mK 0dBm 30dB'];
dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), 0];
dBmin = [dBmin; unique(H0(dB0 == min(dB0))), min(dB0), 0];

% load('Up_0.200_+4dBm.mat');
% dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
% H0 = H0(H0 >= Hmin & H0 <= Hmax);
% % dB0 = medfilt1(dB0,order);
% dB0 = smoothdata(dB0,'loess');
% plot(H0,dB0,'.-');
% lgd = [lgd, 'Up 200mK +4dBm 20dB'];
% dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), 4];
% dBmin = [dBmin; unique(H0(dB0 == min(dB0))), min(dB0), 4];

% load('Up_0.320_+4dBm.mat');
% dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
% H0 = H0(H0 >= Hmin & H0 <= Hmax);
% % dB0 = medfilt1(dB0,order);
% dB0 = smoothdata(dB0,'loess');
% plot(H0,dB0,'s-');
% lgd = [lgd, 'Up 300mK +4dBm 30dB'];
% dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), 4];
% dBmin = [dBmin; unique(H0(dB0 == min(dB0))), min(dB0), 4];

load('Up_0.400_+4dBm.mat');
dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
H0 = H0(H0 >= Hmin & H0 <= Hmax);
% dB0 = medfilt1(dB0,order+5);
dB0 = smoothdata(dB0,'loess');
plot(H0(1:20:end),dB0(1:20:end),'.-');
% plot(H0(1:20:end),smoothdata(dB0(1:20:end)),'.-');
lgd = [lgd, 'Up 400mK +4dBm 30dB'];
dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), 4];
dBmin = [dBmin; unique(H0(dB0 == min(dB0))), min(dB0), 4];

% load('Down_0.250_7dBm.mat');
% dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
% H0 = H0(H0 >= Hmin & H0 <= Hmax);
% % dB0 = medfilt1(dB0,order);
% dB0 = smoothdata(dB0,'loess');
% plot(H0,dB0,'x-');
% lgd = [lgd, 'Down 250mK +7dBm 20dB'];
% dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), 7];
% dBmin = [dBmin; unique(H0(dB0 == min(dB0))), min(dB0), 7];

% load('Down_0.350_+7dBm.mat');
% dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
% H0 = H0(H0 >= Hmin & H0 <= Hmax);
% % dB0 = medfilt1(dB0,order);
% dB0 = smoothdata(dB0,'loess');
% plot(H0,dB0,'x-');
% lgd = [lgd, 'Down 350mK +7dBm 30dB'];
% dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), 7];
% dBmin = [dBmin; unique(H0(dB0 == min(dB0))), min(dB0), 7];

load('Down_0.400_+7dBm.mat');
dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
H0 = H0(H0 >= Hmin & H0 <= Hmax);
% dB0 = medfilt1(dB0,2*order);
dB0 = smoothdata(dB0,'loess');
plot(H0(1:20:end),dB0(1:20:end),'>-');
% plot(H0(1:20:end),smoothdata(dB0(1:20:end)),'>-');
lgd = [lgd, 'Down 400mK +7dBm 30dB'];
dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), 7];
dBmin = [dBmin; unique(H0(dB0 == min(dB0))), min(dB0), 7];

% load('Down_0.350_+10dBm.mat');
% dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
% H0 = H0(H0 >= Hmin & H0 <= Hmax);
% % dB0 = medfilt1(dB0,3*order);
% dB0 = smoothdata(dB0,'loess');
% plot(H0,dB0,'^-');
% lgd = [lgd, 'Down 350mK +10dBm 30dB'];
% dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), 10];
% dBmin = [dBmin; unique(H0(dB0 == min(dB0))), min(dB0), 10];

% load('Up_0.400_+10dBm_1.mat');
% dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
% H0 = H0(H0 >= Hmin & H0 <= Hmax);
% % dB0 = medfilt1(dB0,3*order);
% dB0 = smoothdata(dB0,'loess');
% plot(H0(1:20:end),dB0(1:20:end),'x-');
% % plot(H0(1:20:end),smoothdata(dB0(1:20:end)),'x-');
% lgd = [lgd, 'Down 400mK +10dBm 30dB'];
% dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), 10];
% dBmin = [dBmin; unique(H0(dB0 == min(dB0))), min(dB0), 10];

% load('Up_0.400_+10dBm_2.mat');
% dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
% H0 = H0(H0 >= Hmin & H0 <= Hmax);
% dB0 = medfilt1(dB0,3*order);
% % dB0 = smoothdata(dB0,'loess');
% plot(H0(1:20:end),dB0(1:20:end),'x-');
% lgd = [lgd, 'Down 400mK +10dBm 30dB'];
% dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), 10];
% dBmin = [dBmin; unique(H0(dB0 == min(dB0))), min(dB0), 10];

load('Down_0.400_+15dBm.mat');
dB0 = dB0(H0 >= Hmin & H0 <= Hmax);
H0 = H0(H0 >= Hmin & H0 <= Hmax);
dB0 = medfilt1(dB0, order);
plot(H0,dB0,'<-');
lgd = [lgd, 'Down 400mK +15dBm 30dB'];
dBmax = [dBmax; unique(H0(dB0 == max(dB0))), max(dB0), 15];
dBmin = [dBmin; mean(unique(H0(dB0 == min(dB0)))), min(dB0), 15];

xlabel('Magnetic field (T)');
ylabel('S11 Amplitude');
legend(lgd(~cellfun('isempty',lgd)));
axis([Hmin Hmax -40 5]);

% figure
% plot(dBmax(:,1),dBmax(:,2),'o-');
% xlabel('Magnetic field (T)');
% axis([1.5 3 -40 5]);
% ylabel('Maximum S11');

figure
plot(dBmax(:,3),dBmax(:,2),'s-');
yfit = polyfit(dBmax(:,3),dBmax(:,2),2);
hold on
plot(dBmax(:,3),polyval(yfit, dBmax(:,3)));
xlabel('Power (dBm)');
ylabel('Maximum S11');

clearvars