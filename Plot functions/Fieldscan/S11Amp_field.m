%First two lines for MacOS and last two for Windows
addpath('/Users/yikaiyang/Google Drive/File sharing/Programming scripts/Matlab/Plot functions/Fieldscan')
cd('/Users/yikaiyang/Google Drive/File sharing/PhD projects/LiHoF4 project/Data/Experiment/LiHoF4/SC127/S11_Amp_Field');
% addpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\Fieldscan')
% cd('G:\My Drive\File sharing\PhD projects\LiHoF4 project\Data\Experiment\LiHoF4\SC127\S11_Amp_Field')
order = 10;
figure
hold on
lgd = {''};

% load('Up_0.400_-30dBm.mat');
% dB0 = medfilt1(dB0,order);
% plot(H0,dB0,'>-');
% lgd = [lgd, 'Up 400mK -30dBm'];
% 
% load('Up_0.084_-30dBm.mat');
% dB0 = medfilt1(dB0,order);
% plot(H0,dB0,'o-');
% lgd = [lgd,'Up 84mK -30dBm'];

% load('Down_0.080_-20dBm.mat');
% dB0 = medfilt1(dB0,order);
% plot(H0,dB0,'s-');
% lgd = [lgd,'Down 80mK -20dBm']; 
% 
% load('Up_0.400_-20dBm.mat');
% dB0 = medfilt1(dB0,order);
% plot(H0,dB0,'>-');
% lgd = [lgd, 'Up 400mK -20dBm'];

% load('Up_0.080_-10dBm.mat');
% dB0 = medfilt1(dB0,order);
% plot(H0,dB0,'*-');
% lgd = [lgd,'Up 80mK -10dBm'];
% 
% load('Down_0.400_-10dBm.mat');
% dB0 = medfilt1(dB0,order);
% plot(H0,dB0,'<-');
% lgd = [lgd, 'Down 400mK -10dBm'];

% load('Up_0.150_0dBm.mat');
% dB0 = medfilt1(dB0,order);
% plot(H0,dB0,'d-');
% lgd = [lgd, 'Up 150mK 0dBm'];
% 
% load('Down_0.400_0dBm.mat');
% dB0 = medfilt1(dB0,order);
% plot(H0,dB0,'<-');
% lgd = [lgd, 'Down 400mK 0dBm'];
% 
% load('Up_0.200_4dBm.mat');
% dB0 = medfilt1(dB0,order);
% plot(H0,dB0,'.-');
% lgd = [lgd, 'Up 200mK +4dBm'];
% 
% load('Up_0.400_+4dBm.mat');
% dB0 = medfilt1(dB0,order);
% plot(H0,dB0,'<-');
% lgd = [lgd, 'Up 400mK +4dBm'];
% 
% load('Down_0.250_7dBm.mat');
% dB0 = medfilt1(dB0,order);
% plot(H0,dB0,'x-');
% lgd = [lgd, 'Down 250mK +7dBm'];
% 
% load('Down_0.400_+7dBm.mat');
% dB0 = medfilt1(dB0,2*order);
% plot(H0,dB0,'<-');
% lgd = [lgd, 'Down 400mK +7dBm'];

% load('Down_0.350_+10dBm.mat');
% dB0 = medfilt1(dB0,order);
% plot(H0,dB0,'^-');
% lgd = [lgd, 'Down 350mK +10dBm'];
% 
% load('Up_0.400_+10dBm.mat');
% dB0 = medfilt1(dB0,3*order);
% plot(H0,dB0,'<-');
% lgd = [lgd, 'Down 400mK +10dBm'];

% load('Down_0.400_+15dBm.mat');
% dB0 = medfilt1(dB0,order);
% plot(H0,dB0,'<-');
% lgd = [lgd, 'Down 400mK +15dBm'];

xlabel('Magnetic field (T)');
ylabel('S11 Amplitude');
legend(lgd(~cellfun('isempty',lgd)));
axis([1.5 3 -40 5]);
clearvars