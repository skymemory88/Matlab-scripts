addpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\spec1d--Henrik');
addpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\Fieldscan\functions');
format long g;
scan = importdata('G:\My Drive\File sharing\PhD projects\LiReF4\LiHoF4 project\Data\Experiment\LiHoF4\SC136\2020.02.19\detune_300K_air.dat',',',3);

freq11=scan.data(:,1);
S11=scan.data(:,2);
S11i=scan.data(:,3);
amp11=abs(S11+1i*S11i);
dB_S11=mag2db(amp11);
% psd = db2pow(-amp11)./freq11;    %Calculate power spectral density
% compl_11 = mag2db(1 - sqrt(S11.^2+S11i.^2));

%Option 1: Calculate graphic Quality factors with peak finding and -3dB
[~,Idx11] = min(dB_S11); %Find the peak and its corresponding index
bw11 = range(freq11(dB_S11 <= -3));
f0 = freq11(Idx11);
Q11 = freq11(Idx11)/bw11; %Calculate the quality factor

% Option 2: Calculate Quality factors from custom lorentzian fit
param = [bw11 f0 0 1];
% Param = [1.Bandwidth 2.Resonant frequency 3.Noise floor 4.Scaling factor]
low_bound = [min(freq11) min(freq11) -Inf -Inf];
upr_bound = [max(freq11) max(freq11) Inf Inf];
fdB_S11 = Lorentz_fit(freq11,dB_S11,param, [low_bound upr_bound]);
param = coeffvalues(fdB_S11);
fbw11 = param(1);
ff0 = param(2);
fQ11 = ff0/fbw11;
clearvars param fbw11 fQ11 Idx11 bw11 f0 Q11

% % Option 3: Calculate Quality factors from spec1d lorentzian fit
% zqf = medfilt1(dB_S11,10);
% fix = [1 1 1 1];
% s = spec1d(freq11,zqf,min(zqf)*0.01);
% p = [-1 f0 range(bw11) max(dB_S11)];
% %starting point for the (Lorentzian) fitting parameters(p1: scaling factor ,p2: resonant frequency; p3: FWHM; p4:noise floor)
% [fQ_1, fbck]=fits(s,'lorz', p, fix);
% Q11_fit = abs(fbck.pvals(2)/fbck.pvals(3)/2.0); %Calculate the quality factor
% chi11 = 1/Q11_fit;
% disp(num2str(Q11_fit,'Calculated quality factor for curve 1 from fitting: %3.2f.'));
% clear fbck s p zqf mag bw11 Idx11 pk;

% freq21=scan.data(:,1);
% S21=scan.data(:,4);
% S21i=scan.data(:,5);
% amp21=sqrt(S21.^2+S21i.^2);
% dB_S21=mag2db(amp21);
% %psd = db2pow(dB_S21)./freq2;    %Calculate power density spectrum
% %compl_21 = mag2db(1 - sqrt(S21.^2+S21i.^2));
% bw21 = [];
% j=1;
% for i=1:size(dB_S21)
%     if dB_S21(i) <= -3
%         bw21(j)=freq21(i); %group all frequencies below -3 dB.
%         j=j+1;
%         i=i+1;
%     else
%         i=i+1;
%     end
% end
% 
% clear j i;
% 
% bw21 = nonzeros(bw21);  %remove zero elements if necessary
% [Peak21,Idx21] = min(dB_S21); %Find the peak and its corresponding index
% Q21 = freq21(Idx21)/(max(bw21)-min(bw21)); %Calculate the quality factor
% 
% scan = importdata('transsmission_300K.dat',',',3);
%For import from a second file.

freq11_2=scan.data(:,1);
S11_2=scan.data(:,4);
S11i_2=scan.data(:,5);
amp11_2=sqrt(S11_2.^2+S11i_2.^2);
dB_S11_2=mag2db(amp11_2);
% psd_2 = db2pow(-amp11_2)./freq11_2;    %Calculate power density spectrum

% %Option 1: Calculate graphic Quality factors from peaking finding and -3 dB
[~,Idx11_2] = min(amp11_2); %Find the peak and its corresponding index
bw11_2 = range(freq11_2(dB_S11_2 <= -3));
f0_2 = freq11_2(Idx11_2);
Q11_2 = freq11_2(Idx11_2)/range(bw11_2); %Calculate the quality factor

% Option 2: Calculate Quality factors from custom lorentzian fit
param = [bw11_2 f0_2 0 1]; % Set up starting point of the parameters for the lorentzian fit
% Param = [1.Bandwidth 2.Resonant frequency 3.Noise floor 4.Scaling factor]
fdB_S11_2 = Lorentz_fit(freq11_2,dB_S11_2,param, [low_bound upr_bound]);
param = coeffvalues(fdB_S11_2);
fbw11_2 = param(1);
ff0_2 = param(2);
fQ11_2 = ff0_2/fbw11_2;
clearvars Idx11_2 bw11_2 Q11_2 param fbw11_2 fQ11_2

% %Option 3: Calculate Quality factors from lorentzian fitting
% zqf = medfilt1(dB_S11_2,10);
% s2 = spec1d(freq11_2, zqf, min(zqf)*0.01);
% p2 = [-1 f0 range(bw11_2) max(dB_S11_2)];
% %starting point for the (Lorentzian) fitting parameters(p1: scaling factor ,p2: resonant frequency; p3: FWHM; p4:noise floor?)
% [fQ_2, fbck]=fits(s2,'lorz', p2, fix);
% Q11_fit_2 = abs(fbck.pvals(2)/fbck.pvals(3)/2.0); %Calculate the quality factor
% chi11_2 = 1/Q11_fit;
% disp(num2str(Q11_fit_2,'Calculated quality factor for curve 2 from fitting: %3.2f.'));
% clear fbck s p zqf mag_2 fix mag_2 bw11_2 Idx11_2 pk2;

% freq11_3=scan.data3(:,1);
% S11_3=scan.data3(:,3);
% S11i_3=scan.data3(:,2);
% amp3=sqrt(S11_3.^2+S11i_3.^2);
% dB3=mag2db(amp3);

plot(freq11,dB_S11,'-o','Markersize',1.5);
hold on
plot(fdB_S11,'red');
%plot(freq21,dB_S21,'-x','Markersize',1.5);
%plot(freq11,compl_11,'s', 'Markersize', 1.5);
%plot(freq21,compl_21,'s', 'Markersize', 1.5);

plot(freq11_2,dB_S11_2,'-+','Markersize',1.5);
plot(fdB_S11_2,'red');
% %Mark down the resonant frequencies from numerical calculations with verticle lines
%{
plot(freq11_3,dB3,'s','Markersize',1.5);
vline=zeros(50,2);
vline(:,1)=5.4639E9;
vline(:,2)=0:-60/69:-60;
vline(:,3)=5.7454E9;
vline(:,4)=vline(:,2);
plot(vline(:,1),vline(:,2),'r',vline(:,3),vline(:,4),'r');
plot(vline(:,1),vline(:,2),'r');
%}
title('S11 response at 300K');
% legend(sprintf('Cavity (f_0 = 3.57 GHz), Q_g = %.2f and Q_f = %.2f', Q11, Q11_fit), sprintf('Active (f_0 = 3.40GHz), Q_g = %.2f and Q_f = %.2f', Q11_2, Q11_fit_2));
legend(num2str(ff0/1e9, 'Cavity (f_0 = %3.2f GHz)'),num2str(ff0_2/1e9, 'Active (f_0 = %3.2f GHz)'));
xlabel('Frequency (Hz)');ylabel('S11 (dB)');
hold off
clearvars