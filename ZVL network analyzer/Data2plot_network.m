cd('C:\Users\yiyang\Google Drive\File sharing\PhD projects\Transverse Ising model\Data\Experiment\LiHoF4\18.04.2019');
format long g;
scan = importdata('detuning_300K_2nd_try.dat',',',3);

freq11=scan.data(:,1);
S11=scan.data(:,2);
S11i=scan.data(:,3);
amp11=abs(S11+1i*S11i);
ampdB_S11=mag2db(amp11);
% psd = db2pow(-amp11)./freq11;    %Calculate power spectral density
% compl_11 = mag2db(1 - sqrt(S11.^2+S11i.^2));

%Option 1: Calculate graphic Quality factors with peak finding and bw and -3dB
[pk,Idx11] = min(amp11); %Find the peak and its corresponding index
bw11 = freq11(amp11<=(1-pk/2));
f0 = freq11(Idx11);
Q11 = freq11(Idx11)/range(bw11); %Calculate the quality factor

clear j i;

mag = (1- S11.^2 -S11i.^2)./((1- S11).^2 + S11i.^2); %1-abs(S11)^2/( (1-S11_real)^2+S11_img^2 )
zqf = medfilt1(mag,10);
s = spec1d(freq11,zqf,zqf.*0 + 0.05);
p = [-1 f0 range(bw11) 2];
%starting point for the (Lorentzian) fitting parameters(p1: scaling factor ,p2: resonant frequency; p3: FWHM; p4:noise floor?)
fix = [1 1 1 1];
[fQ, fbck]=fits(s,'lorz',p,fix);
Q11_fit = fbck.pvals(2)/fbck.pvals(3); %Calculate the quality factor
chi11 = 1/Q11_fit;
disp(num2str(Q11_fit,'First calculated quality factor from fitting: %3.2f.'));
clear fQ fbck s p zqf fix mag bw11 f0 Idx11 pk;

% freq21=scan.data(:,1);
% S21=scan.data(:,4);
% S21i=scan.data(:,5);
% amp21=sqrt(S21.^2+S21i.^2);
% ampdB_S21=mag2db(amp21);
% %psd = db2pow(ampdB_S21)./freq2;    %Calculate power density spectrum
% %compl_21 = mag2db(1 - sqrt(S21.^2+S21i.^2));
% bw21 = [];
% j=1;
% for i=1:size(ampdB_S21)
%     if ampdB_S21(i) <= -3
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
% [Peak21,Idx21] = min(ampdB_S21); %Find the peak and its corresponding index
% Q21 = freq21(Idx21)/(max(bw21)-min(bw21)); %Calculate the quality factor
% 
freq11_2=scan.data(:,1);
S11_2=scan.data(:,4);
S11i_2=scan.data(:,5);
amp11_2=sqrt(S11_2.^2+S11i_2.^2);
ampdB_S11_2=mag2db(amp11_2);
% psd_2 = db2pow(-amp11_2)./freq11_2;    %Calculate power density spectrum
% 
[pk2,Idx11_2] = min(amp11_2); %Find the peak and its corresponding index
bw11_2= freq11_2(amp11_2<=(1-pk2/2));
f0_2 = freq11_2(Idx11_2);
Q11_2 = freq11_2(Idx11_2)/range(bw11_2); %Calculate the quality factor

clear j i; 

mag_2 = (1- S11_2.^2 -S11i_2.^2)./((1- S11_2).^2 + S11i_2.^2); %1-abs(S11)^2/( (1-S11_real)^2+S11_img^2 )
zqf = medfilt1(mag_2,10);
s = spec1d(freq11_2,zqf,zqf.*0 + 0.05);
p = [-1 f0_2 range(bw11_2) 2];
%starting point for the (Lorentzian) fitting parameters(p1: scaling factor ,p2: resonant frequency; p3: FWHM; p4:noise floor?)
fix = [1 1 1 1];
[fQ, fbck]=fits(s,'lorz',p,fix);
Q11_fit_2 = fbck.pvals(2)/fbck.pvals(3); %Calculate the quality factor
chi11_2 = 1/Q11_fit;
disp(num2str(Q11_fit_2,'Second calculated quality factor from fitting: %3.2f.'));
clear fQ fbck s p zqf mag_2 fix mag_2 bw11_2 f0_2 Idx11_2 pk2;

% freq11_3=scan.data3(:,1);
% S11_3=scan.data3(:,3);
% S11i_3=scan.data3(:,2);
% amp3=sqrt(S11_3.^2+S11i_3.^2);
% ampdB3=mag2db(amp3);

plot(freq11,ampdB_S11,'-o','Markersize',1.5);
hold on
%plot(freq21,ampdB_S21,'-x','Markersize',1.5);
%plot(freq11,compl_11,'s', 'Markersize', 1.5);
%plot(freq21,compl_21,'s', 'Markersize', 1.5);

plot(freq11_2,ampdB_S11_2,'-+','Markersize',1.5);
%{
plot(freq11_3,ampdB3,'s','Markersize',1.5);
vline=zeros(50,2);
vline(:,1)=5.4639E9;
vline(:,2)=0:-60/69:-60;
vline(:,3)=5.7454E9;
vline(:,4)=vline(:,2);
plot(vline(:,1),vline(:,2),'r',vline(:,3),vline(:,4),'r');
plot(vline(:,1),vline(:,2),'r');
%}
title('S11 response at 300K');
legend(sprintf('Detunned, Q_g = %.2f and Q_f = %.2f', Q11, Q11_fit), sprintf('Critical, Q_g = %.2f and Q_f = %.2f', Q11_2, Q11_fit_2));
xlabel('Frequency (Hz)');ylabel('S11 (dB)');
hold off