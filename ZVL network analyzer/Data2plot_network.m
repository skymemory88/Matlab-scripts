cd('C:\Users\yiyang\Google Drive\File sharing\PhD projects\Transverse Ising model\Data\Experiment\Cavity resonator\D24mm_T5mm_G0.2mm\24.04.2019');
scan = importdata('D24mm_t5mm_g0.2mm_detuning_300K.dat',',',3);

freq11=scan.data(:,1);
S11=scan.data(:,2);
S11i=scan.data(:,3);
amp11=sqrt(S11.^2+S11i.^2);
ampdB_S11=mag2db(amp11);
% psd = db2pow(ampdB_S11)./freq;    %Calculate power density spectrum
% compl_11 = mag2db(1 - sqrt(S11.^2+S11i.^2));
bw11 = [];
j=1;
for i=1:size(ampdB_S11)
    if ampdB_S11(i) <= -3.0
        bw11(j)=freq11(i); %group all frequencies below -3 dB.
        j=j+1;
    end
end
% bw11 = nonzeros(bw11);  %remove zero elements if necessary
[~,Idx11] = min(ampdB_S11); %Find the peak and its corresponding index
Q11 = freq11(Idx11)/(max(bw11)-min(bw11)); %Calculate the quality factor

clear j i;

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
%psd = db2pow(ampdB2)./freq2;    %Calculate power density spectrum
% 
bw11_2 = [];
j=1;
for i=1:size(ampdB_S11_2)
    if ampdB_S11_2(i) <= -3
        bw11_2(j)=freq11_2(i); %group all frequencies below -3 dB.
        j=j+1;
        i=i+1;
    else
        i=i+1;
    end
end
% bw11_2 = nonzeros(bw11_2);  %remove zero elements if necessary
[Peak11_2,Idx11_2] = min(ampdB_S11_2); %Find the peak and its corresponding index
Q11_2 = freq11_2(Idx11_2)/(max(bw11_2)-min(bw11_2)); %Calculate the quality factor

clear j i;
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
legend(sprintf('S11 undercoupled 300K with Q = %.2f', Q11), sprintf('S11 critical 300K with Q = %.2f', Q11_2));
xlabel('Frequency (Hz)');ylabel('S11 (dB)');