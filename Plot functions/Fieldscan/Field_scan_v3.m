function Field_scan
format long;  
% matlabpool open

%addpath('C:\Users\babkevic\Documents\MATLAB\legendflex')

% % curdir = cd;
addpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\Fieldscan\functions');
addpath(genpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\spec1d--Henrik'));
% addpath('/Volumes/GoogleDrive/My Drive/File sharing/Programming scripts/Matlab/Plot functions/Fieldscan/functions')
% addpath(genpath('/Volumes/GoogleDrive/My Drive/File sharing/Programming scripts/Matlab/Plot functions/spec1d--Henrik'));
%The first line is for windows, the second line is for mac OS
%% Define input ===========================================================
% Figure plot options:
plotopt.lnwd = 2;
plotopt.ftsz = 12;
plotopt.mksz = 5;

filepath = 'G:\My Drive\File sharing\PhD projects\LiHoF4 project\Data\Experiment\LiHoF4\SC147\14.11.2019/';
% filepath = '/Users/yikaiyang/Google Drive/File sharing/PhD projects/LiHoF4 project/Data/Experiment/LiHoF4/SC127/05.11.2019/';
%The first line is for windows, the second line is for mac OS
filename = '2019_11_0033';
%% Read ZVL
%Set data range and parameters
opt  = 2;

switch opt
    case 1
        % Test of frequency scans using ZVL as a function of field
        option1(filepath, filename, plotopt)
    case 2
        % Analyse the fieldscan data and display plots
        option2(filepath, filename, plotopt)
    case 3
        % Use routine to analyze the fieldscan data and save to files
        option3(filepath, filename)
end
% matlabpool close

end

%% ------------------------------------------------------------------------

function option1(filepath,filename, plotopt)

run = 1;
out = readdata_v3(filepath,filename,run);

%% Read & Plot data
freq = out.data.ZVLfreq/1e9;
S11 = out.data.ZVLreal + 1i*out.data.ZVLimag;
H = out.data.DCField1;
HH = repmat(H,1,size(freq,2)); %populate the magnetic field to match the dimension of S11 matrix
T1 = out.data.Temperature1;
% T2 = out.data.Temperature2;

S11 = S11';
S11 = S11(:);
dB = mag2db(abs(S11));
freq = freq';
freq = freq(:);
HH = HH';
HH = HH(:); %the third argument is the number of frequency points in each line/segment

%Set data range and parameters
freq_l = min(freq); %set frequency range, l: lower limit, h: higher limit
freq_h = max(freq);
field_l = min(H);  %set field range, l: lower limit, h: higher limit
field_h = max(H);

% Could use "scatteredInterpolant()" to replace "TriScatteredInterp()" as recommended by MATLAB, but it may generate artifacts
FdB = TriScatteredInterp(HH,freq,dB);
% FrS = TriScatteredInterp(HH,freq,real(S11));
% FiS = TriScatteredInterp(HH,freq,imag(S11)); %intrapolate points on a 2D grid

% Plot frequency-field colour map
[xq,yq] = meshgrid(linspace(field_l,field_h,301),linspace(freq_l,freq_h,310)); %set the X and Y range
zq = FdB(xq,yq);

% Plot the interpolated frequency response data in a field scan using color map
% hfig1 = setfig(11);
figure
cmap = pcolor(xq,yq,zq);
set(cmap, 'edgeColor','none')
shading interp;
colorbar
set(gca,'fontsize',plotopt.ftsz)
t(1) = xlabel('Field (T)');
t(2) = ylabel('Frequency (GHz)');
tt(1) = title([filename num2str(run)]);

% Plot the temperature profile against magnetic field
% hfig2 = setfig(12);
figure
plot(H(1:100:end),T1(1:100:end),'o-')
xlabel('DC Magnetic field')
ylabel('Temperature')
title('Magnetic field vs Temperature')
%
% yq = yq.*0 + 3.200; % Plot field scan cut at 3.200 GHz (H vs dB)
% hfig3 = setfig(12);
% zq = FdB(xq,yq);
% h = plot(xq(1,:),zq(1,:));
% xlim([0 9])
% set(gca,'fontsize',plotopt.ftsz)
% t(3) = xlabel('Field (T)');
% t(4) = ylabel('S11 (dB)');
% tt(2) = title([filename num2str(run) ', f = 3.200 GHz']);

% hfig4 = setfig(13); % Plot field scan cut at 3.200 GHz (H vs real(S))
% zq = FrS(xq,yq);
% h = plot(xq(1,:),zq(1,:));
% xlim([0 9])
% set(gca,'fontsize',plotopt.ftsz)
% t(5) = xlabel('Field (T)');
% t(6) = ylabel('Re{S11}');
% tt(3) = title([filename num2str(run) ', f = 3.200 GHz']);
%
% hfig5 = setfig(14); % Plot field scan cut at 3.200 GHz (H vs real(S))
% zq = FiS(xq,yq);
% h = plot(xq(1,:),zq(1,:));
% xlim([0 9])
% set(gca,'fontsize',plotopt.ftsz)
% t(7) = xlabel('Field (T)');
% t(8) = ylabel('Im{S11}');
% tt(4) = title([filename num2str(run) ', f = 3.200 GHz']);
%
%
% % Plot frequency scans at H (f vs dB)
% [xq,yq] = meshgrid(linspace(0,9,501),linspace(3.2,3.7,101));
% xq1 = xq.*0 + 0;
% zq1 = FdB(xq1,yq);
% xq2 = xq.*0 + 2;
% zq2 = FdB(xq2,yq);
% xq3 = xq.*0 + 3.612;
% zq3 = FdB(xq3,yq);
% xq4 = xq.*0 + 3.685;
% zq4 = FdB(xq4,yq);
% xq5 = xq.*0 + 5;
% zq5 = FdB(xq5,yq);
%
% hfig6 = setfig(15);
% h = plot(yq(:,1),zq1(:,1),yq(:,1),zq2(:,1),yq(:,1),zq3(:,1),yq(:,1),zq4(:,1),yq(:,1),zq5(:,1));
% xlim([3.2 3.7])
% set(gca,'fontsize',plotopt.ftsz)
% t(9) = xlabel('Frequency (GHz)');
% t(10) = ylabel('S11 (dB)');
% tt(5) = title([filename num2str(run)]);
% t(11) = legend(h,'0','2','3.612','3.685','5');

% hfig7 = setfig(16); % Plot field-temperature
% h = plot(H,T1,H,T2);
% xlim([0 9])
% set(gca,'fontsize',plotopt.ftsz)
% t(12) = xlabel('Field (T)');
% t(13) = ylabel('Temperature (K)');
% tt(6) = title([filename num2str(run)]);
% t(14) = legend(h,'MC','Sample');


set(t,'fontsize',plotopt.ftsz,'interpreter','latex')
set(tt,'fontsize',plotopt.ftsz,'interpreter','none')

% Save plots
%figname = [filename num2str(run)];
%saveplots(hfig1,[figname '_field-freq_dB_map'])
% saveplots(hfig2,[figname '_field-dB_3.200GHz'])
% saveplots(hfig3,[figname '_field-rS_3.200GHz'])
% saveplots(hfig4,[figname '_field-iS_3.200GHz'])
% saveplots(hfig5,[figname '_freq_dB_H'])
% saveplots(hfig6,[figname '_field_T'])

end
function option2(filepath,filename, plotopt)
Temperature = 0.080;
%Set data range and parameters
order = 4; % set to what order the median filters is applied
FWHM = -3.0; % Define full-width-half-max for calculation of quality factors
base = 2; % Manually define the noise floor

clear freq S11 dB N FdB FrS FiS FTT1 FTT2

% extract data from raw data file
out = readdata_v3(filepath,filename,1);
freq = out.data.ZVLfreq/1e9;
S11 = out.data.ZVLreal + 1i*out.data.ZVLimag;
mag = (1- out.data.ZVLreal.^2 -out.data.ZVLimag.^2)./((1- out.data.ZVLreal).^2 + out.data.ZVLimag.^2);    
H = out.data.DCField1;
T1 = out.data.Temperature1;
HH = repmat(H,1,size(freq,2)); %populate the magnetic field to match the dimension of S11 matrix

freq = freq';
freq = freq(:);
S11 = S11';
S11 = S11(:);
HH = HH';
HH = HH(:);
mag = mag';
mag = mag(:);

freq_l = min(freq); %set frequency range, l: lower limit, h: higher limit
freq_h = max(freq);
field_l = min(H);  %set field range, l: lower limit, h: higher limit
field_h = max(H);

dif = nonzeros(diff(freq));
dif = dif(dif>0);
step = mean(dif);
% step = freq(r,c+1)-freq(r,c); %Scanning step in unit of GHz, default is 1 MHz according to ZVL network analyser
nop = ceil((freq_h-freq_l)/step)+1; %compute how many points pers complete frequency scan.
% lineDiv = (freq_h-freq_l)/step/(out.nop-1); %number of lines per complete frequency scan in the raw data file

clearvars dif

%Plot the temperature vs magnetic field to check the temperature variation
figure
plot(H(1:100:end),T1(1:100:end),'o-')
xlabel('DC Magnetic field')
ylabel('Temperature')
title('Magnetic field vs Temperature')

%Clean up the raw data by removing incomplete scans (step 1) and duplicates
%(step 2)
% step 1: truncate the beginning and end part to keep only complete frequency scans and reshape the matrices into single column vectors
trunc1 = find(freq==freq_l,1,'first'); 
trunc2 = find(freq==freq_h,1,'last'); 
S11_temp = S11(trunc1:trunc2);
dB_temp = mag2db(abs(S11_temp));
freq_temp = freq(trunc1:trunc2);
HH_temp = HH(trunc1:trunc2);
mag_temp = mag(trunc1:trunc2);

% Step 2: remove duplicates
dupl = find(diff(freq_temp) == 0);
freq_temp(dupl+1)=[];
dB_temp(dupl+1)=[];
HH_temp(dupl+1)=[];
mag_temp(dupl+1)=[];

dB = reshape(dB_temp,nop,[]);  %reshape the matrix so that each complete frequency scan occupy one column
freq = reshape(freq_temp,nop,[]);
HH = reshape(HH_temp,nop,[]);

%Find all the resonant peaks
f0 = zeros(size(dB,2),1);
H0 = zeros(size(dB,2),1);
Q0 = zeros(size(dB,2),1);
dB0 = zeros(size(dB,2),1);
%find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
for ii = 1:size(dB,2) %Searching column minima (fixed field) is better than searching row minima (fixed frequency)
[~,idx] = min( dB(:,ii) );
if(length(idx)>1)
    disp(num2str(H0(ii),'multiple minima found at H = %.3f'))
end
f0(ii) = freq(idx,ii);
H0(ii) = HH(idx,ii); 
dB0(ii) = dB(idx,ii);
% Calculate quality factor using f0/FWHM
if isnan(1/range(freq(dB(:,ii) <= FWHM)))
   Q0(ii) = 0;
elseif isempty(range(freq(dB(:,ii) <= FWHM)))
   Q0(ii) = 0;
else
Q0(ii) = freq(idx,ii)/range(freq(dB(:,ii) <= FWHM));
end
Q0(isnan(Q0)) = 0; % Cut out NaN from the array
Q0(isinf(Q0)) = 0; % Cut out inf from the array
Q0(Q0>4000) = 0; % Cut out unrealistic quality factor from noise data
end
%   For noisy data, we need to remove duplicates of minima
[H0,ia,~] = unique(H0,'stable');
f0 = f0(ia);
Q0 = Q0(ia);
dB0 = dB0(ia);

clearvars idx ia ii

% Plot the resonant frequency versus DC magnetic field
figure
f0 = medfilt1(f0,order); % apply median filter to remove some noise
hfig1 = plot(H0,f0,'o','MarkerSize',2);
xlabel('Field (T)');
ylabel('Resonant frequency (GHz)');
title(num2str(Temperature,'Minimal S11 at T = %3.3f K'));
axis([field_l field_h freq_l freq_h]);
% clearvars *_temp;

% Plot the peak amplitude vs. magnetic field
figure
dB0 = medfilt1(dB0, order); % apply median filter to remove some noise
nfig2 = plot(H0, dB0, 'o', 'MarkerSize', 2);
xlabel('Field(T)');
ylabel('S11 amplitute');
title(num2str(Temperature,'Minimal S11 at T = %3.3f K'));

[xq,yq] = meshgrid(linspace(field_l,field_h,301),linspace(freq_l,freq_h,310));

switch 1 %choose data interpolation method and plot the color map of S11 response
case 1 % Option_1 Interpolate data along only the frequency axis.
    interp_dB = zeros(size(yq,1),size(HH,2));
    for ii = 1:size(HH,2)
        interp_dB(:,ii) = interp1(freq(:,ii),dB(:,ii),yq(:,1));
    end
    zq = interp_dB'; % For later use in lorentzian fitting
    HH = repmat(H0,1,size(xq,1));
    
    switch 1 % Choose plotting option
        case 1
            %     plot single-direction interpolated data using pseudo-colormap
            cmap = pcolor(HH(:,1),yq(:,1),interp_dB);
            set(cmap, 'edgeColor','none')
            shading interp;
            colorbar
        case 2
            %     plot single-direction interpolated data using scatter plot
            yy = repmat(yq(:,1),1,length(H0));
            yy = yy(:);
            interp_dB = interp_dB(:);
            HH = HH';
            HH = HH(:);
            smap = scatter3(HH,yy,interp_dB,2,interp_dB,'o','filled','MarkerEdgeColor','none');
            colormap(hsv);
    end
case 2 % Option_2 Interpolate the data along both axis.
%             S11 = S11(:);
    mag = mag_temp;
    freq = freq_temp;
    dB = dB_temp;
    HH = HH_temp; 
%             FdB  = TriScatteredInterp(HH,freq,dB);
    FdB = scatteredInterpolant(HH, freq, dB);
    Fmag  = scatteredInterpolant(HH,freq, mag);
    %     FrS  = TriScatteredInterp(HH,freq,real(S11));
    %     FiS  = TriScatteredInterp(HH,freq,imag(S11));
    %     FTT1 = TriScatteredInterp(HH,freq,TT1);
    %     FTT2 = TriScatteredInterp(HH,freq,TT2);
    figure
    hold on
    box on
    zq = FdB(xq,yq);
    cmap = pcolor(xq,yq,zq);
    set(cmap, 'edgeColor','none')
    shading interp;
    colorbar
    zq = Fmag(xq,yq); % For later use in lorentzian fitting
end

set(gca,'fontsize',plotopt.ftsz)
axis([field_l field_h freq_l freq_h]);
tx(4) = xlabel('Field (T)');
ty(4) = ylabel('Frequency (GHz)');
tt(4) = title(num2str(Temperature,'S11 response at T = %3.3f K'));

%Fit each frequency scan to Lorentzian form to extract the quality factor
%     figure
disp(num2str(Temperature,'Fitting: T = %3.3f K'))

%   Extract frequency cuts
Hx = field_l:range(H)/49:field_h;
switch 2
case 1
%   Use the raw data for fitting    
    Hx = H0';
    ff0 = f0'; 
case 2
%   Use an interpolation of the resonant frequency along the DC magnetic field (Hx) axis
    ff0 = interp1(H0,f0,Hx,'spline','extrap'); 
end

for ii = 1:length(Hx)
xq = xq.*0 + Hx(ii);
zqf = zq(:,1); % zq was defined in the 'switch' section regarding interpolation direction 
%         zqf = medfilt1(zq(:,1), order); %Apply median filter to remove some noises
s = spec1d(yq(:,1), zqf, zqf.*0 + 0.05);
p = [-1 ff0(ii) 0.001 base];
%starting point for the (Lorentzian) fitting parameters(p1: scaling factor ,p2: resonant frequency; p3: FWHM; p4:noise floor?)

fix = [1 1 1 1];
[~, fbck]=fits(s, 'lorz', p, fix);
Qf(ii) = fbck.pvals(2)/fbck.pvals(3); %Calculate the quality factor
chi(ii) = 1/Qf(ii);
disp(num2str(Hx(ii),'Current magnetic field: %3.2f.'));
end

figure
freqPlot = plot(Hx,ff0,'or','MarkerSize',2);
xlabel('Field (T)');
ylabel('Frequency (GHz)');
title(num2str(Temperature,'Resonant frequency from fitted data at T = %3.3f K'));
axis([field_l field_h freq_l freq_h]);

%     yq = yq.*0 + 3.49;
%     zq = Fmag(xq,yq);
%     h1(n) = plot(xq(1,:),zq); Error: the indices on the left side are not compatible with the size of the right side.
%     hk = plot(xq(1,:),Q);
%     plot(xq(1,:),zq);
%         col = setcolours(n/length(T),'jet');
figure
Hx = Hx(Qf >=0);
Qf = Qf(Qf >=0); %Remove unphysical points
Qf = medfilt1(Qf, order);
Qplot1 = plot(Hx, Qf,'o-','MarkerSize',4);
hold on
H0 = H0(Q0 >=0);
Q0 = Q0(Q0 >=0);
Q0 = medfilt1(Q0,10);
Qplot2 = plot(H0, Q0,'s-','MarkerSize',4);
gca;
xlabel('Field (T)');
ylabel('Q factor');
legend('Quality factor from Lorentzian fit', 'Quality factor from FWHM');
title('Quality factor');
%         set(hk(n),'color',col,'linewidth',plotopt.lnwd);

% set(ax,'fontsize',plotopt.ftsz,'xlim',[field_l field_h])
%
% set(ax(2),'ylim',[0.2 0.6])
% tl = legend(h1,num2str(T(:),'%3.1f K'),'location','se');
% set([tx ty tl tt],'fontsize',plotopt.ftsz,'interpreter','latex')
%     tt(2) = title([filename num2str(runs(n)) ', f = 3.437 GHz']);

% saveplots(hfig(2),'7-Temperature_field-S_-16dBm_Q')



end
function option3(filepath, filename)
clear freq S11 dB N FdB FrS FiS FTT1 FTT2

%% Plot data

%Set data range and parameters
Temperature = 0.080;
FWHM = -3.0; % Define full-width-half-max for calculation of quality factors

clear freq S11 dB N FdB FrS FiS FTT1 FTT2

out = readdata_v3(filepath,filename,1);
freq = out.data.ZVLfreq/1e9;
S11 = out.data.ZVLreal + 1i*out.data.ZVLimag;
mag = (1- out.data.ZVLreal.^2 -out.data.ZVLimag.^2)./((1- out.data.ZVLreal).^2 + out.data.ZVLimag.^2);    
H = out.data.DCField1;
T1 = out.data.Temperature1;
HH = repmat(H,1,size(freq,2)); %populate the magnetic field to match the dimension of S11 matrix

% Determing the fieldscan direction
if H(end)-H(1)>0
    direction = 'Up_';
elseif H(end)-H(1)<0
    direction = 'Down_';
else
    error('Could not determing fieldscan direction!');
end
excitation = '_-20dBm_0dB';

%Turn all coordinates into 1D vectors
freq = freq';
freq = freq(:);
S11 = S11';
S11 = S11(:);
HH = HH';
HH = HH(:);

freq_l = min(freq); %set frequency range, l: lower limit, h: higher limit
freq_h = max(freq);

% Use the average of the subsequent difference to caculate frequency scan
% step
dif = nonzeros(diff(freq));
dif = dif(dif>0);
step = mean(dif);
nop = ceil((freq_h-freq_l)/step)+1; %compute how many points pers complete frequency scan.
clearvars dif

%Clean up the raw data by removing incomplete scans (step 1) and duplicates
%(step 2)
% step 1: truncate the beginning and end part to keep only complete frequency scans and reshape the matrices into single column vectors
trunc1 = find(freq==freq_l,1,'first'); 
trunc2 = find(freq==freq_h,1,'last'); 
S11 = S11(trunc1:trunc2);
dB = mag2db(abs(S11));
freq = freq(trunc1:trunc2);
HH = HH(trunc1:trunc2);

% Step 2: remove duplicates
dupl = find(diff(freq) == 0);
freq(dupl+1)=[];
dB(dupl+1)=[];
HH(dupl+1)=[];

dB = reshape(dB,nop,[]);  %reshape the matrix so that each complete frequency scan occupy one column
freq = reshape(freq,nop,[]);
HH = reshape(HH,nop,[]);
%Find all the resonant peaks
HH = HH';
dB = dB';
freq = freq';

f0 = zeros(size(dB,2),1);
H0 = zeros(size(dB,2),1);
Q0 = zeros(size(dB,2),1);
dB0 = zeros(size(dB,2),1);
%find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
for ii = 1:size(dB,2) %Searching column minima (fixed field) is better than searching row minima (fixed frequency)
    [~,idx] = min( dB(:,ii) );
    if(length(idx)>1)
        disp('multiple minima found')
    end
    f0(ii) = freq(idx,ii);
    H0(ii) = HH(idx,ii); 
    dB0(ii) = dB(idx,ii);
    % Calculate quality factor using f0/FWHM
    if isnan(1/range(freq(dB(:,ii) <= FWHM)))
       Q0(ii) = 0;
    elseif isempty(range(freq(dB(:,ii) <= FWHM)))
       Q0(ii) = 0;
    else
    Q0(ii) = freq(idx,ii)/range(freq(dB(:,ii) <= FWHM));
    end
    Q0(isnan(Q0)) = 0; % Cut out NaN from the array
    Q0(isinf(Q0)) = 0; % Cut out inf from the array
    Q0(Q0>4000) = 0; % Cut out unrealistic quality factor from noise data
end
%   For noisy data, we need to remove duplicates of minima
[H0,ia,~] = unique(H0,'stable');
dB0 = dB0(ia);

cd(filepath);
tit=[direction,num2str(Temperature,'%3.3f'), excitation,'.mat'];
save(tit,'H0','dB0');

end
%% --------------------------------------------------------------------------------
function hfig = setfig(nfig)

hfig = figure(nfig);
clf
pos = get(hfig,'position');
set(hfig,'position',[pos(1:2) 700 500])

ha = get(hfig,'currentAxes');
if isempty(ha)
    ha = axes('position',[0.15 0.15 0.7 0.7]);
else
    set(ha,'position',[0.15 0.15 0.7 0.7])
end
end

function saveplots(hfig,figname)
% Function to save the figures to a specified directory and in different
% file formats
% INPUT:
% hfig      figure handle
% figname   file name to save figure to
curdir = cd;
cd('G:\My Drive\File sharing\PhD projects\LiHoF4\Data\Matlab\Cavity resonator\D24mm_T5mm_G0.2mm\12.02.2019')
%cd('/Users/yikaiyang/Google Drive/File sharing/PhD projects/LiHoF4/Data/Matlab/Cavity resonator/D24mm_T5mm_G0.2mm/12.02.2019')

%saveas(figure(hfig),[figname '.fig'],'fig');
%     print(figure(hfig),[figname  '.jpg'],'-djpeg','-r600');
%     print(figure(hfig),[figname  '.png'],'-dpng','-r600');
%print2eps(figname,hfig)
%[~,~] = eps2xxx([figname '.eps'],{'jpeg','pdf'});

disp(['Figure ' figname ' saved to '])
disp(cd)
cd(curdir)
end
