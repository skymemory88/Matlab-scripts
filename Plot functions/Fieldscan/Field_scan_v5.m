format long;
% For macOS
%     addpath('/Volumes/GoogleDrive/My Drive/File sharing/Programming scripts/Matlab/Plot functions/Fieldscan/functions')
%     addpath(genpath('/Volumes/GoogleDrive/My Drive/File sharing/Programming scripts/Matlab/Plot functions/spec1d--Henrik'));
%     location = '/Volumes/GoogleDrive/My Drive/File sharing/PhD program/Research projects/LiHoF4 project/Data/Experiment/LiHoF4/SC200/2021.04.12';

% For windows
addpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\Fieldscan\functions\');
addpath(genpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\spec1d--Henrik\'));
location = ['G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\',...
             'SC239\2021.06.26'];
loadname = '2021_06_0042';
LoadObj = fullfile(location,[loadname, '.dat']);
SaveObj = fullfile(location,[loadname, '_interp', '.mat']);
Options.analysis = 1; % Analysis options

% Figure plot options:
Options.plot = false; % Plot option in analysis part (Option 3)
Options.lnwd = 1.5; % plot linewidth
Options.ftsz = 12; % plot font size
Options.mksz = 3; % plot marker size
Options.order = 4; % 1D filter order
Options.bgdmode = 1; % Background normalization (0: no normalization. 1: normalized with stitched background. 2: normalize with loaded file)
Options.fitfunc = 1; % Pick fitting function from either (1) custom function of (2) spec1d
Options.nData = 1; % Number of dataset from VNA

switch Options.analysis
    case 1
        % Simple color plot of S11 (+ S21)
        [fields,freq,S11,analysis] = option1(LoadObj, Options);
    case 2
        % On-resonance measurement processing (w/ plots)
        [fields,freq,S11,analysis] = option2(LoadObj, Options);
    case 3
        % Off-resonance measurement processing
        [fields,freq,S11,analysis] = option3(LoadObj, Options);
    case 4
        % Temperature scan
        [fields,freq,S11,analysis] = option4(LoadObj, Options);
end
clearvars -except fields freq S11 analysis Options SaveObj

prompt = sprintf('Save the analysis?\n');
answer = input(prompt,'s');
switch lower(answer)
    case {'y','yes'}
        save(SaveObj,'analysis','fields','freq','S11','-v7.3')
%         if Options.analysis == 4 % Add a point to the phase diagram from off-resonance measurement
%             phase(1,1) = Temperature;
%             phase(1,2) = mean(H0(f0==min(f0)));
%             if Options.savefile == true && isfile(dataobj)
%                 save(dataobj,'phase','-append');
%             elseif Options.savefile == true
%                 save(dataobj,'phase','-v7.3');
%             end
%         end
    otherwise
        disp('Aborted saving analysis')
end
clearvars prompt SaveObj

function [xq,yq,zq,analysis] = option1(LoadObj, Options)
out = readdata_v4(LoadObj,Options.nData);
freq = out.data.ZVLfreq/1e9;
S11 = out.data.ZVLreal + 1i*out.data.ZVLimag;
H = out.data.DCField1;
HH = repmat(H,1,size(freq,2)); %populate the magnetic field to match the dimension of S11 matrix
T1 = out.data.Temperature1; % Sample cell temperature
Temperature = mean(T1(H == min(H)));
analysis.temp = T1;
% T2 = out.data.Temperature2; % Mixing chamber temperature
lck1 = out.data.lockin1;
lck2 = out.data.lockin2;
phaseP = false; % Option to plot phase

S11 = S11';
S11 = S11(:);
freq = freq';
freq = freq(:);
HH = HH';
HH = HH(:); %the third argument is the number of frequency points in each line/segment

[rows,~,freq] = find(freq); % Remove nonsensical zeros from the frequency data
HH = HH(rows);
S11 = S11(rows);

%Set data range and parameters
freq_l = min(freq); %set frequency range, l: lower limit, h: higher limit
freq_h = max(freq);
field_l = 2.6;  %set field range, l: lower limit, h: higher limit
field_h = max(H);

%% Code For R&S ZVL-6
% Clean up the raw data by removing incomplete frequency scans
trunc1 = find(freq == min(freq),1,'first'); 
trunc2 = find(freq == min(freq),1,'last')-1; 
freq_temp = freq(trunc1:trunc2);
HH_temp = HH(trunc1:trunc2);
S11_temp = S11(trunc1:trunc2);

% Truncate according to the upper and lower field limit
S11_temp = S11_temp(HH_temp>=field_l & HH_temp<=field_h);
freq_temp = freq_temp(HH_temp>=field_l & HH_temp<=field_h);
HH_temp = HH_temp(HH_temp>=field_l & HH_temp<=field_h);

% Truncate according to the upper and lower frequency limit
S11_temp = S11_temp(freq_temp>=freq_l & freq_temp<=freq_h);
freq_temp = freq_temp(freq_temp>=freq_l & freq_temp<=freq_h);
HH_temp = HH_temp(freq_temp>=freq_l & freq_temp<=freq_h);

% reasign actual lower and upper limits of the frequency and field
freq_l = min(freq_temp);
freq_h = max(freq_temp);
field_l = min(HH_temp);
field_h = max(HH_temp);
[xq,yq] = meshgrid(linspace(field_l,field_h,1501),linspace(freq_l,freq_h,1201)); %set the X and Y range

dif = diff(freq); % frequency increments
resets = find(dif <= 1*(freq_l-freq_h)); % Find the termination point of a complete scans (5% tolerance)
nop = round(mean(diff(resets))); % Set the number of points per frequency scan

mag_temp = abs(S11_temp);
mag = reshape(mag_temp,nop,[]);
freq = reshape(freq_temp,nop,[]);
HH = reshape(HH_temp,nop,[]);
clearvars dif step mag_temp freq_temp HH_temp

if Options.bgdmode ~= 0
    %Find all the resonant peaks
    f0 = zeros(size(mag,2),1);
    H0 = zeros(size(mag,2),1);
    mag0 = zeros(size(mag,2),1);
%     find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
%     [~,bidx] = min( HH(1,:) ); % find the column index for zero-field frequency scan
%     [~,bidx] = max( HH(1,:) ); % find the column index for high-field frequency scan
    [bidx,~,~] = find(mag == min(mag( HH(1,:)==min(HH(1,:)) ),mag( HH(1,:) == max(HH(1,:)) )) ); % pick the finner peak as the first part of the background
    for ii = 1:size(mag,2) %Searching column minima (fixed field)
        [~,idx] = min( mag(:,ii) );
        if(length(idx)>1)
            disp(num2str(H0(ii),'multiple minima found at H = %.3f'))
        end
        f0(ii) = freq(idx,ii);
        H0(ii) = HH(idx,ii);
        mag0(ii) = mag(idx,ii);
    end
    % Stich the background for normalization
    [~,cidx] = min(mag(:,bidx)); % find the resonant peak to be replaced
%     [~,Hpos] = max(mag0); % Method 1. scan of max reflection
%     [~,Hpos] = max(mag(cidx,:)); % Method 2. scan of max reflection at f0 (frequency)
    [~,Hpos] = max(medfilt1(abs(f0-f0(bidx)))); % Method 3. scan with the resonant peak deviates from f0 the furtherest
%     [~,Hpos] = min(abs(H0-1.92)); % Method 4. Manually set the patch slice of the frequency scan
    mag_dif = abs(mag(:,bidx)-mag(:,Hpos));
    [lidx,~] = find(mag_dif(1:cidx) <= 1e-3, 1, 'last'); % left boundary of the resonant peak to be replaced
    [ridx,~] = find(mag_dif(cidx:end) <= 1e-3, 1, 'first'); % right boundary of the resonant peak to be replaced
    if isempty(lidx); lidx = 1; end
    if isempty(ridx); ridx = length(mag_dif);else; ridx = cidx + ridx -1; end
    analysis.Hpos = Hpos; % store field location of the anti-crossing
    analysis.cidx = cidx; % peak center
end

if Options.bgdmode == 1% Construct the background noise by stitching together zero-field-scan and anti-crossing
    bf0 = freq(:,bidx);
    bgd0 = mag(:,bidx); % Base frequency scan
    analysis.bf0 = bf0;
    analysis.bgd0 = bgd0;
    bgd0(lidx:ridx) = mag(lidx:ridx,Hpos);
    figure
    plot(bf0,mag(:,bidx))
    hold on
    plot(freq(:,Hpos),mag(:,Hpos))
    plot(bf0,bgd0);
    xlabel('Frequency (GHz)')
    ylabel('S11')
    legend(num2str(H0(bidx),'B = %.2f T'),num2str(H0(Hpos),'B = %.2f T'),'Stitched')
    [~,trimIdx] = unique(bf0);
    bgd0 = bgd0(trimIdx);
elseif Options.bgdmode == 2% load background data from existing file.
%     backpath = ('/Volumes/GoogleDrive/My Drive/File sharing/PhD program/Research projects/LiHoF4 project/Data/Experiment/Cavity resonator/D24mm_T5mm_G_flex/SuperCoax calibration');
    backpath = ('G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data/Experiment\Cavity resonator\D24mm_T5mm_G_flex\SuperCoax calibration');
    backfile = ('background.mat');
    load(fullfile(backpath,backfile),'bf0','bgd0');
    figure
    plot(bf0,bgd0);
    xlim([freq_l freq_h])
    xlabel('Frequency (GHz)')
    ylabel('S11')
    legend('Background')
    [~,trimIdx] = unique(bf0);
    bgd0 = bgd0(trimIdx)+2.5; % Shift the loaded background to the noise floor of current data
else
    disp('No background normalization')
end

% Renormailize the background
if Options.bgdmode ~= 0
    mag = mag(trimIdx,:); % Remove duplicate points
    freq = freq(trimIdx,:);
    HH = HH(trimIdx,:);
    bgdM = repmat(bgd0,1,size(mag,2));
    mag = mag./bgdM;
end
% Interpolate the 3D data
Fmag = scatteredInterpolant(HH(:),freq(:),mag(:));
zq = Fmag(xq,yq);
clearvars c idx ia lidx ridx widx ii HM trunc1 trunc2 dupl nop bf0

% Plot additional data if there are more than one set of data from VNA
if Options.nData > 1
    S21 = out.data.ZVLreal2 + 1i*out.data.ZVLimag2;
    S21 = S21';
    S21 = S21(:);
    mag2 = abs(S21);
    FaS2 = scatteredInterpolant(HH,freq,mag2);
    zq2 = FaS2(xq,yq);
    
    % Plot the interpolated frequency response data in a field scan using color map
    figure
    cmap = pcolor(xq,yq,mag2db(zq2));
    set(cmap, 'edgeColor','none')
    shading interp;
    colorbar
    set(gca,'fontsize',Options.ftsz)
    xlabel('Field (T)');
    ylabel('Frequency (GHz)');
    xticks(linspace(field_l,field_h,6));
    title(num2str(Temperature,'S21 response at T = %3.3f K'));
end

% Plot the interpolated frequency response data in a field scan using color map
figure
cmap = pcolor(xq,yq,mag2db(zq));
set(cmap, 'edgeColor','none')
shading interp;
caxis([-20 1])
colorbar
set(gca,'fontsize',Options.ftsz)
xlabel('Field (T)');
ylabel('Frequency (GHz)');
xticks(linspace(field_l,field_h,6));
title(num2str(Temperature,'S11 (dB) response at T = %3.3f K'));

% Plot the interpolated imaginary part of the response function in a field scan using color map
if phaseP == true
    FrS = scatteredInterpolant(HH_temp,freq_temp,real(S11_temp)); %intrapolate the real part of S11 response on a 2D grid
    FiS = scatteredInterpolant(HH_temp,freq_temp,imag(S11_temp)); %intrapolate the imaginary part of S11 response on a 2D grid
    zq1 = FrS(xq,yq) + 1i*FiS(xq,yq);
    figure
    cmap1 = pcolor(xq,yq,imag(zq1));
    set(cmap1, 'edgeColor','none')
    shading interp;
    % caxis([-30 0])
    colorbar
    set(gca,'fontsize',Options.ftsz)
    xlabel('Field (T)');
    ylabel('Frequency (GHz)');
    xticks(linspace(field_l,field_h,6));
    title(num2str(Temperature,'Imaginary part of S11 at T = %3.3f K'));
end

% Plot lock-in data if it isn't empty
if any(lck1)
    figure
    plot(H(1:100:end),lck1(1:100:end),'s-')
    xlabel('DC Magnetic field')
    ylabel('Hall resistence')
    title('Hall resistance vs Field')
end

if any(lck2)
    figure
    plot(H(1:100:end),lck2(1:100:end),'x-')
    xlabel('DC Magnetic field')
    ylabel('Sample resistence')
    title('Sample resistance vs Field')
end

% Plot the temperature profile against magnetic field
figure
% plot(H(1:100:end),T1(1:100:end),'o-')
plot(H,T1,'o-')
xlabel('DC Magnetic field')
ylabel('Temperature')
title('Magnetic field vs Temperature')
end
% End of option 1
function [xq,yq,zq,analysis] = option2(LoadObj, Options)
%Set data range and parameters
order = Options.order; % set to what order the median filters is applied
clear freq S11 dB N FdB FrS FiS FTT1 FTT2

% extract data from raw data file
if ~exist('out','var')
    out = readdata_v4(LoadObj, Options.nData);
end

freq = out.data.ZVLfreq/1e9;
S11 = out.data.ZVLreal + 1i*out.data.ZVLimag;  
H = out.data.DCField1;
T1 = out.data.Temperature1;
% T2 = out.data.Temperature2;

Temperature = min(T1(H == min(H))); %Extract the measurement temperature at the lowest field (to avoid magnetoresistance in the thermometer)
analysis.temp = T1;
HH = repmat(H,1,size(freq,2)); %populate the magnetic field to match the dimension of S11 matrix

freq = freq';
freq = freq(:);
S11 = S11';
S11 = S11(:);
HH = HH';
HH = HH(:);

[rows,~,freq] = find(freq); % Remove nonsensical zeros from the frequency data
HH = HH(rows);
S11 = S11(rows);

freq_l = min(freq); % set frequency range, l: lower limit, h: higher limit
freq_h = max(freq);
field_l = min(H);  % set field range, l: lower limit, h: lower limit
field_h = max(H);
field_cut = field_l + 2.5; % Field window for cavity parameter fit

%Plot the temperature vs magnetic field to check the temperature variation
figure
plot(H(1:round(length(T1)/100):end),T1(1:round(length(T1)/100):end),'o-')
xlabel('DC Magnetic field')
ylabel('Temperature')
title('Magnetic field vs Temperature')

%% Code For R&S ZVL-6
% Clean up the raw data by removing incomplete scans (step 1) and duplicates(step 2)
% step 1: truncate the beginning and end part to keep only complete frequency scans and reshape the matrices into single column vectors
trunc1 = find(freq == min(freq),1,'first'); 
trunc2 = find(freq == min(freq),1,'last')-1; 
freq_temp = freq(trunc1:trunc2);
HH_temp = HH(trunc1:trunc2);
S11_temp = S11(trunc1:trunc2);

% Truncate according to the upper and lower field limit
S11_temp = S11_temp(HH_temp>=field_l & HH_temp<=field_h);
freq_temp = freq_temp(HH_temp>=field_l & HH_temp<=field_h);
HH_temp = HH_temp(HH_temp>=field_l & HH_temp<=field_h);

% Truncate according to the upper and lower frequency limit
S11_temp = S11_temp(freq_temp>=freq_l & freq_temp<=freq_h);
freq_temp = freq_temp(freq_temp>=freq_l & freq_temp<=freq_h);
HH_temp = HH_temp(freq_temp>=freq_l & freq_temp<=freq_h);

freq_l = min(freq_temp); %set frequency range, l: lower limit, h: higher limit
freq_h = max(freq_temp);
field_l = min(HH_temp);  %set field range, l: lower limit, h: higher limit
field_h = max(HH_temp);

% Interpolate the data on a 2D grid for the colormap
[xq,yq] = meshgrid(linspace(field_l,field_h,1001),linspace(freq_l,freq_h,801));

dif = diff(freq_temp); % frequency increments
resets = find( dif <= (freq_l-freq_h) ); % Find the termination point of a complete scans
nop = round(mean(diff(resets))); % Set the number of points per frequency scan

mag_temp = abs(S11_temp);
mag = reshape(mag_temp,nop,[]);
freq = reshape(freq_temp,nop,[]);
HH = reshape(HH_temp,nop,[]);
%% Temperary code for Keysight PNA-x N5242B
% S11_temp = S11;
% dB_temp = mag2db(abs(S11_temp));
% dB_temp = dB_temp - max(dB_temp,[],'all'); % Shift the data to compensate changes in RE(impedence) at low temperatures
% freq_temp = freq;
% HH_temp = HH;
% 
% dif = diff(freq); % frequency increments
% resets = find(dif<=0.9*(freq_l-freq_h)); % Find the termination point of a complete scans (10% error)
% nop = round(mean(diff(resets))); % Set the number of points per frequency scan
% dB = reshape(mag2db(abs(S11)),nop,[]);  %reshape the matrix so that each complete frequency scan occupy one column
% dB = dB - max(dB,[],'all'); % Shift the data to compensate changes in RE(impedence) at low temperatures
% freq = reshape(freq,nop,[]);
% HH = reshape(HH,nop,[]);
%% Clean up raw data
clearvars dif step freq_temp S11_temp mag_temp HH_temp dB_temp

if Options.bgdmode == 1% Construct the background noise by stitching together zero-field-scan and anti-crossing
    %Find all the resonant peaks
    f0 = zeros(size(mag,2),1);
    H0 = zeros(size(mag,2),1);
    mag0 = zeros(size(mag,2),1);
    %find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
%     [~,bidx] = min( HH(1,:) ); % find the column index for zero-field frequency scan
%     [~,bidx] = max( HH(1,:) ); % find the column index for high-field frequency scan
    [bidx,~,~] = find(mag == min( mag(unique(HH(1,:)==min(HH(1,:)))), mag(unique(HH(1,:) == max(HH(1,:)))) ) ); % pick the finner peak as the first part of the background
    for ii = 1:size(mag,2) %Searching column minima (fixed field)
        [~,idx] = min( mag(:,ii) );
        if(length(idx)>1)
            disp(num2str(H0(ii),'multiple minima found at H = %.3f'))
        end
        f0(ii) = freq(idx,ii);
        H0(ii) = HH(idx,ii);
        mag0(ii) = mag(idx,ii);
    end 
    [~,cidx] = min(mag(:,bidx)); % find the resonant peak to be replaced
%     [~,Hpos] = max(medfilt1(mag0)); % Method 1. scan of max reflection
%     [~,Hpos] = max(medfilt1(mag(cidx,:))); % Method 2. scan of max reflection at f0 (frequency)
    [~,Hpos] = max(medfilt1(abs(f0-f0(bidx)))); % Method 3. scan with the resonant peak deviates from f0 the furtherest
%     [~,Hpos] = min(abs(H0-3.465)); % Method 4. Manually set the patch slice of the frequency scan
    
    mag_dif = abs(mag(:,bidx)-mag(:,Hpos));
    [lidx,~] = find(mag_dif(1:cidx) <= 1e-3, 1, 'last'); % left boundary of the resonant peak to be replaced
    [ridx,~] = find(mag_dif(cidx:end) <= 1e-3, 1, 'first'); % right boundary of the resonant peak to be replaced
    if isempty(lidx); lidx = 1; end
    if isempty(ridx); ridx = length(mag_dif);else; ridx = cidx + ridx -1; end
    
    bgd0 = mag(:,bidx); % zero-field frequency scan
    bgd0(lidx:ridx) = mag(lidx:ridx,Hpos);
    bf0 = freq(:,bidx);
    figure
    plot(freq(:,bidx),mag(:,bidx))
    hold on
    plot(freq(:,Hpos),mag(:,Hpos))
    plot(freq(:,bidx),bgd0);
    xlabel('Frequency (GHz)')
    ylabel('S11')
    legend('B = 0',num2str(H0(Hpos),'B = %.2f T'),'Stitched')
    disp('Press any key to continue, or Ctrl + c to abort.')
    pause
elseif Options.bgdmode == 2 % load background data from existing file.
    %     backpath = ('/Volumes/GoogleDrive/My Drive/File sharing/PhD program/Research projects/LiHoF4 project/Data/Experiment/Cavity resonator/D24mm_T5mm_G_flex/SuperCoax calibration');
    backpath = ('G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data/Experiment\Cavity resonator\D24mm_T5mm_G_flex\SuperCoax calibration');
    backfile = ('background.mat');
    load(fullfile(backpath,backfile),'bf0','bgd0');
    figure
    plot(bf0,bgd0);
    xlim([freq_l freq_h])
    xlabel('Frequency (GHz)')
    ylabel('S11')
    legend('Background')
    [bf0,trimIdx] = unique(bf0);
    bgd0 = sbgd0(trimIdx);
    bgd0 = interp1(bf0,bgd0,freq(:,1));
    bf0 = freq(:,1);
    hold on
    plot(bf0,bgd0);
    legend('Loaded data','Interpolated data');
    disp('Press any key to continue.')
    pause
else
    disp('No background normalization')
    bf0 = freq(:,1);
    bgd0 = ones(size(freq(:,1),1),1);
end
%Find all the resonant peaks
f0 = zeros(size(mag,2),1);
Q0 = zeros(size(mag,2),1);
FWHM0 = zeros(size(mag,2),1);
%find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
for ii = 1:size(mag,2) %Searching column minima (fixed field)
    [~,idx] = min( mag(:,ii)./bgd0 );
    if(length(idx)>1)
        disp(num2str(H0(ii),'multiple minima found at H = %.3f'))
    end
    f0(ii) = freq(idx,ii);
    H0(ii) = HH(idx,ii); 
    mag0(ii) = mag(idx,ii);
    HM = ( max(mag(:,ii)) + min(mag(:,ii)) )/2; % 
    % Calculate quality factor using f0/FWHM
    if isnan(1/range(freq(mag(:,ii) <= HM)))
       Q0(ii) = 0;
    elseif isempty(range(freq(mag(:,ii) <= HM)))
       Q0(ii) = 0;
    else
       Q0(ii) = freq(idx,ii)/range(freq(mag(:,ii) <= HM));
       FWHM0(ii) = range(freq(mag(:,ii) <= HM));
    end
end
Q0(isinf(Q0)) = NaN; % Cut out inf from the array
[Q0,c] = rmmissing(Q0); % Cut out NaN from the array
H0 = H0(~c); % Remove corresponding elements in H0 array as well
f0 = f0(~c);
mag0 = mag0(~c);
FWHM0 = FWHM0(~c);

% For noisy data, we need to remove duplicates of minima
[H0,ia,~] = unique(H0,'stable');
f0 = f0(ia);
Q0 = Q0(ia);
mag0 = mag0(ia);
FWHM0 = FWHM0(ia);

f0 = medfilt1(f0,order); % apply median filter to remove some noise
f0 = f0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
H0 = H0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
Q0 = Q0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
mag0 = mag0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
FWHM0 = FWHM0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints

% Fit the field dependent resonant frequency data with weak coupling function
spin = -5.5;
hPara = [H0(Hpos), field_l, field_h];
fPara = [(freq_l+freq_h)/2, freq_l, freq_h];
[fitP,~] = wk_cpl_fit(H0,f0,spin,hPara,fPara);
spin = fitP.spin;
wc = fitP.wc;
% wc = 3.645;
gc = fitP.g;
% Brf = fitP.x0; % Level crossing location from perturbative fitting
% Brf = H0(Hpos); % Level crossing location from minimum search
Brf = 3.41;

% Resonant frequency trace with fitting parameters
figure
B = linspace(field_l,field_h,100);
Delt = spin*(B-Brf); % Assume linear relation near the line cross
plot(H0(1:round(length(H0)/200):end),f0(1:round(length(f0)/200):end),'ok','MarkerSize',4);
hold on
wp = wc + Delt./2 + sqrt(Delt.^2 + 4*gc^2)/2;
wm = wc + Delt./2 - sqrt(Delt.^2 + 4*gc^2)/2;
plot(B,wm,'-r',B,wp,'-r','LineWidth',Options.lnwd);
plot(B,Delt + wc,'--k');
% hfig1 = plot(H0, f0, 'o', 'MarkerSize', 2);
xlabel('Field (T)');
ylabel('Resonant frequency (GHz)');
title(num2str(Temperature,'Resonant frequency from minimum search at T = %3.3f K'));
legend(sprintf('gc = %.3f GHz',fitP.g/2));
axis([field_l field_h freq_l freq_h]);

% % Plot frequency scan at line crossing
% figure
% % plot(freq(1:20:end,Hpos),mag(1:20:end,Hpos),'.');
% plot(freq(:,Hpos),(mag(:,Hpos)./bgd0),'.');
% xlabel('Frequency (GHz)');
% ylabel('S11');
% legend(sprintf('Frequency cut at %.2f T',H0(Hpos)));
% title('Frequency scan at line crossing');
clearvars B Delt hPara fPara wp wm fitPara H_res f_res

% Remove duplicate points
[~,trimIdx] = unique(bf0);
bgd0 = bgd0(trimIdx);
mag = mag(trimIdx,:);
freq = freq(trimIdx,:);
HH = HH(trimIdx,:);
mag0 = interp1(H0,mag0,xq(1,:)); % amplitude along resonant frequency trace
[~,Hpos] = max(mag0); % find the line crossing position on field axis along interpolated data

% Renormailize the background
if Options.bgdmode ~= 0
    bgdM = repmat(bgd0,1,size(mag,2));
    mag = mag./bgdM;
end
clearvars c idx ia lidx ridx widx ii HM T1 trunc1 trunc2 dupl nop bf0

switch '2D' %choose data interpolation method and plot the color map of S11 response
    case '1D' % Option_1 1D data interpolation 
        interp_amp = zeros(size(yq,1),size(HH,2));
        parfor ii = 1:size(HH,2)
            interp_amp(:,ii) = interp1(freq(:,ii),mag(:,ii),yq(:,1));
        end
        zq = interp_amp'; % For later use in lorentzian fitting
        HH = repmat(H0,1,size(xq,1));
        
        switch 1 % Choose plotting option
            case 1
                % plot single-direction interpolated data using pseudo-colormap
                cmap0 = pcolor(HH(:,1),yq(:,1),mag2db(interp_amp));
                set(cmap0, 'edgeColor','none')
                shading interp;
                colorbar
            case 2
                % plot single-direction interpolated data using scatter plot
                yy = repmat(yq(:,1),1,length(H0));
                yy = yy(:);
                interp_amp = interp_amp(:);
                HH = HH';
                HH = HH(:);
                cmap0 = scatter3(HH,yy,interp_amp,2,mag2db(interp_amp),'o','filled','MarkerEdgeColor','none');
                colormap(hsv);
        end
    case '2D' % Option_2 2D data interpolation
        Fmag = scatteredInterpolant(HH(:),freq(:),mag(:)); % 2D interpolation of amplitude data
        zq = Fmag(xq,yq);     
        figure
        hold on
        box on   
        cmap0 = pcolor(xq,yq,mag2db(zq));
        set(cmap0, 'edgeColor','none')
        shading interp;
        colorbar
        caxis('auto')
%         caxis([max([mag2db(min(mag0)), -30]) max(mag2db(max(mag0)),0)+1]);
end
clearvars *_temp;
        
set(gca,'fontsize',Options.ftsz)
axis([field_l field_h freq_l freq_h]);
xlabel('Field (T)');
ylabel('Frequency (GHz)');
title(num2str(Temperature,'S11 response at T = %3.3f K'));

% Use interpolated data to extract quality factor
Hx = xq(1,:);
FWHM = interp1(H0,FWHM0,Hx,'pchip','extrap'); % Using spline interpolation to smooth the FWHM
omg = spin.*(Hx-Brf) + wc; % linearized dispersion relation
            
ff0 = double.empty(length(Hx),2,0); 
Qf = double.empty(length(Hx),0);
gamma = double.empty(length(Hx),0);
Gc = double.empty(length(Hx),0);
zq(isnan(zq))=0;

weight = double.empty(size(zq,1),size(zq,2),0);
for jj = 1:length(Hx)
    for ii = 1:length(zq(:,jj))
        weight(ii,jj,1) = abs(1./zq(ii,jj)); % Weight function option 1
    end
end
% weight = 1./gradientweight(zq); % Weight function option 2

%% Step 1: Fit the data for the first time to extract "kpe" and "w0" far from the level crossing
Hc0 = find(Hx >= field_l, 1,'first');
Hc1 = find(Hx <= field_cut, 1, 'last');
kpe = double.empty(length(Hx),0);
kpi = double.empty(length(Hx),0);
parfor ii = Hc0:Hc1
    if mod(ii,20) == 0
        worker = getCurrentTask();
        fprintf('First fit, current magnetic field: %1$3.2f. on core %2$u.\n', Hx(ii), worker.ID);
    end
% for ii = Hc0:Hc1 % for debugging
%     % Selectively plot the fit for visual inspection
%     figWin = [1:round((Hc1-Hc0)/10):Hc1]; % The iteration window over which shows fitting graph
%     if ismember(ii,figWin) % Not useable in parallel mode (parfor)
%         plt = true;
%     end
    plt = false;
    % starting value for param = {'kpe', 'kpi', 'omega', 'Gc', 'gma', 'w0', 'attn'}
    param = [FWHM(ii)  FWHM(ii)  omg(ii)  gc  gc/1000  wc 1/mean(zq(:,ii))];
    bound_l = [0   0  omg(ii)  0   gc/1000  freq_l  0.95]; % lower bound of fitting parameters
    bound_h = [1   1  omg(ii)  1   gc/1000  freq_h  1.05]; % upper bound of fitting parameters
    fit = iptopt(yq(:,ii),zq(:,ii),Hx(ii),param,bound_l,bound_h,weight(:,ii),plt);
    
    param = coeffvalues(fit);
    kpe(ii) = param(1);
    kpi(ii) = param(2);
    Gc(ii) = param(4);
    ff0(ii,:,1) = param(6);
end
%% Step 2: fit the data for the second time with fixed "kpe" and "w0"
kpe0 = mean(kpe);
kpi0 = mean(kpi); % keep constant kpe/kpi ratio
wc = mean(ff0(1:Hc1,1));
Hc2 = find(Hx >= field_cut, 1,'first');
Hc3 = find(Hx <= field_h, 1, 'last');
% parfor ii = Hc2:Hc3
%     if mod(ii,20) == 0
%         worker = getCurrentTask();
%         fprintf('Second fit, current magnetic field: %1$3.2f. on core %2$u.\n', Hx(ii), worker.ID);
%     end
for ii = Hc2:Hc3
    % Selectively plot the fit for visual inspection
    plt = false;
    figWin = Hpos-8:2:Hpos+18; % The iteration window over which shows fitting graph
    if ismember(ii,figWin) % Not useable in parallel mode (parfor)
        plt = true;
    end
    % starting value for param = {'kpe', 'kpi', 'omega', 'Gc', 'gma', 'w0', 'attn'}
    param  =  [kpe0  kpi0  omg(ii)   gc   gc/1000   wc  1/mean(zq(:,ii))];
    bound_l = [ 0    kpi0     0      0    0    wc-0.1   0.95]; % lower bound of fitting parameters
    bound_h = [inf   kpi0    inf     1    1    wc+0.1   1.05]; % upper bound of fitting parameters
    fit = iptopt(yq(:,ii),zq(:,ii),Hx(ii),param,bound_l,bound_h,weight(:,ii),plt);
    
    param = coeffvalues(fit);
    kpe(ii) = param(1);
    kpi(ii) = param(2);
    Gc(ii) = param(4);
    gamma(ii) = param(5);
    
    [~,idx] = findpeaks(-fit(yq(:,ii)));
    if isempty(idx)
        ff0(ii,:) = yq(fit(yq(:,ii)) == min(fit(yq(:,ii))));  % Find the resonant frequency by minimum search
    else
        ff0(ii,:) = yq(idx,ii);  % Find the resonant frequency by minimum search
    end
    Qf(ii) = mean(ff0(ii,:)./(param(1)+param(2))); % Calculate quality factor
end
analysis.kpe = kpe;
analysis.kpi = kpi;
analysis.w0 = ff0;
analysis.Gc = Gc;
analysis.gamma = gamma;

figure
plot(Hx(Hc0:Hc3),kpe(Hc0:Hc3),'-');
ylabel('K_e (GHz)');
hold on
yyaxis right
plot(Hx(Hc0:Hc3),kpi(Hc0:Hc3),'-');
xlabel('Field (T)');
ylabel('K_i (GHz)');
title('Dissipation rates');
legend('External dissipation rate','Internal dissipation rate')

figure
plot(Hx(Hc0:Hc3),medfilt1(kpi(Hc0:Hc3)./kpe(Hc0:Hc3),Options.order),'-');
ylabel('K_i/K_e');
xlabel('Field (T)');
title('Ratio of Dissipation rates');

% Plot the resonant frequency from Lorentzian fit versus DC magnetic field
figure
cmap2 = copyobj(cmap0, gca);
set(cmap2, 'edgeColor','none');
shading interp;
colorbar
caxis([max([mag2db(min(mag0)), -30]) max(mag2db(max(mag0)),0)+1]);
% caxis('auto');
hold on
% plot(Hx(Hc2:round(length(Hx)/200):Hc3),ff0(Hc2:round(length(Hx)/200):Hc3,:),'or','MarkerSize',2,'MarkerFaceColor','red');
plot(Hx(Hc0:round(length(Hx)/200):Hc3),ff0(Hc0:round(length(Hx)/200):Hc3,:),'or','MarkerSize',2,'MarkerFaceColor','red');
xlabel('Field (T)');
ylabel('Frequency (GHz)');
title(num2str(Temperature,'Resonant frequency from fitted data at T = %3.3f K'));
axis([field_l field_h freq_l freq_h]);

%copy over the colormap of the raw data as the background of the next plot
figure
cmap3 = copyobj(cmap0, gca);
set(cmap3, 'edgeColor','none');
shading interp;
colorbar
caxis([max([mag2db(min(mag0)), -30]) max(mag2db(max(mag0)),0)+1]);
% caxis('auto');
axis([field_l field_h freq_l freq_h]);

% Plot the resonant frequency from minimum search versus DC magnetic field
% figure
hold on
f0 = medfilt1(f0,order); % apply median filter to remove some noise
% hfig1 = plot(H0(1:round(length(H0)/200):end),f0(1:round(length(f0)/200):end),'ok','MarkerSize',2,'MarkerFaceColor','black');
plot(H0, f0, 'ok', 'MarkerSize', 2);
xlabel('Field (T)');
ylabel('Resonant frequency (GHz)');
title(num2str(Temperature,'Resonant frequency from minimum search at T = %3.3f K'));
axis([field_l field_h freq_l freq_h]);

% Plot the peak amplitude from minimum search vs. magnetic field
figure
mag0 = medfilt1(mag0, order); % apply median filter to remove some noise
plot(xq(1,:), mag2db(mag0), 'o', 'MarkerSize', 2);
xlabel('Field(T)');
ylabel('S11 (dB)');
title(num2str(Temperature,'Minimal S11 at T = %3.3f K'));

%Plot Quality factor from minimum search vs magnetic field
figure
H0 = H0(Q0 >=0); %Remove unphysical points
Q0 = Q0(Q0 >=0); %Remove unphysical points
Q0 = medfilt1(Q0, order);
Qplot1 = plot(H0(1:round(length(H0)/100):end), Q0(1:round(length(Q0)/100):end),'s-','MarkerSize',2);
% plot(H0, Q0,'s-','MarkerSize',2);
% Plot Quality factor from Lorentzian fit vs magnetic field
hold on
Hx = Hx(Qf >=0);
Qf = Qf(Qf >=0); %Remove unphysical points
Qf = medfilt1(Qf, order);
Qplot2 = plot(Hx(Hc2:round(length(Hx)/100):Hc3), Qf(Hc2:round(length(Qf)/100):Hc3),'o-','MarkerSize',2);
% plot(Hx, Qf,'o-','MarkerSize',2);
gca;
xlabel('Field (T)');
ylabel('Q factor');
legend('Quality factor from FWHM', 'Quality factor from fitting');
title(num2str(Temperature,'Quality factor, T= %.3f'));

if Options.fitfunc == 1
    figure
%     plot(Hx,Gc,'-o')
    plot(Hx(Hc2:Hc3),medfilt1(Gc(Hc2:Hc3),order),'-o')
    xlabel('Field (T)')
    ylabel('Coupling strength (GHz)')
    title(num2str(Brf,'Fitting paramters from fitting, B_0 = %.2f'))
    hold on
    yyaxis right
%     plot(Hx,gamma,'-s')
    plot(Hx(Hc2:Hc3), medfilt1(1./gamma(Hc2:Hc3),order),'-s')
    ylabel('Spin lifetime (ns)')
%     legend(['Gc, \tau_{avg} = ', num2str(1/mean(gamma),'%.2e')])
    legend('Gc','\tau')
%     legend(['\tau_0 = ', num2str(1/mean(gamma),'%.2e')]);
%     figure
%     plot(Hx, B0, '-o')
%     xlabel('Field (T)')
%     ylabel('Line cross field (T)')
end
end
% End of option 2
function [xq,yq,zq,analysis] = option3(LoadObj, Options)
%% Set data range and parameters
clear freq S11 dB N FdB FrS FiS FTT1 FTT2

% extract data from raw data file
if ~exist('out','var')
    out = readdata_v4(LoadObj, Options.nData);
end

freq = out.data.ZVLfreq/1e9;
S11 = out.data.ZVLreal + 1i*out.data.ZVLimag;  
H = out.data.DCField1;
T1 = out.data.Temperature1;
% T2 = out.data.Temperature2;

Temperature = unique(T1(H == min(H))); %Extract the measurement temperature at the lowest field (to avoid magnetoresistance in the thermometer)
HH = repmat(H,1,size(freq,2)); %populate the magnetic field to match the dimension of S11 matrix

freq = freq';
freq = freq(:);
S11 = S11';
S11 = S11(:);
HH = HH';
HH = HH(:);

[rows,~,freq] = find(freq); % Remove nonsensical zeros from the frequency data
HH = HH(rows);
S11 = S11(rows);

freq_l = min(freq); %set frequency range, l: lower limit, h: higher limit
freq_h = max(freq);
field_l = min(H);  % set field range, l: lower limit, h: lower limit
field_h = max(H); 
% Interpolate the data on a 2D grid for the colormap

%Plot the temperature vs magnetic field to check the temperature variation
figure
[H,idx,~] = unique(H,'stable');
T1 = T1(idx);
plot(H(1:round(length(T1)/100):end),T1(1:round(length(T1)/100):end),'o-')
xlabel('DC Magnetic field')
ylabel('Temperature')
title('Magnetic field vs Temperature')
%% Code For R&S ZVL-6
% Clean up the raw data by removing incomplete scans (step 1) and duplicates(step 2)
trunc1 = find(freq == min(freq),1,'first'); 
trunc2 = find(freq == min(freq),1,'last')-1; 
freq_temp = freq(trunc1:trunc2);
HH_temp = HH(trunc1:trunc2);
S11_temp = S11(trunc1:trunc2);

% Truncate according to the upper and lower field limit
S11_temp = S11_temp(HH_temp>=field_l & HH_temp<=field_h);
freq_temp = freq_temp(HH_temp>=field_l & HH_temp<=field_h);
HH_temp = HH_temp(HH_temp>=field_l & HH_temp<=field_h);

% Truncate according to the upper and lower frequency limit
S11_temp = S11_temp(freq_temp>=freq_l & freq_temp<=freq_h);
freq_temp = freq_temp(freq_temp>=freq_l & freq_temp<=freq_h);
HH_temp = HH_temp(freq_temp>=freq_l & freq_temp<=freq_h);

freq_l = min(freq_temp); %set frequency range, l: lower limit, h: higher limit
freq_h = max(freq_temp);
field_l = min(HH_temp);  % set field range, l: lower limit
field_h = max(HH_temp);  % set field range, h: lower limit
[xq,yq] = meshgrid(linspace(field_l,field_h,1001),linspace(freq_l,freq_h,1201));

dif = diff(freq); % frequency increments
resets = find(dif <= (freq_l-freq_h)); % Find the termination point of a complete scans
nop = round(mean(diff(resets))); % Set the number of points per frequency scan

mag_temp = abs(S11_temp);
mag = reshape(mag_temp,nop,[]);
% dB = reshape(mag2db(mag_temp),nop,[]);  
freq = reshape(freq_temp,nop,[]);
HH = reshape(HH_temp,nop,[]);

analysis.temp = interp1(H,T1,xq(1,:));
%% Temperary code for Keysight PNA-x N5242B
% S11_temp = S11;
% dB_temp = mag2db(abs(S11_temp));
% dB_temp = dB_temp - max(dB_temp,[],'all'); % Shift the data to compensate changes in RE(impedence) at low temperatures
% freq_temp = freq;
% HH_temp = HH;
% 
% dif = diff(freq); % frequency increments
% resets = find(dif<=0.9*(freq_l-freq_h)); % Find the termination point of a complete scans (10% error)
% nop = round(mean(diff(resets))); % Set the number of points per frequency scan
% dB = reshape(mag2db(abs(S11)),nop,[]);  %reshape the matrix so that each complete frequency scan occupy one column
% dB = dB - max(dB,[],'all'); % Shift the data to compensate changes in RE(impedence) at low temperatures
% freq = reshape(freq,nop,[]);
% HH = reshape(HH,nop,[]);
%% Clean up raw data
clearvars dif step S11_temp mag_temp freq_temp HH_temp

if Options.bgdmode == 1% Construct the background noise by stitching together zero-field-scan and anti-crossing
    %Find the line crossing by tracing the resonant peaks
    H0 = zeros(size(mag,2),1);
    f0 = zeros(size(mag,2),1);
    mag0 = zeros(size(mag,2),1);
    %find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
%     [~,bidx] = min( HH(1,:) ); % find the column index for zero-field frequency scan
%     [~,bidx] = max( HH(1,:) ); % find the column index for high-field frequency scan
    [bidx,~,~] = find(mag == min(mag( HH(1,:)==min(HH(1,:)) ),mag( HH(1,:) == max(HH(1,:)) )) ); % pick the finner peak as the first part of the background
    for ii = 1:size(mag,2) %Searching column minima (fixed field)
        [~,idx] = min( mag(:,ii) );
        if(length(idx)>1)
            disp(num2str(H0(ii),'multiple minima found at H = %.3f'))
        end
        f0(ii) = freq(idx,ii);
        H0(ii) = HH(idx,ii);
        mag0(ii) = mag(idx,ii);
    end
    [~,cidx] = min(mag(:,bidx)); % find the resonant peak to be replaced
    [~,Hpos] = max(medfilt1(abs(f0-f0(bidx)))); % scan with the resonant peak deviates from f0 the furtherest
    mag_dif = abs(mag(:,bidx)-mag(:,Hpos));
    [lidx,~] = find(mag_dif(1:cidx) <= 3e-2, 1, 'last'); % left boundary of the resonant peak to be replaced
    [ridx,~] = find(mag_dif(cidx:end) <= 3e-2, 1, 'first'); % right boundary of the resonant peak to be replaced
    if isempty(lidx); lidx = 1; end
    if isempty(ridx); ridx = length(mag_dif);else; ridx = cidx + ridx -1; end
    
    bgd0 = mag(:,bidx); % zero-field frequency scan
    bgd0(lidx:ridx) = mag(lidx:ridx,Hpos);
    bf0 = freq(:,bidx);
    [bf0,trimIdx] = unique(bf0);
    bgd0 = bgd0(trimIdx);
    figure
    plot(freq(:,bidx),mag(:,bidx))
    hold on
    plot(freq(:,Hpos),mag(:,Hpos))
    plot(bf0,bgd0);
    xlabel('Frequency (GHz)')
    ylabel('S11')
    legend('B = 0',num2str(H0(Hpos),'B = %.2f T'),'Stitched')
    disp('Press any key to continue, or Ctrl + C to abort.')
    pause
elseif Options.bgdmode == 2% load background data from existing file.
    %     backpath = ('/Volumes/GoogleDrive/My Drive/File sharing/PhD program/Research projects/LiHoF4 project/Data/Experiment/Cavity resonator/D24mm_T5mm_G_flex/SuperCoax calibration');
    backpath = ('G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data/Experiment\Cavity resonator\D24mm_T5mm_G_flex\SuperCoax calibration');
    backfile = ('background.mat');
    load(fullfile(backpath,backfile),'bf0','bgd0');
    figure
    plot(bf0,bgd0);
    xlim([freq_l freq_h])
    xlabel('Frequency (GHz)')
    ylabel('S11')
    legend('Background')
    [bf0,trimIdx] = unique(bf0);
    bgd0 = bgd0(trimIdx);
    bgd0 = interp1(bf0,bgd0,freq(:,1));
    bf0 = freq(:,1);
    hold on
    plot(bf0,bgd0);
    legend('Loaded data','Interpolated data');
    disp('Press any key to continue')
    pause
else
    disp('No background normalization')
    bgd0 = zeros(size(freq(:,1),1),1);
end

% Normalize the data to background noise
if Options.bgdmode ~= 0
    mag = mag(trimIdx,:); % remove duplicate points
    freq = freq(trimIdx,:);
    HH = HH(trimIdx,:);
    bgdM = repmat(bgd0,1,size(mag,2));
    mag = mag./bgdM;
end
FaS = scatteredInterpolant(HH(:),freq(:),mag(:));
zq = FaS(xq,yq);

% Find all the resonant peaks
f0 = zeros(size(zq,2),1);
Q0 = zeros(size(zq,2),1);
% mag0 = zeros(size(zq,2),1);
FWHM = zeros(size(zq,2),1);
%find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
for ii = 1:size(zq,2) %Searching column minima (fixed field)
    [~,idx] = min( zq(:,ii) );
    if(length(idx)>1)
        disp(num2str(H0(ii),'multiple minima found at H = %.3f'))
    end
    H0(ii) = xq(idx,ii); 
    f0(ii) = yq(idx,ii);
%     mag0(ii) = zq(idx,ii);
    HM = ( max(zq(:,ii)) + min(zq(:,ii)) )/2; % 
    % Calculate quality factor using f0/FWHM
    if isnan(1/range(yq(zq(:,ii) <= HM)))
       Q0(ii) = 0;
    elseif isempty(range(yq(zq(:,ii) <= HM)))
       Q0(ii) = 0;
    else
       FWHM(ii) = range(yq(zq(:,ii) <= HM));
       Q0(ii) = yq(idx,ii)/FWHM(ii);
    end
end

% Cleaning data
Q0(isinf(Q0)) = NaN; % Cut out inf from the array
[Q0,c] = rmmissing(Q0); % Cut out NaN from the array
H0 = H0(~c); % Remove corresponding elements in H0 array as well
f0 = f0(~c);
% mag0 = mag0(~c);
FWHM = FWHM(~c);
xq = xq(:,~c);
yq = yq(:,~c);
zq = zq(:,~c);
analysis.temp = analysis.temp(~c);

% For noisy data, we need to remove duplicates of minima
[H0,ia,~] = unique(H0,'stable');
f0 = f0(ia);
Q0 = Q0(ia);
% mag0 = mag0(ia);
FWHM = FWHM(ia);
xq = xq(:,ia);
yq = yq(:,ia);
zq = zq(:,ia);
analysis.temp = analysis.temp(ia);

f0 = medfilt1(f0,Options.order); % apply median filter to remove some noise
f0 = f0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
ff0 = f0(2); % Initial guess of the cavity resonant frequency
H0 = H0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
Q0 = Q0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
% mag0 = mag0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
FWHM = FWHM(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
xq = xq(:,f0 >= freq_l & f0 <= freq_h);
yq = yq(:,f0 >= freq_l & f0 <= freq_h);
zq = zq(:,f0 >= freq_l & f0 <= freq_h);
analysis.temp = analysis.temp(f0 >= freq_l & f0 <= freq_h);
clearvars c idx ia ii HM trunc1 trunc2 dupl nop trimIdx

weight = double.empty(size(zq,1),size(zq,2),0);
for jj = 1:size(zq,2)
    for ii = 1:size(zq,1)
        weight(ii,jj,1) = abs((zq(ii,jj)-max(zq(:,jj)))./zq(ii,jj));
    end
    weight(:,jj) = weight(:,jj)./max(weight(:,jj));
end
switch Options.fitfunc % Pick fitting function from either (1) custom function of (2) spec1d
    case 1 %Option 1: Custom function by Input-output formalism
        % Step 1: Fit the data for the first time to extract "kpe" and "w0" far from the level crossing
        kpe = double.empty(length(H0),0);
        kpi = double.empty(length(H0),0);
        w0 = double.empty(length(H0),0);
        Qf = double.empty(length(H0),0);
        parfor ii = 1:length(H0)
% %         for ii = 1:length(H0)
            plt = false;
%             figWin = Hpos-100:20:Hpos+100; % The iteration window over which shows fitting graph
%             if ismember(ii,figWin) % Not useable in parallel mode (parfor)
%                 plt = true;
%             end
            % Fit using input-output formalism
            param =  [1e-4   1e-4   1e-3    5e-4   ff0]; % starting value for param = {'kpe', 'kpi', 'Re[chi]', 'Im[chi]', 'w0'}
            bound_l = [0      0      0      5e-4   ff0]; % lower bound of fitting parameters
            bound_h = [Inf   Inf    Inf     Inf    ff0]; % upper bound of fitting parameters
            fit = iptoptx(yq(:,ii),zq(:,ii),H0(ii),param,bound_l,bound_h,weight(:,ii),plt);
            
            if mod(ii,20) == 0
                worker = getCurrentTask();
                fprintf('Fitting, current field: %1$3.2f T. Core %2$u.\n', H0(ii), worker.ID);
            end
          
            param = coeffvalues(fit);
            [idx,~,~] = find(fit(yq(:,ii)) == min(fit(yq(:,ii)))); % Find the resonant frequency by (unique) minimum search
            w0(ii) = yq(idx(1)); 
            kpe(ii) = param(1);
            kpi(ii) = param(2);
            Qf(ii) = w0(ii)/(param(1)+param(2));
            xi(ii) = param(3);
            xr(ii) = param(4);
        end
        analysis.kpe = kpe;
        analysis.kpi = kpi;
        analysis.w0 = w0;
        analysis.xr = xr;
        analysis.xi = xi;        
    case 2 %Option 2: use spec1d package to fit the data using Lorentzian form.
        parfor ii = 1:length(H0)
            s = spec1d(yq(:,ii), -zq(:,ii), max(-zq(:,ii)))*0.001; % create spec1d object
            %starting point for the (Lorentzian) fitting parameters
            p = [0.1 ff0(ii) FWHM(ii) min(zq(:,ii))]; % (p1: scaling factor ,p2: resonant frequency; p3: FWHM; p4:noise floor(?) )
            fix = [0 0 0 0]; % Denoting if the fitting parameters are fixed
            [~, fbck] = fits(s, 'lorz', p, fix);
            %             [~, fbck] = fits(s, 'lorz');
            ff0(ii) = fbck.pvals(2); % Retrieve the resonant frequency from fitted data
            Qf(ii) = abs(fbck.pvals(2)/fbck.pvals(3)/2); %Calculate the quality factor
            %             chi(ii) = 1/Qf(ii);
            if mod(ii,20) == 0
                worker = getCurrentTask();
                fprintf('Current magnetic field: %1$3.2f. on core %2$u.\n', H0(ii), worker.ID);
            end
        end
end

% Plot the interpolated frequency response data in a field scan using color map
figure
cmap = pcolor(xq,yq,mag2db(zq));
set(cmap, 'edgeColor','none')
shading interp;
colorbar
set(gca,'fontsize',Options.ftsz)
xlabel('Field (T)');
ylabel('Frequency (GHz)');
xticks(linspace(field_l,field_h,6));
title(num2str(Temperature,'S11 response at T = %3.3f K'));
legend('Experimental data','Minimum search');

figure
cmap1 = copyobj(cmap,gca);
hold on
plot(H0,f0,'.k','MarkerSize',6);
set(cmap1, 'edgeColor','none')
shading interp;
colorbar
set(gca,'fontsize',Options.ftsz)
xlabel('Field (T)');
ylabel('Frequency (GHz)');
xticks(linspace(field_l,field_h,6));
ylim([freq_l,freq_h]);
title(num2str(Temperature,'S11 response at T = %3.3f K'));
legend('Experimental data','Minimum search');

% Plot the interpolated frequency response data in a field scan using color map
figure
cmap2 = copyobj(cmap,gca);
hold on
plot(H0,w0,'.r','MarkerSize',6);
set(cmap2, 'edgeColor','none')
shading interp;
colorbar
set(gca,'fontsize',Options.ftsz)
xlabel('Field (T)');
ylabel('Frequency (GHz)');
xticks(linspace(field_l,field_h,6));
ylim([freq_l,freq_h]);
title(num2str(Temperature,'S11 response at T = %3.3f K'));
legend('Experimental data','|S11| fit');

figure
plot(H0,f0,'ok','MarkerSize',4);
hold on
plot(mean(H0(f0==min(medfilt1(f0,Options.order)))), mean(f0(f0==min(f0))), 'o', 'MarkerFaceColor', 'red','MarkerSize',4);
xlabel('Field (T)');
ylabel('Resonant frequency (GHz)');
title(num2str(Temperature,'Resonant frequency from minimum search at T = %3.3f K'));

figure
plot(H0,medfilt1(1./Q0,Options.order),'.k');
hold on
[~,idx,Qf] = find(Qf);
plot(xq(1,idx),medfilt1(1./Qf,Options.order),'.r');
xlabel('Field (T)');
ylabel('1/Q');
title(num2str(Temperature,'Inverse of Q factor at T = %3.3f K'));
legend('f_0 / FWHM','|S11| fitting');

if Options.fitfunc == 1
    figure
    plot(H0,medfilt1(xr,Options.order),'-');
    ylabel('Susceptibility');
    xlabel('Field (T)');
    title('Re[\chi]')
    
    figure
    plot(H0,medfilt1(xi,Options.order),'-');
    ylabel('Susceptibility');
    xlabel('Field (T)');
    title('Im[\chi]')
    
    figure
    plot(H0,medfilt1(kpe,Options.order),'o');
    ylabel('K_e (GHz)');
    hold on
    yyaxis right
    plot(H0,medfilt1(kpi,Options.order),'*');
    xlabel('Field (T)');
    ylabel('K_i (GHz)');
    title('Dissipation rates');
    legend('External dissipation rate','Internal dissipation rate')
end
end
% End of option 3
function [xq,yq,zq,analysis] = option4(LoadObj, Options)
% extract data from raw data file
out = readdata_v4(LoadObj, Options.nData);
freq = out.data.ZVLfreq/1e9;
S11 = out.data.ZVLreal + 1i*out.data.ZVLimag;  
H = out.data.DCField1;
T1 = out.data.Temperature1;
% T2 = out.data.Temperature2;

Field = mean(H);
analysis.field = Field;
TT = repmat(T1,1,size(freq,2)); %populate the magnetic field to match the dimension of S11 matrix

freq = freq';
freq = freq(:);
S11 = S11';
S11 = S11(:);
dB = mag2db(abs(S11));
TT = TT';
TT = TT(:);

[rows,~,freq] = find(freq); % Remove nonsensical zeros from the frequency data
TT = TT(rows);
S11 = S11(rows);
dB = dB(rows);

freq_l = min(freq); %set frequency range, l: lower limit, h: higher limit
freq_h = max(freq);
% freq_h = 2.503; % Manually set the upper limit of frequency scan when ZVL fails
temp_l = min(T1);  % set field range, l: lower limit, h: higher limit
temp_h = max(T1);

% Could use "scatteredInterpolant()" to replace "TriScatteredInterp()" as recommended by MATLAB, but it may generate artifacts
FdB = TriScatteredInterp(TT,freq,dB);
% FrS = TriScatteredInterp(TT,freq,real(S11));
% FiS = TriScatteredInterp(TT,freq,imag(S11)); %intrapolate points on a 2D grid

% Plot frequency-field colour map
[xq,yq] = meshgrid(linspace(temp_l,temp_h,301),linspace(freq_l,freq_h,310)); %set the X and Y range
zq = FdB(xq,yq);

%Clean up the raw data by removing incomplete scans (step 1) and duplicates(step 2)
dif = nonzeros(diff(freq)); % Extract from raw data the step size in frequency scan
dif = dif(dif>0); % Keep only positive steps
step = min (dif); % Calculate the value of the frequency scan step
nop = ceil(abs(freq_h-freq_l)/step); %compute how many points pers complete frequency scan.

clearvars dif out step
% step 1: truncate the beginning and end part to keep only complete frequency scans and reshape the matrices into single column vectors
trunc1 = find(freq==freq_l,1,'first'); 
trunc2 = find(freq==freq_l,1,'last')-1; 
S11_temp = S11(trunc1:trunc2);
dB_temp = mag2db(abs(S11_temp));
freq_temp = freq(trunc1:trunc2);
TT_temp = TT(trunc1:trunc2);

% Step 2: remove duplicates
dupl = find(diff(freq_temp) == 0);
freq_temp(dupl+1)=[];
dB_temp(dupl+1)=[];
TT_temp(dupl+1)=[];

S11 = reshape(S11_temp,nop,[]);
dB = reshape(dB_temp,nop,[]);  %reshape the matrix so that each complete frequency scan occupy one column
freq = reshape(freq_temp,nop,[]);
TT = reshape(TT_temp,nop,[]);

%Find all the resonant peaks
f0 = zeros(size(dB,2),1);
T0 = zeros(size(dB,2),1);

%find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
for ii = 1:size(dB,2) %Searching column minima (fixed field)
    [~,idx] = min( dB(:,ii) );
    if(length(idx)>1)
        disp(num2str(T0(ii),'multiple minima found at H = %.3f T'))
    end
    f0(ii) = freq(idx,ii);
    T0(ii) = TT(idx,ii); 
end

% For noisy data, we need to remove duplicates of minima
[T0,ia,~] = unique(T0,'stable');
f0 = f0(ia);

f0 = medfilt1(f0); % apply median filter to remove some noise
f0 = f0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
T0 = T0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints

clearvars c idx ia ii HM trunc1 trunc2 dupl nop

% Plot the interpolated frequency response data in a field scan using color map
figure
cmap = pcolor(xq,yq,zq);
hold on
plot(T0,f0,'.r','MarkerSize',6);
set(cmap, 'edgeColor','none')
shading interp;
colorbar
set(gca,'fontsize',Options.ftsz)
xlabel('Field (T)');
ylabel('Frequency (GHz)');
xticks(linspace(temp_l,temp_h,6));
title(num2str(Field,'S11 response at H = %3.3f T'));

figure
plot(T0,f0,'ok','MarkerSize',4);
hold on
plot(mean(T0(f0==min(f0))), mean(f0(f0==min(f0))), 'o', 'MarkerFaceColor', 'red','MarkerSize',4);
xlabel('Temperature (K)');
ylabel('Resonant frequency (GHz)');
title(num2str(Field,'Resonant frequency from minimum search at H = %3.3f T'));
end
% End of option 4