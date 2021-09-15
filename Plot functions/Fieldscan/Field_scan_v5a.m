format long;
clearvars -except mfields freq S11 analysis
% For macOS
% addpath('/Volumes/GoogleDrive/My Drive/File sharing/Programming scripts/Matlab/Plot functions/Fieldscan/functions')
% addpath(genpath('/Volumes/GoogleDrive/My Drive/File sharing/Programming scripts/Matlab/Plot functions/spec1d--Henrik'));

% For windows
addpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\Fieldscan\functions\');
addpath(genpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\spec1d--Henrik\'));

% Process options:
Options.analysis = 2; % Analysis options (1.simple plot, 2.on-resonance fit, 3.off-resonance fit, 4.temp scan)
Options.plot = false; % Plot option in analysis part (Option 3)
% Options.save = 'n'; % Option to save the analysis
Options.lnwd = 1.5; % plot linewidth
Options.ftsz = 12; % plot font size
Options.mksz = 3; % plot marker size
Options.bgdmode = 2; % Background normalization (0: no normalization. 1: normalized with stitched background. 2: normalize with loaded file)
Options.nData = 1; % Number of dataset from VNA
Options.dType = 'sim'; % Input data type: 1. Experimental ('exp'), 2. Simulated ('sim')
Options.phPlot = false; % Option to plot phase in color plots
Options.fitfunc = 1; % Fitting function (only for analysis-3): (1) Input-out or (2) Lorentzian

% Determin file location base on data type
if strcmp(Options.dType, 'exp')
    location = ['G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\',...
        'SC127\SC127_1 (2.5 x1 x 0.5 mm, rectangle)\17.08.2019'];
    % location = ['/Volumes/GoogleDrive/My Drive/File sharing/PhD program/Research projects/LiHoF4 project/',...
        % 'Data/Experiment/LiHoF4/SC200/2021.04.12'];
    loadname = '2019_08_0006';
    LoadObj = fullfile(location,[loadname, '.dat']);
    SaveObj = fullfile(location,[loadname, '_interp', '.mat']);
elseif strcmp(Options.dType, 'sim')
    location = ['G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Simulations\',...
        'Matlab\Susceptibilities\S11 parameters'];
    loadname = sprintf('sim_3.64GHz_RPA_S11_0.180K_33deg');
    LoadObj = fullfile(location,[loadname, '.mat']);
    
    Options.save = 'n'; % Option to save the analysis
    Options.bgdmode = 0; % no background renormalization for simulated data
else
    fprintf(['Unknow data type:',Options.dType])
    return
end

if exist('analysis','var') && exist('freq','var') && exist('S11','var') && exist('mfields','var')
    [~,fname,~] = fileparts(LoadObj);
    if isfield(analysis, 'name') && strcmp(fname, analysis.name) % check if target data is in the workspace
        arguments =  {analysis mfields freq S11 LoadObj Options};
    else
        arguments = {LoadObj Options};
    end
else
    arguments = {LoadObj Options};
end

switch Options.analysis
    case 1
        % Simple color plot of S11 (& S21)
        [mfields,freq,S11,analysis] = option1(arguments{:});
    case 2
        % On-resonance measurement processing (w/ plots)
        [mfields,freq,S11,analysis] = option2(arguments{:});
    case 3
        % Off-resonance measurement processing (w/ plots)
        [mfields,freq,S11,analysis] = option3(arguments{:});
    case 4
        % Temperature scan
        [mfields,freq,S11,analysis] = option4(LoadObj, Options);
end
clearvars -except mfields freq S11 analysis Options SaveObj

if ~isfield(Options,'save')
    prompt = sprintf('Save the analysis?\n');
    answer = input(prompt,'s');
else
    answer = Options.save;
end

switch lower(answer)
    case {'y','yes'}
        fprintf('Saving analysis...\n')
        if isfield(analysis,'bgd0')
            background.f = analysis.bf0;
            background.d = analysis.bgd0;
            analysis = rmfield(analysis,'bf0');
            analysis = rmfield(analysis,'bgd0');
            save(SaveObj,'analysis','background','mfields','freq','S11','-v7.3')
        else
            save(SaveObj,'analysis','mfields','freq','S11','-v7.3')
        end
        %         if Options.analysis == 4 % Add a point to the phase diagram from off-resonance measurement
        %             fprintf('Updating phase diagram\n')
        %             phase(1,1) = Temperature;
        %             phase(1,2) = mean(H0(f0==min(f0)));
        %             if Options.savefile == true && isfile(dataobj)
        %                 save(dataobj,'phase','-append');
        %             elseif Options.savefile == true
        %                 save(dataobj,'phase','-v7.3');
        %             end
        %         end
    otherwise
        fprintf('Aborted saving analysis\n')
end
clearvars prompt SaveObj

function [xq,yq,zq,analysis] = option1(varargin)
clearvars -except varargin

[HH,freq,mag,LoadObj,Options,analysis] = load_data(varargin{:});
Temperature = analysis.temp;

% field_l = 1.6;
% field_h = 5.5;
freq_l = min(freq(:,1));
freq_h = max(freq(:,1));

if exist('field_l','var') && exist('field_h','var')
    % Truncate the data to field limits
    [~,lidx] = find(HH(1,:) >= field_l,1,'first');
    [~,uidx] = find(HH(1,:) <= field_h,1,'last');
    mag = mag(:,lidx:uidx);
    freq = freq(:,lidx:uidx);
    HH = HH(:,lidx:uidx);
end

field_l = min(HH(1,:));  %set field range, l: lower limit, h: higher limit
field_h = max(HH(1,:));
[xq,yq] = meshgrid(linspace(field_l,field_h,1501),linspace(freq_l,freq_h,1201)); %set the X and Y range
%% Background cleanup
clearvars dif step lidx uidx T1 *_temp
if Options.bgdmode ~= 0
    %Find all the resonant peaks
    f0 = zeros(size(mag,2),1);
    H0 = zeros(size(mag,2),1);
    mag0 = zeros(size(mag,2),1);
%     find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
    mag0l = min(mag(:,HH(1,:)==min(HH(1,:)) )); % peak depth at the lowest field
    mag0h = min(mag(:,HH(1,:) == max(HH(1,:)) )); % peak depth at the highest field
    [~,hl] = min([mag0l mag0h]); % pick the finner peak as the first part of the background
    if hl == 1
        [~,bidx] = min( HH(1,:) );
    else
        [~,bidx] = max( HH(1,:) );
    end
    for ii = 1:size(mag,2) %Searching column minima (fixed field)
        [~,idx] = min( mag(:,ii) );
        if(length(idx)>1)
            fprintf(num2str(H0(ii),'multiple minima found at H = %.3f T\n'))
        end
        f0(ii) = freq(idx,ii);
        H0(ii) = HH(idx,ii);
        mag0(ii) = mag(idx,ii);
    end
end

while true
    if Options.bgdmode == 1 % Construct the background noise by stitching together zero-field-scan and anti-crossing
        [~,cidx] = min(mag(:,bidx)); % find the resonant peak to be replaced
        [~,Hpos] = max(mag0); % Method 1. scan of max reflection
%     [~,Hpos] = max(mag(cidx,:)); % Method 2. scan of max reflection at f0 (frequency)
%     [~,Hpos] = max(medfilt1(abs(f0-f0(bidx)))); % Method 3. scan with the resonant peak deviates from f0 the furtherest
%     [~,Hpos] = min(abs(H0-6.258)); % Method 4. Manually set the patch slice of the frequency scan
        mag_dif = abs(mag(:,bidx)-mag(:,Hpos));
        [lidx,~] = find(mag_dif(1:cidx) <= 1e-3, 1, 'last'); % left boundary of the resonant peak to be replaced
        [ridx,~] = find(mag_dif(cidx:end) <= 1e-3, 1, 'first'); % right boundary of the resonant peak to be replaced
        if isempty(lidx); lidx = 1; end
        if isempty(ridx); ridx = length(mag_dif); else; ridx = cidx + ridx -1; end
        analysis.Hpos = Hpos; % store field location of the anti-crossing
        analysis.cidx = cidx; % peak center
        bf0 = freq(:,bidx);
        bgd0 = mag(:,bidx); % Base frequency scan
        bgd0(lidx:ridx) = mag(lidx:ridx,Hpos);
        
        figure
        plot(bf0,mag(:,bidx))
        hold on
        plot(freq(:,Hpos),mag(:,Hpos))
        plot(bf0,bgd0);
        xlabel('Frequency (GHz)')
        ylabel('S11')
        legend(num2str(H0(bidx),'B = %.2f T'),num2str(H0(Hpos),'B = %.2f T'),'Stitched','location','southeast')
        [bf0,trimIdx] = unique(bf0);
        bgd0 = bgd0(trimIdx);
        break
    elseif Options.bgdmode == 2 % load background data from existing file
        fprintf('Loading background normalization data...\n')
        [backpath,backfile,~] = fileparts(LoadObj);
        load(fullfile(backpath,[backfile,'_interp.mat']),'background');
        if exist('background','var')
            bf0 = background.f;
            bgd0 = background.d;
            bgd0 = interp1(bf0,bgd0,freq(:,1));
            bf0 = freq(:,1);
            figure
            plot(bf0,bgd0);
            hold on
            plot(freq(:,bidx),mag(:,bidx))
            xlim([freq_l freq_h])
            xlabel('Frequency (GHz)')
            ylabel('S11')
            legend('Background','Raw data')
            [~,trimIdx] = unique(bf0);
            bgd0 = bgd0(trimIdx); % Shift the loaded background to the noise floor of current data
            break
        else
            fprintf('Background data not found! Reverting to extraction from data.\n')
            Options.bgdmode = 1;
        end
    else
        fprintf('No background normalization\n')
        break
    end
end

% Renormailize the background
if Options.bgdmode ~= 0
    
    mag = mag(trimIdx,:); % Remove duplicate points
    freq = freq(trimIdx,:);
    HH = HH(trimIdx,:);
    bgd0 = medfilt1(bgd0);
    bgdM = repmat(bgd0,1,size(mag,2));
    mag = mag./bgdM;
    if Options.bgdmode == 1
        analysis.bgd0 = bgd0;
        analysis.bf0 = bf0;
    end
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
if Options.phPlot == true
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
if exist('lck1','var')
    figure
    plot(H(1:100:end),lck1(1:100:end),'s-')
    xlabel('DC Magnetic field')
    ylabel('Hall resistence')
    title('Hall resistance vs Field')
end

if exist('lck2','var')
    figure
    plot(H(1:100:end),lck2(1:100:end),'x-')
    xlabel('DC Magnetic field')
    ylabel('Sample resistence')
    title('Sample resistance vs Field')
end

% Plot the temperature profile against magnetic field
if isfield(analysis,'Ts') && isfield(analysis,'Hs')
    figure
    plot(analysis.Hs,analysis.Ts,'o-')
    xlabel('DC Magnetic field')
    ylabel('Temperature')
    title('Magnetic field vs Temperature')
end
end
% End of option 1 (simple color plot)
function [xq,yq,zq,analysis] = option2(varargin)
%Set data range and parameters
clearvars -except varargin
[HH,freq,mag,LoadObj,Options,analysis] = load_data(varargin{:});
Temperature = analysis.temp;

% set desirable field range
field_l = 2.0;
field_h = 5.0;
field_cut = field_l + 0.5; % Field window for cavity parameter fit
freq_l = min(freq(:,1));
freq_h = max(freq(:,1));
strnth = 'strong';

%Plot the temperature vs magnetic field to check the temperature variation
if isfield(analysis,'Ts') && isfield(analysis,'Hs')
    figure
    plot(analysis.Hs,analysis.Ts,'o-')
    xlabel('DC Magnetic field')
    ylabel('Temperature')
    title('Magnetic field vs Temperature')
end

% Truncate the data to field limits
[~,lidx] = find(HH(1,:) >= field_l,1,'first');
[~,uidx] = find(HH(1,:) <= field_h,1,'last');

mag = mag(:,lidx:uidx);
freq = freq(:,lidx:uidx);
HH = HH(:,lidx:uidx);
% set field range according to actual limits
field_l = min(HH(1,:));
field_h = max(HH(1,:));

% Interpolate the data on a 2D grid for the colormap
[xq,yq] = meshgrid(linspace(field_l,field_h,601),linspace(freq_l,freq_h,2001));
Hx = xq(1,:);
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
clearvars dif step *_temp lidx uidx T1

%Find all the resonant peaks
f0 = zeros(size(mag,2),1);
H0 = zeros(size(mag,2),1);
mag0 = zeros(size(mag,2),1);
%find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
mag0l = min(mag(:,HH(1,:)==min(HH(1,:)) )); % peak depth at the lowest field
mag0h = min(mag(:,HH(1,:) == max(HH(1,:)) )); % peak depth at the highest field
[~,hl] = min([mag0l mag0h]); % pick the finner peak as the first part of the background
if hl == 1
    [~,bidx] = min( HH(1,:) );
else
    [~,bidx] = max( HH(1,:) );
end
for ii = 1:size(mag,2) %Searching column minima (fixed field)
    [~,idx] = min( mag(:,ii) );
    if(length(idx)>1)
        fprintf(num2str(H0(ii),'multiple minima found at H = %.3f\n'))
    end
    f0(ii) = freq(idx,ii);
    H0(ii) = HH(idx,ii);
    mag0(ii) = mag(idx,ii);
end 
[~,cidx] = min(mag(:,bidx)); % find the resonant peak to be replaced
[~,Hpos] = max(mag0); % Method 1. scan of max reflection
%     [~,Hpos] = max(medfilt1(mag(cidx,:))); % Method 2. scan of max reflection at f0 (frequency)
%     [~,Hpos] = max(medfilt1(abs(f0-f0(bidx)))); % Method 3. scan with the resonant peak deviates from f0 the furtherest
%     [~,Hpos] = min(abs(H0-3.41)); % Method 4. Manually set the patch slice of the frequency scan

while true
    if Options.bgdmode == 1
        fprintf('Normalizing background...\n')
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
        legend('B = 0',num2str(H0(Hpos),'B = %.2f T'),'Stitched','location','SouthEast')
        fprintf('Checkpoint: press any key to continue, or Ctrl + c to abort.\n')
        pause
        break
    elseif Options.bgdmode == 2 % load background data from existing file.
        fprintf('Loading background normalization data...\n')
        [backpath,backfile,~] = fileparts(LoadObj);
        if isfile(fullfile(backpath,[backfile,'_interp.mat']))
            load(fullfile(backpath,[backfile,'_interp.mat']),'background');
            if exist('background','var')
                bf0 = background.f;
                bgd0 = background.d;
                figure
                plot(bf0,bgd0);
                xlim([freq_l freq_h])
                xlabel('Frequency (GHz)')
                ylabel('S11')
                legend('Background')
                bgd0 = interp1(bf0,bgd0,freq(:,1));
                bf0 = freq(:,1);
                plot(freq(:,HH(1,:)==min(HH(1,:)) ),mag(:,HH(1,:)==min(HH(1,:)) ));
                hold on
                plot(bf0,bgd0);
                legend('Zero field raw data','Loaded background', 'Location', 'SouthEast');
                fprintf('Checkpoint: press any key to continue, or Ctrl + c to abort.\n')
                pause
                break
            else
                fprintf('Background data not found! Reverting to extraction from data.\n')
                Options.bgdmode = 1;
            end
        else
            fprintf('Background data not found! Reverting to extraction from data.\n')
            Options.bgdmode = 1;
        end
    else
        fprintf('No background normalization.\n')
        bf0 = freq(:,1);
        bgd0 = ones(size(freq(:,1),1),1);
        break
    end
end

% Remove duplicate points
[~,trimIdx] = unique(bf0);
bgd0 = medfilt1(bgd0(trimIdx));
mag = medfilt1(mag(trimIdx,:));
freq = freq(trimIdx,:);
HH = HH(trimIdx,:);
if Options.bgdmode == 1
    analysis.bgd0 = bgd0;
    analysis.bf0 = bf0;
end

% Renormailize the background
if Options.bgdmode ~= 0
    bgdM = repmat(bgd0,1,size(mag,2));
    mag = mag./bgdM;
end
clearvars c idx ia lidx ridx widx ii HM trunc1 trunc2 dupl nop bf0

%Find all the resonant peaks
f0 = zeros(size(mag,2),1);
FWHM0 = zeros(size(mag,2),1);
%find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
for ii = 1:size(mag,2) %Searching column minima (fixed field)
    [~,idx] = min( medfilt1(mag(:,ii)) );
    if(length(idx)>1)
        fprintf(num2str(H0(ii),'multiple minima found at H = %.3f'))
    end
    f0(ii) = freq(idx,ii);
    H0(ii) = HH(idx,ii); 
    mag0(ii) = mag(idx,ii);
    HM = ( max(mag(:,ii)) + min(mag(:,ii)) )/2; % 
    left = mag(1:idx,ii);
    right = mag(idx:end,ii);
    [lidx,~] = find(left >= HM, 1, 'last');
    [ridx,~] = find(right >= HM, 1, 'first');
    if lidx == idx+ridx-1
        FWHM0(ii) = Inf;
    elseif isempty(lidx) || isempty(ridx)
        FWHM0(ii) = Inf;
    else
        FWHM0(ii) = freq(idx+ridx-1,ii)-freq(lidx,ii);
    end
end
Q0 = f0./FWHM0;
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

% Discard nonsensical datapoints
f0 = f0(f0 >= freq_l & f0 <= freq_h);
H0 = H0(f0 >= freq_l & f0 <= freq_h);
Q0 = Q0(f0 >= freq_l & f0 <= freq_h);
mag0 = mag0(f0 >= freq_l & f0 <= freq_h);
FWHM0 = FWHM0(f0 >= freq_l & f0 <= freq_h); 

% Reduce data noise
H0 = H0(Q0 >=0); 
Q0 = Q0(Q0 >=0);
f0 = medfilt1(f0);
Q0 = medfilt1(Q0);

% Remove outlier data
[f0,rmidx] = rmoutliers(f0,'mean');
Q0 = Q0(~rmidx);
H0 = H0(~rmidx);
mag0 = mag0(~rmidx);
FWHM0 = FWHM0(~rmidx);

% Plot Resonant frequency trace with fitting parameters
Tplot = figure;
plot(H0,f0,'ok','MarkerSize',Options.mksz);
% plot(H0(1:round(length(H0)/200):end),f0(1:round(length(f0)/200):end),'ok','MarkerSize',Options.mksz);
hold on

tick = 0;
slope = 1;
while true
    switch strnth
        case 'weak'
            % Fit the field dependent resonant frequency data with weak coupling function
            init = [0.1   0.1   slope   f0(1)  H0(Hpos)]; % [g  gamma slope wc x0]
            bound_l = [0   0  -Inf   freq_l  field_l];
            bound_h = [1   1   Inf   freq_h  field_h];
            wt = abs(f0-f0(1))/max(abs(f0-f0(1)));
            yh = 3.6401;
            yl = 3.6357;
            xh = 7.4456;
            xl = 5.0657;
            gradiant = (yh-yl)/(xh-xl); % tilting angle to compensate off-resonance effect on the resonant frequencies
            % gradiant = 0;
            % [fitP,~] = wk_cpl_fit(H0,f0,init,bound_l,bound_h,wt,true);
            [fitP,~] = wk_cpl_fit(H0,f0-( gradiant*(H0-field_l) ),init,bound_l,bound_h,wt,true);
            slopes = fitP.slope;
            wc = fitP.wc;
            % wc = 3.645;
            gcs = fitP.g;
            Brfs = fitP.x0; % Level crossing location from perturbative fitting
            fprintf('Weak coupling fit complete\n')
            break
        case 'strong'
            gradiant = 0;
            [~, lidx] = min(f0);
            [~, uidx] = max(f0(lidx:end));
%             [~, lidx] = max(f0);
%             [~, uidx] = min(f0(lidx:end));
            wt = ones(size(f0));
%% single dispersion split fit-----------------------            
            init = [0.05   -slope   f0(1)  H0(Hpos)]; % [g slope wc x0]
            bound_l = [ 0   -Inf   freq_l   field_l];
            bound_h = [Inf   Inf   freq_h   field_h];
            [lofit,~] = str_cpl_lfit(H0(1:lidx),f0(1:lidx),init,bound_l,bound_h,wt(1:lidx),false);

            init = [lofit.g   lofit.slope   lofit.wc  lofit.x0]; % [g slope wc x0]
            bound_l = [ 0   -Inf   freq_l   field_l];
            bound_h = [Inf   Inf   freq_h   field_h];
            [upfit,~] = str_cpl_ufit(H0(lidx+uidx-1:end),f0(lidx+uidx-1:end),init,bound_l,bound_h,wt(lidx+uidx-1:end),false);

            wc = mean([lofit.wc upfit.wc]);
            Brfs = mean([lofit.x0 upfit.x0]);
            slopes = mean([lofit.slope upfit.slope]);
            Hx1 = linspace(field_l,H0(lidx+uidx-1),501);
            Hx2 = linspace(H0(lidx),field_h,length(Hx1));
            lofit_ex = interp1(H0(1:lidx),lofit(H0(1:lidx)),Hx1,'spline','extrap');
            upfit_ex = interp1(H0(lidx+uidx-1:end),upfit(H0(lidx+uidx-1:end)),Hx2,'spline','extrap');
            analysis.lofit = lofit_ex;
            analysis.upfit = upfit_ex;
%% single dispersion optimization fit-----------------------------------
%             fitting parameter = ['wc', 'slope', 'x0', 'g']
            init = [wc   -slopes   Brfs  0.05];
            bound_l = [freq_l   -Inf   field_l   0];
            bound_h = [freq_h    Inf   field_h   1];
            l_branch = @(param, x) (param(1) + param(2)*(x-param(3))./2 - sqrt((param(2)*(x-param(3))).^2 + 4*param(4)^2)./2);
            u_branch = @(param, x) (param(1) + param(2)*(x-param(3))./2 + sqrt((param(2)*(x-param(3))).^2 + 4*param(4)^2)./2);
            branches = @(param, x) [l_branch(param, x(:,1)), u_branch(param, x(:,2))];

            fitP = lsqcurvefit(branches, init, [Hx1', Hx2'], [lofit_ex', upfit_ex'], bound_l, bound_h); % fit extrapolated data
%             fitP = lsqcurvefit(branches, init, [H0(lidx-300:lidx), H0(lidx+uidx-1:end)],...
%                 [f0(lidx-300:lidx), f0(lidx+uidx-1:end)], bound_l, bound_h); % fit the original data
            
            wc = fitP(1);
            slopes = fitP(2);
            Brfs = fitP(3);
            gcs = fitP(4);
            [~,Hpos] = min(abs(Hx-Brfs));
            omg1 = slopes.*(Hx-Brfs) + wc;
%% double dispersion split fit--------------------------------------
%             init = [0.05   0.01   -slope  -slope   mean(f0)  H0(Hpos)  H0(Hpos)]; % [g slope wc x0]
%             bound_l = [0   0   -Inf   -Inf    freq_l   H0(lidx)   H0(lidx)];
%             bound_h = [1   1    Inf    Inf    freq_h   H0(lidx+uidx)   H0(lidx+uidx)];
%             [upfit,~] = str_cpl_ufit2(H0(lidx+uidx-1:end),f0(lidx+uidx-1:end),init,bound_l,bound_h,wt(lidx+uidx-1:end),false);
%             
%             init = [upfit.g   upfit.g1   upfit.slope  upfit.slope1   upfit.wc  upfit.x0  upfit.x1]; % [g slope wc x0]
%             bound_l = [0   0   -Inf   -Inf    freq_l   H0(lidx)   H0(lidx)];
%             bound_h = [1   1    Inf    Inf    freq_h   H0(lidx+uidx)   H0(lidx+uidx)];
%             [lofit,~] = str_cpl_lfit2(H0(1:lidx),f0(1:lidx),init,bound_l,bound_h,wt(1:lidx),false);
%             
%             wc = mean([upfit.wc lofit.wc]);
%             slopes = [lofit.slope lofit.slope1];
%                         
%             init = [mean([lofit.g upfit.g])   mean([lofit.g1 upfit.g1])   mean(slopes)  mean(slopes)   wc  upfit.x0  upfit.x1]; % [g slope wc x0]
%             bound_l = [0    0    min([lofit.slope upfit.slope])   min([lofit.slope1 upfit.slope1])    freq_l   lofit.x0  lofit.x1];
%             bound_h = [1    1    max([lofit.slope upfit.slope])   max([lofit.slope1 upfit.slope1])    freq_h   lofit.x0  lofit.x1];
%             [upfit,~] = str_cpl_ufit2(H0(lidx+uidx-1:end),f0(lidx+uidx-1:end),init,bound_l,bound_h,wt(lidx+uidx-1:end),false);
%             
%             init = [mean([lofit.g upfit.g])   mean([lofit.g1 upfit.g1])   mean([lofit.slope upfit.slope])  mean([lofit.slope1 upfit.slope1])   wc  upfit.x0  upfit.x1]; % [g slope wc x0]
%             bound_l = [0    0    min([lofit.slope upfit.slope])   min([lofit.slope1 upfit.slope1])    mean([lofit.wc upfit.wc])   lofit.x0  lofit.x1];
%             bound_h = [1    1    max([lofit.slope upfit.slope])   max([lofit.slope1 upfit.slope1])    mean([lofit.wc upfit.wc])   lofit.x0  lofit.x1];
%             [lofit,~] = str_cpl_lfit2(H0(1:lidx),f0(1:lidx),init,bound_l,bound_h,wt(1:lidx),false);
%             
% %             lofit_ex = interp1(H0(1:lidx),lofit(H0(1:lidx)),Hx,'spline','extrap');
% %             upfit_ex = interp1(H0(lidx+uidx-1:end),upfit(H0(lidx+uidx-1:end)),Hx,'spline','extrap');
            
% %             wc = mean([upfit.wc lofit.wc]);
% %             gcs = [upfit.g upfit.g1];
% %             Brfs = [upfit.x0 upfit.x1];
% %             slopes = [upfit.slope upfit.slope1];
% %             [~,Hpos1] = min(abs(Hx-min(upfit.x0,upfit.x1)));
% %             [~,Hpos2] = min(abs(Hx-max(upfit.x0,upfit.x1)));
% %             omg1 = upfit.slope.*(Hx-upfit.x0) + upfit.wc;
% %             omg2 = upfit.slope1.*(Hx-upfit.x1) + upfit.wc;
% %             fprintf('Strong coupling fit complete\n')
            
%             wc = lofit.wc;
%             gcs = [lofit.g lofit.g1];
%             Brfs = [lofit.x0 lofit.x1];
%             slopes = [lofit.slope lofit.slope1];
%             [~,Hpos1] = min(abs(Hx-min(lofit.x0,lofit.x1)));
% %             [~,Hpos2] = min(abs(Hx-max(lofit.x0,lofit.x1)));
%             omg1 = lofit.slope.*(Hx-lofit.x0) + lofit.wc;
%             omg2 = lofit.slope1.*(Hx-lofit.x1) + lofit.wc;
%             fprintf('Strong coupling fit complete\n')
%             analysis.lofit = lofit;
%             analysis.upfit = upfit;
%% double dispersion piecewise fit-----------------------------------
% %             fitting parameter : wc, slope, x0, g, slope1, x1, g1
%             init = [mean(f0)   -slope   H0(Hpos)  0.05   -slope   H0(Hpos)  0.05];
%             bound_l = [freq_l   -Inf   mean(f0)   0   -Inf   field_l   0];
%             bound_h = [freq_h    Inf   mean(f0)   1    Inf   field_h   1];
%             [pcfit,~] = str_cpl_pcfit2(H0,f0,H0(lidx),H0(lidx),init,bound_l,bound_h,wt,false);
% %             wt = ones(size(Hx));
% %             [pcfit,~] = str_cpl_pcfit2(Hx, (Hx <= H0(lidx+uidx-1)).*lofit_ex + (Hx >= H0(lidx)).*upfit_ex,...
% %                 H0(lidx),H0(lidx),init,bound_l,bound_h,wt,true);
%             analysis.pcfit = pcfit;
% 
%             gcs = [pcfit.g  pcfit.g1];
%             slopes = [pcfit.slope   pcfit.slope1];
%             Brfs = [pcfit.x0    pcfit.x1];
%             wc = pcfit.wc;
%             [~,Hpos1] = min(abs(Hx-pcfit.x0));
% %             [~,Hpos2] = min(abs(Hx-pcfit.x1));
%             omg1 = pcfit.slope.*(Hx-pcfit.x0) + pcfit.wc;
%             omg2 = pcfit.slope1.*(Hx-pcfit.x1) + pcfit.wc;
%% -----------------------------------
            figure(Tplot)
            if exist('lofit_ex','var')
                lf = plot(Hx1, lofit_ex,'-r');
                set(lf,'LineWidth',Options.lnwd);
                if exist('upfit_ex','var')
                    uf = plot(Hx2, upfit_ex,'-r');
                    set(uf,'LineWidth',Options.lnwd);
                end
            end
            xlabel('Field (T)');
            ylabel('Resonant frequency (GHz)');
            title(num2str(Temperature,'Resonant frequency from minimum search at T = %3.3f K'));
            axis([field_l field_h freq_l freq_h]);
            break
        otherwise
            prompt = sprintf('Choose coupling strength regime: "strong" or "weak" \n');
            strnth = input(prompt,'s');
            tick = tick + 1;
            if tick >= 5
                clearvars tick
                break
            end
    end
end
% Interpolate data along sampling mesh
FWHM0 = interp1(H0,FWHM0,Hx,'pchip','extrap'); % Using spline interpolation to smooth the FWHM
mag0 = interp1(H0,mag0,Hx,'pchip','extrap'); % amplitude along resonant frequency trace
Q0 = interp1(H0,Q0,Hx,'pchip','extrap');
analysis.mag0 = mag0;
analysis.FWHM0 = FWHM0;
analysis.Q0 = Q0;

% Fit Quality factor curve as a function of field
init = [1e-3  mean(gcs)  1e-3   mean(slopes) + gradiant  Hx(Hpos)];  % [gamma gc kappa slope x0]
bound_l = [0   0   0   -Inf   field_l];
bound_h = [1   1   1    Inf   field_h];
[Qfit0, ~] = Qf_fit(Hx(Hpos-50:Hpos+50), Q0(Hpos-50:Hpos+50), init, bound_l, bound_h, false);
% [Qfit0, ~] = Qf_fit(Hx, Q0, init, bound_l, bound_h, false);
analysis.Qfit0 = Qfit0;% gc = Qfit0.gc;
omg0 = Qfit0.slope.*(Hx-Qfit0.x0) + wc; % linearized dispersion relation
% slope = Qfit0.slope;
% Brfs = Qfit0.x0;
% [~,Hpos] = min(abs(Hx-Qfit0.x0));

figure(Tplot)
hold on
plot(Hx,omg1,'--k',Hx,omg0,'-.r');
if exist('omg2','var')
    plot(Hx,omg2,'--b');
    if exist('omg0','var')
        plot(Hx,omg0,'-.r');
        lgd = {'Exp. data','Theo. data','',sprintf('gc = %.3f GHz', gcs(1)), sprintf('gc = %.3f GHz',gcs(2)),...
            sprintf('Q_0 fit, gc = %.3f GHz', Qfit0.gc)};
    else
        lgd = {'Exp. data','Theo. data','',sprintf('gc = %.3f GHz', gcs(1)), sprintf('gc = %.3f GHz',gcs(2))};
    end
else
    lgd = {'Exp. data','Theo. data','',sprintf('f_0 fit, gc = %.3f GHz', mean(gcs))};
end
legend(lgd,'location', 'southwest')

%Plot Quality factor from minimum search vs magnetic field
Qplot = figure;
plot(Hx, Q0,'s-','MarkerSize',Options.mksz);
% plot(Hx(1:round(length(H0)/20):end), Q0(1:round(length(Q0)/20):end),'s-','MarkerSize',Options.mksz);
hold on
Qp0 = plot(Qfit0,'-k');
set(Qp0,'LineWidth',Options.lnwd);

% Plot frequency scan at line crossing vs. zero-field scan
figure
plot(freq(:,HH(1,:)==min(HH(1,:)) ),mag(:,HH(1,:)==min(HH(1,:)) ),'-');
hold on
[~,Hp] = min(abs(HH(1,:)-mean(Brfs)));
% plot(freq(1:20:end,Hpos),mag(1:20:end,Hpos),'.');
plot(freq(:,Hp),mag(:,Hp),'-');
xlabel('Frequency (GHz)');
ylabel('S11');
legend(sprintf('Frequency cut at %.2f T',min(HH(1,:)) ),sprintf('Frequency cut at %.2f T',HH(1,Hp)), 'Location', 'SouthEast');
clearvars B Delt hPara fPara wp wm fitPara H_res f_res Hp Qparam

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
        Fmag = scatteredInterpolant( HH(:), freq(:), mag(:) ); % 2D interpolation of amplitude data
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
zq(isnan(zq))=0;
clearvars init bound_l bound_h *_temp
        
set(gca,'fontsize',Options.ftsz)
axis([field_l field_h freq_l freq_h]);
xlabel('Field (T)');
ylabel('Frequency (GHz)');
title(num2str(Temperature,'S11 response at T = %3.3f K'));

weight = double.empty(size(zq,1),size(zq,2),0); % Weight function option 1
for jj = 1:length(Hx)
    for ii = 1:length(zq(:,jj))
        weight(ii,jj,1) = abs(1./zq(ii,jj));
    end
end
% weight = 1./gradientweight(zq); % Weight function option 2
% weight = ones(size(zq));  % Weight function option 3
%% Step 1: Fit the data for the first time to extract "kpe" and "w0" far from the level crossing
ff0 = double.empty(length(Hx),2,0); 
Qf = double.empty(length(Hx),0);
Gc1 = double.empty(length(Hx),0);
Gc2 = double.empty(length(Hx),0);
kpe = double.empty(length(Hx),0);
kpi = double.empty(length(Hx),0);
gamma = double.empty(length(Hx),0);

Hc0 = find(Hx >= field_l, 1,'first');
Hc1 = find(Hx <= field_cut, 1, 'last');
for ii = Hc0:Hc1 % fit the largely detunned data
%     % Selectively plot the fit for visual inspection
    plt = false;
%     figWin = [1:round((Hc1-Hc0)/10):Hc1]; % The iteration window over which shows fitting graph
%     if ismember(ii,figWin)
%         plt = true;
%     end
if exist('omg2','var')
    % starting value for param = {'kpe', 'kpi', 'w0', 'Gc1', 'Gc2', 'gma','attn'}
    param = [FWHM0(ii)  FWHM0(ii)  wc  gcs(1)  gcs(2)  mean(gcs)/10  1/mean(zq(:,ii))];
    bound_l = [0   0   freq_l   gcs(1)   gcs(2)   mean(gcs)/10   0.5]; % lower bound of fitting parameters
    bound_h = [1   1   freq_h   gcs(1)   gcs(2)   mean(gcs)/10   1.5]; % upper bound of fitting parameters
    fit = iptopt2(yq(:,ii), zq(:,ii), Hx(ii), omg1(ii), omg2(ii), param, bound_l, bound_h, weight(:,ii), plt);
    Gc1(ii,1) = fit.Gc1;
    Gc2(ii,1) = fit.Gc2;
else
    %   starting value for param = {'kpe', 'kpi', 'omega', 'Gc', 'gma', 'wc', 'attn'}
    param = [FWHM0(ii)  FWHM0(ii)  omg1(ii)  gcs  gcs/10  wc  1/mean(zq(:,ii))];
    bound_l = [0   0   omg1(ii)    gcs    gcs/10   freq_l    0.5]; % lower bound of fitting parameters
    bound_h = [1   1   omg1(ii)    gcs    gcs/10   freq_h    1.5]; % upper bound of fitting parameters
    fit = iptopt(yq(:,ii), zq(:,ii), Hx(ii), param, bound_l, bound_h, weight(:,ii), plt);
    Gc1(ii,1) = fit.Gc;
end
    kpe(ii,1) = fit.kpe;
    kpi(ii,1) = fit.kpi;
    gamma(ii,1) = fit.gma;

    [~,idx] = findpeaks(-fit(yq(:,ii)));
    if ~isempty(idx)
        if length(idx) == 2
            ff0(ii,1,1) = yq(idx(1),ii);  % Find the resonant frequency by minimum search
            ff0(ii,2,1) = yq(idx(2),ii);
        elseif length(idx) > 2
            [~,idxt] = min(zq(idx,ii));
            ff0(ii,:,1) = yq(idx(idxt),ii);
        else
            ff0(ii,:,1) = yq(idx,ii);
        end
    else
        ff0(ii,:,1) = fit.wc;
    end  
    Qf(ii) = mean(ff0(ii,:,1)./(fit.kpe + fit.kpi));
end
%% Step 2: fit the data for the second time with fixed "kpe" and "w0"
wc = mean(rmoutliers(ff0(Hc0:Hc1,1),'median'));
kpe0 = mean(rmoutliers(kpe(Hc0:Hc1),'median'));
kpi0 = mean(rmoutliers(kpi(Hc0:Hc1),'median'));
Hc2 = find(Hx <= field_h, 1, 'last');
for ii = Hc1:Hc2
    % Selectively plot the fit for visual inspection
    plt = false;
    figWin = Hpos-10:2:Hpos+10; % The iteration window over which shows fitting graph
    if ismember(ii,figWin) % Not useable in parallel mode (parfor)
        plt = true;
    end

if exist('omg2','var')
    % starting value for param = {'kpe', 'kpi', 'w0', 'Gc1', 'Gc2', 'gma','attn'}
    param = [kpe0  kpi0  wc  Gc1(ii-1)  Gc2(ii-1)  gamma(ii-1)  1/mean(zq(:,ii))];
    bound_l = [kpe0   kpi0   wc   0   0   0   0.5]; % lower bound of fitting parameters
    bound_h = [kpe0   kpi0   wc   1   1   1   1.5]; % upper bound of fitting parameters
    fit = iptopt2(yq(:,ii), zq(:,ii), Hx(ii), omg1(ii), omg2(ii), param, bound_l, bound_h, weight(:,ii), plt);
    Gc1(ii) = fit.Gc1;
    Gc2(ii) = fit.Gc2;
else
    % starting value for param = {'kpe', 'kpi', 'omega', 'Gc', 'gma', 'wc', 'attn'}
    param  =  [kpe0  kpi0   omg1(ii)   Gc1(ii-1)  gamma(ii-1)  wc  1/mean(zq(:,ii))];
    bound_l = [kpe0  kpi0   omg1(ii)    0     0     wc   0.5]; % lower bound of fitting parameters
    bound_h = [kpe0  kpi0   omg1(ii)    1     1     wc   1.5]; % upper bound of fitting parameters
    fit = iptopt(yq(:,ii), zq(:,ii), Hx(ii), param, bound_l, bound_h, weight(:,ii), plt);             
    Gc1(ii) = fit.Gc;
    omg1(ii) = fit.omega;
end
    kpe(ii,1) = fit.kpe;
    kpi(ii,1) = fit.kpi;
    gamma(ii,1) = fit.gma;
    
    [~,idx] = findpeaks(-fit(yq(:,ii)));
    if ~isempty(idx)
        if length(idx) == 2
            ff0(ii,1) = yq(idx(1),ii);  % Find the resonant frequency by minimum search
            ff0(ii,2) = yq(idx(2),ii);
        elseif length(idx) > 2
            [~,idxt] = min(zq(idx,ii));
            ff0(ii,:) = yq(idx(idxt),ii);
        else
            ff0(ii,:) = yq(idx,ii);
        end
    else
        ff0(ii,1) = yq(fit(yq(:,ii)) == min(fit(yq(:,ii))));  % Find the resonant frequency by minimum search
    end
    Qf(ii) = mean(ff0(ii,:)./(param(1)+param(2))); % Calculate quality factor
end
% 
% % Fit Qf data
% init = [mean(gamma)  mean(gcs)  mean(FWHM0)  mean(slopes)  mean(Brfs)];  % [gamma gc kappa slope x0]
% bound_l = [0  0  0  -Inf   0 ];
% bound_h = [1  1  1   Inf  Inf];
% [Qfit, ~] = Qf_fit(Hx, Qf, init, bound_l, bound_h, false);

analysis.kpe = kpe;
analysis.kpi = kpi;
analysis.wc = ff0;
analysis.w0 = omg1;
analysis.Gc1 = Gc1;
if exist('omg2','var'); analysis.Gc2 = Gc2; analysis.w1 = omg2; end
analysis.gamma = gamma;
analysis.Qf = Qf;
% analysis.Qfit = Qfit;

figure
plot(Hx(Hc0:Hc2),kpe(Hc0:Hc2),'o','MarkerSize',Options.mksz);
ylabel('K_e (GHz)');
hold on
yyaxis right
plot(Hx(Hc0:Hc2),kpi(Hc0:Hc2),'o','MarkerSize',Options.mksz);
xlabel('Field (T)');
ylabel('K_i (GHz)');
title('Dissipation rates');
legend('External dissipation rate','Internal dissipation rate')

figure
plot(Hx(Hc0:Hc2),medfilt1(kpi(Hc0:Hc2)./kpe(Hc0:Hc2)),'-');
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
% plot(Hx(Hc2:round(length(Hx)/200):Hc2),ff0(Hc2:round(length(Hx)/200):Hc2,:),'or','MarkerSize',Options.mksz,'MarkerFaceColor','red');
plot(Hx(Hc0:round(length(Hx)/200):Hc2),ff0(Hc0:round(length(Hx)/200):Hc2,:),'or','MarkerSize',Options.mksz,'MarkerFaceColor','red');
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
hold on
f0 = medfilt1(f0); % apply median filter to remove some noise
% hfig1 = plot(H0(1:round(length(H0)/200):end),f0(1:round(length(f0)/200):end),'ok','MarkerSize',Options.mksz,'MarkerFaceColor','black');
plot(H0, f0, 'ok', 'MarkerSize', Options.mksz); % Plot the resonant frequency from minimum search
xlabel('Field (T)');
ylabel('Resonant frequency (GHz)');
title(num2str(Temperature,'Resonant frequency from minimum search at T = %3.3f K'));
% caxis('auto');
axis([field_l field_h freq_l freq_h]);

% Plot the peak amplitude from minimum search vs. magnetic field
figure
mag0 = medfilt1(mag0); % apply median filter to remove some noise
plot(xq(1,:), mag2db(mag0), 'o', 'MarkerSize', Options.mksz);
xlabel('Field(T)');
ylabel('S11 (dB)');
title(num2str(Temperature,'Minimal S11 at T = %3.3f K'));

% Plot Quality factor from fitting function vs magnetic field
% figure(Qplot)
% hold on
% Hx = Hx(Qf >=0);
% Qf = Qf(Qf >=0); %Remove unphysical points
% Qf = medfilt1(Qf, order);
% plot(Hx, Qf,'ro-','MarkerSize',Options.mksz);
% % plot(Hx(Hc2:round(length(Hx)/20):Hc2), Qf(Hc2:round(length(Qf)/20):Hc2),'o-','MarkerSize',Options.mksz);
% plot(Hx, Qf,'o-','MarkerSize',Options.mksz);
% Qp = plot(Qfit,'-k');
% set(Qp,'LineWidth',Options.lnwd);
gca;
xlabel('Field (T)');
ylabel('Q factor');
legend('Quality factor from FWHM', 'Q0 fit', 'Quality factor from fitting', 'Qf fit');
title(num2str(Temperature,'Quality factor, T= %.2f mK'));

figure
% plot(Hx,Gc,'-o')
plot(Hx(Hc1:Hc2),medfilt1(Gc1(Hc1:Hc2)),'-o')
hold on
if exist('omg2','var')
    plot(Hx(Hc1:Hc2),medfilt1(Gc2(Hc1:Hc2)),'-o')
    lg = {'Gc1','Gc2'};
else
    lg = {'Gc'};
end
xlabel('Field (T)')
ylabel('Coupling strength (GHz)')
title(num2str(Brfs(1),'Fitting paramters from fitting, B_0 = %.2f'))
hold on
yyaxis right
% plot(Hx,gamma,'-s')
plot(Hx(Hc1:Hc2), medfilt1(gamma(Hc1:Hc2)),'-s')
lg{end+1} = '\gamma';
ylabel('Spin linewidth (GHz)')
legend(lg)
% figure
% plot(Hx, B0, '-o')
% xlabel('Field (T)')
% ylabel('Line cross field (T)')
end
% End of option 2 (on-resonance fitting)
function [xq,yq,zq,analysis] = option3(varargin)
%% Set data range and parameters
clear freq S11 dB N FdB FrS FiS FTT1 FTT2

[HH,freq,mag,~,Options,analysis] = load_data(varargin{:});
Temperature = analysis.temp;

% set desirable field range
% field_l = 2.0;
% field_h = 5.0;
freq_l = min(freq(:,1));
freq_h = max(freq(:,1));

if exist('field_l','var') && exist('field_h','var')
    % Truncate the data to field limits
    [~,lidx] = find(HH(1,:) >= field_l,1,'first');
    [~,uidx] = find(HH(1,:) <= field_h,1,'last');
    mag = mag(:,lidx:uidx);
    freq = freq(:,lidx:uidx);
    HH = HH(:,lidx:uidx);
end

%Plot the temperature vs magnetic field to check the temperature variation
if isfield(analysis,'Ts') && isfield(analysis,'Hs')
    figure
    plot(analysis.Hs,analysis.Ts,'o-')
    xlabel('DC Magnetic field')
    ylabel('Temperature')
    title('Magnetic field vs Temperature')
end
%% Clean up raw data
clearvars dif step *_temp

if Options.bgdmode == 1% Construct the background noise by stitching together zero-field-scan and anti-crossing
    %Find the line crossing by tracing the resonant peaks
    H0 = zeros(size(mag,2),1);
    f0 = zeros(size(mag,2),1);
    mag0 = zeros(size(mag,2),1);
    %find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
    mag0l = min(mag(:,HH(1,:)==min(HH(1,:)) )); % peak depth at the lowest field
    mag0h = min(mag(:,HH(1,:) == max(HH(1,:)) )); % peak depth at the highest field
    [~,hl] = min([mag0l mag0h]); % pick the finner peak as the first part of the background
    if hl == 1
        [~,bidx] = min( HH(1,:) );
    else
        [~,bidx] = max( HH(1,:) );
    end
    for ii = 1:size(mag,2) %Searching column minima (fixed field)
        [~,idx] = min( mag(:,ii) );
        if(length(idx)>1)
            fprintf(num2str(H0(ii),'multiple minima found at H = %.3f\n'))
        end
        f0(ii) = freq(idx,ii);
        H0(ii) = HH(idx,ii);
        mag0(ii) = mag(idx,ii);
    end
    [~,cidx] = min(mag(:,bidx)); % find the resonant peak to be replaced
    [~,Hpos] = max(medfilt1(abs(f0-f0(bidx)))); % scan with the resonant peak deviates from f0 the furtherest
    mag_dif = abs(mag(:,bidx)-mag(:,Hpos));
    [lidx,~] = find(mag_dif(1:cidx) <= 1e-3, 1, 'last'); % left boundary of the resonant peak to be replaced
    [ridx,~] = find(mag_dif(cidx:end) <= 1e-3, 1, 'first'); % right boundary of the resonant peak to be replaced
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
    fprintf('Press any key to continue, or Ctrl + C to abort.')
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
    fprintf('Press any key to continue\n')
    pause
else
    fprintf('No background normalization\n')
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
        fprintf(num2str(H0(ii),'multiple minima found at H = %.3f'))
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

f0 = medfilt1(f0); % apply median filter to remove some noise
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
plot(H0,f0,'.k','MarkerSize',Options.mksz);
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
plot(H0,w0,'.r','MarkerSize',Options.mksz);
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
plot(H0,f0,'ok','MarkerSize',Options.mksz);
hold on
plot(mean(H0(f0==min(medfilt1(f0)))), mean(f0(f0==min(f0))), 'o', 'MarkerFaceColor', 'red','MarkerSize',Options.mksz);
xlabel('Field (T)');
ylabel('Resonant frequency (GHz)');
title(num2str(Temperature,'Resonant frequency from minimum search at T = %3.3f K'));

figure
plot(H0,medfilt1(1./Q0),'.k');
hold on
[~,idx,Qf] = find(Qf);
plot(xq(1,idx),medfilt1(1./Qf),'.r');
xlabel('Field (T)');
ylabel('1/Q');
title(num2str(Temperature,'Inverse of Q factor at T = %3.3f K'));
legend('f_0 / FWHM','|S11| fitting');

if Options.fitfunc == 1
    figure
    plot(H0,medfilt1(xr),'-');
    ylabel('Susceptibility');
    xlabel('Field (T)');
    title('Re[\chi]')
    
    figure
    plot(H0,medfilt1(xi),'-');
    ylabel('Susceptibility');
    xlabel('Field (T)');
    title('Im[\chi]')
    
    figure
    plot(H0,medfilt1(kpe),'o');
    ylabel('K_e (GHz)');
    hold on
    yyaxis right
    plot(H0,medfilt1(kpi),'*');
    xlabel('Field (T)');
    ylabel('K_i (GHz)');
    title('Dissipation rates');
    legend('External dissipation rate','Internal dissipation rate')
end
end
% End of option 3 (off-resonance fitting)
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
        fprintf(num2str(T0(ii),'multiple minima found at H = %.3f T\n'))
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
plot(T0,f0,'.r','MarkerSize',Options.mksz);
set(cmap, 'edgeColor','none')
shading interp;
colorbar
set(gca,'fontsize',Options.ftsz)
xlabel('Field (T)');
ylabel('Frequency (GHz)');
xticks(linspace(temp_l,temp_h,6));
title(num2str(Field,'S11 response at H = %3.3f T'));

figure
plot(T0,f0,'ok','MarkerSize',Options.mksz);
hold on
plot(mean(T0(f0==min(f0))), mean(f0(f0==min(f0))), 'o', 'MarkerFaceColor', 'red','MarkerSize',Options.mksz);
xlabel('Temperature (K)');
ylabel('Resonant frequency (GHz)');
title(num2str(Field,'Resonant frequency from minimum search at H = %3.3f T'));
end
% End of option 4 (temperature-scan measurement)
function [HH,freq,mag,LoadObj,Options,analysis] = load_data(varargin)
clearvars -except varargin
if nargin > 2
    fprintf('Data already exist in workspace, skip loading...\n');
    analysis = varargin{1};
    HH = varargin{2};
    freq = varargin{3};
    mag = abs(varargin{4});
    LoadObj = varargin{5};
    Options = varargin{6};
    Options.bgdmode = 0;
    [~,analysis.name,~] = fileparts(LoadObj);
else
    LoadObj = varargin{1};
    Options = varargin{2};
    if Options.dType == 'exp'
        out = readdata_v4(LoadObj,Options.nData);
        freq = out.data.ZVLfreq/1e9;
        S11 = out.data.ZVLreal + 1i*out.data.ZVLimag;
        H = out.data.DCField1;
        HH = repmat(H,1,size(freq,2)); %populate the magnetic field to match the dimension of S11 matrix
        T1 = out.data.Temperature1; % Sample cell temperature
        [~,analysis.name,~] = fileparts(LoadObj);
        analysis.power = mean(out.data.Power1);
%        T2 = out.data.Temperature2; % Mixing chamber temperature
%        lck1 = out.data.lockin1;
%        lck2 = out.data.lockin2;
        
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
        
        %% Code For R&S ZVL/ZNL units
        % Clean up the raw data by removing incomplete scans (step 1) and duplicates(step 2)
        % step 1: truncate the beginning and end part to keep only complete frequency scans and reshape the matrices into single column vectors
        trunc1 = find(freq == min(freq),1,'first');
        trunc2 = find(freq == min(freq),1,'last')-1;
        freq_temp = freq(trunc1:trunc2);
        HH_temp = HH(trunc1:trunc2);
        S11_temp = S11(trunc1:trunc2);
        
        % Truncate according to the upper and lower frequency limit
        S11_temp = S11_temp(freq_temp>=freq_l & freq_temp<=freq_h);
        freq_temp = freq_temp(freq_temp>=freq_l & freq_temp<=freq_h);
        HH_temp = HH_temp(freq_temp>=freq_l & freq_temp<=freq_h);
        mag_temp = abs(S11_temp);
        
        freq_l = min(freq_temp); %set frequency range, l: lower limit, h: higher limit
        freq_h = max(freq_temp);
        
        dif = diff(freq_temp); % frequency increments
        resets = find( dif <= (freq_l-freq_h) ); % Find the termination point of a complete scans
        nop = round(mean(diff(resets))); % Set the number of points per frequency scan
        
        % Rearrange data to ascending field
        if HH(end) > HH(1)
            analysis.direction = 'up';
            mag = reshape(mag_temp,nop,[]);
            freq = reshape(freq_temp,nop,[]);
            HH = reshape(HH_temp,nop,[]);
        else
            analysis.direction = 'down';
            mag = flip(reshape(mag_temp,nop,[]),2);
            freq = flip(reshape(freq_temp,nop,[]),2);
            HH = flip(reshape(HH_temp,nop,[]),2);
        end
        analysis.Ts = T1;
        analysis.Hs = H;
        analysis.temp = min(analysis.temp(analysis.HH == min(analysis.HH))); % Measurement temperature at the lowest field (to avoid magnetoresistance in the thermometer)
        %% Temperary code for Keysight PNA-x N5242B
        % S11_temp = S11;
        % mag_temp = abs(S11_temp);
        % freq_temp = freq;
        % HH_temp = HH;
        % freq_l = min(freq_temp);
        % freq_h = max(freq_temp);
        % field_l = min(HH_temp);
        % % field_l = 1.6;
        % field_h = max(HH_temp);
        % [xq,yq] = meshgrid(linspace(field_l,field_h,1501),linspace(freq_l,freq_h,1201)); %set the X and Y range
        %
        % dif = diff(freq); % frequency increments
        % resets = find(dif<=0.9*(freq_l-freq_h)); % Find the termination point of a complete scans (10% error)
        % nop = round(mean(diff(resets))); % Set the number of points per frequency scan
        % mag = reshape(mag_temp,nop,[]);  %reshape the matrix so that each complete frequency scan occupy one column
        % freq = reshape(freq,nop,[]);
        % HH = reshape(HH,nop,[]);
    else
        fprintf('Loading simulated data...\n')
        load(LoadObj,'mfields','freq','S11','analysis');
        HH = mfields;
        mag = abs(S11);
    end
end
end
% End of data loading function