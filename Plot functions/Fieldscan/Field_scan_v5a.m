format long;
clearvars -except mfields freq S11 analysis

% For windows
addpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\Fieldscan\functions\');
addpath(genpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\spec1d--Henrik\'));

% Process options:
Options.analysis = 4; % Analysis options (1.simple plot, 2.on-resonance fit, 3.off-resonance fit, 4.temp scan)
Options.save = 'n'; % Option to save the analysis
Options.dType = 'exp'; % Input data type: 1. Experimental ('exp'), 2. Simulated ('sim')
Options.lnwd = 1.5; % plot linewidth
Options.ftsz = 12; % plot font size
Options.mksz = 2; % plot marker size
Options.bgdmode = 3; % Background normalization (0: no normalization. 1: Stitched background. 2: Loaded file. 3: S11 fit)
Options.nData = 1; % Number of dataset from VNA
Options.phPlot = false; % Option to plot phase in color plots
Options.fitfunc = 1; % Fitting function (only for analysis-3): (1) Input-output or (2) Lorentzian
Options.filFctr = 0.0127*2.1; % Filling factor (SC239, SC200, SC251)

% Determin file location base on data type
if strcmp(Options.dType, 'exp')
    location = ['G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\',...
        'SC251 (2.49 x 0.92 x 3.45 mm)\Lev (1x4.5x3 mm)\2021.12.18'];
    loadname = '2021_12_0057';
%     location = ['G:\My Drive\File sharing\PhD program\Research projects\Collaborations\JianRui Soh\Data\2021.12.05'];
%     loadname = '2021_12_0018';
    LoadObj = fullfile(location,[loadname,'.dat']);
    SaveObj = fullfile(location,[loadname,'_interp', '.mat']);
elseif strcmp(Options.dType, 'sim')
    location = ['G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Simulations\',...
        'Matlab\Susceptibilities\S11 parameters'];
    loadname = 'S11_LiHoF4_0.200K_3.65GHz_0.02Deg_15.5Deg_1.00e-04';
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
clearvars -except location arguments LoadObj Options SaveObj

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
        [temps,freq,S11,analysis] = option4(LoadObj, Options);
end
analysis.location = location;
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
        if Options.analysis == 4 % Add a point to the phase diagram from off-resonance measurement
            fprintf('Updating phase diagram\n')
            phase(1,1) = Temperature;
            phase(1,2) = mean(H0(f0==min(f0)));
            if Options.savefile == true && isfile(dataobj)
                save(dataobj,'phase','-append');
            elseif Options.savefile == true
                save(dataobj,'phase','-v7.3');
            end
        end
    otherwise
        fprintf('Aborted saving analysis\n')
end
clearvars prompt SaveObj

function [HH,freq,mag,analysis] = option1(varargin)
clearvars -except varargin

[HH,freq,mag,LoadObj,Options,analysis] = load_data(varargin{:});
Temperature = analysis.temp;

freq_l = min(freq(:,1));
freq_h = max(freq(:,1));
field_l = min(HH(1,:));  %set field range, l: lower limit, h: higher limit
field_h = max(HH(1,:));
% field_l = 1.6;
% field_h = 5.5;

%% Background cleanup
clearvars dif step lidx uidx T1 *_temp
if Options.bgdmode ~= 0
    %Find all the resonant peaks
    f0 = zeros(size(mag,2),1);
    H0 = zeros(size(mag,2),1);
    mag0 = zeros(size(mag,2),1);
%     find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
    [~,lidx] = min(HH(1,:));
    [~,uidx] = max(HH(1,:));
    mag0l = min(mag(:,lidx)); % peak depth at the lowest field
    mag0h = min(mag(:,uidx)); % peak depth at the highest field
    [~,hl] = min([mag0l mag0h]); % pick the finner peak as the first part of the background
    if hl == 1
        [~,bidx] = min( HH(1,:) );
        fprintf('Use low field data for background normalization...\n')
    else
        [~,bidx] = max( HH(1,:) );
        fprintf('Use high field data for background normalization...\n')
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
clearvars uidx lidx

while true
    if Options.bgdmode == 1 % Construct the background noise by stitching together zero-field-scan and anti-crossing
        [~,cidx] = min(mag(:,bidx)); % find the resonant peak to be replaced
        [~,Hpos] = max(mag0); % Method 1. scan of max reflection
%         [~,Hpos] = max(mag(cidx,:)); % Method 2. scan of max reflection at f0 (frequency)
%         [~,Hpos] = max(medfilt1(abs(f0-f0(bidx)))); % Method 3. scan with the resonant peak deviates from f0 the furtherest
%         [~,Hpos] = min(abs(H0-14.2)); % Method 4. Manually set the patch slice of the frequency scan
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
        legend(num2str(H0(bidx),'B = %.2f T'),num2str(H0(Hpos),'B = %.2f T'),'Stitched background','location','southeast')
        break
    elseif Options.bgdmode == 2 % load background data from existing file
        fprintf('Loading background normalization data...\n')
        [backpath,backfile,~] = fileparts(LoadObj);
        if isfile(fullfile(backpath,[backfile,'_interp.mat']))
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
                legend('Loaded background','Raw data')
                break
            else
                fprintf('Background data not found! Reverting to extraction from data.\n')
                Options.bgdmode = 1;
            end
        else
            fprintf('Background data not found! Reverting to extraction from data.\n')
            Options.bgdmode = 1;
        end
    elseif Options.bgdmode == 3 % Background normalization through |S11| fitting
        fprintf('Normalization by |S11| fitting...\n')
        % starting value for param = {'kpe', 'kpi', 'omega', 'Gc', 'gma', 'wc', 'attn'}
        param  =  [1e-3   1e-3   freq_h    0   1e-3   f0(bidx)  1/mean(mag(:,bidx))];
        bound_l = [ 0      0     freq_h    0   1e-3    freq_l   0.1]; % lower bound of fitting parameters
        bound_h = [Inf    Inf    freq_h    0   1e-3    freq_h   2.0]; % upper bound of fitting parameters
        fit = iptopt(freq(:,bidx), mag(:,bidx), H0(bidx), param, bound_l, bound_h, 1./mag(:,bidx), false);
        bf0 = freq(:,bidx);
        bgd0 = medfilt1(mag(:,bidx)-fit(freq(:,bidx))+1)*fit.attn;
%         bgd0 = medfilt1(mag(:,bidx)./fit(freq(:,bidx)))*fit.attn;
        figure
        plot(bf0,bgd0);
        xlim([freq_l freq_h])
        xlabel('Frequency (GHz)')
        ylabel('S11')
        legend('Background')
        plot(freq(:,bidx),mag(:,bidx));
        hold on
        plot(bf0,bgd0);
        legend('Zero field raw data','Fitted background', 'Location', 'SouthEast');
        fprintf('Checkpoint: press any key to continue, or Ctrl + c to abort.\n')
        pause
        break
    else
        fprintf('No background normalization\n')
        break
    end
end

% Renormailize the background
if Options.bgdmode ~= 0
    [bf0,trimIdx] = unique(bf0);
    bgd0 = bgd0(trimIdx); % Shift the loaded background to the noise floor of current data
    mag = mag(trimIdx,:); % Remove duplicate points
    freq = freq(trimIdx,:);
    HH = HH(trimIdx,:);
    bgd0 = medfilt1(bgd0);
    bgdM = repmat(bgd0,1,size(mag,2));
    fprintf('Extracting dissipation rates by fitting...\n')
    % starting value for param = {'kpe', 'kpi', 'omega', 'Gc', 'gma', 'wc', 'attn'}
    param  =  [1e-3   1e-3   freq_h   0   1e-3   (freq_l+freq_h)/2  1/mean(bgd0)];
    bound_l = [ 0     0    freq_h     0   1e-3    freq_l   0.8]; % lower bound of fitting parameters
    bound_h = [Inf   Inf   freq_h     0   1e-3    freq_h   1.2]; % upper bound of fitting parameters
    fit = iptopt(freq(:,bidx), mag(:,bidx)./bgd0, H0(bidx), param, bound_l, bound_h, bgd0./mag(:,bidx), false);
%     bgd0 = mag(:,bidx)./fit(freq(:,bidx));  % further normalization through fitting
    mag = mag./bgdM;
    analysis.kpe0 = fit.kpe;
    analysis.kpi0 = fit.kpi;
    analysis.wc0 = fit.wc;
    if Options.bgdmode ~= 0
        analysis.bgd0 = bgd0;
        analysis.bf0 = bf0;
    end
end
clearvars mag0l mag0h idx ia trimIdx ii bgd0 bf0

% Plot additional data if there are more than one set of data from VNA
if Options.nData > 1
    S21 = out.data.ZVLreal2 + 1i*out.data.ZVLimag2;
    S21 = S21';
    S21 = S21(:);
    mag2 = abs(S21);
    
    % Plot the interpolated frequency response data in a field scan using color map
    figure
    cmap = pcolor(HH,freq,mag2db(mag2));
    set(cmap, 'edgeColor','none')
    shading interp;
    colorbar
    set(gca,'fontsize',Options.ftsz)
    xlabel('Magnetic Field (T)');
    ylabel('Frequency (GHz)');
    xticks(linspace(field_l,field_h,6));
    title(num2str(Temperature,'S21 response at T = %3.3f K'));
end

% Plot the interpolated frequency response data in a field scan using color map
figure
cmap = pcolor(HH,freq,mag2db(mag));
set(cmap, 'edgeColor','none')
shading interp;
caxis([-20 1])
colorbar
set(gca,'fontsize',Options.ftsz)
xlabel('Magnetic Field (T)');
ylabel('Frequency (GHz)');
xticks(linspace(field_l,field_h,6));
title(num2str(Temperature,'S11 (dB) response at T = %3.3f K'));

% Plot the imaginary part of the response function in a field scan using color map
if Options.phPlot == true
    figure
    cmap1 = pcolor(HH_temp,freq_temp,imag(S11_temp));
    set(cmap1, 'edgeColor','none')
    shading interp;
    % caxis([-30 0])
    colorbar
    set(gca,'fontsize',Options.ftsz)
    xlabel('Magnetic Field (T)');
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
    plot(analysis.Hs,medfilt1(analysis.Ts),'o-')
    xlabel('DC Magnetic field')
    ylabel('Temperature')
    title('Magnetic field vs Temperature')
end
end
% End of option 1 (simple color plot)
function [HH,freq,mag,analysis] = option2(varargin)
%Set data range and parameters
clearvars -except varargin
[HH,freq,mag,LoadObj,Options,analysis] = load_data(varargin{:});
Temperature = analysis.temp;

% set desirable frequency range
% freq_l = 3.6;
freq_l = min(freq(:,1));
freq_h = max(freq(:,1));
if freq_l >= max(freq(:,1)) || freq_h <= min(freq(:,1))
    freq_l = min(freq(:,1));
    freq_h = max(freq(:,1));
    fprintf('Manual frequency range out of data range! resort back to default setting.\n')
end

% set desirable field range
field_l = 0; % manually set the field limit
field_h = 1;
if field_l >= max(HH(1,:)) || field_h <= min(HH(1,:))
    field_l = min(HH(1,:));
    field_h = max(HH(1,:));
    fprintf('Manual field range out of data range! resort back to default setting.\n')
end

strnth = 'weak'; % 'strong' or 'weak'
sgn = +1; % sign of the gradiant of the dispersion relations
slope = 0.1; % linear spin resonance slope
aux = false; % auxiliery spin resonance for fitting
    Brf2 = 13.6; % anticrossing location of auxiliary spin resonance

%Plot the temperature vs magnetic field to check the temperature variation
if isfield(analysis,'Ts') && isfield(analysis,'Hs')
    figure
    plot(analysis.Hs,medfilt1(analysis.Ts),'o-')
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

% Truncate the data to frequency limits
[lidx,~] = find(freq(:,1) >= freq_l,1,'first');
[uidx,~] = find(freq(:,1) <= freq_h,1,'last');
mag = mag(lidx:uidx,:);
freq = freq(lidx:uidx,:);
HH = HH(lidx:uidx,:);

% set field and frequency ranges according to actual limits
field_l = min(HH(1,:));
field_h = max(HH(1,:));
freq_l = min(freq(:,1));
freq_h = max(freq(:,1));
analysis.field_l = field_l;
analysis.field_h = field_h;

%% Clean up raw data
clearvars dif step *_temp lidx uidx T1

mag0l = min(mag(:,HH(1,:)==min(HH(1,:)) )); % peak depth at the lowest field
mag0h = min(mag(:,HH(1,:) == max(HH(1,:)) )); % peak depth at the highest field
[~,hl] = min([mag0l mag0h]); % pick the finner peak as the first part of the background
if hl == 1
    [~,bidx] = min( HH(1,:) );
else
    [~,bidx] = max( HH(1,:) );
end
mag = medfilt2(mag,'symmetric');
[~,fidx] = min(mag(:,bidx));
while true
    if Options.bgdmode == 1
        fprintf('Normalizing background by stiching...\n')
        [~,cidx] = min(mag(:,bidx)); % find the resonant peak to be replaced
        % [~,Hp0] = max(mag0); % Method 1. scan of max reflection
        [~,Hp0] = max(mag(cidx,:)); % Method 2. scan of max reflection at f0 (frequency)
        % [~,Hp0] = max(medfilt1(abs(f0-f0(bidx)))); %clear Method 3. scan with the resonant peak deviates from f0 the furtherest
        % [~,Hp0] = min(abs(H0-9.64)); % Method 4. Manually set the patch slice of the frequency scan
        mag_dif = abs(mag(:,bidx)-mag(:,Hp0));
        [lidx,~] = find(mag_dif(1:cidx) <= 1e-3, 1, 'last'); % left boundary of the resonant peak to be replaced
        [ridx,~] = find(mag_dif(cidx:end) <= 1e-3, 1, 'first'); % right boundary of the resonant peak to be replaced
        if isempty(lidx); lidx = 1; end
        if isempty(ridx); ridx = length(mag_dif);else; ridx = cidx + ridx -1; end
        bgd0 = mag(:,bidx); % zero-field frequency scan
        bgd0(lidx:ridx) = mag(lidx:ridx,Hp0);
        bf0 = freq(:,bidx);
        figure
        plot(freq(:,bidx),mag(:,bidx))
        hold on
        plot(freq(:,Hp0),mag(:,Hp0))
        plot(freq(:,bidx),bgd0);
        xlabel('Frequency (GHz)')
        ylabel('S11')
        legend('B = 0',num2str(HH(Hp0,bidx),'B = %.2f T'),'Stitched','location','SouthEast')
        fprintf('Checkpoint: press any key to continue, or Ctrl + c to abort.\n')
        pause
        fig_norm = figure;
        plot( freq(:,bidx), mag(:,bidx)./bgd0);
        break 
    elseif Options.bgdmode == 2 % load background data from existing file.
        fprintf('Loading background data...\n')
        [backpath,backfile,~] = fileparts(LoadObj);
        if isfile(fullfile(backpath,[backfile,'_interp.mat']))
            load(fullfile(backpath,[backfile,'_interp.mat']),'background');
            temp_load = load(fullfile(backpath,[backfile,'_interp.mat']),'analysis');
            analysis.kpe0 = temp_load.analysis.kpe0;
            analysis.kpi0 = temp_load.analysis.kpi0;
            analysis.wc0 = temp_load.analysis.wc0;            
            if exist('background','var')
                bf0 = freq(:,bidx); 
                bgd0 = interp1(background.f, background.d, bf0);  
                fig_norm = figure;
                plot( freq(:,bidx), mag(:,bidx)./bgd0);
                break
            else
                fprintf('Background data not found! Please choose normalization mode\n')
                prompt = sprintf('0. no normalization; 1. Stiched background; 3. Normalization from |S11| fit\n');
                Options.bgdmode = lower(input(prompt));
            end
        else
            fprintf('Background data not found! Please choose normalization mode\n')
            prompt = sprintf('0. no normalization; 1. Stiched background; 3. Normalization from |S11| fit\n');
            Options.bgdmode = lower(input(prompt));
        end
    elseif Options.bgdmode == 3
        % further normalization through |S11| fitting
        fprintf('Normalization by |S11| fitting...\n')
        % starting value for param = {'kpe', 'kpi', 'omega', 'Gc', 'gma', 'wc', 'attn'}
        param  =  [1e-3   1e-3   freq_h    0   1e-3    freq(fidx,bidx)   mean(mag(:,bidx))];
        bound_l = [ 0      0     freq_h    0   1e-3    freq_l   0.5]; % lower bound of fitting parameters
        bound_h = [Inf    Inf    freq_h    0   1e-3    freq_h   1.5]; % upper bound of fitting parameters
        fit = iptopt(freq(:,bidx), mag(:,bidx), HH(1,bidx), param, bound_l, bound_h, 1./mag(:,bidx), true);
        bf0 = freq(:,bidx);
%         bgd0 = medfilt1(mag(:,bidx)-fit(freq(:,bidx))+1)*fit.attn;
        bgd0 = medfilt1(mag(:,bidx)./fit(freq(:,bidx)))*fit.attn;
        figure;        
        plot(freq(:,bidx),mag(:,bidx));
        hold on
        plot(bf0,bgd0);
        xlim([freq_l freq_h])
        xlabel('Frequency (GHz)')
        ylabel('S11')
        legend('Background')
        legend('Zero field raw data','Fitted background', 'Location', 'SouthEast');
        fprintf('Checkpoint: press any key to continue, or Ctrl + c to abort.\n')
        pause
        fig_norm = figure;
        plot(freq(:,HH(1,:)==min(HH(1,:)) ),mag(:,HH(1,:)==min(HH(1,:)) )./bgd0)
        break
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

% Renormailize the background
if Options.bgdmode ~= 0
    analysis.bgd0 = bgd0;
    analysis.bf0 = bf0;
    bgdM = repmat(bgd0,1,size(mag,2));
    mag = mag./bgdM;
end

%Find all the resonant peaks
f0 = zeros(size(mag,2),1);
H0 = zeros(size(mag,2),1);
mag0 = zeros(size(mag,2),1);
%find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
for ii = 1:size(mag,2) %Searching column minima (fixed field)
    [~,idx] = min( mag(:,ii) );
    if(length(idx)>1)
        fprintf(num2str(H0(ii),'multiple minima found at H = %.3f\n'))
    end
    f0(ii) = freq(idx,ii);
    H0(ii) = HH(idx,ii);
    mag0(ii) = mag(idx,ii);
end 
% further normalization through |S11| fitting
% starting value for param = {'kpe', 'kpi', 'omega', 'Gc', 'gma', 'wc', 'attn'}
param  =  [1e-2   1e-2   freq_h    0   1e-3   f0(bidx)  mean(mag(:,bidx))];
bound_l = [ 0      0     freq_h    0   1e-3    freq_l   0.5]; % lower bound of fitting parameters
bound_h = [Inf    Inf    freq_h    0   1e-3    freq_h   1.5]; % upper bound of fitting parameters
fit = iptopt(freq(:,bidx), mag(:,bidx), H0(bidx), param, bound_l, bound_h, 1./mag(:,bidx), false);
bgd0 = mag(:,bidx)./fit(freq(:,bidx))*fit.attn;
% bgd0 = medfilt1(mag(:,bidx)./fit(freq(:,bidx)))*fit.attn;
if Options.bgdmode ~= 0
figure(fig_norm)
hold on
plot(freq(:,bidx),mag(:,bidx)./bgd0);
legend(num2str(H0(bidx),'Normalized with Loaded background at B = %.2f T'),...
    num2str(H0(bidx),'Further normalized with fitted background at B = %.2f T'),'Location', 'SouthEast');
end
prompt = sprintf('Implement additional normalization by fitting? \n');
add_fit = lower(input(prompt,'s'));
if strcmp(add_fit, 'yes') || strcmp(add_fit,'y') || num2str(add_fit) == 1
    fprintf('Further normalization by fitting...\n')
    bgdM = repmat(bgd0,1,size(mag,2));
    mag = mag./bgdM;
end
clearvars c idx ia lidx ridx widx ii HM trunc1 trunc2 dupl nop bf0 temp_load

%Find all the resonant peaks
f0 = zeros(size(mag,2),2);
mag0 = zeros(size(mag,2),2);
FWHM0 = zeros(size(mag,2),2);
%find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
for ii = 1:size(mag,2) %Searching column minima (fixed field)
    switch strnth
        case 'weak'
            [~,idx] = min(medfilt1(mag(:,ii)));
        case 'strong'
            [~,idx] = findpeaks(-medfilt1(mag(:,ii)),'MinPeakProminence',0.08,'MinPeakDistance',100);
            if isempty(idx); [~,idx] = min(mag(:,ii)); end
    end
    H0(ii) = HH(idx(1),ii);
    for jj = 1:length(idx)
        f0(ii,jj) = freq(idx(jj),ii);
        mag0(ii,jj) = mag(idx(jj),ii);
%         HM = ( max(mag(:,ii)) + mag0(ii,jj) )/2; % Half max point
%         HM = db2mag(-3); % find -3dB level
        HM = max(mag(:,ii)) - (max(mag(:,ii)) - mag0(ii,jj))*0.3; % 70% of the height
        left = mag(1:idx(jj),ii);
        right = mag(idx(jj):end,ii);
        [lidx,~] = find(left >= HM, 1, 'last');
        [ridx,~] = find(right >= HM, 1, 'first');
        if lidx == idx(jj)+ridx-1
            FWHM0(ii,jj) = Inf;
        elseif isempty(lidx) || isempty(ridx)
            FWHM0(ii,jj) = Inf;
        else
            FWHM0(ii,jj) = freq(idx(jj)+ridx-1,ii)-freq(lidx,ii);
        end
    end
end
% Reduce data noise
f0 = medfilt1(f0);
FWHM0 = medfilt1(FWHM0);
Q0 = f0./FWHM0;
Q0(isnan(Q0)) = 0;

% % Remove outlier data
% [f0,rmidx] = rmoutliers(f0,'mean',1);
% Q0 = Q0(~rmidx);
% H0 = H0(~rmidx);
% mag0 = mag0(~rmidx,:);
% FWHM0 = FWHM0(~rmidx);

analysis.mag0 = mag0;
analysis.FWHM0 = FWHM0;
analysis.Q0 = Q0;

trigger = 0;
wc = f0(1,1);
[~,Hp0] = min(Q0(:,1)); % Option 1: Use the minimum of Q0 to identify anticrossing location
% [~,Hp0] = max(mag0(:,1)); % Option 2: Use the maximum of mag0 to identify anticrossing location
% [~,Hp0] = min(abs(H0-0.33)); % Method 3. Manually set position of linecrossing
while true
    switch strnth
        case 'weak'
            % Fit the field dependent resonant frequency data with weak coupling function
% %             [~,idx] = max(f0);
            [~,idx] = findpeaks(f0(:,1)); % location anticrossing by inflection point
            [pfit,~,mu] = polyfit(H0(1:idx(1)),f0(1:idx(1),1),1); % use the segment before anticrossing
%             [pfit,~,mu] = polyfit(H0,f0(:,1),1); % use the segment over the full field range
            detrend = polyval(pfit,H0,[],mu) - pfit(2);
            init = [1e-2  1e-2   sgn*slope   wc  H0(Hp0)]; % [g  gamma slope wc x0]
            bound_l = [0   0  -Inf   freq_l  field_l];
            bound_h = [1   1   Inf   freq_h  field_h];
            wt = abs(f0(:,1)-f0(1,1))/max(abs(f0(:,1)-f0(1,1)));
%             wt = ones(size(f0));
            [fitP,~] = wk_cpl_fit(H0,f0(:,1)-detrend,init,bound_l,bound_h,wt,true);
%             [fitP,~] = wk_cpl_fit(H0,f0,init,bound_l,bound_h,wt,true);
            analysis.Ffit = fitP;
            slopes = fitP.slope;
            gma0 = fitP.gamma;
            wc = fitP.wc;
            % wc = 3.642;
            gcs = fitP.g;
            Brf = fitP.x0; % Level crossing location from perturbative fitting
            omg1 = (slopes + pfit(1)).*(H0-Brf) + wc;
            if aux == true % auxillary spin resonance
                omg2 = (slopes + pfit(1)).*(H0-Brf2) + wc; % create an auxillary spin dispersion
                gcs = [gcs gcs];
                fprintf('Auxiliary spin resonance generated...\n');
            end
            fprintf('Weak coupling fit complete\n')
            
            Tplot = figure;
            plot(H0,f0(:,1),'ok','MarkerSize',Options.mksz);
            hold on
            plot(H0,fitP(H0) + detrend);
            lgd = {'Exp. data','Theo. data'};
            xlabel('Magnetic Field (T)');
            ylabel('Resonant frequency (GHz)');
            title(num2str(Temperature,'Resonant frequency from minimum search at T = %3.3f K'));
            axis([field_l field_h freq_l freq_h]);
            break
        case 'strong'
            % split the resonant frequencies into upper and lower branches
            if any(f0(:,2)) % search if there is any overlapping region in avoide crossing pattern
                [bkr1,~] = find(f0(:,2) ~= 0,1,'first');
                [bkr2,~] = find(f0(bkr1:end,2) == 0,1,'first');
                bkr2 = bkr1 + bkr2 - 2;
                [~,Hp0] = min(abs(mag0(bkr1:bkr2,1)-mag0(bkr1:bkr2,2)));
                Hp0 = Hp0 + bkr1;
                if sgn > 0 % for dispersion with positive slope
                    upB = cat(1,f0(1:bkr1-1,1),f0(bkr1:bkr2,2)); % upper branch
                    upH = H0(1:bkr2);
                    loB = f0(bkr1:end,1); % lower branch
                    loH = H0(bkr1:end);
                    Hx1 = linspace(min(upH),max(upH),501); % for trace extrapolation
                    Hx2 = linspace(min(loH),max(loH),length(Hx1));
                else % for dispersion with negative slope
                    upB = cat(1,f0(bkr1:bkr2,2),f0(bkr2+1:end,1)); % upper branch
                    upH = H0(bkr1:end);
                    loB = f0(1:bkr2,1); % lower branch
                    loH = H0(1:bkr2);
                    Hx1 = linspace(min(loH),max(loH),501); % for trace extrapolation
                    Hx2 = linspace(min(upH),max(upH),length(Hx1));
                end
            else % if no overlap found
                if sgn > 0 % for dispersion with positive slope
                    [~, lidx] = max(f0(:,1));
                    [~, ridx] = min(f0(lidx:end,1));
                    upB = f0(1:lidx,1);
                    upH = H0(1:lidx);
                    loB = f0(lidx+ridx-1:end,1);
                    loH = H0(lidx+ridx-1:end);
                else
                    [~, lidx] = min(f0(:,1));
                    [~, ridx] = max(f0(lidx:end,1));
                    loB = f0(1:lidx,1);
                    loH = H0(1:lidx);
                    upB = f0(lidx+ridx-1:end,1);
                    upH = H0(lidx+ridx-1:end);
                end
                Hx1 = linspace(field_l,H0(lidx+ridx-1),501);
                Hx2 = linspace(H0(lidx),field_h,length(Hx1));
            end
            % Plot Resonant frequency trace with fitting parameters
            Tplot = figure;
            plot(upH,upB,'ok',loH,loB,'ok','MarkerSize',Options.mksz);
            % plot(H0(1:round(length(H0)/200):end),f0(1:round(length(f0)/200):end),'ok','MarkerSize',Options.mksz);
            hold on
            if aux == false % auxiliary spin resonance option
                %% single dispersion split fit-----------------------
                init = [0.05   sgn*slope   wc  H0(Hp0)]; % [g slope wc x0]
                bound_l = [ 0   -Inf   freq_l   field_l];
                bound_h = [Inf   Inf   freq_h   field_h];
                wt = ones(size(loB));
                [lofit,~] = str_cpl_fit(loH,loB,-1,init,bound_l,bound_h,wt,false);
                
                init = [lofit.g   lofit.slope   lofit.wc  lofit.x0]; % [g slope wc x0]
                bound_l = [ 0   -Inf   freq_l   field_l];
                bound_h = [Inf   Inf   freq_h   field_h];
                wt = ones(size(upB));
                [upfit,~] = str_cpl_fit(upH,upB,+1,init,bound_l,bound_h,wt,false);
                wc = mean([lofit.wc upfit.wc]);
                Brf = mean([lofit.x0 upfit.x0]);
                slopes = mean([lofit.slope upfit.slope]);
                %% single dispersion optimization fit-----------------------------------
                %             fitting parameter = ['wc', 'slope', 'x0', 'g']
                init = [wc   slopes   Brf  0.05];
                bound_l = [freq_l   -Inf   field_l   0];
                bound_h = [freq_h    Inf   field_h   1];
                if sgn > 0
                    brch1 = interp1(upH,upfit(upH),Hx1,'spline','extrap');
                    brch2 = interp1(loH,lofit(loH),Hx2,'spline','extrap');
                    front = @(param, x) (param(1) + param(2)*(x-param(3))./2 + sqrt((param(2)*(x-param(3))).^2 + 4*param(4)^2)./2);
                    back = @(param, x) (param(1) + param(2)*(x-param(3))./2 - sqrt((param(2)*(x-param(3))).^2 + 4*param(4)^2)./2);
                else
                    brch1 = interp1(loH,lofit(loH),Hx1,'spline','extrap');
                    brch2 = interp1(upH,upfit(upH),Hx2,'spline','extrap');
                    front = @(param, x) (param(1) + param(2)*(x-param(3))./2 - sqrt((param(2)*(x-param(3))).^2 + 4*param(4)^2)./2);
                    back = @(param, x) (param(1) + param(2)*(x-param(3))./2 + sqrt((param(2)*(x-param(3))).^2 + 4*param(4)^2)./2);
                end
                branches = @(param, x) [front(param, x(:,1)), back(param, x(:,2))];
                fitP = lsqcurvefit(branches, init, [Hx1', Hx2'], [brch1', brch2'], bound_l, bound_h); % fit extrapolated data
%                fitP = lsqcurvefit(branches, init, [H0(lidx-300:lidx), H0(lidx+uidx-1:end)],...
%                    [f0(lidx-300:lidx), f0(lidx+uidx-1:end)], bound_l, bound_h); % fit the original data

%                 wc = fitP(1);
                slopes = fitP(2);
                Brf = fitP(3);
                gcs = fitP(4);
                omg1 = slopes.*(H0-Brf) + wc;
                
                analysis.wfit1 = brch1;
                analysis.wfit2 = brch2;
            else
                %% double dispersion split fit--------------------------------------
                wt = ones(size(H0));
                init = [0.05   0.01   -slope  -slope   mean(f0)  H0(Hp0)  H0(Hp0)]; % [g slope wc x0]
                bound_l = [0   0   -Inf   -Inf    freq_l   H0(lidx)   H0(lidx)];
                bound_h = [1   1    Inf    Inf    freq_h   H0(lidx+uidx)   H0(lidx+uidx)];
                [upfit,~] = str_cpl_ufit2(H0(lidx+uidx-1:end),f0(lidx+uidx-1:end),init,bound_l,bound_h,wt(lidx+uidx-1:end),false);
                
                init = [upfit.g   upfit.g1   upfit.slope  upfit.slope1   upfit.wc  upfit.x0  upfit.x1]; % [g slope wc x0]
                bound_l = [0   0   -Inf   -Inf    freq_l   H0(lidx)   H0(lidx)];
                bound_h = [1   1    Inf    Inf    freq_h   H0(lidx+uidx)   H0(lidx+uidx)];
                [lofit,~] = str_cpl_lfit2(H0(1:lidx),f0(1:lidx),init,bound_l,bound_h,wt(1:lidx),false);
                
                wc = mean([upfit.wc lofit.wc]);
                slopes = [lofit.slope lofit.slope1];
                
                init = [mean([lofit.g upfit.g])   mean([lofit.g1 upfit.g1])   mean(slopes)  mean(slopes)   wc  upfit.x0  upfit.x1]; % [g slope wc x0]
                bound_l = [0    0    min([lofit.slope upfit.slope])   min([lofit.slope1 upfit.slope1])    freq_l   lofit.x0  lofit.x1];
                bound_h = [1    1    max([lofit.slope upfit.slope])   max([lofit.slope1 upfit.slope1])    freq_h   lofit.x0  lofit.x1];
                [upfit,~] = str_cpl_ufit2(H0(lidx+uidx-1:end),f0(lidx+uidx-1:end),init,bound_l,bound_h,wt(lidx+uidx-1:end),false);
                
                init = [mean([lofit.g upfit.g])   mean([lofit.g1 upfit.g1])   mean([lofit.slope upfit.slope])  mean([lofit.slope1 upfit.slope1])   wc  upfit.x0  upfit.x1]; % [g slope wc x0]
                bound_l = [0    0    min([lofit.slope upfit.slope])   min([lofit.slope1 upfit.slope1])    mean([lofit.wc upfit.wc])   lofit.x0  lofit.x1];
                bound_h = [1    1    max([lofit.slope upfit.slope])   max([lofit.slope1 upfit.slope1])    mean([lofit.wc upfit.wc])   lofit.x0  lofit.x1];
                [lofit,~] = str_cpl_lfit2(H0(1:lidx),f0(1:lidx),init,bound_l,bound_h,wt(1:lidx),false);
                
%             lofit_ex = interp1(H0(1:lidx),lofit(H0(1:lidx)),Hx,'spline','extrap');
%             upfit_ex = interp1(H0(lidx+uidx-1:end),upfit(H0(lidx+uidx-1:end)),Hx,'spline','extrap');

%             wc = mean([upfit.wc lofit.wc]);
%             gcs = [upfit.g upfit.g1];
%             Brf = [upfit.x0 upfit.x1];
%             slopes = [upfit.slope upfit.slope1];
%             [~,Hpos1] = min(abs(Hx-min(upfit.x0,upfit.x1)));
%             [~,Hpos2] = min(abs(Hx-max(upfit.x0,upfit.x1)));
%             omg1 = upfit.slope.*(Hx-upfit.x0) + upfit.wc;
%             omg2 = upfit.slope1.*(Hx-upfit.x1) + upfit.wc;
%             fprintf('Strong coupling fit complete\n')
                
                wc = lofit.wc;
                gcs = [lofit.g lofit.g1];
                Brf = [lofit.x0 lofit.x1];
                slopes = [lofit.slope lofit.slope1];
                [~,Hpos1] = min(abs(Hx-min(lofit.x0,lofit.x1)));
                %             [~,Hpos2] = min(abs(Hx-max(lofit.x0,lofit.x1)));
                omg1 = lofit.slope.*(Hx-lofit.x0) + lofit.wc;
                omg2 = lofit.slope1.*(Hx-lofit.x1) + lofit.wc;
                fprintf('Strong coupling fit complete\n')
                analysis.lofit = lofit;
                analysis.upfit = upfit;
                % double dispersion piecewise fit-----------------------------------
                %             fitting parameter : wc, slope, x0, g, slope1, x1, g1
                init = [mean(f0)   -slope   H0(Hpos)  0.05   -slope   H0(Hpos)  0.05];
                bound_l = [freq_l   -Inf   mean(f0)   0   -Inf   field_l   0];
                bound_h = [freq_h    Inf   mean(f0)   1    Inf   field_h   1];
                [pcfit,~] = str_cpl_pcfit2(H0,f0,H0(lidx),H0(lidx),init,bound_l,bound_h,wt,false);
                %             wt = ones(size(Hx));
                %             [pcfit,~] = str_cpl_pcfit2(Hx, (Hx <= H0(lidx+uidx-1)).*lofit_ex + (Hx >= H0(lidx)).*upfit_ex,...
                %                 H0(lidx),H0(lidx),init,bound_l,bound_h,wt,true);
                analysis.pcfit = pcfit;
                
                gcs = [pcfit.g  pcfit.g1];
                slopes = [pcfit.slope   pcfit.slope1];
                Brf = [pcfit.x0    pcfit.x1];
                wc = pcfit.wc;
                [~,Hpos1] = min(abs(Hx-pcfit.x0));
                % [~,Hpos2] = min(abs(Hx-pcfit.x1));
                omg1 = pcfit.slope.*(Hx-pcfit.x0) + pcfit.wc;
                omg2 = pcfit.slope1.*(Hx-pcfit.x1) + pcfit.wc;
            end
%% -----------------------------------
            figure(Tplot)
            if exist('brch1','var')
                lf = plot(Hx1, brch1,'-r');
                set(lf,'LineWidth',Options.lnwd);
                if exist('brch2','var')
                    uf = plot(Hx2, brch2,'-r');
                    set(uf,'LineWidth',Options.lnwd);
                end
            end
            lgd = {'Exp. data','Theo. data',''};
            xlabel('Magnetic Field (T)');
            ylabel('Resonant frequency (GHz)');
            title(num2str(Temperature,'Resonant frequency from minimum search at T = %3.3f K'));
            axis([field_l field_h freq_l freq_h]);
            break
        otherwise
            prompt = sprintf('Choose coupling strength regime: "strong" or "weak": \n');
            strnth = input(prompt,'s');
            sprintf('Choose the sign of the spin resonance (+1 or -1): \n');
            sgn = input();
            trigger = trigger + 1;
            if trigger >= 5
                clearvars tick
                break
            end
    end
end
% Fit Quality factor curve as a function of field
% [~,Hp0] = min(abs(H0-Brf)); % update anticrossing location with fitted value
init = [max(nonzeros(FWHM0))  mean(gcs)  1e-3   mean(slopes)  H0(Hp0)  wc];  % [gamma gc kappa slope x0 wc]
bound_l = [ 0    0    0   min([0 sgn*Inf])   field_l   wc];
bound_h = [Inf   1    1   max([0 sgn*Inf])   field_h   wc];
% fpts = min([length(H0(1:Hp0))-1 length(H0(Hp0:end))-1 20]);
% [Qfit0, ~] = Qf_fit(H0(Hp0-fpts:Hp0+fpts), Q0(Hp0-fpts:Hp0+fpts), init, bound_l, bound_h, false);
[Qfit0, ~] = Qf_fit(H0, Q0(:,1), init, bound_l, bound_h, false);
analysis.Qfit0 = Qfit0;

if ~exist('gma0','var'); gma0 = Qfit0.gamma; end
omgQ = Qfit0.slope.*(H0-Qfit0.x0) + Qfit0.wc; % linearized dispersion relation

figure(Tplot)
hold on
plot(H0,omg1,'--k');
switch strnth
    case 'weak'
        lgd{end+1} = sprintf('f0 fit, r = %.3f GHz', fitP.gamma);
    case 'strong'
        lgd{end+1} = sprintf('f0 fit, gc = %.3f GHz', gcs(1));
end
if exist('omgQ','var')
    plot(H0,omgQ,'-.r');
    lgd{end+1} = sprintf('Q0 fit, r = %.3f GHz', Qfit0.gamma);
end
if exist('omg2','var')
    plot(H0,omg2,'--b');
    lgd{end+1} = sprintf('f0 fit, gc = %.3f GHz',gcs(2)); 
end
legend(lgd,'location', 'southwest')

%Plot Quality factor from minimum search vs magnetic field
Qplot = figure;
plot(H0, Q0(:,1),'s-','MarkerSize',Options.mksz);
% plot(H0(1:round(length(H0)/20):end), Q0(1:round(length(Q0)/20):end),'s-','MarkerSize',Options.mksz);
hold on
Qp0 = plot(Qfit0,'-k');
set(Qp0,'LineWidth',Options.lnwd);
legend('Exp. data', sprintf('B0 = %1$.2f T\n r = %2$.2e GHz\n Gc = %3$.2e GHz',...
    Qfit0.x0, Qfit0.gamma, Qfit0.gc),'location','SouthWest')

prompt = 'Choose initial guess for spin resonance (1. omg1, 2. omgQ): ';
if input(prompt) == 2
    omg1 = omgQ;
    [~,Hp0] = min(abs(H0-Qfit0.x0));
    gcs = Qfit0.gc;
    wc = Qfit0.wc;
end

clearvars B Delt hPara fPara wp wm fitPara H_res f_res Qparam

figure
hold on
box on
cmap0 = pcolor(HH,freq,mag2db(mag));
set(cmap0, 'edgeColor','none')
shading interp;
colorbar
caxis('auto')
% caxis([max([mag2db(min(mag0)), -30]) max(mag2db(max(mag0)),0)+1]);
clearvars init bound_l bound_h *_temp

set(gca,'fontsize',Options.ftsz)
axis([field_l field_h freq_l freq_h]);
xlabel('Magnetic Field (T)');
ylabel('Frequency (GHz)');
title(num2str(Temperature,'S11 response at T = %3.3f K'));
        
% Plot frequency scan at line crossing vs. zero-field scan with fittings
Cplot = figure;
[~,Hl] = min(H0);
plot(freq(1:10:end,Hl),mag(1:10:end,Hl),'o');
hold on
plot(freq(1:10:end,Hp0),mag(1:10:end,Hp0),'or');
xlabel('Frequency (GHz)');
ylabel('S11');
clgd = {sprintf('Frequency cut at %.2f T',min(H0) ), sprintf('Frequency cut at %.2f T',H0(Hp0))};

%% Step 1: Fit the data for the first time to extract "kpe" and "w0" far from the level crossing
ff0 = double.empty(length(H0),2,0); 
Qf = double.empty(length(H0),0);
Gc1 = double.empty(length(H0),0);
if exist('omg2','var'); Gc2 = double.empty(length(H0),0); end

kpe = double.empty(length(H0),0);
kpe(1,1) = FWHM0(1)/10;
kpi = double.empty(length(H0),0);
kpi(1,1) = FWHM0(1)/10;

gamma = double.empty(length(H0),0);
gamma_ci = double.empty(length(H0),2,0);
Gc1_ci = double.empty(length(H0),2,0);

wcs = double.empty(length(H0),0);
wcs_ci = double.empty(length(H0),2,0);
omg1_ci = double.empty(length(H0),2,0);
if exist('omg2','var')
    omg2_ci = double.empty(length(H0),2,0);
    Gc2_ci = double.empty(length(H0),2,0);
end
attn = double.empty(length(H0),0);

[~,Hc0] = min(H0);
Hc1 = Hc0 + 1;
weight(:,:,1) = abs(1./mag(:,:));      % Weight function option 1
% weight = 1./gradientweight(mag);       % Weight function option 2
% weight(isinf(weight)) = 1e4;           % Remove infinities in weight function
% weight = ones(size(mag));                % Weight function option 3 (uniform)
if exist('omg2','var')
    % starting value for param = {'kpe', 'kpi', 'wc', 'w1', 'Gc1', 'w2', 'Gc2', 'gma','attn'}
    param  =  [kpe(Hc0)  kpi(Hc0)    wc    omg1(Hc0)  gcs(1)   omg2(Hc0)   gcs(2)    gma0   1/mean(mag(:,Hc0))];
    bound_l = [   0        0      freq_l   omg1(Hc0)  gcs(1)   omg2(Hc0)   gcs(2)    gma0    0.8]; % lower bound of fitting parameters
    bound_h = [  Inf      Inf     freq_h   omg1(Hc0)  gcs(1)   omg2(Hc0)   gcs(2)    gma0    1.2]; % upper bound of fitting parameters
    fit = iptopt2(freq(:,Hc0), mag(:,Hc0), H0(Hc0), param, bound_l, bound_h, weight(:,Hc0), false);
else
    % starting value for param = {'kpe', 'kpi', 'omega', 'Gc', 'gma', 'wc', 'attn'}
    param = [kpe(Hc0)  kpi(Hc0)  omg1(Hc0)   gcs(1)    gma0     wc     1/mean(mag(:,Hc0))];
    bound_l = [ 0         0      omg1(Hc0)   gcs(1)    gma0   freq_l   0.8]; % lower bound of fitting parameters
    bound_h = [Inf       Inf     omg1(Hc0)   gcs(1)    gma0   freq_h   1.2]; % upper bound of fitting parameters
    fit = iptopt(freq(:,Hc0), mag(:,Hc0), H0(Hc0), param, bound_l, bound_h, weight(:,Hc0), false);
end
conf_intv = confint(fit,0.95);
kpe(Hc0,1) = fit.kpe;
kpi(Hc0,1) = fit.kpi;
wcs(Hc0,1) = fit.wc;
wcs_ci(Hc0,:,1) = conf_intv(:,6);
attn(Hc0,1) = fit.attn;

[~,idx] = findpeaks(-fit(freq(:,Hc0))); % Find the resonant frequencies
if ~isempty(idx)
    if length(idx) == 2
        ff0(Hc0,1,1) = freq(idx(1),Hc0);  
        ff0(Hc0,2,1) = freq(idx(2),Hc0);
    elseif length(idx) > 2
        [~,idxt] = sort(mag(:,Hc0));
        ff0(Hc0,1,1) = freq(idxt(1),Hc0);
        ff0(Hc0,2,1) = freq(idxt(2),Hc0);
    else
        ff0(Hc0,:,1) = freq(idx,Hc0);
    end
else
    ff0(Hc0,:,1) = fit.wc;
end
Qf(Hc0) = mean(ff0(Hc0,:,1)./(fit.kpe + fit.kpi));
figure(Cplot)
plot(freq(:,Hc0),fit(freq(:,Hc0)),'-k');
clgd{end+1} = '';

%% Step 2: fit the data for the second time with fixed "kpe" and "w0"
% Use fit far away from anticrossing for dissipation rates
kpe0 = kpe(Hc0);
kpi0 = kpi(Hc0);
wc0 = wcs(Hc0);
attn0 = attn(Hc0);
gamma(1,1) = gma0;
for ii = Hc1:length(H0)
    % Selectively plot the fit for visual inspection
    plt = false;
%     figWin = Hp0-10:2:Hp0+10; % The iteration window over which shows fitting graph
%     if ismember(ii,figWin) % Not useable in parallel mode (parfor)
%         plt = true;
%     end
    if exist('omg2','var')
        % starting value for param = {'kpe', 'kpi', 'wc', 'w1', 'Gc1', 'w2', 'Gc2', 'gma','attn'}
        param  =  [kpe0   kpi0   wc0  omg1(ii)  gcs(1)  omg2(ii)  gcs(2)  gamma(ii-1)  attn0];
        bound_l = [kpe0   kpi0   freq_l   -Inf    0   -Inf    0     0    0.8]; % lower bound of fitting parameters
        bound_h = [kpe0   kpi0   freq_h    Inf   Inf   Inf   Inf   Inf   1.2]; % upper bound of fitting parameters
        fit = iptopt2(freq(:,ii), mag(:,ii), H0(ii), param, bound_l, bound_h, weight(:,ii), plt);
        conf_intv = confint(fit,0.95);
        wcs_ci(ii,:,1) = conf_intv(:,3);
        omg1(ii) = fit.w1;
        omg1_ci(ii,:,1) = conf_intv(:,4);
        Gc1(ii) = fit.Gc1;
        Gc1_ci(ii,:,1) = conf_intv(:,5);
        omg2(ii) = fit.w2;
        omg2_ci(ii,:,1) = conf_intv(:,6);
        Gc2(ii) = fit.Gc2;
        Gc2_ci(ii,:,1) = conf_intv(:,7);
        gamma_ci(ii,:,1) = conf_intv(:,8);
    else
        % starting value for param = {'kpe', 'kpi', 'omega', 'Gc', 'gma', 'wc', 'attn'}
        param  =  [kpe0   kpi0   omg1(ii)  gcs(1)   gamma(ii-1)    wc0    attn0];
        bound_l = [kpe0   kpi0    -Inf    0     0    freq_l  0.8]; % lower bound of fitting parameters
        bound_h = [kpe0   kpi0     Inf   Inf   Inf   freq_h  1.2]; % upper bound of fitting parameters
        fit = iptopt(freq(:,ii), mag(:,ii), H0(ii), param, bound_l, bound_h, weight(:,ii), plt);
        conf_intv = confint(fit,0.95);
        omg1(ii) = fit.omega;
        omg1_ci(ii,:,1) = conf_intv(:,3);
        gamma_ci(ii,:,1) = conf_intv(:,5);
        wcs_ci(ii,:,1) = conf_intv(:,6);
        Gc1(ii) = fit.Gc;
        Gc1_ci(ii,:,1) = conf_intv(:,4);
    end
    kpe(ii,1) = fit.kpe;
    kpi(ii,1) = fit.kpi;
    gamma(ii,1) = fit.gma;
    wcs(ii,1) = fit.wc;
    attn(ii,1) = fit.attn;
    
    if ii == Hc0 || ii == Hp0
        figure(Cplot)
        plot(freq(:,ii),fit(freq(:,ii)),'-k');
        clgd{end+1} = sprintf('Fitting: G_c = %1$.2e GHz, r = %2$.2e GHz.',Gc1(ii),gamma(ii,1));
        legend(clgd,'location','SouthEast')
    end
    
    switch strnth
        case 'weak' % Anticipate no branching in weak coupling regime
            [~,idx] = min(fit(freq(:,ii)));
            ff0(ii,:,:) = freq(idx,ii);
        case 'strong'
            [~,idx] = findpeaks(-fit(freq(:,ii))); % Find the resonant frequency with peak finder
            if ~isempty(idx)
                if length(idx) == 2
                    ff0(ii,1,1) = freq(idx(1),ii);
                    ff0(ii,2,1) = freq(idx(2),ii);
                elseif length(idx) > 2
                    [~,idxt] = min(mag(idx,ii));
                    ff0(ii,:,1) = freq(idx(idxt),ii);
                else
                    ff0(ii,:,1) = freq(idx,ii);
                end
            else
%                ff0(ii,1,1) = min(fit(freq(:,ii)));  % alternatively by minimum search
                ff0(ii,1,1) = ff0(ii-1,1,1);
            end
    end
    Qf(ii) = mean(ff0(ii,:,1)./(kpe0+kpi0)); % Calculate quality factor
end
analysis.attn = attn;
analysis.kpe = kpe;
analysis.kpi = kpi;
analysis.wc = wcs;
analysis.wc_ci = wcs_ci;
analysis.wr = ff0;
analysis.w0 = omg1;
analysis.w0_ci = omg1_ci;
if exist('omg2','var'); analysis.wa = omg2; analysis.wa_ci = omg2_ci; end
analysis.gma = gamma;
analysis.gma_ci = gamma_ci;
analysis.Qf = Qf;
% analysis.Qfit = Qfit;
analysis.Gc1 = Gc1;
analysis.Gc1_ci = Gc1_ci;
if exist('omg2','var')
    analysis.Gc2 = Gc2;
    analysis.w1 = omg2;
    analysis.Gc2_ci = Gc2_ci;
end

% figure
% plot(H0(Hc0:end),kpe(Hc0:end),'o','MarkerSize',Options.mksz);
% ylabel('K_e (GHz)');
% hold on
% yyaxis right
% plot(H0(Hc0:end),kpi(Hc0:end),'o','MarkerSize',Options.mksz);
% xlabel('Magnetic Field (T)');
% ylabel('K_i (GHz)');
% title('Dissipation rates');
% legend('External dissipation rate','Internal dissipation rate')

% figure
% plot(H0(Hc0:end),medfilt1(kpi(Hc0:end)./kpe(Hc0:end)),'-');
% ylabel('K_i/K_e');
% xlabel('Magnetic Field (T)');
% title('Ratio of Dissipation rates');

% Plot the resonant frequency from Lorentzian fit versus DC magnetic field
figure
cmap2 = copyobj(cmap0, gca);
set(cmap2, 'edgeColor','none');
shading interp;
colorbar
% caxis([max([mag2db(min(mag0)), -30]) max(mag2db(max(mag0)),0)+1]);
caxis('auto');
hold on
plot(H0(Hc0:end),ff0(Hc0:end,:),'or','MarkerSize',Options.mksz,'MarkerFaceColor','red');
plot(H0,omg1,'.k');
plot(H0,omg1_ci,'-r');
xlabel('Magnetic Field (T)');
ylabel('Frequency (GHz)');
title(num2str(Temperature,'Data with fitted paramters at T = %3.3f K'));
axis([field_l field_h freq_l freq_h]);

%copy over the colormap of the raw data as the background of the next plot
% figure
% cmap3 = copyobj(cmap0, gca);
% set(cmap3, 'edgeColor','none');
% shading interp;
% colorbar
% % caxis([max([mag2db(min(mag0)), -30]) max(mag2db(max(mag0)),0)+1]);
% caxis('auto')
% hold on
% f0 = medfilt1(f0); % apply median filter to remove some noise
% % hfig1 = plot(H0(1:round(length(H0)/200):end),f0(1:round(length(f0)/200):end),'ok','MarkerSize',Options.mksz,'MarkerFaceColor','black');
% plot(H0, f0, 'ok', 'MarkerSize', Options.mksz); % Plot the resonant frequency from minimum search
% xlabel('Magnetic Field (T)');
% ylabel('Resonant frequency (GHz)');
% title(num2str(Temperature,'Resonant frequency from minimum search at T = %3.3f K'));
% % caxis('auto');
% axis([field_l field_h freq_l freq_h]);

% Plot the peak amplitude from minimum search vs. magnetic field
figure
mag0 = medfilt1(mag0); % apply median filter to remove some noise
plot(H0, mag2db(mag0), 'o', 'MarkerSize', Options.mksz);
xlabel('Field(T)');
ylabel('S11 (dB)');
title(num2str(Temperature,'Minimal S11 at T = %3.3f K'));

% Plot Quality factor from fitting function vs magnetic field
% init = [max(gamma)  mean(gcs)  kpe0   mean(slopes)  H0(Hp0)  wc];  % [gamma gc kappa slope x0 wc]
% bound_l = [0   0   0   -Inf   field_l   freq_l];
% bound_h = [1   1   1    Inf   field_h   freq_h];
% [Qfit, ~] = Qf_fit(H0(Hp0-15:Hp0+15), Qf(Hp0-15:Hp0+15), init, bound_l, bound_h, false);
figure(Qplot)
% yyaxis right
% hold on
% H0 = H0(Qf >=0);
% Qf = Qf(Qf >=0); %Remove unphysical points
% Qf = medfilt1(Qf);
% plot(H0, Qf,'ro-','MarkerSize',Options.mksz);
% % plot(H0(Hc2:round(length(H0)/20):Hc2), Qf(Hc2:round(length(Qf)/20):Hc2),'o-','MarkerSize',Options.mksz);
% % plot(H0, Qf,'o','MarkerSize',Options.mksz);
% Qp = plot(Qfit,'-k');
% set(Qp,'LineWidth',Options.lnwd);
% gca;
xlabel('Magnetic Field (T)');
ylabel('Q factor');
% legend('Quality factor from FWHM', 'Q0 fit', 'Quality factor from fitting', 'Qf fit');
title(num2str(Temperature,'Quality factor, T= %.2f K'));

gplot = figure;
errorbar(H0(Hc1:end),Gc1(Hc1:end),abs(Gc1_ci(Hc1:end,1)-Gc1_ci(Hc1:end,2)),'o');
if exist('omg2','var')
    figure(gplot)
    hold on
%     plot(H0(Hc1:end),medfilt1(Gc2(Hc1:end)),'o');
    plot(H0(Hc1:end),Gc2(Hc1:end),'o');
%     errorbar(H0(Hc1:end),Gc2(Hc1:end),abs(Gc2_ci(Hc1:end,1)-Gc2_ci(Hc1:end,2)),'o');
    lg = {'Gc1','Gc2'};
else
    lg = {'Gc'};
end
xlabel('Magnetic Field (T)')
ylabel('Coupling strength (GHz)')
title(num2str(Temperature,'Fitting paramters from input-output formalism at T = %.2f K'))
hold on
yyaxis right
% plot(H0(Hc1:end), medfilt1(gamma(Hc1:end)),'s')
errorbar(H0(Hc1:end), gamma(Hc1:end),abs(gamma_ci(Hc1:end,1)-gamma_ci(Hc1:end,2)),'s')
lg{end+1} = '\gamma';
ylabel('Spin linewidth (GHz)')
legend(lg)
% figure
% plot(Hx, B0, '-o')
% xlabel('Magnetic Field (T)')
% ylabel('Line cross field (T)')
end
% End of option 2 (on-resonance fitting)
function [HH,freq,mag,analysis] = option3(varargin)
%% Set data range and parameters
[HH,freq,mag,~,Options,analysis] = load_data(varargin{:});
Temperature = analysis.temp;

% field_l = 2.0; % Manually set field range
% field_h = 5.0;
field_l = min(HH(1,:));
field_h = max(HH(1,:));
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
    plot(analysis.Hs,medfilt1(analysis.Ts),'o-')
    xlabel('DC Magnetic field')
    ylabel('Temperature')
    title('Magnetic field vs Temperature')
end
%% Clean up raw data
clearvars dif step *_temp

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

if Options.bgdmode == 1 % Construct the background noise by stitching together zero-field-scan and anti-crossing
    [~,cidx] = min(mag(:,bidx)); % find the resonant peak to be replaced
%     [~,Hpos] = max(mag0); % Method 1. scan of max reflection
    [~,Hpos] = max(medfilt1(abs(f0-f0(bidx)))); % Method 2. scan with the resonant peak deviates from f0 the furtherest
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
    fprintf('Press any key to continue, or Ctrl + C to abort.\n')
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
    bgd0 = bgd0(trimIdx);
    bgd0 = interp1(bf0,bgd0,freq(:,1));
    bf0 = freq(:,1);
    hold on
    plot(bf0,bgd0);
    legend('Loaded data','Interpolated data');
    fprintf('Press any key to continue\n')
    pause
elseif Options.bgdmode == 3 % Background normalization through |S11| fitting
        % further normalization through |S11| fitting
        fprintf('Normalization by |S11| fitting...\n')
        % starting value for param = {'kpe', 'kpi', 'omega', 'Gc', 'gma', 'wc', 'attn'}
        param  =  [1e-3   1e-3   freq_h    0   1e-3   f0(bidx)  1/mean(mag(:,bidx))];
        bound_l = [ 0      0     freq_h    0   1e-3    freq_l   0.1]; % lower bound of fitting parameters
        bound_h = [Inf    Inf    freq_h    0   1e-3    freq_h   2.0]; % upper bound of fitting parameters
        fit = iptopt(freq(:,bidx), mag(:,bidx), H0(bidx), param, bound_l, bound_h, 1./mag(:,bidx), true);
        bf0 = freq(:,bidx);
        bgd0 = medfilt1(mag(:,bidx)-fit(freq(:,bidx))+1)*fit.attn;
%         bgd0 = medfilt1(mag(:,bidx)./fit(freq(:,bidx)))*fit.attn;
        [bf0,trimIdx] = unique(bf0);
        bgd0 = bgd0(trimIdx);
        figure
        plot(freq(:,bidx),mag(:,bidx));
        hold on
        plot(bf0,bgd0);
        xlim([freq_l freq_h])
        xlabel('Frequency (GHz)')
        ylabel('S11')
        legend('Background')
        legend('Zero field raw data','Fitted background', 'Location', 'SouthEast');
        fprintf('Checkpoint: press any key to continue, or Ctrl + c to abort.\n')
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

% Find all the resonant peaks
f0 = zeros(size(mag,2),1);
Q0 = zeros(size(mag,2),1);
% mag0 = zeros(size(zq,2),1);
FWHM = zeros(size(mag,2),1);
%find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
for ii = 1:size(mag,2) %Searching column minima (fixed field)
    [~,idx] = min( mag(:,ii) );
    if(length(idx)>1)
        fprintf(num2str(H0(ii),'multiple minima found at H = %.3f'))
    end
    H0(ii) = HH(idx,ii); 
    f0(ii) = freq(idx,ii);
%     mag0(ii) = zq(idx,ii);
    HM = ( max(mag(:,ii)) + min(mag(:,ii)) )/2; % 
    % Calculate quality factor using f0/FWHM
    if isnan(1/range(freq(mag(:,ii) <= HM)))
       Q0(ii) = 0;
    elseif isempty(range(freq(mag(:,ii) <= HM)))
       Q0(ii) = 0;
    else
       FWHM(ii) = range(freq(mag(:,ii) <= HM));
       Q0(ii) = freq(idx,ii)/FWHM(ii);
    end
end

% Cleaning data
Q0(isinf(Q0)) = NaN; % Cut out inf from the array
[Q0,c] = rmmissing(Q0); % Cut out NaN from the array
H0 = H0(~c); % Remove corresponding elements in H0 array as well
f0 = f0(~c);
% mag0 = mag0(~c);
FWHM = FWHM(~c);
HH = HH(:,~c);
freq = freq(:,~c);
mag = mag(:,~c);
% analysis.temp = analysis.temp(~c);

% For noisy data, we need to remove duplicates of minima
[H0,ia,~] = unique(H0,'stable');
f0 = f0(ia);
Q0 = Q0(ia);
% mag0 = mag0(ia);
FWHM = FWHM(ia);
HH = HH(:,ia);
freq = freq(:,ia);
mag = mag(:,ia);
% analysis.temp = analysis.temp(ia);

f0 = medfilt1(f0); % apply median filter to remove some noise
f0 = f0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
f_init = f0(bidx); % Initial guess of the cavity resonant frequency
H0 = H0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
Q0 = Q0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
% mag0 = mag0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
FWHM = FWHM(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
HH = HH(:,f0 >= freq_l & f0 <= freq_h);
freq = freq(:,f0 >= freq_l & f0 <= freq_h);
mag = mag(:,f0 >= freq_l & f0 <= freq_h);
% analysis.temp = analysis.temp(f0 >= freq_l & f0 <= freq_h);
clearvars c idx ia ii HM trunc1 trunc2 dupl nop trimIdx

weight = double.empty(size(mag,1),size(mag,2),0);
for jj = 1:size(mag,2)
    for ii = 1:size(mag,1)
        weight(ii,jj,1) = abs((mag(ii,jj)-max(mag(:,jj)))./mag(ii,jj));
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
        wc = double.empty(length(H0),0);
        fill = Options.filFctr;
        [~,Hc0] = min(H0);
        % starting value for param = {'kpe', 'kpi', 'omega', 'Gc', 'gma', 'wc', 'attn'}
        param =  [1e-4   1e-4    0     0    f_init]; % param = {'kpe', 'kpi', 'Re[chi]', 'Im[chi]', 'wc'}
        bound_l = [0      0      0     0    freq_l]; % lower bound of fitting parameters
        bound_h = [Inf   Inf     0     0    freq_h]; % upper bound of fitting parameters
        fit0 = iptoptx(freq(:,Hc0),mag(:,Hc0),H0(Hc0),fill,param,bound_l,bound_h,weight(:,Hc0),false);
        kpe(Hc0) = fit0.kpe;
        kpi(Hc0) = fit0.kpi;
        w0(Hc0) = fit0.wc;
        kpe0 = fit0.kpe;
        kpi0 = fit0.kpi;
        parfor ii = Hc0+1:length(H0)
% %         for ii = 1:length(H0)
            plt = false;
%             figWin = Hpos-100:20:Hpos+100; % The iteration window over which shows fitting graph
%             if ismember(ii,figWin) % Not useable in parallel mode (parfor)
%                 plt = true;
%             end
            % Fit using input-output formalism
            param =  [ kpe0  kpi0    2e-4   2e-4   f_init]; % param = {'kpe', 'kpi', 'Re[chi]', 'Im[chi]', 'wc'}
            bound_l = [kpe0  kpi0      0      0    f_init]; % lower bound of fitting parameters
            bound_h = [kpe0  kpi0    Inf    Inf    f_init]; % upper bound of fitting parameters
            fit = iptoptx(freq(:,ii),mag(:,ii),H0(ii),fill,param,bound_l,bound_h,weight(:,ii),plt);
            
            if mod(ii,20) == 0
                worker = getCurrentTask();
                fprintf('Fitting, current field: %1$3.2f T. Core %2$u.\n', H0(ii), worker.ID);
            end
          
            param = coeffvalues(fit);
            [idx,~,~] = find(fit(freq(:,ii)) == min(fit(freq(:,ii)))); % Find the resonant frequency by (unique) minimum search
            w0(ii) = freq(idx(1)); 
            kpe(ii) = param(1);
            kpi(ii) = param(2);
            xr(ii) = param(3);
            xi(ii) = param(4);   
            wc(ii) = param(5);
            Qf(ii) = w0(ii)/(param(1)+param(2));
        end
        analysis.kpe = kpe;
        analysis.kpi = kpi;
        analysis.w0 = w0;
        analysis.xr = xr;
        analysis.xi = xi;        
    case 2 %Option 2: use spec1d package to fit the data using Lorentzian form.
        parfor ii = 1:length(H0)
            s = spec1d(freq(:,ii), -mag(:,ii), max(-mag(:,ii)))*0.001; % create spec1d object
            %starting point for the (Lorentzian) fitting parameters
            p = [0.1 wc(ii) FWHM(ii) min(mag(:,ii))]; % (p1: scaling factor ,p2: resonant frequency; p3: FWHM; p4:noise floor(?) )
            fix = [0 0 0 0]; % Denoting if the fitting parameters are fixed
            [~, fbck] = fits(s, 'lorz', p, fix);
%             [~, fbck] = fits(s, 'lorz');
            wc(ii) = fbck.pvals(2); % Retrieve the resonant frequency from fitted data
            Qf(ii) = abs(fbck.pvals(2)/fbck.pvals(3)/2); %Calculate the quality factor
%             chi(ii) = 1/Qf(ii);
            if mod(ii,20) == 0
                worker = getCurrentTask();
                fprintf('Current magnetic field: %1$3.2f. on core %2$u.\n', H0(ii), worker.ID);
            end
        end
end
analysis.Qf = Qf;
analysis.wc = wc;

% Plot the interpolated frequency response data in a field scan using color map
figure
cmap = pcolor(HH,freq,mag2db(mag));
set(cmap, 'edgeColor','none')
shading interp;
colorbar
set(gca,'fontsize',Options.ftsz)
xlabel('Magnetic Field (T)');
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
xlabel('Magnetic Field (T)');
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
xlabel('Magnetic Field (T)');
ylabel('Frequency (GHz)');
xticks(linspace(field_l,field_h,6));
ylim([freq_l,freq_h]);
title(num2str(Temperature,'S11 response at T = %3.3f K'));
legend('Experimental data','|S11| fit');

figure
plot(H0,f0,'ok','MarkerSize',Options.mksz);
hold on
plot(mean(H0(f0==min(medfilt1(f0)))), mean(f0(f0==min(f0))), 'o', 'MarkerFaceColor', 'red','MarkerSize',Options.mksz);
xlabel('Magnetic Field (T)');
ylabel('Resonant frequency (GHz)');
title(num2str(Temperature,'Resonant frequency from minimum search at T = %3.3f K'));

figure
plot(H0,medfilt1(1./Q0),'.k');
yyaxis right
hold on
[~,idx,Qf] = find(Qf);
plot(HH(1,idx),medfilt1(1./Qf),'.r');
xlabel('Magnetic Field (T)');
ylabel('1/Q');
title(num2str(Temperature,'Inverse of Q factor at T = %3.3f K'));
legend('f_0 / FWHM','|S11| fitting');

if Options.fitfunc == 1
    figure
    plot(H0,medfilt1(xr),'-');
    ylabel('Susceptibility');
    xlabel('Magnetic Field (T)');
    title('Re[\chi]')
    
    figure
    plot(H0,medfilt1(xi),'-');
    ylabel('Susceptibility');
    xlabel('Magnetic Field (T)');
    title('Im[\chi]')
    
%     figure
%     plot(H0,medfilt1(kpe),'o');
%     ylabel('K_e (GHz)');
%     hold on
%     yyaxis right
%     plot(H0,medfilt1(kpi),'*');
%     xlabel('Magnetic Field (T)');
%     ylabel('K_i (GHz)');
%     title('Dissipation rates');
%     legend('External dissipation rate','Internal dissipation rate')
end
end
% End of option 3 (off-resonance fitting)
function [TT,freq,mag,analysis] = option4(varargin)
% extract data from raw data file
[~,freq,mag,~,Options,analysis] = load_data(varargin{:});
T1 = analysis.Ts;
TT = repmat(T1,size(freq,1),1);
H0 = mean(analysis.Hs);

freq_l = min(freq(:,1));
freq_h = max(freq(:,1));
temp_l = min(T1);  % set field range, l: lower limit, h: higher limit
temp_h = max(T1);

clearvars dif step lidx uidx T1 *_temp

% Find all the resonant peaks
f0 = zeros(size(mag,2),1);
T0 = zeros(size(mag,2),1);
mag0 = zeros(size(mag,2),1);
% find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
[~,lidx] = min(TT(1,:));
[~,uidx] = max(TT(1,:));
mag0l = min(mag(:,lidx)); % peak depth at the lowest field
mag0h = min(mag(:,uidx)); % peak depth at the highest field
[~,hl] = min([mag0l mag0h]); % pick the finner peak as the first part of the background
if hl == 1
    [~,tidx] = min( TT(1,:) );
    fprintf('Use low field data for background normalization...\n')
else
    [~,tidx] = max( TT(1,:) );
    fprintf('Use high field data for background normalization...\n')
end
for ii = 1:size(mag,2) %Searching column minima (fixed field)
    [~,idx] = min( mag(:,ii) );
    if(length(idx)>1)
        fprintf(num2str(T0(ii),'multiple minima found at H = %.3f T\n'))
    end
    f0(ii) = freq(idx,ii);
    T0(ii) = TT(idx,ii);
    mag0(ii) = mag(idx,ii);
end
clearvars uidx lidx

while true
    if Options.bgdmode == 1 % Construct the background noise by stitching together zero-field-scan and anti-crossing
        [~,cidx] = min(mag(:,tidx)); % find the resonant peak to be replaced
        [~,Tpos] = max(mag0); % Method 1. scan of max reflection
%         [~,Hpos] = max(mag(cidx,:)); % Method 2. scan of max reflection at f0 (frequency)
%         [~,Hpos] = max(medfilt1(abs(f0-f0(bidx)))); % Method 3. scan with the resonant peak deviates from f0 the furtherest
%         [~,Hpos] = min(abs(H0-14.2)); % Method 4. Manually set the patch slice of the frequency scan
        mag_dif = abs(mag(:,tidx)-mag(:,Tpos));
        [lidx,~] = find(mag_dif(1:cidx) <= 1e-3, 1, 'last'); % left boundary of the resonant peak to be replaced
        [ridx,~] = find(mag_dif(cidx:end) <= 1e-3, 1, 'first'); % right boundary of the resonant peak to be replaced
        if isempty(lidx); lidx = 1; end
        if isempty(ridx); ridx = length(mag_dif); else; ridx = cidx + ridx -1; end
        analysis.Tpos = Tpos; % store field location of the anti-crossing
        analysis.cidx = cidx; % peak center
        bf0 = freq(:,tidx);
        bgd0 = mag(:,tidx); % Base frequency scan
        bgd0(lidx:ridx) = mag(lidx:ridx,Tpos);

        figure
        plot(bf0,mag(:,tidx))
        hold on
        plot(freq(:,Tpos),mag(:,Tpos))
        plot(bf0,bgd0);
        xlabel('Frequency (GHz)')
        ylabel('S11')
        legend(num2str(T0(tidx),'B = %.2f T'),num2str(T0(Tpos),'T = %.2f K'),'Stitched background','location','southeast')
        break
    elseif Options.bgdmode == 2 % load background data from existing file
        fprintf('Loading background normalization data...\n')
        [backpath,backfile,~] = fileparts(LoadObj);
        if isfile(fullfile(backpath,[backfile,'_interp.mat']))
            load(fullfile(backpath,[backfile,'_interp.mat']),'background');
            if exist('background','var')
                bf0 = background.f;
                bgd0 = background.d;          
                bgd0 = interp1(bf0,bgd0,freq(:,1));
                bf0 = freq(:,1);
                
                figure
                plot(bf0,bgd0);
                hold on
                plot(freq(:,tidx),mag(:,tidx))
                xlim([freq_l freq_h])
                xlabel('Frequency (GHz)')
                ylabel('S11')
                legend('Loaded background','Raw data')
                break
            else
                fprintf('Background data not found! Please choose normalization mode:\n')
                prompt = sprintf('0. no normalization; 1. Stiched background; 3. Normalization from |S11| fit\n');
                Options.bgdmode = lower(input(prompt));
            end
        else
            fprintf('Background data not found! Please choose normalization mode:\n')
            prompt = sprintf('0. no normalization; 1. Stiched background; 3. Normalization from |S11| fit\n');
            Options.bgdmode = lower(input(prompt));
        end
    elseif Options.bgdmode == 3 % Background normalization through |S11| fitting
        fprintf('Normalization by |S11| fitting...\n')
        % starting value for param = {'kpe', 'kpi', 'omega', 'Gc', 'gma', 'wc', 'attn'}
        param  =  [1e-3   1e-3   freq_h    0   1e-3   f0(tidx)  1/mean(mag(:,tidx))];
        bound_l = [ 0      0     freq_h    0   1e-3    freq_l   0.1]; % lower bound of fitting parameters
        bound_h = [Inf    Inf    freq_h    0   1e-3    freq_h   2.0]; % upper bound of fitting parameters
        fit = iptopt(freq(:,tidx), mag(:,tidx), H0, param, bound_l, bound_h, 1./mag(:,tidx), false);
        bf0 = freq(:,tidx);
        bgd0 = medfilt1(mag(:,tidx)-fit(freq(:,tidx))+1)*fit.attn;
%         bgd0 = medfilt1(mag(:,bidx)./fit(freq(:,bidx)))*fit.attn;
        
        figure
        plot(freq(:,tidx),mag(:,tidx));
        hold on
        xlim([freq_l freq_h])
        plot(bf0,bgd0);
        xlabel('Frequency (GHz)')
        ylabel('S11')
        legend('Zero field raw data','Fitted background', 'Location', 'SouthEast');
        fprintf('Checkpoint: press any key to continue, or Ctrl + c to abort.\n')
        pause
        break
    else
        fprintf('No background normalization\n')
        break
    end
end

% Renormailize the background
if Options.bgdmode ~= 0
    [bf0,trimIdx] = unique(bf0);
    bgd0 = bgd0(trimIdx); % Shift the loaded background to the noise floor of current data
    mag = mag(trimIdx,:); % Remove duplicate points
    freq = freq(trimIdx,:);
    TT = TT(trimIdx,:);
    bgd0 = medfilt1(bgd0);
    bgdM = repmat(bgd0,1,size(mag,2));
    fprintf('Extracting dissipation rates by fitting...\n')
    % starting value for param = {'kpe', 'kpi', 'omega', 'Gc', 'gma', 'wc', 'attn'}
    param  =  [1e-3  1e-3  freq_h    0   1e-3   (freq_l+freq_h)/2  1/mean(bgd0)];
    bound_l = [ 0     0    freq_h    0   1e-3    freq_l   0.8]; % lower bound of fitting parameters
    bound_h = [Inf   Inf   freq_h    0   1e-3    freq_h   1.2]; % upper bound of fitting parameters
    fit = iptopt(freq(:,tidx), mag(:,tidx)./bgd0, H0, param, bound_l, bound_h, bgd0./mag(:,tidx), false);
%     bgd0 = mag(:,bidx)./fit(freq(:,bidx));  % further normalization through fitting
    mag = mag./bgdM;
    analysis.kpe0 = fit.kpe;
    analysis.kpi0 = fit.kpi;
    analysis.wc0 = fit.wc;
    if Options.bgdmode ~= 0
        analysis.bgd0 = bgd0;
        analysis.bf0 = bf0;
    end
end
clearvars mag0l mag0h idx ia trimIdx ii bgd0 bf0

% Plot additional data if there are more than one set of data from VNA
if Options.nData > 1
    S21 = out.data.ZVLreal2 + 1i*out.data.ZVLimag2;
    S21 = S21';
    S21 = S21(:);
    mag2 = abs(S21);
    
    % Plot the interpolated frequency response data in a field scan using color map
    figure
    cmap = pcolor(TT,freq,mag2db(mag2));
    set(cmap, 'edgeColor','none')
    shading interp;
    colorbar
    set(gca,'fontsize',Options.ftsz)
    xlabel('Temperature (K)');
    ylabel('Frequency (GHz)');
    xticks(linspace(field_l,field_h,6));
    title(num2str(Temperature,'S21 response at T = %3.3f K'));
end

%Find all the resonant peaks
f0 = zeros(size(mag,2),1);
Q0 = zeros(size(mag,2),1);
FWHM = zeros(size(mag,2),1);
%find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
for ii = 1:size(mag,2) %Searching column minima (fixed field)
    [~,idx] = min( mag(:,ii) );
    if(length(idx)>1)
        fprintf(num2str(T0(ii),'multiple minima found at H = %.3f'))
    end
    T0(ii) = TT(idx,ii); 
    f0(ii) = freq(idx,ii);
%     mag0(ii) = zq(idx,ii);
    HM = ( max(mag(:,ii)) + min(mag(:,ii)) )/2; % 
    % Calculate quality factor using f0/FWHM
    if isnan(1/range(freq(mag(:,ii) <= HM)))
       Q0(ii) = 0;
    elseif isempty(range(freq(mag(:,ii) <= HM)))
       Q0(ii) = 0;
    else
       FWHM(ii) = range(freq(mag(:,ii) <= HM));
       Q0(ii) = freq(idx,ii)/FWHM(ii);
    end
end

% For noisy data, we need to remove duplicates of minima
[T0,ia,~] = unique(T0,'stable');
f0 = f0(ia);
Q0 = Q0(ia);

f0 = medfilt1(f0); % apply median filter to remove some noise
f0 = f0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
T0 = T0(f0 >= freq_l & f0 <= freq_h);
Q0 = Q0(f0 >= freq_l & f0 <= freq_h);

weight = double.empty(size(mag,1),size(mag,2),0);
for jj = 1:size(mag,2)
    for ii = 1:size(mag,1)
        weight(ii,jj,1) = abs((mag(ii,jj)-max(mag(:,jj)))./mag(ii,jj));
    end
    weight(:,jj) = weight(:,jj)./max(weight(:,jj));
end
switch Options.fitfunc % Pick fitting function from either (1) custom function of (2) spec1d
    case 1 %Option 1: Custom function by Input-output formalism
        % Step 1: Fit the data for the first time to extract "kpe" and "w0" far from the level crossing
        kpe = double.empty(length(T0),0);
        kpi = double.empty(length(T0),0);
        w0 = double.empty(length(T0),0);
        Qf = double.empty(length(T0),0);
        wc = double.empty(length(T0),0);
        xr = double.empty(length(T0),0);
        xi = double.empty(length(T0),0);
        fill = Options.filFctr;
        [~,Tc0] = min(T0);
        f_init = min(f0(Tc0));
        % starting value for param = {'kpe', 'kpi', 'omega', 'Gc', 'gma', 'wc', 'attn'}
        param =  [1e-4   1e-4    0     0    f_init]; % param = {'kpe', 'kpi', 'Re[chi]', 'Im[chi]', 'wc'}
        bound_l = [0      0      0     0    freq_l]; % lower bound of fitting parameters
        bound_h = [Inf   Inf     0     0    freq_h]; % upper bound of fitting parameters
        fit0 = iptoptx(freq(:,Tc0),mag(:,Tc0),H0,fill,param,bound_l,bound_h,weight(:,Tc0),false);
        kpe(Tc0) = fit0.kpe;
        kpi(Tc0) = fit0.kpi;
        w0(Tc0) = fit0.wc;
        f_init = fit0.wc;
        kpe0 = fit0.kpe;
        kpi0 = fit0.kpi;
        parfor ii = Tc0+1:length(T0)
%         for ii = 1:length(T0)
            plt = false;
%             figWin = 10:10:100; % The iteration window over which shows fitting graph
%             if ismember(ii,figWin) % Not useable in parallel mode (parfor)
%                 plt = true;
%             end
            % Fit using input-output formalism
            param =  [ kpe0  kpi0    2e-4   2e-4   f_init]; % param = {'kpe', 'kpi', 'Re[chi]', 'Im[chi]', 'wc'}
            bound_l = [kpe0  kpi0     0      0     f_init]; % lower bound of fitting parameters
            bound_h = [kpe0  kpi0    Inf    Inf    f_init]; % upper bound of fitting parameters
            fit = iptoptx(freq(:,ii),mag(:,ii),H0,fill,param,bound_l,bound_h,weight(:,ii),plt);
            
            if mod(ii,20) == 0
                worker = getCurrentTask();
                fprintf('Fitting, current temperature: %1$3.2f K. Core %2$u.\n', T0(ii), worker.ID);
            end
          
            param = coeffvalues(fit);
            [idx,~,~] = find(fit(freq(:,ii)) == min(fit(freq(:,ii)))); % Find the resonant frequency by (unique) minimum search
            w0(ii) = freq(idx(1)); 
            kpe(ii) = param(1);
            kpi(ii) = param(2);
            xr(ii,1) = param(3);
            xi(ii,1) = param(4);   
            wc(ii) = param(5);
            Qf(ii) = w0(ii)/(param(1)+param(2));
        end
        analysis.kpe = kpe;
        analysis.kpi = kpi;
        analysis.w0 = w0;
        analysis.xr = xr;
        analysis.xi = xi;        
    case 2 %Option 2: use spec1d package to fit the data using Lorentzian form.
        parfor ii = 1:length(T0)
            s = spec1d(freq(:,ii), -mag(:,ii), max(-mag(:,ii)))*0.001; % create spec1d object
            %starting point for the (Lorentzian) fitting parameters
            p = [0.1 wc(ii) FWHM(ii) min(mag(:,ii))]; % (p1: scaling factor ,p2: resonant frequency; p3: FWHM; p4:noise floor(?) )
            fix = [0 0 0 0]; % Denoting if the fitting parameters are fixed
            [~, fbck] = fits(s, 'lorz', p, fix);
%             [~, fbck] = fits(s, 'lorz');
            wc(ii) = fbck.pvals(2); % Retrieve the resonant frequency from fitted data
            Qf(ii) = abs(fbck.pvals(2)/fbck.pvals(3)/2); %Calculate the quality factor
%             chi(ii) = 1/Qf(ii);
            if mod(ii,20) == 0
                worker = getCurrentTask();
                fprintf('Current Temperature: %1$3.2f K. on core %2$u.\n', T0(ii), worker.ID);
            end
        end
end
analysis.Qf = Qf;
analysis.wc = wc;

% Plot the interpolated frequency response data in a field scan using color map

% Plot the interpolated frequency response data in a field scan using color map
figure
cmap = pcolor(TT,freq,mag2db(mag));
hold on
plot(T0,f0,'.k','MarkerSize',Options.mksz);
set(cmap, 'edgeColor','none')
shading interp;
colorbar
set(gca,'fontsize',Options.ftsz)
xlabel('Temperature (K)');
ylabel('Frequency (GHz)');
xticks(linspace(temp_l,temp_h,6));
title(num2str(H0,'Resonant frequency from minimum search at H = %3.3f T'));
legend('Experimental data','Minimum search');

% Plot the interpolated frequency response data in a field scan using color map
figure
cmap1 = copyobj(cmap,gca);
hold on
plot(T0,w0,'.r','MarkerSize',Options.mksz);
set(cmap1, 'edgeColor','none')
shading interp;
colorbar
set(gca,'fontsize',Options.ftsz)
xlabel('Temperature (K)');
ylabel('Frequency (GHz)');
xticks(linspace(temp_l,temp_h,6));
title(num2str(H0,'S11 response at H = %3.3f T'));
legend('Experimental data','|S11| fit');

figure
plot(T0,medfilt1(1./Q0),'.k');
yyaxis right
hold on
[~,idx,Qf] = find(Qf);
plot(TT(1,idx),medfilt1(1./Qf),'.r');
xlabel('Temperature (K)');
ylabel('1/Q');
title(num2str(H0,'Inverse of Q factor at H = %3.3f T'));
legend('f_0 / FWHM','|S11| fitting');

if Options.fitfunc == 1
    figure
    plot(T0,xr,'o');
%     plot(H0,medfilt1(xr),'-');
    ylabel('Susceptibility');
    xlabel('Temperature (K)');
    title('Re[\chi]')
    
    figure
    plot(T0,xi,'o');
%     plot(T0,medfilt1(xi),'-');
    ylabel('Susceptibility');
    xlabel('Temperature (K)');
    title('Im[\chi]')
    
%     figure
%     plot(T0,medfilt1(kpe),'o');
%     ylabel('K_e (GHz)');
%     hold on
%     yyaxis right
%     plot(T0,medfilt1(kpi),'*');
%     xlabel('Temperature (K)');
%     ylabel('K_i (GHz)');
%     title('Dissipation rates');
%     legend('External dissipation rate','Internal dissipation rate')
end
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
    while true
        switch Options.dType
            case 'exp'
                out = readdata_v4(LoadObj,Options.nData);
                freq = out.data.ZVLfreq/1e9;
                S11 = out.data.ZVLreal + 1i*out.data.ZVLimag;
                H = out.data.DCField1;
                HH = repmat(H,1,size(freq,2)); %populate the magnetic field to match the dimension of S11 matrix
                T1 = out.data.Temperature1; % Sample cell temperature
                T1 = repmat(T1,1,size(freq,2));
                [~,analysis.name,~] = fileparts(LoadObj);
                analysis.power = mean(out.data.Power1);
%                T2 = out.data.Temperature2; % Mixing chamber temperature
%                lck1 = out.data.lockin1;
%                lck2 = out.data.lockin2;
                
                S11 = S11';
                S11 = S11(:);
                freq = freq';
                freq = freq(:);
                HH = HH';
                HH = HH(:); %the third argument is the number of frequency points in each line/segment
                T1 = T1';
                T1 = T1(:);
                
                [rows,~,freq] = find(freq); % Remove nonsensical zeros from the frequency data
                HH = HH(rows);
                S11 = S11(rows);
                T1 = T1(rows);
                
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
                T1 = T1(trunc1:trunc2);
                
                % Truncate according to the upper and lower frequency limit
                S11_temp = S11_temp(freq_temp>=freq_l & freq_temp<=freq_h);
                freq_temp = freq_temp(freq_temp>=freq_l & freq_temp<=freq_h);
                HH_temp = HH_temp(freq_temp>=freq_l & freq_temp<=freq_h);
                mag_temp = abs(S11_temp);
                T1 = T1(freq_temp>=freq_l & freq_temp<=freq_h);
                
                freq_l = min(freq_temp); %set frequency range, l: lower limit, h: higher limit
                freq_h = max(freq_temp);               
                dif = diff(freq_temp); % frequency increments
                resets = find( dif <= (freq_l-freq_h) ); % Find the termination point of a complete scans
                nop = round(mean(diff(resets))); % Set the number of points per frequency scan
                
                % Rearrange data to ascending field
                mag = reshape(mag_temp,nop,[]);
                freq = reshape(freq_temp,nop,[]);
                HH = reshape(HH_temp,nop,[]);
                T1 = reshape(T1,nop,[]);
                
                % Determine field/temperature scan direction and sort the data in ascending order
                if Options.analysis == 4
                    if T1(10) > T1
                        analysis.direction = 'up';
                    else
                        analysis.direction = 'down';
                    end
                    [T1,sIdx] = sort(T1(1,:));
                    HH = HH(:,sIdx);
                else
                   if H(10) > H(1)
                       analysis.direction = 'up';
                    else
                       analysis.direction = 'down';
                   end
                   [HH,sIdx] = sort(HH(1,:));
                   T1 = T1(:,sIdx);
                end
                freq = freq(:,sIdx);
                mag = mag(:,sIdx);
                
                analysis.Ts = T1(1,:);
                analysis.Hs = HH(1,:);
                analysis.temp = min(analysis.Ts(analysis.Hs == min(analysis.Hs))); % Measurement temperature at the lowest field (to avoid magnetoresistance in the thermometer)
                break
            case 'sim'
                fprintf('Loading simulated data...\n')
                load(LoadObj,'mfields','freq','S11','analysis');
                HH = mfields;
                mag = abs(S11);
                break
            otherwise
                Options.dType = input('Unknown data type!\n','s');
        end
    end
end
end
% End of data loading function