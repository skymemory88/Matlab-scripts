%% 2D |S11|/|S21| data analysis and plot
%  input: raw data (column format)
%  output: 2D color plot and 1D fits
%  Options: 
%       1: simple plot without analysis.
%       2: On-resonance field scan data with input-output formalism fitting.
%       3: Off-resonance field scan data with input-output formalism fitting.
%       4: temperature scan data with input-output formalism fitting.

format long;
clearvars -except continu_var freq S11 analysis

addpath('G:\.shortcut-targets-by-id\1CapZB_um4grXCRbK6t_9FxyzYQn8ecQE\File sharing\Programming scripts\Matlab\Plot functions\Fieldscan\functions\');
addpath(genpath('G:\.shortcut-targets-by-id\1CapZB_um4grXCRbK6t_9FxyzYQn8ecQE\File sharing\Programming scripts\Matlab\Plot functions\spec1d--Henrik\'));

% Process options:
Options.analysis = 3; % Analysis options (1.simple plot, 2.on-resonance fit, 3.off-resonance fit, 4.temp scan)
% Options.save = 'n'; % Option to save the analysis
Options.dType = 'raw'; % Input data type: 1. Experimental ('raw'), 2. Simulated ('sim'), 3. pre-processed exp. data ('proc')
Options.lnwd = 1.5; % plot linewidth
Options.ftsz = 12; % plot font size
Options.mksz = 2; % plot marker size
Options.bgdmode = 0; % Background normalization (0: none. 1: Stitch background. 2: Load background. 3: S11 fit. 4: Manual offset)
Options.offset = 0.915; % background height (only for bgdmode 4)
Options.nData = 1; % Number of dataset from VNA
Options.phPlot = false; % Option to plot phase in color plots
Options.fitfunc = 1; % Fitting function (only for analysis-3): (1) Input-output or (2) Lorentzian
Options.fill = 0.01005*1.06; % Filling factor * scaling factor (1.06 for SC239, SC200, SC251)
Options.cmap = unique([[0 0 0];[zeros(20,1),linspace(0,0.5,20)',linspace(0.4,1,20)'];...
              [linspace(0,1,36)',0.5*ones(36,1),linspace(1,0,36)'];...
              [ones(30,1),linspace(0.5,1,30)',linspace(0,1,30)']],'rows');
Options.cmap = flip(Options.cmap,1);

Options.fileloc = ['G:\.shortcut-targets-by-id\1CapZB_um4grXCRbK6t_9FxyzYQn8ecQE\',...
    'File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4',...
    '\SC239\2022.09.01'];
% Options.fileloc = ['G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\',...
%     'SC107 (4x5x2mm)\19.05.2019'];
loadname = '2022_09_0002';

% Determin file location base on data type
switch Options.dType
    case 'raw'
        LoadObj = fullfile(Options.fileloc,[loadname,'.dat']);
        SaveObj = fullfile(Options.fileloc,[loadname,'_interp', '.mat']);
    case 'sim'
        Options.fileloc = ['G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Simulations\',...
            'Matlab\Susceptibilities\S11 parameters'];
        loadname = 'S11_LiHoF4_0.200K_3.65GHz_0.02Deg_15.5Deg_1.00e-04';
        LoadObj = fullfile(Options.fileloc,[loadname, '.mat']);
        Options.save = 'n'; % Option to save the analysis
        Options.bgdmode = 0; % no background renormalization for simulated data
    case 'proc'
        LoadObj = fullfile(Options.fileloc,[loadname,'_interp.mat']);
    otherwise
        fprintf(['Unknown data type: ',Options.dType,'\n'])
end

if exist('analysis','var') && exist('freq','var') && exist('S11','var') && exist('continu_var','var')
    [~,fname,~] = fileparts(LoadObj);
    if isfield(analysis, 'name') && strcmp(fname, analysis.name) % check if target data is in the workspace
        arguments =  {analysis continu_var freq S11 LoadObj Options};
    else
        arguments = {LoadObj Options};
    end
else
    arguments = {LoadObj Options};
end
clearvars -except arguments LoadObj Options SaveObj

switch Options.analysis
    case 1
        % Simple color plot of field sweep S11 (& S21)
        [continu_var,freq,S11,analysis] = option1(arguments{:});
    case 2
        % On-resonance field sweep measurement processing (w/ plots)
        [continu_var,freq,S11,analysis] = option2(arguments{:});
    case 3
        % Off-resonance field sweep measurement processing (w/ plots)
        [continu_var,freq,S11,analysis] = option3(arguments{:});
    case 4
        % Temperature scan
        [continu_var,freq,S11,analysis] = option4(LoadObj, Options);
end
analysis.location = Options.fileloc;
clearvars -except temps continu_var freq S11 analysis Options SaveObj

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
            save(SaveObj,'analysis','background','continu_var','freq','S11','-v7.3')
        else
            save(SaveObj,'analysis','continu_var','freq','S11','-v7.3')
        end
%         if Options.analysis == 4 % Add a point to the phase diagram from off-resonance measurement
%             fprintf('Updating phase diagram\n')
%             phase(1,1) = analysis.temp;
%             phase(1,2) = mean(analysis.Hs(analysis.w0==min(analysis.w0)));
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

function [HH,freq,mag,analysis] = option1(varargin)
clearvars -except varargin

[HH,freq,mag,LoadObj,Options,analysis] = load_data(varargin{:});
Temperature = analysis.temp;

freq_l = min(freq(:,1));
freq_h = max(freq(:,1));
field_l = min(HH(1,:));  %set field range, l: lower limit, h: higher limit
field_h = max(HH(1,:));
% field_l = 12.2; % manually set the field limit
% field_h = 17;
if field_l ~= min(HH(1,:)) || field_h ~= max(HH(1,:))
    [~,idxl,~] = find(HH(1,:)>=field_l, 1, 'first');
    [~,idxh,~] = find(HH(1,:)<=field_h, 1, 'last');
    
    HH = HH(:,idxl:idxh); % truncate the data to the selected field range
    mag = mag(:,idxl:idxh);
    freq = freq(:,idxl:idxh);
    
    [idxl,~,~] = find(freq(:,1)>=freq_l);
    [idxh,~,~] = find(freq(:,1)<=freq_h);
    idxl = min(idxl);
    idxh = max(idxh);
    
    HH = HH(idxl:idxh,:); % truncate the data to the selected frequency range
    mag = mag(idxl:idxh,:);
    freq = freq(idxl:idxh,:);
end

%% Background cleanup
clearvars dif step lidx uidx T1 *_temp
if Options.bgdmode ~= 0
    %Find all the resonant peaks
    f0 = zeros(size(mag,2),1);
    B0 = zeros(size(mag,2),1);
    mag0 = zeros(size(mag,2),1);
%     find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
    [~,lidx] = min(HH(1,:));
    [~,uidx] = max(HH(1,:));
    mag0l = min(mag(:,lidx)); % peak depth at the lowest field
    mag0h = min(mag(:,uidx)); % peak depth at the highest field
    [~,hl] = min([mag0l mag0h]); % pick the deeper peak as the first part of the background
    if hl == 2
        [~,bidx] = min( HH(1,:) );
        fprintf('Use low field data for background normalization...\n')
    else
        [~,bidx] = max( HH(1,:) );
        fprintf('Use high field data for background normalization...\n')
    end
    for ii = 1:size(mag,2) %Searching column minima (fixed field)
        [~,idx] = min( mag(:,ii) );
        if(length(idx)>1)
            fprintf(num2str(B0(ii),'multiple minima found at H = %.3f T\n'))
        end
        f0(ii) = freq(idx,ii);
        B0(ii) = HH(idx,ii);
        mag0(ii) = mag(idx,ii);
    end
end
clearvars uidx lidx

while true
    if Options.bgdmode == 1 % Construct the background noise by stitching together zero-field-scan and anti-crossing
        [~,cidx] = min(mag(:,bidx)); % find the resonant peak to be replaced
%         [~,Hpos] = max(mag0); % Method 1. scan of max reflection
%         [~,Hpos] = max(mag(cidx,:)); % Method 2. scan of max reflection at f0 (frequency)
        [~,Hpos] = max(medfilt1(abs(f0-f0(bidx)))); % Method 3. scan with the resonant peak deviates from f0 the furtherest
%         [~,Hpos] = min(abs(B0-3.12)); % Method 4. Manually set the patch slice of the frequency scan
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
        legend(num2str(B0(bidx),'B = %.2f T'),num2str(B0(Hpos),'B = %.2f T'),'Stitched background','location','southeast')
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
        fit = iptopt0(freq(:,bidx), mag(:,bidx), B0(bidx), [param; bound_l; bound_h], 1./mag(:,bidx), false);
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
    fit = iptopt0(freq(:,bidx), mag(:,bidx)./bgd0, B0(bidx), [param; bound_l; bound_h], bgd0./mag(:,bidx), false);
%     bgd0 = mag(:,bidx)./fit(freq(:,bidx));  % further normalization through fitting
    mag = mag./bgdM;
    analysis.kpe0 = fit.kpe;
    analysis.kpi0 = fit.kpi;
    analysis.wc = fit.wc;
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
    colmap = pcolor(HH,freq,mag2db(mag2));
    set(colmap, 'edgeColor','none')
    cmap = unique([[0 0 0];[zeros(20,1),linspace(0,0.5,20)',linspace(0.4,1,20)'];...
                  [linspace(0,1,36)',0.5*ones(36,1),linspace(1,0,36)'];...
                  [ones(30,1),linspace(0.5,1,30)',linspace(0,1,30)']],'rows');
    cmap = flip(cmap,1);
    colormap(cmap)
    caxis([max([1.2*mag2db(min(mag2,[],'all')), -30]) max(mag2db(max(mag2,[],'all')),0)+1]);
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
colmap = pcolor(HH,freq,mag2db(mag));
set(colmap, 'edgeColor','none')
cmap = unique([[0 0 0];[zeros(20,1),linspace(0,0.5,20)',linspace(0.4,1,20)'];...
              [linspace(0,1,36)',0.5*ones(36,1),linspace(1,0,36)'];...
              [ones(30,1),linspace(0.5,1,30)',linspace(0,1,30)']],'rows');
cmap = flip(cmap,1);
colormap(cmap)
caxis([max([1.2*mag2db(min(mag,[],'all')), -30]) max(mag2db(max(mag,[],'all')),0)+1]);
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
    colmap1 = pcolor(HH_temp,freq_temp,imag(S11_temp));
    set(colmap1, 'edgeColor','none')
    shading interp;
    caxis([max([1.2*mag2db(min(imag(S11_temp),[],1)), -30]) max(mag2db(max(imag(S11_temp),[],1)),0)+1]);
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
% freq_l = 3.62;
freq_l = min(freq(:,1));
freq_h = max(freq(:,1));
if freq_l >= max(freq(:,1)) || freq_h <= min(freq(:,1))
    freq_l = min(freq(:,1));
    freq_h = max(freq(:,1));
    fprintf('Manual frequency range out of data range! resort back to default setting.\n')
end

% set input-output fit parameters
ioFit.mode = 1; % fitting function choice (1-6, "0" reserved for empty cavity fitting)
ioFit.tarIdx = 1; % choice of spin resonance to fit (1-7)
init_guess = 0.24; % temperature index of the files for initial guesses
ioFit.init = true; % perform an initial fit to obtain cavity parameters
ioFit.total = 7; % total number of transitions

% set desirable field range
field_l = 0; % manually set the field limit
field_h = 4.6;
% mB0 = 6.27;
strnth = 'strong'; % 'strong' or 'weak'

if field_l >= max(HH(1,:)) || field_h <= min(HH(1,:))
    field_l = min(HH(1,:));
    field_h = max(HH(1,:));
    fprintf('Manual field range out of data range! resort back to default setting.\n')
end

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

mag0l = min(mag(:,HH(1,:) == min(HH(1,:)) )); % peak depth at the lowest field
mag0h = min(mag(:,HH(1,:) == max(HH(1,:)) )); % peak depth at the highest field
[~,hl] = min([mag0l mag0h]); % select the finner peak as the first part of the background
if hl == 2
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
        
        if exist('mB0','var') % Prioritize manual anti-crossing location
            [~,Hp0] = min(abs(HH(1,:)-mB0)); % Method 4. Manually set the patch slice of the frequency scan            
        else
%             [~,Hp0] = max(mag(fidx,:)); % Method 1. scan of max reflection
            [~,Hp0] = max(mag(cidx,:)); % Method 2. scan of max reflection at f0 (frequency)
%             [~,Hp0] = max(medfilt1(abs(f0-f0(bidx)))); % Method 3. scan with the resonant peak deviates from f0 the furtherest
        end

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
%             analysis.wc = temp_load.analysis.wc;            
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
        % Normalization through |S11| fitting
        fprintf('Normalization by |S11| fitting...\n')
        % starting value for param = {'kpe', 'kpi', 'omega', 'Gc', 'gma', 'wc', 'attn'}
        param  =  [1e-3   1e-3   freq_h    0   1e-3    freq(fidx,bidx)   mean(mag(:,bidx))];
        bound_l = [ 0      0     freq_h    0   1e-3    freq_l   0.5]; % lower bound of fitting parameters
        bound_h = [Inf    Inf    freq_h    0   1e-3    freq_h   1.5]; % upper bound of fitting parameters
        fit = iptopt0(freq(:,bidx), mag(:,bidx), HH(1,bidx), [param; bound_l; bound_h], 1./mag(:,bidx), true);
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
        plot(freq(:,bidx),mag(:,bidx)./bgd0)
        break
    else
        fprintf('No background normalization.\n')
        bf0 = freq(:,1);
        bgd0 = ones(size(freq(:,1),1),1);
        break
    end
end

% Remove duplicate points
[bf0,trimIdx] = unique(bf0);
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
B0 = zeros(size(mag,2),1);
mag0 = zeros(size(mag,2),1);
%find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
for ii = 1:size(mag,2) %Searching column minima (fixed field)
    [~,idx] = min( mag(:,ii) );
    if(length(idx)>1)
        fprintf(num2str(B0(ii),'multiple minima found at H = %.3f\n'))
    end
    f0(ii) = freq(idx,ii);
    B0(ii) = HH(idx,ii);
    mag0(ii) = mag(idx,ii);
end 
% further normalization through |S11| fitting
% starting value for param = {'kpe', 'kpi', 'omega', 'Gc', 'gma', 'wc', 'attn'}
param  =  [1e-2   1e-2   freq_h    0   1e-3   f0(bidx)  mean(mag(:,bidx))];
bound_l = [ 0      0     freq_h    0   1e-3    freq_l   0.5]; % lower bound of fitting parameters
bound_h = [Inf    Inf    freq_h    0   1e-3    freq_h   1.5]; % upper bound of fitting parameters
fit = iptopt0(freq(:,bidx), mag(:,bidx), B0(bidx), [param; bound_l; bound_h], 1./mag(:,bidx), false);
bgd0 = mag(:,bidx)./fit(freq(:,bidx))*fit.attn;
% bgd0 = medfilt1(mag(:,bidx)./fit(freq(:,bidx)))*fit.attn;

if Options.bgdmode ~= 0
figure(fig_norm)
hold on
plot(freq(:,bidx),mag(:,bidx)./bgd0);
legend(num2str(B0(bidx),'Normalized with Loaded background at B = %.2f T'),...
    num2str(B0(bidx),'Further normalized with fitted background at B = %.2f T'),'Location', 'SouthEast');
end

prompt = sprintf('Implement additional normalization by fitting? \n');
add_fit = lower(input(prompt,'s'));
if strcmp(add_fit, 'yes') || strcmp(add_fit,'y') || str2double(add_fit) == 1
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
    B0(ii) = HH(idx(1),ii);
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
% B0 = B0(~rmidx);
% mag0 = mag0(~rmidx,:);
% FWHM0 = FWHM0(~rmidx);

analysis.mag0 = mag0;
analysis.FWHM0 = FWHM0;
analysis.Q0 = Q0;

slope = 0.1; % linear spin resonance slope
trigger = 0;
if ioFit.tarIdx > 1
    sgn = +1; % sign of the gradiant of the dispersion relations
else
    sgn = -1;
end

if exist('mB0','var') % Prioritize manual anti-crossing location
    [~,Hp0] = min(abs(B0-mB0)); % Method 3. Manually set position of linecrossing
else
    [~,Hp0] = max(mag0(:,1)); % Option 1: Use the maximum of mag0 to identify anticrossing location
    % [~,Hp0] = min(Q0(:,1)); % Option 2: Use the minimum of Q0 to identify anticrossing location
end

wc = f0(1,1);
gc0 = zeros(1,length(B0));
while true
    switch strnth
        case 'weak'
            % Fit the field dependent resonant frequency data with weak coupling function
%             idx = Hp0; % method 1: use consistent value with the background fitting
%             [~,idx] = max(f0); % method 2: maximal reflection as anticrossing point
%             [~,idx] = findpeaks(f0(:,1)); % method 3: locate the anticrossing by inflection point
%             [~,idx] = min(abs(B0-11.2)); % method 4: manually set the cut-off location
%             [pfit,~,mu] = polyfit(B0(1:idx(1)),f0(1:idx(1),1),1); % use the segment before anticrossing
            [pfit,~,mu] = polyfit(B0,f0(:,1),1); % use the segment over the full field range
            detrend = polyval(pfit,B0,[],mu) - pfit(2);
            init = [1e-2  1e-2   sgn*slope   wc   B0(Hp0)]; % [g  gamma slope wc x0]
            bound_l = [ 0     0   -Inf   freq_l  field_l];
            bound_h = [Inf   Inf   Inf   freq_h  field_h];
            wt = abs(f0(:,1)-f0(1,1))/max(abs(f0(:,1)-f0(1,1)));
%             wt = ones(size(f0));
            [fitP,~] = wk_cpl_fit(B0,f0(:,1)-detrend,init,bound_l,bound_h,wt,true);
%             [fitP,~] = wk_cpl_fit(B0,f0,init,bound_l,bound_h,wt,true);
            analysis.Ffit = fitP;
            slopes = fitP.slope;
            gma0 = fitP.gamma;
            wc = fitP.wc;
            % wc = 3.642;
            gc0(:) = fitP.g;
            if ~exist('mB0','var')
                Brf = fitP.x0; % Level crossing location from perturbative fitting
            else
                Brf = B0(Hp0);
            end
            omg = (slopes + pfit(1)).*(B0-Brf) + wc;
            fprintf('Weak coupling fit complete\n')
            
            Tplot = figure;
            plot(B0,f0(:,1),'ok','MarkerSize',Options.mksz);
            hold on
            plot(B0,fitP(B0) + detrend);
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
                    upH = B0(1:bkr2);
                    loB = f0(bkr1:end,1); % lower branch
                    loH = B0(bkr1:end);
                    Hx1 = linspace(min(upH),max(upH),501); % for trace extrapolation
                    Hx2 = linspace(min(loH),max(loH),length(Hx1));
                else % for dispersion with negative slope
                    upB = cat(1,f0(bkr1:bkr2,2),f0(bkr2+1:end,1)); % upper branch
                    upH = B0(bkr1:end);
                    loB = f0(1:bkr2,1); % lower branch
                    loH = B0(1:bkr2);
                    Hx1 = linspace(min(loH),max(loH),501); % for trace extrapolation
                    Hx2 = linspace(min(upH),max(upH),length(Hx1));
                end
            else % if no overlap found
                if sgn > 0 % for dispersion with positive slope
                    [~, lidx] = max(f0(:,1));
                    [~, ridx] = min(f0(lidx:end,1));
                    upB = f0(1:lidx,1);
                    upH = B0(1:lidx);
                    loB = f0(lidx+ridx-1:end,1);
                    loH = B0(lidx+ridx-1:end);
                else
                    [~, lidx] = min(f0(:,1));
                    [~, ridx] = max(f0(lidx:end,1));
                    loB = f0(1:lidx,1);
                    loH = B0(1:lidx);
                    upB = f0(lidx+ridx-1:end,1);
                    upH = B0(lidx+ridx-1:end);
                end
                Hx1 = linspace(field_l,B0(lidx+ridx-1),501);
                Hx2 = linspace(B0(lidx),field_h,length(Hx1));
            end
            % Plot Resonant frequency trace with fitting parameters
            Tplot = figure;
            plot(upH,upB,'ok',loH,loB,'ok','MarkerSize',Options.mksz);
            % plot(B0(1:round(length(B0)/200):end),f0(1:round(length(f0)/200):end),'ok','MarkerSize',Options.mksz);
            hold on
            init = [0.05   sgn*slope   wc  B0(Hp0)]; % [g slope wc x0]
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
            
            % fitting parameter = ['wc', 'slope', 'x0', 'g']
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
%                fitP = lsqcurvefit(branches, init, [B0(lidx-300:lidx), B0(lidx+uidx-1:end)],...
%                    [f0(lidx-300:lidx), f0(lidx+uidx-1:end)], bound_l, bound_h); % fit the original data

%                 wc = fitP(1);
            slopes = fitP(2);
            Brf = fitP(3);
            gc0(:) = fitP(4);
            omg = slopes.*(B0-Brf) + wc;
            analysis.wfit1 = brch1;
            analysis.wfit2 = brch2;
            detrend = zeros(1,length(B0));
            
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
% [~,Hp0] = min(abs(B0-Brf)); % update anticrossing location with fitted value
init = [max(nonzeros(FWHM0))  gc0(1)  1e-3   mean(slopes)  B0(Hp0)  wc];  % [gamma gc kappa slope x0 wc]
bound_l = [ 0     0     0    min([0 sgn*Inf])   field_l   wc];
bound_h = [Inf   Inf   Inf   max([0 sgn*Inf])   field_h   wc];
% fpts = min([length(B0(1:Hp0))-1 length(B0(Hp0:end))-1 20]);
% [Qfit0, ~] = Qf_fit(B0(Hp0-fpts:Hp0+fpts), Q0(Hp0-fpts:Hp0+fpts), init, bound_l, bound_h, false);
[Qfit0, ~] = Qf_fit(B0, Q0(:,1), init, bound_l, bound_h, false);
analysis.Qfit0 = Qfit0;

%Plot Quality factor from minimum search vs magnetic field
Qplot = figure;
plot(B0, Q0(:,1),'s-','MarkerSize',Options.mksz);
% plot(B0(1:round(length(B0)/20):end), Q0(1:round(length(Q0)/20):end),'s-','MarkerSize',Options.mksz);
hold on
Qp0 = plot(Qfit0,'-k');
set(Qp0,'LineWidth',Options.lnwd);
legend('Exp. data', sprintf('B0 = %1$.2f T\n r = %2$.2e GHz\n Gc = %3$.2e GHz',...
    Qfit0.x0, Qfit0.gamma, Qfit0.gc),'location','SouthWest')

omgQ = Qfit0.slope.*(B0-Qfit0.x0) + Qfit0.wc; % linearized dispersion relation

figure(Tplot)
hold on
plot(B0,omg,'--k');
switch strnth
    case 'weak'
        lgd{end+1} = sprintf('f0 fit, r = %.3f GHz', fitP.gamma);
    case 'strong'
        lgd{end+1} = sprintf('f0 fit, gc = %.3f GHz', gc0(1));
end
if exist('omgQ','var')
    plot(B0,omgQ,'-.r');
    lgd{end+1} = sprintf('Q0 fit, r = %.3f GHz', Qfit0.gamma);
end

% load and erxtrapolate spin resonant frequencies as auxiliary field
if ioFit.mode > 0
    aux_files = fit_list(init_guess);
    aux_param = spinRes_load(Options.fileloc, aux_files, B0);
    omg_temp = double.empty(length(B0),0);
    for ii = 1:length(B0)
        omg_temp(ii,1) = aux_param(ii).f0(ioFit.tarIdx); % initial guess of spin resonant frequencies
        analysis.auxW(ii,:) = aux_param(ii).f0(:); % auxiliary resonant frequencies
        analysis.auxG(ii,:) = aux_param(ii).gc(:); % coupling strength of auxiliary resonant frequencies
        analysis.auxR(ii,:) = aux_param(ii).gma(:); % spin linewidth of auxiliary resonant frequencies
    end
    analysis.w0_init = omg_temp;
    plot(B0,omg_temp,'-.b');
    lgd{end+1} = sprintf('Previous fit');
    legend(lgd,'location', 'southwest')
else
    if exist('detrend','var')
        aux_param = detrend;
    else
        aux_param = zeros(1,length(B0));
    end
end


prompt = 'Choose initial guess for spin resonance (1. Perturbation fit, 2. Q-factor fit, 3. Previous fit): ';
init_opt = input(prompt);
if init_opt == 1
    gma0 = 2e-2; % arbitrary guess of spin linewidth (GHz)
    analysis.w0_init = omg;
%     for ii = 1:length(B0)
%         gc0(ii) = aux_param(ii).gc(ioFit.tarIdx); % manually set initial values of coupling strength
%     end
    gc0(:) = 1e-2; % [GHz]
elseif init_opt == 2
    gma0 = Qfit0.gamma; 
    omg = omgQ;
    analysis.w0_init = omgQ;
    [~,Hp0] = min(abs(B0-Qfit0.x0));
    gc0(:) = Qfit0.gc;
    wc = Qfit0.wc;
elseif init_opt == 3 && ioFit.mode > 0
    omg = omg_temp;
    gma0 = aux_param(bidx).gma(ioFit.tarIdx);
    for ii = 1:length(B0)
        gc0(ii) = aux_param(ii).gc(ioFit.tarIdx);
        analysis.Gci(ii) = aux_param(ii).gc(ioFit.tarIdx);
    end
else
    % makes no change
end
clearvars B Delt hPara fPara wp wm fitPara H_res f_res Qparam omg_temp

figure
hold on
box on
cmap0 = pcolor(HH,freq,mag2db(mag));
set(cmap0, 'edgeColor','none')
shading interp;
colorbar
% caxis('auto')
caxis([max([1.2*mag2db(min(mag,[],'all')), -30]) max(mag2db(max(mag,[],'all')),0)+1]);
clearvars init bound_l bound_h *_temp

set(gca,'fontsize',Options.ftsz)
axis([field_l field_h freq_l freq_h]);
xlabel('Magnetic Field (T)');
ylabel('Frequency (GHz)');
title(num2str(Temperature,'S11 response at T = %3.3f K'));
        
% Plot frequency scan at line crossing vs. zero-field scan with fittings
Cplot = figure;
% [~,bidx] = min(B0);
plot(freq(1:10:end,bidx),mag(1:10:end,bidx),'o');
hold on
plot(freq(1:10:end,Hp0),mag(1:10:end,Hp0),'or');
xlabel('Frequency (GHz)');
ylabel('S11');
clgd = {sprintf('Frequency cut at %.2f T',min(B0) ), sprintf('Frequency cut at %.2f T',B0(Hp0))};

%% Step 1: Fit the data for the first time to extract "kpe" and "w0" far from the level crossing
ff0 = double.empty(length(B0),2,0); 
Qf = double.empty(length(B0),0);
Gc = double.empty(length(B0),0);
Gc_ci = zeros(length(B0),2);
Gc1 = double.empty(length(B0),0);
w1 = double.empty(length(B0),0);
w2 = double.empty(length(B0),0);

xr = double.empty(length(B0),0);
xi = double.empty(length(B0),0);
gamma = double.empty(length(B0),0);
gamma_ci = zeros(length(B0),2);
omg_ci = zeros(length(B0),2);
attn = double.empty(length(B0),0);
fits = double.empty(size(mag,1),size(mag,2),0);

weight(:,:,1) = abs(1./mag(:,:));      % Weight function option 1
% weight = 1./gradientweight(mag);       % Weight function option 2
% weight(isinf(weight)) = 1e4;           % Remove infinities in weight function
% weight = ones(size(mag));                % Weight function option 3 (uniform)
param  =  [omg(bidx)   gc0(1)    gma0   FWHM0(1)/10    FWHM0(1)/10    wc   1/mean(mag(:,bidx))];  % {omega, Gc, gma, kpe, kpi, wc, attn}
bound_l = [omg(bidx)   gc0(1)    gma0    0     0    freq_l   0.8]; % lower bound of fitting parameters
bound_h = [omg(bidx)   gc0(1)    gma0   Inf   Inf   freq_h   1.2]; % upper bound of fitting parameters
fit0 = iptopt(freq(:,bidx), mag(:,bidx), [0 0 0 1 B0(bidx) false], [param; bound_l; bound_h], [], weight(:,bidx), ioFit);
ioFit.init = false;

% analysis.kpe = fit0.kpe;
% analysis.kpi = fit0.kpi;
conf_intv = confint(fit0,0.95);
omg(bidx) = fit0.omega;
omg_ci(bidx,:) = conf_intv(:,1);
attn(bidx,1) = fit0.attn;
gamma(1,1) = gma0;
fits(:,bidx,1) = fit0(freq(:,bidx));

[~,idx] = findpeaks(-fit0(freq(:,bidx))); % Find the resonant frequencies
if ~isempty(idx)
    if length(idx) == 2
        ff0(bidx,1,1) = freq(idx(1),bidx);  
        ff0(bidx,2,1) = freq(idx(2),bidx);
    elseif length(idx) > 2
        [~,idxt] = sort(mag(:,bidx));
        ff0(bidx,1,1) = freq(idxt(1),bidx);
        ff0(bidx,2,1) = freq(idxt(2),bidx);
    else
        ff0(bidx,:,1) = freq(idx,bidx);
    end
else
    ff0(bidx,:,1) = fit0.wc;
end
Qf(bidx) = mean(ff0(bidx,:,1)./(fit0.kpe + fit0.kpi));
figure(Cplot)
plot(freq(:,bidx),fit0(freq(:,bidx)),'-k');
% clgd{end+1} = '';

%% Step 2: fit the data for the second time with fixed "kpe" and "w0"
% Construct fixed coefficients for the second fit
ff0_1 = ff0(:,1);
ff0_2 = ff0(:,2);

coef = [fit0.kpe  fit0.kpi  fit0.wc  fit0.attn]; % using the fit to data away from anti-crossings for cavity paramters
% coef = [analysis.kpe0  analysis.kpi0  fit0.wc  fit0.attn]; % load cavity parameters
plt = zeros(1,length(B0)); % fit inspection option
% plt(Hp0-10:3:Hp0+10) = 1; % Selectively plot the fit for visual inspection
% for ii = 1:length(B0) % for debugging
parfor ii = 1:length(B0)
    if mod(ii,20) == 0
        worker = getCurrentTask();
        fprintf('Core %2$u: Fitting, current field: %1$3.2f T. \n', B0(ii), worker.ID);
    end
    param  =  [omg(ii)  gc0(ii)   gma0     0       0    gc0(ii)   gc0(ii)]; % {'omega', 'Gc', 'gma', 'xr', 'xi', 'Gc1', 'Gc2'}
    bound_l = [ -Inf       0       0     -Inf      0       0        0  ]; % lower bound of fitting parameters
    bound_h = [  Inf      Inf     Inf     Inf      0      Inf      Inf ]; % upper bound of fitting parameters
    fit = iptopt(freq(:,ii), mag(:,ii), [coef B0(ii) plt(ii)], [param; bound_l; bound_h], aux_param(ii), weight(:,ii), ioFit);
    
    conf_intv = confint(fit,0.95);
    omg(ii) = fit.omega;
    omg_ci(ii,:) = conf_intv(:,1);  
    Gc(ii) = fit.Gc;
    Gc_ci(ii,:) = conf_intv(:,2);
    if any(strcmp(fieldnames(fit),'Gc1')); Gc1(ii) = fit.Gc1; end
    if any(strcmp(fieldnames(fit),'w1')); w1(ii) = fit.w1; end
    if any(strcmp(fieldnames(fit),'w2')); w2(ii) = fit.w2; end
    if any(strcmp(fieldnames(fit),'xr')); xr(ii,1) = fit.xr; end
    if any(strcmp(fieldnames(fit),'xi')); xi(ii,1) = fit.xi; end
    gamma(ii,1) = fit.gma;
    gamma_ci(ii,:) = conf_intv(:,3);
    attn(ii,1) = fit.attn;
    fits(:,ii) = fit(freq(:,ii));
    
    switch strnth
        case 'weak' % Anticipate no branching in weak coupling regime
            [~,idx] = min(fit(freq(:,ii)));
            ff0_1(ii) = freq(idx,ii);
            ff0_2(ii) = freq(idx,ii);
        case 'strong'
            [~,idx] = findpeaks(-fit(freq(:,ii))); % Find the resonant frequency with peak finder
            if ~isempty(idx)
                if length(idx) == 2
                    ff0_1(ii) = freq(idx(1),ii);
                    ff0_2(ii) = freq(idx(2),ii);
                elseif length(idx) > 2
                    [~,idxt] = min(mag(idx,ii));
                    ff0_1(ii) = freq(idx(idxt),ii);
                    ff0_2(ii) = freq(idx(idxt),ii);
                else
                    ff0_1(ii) = freq(idx,ii);
                    ff0_2(ii) = freq(idx,ii);
                end
            else
                ff0_1(ii) = min(fit(freq(:,ii)));  % alternatively by minimum search
                ff0_2(ii) = min(fit(freq(:,ii)));
            end
    end
    qf1 = ff0_1(ii)./(fit0.kpe + fit0.kpi);
    qf2 = ff0_2(ii)./(fit0.kpe + fit0.kpi);
    Qf(ii) = mean([qf1 qf2]); % Calculate quality factor
end
ff0(:,1) = ff0_1(:);
ff0(:,2) = ff0_2(:);
figure(Cplot)
plot(freq(:,Hp0),fits(:,Hp0),'-k');
clgd{end+1} = sprintf('Fitting: G_c = %1$.2e GHz, r = %2$.2e GHz.',Gc(Hp0),gamma(Hp0,1));
legend(clgd,'location','SouthEast')

% mu0 = 4*pi*10^-7; % Vacuum permeability ([H/m])
% hbar = 1.055E-34; % Reduced Planck constant [J.s]
% meV2J = 1.602217e-22; % [J/meV]
% rho = 4e30/(5.17*5.17*10.75); % Holmium (magnetic moment) number density [m^-3]
% f2E = hbar*2*pi*10^9/meV2J; % [meV/GHz]
% g0 = sqrt(mu0*wc*10^9*2*pi*rho*Options.fill /hbar/2); % susceptibility prefactor [T.(J.s)^-1]
% g0 = g0 * meV2J * f2E * 10^-9;

if ~isempty(xr); analysis.xr = xr; end % contribution from off-resonance spin resonance
if ~isempty(xi); analysis.xi = xi; end % contribution from off-resonance spin resonance
analysis.B0 = B0;
analysis.wr = ff0; % CMP resonance frequencies
analysis.w0 = omg; % spin resonance frequencies
analysis.w0_ci = omg_ci;
analysis.gma = gamma; % spin linewidth
analysis.gma_ci = gamma_ci;
analysis.Qf = Qf; % CMP quality factor
% analysis.Qfit = Qfit;
analysis.Gc = Gc; % on-resonance coupling strength
analysis.Gc_ci = Gc_ci;
if ~isempty(Gc1); analysis.aGc = Gc1; end % on-resonance coupling strength
if ~isempty(w1); analysis.aW1 = w1; end % on-resonance coupling strength
if ~isempty(w2); analysis.aW2 = w2; end % on-resonance coupling strength
analysis.attn = attn; % attenuation rate
analysis.fits = fits; % fits to individual |S11| slice

% Plot the resonant frequency from Lorentzian fit versus DC magnetic field
figure
cmap2 = copyobj(cmap0, gca);
set(cmap2, 'edgeColor','none');
shading interp;
colorbar
% caxis([max([mag2db(min(mag0)), -30]) max(mag2db(max(mag0)),0)+1]);
caxis('auto');
hold on
plot(B0(bidx:end),ff0(bidx:end,:),'or','MarkerSize',Options.mksz,'MarkerFaceColor','red');
plot(B0,omg,'.k');
plot(B0,omg_ci(:,1),'-r',B0,omg_ci(:,2),'-r');
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
% % hfig1 = plot(B0(1:round(length(B0)/200):end),f0(1:round(length(f0)/200):end),'ok','MarkerSize',Options.mksz,'MarkerFaceColor','black');
% plot(B0, f0, 'ok', 'MarkerSize', Options.mksz); % Plot the resonant frequency from minimum search
% xlabel('Magnetic Field (T)');
% ylabel('Resonant frequency (GHz)');
% title(num2str(Temperature,'Resonant frequency from minimum search at T = %3.3f K'));
% % caxis('auto');
% axis([field_l field_h freq_l freq_h]);

% Plot the peak amplitude from minimum search vs. magnetic field
figure
mag0 = medfilt1(mag0); % apply median filter to remove some noise
plot(B0, mag2db(mag0), 'o', 'MarkerSize', Options.mksz);
xlabel('Field(T)');
ylabel('S11 (dB)');
title(num2str(Temperature,'Minimal S11 at T = %3.3f K'));

% Plot Quality factor from fitting function vs magnetic field
% init = [max(gamma)  mean(gc0)  kpe0   mean(slopes)  B0(Hp0)  wc];  % [gamma gc kappa slope x0 wc]
% bound_l = [0   0   0   -Inf   field_l   freq_l];
% bound_h = [1   1   1    Inf   field_h   freq_h];
% [Qfit, ~] = Qf_fit(B0(Hp0-15:Hp0+15), Qf(Hp0-15:Hp0+15), init, bound_l, bound_h, false);
figure(Qplot)
% yyaxis right
% hold on
% B0 = B0(Qf >=0);
% Qf = Qf(Qf >=0); %Remove unphysical points
% Qf = medfilt1(Qf);
% plot(B0, Qf,'ro-','MarkerSize',Options.mksz);
% % plot(B0(Hc2:round(length(B0)/20):Hc2), Qf(Hc2:round(length(Qf)/20):Hc2),'o-','MarkerSize',Options.mksz);
% % plot(B0, Qf,'o','MarkerSize',Options.mksz);
% Qp = plot(Qfit,'-k');
% set(Qp,'LineWidth',Options.lnwd);
% gca;
xlabel('Magnetic Field (T)');
ylabel('Q factor');
% legend('Quality factor from FWHM', 'Q0 fit', 'Quality factor from fitting', 'Qf fit');
title(num2str(Temperature,'Quality factor, T= %.2f K'));

figure;
errorbar(B0,Gc,abs(Gc_ci(:,1)-Gc_ci(:,2)),'o');
ylim([0 1])
lg = {'Gc'};
xlabel('Magnetic Field (T)')
ylabel('Coupling strength (GHz)')
title(num2str(Temperature,'Fitting paramters from input-output formalism at T = %.2f K'))
hold on
yyaxis right
% plot(B0(Hc1:end), medfilt1(gamma(Hc1:end)),'s')
errorbar(B0, gamma,abs(gamma_ci(:,1)-gamma_ci(:,2)),'s')
% errorbar(B0(Hc1:end), gamma(Hc1:end),abs(gamma_ci(Hc1:end,1)-gamma_ci(Hc1:end,2)),'s')
lg{end+1} = '\gamma';
ylabel('Spin linewidth (GHz)')
legend(lg)
ylim([0 1])
end
% End of option 2 (on-resonance fitting)
function [HH,freq,mag,analysis] = option3(varargin)
%% Set data range and parameters
[HH,freq,mag,LoadObj,Options,analysis] = load_data(varargin{:});
Temperature = analysis.temp;

% field_l = 2.0; % Manually set field range
% field_h = 5.0;
% mB0 = 4.786;
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
B0 = zeros(size(mag,2),1);
f0 = zeros(size(mag,2),1);
mag0 = zeros(size(mag,2),1);
%find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
mag0l = min(mag(:,HH(1,:)==min(HH(1,:)) )); % peak depth at the lowest field
mag0h = min(mag(:,HH(1,:) == max(HH(1,:)) )); % peak depth at the highest field
[~,hl] = min([mag0l mag0h]); % pick the finner peak as the first part of the background
if hl == 2
    [~,bidx] = min( HH(1,:) );
else
    [~,bidx] = max( HH(1,:) );
end
for ii = 1:size(mag,2) %Searching column minima (fixed field)
    [~,idx] = min( mag(:,ii) );
    if(length(idx)>1)
        fprintf(num2str(B0(ii),'multiple minima found at H = %.3f\n'))
    end
    f0(ii) = freq(idx,ii);
    B0(ii) = HH(idx,ii);
    mag0(ii) = mag(idx,ii);
end
[~,fidx] = min(mag(:,bidx));

while true
    if Options.bgdmode == 1
        fprintf('Normalizing background by stiching...\n')
        [~,cidx] = min(mag(:,bidx)); % find the resonant peak to be replaced
        
        if exist('mB0','var') % Manually set the anti-crossing location
            [~,Hp0] = min(abs(HH(1,:)-mB0)); % Method 4. Manually set the patch slice of the frequency scan            
        else
%             [~,Hp0] = max(mag(fidx,:)); % Method 1. scan of max reflection
%             [~,Hp0] = max(mag(cidx,:)); % Method 2. scan of max reflection at f0 (frequency)
            [~,Hp0] = max((abs(f0-f0(bidx)))); % Method 3. scan with the resonant peak deviates from f0 the furtherest
        end

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
        legend(sprintf('B = %.2f T', HH(1,bidx)), sprintf('B = %.2f T',HH(1,Hp0)),'Stitched','location','SouthEast')
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
        % Normalization through |S11| fitting
        fprintf('Normalization by |S11| fitting...\n')
        % starting value for param = {'kpe', 'kpi', 'omega', 'Gc', 'gma', 'wc', 'attn'}
        param  =  [1e-3   1e-3   freq_h    0   1e-3    freq(fidx,bidx)   mean(mag(:,bidx))];
        bound_l = [ 0      0     freq_h    0   1e-3    freq_l   0.5]; % lower bound of fitting parameters
        bound_h = [Inf    Inf    freq_h    0   1e-3    freq_h   1.5]; % upper bound of fitting parameters
        fit = iptopt0(freq(:,bidx), mag(:,bidx), HH(1,bidx), [param; bound_l; bound_h], 1./mag(:,bidx), false);
        bf0 = freq(:,bidx);
%         bgd0 = medfilt1(mag(:,bidx)-fit(freq(:,bidx))+1)*fit.attn;
        bgd0 = medfilt1(fit(freq(:,bidx))./mag(:,bidx))*fit.attn;
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
        plot(freq(:,bidx),mag(:,bidx)./bgd0)
        break
    else
        fprintf('No background normalization.\n')
        bf0 = freq(:,1);
        bgd0 = ones(size(freq(:,1),1),1);
        break
    end
end

% Renormailize the background
if Options.bgdmode ==(1||2||3)
    figure(fig_norm)
    xlabel('Frequency (GHz)')
    ylabel('|S11|')
    title(sprintf('Background normalized data at %.2f T', HH(1,bidx)))

    [bf0,trimIdx] = unique(bf0);
    bgd0 = bgd0(trimIdx); % Shift the loaded background to the noise floor of current data
    mag = mag(trimIdx,:); % Remove duplicate points
    freq = freq(trimIdx,:);
    HH = HH(trimIdx,:);
    bgd0 = medfilt1(bgd0);
    bgdM = repmat(bgd0,1,size(mag,2));
    fprintf('Extracting dissipation rates by fitting...\n')
    % starting value for param = {'kpe', 'kpi', 'omega', 'Gc', 'gma', 'wc', 'attn'}
    param  =  [1e-4   1e-4   freq_h    0   1e-3   (freq_l+freq_h)/2  1/mean(bgd0)];
    bound_l = [ 0      0     freq_h    0   1e-3    freq_l   0.8]; % lower bound of fitting parameters
    bound_h = [Inf    Inf    freq_h    0   1e-3    freq_h   1.2]; % upper bound of fitting parameters
    fit = iptopt0(freq(:,bidx), mag(:,bidx)./bgd0, B0(bidx), [param; bound_l; bound_h], mag(:,bidx)./bgd0, false);
%     bgd0 = mag(:,bidx)./fit(freq(:,bidx));  % further normalization through fitting
    mag = mag./bgdM;
    analysis.kpe0 = fit.kpe;
    analysis.kpi0 = fit.kpi;
    analysis.wc = fit.wc;
    if Options.bgdmode ~= 0
        analysis.bgd0 = bgd0;
        analysis.bf0 = bf0;
    end
elseif Options.bgdmode == 4
    mag = mag./Options.offset;
end
clearvars mag0l mag0h idx ia trimIdx ii bgd0 bf0

% Find all the resonant peaks
f0 = zeros(size(mag,2),1);
Q0 = zeros(size(mag,2),1);
% mag0 = zeros(size(zq,2),1);
FWHM = zeros(size(mag,2),1);
%find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
for ii = 1:size(mag,2) %Searching column minima (fixed field)
    [~,idx] = min( mag(:,ii) );
    if(length(idx)>1)
        fprintf(num2str(B0(ii),'multiple minima found at H = %.3f'))
    end
    B0(ii) = HH(idx,ii); 
    f0(ii) = freq(idx,ii);
%     mag0(ii) = zq(idx,ii);
    HM = ( max(mag(:,ii)) + min(mag(:,ii)) )*0.7; % equivalent to Full-Width-Half-Max in power spectrum
%     HM = ( max(mag(:,ii)) + min(mag(:,ii)) )/2; % Full-Width-Half-Max
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
B0 = B0(~c); % Remove corresponding elements in B0 array as well
f0 = f0(~c);
% mag0 = mag0(~c);
FWHM = FWHM(~c);
HH = HH(:,~c);
freq = freq(:,~c);
mag = mag(:,~c);
% analysis.temp = analysis.temp(~c);

% For noisy data, we need to remove duplicates of minima
[B0,ia,~] = unique(B0,'stable');
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
B0 = B0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
Q0 = Q0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
% mag0 = mag0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
FWHM = FWHM(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
HH = HH(:,f0 >= freq_l & f0 <= freq_h);
freq = freq(:,f0 >= freq_l & f0 <= freq_h);
mag = mag(:,f0 >= freq_l & f0 <= freq_h);
% analysis.temp = analysis.temp(f0 >= freq_l & f0 <= freq_h);
clearvars c idx ia ii HM trunc1 trunc2 dupl nop trimIdx

f_init = f0(B0 == min(B0)); % Initial guess of the cavity resonant frequency
weight(:,:,1) = abs(1./mag(:,:));      % Weight function option 1
% weight = 1./gradientweight(mag);       % Weight function option 2
% weight = ones(size(mag));                % Weight function option 3 (uniform)
% weight = double.empty(size(mag,1),size(mag,2),0); % Weight function option 4
% for jj = 1:size(mag,2)
%     for ii = 1:size(mag,1)
%         weight(ii,jj,1) = abs((mag(ii,jj)-max(mag(:,jj)))./mag(ii,jj));
%     end
%     weight(:,jj) = weight(:,jj)./max(weight(:,jj));
% end
weight(isinf(weight)) = 1e4;           % Remove infinities in weight function

switch Options.fitfunc % Pick fitting function from either (1) custom function of (2) spec1d
    case 1 %Option 1: Custom function by Input-output formalism
        % Step 1: Fit the data for the first time to extract "kpe" and "w0" far from the level crossing
        kpe = double.empty(length(B0),0);
        kpi = double.empty(length(B0),0);
        w0 = double.empty(length(B0),0);
        Qf = double.empty(length(B0),0);
        wc = double.empty(length(B0),0);
        [~,Hc0] = min(B0);
%         Hc0 = bidx;
        % Fit the zero-field data for cavity properties
        param  =  [1e-4    1e-4    0     0    f_init]; % param = {'kpe', 'kpi', 'Re[chi]', 'Im[chi]', 'wc'}
        bound_l = [ 0       0      0     0    freq_l]; % lower bound of fitting parameters
        bound_h = [Inf     Inf     0     0    freq_h]; % upper bound of fitting parameters
        fit0 = iptoptx(freq(:,Hc0), mag(:,Hc0), B0(Hc0), param, bound_l, bound_h, weight(:,Hc0), true);
        kpe(Hc0) = fit0.kpe;
        kpi(Hc0) = fit0.kpi;
        w0(Hc0) = fit0.wc;
        kpe0 = fit0.kpe;
        kpi0 = fit0.kpi;
        mu0 = 4*pi*1e-7; % Vacuum permeability ([H/m])
        meV2J = 1.60218e-22; % [J/meV]
        rho = 4 / (5.175e-10 * 5.175e-10 * 10.75e-10); % Holmium (magnetic moment) number density [m^-3]
        g = sqrt(mu0 * 2*pi * f_init*1e9 * rho/2) * Options.fill;  % susceptibility prefactor [T^2/J. rad/s]^1/2
        g2 = g^2  * 2*pi * 1e-9; % coupling prefactor [T^2/J. GHz]
%         analysis.f_init = f_init; % for debugging

%         for ii = 1:length(B0)
            plt = false;
%             figWin = Hpos-100:20:Hpos+100; % The iteration window over which shows fitting graph
% %             figWin = length(B0);
%             if ismember(ii,figWin) % Not useable in parallel mode (parfor)
%                 plt = true;
%             end
        parfor ii = Hc0+1:length(B0)
            if mod(ii,20) == 0
                worker = getCurrentTask();
                fprintf('Fitting, current field: %1$3.2f T. Core %2$u.\n', B0(ii), worker.ID);
            end
            % Fit using input-output formalism
            param  =  [kpe0   kpi0     0      0    f_init]; % param = {'kpe', 'kpi', 'Re[chi]', 'Im[chi]', 'wc'}
            bound_l = [kpe0   kpi0   -Inf   -Inf  f_init]; % lower bound of fitting parameters
            bound_h = [kpe0   kpi0    Inf    Inf   f_init]; % upper bound of fitting parameters
            fit = iptoptx(freq(:,ii), mag(:,ii), B0(ii), param, bound_l, bound_h, weight(:,ii), plt);
                     
            [idx,~,~] = find(fit(freq(:,ii)) == min(fit(freq(:,ii)))); % Find the resonant frequency by (unique) minimum search
            w0(ii) = freq(idx(1)); 
            kpe(ii) = fit.kpe;
            kpi(ii) = fit.kpi;
            xr(ii) = fit.xr / meV2J / g2; % [meV/T^2]
            xi(ii) = fit.xi / meV2J / g2; % [meV/T^2]  
            wc(ii) = fit.wc;
            Qf(ii) = w0(ii)/(fit.kpe + fit.kpi);
        end
        analysis.kpe = kpe;
        analysis.kpi = kpi;
        analysis.w0 = w0;
        analysis.xr = xr;
        analysis.xi = xi;        
    case 2 %Option 2: use spec1d package to fit the data using Lorentzian form.
        parfor ii = 1:length(B0)
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
                fprintf('Current magnetic field: %1$3.2f. on core %2$u.\n', B0(ii), worker.ID);
            end
        end
end
analysis.Qf = Qf;
analysis.wc = wc;

% Plot the interpolated frequency response data in a field scan using color map
figure
colmap = pcolor(HH,freq,mag2db(mag));
set(colmap, 'edgeColor','none')
colormap(Options.cmap)
caxis([max([1.2*mag2db(min(mag,[],'all')), -30]) max(mag2db(max(mag,[],'all')),0)+1]);
shading interp;
colorbar
set(gca,'fontsize',Options.ftsz)
xlabel('Magnetic Field (T)');
ylabel('Frequency (GHz)');
xlim([field_l,field_h]);
xticks(linspace(field_l,field_h,6));
title(num2str(Temperature,'S11 response at T = %3.3f K'));
legend('Exp. data','Minimum search');

figure
colmap1 = copyobj(colmap,gca);
hold on
plot(B0,f0,'.k','MarkerSize',Options.mksz);
colormap(Options.cmap)
set(colmap1, 'edgeColor','none')
shading interp;
colorbar
set(gca,'fontsize',Options.ftsz)
xlabel('Magnetic Field (T)');
ylabel('Frequency (GHz)');
xlim([field_l,field_h]);
xticks(linspace(field_l,field_h,6));
ylim([freq_l,freq_h]);
title(num2str(Temperature,'S11 response at T = %3.3f K'));
legend('Exp. data','Minimum search');

% Plot the interpolated frequency response data in a field scan using color map
figure
cmap2 = copyobj(colmap,gca);
hold on
plot(B0,w0,'.r','MarkerSize',Options.mksz);
colormap(Options.cmap)
set(cmap2, 'edgeColor','none')
shading interp;
colorbar
set(gca,'fontsize',Options.ftsz)
xlabel('Magnetic Field (T)');
ylabel('Frequency (GHz)');
xlim([field_l,field_h]);
xticks(linspace(field_l,field_h,6));
ylim([freq_l,freq_h]);
title(num2str(Temperature,'S11 response at T = %3.3f K'));
legend('Exp. data','|S11| fit');

figure
plot(B0,f0,'ok','MarkerSize',Options.mksz);
hold on
plot(mean(B0(f0==min(medfilt1(f0)))), mean(f0(f0==min(f0))), 'o', 'MarkerFaceColor', 'red','MarkerSize',Options.mksz);
xlim([field_l,field_h]);
xlabel('Magnetic Field (T)');
ylabel('Resonant frequency (GHz)');
title(num2str(Temperature,'Resonant frequency from minimum search at T = %3.3f K'));

figure
plot(B0,medfilt1(1./Q0),'.k');
yyaxis right
hold on
[~,idx,Qf] = find(Qf);
plot(HH(1,idx),medfilt1(1./Qf),'.r');
xlim([field_l,field_h]);
xlabel('Magnetic Field (T)');
ylabel('1/Q');
title(num2str(Temperature,'Inverse of Q factor at T = %3.3f K'));
legend('f_0 / FWHM','|S11| fitting');

if Options.fitfunc == 1
    figure
    plot(B0,medfilt1(xr),'-');
    ylabel('Susceptibility');
    xlabel('Magnetic Field (T)');
    title('Re[\chi]')
    
    figure
    plot(B0,medfilt1(xi),'-');
    ylabel('Susceptibility');
    xlabel('Magnetic Field (T)');
    title('Im[\chi]')
    
%     figure
%     plot(B0,medfilt1(kpe),'o');
%     ylabel('K_e (GHz)');
%     hold on
%     yyaxis right
%     plot(B0,medfilt1(kpi),'*');
%     xlabel('Magnetic Field (T)');
%     ylabel('K_i (GHz)');
%     title('Dissipation rates');
%     legend('External dissipation rate','Internal dissipation rate')
end
end
% End of option 3 (off-resonance fitting)
function [TT,freq,mag,analysis] = option4(varargin)
% extract data from raw data file
[~,freq,mag,LoadObj,Options,analysis] = load_data(varargin{:});
T1 = analysis.Ts;
TT = repmat(T1,size(freq,1),1);
B0 = mean(analysis.Hs);

freq_l = min(freq(:,1));
freq_h = max(freq(:,1));
temp_l = min(T1);  % set field range, l: lower limit, h: higher limit
temp_h = max(T1);

ioFit.mode = 1; % input-output fitting mode (1-6, "0" reserved for empty cavity fitting)
ioFit.tarIdx = 4; % choice of spin resonance to fit (0-6)
ioFit.init = true; % perform an initial fit to obtain cavity parameters
ioFit.total = 1;
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
    if Options.bgdmode ~= 0; fprintf('Low field data for background normalization...\n'); end
else
    [~,tidx] = max( TT(1,:) );
    if Options.bgdmode ~= 0; fprintf('Use high field data for background normalization...\n'); end
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
%         [~,Hpos] = min(abs(B0-14.2)); % Method 4. Manually set the patch slice of the frequency scan
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
        fit = iptopt0(freq(:,tidx), mag(:,tidx), B0, [param; bound_l; bound_h], 1./mag(:,tidx), false);
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
    [bf0,trimIdx] = unique(bf0); % Remove duplicate points
    bgd0 = bgd0(trimIdx);
    mag = mag(trimIdx,:);
    freq = freq(trimIdx,:);
    TT = TT(trimIdx,:);
    bgd0 = medfilt1(bgd0);
    bgdM = repmat(bgd0,1,size(mag,2));
%     bgd0 = mag(:,bidx)./fit(freq(:,bidx));  % further normalization through fitting
    mag = mag./bgdM;
    analysis.bgd0 = bgd0;
    analysis.bf0 = bf0;
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
    colmap = pcolor(TT,freq,mag2db(mag2));
    set(colmap, 'edgeColor','none')
    colormap(Options.cmap)
    caxis([max([1.2*mag2db(min(mag2,[],'all')), -30]) max(mag2db(max(mag2,[],'all')),0)+1]);
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

f0 = medfilt1(f0); % apply median filter to remove some noise
f0 = f0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
T0 = T0(f0 >= freq_l & f0 <= freq_h);
Q0 = Q0(f0 >= freq_l & f0 <= freq_h);
analysis.Q0 = Q0;

weight = double.empty(size(mag,1),size(mag,2),0);
for jj = 1:size(mag,2)
    for ii = 1:size(mag,1)
        weight(ii,jj,1) = abs((mag(ii,jj)-max(mag(:,jj)))./mag(ii,jj));
    end
    weight(:,jj) = weight(:,jj)./max(weight(:,jj));
end
weight(isinf(weight)) = 1;
weight(isnan(weight)) = 0;
switch Options.fitfunc % Pick fitting function from either (1) custom function of (2) spec1d
    case 1 %Option 1: Custom function by Input-output formalism
        % Step 1: First fit to extract "kpe", "kpi" and "w0" far from the anti-crossing
        fprintf('Extracting dissipation rates by fitting...\n')
        % Fit using input-output formalism {'omega', 'Gc', 'gma', 'attn', 'kpe', 'kpi', 'wc'}
        param   = [freq_l    0    1e-4   1/mean(mag(:,tidx))    1e-4    1e-4    (freq_l+freq_h)/2];
        bound_l = [freq_l    0    1e-4   0.8      0      0     freq_l]; % lower bound of fitting parameters
        bound_h = [freq_l    0    1e-4   1.2     Inf    Inf    freq_h]; % upper bound of fitting parameters
        fit0 = iptopt(freq(:,tidx), mag(:,tidx), [0 0 0 T0(tidx) 1 false], [param; bound_l; bound_h], [], ones(size(mag(:,tidx))), ioFit);
        ioFit.init = false; % reset and free the fitting function choice
        analysis.kpe = fit0.kpe;
        analysis.kpi = fit0.kpi;
        analysis.wc = fit0.wc;
        
        % Construct fixed coefficients for the second fit
        mu0 = 4*pi*10^-7; % Vacuum permeability ([H/m])
        hbar = 1.055E-34; % Reduced Planck constant [J.s]
        meV2J = 1.602217e-22; % [J/meV]
        rho = 4e30/(5.17*5.17*10.75); % Holmium (magnetic moment) number density [m^-3]
        f2E = hbar*2*pi*10^9/meV2J; % [meV/GHz]
        g0 = sqrt(mu0*fit0.wc*10^9*2*pi*rho*Options.fill /hbar/2); % susceptibility prefactor [T.(J.s)^-1]
        g0 = g0 * meV2J * f2E * 10^-9;
        coef = [fit0.kpe  fit0.kpi  fit0.wc  fit0.attn];
        % Step 2: Fit the rest of the data using "kpe", "kpi" and "w0" from initial fit
        wr = double.empty(length(T0),0);
        Qf = double.empty(length(T0),0);
        omg = double.empty(length(T0),0);
        xr = double.empty(length(T0),0);
        xi = double.empty(length(T0),0);
        Gc = double.empty(length(T0),0);
        gma = double.empty(length(T0),0);
        plt = zeros(length(T0),1);
%         plt(10:80:length(T0)) = 1;
%         for ii = 1:length(T0)
        parfor ii = 1:length(T0)
            if mod(ii,20) == 0
                worker = getCurrentTask();
                fprintf('Fitting, current temperature: %1$3.2f K. Core %2$u.\n', T0(ii), worker.ID);
            end
            % Fit using input-output formalism {'omega', 'Gc', 'gma', 'attn', 'xr', 'xi', 'Gc1'}
            param  =  [freq_l   1e-3   1e-4   1/mean(mag(:,ii))     0     0     0];
            bound_l = [-Inf      0      0      0.8      0     0    0]; % lower bound of fitting parameters
            bound_h = [ Inf     Inf    Inf     1.2     Inf   Inf   0]; % upper bound of fitting parameters
            fit = iptopt(freq(:,ii), mag(:,ii), [coef T0(ii) plt(ii)], [param; bound_l; bound_h], [], weight(:,ii), ioFit);

          
            [idx,~,~] = find(fit(freq(:,ii)) == min(fit(freq(:,ii)))); % Find the resonant frequency by (unique) minimum search
            wr(ii) = freq(idx(1)); 
            if any(strcmp(fieldnames(fit),'xr')); xr(ii,1) = fit.xr; end
            if any(strcmp(fieldnames(fit),'xi')); xi(ii,1) = fit.xi; end
            omg(ii) = fit.omega;   
            Gc(ii) = fit.Gc;
            gma(ii) = fit.gma; 
        end
        analysis.xr = xr;
        analysis.xi = xi;  
        analysis.g = Gc;
        analysis.gma = gma;
        analysis.wr = wr;
        analysis.Qf(:) = wr(:)./(coef(1) + coef(2) + g0^2.*xr(:));
        analysis.w0 = omg;
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
% Plot the interpolated frequency response data in a field scan using color map
figure
colmap = pcolor(TT,freq,mag2db(mag));
caxis([max([1.2*mag2db(min(mag,[],'all')), -30]) max(mag2db(max(mag,[],'all')),0)+1]);
hold on
plot(T0,f0,'.k','MarkerSize',Options.mksz);
set(colmap, 'edgeColor','none')
colormap(Options.cmap)
shading interp;
colorbar
set(gca,'fontsize',Options.ftsz)
xlabel('Temperature (K)');
ylabel('Frequency (GHz)');
xticks(linspace(temp_l,temp_h,6));
title(num2str(B0,'S11 response and resonant frequqncies at H = %3.3f T'));
legend('Experimental data','Minimum search');

% Plot the interpolated frequency response data in a field scan using color map
figure
colmap1 = copyobj(colmap,gca);
hold on
plot(T0,wr,'.r','MarkerSize',Options.mksz);
set(colmap1, 'edgeColor','none')
colormap(Options.cmap)
shading interp;
colorbar
set(gca,'fontsize',Options.ftsz)
xlabel('Temperature (K)');
ylabel('Frequency (GHz)');
xlim([temp_l temp_h]);
xticks(linspace(temp_l,temp_h,6));
title(num2str(B0,'S11 response and resonant frequencies at H = %3.3f T'));
legend('Experimental data','|S11| fit');

figure
plot(T0,medfilt1(1./Q0),'.k');
yyaxis right
hold on
[~,idx,Qf] = find(Qf);
plot(TT(1,idx),medfilt1(1./Qf),'.r');
xlabel('Temperature (K)');
ylabel('1/Q');
title(num2str(B0,'Inverse of Q factor at H = %3.3f T'));
legend('f_0 / FWHM','|S11| fitting');

if Options.fitfunc == 1
    figure
    plot(T0,xr,'o');
%     plot(B0,medfilt1(xr),'-');
    ylabel('Susceptibility');
    xlabel('Temperature (K)');
    title('Re[\chi_{off}]')
    
    figure
    plot(T0,xi,'o');
%     plot(T0,medfilt1(xi),'-');
    ylabel('Susceptibility');
    xlabel('Temperature (K)');
    title('Im[\chi_{off}]')
    
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
            case 'raw'
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
                
                % For data taken with ZVL/ZNL/ZNB 
                % Clean up the raw data by removing incomplete scans (step 1) and duplicates(step 2)
                % step 1: truncate the beginning and end part to keep only complete frequency scans and reshape the matrices into single column vectors
                trunc1 = find(freq == freq_l,1,'first');
                trunc2 = find(freq == freq_h,1,'last');
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
                             
%                 %% for data taken with PNA-X
%                 dif = diff(freq); % frequency increments
%                 resets = find(dif<=0.9*(freq_l-freq_h)); % Find the termination point of a complete scans (10% error)
%                 nop = round(mean(diff(resets))); % Set the number of points per frequency scan
%                 mag = reshape(abs(S11./max(S11)),nop,[]);  %reshape the matrix so that each complete frequency scan occupy one column
%                 freq = reshape(freq,nop,[]);
%                 HH = reshape(HH,nop,[]);
%                 T1 = reshape(T1,nop,[]);

                % Determine field/temperature scan direction and sort the data in ascending order
                if Options.analysis == 4
                    if T1(1,10) > T1
                        analysis.direction = 'up';
                    else
                        analysis.direction = 'down';
                    end
                    [~,sIdx] = sort(T1(1,:));
                    T1 = T1(:,sIdx);
                    HH = HH(:,sIdx);
                else
                   if H(10) > H(1)
                       analysis.direction = 'up';
                    else
                       analysis.direction = 'down';
                   end
                   [~,sIdx] = sort(HH(1,:));
                   T1 = T1(:,sIdx);
                   HH = HH(:,sIdx);
                end
                freq = freq(:,sIdx);
                mag = mag(:,sIdx);
                
                analysis.Ts = T1(1,:);
                analysis.Hs = HH(1,:);
                analysis.temp = min(analysis.Ts(analysis.Hs == min(analysis.Hs))); % Measurement temperature at the lowest field (to avoid magnetoresistance in the thermometer)
                break
            case 'sim'
                fprintf('Loading simulated data...\n')
                load(LoadObj,'continu_var','freq','S11','analysis');
                HH = continu_var;
                mag = abs(S11);
                break
            case 'proc'
                load(LoadObj,'continu_var', 'freq', 'S11', 'analysis');
                HH = continu_var;
                mag = S11;
                break
            otherwise
                Options.dType = input('Unknown data type!\n','s');
        end
    end
end
end
% End of data loading function
function aux_param = spinRes_load(varargin)
% input argument: location,  files, query B field (optional)
% output: resonant frequencies, coupling strength, spin linewidth
file_loc = varargin{1};
file_names = varargin{2};
aux_param = struct;
for ii = 1:length(file_names)
    if strlength(file_names(ii)) ~= 0
        aux_file = fullfile(file_loc, strcat(file_names(ii),".mat"));
        load(aux_file,'H0','w0','gc','gma');
        if nargin == 3
            B0 = varargin{3};
        elseif nargin == 2
            B0 = H0;
        end
        [w0, rmidx] = rmoutliers(w0); % remove outliers (default 3 MAD)
        H0 = H0(~rmidx);
        gc = gc(~rmidx);
        gma = gma(~rmidx);
        [pfit,~,mu] = polyfit(H0,w0,1);
        f0 = polyval(pfit,B0,[],mu);
        Gc = interp1(H0,gc,B0,'nearest','extrap');
        gma = interp1(H0,gma,B0,'nearest','extrap');
        for kk = 1:length(B0)
            aux_param(kk).f0(ii) = f0(kk);
            aux_param(kk).gc(ii) = Gc(kk);
            aux_param(kk).gma(ii) = gma(kk);
%             aux_param(kk).gc(ii) = mean(Gc);
%             aux_param(kk).gma(ii) = mean(gma);
        end
    else
        for kk = 1:length(B0)
            aux_param(kk).f0(ii) = 0;
            aux_param(kk).gc(ii) = 0;
            aux_param(kk).gma(ii) = 1e-5;
        end
    end
end
end