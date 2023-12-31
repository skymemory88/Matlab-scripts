%% Set data range and parameters
location = '';
filename = '';
loadobj = fullfile(location,loadname);
Options.dType = 'exp';
Options.nData = 1;
Options.bgdmode = 0; % background normalization (0: no normalization, 1: stitched background, 2: loaded background)
argmnt = {loadobj, Options};

[HH,freq,mag,~,Options,analysis] = load_data(argmnt{:});
Temperature = analysis.temp;

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
    plot(analysis.Hs,analysis.Ts,'o-')
    xlabel('DC Magnetic field')
    ylabel('Temperature')
    title('Magnetic field vs Temperature')
end
%% Clean up raw data
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
    fprintf('Press any key to continue, or Ctrl + C to abort.\n')
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
ff0 = f0(2); % Initial guess of the cavity resonant frequency
H0 = H0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
Q0 = Q0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
% mag0 = mag0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
FWHM = FWHM(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
HH = HH(:,f0 >= freq_l & f0 <= freq_h);
freq = freq(:,f0 >= freq_l & f0 <= freq_h);
mag = mag(:,f0 >= freq_l & f0 <= freq_h);
% analysis.temp = analysis.temp(f0 >= freq_l & f0 <= freq_h);
clearvars idx ia ii dupl nop trimIdx

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
        kpe = double.empty(length(H0),0); % external dissipations rate
        kpi = double.empty(length(H0),0); % internal dissipations rate
        w0 = double.empty(length(H0),0); % Spin resonance
        Qf = double.empty(length(H0),0); % Quality factor
        xr = double.empty(length(H0),0); % Re[chi]
        xi = double.empty(length(H0),0);% Im[chi]
        for ii = 1:length(H0)
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
            fit = iptoptx(freq(:,ii),mag(:,ii),H0(ii),param,bound_l,bound_h,weight(:,ii),plt);
            
            if mod(ii,20) == 0
                worker = getCurrentTask();
                fprintf('Fitting, current field: %1$3.2f T. Core %2$u.\n', H0(ii), worker.ID);
            end
          
            param = coeffvalues(fit);
            [idx,~,~] = find(fit(freq(:,ii)) == min(fit(freq(:,ii)))); % Find the resonant frequency by (unique) minimum search
            w0(ii) = freq(idx(1)); 
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
    case 2 %Option 2: use lorentizan form to fit the data
    return
end

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
    
    figure
    plot(H0,medfilt1(kpe),'o');
    ylabel('K_e (GHz)');
    hold on
    yyaxis right
    plot(H0,medfilt1(kpi),'*');
    xlabel('Magnetic Field (T)');
    ylabel('K_i (GHz)');
    title('Dissipation rates');
    legend('External dissipation rate','Internal dissipation rate')
end
end