function Field_scan_v5
    format long;  
%     addpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\Fieldscan\functions\');
%     addpath(genpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\spec1d--Henrik\'));
    addpath('/Volumes/GoogleDrive/My Drive/File sharing/Programming scripts/Matlab/Plot functions/Fieldscan/functions')
    addpath(genpath('/Volumes/GoogleDrive/My Drive/File sharing/Programming scripts/Matlab/Plot functions/spec1d--Henrik'));
    %The first line is for windows, the second line is for mac OS

    % Figure plot options:
    Options.plot = false; % Plot option in analysis part (Option 3)
    Options.lnwd = 1.5;
    Options.ftsz = 12;
    Options.mksz = 3;
    Options.fitfunc = 1; % Pick fitting function from either (1) custom function of (2) spec1d
    
%     loadpath = 'G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\SC200\2021.02.05';
    loadpath = '/Volumes/GoogleDrive/My Drive/File sharing/PhD program/Research projects/LiHoF4 project/Data/Experiment/LiHoF4/SC200/2021.02.10';
    %The first line is for windows, the second line is for mac OS
    loadname = '2021_02_0013.dat';
    opt = 2;% Analysis options
    nZVL = 1; % Number of dataset from ZVL
    fileobj = fullfile(loadpath,loadname);

    % Path and file name to save
    savepath = 'G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\SC200';
    savename = 'Phase Diagram.mat';
    dataobj = fullfile(savepath,savename);

    switch opt
        case 1
            % Simple color plot of S11 (+ S21)
            option1(fileobj, Options, nZVL);
        case 2
            % On-resonance measurement processing (w/ plots)
            option2(fileobj, Options, nZVL);
        case 3
            % Data fitting and file saving (w/o plots)
            option3(fileobj, nZVL);
        case 4
            % Off-resonance measurement processing
            option4(fileobj, dataobj, Options, nZVL);
        case 5
            % Temperature scan
            option5(fileobj, Options, nZVL);
    end
end

function option1(fileobj, Options, nZVL)

out = readdata_v4(fileobj,nZVL);
freq = out.data.ZVLfreq/1e9;
S11 = out.data.ZVLreal + 1i*out.data.ZVLimag;
H = out.data.DCField1;
HH = repmat(H,1,size(freq,2)); %populate the magnetic field to match the dimension of S11 matrix
T1 = out.data.Temperature1; % Sample cell temperature
% T2 = out.data.Temperature2; % Mixing chamber temperature
lck1 = out.data.lockin1;
lck2 = out.data.lockin2;

Temperature = mean(T1(H == min(H)));

S11 = S11';
S11 = S11(:);
freq = freq';
freq = freq(:);
HH = HH';
HH = HH(:); %the third argument is the number of frequency points in each line/segment

[rows,~,freq] = find(freq); % Remove nonsensical zeros from the frequency data
HH = HH(rows);
S11 = S11(rows);
dB = mag2db(abs(S11));

%Set data range and parameters
freq_l = min(freq); %set frequency range, l: lower limit, h: higher limit
freq_h = max(freq);
field_l = min(H);  %set field range, l: lower limit, h: higher limit
field_h = max(H);

% Could use "scatteredInterpolant()" to replace "TriScatteredInterp()" as recommended by MATLAB, but it may generate artifacts
FdB = TriScatteredInterp(HH,freq,dB);
FrS = TriScatteredInterp(HH,freq,real(S11)); %intrapolate the real part of S11 response on a 2D grid
FiS = TriScatteredInterp(HH,freq,imag(S11)); %intrapolate the imaginary part of S11 response on a 2D grid

% Plot frequency-field colour map
[xq,yq] = meshgrid(linspace(field_l,field_h,1000),linspace(freq_l,freq_h,1000)); %set the X and Y range
zq = FdB(xq,yq);
zq1 = FrS(xq,yq)+1i*FiS(xq,yq);

% Plot additional data if there are more than one set of data from VNA
if nZVL > 1
    S21 = out.data.ZVLreal2 + 1i*out.data.ZVLimag2;
    S21 = S21';
    S21 = S21(:);
    dB2 = mag2db(abs(S21));
    FdB2 = TriScatteredInterp(HH,freq,dB2);
    zq2 = FdB2(xq,yq);
    
    % Plot the interpolated frequency response data in a field scan using color map
    figure
    cmap = pcolor(xq,yq,zq2);
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
cmap = pcolor(xq,yq,zq);
set(cmap, 'edgeColor','none')
shading interp;
% caxis([-30 0])
colorbar
set(gca,'fontsize',Options.ftsz)
xlabel('Field (T)');
ylabel('Frequency (GHz)');
xticks(linspace(field_l,field_h,6));
title(num2str(Temperature,'S11 response at T = %3.3f K'));

% Plot the interpolated imaginary part of the response function in a field scan using color map
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
function [out] = option2(fileobj, Options, nZVL)
%Set data range and parameters
order = 4; % set to what order the median filters is applied
clear freq S11 dB N FdB FrS FiS FTT1 FTT2

% extract data from raw data file
if ~exist('out','var')
    out = readdata_v4(fileobj, nZVL);
end

freq = out.data.ZVLfreq/1e9;
S11 = out.data.ZVLreal + 1i*out.data.ZVLimag;  
H = out.data.DCField1;
T1 = out.data.Temperature1;
% T2 = out.data.Temperature2;

Temperature = mean(T1(H == min(H))); %Extract the measurement temperature at the lowest field (to avoid magnetoresistance in the thermometer)
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
field_l = min(H);  % set field range, l: lower limit
field_h = max(H);  % set field range, h: lower limit
% field_l = 5;  % Manually set field range, l: lower limit
% field_h = 6;  % Manually set field range, h: lower limit
Hcut0 = field_l; % Field window for cavity parameter fit
Hcut1 = 2; % Field window for cavity parameter fit
Hcut2 = 2; % Field window for line-crossing fit
Hcut3 = field_h; % Field window for line-crossing fit
% Interpolate the data on a 2D grid for the colormap
[xq,yq] = meshgrid(linspace(field_l,field_h,501),linspace(freq_l,freq_h,801));

%Plot the temperature vs magnetic field to check the temperature variation
figure
plot(H(1:round(length(T1)/100):end),T1(1:round(length(T1)/100):end),'o-')
xlabel('DC Magnetic field')
ylabel('Temperature')
title('Magnetic field vs Temperature')

%% Code For R&S ZVL-6
% Clean up the raw data by removing incomplete scans (step 1) and duplicates(step 2)
% step 1: truncate the beginning and end part to keep only complete frequency scans and reshape the matrices into single column vectors
trunc1 = find(freq==freq_l,1,'first'); 
trunc2 = find(freq==freq_l,1,'last')-1; 
trunc3 = find(HH>=field_l,1,'first'); % Truncate the data according to the set field range
trunc4 = find(HH<=field_h,1,'last'); 
S11_temp = S11(max(trunc1,trunc3):min(trunc2,trunc4));
freq_temp = freq(max(trunc1,trunc3):min(trunc2,trunc4));
HH_temp = HH(max(trunc1,trunc3):min(trunc2,trunc4));

% % Step 2: remove duplicates (Not for data since 2019)
% dupl = find(diff(freq_temp) == 0.0);
% freq_temp(dupl+1)=[];
% S11_temp(dupl+1)=[];
% HH_temp(dupl+1)=[];

dif = diff(freq); % frequency increments
resets = find(dif<=0.9*(freq_l-freq_h)); % Find the termination point of a complete scans (10% error)
nop = round(mean(diff(resets))); % Set the number of points per frequency scan

% % Alternative way of finding number of points per complete frequency scan
% dif = nonzeros(diff(freq)); % Extract from raw data the step size in frequency scan
% dif = dif(dif>0); % Keep only positive steps
% step = min (dif); % Calculate the value of the frequency scan step
% nop = ceil(abs(freq_h-freq_l)/step); %compute how many points pers complete frequency scan.

mag_temp = abs(S11_temp);
dB_temp = mag2db(mag_temp);
mag = reshape(mag_temp,nop,[]);
dB = reshape(dB_temp,nop,[]);  
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
clearvars dif step

%Find all the resonant peaks
f0 = zeros(size(dB,2),1);
H0 = zeros(size(dB,2),1);
Q0 = zeros(size(dB,2),1);
dB0 = zeros(size(dB,2),1);
FWHM0 = zeros(size(dB,2),1);
%find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
[~,bidx] = min( HH(1,:) ); % find the column index for zero-field frequency scan
for ii = 1:size(dB,2) %Searching column minima (fixed field)
    [~,idx] = min( dB(:,ii) );
    if(length(idx)>1)
        disp(num2str(H0(ii),'multiple minima found at H = %.3f'))
    end
    f0(ii) = freq(idx,ii);
    H0(ii) = HH(idx,ii); 
    dB0(ii) = dB(idx,ii);
%     HM = dB(idx,ii)*0.3;
    HM = -abs(dB(idx,ii)-max(dB(:,ii)))*0.3;
    % Calculate quality factor using f0/FWHM
    if isnan(1/range(freq(dB(:,ii) <= HM)))
       Q0(ii) = 0;
    elseif isempty(range(freq(dB(:,ii) <= HM)))
       Q0(ii) = 0;
    else
       Q0(ii) = freq(idx,ii)/range(freq(dB(:,ii) <= HM));
       FWHM0(ii) = range(freq(dB(:,ii) <= HM));
    end
    if ii == bidx
        lidx = find(dB(1:idx,ii) >= HM,1,'last'); % left stitching point for background
        ridx = idx+ find(dB(idx:end,ii) >= HM,1,'first'); % right stitching point for background
        widx = ridx-lidx;
    end
end
Q0(isinf(Q0)) = NaN; % Cut out inf from the array
[Q0,c] = rmmissing(Q0); % Cut out NaN from the array
H0 = H0(~c); % Remove corresponding elements in H0 array as well
f0 = f0(~c);
dB0 = dB0(~c);
FWHM0 = FWHM0(~c);

% For noisy data, we need to remove duplicates of minima
[H0,ia,~] = unique(H0,'stable');
f0 = f0(ia);
Q0 = Q0(ia);
dB0 = dB0(ia);
FWHM0 = FWHM0(ia);

f0 = medfilt1(f0,order); % apply median filter to remove some noise
f0 = f0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
H0 = H0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
Q0 = Q0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
dB0 = dB0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
FWHM0 = FWHM0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints

[~,Hpos] = max(dB0); % find the line crossing position on field axis
% Hpos = 61; % Manually set the line crossing position (index)

% Fit the field dependent resonant frequency data with weak coupling function
spin = 8; % Electronic spin
hPara = [H0(Hpos), field_l, field_h];
fPara = [(freq_l+freq_h)/2, freq_l, freq_h];
[fitP,~] = wk_cpl_fit(H0,f0,spin,hPara,fPara);
% Brf = H0(Hpos); % Level crossing location from minimum search
Brf = fitP.x0; % Level crossing location from perturbative fitting
gc = fitP.g;

% Resonant frequency trace with fitting parameters
B = linspace(field_l,field_h,100);
Delt = -spin*(B-Brf); % Assume linear relation near the line cross
plot(H0(1:round(length(H0)/200):end),f0(1:round(length(f0)/200):end),'ok','MarkerSize',4);
hold on
wp = fitP.wc + Delt./2 + sqrt(Delt.^2+4*fitP.g^2)/2;
wm = fitP.wc + Delt./2 - sqrt(Delt.^2+4*fitP.g^2)/2;
plot(B,wm,'-r',B,wp,'-r','LineWidth',Options.lnwd);
hold on
plot(B,Delt+fitP.wc,'--k');
% hfig1 = plot(H0, f0, 'o', 'MarkerSize', 2);
xlabel('Field (T)');
ylabel('Resonant frequency (GHz)');
title(num2str(Temperature,'Resonant frequency from minimum search at T = %3.3f K'));
legend(sprintf('gc = %.3f GHz',fitP.g/2));
axis([field_l field_h freq_l freq_h]);

% Plot frequency scan at line crossing
figure
% plot(freq(1:10:end,Hpos),dB(1:10:end,Hpos),'-o');
plot(freq(:,Hpos),dB(:,Hpos),'-o');
xlabel('Frequency (GHz)');
ylabel('S11 (dB)');
legend(sprintf('Frequency cut at %.2f T',H0(Hpos)));
title('Frequency scan at line crossing');
clearvars B Delt hPara fPara wp wm spin fitPara H_res f_res

% Construct the background noise by stitching together zero-field-scan and anti-crossing
bgd0 = mag(:,bidx); % zero-field frequency scan
bgd1 = mag(lidx-3*widx:ridx+3*widx,Hpos); % center segment of frequency scan slightly away from the anti-crossing
bgd0(lidx-3*widx:ridx+3*widx) = bgd1; % substitute the center segment of zero-field frequency scan with that from anti-crossing
figure
plot(freq(:,bidx),mag(:,bidx))
hold on
plot(freq(:,Hpos),mag(:,Hpos))
plot(freq(:,bidx),bgd0);
xlabel('Frequency (GHz)')
ylabel('S11 (dB)')
legend('B = 0','B = Br','Stitched')
[bf0,trimIdx] = unique(freq(:,bidx));
bgd0 = interp1(bf0,bgd0(trimIdx),yq(:,1));
dB0 = interp1(H0,dB0,xq(1,:));
[~,Hpos] = max(dB0); % find the line crossing position on field axis along interpolated data
clearvars c idx ia ii HM T1 trunc1 trunc2 dupl nop bf0

switch 2 %choose data interpolation method and plot the color map of S11 response
    case 1 % Option_1 Interpolate data along only the frequency axis.
        interp_dB = zeros(size(yq,1),size(HH,2));
        parfor ii = 1:size(HH,2)
            interp_dB(:,ii) = interp1(freq(:,ii),dB(:,ii),yq(:,1));
        end
        zq = interp_dB'; % For later use in lorentzian fitting
        HH = repmat(H0,1,size(xq,1));
        
        switch 1 % Choose plotting option
            case 1
                %     plot single-direction interpolated data using pseudo-colormap
                cmap0 = pcolor(HH(:,1),yq(:,1),interp_dB);
                set(cmap0, 'edgeColor','none')
                shading interp;
                colorbar
            case 2
                %     plot single-direction interpolated data using scatter plot
                yy = repmat(yq(:,1),1,length(H0));
                yy = yy(:);
                interp_dB = interp_dB(:);
                HH = HH';
                HH = HH(:);
                cmap0 = scatter3(HH,yy,interp_dB,2,interp_dB,'o','filled','MarkerEdgeColor','none');
                colormap(hsv);
        end
    case 2 % Option_2 Interpolate the data along both axis.
        %    FdB  = TriScatteredInterp(HH,freq,dB);
        FdB = scatteredInterpolant(HH_temp, freq_temp, mag_temp); % 2D interpolation of amplitude data
%         FdB = scatteredInterpolant(HH_temp, freq_temp, dB_temp); % 2D interpolation of return loss data
        %     FrS  = TriScatteredInterp(HH,freq,real(S11));
        %     FiS  = TriScatteredInterp(HH,freq,imag(S11));
        %     FTT1 = TriScatteredInterp(HH,freq,TT1);
        %     FTT2 = TriScatteredInterp(HH,freq,TT2);
        figure
        hold on
        box on               
        zq = FdB(xq,yq);
        
        % Normalize the data to background noise
        bgdM = repmat(bgd0,1,size(zq,2));
        zq = zq - bgdM + 1;
%         zq = zq - 0.3; % Manually adjust the floor
        
        cmap0 = pcolor(xq,yq,mag2db(zq));
        set(cmap0, 'edgeColor','none')
        shading interp;
        colorbar
%         caxis('auto')
        caxis([max([min(dB0), -30]) max(max(dB0),0)+1]);
end
% clearvars *_temp;
        
set(gca,'fontsize',Options.ftsz)
axis([field_l field_h freq_l freq_h]);
xlabel('Field (T)');
ylabel('Frequency (GHz)');
title(num2str(Temperature,'S11 response at T = %3.3f K'));

%Fit each frequency scan to Lorentzian form to extract the quality factor
%     figure
disp(num2str(Temperature,'Fitting: T = %3.3f K'))

switch 2 % Use interpolated data or raw data to extract quality factor
    case 1 % Use raw data
        Hx = H0;
        yq = freq;
        zq = mag;
        ff0 = f0;
        FWHM = FWHM0;
    case 2 % Use interpolated data
        Hx = xq(1,:);
        ff0 = interp1(H0,f0,Hx,'pchip','extrap'); % Using spline interpolation to smooth the resonant frequency trace
        FWHM = interp1(H0,FWHM0,Hx,'pchip','extrap');
end

Qf = double.empty(length(Hx),0);
gamma = double.empty(length(Hx),0);
Gc = double.empty(length(Hx),0);
% B0 = double.empty(length(Hx),0);
zq(isnan(zq))=0;
weight = double.empty(size(zq,1),size(zq,2),0);
for jj = 1:length(Hx)
    for ii = 1:length(zq(:,jj))
        weight(ii,jj,1) = abs(zq(ii,jj)-max(zq(:,jj)));
    end
end
switch Options.fitfunc % Pick fitting function from either (1) custom function of (2) spec1d
    case 1 %Option 1: Custom function by Input-output formalism
        % Step 1: Fit the data for the first time to extract "kpe" and "w0" far from the level crossing
        Hc0 = find(Hx >= Hcut0, 1,'first');
        Hc1 = find(Hx <= Hcut1, 1, 'last');
        kpe = double.empty(Hc1-Hc0,0);
        kpi = double.empty(Hc1-Hc0,0);
        w0 = double.empty(Hc1-Hc0,0);
%         B0 = double.empty(Hc1-Hc0,0);
        parfor ii = Hc0:Hc1
            plt = false;
%             % Selectively plot the fit for visual inspection
%             figWin = [1 10 50 Hc1]; % The iteration window over which shows fitting graph
%             if ismember(ii,figWin) % Not useable in parallel mode (parfor)
%                 plt = true;
%             end
            % Fit using input-output formalism
            param = [FWHM(ii) FWHM(ii)/10  ff0(ii)  0   1e3]; % starting value for param = {'kpe', 'kpi', 'w0', 'Gc', 'gma'}
%             Set up boundaries for the fitting parameters
            bound_l = [ 0   0   0   0  1]; % lower bound of fitting parameters
            bound_h = [inf inf inf  0  1]; % upper bound of fitting parameters
            fit = iptopt_0(yq(:,ii),zq(:,ii),Hx(ii),Brf,param,bound_l,bound_h,weight(:,ii),plt);
 
            if mod(ii,20) == 0
                worker = getCurrentTask();
                fprintf('First fit, current magnetic field: %1$3.2f. on core %2$u.\n', Hx(ii), worker.ID);
            end

            param = coeffvalues(fit);
            kpe(ii) = param(1);
            kpi(ii) = param(2);
            w0(ii) = param(3);
        end
%% Step 2: fit the data for the second time with fixed "kpe" and "w0"
        kpe0 = mean(medfilt1(kpe));
%         kemx = max(medfilt1(kpe));
%         kemn = min(medfilt1(kpe));
        kpi0 = mean(medfilt1(kpi));
%         kimx = max(medfilt1(kpi));
%         kimn = min(medfilt1(kpi));
        w0 = w0(end);
        Hc2 = find(Hx >= Hcut2, 1,'first');
        Hc3 = find(Hx <= Hcut3, 1, 'last');
        for ii = Hc2:Hc3
            % Selectively plot the fit for visual inspection
            plt = false;
            figWin = Hpos-30:6:Hpos+30; % The iteration window over which shows fitting graph
            if ismember(ii,figWin) % Not useable in parallel mode (parfor)
                plt = true;
            end
%             freq_r = range(yq(:,ii)); % Frequency range of the sweep
%             dw = freq_r*1e-2;
            dw = 0;
            param = [kpe0 kpi0 w0 gc gc*1e-2]; % starting value for Param = {'kpe', 'kpi', 'w0', 'Gc', 'gma'}
%             Set up boundaries for the fitting parameters
            bound_l = [kpe0  kpi0  w0-dw   0   0]; % lower bound of fitting parameters
            bound_h = [kpe0  kpi0  w0+dw  inf inf]; % upper bound of fitting parameters            
            fit = iptopt_0(yq(:,ii),zq(:,ii),Hx(ii),Brf,param,bound_l,bound_h,weight(:,ii),plt);
            
%             if mod(ii,20) == 0
%                 worker = getCurrentTask();
%                 fprintf('Second fit, current magnetic field: %1$3.2f. on core %2$u.\n', Hx(ii), worker.ID);
%             end
            
            param = coeffvalues(fit);
            kpe(ii) = param(1);
            kpi(ii) = param(2);
            ff0(ii) = yq(fit(yq(:,ii)) == min(fit(yq(:,ii)))); % Find the resonant frequency by minimum search
%             ff0(ii) = param(3);
%             Qf(ii) = ff0(ii)/(param(1)+param(2));
            Qf(ii) = param(3)/(param(1)+param(2));
            Gc(ii) = param(4);
            gamma(ii) = param(5);
        end
%         gama = mean(gamma); % spin lifetime^-1 as fixed parameter
        
        figure
        plot(Hx,kpe,'-');
        ylabel('K_e (GHz)');
        hold on
        yyaxis right
        plot(Hx,kpi,'-');
        xlabel('Field (T)');
        ylabel('K_i (GHz)');
        title('Dissipation rates');
        legend('External dissipation rate','Internal dissipation rate')
                
        figure
        plot(Hx,kpi./kpe,'-');
        ylabel('K_i/K_e');
        xlabel('Field (T)');
        title('Ratio of Dissipation rates');
        
    case 2 %Option 2: spec1d package
        parfor ii = 1:length(Hx)
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
                fprintf('Current magnetic field: %1$3.2f. on core %2$u.\n', Hx(ii), worker.ID);
            end
        end
end

% Plot the resonant frequency from Lorentzian fit versus DC magnetic field
figure
cmap2 = copyobj(cmap0, gca);
set(cmap2, 'edgeColor','none');
shading interp;
colorbar
caxis([max([min(dB0), -30]) max(max(dB0),0)+1]);
% caxis('auto');
hold on
plot(Hx(Hc2:round(length(Hx)/200):Hc3),ff0(Hc2:round(length(Hx)/200):Hc3),'or','MarkerSize',2,'MarkerFaceColor','red');
% plot(Hx(1:length(Hx)/100:end),ff0_2(1:length(Hx)/100:end),'sk','MarkerSize',2,'MarkerFaceColor','black');
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
caxis([max([min(dB0), -30]) max(max(dB0),0)]);
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
dB0 = medfilt1(dB0, order); % apply median filter to remove some noise
plot(xq(1,:), dB0, 'o', 'MarkerSize', 2);
xlabel('Field(T)');
ylabel('S11 amplitute');
title(num2str(Temperature,'Minimal S11 at T = %3.3f K'));

%Plot Quality factor from minimum search vs magnetic field
figure
H0 = H0(Q0 >=0); %Remove unphysical points
Q0 = Q0(Q0 >=0); %Remove unphysical points
Q0 = medfilt1(Q0, order);
Qplot1 = plot(H0(1:round(length(H0)/100):end), Q0(1:round(length(Q0)/100):end),'s-','MarkerSize',2);
% plot(H0, Q0,'s-','MarkerSize',2);
%Plot Quality factor from Lorentzian fit vs magnetic field
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
    ylabel('Coupling strength')
    title(num2str(Brf,'Fitting paramters from fitting, B_0 = %.2f'))
    hold on
    yyaxis right
%     plot(Hx,gamma,'-s')
    plot(Hx(Hc2:Hc3), medfilt1(1./gamma(Hc2:Hc3),order),'-s')
    ylabel('Spin lifetime')
%     legend(['Gc, \tau_{avg} = ', num2str(1/gama,'%.2e')])
    legend('Gc','\tau')
%     legend(['\tau_0 = ', num2str(1/gama,'%.2e')]);
%     figure
%     plot(Hx, B0, '-o')
%     xlabel('Field (T)')
%     ylabel('Line cross field (T)')
end
end
% End of option 2
function option3(fileobj, nZVL)
%% Plot data
%Set data range and parameters
clear freq S11 dB N FdB FrS FiS FTT1 FTT2
order = 4; % set to what order the median filters is applied

% extract data from raw data file
out = readdata_v4(fileobj, nZVL);
freq = out.data.ZVLfreq/1e9;
S11 = out.data.ZVLreal + 1i*out.data.ZVLimag;  
H = out.data.DCField1;
HH = repmat(H,1,size(freq,2)); %populate the magnetic field to match the dimension of S11 matrix
T1 = out.data.Temperature1;
Temperature = T1(H==min(H));

% Determing the fieldscan direction
if H(end)-H(1)>0
    direction = 'Up_';
elseif H(end)-H(1)<0
    direction = 'Down_';
else
    error('Could not determing fieldscan direction!');
end
excitation = num2str(out.Power1,'_%udBm_0dB');

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
field_l = min(H);  %set field range, l: lower limit, h: higher limit
field_h = max(H);

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
dB = mag2db(abs(S11(trunc1:trunc2)));
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
f0 = zeros(size(dB,2),1);
H0 = zeros(size(dB,2),1);
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
    HM = dB(idx,ii)*0.3; % Define full-width-half-max position
    % Calculate quality factor using f0/FWHM
    if isnan(1/range(freq(dB(:,ii) <= HM)))
       Q0(ii) = 0;
    elseif isempty(range(freq(dB(:,ii) <= HM)))
       Q0(ii) = 0;
    else
    Q0(ii) = freq(idx,ii)/range(freq(dB(:,ii) <= HM));
    end
end
Q0(isinf(Q0)) = NaN; % Cut out inf from the array
[Q0,c] = rmmissing(Q0); % Cut out NaN from the array
H0 = H0(~c); % Remove corresponding elements in H0 array as well
f0 = f0(~c);
dB0 = dB0(~c);

%   For noisy data, we need to remove duplicates of minima
[H0,ia,~] = unique(H0,'stable');
f0 = f0(ia);
Q0 = Q0(ia);
dB0 = dB0(ia);

f0 = medfilt1(f0,order); % apply median filter to remove some noise
f0 = f0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
H0 = H0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
dB0 = medfilt1(dB0); % apply median filter to remove some noise
dB0 = dB0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints

[~,Hpos] = max(dB0); % find the line crossing position on field axis
hPara = [H0(Hpos), field_l, field_h];
fPara = [(freq_l+freq_h)/2, freq_l, freq_h];
[fitPara,~] = wk_cpl_fit(H0,f0,hPara,fPara);

cd(fileparts(fileobj));
tit=[direction,num2str(Temperature,'%3.3f'), excitation,'.mat'];
save(tit,'H0','f0','dB0','Q0','hPara','fPara','fitPara');
end
% End of option 3
function option4(fileobj, dataobj, Options, nZVL)
% extract data from raw data file
out = readdata_v4(fileobj, nZVL);
freq = out.data.ZVLfreq/1e9;
S11 = out.data.ZVLreal + 1i*out.data.ZVLimag;  
H = out.data.DCField1;
T1 = out.data.Temperature1;
% T2 = out.data.Temperature2;

Temperature = mean(T1(H == min(H))); %Extract the measurement temperature at the lowest field (to avoid magnetoresistance in the thermometer)
HH = repmat(H,1,size(freq,2)); %populate the magnetic field to match the dimension of S11 matrix

freq = freq';
freq = freq(:);
S11 = S11';
S11 = S11(:);
dB = mag2db(abs(S11));
HH = HH';
HH = HH(:);

[rows,~,freq] = find(freq); % Remove nonsensical zeros from the frequency data
HH = HH(rows);
S11 = S11(rows);
dB = dB(rows);
dB = dB - max(dB,[],'all'); % Shift the data to compensate changes in RE(impedence) at low temperatures

freq_l = min(freq); %set frequency range, l: lower limit, h: higher limit
freq_h = max(freq);
% freq_h = 2.503; % Manually set the upper limit of frequency scan when ZVL fails
field_l = min(H);  % set field range, l: lower limit, h: higher limit
field_h = max(H);

% Could use "scatteredInterpolant()" to replace "TriScatteredInterp()" as recommended by MATLAB, but it may generate artifacts
FdB = TriScatteredInterp(HH,freq,dB);
% FrS = TriScatteredInterp(HH,freq,real(S11));
% FiS = TriScatteredInterp(HH,freq,imag(S11)); %intrapolate points on a 2D grid

% Plot frequency-field colour map
[xq,yq] = meshgrid(linspace(field_l,field_h,301),linspace(freq_l,freq_h,601)); %set the X and Y range
zq = FdB(xq,yq);

%Clean up the raw data by removing incomplete scans (step 1) and duplicates(step 2)
% For R&S ZVL-6
% step 1: truncate the beginning and end part to keep only complete frequency scans and reshape the matrices into single column vectors
trunc1 = find(freq==freq_l,1,'first'); 
trunc2 = find(freq==freq_l,1,'last')-1; 
S11_temp = S11(trunc1:trunc2);
freq_temp = freq(trunc1:trunc2);
HH_temp = HH(trunc1:trunc2);

% Step 2: remove duplicates
dupl = find(diff(freq_temp) == 0);
S11_temp(dupl+1)=[];
freq_temp(dupl+1)=[];
dB_temp = mag2db(abs(S11_temp));
HH_temp(dupl+1)=[];

dif = nonzeros(diff(freq_temp)); % Extract from raw data the step size in frequency scan
nop = find(dif<0,1,'first'); % Calculate the number of points per complete scan
clearvars dif out step

% S11 = reshape(S11_temp,nop,[]);
dB = reshape(dB_temp,nop,[]); % reshape the matrix so that each complete frequency scan occupy one column
dB = dB - max(dB,[],'all'); % Shift the data to compensate changes in RE(impedence) at low temperatures
freq = reshape(freq_temp,nop,[]);
HH = reshape(HH_temp,nop,[]);

% % Temperary code for Keysight PNA-x N5242B
% nop = find(freq==freq_h,1)-find(freq==freq_l,1)+1;
% S11 = reshape(S11,nop,[]);
% dB = reshape(dB,nop,[]);  %reshape the matrix so that each complete frequency scan occupy one column
% freq = reshape(freq,nop,[]);
% HH = reshape(HH,nop,[]);

%Find all the resonant peaks
f0 = zeros(size(dB,2),1); % Resonant frequency from minimum search
% ff0 = zeros(size(dB,2),1); % Resonant frequency from fitting 
H0 = zeros(size(dB,2),1); % External magnetic field
FWHM = zeros(size(dB,2),1); % Full Width Half Max
Q0 = double.empty(0,length(H0)); % Quality factor from FWHM

%find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
for ii = 1:size(dB,2) %Searching column minima (fixed field)
    [~,idx] = min( dB(:,ii) );
    if(length(idx)>1)
        disp(num2str(H0(ii),'multiple minima found at H = %.3f'))
    end
    f0(ii) = freq(idx,ii);
    H0(ii) = HH(idx,ii); 
    HM = -abs(dB(idx,ii))*0.3;
    % Calculate quality factor using f0/FWHM
    if isnan(1/range(freq(dB(:,ii) <= HM)))
        Q0(ii) = 0;
    elseif isempty(range(freq(dB(:,ii) <= HM)))
        Q0(ii) = 0;
    else
        Q0(ii) = freq(idx,ii)/range(freq(dB(:,ii) <= HM));
        FWHM(ii) = range(freq(dB(:,ii) <= HM));
    end
end

% For noisy data, we need to remove duplicates of minima
[H0,ia,~] = unique(H0,'stable');
f0 = f0(ia);

f0 = medfilt1(f0); % apply median filter to remove some noise
f0 = f0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
H0 = H0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints

clearvars c idx ia ii HM trunc1 trunc2 dupl nop

switch Options.fitfunc % Pick Lorentzian fit function from either custom function of spec1d
    case 1 %Option 1: Custom function
        Qf = double.empty(0,length(H0)); % Quality factor from fitting
        Gc = double.empty(length(H0),0); % Coupling strength
        gamma = double.empty(length(H0),0); % Spin level linewidth
        parfor ii = 1:length(H0)
            % fitting by Input-output formalism
            % Param = {'kpe', 'w0', 'Gc', 'Br', 'gma'}
            param = [FWHM(ii) f0(ii) 0.01 field_l*0.1 0.1]; % Fitting parameter starting point
            bound_l = [1e-5 f0(ii) 0 field_l*0.1 0];
            bound_h = [1e-2 f0(ii) 1 field_h*0.1 1];
            % Set up boundaries for the fitting parameters
            fit = iptopt_0(freq(:,ii),-dB(:,ii),H0(ii),param,bound_l,bound_h);
            if mod(ii,20) == 0
                worker = getCurrentTask();
                fprintf('Current magnetic field: %1$3.2f. on core %2$u.\n', H0(ii), worker.ID);
            end           
            param = coeffvalues(fit);
%             ff0(ii) = param(2);
            Qf(ii) = param(2)/param(1);
            Gc(ii) = param(3);
            gamma(ii) = param(5);
        end
    case 2 %Option 2: spec1d package
        Qf = double.empty(0,size(xq,2)); % Quality factor from fitting
        zq(isnan(zq))=0;
        parfor ii = 1:size(xq,2)
            s = spec1d(yq(:,ii), db2mag(zq(:,ii)), max(db2mag(zq(:,ii))))*0.001; % create spec1d object
%             %starting point for the (Lorentzian) fitting parameters
%             p = [0.1 ff0(ii) FWHM(ii) min(zq(:,ii))]; % (p1: scaling factor ,p2: resonant frequency; p3: FWHM; p4:noise floor(?) )
%             fix = [0 0 0 0]; % Denoting if the fitting parameters are fixed
%             [~, fbck] = fits(s, 'lorz', p, fix);
            [~, fbck] = fits(s, 'lorz');
%             ff0(ii) = fbck.pvals(2); % Retrieve the resonant frequency from fitted data
            Qf(ii) = abs(fbck.pvals(2)/fbck.pvals(3)/2); %Calculate the quality factor
%             chi(ii) = 1/Qf(ii);
            if mod(ii,20) == 0
                worker = getCurrentTask();
                fprintf('Current magnetic field: %1$3.2f. on core %2$u.\n', xq(1,ii), worker.ID);
            end
        end
end
        
% Plot the interpolated frequency response data in a field scan using color map
figure
cmap = pcolor(xq,yq,zq);
hold on
plot(H0,f0,'.r','MarkerSize',6);
set(cmap, 'edgeColor','none')
shading interp;
colorbar
set(gca,'fontsize',Options.ftsz)
xlabel('Field (T)');
ylabel('Frequency (GHz)');
xticks(linspace(field_l,field_h,6));
title(num2str(Temperature,'S11 response at T = %3.3f K'));

% Plot the temperature profile against magnetic field
figure
plot(H(1:100:end),T1(1:100:end),'o-')
xlabel('DC Magnetic field')
ylabel('Temperature')
title('Magnetic field vs Temperature')

figure
plot(H0,f0,'ok','MarkerSize',4);
hold on
plot(mean(H0(f0==min(f0))), mean(f0(f0==min(f0))), 'o', 'MarkerFaceColor', 'red','MarkerSize',4);
xlabel('Field (T)');
ylabel('Resonant frequency (GHz)');
title(num2str(Temperature,'Resonant frequency from minimum search at T = %3.3f K'));

[~,idx,Q0] = find(Q0);
figure
plot(H0(idx),Q0,'-k');
hold on
[~,idx,Qf] = find(Qf);
plot(xq(1,idx),Qf,'-r');
xlabel('Field (T)');
ylabel('Q');
title(num2str(Temperature,'Inverse of Q factor at T = %3.3f K'));
legend('Minimum search','Lorentzian fit');

if Options.fitfunc == 1
    figure
    plot(H0,Gc,'-o')
    xlabel('Field (T)')
    ylabel('Coupling strenght')
    title('Fitting paramters from citting')
    hold on
    yyaxis right
    plot(H0, 1./gamma,'-s')
    ylabel('Spin lifetime')
    legend('Gc','\tau')
end

% Save the data
phase(1,1) = Temperature;
phase(1,2) = mean(H0(f0==min(f0)));

if isfile(dataobj)
    save(dataobj,'phase','-append');
else
    save(dataobj,'phase','-v7.3');
end

end
% End of option 4
function option5(fileobj, Options, nZVL)
% extract data from raw data file
out = readdata_v4(fileobj, nZVL);
freq = out.data.ZVLfreq/1e9;
S11 = out.data.ZVLreal + 1i*out.data.ZVLimag;  
H = out.data.DCField1;
T1 = out.data.Temperature1;
% T2 = out.data.Temperature2;

Field = mean(H);
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
% End of option 5
function saveplots(hfig,figname)
    % Function to save the figures to a specified directory and in different
    % file formats
    % INPUT:
    % hfig      figure handle
    % figname   file name to save figure to
    curdir = cd;
    cd('G:\My Drive\File sharing\PhD projects\LiHoF4\Data\Matlab\Cavity resonator\D24mm_T5mm_G0.2mm\12.02.2019')
    %cd('/Users/yikaiyang/Google Drive/File sharing/PhD projects/LiHoF4/Data/Matlab/Cavity resonator/D24mm_T5mm_G0.2mm/12.02.2019')
    
    saveas(figure(hfig),[figname '.fig'],'fig');
    print(figure(hfig),[figname  '.jpg'],'-djpeg','-r600');
    print(figure(hfig),[figname  '.png'],'-dpng','-r600');
    print2eps(figname,hfig)
    [~,~] = eps2xxx([figname '.eps'],{'jpeg','pdf'});
    
    disp(['Figure ' figname ' saved to '])
    disp(cd)
    cd(curdir)
end