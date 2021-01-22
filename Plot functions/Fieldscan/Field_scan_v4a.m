function Field_scan_v4a
    format long;  
    addpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\Fieldscan\functions\');
    addpath(genpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\spec1d--Henrik\'));
%     addpath('/Volumes/GoogleDrive/My Drive/File sharing/Programming scripts/Matlab/Plot functions/Fieldscan/functions')
%     addpath(genpath('/Volumes/GoogleDrive/My Drive/File sharing/Programming scripts/Matlab/Plot functions/spec1d--Henrik'));
    %The first line is for windows, the second line is for mac OS

    % Figure plot options:
    Options.plot = false; % Plot option in analysis part (Option 3)
    Options.lnwd = 2;
    Options.ftsz = 12;
    Options.mksz = 3;
    Options.fitfunc = 1; % Pick fitting function from either (1) custom function of (2) spec1d
    
    loadpath = 'G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\SC107 (4x5x2mm)\19.05.2019';
%     loadpath = '/Volumes/GoogleDrive/My Drive/File sharing/PhD program/Research projects/LiHoF4 project/Data/Experiment/LiHoF4/SC199/2020.11.05';
    %The first line is for windows, the second line is for mac OS
    loadname = '2019_05_0025.dat';
    opt = 4;% Analysis options
    nZVL = 1; % Number of dataset from ZVL
    fileobj = fullfile(loadpath,loadname);

    % Path and file name to save
    savepath = 'G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\SC200';
    savename = 'Phase Diagram.mat';
    dataobj = fullfile(savepath,savename);

    switch opt
        case 1
            % Simple color plot of S11 (+ S21)
            option1(fileobj, Options, nZVL)
        case 2
            % On-resonance measurement processing (w/ plots)
            option2(fileobj, Options, nZVL)
        case 3
            % Data fitting and file saving (w/o plots)
            option3(fileobj, nZVL)
        case 4
            % Off-resonance measurement processing
            option4(fileobj, dataobj, Options, nZVL)
        case 5
            % Temperature scan
            option5(fileobj, Options, nZVL)
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

function option2(fileobj, Options, nZVL)
%Set data range and parameters
order = 4; % set to what order the median filters is applied
clear freq S11 dB N FdB FrS FiS FTT1 FTT2

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
HH = HH';
HH = HH(:);

[rows,~,freq] = find(freq); % Remove nonsensical zeros from the frequency data
HH = HH(rows);
S11 = S11(rows);

freq_l = min(freq); %set frequency range, l: lower limit, h: higher limit
freq_h = max(freq);
field_l = min(H);  % set field range, l: lower limit, h: higher limit
field_h = max(H);

%Plot the temperature vs magnetic field to check the temperature variation
figure
plot(H(1:length(T1)/100:end),T1(1:length(T1)/100:end),'o-')
xlabel('DC Magnetic field')
ylabel('Temperature')
title('Magnetic field vs Temperature')

%% Code For R&S ZVL-6
% Clean up the raw data by removing incomplete scans (step 1) and duplicates(step 2)
% step 1: truncate the beginning and end part to keep only complete frequency scans and reshape the matrices into single column vectors
trunc1 = find(freq==freq_l,1,'first'); 
trunc2 = find(freq==freq_l,1,'last')-1; 
S11_temp = S11(trunc1:trunc2);
freq_temp = freq(trunc1:trunc2);
HH_temp = HH(trunc1:trunc2);

% Step 2: remove duplicates
dupl = find(diff(freq_temp) == 0.0);
freq_temp(dupl+1)=[];
S11_temp(dupl+1)=[];
HH_temp(dupl+1)=[];
dB_temp = mag2db(abs(S11_temp));

dif = diff(freq); % frequency increments
resets = find(dif<=0.9*(freq_l-freq_h)); % Find the termination point of a complete scans (10% error)
nop = mean(diff(resets)); % Set the number of points per frequency scan

% % Alternative way of finding number of points per complete frequency scan
% dif = nonzeros(diff(freq)); % Extract from raw data the step size in frequency scan
% dif = dif(dif>0); % Keep only positive steps
% step = min (dif); % Calculate the value of the frequency scan step
% nop = ceil(abs(freq_h-freq_l)/step); %compute how many points pers complete frequency scan.

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
% nop = mean(diff(resets)); % Set the number of points per frequency scan
% dB = reshape(mag2db(abs(S11)),nop,[]);  %reshape the matrix so that each complete frequency scan occupy one column
% dB = dB - max(dB,[],'all'); % Shift the data to compensate changes in RE(impedence) at low temperatures
% freq = reshape(freq,nop,[]);
% HH = reshape(HH,nop,[]);
%% end of temporary code for PNA-x N5242B
clearvars dif out step

%Find all the resonant peaks
f0 = zeros(size(dB,2),1);
H0 = zeros(size(dB,2),1);
Q0 = zeros(size(dB,2),1);
dB0 = zeros(size(dB,2),1);
FWHM = zeros(size(dB,2),1);
%find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
for ii = 1:size(dB,2) %Searching column minima (fixed field)
    [~,idx] = min( dB(:,ii) );
    if(length(idx)>1)
        disp(num2str(H0(ii),'multiple minima found at H = %.3f'))
    end
    f0(ii) = freq(idx,ii);
    H0(ii) = HH(idx,ii); 
    dB0(ii) = dB(idx,ii);
%     HM = dB(idx,ii)*0.3;
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
Q0(isinf(Q0)) = NaN; % Cut out inf from the array
[Q0,c] = rmmissing(Q0); % Cut out NaN from the array
H0 = H0(~c); % Remove corresponding elements in H0 array as well
f0 = f0(~c);
dB0 = dB0(~c);
FWHM0 = FWHM(~c);

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

clearvars c idx ia ii HM T1 trunc1 trunc2 dupl nop

% Fit the field dependent resonant frequency data with weak coupling function
hPara = [H0(Hpos), field_l, field_h];
fPara = [(freq_l+freq_h)/2, freq_l, freq_h];
[fitP,~] = wk_cpl_fit(H0,f0,hPara,fPara);
Br = fitP.x0;

% Coupling curve with fit parameters
B = linspace(field_l,field_h,100);
spin = 7/2;
Delt = -spin*(B-fitP.x0);
plot(H0(1:round(length(H0)/200):end),f0(1:round(length(f0)/200):end),'ok','MarkerSize',4);
hold on
wp = fitP.wc + Delt./2 + sqrt(Delt.^2+4*fitP.g^2)/2;
wm = fitP.wc + Delt./2 - sqrt(Delt.^2+4*fitP.g^2)/2;
plot(B,wm,'-r',B,wp,'-r','LineWidth',Options.lnwd);
% hfig1 = plot(H0, f0, 'o', 'MarkerSize', 2);
xlabel('Field (T)');
ylabel('Resonant frequency (GHz)');
title(num2str(Temperature,'Resonant frequency from minimum search at T = %3.3f K'));
legend(sprintf('gc = %.3f GHz',fitP.g));
axis([field_l field_h freq_l freq_h]);

% Plot frequency scan at line crossing
figure
% plot(freq(1:10:end,Hpos),dB(1:10:end,Hpos),'-o');
plot(freq(:,Hpos),dB(:,Hpos),'-o');
xlabel('Frequency (GHz)');
ylabel('S11 (dB)');
legend(sprintf('Frequency cut at %.2f T',H0(Hpos)));
title('Frequency scan at line crossing');

clearvars B Delt hPara fPara wp wm spin fitPara H_res f_res gc Hpos

% Interpolate the data on a 2D grid for the colormap
[xq,yq] = meshgrid(linspace(field_l,field_h,301),linspace(freq_l,freq_h,601));

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
                cmap = scatter3(HH,yy,interp_dB,2,interp_dB,'o','filled','MarkerEdgeColor','none');
                colormap(hsv);
        end
    case 2 % Option_2 Interpolate the data along both axis.
        %    FdB  = TriScatteredInterp(HH,freq,dB);
        FdB = scatteredInterpolant(HH_temp, freq_temp, dB_temp);
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
        caxis('auto')
        %     caxis([max([min(dB0), -30]) 5]);
end
clearvars *_temp;

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
        zq = dB;
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
zq(isnan(zq))=0;
switch Options.fitfunc % Pick fitting function from either (1) custom function of (2) spec1d
    case 1 %Option 1: Custom function by Input-output formalism
        parfor ii = 1:length(Hx)
            % Fitting parameter starting point
            param = [FWHM(ii) ff0(ii) 1e-3 Br 1e-3]; % Param = {'kpe', 'w0', 'Gc', 'Br', 'gma'}
            % Set up boundaries for the fitting parameters
            bound_l = [0 ff0(ii) 0 field_l 0]; % lower bound of fitting parameters
            bound_h = [inf ff0(ii) 5 field_h 5]; % upper bound of fitting parameters
            fit = iptopt_0(yq(:,ii),zq(:,ii),Hx(ii),param,bound_l,bound_h);

%             param = [FWHM(ii) ff0(ii) 1e-3 1e-3 Br 1e-3]; % Param = {'kpe', 'w0', 'x1', 'x2', 'Br', 'gma'}
%             bound_l = [0 freq_l 0 0 field_l 0]; % lower bound of fitting parameters
%             bound_h = [inf freq_h 5 5 field_h 0.1]; % upper bound of fitting parameters
%             fit = iptopt(yq(:,ii),zq(:,ii),Hx(ii),param,bound_l,bound_h);

            if mod(ii,20) == 0
                worker = getCurrentTask();
                fprintf('Current magnetic field: %1$3.2f. on core %2$u.\n', Hx(ii), worker.ID);
            end
            param = coeffvalues(fit);
            ff0(ii) = param(2);
            Qf(ii) = param(2)/param(1);
            Gc(ii) = param(3);
            gamma(ii) = param(5);
        end
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
% figure
plot(Hx(1:round(length(Hx)/200):end),ff0(1:round(length(Hx)/200):end),'or','MarkerSize',2,'MarkerFaceColor','red');
% plot(Hx(1:length(Hx)/100:end),ff0_2(1:length(Hx)/100:end),'sk','MarkerSize',2,'MarkerFaceColor','black');
xlabel('Field (T)');
ylabel('Frequency (GHz)');
title(num2str(Temperature,'Resonant frequency from fitted data at T = %3.3f K'));
axis([field_l field_h freq_l freq_h]);

%copy over the colormap of the raw data as the background of the next plot
figure
cmap2 = copyobj(cmap, gca);
set(cmap2, 'edgeColor','none');
shading interp;
colorbar
% caxis([max([min(dB0), -30]) 5]);
caxis('auto');
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
plot(H0, dB0, 'o', 'MarkerSize', 2);
xlabel('Field(T)');
ylabel('S11 amplitute');
title(num2str(Temperature,'Minimal S11 at T = %3.3f K'));

%Plot Quality factor from minimum search vs magnetic field
figure
H0 = H0(Q0 >=0);
Q0 = Q0(Q0 >=0);
Q0 = medfilt1(Q0, 10);
% Qplot2 = plot(H0(1:round(length(H0)/100):end), Q0(1:round(length(Q0)/100):end),'s-','MarkerSize',2);
plot(H0, Q0,'s-','MarkerSize',2);
%Plot Quality factor from Lorentzian fit vs magnetic field
hold on
Hx = Hx(Qf >=0);
Qf = Qf(Qf >=0); %Remove unphysical points
Qf = medfilt1(Qf, order);
% Qplot1 = plot(Hx(1:round(length(Hx)/100):end), Qf(1:round(length(Qf)/100):end),'o-','MarkerSize',2);
plot(Hx, Qf,'o-','MarkerSize',2);
gca;
xlabel('Field (T)');
ylabel('Q factor');
legend('Quality factor from FWHM', 'Quality factor from Lorentzian fit');
title(num2str(Temperature,'Quality factor, T= %.3f'));

if Options.fitfunc == 1
    figure
    plot(Hx,Gc,'-o')
    xlabel('Field (T)')
    ylabel('Coupling strenght')
    title('Fitting paramters from citting')
    hold on
    yyaxis right
    plot(Hx, gamma,'-s')
    ylabel('Spin lifetime')
    legend('Ge','gamma')
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
save(tit,'H0','f0','dB0','hPara','fPara','fitPara');
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
            bound_l = [1e-4 freq_l 0 0 0];
            bound_h = [1e-2 freq_h 1 field_h*10 1];
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
            s = spec1d(yq(:,ii), -zq(:,ii), max(-zq(:,ii)))*0.001; % create spec1d object
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
xlabel('Field (T)');
ylabel('1/Q');
hold on
[~,idx,Qf] = find(Qf);
plot(xq(1,idx),Qf,'-r');
xlabel('Field (T)');
ylabel('1/Q');
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
    plot(H0, gamma,'-s')
    ylabel('Spin lifetime')
    legend('Ge','gamma')
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