function Field_scan_v4
    format long;  
    addpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\Fieldscan\functions\');
    addpath(genpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\spec1d--Henrik\'));
%     addpath('/Volumes/GoogleDrive/My Drive/File sharing/Programming scripts/Matlab/Plot functions/Fieldscan/functions')
%     addpath(genpath('/Volumes/GoogleDrive/My Drive/File sharing/Programming scripts/Matlab/Plot functions/spec1d--Henrik'));
    %The first line is for windows, the second line is for mac OS

    % Figure plot options:
    plotopt.plot = false;
    plotopt.lnwd = 2;
    plotopt.ftsz = 12;
    plotopt.mksz = 3;

    filepath = 'G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\SC200\2020.10.09';
%     filepath = '/Volumes/GoogleDrive/My Drive/File sharing/PhD program/Research projects/LiHoF4 project/Data/Experiment/LiHoF4/SC127/SC127_2 (2.5 x 1 x 0.5 mm, triangle)/25.07.2020';
    %The first line is for windows, the second line is for mac OS
    filename = '2020_10_0009.dat';
    nZVL = 1; % Number of dataset from ZVL
    fileobj = fullfile(filepath,filename);
    
    datapath = 'G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\SC200';
    dataname = 'Phase Diagram.mat';
    dataobj = fullfile(datapath,dataname);
    
% Operation options
    opt  = 4;

    switch opt
        case 1
            % Simple color plot of S11 (+ S21)
            option1(fileobj, plotopt, nZVL)
        case 2
            % Color plot + Data fitting
            option2(fileobj, plotopt, nZVL)
        case 3
            % Data fitting and file saving (w/o plots)
            option3(fileobj, plotopt, nZVL)
        case 4
            % Off-resonance measurement processing
            option4(fileobj, dataobj, plotopt, nZVL)
    end
end

function option1(fileobj, plotopt, nZVL)

out = readdata_v4(fileobj,nZVL);
freq = out.data.ZVLfreq/1e9;
S11 = out.data.ZVLreal + 1i*out.data.ZVLimag;
H = out.data.DCField1;
HH = repmat(H,1,size(freq,2)); %populate the magnetic field to match the dimension of S11 matrix
T1 = out.data.Temperature1;
% T2 = out.data.Temperature2;
lck1 = out.data.lockin1;
lck2 = out.data.lockin2;

Temperature = mean(T1(H == min(H)));

S11 = S11';
S11 = S11(:);
dB = mag2db(abs(S11));
freq = freq';
freq = freq(:);
HH = HH';
HH = HH(:); %the third argument is the number of frequency points in each line/segment

[rows,~,freq] = find(freq); % Remove nonsensical zeros from the frequency data
HH = HH(rows);
dB = dB(rows);

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
    set(gca,'fontsize',plotopt.ftsz)
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
set(gca,'fontsize',plotopt.ftsz)
xlabel('Field (T)');
ylabel('Frequency (GHz)');
xticks(linspace(field_l,field_h,6));
title(num2str(Temperature,'S11 response at T = %3.3f K'));

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
plot(H(1:100:end),T1(1:100:end),'o-')
xlabel('DC Magnetic field')
ylabel('Temperature')
title('Magnetic field vs Temperature')
end

function option2(fileobj, plotopt, nZVL)
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
% freq_h = 2.503; % Manually set the upper limit of frequency scan when ZVL fails
field_l = min(H);  % set field range, l: lower limit, h: higher limit
field_h = max(H);

dif = nonzeros(diff(freq)); % Extract from raw data the step size in frequency scan
dif = dif(dif>0); % Keep only positive steps
step = min (dif); % Calculate the value of the frequency scan step
nop = ceil(abs(freq_h-freq_l)/step); %compute how many points pers complete frequency scan.

clearvars dif out step

%Plot the temperature vs magnetic field to check the temperature variation
figure
plot(H(1:length(T1)/100:end),T1(1:length(T1)/100:end),'o-')
xlabel('DC Magnetic field')
ylabel('Temperature')
title('Magnetic field vs Temperature')

%Clean up the raw data by removing incomplete scans (step 1) and duplicates(step 2)
% step 1: truncate the beginning and end part to keep only complete frequency scans and reshape the matrices into single column vectors
trunc1 = find(freq==freq_l,1,'first'); 
trunc2 = find(freq==freq_l,1,'last')-1; 
S11_temp = S11(trunc1:trunc2);
dB_temp = mag2db(abs(S11_temp));
freq_temp = freq(trunc1:trunc2);
HH_temp = HH(trunc1:trunc2);

% Step 2: remove duplicates
dupl = find(diff(freq_temp) == 0);
freq_temp(dupl+1)=[];
dB_temp(dupl+1)=[];
HH_temp(dupl+1)=[];

dB = reshape(dB_temp,nop,[]);  %reshape the matrix so that each complete frequency scan occupy one column
freq = reshape(freq_temp,nop,[]);
HH = reshape(HH_temp,nop,[]);

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
    HM = dB(idx,ii)*0.3;
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
FWHM = FWHM(~c);

% For noisy data, we need to remove duplicates of minima
[H0,ia,~] = unique(H0,'stable');
f0 = f0(ia);
Q0 = Q0(ia);
dB0 = dB0(ia);
FWHM = FWHM(ia);
[~,Hpos] = max(dB0); % find the line crossing position on field axis

f0 = medfilt1(f0,order); % apply median filter to remove some noise
f0 = f0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
H0 = H0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
Q0 = Q0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
dB0 = dB0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints

clearvars c idx ia ii HM T1 trunc1 trunc2 dupl nop

% Fit the field dependent resonant frequency data with weak coupling function
hPara = [H0(Hpos), field_l, field_h];
fPara = [(freq_l+freq_h)/2, freq_l, freq_h];
[fitP,~] = wk_cpl_fit(H0,f0,hPara,fPara);

% Coupling curve with fit parameters
B = linspace(field_l,field_h,100);
spin = 7/2;
Delt = -spin*(B-fitP.x0);
plot(H0(1:round(length(H0)/200):end),f0(1:round(length(f0)/200):end),'ok','MarkerSize',4);
hold on
wp = fitP.wc + Delt./2 + sqrt(Delt.^2+4*fitP.g^2)/2;
wm = fitP.wc + Delt./2 - sqrt(Delt.^2+4*fitP.g^2)/2;
plot(B,wm,'-r',B,wp,'-r','LineWidth',plotopt.lnwd);
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

clearvars B Delt hPara fPara wp wm spin fitPara H_res f_res Hpos gc

% Interpolate the data on a 2D grid for the colormap
[xq,yq] = meshgrid(linspace(field_l,field_h,301),linspace(freq_l,freq_h,310));

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
%             FdB  = TriScatteredInterp(HH,freq,dB);
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
    caxis([max([min(dB0), -30]) 5]);
end
    clearvars *_temp;

set(gca,'fontsize',plotopt.ftsz)
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
        Qf = double.empty(length(Hx),0);
    case 2 % Use interpolated data
        Hx = xq(1,:);
        ff0 = interp1(H0,f0,Hx,'pchip','extrap'); % Using spline interpolation to smooth the resonant frequency trace
        Qf = double.empty(length(Hx),0);
end

switch 2 % Pick Lorentzian fit function from either custom function of spec1d
    case 1 %Option 1: Custom function
%        parfor ii = 1:length(Hx)
            % fitting by Input-output formalism
            % Fitting parameter starting point
            param = [FWHM(ii) ff0(ii) 1e-3 Br 1e-3]; % Param = {'kpe', 'w0', 'Gc', 'Br', 'gma'}
            bound_l = [0 freq_l 0 Br*0.5 0]; % lower bound of fitting parameters
            bound_h = [inf freq_h 5 Br*1.5 0.1]; % upper bound of fitting parameters
            % Set up boundaries for the fitting parameters
            fit = iptopt(yq(:,ii),zq(:,ii),Hx(ii),param,bound_l,bound_h);
            if mod(ii,20) == 0
                worker = getCurrentTask();
                fprintf('Current magnetic field: %1$3.2f. on core %2$u.\n', Hx(ii), worker.ID);
            end
            param = coeffvalues(fit);
            ff0(ii) = param(2);
            Qf(ii) = param(2)/param(1);
        end
    case 2 %Option 2: spec1d package
        parfor ii = 1:length(Hx)
            s = spec1d(yq(:,ii), -zq(:,ii), max(-zq(:,ii)))*0.001; % create spec1d object
            %starting point for the (Lorentzian) fitting parameters
            p = [0.1 ff0(ii) FWHM(ii) min(zq(:,ii))]; % (p1: scaling factor ,p2: resonant frequency; p3: FWHM; p4:noise floor(?) )
            fix = [0 0 0 0]; % Denoting if the fitting parameters are fixed
            [~, fbck] = fits(s, 'lorz', p, fix);
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
caxis([max([min(dB0), -30]) 5]);
axis([field_l field_h freq_l freq_h]);

% Plot the resonant frequency from minimum search versus DC magnetic field
% figure
hold on
f0 = medfilt1(f0,order); % apply median filter to remove some noise
% hfig1 = plot(H0(1:round(length(H0)/200):end),f0(1:round(length(f0)/200):end),'ok','MarkerSize',2,'MarkerFaceColor','black');
plot(H0, f0, 'o', 'MarkerSize', 2);
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
legend('Quality factor from Lorentzian fit', 'Quality factor from FWHM');
title(num2str(Temperature,'Quality factor, T= %.3f'));
end

function option3(fileobj, plotopt, nZVL)
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
excitation = '_-20dBm_0dB';

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

figure
f0 = medfilt1(f0,order); % apply median filter to remove some noise
f0 = f0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
H0 = H0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
plot(H0(1:round(length(H0)/100):end),f0(1:round(length(f0)/100):end),'ok','MarkerSize',2,'MarkerFaceColor','black');
% hfig1 = plot(H0, f0, 'o', 'MarkerSize', 2);
xlabel('Field (T)');
ylabel('Resonant frequency (GHz)');
title(num2str(Temperature,'Resonant frequency from minimum search at T = %3.3f K'));
axis([field_l field_h freq_l freq_h]);

figure
dB0 = medfilt1(dB0); % apply median filter to remove some noise
dB0 = dB0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
plot(H0(1:length(H0)/100:end), dB0(1:length(dB0)/100:end), 'o', 'MarkerSize', 2);
xlabel('Field(T)');
ylabel('S11 amplitute');
title(num2str(Temperature,'Minimal S11 at T = %3.3f K'));
set(gca,'fontsize',plotopt.ftsz)

[~,Hpos] = max(dB0); % find the line crossing position on field axis
hPara = [H0(Hpos), field_l, field_h];
fPara = [(freq_l+freq_h)/2, freq_l, freq_h];
[fitPara,~] = wk_cpl_fit(H0,f0,hPara,fPara);

cd(fileparts(fileobj));
tit=[direction,num2str(Temperature,'%3.3f'), excitation,'.mat'];
save(tit,'H0','f0','dB0','hPara','fPara','fitPara');
end

function option4(fileobj, dataobj, plotopt, nZVL)
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
[xq,yq] = meshgrid(linspace(field_l,field_h,301),linspace(freq_l,freq_h,310)); %set the X and Y range
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
HH_temp = HH(trunc1:trunc2);

% Step 2: remove duplicates
dupl = find(diff(freq_temp) == 0);
freq_temp(dupl+1)=[];
dB_temp(dupl+1)=[];
HH_temp(dupl+1)=[];

dB = reshape(dB_temp,nop,[]);  %reshape the matrix so that each complete frequency scan occupy one column
freq = reshape(freq_temp,nop,[]);
HH = reshape(HH_temp,nop,[]);

%Find all the resonant peaks
f0 = zeros(size(dB,2),1);
H0 = zeros(size(dB,2),1);

%find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
for ii = 1:size(dB,2) %Searching column minima (fixed field)
    [~,idx] = min( dB(:,ii) );
    if(length(idx)>1)
        disp(num2str(H0(ii),'multiple minima found at H = %.3f'))
    end
    f0(ii) = freq(idx,ii);
    H0(ii) = HH(idx,ii); 
end

% For noisy data, we need to remove duplicates of minima
[H0,ia,~] = unique(H0,'stable');
f0 = f0(ia);

f0 = medfilt1(f0); % apply median filter to remove some noise
f0 = f0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints
H0 = H0(f0 >= freq_l & f0 <= freq_h); % Discard nonsensical datapoints

clearvars c idx ia ii HM trunc1 trunc2 dupl nop

% Plot the interpolated frequency response data in a field scan using color map
figure
cmap = pcolor(xq,yq,zq);
hold on
plot(H0(1:round(length(H0)/200):end),f0(1:round(length(f0)/200):end),'.r','MarkerSize',6);
set(cmap, 'edgeColor','none')
shading interp;
colorbar
set(gca,'fontsize',plotopt.ftsz)
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
plot(H0(1:round(length(H0)/200):end),f0(1:round(length(f0)/200):end),'ok','MarkerSize',4);
hold on
plot(mean(H0(f0==min(f0))), mean(f0(f0==min(f0))), 'o', 'MarkerFaceColor', 'red','MarkerSize',4);
xlabel('Field (T)');
ylabel('Resonant frequency (GHz)');
title(num2str(Temperature,'Resonant frequency from minimum search at T = %3.3f K'));

% Save the data
phase(1,1) = Temperature;
phase(1,2) = mean(H0(f0==min(f0)));
if isfile(dataobj)
    save(dataobj,'phase','-append');
else
    save(dataobj,'phase','-v7.3');
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
    
    saveas(figure(hfig),[figname '.fig'],'fig');
    print(figure(hfig),[figname  '.jpg'],'-djpeg','-r600');
    print(figure(hfig),[figname  '.png'],'-dpng','-r600');
    print2eps(figname,hfig)
    [~,~] = eps2xxx([figname '.eps'],{'jpeg','pdf'});
    
    disp(['Figure ' figname ' saved to '])
    disp(cd)
    cd(curdir)
end