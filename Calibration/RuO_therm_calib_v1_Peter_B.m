function RuO_thermometer_calibration_v1
% Function to use CMN and SRD calibration thermometers to calibrate the
% Kelvinox 400 dilution fridge thermometers

%% Load data ==============================================================
% Calibration data folder
pathMIDS = 'C:\Users\Peter\Documents\Projects\LiReF4\20130202 - Kelvinox 400 thermometer calibration';
% Fridge log file folder
pathKelvinox = 'C:\Users\Peter\Documents\Projects\LiReF4\20130202 - Kelvinox 400 thermometer calibration';

opt = 1;

switch opt
    case 1
        % 1st run
        fileMIDS = '20130202_MIDS202-221-warmscan.dat';
        fileKelvinox = 'log 130202 114754.csv';
        dataMIDS = load_MIDS(pathMIDS,fileMIDS);
        dataKelvinox = load_kelvinox(pathKelvinox,fileKelvinox);
        % Interpolate the temperature over serial date number range specified:
        rtimeserial = [min(dataMIDS.timeserial)  7.352681398670982e5];   
        
    case 2
        % 2nd run
        fileMIDS = '20130203_MIDS202-221-warmscan.dat';
        fileKelvinox = 'log 130202 114754.csv';
        dataMIDS = load_MIDS(pathMIDS,fileMIDS);
        dataKelvinox = load_kelvinox(pathKelvinox,fileKelvinox);
        rtimeserial = [7.352685622119815e5 7.352692350230415e5];
        
    case 3
        % 3rd run
        fileMIDS = '20130204_MIDS202-221-warmscan.dat';
        fileKelvinox = 'log 130202 114754.csv';
        dataMIDS = load_MIDS(pathMIDS,fileMIDS);
        dataKelvinox = load_kelvinox(pathKelvinox,fileKelvinox);
        rtimeserial = [7.352697657061547e5 7.352703059210526e5];
        
    case 4
        % 4th run - going up in steps and doing slow scans close to Tc
        fileMIDS = '20130205_MIDS202-221-steps.dat';
        fileKelvinox = 'log 130205 145805.csv';
        dataMIDS = load_MIDS(pathMIDS,fileMIDS);
        dataKelvinox = load_kelvinox(pathKelvinox,fileKelvinox);
        rtimeserial = [min(dataMIDS.timeserial) max(dataMIDS.timeserial)];
        
    case 5
        % 5th run - slow scan at 0.001 K/min from base to 0.6K
        fileMIDS = '20130205_MIDS202-221-warmscan.dat';
        fileKelvinox = 'log 130205 145805.csv';
        dataMIDS = load_MIDS(pathMIDS,fileMIDS);
        dataKelvinox = load_kelvinox(pathKelvinox,fileKelvinox);
        rtimeserial = [min(dataMIDS.timeserial) max(dataMIDS.timeserial)];      
        
    case 6
        % 6th run - slow heater ramp to pass through lowest two transitions
        fileMIDS = '20130206_MIDS202-221-1.dat';
        fileKelvinox = 'log 130205 145805.csv';
        dataMIDS = load_MIDS(pathMIDS,fileMIDS);
        dataKelvinox = load_kelvinox(pathKelvinox,fileKelvinox);
        rtimeserial = [min(dataMIDS.timeserial) 7.352716451233947e+05];
        
    case 7
        % Measurements of the higher-temperature transition points
        fileMIDS = '20130212_MIDS202-221-1.dat';
        fileKelvinox = 'log 130212 120437.csv';
        dataMIDS = load_MIDS(pathMIDS,fileMIDS);
        dataKelvinox = load_kelvinox(pathKelvinox,fileKelvinox);
        rtimeserial = [min(dataMIDS.timeserial) max(dataMIDS.timeserial)];
        
    case 8
        % Measurements of the higher-temperature transition points
        fileMIDS = '20130212_MIDS202-221-2.dat';
        fileKelvinox = 'log 130212 120437.csv';
        dataMIDS = load_MIDS(pathMIDS,fileMIDS);
        dataKelvinox = load_kelvinox(pathKelvinox,fileKelvinox);
        rtimeserial = [min(dataMIDS.timeserial) max(dataMIDS.timeserial)];
        
    case 9
        % Look at the 3.4 K phase transition
        fileMIDS = '20130212_MIDS202-221-2.dat';
        fileKelvinox = 'log 130212 120437.csv';
        dataMIDS = load_MIDS(pathMIDS,fileMIDS);
        dataKelvinox = load_kelvinox(pathKelvinox,fileKelvinox);
        rtimeserial = [7.352776669139069 7.352776795554093]*1e5;
        
    case 10
        % Look at the 4.9 K phase transition
        fileMIDS = '20130212_MIDS202-221-2.dat';
        fileKelvinox = 'log 130212 120437.csv';
        dataMIDS = load_MIDS(pathMIDS,fileMIDS);
        dataKelvinox = load_kelvinox(pathKelvinox,fileKelvinox);
        rtimeserial = [7.352776906919710 7.352777076978016]*1e5;
        
    case 11
        % Look at the 1.1 K phase transition
        fileMIDS = '20130212_MIDS202-221-3.dat';
        fileKelvinox = 'log 130212 120437.csv';
        dataMIDS = load_MIDS(pathMIDS,fileMIDS);
        dataKelvinox = load_kelvinox(pathKelvinox,fileKelvinox);
        rtimeserial = [7.352780006285468 7.352780938763814]*1e5;

    case 12
        % Look at the 0.8 K phase transition
        fileMIDS = '20130212_MIDS202-221-2.dat';
        fileKelvinox = 'log 130212 120437.csv';
        dataMIDS = load_MIDS(pathMIDS,fileMIDS);
        dataKelvinox = load_kelvinox(pathKelvinox,fileKelvinox);
        rtimeserial = [7.352775927203275 7.352776025024425]*1e5;

end





%% Analyse results ========================================================
% Plot mixing
hfig1 = figure(1); 
% set(hfig1,'position',[-1591        -171         784         808])
clf

subplot(2,1,1)
plot(dataKelvinox.timeserial,dataKelvinox.diag.MC.T)
xlim(rtimeserial)
ylabel('M.C. (K)','interpreter','latex')

subplot(2,1,2)
plot(dataMIDS.timeserial,dataMIDS.VA,dataMIDS.timeserial,dataMIDS.VB)
xlim(rtimeserial)
xlabel('time (arb.)','interpreter','latex')
ylabel('V (V)','interpreter','latex')

dataMIDS = cutMIDS(dataMIDS,rtimeserial);
dataKelvinox = cutKelvinox(dataKelvinox,rtimeserial);

% set(get(hfig1,'children'),'xlim',rtimeserial)
% set(get(hfig2,'children'),'xlim',rtimeserial)

Tin = interp1(dataKelvinox.timeserial,dataKelvinox.diag.MC.T,dataMIDS.timeserial);
Rin = interp1(dataKelvinox.timeserial,dataKelvinox.diag.MC.R,dataMIDS.timeserial);

% Plot figures using time as the reference point ==========================

data_fixed = calibrationdata_SRD1000;

hfig3 = figure(4);
set(hfig3, ... %'position',[-791  -171   784   808],...
            'tag','cal_kelvinox')
clf

subplot(3,1,1)
plot(Tin,dataMIDS.VA)
ylabel('Channel A V (V)','interpreter','latex')

% Draw fixed points
yb = ylim;
xb = xlim;
hold on
for n=1:13
    l(n) = line([data_fixed.Tc(n) data_fixed.Tc(n)],yb,'color','r');
    t(n) = text(data_fixed.Tc(n),yb(1) + diff(yb)*0.1,num2str(data_fixed.Tc(n)),'color','r','rotation',90);
end
hold off
xlim(xb)

subplot(3,1,2)
plot(Tin,dataMIDS.VB)
ylabel('Channel B V (V)','interpreter','latex')

subplot(3,1,3)
plot(Tin,Rin)
ylabel('M.C. R (Ohms)','interpreter','latex')
xlabel('Uncalibrated T (K)','interpreter','latex')

% Add data to figure window
save2fig.Tin = Tin;
save2fig.Rin = Rin;
save2fig.VA = dataMIDS.VA;
save2fig.VB = dataMIDS.VB;
save2fig.timeserial = dataMIDS.timeserial;

guidata(hfig3,save2fig)

% Save data to a .mat file
save([num2str(opt) ' ' fileMIDS(1:(end-4))],'save2fig')
disp(['...' num2str(opt) ' ' fileMIDS(1:(end-4)) '.mat saved'])


% % Convert into time elapsed since the start of the measurement (s)
% vec = datevec(rtimeserial);                 % [YY MM DD HR MM SS]
% vec = vec - repmat(vec(1,:),size(vec,1),1); % time in [YY MM DD HR MM SS] elapsed
% % convert into [0 0 SS SS SS SS]
% vec = (vec - repmat(vec(1,:),size(vec,1),1)).*repmat([0 0 86400 3600 60 1],size(vec,1),1);
% timesec = sum(vec,2);

end


function data = load_MIDS(path,file)
% Load the calibration data measured by the CMN and SRD thermometers

curdir = cd;
cd(path)

disp(file)
disp('... loading calibrated thermometers')

% Import the file
newData1 = importdata(file);

data.timeserial = datenum(newData1.textdata(:,1));    % convert date-time
data.timeelapsed = newData1.data(:,1);               % elapsed time (s)

data.VA = newData1.data(:,3);
data.VB = newData1.data(:,5);

% Statistics:
disp(char(['Measured: ' newData1.textdata(1) ' ... ' newData1.textdata(end)]))

cd(curdir)
end


function data = calibrationdata_SRD1000
% Calibration data of SRD 1000 as given in HDL user guide

% Reference point compounds
data.element = { 'W',...
                 'Be',...
                 'Ir80Rh20',...
                 'Ir92Rh08',...
                 'Ir',...
                 'AuAl2',...
                 'AuIn2',...
                 'Cd',...
                 'Zn',...
                 'Al',...
                 'In',...
                 'V',...
                 'Pb'};

% Temperature (K) at which sample is sufficiently in the s.c. state            
data.Tsc = [13.5...
            19.5...
            30 ...
            64 ...
            90 ...
            150 ...
            204 ...
            512 ...
            840 ...
            1170 ...
            3400 ...
            4800 ...
            7170]*1e-3;

% Temperature (K) at which sample is sufficiently in the normal state  
data.Tnc = [18.0 ...
            21.5 ...
            38 ...
            68 ...
            95 ...
            158 ...
            212 ...
            518 ...
            855 ...
            1190 ...
            3430 ...
            5000 ...
            7230]*1e-3;

% Reference temperature (K) of the s.c. phase transition        
data.Tc = [ 14.14 ...
            20.60 ...
            35.24 ...
            65.97 ...
            92.73 ...
            155.28 ...
            207.68 ...
            515.57 ...
            848.66 ...
            1182 ...
            3418 ...
            4925 ...
            7205]*1e-3;

% Width of the s.c. phase transition (K)      
data.Wc = [ 0.5 ...
            0.3 ...
            1.0 ...
            1.1 ...
            1.3 ...
            0.7 ...
            0.9 ...
            1.0 ...
            2.1 ...
            1.2 ...
            3 ...
            23 ...
            6]*1e-3;     

% Estimate of the uncertainty of determination of Tc (K)    
data.eTc = [0.11 ...
            0.10 ...
            0.26 ...
            0.22 ...
            0.28 ...
            0.18 ...
            0.29 ...
            1.00 ...
            0.48 ...
            1.6 ...
            6.2 ...
            8.4 ...
            8.9]*1e-3;
end


function data = load_kelvinox(path,file)
% Load the .csv file which logs the Kelvinox dilution fridge state

curdir = cd;
cd(path)

disp(file)
disp('... loading dilution fridge state')

% Read worksheet
[ndata, text, ~] = xlsread([path filesep file]);

disp(['File stored at: ' char(text(1))])

data.variables = text(2,1:64);              % headers of parameters recorded
data.recorded = ndata(2:end,1:64);          % matrix containing all data


%% Parameterise the data measured =========================================
data.npoint = data.recorded(:,2);


% Convert epoch into matlab timestamp -------------------------------------
% offset (serial date number for 1/1/1970)
dnOffset = datenum('01-Jan-1970 01:00:00');         % GMT -> CET
% assuming it's read in as a string originally
tstamp = num2str(data.recorded(:,3));
% convert to a number, dived by number of seconds
% per day, and add offset
dnNow = str2num(tstamp)/(24*60*60) + dnOffset;

data.timeserial = dnNow;
% -------------------------------------------------------------------------

% Pressure guages (Bar)
data.pressure.P2 = data.recorded(:,4);              % Bar
data.pressure.P1 = data.recorded(:,5);              % Bar
data.pressure.P5 = data.recorded(:,6);              % Bar
data.pressure.P3 = data.recorded(:,7)*1e-3;         % Bar
data.pressure.P4 = data.recorded(:,8)*1e-3;         % Bar

% 1K plate readings
data.diag.Kplate.t = data.recorded(:,9);            % time (s)
data.diag.Kplate.T = data.recorded(:,10);           % temperature (K)
data.diag.Kplate.R = data.recorded(:,11);           % resistance (Ohms)

% Sorb readings
data.diag.Sorb.t = data.recorded(:,12);             % time (s)
data.diag.Sorb.T = data.recorded(:,13);             % temperature (K)
data.diag.Sorb.R = data.recorded(:,14);             % resistance (Ohms)

% Cold plate readings
data.diag.Coldplate.t = data.recorded(:,15);        % time (s)
data.diag.Coldplate.T = data.recorded(:,16);        % temperature (K)
data.diag.Coldplate.R = data.recorded(:,17);        % resistance (Ohms)

% Still readings
data.diag.Still.t = data.recorded(:,18);            % time (s)
data.diag.Still.T = data.recorded(:,19);            % temperature (K)
data.diag.Still.R = data.recorded(:,20);            % resistance (Ohms)

% Mixing chamber readings
data.diag.MC.t = data.recorded(:,21);               % time (s)
data.diag.MC.T = data.recorded(:,22);               % temperature (K)
data.diag.MC.R = data.recorded(:,23);               % resistance (Ohms)

% Heaters
data.heaters.still = data.recorded(:,57);           % Still heater (W)
data.heaters.chamber = data.recorded(:,58);         % Chamber heater (W)
data.heaters.ivc = data.recorded(:,59);             % IVC heater (W)

% State of the turbo pump
data.turbo.current = data.recorded(:,60);           % current (A)
data.turbo.power = data.recorded(:,61);             % power (W)
data.turbo.speed = data.recorded(:,62);             % speed (Hz)
data.turbo.motor = data.recorded(:,63);             % motor (C)
data.turbo.bottom = data.recorded(:,64);            % bottom (C)    

cd(curdir)
end


function data = cutMIDS(data,rtimeserial)
% Function to trim data to within the bounds of the serial time
% MIDS data structure
%
% dataMIDS = 
% 
%      timeserial: [Xx1 double]
%     timeelapsed: [Xx1 double]
%              VA: [Xx1 double]
%              VB: [Xx1 double]

ind = data.timeserial > rtimeserial(1) & data.timeserial < rtimeserial(2);

data.timeserial = data.timeserial(ind);
data.timeelapsed = data.timeelapsed(ind);
data.VA = data.VA(ind);
data.VB = data.VB(ind);


end


function data = cutKelvinox(data,rtimeserial)
% Function to trim data to within the bounds of the serial time
% Kelvinox data structure
%
% dataKelvinox = 
% 
%      variables: {1x64 cell}
%       recorded: [Xx64 double]
%         npoint: [Xx1 double]
%     timeserial: [Xx1 double]
%       pressure: [1x1 struct]
%           diag: [1x1 struct]
%        heaters: [1x1 struct]
%          turbo: [1x1 struct]

ind = data.timeserial > rtimeserial(1) & data.timeserial < rtimeserial(2);

data.npoint = data.npoint(ind);
data.timeserial = data.timeserial(ind);
data.pressure.P1 = data.pressure.P1(ind);
data.pressure.P2 = data.pressure.P2(ind);
data.pressure.P3 = data.pressure.P3(ind);
data.pressure.P4 = data.pressure.P4(ind);
data.pressure.P5 = data.pressure.P5(ind);

data.diag.Kplate.t = data.diag.Kplate.t(ind);
data.diag.Kplate.T = data.diag.Kplate.T(ind);
data.diag.Kplate.R = data.diag.Kplate.R(ind);

data.diag.Sorb.t = data.diag.Sorb.t(ind);
data.diag.Sorb.T = data.diag.Sorb.T(ind);
data.diag.Sorb.R = data.diag.Sorb.R(ind);

data.diag.Coldplate.t = data.diag.Coldplate.t(ind);
data.diag.Coldplate.T = data.diag.Coldplate.T(ind);
data.diag.Coldplate.R = data.diag.Coldplate.R(ind);

data.diag.Still.t = data.diag.Still.t(ind);
data.diag.Still.T = data.diag.Still.T(ind);
data.diag.Still.R = data.diag.Still.R(ind);

data.diag.MC.t = data.diag.MC.t(ind);
data.diag.MC.T = data.diag.MC.T(ind);
data.diag.MC.R = data.diag.MC.R(ind);

data.heaters.still = data.heaters.still(ind);
data.heaters.chamber = data.heaters.chamber(ind);
data.heaters.ivc = data.heaters.ivc(ind);

data.turbo.current = data.turbo.current(ind);
data.turbo.power = data.turbo.power(ind);
data.turbo.speed = data.turbo.speed(ind);
data.turbo.motor = data.turbo.motor(ind);
data.turbo.bottom = data.turbo.bottom(ind);

end











