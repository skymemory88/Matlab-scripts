function Field_scan

% matlabpool open

%addpath('C:\Users\babkevic\Documents\MATLAB\legendflex')

curdir = cd;
cd('C:\Users\yiyang\Google Drive\File sharing\Programming scripts\Matlab\Fieldscan\functions')
%cd('/Users/yikaiyang/Google Drive/File sharing/PhD projects/LiHoF4/Data/Matlab/Cavity resonator/D24mm_T5mm_G0.2mm/12.02.2019/functions')
%The first line is for windows, the second line is for mac OS
%% Define input ===========================================================
% Figure plot options:
plotopt.col(1,:) = [0.2 0.2 0.7];
plotopt.col(2,:) = [0.7 0.2 0.7];
plotopt.col(3,:) = [0.2 0.7 0.2];
plotopt.col(4,:) = [0.7 0.2 0.2];
plotopt.col(5,:) = [0.2 0.7 0.7];
plotopt.col(6,:) = [0.7 0.7 0.2];
plotopt.col(7,:) = [0.7 0.7 0.7];
plotopt.col(8,:) = [0.5 0.9 0.5];
plotopt.col(9,:) = [0.9 0.5 0.5];
plotopt.lnwd = 2;
plotopt.ftsz = 12;
plotopt.mksz = 5;


%% Read ZVL

opt  = 7;

switch opt
    case 1
        % Test of frequency scans using ZVL as a function of field, ramp
        % rate 0.1 T/min and -10 dBm input power
        option1(plotopt)
    case 2
        % Test at 600 mK using different input powers
        option2(plotopt)
    case 3
        % Examine the anomaly at 2 dBm
        option3(plotopt)
    case 4
        % Field scans at different temperatures
        option4(plotopt)
        
        
    % Coil measurements ===================================================
    case 5        
        % Field sweeps at constant temperature
        option5(plotopt)        
    case 6
        % temperature sweeps in 0 and 0.5 T
        option6(plotopt)
    case 7
        % Analyse the fieldscan (case 5) data to extract quality factor value
        option7(plotopt)
        
    % Measurements using snake 1.7 GHz resonator ==========================
    case 8
        % Field scans at different temperatures
        option8(plotopt)
    case 9
        % Temperature scans at different fields
        option9(plotopt)
        
    % Measurements using snake 4.4 GHz resonator ==========================
    case 10
        % Field scans at different temperatures
        option10(plotopt)
        
    % Measurements using 5.5 GHz resonator ==========================
    case 11
        % Field scans at different temperatures
        option11(plotopt)
        
    % Measurements using 3.925 GHz resonator ==========================
    case 12
        % Field scans at different temperatures
        option12(plotopt)
        
    case 999
        % Testing colour plots
        option999(plotopt)
end





% matlabpool close

end

%% ------------------------------------------------------------------------

function option1(plotopt)

filepath = 'G:\My Drive\File sharing\PhD projects\Transverse Ising model\Data\Experiment\Cavity resonator\D24mm_T5mm_G0.2mm\26.04.2019';
%filepath = '/Users/yikaiyang/Google Drive/File sharing/PhD projects/LiHoF4/Data/Matlab/Cavity resonator/D24mm_T5mm_G0.2mm/12.02.2019';
%The first line is for windows, the second line is for mac OS
filename = '2019_04_0013';
run = 12;
nop = 31; % Number of points per segment for the frequency scan

out = readdata_v3(filepath,filename,run);


%% Plot data
freq = out.data.ZVLfreq;
S11 = out.data.ZVLreal + 1i*out.data.ZVLimag;
dB = 20*log10(abs(S11));

freq = freq(:)/1e9;
S11 = S11(:);
dB = dB(:);


H = out.data.DCField1;
T1 = out.data.Temperature1;
T2 = out.data.Temperature2;
HH = repmat(H,1,nop); HH = HH(:); %the third argument is the number of frequency points in each line/segment

FdB = TriScatteredInterp(HH,freq,dB);
% FrS = TriScatteredInterp(HH,freq,real(S11));
% FiS = TriScatteredInterp(HH,freq,imag(S11)); %intrapolate points on a 2D grid

% Plot frequency-field colour map
freq_h = 3.520;
freq_l = 3.540; %set frequency range, l: lower limit, h: higher limit
field_h = 9.0;
field_l = 0.0;  %set field range, l: lower limit, h: higher limit
[xq,yq] = meshgrid(linspace(field_l,field_h,501),linspace(freq_l,freq_h,310)); %set the X and Y range
zq = FdB(xq,yq);

hfig1 = setfig(11);
hp = pcolor(xq,yq,zq);
set(hp, 'edgeColor','none')
colorbar
set(gca,'fontsize',plotopt.ftsz)
t(1) = xlabel('Field (T)');
t(2) = ylabel('Frequency (GHz)');
tt(1) = title([filename num2str(run)]);

% 
% yq = yq.*0 + 3.200; % Plot field scan cut at 3.200 GHz (H vs dB)
% hfig2 = setfig(12);
% zq = FdB(xq,yq);
% h = plot(xq(1,:),zq(1,:));
% xlim([0 9])
% set(gca,'fontsize',plotopt.ftsz)
% t(3) = xlabel('Field (T)');
% t(4) = ylabel('S11 (dB)');
% tt(2) = title([filename num2str(run) ', f = 3.200 GHz']);

% hfig3 = setfig(13); % Plot field scan cut at 3.200 GHz (H vs real(S))
% zq = FrS(xq,yq);
% h = plot(xq(1,:),zq(1,:));
% xlim([0 9])
% set(gca,'fontsize',plotopt.ftsz)
% t(5) = xlabel('Field (T)');
% t(6) = ylabel('Re{S11}');
% tt(3) = title([filename num2str(run) ', f = 3.200 GHz']);
% 
% hfig4 = setfig(14); % Plot field scan cut at 3.200 GHz (H vs real(S))
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
% hfig5 = setfig(15);
% h = plot(yq(:,1),zq1(:,1),yq(:,1),zq2(:,1),yq(:,1),zq3(:,1),yq(:,1),zq4(:,1),yq(:,1),zq5(:,1));
% xlim([3.2 3.7])
% set(gca,'fontsize',plotopt.ftsz)
% t(9) = xlabel('Frequency (GHz)');
% t(10) = ylabel('S11 (dB)');
% tt(5) = title([filename num2str(run)]);
% t(11) = legend(h,'0','2','3.612','3.685','5');

% hfig6 = setfig(16); % Plot field-temperature
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

function option2(plotopt)
filepath = 'V:\2014\11';
filename = '2014_11_';
runs =  [ 53 56 59 62 65  68  71  74  77];
power = [-10 -7 -4 -1  2 -16 -19 -22 -25];

out = readdata_v3(filepath,filename,runs);


%% Plot data
hfig(1) = setfig(21);
hold on
box on
hfig(2) = setfig(22);
hold on
box on
hfig(3) = setfig(23);
hold on
box on
hfig(4) = setfig(24);
hold on
box on

for n = 1:length(out)
    freq = out(n).data.ZVLfreq;
    S11 = out(n).data.ZVLreal + 1i*out(n).data.ZVLimag;
    dB = 20*log10(abs(S11));
    N = size(freq,2);

    freq = freq(:)/1e9;
    S11 = S11(:);
    dB = dB(:);

    H = out(n).data.DCField1;
    T1 = out(n).data.Temperature1;
    T2 = out(n).data.Temperature2;
    HH = repmat(H,1,N); HH = HH(:);
    TT1 = repmat(T1,1,N); TT1 = TT1(:);
    TT2 = repmat(T2,1,N); TT2 = TT2(:);

    FdB = TriScatteredInterp(HH,freq,dB);
    FrS = TriScatteredInterp(HH,freq,real(S11));
    FiS = TriScatteredInterp(HH,freq,imag(S11));
    FTT1 = TriScatteredInterp(HH,freq,TT1);
    FTT2 = TriScatteredInterp(HH,freq,TT2);
    [xq,yq] = meshgrid(linspace(2.5,4.5,501),linspace(3.4,3.5,101));    
    zq = FdB(xq,yq);
    
    % Plot field scan cut at 3.200 GHz (H vs dB)
    yq = yq.*0 + 3.200;

    figure(hfig(1))
    zq = FdB(xq,yq);
    h1(n) = plot(xq(1,:),zq(1,:));     
    ax(1) = gca;
    t(1) = xlabel('Field (T)');
    t(2) = ylabel('S11 (K)');    
    set(h1(n),'color',plotopt.col(n,:),'linewidth',plotopt.lnwd);
    
    figure(hfig(2))
    zq = FrS(xq,yq);
    h2(n) = plot(xq(1,:),zq(1,:));    
    ax(2) = gca;
    t(3) = xlabel('Field (T)');
    t(4) = ylabel('Re{S11} (K)');
    set(h2(n),'color',plotopt.col(n,:),'linewidth',plotopt.lnwd);
    
    figure(hfig(3))
    zq = FiS(xq,yq);
    h3(n) = plot(xq(1,:),zq(1,:));  
    ax(3) = gca;
    t(5) = xlabel('Field (T)');
    t(6) = ylabel('Im{S11} (K)');
    set(h3(n),'color',plotopt.col(n,:),'linewidth',plotopt.lnwd);
    
    figure(hfig(4))
    zq1 = FTT1(xq,yq);
    zq2 = FTT2(xq,yq);
    h4(n,:) = plot(xq(1,:),zq1(1,:),xq(1,:),zq2(1,:));   
    ax(4) = gca;
    t(7) = xlabel('Field (T)');
    t(8) = ylabel('Temperature (K)');
    set(h4(n,:),'color',plotopt.col(n,:),'linewidth',plotopt.lnwd);

end


set(ax,'fontsize',plotopt.ftsz,'xlim',[2.5 4.5])
t(9) = legend(h1,num2str(power'));    
set(t,'fontsize',plotopt.ftsz,'interpreter','latex')
%     tt(2) = title([filename num2str(runs(n)) ', f = 3.200 GHz']);

saveplots(hfig(1),'2-Temperature=0.6K_field-S_powdep_3.200GHz')
saveplots(hfig(2),'2-Temperature=0.6K_field-rS_powdep_3.200GHz')
saveplots(hfig(3),'2-Temperature=0.6K_field-iS_powdep_3.200GHz')
saveplots(hfig(4),'2-Temperature=0.6K_field-T_powdep_3.200GHz')
end

function option3(plotopt)
filepath = 'V:\2014\11';
filename = '2014_11_';
runs =  [65 83 68 85];
% power=  [  2  2 -16 -16];

% out(1) = readdata_v3(filepath,filename,65);
% out(2) = readdata_v3(filepath,filename,85);

out = readdata_v3(filepath,filename,runs);

%% Plot data
hfig(1) = setfig(31);
hold on
box on
hfig(2) = setfig(32);
hold on
box on
hfig(3) = setfig(33);
hold on
box on
hfig(4) = setfig(34);
hold on
box on


[xq,yq] = meshgrid(linspace(2.5,4.5,301),linspace(3.4,3.5,101));  
for n = 1:length(out)
    clear freq S11 dB N FdB FrS FiS FTT1 FTT2
    freq = out(n).data.ZVLfreq;
    S11 = out(n).data.ZVLreal + 1i*out(n).data.ZVLimag;
    dB = 20*log10(abs(S11));
    N = size(freq,2);

    freq = freq(:)/1e9;
    S11 = S11(:);
    dB = dB(:);

    H = out(n).data.DCField1;
    T1 = out(n).data.Temperature1;
    T2 = out(n).data.Temperature2;
    HH = repmat(H,1,N); HH = HH(:);
    TT1 = repmat(T1,1,N); TT1 = TT1(:);
    TT2 = repmat(T2,1,N); TT2 = TT2(:);

    FdB  = TriScatteredInterp(HH,freq,dB);
    FrS  = TriScatteredInterp(HH,freq,real(S11));
    FiS  = TriScatteredInterp(HH,freq,imag(S11));
    FTT1 = TriScatteredInterp(HH,freq,TT1);
    FTT2 = TriScatteredInterp(HH,freq,TT2);
    
    % Plot field scan cut at 3.437 GHz (H vs dB)
    yq = yq.*0 + 3.437;

    figure(hfig(1))
    zq = FdB(xq,yq);
    h1(n) = plot(xq(1,:),zq(1,:));     
    ax(1) = gca;
    t(1) = xlabel('Field (T)');
    t(2) = ylabel('S11 (K)');    
    set(h1(n),'color',plotopt.col(n,:),'linewidth',plotopt.lnwd);
    
    figure(hfig(2))
    zq = FrS(xq,yq);
    h2(n) = plot(xq(1,:),zq(1,:));    
    ax(2) = gca;
    t(3) = xlabel('Field (T)');
    t(4) = ylabel('Re{S11} (K)');
    set(h2(n),'color',plotopt.col(n,:),'linewidth',plotopt.lnwd);
    
    figure(hfig(3))
    zq = FiS(xq,yq);
    h3(n) = plot(xq(1,:),zq(1,:));  
    ax(3) = gca;
    t(5) = xlabel('Field (T)');
    t(6) = ylabel('Im{S11} (K)');
    set(h3(n),'color',plotopt.col(n,:),'linewidth',plotopt.lnwd);
    
    figure(hfig(4))
    zq1 = FTT1(xq,yq);
    zq2 = FTT2(xq,yq);
    h4(n,:) = plot(xq(1,:),zq1(1,:),xq(1,:),zq2(1,:));   
    ax(4) = gca;
    t(7) = xlabel('Field (T)');
    t(8) = ylabel('Temperature (K)');
    set(h4(n,:),'color',plotopt.col(n,:),'linewidth',plotopt.lnwd);

end


set(ax,'fontsize',plotopt.ftsz,'xlim',[2.5 4.5])
t(9) = legend(h1,'2 dBm (0.6 K)','2 dBm (0.6 K)','-16 dBm (0.6 K)','-16 dBm (0.7 K)');    
set(t,'fontsize',plotopt.ftsz,'interpreter','latex')
%     tt(2) = title([filename num2str(runs(n)) ', f = 3.437 GHz']);

saveplots(hfig(1),'3-Temperature=0.6K_field-S_2dBm_3.437GHz')
saveplots(hfig(2),'3-Temperature=0.6K_field-rS_2dBm_3.437GHz')
saveplots(hfig(3),'3-Temperature=0.6K_field-iS_2dBm_3.437GHz')
saveplots(hfig(4),'3-Temperature=0.6K_field-T_2dBm_3.437GHz')
end

function option4(plotopt)
filepath = 'V:\2014\11';
filename = '2014_11_';

ll = 1;
runs{ll} =  90:92;     T(ll) = 0.1;  ll=ll+1;
runs{ll} =  197:199;   T(ll) = 0.2;  ll=ll+1;
runs{ll} =  95:97;     T(ll) = 0.3;  ll=ll+1;
runs{ll} =  202:204;   T(ll) = 0.4;  ll=ll+1;
runs{ll} =  207:209;   T(ll) = 0.6;  ll=ll+1;
runs{ll} =  111:113;   T(ll) = 1.0;  ll=ll+1;
runs{ll} =  137:139;   T(ll) = 1.2;  ll=ll+1;
runs{ll} =  116:118;   T(ll) = 1.5;  ll=ll+1;
runs{ll} =  142:144;   T(ll) = 1.7;  ll=ll+1;
runs{ll} =  147:149;   T(ll) = 2.0;  ll=ll+1;
runs{ll} =  152:154;   T(ll) = 2.5;  ll=ll+1;
runs{ll} =  126:134;   T(ll) = 4.0;  ll=ll+1;

    
%% Plot data
hfig(1) = setfig(41);
hold on
box on
hfig(2) = setfig(42);
hold on
box on
hfig(3) = setfig(43);
hold on
box on
hfig(4) = setfig(44);
hold on
box on

g = gausswin(20); % <-- this value determines the width of the smoothing window
g = g/sum(g);
   
[xq,yq] = meshgrid(linspace(0,9,901),linspace(3.4,3.5,101));  
for n = 1:length(T)
    tout = readdata_v3(filepath,filename,runs{n});
    out(n) = mergeout_v1(tout);
    
    col = setcolours(n/length(T),'jet');
    
    clear freq S11 dB N FdB FrS FiS FTT1 FTT2
    freq = out(n).data.ZVLfreq;
    S11 = out(n).data.ZVLreal + 1i*out(n).data.ZVLimag;
    dB = 20*log10(abs(S11));
    N = size(freq,2);

    freq = freq(:)/1e9;
    S11 = S11(:);
    dB = dB(:);

    H = out(n).data.DCField1;
    T1 = out(n).data.Temperature1;
    T2 = out(n).data.Temperature2;
    HH = repmat(H,1,N); HH = HH(:);
    TT1 = repmat(T1,1,N); TT1 = TT1(:);
    TT2 = repmat(T2,1,N); TT2 = TT2(:);

    FdB  = TriScatteredInterp(HH,freq,dB);
    FrS  = TriScatteredInterp(HH,freq,real(S11));
    FiS  = TriScatteredInterp(HH,freq,imag(S11));
    FTT1 = TriScatteredInterp(HH,freq,TT1);
    FTT2 = TriScatteredInterp(HH,freq,TT2);
    
    % Plot field scan cut at 3.437 GHz (H vs dB)
    yq = yq.*0 + 3.437;

    figure(hfig(1))
    zq = FdB(xq,yq);
    zqf = medfilt1(zq(1,:),10);

    h1(n) = plot(xq(1,:),zqf);     
    ax(1) = gca;
    t(1) = xlabel('Field (T)');
    t(2) = ylabel('S11 (arb.)');    
    set(h1(n),'color',col,'linewidth',plotopt.lnwd);
    
    figure(hfig(2))
    zq = FrS(xq,yq);
    zqf = medfilt1(zq(1,:),10);
    h2(n) = plot(xq(1,:),zqf);    
    ax(2) = gca;
    t(3) = xlabel('Field (T)');
    t(4) = ylabel('Re{S11} (arb.)');
    set(h2(n),'color',col,'linewidth',plotopt.lnwd);
    
    figure(hfig(3))
    zq = FiS(xq,yq);
	zqf = medfilt1(zq(1,:),10);
    h3(n) = plot(xq(1,:),zqf);  
    ax(3) = gca;
    t(5) = xlabel('Field (T)');
    t(6) = ylabel('Im{S11} (arb.)');
    set(h3(n),'color',col,'linewidth',plotopt.lnwd);
    
    figure(hfig(4))
    zq1 = FTT1(xq,yq);
    zq2 = FTT2(xq,yq);
    h4(n,:) = plot(xq(1,:),zq1(1,:),xq(1,:),zq2(1,:));   
    ax(4) = gca;
    t(7) = xlabel('Field (T)');
    t(8) = ylabel('Temperature (K)');
    set(h4(n,:),'color',col,'linewidth',plotopt.lnwd);

end


set(ax,'fontsize',plotopt.ftsz,'xlim',[0 9])
set(ax(2),'ylim',[0.2 0.6])
t(9) = legend(h1,num2str(T(:),'%3.1f K'),'location','EastOutside');    
set(t,'fontsize',plotopt.ftsz,'interpreter','latex')
%     tt(2) = title([filename num2str(runs(n)) ', f = 3.437 GHz']);

saveplots(hfig(1),'4-Temperature=field-S_-16dBm_3.437GHz')
saveplots(hfig(2),'4-Temperature=field-rS_-16dBm_3.437GHz')
saveplots(hfig(3),'4-Temperature=field-iS_-16dBm_3.437GHz')
saveplots(hfig(4),'4-Temperature=field-T_-16dBm_3.437GHz')
end

function option5(plotopt)
filepath = 'G:\My Drive\File sharing\PhD projects\Transverse Ising model\Data\Experiment\Cavity resonator\D24mm_T5mm_G0.2mm\22.04.2019';
filename = '2019_04_0009';

ll = 1;
runs{ll} =  283:285;    T(ll) = 0.1;  ll=ll+1;
runs{ll} =  288:290;    T(ll) = 0.3;  ll=ll+1;
runs{ll} =  293:295;    T(ll) = 0.5;  ll=ll+1;
runs{ll} =  2803:2805;  T(ll) = 0.7;  ll=ll+1;
runs{ll} =  264:266;    T(ll) = 0.8;  ll=ll+1;
runs{ll} =  243:245;    T(ll) = 1.0;  ll=ll+1;
runs{ll} =  248:250;    T(ll) = 1.2;  ll=ll+1;
runs{ll} =  2813:2815;  T(ll) = 1.3;  ll=ll+1;
runs{ll} =  253:255;    T(ll) = 1.4;  ll=ll+1;
runs{ll} =  2818:2820;  T(ll) = 1.45;  ll=ll+1;
runs{ll} =  2823:2825;  T(ll) = 1.5;  ll=ll+1;
runs{ll} =  2828:2830;  T(ll) = 1.55;  ll=ll+1;
%runs{ll} =  259:261;   T(ll) = 1.6;  ll=ll+1;
runs{ll} =  2833:2835;  T(ll) = 1.6;  ll=ll+1;
runs{ll} =  2838:2840;  T(ll) = 1.7;  ll=ll+1;
runs{ll} =  2843:2845;  T(ll) = 1.8;  ll=ll+1;
runs{ll} =  2848:2850;  T(ll) = 3.0;  ll=ll+1;
runs{ll} =  271:273;    T(ll) = 4.0;  ll=ll+1;

    
%% Plot data
hfig(1) = setfig(51);
hold on
box on
hfig(2) = setfig(52);
hold on
box on
hfig(3) = setfig(53);
hold on
box on
hfig(4) = setfig(54);
hold on
box on

fcut = 0.4674;
fac = 10;


for n = 1:length(T)
    [xq,yq] = meshgrid(linspace(0,0.2,951),linspace(0.46,0.475,201)*fac);  
    
    tout = readdata_v3(filepath,filename,runs{n});
    out(n) = mergeout_v1(tout);
    
    col = setcolours(n/length(T),'jet');
    
    clear freq S11 dB N FdB FrS FiS FTT1 FTT2
    freq = out(n).data.ZVLfreq;
    S11 = out(n).data.ZVLreal + 1i*out(n).data.ZVLimag;
    
    
    dB = 20*log10(abs(S11));
    N = size(freq,2);

    freq = freq(:)/1e9;
    S11 = S11(:);
    dB = dB(:);

    H = out(n).data.DCField1;
    T1 = out(n).data.Temperature1;
    T2 = out(n).data.Temperature2;
    HH = repmat(H,1,N); HH = HH(:);
    TT1 = repmat(T1,1,N); TT1 = TT1(:);
    TT2 = repmat(T2,1,N); TT2 = TT2(:);

    FdB  = TriScatteredInterp(HH,freq*fac,dB);
    FrS  = TriScatteredInterp(HH,freq*fac,real(S11));
    FiS  = TriScatteredInterp(HH,freq*fac,imag(S11));
    FTT1 = TriScatteredInterp(HH,freq*fac,TT1);
    FTT2 = TriScatteredInterp(HH,freq*fac,TT2);
    
    % Plot colour map of freq-field scans
    zq = FdB(xq,yq);
    hfig_map = setfig(0);    
    hp = pcolor(xq,yq/fac,zq);
    set(hp, 'edgeColor','none')
    colorbar
    set(gca,'fontsize',plotopt.ftsz)
    txt = xlabel('Field (T)');
    tyt = ylabel('Frequency (GHz)');
    ttt = title(num2str(T(n),'$T = %3.2f$ K'));
    set([txt tyt ttt],'fontsize',plotopt.ftsz,'interpreter','latex')
    saveplots(hfig_map,['5-Temperature=' num2str(T(n)) 'K_field-S_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
    
    % Plot field scan cut at 0.467 GHz (H vs dB)
    yq = yq.*0 + fcut*fac;

    figure(hfig(1))
    zq = FdB(xq,yq);
    zqf = medfilt1(zq(1,:),10);
    
    h1(n) = plot(xq(1,:),zqf);     
    ax(1) = gca;
    tx(1) = xlabel('Field (T)');
    ty(1) = ylabel('S11 (arb.)');    
    tt(1) = title([num2str(T(n),'$T = %3.2f$ K, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h1(n),'color',col,'linewidth',plotopt.lnwd);
    
    figure(hfig(2))
    zq = FrS(xq,yq);    
    zqf = medfilt1(zq(1,:),10);
    
    h2(n) = plot(xq(1,:),zqf(1,:));    
    ax(2) = gca;
    tx(2) = xlabel('Field (T)');
    ty(2) = ylabel('Re{S11} (arb.)');
    tt(2) = title([num2str(T(n),'$T = %3.2f$ K, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h2(n),'color',col,'linewidth',plotopt.lnwd);
    
    figure(hfig(3))
    zq = FiS(xq,yq);
    zqf = medfilt1(zq(1,:),10);
    
    h3(n) = plot(xq(1,:),zqf(1,:));  
    ax(3) = gca;
    tx(3) = xlabel('Field (T)');
    ty(3) = ylabel('Im{S11} (arb.)');
    tt(3) = title([num2str(T(n),'$T = %3.2f$ K, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h3(n),'color',col,'linewidth',plotopt.lnwd);
    
    figure(hfig(4))
    zq1 = FTT1(xq,yq);
    zqf1 = medfilt1(zq1(1,:),10);
    zq2 = FTT2(xq,yq);
    zqf2 = medfilt1(zq2(1,:),10);
    h4(n,:) = plot(xq(1,:),zqf1(1,:),xq(1,:),zqf2(1,:));   
    ax(4) = gca;
    tx(4) = xlabel('Field (T)');
    ty(4) = ylabel('Temperature (K)');
    tt(4) = title([num2str(T(n),'$T = %3.2f$ K, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h4(n,:),'color',col,'linewidth',plotopt.lnwd);

    T_MC(n) = mean(zq1(~isnan(zq1)));
    T_sample(n) = mean(zq2(~isnan(zq2)));
end


set(ax,'fontsize',plotopt.ftsz,'xlim',[-0.5 9])
set(ax(1),'ylim',[-13.5 -10.5])
set(ax(2),'ylim',[-0.02 0.2])
set(ax(3),'ylim',[0.19 0.3])
set(ax(4),'ylim',[0 4.5])

for nn = 1:length(T_sample)
    lbl{nn} = num2str(T_sample(nn),'%3.2f K');
end

% % hl(4) = legend(h1, lbl,'location','sw');
% 
% [dum h] = legendflex(h1, lbl, 'ref', hfig(1), ...
%                        'anchor', {'s','s'}, ...
%                        'buffer',[80 80], ...
%                        'nrow',6, ...
%                        'fontsize',plotopt.ftsz,...
%                        'box','on',...
%                        'interpreter','latex');
                   
tl = legend(h1,num2str(T(:),'%3.1f K'),'location','EastOutside');    
set([tx ty tl tt],'fontsize',plotopt.ftsz,'interpreter','latex')


saveplots(hfig(1),['5-Temperature_field-S_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
saveplots(hfig(2),['5-Temperature_field-rS_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
saveplots(hfig(3),['5-Temperature_field-iS_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
saveplots(hfig(4),['5-Temperature_field-T_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])



end

function option6(plotopt)
filepath = 'V:\2014\11';
filename = '2014_11_';

ll = 1;
% runs{ll} =  269;   Hx(ll) = 0;      ll=ll+1;
% runs{ll} =  276;   Hx(ll) = 0.5;    ll=ll+1;
runs{ll} =  2853;   Hx(ll) = 0;         ll=ll+1;
runs{ll} =  2856;   Hx(ll) = 0.5;       ll=ll+1;
runs{ll} =  2859;   Hx(ll) = 1;         ll=ll+1;
runs{ll} =  2862;   Hx(ll) = 1.5;       ll=ll+1;
runs{ll} =  2865;   Hx(ll) = 2;         ll=ll+1;
runs{ll} =  2868;   Hx(ll) = 2.5;       ll=ll+1;
runs{ll} =  2871;   Hx(ll) = 3;         ll=ll+1;
runs{ll} =  2874;   Hx(ll) = 5;         ll=ll+1;

    
%% Plot data
hfig(1) = setfig(61);
hold on
box on
hfig(2) = setfig(62);
hold on
box on
hfig(3) = setfig(63);
hold on
box on

fcut = 0.4674;
fac = 10;

for n = 1:length(Hx)
    [xq,yq] = meshgrid(linspace(0.8,4,501),linspace(0.46,0.475,201)*fac);  
    
    tout = readdata_v3(filepath,filename,runs{n});
    out(n) = mergeout_v1(tout);
    
    col = setcolours(n/length(Hx),'jet');
    
    clear freq S11 dB N FdB FrS FiS FTT1 FTT2
    freq = out(n).data.ZVLfreq;
    S11 = out(n).data.ZVLreal + 1i*out(n).data.ZVLimag;
    dB = 20*log10(abs(S11));
    N = size(freq,2);

    freq = freq(:)/1e9;
    S11 = S11(:);
    dB = dB(:);

    H = out(n).data.DCField1;
    T1 = out(n).data.Temperature1;
    T2 = out(n).data.Temperature2;
    HH = repmat(H,1,N); HH = HH(:);
    TT1 = repmat(T1,1,N); TT1 = TT1(:);
    TT2 = repmat(T2,1,N); TT2 = TT2(:);

    FdB  = TriScatteredInterp(TT2,freq*fac,dB);
    FrS  = TriScatteredInterp(TT2,freq*fac,real(S11));
    FiS  = TriScatteredInterp(TT2,freq*fac,imag(S11));
%     FTT1 = TriScatteredInterp(HH,freq*fac,TT1);
%     FTT2 = TriScatteredInterp(HH,freq*fac,TT2);
    
    % Plot colour map of freq-field scans
    zq = FdB(xq,yq);
    hfig_map = setfig(1545);
    hp = pcolor(xq,yq/fac,zq);
    set(hp, 'edgeColor','none')
    colorbar
    set(gca,'fontsize',plotopt.ftsz)
    txt = xlabel('Temperature (K)');
    tyt = ylabel('Frequency (GHz)');
    ttt = title(num2str(Hx(n),'$H = %3.2f$ T'));    
    set([txt tyt ttt],'fontsize',plotopt.ftsz,'interpreter','latex')
    saveplots(hfig_map,['6-Field=' num2str(Hx(n)) 'T_temperature-S_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
    
    % Plot field scan cut at 0.4674 GHz (H vs dB)
    yq = yq.*0 + fcut*fac;

    figure(hfig(1))
    zq = FdB(xq,yq);
    zqf = medfilt1(zq(1,:),10);
    
    h1(n) = plot(xq(1,:),zqf);     
    ax(1) = gca;
    tx(1) = xlabel('Temperature(sample) (K)');
    ty(1) = ylabel('S11 (arb.)');    
    tt(1) = title([num2str(Hx(n),'$H = %3.2f$ T, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h1(n),'color',col,'linewidth',plotopt.lnwd);
    
    figure(hfig(2))
    zq = FrS(xq,yq);   
    zqf = medfilt1(zq(1,:),10);
    
    h2(n) = plot(xq(1,:),zqf(1,:));    
    ax(2) = gca;
    tx(2) = xlabel('Temperature(sample) (K)');
    ty(2) = ylabel('Re{S11} (arb.)');
    tt(2) = title([num2str(Hx(n),'$H = %3.2f$ T, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h2(n),'color',col,'linewidth',plotopt.lnwd);
    
    figure(hfig(3))
    zq = FiS(xq,yq);
    zqf = medfilt1(zq(1,:),10);
    
    h3(n) = plot(xq(1,:),zqf(1,:));  
    ax(3) = gca;
    tx(3) = xlabel('Temperature(sample) (K)');
    ty(3) = ylabel('Im{S11} (arb.)');
    tt(3) = title([num2str(Hx(n),'$H = %3.2f$ T, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h3(n),'color',col,'linewidth',plotopt.lnwd);
    

end


set(ax,'fontsize',plotopt.ftsz,'xlim',[0.8 2])
% set(ax(2),'ylim',[0.2 0.6])
tl = legend(h1,num2str(Hx(:),'%3.1f T'),'location','nw');    
set([tx ty tl tt],'fontsize',plotopt.ftsz,'interpreter','latex')
%     tt(2) = title([filename num2str(runs(n)) ', f = 3.437 GHz']);

saveplots(hfig(1),['6-Field_temperature-S_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
saveplots(hfig(2),['6-Field_temperature-rS_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
saveplots(hfig(3),['6-Field_temperature-iS_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])

end

function option7(plotopt)


filepath = 'C:\Users\yiyang\Google Drive\File sharing\PhD projects\Transverse Ising model\Data\Experiment\LiHoF4\22.04.2019';
filename = '2019_04_0009';

ll = 1;
runs{ll} =  1;   T(ll) = 0.3;
% runs{ll} =  264:264;   T(ll) = 0.3;  ll=ll+1;
% runs{ll} =  243:245;   T(ll) = 1.0;  ll=ll+1;
% runs{ll} =  248:250;   T(ll) = 1.2;  ll=ll+1;
% runs{ll} =  253:255;   T(ll) = 1.4;  ll=ll+1;
% runs{ll} =  259:261;   T(ll) = 1.6;  ll=ll+1;
% runs{ll} =  271:273;   T(ll) = 4.0;  ll=ll+1;

    
%% Plot data
% hfig(1) = setfig(71);
% hold on
% box on

fac = 1.0;
freq_h = 3.520;
freq_l = 3.490; %set frequency range, l: lower limit, h: higher limit
field_h = 9.0;
field_l = 0.0;  %set field range, l: lower limit, h: higher limit

for n = 1:length(T)
    [xq,yq] = meshgrid(linspace(field_l,field_h,901),linspace(freq_l,freq_h,310*fac));
    
    tout = readdata_v3(filepath,filename,runs{n});
    out(n) = mergeout_v1(tout);
    
    col = setcolours(n/length(T),'jet');
    
    clear freq S11 dB N FdB FrS FiS FTT1 FTT2
    freq = out(n).data.ZVLfreq;
    S11 = out(n).data.ZVLreal + 1i*out(n).data.ZVLimag;
    H = out(n).data.DCField1;
    %     T1 = out(n).data.Temperature1;
    dB = mag2db(abs(S11));
    f0 = zeros(length(dB),1);
    %Collect all the resonant peaks
    for i = 1:length(dB)
        [~,idx] = min(dB(i,:));
        f0(i,1) = freq(i,idx)/1e9;
    end
      
    %     N = length(yq(:,1));  % Option_1 Number of elements of the interpolated frequency data
    N = size(freq,2);   % Option_2 Number of elements of the measured frequency data
    
    freq = freq(:)/1e9;
    %     S11 = S11(:);
    mag = (1- out(n).data.ZVLreal.^2 -out(n).data.ZVLimag.^2)./((1- out(n).data.ZVLreal).^2 + out(n).data.ZVLimag.^2);
    mag = mag(:); %Unsure what these two lines are meant to do
   
%     % Option_1 Interpolate data along only the frequency axis.
%     interp_dB = zeros(length(H),N);
%     for i = 1:length(H)
%         interp_dB(i,:) = interp1(freq(i,:),dB(i,:),yq(:,1));
%     end
%     interp_dB = interp_dB(:);
%     yy = repmat(yq(:,1),1,length(H));
%     yy = yy(:);
     
    %     T2 = out(n).data.Temperature2;    %If there are two sets of temperature data
    %     TT1 = repmat(T1,1,N);  %populate the temperature to match the dimension of S11 matrix
    HH = repmat(H,1,N); %populate the magnetic field to match the dimension of S11 matrix
%     dB = dB(:);
    HH = HH(:);
    %     TT1 = TT1(:);
    %     TT2 = repmat(T2,1,N); TT2 = TT2(:);
    
    % Option_2 Interpolate the data along all (two) axis.
%     FdB  = scatteredInterpolant(HH,freq,dB);
    Fmag  = scatteredInterpolant(HH,freq*fac,mag);
    %     FrS  = scatteredInterpolant(HH,freq*fac,real(S11));
    %     FiS  = scatteredInterpolant(HH,freq*fac,imag(S11));
    %     FTT1 = scatteredInterpolant(HH,freq*fac,TT1);
    %     FTT2 = scatteredInterpolant(HH,freq*fac,TT2);
    
    % option_1 plot surface interpolated data using pseudo-colormap
%     zq = FdB(xq,yq);
%     %     hfig_map(n) = setfig(n);
%     hp = pcolor(xq,yq/fac,zq);
%     set(hp, 'edgeColor','none')
%     shading interp;
    
    %     option_2.1 plot single-direction interpolated data using pseudo-colormap
%     hp = pcolor(H,yq(:,1),interp_dB');
%     set(hp, 'edgeColor','none')
   
    %option_2.2 plot single-direction interpolated data using scatter plot
%     C = linspace(min(interp_dB),max(interp_dB),length(HH));
%     hp = scatter3(HH,yy,interp_dB,2,C,'o','filled','MarkerEdgeColor','none');   
%     colormap(hsv);

%     colorbar
%     set(gca,'fontsize',plotopt.ftsz)
%     axis([field_l field_h freq_l freq_h]);
%     tx(n+4) = xlabel('Field (T)');
%     ty(n+4) = ylabel('Frequency (GHz)');
%     tt(n+4) = title(num2str(T(n),'S11 response at T = %3.3f K'));
    
    figure
    axis([field_l field_h freq_l freq_h]);
    plot(H,f0,'-o','MarkerSize',3);
    tx(n+5) = xlabel('Field (T)');
    ty(n+5) = ylabel('Resonant frequency (GHz)');
    tt(n+5) = title(num2str(T(n),'Resonant frequency at T = %3.3f K'));
    
    figure
    disp(num2str(T(n),'Fitting: T = %3.2f K'))
    % Extract frequency cuts
    ll = 1;
    Hx = 0.1:0.1:8.9;
    for m = Hx
        xq = xq.*0 + m;
        zq = Fmag(xq,yq);
        zqf = medfilt1(zq(:,1),10);
        s = spec1d(yq(:,1)/fac,zqf,zqf.*0 + 0.05);
        p = [-1 0.467 0.001 1.2];
        fix = [1 1 1 1];
        [fQ fbck]=fits(s,'lorz',p,fix);
        Q = fbck.pvals(2)/fbck.pvals(3);
        chi(ll) = 1/Q;
        ll = ll + 1;
    end
    
    yq = yq.*0 + 0.349;
    zq = Fmag(xq,yq);
    figure(hfig(1))
    h1(n) = plot(xq(1,:),zq);
    ax(1) = gca;
    tx(1) = xlabel('Field (T)');
    ty(1) = ylabel('Quality factor');
    set(h1(n),'color',col,'linewidth',plotopt.lnwd);
end

set(ax,'fontsize',plotopt.ftsz,'xlim',[field_l field_h])
% set(ax(2),'ylim',[0.2 0.6])
% tl = legend(h1,num2str(T(:),'%3.1f K'),'location','se');
% set([tx ty tl tt],'fontsize',plotopt.ftsz,'interpreter','latex')
%     tt(2) = title([filename num2str(runs(n)) ', f = 3.437 GHz']);

% saveplots(hfig(1),'7-Temperature_field-S_-16dBm_Q')



end

function option8(plotopt)

filepath = 'V:\2014\11';
filename = '2014_11_';

year = '2014';

ll = 1;
runs{ll} =  2918:2920;   T(ll) = 0.1;   mon{ll} = '11';     ll=ll+1; 
runs{ll} =  2933:2935;   T(ll) = 0.3;   mon{ll} = '11';     ll=ll+1;
runs{ll} =  2943:2945;   T(ll) = 0.6;   mon{ll} = '11';     ll=ll+1;
runs{ll} =  2880:2882;   T(ll) = 0.85;  mon{ll} = '11';     ll=ll+1;
runs{ll} =  2959:2961;   T(ll) = 1;     mon{ll} = '11';     ll=ll+1;
runs{ll} =  2964:2966;   T(ll) = 1.25;  mon{ll} = '11';     ll=ll+1;
runs{ll} =  3015:3017;   T(ll) = 1.35;  mon{ll} = '11';     ll=ll+1;
runs{ll} =  2969:2971;   T(ll) = 1.45;  mon{ll} = '11';     ll=ll+1;
runs{ll} =  2974:2976;   T(ll) = 1.5;   mon{ll} = '11';     ll=ll+1;
runs{ll} =  2979:2981;   T(ll) = 1.55;  mon{ll} = '11';     ll=ll+1;
runs{ll} =  2984:2986;   T(ll) = 1.6;   mon{ll} = '11';     ll=ll+1;
runs{ll} =  2989:2991;   T(ll) = 1.8;   mon{ll} = '11';     ll=ll+1;
runs{ll} =  3020:3022;   T(ll) = 2.2;   mon{ll} = '11';     ll=ll+1;
runs{ll} =  3:5;         T(ll) = 2.6;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  2994:2996;   T(ll) = 3.0;   mon{ll} = '11';     ll=ll+1;
runs{ll} =  8:10;        T(ll) = 4.0;   mon{ll} = '12';     ll=ll+1;
    
%% Plot data
hfig(1) = setfig(81);
hold on
box on
hfig(2) = setfig(82);
hold on
box on
hfig(3) = setfig(83);
hold on
box on
hfig(4) = setfig(84);
hold on
box on

fcut = 1.682;

fac = 10;

for n = 1:length(T)
    [xq,yq] = meshgrid(linspace(-0.5,9,951),linspace(1.63,1.73,201)*fac);  
    
    tout = readdata_v3(['V:\' year '\' mon{n}],[year '_' mon{n} '_'],runs{n});
    out(n) = mergeout_v1(tout);
    
    col = setcolours(n/length(T),'jet');
    
    clear freq S11 dB N FdB FrS FiS FTT1 FTT2
    freq = out(n).data.ZVLfreq;
    S11 = out(n).data.ZVLreal + 1i*out(n).data.ZVLimag;
        
    dB = 20*log10(abs(S11));
    N = size(freq,2);

    freq = freq(:)/1e9;
    S11 = S11(:);
    dB = dB(:);

    H = out(n).data.DCField1;
    T1 = out(n).data.Temperature1;
    T2 = out(n).data.Temperature2;
    HH = repmat(H,1,N); HH = HH(:);
    TT1 = repmat(T1,1,N); TT1 = TT1(:);
    TT2 = repmat(T2,1,N); TT2 = TT2(:);

    FdB  = TriScatteredInterp(HH,freq*fac,dB);
    FrS  = TriScatteredInterp(HH,freq*fac,real(S11));
    FiS  = TriScatteredInterp(HH,freq*fac,imag(S11));
    FTT1 = TriScatteredInterp(HH,freq*fac,TT1);
    FTT2 = TriScatteredInterp(HH,freq*fac,TT2);
    
    % Plot colour map of freq-field scans
    zq = FdB(xq,yq);
    hfig_map = setfig(9787);    
    hp = pcolor(xq,yq/fac,zq);
    set(hp, 'edgeColor','none')
    colorbar
    set(gca,'fontsize',plotopt.ftsz)
    txt = xlabel('Field (T)');
    tyt = ylabel('Frequency (GHz)');
    ttt = title(num2str(T(n),'$T = %3.2f$ K'));
    set([txt tyt ttt],'fontsize',plotopt.ftsz,'interpreter','latex')
    saveplots(hfig_map,...
        ['8-Temperature=' num2str(T(n)) 'K_field-S_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
    
    % Plot field scan cut at 0.467 GHz (H vs dB)
    yq = yq.*0 + fcut*fac;

    figure(hfig(1))
    zq = FdB(xq,yq);
    zqf = medfilt1(zq(1,:),10);
    
    h1(n) = plot(xq(1,:),zqf);     
    ax(1) = gca;
    tx(1) = xlabel('Field (T)');
    ty(1) = ylabel('S11 (arb.)');    
    tt(1) = title([num2str(T(n),'$T = %3.2f$ K, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h1(n),'color',col,'linewidth',plotopt.lnwd);
    
    figure(hfig(2))
    zq = FrS(xq,yq);    
    zqf = medfilt1(zq(1,:),10);
    
    h2(n) = plot(xq(1,:),zqf(1,:));    
    ax(2) = gca;
    tx(2) = xlabel('Field (T)');
    ty(2) = ylabel('Re{S11} (arb.)');
    tt(2) = title([num2str(T(n),'$T = %3.2f$ K, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h2(n),'color',col,'linewidth',plotopt.lnwd);
    
    figure(hfig(3))
    zq = FiS(xq,yq);
    zqf = medfilt1(zq(1,:),10);
    
    h3(n) = plot(xq(1,:),zqf(1,:));  
    ax(3) = gca;
    tx(3) = xlabel('Field (T)');
    ty(3) = ylabel('Im{S11} (arb.)');
    tt(3) = title([num2str(T(n),'$T = %3.2f$ K, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h3(n),'color',col,'linewidth',plotopt.lnwd);
    
    figure(hfig(4))
    zq1 = FTT1(xq,yq);
    zqf1 = medfilt1(zq1(1,:),10);
    zq2 = FTT2(xq,yq);
    zqf2 = medfilt1(zq2(1,:),10);
    h4(n,:) = plot(xq(1,:),zqf1(1,:),xq(1,:),zqf2(1,:));   
    ax(4) = gca;
    tx(4) = xlabel('Field (T)');
    ty(4) = ylabel('Temperature (K)');
    tt(4) = title([num2str(T(n),'$T = %3.2f$ K, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h4(n,:),'color',col,'linewidth',plotopt.lnwd);   
    T_MC(n) = mean(zq1(~isnan(zq1)));
    T_sample(n) = mean(zq2(~isnan(zq2)));

end


set(ax,'fontsize',plotopt.ftsz,'xlim',[-0.5 9])
% set(ax(1),'ylim',[-10 -6])
% set(ax(2),'ylim',[0.25 0.45])
% set(ax(3),'ylim',[0.17 0.25])
% set(ax(4),'ylim',[0 4.5])
% tl = legend(h1,num2str(T_sample(:),'%3.2f K'),'location','se');    

for nn = 1:length(T_sample)
    lbl{nn} = num2str(T_sample(nn),'%3.2f K');
end

% % hl(4) = legend(h1, lbl,'location','sw');
% 
[dum h] = legendflex(h1, lbl, 'ref', hfig(1), ...
                       'anchor', {'s','s'}, ...
                       'buffer',[80 80], ...
                       'nrow',6, ...
                       'fontsize',plotopt.ftsz,...
                       'box','on',...
                       'interpreter','latex');
                   
                   
set([tx ty tt],'fontsize',plotopt.ftsz,'interpreter','latex')
% set([tl],'fontsize',plotopt.ftsz-2,'interpreter','latex')

saveplots(hfig(1),['8-Temperature_field-S_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
saveplots(hfig(2),['8-Temperature_field-rS_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
saveplots(hfig(3),['8-Temperature_field-iS_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
saveplots(hfig(4),['8-Temperature_field-T_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])


end

function option9(plotopt)
filepath = 'V:\2014\11';
filename = '2014_11_';

ll = 1;
runs{ll} =  2999;   Hx(ll) = 0;      ll=ll+1;
runs{ll} =  3002;   Hx(ll) = 1;      ll=ll+1;
runs{ll} =  3005;   Hx(ll) = 2;      ll=ll+1;
runs{ll} =  3008;   Hx(ll) = 2.5;    ll=ll+1;
runs{ll} =  3011;   Hx(ll) = 3.0;    ll=ll+1;

%% Plot data
hfig(1) = setfig(91);
hold on
box on
hfig(2) = setfig(92);
hold on
box on
hfig(3) = setfig(93);
hold on
box on

fcut = 1.682;

fac = 10;

for n = 1:length(Hx)
    [xq,yq] = meshgrid(linspace(0.8,2,501),linspace(1.63,1.73,201)*fac);  
    
    tout = readdata_v3(filepath,filename,runs{n});
    out(n) = mergeout_v1(tout);
    
    col = setcolours(n/length(Hx),'jet');
    
    clear freq S11 dB N FdB FrS FiS FTT1 FTT2
    freq = out(n).data.ZVLfreq;
    S11 = out(n).data.ZVLreal + 1i*out(n).data.ZVLimag;
    dB = 20*log10(abs(S11));
    N = size(freq,2);

    freq = freq(:)/1e9;
    S11 = S11(:);
    dB = dB(:);

    H = out(n).data.DCField1;
    T1 = out(n).data.Temperature1;
    T2 = out(n).data.Temperature2;
    HH = repmat(H,1,N); HH = HH(:);
    TT1 = repmat(T1,1,N); TT1 = TT1(:);
    TT2 = repmat(T2,1,N); TT2 = TT2(:);

    FdB  = TriScatteredInterp(TT2,freq*fac,dB);
    FrS  = TriScatteredInterp(TT2,freq*fac,real(S11));
    FiS  = TriScatteredInterp(TT2,freq*fac,imag(S11));
    
    % Plot colour map of freq-field scans
    zq = FdB(xq,yq);
    hfig_map = setfig(3243);
    hp = pcolor(xq,yq/fac,zq);
    set(hp, 'edgeColor','none')
    colorbar
    set(gca,'fontsize',plotopt.ftsz)
    txt = xlabel('Temperature(sample) (K)');
    tyt = ylabel('Frequency (GHz)');
    ttt = title(num2str(Hx(n),'$H = %3.2f$ T'));
    set([txt tyt ttt],'fontsize',plotopt.ftsz,'interpreter','latex')
    saveplots(hfig_map,['6-Field=' num2str(Hx(n)) 'T_temperature-S_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])

    
    % Plot field scan cut (H vs dB)
    yq = yq.*0 + fcut*fac;

    figure(hfig(1))
    zq = FdB(xq,yq);
    zqf = medfilt1(zq(1,:),10);
    
    h1(n) = plot(xq(1,:),zqf);     
    ax(1) = gca;
    tx(1) = xlabel('Temperature(sample) (K)');
    ty(1) = ylabel('S11 (arb.)');   
    tt(1) = title([num2str(Hx(n),'$H = %3.2f$ T, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h1(n),'color',col,'linewidth',plotopt.lnwd);
    
    figure(hfig(2))
    zq = FrS(xq,yq);
    zqf = medfilt1(zq(1,:),10);
    
    h2(n) = plot(xq(1,:),zqf(1,:));    
    ax(2) = gca;
    tx(2) = xlabel('Temperature(sample) (K)');
    ty(2) = ylabel('Re{S11} (arb.)');
    tt(2) = title([num2str(Hx(n),'$H = %3.2f$ T, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h2(n),'color',col,'linewidth',plotopt.lnwd);
    
    figure(hfig(3))
    zq = FiS(xq,yq);
    zqf = medfilt1(zq(1,:),10);
    
    h3(n) = plot(xq(1,:),zqf(1,:));  
    ax(3) = gca;
    tx(3) = xlabel('Temperature(sample) (K)');
    ty(3) = ylabel('Im{S11} (arb.)');
    tt(3) = title([num2str(Hx(n),'$H = %3.2f$ T, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h3(n),'color',col,'linewidth',plotopt.lnwd);
    
end


set(ax,'fontsize',plotopt.ftsz,'xlim',[0.8 2])
% set(ax(2),'ylim',[0.2 0.6])
tl = legend(h1,num2str(Hx(:),'%3.1f T'),'location','sw');    
set([tx ty tl tt],'fontsize',plotopt.ftsz,'interpreter','latex')

saveplots(hfig(1),['9-Field_temperature-S_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
saveplots(hfig(2),['9-Field_temperature-rS_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
saveplots(hfig(3),['9-Field_temperature-iS_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])


end

function option10(plotopt)

filepath = 'V:\2014\11';
filename = '2014_12_';

year = '2014';

ll = 1;
runs{ll} =  76:78;   T(ll) = 0.01;  mon{ll} = '12';     ll=ll+1;
runs{ll} =  81:83;   T(ll) = 0.1;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  86:88;   T(ll) = 0.2;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  91:93;   T(ll) = 0.3;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  96:98;   T(ll) = 0.4;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  101:103; T(ll) = 0.5;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  106:108; T(ll) = 0.6;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  111:113; T(ll) = 0.7;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  23:25;   T(ll) = 0.9;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  28:30;   T(ll) = 1;     mon{ll} = '12';     ll=ll+1;
runs{ll} =  58:60;   T(ll) = 1.1;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  33:35;   T(ll) = 1.2;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  63:65;   T(ll) = 1.3;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  38:40;   T(ll) = 1.4;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  68:70;   T(ll) = 1.5;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  43:45;   T(ll) = 1.6;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  48:50;   T(ll) = 1.8;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  53:55;   T(ll) = 3.0;   mon{ll} = '12';     ll=ll+1;


%% Plot data
hfig(1) = setfig(81);
hold on
box on
hfig(2) = setfig(82);
hold on
box on
hfig(3) = setfig(83);
hold on
box on
hfig(4) = setfig(84);
hold on
box on

fcut = 4.449;

fac = 10;

for n = 1:length(T)
    [xq,yq] = meshgrid(linspace(-0.5,9,951),linspace(4.4,4.5,201)*fac);  
    
    tout = readdata_v3(['V:\' year '\' mon{n}],[year '_' mon{n} '_'],runs{n});
    out(n) = mergeout_v1(tout);
    
    col = setcolours(n/length(T),'jet');
    
    clear freq S11 dB N FdB FrS FiS FTT1 FTT2
    freq = out(n).data.ZVLfreq;
    S11 = out(n).data.ZVLreal + 1i*out(n).data.ZVLimag;
        
    dB = 20*log10(abs(S11));
    N = size(freq,2);

    freq = freq(:)/1e9;
    S11 = S11(:);
    dB = dB(:);

    H = out(n).data.DCField1;
    T1 = out(n).data.Temperature1;
    T2 = out(n).data.Temperature2;
    HH = repmat(H,1,N); HH = HH(:);
    TT1 = repmat(T1,1,N); TT1 = TT1(:);
    TT2 = repmat(T2,1,N); TT2 = TT2(:);

    FdB  = TriScatteredInterp(HH,freq*fac,dB);
    FrS  = TriScatteredInterp(HH,freq*fac,real(S11));
    FiS  = TriScatteredInterp(HH,freq*fac,imag(S11));
    FTT1 = TriScatteredInterp(HH,freq*fac,TT1);
    FTT2 = TriScatteredInterp(HH,freq*fac,TT2);
    
    % Plot colour map of freq-field scans
    zq = FdB(xq,yq);
    hfig_map = setfig(9787);    
    hp = pcolor(xq,yq/fac,zq);
    set(hp, 'edgeColor','none')
    colorbar
    set(gca,'fontsize',plotopt.ftsz)
    txt = xlabel('Field (T)');
    tyt = ylabel('Frequency (GHz)');
    ttt = title(num2str(T(n),'$T = %3.2f$ K'));
    set([txt tyt ttt],'fontsize',plotopt.ftsz,'interpreter','latex')
    saveplots(hfig_map,...
        ['10-Temperature=' num2str(T(n)) 'K_field-S_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
    
    % Plot field scan cut at 0.467 GHz (H vs dB)
    yq = yq.*0 + fcut*fac;

    figure(hfig(1))
    zq = FdB(xq,yq);
    zqf = medfilt1(zq(1,:),10);
    
    h1(n) = plot(xq(1,:),zqf);     
    ax(1) = gca;
    tx(1) = xlabel('Field (T)');
    ty(1) = ylabel('S11 (arb.)');    
    tt(1) = title([num2str(T(n),'$T = %3.2f$ K, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h1(n),'color',col,'linewidth',plotopt.lnwd);
    
    figure(hfig(2))
    zq = FrS(xq,yq);    
    zqf = medfilt1(zq(1,:),10);
    
    h2(n) = plot(xq(1,:),zqf(1,:));    
    ax(2) = gca;
    tx(2) = xlabel('Field (T)');
    ty(2) = ylabel('Re{S11} (arb.)');
    tt(2) = title([num2str(T(n),'$T = %3.2f$ K, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h2(n),'color',col,'linewidth',plotopt.lnwd);
    
    figure(hfig(3))
    zq = FiS(xq,yq);
    zqf = medfilt1(zq(1,:),10);
    
    h3(n) = plot(xq(1,:),zqf(1,:));  
    ax(3) = gca;
    tx(3) = xlabel('Field (T)');
    ty(3) = ylabel('Im{S11} (arb.)');
    tt(3) = title([num2str(T(n),'$T = %3.2f$ K, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h3(n),'color',col,'linewidth',plotopt.lnwd);
    
    figure(hfig(4))
    zq1 = FTT1(xq,yq);
    zqf1 = medfilt1(zq1(1,:),10);
    zq2 = FTT2(xq,yq);
    zqf2 = medfilt1(zq2(1,:),10);
    h4(n,:) = plot(xq(1,:),zqf1(1,:),xq(1,:),zqf2(1,:));   
    ax(4) = gca;
    tx(4) = xlabel('Field (T)');
    ty(4) = ylabel('Temperature (K)');
    tt(4) = title([num2str(T(n),'$T = %3.2f$ K, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h4(n,:),'color',col,'linewidth',plotopt.lnwd);   
    T_MC(n) = mean(zq1(~isnan(zq1)));
    T_sample(n) = mean(zq2(~isnan(zq2)));

end


set(ax,'fontsize',plotopt.ftsz,'xlim',[-0.5 9])
% set(ax(1),'ylim',[-10 -6])
% set(ax(2),'ylim',[0.25 0.45])
% set(ax(3),'ylim',[0.17 0.25])
% set(ax(4),'ylim',[0 4.5])
% tl = legend(h1,num2str(T_sample(:),'%3.2f K'),'location','se');    

for nn = 1:length(T_sample)
    lbl{nn} = num2str(T_sample(nn),'%3.2f K');
end

% % hl(4) = legend(h1, lbl,'location','sw');
% 
[dum1 h] = legendflex(h1, lbl, 'ref', hfig(1), ...
                       'anchor', {'s','s'}, ...
                       'buffer',[80 80], ...
                       'nrow',6, ...
                       'fontsize',plotopt.ftsz,...
                       'box','on',...
                       'xscale',0.3,...
                       'interpreter','none');
                   
[dum2 h] = legendflex(h1, lbl, 'ref', hfig(2), ...
                       'anchor', {'n','n'}, ...
                       'buffer',[140 -80], ...
                       'nrow',6, ...
                       'fontsize',plotopt.ftsz,...
                       'box','on',...
                       'xscale',0.3,...
                       'interpreter','none');

[dum3 h] = legendflex(h1, lbl, 'ref', hfig(3), ...
                       'anchor', {'s','s'}, ...
                       'buffer',[80 80], ...
                       'nrow',6, ...
                       'fontsize',plotopt.ftsz,...
                       'box','on',...
                       'xscale',0.3,...
                       'interpreter','none');
                   
set([tx ty tt],'fontsize',plotopt.ftsz,'interpreter','latex')
% set([tl],'fontsize',plotopt.ftsz-2,'interpreter','latex')

saveplots(hfig(1),['10-Temperature_field-S_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
saveplots(hfig(2),['10-Temperature_field-rS_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
saveplots(hfig(3),['10-Temperature_field-iS_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
saveplots(hfig(4),['10-Temperature_field-T_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])


end

function option11(plotopt)

filepath = 'V:\2014\11';
filename = '2014_12_';

year = '2014';

ll = 1;
runs{ll} =  183:185;   T(ll) = 0.01;  mon{ll} = '12';     ll=ll+1;
runs{ll} =  188:190;   T(ll) = 0.1;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  193:195;   T(ll) = 0.2;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  198:200;   T(ll) = 0.3;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  203:205;   T(ll) = 0.4;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  208:210;   T(ll) = 0.5;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  213:215;   T(ll) = 0.6;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  218:220;   T(ll) = 0.7;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  120:122;   T(ll) = 0.9;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  125:127;   T(ll) = 1;     mon{ll} = '12';     ll=ll+1;
runs{ll} =  130:132;   T(ll) = 1.1;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  135:137;   T(ll) = 1.2;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  140:142;   T(ll) = 1.3;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  145:147;   T(ll) = 1.4;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  150:152;   T(ll) = 1.5;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  155:157;   T(ll) = 1.6;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  160:162;   T(ll) = 1.7;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  165:167;   T(ll) = 1.8;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  170:172;   T(ll) = 2.4;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  175:177;   T(ll) = 3.0;   mon{ll} = '12';     ll=ll+1;


%% Plot data
hfig(1) = setfig(111);
hold on
box on
hfig(2) = setfig(112);
hold on
box on
hfig(3) = setfig(113);
hold on
box on
hfig(4) = setfig(114);
hold on
box on

fcut = 5.601;

fac = 10;

for n = 1:length(T)
    [xq,yq] = meshgrid(linspace(-0.5,9,951),linspace(5.55,5.65,201)*fac);  
    
    tout = readdata_v3(['V:\' year '\' mon{n}],[year '_' mon{n} '_'],runs{n});
    out(n) = mergeout_v1(tout);
    
    col = setcolours(n/length(T),'jet');
    
    clear freq S11 dB N FdB FrS FiS FTT1 FTT2
    freq = out(n).data.ZVLfreq;
    S11 = out(n).data.ZVLreal + 1i*out(n).data.ZVLimag;
        
    dB = 20*log10(abs(S11));
    N = size(freq,2);

    freq = freq(:)/1e9;
    S11 = S11(:);
    dB = dB(:);

    H = out(n).data.DCField1;
    T1 = out(n).data.Temperature1;
    T2 = out(n).data.Temperature2;
    HH = repmat(H,1,N); HH = HH(:);
    TT1 = repmat(T1,1,N); TT1 = TT1(:);
    TT2 = repmat(T2,1,N); TT2 = TT2(:);

    FdB  = TriScatteredInterp(HH,freq*fac,dB);
    FrS  = TriScatteredInterp(HH,freq*fac,real(S11));
    FiS  = TriScatteredInterp(HH,freq*fac,imag(S11));
    FTT1 = TriScatteredInterp(HH,freq*fac,TT1);
    FTT2 = TriScatteredInterp(HH,freq*fac,TT2);
    
    % Plot colour map of freq-field scans
    zq = FdB(xq,yq);
    hfig_map = setfig(9787);    
    hp = pcolor(xq,yq/fac,zq);
    set(hp, 'edgeColor','none')
    colorbar
    set(gca,'fontsize',plotopt.ftsz)
    txt = xlabel('Field (T)');
    tyt = ylabel('Frequency (GHz)');
    ttt = title(num2str(T(n),'$T = %3.2f$ K'));
    set([txt tyt ttt],'fontsize',plotopt.ftsz,'interpreter','latex')
    saveplots(hfig_map,...
        ['11-Temperature=' num2str(T(n)) 'K_field-S_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
    
    % Plot field scan cut at 0.467 GHz (H vs dB)
    yq = yq.*0 + fcut*fac;

    figure(hfig(1))
    zq = FdB(xq,yq);
    zqf = medfilt1(zq(1,:),10);
    
    h1(n) = plot(xq(1,:),zqf);     
    ax(1) = gca;
    tx(1) = xlabel('Field (T)');
    ty(1) = ylabel('S11 (arb.)');    
    tt(1) = title([num2str(T(n),'$T = %3.2f$ K, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h1(n),'color',col,'linewidth',plotopt.lnwd);
    
    figure(hfig(2))
    zq = FrS(xq,yq);    
    zqf = medfilt1(zq(1,:),10);
    
    h2(n) = plot(xq(1,:),zqf(1,:));    
    ax(2) = gca;
    tx(2) = xlabel('Field (T)');
    ty(2) = ylabel('Re{S11} (arb.)');
    tt(2) = title([num2str(T(n),'$T = %3.2f$ K, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h2(n),'color',col,'linewidth',plotopt.lnwd);
    
    figure(hfig(3))
    zq = FiS(xq,yq);
    zqf = medfilt1(zq(1,:),10);
    
    h3(n) = plot(xq(1,:),zqf(1,:));  
    ax(3) = gca;
    tx(3) = xlabel('Field (T)');
    ty(3) = ylabel('Im{S11} (arb.)');
    tt(3) = title([num2str(T(n),'$T = %3.2f$ K, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h3(n),'color',col,'linewidth',plotopt.lnwd);
    
    figure(hfig(4))
    zq1 = FTT1(xq,yq);
    zqf1 = medfilt1(zq1(1,:),10);
    zq2 = FTT2(xq,yq);
    zqf2 = medfilt1(zq2(1,:),10);
    h4(n,:) = plot(xq(1,:),zqf1(1,:),xq(1,:),zqf2(1,:));   
    ax(4) = gca;
    tx(4) = xlabel('Field (T)');
    ty(4) = ylabel('Temperature (K)');
    tt(4) = title([num2str(T(n),'$T = %3.2f$ K, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h4(n,:),'color',col,'linewidth',plotopt.lnwd);   
    T_MC(n) = mean(zq1(~isnan(zq1)));
    T_sample(n) = mean(zq2(~isnan(zq2)));

end


set(ax,'fontsize',plotopt.ftsz,'xlim',[-0.5 9])
% set(ax(1),'ylim',[-10 -6])
% set(ax(2),'ylim',[0.25 0.45])
set(ax(3),'ylim',[-0.6 0])
% set(ax(4),'ylim',[0 4.5])
% tl = legend(h1,num2str(T_sample(:),'%3.2f K'),'location','se');    

for nn = 1:length(T_sample)
    lbl{nn} = num2str(T_sample(nn),'%3.2f K');
end

set([tx ty tt],'fontsize',plotopt.ftsz,'interpreter','latex')

figure(hfig(1))
[dum1 h] = legendflex(h1, lbl, 'ref', hfig(1), ...
                       'anchor', {'n','n'}, ...
                       'buffer',[135 -80], ...
                       'nrow',7, ...
                       'fontsize',plotopt.ftsz,...
                       'box','on',...
                       'xscale',0.3,...
                       'interpreter','none');

figure(hfig(2))                   
[dum2 h] = legendflex(h1, lbl, 'ref', hfig(2), ...
                       'anchor', {'s','s'}, ...
                       'buffer',[80 80], ...
                       'nrow',7, ...
                       'fontsize',plotopt.ftsz,...
                       'box','on',...
                       'xscale',0.3,...
                       'interpreter','none');

figure(hfig(3))                   
[dum3 h] = legendflex(h1, lbl, 'ref', hfig(3), ...
                       'anchor', {'s','s'}, ...
                       'buffer',[140 80], ...
                       'nrow',7, ...
                       'fontsize',plotopt.ftsz,...
                       'box','on',...
                       'xscale',0.3,...
                       'interpreter','none');
                   

saveplots(hfig(1),['11-Temperature_field-S_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
saveplots(hfig(2),['11-Temperature_field-rS_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
saveplots(hfig(3),['11-Temperature_field-iS_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
saveplots(hfig(4),['11-Temperature_field-T_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])


end

function option12(plotopt)

filepath = 'G:\My Drive\File sharing\PhD projects\LiHoF4\Data\Matlab\Cavity resonator\D24mm_T5mm_G0.2mm\06.02.2019';
filename = '2019_02_0026';

year = '2019';

ll = 1;
runs{ll} =  280:282;   T(ll) = 0.01;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  285:287;   T(ll) = 0.2;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  290:292;   T(ll) = 0.3;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  295:297;   T(ll) = 0.4;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  300:302;   T(ll) = 0.5;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  305:307;   T(ll) = 0.6;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  310:312;   T(ll) = 0.7;   mon{ll} = '12';     ll=ll+1;

runs{ll} =  240:242;   T(ll) = 0.8;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  245:247;   T(ll) = 1;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  250:252;   T(ll) = 1.2;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  255:257;   T(ll) = 1.4;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  260:262;   T(ll) = 1.5;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  265:267;   T(ll) = 1.6;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  270:272;   T(ll) = 2;   mon{ll} = '12';     ll=ll+1;
runs{ll} =  275:277;   T(ll) = 3;   mon{ll} = '12';     ll=ll+1;


%% Plot data
hfig(1) = setfig(121);
hold on
box on
hfig(2) = setfig(122);
hold on
box on
hfig(3) = setfig(123);
hold on
box on
hfig(4) = setfig(124);
hold on
box on

set(hfig,'visible','off')

fcut = 3.923;

fac = 10;

for n = 1:length(T)
    [xq,yq] = meshgrid(linspace(-0.5,9,951),linspace(3.875,3.975,201)*fac);  
    
    tout = readdata_v3(filepath,filename,runs{n});
    out(n) = mergeout_v1(tout);
    
    col = setcolours(n/length(T),'jet');
    
    clear freq S11 dB N FdB FrS FiS FTT1 FTT2
    freq = out(n).data.ZVLfreq;
    S11 = out(n).data.ZVLreal + 1i*out(n).data.ZVLimag;
        
    dB = 20*log10(abs(S11));
    N = size(freq,2);

    freq = freq(:)/1e9;
    S11 = S11(:);
    dB = dB(:);

    H = out(n).data.DCField1;
    T1 = out(n).data.Temperature1;
    T2 = out(n).data.Temperature2;
    HH = repmat(H,1,N); HH = HH(:);
    TT1 = repmat(T1,1,N); TT1 = TT1(:);
    TT2 = repmat(T2,1,N); TT2 = TT2(:);

    FdB  = TriScatteredInterp(HH,freq*fac,dB);
    FrS  = TriScatteredInterp(HH,freq*fac,real(S11));
    FiS  = TriScatteredInterp(HH,freq*fac,imag(S11));
    FTT1 = TriScatteredInterp(HH,freq*fac,TT1);
    FTT2 = TriScatteredInterp(HH,freq*fac,TT2);
    
    % Plot colour map of freq-field scans
    zq = FdB(xq,yq);
    hfig_map = setfig(9787);    
    hp = pcolor(xq,yq/fac,zq);
    set(hp, 'edgeColor','none')
    colorbar
    set(gca,'fontsize',plotopt.ftsz)
    txt = xlabel('Field (T)');
    tyt = ylabel('Frequency (GHz)');
    ttt = title(num2str(T(n),'$T = %3.2f$ K'));
    set([txt tyt ttt],'fontsize',plotopt.ftsz,'interpreter','latex')
    saveplots(hfig_map,...
        ['12-Temperature=' num2str(T(n)) 'K_field-S_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
    
    
    % Plot field scan cut at 0.467 GHz (H vs dB)
    yq = yq.*0 + fcut*fac;

    figure(hfig(1))
    zq = FdB(xq,yq);
    zqf = medfilt1(zq(1,:),10);
    
    h1(n) = plot(xq(1,:),zqf);     
    ax(1) = gca;
    tx(1) = xlabel('Field (T)');
    ty(1) = ylabel('S11 (arb.)');    
    tt(1) = title([num2str(T(n),'$T = %3.2f$ K, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h1(n),'color',col,'linewidth',plotopt.lnwd);
    set(hfig,'visible','off')
    
    figure(hfig(2))
    zq = FrS(xq,yq);    
    zqf = medfilt1(zq(1,:),10);
    
    h2(n) = plot(xq(1,:),zqf(1,:));    
    ax(2) = gca;
    tx(2) = xlabel('Field (T)');
    ty(2) = ylabel('Re{S11} (arb.)');
    tt(2) = title([num2str(T(n),'$T = %3.2f$ K, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h2(n),'color',col,'linewidth',plotopt.lnwd);
    set(hfig,'visible','off')
    
    figure(hfig(3))
    zq = FiS(xq,yq);
    zqf = medfilt1(zq(1,:),10);
    
    h3(n) = plot(xq(1,:),zqf(1,:));  
    ax(3) = gca;
    tx(3) = xlabel('Field (T)');
    ty(3) = ylabel('Im{S11} (arb.)');
    tt(3) = title([num2str(T(n),'$T = %3.2f$ K, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h3(n),'color',col,'linewidth',plotopt.lnwd);
    set(hfig,'visible','off')
    
    figure(hfig(4))
    zq1 = FTT1(xq,yq);
    zqf1 = medfilt1(zq1(1,:),10);
    zq2 = FTT2(xq,yq);
    zqf2 = medfilt1(zq2(1,:),10);
    h4(n,:) = plot(xq(1,:),zqf1(1,:),xq(1,:),zqf2(1,:));   
    ax(4) = gca;
    tx(4) = xlabel('Field (T)');
    ty(4) = ylabel('Temperature (K)');
    tt(4) = title([num2str(T(n),'$T = %3.2f$ K, ') num2str(fcut,'%3.4f') ' GHz']);
    set(h4(n,:),'color',col,'linewidth',plotopt.lnwd);   
    T_MC(n) = mean(zq1(~isnan(zq1)));
    T_sample(n) = mean(zq2(~isnan(zq2)));
    set(hfig,'visible','off')
end

set(hfig,'visible','on')

set(ax,'fontsize',plotopt.ftsz,'xlim',[-0.5 9])
% set(ax(1),'ylim',[-10 -6])
% set(ax(2),'ylim',[0.25 0.45])
set(ax(3),'ylim',[-1 -0.3])
% set(ax(4),'ylim',[0 4.5])
% tl = legend(h1,num2str(T_sample(:),'%3.2f K'),'location','se');    

for nn = 1:length(T_sample)
    lbl{nn} = num2str(T_sample(nn),'%3.2f K');
end

set([tx ty tt],'fontsize',plotopt.ftsz,'interpreter','latex')

figure(hfig(1))
[dum1 h] = legendflex(h1, lbl, 'ref', hfig(1), ...
                       'anchor', {'n','n'}, ...
                       'buffer',[100 -80], ...
                       'nrow',4, ...
                       'fontsize',plotopt.ftsz,...
                       'box','on',...
                       'xscale',0.3,...
                       'interpreter','none');

figure(hfig(2))                   
[dum2 h] = legendflex(h1, lbl, 'ref', hfig(2), ...
                       'anchor', {'s','s'}, ...
                       'buffer',[100 80], ...
                       'nrow',4, ...
                       'fontsize',plotopt.ftsz,...
                       'box','on',...
                       'xscale',0.3,...
                       'interpreter','none');

figure(hfig(3))                   
[dum3 h] = legendflex(h1, lbl, 'ref', hfig(3), ...
                       'anchor', {'s','s'}, ...
                       'buffer',[100 80], ...
                       'nrow',4, ...
                       'fontsize',plotopt.ftsz,...
                       'box','on',...
                       'xscale',0.3,...
                       'interpreter','none');
                   

saveplots(hfig(1),['12-Temperature_field-S_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
saveplots(hfig(2),['12-Temperature_field-rS_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
saveplots(hfig(3),['12-Temperature_field-iS_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])
saveplots(hfig(4),['12-Temperature_field-T_-16dBm_' num2str(fcut,'%3.4f') 'GHz'])


end



function option999(plotopt)

filepath = 'V:\2014\11';
filename = '2014_11_';
runs =  [ 283:285];

tout = readdata_v3(filepath,filename,runs);
out = mergeout_v1(tout);

%% Plot data
% hfig(1) = setfig(51);
% hold on
% box on
% hfig(2) = setfig(52);
% hold on
% box on
% hfig(3) = setfig(53);
% hold on
% box on
% hfig(4) = setfig(54);
% hold on
% box on

for n = 1:length(out)
    freq = out(n).data.ZVLfreq;
    S11 = out(n).data.ZVLreal + 1i*out(n).data.ZVLimag;
    dB = 20*log10(abs(S11));
    N = size(freq,2);

    [a,b] = histc(freq(:,1),unique(freq(:,1)));
    ind = find(diff(freq(:,1))<0);
    steps_freq = max(b);            

    H = out(n).data.DCField1;       HH = repmat(H,1,N); 
    T1 = out(n).data.Temperature1;  TT1 = repmat(T1,1,N); 
    T2 = out(n).data.Temperature2;  TT2 = repmat(T2,1,N); 
                
    % Remove unfinished initial sweep of frequency
    freq((ind(end)+1):end,:) = [];
    dB(ind(end)+1:end,:) = [];
    S11(ind(end)+1:end,:) = [];
    HH(ind(end)+1:end,:) = [];
    TT1(ind(end)+1:end,:) = [];
    TT2(ind(end)+1:end,:) = [];
    
    freq(1:ind(1),:) = [];
    dB(1:ind(1),:) = [];
    S11(1:ind(1),:) = [];
    HH(1:ind(1),:) = [];
    TT1(1:ind(1),:) = [];
    TT2(1:ind(1),:) = [];        
    
    NL = size(freq,1);
       
    for mm=1:NL/steps_freq
        ind2 = (steps_freq*(mm-1) + 1):(steps_freq*(mm));
        t1 = freq(ind2,:);
        t2 = dB(ind2,:);
        t3 = S11(ind2,:);
        t4 = HH(ind2,:);
        t5 = TT1(ind2,:);
        t6 = TT2(ind2,:);
        
        t1 = t1(:)'; t2 = t2(:)'; t3 = t3(:)'; t4 = t4(:)'; t5 = t5(:)'; t6 = t6(:)';
        
        [dum inds] = sort(t1);
        t1 = t1(inds); t2 = t2(inds); t3 = t3(inds); t4 = t4(inds); t5 = t5(inds); t6 = t6(inds);
        
        freq_s(mm,:) = t1/1e9;
        dB_s(mm,:) = t2;
        S11_s(mm,:) = t3;
        HH_s(mm,:) = t4;
        TT1_s(mm,:) = t5;
        TT2_s(mm,:) = t6;
    end
   
    meth = 'nearest';
    FdB = TriScatteredInterp(HH_s(:),freq_s(:)*10,dB_s(:),meth);
    FrS = TriScatteredInterp(HH_s(:),freq_s(:)*10,real(S11_s(:)),meth);
    FiS = TriScatteredInterp(HH_s(:),freq_s(:)*10,imag(S11_s(:)),meth);
    FTT1 = TriScatteredInterp(HH_s(:),freq_s(:)*10,TT1_s(:),meth);
    FTT2 = TriScatteredInterp(HH_s(:),freq_s(:)*10,TT2_s(:),meth);
    [xq,yq] = meshgrid(linspace(0,9,901),linspace(4.6,4.75,101));    
    zq = FdB(xq,yq);
    
    
    hfig1 = setfig(40);
    hp = pcolor(HH_s,freq_s,dB_s);
    set(hp, 'edgeColor','none')
    colorbar
    set(gca,'fontsize',plotopt.ftsz)
    t(1) = xlabel('Field (T)');
    t(2) = ylabel('Frequency (GHz)');
    
    % Plot frequency-field colour map
    hfig1 = setfig(50);
    hp = pcolor(xq,yq,zq);
    set(hp, 'edgeColor','none')
    colorbar
    set(gca,'fontsize',plotopt.ftsz)
    t(1) = xlabel('Field (T)');
    t(2) = ylabel('Frequency (GHz)');
%     tt(1) = title([filename num2str(run)]);

    hfig1 = setfig(60);
%     F = TriScatteredInterp(x,y,z);
    z2i=FdB(xq,yq);
    tri = delaunay(xq,yq);
    ph = trisurf(tri,xq,yq,z2i);


    % Plot field scan cut at 0.467 GHz (H vs dB)
    yq = yq.*0 + 0.467;

    figure(hfig(1))
    zq = FdB(xq,yq);
    h1(n) = plot(xq(1,:),zq(1,:));     
    ax(1) = gca;
    t(1) = xlabel('Field (T)');
    t(2) = ylabel('S11 (K)');    
    set(h1(n),'color',plotopt.col(n,:),'linewidth',plotopt.lnwd);
    
    figure(hfig(2))
    zq = FrS(xq,yq);
    h2(n) = plot(xq(1,:),zq(1,:));    
    ax(2) = gca;
    t(3) = xlabel('Field (T)');
    t(4) = ylabel('Re{S11} (K)');
    set(h2(n),'color',plotopt.col(n,:),'linewidth',plotopt.lnwd);
    
    figure(hfig(3))
    zq = FiS(xq,yq);
    h3(n) = plot(xq(1,:),zq(1,:));  
    ax(3) = gca;
    t(5) = xlabel('Field (T)');
    t(6) = ylabel('Im{S11} (K)');
    set(h3(n),'color',plotopt.col(n,:),'linewidth',plotopt.lnwd);
    
    figure(hfig(4))
    zq1 = FTT1(xq,yq);
    zq2 = FTT2(xq,yq);
    h4(n,:) = plot(xq(1,:),zq1(1,:),xq(1,:),zq2(1,:));   
    ax(4) = gca;
    t(7) = xlabel('Field (T)');
    t(8) = ylabel('Temperature (K)');
    set(h4(n,:),'color',plotopt.col(n,:),'linewidth',plotopt.lnwd);

end


set(ax,'fontsize',plotopt.ftsz,'xlim',[2.5 4.5])  
set(t,'fontsize',plotopt.ftsz,'interpreter','latex')
%     tt(2) = title([filename num2str(runs(n)) ', f = 3.437 GHz']);

saveplots(hfig(1),'5-Temperature=0.6K_field-S_powdep_3.437GHz')
saveplots(hfig(2),'5-Temperature=0.6K_field-rS_powdep_3.437GHz')
saveplots(hfig(3),'5-Temperature=0.6K_field-iS_powdep_3.437GHz')
saveplots(hfig(4),'5-Temperature=0.6K_field-T_powdep_3.437GHz')


end

%% ------------------------------------------------------------------------

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


