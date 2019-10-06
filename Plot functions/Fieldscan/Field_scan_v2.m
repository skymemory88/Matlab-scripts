function Field_scan
format long;  
% matlabpool open

%addpath('C:\Users\babkevic\Documents\MATLAB\legendflex')

% % curdir = cd;
addpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\Fieldscan\functions');
addpath(genpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\spec1d--Henrik'));

% addpath('/Users/yikaiyang/Google Drive/File sharing/Programming scripts/Matlab/Fieldscan/functions')
% addpath(genpath('/Users/yikaiyang/Google Drive/File sharing/Programming scripts/Matlab/spec1d--Henrik'));
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
%Set data range and parameters
opt  = 2;

switch opt
    case 1
        % Test of frequency scans using ZVL as a function of field
        option1(plotopt)
    case 2
        % Analyse the fieldscan (case 5) data to extract quality factor value
        option2(plotopt)
end
% matlabpool close

end

%% ------------------------------------------------------------------------

function option1(plotopt)

filepath = 'G:\My Drive\File sharing\PhD projects\LiHoF4 project\Data\Experiment\LiHoF4\SC127\05.09.2019';
% filepath = '/Users/yikaiyang/Google Drive/File sharing/PhD projects/LiHoF4 project/Data/Experiment/LiHoF4/SC127/03.09.2019/';
%The first line is for windows, the second line is for mac OS
filename = '2019_09_0010';
run = 1;
out = readdata_v3(filepath,filename,run);

%Set data range and parameters
freq_l = 3.630; %set frequency range, l: lower limit, h: higher limit
freq_h = 3.690;
field_l = 1.5;  %set field range, l: lower limit, h: higher limit
field_h = 4.0;
nop = 31; % Number of points per segment for the frequency scan according to ZVL network analyzer
step = 1e-4; %Scanning step in unit of GHz,default is 1 MHz according to ZVL network analyser
lineDiv = floor((freq_h-freq_l)/((nop - 1) * step)); %compute how many lines per complete frequency scan.

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

% Could use "scatteredInterpolant()" to replace "TriScatteredInterp()" as recommended by MATLAB, but it may generate artifacts
FdB = TriScatteredInterp(HH,freq,dB);
% FrS = TriScatteredInterp(HH,freq,real(S11));
% FiS = TriScatteredInterp(HH,freq,imag(S11)); %intrapolate points on a 2D grid

% Plot frequency-field colour map
[xq,yq] = meshgrid(linspace(field_l,field_h,501),linspace(freq_l,freq_h,310)); %set the X and Y range
zq = FdB(xq,yq);

% Plot the interpolated frequency response data in a field scan using color map
% hfig1 = setfig(11);
figure
cmap = pcolor(xq,yq,zq);
set(cmap, 'edgeColor','none')
shading interp;
colorbar
set(gca,'fontsize',plotopt.ftsz)
t(1) = xlabel('Field (T)');
t(2) = ylabel('Frequency (GHz)');
tt(1) = title([filename num2str(run)]);

% Plot the temperature profile against magnetic field
% hfig2 = setfig(12);
figure
plot(H(1:100:end),T1(1:100:end),'o-')
xlabel('DC Magnetic field')
ylabel('Temperature')
title('Magnetic field vs Temperature')
%
% yq = yq.*0 + 3.200; % Plot field scan cut at 3.200 GHz (H vs dB)
% hfig3 = setfig(12);
% zq = FdB(xq,yq);
% h = plot(xq(1,:),zq(1,:));
% xlim([0 9])
% set(gca,'fontsize',plotopt.ftsz)
% t(3) = xlabel('Field (T)');
% t(4) = ylabel('S11 (dB)');
% tt(2) = title([filename num2str(run) ', f = 3.200 GHz']);

% hfig4 = setfig(13); % Plot field scan cut at 3.200 GHz (H vs real(S))
% zq = FrS(xq,yq);
% h = plot(xq(1,:),zq(1,:));
% xlim([0 9])
% set(gca,'fontsize',plotopt.ftsz)
% t(5) = xlabel('Field (T)');
% t(6) = ylabel('Re{S11}');
% tt(3) = title([filename num2str(run) ', f = 3.200 GHz']);
%
% hfig5 = setfig(14); % Plot field scan cut at 3.200 GHz (H vs real(S))
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
% hfig6 = setfig(15);
% h = plot(yq(:,1),zq1(:,1),yq(:,1),zq2(:,1),yq(:,1),zq3(:,1),yq(:,1),zq4(:,1),yq(:,1),zq5(:,1));
% xlim([3.2 3.7])
% set(gca,'fontsize',plotopt.ftsz)
% t(9) = xlabel('Frequency (GHz)');
% t(10) = ylabel('S11 (dB)');
% tt(5) = title([filename num2str(run)]);
% t(11) = legend(h,'0','2','3.612','3.685','5');

% hfig7 = setfig(16); % Plot field-temperature
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


filepath = 'G:\My Drive\File sharing\PhD projects\LiHoF4 project\Data\Experiment\LiHoF4\Bulk(6x5x4.5mm)\05.05.2019';
% filepath = '/Users/yikaiyang/Google Drive/File sharing/PhD projects/LiHoF4 project/Data/Experiment/LiHoF4/SC107/19.05.2019/';
%The first line is for windows, the second line is for mac OS
filename = '2019_05_0006';

ii = 1;
runs{ii} =  1;   T(ii) = 0.085;
%% Plot data
%Set data range and parameters
fac = 1e9; %Swtich between "Hz" and "GHz"
freq_l = 3.400; %set frequency range, l: lower limit, h: higher limit
freq_h = 3.580;
field_l = 0.0;  %set field range, l: lower limit, h: higher limit
field_h = 9.0;
nop = 31; % Number of points per segment of the frequency scan according to ZVL network analyzer
step = 1e-4; %Scanning step in unit of GHz, default is 1 MHz according to ZVL network analyser
lineDiv = floor((freq_h-freq_l)/((nop - 1) * step)); %compute how many lines per complete frequency scan.

for n = 1:length(T)
    [xq,yq] = meshgrid(linspace(field_l,field_h,501),linspace(freq_l,freq_h,310));
    
    tout = readdata_v3(filepath,filename,runs{n});
    out(n) = mergeout_v1(tout);

    clear freq S11 dB N FdB FrS FiS FTT1 FTT2
    freq = out(n).data.ZVLfreq;
    S11 = out(n).data.ZVLreal + 1i*out(n).data.ZVLimag;
    H = out(n).data.DCField1;
    %     T1 = out(n).data.Temperature1;
    N = size(freq,2);   % Option_1 Number of elements of the measured frequency data
    %     N = length(yq(:,1));  % Option_2 Number of elements of the interpolated frequency data
    dB = mag2db(abs(S11));
    HH = repmat(H,1,N); %populate the magnetic field to match the dimension of S11 matrix
    %reshape the matrices into single column vectors
     
    %Find all the resonant peaks
    HH_temp = HH';
    dB_temp = dB';
    freq_temp = freq'/fac;
    
    trunc = mod(numel(dB),(lineDiv*nop)); %compute the number of complete frequency scans
   
    dB_temp = dB_temp(1:end-trunc); %discard the (last) incomplete frequency scan
    freq_temp = freq_temp(1:end-trunc);
    HH_temp = HH_temp(1:end-trunc);
    
    dB_temp = reshape(dB_temp,lineDiv*nop,[]);  %reshape the matrix so that each complete frequency scan occupy one column
    freq_temp = reshape(freq_temp,lineDiv*nop,[]);
    HH_temp = reshape(HH_temp,lineDiv*nop,[]);
    
    f0 = zeros(length(dB_temp),1);
    H0 = zeros(length(dB_temp),1);
    
%find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
    for i = 1:size(dB_temp,2) %Searching column minima (fixed field) is better than searching row minima (fixed frequency)
        [~,idx] = min( dB_temp(:,i) );
%         [value,idx] = min( dB_temp(:,i) );
%         H0(i) = HH(dB==value);    %use value matching to find the corresponding magnetic field from the original dataset (slow method)
        f0(i) = freq_temp(idx,i);
        H0(i) = HH_temp(idx,i);
    end

    %plot the resonant frequency versus DC magnetic field
    figure
    hfig(1) = plot(H0,f0,'o','MarkerSize',2);
    xlabel('Field (T)');
    ylabel('Resonant frequency (GHz)');
    title(num2str(T(n),'Resonant frequency from raw data at T = %3.3f K'));
    axis([field_l field_h freq_l freq_h]);
    clear dB_temp HH_temp freq_temp;
    
    freq = freq(:)/fac;
    dB = dB(:);
    HH = HH(:);
%     S11 = S11(:);
    mag = (1- out(n).data.ZVLreal.^2 -out(n).data.ZVLimag.^2)./((1- out(n).data.ZVLreal).^2 + out(n).data.ZVLimag.^2);
    mag = mag(:); %Unsure what these two lines are meant to do
    
    switch 2 %choose data interpolation method
        case 1
            % Option_1 Interpolate data along only the frequency axis.
            interp_dB = zeros(length(H),N);
            for i = 1:length(H)
                interp_dB(i,:) = interp1(freq(i,:),dB(i,:),yq(:,1));
            end
            interp_dB = interp_dB(:);
            yy = repmat(yq(:,1),1,length(H));
            yy = yy(:);

            T2 = out(n).data.Temperature2;    %If there are two sets of temperature data
            TT1 = repmat(T1,1,N);  %populate the temperature to match the dimension of S11 matrix
            TT1 = TT1(:);
            TT2 = repmat(T2,1,N); TT2 = TT2(:);
        case 2
            % Option_2 Interpolate the data along all (two) axis.
            FdB  = TriScatteredInterp(HH,freq,dB);
            Fmag  = scatteredInterpolant(HH,freq,mag);
            %     FrS  = TriScatteredInterp(HH,freq*fac,real(S11));
            %     FiS  = TriScatteredInterp(HH,freq*fac,imag(S11));
            %     FTT1 = TriScatteredInterp(HH,freq*fac,TT1);
            %     FTT2 = TriScatteredInterp(HH,freq*fac,TT2);
    end
    
    %plot the color map of S11 response
    switch 1 %choose colormap method
        case 1
            %option_1 plot surface interpolated data using pseudo-colormap
            figure
            hold on
            box on
            zq = FdB(xq,yq);
            cmap = pcolor(xq,yq,zq);
        case 2
    %     option_2.1 plot single-direction interpolated data using pseudo-colormap
            cmap = pcolor(H,yq(:,1),interp_dB');
        case 3
    %     option_2.2 plot single-direction interpolated data using scatter plot
            C = linspace(min(interp_dB),max(interp_dB),length(HH));
            cmap = scatter3(HH,yy,interp_dB,2,C,'o','filled','MarkerEdgeColor','none');
            colormap(hsv);
    end

    set(cmap, 'edgeColor','none')
    shading interp;
    colorbar
    set(gca,'fontsize',plotopt.ftsz)
    axis([field_l field_h freq_l freq_h]);
    tx(n+4) = xlabel('Field (T)');
    ty(n+4) = ylabel('Frequency (GHz)');
    tt(n+4) = title(num2str(T(n),'S11 response at T = %3.3f K'));
    
    %Fit each frequency scan to Lorentzian form to extract the quality factor
%     figure
    disp(num2str(T(n),'Fitting: T = %3.3f K'))
    % Extract frequency cuts
    Hx = field_l:0.1:field_h;
    
%   For noisy data, we need to remove duplicates of minima
    [H0,ia,~] = unique(H0,'last');
    f0 = f0(ia);
    ff0 = interp1(H0,f0,Hx,'spline','extrap'); %Create an interpolation of the resonant frequency along the DC magnetic field (Hx) axis
%     Hx = H0;
%     ff0 = f0; %Alernatively use the raw data for fitting
    ii = 1;
    for m = Hx
        xq = xq.*0 + m;
        zq = Fmag(xq,yq);
        zqf = medfilt1(zq(:,1),10);
        s = spec1d(yq(:,1),zqf,zqf.*0 + 0.05);
        p = [-1 ff0(ii) 0.001 2];
        %starting point for the (Lorentzian) fitting parameters(p1: scaling factor ,p2: resonant frequency; p3: FWHM; p4:noise floor?)
        
        fix = [1 1 1 1];
        [fQ, fbck]=fits(s,'lorz',p,fix);
        Q(ii) = fbck.pvals(2)/fbck.pvals(3); %Calculate the quality factor
        chi(ii) = 1/Q(ii);
        ii = ii + 1;
        disp(num2str(m,'Current magnetic field: %3.1f.'));
    end
        figure
        freqPlot = plot(Hx,ff0,'or','MarkerSize',2)
        xlabel('Field (T)');
        ylabel('Frequency (GHz)');
        title(num2str(T(n),'Resonant frequency from fitted data at T = %3.3f K'));
        axis([field_l field_h freq_l freq_h]);
        
    %     yq = yq.*0 + 3.49;
    %     zq = Fmag(xq,yq);
    %     h1(n) = plot(xq(1,:),zq); Error: the indices on the left side are not compatible with the size of the right side.
    %     hk = plot(xq(1,:),Q);
    %     plot(xq(1,:),zq);
%         col = setcolours(n/length(T),'jet');
        figure
        Qplot = plot(Hx,Q,'o-','MarkerSize',4);
        gca;
        xlabel('Field (T)');
        ylabel('Q factor');
        title('Quality factor');
%         set(hk(n),'color',col,'linewidth',plotopt.lnwd);
end
%
% set(ax,'fontsize',plotopt.ftsz,'xlim',[field_l field_h])
%
% set(ax(2),'ylim',[0.2 0.6])
% tl = legend(h1,num2str(T(:),'%3.1f K'),'location','se');
% set([tx ty tl tt],'fontsize',plotopt.ftsz,'interpreter','latex')
%     tt(2) = title([filename num2str(runs(n)) ', f = 3.437 GHz']);

% saveplots(hfig(2),'7-Temperature_field-S_-16dBm_Q')



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


