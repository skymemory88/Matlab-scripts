clear
meV2GHz = 241.7681; % meV to GHz

% Temperature (or range of temperatures)
sample.temp = 0.1; % [K]

sample.Erbium = 1.0;                  % Erbium proportion 
sample.ErHyp = 0.23;                  % Erbium hyperfine isotop proportion
% sample.h4 = 0.11e-3;                  % ih4 anisotropy (meV)
sample.h4 = 0;
sample.Holmium = 1.0-sample.Erbium;   % holmium proportion 
sample.HoHyp = 1;                     % Holmium hyperfine isotop proportion

sample.alpha = 0; % sha pe is a sphere (1) or needle (0, default)

freq = linspace(0,50,201);
sample.omega = freq/meV2GHz;
% sample.omega = linspace(-0.1,0.8,201); % meV E-cut
% sample.omega = linspace(-0.01,0.16,201);

% sample.epsilon = 1.68e-4; % Spin linewidth [meV] for Ho
sample.epsilon = 0.005;  % Spin linewidth [meV] for Er

% sample.qvec = 1:0.02:2;
sample.qvec = 0;

sample.h = sample.qvec;
sample.k = 0;
sample.l = 0;

sample.hvec = linspace(0,0.6,151); 
% sample.hvec = 4.2; % for 300 mK
sample.fields=[sample.hvec
               zeros(size(sample.hvec))
               zeros(size(sample.hvec))];

sample.pop = [0.5 0.5]; % Domain population
% sample.pop = [1 0];
calc1 = sim_RPA_4D_v5(sample);

if length(sample.qvec) > length(sample.hvec)
    for ii = 1:length(sample.hvec)
        figure
        sanePColor(sample.qvec,sample.omega,log(squeeze(real(calc1.chitot(:,ii,3,3,:,:)))'));
        set(gca, 'xdir', 'reverse' )
        ylabel('Energy (meV)')
        xlabel('[h 0 0]')
        legend(['Re(\chi),', num2str(sample.hvec(ii),'B = %.2f T')])
        figure
        sanePColor(sample.qvec,sample.omega,log(squeeze(imag(calc1.chitot(:,ii,3,3,:,:)))'));
        set(gca, 'xdir', 'reverse' )
        ylabel('Energy (meV)')
        xlabel('[h 0 0]')
        legend(['Im(\chi), ', num2str(sample.hvec(ii),'B = %.2f T')])
    end
else
    for ii = 1:length(sample.qvec)
        figure
        sanePColor(sample.hvec,sample.omega,log(squeeze(real(calc1.chitot(:,:,3,3,ii,:)))'));
        ylabel('Energy (meV)')
        xlabel('Magnetic field (T)')
        legend(['Re(\chi), ',sprintf('q = [%1$u %2$u %3$u]',sample.h,sample.k,sample.l)])
        figure
        sanePColor(sample.hvec,sample.omega,log(squeeze(imag(calc1.chitot(:,:,3,3,ii,:)))'));
        ylabel('Energy (meV)')
        xlabel('Magnetic field (T)')
        legend(['Im(\chi), ',sprintf('q = [%1$u %2$u %3$u]',sample.h,sample.k,sample.l)])
    end
end

% calc1.SS = squeeze(real(calc1.chitot(:,:,3,3,:,:)));
% calc1.SS = squeeze(imag(calc1.chitot(:,:,3,3,:,:)));

% calc1 = extract_peaks(calc1,'qvec','omega','SS',1,true,1,1,true);
% calc1 = extract_peaks(calc1,'hvec','omega','SS',1,true,1,1,true);

% calc1 = extract_peaks(calc1,'qvec','omega','SS',1,true,10);
% calc1 = extract_peaks(calc1,'hvec','omega','SS',1,true,10);

%% Look how chi'' changes at H = Hc when Er(Hyp) = 0, 1 or 0.23
% sample.ErHyp = 0;
% calc2 = sim_RPA_4D_v5(sample);
% 
% figure(654655)
% clf
% sanePColor(sample.qvec,sample.omega,log(squeeze(imag(calc2.chitot(:,:,3,3,:,:)))'));
% % sanePColor(sample.hvec,sample.omega,log(squeeze(imag(calc2.chitot(:,:,3,3,:,:)))'));
% 
% % calc2.SS = squeeze(imag(calc2.chitot(:,:,3,3,:,:)))';
% % calc2 = extract_peaks(calc2,'qvec','omega','SS',2,true,10);
% % calc2 = extract_peaks(calc2,'hvec','omega','SS',2,true,10);
% % figure(212121)
% % clf
% % hold on
% % h1 = plot(calc1.extpeaks.xpar,calc1.extpeaks.sloc,'-o');
% % h2 = plot(calc2.extpeaks.xpar,calc2.extpeaks.sloc,'-s');
% 
% plot(sample.omega*meV2GHz,squeeze(imag(calc1.chitot(:,:,3,3,:,:)))')
% 
% %% Examine the excitations at the elastic line
% 
% clear sample
% % sample.omega=linspace(-0.02,0.02,201);     % meV E-cut
% sample.omega = linspace(-0.01,0.16,201);
% 
% % sample.qvec = linspace(0.4,1.6,31);
% sample.qvec = 0:-0.02:-1;
% sample.h = sample.qvec;
% sample.k = 0;
% sample.l = 0;
% 
% % sample.epsilon = 0.0001; % Original code
% sample.epsilon = 0.005;
% 
% sample.hvec = 0.36; 
% sample.fields=[
%     zeros(size(sample.hvec))
%     zeros(size(sample.hvec))
%     sample.hvec
%     ];
% 
% % sample.A_Er = 0.0005;   % hyperfine coupling energy in meV for Erbium (Phys. Rev. B 2, 2298 - 2301 (1970))
% % sample.ErHyp = 0.23;
% % sample.exEr = 0.0000;
% sample.pop = [0.5 0.5];
% % sample.pop = [1 0];
% calcE = sim_RPA_4D_v5(sample);
% 
% figure
% sanePColor(sample.qvec,sample.omega,squeeze(calcE.Stot)')
% % set(gca, 'xdir', 'reverse' )
% % sanePColor(sample.hvec,sample.omega,squeeze(calcE.Stot)')
% caxis([0 2000])
