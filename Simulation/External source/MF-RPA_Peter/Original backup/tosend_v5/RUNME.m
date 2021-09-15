
meV2GHz = 241.7681; % meV -> GHz

sample.omega=linspace(-0.05,0.25,201);     % meV E-cut
sample.qvec = 0;
sample.h = 0;
sample.k = 0;
sample.l = sample.qvec;

sample.hvec = linspace(0.5,0.6,3); 
sample.fields=[
    zeros(size(sample.hvec))
    zeros(size(sample.hvec))
    sample.hvec
    ];

sample.pop = [0.5 0.5];
calc1 = sim_RPA_4D_v5(sample);

figure(654654)
clf
sanePColor(sample.hvec,sample.omega,log(squeeze(imag(calc1.chitot(:,:,3,3,:,:)))'));

calc1.SS = squeeze(imag(calc1.chitot(:,:,3,3,:,:)));
calc1 = extract_peaks(calc1,'hvec','omega','SS',1,true,10);


%% Look how chi'' changes at H = Hc when Er(Hyp) = 0, 1 or 0.23
sample.ErHyp = 0;
calc2 = sim_RPA_4D_v5(sample);

figure(654655)
clf
sanePColor(sample.hvec,sample.omega,log(squeeze(imag(calc2.chitot(:,:,3,3,:,:)))'));

calc2.SS = squeeze(imag(calc2.chitot(:,:,3,3,:,:)))';
calc2 = extract_peaks(calc2,'hvec','omega','SS',2,true,10);

figure(212121)
clf
hold on

h1 = plot(calc1.extpeaks.xpar,calc1.extpeaks.sloc,'-o');
h2 = plot(calc2.extpeaks.xpar,calc2.extpeaks.sloc,'-s');

plot(sample.omega*meV2GHz,squeeze(imag(calc1.chitot(:,:,3,3,:,:)))')

%% Examine the excitations at the elastic line

clear sample
sample.omega=linspace(-0.02,0.02,201);     % meV E-cut
sample.qvec = linspace(0.4,1.6,31);
sample.h = sample.qvec;
sample.k = 0;
sample.l = 0;
sample.epsilon = 0.0001;

sample.hvec = 0.53; 
sample.fields=[
    zeros(size(sample.hvec))
    zeros(size(sample.hvec))
    sample.hvec
    ];

% sample.A_Er = 0.0005;   % hyperfine coupling energy in meV for Erbium (Phys. Rev. B 2, 2298 - 2301 (1970))
% sample.ErHyp = 0.23;
% sample.exEr = 0.0000;
sample.pop = [0.5 0.5];
calcE = sim_RPA_4D_v5(sample);

figure
sanePColor(sample.qvec,sample.omega,squeeze(calcE.Stot)')
caxis([0 2000])

