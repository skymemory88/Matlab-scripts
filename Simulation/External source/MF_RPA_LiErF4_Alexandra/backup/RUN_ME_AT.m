% Function to simulate the inelastic neutron scattering measurements on LiErF4 using MF-RPA
% Requires spec1d to run
% RUN_ME tested on Matlab 2020a on April 1 by Alexandra Turrini
% I've solved the problem of calling the Hamiltonian in remf, chi0, and
% chiq by making them into functions accessible in all, but there are many
% solutions. It works for single samples, not combined RE ions like the
% earlier versions

%% Adds Spectra and My Current Folder to the Database
% if contains(pwd,'Alexandra')
%     good_path = 'C:\Users\Alexandra\Dropbox\';
% elseif contains(pwd,'alext')
%     good_path = 'C:\Users\alext\Dropbox\';
% else
%     good_path = 'C:\Users\turrini_a\Dropbox\';
% end
% addpath(genpath([good_path,'Spectra'])) %% Call your spec1d here

addpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\spec1d--Henrik')

% Function to simulate the LET field experiment using RPA
plotopt.col  = [0.2 0.2 0.7];
plotopt.col2 = abs(0.8-plotopt.col);
plotopt.lnwd = 2;
plotopt.ftsz = 12;
plotopt.mksz = 5;

global strategies; % global convergence strategies switches
strategies.powerlaw=true; % Use a starting guess with the assumption the variables follow a power law.
strategies.accelerator=0.0; % accelerate. Beware: it often doesn't accelerate actually. Can speed up if the system has to break a symetry
strategies.damping=0.1; % damping factor. May reduce exponential fit efficiency
strategies.expfit=true; % turn on fitting to an exponential.
strategies.expfit_period=50; % period between exponential fit.
strategies.expfit_deltaN=20; % The exponential fit points distance. Three points used: N, N-deltaN, N-2deltaN.

%Define Compounds
run LiReF4_Compound_Properties.m

% Define Domains
momente_er_d1=[1 0 0;-1 0 0;-1 0 0;1 0 0];
momente_er_d2=[0 1 0;0 1 0;0 -1 0;0 -1 0];
momente_ho=repmat([0 0 1]*2.6,[4,1]);

%Define Dipole Range;
dip_range = 100;

% Define Demagnetization
withdemagn=1;
alpha=0; % shape is a sphere (1), needle (0).

% Define temperature
t = 0.2;

% Define field
h = 0.53;
% h = [0:0.1:0.5 0.53 0.6:0.1:1];
hvec = [zeros(size(h))',zeros(size(h))',h'];

if length(t)>1
    runfield = false;
    hist_d1 = zeros(length(t),4,3);
    hist_d2 = zeros(length(t),4,3);
else
    runfield = true;
    hist_d1 = zeros(length(h),4,3);
    hist_d2 = zeros(length(h),4,3);

end

%Define S(q,omega)
omega=[-0.1:0.001:0.252];%f2E = 1/241.8;  % GHz to meV
epsilon=0.005; % electronuclear spin line width
chi0_1 = zeros(3,3,size(LiErF4.tau,1),length(omega),length(h));
chi0_2 = zeros(3,3,size(LiErF4.tau,1),length(omega),length(h));

qh = -1.25:0.1:0.15;
q = qh'*[1 0 0]+ones(size(qh'))*[0 0 0];

chi_1 = zeros(3,3,length(omega),size(q,1),length(h));
chi_2 = zeros(3,3,length(omega),size(q,1),length(h));
    
%% Fixed-temperature field scan
for nh=1:size(hvec,1)
       
    %% Domain 1 ===========================================================
    
    [momente_er_1,momente_er_1_hf,eigenvals,wavefunctions] = remf(LiErF4,hvec(nh,:),t,momente_er_d1,dip_range,withdemagn,alpha);
    hist_d1(nh,:,:) = momente_er_1;

    [chi0_1(:,:,:,:,nh),eigenvals_chi0,wavefunctions_chi0] = chi0_w(LiErF4,hvec(nh,:),t,momente_er_1,momente_er_1_hf,0,dip_range,withdemagn,alpha,omega,epsilon);
    srpa_1 = gen_scattering_crosssec_chi0(LiErF4,squeeze(chi0_1(:,:,:,:,nh)),q,omega,t,dip_range,withdemagn,alpha);
     
   %% Domain 2 ===========================================================
    
    [momente_er_2,momente_er_2_hf] = remf(LiErF4,hvec(nh,:),t,momente_er_d2,dip_range,withdemagn,alpha);
        hist_d2(nh,:,:) = momente_er_2; % records moment strength and vector at each field
    
    chi0_2(:,:,:,:,nh) = chi0_w(LiErF4,hvec(nh,:),t,momente_er_2,momente_er_2_hf,0,dip_range,withdemagn,alpha,omega,epsilon);
    srpa_2 = gen_scattering_crosssec_chi0(LiErF4,squeeze(chi0_2(:,:,:,:,nh)),q,omega,t,dip_range,withdemagn,alpha);
    
    
        % Total contribution to the scattering cross-section from both domains
    for nq = 1:size(q,1)
       srpa_total(nh,nq) = srpa_1(nq) + srpa_2(nq);
%        srpa_total(nh,nq) = srpa_1(nq);
%        srpa_total(nh,nq) = combine(srpa_1(nq),srpa_2(nq));
    end

end

[x_HcT,y_HcT,e_HcT] = extract(srpa_total(1,:));

figure;
p = pcolor(qh,omega,y_HcT);
 set(p,'EdgeColor','none')
cb = colorbar;
cb.TickLabelInterpreter = 'latex';
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\chi_q^{\prime \prime}$ at 0.53 T, 200 mK';
set(gca,'TickLabelInterpreter','latex','FontSize',14)
% legend([p2;p3;p4;p5],'E$_{remf}$ D1','E$_{\chi^0}$ D1','E$_{remf}$ D2','E$_{\chi^0}$ D2','Interpreter','latex','Box','off')
xlabel('$(h, 0, 0)$ (r.l.u.)','Interpreter','latex')
ylabel('Energy (meV)','Interpreter','latex')
% caxis([0 50])

% % 