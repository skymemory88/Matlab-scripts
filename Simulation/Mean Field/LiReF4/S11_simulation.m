%% Strong/weak coupling simulation
% B = linspace(0,9,100);
% Br = 3.72;
% m = 7/2;
% hbar = 1;
% Delt = -m*(B-Br)/hbar;
% wc = 3.674;
% gc = 0.055;
% gamma = 0.1;
% w = wc - gc^2.*Delt./(Delt.^2+gamma^2);
% hp1 = plot(B,w,'-k','LineWidth',3);
% hold on
% wp = wc + Delt./2 + sqrt(Delt.^2+4*gc^2)/2;
% wm = wc + Delt./2 - sqrt(Delt.^2+4*gc^2)/2;
% hp2 = plot(B,wm,'-r',B,wp,'-r','LineWidth',3);
%% Bare cavity S11 Simulation
% clearvars
% kappa = 0.2;
% gamma = -2.0*kappa^2;
% w0 = 3.5;
% w = linspace(3,4,101);
% S11 = 1+gamma./(1i.*(w-w0)-gamma/2);
% 
% S11_real = real(S11);
% S11_img = imag(S11);
% % S11_2 = S11.^2;
% 
% figure
% plot(w,S11_real,'-o');
% title('Real part of S11 respons')
% xlabel('Frequency (GHz)');
% ylabel('S11 response');
% 
% figure
% plot(w,S11_img,'-o');
% title('Imaginary part of S11 respons')
% xlabel('Frequency (GHz)');
% ylabel('S11 response');
%% 2D color plot of S11
clearvars -except delta
% freq_l = 3.685; % for use in mathematica
% freq_h = 3.735;
% w0 = 3.71;
% filFctr = 5E-3;

temp = 0.150; % temperature(s)
theta = 0; % Angle (in degrees) in the transverse field plane

%frequency parameters
freq_l = 3.49;
freq_h = 3.51;
freq_pts = 501; % number of points along frequency axis
freq = linspace(freq_l,freq_h,freq_pts); % frequency range

% field parameters
field_l = 0;
field_h = 9;
field_pts = 501; % number of points along field axis
field = linspace(field_l,field_h,field_pts); % sampling points along field axis

w0 = 3.50; % Resonant frequency for bare cavity
filFctr = 1E-4; % Calculated from COMSOL


f2E = 1/241.8; % Convert from GHz to meV
kB = 0.08617; % [meV/K]
kappa = 0.09; % Define coupling strength between the driving field and the cavity field
gamma = -2*kappa^2;
Gamma = 0.4*gamma; % Coupling strength between the cavity field and the spin system
% Gamma = -0.1*0.05^2; % fix Gamma for checkpoint

Option = 1;
Options.Elevel = false;
Options.Ediff = false;
Options.savedata = false;

path = genpath('G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\Mean Field\LiReF4\output');
addpath(path);
phi = 0; % field angle in ab-plane

if Option == 1
%% Option 1: Perturbative treatment of resonant frequency
beta = 1/(temp*kB);
g = 0.001; % Coupling strength measured against the ground state energy
% g = 0.01*w0;

%         cd('G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\External source\Mean Field --Peter & Ivan\output\A=1.0_angles_Yikai');
addpath(genpath('G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\Mean Field\LiReF4\output'));
%         cd('/Volumes/GoogleDrive/My Drive/File sharing/Programming scripts/Matlab/Simulation/External source/Mean Field --Peter & Ivan/output/A=1.0_angles_Peter')
lname=['Hscan_LiHoF4_', sprintf('%1$3.3fK_%2$uDeg',temp,theta),'.mat'];

load(lname,'-mat','eee','fff');
fields = vecnorm(fff);
E2f = 241.8; % Convert Energy to frequency
E(:,:)=squeeze(eee)*E2f;
%E(:,:) = eee(:,1,:)*E2f;
Ediff = double.empty(7,size(E,1),0);
bzF = double.empty(size(Ediff,1),size(fields,2),0);
for ii = 1:7
    Ediff(ii,:,1)=E(:,ii+1)-E(:,ii);
    bzF(ii,:,1) = exp(-Ediff(ii,:)*beta); % Boltzmann factor for each transition line
end

%Optional: Plot the lowest eight energy levels
if Options.Elevel == true
    color = ["black","red","blue","magenta","green","yellow","cyan"];
    theta = [0]; % Angle (in degrees) in the transverse field plane
    marker = [":","-.","--","-"];
    figure
    hold on
    plot(fields,E(:,1:8),'Color',[35 107 142]/255,'linewidth',2);
    set(gca,'fontsize',15);
    % text(8,(j-1)*D+0.02,'0.1');
    %         set(fig1,'position',[100 100 600 300])
    xlabel('Field (T)','FontSize',15)
    ylabel('Energy (GHz)','FontSize',15)
    grid off
    box on;
    xlim([0 9]);
    %         ylim([2 5]);
    tit1='Energy levels';
    title(tit1,'FontSize',15)
    legend(num2str(temp*1000,'T = %u mK,  A = A_{th}'))
end

% Optional: Energy difference bewteen neighbour levels
if Options.Ediff == true
    figure 
    hold on
    figs(:, iter2) = plot(fields, Ediff, 'Marker', 'none', 'LineStyle',marker(1), 'Color',color(iter2),'LineWidth',2);
    hold on
    lg(iter, iter2) = [num2str(temp(iter)*1000,"T = %u mK, A = A_{th}"),num2str(theta(iter2),", theta = %u")];
    set(gca,'fontsize',15,'Xtick',0:1:max(fields));
    xlabel('Field (T)','FontSize',15);
    ylabel('Energy difference (GHz)','FontSize',15);
    grid off
    box on;
    xlim([min(fields) max(fields)]);
    ylim([0.9*min(Ediff,[],'All') 1.1*max(Ediff,[],'All')]);
    title('Difference between energy levels','FontSize',15);
end
% Calculate perturbed discpersion relations

w = double.empty(size(Ediff,1)+1,size(field,2),0); % eigen-frequencys
% boltFac = double.empty(size(Ediff,1),size(field,2),0); % Boltzmann factor 
temp_Ediff = double.empty(size(Ediff,1),size(field,2),0);
w(1,:,1) = w0;
for ii = 2:size(Ediff,1)+1
    w(ii,:) = interp1(fields,Ediff(ii-1,:),field);
    temp_Ediff(ii-1,:,1) = interp1(fields,Ediff(ii-1,:),field);
%     boltFac(ii-1,:,1) = interp1(fields,bzF(ii-1,:),field);
end
Ediff = temp_Ediff;
clearvars temp_Ediff

% Construct interaction Hamiltonian along magnetic field axis
ww = double.empty(0,size(Ediff,1)+1,size(field,2));
for jj = 1:size(field,2)
    % off diagonal elements (coupling strength)
    gs = g./abs(w0-w(2:end,jj));
%    gs = g./abs(w0-w(2:end,jj)).*(boltFac(:,jj)/max(boltFac(:,jj)));
%    gs = g/abs(w0-w(2,jj))*ones(1,size(Ediff,1));
%    gs = g/abs(w0-w(2,jj))*ones(1,size(Ediff,1))*(boltFac(:,jj)/max(boltFac(:,jj)));

    % check point: use to turn off interactions bewteen cavity and certain spin levels
%    gs(1:2) = 0; 
%    gs(1:end-1) = 0; 

    Ham = diag(w(:,jj));
    Ham(1,2:end) = gs;
    Ham(2:end,1) = gs';
    ww(1,:,jj) = eig(Ham);
end
% pick out the frequencies along the field axis that is closes to bare frequency of the cavity
wn = abs(squeeze(ww)-w0);
w0 = ww(wn == min(wn,[],1));
clearvars wn ww

figure
plot(field,w0,'--k','linewidth',1);
hold on
plot(field,Ediff,'-r');
xlim([field_l field_h])
ylim([freq_l freq_h])

w0 = repmat(w0,1,size(freq,2)); % expand the resonant frequency to a matrix for calculations in the next step
[freq,field] = meshgrid(freq,field);
elseif Option == 2
%% Option 2 Load existing susceptibilities and interpolate
    filename = ['Hscan_LiHoF4_',sprintf('%1$3.3fK_%2$uDeg.mat',temp,phi)];
    load(filename,'-mat','eee','vvv','fff'); % loads variables "Energy" and "EigenVector", and "Fields"
    E = squeeze(eee);
    fields = vecnorm(fff);
    V = vvv;
    E2f = 241.8; % Convert Energy (meV) to frequency (GHz)
    Ediff = double.empty(0,size(E,1));
    for i=1:7
        Ediff(i,:) = (E(:,i+1)-E(:,i))*E2f;
    end
    
    % interpolate the dispersion relation along field axis
    temp_Ediff = double.empty(size(Ediff,1),size(field,2),0);
    for ii = 1:size(Ediff,1)
        temp_Ediff(ii,:,1) = interp1(fields,Ediff(ii,:),field);
    end
    Ediff = temp_Ediff;
    clearvars temp_Ediff fff vvv eee
    
    filename = strcat('LiHoF4_x1z_x2z_',sprintf('%1$3.3fK_%2$uDeg.mat',temp,phi));
    load(filename,'-mat','fields','freq_total','x1z','x2z');
    
    [freq,field] = meshgrid(freq,field);
    x1 = interp2(fields,freq_total,x1z,field,freq);
    x2 = interp2(fields,freq_total,x2z,field,freq);
    w0 = w0./(sqrt(1+filFctr.*(x1+1i*x2)));
% w0 = w0./(sqrt(1+filFctr.*x1)); % Use only real part of the susceptibility
% w0 = w0./(sqrt(1+filFctr.*1i*x2)); % Use only imaginary part of the susceptibility
elseif Option == 3
%% Option 3: Calculate susceptabilities for resonant frequency shift (takes long time)
    filename = ['LHF_',num2str(temp,'%.3f.mat')];
    load(filename,'-mat','eee','vvv','fff'); % loads variables "Energy" and "EigenVector", and "Fields"
    if ~exist ('x1','var') && ~exist('x2','var')
        gama = 0.00005; % define lifetime (meV) for hyperfine levels
        [fields, Ediff, x1, x2] = linear_response(eee,fff,ttt,vvv,gama); % Calculate susceptibilities
    end
    x1 = x1';
    x2 = x2';
    [freq,field] = meshgrid(freq,field);
    w0 = w0./(sqrt(1+filFctr.*(x1+1i*x2)));
% w0 = w0./(sqrt(1+filFctr.*x1)); % Use only real part of the susceptibility
% w0 = w0./(sqrt(1+filFctr.*1i*x2)); % Use only imaginary part of the susceptibility
end

spins = double.empty(size(Ediff,2),0); % A container for the elements of spin term summation in the denominator of S11
spin_term = double.empty(size(Ediff,2),size(freq,2),0); % A container for the spin term summation in the denominator of S11
popul = double.empty(size(Ediff,1),size(Ediff,2),0); % A container for the population factor of the transition levels
for ii = 1:size(field,1)
    for kk = 1:size(Ediff,1)
       popul(kk,ii,1) = exp(-Ediff(kk,ii).*f2E./(kB*temp))./exp(-min(Ediff(kk,ii)).*f2E./(kB*temp)); % occupation factor
       spins(kk) = 2*Gamma^2.*sqrt(popul(kk,ii,1))./(1i.*(Ediff(kk,ii)-w0(ii,1)) + 3*abs(gamma)); % calculate the spin term for each transition level
%        spins(kk) = 2*Gamma^2.*sqrt(popul(kk,ii,1))./(1i.*(Ediff(kk,ii)-w0(kk,ii)) + 3*abs(gamma)); % calculate the spin term for each transition level
    end
    spin_term(ii,:,1) = sum(spins(:)); % Sum over all levels
end
S11_2 = 1 + gamma./(1i.*(w0-freq) -gamma + spin_term); % Calculate the S11 response
close all

figure
map = pcolor(field,freq,log(abs(S11_2))); % color plot of the S11 response
map.EdgeColor = 'none';
colorbar
% xlim([3 4.5])
% ylim([3.69 3.73])
xlabel('Magnetic field (T)');
ylabel('Frequency (GHz)');
title(num2str(temp*1000,'Simulated data of S11 at %umK'));
caxis([-7 0])
% hold on
% levels = plot(fields, Ediff, '-k'); % plot the transition levels on top of the S11 color map
clearvars -except delta