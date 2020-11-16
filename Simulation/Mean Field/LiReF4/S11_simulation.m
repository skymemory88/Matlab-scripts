%% Strong/weak coupling simulation
% B = linspace(0,9,100);
% Br = 3.72;
% m = 7/2;
% hbar = 1;
% Delt = -m*(B-Br)/hbar;
% wc = 3.674;
% gc = 0.055;
% gamma_e = 0.1;
% w = wc - gc^2.*Delt./(Delt.^2+gamma_e^2);
% hp1 = plot(B,w,'-k','LineWidth',3);
% hold on
% wp = wc + Delt./2 + sqrt(Delt.^2+4*gc^2)/2;
% wm = wc + Delt./2 - sqrt(Delt.^2+4*gc^2)/2;
% hp2 = plot(B,wm,'-r',B,wp,'-r','LineWidth',3);
%% Bare cavity S11 Simulation
% clearvars
% kappa = 0.2;
% gamma_e = -2.0*kappa^2;
% w0 = 3.5;
% w = linspace(3,4,101);
% S11 = 1+gamma_e./(1i.*(w-w0)-gamma_e/2);
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

temp = 1.2; % temperature(s)
theta = 0.0; % angular deviation (in degrees) of the field angle from c-axis direction
phi = 0.0; % angle (in degrees) in ab-plane rotation
% alpha = linspace(0,pi,18); % Phase angle between coherent and dissipative couplings
alpha = 0;

% field parameters
field_l = 0;
field_h = 3;
field_pts = 601; % number of points along field axis
field = linspace(field_l,field_h,field_pts); % sampling points along field axis

%frequency parameters
w0 = 3.51; % Resonant frequency for bare cavity
freq_l = 3.435;
freq_h = 3.545;
freq_pts = 601; % number of points along frequency axis
freq = linspace(freq_l,freq_h,freq_pts); % frequency range

f2E = 1/241.8; % Convert from GHz to meV
kB = 0.08617; % [meV/K]
filFctr = 0.01; % Calculated from COMSOL
gama = 7e-4; % Spin state lifetime (meV)
gamma_i = -0.01; % internal dissipation rate
% gamma_i = -0.01; % internal dissipation rate
gamma_e = 1.25*gamma_i; % external dissipation rate
Gamma = 1*gamma_i; % Coupling strength between the cavity field and the spin system
% Gamma = 0.001*0.004; % fix Gamma for checkpoint

Option = 2; % Analysis options
Options.Ediff = false; % Plot energy levels on top
Options.x1x2 = false;
Options.savedata = false; % Save results to a file

if Option == 1
%% Option 1: Perturbative treatment of resonant frequency
beta = 1/(temp*kB);
g = 0.001; % Coupling strength measured against the ground state energy
% g = 0.01*w0;

%         cd('G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\External source\Mean Field --Peter & Ivan\output\A=1.0_angles_Yikai');
location = 'G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\Mean Field\LiReF4\output';
%         cd('/Volumes/GoogleDrive/My Drive/File sharing/Programming scripts/Matlab/Simulation/External source/Mean Field --Peter & Ivan/output/A=1.0_angles_Peter')
lname=['Hscan_LiHoF4_', sprintf('%1$3.3fK_%2$uDeg',temp,theta),'.mat'];
file = fullfile(location,lname);
load(file,'-mat','eee','fff');
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
    location = 'G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\Mean Field\LiReF4\output\without Hz_I';
    filename = ['Hscan_LiHoF4_',sprintf('%1$3.3fK_%2$.1fDeg_%3$.1fDeg.mat',temp,theta,phi)];
    file = fullfile(location,filename);
    load(file,'-mat','eee','vvv','fff'); % loads variables "Energy" and "EigenVector", and "Fields"
    E = squeeze(eee);
    fields = vecnorm(fff);
    V = vvv;
    E2f = 241.8; % Convert Energy (meV) to frequency (GHz)
    %% Calculate the expecation value of the spin moment
    %     gLande_Ho=1.25;
    %     ELEf = gLande_Ho * 0.05788;     % Lande factor * Bohr magneton (meV T^-1)
    %     NUCf = 4.173 * 3.15245e-5;   % Nuclear Lande factor, mu/mu_N = 4.173
    %     % NUCf = 4.732 * 3.1519e-5;   % Original code
    %     J=8; % Electronic moment for Ho3+
    %     I=3.5; % Nuclear moment for Ho3+
    %     %Initiate J operators
    %     Jz=diag(J:-1:-J); % Jz = -J, -J+1,...,J-1,J
    %     Jzh=kron(Jz,eye(2*I+1)); % Expand Jz space to include nuclear degree of freedom
    % %     Iz=diag(I:-1:-I); %Iz = -I, -I+1,...,I-1,I
    % %     Izh=kron(eye(2*J+1),Iz); % Expand Hilbert space
    %     Jz_exp = double.empty(0,length(fields(1,:)));
    % %     JIz_exp = double.empty(0,length(fields(1,:)));
    %     for kk = 1:length(fields(1,:)) % calculate susceptibility for all fields
    %         v = squeeze(squeeze(V(kk,:,:,:))); % Obtain the corresponding eigen vectors
    %         JzhT = Jzh * ELEf;
    % %         IzhT = Izh * NUCf;
    %         tz = v'  * JzhT * v;
    % %         ttz  = v'  * (JzhT+IzhT) * v;
    % %         Jz_exp(1,kk) = sqrt(sum(sum((tz) .* (tz.'))));
    %         Jz_exp(1,kk) = real(sum(diag(tz)));
    % %         JIz_exp(1,kk) = sqrt(sum(sum((ttz) .* (ttz.'))));
    % %         JIz_exp(1,kk) = real(sum(diag(ttz)));
    %     end
    %     Jz_exp = interp1(fields,Jz_exp,field);
    % %     JIz_exp = interp1(fields,JIz_exp,field);
    %% interpolate the dispersion relation along field axis
    if Options.Ediff == true
        Ediff = double.empty(0,size(E,1));
        for i=1:7 % upper limit = size(E,2)-1
            Ediff(i,:) = (E(:,i+1)-E(:,i))*E2f;
        end
        temp_Ediff = double.empty(size(Ediff,1),size(field,2),0);
        for ii = 1:size(Ediff,1)
            temp_Ediff(ii,:,1) = interp1(fields,Ediff(ii,:),field);
        end
        Ediff = temp_Ediff;
    end
    clearvars temp_Ediff fff vvv eee
    
    % Load the susceptibilities from MF-linear response calculations
    filename = strcat('LiHoF4_x1z_x2z_',sprintf('%1$3.3fK_%2$.1fDeg_%3$.1fDeg_%4$.2e.mat',temp,theta,phi,gama));
    file = fullfile(location,filename);
    load(file,'-mat','fields','freq_total','x1z','x2z');
    
    if Options.x1x2 == true
        figure;
        hp1 = pcolor(fields(1,:),freq_total,x1z);
        set(hp1, 'edgeColor','none')
        caxis([-23 2]);
        colorbar
        xlabel('Magnetic field (T)')
        ylabel('Frequency (GHz)')
        title({'Real part of $\chi$ in z direction'})
        
        figure
        hp2 = pcolor(fields(1,:),freq_total,log(x2z));
        set(hp2, 'edgeColor','none')
        colorbar
        xlabel('Magnetic field (T)')
        ylabel('Frequency (GHz)')
        title({'Imaginary part of $\chi$ (log scale)', 'in x direction'})
    end
    
    [freq,field] = meshgrid(freq,field);
    x1 = interp2(fields,freq_total,x1z,field,freq);
    x2 = interp2(fields,freq_total,x2z,field,freq);
    w0 = w0./(sqrt(1+filFctr.*(x1+1i*x2)));
%     w0 = w0./(sqrt(1+filFctr.*x1)); % Use only real part of the susceptibility
%     w0 = w0./(sqrt(1+filFctr.*1i*x2)); % Use only imaginary part of the susceptibility
%     w0 = w0./(sqrt(1+filFctr.*(0.0+0.0))); % Equivalent to empty cavity
elseif Option == 3
%% Option 3: Calculate susceptabilities for resonant frequency shift (takes long time)
    location = 'G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\Mean Field\LiReF4\output';
    filename = ['LHF_',num2str(temp,'%.3f.mat')];
    file = fullfile(location,filename);
    load(file,'-mat','eee','vvv','fff'); % loads variables "Energy" and "EigenVector", and "Fields"
    if ~exist ('x1','var') && ~exist('x2','var')
        life = 1/0.00005; % define lifetime (meV-1) for hyperfine levels
        [fields, Ediff, x1, x2] = linear_response(eee,fff,ttt,vvv,1/life); % Calculate susceptibilities
    end
    x1 = x1';
    x2 = x2';
    [freq,field] = meshgrid(freq,field);
    w0 = w0./(sqrt(1+filFctr.*(x1+1i*x2)));
%     w0 = w0./(sqrt(1+filFctr.*x1)); % Use only real part of the susceptibility
    % w0 = w0./(sqrt(1+filFctr.*1i*x2)); % Use only imaginary part of the susceptibility
end
%   % Calculate the spin terms in the denominator
%   life = 1/(200*abs(gamma_e)); % define spin level lifetime (meV-1)
%   % life = 0.0001;
%   % life = 0.00005;
%   spins = double.empty(size(Ediff,2),0); % A container for the elements of spin term summation in the denominator of S11
%   spin_term = double.empty(size(Ediff,2),size(freq,2),0); % A container for the spin term summation in the denominator of S11
%   popul = double.empty(size(Ediff,1),size(Ediff,2),0); % A container for the population factor of the transition levels
%   for ii = 1:size(field,1)
%       for kk = 1:size(Ediff,1)
%    %        popul(kk,ii,1) = exp(-Ediff(kk,ii).*f2E./(kB*temp))./exp(-min(Ediff(kk,ii)).*f2E./(kB*temp)); % occupation factor
%    %        spins(kk) = JIz_exp(1,ii)*2*Gamma^2.*sqrt(popul(kk,ii,1))./(1i.*(Ediff(kk,ii)-w0(ii,1)) + 1/life); % calculate the spin term for each transition level
%           spins(kk) = 2*Gamma^2.*Jz_exp(1,ii)./(1i.*(Ediff(kk,ii)-w0(ii,1)) + 1/life);
%    %        spins(kk) = 2*Gamma^2.*sqrt(popul(kk,ii,1))./(1i.*(Ediff(kk,ii)-w0(kk,ii)) + 3*abs(gamma_e)); % calculate the spin term for each transition level
%       end
%        spin_term(ii,:,1) = sum(spins(:)); % Sum over all levels
%   end
% S11_2 = 1 + 2*gamma_e./(1i.*(w0-freq) - gamma_i - gamma_e + spin_term); % Calculate the S11 response using spin turms
% S11_2 = 1 + 2*gamma_e./(1i.*(w0-freq) - gamma_i - gamma_e + Gamma.*x2); % Calculate the S11 response using susceptibilites
for ii = 1:length(alpha)
%     S11_2 = 1 + 2*gamma_e./(1i.*(w0-freq) - gamma_i - gamma_e + (1+exp(1i*alpha(ii)))*Gamma^2.*x2);
    S11_2 = 1 + 2*gamma_e./(1i.*(w0-freq) - gamma_i - gamma_e + Gamma^2.*x2);
%     close all
    figure
    map = pcolor(field,freq,log(abs(S11_2))); % color plot of the S11 response
    map.EdgeColor = 'none';
    colorbar
    xlim([field_l field_h])
    ylim([freq_l freq_h])
    xlabel('Magnetic field (T)');
    ylabel('Frequency (GHz)');
    title(num2str(temp*1000,'Simulated data of S11 at %umK'));
    caxis([-10 0])
    if Options.Ediff == true
        hold on
        yyaxis right
        levels = plot(field(:,1), Ediff(:,:), '--k'); % plot the transition levels on top of the S11 color map
    end
    xlim([field_l field_h])
    ylim([freq_l freq_h])
end
clearvars -except delta