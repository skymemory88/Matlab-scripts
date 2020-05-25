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
% % use 'dealta' to move the frequency window and central freqeuncy by steps
% clearvars -except delta freq_l freq_h w0
% delta = -0.02;
% freq_l = freq_l+delta;
% freq_h = freq_h+delta;
% w0 = w0+delta; % Resonant frequency of the bare cavity
% freq_l = freq_l;
% freq_h = freq_h;
% w0 = w0; % Resonant frequency of the bare cavity

clearvars -except delta
freq_l = 3.69;
freq_h = 3.73;
freq_pts = 501;
freq = linspace(freq_l,freq_h,freq_pts);
filFctr = 1E-5;
w0 = 3.71; % Resonant frequency for bare cavity
f2E = 1/241.8; % Convert from GHz to meV
kB = 0.08617; % [meV/K]
kappa = 0.06; % Define coupling strength between the driving field and the cavity field
gamma = -2*kappa^2;
Gamma = 1.3*gamma; % Coupling strength between the cavity field and the spin system
% Gamma = -2.0*0.5*0.08^2; % fix Gamma for checkpoint

Option = 2;

path = 'G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\External source\Mean Field --Peter & Ivan\output';
cd(path);
temp = 0.150;
phi = 0; % field angle in ab-plane

if Option == 1
%% Option 1: Perturbative treatment of resonant frequency (working progress)
    load(num2str(temp*1000,'sim_%umK_trans.mat')); % Load the energy levels from mean field calculations
    w_up = 0.5.*(Ediff(:,:) + w0)+0.5.*sqrt((Ediff(:,:)-w0).^2 + 4*Gamma^2); %Upper perturbed resonant frequencies of the cavity by different transition lines
    w_low = 0.5.*(Ediff(:,:) + w0)-0.5.*sqrt((Ediff(:,:)-w0).^2 + 4*Gamma^2); % Lower perturbed resonant frequencies of the cavity by different transition lines

    wh = double.empty(size(w,1),0);
    wl = double.empty(size(w,1),0);
    ww = double.empty(size(w,1),0);
    for ii = 1:length(fields)
        [~,idx] = min(Ediff(:,ii)-w0);
        wh(ii) = 0.5.*(Ediff(idx,ii) + w0)+0.5.*sqrt((Ediff(idx,ii)-w0).^2 + 4*Gamma^2);
        wl(ii) = 0.5.*(Ediff(idx,ii) + w0)-0.5.*sqrt((Ediff(idx,ii)-w0).^2 + 4*Gamma^2);
        ww(ii) = (19*wh(ii)+wl(ii))/20; % Take a weighted average between the upper and lower levels
        ww(ii) = mean(w_up(:,ii));
    end
    ww = w0 - Gamma^2*(Ediff(idx,:)-w0)./((Ediff(idx,:)-w0).^2 + Gamma^2);
    ww = repmat(ww',[1,501]);
    figure
    plot(fields, w_up, '-k', fields, w_low, '-r'); % check point
    figure
    plot(fields, ws,'-k'); % check point
    figure
    plot(fields, ww,'-k'); % check point

    ws1 = -0.35.*B0+5; % make two fictitious dispersion relations for checkpoints
    ws2 = +0.75.*B0;
    plot(y,ws1,y,ws2);
    xlim([0 9])
    ylim([1 5])
elseif Option == 2
%% Option 2 Load existing susceptibilities and interpolate
    filename = ['Hscan_LiHoF4_',sprintf('%1$3.3fK_%2$uDeg.mat',temp,phi)];
    load(filename,'-mat','eee','vvv','fff'); % loads variables "Energy" and "EigenVector", and "Fields"
    E = squeeze(eee);
    V = vvv;
    field = vecnorm(fff);
    [w,B0] = meshgrid(freq,field);
    E2f = 241.8; % Convert Energy (meV) to frequency (GHz)
    Ediff = double.empty(0,size(E,1));
    for i=1:7
        Ediff(i,:) = (E(:,i+1)-E(:,i))*E2f;
    end
    filename = strcat('LiHoF4_x1z_x2z_',sprintf('%1$3.3fK_%2$uDeg.mat',temp,phi));
    load(filename,'-mat','fields','freq_total','x1z','x2z');
    x1 = interp2(fields,freq_total,x1z,B0,w);
    x2 = interp2(fields,freq_total,x2z,B0,w);
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
    [w,B0] = meshgrid(freq,field);
end

w0 = w0./(sqrt(1+filFctr.*(x1+1i*x2)));
Ediff = Ediff(1:7,:); % Select the levels to include

spins = double.empty(size(Ediff,2),0); % A container for the elements of spin term summation in the denominator of S11
spin_term = double.empty(size(Ediff,2),size(w,2),0); % A container for the spin term summation in the denominator of S11
popul = double.empty(size(Ediff,1),size(Ediff,2),0); % A container for the population factor of the transition levels
for ii = 1:size(Ediff,2)
    for kk = 1:size(Ediff,1)
       popul(kk,ii,1) = exp(-Ediff(kk,ii).*f2E./(kB*temp))./exp(-min(Ediff(kk,ii)).*f2E./(kB*temp)); % occupation factor
       spins(kk) = 2*Gamma^2.*sqrt(popul(kk,ii,1))./(1i.*(Ediff(kk,ii)-w0(ii,1)) + 3*abs(gamma)); % calculate the spin term for each transition level
    end
    spin_term(ii,:,1) = sum(spins(:)); % Sum over all levels
end
S11_2 = 1 + gamma./(1i.*(w0-w) -gamma + spin_term); % Calculate the S11 response
close all

figure
map = pcolor(B0,w,log(abs(S11_2))); % color plot of the S11 response
map.EdgeColor = 'none';
colorbar
xlim([0 17])
% ylim([3.69 3.73])
xlabel('Magnetic field (T)');
ylabel('Frequency (GHz)');
title(num2str(temp*1000,'Simulated data of S11 at %umK'));
caxis([-7 0])
% hold on
% levels = plot(fields, Ediff, '-k'); % plot the transition levels on top of the S11 color map
clearvars -except delta