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
clearvars
cd('G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\External source\Mean Field --Peter & Ivan\output');
temp = 0.150;
load(num2str(temp*1000,'sim_%umK_trans.mat')); % Load the energy levels from mean field calculations

Ediff = Ediff(1:7,:); % Select the levels to include
f2E = 1/241.8; % Convert from GHz to meV
kB = 0.08617; % [meV/K]
kappa = 0.05; % Define coupling strength between the driving field and the cavity field
gamma = -2.0*kappa^2;
Gamma = 2.2*gamma; % Coupling strength between the cavity field and the spin system
% Gamma = -2.0*0.5*0.08^2; % fix Gamma for checkpoint
w0 = 3.85; % Resonant frequency of the bare cavity
x = linspace(3.83,3.88,501); % calculation range of the frequency
[w,B0] = meshgrid(x,fields);

% ws1 = -0.35.*B0+5; % make two pseudo-dispersion relations for checkpoints
% ws2 = +0.75.*B0;
% plot(y,ws1,y,ws2);
% xlim([0 9])
% ylim([1 5])

spins = double.empty(size(Ediff,2),0); % A container for the elements of spin term summation in the denominator of S11
spin_term = double.empty(size(Ediff,2),size(w,2),0); % A container for the spin term summation in the denominator of S11
popul = double.empty(size(Ediff,1),size(Ediff,2),0); % A container for the population factor of the transition levels
for ii = 1:size(Ediff,2)
    for kk = 1:size(Ediff,1)
       popul(kk,ii,1) = exp(-Ediff(kk,ii).*f2E./(kB*temp)); 
       spins(kk) = 2*Gamma^2.*sqrt(popul(kk,ii,1))./(1i.*(Ediff(kk,ii)-w0)+ 4*abs(gamma)); % calculate the spin term for each transition level
    end
    spin_term(ii,:,1) = sum(spins(:)); % Sum over all levels
end
S11_2 = 1 - 2*kappa^2./(1i.*(w0-w) + 1.9*kappa^2 + spin_term); 

map = pcolor(B0,w,log(abs(S11_2))); % color plot of the S11 response
map.EdgeColor = 'none';
colorbar
xlim([0 9])
% ylim([3.69 3.73])
xlabel('Magnetic field (T)');
ylabel('Frequency (GHz)');
title(num2str(temp*1000,'Simulated data of S11 at %umK'));
caxis([-7 0])
% hold on
% levels = plot(fields, Ediff, '-k'); % plot the transition levels on top of the S11 color map
