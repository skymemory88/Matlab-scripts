function MF_linear_response
clearvars;
Options.RPA = false;
Options.plotting = false; % Decide whether or not to plot the data at the end
Options.saving = true;

% Temperatures = [0.08 0.1 0.12 0.15 0.2 0.24];
Temperatures = 1.7;
theta = 0.0;
phi = 0.0;
    for ii = 1:length(Temperatures)
        location = 'G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\Mean Field\LiReF4\output\without Hz_I';
        filename = strcat('Hscan_LiHoF4_', sprintf('%1$3.3fK_%2$.1fDeg_%3$.1fDeg', Temperatures(ii), theta, phi),'.mat');
        file = fullfile(location,filename);
        load(file,'-mat','eee','fff','ttt','vvv'); % loads variables "fields", "temp", "E" and "V" 
        % which are eigenstates and eigenvalues calculated in the mean-field model 
        % as a function of transverse field and temperature
        
        fprintf('Calculating for T = %.3f.\n', Temperatures(ii));
%         [fields, freq_total, rechi, imchi] = linear_response(eee,fff,ttt,vvv);
        [fields, freq_total, rechi, imchi, ~] = linear_response(eee,fff,ttt,vvv);
        if Options.RPA == true
            [fields, freq_total, rechi, imchi, ~] = RPA(fields, freq_total, rechi, imchi);
        end
        if Options.saving == true
            save_file(fields,freq_total,rechi,imchi,ttt,theta,phi,location);
        end
    end
%% Plot the susceptibilities
if Options.plotting == true
% Color plot of the imaginary part of the susceptibility of x component
    figure (1);
    clf
    hp1 = pcolor(fields(1,:),freq_total,squeeze(log(imchi.x)));
    set(hp1, 'edgeColor','none')
    caxis([-23 2]);
    colorbar
    xlabel('Magnetic field (T)')
    ylabel('Frequency (GHz)')
    title({'Imaginary part of Susceptibility (log scale)', 'in x direction'})
    
% Color plot of the imaginary part of the susceptibility of y component
    hfig1 = figure (2);
    clf
    hp1 = pcolor(fields(1,:),freq_total,squeeze(log(imchi.y)));
    set(hp1, 'edgeColor','none')
    caxis([-23 2]);
    colorbar
    xlabel('Magnetic field (T)')
    ylabel('Frequency (GHz)')
    title({'Imaginary part of Susceptibility (log scale)', 'in y direction'})
    
% Color plot the imaginary part of the susceptibilities of z component
    figure (3);
    clf
    hp2 = pcolor(fields(1,:),freq_total,squeeze(log(imchi.z)));
    set(hp2, 'edgeColor','none')
    caxis([-23 2]);
    colorbar
    xlabel('Magnetic field (T)')
    ylabel('Frequency (GHz)')
    title({'Imaginary part of Susceptibility (log scale)', 'in z direction'})

% Plot the real part of the susceptibility of the z component
    figure (4);
    clf
    set(hfig1,'position',[50 100 600 400])
    hp3 = pcolor(fields(1,:),freq_total,squeeze(rechi.z));
    set(hp3, 'edgeColor','none')
    caxis([-23 2]);
    colorbar
    xlabel('Magnetic field (T)')
    ylabel('Frequency (GHz)')
    title({'Real part of Susceptibility', 'in z direction'})
    
% % Plot the expectation value of <Jz+Iz>
%     figure (5);
%     clf
%     hp4 = plot(fields(1,:),squeeze(Jz_exp),'-o');
%     xlabel('Magnetic field (T)')
%     ylabel('<J_z + I_z>')
%     title({'Expectation value of J_z + I_z', 'in z direction'})
end
%% Save the susceptibilities
end
function save_file(fields,freq_total,rechi,imchi,ttt,theta,phi,location)

% x1x = squeeze(rechi1.x);
% x2x = squeeze(imchi.x);
% save(strcat('LiHoF4_x1x_x2x_',sprintf('%1$3.3fK_%2$uDeg', ttt, theta),'fields','freq_total','x1x','x2x'));
%
% x1y = squeeze(rechi1.y);
% x2y = squeeze(imchi.y);
% save(strcat('LiHoF4_x1y_x2y_',sprintf('%1$3.3fK_%2$uDeg', ttt, theta),'fields','freq_total','x1y','x2y'));

x1z = squeeze(rechi.z);
x2z = squeeze(imchi.z);
sname = strcat('LiHoF4_x1z_x2z_',sprintf('%1$3.3fK_%2$.1fDeg_%3$.1fDeg.mat', ttt, theta, phi));
savefile = fullfile(location,sname);
save(savefile,'fields','freq_total','x1z','x2z','-v7.3');
end

function [fields, freq_total, rechi, imchi, JIz_exp]=linear_response(eee,fff,ttt,vvv)
% Calculation of susceptibilities
E = eee;
V = vvv;
fields = vecnorm(fff);
freq_total = (1:0.01:5);

imchix = double.empty(length(freq_total(1,:)),length(fields(1,:)),0);
imchiy = double.empty(length(freq_total(1,:)),length(fields(1,:)),0);
imchiz = double.empty(length(freq_total(1,:)),length(fields(1,:)),0);

rechi1x = double.empty(length(freq_total(1,:)),length(fields(1,:)),0);
rechi1y = double.empty(length(freq_total(1,:)),length(fields(1,:)),0);
rechi1z = double.empty(length(freq_total(1,:)),length(fields(1,:)),0);

J=8;
I=3.5;
gLande_Ho=1.25;
ELEf = gLande_Ho * 0.05788;     % Lande factor * Bohr magneton (meV T^-1)
NUCf = 4.173 * 3.15245e-5;   % Nuclear Lande factor, mu/mu_N = 4.173
% NUCf = 4.732 * 3.1519e-5;   % Original code
gama = 0.02; % define lifetime (meV)
% gama = 0.00005; 
    
%Initiate J operators
Jz=diag(J:-1:-J); % Jz = -J, -J+1,...,J-1,J
Jzh=kron(Jz,eye(2*I+1)); % Expand Jz space to include nuclear degree of freedom
Jp=diag(sqrt((J-((J-1):-1:-J) ).*(J+1+( (J-1):-1:-J) )),1); % electronic spin ladder operator
Jm=Jp'; % electronic spin ladder operator
Jph=kron(Jp,eye(2*I+1)); % Expand Hilbert space
Jmh=kron(Jm,eye(2*I+1)); % Expand Hilbert space
Jxh=(Jph+Jmh)/2;
Jyh=(Jph-Jmh)/2i;
%tensor product of cristal field to include nuclear moments
%Initiate I operators
Iz=diag(I:-1:-I); %Iz = -I, -I+1,...,I-1,I
Izh=kron(eye(2*J+1),Iz); % Expand Hilbert space
Ip=diag(sqrt((I-((I-1):-1:-I)).*(I+1+((I-1):-1:-I))),1); % Nuclear spin ladder operator
Im=Ip'; % Nuclear spin ladder operator
Iph=kron(eye(2*J+1),Ip); % Expand to match the Hilbert space
Imh=kron(eye(2*J+1),Im); % Expand to match the Hilbert space
Ixh=(Iph+Imh)/2;
Iyh=(Iph-Imh)/2i;

% Single out <Jz+Iz> calculations
% Jz_exp = double.empty(0,length(fields(1,:)));
JIz_exp = double.empty(0,length(fields(1,:)));
for kk = 1:length(fields(1,:)) % calculate susceptibility for all fields
    v = squeeze(squeeze(V(kk,:,:,:))); % Obtain the corresponding eigen vectors
    JzhT = Jzh * ELEf;
    IzhT = Izh * NUCf;
%     tz = v'  * JzhT * v;
    ttz  = v'  * (JzhT+IzhT) * v;
%     Jz_exp(1,kk) = sqrt(sum(sum((tz) .* (tz.'))));
    JIz_exp(1,kk) = sqrt(sum(sum((ttz) .* (ttz.'))));
end

for m = 1:length(freq_total(1,:)) %calculate susceptibility for all frequencies
    freq = freq_total (m);
    f2E = 1/241.8;  % GHz to meV
    omega = freq*f2E;   % define frequency sweep range (meV)
%     for k = 1:length(fields(1,:)) % for debugging: calculate susceptibility for all fields
    parfor k = 1:length(fields(1,:)) % calculate susceptibility for all fields
        v = squeeze ( squeeze(V(k,:,:,:)) ); % Obtain the corresponding eigen vectors
        en = squeeze ( squeeze(E(k,:,:)) ); % Obtain the corresponding eigen energies in meV
%       N = length(en);
%       chi_t = zeros(1,N^2);
%       ll = 1;
        beta = 11.6/ttt; %[meV^-1]
        z=sum(exp(-beta*en));
        zz=exp(-beta*en)/z;
        [n,np]=meshgrid(zz,zz);
        NN=n-np;
        [ee,eep]=meshgrid(en,en);
        EE = eep-ee-omega;
        gamma = ones(size(EE))*gama;
        G = gamma ./ (EE.^2 + gamma.^2); 
        G1 = EE ./ (EE.^2 + gamma.^2);  

        JxhT = Jxh * ELEf;
        IxhT = Ixh * NUCf;

        JyhT = Jyh * ELEf;
        IyhT = Iyh * NUCf;

        JzhT = Jzh * ELEf;
        IzhT = Izh * NUCf;

% Calculate susceptibilities along x-axis
        ttx  = v'  * (JxhT+IxhT) * v; 
        chi_tx  = (ttx) .* (ttx.') .* NN .* G;
        chi_t1x = (ttx) .* (ttx.') .* NN .* G1;
        sss=sum(sum(chi_tx)); 
        sss1=sum(sum(chi_t1x));
        imchix  (m,k,1) =  real(sss)   ;
        rechi1x (m,k,1) =  real(sss1)  ; 

% Calculate susceptibilities along y-axis
        tty  = v'  * (JyhT+IyhT) * v; 
        chi_ty  = (tty) .* (tty.') .* NN .* G;
        chi_t1y = (tty) .* (tty.') .* NN .* G1;
        sss=sum(sum(chi_ty));
        sss1=sum(sum(chi_t1y));
        imchiy  (m,k,1) =  real(sss)   ;
        rechi1y (m,k,1) =  real(sss1)  ;   

% Calculate susceptibilities along c-axis
        ttz  = v'  * (JzhT+IzhT) * v;
        chi_tz  = (ttz) .* (ttz.') .* NN .* G;
        chi_t1z = (ttz) .* (ttz.') .* NN .* G1;
        sss=sum(sum(chi_tz));
        sss1=sum(sum(chi_t1z));
        imchiz  (m,k,1) =  real(sss)   ;
        rechi1z (m,k,1) =  real(sss1)  ;
    end
end
imchi.x = imchix;
imchi.y = imchiy;
imchi.z = imchiz;

rechi.x = rechi1x;
rechi.y = rechi1y;
rechi.z = rechi1z;
end