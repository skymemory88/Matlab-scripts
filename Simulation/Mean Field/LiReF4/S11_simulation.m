function sim = S11_simulation(mion,scanMode,dscrt_var,theta,phi,gama)
%% Strong/weak coupling simulation
% B = linspace(0,9,100);
% Br = 3.72;
% m = 7/2;
% hbar = 1;
% Delt = -m*(B-Br)/hbar;
% wc = 3.674;
% gc = 0.055;
% kappa_e = 0.1;
% w = wc - gc^2.*Delt./(Delt.^2+kappa_e^2);
% hp1 = plot(B,w,'-k','LineWidth',3);
% hold on
% wp = wc + Delt./2 + sqrt(Delt.^2+4*gc^2)/2;
% wm = wc + Delt./2 - sqrt(Delt.^2+4*gc^2)/2;
% hp2 = plot(B,wm,'-r',B,wp,'-r','LineWidth',3);
%% Bare cavity S11 Simulation
% clearvars
% kappa = 0.2;
% kappa_e = -2.0*kappa^2;
% wc = 3.5;
% w = linspace(3,4,101);
% S11 = 1+kappa_e./(1i.*(w-wc)-kappa_e/2);
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
clearvars -except dscrt_var scanMode theta phi gama RPA_opt dE mion

mu0 = 4*pi*10^-7; % Vacuum permeability ([H/m])
hbar = 1.055E-34; % Reduced Planck constant [J.s]
meV2J = 1.602217e-22; % [J/meV]
f2E = hbar*2*pi*10^9/meV2J; % [meV/GHz]
kB = 0.08617; % [meV/K]
rho = 4e30/(5.17*5.17*10.75); % Holmium (magnetic moment) number density [m^-3]

% frequency parameter setup
dE = -0.0;
% wc = 4.7565 +dE; % Fundamental mode of the cavity
wc = 3.643+dE;
% wc = 3.3+dE;
% w2 = 4.2; % Second mode of the cavity
% freq_l = 4.715;
% freq_h = 4.775;
freq_l = 3.52;
freq_h = 3.74;
% freq_l = wc-0.1;
% freq_h = wc+0.1;
freq_pts = 2001; % number of points along frequency axis
freq = linspace(freq_l,freq_h,freq_pts); % Only applies when calculating from scratch

Options.scanMode = scanMode; % continuous variable choice: 1. Field, 2. Temperature
    cVar_l = 0.0; % continuous variable setup
    cVar_h = 17.0;
    cVar_np = 1201; % number of points along field axis
    cVar = linspace(cVar_l,cVar_h,cVar_np); % Only applies when calculating from scratch
    
Options.simType = 2; % Analysis options (1) Perturbation (2) Load MF/RPA susceptibilities (3) MF/RPA calculation
Options.Ediff = false; % Plot energy levels on top
    hyp = 1.0; % hyperfine isotope proportionality
Options.RPA = true; % Use RPA susceptibilities
Options.noise = false; % Add white noises to the backgroud
    nlevel = 30; % signal-to-noise level (dB)
Options.x1x2 = false; % Plot the matrix elements of the susceptibilities
Options.trace = false; % Calculate the trace of resonant frequency along field axis
Options.Q_1 = false; % Calculate 1/Q plot along the field axis
Options.plot = true; % Option to plot data
Options.savedata = false; % Save results to a file
Options.savegif = false; % Save a series of color plots as a gif
    saveloc = 'G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Simulations\Matlab\Susceptibilities\S11 parameters';
    location = ['G:\My Drive\File sharing\PhD program\Research projects\Li',mion,...
            'F4 project\Data\Simulations\Matlab\Susceptibilities\with Hz_I'];
        
switch Options.scanMode
    case 'field'
        file_part = sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_hp=%4$.2f.mat',dscrt_var,theta,phi,hyp);
        file_part2 = sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_%4$.2e.mat', dscrt_var, theta, phi, gama);
        filename = ['Hscan_LiHoF4_',file_part];
    case 'temp'
        file_part = sprintf('%1$3.3fT_%2$.2fDg_%3$.1fDg_hp=%4$.2f.mat',dscrt_var,theta,phi,hyp);
        file_part2 = sprintf('%1$3.3fT_%2$.2fDg_%3$.1fDg_%4$.2e.mat', dscrt_var, theta, phi, gama);
        filename = ['Tscan_LiHoF4_',file_part];
end
MF_file = fullfile(location,filename); % Mean field data to load

if Options.RPA ==true
    
    filename = strcat('RPA_LiHoF4_',file_part2);
else
    filename = strcat('chi0_LiHoF4_',file_part2);
end
chi_file = fullfile(location,filename); % susceptibility data to load

alpha = 0.0; % Phase angle (in radians)
% alpha = [0 pi/6 pi/4 pi/2 pi/3 pi]; % Phase angle (in radians) between coherent and dissipative couplings (can be an array)
% alpha = linspace(0,pi,20);

chi_labl = ["\chi_{xx}","\chi_{xy}","\chi_{xz}"
            "\chi_{yz}","\chi_{yy}","\chi_{yz}"
            "\chi_{zx}","\chi_{zy}","\chi_{zz}"]; % element label of susceptibility tensor
chi_elem = [0 0 0
            0 0 0
            0 0 1]; % tensorial element of susceptibility to use
chi_idx = find(chi_elem); % index of chi element choice.

if Options.RPA == true
%     scale = 2.1; % MF-RPA scaling factor
    scale = 2.1; % MF-RPA scaling factor
else
%     scale = 2.3; % MF scaling factor
    scale = 2.4; % MF scaling factor
end
% Filling factor:  %SC108: 0.112, SC200: 0.0127
filFctr = 0.0127*scale; % SC200
% filFctr = 0.112*1.5; % SC108
% filFctr_1 = 0.004;
% filFctr_2 = 0.15; % Filling factor for the second mode
% filFctr = [0.005 0.01 0.015 0.02 0.025 0.03 0.035 0.04];

% Cavity loss rates
% kappa_i = 7.69e-4; % SC200: 2020_10_0009.dat
% kappa_e = 9.34e-4; % SC200: 2020_10_0009.dat
% kappa_i = 2.64e-4; % SC200: 2020_11_0031.dat
% kappa_e = 7.02e-4; % SC200: 2020_11_0031.dat
% kappa_i = 1.07e-4; % SC200: 2020_11_0032.dat
% kappa_e = 8.16e-4; % SC200: 2020_11_0032.dat

% kappa_i = 5.74e-4; % SC239: 2021_06_0040.dat
% kappa_e = 3.66e-4; % SC239: 2021_06_0040.dat
% kappa_i = 3.81e-4; % SC239: 2021_07_0003.dat
% kappa_e = 5.56e-4; % SC239: 2021_07_0003.dat

% kappa_i = 1.77e-4; % SC127_1
% kappa_e = 6.92e-4; % SC127_1
% kappa_i = 1.15e-4; % SC108
% kappa_e = 8.18e-4; % SC108
kappa_i = 2.91e-4;
kappa_e = 4.02e-4;

% coupling strength
gw0 = sqrt(mu0*wc*10^9*2*pi*rho*filFctr/hbar/2); % susceptibility prefactor [T.(J.s)^-1]
gw0 = gw0 * meV2J * f2E * 10^-9; % [T]
% gw0 = filFctr^2*wc^2; % Phys.Rev.Appl. 2, 054002 (2014)


if Options.savegif == true
%     im_t = numel(alpha); % a series with varying phase angle
    im_t = numel(filFctr); % a series with varying filling factor
    im_idx = 1;
    S11_frame = []*im_t;
end

while true
    switch Options.simType
        case 1 % Option 1: Perturbative treatment of resonant frequency
            beta = 1/(analysis.temp*kB);
            g = 0.001; % Coupling strength measured against the ground state energy
            % g = 0.01*wc;
            load(MF_file,'-mat','eee','fff');
            continu_var = vecnorm(fff);
            En(:,:)=squeeze(eee)/f2E;
            %En(:,:) = eee(:,1,:)/f2E;
            Ediff = double.empty(7,size(En,1),0);
            bzF = double.empty(size(Ediff,1),size(continu_var,2),0);
            for ii = 1:7
                Ediff(ii,:,1)=En(:,ii+1)-En(:,ii);
                bzF(ii,:,1) = exp(-Ediff(ii,:)*beta); % Boltzmann factor for each transition line
            end
            
            %Optional: Plot the lowest eight energy levels
            if Options.Elevel == true
                color = ["black","red","blue","magenta","green","yellow","cyan"];
                theta = [0]; % Angle (in degrees) in the transverse field plane
                marker = [":","-.","--","-"];
                figure
                hold on
                plot(continu_var,En(:,1:8),'Color',[35 107 142]/255,'linewidth',2);
                set(gca,'fontsize',15);
                % text(8,(j-1)*D+0.02,'0.1');
                %         set(fig1,'position',[100 100 600 300])
                xlabel('Field (T)','FontSize',15)
                ylabel('Energy (GHz)','FontSize',15)
                grid off
                box on;
                xlim([0 9]);
    %            ylim([2 5]);
                tit1='Energy levels';
                title(tit1,'FontSize',15)
                legend(num2str(analysis.temp*1000,'T = %u mK,  A = A_{th}'))
            end
            
            % Optional: Energy difference bewteen neighbour levels
            if Options.Ediff == true
                figure
                hold on
                figs(:, iter2) = plot(continu_var, Ediff, 'Marker', 'none', 'LineStyle',marker(1), 'Color',color(iter2),'LineWidth',1.5);
                hold on
    %             lg(iter, iter2) = [num2str(analysis.temp(iter)*1000,"T = %u mK, A = A_{th}"),num2str(theta(iter2),", theta = %u")];
                set(gca,'fontsize',15,'Xtick',0:1:max(continu_var));
                xlabel('Field (T)','FontSize',15);
                ylabel('Energy difference (GHz)','FontSize',15);
                grid off
                box on;
    %            xlim([min(continu_var) max(continu_var)]);
    %            ylim([0.9*min(Ediff,[],'All') 1.1*max(Ediff,[],'All')]);
                title('Difference between energy levels','FontSize',15);
            end
            
            % Calculate perturbed discpersion relations
            w = double.empty(size(Ediff,1)+1,size(cVar,2),0); % eigen-frequencys
            % boltFac = double.empty(size(Ediff,1),size(field,2),0); % Boltzmann factor
            temp_Ediff = double.empty(size(Ediff,1),size(cVar,2),0);
            w(1,:,1) = wc;
            for ii = 2:size(Ediff,1)+1
                w(ii,:) = interp1(continu_var,Ediff(ii-1,:),cVar);
                temp_Ediff(ii-1,:,1) = interp1(continu_var,Ediff(ii-1,:),cVar);
%                boltFac(ii-1,:,1) = interp1(continu_var,bzF(ii-1,:),field);
            end
            Ediff = temp_Ediff;
            clearvars temp_Ediff

            % Construct interaction Hamiltonian along magnetic field axis
            ww = double.empty(0,size(Ediff,1)+1,size(cVar,2));
            for jj = 1:size(cVar,2)
                % off diagonal elements (coupling strength)
                gs = g./abs(wc-w(2:end,jj));
%                 gs = g./abs(wc-w(2:end,jj)).*(boltFac(:,jj)/max(boltFac(:,jj)));
%                 gs = g/abs(wc-w(2,jj))*ones(1,size(Ediff,1));
%                 gs = g/abs(wc-w(2,jj))*ones(1,size(Ediff,1))*(boltFac(:,jj)/max(boltFac(:,jj)));

                % check point: use to turn off interactions bewteen cavity and certain spin levels
%                 gs(1:2) = 0;
%                 gs(1:end-1) = 0;

                Ham = diag(w(:,jj));
                Ham(1,2:end) = gs;
                Ham(2:end,1) = gs';
                ww(1,:,jj) = eig(Ham);
            end
            % pick out the frequencies along the field axis that is closes to bare frequency of the cavity
            wn = abs(squeeze(ww)-wc);
            wc = ww(wn == min(wn,[],1));
            clearvars wn ww

            figure
            plot(cVar,wc,'--k','linewidth',1);
            hold on
            plot(cVar,Ediff,'-r');
            xlim([cVar_l cVar_h])
            ylim([freq_l freq_h])

            wc = repmat(wc,1,size(freq,2)); % expand the resonant frequency to a matrix for calculations in the next step
            [cVar,freq] = meshgrid(cVar,freq);
            break;
        case 2 % Option 2 Load existing susceptibilities and interpolate
            if isfile(MF_file) && isfile(chi_file)
                load(MF_file,'-mat','eee'); % loads eigenEnergy
                En = squeeze(eee);
                switch Options.scanMode
                    case 'field'
                        load(chi_file,'-mat','fields','freq_total','chi');
                        continu_var = fields;                
                        [cVar,freq] = meshgrid(continu_var(1,:),freq_total);
                        sim.freq = freq;
                        sim.field = cVar;
                        xlab = 'Magnetic Field (T)';
                        titl = num2str(dscrt_var*1000, '%u mK');
                    case 'temp'
                        load(chi_file,'-mat','temp','freq_total','chi');
                        continu_var = temp;                
                        [cVar,freq] = meshgrid(continu_var(1,:),freq_total);
                        sim.freq = freq;
                        sim.temp = cVar;
                        xlab = 'Temperature (K)';
                        titl = num2str(dscrt_var, '%.2f T');
                end
                chi = reshape(chi,[],size(chi,3),size(chi,4));
                
                chi_elem = find(chi_elem);
                x1 = real(squeeze(chi(chi_elem,:,:)));
                x2 = imag(squeeze(chi(chi_elem,:,:)));
                %% Calculate excitation spectrum
                if Options.Ediff == true
                    Ediff = double.empty(0,size(En,1));
                    for i=1:7 % upper limit = size(En,2)-1
                        Ediff(i,:) = (En(:,i+1)-En(:,i))/f2E;
                    end
                end
                clearvars fff vvv eee
                
                if Options.x1x2 == true
                    figure;
                    hp1 = pcolor(continu_var(1,:),freq,x1);
                    set(hp1, 'edgeColor','none')
%                     caxis([-23 2]);
                    colorbar
                    xlabel(xlab)
                    ylabel('Frequency (GHz)')
                    if Options.RPA == true
                        title({'Real part of \chi_{zz} with RPA'})
                    else
                        title({'Real part of \chi_{zz}'})
                    end
                    
                    figure
                    hp2 = pcolor(continu_var(1,:),freq,(x2));
                    set(hp2, 'edgeColor','none')
                    colorbar
                    xlabel(xlab)
                    ylabel('Frequency (GHz)')
                    if Options.RPA == true
                        title({strcat("Imaginary part of ", chi_labl(chi_idx), " with RPA")})
                    else
                        title({strcat("Imaginary part of ", chi_labl(chi_idx))})
                    end
                end
                break;
            else % if file not exist, call calculation function first
                prompt = sprintf('Not all required data file exists, run calculation and generate the files?\n');
                answer = input(prompt,'s');
                switch lower(answer)
                    case {'y','yes'}
                        Options.simType = 3;
                    case {'n','no'}
                        return
                    otherwise
                        disp('Input error! Does not recognize your answer')
                        return
                end
            end
        case 3 % Option 3: Calculate susceptabilities for resonant frequency shift (takes long time)
            LiReF4_MF_Yikai('Ho',dscrt_var,cVar,theta,phi);
            MF_RPA_Yikai('Ho',dscrt_var,freq,theta,phi,gama,1,Options.RPA);
            
            load(chi_file,'-mat','chi');
            chi = reshape(chi,[],size(chi,3),size(chi,4));
            
            chi_elem = find(chi_elem);
            x1 = real(squeeze(chi(chi_elem,:,:)));
            x2 = imag(squeeze(chi(chi_elem,:,:)));
            sim.freq = freq;
            sim.field = cVar;
            %% Calculate the expecation value of the spin moment
%             load(MF_file,'-mat','eee','vvv','fff','ion'); % loads variables "Energy" and "EigenVector", and "Fields"
%             V = vvv;
%             field = vecnorm(fff);
%             %Initiate J operators
%             J=8; % Electronic moment for Ho3+
%             I=3.5; % Nuclear moment for Ho3+
%             Jz=diag(J:-1:-J); % Jz = -J, -J+1,...,J-1,J
%             Jzh=kron(Jz,eye(2*I+1)); % Expand Jz space to include nuclear degree of freedom
%             %     Iz=diag(I:-1:-I); %Iz = -I, -I+1,...,I-1,I
%             %     Izh=kron(eye(2*J+1),Iz); % Expand Hilbert space
%             Jz_exp = double.empty(0,length(continu_var(1,:)));
%             %     JIz_exp = double.empty(0,length(continu_var(1,:)));
%             for kk = 1:length(continu_var(1,:)) % calculate susceptibility for all continu_var
%                 v = squeeze(squeeze(V(kk,:,:,:))); % Obtain the corresponding eigen vectors
%                 JzhT = Jzh * ELEf;
%                 %         IzhT = Izh * NUCf;
%                 tz = v'  * JzhT * v;
%                 %         ttz  = v'  * (JzhT+IzhT) * v;
%                 %         Jz_exp(1,kk) = sqrt(sum(sum((tz) .* (tz.'))));
%                 Jz_exp(1,kk) = real(sum(diag(tz)));
%                 %         JIz_exp(1,kk) = sqrt(sum(sum((ttz) .* (ttz.'))));
%                 %         JIz_exp(1,kk) = real(sum(diag(ttz)));
%             end
%             Jz_exp = interp1(continu_var,Jz_exp,field);
%             % JIz_exp = interp1(continu_var,JIz_exp,field);
            break;
    end
end
%% Calculate the spin terms in the denominator without susceptibilities
%   load(MF_file,'-mat','eee','vvv'); % loads variables "Energy" and "EigenVector", and "Fields"
%   Ediff = double.empty(0,size(En,1));
%   for i=1:7 % upper limit = size(En,2)-1
%      Ediff(i,:) = (En(:,i+1)-En(:,i))/f2E;
%   end
%   Gamma = 150*kappa_i; % Coupling strength between the cavity field and the spin system
%   life = 1/(200*abs(kappa_e)); % define spin level lifetime (meV-1)
%   % life = 0.0001;
%   % life = 0.00005;
%   spins = double.empty(size(Ediff,2),0); % A container for the elements of spin term summation in the denominator of S11
%   spin_term = double.empty(size(Ediff,2),size(freq,2),0); % A container for the spin term summation in the denominator of S11
%   popul = double.empty(size(Ediff,1),size(Ediff,2),0); % A container for the population factor of the transition levels
%   for ii = 1:size(field,1)
%       for kk = 1:size(Ediff,1)
%    %        popul(kk,ii,1) = exp(-Ediff(kk,ii).*f2E./(kB*analysis.temp))./exp(-min(Ediff(kk,ii)).*f2E./(kB*analysis.temp)); % occupation factor
%    %        spins(kk) = JIz_exp(1,ii)*2*Gamma^2.*sqrt(popul(kk,ii,1))./(1i.*(Ediff(kk,ii)-wc(ii,1)) + 1/life); % calculate the spin term for each transition level
%           spins(kk) = 2*Gamma^2.*Jz_exp(1,ii)./(1i.*(Ediff(kk,ii)-wc(ii,1)) + 1/life);
%    %        spins(kk) = 2*Gamma^2.*sqrt(popul(kk,ii,1))./(1i.*(Ediff(kk,ii)-wc(kk,ii)) + 3*abs(kappa_e)); % calculate the spin term for each transition level
%       end
%        spin_term(ii,:,1) = sum(spins(:)); % Sum over all levels
%   end
% S11_2 = 1 + 2*kappa_e./(1i.*(wc-freq) - kappa_i - kappa_e + spin_term); % Calculate the S11 response using spin turms
% S11_2 = 1 + 2*kappa_e./(1i.*(wc-freq) - kappa_i - kappa_e + Gamma.*x2); % Calculate the S11 response using susceptibilites
%% Simulat S11 using input-output formalism
for ii = 1:length(alpha) % varying phase angle
    chi = (x1+1i*x2); % [meV.T^-2]
    S11 = 1 + 2*kappa_e./(1i.*(freq-wc) - (kappa_i + kappa_e) + gw0^2 * 1i*chi * exp(1i*alpha(ii)));
%     S21_1 = kappa_e./(1i.*(freq-wc) - (kappa_i + kappa_e) + filFctr_1*1i*chi*exp(1i*alpha(ii)));
%     S21_2 = kappa_e_2./(1i.*(freq-w2) - (kappa_i_2 + kappa_e_2) + 1i*filFctr_2*chi*exp(1i*(alpha(ii))));
%     S21 = abs(S21_1 + S21_2);
    if Options.noise == true
        fprintf('Adding white noise (%d dB) to the background.\n',nlevel);
        S11 = awgn(S11,nlevel,'measured');
%         S11 = S11 + randn(size(S11))*0.01;
    end
    sim.S11{ii} = S11;
    if Options.plot == true
        s_para = figure;
        map = pcolor(cVar,freq,mag2db(abs(S11))); % color plot of the S11 response
%         map = pcolor(field,freq,mag2db(real(S11))); % color plot of the S11 response        
%         map = pcolor(field,freq,mag2db(S21)); % color plot of the S21 response
        map.EdgeColor = 'none';
        colorbar
        xlim([cVar_l cVar_h])
        ylim([freq_l freq_h])
        xlabel(xlab);
        ylabel('Frequency (GHz)');
        if Options.RPA == true
            title(['MF-RPA simulation of S11 at ', titl])
        else
            title(['MF simulation of S11 at ', titl])
        end
        caxis('auto')
%         caxis([-30 2])
        if Options.Ediff == true
            hold on
            plot(cVar(1,:), Ediff(:,:), ':k', 'linewidth', 1.3); % plot the transition levels on top of the S11 color map
        end
        xlim([cVar_l cVar_h])
        ylim([freq_l freq_h])
    end

   if Options.trace
       f0 = zeros(size(cVar,2),2);
       mag11 = abs(S11);
       for jj = 1:size(mag11,2)-1
           f = freq(:,jj);
           f0(jj) = f(mag11(:,jj)==min(mag11(:,jj)));
       end
       sim.f0 = f0;
       if Options.plot == true
           figure
           plot(cVar(1,:),f0,'-k','lineWidth',1.5);
           xlim([cVar_l cVar_h])
           ylim([freq_l freq_h])
           xlabel(xlab)
           ylabel('Frequency (GHz)')
       end
    end
    if Options.Q_1 == true
        Q0 = zeros(size(cVar,2),2);
        f0 = zeros(size(cVar,2),2);
        mag11 = abs(S11);
        for jj = 1:size(mag11,2)
            [~,fidx] = min(mag11(:,jj));
            f0(jj) = freq(fidx,jj);
            Q0(jj) = f0(jj)/(kappa_e+kappa_i+imag(chi(fidx,jj)));
        end
        sim.Qf = Q0;       
        mag = abs(S11)';
        FWHM = zeros(size(mag,1),1);
        %find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
        for jj = 1:size(cVar,2) %Searching column minima (fixed field)
            [~,idx] = min( mag(jj,:) );
            HM = ( max(mag(jj,:)) + min(mag(jj,:)) )/2; %
            % Calculate quality factor using f0/FWHM
            left = mag(jj,1:idx);
            right = mag(jj,idx:end);
            [chi_idx,~] = find(left >= HM, 1, 'last');
            [ridx,~] = find(right >= HM, 1, 'first');
            if chi_idx == idx+ridx-1
                FWHM(jj) = Inf;
            else
                FWHM(jj) = freq(idx+ridx-1,1)-freq(chi_idx,1);
            end
        end
        sim.FWHM = FWHM;
        if Options.plot == true
            figure
            plot(cVar(1,:),1./Q0,'-')
            xlabel(xlab)
            ylabel('1/Q')
            figure
            plot(cVar(1,:),FWHM,'-')
            xlabel(xlab)
            ylabel('FWHM (MHz)')
        end
    end
    if Options.savegif == true
        drawnow
        frame = getframe(s_para);
        S11_frame{im_idx} = frame2im(frame);
        im_idx = im_idx + 1;
    end
end

if Options.savegif == true
    filename = strcat('S11_LiHoF4_', file_part, sprintf('alpha=%1$d-%2$d.gif',min(alpha)*180/pi,max(alpha)*180/pi)); % Varying phase angle
%     filename = strcat('S11_LiHoF4_',sprintf('%1$3.3fK_%2$.1fDeg_%3$.1fDeg_%4$.2e',analysis.temp,theta,phi,gama),...
%         sprintf('FF=%1$.2f-%2$.2f.gif',min(filFctr),max(filFctr))); % Varying filling factor
    fileobj = fullfile(saveloc,filename);
    for ii = 1:length(S11_frame)
        [img,cmp] = rgb2ind(S11_frame{ii},256);
        if ii == 1
            imwrite(img,cmp,fileobj,'gif','LoopCount',inf,'DelayTime',1);
        else
            imwrite(img,cmp,fileobj,'gif','WriteMode','append','DelayTime',1);
        end
    end
end

if Options.savedata == true
    switch Options.scanMode
        case 'field'
            file_part = sprintf('%1$3.3fK_%2$.2fGHz_%3$.2fDeg_%4$.1fDeg_%5$.2e',dscrt_var,wc,theta,phi,gama);
            analysis.temp = dscrt_var; % Simulation temperature(s)
        case 'temp'
            file_part = sprintf('%1$3.3fT_%2$.2fDg_%3$.1fDg_hp=%4$.2f.mat',dscrt_var,theta,phi,hyp);
            analysis.field = dscrt_var; % Simulation field(s)
    end
    savename = strcat('S11_LiHoF4_',file_part);
    
    % save the simulation parameters
    analysis.scanMode = Options.scanMode;
    analysis.wc0 = wc;
    analysis.alpha = alpha;
    analysis.kpi0 = kappa_i;
    analysis.kpe0 = kappa_e;
    analysis.gamma0 = gama;
    analysis.fil = filFctr;
    analysis.dType = 'sim';
    analysis.RPA = RPA_opt;
    analysis.name = savename;
    if length(alpha) == 1
        S11 = sim.S11{:};
    else
        S11 = sim.S11{1}; % take only the first struct element for quick plots
    end

    saveObj = fullfile(saveloc,savename);
    save([saveObj,'.mat'],'cVar','freq','S11','analysis','-v7.3');
end
clearvars -except Options sim
end