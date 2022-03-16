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

mu0 = 4*pi*1e-7; % Vacuum permeability ([H/m])
hbar = 1.05457E-34; % Reduced Planck constant [J.s/rad]
meV2J = 1.60218e-22; % [J/meV]
Gh2mV = hbar*2*pi*10^9/meV2J; % [meV/GHz]
kB = 8.61733e-2; % [meV/K]
rho = 4/(5.175e-10*5.175e-10*10.75e-10); % magnetic moment number density in LiReF4 [m^-3]

% frequency parameter setup
dE = -0.0;
% wc = 4.7565 +dE; % Fundamental mode of the cavity
fc = 3.643+dE; % cavity resonance frequency GHz
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
    cVar_np = 1601; % number of points along field axis
    cVar = linspace(cVar_l,cVar_h,cVar_np); % Only applies when calculating from scratch
    
Options.simType = 3; % Analysis options (1) Perturbation (2) Load MF/RPA susceptibilities (3) MF/RPA calculation
Options.nZee = false;
Options.RPA = true; % Use RPA susceptibilities
Options.noise = false; % Add white noises to the backgroud
    sig2n = 30; % signal-to-noise level (dB)
Options.plot = false; % Option to plot data
    Options.x1x2 = false; % Plot the matrix elements of the susceptibilities
    Options.trace = true; % Calculate the trace of resonant frequency along field axis
    Options.Q_1 = false; % Calculate 1/Q plot along the field axis
    Options.Elevel = false; % Plot energy eigenstates
    Options.Ediff = false; % Plot energy levels on top
        nLevel = 7; % number of excitation spectrum (max 135);
        ndiff = 1; % order of excitations (min 1)
        nMax = nLevel + ndiff; % Hilbert space dimension
        hyp = 1.0; % hyperfine isotope proportionality
Options.savedata = false; % Save results to a file
Options.savegif = false; % Save a series of color plots as a gif

if Options.nZee == true
    nZeePath = 'with Hz_I';
else
    nZeePath = 'without Hz_I';
end
location = ['G:\My Drive\File sharing\PhD program\Research projects\Li',mion,...
    'F4 project\Data\Simulations\Matlab\Susceptibilities\'];
saveloc = [location,'\S11 parameters\', nZeePath];

switch Options.scanMode
    case 'field'
        file_part = sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_hp=%4$.2f.mat',dscrt_var, theta, phi, hyp);
        file_part2 = sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_%4$.2e.mat', dscrt_var, theta, phi, gama);
        filename = ['Hscan_LiHoF4_',file_part];
    case 'temp'
        file_part = sprintf('%1$3.3fT_%2$.2fDg_%3$.1fDg_hp=%4$.2f.mat',dscrt_var, theta, phi, hyp);
        file_part2 = sprintf('%1$3.3fT_%2$.2fDg_%3$.1fDg_%4$.2e.mat', dscrt_var, theta, phi, gama);
        filename = ['Tscan_LiHoF4_',file_part];
end
MF_file = fullfile([location, nZeePath],filename); % Mean field data to load

if Options.RPA ==true
    
    filename = strcat('RPA_LiHoF4_',file_part2);
else
    filename = strcat('chi0_LiHoF4_',file_part2);
end
chi_file = fullfile([location, nZeePath],filename); % susceptibility data to load

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
    scale = 0.4;
else
    scale = 0.5;
end
% Filling factor:  %SC108: 0.112, SC200: 0.0127
filFctr = 0.0127*scale; % SC200
% filFctr = 0.112*scale; % SC108
% filFctr = 0.0127;
% filFctr_2 = 0.15; % Filling factor for the second mode
% filFctr = [0.005 0.01 0.015 0.02 0.025 0.03 0.035 0.04];

% coupling strength
gw0 = sqrt(mu0*2*pi*fc*10^9*rho*filFctr/hbar/2) / (2*pi) * 1e-9; % susceptibility prefactor [T/J. GHz]
% gw0 = filFctr^2*wc^2; % Phys.Rev.Appl. 2, 054002 (2014)

% Cavity loss rates (GHz)
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
kappa_i = 8.43e-4;
kappa_e = 4.02e-4;

if Options.savegif == true
%     im_t = numel(alpha); % a series with varying phase angle
    im_t = numel(filFctr); % a series with varying filling factor
    im_idx = 1;
    S11_frame = []*im_t;
end

while true
    switch Options.simType
        case 1 % Option 1: Perturbative treatment of resonant frequency
            load(MF_file,'-mat','eee','fff');
            continu_var = vecnorm(fff);
            switch Options.scanMode
                case 'field'
                    sim.temp = dscrt_var; % Simulation temperature(s)
                    sim.freq = freq;
                    sim.field = cVar;
                case 'temp'
                    sim.temp = continu_var; % Simulation field(s)
                    sim.freq = freq;
                    sim.field = cVar;
            end
            beta = Gh2mV/(sim.temp.*kB);
            g = 0.01; % Coupling strength measured against the ground state energy
            % g = 0.01*wc;
            
            eee = eee - min(eee,[],2);
            En(:,:) = squeeze(eee)/Gh2mV;
            Ediff = double.empty(7,size(En,1),0);
            bzF = double.empty(size(Ediff,1),size(continu_var,2),0);
            Zn = sum(exp(-En.*beta),2); % Partition function
            for ii = 1:nLevel
                Ediff(ii,:,1) = En(:,ii+1)-En(:,ii);
                bzF(ii,:,1) = (exp(-En(:,ii).*beta) - exp(-En(:,ii+1).*beta))./Zn; % Boltzmann factor for each transition line
            end
                      
            % Calculate perturbed discpersion relations
            w = double.empty(size(Ediff,1)+1,size(cVar,2),0); % eigen-frequencys
            boltFac = double.empty(size(Ediff,1),size(cVar,2),0); % Boltzmann factor
            temp_Ediff = double.empty(size(Ediff,1),size(cVar,2),0);
            w(1,:,1) = fc;
            for ii = 2:size(Ediff,1)+1
                w(ii,:) = interp1(continu_var,Ediff(ii-1,:),cVar);
                temp_Ediff(ii-1,:,1) = interp1(continu_var,Ediff(ii-1,:),cVar);
                boltFac(ii-1,:,1) = interp1(continu_var,bzF(ii-1,:),cVar);
            end
            Ediff = temp_Ediff;
            clearvars temp_Ediff

            % Construct interaction Hamiltonian along magnetic field axis
            ww = double.empty(0,size(Ediff,1)+1,size(cVar,2));
            for jj = 1:size(cVar,2)
                % off diagonal elements (coupling strength)
%                 gs = g./abs(wc-w(2:end,jj));
                gs = g./abs(fc-w(2:end,jj)).*(boltFac(:,jj));
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
            % pick out the frequencies along the field axis that is closest to bare frequency of the cavity
            wn = abs(squeeze(ww) - fc);
            fc = ww(wn == min(wn,[],1));
            sim.wr = fc;
            clearvars wn ww

            figure
            plot(cVar,fc,'.k','linewidth',1);
            hold on
            plot(cVar,Ediff,'-r');
            xlim([cVar_l cVar_h])
            ylim([freq_l freq_h])
            break;
        case 2 % Option 2 Load existing susceptibilities and interpolate
            if isfile(MF_file) && isfile(chi_file)
                load(MF_file,'-mat','eee'); % loads eigenEnergy
                En = squeeze(eee - min(eee,[],2)); % normalize to the ground state energy [meV]
                switch Options.scanMode
                    case 'field'
                        load(chi_file,'-mat','fields','freq_total','chi','unit');
                        continu_var = fields;                
                        [cVar,freq] = meshgrid(continu_var(1,:),freq_total);
                        sim.freq = freq;
                        sim.field = cVar;
                        xlab = 'Magnetic Field (T)';
                        titl = num2str(dscrt_var*1000, '%u mK');
                    case 'temp'
                        load(chi_file,'-mat','temp','freq_total','chi','unit');
                        continu_var = temp;                
                        [cVar,freq] = meshgrid(continu_var(1,:),freq_total);
                        sim.freq = freq;
                        sim.temp = cVar;
                        xlab = 'Temperature (K)';
                        titl = num2str(dscrt_var, '%.2f T');
                end
                chi = reshape(chi,[],size(chi,3),size(chi,4));
                %% Calculate excitation spectrum
                if Options.Ediff == true
                    Ediff = diff(En,1,2)/Gh2mV; % [GHz]
                    Ediff = Ediff(:,1:nLevel);
                end
                clearvars fff vvv eee
                break;
            else % if file not exist, call calculation function first
                prompt = sprintf('Not all required data file exists, run calculation and generate the files?\n');
                answer = input(prompt,'s');
                switch lower(answer)
                    case {'y','yes'}
                        Options.simType = 4;
                    case {'n','no'}
                        return
                    otherwise
                        disp('Input error! Does not recognize your answer')
                        return
                end
            end
        case 3 % Option 3: on-resonance cavity magnon-polariton interaction only
            if ~isfile(MF_file)
                prompt = sprintf('Required data not found, run calculation and generate the files?\n');
                answer = input(prompt,'s');
                switch lower(answer)
                    case {'y','yes'}
                        LiReF4_MF_Yikai('Ho',dscrt_var,cVar,theta,phi);
                        load(MF_file,'-mat','vvv','ttt','fff','eee'); % loads "eigenfunction" "temperatures", "field", "Energy"
                    case {'n','no'}
                        return
                    otherwise
                        disp('Input error! Does not recognize your answer')
                        return
                end
            else
                load(MF_file,'-mat','vvv','ttt','fff', 'eee'); % loads "eigenfunction" "temperatures", "field", "Energy"
            end
            eee = eee - min(eee,[],2); % normalize the energy levels to the ground state
            Ediff = diff(eee,1,2);
            
            switch Options.scanMode
                case 'field'
                    continu_var = vecnorm(fff);
                    sim.temp = dscrt_var; % Simulation temperature(s)
                    sim.field = continu_var; % Simulation field(s)
                    xlab = 'Magnetic field (T)';
                    titl = num2str(dscrt_var, '%.3f K');
                case 'temp'
                    continu_var = ttt;
                    sim.temp = continu_var; % Simulation field(s)
                    sim.field = dscrt_var; % Simulation field(s)
                    xlab = 'Temperature (K)';
                    titl = num2str(dscrt_var, '%.3f T');
            end
            sim.freq = freq;
                       
            %Initiate J operators            
            %        Er       Ho      Yb      Tm      Gd      Y
            J_tab = [15/2;    8;      7/2;    6;      7/2;    1];
            L_tab = [6;       6;       3;     5;      0;      1];
            S_tab = [3/2;     2;      1/2;    1;      7/2;    1];
            I_tab = [3;      3.5;      0;     0;      1.50;   0];
            nLande = [-0.1611; 1.192; -0.2592; -0.462; -0.2265; 0]; % https://easyspin.org/documentation/isotopetable.html
                        
            elem_idx = 2; % LiReF4, Re = Ho
            J = J_tab(elem_idx); % Electronic moment for Ho3+
            I = I_tab(elem_idx); % Nuclear moment for Ho3+
            muB = 9.27401e-24; % Bohr magneton [J/T]
            muN = 5.05078e-27; % Nuclear magneton [J/T]
            ELEf = gLande(L_tab(elem_idx),S_tab(elem_idx)) * muB; % Lande factor * Bohr magneton [J/T]
            NUCf = nLande(elem_idx) * muN; % [J/T]
            Jz=diag(J:-1:-J); % Jz = -J, -J+1,...,J-1,J
            Jzh=kron(Jz,eye(2*I+1)); % Expand Jz space to include nuclear degree of freedom
            Iz=diag(I:-1:-I); % Iz = -I, -I+1,...,I-1,I
            Izh=kron(eye(2*J+1),Iz); % Expand Hilbert space
            JzhT = Jzh * ELEf; % electronic spin moment [J/T]
            IzhT = Izh * NUCf; % nuclear spin moment [J/T]
                                   
            spins = double.empty(length(continu_var), nLevel, size(freq,2),0); % A container for the elements of spin term summation in the denominator of S11
            Jz_exp = double.empty(length(continu_var), 0);
            JIz_exp = double.empty(length(continu_var), 0);
            Gc2 = double.empty(length(continu_var), nMax, nMax, 0);
            w0 = double.empty(length(continu_var), nLevel,0);
%             Gc2 = double.empty(length(continu_var),0);
            for ii = 1:length(continu_var) % calculate susceptibility for all continu_var
                En = squeeze(eee(ii,:));
                if sim.temp ~= 0
                    beta = 1/(kB*sim.temp); % [1/meV]
                    Z = sum(exp(-beta*En)); % Partition function
                    zn = exp(-beta*En)/Z; % Boltzmann weight
                else
                    zn = zeros(size(En));
                    zn(1) = 1; % set the ground state occupation to unity
                end
                zn = zn(1:nMax);
                En = En(1:nMax);
                [n,np] = meshgrid(zn,zn);
                [e,ee] = meshgrid(En,En);
                NN = n-np; % population difference among levels
                dE = e-ee; % energy gaps

                v = squeeze(vvv(ii,:,:)); % Obtain the corresponding eigen vectors
                tz = v' * JzhT * v; % [J/T]
                ttz  = v' * (JzhT+IzhT) * v; % [J/T]
                
                Jz_exp(ii,1) = zn*diag(tz(1:nMax,1:nMax));
                JIz_exp(ii,1) = zn*diag(ttz(1:nMax,1:nMax));
                Gc2(ii,:,:,1) = gw0^2 * ttz(1:nMax,1:nMax) .* ttz(1:nMax,1:nMax).' .* NN; % [GHz^2]
%                 Gc2(ii,:,:) = (squeeze(Gc2(ii,:,:)) - squeeze(Gc2(ii,:,:)).')/2; % take symmetric average to reduce numerical error
                for kk = 1:nLevel
                    w0(ii,kk,1) = dE(kk,kk+ndiff)/Gh2mV;
                    spins(ii,kk,:,1) = Gc2(ii,kk,kk+ndiff)./(1i*(w0(ii,kk) - freq) - gama); % calculate the spin term for each transition level
%                     spins(ii,kk) = Gc2(ii)./(Ediff(ii,kk) - wc(1,ii) - 1i*gama);
%                     spin_term(ii,:) = spin_term(ii,:) + squeeze(spins(ii,kk,:,1))';
                end
            end
            sim.Jz = Jz_exp;
            sim.JIz = JIz_exp;
            sim.omega = w0;
            sim.Gc2 = Gc2;
            spin_term = squeeze(sum(spins,2)); % Sum over all levels [GHz]
            break;
        case 4 % calculate MF-RPA susceptibilities from scratch (takes long time)
            LiReF4_MF_Yikai('Ho',dscrt_var,cVar,theta,phi);
            MF_RPA_Yikai('Ho',Options.scanMode,dscrt_var,freq,theta,phi,gama,1,Options.RPA);
            
            load(chi_file,'-mat','chi','unit');
            chi = reshape(chi,[],size(chi,3),size(chi,4));
            sim.freq = freq;
            sim.field = cVar;
    end
end
%% Simulat S11 using input-output formalism
switch Options.simType
    case {2,4}
        sim.unit = unit;
        for ii = 1:length(alpha) % varying phase angle
            chi_elem = find(chi_elem); % select the susceptibility tensorial component
            if strcmp(unit,'J')
                chi = squeeze(chi(chi_elem,:,:)); % [J/T^2]
            elseif strcmp(unit,'GHz')
                chi = squeeze(chi(chi_elem,:,:))*Gh2mV*meV2J; % [GHz/T^2] --> [J/T^2]
            else
                chi = squeeze(chi(chi_elem,:,:))*meV2J; % [meV/T^2] --> [J/T^2]
            end
            gw2 = gw0^2*Gh2mV*meV2J; % [GHz*T^2/J]

            S11 = 1 + 2*kappa_e./(1i.*(freq-fc) - (kappa_i + kappa_e) + gw2*1i*chi * exp(1i*alpha(ii)));
%             S21_1 = kappa_e./(1i.*(freq-wc) - (kappa_i + kappa_e) + filFctr_1*1i*chi*exp(1i*alpha(ii)));
%             S21_2 = kappa_e_2./(1i.*(freq-w2) - (kappa_i_2 + kappa_e_2) + 1i*filFctr_2*chi*exp(1i*(alpha(ii))));
%             S21 = abs(S21_1 + S21_2);
            sim.S11{ii} = S11;
            if Options.noise == true
                fprintf('Adding white noise (%d dB) to the background.\n',sig2n);
                S11 = awgn(S11,sig2n,'measured');
                %         S11 = S11 + randn(size(S11))*0.01;
            end
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
                        
                if Options.x1x2 == true
                    x1 = real(chi/meV2J); % [meV/T^2]
                    x2 = imag(chi/meV2J); % [meV/T^2]
                    figure;
                    hp1 = pcolor(continu_var(1,:), freq, x1);
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
                    hp2 = pcolor(continu_var(1,:), freq, x2);
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
                    Q0(jj) = f0(jj)/(kappa_e + kappa_i + imag(gw2*chi(fidx,jj)));
                end
                sim.Qf = Q0;
                mag = abs(S11)';
                FWHM = zeros(size(mag,1),1);
                %find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
                for jj = 1:size(cVar,2) %Searching column minima (fixed field)
                    [~,idx] = min( mag(jj,:) );
                    HM = max(mag(jj,:))*0.7; % 70% of the magnitude as Half max of power
%                     HM = ( max(mag(jj,:)) + min(mag(jj,:)) )/2; % midpoint of the magnitude
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
    case 3
        S11 = 1 + 2*kappa_e./(1i.*(freq-fc) - (kappa_i + kappa_e) + spin_term); % Calculate the S11 response using spin terms
        S11 = S11';
        sim.S11 = S11;
        [cVar,freq] = meshgrid(continu_var,freq); % Calculate the spin terms in the denominator without susceptibilities
        cVar_l = min(continu_var);
        cVar_h = max(continu_var);
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
end

% Optional: plot the simulated |S11| data
if Options.plot == true
    s_para = figure;
    map = pcolor(cVar,freq,mag2db(abs(S11))); % color plot of the S11 response
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
%     caxis([-30 2])
    if Options.Ediff == true
        hold on
        plot(cVar(1,:), Ediff(:,1:nLevel), ':k', 'linewidth', 1.3); % plot the transition levels on top of the S11 color map
    end
    xlim([cVar_l cVar_h])
    ylim([freq_l freq_h])
end

% Optional: Plot the lowest eight energy levels
if Options.Elevel == true
%     color = ["black","red","blue","magenta","green","yellow","cyan"];
    theta = [0]; % Angle (in degrees) in the transverse field plane
    marker = [":","-.","--","-"];
    figure
    hold on
    plot(continu_var,En(:,1:nLevel+1),'LineStyle',marker(1),'Color',[35 107 142]/255,'linewidth',2);
    set(gca,'fontsize',10);
    % text(8,(j-1)*D+0.02,'0.1');
%     set(fig1,'position',[100 100 600 300])
    xlabel('Field (T)','FontSize',10)
    ylabel('Energy (GHz)','FontSize',10)
    grid off
    box on;
%     xlim([0 9]);
%     ylim([2 5]);
    tit1='Energy levels';
    title(tit1,'FontSize',10)
    legend(num2str(sim.temp*1000,'T = %u mK,  A = A_{th}'))
end

% Optional: Energy difference bewteen neighbour levels
if Options.Ediff == true
    marker = [":","-.","--","-"];
    figure
    hold on
    plot(continu_var, Ediff(:,1:nLevel), 'Marker', 'none', 'LineStyle',marker(1), 'Color','k','LineWidth',1.5);
    hold on
%     lg(iter, iter2) = [num2str(sim.temp(iter)*1000,"T = %u mK, A = A_{th}"),num2str(theta(iter2),", theta = %u")];
    set(gca,'fontsize',10,'Xtick',0:2:max(continu_var));
    xlabel('Field (T)','FontSize',10);
    ylabel('Energy difference (GHz)','FontSize',10);
    grid off
    box on;
%     xlim([min(continu_var) max(continu_var)]);
%     ylim([0.9*min(Ediff,[],'All') 1.1*max(Ediff,[],'All')]);
    title('Difference between energy levels','FontSize',10);
end

if Options.savegif == true
    filename = strcat('S11_LiHoF4_', file_part, sprintf('alpha=%1$d-%2$d.gif',min(alpha)*180/pi,max(alpha)*180/pi)); % Varying phase angle
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
            file_part = sprintf('%1$3.3fK_%2$.2fGHz_%3$.2fDeg_%4$.1fDeg_%5$.2e',dscrt_var,fc,theta,phi,gama);
            sim.temp = dscrt_var; % Simulation temperature(s)
        case 'temp'
            file_part = sprintf('%1$3.3fT_%2$.2fDg_%3$.1fDg_hp=%4$.2f.mat',dscrt_var,theta,phi,hyp);
            sim.field = dscrt_var; % Simulation field(s)
    end
    savename = strcat('S11_LiHoF4_',file_part);
    
    % save the simulation parameters
    analysis.scanMode = Options.scanMode;
    analysis.wc0 = fc;
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