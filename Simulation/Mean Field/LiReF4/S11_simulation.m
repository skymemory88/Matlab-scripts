function sim = S11_simulation(ion,temp,theta,phi,gama,RPA)
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
clearvars -except temp theta phi gama RPA dE ion

sim.temp = temp; % Simulation temperature(s)
mu0 = 4*pi*10^-7; % Vacuum permeability ([H/m])
hbar = 1.055E-34; % Reduced Planck constant [J.s]
meV2J = 1.602217e-22; % [J/meV]
rho = 4e30/(5.17*5.17*10.75); % Holmium (magnetic moment) number density [m^-3]
f2E = hbar*2*pi*10^9/meV2J; % [meV/GHz]
kB = 0.08617; % [meV/K]

location = ['G:\My Drive\File sharing\PhD program\Research projects\Li',ion,...
            'F4 project\Data\Simulations\Matlab\Susceptibilities\with Hz_I'];
% location = ['/Volumes/GoogleDrive/My Drive/File sharing\PhD program\Research projects\Li',ion,...
%             'F4 project\Data\Simulations results\Matlab\Susceptibilities\with Hz_I'];

Option = 2; % Analysis options (1) Perturbation (2) Load MF/RPA susceptibilities (3) MF/RPA calculation
Options.Ediff = false; % Plot energy levels on top
    hyp = 1.0; % hyperfine isotope proportionality
Options.RPA = RPA; % Use RPA susceptibilities
Options.x1x2 = false; % Plot the matrix elements of the susceptibilities
Options.trace = false; % Calculate the trace of resonant frequency along field axis
Options.Q_1 = true; % Calculate 1/Q plot along the field axis
Options.plot = false; % Option to plot data
Options.savedata = false; % Save results to a file
Options.savegif = false; % Save a series of color plots as a gif
    saveloc = 'G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Simulations\Matlab\Susceptibilities\S11 parameters';
    
alpha = 0.0; % Phase angle (in radians)
% alpha = [0 pi/6 pi/4 pi/2 pi/3 pi]; % Phase angle (in radians) between coherent and dissipative couplings (can be an array)
% alpha = linspace(0,pi,20);
sim.alpha = alpha;

if Options.RPA == true
    scale = 2.1; % MF-RPA scaling factor
else
    scale = 2.3; % MF scaling facto
end
% Filling factor:  %SC108: 0.112, SC200: 0.0127
filFctr = 0.0127*scale; % SC200
% filFctr = 0.112*1.5; % SC108
% filFctr_1 = 0.004;
% filFctr_2 = 0.15; % Filling factor for the second mode
% filFctr = [0.005 0.01 0.015 0.02 0.025 0.03 0.035 0.04];

% Cavity loss rates
% gamma_i = 7.69e-4; % SC200: 2020_10_0009.dat
% gamma_e = 9.34e-4; % SC200: 2020_10_0009.dat
% gamma_i = 4.99e-4; % SC200: 2020_11_0032.dat
% gamma_e = 4.85e-4; % SC200: 2020_11_0032.dat
% gamma_i = 5.74e-4; % SC239: 2021_06_0040.dat
% gamma_e = 3.66e-4; % SC239: 2021_06_0040.dat
% gamma_i = 1.77e-4; % SC127_1
% gamma_e = 6.92e-4; % SC127_1
% gamma_e = 6.19e-4; % SC108
% gamma_i = 3.60e-4; % SC108
% gamma_i = 4.00e-4; % SC107
% gamma_e = 6.00e-4; % SC107
gamma_i = 4e-4;
gamma_e = 1e-3;

% field parameters
field_l = 5;
field_h = 7.5;
field_pts = 501; % number of points along field axis
field = linspace(field_l,field_h,field_pts); % sampling points along field axis
sim.field = field;

dE = -0.0;
% frequency parameters
% w0 = 4.73532; % Fundamental mode of the cavitya
w0 = 3.642+dE;
% w2 = 4.2; % Second mode of the cavity
gw0 = sqrt(mu0*w0*10^9*2*pi*rho*filFctr/hbar/2); % susceptibility prefactor [T.(J.s)^-1]
gw0 = gw0 * meV2J * f2E * 10^-9; % [T]
% gw0 = filFctr^2*w0^2; % Phys.Rev.Appl. 2, 054002 (2014)

% freq_l = 4.710;
% freq_h = 4.745;
freq_l = 3.55+dE;
freq_h = 3.73+dE;
freq_pts = 2001; % number of points along frequency axis
freq = linspace(freq_l,freq_h,freq_pts); % frequency range
sim.freq = freq;

if Options.savegif == true
%     im_t = numel(alpha); % a series with varying phase angle
    im_t = numel(filFctr); % a series with varying filling factor
    im_idx = 1;
    S11_frame = []*im_t;
end

switch Option
    case 1 % Option 1: Perturbative treatment of resonant frequency
        beta = 1/(sim.temp*kB);
        g = 0.001; % Coupling strength measured against the ground state energy
        % g = 0.01*w0;
        lname=['Hscan_LiHoF4_', sprintf('%1$3.3fK_%2$uDeg',sim.temp,theta),'.mat'];
        file = fullfile(location,lname);
        load(file,'-mat','eee','fff');
        fields = vecnorm(fff);
        E(:,:)=squeeze(eee)/f2E;
        %E(:,:) = eee(:,1,:)/f2E;
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
%            ylim([2 5]);
            tit1='Energy levels';
            title(tit1,'FontSize',15)
            legend(num2str(sim.temp*1000,'T = %u mK,  A = A_{th}'))
        end
        
        % Optional: Energy difference bewteen neighbour levels
        if Options.Ediff == true
            figure
            hold on
            figs(:, iter2) = plot(fields, Ediff, 'Marker', 'none', 'LineStyle',marker(1), 'Color',color(iter2),'LineWidth',1.5);
            hold on
            lg(iter, iter2) = [num2str(sim.temp(iter)*1000,"T = %u mK, A = A_{th}"),num2str(theta(iter2),", theta = %u")];
            set(gca,'fontsize',15,'Xtick',0:1:max(fields));
            xlabel('Field (T)','FontSize',15);
            ylabel('Energy difference (GHz)','FontSize',15);
            grid off
            box on;
%            xlim([min(fields) max(fields)]);
%            ylim([0.9*min(Ediff,[],'All') 1.1*max(Ediff,[],'All')]);
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
    case 2 % Option 2 Load existing susceptibilities and interpolate
        if Options.Ediff == true
            filename = ['Hscan_LiHoF4_',sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_hp=%4$.2f.mat',sim.temp,theta,phi,hyp)];
            file = fullfile(location,filename);
            if isfile(file)
                load(file,'-mat','eee','vvv','fff'); % loads variables "Energy" and "EigenVector", and "Fields"
            else % if file not exist, call calculation function first
                prompt = sprintf('%s doesn`t exist, run calculation and generate the file?\n',filename);
                answer = input(prompt,'s');
                switch lower(answer)
                    case {'y','yes'}
                        LiReF4_MF_Yikai(sim.temp,field,theta,phi);
                        load(file,'-mat','eee','vvv','fff'); % loads variables "Energy" and "EigenVector", and "Fields"
                    case {'n','no'}
                        return
                    otherwise
                        disp('Input error! Does not recognize your answer')
                        return
                end
            end
            E = squeeze(eee);
            fields = vecnorm(fff);
            %% Calculate the expecation value of the spin moment
            %     V = vvv;
            %     muB = 9.274e-24/meV2J; % [meV/T]
            %     muN = 3.15245e-5; % [meV/T]
            %     gLande_Ho = 1.25;
            %     ELEf = gLande_Ho * muB; 
            %     NUCf = 1.668 * muN; % gN = 1.668. https://easyspin.org/documentation/isotopetable.html
            %     % NUCf = 4.173 * 3.15245e-5;   % Nuclear Lande factor, mu/mu_N = 4.173
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
            Ediff = double.empty(0,size(E,1));
            for i=1:7 % upper limit = size(E,2)-1
                Ediff(i,:) = (E(:,i+1)-E(:,i))/f2E;
            end
            temp_Ediff = double.empty(size(Ediff,1),size(field,2),0);
            for ii = 1:size(Ediff,1)
                temp_Ediff(ii,:,1) = interp1(fields,Ediff(ii,:),field);
            end
            Ediff = temp_Ediff;
        end
        clearvars temp_Ediff fff vvv eee
        
        % Load the susceptibilities from MF-linear response calculations
        if Options.RPA ==true
            filename = strcat('RPA_LiHoF4_x1z_x2z_',sprintf('%1$3.3fK_%2$.2fDeg_%3$.1fDeg_%4$.2e.mat',sim.temp,theta,phi,gama));
%            filename = 'chi125mK_17e-5meV.mat';
        else
            filename = strcat('LiHoF4_x1z_x2z_',sprintf('%1$3.3fK_%2$.2fDeg_%3$.1fDeg_%4$.2e.mat',sim.temp,theta,phi,gama));
        end
        file = fullfile(location,filename);
        load(file,'-mat','fields','freq_total','x1z','x2z');
%         load(file,'-mat','fields','freq_total','x1x','x2x');
%         load(file,'-mat','fields','freq_total','x1y','x2y');
        
        if Options.x1x2 == true
            figure;
            hp1 = pcolor(fields(1,:),freq_total,x1z);
            set(hp1, 'edgeColor','none')
%            caxis([-23 2]);
            colorbar
            xlabel('Magnetic field (T)')
            ylabel('Frequency (GHz)')
            if Options.RPA == true
                title({'Real part of \chi_{zz} with RPA'})
            else
                title({'Real part of \chi_{zz}'})
            end
            
            figure
            hp2 = pcolor(fields(1,:),freq_total,(x2z));
            set(hp2, 'edgeColor','none')
            colorbar
            xlabel('Magnetic field (T)')
            ylabel('Frequency (GHz)')
            if Options.RPA == true
                title({'Imaginary part of \chi_{zz} with RPA'})
            else
                title({'Imaginary part of \chi_{zz}'})
            end
        end
        
        [freq,field] = meshgrid(freq,field);
        x1 = interp2(fields,freq_total,x1z,field,freq);
        x2 = interp2(fields,freq_total,x2z,field,freq);
    case 3 % Option 3: Calculate susceptabilities for resonant frequency shift (takes long time)
        filename = ['LHF_',num2str(sim.temp,'%.3f.mat')];
        file = fullfile(location,filename);
        load(file,'-mat','eee','vvv','fff'); % loads variables "Energy" and "EigenVector", and "Fields"
        if ~exist ('x1','var') && ~exist('x2','var')
            life = 1/0.00005; % define lifetime (meV-1) for hyperfine levels
            [field, Ediff, x1, x2] = linear_response(eee,fff,ttt,vvv,1/life); % Calculate susceptibilities
        end
        x1 = x1';
        x2 = x2';
        [freq,field] = meshgrid(freq,field);
%        w0 = w0./(sqrt(1+filFctr.*x1)); % Use only real part of the susceptibility
%        w0 = w0./(sqrt(1+filFctr.*1i*x2)); % Use only imaginary part of the susceptibility
%        w0 = w0./(sqrt(1+filFctr.*(x1+1i*x2)));
end
%% Calculate the spin terms in the denominator without susceptibilities
%   Gamma = 150*gamma_i; % Coupling strength between the cavity field and the spin system
%   life = 1/(200*abs(gamma_e)); % define spin level lifetime (meV-1)
%   % life = 0.0001;
%   % life = 0.00005;
%   spins = double.empty(size(Ediff,2),0); % A container for the elements of spin term summation in the denominator of S11
%   spin_term = double.empty(size(Ediff,2),size(freq,2),0); % A container for the spin term summation in the denominator of S11
%   popul = double.empty(size(Ediff,1),size(Ediff,2),0); % A container for the population factor of the transition levels
%   for ii = 1:size(field,1)
%       for kk = 1:size(Ediff,1)
%    %        popul(kk,ii,1) = exp(-Ediff(kk,ii).*f2E./(kB*sim.temp))./exp(-min(Ediff(kk,ii)).*f2E./(kB*sim.temp)); % occupation factor
%    %        spins(kk) = JIz_exp(1,ii)*2*Gamma^2.*sqrt(popul(kk,ii,1))./(1i.*(Ediff(kk,ii)-w0(ii,1)) + 1/life); % calculate the spin term for each transition level
%           spins(kk) = 2*Gamma^2.*Jz_exp(1,ii)./(1i.*(Ediff(kk,ii)-w0(ii,1)) + 1/life);
%    %        spins(kk) = 2*Gamma^2.*sqrt(popul(kk,ii,1))./(1i.*(Ediff(kk,ii)-w0(kk,ii)) + 3*abs(gamma_e)); % calculate the spin term for each transition level
%       end
%        spin_term(ii,:,1) = sum(spins(:)); % Sum over all levels
%   end
% S11_2 = 1 + 2*gamma_e./(1i.*(w0-freq) - gamma_i - gamma_e + spin_term); % Calculate the S11 response using spin turms
% S11_2 = 1 + 2*gamma_e./(1i.*(w0-freq) - gamma_i - gamma_e + Gamma.*x2); % Calculate the S11 response using susceptibilites
%% Simulat S11 using input-output formalism
for ii = 1:length(alpha) % varying phase angle
% for ii = 1:length(filFctr) % varying filling factor
%     chi = filFctr(ii)*(x1+1i*x2); % Coupling strength between the spin system and cavity field
    chi = (x1+1i*x2); % [meV.T^-2]
    S11 = 1 + 2*gamma_e./(1i.*(freq-w0) - (gamma_i + gamma_e) + gw0^2 * 1i*chi * exp(1i*alpha(ii)));
%     S21_1 = 2*gamma_e./(1i.*(freq-w0) - (gamma_i + gamma_e) + filFctr_1*1i*chi*exp(1i*alpha(ii)));
%     S21_2 = 2*gamma_e_2./(1i.*(freq-w2) - (gamma_i_2 + gamma_e_2) + 1i*filFctr_2*chi*exp(1i*(alpha(ii))));
%     S21 = abs(S21_1 + S21_2);
    sim.S11{ii} = S11;
    if Options.plot == true
        s_para = figure;
        map = pcolor(field,freq,mag2db(abs(S11))); % color plot of the S11 response
%         map = pcolor(field,freq,mag2db(real(S11))); % color plot of the S11 response        
        %     map = pcolor(field,freq,mag2db(S21)); % color plot of the S21 response
        map.EdgeColor = 'none';
        colorbar
        xlim([field_l field_h])
        ylim([freq_l freq_h])
        xlabel('Magnetic field (T)');
        ylabel('Frequency (GHz)');
        if Options.RPA == true
            title(num2str(sim.temp*1000,'MF-RPA simulation of S11 at %umK'))
        else
            title(num2str(sim.temp*1000,'MF simulation of S11 at %umK'))
        end
        caxis('auto')
%         caxis([-30 2])
        if Options.Ediff == true
            hold on
            plot(field(:,1), Ediff(:,:), ':k', 'linewidth', 1.3); % plot the transition levels on top of the S11 color map
        end
        xlim([field_l field_h])
        ylim([freq_l freq_h])
    end

   if Options.trace
       f0 = zeros(size(field,1),1);
       mag11 = abs(S11);
       for jj = 1:size(mag11,1)-1
           f = freq(jj,:);
           f0(jj) = f(mag11(jj,:)==min(mag11(jj,:)));
       end
       sim.f0 = f0;
       if Options.plot == true
           figure
           plot(field(:,1),f0,'-k','lineWidth',1.5)
           xlim([field_l field_h])
           ylim([freq_l freq_h])
           xlabel('Magnetic field (T)')
           ylabel('Frequency (GHz)')
       end
    end
    if Options.Q_1 == true
        Q0 = zeros(size(field,1),1);
        f0 = zeros(size(field,1),1);
        mag11 = abs(S11);
        for jj = 1:size(mag11,1)
            [~,fidx] = min(mag11(jj,:));
            f0(jj) = freq(jj,fidx);
            Q0(jj) = f0(jj)/(gamma_e+gamma_i+imag(chi(jj,fidx)));
        end
        sim.Qf = Q0;       
        mag = abs(S11)';
        FWHM = zeros(size(mag,2),1);
        %find the indices to the minima (resonant frequency) of each complete frequency scan until the end of the data
        for jj = 1:length(sim.field) %Searching column minima (fixed field)
            [~,idx] = min( mag(:,jj) );
            HM = ( max(mag(:,jj)) + min(mag(:,jj)) )/2; %
            % Calculate quality factor using f0/FWHM
            left = mag(1:idx,jj);
            right = mag(idx:end,jj);
            [lidx,~] = find(left >= HM, 1, 'last');
            [ridx,~] = find(right >= HM, 1, 'first');
            if lidx == idx+ridx-1
                FWHM(jj) = Inf;
            else
                FWHM(jj) = sim.freq(idx+ridx-1)-sim.freq(lidx);
            end
        end
        sim.FWHM = FWHM;
        if Options.plot == true
            figure
            plot(field(:,1),1./Q0,'-')
            xlabel('Magnetic field (T)')
            ylabel('1/Q')
            figure
            plot(field(:,1),FWHM,'-')
            xlabel('Magnetic field (T)')
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
    filename = strcat('S11_LiHoF4_',sprintf('%1$3.3fK_%2$.2fDeg_%3$.1fDeg_%4$.2e',sim.temp,theta,phi,gama),...
        sprintf('alpha=%1$d-%2$d.gif',min(alpha)*180/pi,max(alpha)*180/pi)); % Varying phase angle
%     filename = strcat('S11_LiHoF4_',sprintf('%1$3.3fK_%2$.1fDeg_%3$.1fDeg_%4$.2e',sim.temp,theta,phi,gama),...
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

clearvars -except Options saveloc sim.temp phi sim
if Options.savedata == true
    savename = sprintf('sim_S11_%1$.3fK_%2$udeg.mat',sim.temp,phi);
    saveObj = fullfile(saveloc,savename);
    save(saveObj,'-v7.3');
end
end