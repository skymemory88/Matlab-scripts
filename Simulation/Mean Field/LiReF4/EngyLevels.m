function [results, continu_var] = EngyLevels(varargin)
clearvars -except varargin

% default options
Options.ion = 'Ho'; % support 'Er' and 'Ho'
Options.hyp = 0.00; % Hyperfine isotope proportion
Options.nZee = false; % nuclear Zeeman interaction
Options.Enorm = true; % Normalize the energy states to the ground state
Options.Elevel = false; % eigen-energies
Options.Ediff = false; % excitation spectrum
Options.deltaI2 = false; % Transitions between the next nearest neightbouring levels
Options.Js = true; % spin expectation values
Options.Mx = false; % off diagonal element of spin matrix
Options.RPA = false; % not implemented yet
    omega = linspace(0,35,3001); % frequency range (GHz) for RPA calculation
Options.popul = false; % thermal populations of each states
Options.Espan = false; % energy span of magnon modes
    f_cav = 3.645; % Resonant frequency of the cavity (GHz)
%     filFctr = 0.01005; % filling factor
    filFctr = 0.01005*1.04; % filling factor
Options.eUnit = 'meV'; % unit choice: 'eV', 'meV', 'GHz'
Options.ScanMode = 'field'; % temperature/field scan option
Options.plot = true; % plot the figure
Options.savedata = false;
Options.savegif = false;

if nargin == 4
    dscrt_var = varargin{1};
    theta = varargin{2};
    phi = varargin{3};
    N_level = varargin{4};
elseif nargin == 5
    dscrt_var = varargin{1};
    theta = varargin{2};
    phi = varargin{3};
    N_level = varargin{4};
    opt = varargin{5};
    custom = fieldnames(opt); % options to change
    change = intersect(fieldnames(Options), custom);
    for ii = 1:length(change)
        Options.(change{ii}) = opt.(change{ii});
%         Options = setfield(Options, change{ii}, opt.change(ii));
    end
    Options.plot = false;
elseif nargin == 6
    Options.ion = varargin{1};
    dscrt_var = varargin{2};
    theta = varargin{3};
    phi = varargin{4};
    N_level = varargin{5};
    Options.plot = varargin{6};
end

muB = 9.27401e-24; %[J/T]
muN = 5.05078e-27; % Nuclear magneton [J/T]
mu0 = 4*pi*1e-7; % [H/m]
kB = 8.61733e-2; % [meV/K]
hbar = 1.05457E-34; % Reduced Planck constant [J.s]
J2meV = 6.24151e+21; % [meV/J]
Gh2mV = hbar*2*pi*10^9*J2meV; % [meV/GHz]
elem_idx = find(strcmp(Options.ion,[{'Er'};{'Ho'};{'Yb'};{'Tm'};{'Gd'};{'Y'}]));

if ispc
    Options.location = ['G:\.shortcut-targets-by-id\1CapZB_um4grXCRbK6t_9FxyzYQn8ecQE\File sharing',...
        '\PhD program\Research projects\Li', Options.ion, 'F4 project\Data\Simulations\Matlab\Susceptibilities\'];
elseif ismac
    Options.location = ['/Users/yikaiyang/Library/CloudStorage/GoogleDrive-yikai.yang@epfl.ch/My Drive/'...
        'File sharing/PhD program/Research projects/LiHoF4 project/Data/Simulations/MATLAB/Susceptibilities/'];
end

if Options.nZee == true && Options.hyp ~= 0
    Options.location = [Options.location, 'Hz_I=1'];
else
    Options.location = [Options.location, 'Hz_I=0'];
end

lineSty = [":","-.","--","-"];
% color = ["black","red","blue","magenta","green","yellow","cyan"];
color = {[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.9290, 0.6940, 0.1250]; ...
    [0.4940, 0.1840, 0.5560]; [0.4660, 0.6740, 0.1880];[0.3010, 0.7450, 0.9330]};
if strcmp(Options.eUnit,'meV')
    eUnit = 1; % use default meV unit
    f2E = Gh2mV; % GHz to meV conversion
    ylab = 'Energy (meV)';
    ylab1 = 'Energy gap (meV)';
    ylab2 = 'Freqeuency range (meV)';
elseif strcmp(Options.eUnit, 'eV')
    eUnit = 1e-3; % meV to eV conversion
    f2E = Gh2mV*1e-3; % GHz to eV conversion
    ylab = 'Energy (eV)';
    ylab1 = 'Energy gap (eV)';
    ylab2 = 'Freqeuency range (eV)';
elseif strcmp(Options.eUnit,'GHz')
    eUnit = 1/Gh2mV; % meV to GHz conversion
    f2E = 1;
    ylab = 'Energy (GHz)';
    ylab1 = 'Energy gap (GHz)';
    ylab2 = 'Freqeuency range (GHz)';
end

lg = strings(1,numel(theta),numel(phi)); % Create an empty array for legends
figs_E = gobjects(N_level,numel(theta),numel(phi));
figs_dE = gobjects(N_level-1,numel(theta),numel(phi));
if Options.deltaI2; figs_dE2 = gobjects(N_level-2,numel(theta),numel(phi)); end
if Options.savegif == true
    im_t = numel(dscrt_var)*numel(theta)*numel(phi);
    im_idx = 1;
    if Options.Elevel == true
        Elevel_frame = []*im_t;
    end
    if Options.Ediff == true
        Ediff_frame = []*im_t;
    end
end
results.En = {};
results.Ediff = {};
results.Js = {};
results.Jmx = {};
results.JImx = {};
for iter = 1:numel(dscrt_var)
    for iter2 = 1:numel(theta)
        for iter3 = 1:numel(phi)
            switch Options.ScanMode % 1. Field plot. 2. Temperature plot
                case 'field'
                    lname = strcat(['Hscan_Li', Options.ion, 'F4_'], sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_hp=%4$.2f',...
                        dscrt_var(iter), theta(iter2), phi(iter3), Options.hyp),'.mat');
                    file = fullfile(Options.location,lname);
                    load(file,'-mat','vvv','ttt','eee','fff','ion');
                    continu_var = vecnorm(fff,2,1); % choose the continuous variable to be field
                    temperature = repmat(ttt, length(continu_var), 1);
                    xlab = 'Magnetic Field (T)';
                case 'temp'
                    lname = strcat(['Tscan_Li', Options.ion, 'F4_'], sprintf('%1$3.3fT_%2$.2fDg_%3$.1fDg_hp=%4$.2f',...
                        dscrt_var(iter), theta(iter2), phi(iter3), Options.hyp),'.mat');
                    file = fullfile(Options.location,lname);
                    load(file,'-mat','vvv','ttt','eee','fff','ion');
                    continu_var = ttt; % choose the continuous variable to be temperature
                    temperature = ttt;
                    xlab = 'Temperature (K)';
            end

            if Options.Enorm == true; eee = eee - min(eee,[],2); end % Normalize the energy amplitude to the lowest eigen-energy
%             if Options.Enorm == true
%                 for ii = 1:N_level
%                     eee(:,ii) = eee(:,ii) - E0/E2f; % Normalize the energy amplitude to reference point E0 (manual)
%                 end
%             end
            eigenE = squeeze(eee);
            eigenW = vvv;
            results.En{iter,iter2,iter3} = eigenE(:,1:N_level); % Energy eigenstate
            results.wav = vvv; % needs to implementing truncation by N_level
            Ediff = double.empty(0,size(eigenE,1));
            if Options.deltaI2; Ediff2 = double.empty(0,size(eigenE,1)); end
%% Plot eigen-energies
            if Options.Elevel == true && Options.plot == true
                fig_E = figure;
                figs_E(:, iter2, iter3) = plot(continu_var, eigenE(:,1:N_level) .* eUnit, 'Color', color{1}, 'linewidth', 2);
%                 figs_E(:, iter2, iter3) = plot(fields, eigenE(:,1:N_level), 'Color', 'b', 'linewidth', 2);
%                 figs_E(:, iter2, iter3) = plot(fields, eigenE(:,1:N_level), 'Color', [35 107 142]/255, 'linewidth', 2);
                hold on
                set(gca,'fontsize',14);
%                text(8,(j-1)*D+0.02,'0.1');
%                set(fig1,'position',[100 100 600 300])
                xlabel(xlab,'FontSize',14)
                ylabel(ylab,'FontSize',14)
                grid off
                box on;
                tit1='Energy levels';
                title(tit1,'FontSize',14)
                legend(num2str(dscrt_var(iter)*1000,'T = %u mK,  A = A_{th}'))
                if Options.savegif == true
                    drawnow
                    frame = getframe(fig_E);
                    Elevel_frame{im_idx} = frame2im(frame);
                end
            end
%% Plot excitation spectrum among neighbour levels
            if Options.Ediff == true
                for i = 1:N_level-2 % Up until the second from top level
                    Ediff(i,:) = eigenE(:,i+1)-eigenE(:,i); % Transition between the nearest neighbouring levels
                    if Options.deltaI2
                        Ediff2(i,:) = eigenE(:,i+2)-eigenE(:,i); % Transition between the next nearest neighbouring levels
                    end
                end
                Ediff(N_level-1,:) = eigenE(:,N_level)-eigenE(:,N_level-1); % Transition between the top two levels
                Ediff = Ediff * eUnit; % convert to set energy unit
                results.Ediff{iter,iter2,iter3} = Ediff;
                if exist('Ediff2','var')
                    Ediff2 = Ediff2 * eUnit; % convert to set energy unit
                    results.Ediff2{iter,iter2,iter3} = Ediff2;
                end
                % frequency = { '1.682 GHz', '3.436 GHz', '3.924 GHz', '4.449 GHz', '5.604 GHz' }; % frequency reference lines

                % Plot transitions between the nearest levels

                if Options.plot == true
%                     figure;
                    hold on
                    figs_dE(1, iter2, iter3) = plot(continu_var, Ediff(1,:), 'Marker', 'none', 'LineStyle', lineSty(1),...
                        'Color', 'r', 'LineWidth', 2);
                    if N_level-1 > 1
                        figs_dE(2:end, iter2, iter3) = plot(continu_var, Ediff(2:end,:), 'Marker', 'none', 'LineStyle',...
                            lineSty(1), 'Color', color{iter3},'LineWidth', 2);
                        % figs_dE(2:end, iter2, iter3) = plot(continu_var, Ediff(2:end,:), 'Marker', 'none',...
                        %     'LineStyle', lineSty(1), 'Color', 'w', 'LineWidth', 2);
                    end

                    % Plot transitions between the next nearest levels
                    if Options.deltaI2
                        figs_dE2(1, iter2, iter3) = plot(continu_var, Ediff2(1,:), 'Marker', 'none', 'LineStyle', lineSty(1),...
                            'Color', color{iter3},'LineWidth',2);
                        if N_level-1 > 2
                            figs_dE2(2:end, iter2, iter3) = plot(continu_var, Ediff2(2:end,:), 'Marker', 'none', 'LineStyle', lineSty(1),...
                                'Color', color{iter3},'LineWidth',2);
                        end

%                         hpcfr1 = plot(fields, Ediff,'Color','blue','linewidth',2);
% %                         hpcfr1 = plot(fields, Ediff(2:end,:),'Color','black','linewidth',2);
% %                         hpcfr1 = plot(fields, Ediff,'Color',[35 107 142]/255,'linewidth',2,'LineStyle','-');
%                         legend(num2str(dscrt_var(iter),'T = %.2f K, A = A_{th}'),num2str(dscrt_var(iter),'T = %.2f K, A = A_{th}'))

%                         % Plot the resonant frequency of a bare cavity
%                         plot([0 max(fields)],[f_cav f_cav],'-r','LineWidth',1.5);
%                         yyaxis(gca,'left');
%                         xticks(linspace(0.1,max(fields),10))
%                         text(8,(j-1)*D+0.02,'0.1');
                    end
                    set(gca,'FontSize',14,'Xtick',0:1:max(continu_var));
                    set(gca,'XTickLabelRotation',0)
                    currYlim = get(gca,'ylim');
                    xlabel(xlab,'FontSize',15);
                    ylabel(ylab1,'FontSize',15);
                    grid off
                    box on;
                    xlim([min(continu_var) max(continu_var)]);
                    xticks('auto')
                    ylim(currYlim);
%                     ylim([0.9*min(Ediff,[],'All') 1.1*max(Ediff,[],'All')]);
%                     ylim([1 5]);
%                     ylim([3.59 3.71]) % adjust y-axis range to fit experimental data
                    tit4 = 'Energy gaps at ';
                    title([tit4 sprintf('T = %.3f K',dscrt_var(iter))],'FontSize',15)
                    if Options.savegif == true
                        drawnow
                        frame = getframe(fig_dE);
                        Ediff_frame{im_idx} = frame2im(frame);
                    end
%                     lgd = legend(squeeze(figs_dE(1,:,:)), lg);
%                     lgd.FontSize = 12;
                end
            end
%% Add frequency span of the transitions at each field
            if Options.Espan == true
                hold on
                eRange = double.empty(0,length(Ediff));
                for i=1:size(Ediff,2)
                    eRange(i) = (max(Ediff(:,i))-min(Ediff(:,i)));
                end
                yyaxis right;
                plot(continu_var, eRange,'k-');
                plot(continu_var(eRange == min(eRange)),min(eRange),'o','MarkerSize',4,'MarkerFaceColor','none');
                set(gca,'fontsize',15);
                xlabel(xlab,'FontSize',15)
                ylabel(ylab2,'FontSize',15)
                grid off
                box on;
                xlim([min(continu_var) max(continu_var)]);
                ylim([0.9*min(eRange) 1.2*max(eRange)]);
                tit4='Frequency range the energy levels cover';
                title(tit4,'FontSize',15)
            end
            lg(iter, iter2, iter3) = ['\phi', num2str(phi(iter3)," = %.1f"), char(176)];
            if Options.savegif == true
                im_idx = im_idx + 1;
            end
%% Plot expectation values of spin moment
            if Options.Js == true
                Jx = cell(1,4); % create one cell for each of the four ions in the same unit cell
                Jy = cell(1,4);
                Jz = cell(1,4);
                Jnorm = cell(1,4);
                
                for ii = 1:4
                    Jx{ii} = double.empty(length(continu_var),0);
                    Jy{ii} = double.empty(length(continu_var),0);
                    Jz{ii} = double.empty(length(continu_var),0);
                    for jj = 1:length(continu_var)
                        if Options.hyp ~= 0
                        Jx{ii}(jj,1) = ion.Js_hyp(ii,1,jj,elem_idx);
                        Jy{ii}(jj,1) = ion.Js_hyp(ii,2,jj,elem_idx);
                        Jz{ii}(jj,1) = ion.Js_hyp(ii,3,jj,elem_idx);
                        else
                        Jx{ii}(jj,1) = ion.Js(ii,1,jj,elem_idx);
                        Jy{ii}(jj,1) = ion.Js(ii,2,jj,elem_idx);
                        Jz{ii}(jj,1) = ion.Js(ii,3,jj,elem_idx);
                        end
                    end
                    Jnorm{ii} = vecnorm([Jx{ii} Jy{ii} Jz{ii}], 2, 2);
                end
                                
                Jmx = double.empty(4,length(continu_var),0);
                Jmy = double.empty(4,length(continu_var),0);
                Jmz = double.empty(4,length(continu_var),0);
                for ii = 1:length(fff(1,:))
                    if Options.hyp ~= 0
                        Jmx(:,ii,1) = ion.Jmom_hyp(1,ii,elem_idx);
                        Jmy(:,ii,1) = ion.Jmom_hyp(2,ii,elem_idx);
                        Jmz(:,ii,1) = ion.Jmom_hyp(3,ii,elem_idx);
                    else
                        Jmx(:,ii,1) = ion.Jmom(1,ii,elem_idx);
                        Jmy(:,ii,1) = ion.Jmom(2,ii,elem_idx);
                        Jmz(:,ii,1) = ion.Jmom(3,ii,elem_idx);
                    end
                end
                results.Js{iter,iter2,iter3} = cat(3,Jmx,Jmy,Jmz);
                
                Jx_av = rms(Jmx); % RMS over the unit cell
                Jy_av = rms(Jmy);
                Jz_av = rms(Jmz);
                
                % Jx_av = mean(Jmx); % average over the unit cell
                % Jy_av = mean(Jmy);
                % Jz_av = mean(Jmz);

                if Options.plot == true
                    % figure;
                    % box on
                    % hold on
                    % plot(fields, mean(cell2mat(Jnorm), 2), 'LineStyle', '-', 'LineWidth', 2); % plot <J> averaged over the four ions in the unit cell
                    % legend(strcat(Options.ion,'-1'),strcat(Options.ion,'-2'),strcat(Options.ion,'-3'),strcat(Options.ion,'-4'))
                    % xlabel(xlab)
                    % ylabel('<J>')

                    % figure
                    box on
                    hold on
                    % plot(continu_var(1:30:end), Jx_av(1:30:end), '-o', 'LineWidth', 1, 'Color', 'b')
                    % plot(continu_var(1:30:end), Jy_av(1:30:end), '-s', 'LineWidth', 1, 'Color', 'b')
                    % plot(continu_var(1:15:end), Jz_av(1:15:end), '-^', 'LineWidth', 1, 'Color', 'k')
                    plot(continu_var, Jx_av, 'LineStyle', '-', 'LineWidth', 2, 'Color', color{1})
                    plot(continu_var, Jy_av, 'LineStyle', '-', 'LineWidth', 2, 'Color', color{2})
                    plot(continu_var, Jz_av, 'LineStyle', '-', 'LineWidth', 2, 'Color', color{3})
                    legend('\langle J_x \rangle', '\langle J_y \rangle', '\langle J_z \rangle')
                    xlabel(xlab)
                    ylabel("$\langle J_\alpha \rangle$", 'Interpreter', 'latex')
                end
            end
%% Plot the matrix elements of the CMP coupling strength
            if Options.Mx == true
%                         Er      Ho      Yb      Tm      Gd      Y
                J_tab = [15/2;    8;      7/2;    6;      7/2;    1];
                L_tab = [6;       6;       3;     5;      0;      1];
                S_tab = [3/2;     2;      1/2;    1;      7/2;    1];
                I_tab = [3;      3.5;      0;     0;      1.50;   0];
                nLande = [-0.1611; 1.192; -0.2592; -0.462; -0.2265; 0]; % https://easyspin.org/documentation/isotopetable.html

                % Ions' lattice parameters
                lattice = [{[5.162 0 0; 0 5.162 0; 0 0 10.70]}      %Er
                           {[5.175 0 0; 0 5.175 0; 0 0 10.75]}      %Ho
                           {[5.132 0 0; 0 5.132 0; 0 0 10.59]}      %Yb
                           {[5.150 0 0; 0 5.150 0; 0 0 10.64]}      %Tm
                           {[5.162 0 0; 0 5.162 0; 0 0 10.70]}      %Gd
%                            {[5.132 0 0; 0 5.132 0; 0 0 10.59]}     %Gd
                           {[5.132 0 0; 0 5.132 0; 0 0 10.59]}];    %Y
                
                lattice = lattice{elem_idx};
                ionJ = J_tab(elem_idx); % Electronic moment for Ho3+
                ionI = I_tab(elem_idx); % Nuclear moment for Ho3+

                ELEf = gLande(L_tab(elem_idx),S_tab(elem_idx)) * muB; % Lande factor * Bohr magneton (J/T)
                NUCf = nLande(elem_idx) * muN; % (J/T)
                rho = 4e30 / det(lattice); % magnetic ion number density [m^-3]
                gw0 = sqrt(mu0 * 2*pi * f_cav*1e9 * rho/2) * filFctr; % susceptibility prefactor [T^2/J. rad/s]^1/2
                gw2 = gw0^2 * 2*pi * 1e-9; % [T^2/J. GHz]

                % Initiate ionJ operators
                Jz = diag(ionJ:-1:-ionJ); % Jz = -J, -J+1,...,J-1,J
                JhT.z = kron(Jz,eye(2*ionI+1)); % Expand Jz space to include nuclear degree of freedom
                Jp = diag(sqrt((ionJ-((ionJ-1):-1:-ionJ) ).*(ionJ+1+( (ionJ-1):-1:-ionJ) )),1); % electronic spin ladder operator
                Jm = Jp'; % electronic spin ladder operator
                Jph = kron(Jp,eye(2*ionI+1)); % Expand to match the dimension of Hilbert space
                Jmh = kron(Jm,eye(2*ionI+1));
                JhT.x = (Jph+Jmh)/2;
                JhT.y  =(Jph-Jmh)/2i;
                
                % Initiate I operators
                Iz = diag(ionI:-1:-ionI); %Iz = -I, -I+1,...,I-1,I
                IhT.z = kron(eye(2*ionJ+1),Iz); % Expand Hilbert space
                Ip = diag(sqrt((ionI-((ionI-1):-1:-ionI)).*(ionI+1+((ionI-1):-1:-ionI))),1); % Nuclear spin ladder operator
                Im = Ip'; % Nuclear spin ladder operator
                Iph = kron(eye(2*ionJ+1),Ip); % Expand to match the dimension of Hilbert space
                Imh = kron(eye(2*ionJ+1),Im);
                IhT.x = (Iph+Imh)/2;
                IhT.y = (Iph-Imh)/2i;

                tx = double.empty(size(eigenW,2),size(eigenW,2),size(continu_var,2),0); % Expectation value of Jx pseudo-spin
                ttx = double.empty(size(eigenW,2),size(eigenW,2),size(continu_var,2),0); % Expectation value of Jx-Ix pseudo-spin
                ty = double.empty(size(eigenW,2),size(eigenW,2),size(continu_var,2),0); % Expectation value of Jy pseudo-spin
                tty = double.empty(size(eigenW,2),size(eigenW,2),size(continu_var,2),0); % Expectation value of Jy-Iy pseudo-spin
                tz = double.empty(size(eigenW,2),size(eigenW,2),size(continu_var,2),0); % Expectation value of Jz pseudo-spin
                ttz = double.empty(size(eigenW,2),size(eigenW,2),size(continu_var,2),0); % Expectation value of Jz-Iz pseudo-spin
                Gc2 = double.empty(size(eigenW,2),size(eigenW,2),size(continu_var,2),0); % Coupling strength
                dP = double.empty(size(eigenW,2),size(eigenW,2),size(continu_var,2),0); % Coupling strength
                for kk = 1:size(continu_var,2) % calculate susceptibility for all fields
                    en = squeeze(eigenE(kk,:)); % Obtain the corresponding eigen energies [meV]
                    if temperature(kk) ~= 0
                        beta = 1/(kB*temperature(kk)); % [1/meV]
                        Z = sum(exp(-beta*en));
                        zn = exp(-beta*en)/Z;
                        [n,np] = meshgrid(zn,zn);
                        NN = n-np;
                    else
                        Z = zeros(size(en));
                        Z(1) = 1;
                        [n,np] = meshgrid(Z,Z);
                        NN = n-np;
                    end

                    v = squeeze(eigenW(kk,:,:)); % Obtain the corresponding eigen vectors
                    JxhT = JhT.x * ELEf;
                    IxhT = IhT.x * NUCf;
                    tx(:,:,kk,1) = v' * JxhT * v;
                    ttx(:,:,kk,1)  = v' * (JxhT+IxhT) * v;
                    
                    JyhT = JhT.y * ELEf;
                    IyhT = IhT.y * NUCf;
                    ty(:,:,kk,1) = v' * JyhT * v;
                    tty(:,:,kk,1)  = v' * (JyhT+IyhT) * v;
                    
                    JzhT = JhT.z * ELEf;
                    IzhT = IhT.z * NUCf;
                    tz(:,:,kk,1) = v' * JzhT * v;
                    ttz(:,:,kk,1)  = v' * (JzhT+IzhT) * v;
                    
                    Gc2(:,:,kk,1) = gw2 * ttz(:,:,kk,1) .* ttz(:,:,kk,1).' .* NN * J2meV / Gh2mV; % [meV^2]
                    dP(:,:,kk,1) = NN;
                end
                results.Jmx{iter,iter2,iter3} = tx;
                results.Jmy{iter,iter2,iter3} = ty;
                results.Jmz{iter,iter2,iter3} = tz;
                
                results.JIx{iter,iter2,iter3} = ttx;
                results.JIy{iter,iter2,iter3} = tty;
                results.JIz{iter,iter2,iter3} = ttz;
                
                results.pop{iter,iter2,iter3} = dP;
                results.Gc2{iter,iter2,iter3} = Gc2;

                if Options.plot == true
                    figure;
                    hold on
                    for ii = 1:N_level-1
                        plot(continu_var, sqrt(squeeze(abs(Gc2(ii,ii+1,:)))),'.');
%                         plot(fields, sqrt(abs(squeeze(Gc2(ii+1,ii,:)))),'o'); % symmetrical part of the matrix
                    end
                    xlabel(xlab)
                    ylabel(ylab)
                    title('Nearest level excitation coupling strength')
                end
            end
%% Plot the RPA excitation modes at T=0 (only works for the case of exchange anisotropy)
            if Options.RPA == true
                ionJ = ion.J(elem_idx); % Electronic moment for Ho3+
                ionI = ion.I(elem_idx);
                Jz = diag(ionJ:-1:-ionJ); % Jz = -J, -J+1,...,J-1,J
                if ionI ~= 0
                    Iz = diag(ionI:-1:-ionI); % Iz = -I, -I+1,...,I-1,I
                    Jz = kron(Jz, eye(2*ionI+1));
                    Iz = kron(eye(2*ionJ+1), Iz); % Expand Hilbert space
                end
                
                ELEf = ion.gLande(elem_idx) * muB * J2meV; % Lande factor * Bohr magneton [meV/T]
                NUCf = ion.nLande(elem_idx) * muN * J2meV; % nuclear Lande * magneton [meV/T]
                gfac = mu0/4/pi * (ion.gLande(elem_idx) * muB)^2 * J2meV * 10^30; % mu0/(4pi).(gL.muB)^2 [meV.Ang^3]
                demagn = zeros(3,3,4,4); % demagnetization field
                demagn_t = ellipsoid_demagn(ion.alpha);
                demagn(1,1,:,:) = demagn_t(1,1);
                demagn(2,2,:,:) = demagn_t(2,2);
                demagn(3,3,:,:) = demagn_t(3,3);

                lattice = ion.abc{elem_idx}; % retrieve lattice constants
                Vc = sum( lattice(1,:) .* cross(lattice(2,:), lattice(3,:)) ); % Volume of unit cell [Ang^-3]
                eins = zeros(3,3,4,4);
                eins(1,1,:,:) = 1; eins(2,2,:,:) = 1; eins(3,3,:,:) = 1;
                Dip = gfac * (MF_dipole([0 0 0], ion.dpRng, lattice, ion.tau) + 4*pi/Vc*(eins/3 - ion.demag*demagn)); % [meV] dipole + Lorenz term
                Jex = exchange([0 0 0], abs(ion.ex(elem_idx)), lattice, ion.tau); % [meV] AFM exchange interaction
                
                Dip = sum(sum(Dip,4),3)/size(ion.tau,1); % average over the unit cell
                Jex = sum(sum(Jex,4),3)/size(ion.tau,1);
                Jij = diag(ion.renorm(elem_idx,:)) .* (Dip - Jex) .* eUnit;

                dE = double.empty(length(continu_var), N_level-1, 0);
                deno = double.empty(length(continu_var), length(omega), 0);
                poles = zeros(length(continu_var), N_level-1);
                for kk = 1:length(continu_var)
                    en = squeeze(eigenE(kk,:));
%                     if temperature(kk) > 0            
%                         beta = 1/(kB*temperature(kk)); % [meV^-1]
%                         Z = sum(exp(-beta*en));
%                         zn = exp(-beta*en)/Z;
%                     else
%                         zn = ones(1,N_level);
%                         zn(1) = 0;
%                     end
                    vv = squeeze(eigenW(kk,:,:));
                    jzz = vv' * (Jz + NUCf/ELEf .* Iz) * vv;
                    jmn = jzz(2:N_level,1); % off diagonal elements between <n| and |0>
%                     jmn = (zn(2:N_level)-zn(1))' .* jzz(2:N_level,1); % off diagonal elements between <n| and |0>
                    J0 = [0 0 1] * Jij * [0 0 1]'; % Ising element of q=0 interaction tensor

                    dE0 = en - min(en,[],2); % MF excitation modes at 0 K
                    dE0 = dE0(2:N_level)' .* eUnit;
                    omega = omega .* f2E;
                    parfor nw = 1:length(omega) % search for RPA poles
                        dE2 = (dE0.^2 - omega(nw)^2);
                        deno(kk,nw,1) = prod(dE2)*(1 - J0 * (2*dE0') * (abs(jmn).^2 ./ dE2));
                    end
                    pks = squeeze(mag2db(1./abs(deno(kk,:,1)))); % search for poles
                    [~,loc] = findpeaks(pks, 'SortStr', 'none');
                    poles(kk,1:length(loc)) = flip(omega(loc));
                    for ii = 2:size(poles,2)
                        dE(:,ii-1,1) = poles(:,ii-1)-poles(:,ii);
                    end
                    results.freq = omega;
                    results.poles = poles;
                    results.deno = squeeze(deno);
                    results.dE = dE;
                end

                if Options.plot == true
                    figure;
                    plot(continu_var, poles);
                    xlabel('Magnetic Field (T)')
                    ylabel(ylab1)
                end
            end
%% Plot the occupation number of each state
            if Options.popul == true
                En = eee - min(eee,[],2);
                beta = 11.6/dscrt_var(iter); % 1/(kB T) [meV]^-1
                Pop = cell(N_level,1); % thermal population
                Z = double.empty(size(En,1),0);
                for ii = 1:size(En,1)
                    Z(ii,1) = sum(exp(-squeeze(En(ii,:))*beta)); % partition function at each field
                end
                for jj = 1:N_level
                    Pop{jj} = exp(-En(:,jj)*beta) ./ Z; % Thermal population of each level
                end
                results.popul = Pop;

                if Options.plot == true
                    pFig = figure;
                    hold on
                    box on
                    pRat = figure;
                    hold on
                    box on
                    if jj > 1
                        for jj = 2:N_level
                            figure(pFig)
                            plot(continu_var, Pop{jj}(:)./Pop{1}(:),'-','LineWidth',2)
                            figure(pRat)
                            plot(continu_var, Pop{jj}(:)./Pop{jj-1}(:),'-','LineWidth',2);
                        end
                    end
                    figure(pFig)
                    xlim([0 max(continu_var)])
                    xticks('auto')
                    set(gca,'FontSize', 14)
                    xlabel(xlab);
                    ylabel('P_{|i\rangle}/P_{|0\rangle}')
                    legend(["P_{|1\rangle}/P_{|0\rangle}", "P_{|2\rangle}/P_{|0\rangle}", "P_{|3\rangle}/P_{|0\rangle}",...
                    	"P_{|4\rangle}/P_{|0\rangle}", "P_{|5\rangle}/P_{|0\rangle}", "P_{|6\rangle}/P_{|0\rangle}",...
                    	"P_{|7\rangle}/P_{|0\rangle}"]);

                    figure(pRat)
                    xlim([0 max(continu_var)])
                    xticks('auto')
                    set(gca,'FontSize', 14)
                    xlabel(xlab);
                    ylabel('P_{|i+1\rangle}/P_{|i\rangle}')
                    legend(["P_{|1\rangle}/P_{|0\rangle}", "P_{|2\rangle}/P_{|1\rangle}", "P_{|3\rangle}/P_{|2\rangle}",...
                    	"P_{|4\rangle}/P_{|3\rangle}", "P_{|5\rangle}/P_{|4\rangle}", "P_{|6\rangle}/P_{|5\rangle}",...
                    	"P_{|7\rangle}/P_{|6\rangle}"]);
                end
            end
        end
    end
end

if Options.savedata == true
    lname = strcat('sim_',num2str(dscrt_var*1000,'%u'),'mK_trans.mat','-v7.3');
    save(lname,'continu_var','Ediff');
end
if Options.savegif == true && Options.plot == true
    if Options.Elevel == true
        lname = strcat('sim_',num2str(min(dscrt_var)*1000,'%u'),'-',num2str(max(dscrt_var)*1000,'%u'),'mK_Elevel.gif');
        fileobj = fullfile(Options.location,lname);
        for ii = 1:length(Elevel_frame)
            [img,cmp] = rgb2ind(Elevel_frame{ii},256);
            if ii == 1
                imwrite(img,cmp,fileobj,'gif','LoopCount',inf,'DelayTime',1);
            else
                imwrite(img,cmp,fileobj,'gif','WriteMode','append','DelayTime',1);
            end
        end
    end
    if Options.Ediff == true
        lname = strcat('sim_',num2str(min(dscrt_var)*1000,'%u'),'-',num2str(max(dscrt_var)*1000,'%u'),'mK_Ediff.gif');
        fileobj = fullfile(Options.location,lname);
        for ii = 1:length(Ediff_frame)
            [img,cmp] = rgb2ind(Ediff_frame{ii},256);
            if ii == 1
                imwrite(img,cmp,fileobj,'gif','LoopCount',inf,'DelayTime',1);
            else
                imwrite(img,cmp,fileobj,'gif','WriteMode','append','DelayTime',1);
            end
        end
    end
end
end
