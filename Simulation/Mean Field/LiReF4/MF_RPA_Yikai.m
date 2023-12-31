function  MF_RPA_Yikai(mion, scanMode, dscrt_var, freq_total, theta, phi, gama, hyp, RPA_mode)
% Current version assumes complete symmetrical equivalence among the four spin moments per unit cell
% mion: Magnetic ion type: 'Er', 'Ho'
% scanMode: 'field' or 'temp' scan
% dscrt_var: discrete variable (temperature or field)
% freq_total frequency range [GHz]
% theta: misalignment angle from c-axis [degree]
% phi: inplane angle in ab-plane between the magnetic field and a/b-axis [degree]
% gama: linewidth of hyperfine levels [meV].
% hyp: Hyperfine isotope proportion (0~1)
% RPA_Mode: RPA option (true/false)

Options.plotting = false; % Decide whether or not to plot the data at thes end
Options.unit = 'meV'; % Energy unit choice: J, GHz, meV (default)
Options.saving = true; % Options to save the susceptibility tensors
Options.scanMode = scanMode; % 1. Field plot with RPA; 2. Temp plot with RPA; 3. wavevector plot with RPA
Options.nZee = true;
Options.RPA = RPA_mode; % Apply random phase approximation (RPA) correction
Options.Kplot = false; % k-dependent plot for RPA susceptibilities

if Options.RPA == false
    Options.Kplot = false;
end
if Options.Kplot == true
    % qz = (6.9:0.0025:7.5)';
    qx = linspace(0,2,301)';
    qy = zeros(length(qx),1);
    qz = zeros(length(qx),1);
    qvec = [qx qy qz];
    cVar0 = [4 5]; % selected points of the continuous variable for k-plot
else
    qx = 0.0;
    %         qx = [0.01 0.1 0.3 0.6 1]';
    qy = zeros(size(qx,1),1);
    qz = zeros(size(qx,1),1);
    qvec = [qx qy qz];
end
clearvars -except dscrt_var theta phi gama Options mion freq_total hyp qvec cVar0 nZee_path

if strcmp(pathsep, ':')
    platform = 'Unix';
else
    platform = 'Win';
end
% Declare physical constants as global for consistency
const.hbar = 1.05457E-34; % Reduced Planck constant [J.s]
const.J2meV = 6.24151e+21; % convert Joule to meV
const.Gh2mV = const.hbar * 2*pi * 1e9 * const.J2meV; % convert GHz to meV
const.muB = 9.27401e-24; % Bohr magneton [J/T]
const.muN = 5.05078e-27; % Nuclear magneton [J/T]
const.mu0 = 4*pi*1e-7; % [H/m]
const.dpRng = 100; % dipole summation range (number of unit cell)

% Select the final output unit (default: meV)
switch Options.unit
    case 'GHz'
    ConvUnit = 1/const.Gh2mV;
    case 'J'
    ConvUnit = 1/const.J2meV;
    otherwise
    ConvUnit = 1;
end

for ii = 1:length(dscrt_var) 
    ipt = false;
    counter = 0;
    if Options.nZee == true
        nZee_path = 'Hz_I=1';
    else
        nZee_path = 'Hz_I=0';
    end
    switch platform
        case 'Win'
            Options.location = ['G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\',...
                'Simulations\MATLAB\Susceptibilities\', nZee_path, '\'];
        case 'Unix'
            Options.location = ['/Users/yikaiyang/Library/CloudStorage/GoogleDrive-yikai.yang@epfl.ch/My Drive/'...
                'File sharing/PhD program/Research projects/LiHoF4 project/Data/',...
                'Simulations/MATLAB/Susceptibilities/', nZee_path, '/'];
    end
    while ~ipt
        switch Options.scanMode % 1. Field plot with RPA. 2. wavevector plot with RPA
            case 'field'
                filename = strcat(['Hscan_Li',mion,'F4_'],...
                    sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_hp=%4$.2f', dscrt_var(ii), theta, phi, hyp),'.mat');
                file = fullfile(Options.location, filename);
                load(file,'-mat','eee','fff','ttt','vvv','ion'); % loads variables "fields", "temp", "E" (eigenvalues) and "V" (eigenstates)
                file_part = sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_%4$.2e_%5$.3f-%6$.3fGHz'...
                    , dscrt_var(ii), theta, phi, gama, min(freq_total), max(freq_total));
                cVar = vecnorm(fff,2,1); % choose the continuous variable to be field
                ttt = repelem(ttt, length(cVar));
                ipt = true;
                fprintf('Calculating for T = %.3f K.\n', dscrt_var(ii));
            case 'temp'
                filename = strcat(['Tscan_Li',mion,'F4_'],...
                    sprintf('%1$3.3fT_%2$.2fDg_%3$.1fDg_hp=%4$.2f', dscrt_var(ii), theta, phi, hyp),'.mat');
                file = fullfile(Options.location,filename);
                load(file,'-mat','eee','ttt','vvv','ion'); % loads variables "fields", "temp", "E" (eigenvalues) and "V" (eigenstates)
                file_part = sprintf('%1$3.3fT_%2$.2fDg_%3$.1fDg_%4$.2e_%5$.3f-%6$.3fGHz'...
                    , dscrt_var(ii), theta, phi, gama, min(freq_total), max(freq_total));
                cVar = ttt; % choose the continuous variable to be temperature
                ipt = true;
                fprintf('Calculating for B = %.3f T.\n', dscrt_var(ii));
            otherwise
                prompt = sprintf('Unknown Scan Mode! Please select a mode between tempscan and fieldscan.\n');
                answer = input(prompt,'s');
                Options.scanMode = lower(answer);
                counter = counter + 1; % input request counter
                if counter >= 3
                    return % terminate the program after three failed requests
                end
        end
    end
    const.elem = find(ion.prop);    
    const.ELEf = ion.gLande(const.elem) * const.muB * const.J2meV; % Lande factor * Bohr magneton (meV/T)
    const.NUCf = ion.nLande(const.elem) * const.muN * const.J2meV; % [meV/T]
    const.gfac = ion.gLande(const.elem)^2 * (const.muB)^2 * (const.mu0/4/pi) * const.J2meV * 10^30; % (gL.muB)^2.mu0/4pi [meV.Ang^3]

    eee = eee - min(eee,[],2); % Normalize the eigen-energies to the ground state
    save_name = strcat('chi_Li',mion,'F4_',file_part, '.mat');
    
    eigenE = eee; % load eigenenergy [meV]
%     eigenE(:,9:end) = 0; % truncate the Hilbert space
    eigenW = vvv; % load eigenstates
%     eigenW(:,9:end,9:end) = 0; % truncate the Hilbert space
    if Options.RPA == true
        if Options.Kplot == true
            bidx = int16.empty(0,length(cVar0));
            for jj = 1:length(cVar0)
                [~,bidx(jj)] = min(abs( cVar - cVar0(jj) ));
            end
            eigenW = vvv(bidx,:,:); % eigen-functions
            eigenE = eee(bidx,:); % eigen-energies [meV]
            cVar = cVar(bidx); % magnetic field or temperature
        end
        [cVar, freq_total, chi0, ~] = linear_response(ion, eigenE, cVar, freq_total, ttt, eigenW, gama, const, Options);
        [~, ~, ~, chiq, ~] = RPA(qvec, cVar, freq_total, ion, chi0, const); % Electronic susceptibilitie
%         chiq = const.ELEf^2 * chiq .* ConvUnit; % [J/T^2 or GHz/T^2 or meV/T^2]
        chiq = chiq .* ConvUnit; % [J/T^2 or GHz/T^2 or meV/T^2]
    else
        [cVar, freq_total, chi0, ~] = linear_response(ion, eigenE, cVar, freq_total, ttt, eigenW, gama, const, Options);
    end
%     chi0 = const.ELEf^2 * chi0 .* ConvUnit; % [J/T^2 or GHz/T^2 or meV/T^2]
    chi0 = chi0 .* ConvUnit; % [J/T^2 or GHz/T^2 or meV/T^2]

    if Options.saving == true % Save the susceptibilities
        if Options.RPA == true
            savefile1 = fullfile(Options.location,save_name);
            save_vars = {dscrt_var, cVar, freq_total, ion, chi0, gama, qvec, chiq};
            save_file(save_vars, savefile1, Options); % save the data w. RPA corrections
        else
            savefile2 = fullfile(Options.location,save_name);
            save_vars = {dscrt_var, cVar, freq_total, ion, chi0, gama};
            save_file(save_vars, savefile2, Options); % Save the data w/o RPA corrections
        end
    end
    
    if Options.plotting == true % Plot the susceptibilities
        plot_var = {cVar, freq_total, chi0, gama, [0,0,0], dscrt_var(ii)};
        if Options.Kplot == true
            plot_var{3} = chiq;
            plot_var{5} = qvec;
            figs(plot_var, Options, const, "\chi_{RPA}", Options.Kplot);
            chi0_p = repmat(chi0,1,1,1,1,length(qvec));
%             chi0_p = permute(chi0_p,[1 2 3 5 4]); % Add the additional dimension for the plotting purpose
            plot_var{3} = chi0_p;
            figs(plot_var, Options, const, "chi_0", Options.Kplot);
        else
            figs(plot_var, Options, const, "\chi_0", Options.Kplot);
            if Options.RPA == true
                plot_var{3} = chiq;
                figs(plot_var, Options, const, "\chi_{RPA}", Options.Kplot);
                plot_var{3} = chiq-chi0;
                figs(plot_var, Options, const, "\chi_{RPA}-\chi_0", Options.Kplot);
            end
        end
    end
end
clearvars chi0 chi_p chiq
end

function figs(input_var, Options, const, fig_tit, Qplot)
continu_var = input_var{1};
freq_total = input_var{2};
chi = input_var{3};
gama = input_var{4};
qvec = input_var{5};
dscrt_var = input_var{6};

% custom colormap scale
cmap = unique([[0 0 0];[zeros(20,1),linspace(0,0.5,20)',linspace(0.4,1,20)'];...
[linspace(0,1,36)',0.5*ones(36,1),linspace(1,0,36)'];...
[ones(30,1),linspace(0.5,1,30)',linspace(0,1,30)']],'rows');
cmap = flip(cmap,1);

% convert frequency axis to appropriate choice
if strcmp(Options.unit, 'J')
    freq_total = freq_total*const.Gh2mV/const.J2meV;
    ylab = "Energy (J)"; % yaxis label
elseif strcmp(Options.unit, 'meV')
    freq_total = freq_total*const.Gh2mV;
    ylab = "Energy (meV)";
else
    ylab = "Frequency (GHz)";
end
 
pos0 = [100 300 600 400]; % initial figure position
pos_inc = [150 0 0 0];
indx = [strcat(fig_tit,"^{xx}") strcat(fig_tit,"^{yy}") strcat(fig_tit,"^{zz}")]; % index for figure legend
if Qplot == true
    qv = qvec(:,1); % use qx
    for nb = 1:size(continu_var,1)
        switch Options.scanMode
            case 'field'
                file_part1 = sprintf('T = %.3f K. and ', dscrt_var);
                file_part2 = sprintf('B = %1$.2f T.', continu_var(nb));
            case 'temp'
                file_part1 = sprintf('B = %.3f T. and ', dscrt_var);
                file_part2 = sprintf('B = %$.2f T.', continu_var(nb));
            otherwise
                disp("Unknow scan mode!")
                return
        end
        for jj = 3 % plot only chi_zz
%         for jj = 1:3 % plot all three diagonal (x, y, z) elements
            % Color plot the imaginary part of the susceptibilities of z component
            fig(1) = figure;
            set(fig(1),'position',pos0 + (jj-1)*pos_inc);
            hp1 = pcolor(qv,freq_total,real(squeeze(chi(jj,jj,:,nb,:))));
            set(hp1, 'edgeColor','none')
            set(gca, 'xdir', 'reverse' )
            colormap(cmap)
%             caxis([0 100]);
            caxis('auto');
            colorbar
            legend(['\gamma =' num2str(gama,'%.2e meV')]);
            title(['Re[', char(indx(jj)), '] at ', file_part1, file_part2])
            xlabel(sprintf('Q = [h, 0, 0]'))
            ylabel(ylab)
      
            % Plot the real part of the susceptibility of the z component
            fig(2) = figure;
            set(fig(2),'position',pos0 + jj*pos_inc);
            hp2 = pcolor(qv,freq_total,imag(squeeze(chi(jj,jj,:,nb,:))));
            set(hp2, 'edgeColor','none')
            set(gca, 'xdir', 'reverse' )
            colormap(cmap)
%             caxis([0 100]);
            caxis('auto');
            colorbar
            legend(['\gamma =' num2str(gama,'%.2e meV')]);
            title(['Im[', char(indx(jj)), '] at ', file_part1, file_part2])
            xlabel(sprintf('Q = [h, 0, 0]'))
            ylabel('Frequency (GHz)')
            ylabel(ylab)
        end
    end
else
    for nq = 1:size(qvec,1)
        switch Options.scanMode
            case 'field'
                xlab = 'Magnetic Field (T)';
                file_part1 = sprintf('T = %.3f K and ', dscrt_var);
                file_part2 = sprintf('Q = [%1$.2f %2$.2f %3$.2f]', qvec(nq,1), qvec(nq,2), qvec(nq,3));
            case 'temp'
                xlab = 'Temperature (K)';
                file_part1 = sprintf('T = %.3f K and ', dscrt_var);
                file_part2 = sprintf('Q = [%1$.2f %2$.2f %3$.2f]', qvec(nq,1), qvec(nq,2), qvec(nq,3));
            otherwise
                disp("Unknow scan mode!")
                return
        end
        % Color plot the imaginary part of the susceptibilities
        for jj = 3 % plot only chi_zz component
%             for jj = 1:3
            fig(1) = figure;
            set(fig(1),'position',pos0 + (jj-1)*pos_inc);
%             hp2 = pcolor(fields(1,:),freq_total,squeeze(imag(chi(jj,jj,:,:,nq))));
            hp2 = pcolor(continu_var(1,:),freq_total,mag2db(abs(squeeze(imag(chi(jj,jj,:,:,nq))))));
            set(hp2, 'edgeColor','none')
            colormap(cmap)
%             caxis([0 5]);
            caxis('auto');
            colorbar
            legend(['\gamma =' num2str(gama,'%.2e meV')]);
            title(['Im[', char(indx(jj)), '] at ', file_part1, file_part2])
            xlabel(xlab)
            ylabel(ylab)
            
            % Plot the real part of the susceptibility of the z component
            fig(2) = figure;
            set(fig(2),'position',pos0 + jj*pos_inc);
%             hp3 = pcolor(fields(1,:),freq_total,squeeze(real(chi(jj,jj,:,:,nq))));
            hp3 = pcolor(continu_var(1,:),freq_total,mag2db(abs(squeeze(real(chi(jj,jj,:,:,nq))))));
            set(hp3, 'edgeColor','none')
            colormap(cmap)
%             caxis([0 5]);
            caxis('auto');
            colorbar
            legend(['\gamma =' num2str(gama,'%.2e meV')]);
            title(['Re[', char(indx(jj)), '] at ', file_part1, file_part2])
            xlabel(xlab)
            ylabel(ylab)
        end
    end
end
end

function save_file(save_vars, savefile, opt)
switch opt.scanMode
    case 'field'
        temp = save_vars{1};
        fields = save_vars{2};
    case 'temp'
        fields = save_vars{1};
        temp = save_vars{2};
end
freq_total = save_vars{3};
ion = save_vars{4};
chi = save_vars{5};
gama = save_vars{6};
unit = opt.unit;
if length(save_vars) == 8
    qvec = save_vars{7};
    chiq = save_vars{8};
    save(savefile,'temp','fields','freq_total','ion','chi','gama','qvec','chiq','unit','-v7.3');
else
    save(savefile,'temp','fields','freq_total','ion','chi','gama','unit','-v7.3');
end
end

function [cVar, freq_total, chi, JIz_exp] = linear_response(ion, eigenE, cVar, freq_total,...
    temperatures, eigenW, gma, const, Options)
kB = 8.61733e-2; % Boltzmann constant [meV/K]
gamma = ones(size(eigenE,2))*gma; % spin linewidth [meV]

% Declare susceptibility tensor
chi0 = zeros(3,3,length(freq_total(1,:)),size(cVar,2));

%Initiate ionJ operators
ionJ = ion.J(const.elem);
Jz = diag(ionJ:-1:-ionJ); % Jz = -J, -J+1, ... ,J-1,J
Jp = diag(sqrt((ionJ-((ionJ-1):-1:-ionJ) ).*(ionJ+1+( (ionJ-1):-1:-ionJ) )),1); % electronic spin ladder operator
% Jp = diag(sqrt(ionJ*(ionJ+1) - ((ionJ-1):-1:-ionJ).*(ionJ:-1:(-ionJ+1))),1); % alternative format of spin ladder operator
Jm = Jp'; % electronic spin ladder operator

if ion.hyp(const.elem) > 0 % hyperfine interaction option
    %Initiate I operators
    ionI = ion.I(const.elem);
    Iz = diag(ionI:-1:-ionI); % Iz = -I, -I+1,...,I-1,I
%     Iz = diag([ones(1,4)*0.5 ones(1,4)*(-0.5)]); % projection to Ising space
    IhT.z = kron(eye(2*ionJ+1), Iz); % Expand Hilbert space
    Ip = diag(sqrt((ionI-((ionI-1):-1:-ionI)).*(ionI+1+((ionI-1):-1:-ionI))),1); % Nuclear spin ladder operator
    % Ip = diag(sqrt(ionI*(ionI+1) - ((ionI-1):-1:-ionI).*(ionI:-1:(-ionI+1))),1); % alternative format of spin ladder operator
    Im = Ip'; % Nuclear spin ladder operator
    Iph = kron(eye(2*ionJ+1), Ip); % Expand to match the dimension of Hilbert space
    Imh = kron(eye(2*ionJ+1), Im);
    IhT.x = (Iph+Imh)/2;
    IhT.y = (Iph-Imh)/2i;
    
    Jph = kron(Jp, eye(2*ionI+1)); % Expand to match the dimension of Hilbert space
    Jmh = kron(Jm, eye(2*ionI+1));
    JhT.x = (Jph+Jmh)/2; % Expand J space to include nuclear degree of freedom
    JhT.y = (Jph-Jmh)/2i;
    JhT.z = kron(Jz, eye(2*ionI+1));
else
    JhT.x = (Jp+Jm)/2;
    JhT.y = (Jp-Jm)/2i;
    JhT.z = Jz;

    IhT.x = 0;
    IhT.y = 0;
    IhT.z = 0;
end
JhT.x = JhT.x + ion.hyp(const.elem) * const.NUCf/const.ELEf * IhT.x; % Hybridized electronuclear spin operator
JhT.y = JhT.y + ion.hyp(const.elem) * const.NUCf/const.ELEf * IhT.y;
JhT.z = JhT.z + ion.hyp(const.elem) * const.NUCf/const.ELEf * IhT.z;

% Single out <Jz+Iz> calculations
JIz_exp = double.empty(size(cVar,2), size(JhT.z,1),0); % Expectation value of J-I pseudo-spin
for ii = 1:size(cVar,2) % calculate susceptibility for all fields
    v = squeeze(eigenW(ii,:,:)); % Obtain the corresponding eigen vectors
    ttz  = v'  * JhT.z * v;
    JIz_exp(ii,:,1) = real(diag(ttz));
end

% If exists, load spin-linewidths extracted from experiments for field scan simulation
if temperatures(1) == temperatures(end) % for fixed temperature calculations
    if Options.nZee && ismember(temperatures(1), [0.15 0.18 0.22 0.30])
        gma_load = load([Options.location, sprintf('Exp_fit_gamma_%umK.mat', 1000*temperatures(1))],'gma');
        gma_load = flip(gma_load.gma)*const.Gh2mV; % [meV]
        for kk = 1:length(gma_load)
            gamma(kk,kk+1) = gma_load(kk);
            gamma(kk+1,kk) = gma_load(kk);
        end
    elseif ~Options.nZee && ismember(temperatures(1), [0.1 0.13 0.15 0.18 0.24 0.25])
        gma_load = load([Options.location, sprintf('Exp_fit_gamma_%umK.mat', 1000*temperatures(1))],'gma');
        gma_load = flip(gma_load.gma)*const.Gh2mV; % [meV]
        for kk = 1:length(gma_load)
            gamma(kk,kk+1) = gma_load(kk);
            gamma(kk+1,kk) = gma_load(kk);
        end
    end
end


for nf = 1:length(freq_total(1,:)) %calculate susceptibility for all frequencies
    freq = freq_total (nf);
    omega = freq*const.Gh2mV;   % GHz to meV
%     for jj = 1:size(continu_var,2) % for debugging: calculate susceptibility for all fields
    parfor jj = 1:size(cVar,2) % calculate susceptibility for all fields
        v = squeeze(eigenW(jj,:,:)); % Obtain the corresponding eigen vectors
        en = squeeze(eigenE(jj,:)); % Obtain the corresponding eigen energies [meV]
        if temperatures(jj) ~= 0
            beta = 1/(kB*temperatures(jj)); % [meV^-1]
            Z = sum(exp(-beta*en));
            zn = exp(-beta*en)/Z;
            [n,np] = meshgrid(zn,zn);
            NN = n-np;
        else
            Z = zeros(size(en));
            Z(1) = 1; % Occupy the ground state with unity probability
            [n,np] = meshgrid(Z,Z);
            NN = n-np;
        end
%         % pseudo-statistics of occupation number for debugging/testing
%         Z = 2.^(0:size(en,1)-1); % geometric sequence
%         Z = Z./sum(Z);
%         [n, np] = meshgrid(Z,Z);
%         NN = n - np;
        
        [ee,eep] = meshgrid(en,en);
        EE = (eep-ee-omega); % [meV]

        % computer separatly the electronic and nuclear spin susceptibilities 
        deno0 = 1 ./ (EE - 1i*gamma); % [meV^-1]
        chi0(:,:,nf,jj) = chi_Mx(JhT, JhT, v, NN, deno0); % Electornic spin operators
    end
end
chi = chi0; % Electronic susceptibilities
end

function chi0 = chi_Mx(varargin)
% Calculation of matrix elements for the real and imaginary susceptibilities
% opr1: operator 1
% opr2: operator 2
% wav: eigenfunction of the hamiltonian
% NN: population difference between two eigenstates
% deno1: denominator for real part of the susceptibility
% deno2: denominator for imaginary part of the susceptibility
opr1 = varargin{1};
opr2 = varargin{2};
wav = varargin{3};
NN = varargin{4};

ttx1  = wav'  * opr1.x * wav;
tty1  = wav'  * opr1.y * wav;
ttz1  = wav'  * opr1.z * wav;

ttx2  = wav'  * opr2.x * wav;
tty2  = wav'  * opr2.y * wav;
ttz2  = wav'  * opr2.z * wav;

if nargin == 6 % Compute separately the real and imaginary parts
    deno1 = varargin{5};
    deno2 = varargin{6};
    % Calculate susceptibilities along a-axis
    chi_tx  = (ttx1) .* (ttx2.') .* NN .* deno1; % Imaginary part
    chi_t1x = (ttx1) .* (ttx2.') .* NN .* deno2; % Real part
    xx = sum(sum(chi_tx));
    xx1 = sum(sum(chi_t1x));
    
    % Calculate susceptibilities along b-axis
    chi_ty  = (tty1) .* (tty2.') .* NN .* deno1;
    chi_t1y = (tty1) .* (tty2.') .* NN .* deno2;
    yy = sum(sum(chi_ty));
    yy1 = sum(sum(chi_t1y));
    
    % Calculate susceptibilities along c-axis
    chi_tz  = (ttz1) .* (ttz2.') .* NN .* deno1;
    chi_t1z = (ttz1) .* (ttz2.') .* NN .* deno2;
    zz = sum(sum(chi_tz));
    zz1 = sum(sum(chi_t1z));
    
    % Calculate susceptibilities along ab-axis
    chi_txy  = (ttx1) .* (tty2.') .* NN .* deno1;
    chi_t1xy = (ttx1) .* (tty2.') .* NN .* deno2;
    xy = sum(sum(chi_txy));
    xy1 = sum(sum(chi_t1xy));
    
    % Calculate susceptibilities along ac-axis
    chi_txz  = (ttx1) .* (ttz2.') .* NN .* deno1;
    chi_t1xz = (ttx1) .* (ttz2.') .* NN .* deno2;
    xz = sum(sum(chi_txz));
    xz1 = sum(sum(chi_t1xz));
    
    % Calculate susceptibilities along ba-axis
    chi_tyx  = (tty1) .* (ttx2.') .* NN .* deno1;
    chi_t1yx = (tty1) .* (ttx2.') .* NN .* deno2;
    yx = sum(sum(chi_tyx));
    yx1 = sum(sum(chi_t1yx));
    
    % Calculate susceptibilities along bc-axis
    chi_tyz  = (tty1) .* (ttz2.') .* NN .* deno1;
    chi_t1yz = (tty1) .* (ttz2.') .* NN .* deno2;
    yz = sum(sum(chi_tyz));
    yz1 = sum(sum(chi_t1yz));
    
    % Calculate susceptibilities along ca-axis
    chi_tzx  = (ttz1) .* (ttx2.') .* NN .* deno1;
    chi_t1zx = (ttz1) .* (ttx2.') .* NN .* deno2;
    zx = sum(sum(chi_tzx));
    zx1 = sum(sum(chi_t1zx));
    
    % Calculate susceptibilities along cb-axis
    chi_tzy  = (ttz1) .* (tty2.') .* NN .* deno1;
    chi_t1zy = (ttz1) .* (tty2.') .* NN .* deno2;
    zy = sum(sum(chi_tzy));
    zy1 = sum(sum(chi_t1zy));
    
    xr = [xx xy xz
          yx yy yz
          zx zy zz]; % Real part of susceptibility
    xi = [xx1 xy1 xz1
          yx1 yy1 yz1
          zx1 zy1 zz1]; % Imaginary part of susceptibility
    chi0 = xr + 1i.*xi;
      
elseif nargin == 5 % Compute the complex matrix element
    deno0 = varargin{5};
    % Calculate susceptibilities along a-axis
    chi_tx  = (ttx1) .* (ttx2.') .* NN .* deno0;
    xx = sum(sum(chi_tx));
    
    % Calculate susceptibilities along b-axis
    chi_ty  = (tty1) .* (tty2.') .* NN .* deno0;
    yy = sum(sum(chi_ty));
    
    % Calculate susceptibilities along c-axis
    chi_tz  = (ttz1) .* (ttz2.') .* NN .* deno0;
    zz = sum(sum(chi_tz));
    
    % Calculate susceptibilities along ab-axis
    chi_txy  = (ttx1) .* (tty2.') .* NN .* deno0;
    xy = sum(sum(chi_txy));
    
    % Calculate susceptibilities along ac-axis
    chi_txz  = (ttx1) .* (ttz2.') .* NN .* deno0;
    xz = sum(sum(chi_txz));
    
    % Calculate susceptibilities along ba-axis
    chi_tyx  = (tty1) .* (ttx2.') .* NN .* deno0;
    yx = sum(sum(chi_tyx));
    
    % Calculate susceptibilities along bc-axis
    chi_tyz  = (tty1) .* (ttz2.') .* NN .* deno0;
    yz = sum(sum(chi_tyz));
    
    % Calculate susceptibilities along ca-axis
    chi_tzx  = (ttz1) .* (ttx2.') .* NN .* deno0;
    zx = sum(sum(chi_tzx));
    
    % Calculate susceptibilities along cb-axis
    chi_tzy  = (ttz1) .* (tty2.') .* NN .* deno0;
    zy = sum(sum(chi_tzy));
       
    chi0 = [xx xy xz
            yx yy yz
            zx zy zz];
end
end

function [qvec, cVar, freq_total, chiq, RPA_deno] = RPA(qvec, cVar, freq_total, ion, chi0, const)
unitN = 4; % Number of magnetic atoms in unit cell
lattice = ion.abc{const.elem};
Vc = sum( lattice(1,:) .* cross(lattice(2,:), lattice(3,:)) ); % Volume of unit cell [Ang^-3]
eins = zeros(3,3,4,4);
eins(1,1,:,:) = 1; eins(2,2,:,:) = 1; eins(3,3,:,:) = 1;
chiq = zeros(3, 3, length(freq_total(1,:)), size(cVar,2), size(qvec,1));
D = zeros(3, 3, unitN, unitN, size(qvec,1));
parfor jj = 1:size(qvec,1)
    lorz_on = 1; % keep on the Lorentz term
    if abs(real(sum(exp(1i*2*pi*qvec(jj,:)*ion.tau'))/size(ion.tau,1))-1) > 1e-10
        lorz_on = 0;  % turn off Lorentz term at finite q = [h k l]
    end
%     D(:,:,:,:,jj) = const.gfac * (dipole_direct(qvec(jj,:), const.dpRng, lattice) + lorz_on*4*pi*eins/3/Vc)...
%         + exchange(qvec(jj,:), ion.ex(2), lattice, ion.tau); % [meV] dipole interaction
    D(:,:,:,:,jj) = const.gfac * (MF_dipole(qvec(jj,:), const.dpRng, lattice, ion.tau) + lorz_on*4*pi*eins/3/Vc)...
         - exchange(qvec(jj,:), abs(ion.ex(2)), lattice, ion.tau); % [meV]
end

RPA_deno = zeros(3, 3, size(freq_total,1), size(cVar,2)); % RPA correction factor (denominator)
for nb = 1:size(cVar,2) % field/temperature iterator
    for nq = 1:size(qvec,1) % q vector iterator
        Jav = squeeze( sum(sum(D(:,:,:,:,nq),4),3)/unitN ); % [meV] average over the unit cell
        Jq = diag(ion.renorm(const.elem,:)) .* Jav; % [meV]
        parfor nf = 1:length(freq_total(1,:))
%         for nf = 1:length(freq_total(1,:)) % for debugging
            chi_mf = squeeze(chi0(:,:,nf,nb));
            MM = chi_mf * Jq; % [meV^-1 * meV], non-commuting operation for matrices
            deno = squeeze(eye(size(MM))- MM);
            chiq(:,:,nf,nb,nq) = deno\chi_mf;
            RPA_deno(:,:,nf,nb,nq) = det(deno); % RPA denominator, save for pole analysis
        end
    end
end
% subs = ["xx", "yy", "zz"];
% for ii = 1:3
%     figure
%     hp0 = pcolor(continu_var,freq_total,squeeze((Mx(ii,ii,:,:))));
%     set(hp0, 'edgeColor','none')
%     xlabel('Magnetic field (T)')
%     ylabel('Frequency (GHz)')
%     
%     figure
%     hp0 = pcolor(continu_var,freq_total,squeeze(abs(deno(ii,ii,:,:))));
%     set(hp0, 'edgeColor','none')
%     xlabel('Magnetic field (T)')
%     ylabel('Frequency (GHz)')
%     for nq = 1:size(qvec,1)
%         figure
%         hp0 = pcolor(continu_var,freq_total,real(squeeze(chiq(ii,ii,:,:,nq))));
%         set(hp0, 'edgeColor','none');
%         colorbar
%         xlabel('Magnetic field (T)')
%         ylabel('Frequency (GHz)')
%         title(['Real part of \chi_',char(subs(ii)),'^{RPA}'])
%         figure
%         hp0 = pcolor(continu_var,freq_total,imag(squeeze(chiq(ii,ii,:,:,nq))));
%         set(hp0, 'edgeColor','none');
%         colorbar
%         xlabel('Magnetic field (T)')
%         ylabel('Frequency (GHz)')
%         title(['Imaginary part of \chi_',char(subs(ii)),'^{RPA}'])
%     end
% end
end