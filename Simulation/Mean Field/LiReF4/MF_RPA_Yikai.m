function  MF_RPA_Yikai(mion,scanMode,dscrt_var,freq_total,theta,phi,gama,hyp,RPA_mode)
% Current version assumes complete symmetrical equivalence among the four spin moments per unit cell
% mion: Magnetic ion type: 'Er', 'Ho'
% dscrt_var: discrete variable, can be temperature of field.
% freq_total frequency range in GHz
% theta: misalignment angle from c-axis
% phi: inplane angle in ab-plane between the magnetic field and a/b-axis [degree]
% gama: linewidth of hyperfine levels [meV].
% hyp: Hyperfine isotope proportion
% RPA_Mode: whether or not to apply RPA

Options.plotting = false; % Decide whether or not to plot the data at the end
Options.meV = false; % Energy unit choice: meV or GHz (default)
Options.saving = true; % Options to save the susceptibility tensors
Options.scanMode = scanMode; % 1. Field plot with RPA; 2. Temp plot with RPA; 2. wavevector plot with RPA
Options.nZee = false;
    if Options.nZee == true
        nZee_path = 'with Hz_I';
    else
        nZee_path = 'without Hz_I';
    end
Options.RPA = RPA_mode; % Apply random phase approximation (RPA) correction
    Qplot = false; % k-dependent plot for RPA susceptibilities
    if Qplot == true
        % qz = (6.9:0.0025:7.5)';
        qz = (0:0.01:0.5)';
        qx = zeros(size(qz,1),1);
        qy = zeros(size(qz,1),1);
        qvec = [qx qy qz];
        contnu_var0 = 0; % fixed point of the continuous variable for k-plot
    else
        qz = 0.0;
    %         qz = 7.305; % Wavelength = 0.1369 m, frequency = 2.19 GHz
    %         qz = 12.175; % Wavelength = 0.0821 m, frequency = 3.65 GHz
    %         qx = [0.01 0.1 0.3 0.6 1]';
        qy = zeros(size(qz,1),1);
        qx = zeros(size(qz,1),1);
        qvec = [qx qy qz];
    end
if Options.RPA == false
    Qplot = flase;
end
clearvars -except dscrt_var theta phi gama Options mion freq_total hyp qvec contnu_var0 Qplot nZee_path

global dip_range ex_range muB J2meV mu0 f2E elem_idx ELEf NUCf
% Declare physical constants as global for consistency
hbar = 1.055E-34; % Reduced Planck constant [J.s]
muN = 3.15245e-5; % Nuclear magneton [meV/T]
J2meV = 6.24151e+21; % Joule to meV
f2E = hbar*2*pi*10^9*J2meV;% GHz to meV
muB = 9.274e-24; %[J/T]
mu0 = 4e-7*pi; % [H/m]
dip_range = 100; % dipole summation range (number of unit cell)
ex_range = 1; % Exchange interaction range (number of unit cell)

for ii = 1:length(dscrt_var) 
    ipt = false;
    cntr = 0;
    while ~ipt
        switch Options.scanMode % 1. Field plot with RPA. 2. wavevector plot with RPA
            case 'field'
                location = ['G:\My Drive\File sharing\PhD program\Research projects\Li',mion,...
                            'F4 project\Data\Simulations\Matlab\Susceptibilities\', nZee_path];
                filename = strcat(['Hscan_Li',mion,'F4_'], sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_hp=%4$.2f', dscrt_var(ii), theta, phi, hyp),'.mat');
                file = fullfile(location,filename);
                load(file,'-mat','eee','fff','ttt','vvv','ion'); % loads variables "fields", "temp", "E" (eigenvalues) and "V" (eigenstates)
                file_part = sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_%4$.2e.mat', dscrt_var(ii), theta, phi, gama);
                continu_var = vecnorm(fff,2,1); % choose the continuous variable to be field
                % ttt = 0; % To account for only the lowest eigen-state contributiont
                % continu_var = continu_var(continu_var <= 5); % Truncate the field range
                ttt = repelem(ttt,length(continu_var));
                ipt = true;
                fprintf('Calculating for T = %.3f K.\n', dscrt_var(ii));
            case 'temp'
                location = ['G:\My Drive\File sharing\PhD program\Research projects\Li',mion,...
                    'F4 project\Data\Simulations\Matlab\Susceptibilities\', nZee_path];
                filename = strcat(['Tscan_Li',mion,'F4_'], sprintf('%1$3.3fT_%2$.2fDg_%3$.1fDg_hp=%4$.2f', dscrt_var(ii), theta, phi, hyp),'.mat');
                file = fullfile(location,filename);
                load(file,'-mat','eee','ttt','vvv','ion'); % loads variables "fields", "temp", "E" (eigenvalues) and "V" (eigenstates)
                file_part = sprintf('%1$3.3fT_%2$.2fDg_%3$.1fDg_%4$.2e.mat', dscrt_var(ii), theta, phi, gama);
                continu_var = ttt; % choose the continuous variable to be temperature
                ipt = true;
                fprintf('Calculating for B = %.3f T.\n', dscrt_var(ii));
            otherwise
                prompt = sprintf('Unknown Scan Mode! Please select a mode between tempscan and fieldscan.\n');
                answer = input(prompt,'s');
                Options.scanMode = lower(answer);
                cntr = cntr + 1; % input request counter
                if cntr >= 3
                    return % terminate the program after three failed requests
                end
        end
    end
    eee = eee - min(eee,[],2); % Normalize the energies against the ground state energy
    save_name1 = strcat('chi0_Li',mion,'F4_',file_part);
    save_name2 = strcat('RPA_Li',mion,'F4_',file_part);
    
    elem_idx = find(ion.prop);  
    ELEf = ion.gLande(elem_idx) * muB * J2meV; % Lande factor * Bohr magneton (meV/T)
    NUCf = ion.nLande(elem_idx) * muN;
    eigenE = eee; % load eigenenergy
%     eigenE(:,9:end) = 0; % truncate the Hilbert space
    eigenW = vvv; % load eigenstates
%     eigenW(:,9:end,9:end) = 0; % truncate the Hilbert space
    if Options.RPA == true
        if Qplot == true
            bidx = int16.empty(0,length(contnu_var0));
            for jj = 1:length(contnu_var0)
                [~,bidx(jj)] = min(abs(vecnorm(continu_var,1,1)-contnu_var0(jj)));
            end
            eigenW = vvv(bidx,:,:);
            eigenE = eee(bidx,:);
            continu_var = continu_var(bidx);
        end
        [continu_var, freq_total, chi, ~] = linear_response(ion,eigenE,continu_var,freq_total,ttt,eigenW,gama);
        [continu_var, freq_total, chiq_J ] = RPA(qvec, continu_var, freq_total, ion, chi.ionJ); % Electronic susceptibilities
        [continu_var, freq_total, chiq_IJ] = RPA(qvec, continu_var, freq_total, ion, chi.IJ); % cross term
        [continu_var, freq_total, chiq_I ] = RPA(qvec, continu_var, freq_total, ion, chi.ionI); % Nuclear susceptibilities

        chiq = ELEf^2*chiq_J + 2*ELEf*NUCf*chiq_IJ + NUCf^2*chiq_I; % Full electronuclear susceptibility expansion
%         chiq = ELEf^2*chiq.ionJ + 2*ELEf*NUCf*chiq.IJ; % Electron susceptibility + the cross term
%         chiq = ELEf^2*chiq.ionJ + NUCf^2*chiq.ionI; % Electron susceptibility + nuclear susceptibilities
%         chiq = ELEf^2*chiq.ionJ; % Electron susceptibility only
%         chiq = NUCf^2*chiq.ionI; % Nuclear susceptibility only
    else
        [continu_var, freq_total, chi, ~] = linear_response(ion,eigenE,continu_var,freq_total,ttt,eigenW,gama);
    end
    chi0 = ELEf^2*chi.ionJ + 2*ELEf*NUCf*chi.IJ + NUCf^2*chi.ionI; % Full electronuclear susceptibility expansion
%     chi0 = ELEf^2*chi0.ionJ + 2*ELEf*NUCf*chi0.IJ; % Electron susceptibility + the cross term
%     chi0 = ELEf^2*chi0.ionJ + NUCf^2*chi0.ionI; % Electron susceptibility + nuclear susceptibilities
%     chi0 = ELEf^2*chi0.ionJ; % Electron susceptibility only
%     chi0 = NUCf^2*chi0.ionI; % Nuclear susceptibility only

if Options.saving == true % Save the susceptibilities
    savefile1 = fullfile(location,save_name1);
    save_vars = {dscrt_var,continu_var,freq_total,chi0,gama,qvec};
    save_file(save_vars, savefile1,Options); % Save the data w/o RPA corrections
    if Options.RPA == true
        savefile2 = fullfile(location,save_name2);
        save_vars{4} = chiq;
        save_file(save_vars,savefile2,Options); % save the data w. RPA corrections
    end
end

if Options.plotting == true % Plot the susceptibilities
    plot_var = {continu_var, freq_total, chi0, gama, [0,0,0], dscrt_var(ii)};
    if Qplot == true
        plot_var{3} = chiq;
        plot_var{5} = qvec;
        figs(plot_var,Options,"\chi_{RPA}",Qplot);
        chi0_p = repmat(chi0,1,1,1,length(qvec));
        chi0_p = permute(chi0_p,[1 2 3 5 4]); % Add the additional dimension for the plotting purpose
        plot_var{3} = chi0_p;
        figs(plot_var,Options,"chi_0",Qplot);
    else
        figs(plot_var,Options,"\chi_0",Qplot);
        if Options.RPA == true
            plot_var{3} = chiq;
            figs(plot_var,Options,"\chi_{RPA}",Qplot);
            plot_var{3} = chiq-chi0;
            figs(plot_var,Options,"\chi_{RPA}-\chi_0",Qplot);
        end
    end
end
end
clearvars chi0 chi_p chiq
end

function figs(input_var,Options,fig_tit,Qplot)
continu_var = input_var{1};
freq_total = input_var{2};
chi = input_var{3};
gama = input_var{4};
qvec = input_var{5};
dscrt_var = input_var{6};
global f2E
if Options.meV == true % unit choice between meV and GHz
    freq_total = freq_total*f2E;
end
pos0 = [100 300 600 400]; % initial figure position
pos_inc = [150 0 0 0];
indx = [strcat(fig_tit,"^{xx}") strcat(fig_tit,"^{yy}") strcat(fig_tit,"^{zz}")]; % index for figure legend
if Qplot == true
    qvec = qvec(:,3); % use qz
    for ii = 1:size(continu_var,2)
        switch Options.scanMode
            case 'field'
                file_part1 = sprintf('T = %.3f K. and ', dscrt_var);
                file_part2 = sprintf('B = [%1$.2f %2$.2f %3$.2f] T.',continu_var(1,ii),continu_var(2,ii),continu_var(3,ii));
            case 'temp'
                file_part1 = sprintf('B = %.3f T. and ', dscrt_var);
                file_part2 = sprintf('B = [%1$.2f %2$.2f %3$.2f] T.',continu_var(1,ii),continu_var(2,ii),continu_var(3,ii));
            otherwise
                disp("Unknow scan mode!")
                return
        end
        for jj = 3 % plot only chi_zz component
            % for jj = 1:3
            % Color plot the imaginary part of the susceptibilities of z component
            fig(1) = figure;
            set(fig(1),'position',pos0 + (jj-1)*pos_inc);
            hp1 = pcolor(qvec,freq_total,real(squeeze(chi(jj,jj,:,ii,:))));
            set(hp1, 'edgeColor','none')
            set(gca, 'xdir', 'reverse' )
            caxis([0 100]);
            colorbar
            legend(['\gamma =' num2str(gama,'%.2e meV')]);
            title(['Re[', char(indx(jj)), '] at ', file_part1, file_part2])
            xlabel(sprintf('Q = [h, 0, 0]'))
            ylabel('Frequency (GHz)')
            if Options.meV == true
                ylabel('Energy (meV)')
            end         
            % Plot the real part of the susceptibility of the z component
            fig(2) = figure;
            set(fig(2),'position',pos0 + jj*pos_inc);
            hp2 = pcolor(qvec,freq_total,imag(squeeze(chi(jj,jj,:,ii,:))));
            set(hp2, 'edgeColor','none')
            set(gca, 'xdir', 'reverse' )
            caxis([0 100]);
            colorbar
            legend(['\gamma =' num2str(gama,'%.2e meV')]);
            title(['Im[', char(indx(jj)), '] at ', file_part1, file_part1])
            xlabel(sprintf('Q = [h, 0, 0]'))
            ylabel('Frequency (GHz)')
            if Options.meV == true
                ylabel('Energy (meV)')
            end
        end
    end
else
    for ii = 1:size(qvec,1)
        switch Options.scanMode
            case 'field'
                xlab = 'Magnetic Field (T)';
                file_part1 = sprintf('T = %.3f K and ', dscrt_var);
                file_part2 = sprintf('Q = [%1$.2f %2$.2f %3$.2f]',qvec(ii,1),qvec(ii,2),qvec(ii,3));
            case 'temp'
                xlab = 'Temperature (K)';
                file_part1 = sprintf('T = %.3f K and ', dscrt_var);
                file_part2 = sprintf('Q = [%1$.2f %2$.2f %3$.2f]',qvec(ii,1),qvec(ii,2),qvec(ii,3));
            otherwise
                disp("Unknow scan mode!")
                return
        end
        % Color plot the imaginary part of the susceptibilities
        for jj = 3 % plot only chi_zz component
%             for jj = 1:3
            fig(1) = figure;
            set(fig(1),'position',pos0 + (jj-1)*pos_inc);
%             hp2 = pcolor(fields(1,:),freq_total,squeeze(imag(chi(jj,jj,:,:,ii))));
            hp2 = pcolor(continu_var(1,:),freq_total,mag2db(abs(squeeze(imag(chi(jj,jj,:,:,ii))))));
            set(hp2, 'edgeColor','none')
            caxis([0 5]);
            colorbar
            legend(['\gamma =' num2str(gama,'%.2e meV')]);
            title(['Im[', char(indx(jj)), '] at ', file_part1, file_part2])
            xlabel(xlab)
            ylabel('Frequency (GHz)')
            if Options.meV == true
                ylabel('Energy (meV)')
            end            
            % Plot the real part of the susceptibility of the z component
            fig(2) = figure;
            set(fig(2),'position',pos0 + jj*pos_inc);
%             hp3 = pcolor(fields(1,:),freq_total,squeeze(real(chi(jj,jj,:,:,ii))));
            hp3 = pcolor(continu_var(1,:),freq_total,mag2db(abs(squeeze(real(chi(jj,jj,:,:,ii))))));
            set(hp3, 'edgeColor','none')
            caxis([0 5]);
            colorbar
            legend(['\gamma =' num2str(gama,'%.2e meV')]);
            title(['Re[', char(indx(jj)), '] at ', file_part1, file_part2])
            xlabel(xlab)
            ylabel('Frequency (GHz)')
            if Options.meV == true
                ylabel('Energy (meV)')
            end
        end
    end
end
end

function save_file(save_vars,savefile,opt)
switch opt.scanMode
    case 'field'
        temp = save_vars{1};
        fields = save_vars{2};
    case 'temp'
        fields = save_vars{1};
        temp = save_vars{2};
end
freq_total = save_vars{3};
chi = save_vars{4};
gama = save_vars{5};
qvec = save_vars{6};
save(savefile,'temp','fields','freq_total','chi','gama','qvec','-v7.3');
end

function [continu_var, freq_total, chi, JIz_exp] = linear_response(ion,eigenE,continu_var,freq_total,temperature,eigenW,gma)
global ELEf NUCf f2E elem_idx
% Declare susceptibility tensor
chi0_J = zeros(3,3,length(freq_total(1,:)),size(continu_var,2));
chi0_I = zeros(3,3,length(freq_total(1,:)),size(continu_var,2));
chi0_IJ = zeros(3,3,length(freq_total(1,:)),size(continu_var,2));
ionJ = ion.J(elem_idx);
ionI = ion.I(elem_idx);
%Initiate ionJ operators
Jz=diag(ionJ:-1:-ionJ); % Jz = -J, -J+1,...,J-1,J
JhT.z=kron(Jz,eye(2*ionI+1)); % Expand Jz space to include nuclear degree of freedom
Jp=diag(sqrt((ionJ-((ionJ-1):-1:-ionJ) ).*(ionJ+1+( (ionJ-1):-1:-ionJ) )),1); % electronic spin ladder operator
Jm=Jp'; % electronic spin ladder operator
Jph=kron(Jp,eye(2*ionI+1)); % Expand to match the dimension of Hilbert space
Jmh=kron(Jm,eye(2*ionI+1));
JhT.x=(Jph+Jmh)/2;
JhT.y=(Jph-Jmh)/2i;

%Initiate I operators
Iz=diag(ionI:-1:-ionI); %Iz = -I, -I+1,...,I-1,I
IhT.z=kron(eye(2*ionJ+1),Iz); % Expand Hilbert space
Ip=diag(sqrt((ionI-((ionI-1):-1:-ionI)).*(ionI+1+((ionI-1):-1:-ionI))),1); % Nuclear spin ladder operator
Im=Ip'; % Nuclear spin ladder operator
Iph=kron(eye(2*ionJ+1),Ip); % Expand to match the dimension of Hilbert space
Imh=kron(eye(2*ionJ+1),Im);
IhT.x=(Iph+Imh)/2;
IhT.y=(Iph-Imh)/2i;
        
% IJ_hT.x = JhT.x+NUCf/ELEf*IhT.x; % Hybridized electronuclear spin operator
% IJ_hT.y = JhT.y+NUCf/ELEf*IhT.y;
% IJ_hT.z = JhT.z+NUCf/ELEf*IhT.z;

% Single out <Jz+Iz> calculations
JIz_exp = double.empty(size(continu_var,2),size(JhT.z,1),0); % Expectation value of J-I pseudo-spin
for kk = 1:size(continu_var,2) % calculate susceptibility for all fields
    v = squeeze(eigenW(kk,:,:)); % Obtain the corresponding eigen vectors
    JzhT = JhT.z * ELEf;
    IzhT = IhT.z * NUCf;
%     JzhT = JhT.z;
%     IzhT = 0;
    ttz  = v'  * (JzhT+IzhT) * v;
    JIz_exp(kk,:,1) = real(diag(ttz));
end

for m = 1:length(freq_total(1,:)) %calculate susceptibility for all frequencies
    freq = freq_total (m);
    omega = freq*f2E;   % define frequency sweep range (meV)
%     for k = 1:size(continu_var,2) % for debugging: calculate susceptibility for all fields
    parfor k = 1:size(continu_var,2) % calculate susceptibility for all fields
        v = squeeze(eigenW(k,:,:)); % Obtain the corresponding eigen vectors
        en = squeeze(eigenE(k,:)); % Obtain the corresponding eigen energies [meV]
        if temperature(k) ~= 0
            beta = 11.6/temperature(k); %[meV^-1]
            zn = sum(exp(-beta*en));
            Z = exp(-beta*en)/zn;
%             Z = exp(-beta*en);
%             for nn = 1:8 % Skew the population of the lowest 8 states
%                Z(nn) = Z(nn) - (9-nn)/2e2 + nn/2e2;
%             end
%             Z = Z/sum(Z);
            [n,np] = meshgrid(Z,Z);
            NN = n-np;
        else
            zn = zeros(size(en));
            zn(1) = 1;
            [n,np] = meshgrid(zn,zn);
            NN = n-np;
        end
        [ee,eep] = meshgrid(en,en);
        EE = eep-ee-omega;
        gamma = ones(size(EE))*gma;

%         deno1 = EE ./ (EE.^2 + gamma.^2);
%         deno2 = gamma ./ (EE.^2 + gamma.^2);
%         chi0_J(:,:,m,k) = chi_Mx(JhT,JhT,v,NN,deno1,deno2); % Electornic spin operators
%         chi0_IJ(:,:,m,k) = chi_Mx(IhT,JhT,v,NN,deno1,deno2); % Electro-Nuclear cross term
%         chi0_I(:,:,m,k) = chi_Mx(IhT,IhT,v,NN,deno1,deno2); % Nuclear spin operators        
% %         chi0_J(:,:,m,k) = chi_Mx(IJ_hT,IJ_hT,v,NN,deno1,deno2); % hyperdized electronuclear operator

        EE = eep-ee-omega;
        deno0 = 1 ./ (EE - 1i*gamma);
        chi0_J(:,:,m,k) = chi_Mx(JhT,JhT,v,NN,deno0); % Electornic spin operators
        chi0_IJ(:,:,m,k) = chi_Mx(JhT,IhT,v,NN,deno0); % Electro-Nuclear cross term
        chi0_I(:,:,m,k) = chi_Mx(IhT,IhT,v,NN,deno0); % Nuclear spin operators
%         chi0_J(:,:,m,k) = chi_Mx(IJ_hT,IJ_hT,v,NN,deno0); % hyperdized electronuclear operator
    end
end

chi.ionJ = chi0_J; % Electronic susceptibilities
chi.ionI = chi0_I; % Nuclear susceptibilities
chi.IJ = chi0_IJ; % Electronuclear Cross term
end

function chi0 = chi_Mx(varargin)
% Calculation of matrix elements for the real and imaginary susceptibilities
% opr1: operator 1
% opr2: operator 2
% wav: eigenfunction of the hamiltonian
% NN: population differene between two eigenstates
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

function [continu_var, freq_total, chiq] = RPA(qvec, continu_var, freq_total, ion, chi0)

global muB mu0 J2meV dip_range ex_range elem_idx
N = 4; % Number of magnetic atoms in unit cell
lattice = ion.abc{elem_idx};
gfac = (muB)^2*(mu0/4/pi)*J2meV*10^30;

chiq = zeros(3,3,length(freq_total(1,:)),size(continu_var,2),size(qvec,1));
D = zeros(3,3,N,N,size(qvec,1));
for jj=1:size(qvec,1)
    D(:,:,:,:,jj) = ion.gLande(elem_idx)^2*(gfac*dipole_direct(qvec(jj,:),dip_range,lattice))...
        + exchange(qvec(jj,:),ion.ex(2),lattice,ex_range);
end

deno = zeros(3,3,size(freq_total,1),size(continu_var,2)); % RPA correction factor (denominator)
for ii = 1:size(continu_var,2)
    for nq = 1:size(qvec,1)
        Jq = sum(sum(D(:,:,:,:,nq),4),3)/4; % average over the four sites in the unit cell and avoid double counting
        parfor kk = 1:length(freq_total(1,:))
            MM = chi0(:,:,kk,ii).*Jq;
            deno(:,:,kk,ii) = (ones(size(MM))- MM); % Suppress parfor warning for this line
            chiq(:,:,kk,ii,nq) = squeeze(deno(:,:,kk,ii)).\chi0(:,:,kk,ii);
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