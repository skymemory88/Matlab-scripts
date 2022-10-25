function varargout = linear_response(varargin)
Options.unit = 'GHz';
Options.RPA = true;
filepath = ['G:\.shortcut-targets-by-id\1CapZB_um4grXCRbK6t_9FxyzYQn8ecQE\File sharing',...
        '\PhD program\Research projects\LiHoF4 project\Data\Simulations\Matlab\Susceptibilities',...
        '\toy_model\'];

if nargin == 4
    temp = varargin{1};
    gama = varargin{2};
    qvec = varargin{3};
    Options.saving = varargin{4};

    load_name = sprintf('%1$3.3fK_%2$.2e_*GHz', temp, gama);
    load_name = strcat('chi0_toy_', load_name, '.mat');
    chi_file = dir( fullfile(filepath, load_name) );
    load([filepath chi_file.name],'-mat','ion','fields','freq_total','chi0');

    [chiq, RPA_deno] = RPA_chi(ion, freq_total, fields, qvec, chi0); % calculate RPA susceptibility
    varargout{1} = fields;
    varargout{2} = freq_total;
    varargout{3} = chiq;
elseif nargin == 5
    temp = varargin{1};
    freq_total = varargin{2};
    gama = varargin{3};
    hyp = varargin{4};
    Options.saving = varargin{5};

    load_name = sprintf('%1$3.3fK_*T_hp=%2$.2f', temp, hyp);
    load_name = strcat('Hscan_toy_', load_name, '.mat');
    MF_file = dir( fullfile(filepath, load_name) );
    load([filepath MF_file.name],'ion','-mat','eee','fff','vvv');
    fields = fff;
    En = eee-min(eee,[],2); % normalize the eigen-energies
    wav = vvv;

    [chi0] = MF_chi(ion, freq_total, fields, En, temp, wav, gama, Options);
    if Options.RPA == true
        [chiq, poles, RPA_deno] = RPA_chi(ion, freq_total, fields, [0 0 0], chi0); % calculate RPA susceptibility
    end
    varargout{1} = fields;
    varargout{2} = freq_total;
    varargout{3} = chi0;
    if Options.saving == true % Save the susceptibilities
        file_part = sprintf('%1$3.3fK_%2$.2e_%3$.3f-%4$.3fGHz', temp,...
            gama, min(freq_total), max(freq_total));
        save_name = strcat('chi0_toy_', file_part, '.mat');
        savefile = fullfile(filepath, save_name);
        save(savefile, 'ion', 'temp', 'fields', 'freq_total', 'chi0',...
            'gama', 'Options', '-v7.3');
    end
end

if exist('chiq','var')
    varargout{4} = chiq;
    varargout{5} = poles;
    varargout{6} = RPA_deno;
    if Options.saving == true % Save the susceptibilities
        file_part = sprintf('%1$3.3fK_%2$.2e_%3$.3f-%4$.3fGHz', temp,...
            gama, min(freq_total), max(freq_total));
        save_name = strcat('RPA_toy_', file_part, '.mat');
        savefile = fullfile(filepath, save_name);
        save(savefile, 'temp', 'fields', 'freq_total', 'chi0',...
            'gama', 'chiq', 'poles', 'RPA_deno', 'Options', '-v7.3');
    end
end
end

function [chi0] = MF_chi(ion, freq_total, field, E, T, V, gama, Options)
fprintf('Computing MF susceptibility at T = %.2f K\n', T);
muB = 9.27401e-24; % Bohr magneton [J/T]
muN = 5.05078e-27; % Nuclear magneton [J/T]
J2meV = 6.24151e+21; % Joule to meV
switch Options.unit
    case 'meV'
        hbar = 1.055E-34; % Reduced Planck constant [J.s]
        J2meV = 6.24151e+21; % [meV/J]
        f2E = hbar*2*pi*10^9*J2meV; % [meV/GHz]
    case 'GHz'
        f2E = 1;
end
ELEf = ion.gLande * muB * J2meV; % Lande factor * Bohr magneton (meV/T)
NUCf = ion.nLande * muN * J2meV; % [meV/T]

%Initiate ionJ operators
ionJ = ion.J;
Jz = diag(ionJ:-1:-ionJ); % Jz = -J, -J+1,...,J-1,J
Jp = diag(sqrt((ionJ-((ionJ-1):-1:-ionJ) ).*(ionJ+1+( (ionJ-1):-1:-ionJ) )),1); % electronic spin ladder operator
Jm = Jp'; % electronic spin ladder operator

if any(ion.hyp) % hyperfine interaction option
    %Initiate I operators
    ionI = ion.I;
    Iz = diag(ionI:-1:-ionI); %Iz = -I, -I+1,...,I-1,I
    IhT.z = kron(eye(2*ionJ+1),Iz); % Expand Hilbert space
    Ip = diag(sqrt((ionI-((ionI-1):-1:-ionI)).*(ionI+1+((ionI-1):-1:-ionI))),1); % Nuclear spin ladder operator
    Im = Ip'; % Nuclear spin ladder operator
    Iph = kron(eye(2*ionJ+1),Ip); % Expand to match the dimension of Hilbert space
    Imh = kron(eye(2*ionJ+1),Im);
    IhT.x = (Iph+Imh)/2;
    IhT.y = (Iph-Imh)/2i;

    JhT.z = kron(Jz,eye(2*ionI+1)); % Expand Jz space to include nuclear degree of freedom
    Jph = kron(Jp,eye(2*ionI+1)); % Expand to match the dimension of Hilbert space
    Jmh = kron(Jm,eye(2*ionI+1));
    JhT.x = (Jph+Jmh)/2;
    JhT.y = (Jph-Jmh)/2i;
else
    JhT.x = (Jp+Jm)/2;
    JhT.y = (Jp-Jm)/2i;
    JhT.z = Jz;

    IhT.x = 0;
    IhT.y = 0;
    IhT.z = 0;
end

chi0 = zeros(3,3,length(freq_total(1,:)),size(field,1));
for ii = 1:size(freq_total,2) %calculate susceptibility for all frequencies
    freq = freq_total(ii);
    omega = freq*f2E; % define frequency sweep range (meV)
%     parfor jj = 1:size(field,1) % calculate susceptibility for all fields
    for jj = 1:size(field,1) % for debugging
        wav = squeeze ( squeeze(V(jj,:,:,:)) ); % Obtain the corresponding eigen vectors
        ee = squeeze ( squeeze(E(jj,:,:)) ); % Obtain the corresponding eigen energies in meV
        beta = 1/(T*8.61733E-2); %[meV^-1]
        if T == 0
            zz = zeros(size(ee));
            zz(1) = 1;
        elseif T > 0
            z = sum(exp(-beta*ee));
            zz = exp(-beta*ee)/z;
        end
        [pn, pm] = meshgrid(zz,zz);
        pmn = pn - pm; % population factor
%         pmn = ones(size(pn)); % for debugging
        [en, em] = meshgrid(ee,ee);
        dE = em-en-omega;
        gamma = ones(size(dE))*gama;
        deno0 = 1 ./ (dE - 1i*gamma); % [meV^-1]

        ttx1  = wav'  * (JhT.x + NUCf/ELEf * IhT.x) * wav;
        tty1  = wav'  * (JhT.y + NUCf/ELEf * IhT.y) * wav;
        ttz1  = wav'  * (JhT.z + NUCf/ELEf * IhT.z) * wav;

        % Calculate susceptibilities along a-axis
        chi_tx  = (ttx1) .* (ttx1.') .* pmn .* deno0;
        xx = sum(sum(chi_tx));

        % Calculate susceptibilities along b-axis
        chi_ty  = (tty1) .* (tty1.') .* pmn .* deno0;
        yy = sum(sum(chi_ty));

        % Calculate susceptibilities along c-axis
        chi_tz  = (ttz1) .* (ttz1.') .* pmn .* deno0;
        zz = sum(sum(chi_tz));

        % Calculate susceptibilities along ab-axis
        chi_txy  = (ttx1) .* (tty1.') .* pmn .* deno0;
        xy = sum(sum(chi_txy));

        % Calculate susceptibilities along ac-axis
        chi_txz  = (ttx1) .* (ttz1.') .* pmn .* deno0;
        xz = sum(sum(chi_txz));

        % Calculate susceptibilities along ba-axis
        chi_tyx  = (tty1) .* (ttx1.') .* pmn .* deno0;
        yx = sum(sum(chi_tyx));

        % Calculate susceptibilities along bc-axis
        chi_tyz  = (tty1) .* (ttz1.') .* pmn .* deno0;
        yz = sum(sum(chi_tyz));

        % Calculate susceptibilities along ca-axis
        chi_tzx  = (ttz1) .* (ttx1.') .* pmn .* deno0;
        zx = sum(sum(chi_tzx));

        % Calculate susceptibilities along cb-axis
        chi_tzy  = (ttz1) .* (tty1.') .* pmn .* deno0;
        zy = sum(sum(chi_tzy));

        chi0(:,:,ii,jj) = [xx xy xz
                           yx yy yz
                           zx zy zz];
    end
end
end

function [chiq, pole, RPA_deno] = RPA_chi(ion, freq_total, field, qvec, chi0)
unitN = size(ion.tau,1); % Number of magnetic atoms in one unit cell
muB = 9.274e-24; % [J/T]
mu0 = 4e-7*pi; % [H/m]
J2meV = 6.24151e+21; % [mev/J]
% gfac = ion.gLande^2 * muB^2 * mu0 * J2meV / (4*pi); % (gL * const.muB)^2 * const.mu0/(4pi) [meV.m^3]
gfac = mu0/4/pi * muB^2 * J2meV; % [meV.m^3] removed duplicate Lande factor

Jij = zeros(3,3,unitN,unitN,size(qvec,1));
chiq = zeros(3,3,length(freq_total(1,:)),size(field,1),size(qvec,1));
deno = zeros(3,3,length(freq_total(1,:)),size(field,1),size(qvec,1));
switch ion.int
    case 'dipole'
        for jj = 1:size(qvec,1)
            Jij(:,:,:,:,jj) = gfac * 1e30 * MF_dipole(qvec(jj,:), ion.dpRng, ion.abc, ion.tau); % [meV] dipole only
        end
    case 'exchange'
        for jj = 1:size(qvec,1)
            Jij(:,:,:,:,jj) = MF_exchange(qvec(jj,:), ion.ex, ion.abc, ion.tau); % [meV] exchange only
        end
    case 'both'
        for jj = 1:size(qvec,1)
            Jij(:,:,:,:,jj) = gfac*10^30*MF_dipole(qvec(jj,:), dpRng, ion.abc, ion.tau)...
                + MF_exchange(qvec(jj,:), ion.ex, ion.abc, ion.tau); % [meV]
        end
end

Jex = zeros(3,3,size(qvec,1));
for ionn = 1:unitN
    for ionm = 1:ionn
        Jex = Jex + Jij(:,:,ionn,ionm,:); % average over all sites within one unit cell
    end
end

for ii = 1:size(field,1)
    if ismember(ii, [1:round(size(field,1)/5):size(field,1)])
        fprintf('Computing RPA susceptibility at B = (%.1f %.1f %.1f) T\n',...
            field(ii,1), field(ii,2), field(ii,3));
    end
    for nq = 1:size(qvec,1)
        Jq = squeeze(Jex(:,:,nq) .* ion.mat);  % truncate interaction matrix (Ising, isotropic, diagonal)
        Jq = abs(Jq);
        parfor kk = 1:length(freq_total(1,:))
%         for kk = 1:length(freq_total(1,:))
            chi_mf = squeeze(chi0(:,:,kk,ii));
            MM = chi_mf * Jq; % [meV*meV^-1]
            deno(:,:,kk,ii,nq) = eye(size(MM)) - MM; % Suppress parfor warning for this line
            chiq(:,:,kk,ii,nq) = deno(:,:,kk,ii,nq)\chi_mf;
        end
    end
end
chiq = squeeze(chiq);
RPA_deno = squeeze(deno);

chizq = squeeze(mag2db(imag(chiq(3,3,:,:))));
[~,p0] = findpeaks(chizq(:,1)); % determine the maximum number of poles
pole = zeros(size(fields,1), length(p0));
pole(1,:) = p0;
for ii = 2:size(fields,1)
    [~,loc] = findpeaks(chizq(:,ii));
    for jj = 1:length(loc)
        pole(ii,jj) = freq_total(loc(jj));
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