function varargout = linear_response(varargin)
% input: temperature, linewidth, frequency (range), hyperfine, saving option
% output: field (range), frequency (range), chi0, chiq, RPA denominator

Options.unit = 'meV'; % frequency (NOT energy) unit
% Options.unit = 'GHz'; % frequency (NOT energy) unit
Options.RPA = true;
filepath = ['G:\.shortcut-targets-by-id\1CapZB_um4grXCRbK6t_9FxyzYQn8ecQE\File sharing',...
        '\PhD program\Research projects\LiHoF4 project\Data\Simulations\Matlab\Susceptibilities',...
        '\toy_model\'];
temp = varargin{1};
const.gama = varargin{2};
freq_total = varargin{3};
hyp = varargin{4};
Options.saving = varargin{5};

load_name = sprintf('%1$3.3fK_*T_hp=%2$.2f', temp, hyp);
load_name = strcat('Hscan_toy_', load_name, '.mat');
%     load_name = strcat('Hscan_LHF-Ising_', load_name, '.mat');
MF_file = dir( fullfile(filepath, load_name) );
load([filepath MF_file(1).name],'ion','-mat','eee','fff','vvv');
% truncation of interaction dimension
if ~isfield(ion, 'mat')
    ion.mat = [1 1 1;
               1 1 1;
               1 1 1]; % default: Isotropic
%     ion.mat = [1 0 0;
%                0 1 0;
%                0 0 1]; % diagonal
%     ion.mat = [0 0 0;
%                0 0 0;
%                0 0 1]; % Ising
end
eee = squeeze(eee);
fields = fff;
En = eee - min(eee,[],2); % normalize the eigen-energies
wav = vvv;

const.hbar = 1.055E-34; % Reduced Planck constant [J.s]
const.mu0 = 4e-7*pi; % [H/m]
const.muB = 9.27401e-24; % Bohr magneton [J/T]
const.muN = 5.05078e-27; % Nuclear magneton [J/T]
const.J2meV = 6.24151e+21; % Joule to meV
const.kB = 8.61733e-2; % Boltzmann constant [meV/K]
switch Options.unit
    case 'GHz'
        const.f2E = const.hbar * 2*pi*10^9 * const.J2meV; % [meV/GHz]
    case 'meV'
        const.f2E = 1;
end

[chi] = MF_chi(ion, freq_total, fields, En, temp, wav, const);
if Options.RPA == true
    [chiq, poles] = RPA_chi(ion, En, freq_total, fields, [0 0 0], chi, const); % calculate RPA susceptibility
end
varargout{1} = fields;
varargout{2} = freq_total;
varargout{3} = chi;
if Options.saving == true && Options.RPA == false % Save the susceptibilities
    fields = fields';
    file_part = sprintf('%1$3.3fK_%2$.2e_%3$.3f-%4$.3fGHz_hyp=%5$.2f', temp,...
        const.gama, min(freq_total), max(freq_total), hyp);
    save_name = strcat('chi0_toy_', file_part, '.mat');
    %         save_name = strcat('chi0_LHF-Ising_', file_part, '.mat');
    savefile = fullfile(filepath, save_name);
    save(savefile, 'ion', 'temp', 'fields', 'freq_total', 'chi',...
        'const', 'Options', '-v7.3');
end

if exist('chiq','var')
    chi0 = chi;
    chi = chiq;
    varargout{4} = chiq;
    varargout{5} = poles;
    if Options.saving == true % Save the susceptibilities
        file_part = sprintf('%1$3.3fK_%2$.2e_%3$.3f-%4$.3fGHz', temp,...
            const.gama, min(freq_total), max(freq_total));
        save_name = strcat('RPA_toy_', file_part, '.mat');
%         save_name = strcat('RPA_LHF-Ising_', file_part, '.mat');
        savefile = fullfile(filepath, save_name);
        save(savefile, 'temp', 'fields', 'freq_total', 'chi0',...
            'const', 'chi', 'poles', 'Options', '-v7.3');
    end
end
end

function [chi0] = MF_chi(ion, freq_total, field, E, T, V, const)
fprintf('Computing MF susceptibility at T = %.3f K\n', T);
ELEf = ion.gLande * const.muB * const.J2meV; % Lande factor * Bohr magneton (meV/T)
% NUCf = ELEf; % [meV/T] for debugging
% NUCf = ion.nLande * muN * J2meV; % [meV/T]
NUCf = 0; % excluding nuclear contribution

%Initiate ionJ operators
Jz = diag(ion.J:-1:-ion.J); % Jz = -J, -J+1, ..., J-1, J
Jp = diag(sqrt( (ion.J-((ion.J-1):-1:-ion.J) ).*(ion.J+1+( (ion.J-1):-1:-ion.J) )), 1); % electronic spin ladder operator
Jm = Jp'; % electronic spin ladder operator
if any(ion.hyp) % hyperfine interaction option
    %Initiate I operators
    Iz = diag(ion.I:-1:-ion.I); %Iz = -I, -I+1, ..., I-1, I
    IhT.z = kron(eye(2*ion.J+1), Iz); % Expand Hilbert space
    Ip = diag(sqrt( (ion.I-((ion.I-1):-1:-ion.I) ).*(ion.I+1+((ion.I-1):-1:-ion.I))), 1); % Nuclear spin ladder operator
    Im = Ip'; % Nuclear spin ladder operator
    Iph = kron(eye(2*ion.J+1), Ip); % Expand to match the dimension of Hilbert space
    Imh = kron(eye(2*ion.J+1), Im);
    IhT.x = (Iph+Imh)/2;
    IhT.y = (Iph-Imh)/2i;

    JhT.z = kron(Jz, eye(2*ion.I+1)); % Expand Jz space to include nuclear degree of freedom
    Jph = kron(Jp, eye(2*ion.I+1)); % Expand to match the dimension of Hilbert space
    Jmh = kron(Jm, eye(2*ion.I+1));
    JhT.x = (Jph + Jmh)/2;
    JhT.y = (Jph - Jmh)/2i;
else
    JhT.x = (Jp + Jm)/2;
    JhT.y = (Jp - Jm)/2i;
    JhT.z = Jz;

    IhT.x = 0;
    IhT.y = 0;
    IhT.z = 0;
end

if isfield(ion, 'Czz') % only for projected effective spins
    Cxx = 2*real(ion.Cxx); % factor of 2 due to Pauli matrices' elements being 1 instead of 1/2
    Cyy = 2*real(ion.Cyy);
    Czz = 2*real(ion.Czz);
else
    Cxx = ones(size(field,1),1);
    Cyy = Cxx;
    Czz = Cyy;
end

chi0 = zeros(3,3,length(freq_total(1,:)),size(field,1));
for nb = 1:size(field,1) % calculate susceptibility for all fields
    wav = squeeze(V(nb,:,:,:)); % Obtain the corresponding eigen vectors
    ee = squeeze(E(nb,:,:)); % Obtain the corresponding eigen energies in meV
    beta = 1/(const.kB*T); %[meV^-1]
    if T == 0
        zz = zeros(size(ee));
        zz(1) = 1;
    elseif T > 0
        z = sum(exp(-beta*ee));
        zz = exp(-beta*ee)/z;
    end
    [pn, pm] = meshgrid(zz,zz);
    pmn = pn - pm; % population factor
%     pmn = ones(size(pn)); % for debugging
    [en, em] = meshgrid(ee, ee);
    parfor nw = 1:size(freq_total,2) %calculate susceptibility for all frequencies
        freq = freq_total(nw);
        dE = em - en - freq * const.f2E; % [meV]
        gamma = ones(size(dE)) * const.gama;
        deno0 = 1 ./ (dE - 1i*gamma); % [meV^-1]

        ttx1  = wav'  * (Cxx(nb) * JhT.x + NUCf/ELEf * IhT.x) * wav;
        tty1  = wav'  * (Cyy(nb) * JhT.y + NUCf/ELEf * IhT.y) * wav;
        ttz1  = wav'  * (Czz(nb) * JhT.z + NUCf/ELEf * IhT.z) * wav;

%         ttx1  = wav'  * (Cxx(ii) * JhT.x) * wav;
%         tty1  = wav'  * (Cyy(ii) * JhT.y) * wav;
%         ttz1  = wav'  * (Czz(ii) * JhT.z) * wav;
        
%         ttx1  = wav'  * (IhT.x) * wav;
%         tty1  = wav'  * (IhT.y) * wav;
%         ttz1  = wav'  * (IhT.z) * wav;

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

        chi0(:,:,nw,nb) = [xx xy xz
                           yx yy yz
                           zx zy zz];
    end
end
end

function [chiq, poles] = RPA_chi(ion, En, freq_total, field, qvec, chi0, const)
unitN = size(ion.tau,1); % Number of magnetic atoms in one unit cell
Vc = sum( ion.abc(1,:) .* cross(ion.abc(2,:), ion.abc(3,:)) ); % Volume of unit cell [Ang^-3]
gfac = const.mu0/4/pi * (ion.gLande * const.muB)^2 * const.J2meV * 1e30; % [meV*A^3]
eins = zeros(3,3,4,4);
eins(1,1,:,:) = 1; eins(2,2,:,:) = 1; eins(3,3,:,:) = 1;
Jij = zeros(3, 3, unitN, unitN, size(qvec,1));
switch ion.int
    case 'dipole'
        for nq = 1:size(qvec,1)
            Jij(:,:,:,:,nq) = gfac * (MF_dipole(qvec(nq,:), ion.dpRng, ion.abc, ion.tau) + 4*pi*eins/3/Vc); % [meV] dipole only
        end
    case 'exchange'
        for nq = 1:size(qvec,1)
            Jij(:,:,:,:,nq) = MF_exchange(qvec(nq,:), ion.ex, ion.abc, ion.tau); % [meV] exchange only
%             Jij(:,:,:,:,nq) = exchange(qvec(nq,:), abs(ion.ex), ion.abc, ion.tau);
        end
    case 'both'
        for nq = 1:size(qvec,1)
%             Jij(:,:,:,:,nq) = gfac * (MF_dipole(qvec(nq,:), ion.dpRng, ion.abc, ion.tau) + 4*pi*eins/3/Vc)...
%                 - MF_exchange(qvec(nq,:), ion.ex, ion.abc, ion.tau); % [meV]
            Jij(:,:,:,:,nq) = gfac * (MF_dipole(qvec(nq,:), ion.dpRng, ion.abc, ion.tau) + 4*pi*eins/3/Vc)...
                - exchange(qvec(nq,:), abs(ion.ex), ion.abc, ion.tau); % [meV]
        end
end
J = sum(sum(Jij,4),3)/unitN; % average over the unit cell
J = abs(J);

chiq = zeros(3, 3, length(freq_total(1,:)), size(field,1), size(qvec,1));
pks = zeros(3, 3, length(freq_total(1,:)), size(field,1), size(qvec,1));
poles = zeros(size(field,1), factorial(size(En,2)-1));
for nb = 1:size(field,1)
    if ismember(nb, 1:round(size(field,1)/5):size(field,1))
        fprintf('Computing RPA susceptibility at B = (%.1f %.1f %.1f) T\n',...
            field(nb,1), field(nb,2), field(nb,3));
    end
%     ee = squeeze(En(ii,:,:)); % Obtain the corresponding eigen energies in meV
%     [en, em] = meshgrid(ee, ee);
%     dE = (em - en - 1i*const.gama)./const.f2E; % [GHz]
    for nq = 1:size(qvec,1)
        Jq = squeeze(J(:,:,nq)) .* ion.mat; % truncate the interaction tensor
        parfor nw = 1:length(freq_total(1,:))
%             for kk = 1:length(freq_total(1,:))
            chi_mf = squeeze(chi0(:,:,nw,nb));
            MM = chi_mf * Jq; % [meV * meV^-1]
            deno = squeeze(eye(size(MM)) - MM);
            chiq(:,:,nw,nb,nq) = deno\chi_mf; % RPA susceptibility
            pks(:,:,nw,nb,nq) = inv(det(deno)); % take the inverse of the denominator tensor
%             pks(:,:,kk,ii,nq) = inv(det(deno*prod(dE,'all'))); % take the inverse of the denominator tensor
        end
        [~,loc] = findpeaks(squeeze(abs(pks(3,3,:,nb))), 'SortStr', 'none'); % search for poles along z
        poles(nb, 1:length(loc)) = flip(freq_total(loc));
    end
end
chiq = squeeze(chiq);
poles = squeeze(poles);
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