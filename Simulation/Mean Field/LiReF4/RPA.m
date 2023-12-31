function [qvec, fields, freq_total, chi, deno] = RPA(varargin)
Options.saving = false; % saving option
Options.plot = true; % plot option

const.muB = 9.27401e-24; % Bohr magneton [J/T]
const.muN = 5.05078e-27; % Nuclear magneton [J/T]
const.mu0 = 4*pi*1e-7; % [H/m]
const.hbar = 1.05457E-34; % Reduced Planck constant [J.s]
const.J2meV = 6.24151e+21; % Joule to meV
const.Gh2mV = const.hbar * 2*pi * 1e9 * const.J2meV; % convert GHz to meV
const.kB = 8.61733e-2; % Boltzmann constant [meV/K];

mion = varargin{1};
temp = varargin{2};
theta = varargin{3};
phi = varargin{4};
gama = varargin{5};
Options.hyperfine = varargin{6}; % hyperfine interaction option
Options.nZee = varargin{6}; % nuclear Zeeman interaction option

if Options.nZee == true
    nZee_path = 'Hz_I=1';
else
    nZee_path = 'Hz_I=0';
end

if strcmp(pathsep, ':')
    platform = 'Unix';
else
    platform = 'Win';
end

switch platform
    case 'Win'
        location = ['G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Simulations\Matlab',...
            '\Susceptibilities\',nZee_path];
    case 'Unix'
        Options.location = ['/Users/yikaiyang/Library/CloudStorage/GoogleDrive-yikai.yang@epfl.ch/My Drive/'...
            'File sharing/PhD program/Research projects/LiHoF4 project/Data/',...
            'Simulations/MATLAB/Susceptibilities/', nZee_path];
end

% load the MF and MF susceptibility files
MF_name = strcat(['Hscan_Li',mion,'F4_'],...
    sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_hp=%4$.2f', temp, theta, phi, Options.hyperfine),'.mat');
MF_file = fullfile(location, MF_name);
load(MF_file,'-mat','ion'); % loads variables "fields", "temp", "E" (eigenvalues) and "V" (eigenstates)
unitN = size(ion.tau,1); % moment per unit cell

chi_name = sprintf('chi0_LiHoF4_%1$3.3fK_%2$.2fDg_%3$.1fDg_%4$.2e_*GHz.mat', temp, theta, phi, gama);
chi_file = dir( fullfile(location, chi_name) ); % susceptibility data to load
load([location chi_file.name],'-mat','fields','freq_total','chi','unit');

[~, elem] = ismember(mion, ion.name); % find appropriate element index
if temp > 0; beta = 1/const.kB/temp; else; beta = 1; end % [meV^-1]
lattice = ion.abc{elem}; % retrieve lattice constants
Vc = sum( lattice(1,:) .* cross(lattice(2,:), lattice(3,:)) ); % Volume of unit cell [Ang^-3]
gfac = const.mu0/4/pi * (ion.gLande(elem) * const.muB)^2 * const.J2meV * 10^30; % mu0/(4pi).(gL.muB)^2 [meV.Ang^3]
ELEf = ion.gLande(elem) * const.muB * const.J2meV; % Lande factor * Bohr magneton (meV/T)
chi0 = chi./ELEf^2; % [1/meV]

if nargin == 8 % for q-plot at constant field
    qvec = varargin{7};
    B0 = varargin{8};
    bidx = int32.empty(length(B0),0);
    for ii = 1:length(B0)
        [~, bidx(ii)] = min(abs(fields - B0(ii)));
        B0(ii) = fields(bidx(ii));
    end
else
    qvec = [0 0 0];
    bidx = 1:size(fields,2);
end

chiq = zeros(3, 3, length(freq_total), length(bidx), size(qvec,1));
deno = zeros(3, 3, 4, length(freq_total), length(bidx), size(qvec,1)); % RPA correction factor (denominator)
eins = zeros(3,3,4,4);
eins(1,1,:,:) = 1; eins(2,2,:,:) = 1; eins(3,3,:,:) = 1;
D = zeros(3, 3, unitN, unitN, size(qvec,1)); % dipole interaction container
J = zeros(3, 3, unitN, unitN, size(qvec,1)); % dipole interaction container
parfor jj = 1:size(qvec,1)
% for jj = 1:size(qvec,1) % for debugging
%     D(:,:,:,:,jj) = gfac * MF_dipole(qvec(jj,:), ion.dpRng, lattice, ion.tau); % [meV] dipole interaction summation
    D(:,:,:,:,jj) = gfac * (MF_dipole(qvec(jj,:), ion.dpRng, lattice, ion.tau) + 4*pi*eins/3/Vc); % dipole + Lorenz term
    J(:,:,:,:,jj) = exchange(qvec(jj,:), abs(ion.ex(2)), lattice, ion.tau); % [meV] AFM exchange interaction
end
% Jq = zeros(3,3,unitN,unitN,size(qvec,1));
% Jq(3,3,:,:,:) = repmat(diag(ion.renorm(elem,:)),1,1,4,4) .*...
%                   (D(3,3,:,:,:) - J(3,3,:,:,:)); % truncate to Ising interaction
Jq = repmat(diag(ion.renorm(elem,:)),1,1,4,4) .* (D - J);

counter = 1;
chii = zeros(3,3,4); % susceptibility for individual site within the unit cell
screw = rotz(0); % rotate from site to site in the unit cell
for nb = bidx % field iterator
    for nf = 1:size(freq_total,2) 
        chii(:,:,1) = beta*squeeze(chi0(:,:,nf,nb));
        chii(:,:,2) = beta*screw*squeeze(chii(:,:,1));
        chii(:,:,3) = beta*screw*squeeze(chii(:,:,2));
        chii(:,:,4) = beta*screw*squeeze(chii(:,:,3));
        parfor nq = 1:size(qvec,1)
%         for nq = 1:size(qvec,1)
          M = zeros(3*unitN); % Block matrix: 3 x 3 for each ion in the unit cell  
          for ionn = 1:unitN  % Iterate through all ions
            for ionm = 1:unitN
              M((ionn-1)*3+(1:3),(ionm-1)*3+(1:3)) = chii(:,:,ionn)... % use the same chi0 for all ions
                  *squeeze(Jq(:,:,ionn,ionm,nq)); % matrix product, non-commuting
            end
            deno(:,:,ionn,nf,nb,nq) = (eye(3)-M((ionn-1)*3+(1:3),3));
          end
          chiq(:,:,nf,counter,nq) = 1/4*[eye(3) eye(3) eye(3) eye(3)]*...
            ((eye(size(M))-M)\([chii(:,:,1);chii(:,:,2);chii(:,:,3);chii(:,:,4)]));
        end
    end
    counter = counter + 1;
end
chi = ELEf^2 .* chiq; % [meV/T^2]

if Options.saving == true % Save the susceptibilities
    RPA_deno = squeeze(deno);
    file_part = sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_%4$.2e_%5$.3f-%6$.3fGHz', temp,...
        theta, phi, gama, min(freq_total), max(freq_total));
    save_name = strcat('RPA_Li',mion,'F4_',file_part, '.mat');
    savefile = fullfile(location,save_name);
    save(savefile, 'temp', 'fields', 'freq_total', 'chi', 'gama', 'qvec', 'unit', 'RPA_deno', 'const', '-v7.3');
end

if Options.plot == true % incomplete
    subs = ["xx", "yy", "zz"];
    for ii = 1:3
        figure
        hp0 = pcolor(continu_var,freq_total,squeeze((Mx(ii,ii,:,:))));
        set(hp0, 'edgeColor','none')
        xlabel('Magnetic field (T)')
        ylabel('Frequency (GHz)')
        
        figure
        hp0 = pcolor(continu_var,freq_total,squeeze(abs(RPA_deno(ii,ii,:,:))));
        set(hp0, 'edgeColor','none')
        xlabel('Magnetic field (T)')
        ylabel('Frequency (GHz)')
        for nq = 1:size(qvec,1)
            figure
            hp0 = pcolor(continu_var,freq_total,real(squeeze(chiq(ii,ii,:,:,nq))));
            set(hp0, 'edgeColor','none');
            colorbar
            xlabel('Magnetic field (T)')
            ylabel('Frequency (GHz)')
            title(['Real part of \chi_',char(subs(ii)),'^{RPA}'])
            figure
            hp0 = pcolor(continu_var,freq_total,imag(squeeze(chiq(ii,ii,:,:,nq))));
            set(hp0, 'edgeColor','none');
            colorbar
            xlabel('Magnetic field (T)')
            ylabel('Frequency (GHz)')
            title(['Imaginary part of \chi_',char(subs(ii)),'^{RPA}'])
        end
    end
end

end