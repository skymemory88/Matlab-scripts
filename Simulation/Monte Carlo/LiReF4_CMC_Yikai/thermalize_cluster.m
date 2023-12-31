function [clusterSize, dE, accpRate, spin, coef, en0] = thermalize_cluster(ion, const, params, temp, en0, hamI, coef, spin)
pos = params.pos;
accpRate = 0;
origin = randi(size(pos,1)); % randomly pick a lattice site
nuc = pos(origin,:); % pick the nucleat site for the cluster
spin0 = spin(randi(size(pos,1)),:);

posT = pos - nuc; % shift the lattice center
cluster = zeros(size(posT,1),1); % cluster indices
for ii = 1:size(posT,1)
    prob = rand();
    if exp( -2 * dot(spin0, posT(ii,:)) * dot(spin(ii,:), posT(ii,:)) / const.kB / temp) > prob
        cluster(ii) = ii;
    end
end
sites = nonzeros(cluster); % cluster indices
clusterSize = length(sites); % cluster size

dipE0 = dipSum(const, ion, pos, spin); % energy before update
spin_new = spin; % copy the current spin config
spin_new(sites, 3) = -spin_new(sites, 3);
wav_new = coef;
En = en0;

dipE_new = dipSum(const, ion, pos, spin_new);
dEd = dipE_new - dipE0; % global change of dipole energy

% flip Jz on all sites in the cluster
dEi = zeros(length(sites),1); % local energy change
for ii = 1:length(sites)
    % compute change of the single-ion hamiltonian
    wav_new(:,sites(ii)) = flip(wav_new(:, sites(ii)));
    En(sites(ii)) = real(wav_new(:,sites(ii))' * hamI * wav_new(:, sites(ii))); % single-ion energy after update
    dEi(ii) = En(sites(ii)) - en0(ii); % change in single-ion energy
end
dE = sum(dEi) + dEd; % global energy change

% Glauber algorithm
prob = 0;
crit = rand; % update criteria
if temp > 0
    prob = 1 / ( 1 + exp(2*dE/const.kB/const.J2meV/temp) ); % probability critirion
end

if dE < 0 || prob > crit
    spin = spin_new;
    coef = wav_new;
    en0 = En;
    accpRate = 1; % acceptance rate
end
end