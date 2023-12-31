function [clusterSize, dE, accpRate, Esi, Edip, coef, eSpin] = thermalize_cluster_a(ion, const, params, tt, Esi, Edip, hamI, coef, eSpin, nSpin)
accpRate = 0;
cntr = randi(size(params.pos,1)); % randomly pick a lattice site
origin = params.pos(cntr,:); % pick the nucleat site for the cluster
spin0 = eSpin(randi(size(params.pos,1)),:);

posT = params.pos - origin; % shift the lattice center
cIndices = zeros(size(posT,1),1); % cluster indices
gfac = const.J2meV * const.mu0 / 4 / pi * (ion.gLande(ion.idx) * const.muB)^2;  % prefactor
beta = 1/ (const.J2meV * const.kB * params.temp(tt)); % [meV]
% construct the cluster
for ii = 1:size(posT,1)
    crit = rand();
    bond = gfac * (3*dot(spin0, posT(ii,:)*10^-10) * dot(eSpin(ii,:), posT(ii,:)*10^-10)...,
        /vecnorm(posT(ii,:) * 10^-10)^5 - dot(spin0, eSpin(ii,:)) / vecnorm(posT(ii,:) * 10^-10)^3); % [meV]
    if 1 - exp( -bond * beta ) > crit
        cIndices(ii) = ii; % store the index of the site added to the cluster
    end
end
cluster = nonzeros(cIndices); % cluster indices
clusterSize = length(cluster); % cluster size

spin_new = eSpin; % copy the current spin config
spin_new(cluster, 3) = -spin_new(cluster, 3);
wav_new = coef;
en = Esi;

dipE_old = dipSum(const, ion, params.pos, eSpin);
dipE_new = dipSum(const, ion, params.pos, spin_new);
dEd = dipE_new - dipE_old; % global change of dipole energy

% flip Jz on all sites in the cluster
dEi = zeros(length(cluster),1); % local energy change
for ii = 1:length(cluster)    
    % compute change of the single-ion hamiltonian
    wav_new(:,cluster(ii)) = flip(wav_new(:, cluster(ii)));
    en(cluster(ii)) = real(wav_new(:,cluster(ii))' * hamI * wav_new(:, cluster(ii))); % single-ion energy after update
    % thermalize nuclear spins to the original configurations
    if params.hyp
        en(cluster(ii)) = en(cluster(ii)) + ion.A(ion.idx) * sum(nSpin(cluster(ii),:) .* spin_new(cluster(ii),:)); % include hyperfine energy (meV)
    end
    dEi(ii) = en(cluster(ii)) - Esi(cluster(ii)); % change in single-ion energy
end
dE = sum(dEi) + dEd; % global energy change

% Glauber algorithm
prob = 0;
crit = rand; % update criteria
if params.temp(tt) > 0
    prob = 1 / ( 1 + exp(dE * beta) ); % probability critirion
end

if prob > crit
    eSpin = spin_new;
    coef = wav_new;
    Esi = en;
    for ii = 1:length(cluster)
        Edip(cluster(ii)) = CntrDip(const, ion, params.pos, eSpin, cluster(ii), eSpin(cluster(ii),:));
    end
    accpRate = 1; % acceptance rate
end
end