function [clusterSize, dE, accpRate, Esi, coef, eSpin] = thermalize_cluster_a(ion, const, params, tt, Esi, hamI, coef, eSpin, nSpin)

cntr = randi(size(params.pos,1)); % randomly pick a lattice site
origin = params.pos(cntr,:); % pick the nucleat site for the cluster
spin0 = eSpin(cntr,:); % spin of the nucleat site

dVec = (params.pos - origin) * 10^-10; % [m] displacement vector of all the sites from the nucleat site
cIndices = zeros(size(dVec,1),1); % cluster indices
gfac = const.J2meV * const.mu0 / 4 / pi * (ion.gLande(ion.idx) * const.muB)^2;  % prefactor
beta = 1/ (const.J2meV * const.kB * params.temp(tt)); % [meV]
% construct the cluster
for ii = 1:size(params.pos,1)
    crit = rand();
    bond = gfac * (dot(spin0, eSpin(ii,:)) / vecnorm(dVec(ii,:))^3 -...
        3*dot(spin0, dVec(ii,:)) * dot(eSpin(ii,:), dVec(ii,:))/vecnorm(dVec(ii,:))^5 ); % [meV]
    % bond = gfac * dot(spin0, eSpin(ii,:)) / vecnorm(dVec(ii,:))^3; % [meV]
    if 1 - exp( 2 * bond * beta ) > crit % automatically excludes the center spin
        cIndices(ii) = ii; % store the index of the site added to the cluster
    end
end
cluster = nonzeros(cIndices); % cluster indices
clusterSize = length(cluster); % cluster size

spin_new = eSpin; % copy the current spin config
spin_new(cluster, 3) = -spin_new(cluster, 3);
wav_new = coef;
Esi_new = Esi;

dipE_old = dipSum(const, ion, params.pos, eSpin);
dipE_new = dipSum(const, ion, params.pos, spin_new);
dEd = dipE_new - dipE_old; % global change of dipole energy

% flip Jz on all sites in the cluster
dEi = zeros(length(cluster),1); % local energy change
for ii = 1:length(cluster)    
    % compute change of the single-ion hamiltonian
    wav_new(:,cluster(ii)) = flip(wav_new(:, cluster(ii)));
    Esi_new(cluster(ii)) = real(wav_new(:,cluster(ii))' * hamI * wav_new(:, cluster(ii))); % single-ion energy after update
    dEi(ii) = Esi_new(cluster(ii)) - Esi(cluster(ii)); % change in single-ion energy

    % compute exchange interaction energy change
    if ion.ex(ion.idx)
        Exch_old = exchange(params, ion, eSpin, cluster(ii), eSpin(cluster(ii),:));
        Exch_new = exchange(params, ion, spin_new, cluster(ii), spin_new(cluster(ii),:));
        dEi(ii) = dEi(ii) + Exch_new - Exch_old;
    end

    % thermalize nuclear spins to the original configurations
    if params.hyp
        Ehyp_new = ion.A(ion.idx) * dot(nSpin(cluster(ii),:), spin_new(cluster(ii),:));
        Ehyp_old = ion.A(ion.idx) * dot(nSpin(cluster(ii),:), eSpin(cluster(ii),:));
        dEi(ii) = dEi(ii) + Ehyp_new - Ehyp_old; % include hyperfine energy (meV)
    end
end
dE = sum(dEi) + dEd; % global energy change

% flip the cluster without acceptance criteria
eSpin = spin_new;
coef = wav_new;
Esi = Esi_new;
accpRate = 1; % acceptance rate

% % Alternative: Glauber algorithm for flipping the cluster
% accpRate = 0;
% prob = 0;
% crit = rand; % update criteria
% if params.temp(tt) > 0
%     prob = 1 / ( 1 + exp(dE * beta) ); % probability critirion
% end
% 
% if prob > crit || dE <= 0
%     eSpin = spin_new;
%     coef = wav_new;
%     Esi = Esi_new;
%     accpRate = 1; % acceptance rate
% end
end