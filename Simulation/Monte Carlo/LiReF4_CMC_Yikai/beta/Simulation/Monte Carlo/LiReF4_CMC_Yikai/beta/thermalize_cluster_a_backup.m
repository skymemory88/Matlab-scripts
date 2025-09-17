function [clusterSize, dE, accpRate, Esi, coef, eSpin] = thermalize_cluster_a(ion, const, params, beta, Esi, hamI, coef, basis, eSpin, nSpin, tt, ff)
% Cluster Monte Carlo algorithm for quantum spins using unitary transformations
% Works with wavefunctions on the Bloch sphere rather than classical spin reflections

% Step 1: Generate random unitary transformation
% For spin-J system, we need a (2J+1) x (2J+1) unitary matrix
dim = size(hamI, 1); % Hilbert space dimension = 2J+1

% Generate random unitary matrix using QR decomposition of random complex matrix
% This gives a uniformly distributed unitary matrix (Haar measure)
Z = randn(dim, dim) + 1i*randn(dim, dim);
[Q, R] = qr(Z);
% Ensure proper phase convention
for j = 1:dim
    if real(R(j,j)) < 0
        Q(:,j) = -Q(:,j);
    end
end
U = Q; % Random unitary transformation

% Step 2: Initialize cluster construction
N = size(params.pos,1);
cntr = randi(N); % Randomly pick a seed site
cluster = cntr; % Initialize cluster with seed
stack = cntr; % Stack for sites to process
inCluster = false(N,1); % Track which sites are in cluster
inCluster(cntr) = true;

% Get dipole prefactor
gfac = const.J2meV * const.mu0 / 4 / pi * (ion.gLande(ion.idx) * const.muB)^2;

% Step 3: Build cluster using modified Wolff algorithm for quantum spins
% For quantum spins, we need to compute the change in expectation values
% when applying the unitary transformation

% Pre-compute the transformed basis for efficiency
basis_U = basis * U; % Transformed basis states

while ~isempty(stack)
    % Get current site from stack
    current = stack(1);
    stack(1) = [];

    % Get neighbors of current site
    neighbors = params.nList{current};
    if isempty(neighbors)
        continue;
    end

    % Get distances to neighbors
    r = params.nDists{current}; % Already in meters

    % Current wavefunction and its transformation
    wav_current = coef(:, current);
    wav_current_new = U * wav_current;

    % Compute current spin and transformed spin
    spx_old = real(wav_current' * basis' * ion.Jx * basis * wav_current);
    spy_old = real(wav_current' * basis' * ion.Jy * basis * wav_current);
    spz_old = real(wav_current' * basis' * ion.Jz * basis * wav_current);

    spx_new = real(wav_current_new' * basis_U' * ion.Jx * basis_U * wav_current_new);
    spy_new = real(wav_current_new' * basis_U' * ion.Jy * basis_U * wav_current_new);
    spz_new = real(wav_current_new' * basis_U' * ion.Jz * basis_U * wav_current_new);

    s_current = [spx_old, spy_old, spz_old];
    s_current_new = [spx_new, spy_new, spz_new];
    delta_s = s_current_new - s_current;

    % Check each neighbor for inclusion
    for jj = 1:length(neighbors)
        neighbor = neighbors(jj);

        % Skip if already in cluster
        if inCluster(neighbor)
            continue;
        end

        % Get neighbor spin
        s_neighbor = eSpin(neighbor,:);

        % Calculate coupling strength (isotropic part only for cluster building)
        J_ij = gfac / r(jj)^3;

        % Bond probability for quantum spins
        % P = 1 - exp[-β * J_ij * Δs_i · s_j]
        bond_energy = J_ij * dot(delta_s, s_neighbor);

        % Compute probability
        if bond_energy > 0
            prob = 1 - exp(-beta * bond_energy);
        else
            prob = 0; % Don't add if energy would increase
        end

        % Add to cluster with probability
        if rand() < prob
            cluster = [cluster, neighbor];
            stack = [stack, neighbor];
            inCluster(neighbor) = true;
        end
    end
end

clusterSize = length(cluster);

% Step 4: Calculate energy change for the entire cluster
% Apply unitary transformation to all sites in the cluster

% Create new configuration
coef_new = coef;
eSpin_new = eSpin;
Esi_new = Esi;

% Transform all wavefunctions in the cluster
for ii = 1:clusterSize
    site = cluster(ii);

    % Apply unitary transformation to wavefunction
    coef_new(:, site) = U * coef(:, site);

    % Compute new spin expectation values
    wav_new = coef_new(:, site);
    spx = real(wav_new' * basis' * ion.Jx * basis * wav_new);
    spy = real(wav_new' * basis' * ion.Jy * basis * wav_new);
    spz = real(wav_new' * basis' * ion.Jz * basis * wav_new);
    eSpin_new(site, :) = [spx, spy, spz];

    % New single-ion energy
    Esi_new(site) = real(wav_new' * hamI * wav_new);

    % Add hyperfine contribution if applicable
    if params.hyp && ismember(site, params.isoIdx)
        Esi_new(site) = Esi_new(site) + ion.A(ion.idx) * dot(nSpin(site,:), eSpin_new(site,:));
    end
end

% Step 5: Calculate total energy change
% 5a: Single-ion energy change
dE_si = sum(Esi_new(cluster)) - sum(Esi(cluster));

% 5b: Dipole interaction energy change using GPU with trial configuration
if params.useGPU && ~isempty(params.gpuState)
    % Current configuration energy
    E_dip_old = dipSum_gpu_a(params.gpuState, tt, ff);

    % Trial configuration energy - pass new spins as 4th argument
    E_dip_new = dipSum_gpu_a(params.gpuState, tt, ff, eSpin_new);

    dE_dip = E_dip_new - E_dip_old;
else
    % CPU fallback
    E_dip_old = dipSum_a(params, gfac, eSpin);
    E_dip_new = dipSum_a(params, gfac, eSpin_new);
    dE_dip = E_dip_new - E_dip_old;
end

% 5c: Exchange interaction energy change (if applicable)
dE_ex = 0;
if ion.ex(ion.idx) > 0
    for ii = 1:clusterSize
        site = cluster(ii);
        neighbors = params.nList{site};

        for jj = 1:length(neighbors)
            neighbor = neighbors(jj);
            % Only count each pair once
            if neighbor > site
                % Old configuration
                E_ex_old = ion.ex(ion.idx) * dot(eSpin(site,:), eSpin(neighbor,:));
                % New configuration
                E_ex_new = ion.ex(ion.idx) * dot(eSpin_new(site,:), eSpin_new(neighbor,:));
                dE_ex = dE_ex + (E_ex_new - E_ex_old);
            end
        end
    end
end

% Total energy change
dE = dE_si + dE_dip + dE_ex;

% Step 6: Accept or reject cluster flip (Metropolis criterion)
prob = exp(-beta * dE);
crit = rand();

if prob >= crit || dE <= 0
    % Accept the cluster flip
    coef = coef_new;
    eSpin = eSpin_new;
    Esi = Esi_new;
    accpRate = 1;

    % Update GPU state if using GPU
    if params.useGPU && ~isempty(params.gpuState)
        params.gpuState.updateSpins(eSpin, tt, ff);
    end
else
    % Reject the cluster flip
    accpRate = 0;
    dE = 0; % No energy change
end

% Optional: Print diagnostics
if isfield(params, 'verbose') && params.verbose && clusterSize > 1
    fprintf('Cluster size: %d, ΔE: %.4f meV, Accepted: %d\n', clusterSize, dE, accpRate);
end

end