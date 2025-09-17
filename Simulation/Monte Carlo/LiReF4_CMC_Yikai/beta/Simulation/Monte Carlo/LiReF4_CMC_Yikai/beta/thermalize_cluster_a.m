function [cSize, dE, aRate, E_si, coef, spin_e] =...
    thermalize_cluster_a(ion, const, params, beta, E_si, hamI, coef, basis, spin_e, nSpin, tt, ff)
%THERMALIZE_CLUSTER_A  Long-range Wolff/LB cluster with quantum unitary transformations (LB tables if available)

% --- aliases from your suite ---
r = params.pos; % [angstrom] N×3 positions
D = const.gfac; % dipolar prefactor
N = size(r,1);


if N < 2
    cSize = 0; dE = 0; aRate = 0; return;
end

% --- Generate random unitary transformation for quantum spins ---
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

% --- init cluster state ---
i0 = randi(N);
inC = false(N,1);
stack = i0;
inC(i0) = true;
cSize = 1;


useLB = isfield(params,'LB') && ~isempty(params.LB);
if useLB
    LB = params.LB;
    scale = 2 * beta * D;     % scales the unit cumulative S_unit
    if scale <= 0
        useLB = false;        % high-T degenerate; fall back
    end
end

% Build cluster
if useLB
    % LB growth using precomputed offsets + S_unit
    Off   = LB.Offsets;      % (N-1)×3
    invr3 = LB.inv_r3;       % (N-1)×1
    Sunit = LB.S_unit;       % (N-1)×1
    box   = LB.box;
    posmap= LB.posmap;

    while ~isempty(stack)
        i = stack(end); stack(end) = [];

        z = 0;  % number of candidates already considered for this i
        while true
            u = rand();
            if z > 0
                base = Sunit(z);                             % S_unit(z)
            else
                base = 0;                                    % Start from zero for first candidate
            end
            tau_unit = (-0.5*log(1 - u))/scale + base;       % scaled target
            w = lower_bound_unit(Sunit, tau_unit);           % next candidate
            if isempty(w), break; end
            z = w;

            % map offset → absolute neighbor j of i
            rj = wrap_pos(r(i,:) + Off(w,:), box);
            key = hash_pos(rj, LB.keyDigits);
            if isKey(posmap, key)
                j = posmap(key);
            else
                j = find_site_by_position(rj, r, box);       % robust fallback
                if isempty(j), continue; end
            end
            if j==i || inC(j), continue; end

            % Luijten-Blöte bond probability for dipolar interactions:
            Jij = D * invr3(w);
            
            % Current spin expectation values
            si = spin_e(i,:); sj = spin_e(j,:);
            
            % Standard Luijten-Blöte bond probability using current configuration
            si_dot_sj = dot(si, sj);
            % For dipolar systems, use absolute value - bonds can form regardless of alignment
            % The reflection will handle the proper physics
            if abs(si_dot_sj) > 1e-6  % avoid numerical issues
                p_add = 1 - exp(-2 * beta * Jij * abs(si_dot_sj));
            else
                p_add = 0;
            end
            
            
            % Ensure probability is in valid range
            p_add = max(0, min(1, p_add));

            if rand() < p_add
                inC(j) = true;
                stack(end+1) = j; %#ok<AGROW>
                cSize = cSize + 1;
            end
        end
    end

else
    % fallback growth
    while ~isempty(stack)
        i = stack(end); stack(end) = [];

        cand = find(~inC); cand(cand==i) = [];
        if isempty(cand), continue; end

        rij = r(cand,:) - r(i,:);
        if isfield(params,'box')
            rij = mi_wrap_rows(rij, params.box);
        end
        rnorm = sqrt(sum(rij.^2,2));
        keep = rnorm > 0;
        cand = cand(keep); rij = rij(keep,:); rnorm = rnorm(keep);
        if isempty(cand), continue; end

        Jiso = D ./ (rnorm.^3);
        [Jiso, ord] = sort(Jiso, 'descend');
        cand = cand(ord);

        S0 = cumsum(2*beta*Jiso);

        z = 0;
        while true
            u = rand();
            if z > 0
                base = S0(z);
            else
                base = 0;
            end
            tau  = -0.5*log(1 - u) + base;
            w    = lower_bound(S0, tau);
            if isempty(w), break; end
            z = w;

            j = cand(w);
            
            % Current spin expectation values
            si = spin_e(i,:); sj = spin_e(j,:);
            
            % Standard Luijten-Blöte bond probability using current configuration
            si_dot_sj = dot(si, sj);
            % For dipolar systems, use absolute value - bonds can form regardless of alignment
            % The reflection will handle the proper physics
            if abs(si_dot_sj) > 1e-6  % avoid numerical issues
                p_add = 1 - exp(-2 * beta * Jiso(w) * abs(si_dot_sj));
            else
                p_add = 0;
            end
            
            
            % Ensure probability is in valid range
            p_add = max(0, min(1, p_add));

            if rand() < p_add
                inC(j) = true;
                stack(end+1) = j;
                cSize = cSize + 1;
            end
        end
    end
end

% Calculate energy change for the entire cluster using quantum approach
% Apply unitary transformation to all sites in the cluster

cluster = find(inC);
cSize = length(cluster);

% Create new configuration
coef_new = coef;
spin_e_new = spin_e;
E_si_new = E_si;

% Transform all wavefunctions in the cluster
for ii = 1:cSize
    site = cluster(ii);
    
    % Apply unitary transformation to wavefunction
    coef_new(:, site) = U * coef(:, site);
    
    % Compute new spin expectation values
    wav_new = coef_new(:, site);
    spx = real(wav_new' * basis' * ion.Jx * basis * wav_new);
    spy = real(wav_new' * basis' * ion.Jy * basis * wav_new);
    spz = real(wav_new' * basis' * ion.Jz * basis * wav_new);
    spin_e_new(site, :) = [spx, spy, spz];
    
    % New single-ion energy
    E_si_new(site) = real(wav_new' * hamI * wav_new);
    
    % Add hyperfine contribution if applicable
    if params.hyp && ismember(site, params.isoIdx)
        E_si_new(site) = E_si_new(site) + ion.A(ion.idx) * dot(nSpin(site,:), spin_e_new(site,:));
    end
end

% Calculate total energy change
% Single-ion energy change
dE_si = sum(E_si_new(cluster)) - sum(E_si(cluster));

% Dipole interaction energy change using GPU with trial configuration
if params.useGPU && ~isempty(params.gpuState)
    % Current configuration energy
    E_dip_old = dipSum_gpu_a(params.gpuState, tt, ff);
    
    % Trial configuration energy - pass new spins as 4th argument
    E_dip_new = dipSum_gpu_a(params.gpuState, tt, ff, spin_e_new);
    
    dE_dip = E_dip_new - E_dip_old;
else
    % CPU fallback using dipole sum function
    E_dip_old = dipSum_a(params, D, spin_e);
    E_dip_new = dipSum_a(params, D, spin_e_new);
    dE_dip = E_dip_new - E_dip_old;
end

% Exchange interaction energy change (if applicable)
dE_ex = 0;
if ion.ex(ion.idx) > 0
    for ii = 1:cSize
        site = cluster(ii);
        neighbors = params.nList{site};
        
        for jj = 1:length(neighbors)
            neighbor = neighbors(jj);
            % Only count each pair once
            if neighbor > site
                % Old configuration
                E_ex_old = ion.ex(ion.idx) * dot(spin_e(site,:), spin_e(neighbor,:));
                % New configuration
                E_ex_new = ion.ex(ion.idx) * dot(spin_e_new(site,:), spin_e_new(neighbor,:));
                dE_ex = dE_ex + (E_ex_new - E_ex_old);
            end
        end
    end
end

% Total energy change
dE = dE_si + dE_dip + dE_ex;

% Accept or reject cluster flip (Metropolis criterion)
prob = exp(-beta * dE);
crit = rand();

if prob >= crit || dE <= 0
    % Accept the cluster flip
    coef = coef_new;
    spin_e = spin_e_new;
    E_si = E_si_new;
    aRate = 1;
    
    % Update GPU state if using GPU
    if params.useGPU && ~isempty(params.gpuState)
        params.gpuState.updateSpins(spin_e, tt, ff);
    end
else
    % Reject the cluster flip
    aRate = 0;
    dE = 0; % No energy change
end

% Optional: Print diagnostics
if isfield(params, 'verbose') && params.verbose && cSize > 1
    fprintf('Cluster size: %d, ΔE: %.4f meV, Accepted: %d\n', cSize, dE, aRate);
end
end

%% helper functions
function w = lower_bound(S0, tau)
    if tau >= S0(end), w = []; return; end
    lo=1; hi=numel(S0);
    while lo<hi
        mid = floor((lo+hi)/2);
        if S0(mid) > tau, hi=mid; else, lo=mid+1; end
    end
    w = lo;
end

function w = lower_bound_unit(Sunit, tau_unit)
    if tau_unit >= Sunit(end), w = []; return; end
    lo=1; hi=numel(Sunit);
    while lo<hi
        mid = floor((lo+hi)/2);
        if Sunit(mid) > tau_unit, hi=mid; else, lo=mid+1; end
    end
    w = lo;
end


function rij = mi_wrap(dr, box)
    rij = dr;
    for a = 1:3
        L = box(a);
        if L > 0
            rij(a) = rij(a) - L*round(rij(a)/L);
        end
    end
end

function rij = mi_wrap_rows(dr, box)
    rij = dr;
    for k = 1:size(dr,1)
        rij(k,:) = mi_wrap(dr(k,:), box);
    end
end

function rp = wrap_pos(r, box)
    rp = r;
    for a=1:3
        L = box(a);
        if L>0, rp(a) = r(a) - L*floor((r(a)/L) + 0.5); end
    end
end

function key = hash_pos(ri, keyDigits)
    fmt = sprintf('%%.%df,%%.%df,%%.%df', keyDigits, keyDigits, keyDigits);
    key = sprintf(fmt, ri(1), ri(2), ri(3));
end

function j = find_site_by_position(rj, r, box)
% Robust O(N) mapping by nearest wrapped position (fallback)
    N = size(r,1);
    best = 0; bestd2 = inf;
    for k = 1:N
        dk = mi_wrap(rj - r(k,:), box);
        d2 = dot(dk,dk);
        if d2 < bestd2
            bestd2 = d2; best = k;
        end
    end
    if bestd2 < 1e-10, j = best; else, j = []; end
end