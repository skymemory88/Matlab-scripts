function [Esi, dE, coef, eSpin, accpRate] = thermalize_GPU(const, ion, params, Esi, hamI, coef, basis, eSpin, nSpin, beta)
numSites = size(params.pos,1);
change = gpuArray.zeros(numSites,1);
site = gpuArray.randperm(numSites); % Random update of the whole lattice on GPU
nnD = gpuArray(min(nonzeros(vecnorm(ion.tau * ion.abc{ion.idx}, 2, 2))));  % Nearest neighbour distance for exchange interaction
Jex = gpuArray(ion.ex(ion.idx));
gfac = gpuArray(const.gfac);
pos = gpuArray(params.pos);

% Random rotations on the Bloch sphere, vectorized
coords_new = randSph(numSites, 'Cartesian');
[alp_new, bet_new, ~] = cart2sph(coords_new(:,1), coords_new(:,2), coords_new(:,3));
alp_new = alp_new + pi; % Adjust phase
bet_new = bet_new + pi/2; % Adjust elevation
coef_new = [cos(bet_new'/2); sin(bet_new'/2) .* exp(1i * alp_new')];
wav = gpuArray(coef_new);

% Vectorized computation for spin components and energy
eSpin_temp = arrayfun(@(ii) newSpin_GPU(ion, basis, wav(:,ii), 'electron'), 1:numSites, 'UniformOutput', false);
spin_new = gpuArray(cat(1, eSpin_temp{:})); % concatenate and flatten the results
Esi_new = real(sum(conj(wav) .* (hamI * wav), 1))'; % Vectorized energy calculation

if params.hyp
    % Hyperfine interactions, conditional on sites
    isoSite = ismember(site, params.isoIdx);
    Esi_new(isoSite) = Esi_new(isoSite) + ion.A(ion.idx) .* dot(nSpin(isoSite, :), spin_new(isoSite, :), 2);
end

% Energy changes
dE_si = Esi_new - Esi(site); 
Eint_old = arrayfun(@(ii) CntrDip_GPU(gfac, pos, eSpin, site(ii), eSpin(ii, :)), 1:numSites);
Eint_new = arrayfun(@(ii) CntrDip_GPU(gfac, pos, eSpin, site(ii), spin_new(ii, :)), 1:numSites);
if ion.ex(ion.idx) % check exchange interaction strength
    Eint_old = Eint_old + arrayfun(@(jj) exchange_GPU(pos, Jex, nnD, eSpin, site(jj), eSpin(jj, :)), 1:numSites);
    Eint_new = Eint_new + arrayfun(@(jj) exchange_GPU(pos, Jex, nnD, eSpin, site(jj), spin_new(jj, :)), 1:numSites);
end
dE_int = Eint_new - Eint_old;
dEi = dE_si + dE_int';

% Glauber dynamics in a vectorized form
prob = 1 ./ (1 + exp(dEi * beta));
updates = prob >= rand("gpuArray");

% Apply updates
Esi(site(updates)) = Esi_new(updates);
coef(:, site(updates)) = wav(:, updates);
eSpin(site(updates), :) = spin_new(updates, :);
change(updates) = 1;

% Compute outputs
dE = sum(dEi .* change);
accpRate = sum(change) / numSites;
end