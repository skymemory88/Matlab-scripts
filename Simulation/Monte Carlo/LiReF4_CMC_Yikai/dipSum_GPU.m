function E_tot = dipSum_GPU(gfac, pos, spins, cutoff)
% Calculate the total magnetic dipole-dipole interaction energy
% for all pairs in the ensemble of spins without double counting on GPU.

% Convert input arrays to GPU arrays if they are not already
if ~isa(pos, 'gpuArray'); pos = gpuArray(pos); end
if ~isa(spins, 'gpuArray'); spins = gpuArray(spins); end

% Calculate all pairwise distances and difference vectors
[pos_i, pos_j] = ndgrid(1:size(pos, 1), 1:size(pos, 1));
r_vecs = (pos(pos_i, :) - pos(pos_j, :)) * 10^-10; % [Angstrom --> meter]
r_ij = vecnorm(r_vecs,2,2); % [m]

ss_dots = sum(spins(pos_i, :) .* spins(pos_j, :), 3);
sr_dot_is = sum(r_vecs .* spins(pos_i, :), 3);
sr_dot_js = sum(r_vecs .* spins(pos_j, :), 3);

% Apply the cutoff mask
if nargin > 4
    cutoff_mask = (r_ij <= cutoff) & (r_ij > 0);
else
    cutoff_mask = r_ij > 0; % avoid self interaction
end
r_ij = r_ij(cutoff_mask);
rij3 = r_ij .^ 3;
rij5 = r_ij .^ 5;

ss_dots = ss_dots(cutoff_mask);
sr_dot_is = sr_dot_is(cutoff_mask);
sr_dot_js = sr_dot_js(cutoff_mask);

% Calculate interaction energy (meV) for all pairs
E_dips = gfac * (ss_dots ./ rij3 - 3 * (sr_dot_is .* sr_dot_js) ./ rij5);

% Sum up all interaction energies
E_tot = gather(sum(E_dips(:))); % [meV]
end
