function E_dip = CntrDip_GPU(gfac, pos, spin, idx, spin0, cutoff)
% Ensure inputs are on the GPU
if ~isa(pos, 'gpuArray'); pos = gpuArray(pos); end
if ~isa(spin, 'gpuArray'); spin = gpuArray(spin); end
if ~isa(spin0, 'gpuArray'); spin0 = gpuArray(spin0);end

% Calculate position vectors and distances
r_vec = (pos - pos(idx, :)) * 1e-10; % angstrom -> meter
% r = vecnorm(r_vec, 2, 2); % [m]
dist = sqrt(sum(r_vec .^ 2, 2)); % [m]


% Apply isotropic cutoff
if nargin > 6
    mask = dist <= cutoff;
else
    mask = true(size(dist));
end

% Calculate dot products using repmat for implicit expansion
spin0_exd = repmat(spin0, size(spin, 1), 1);
spin_dot = dot(spin, spin0_exd, 2);
rm_dot1 = sum(r_vec .* spin0_exd, 2);
rm_dot2 = sum(r_vec .* spin, 2);

% Remove the self-interaction term
spin_dot(idx) = [];
mask(idx) = [];
dist(idx) = [];
rm_dot1(idx) = [];
rm_dot2(idx) = [];

% Calculate interaction energy (meV)
E_dip = gfac * (spin_dot ./ dist.^3 - 3 * (rm_dot1 .* rm_dot2) ./ dist.^5);
E_dip = sum(E_dip .* mask); % apply cutoff criteria and avoid double counting
end
