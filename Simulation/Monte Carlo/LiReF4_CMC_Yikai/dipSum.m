function E_dip = dipSum(gfac, pos, spins, cutoff)
% Calculate the total magnetic dipole-dipole interaction energy
% for all pairs in the ensemble of spins without double counting.
% Initialize total interaction energy (SI units)
E_dip = 0; % Initialize the total energy
spins = squeeze(spins);
pos = pos * 10^-10; % Angstrom -> meter

% Check for cutoff argument
if nargin > 4
    cutoff = abs(cutoff); % Ensure cutoff is a positive scalar
else
    cutoff = 2*max(vecnorm(pos,2,2)); % set cutoff beyound simulation boundary
end

% % Calculate all pairwise distances
% [pos_i, pos_j] = ndgrid(1:numSpins, 1:numSpins);
% r_vecs = (pos(pos_i, :) - pos(pos_j, :)) * scaling_factor; % [m]
% r_ij = norm(r_vecs); % [m]

% Iterate over spin pairs and calculate interaction energy
for ii = 1:size(spins, 1)
    for jj = ii+1:size(spins, 1)
        % Check for self-interaction
        if ii == jj
            continue;
        end

        % Calculate position vector and distance between spin i and j
        r_vec = pos(ii, :) - pos(jj, :); % [Angstrom]
        r_ij = norm(r_vec); % [m]

        % Apply cutoff if provided
        if nargin > 4 && r_ij > cutoff
            continue; % Skip if beyond cutoff distance
        end

        % Calculate dot products
        ss_dot = dot(spins(ii, :), spins(jj, :));
        sr_dot_i = dot(r_vec, spins(ii, :));
        sr_dot_j = dot(r_vec, spins(jj, :));

        % Calculate interaction energy (meV) for the pair
        E_ij = gfac * (ss_dot / r_ij^3 - 3 * (sr_dot_i * sr_dot_j) / r_ij^5);

        % Add to total E_tot
        E_dip = E_dip + E_ij;
    end
end
end