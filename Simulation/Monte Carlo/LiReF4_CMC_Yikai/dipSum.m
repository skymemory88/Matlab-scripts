function E_tot = dipSum(const, ion, pos, spins, cutoff)
% Calculate the total magnetic dipole-dipole interaction energy
% for all pairs in the ensemble of spins without double counting.
% Initialize total interaction energy (SI units)
gfac = const.mu0 / 4 / pi * (ion.gLande(ion.idx) * const.muB)^2;  % prefactor
E_tot = 0; % Initialize the total energy

% Check for cutoff argument
if nargin > 4
    cutoff = abs(cutoff); % Ensure cutoff is a positive scalar
else
    cutoff = 2*max(vecnorm(pos,2,2)); % set cutoff beyound simulation boundary
end

% Iterate over spin pairs and calculate interaction energy
for ii = 1:size(spins, 1)
    for jj = ii+1:size(spins, 1)
        % Check for self-interaction
        if ii == jj
            continue;
        end

        % Calculate position vector and distance between spin i and j
        r_vec = pos(ii, :) - pos(jj, :); % [Angstrom]
        r_vec = r_vec * 10^-10; % Angstrom -> meter
        r = norm(r_vec); % [m]

        % Apply cutoff if provided
        if nargin > 4 && r > cutoff
            continue; % Skip if beyond cutoff distance
        end

        % Calculate dot products
        spin_dot = dot(spins(ii, :), spins(jj, :));
        r_dot_i = dot(r_vec, spins(ii, :));
        r_dot_j = dot(r_vec, spins(jj, :));

        % Calculate interaction energy (meV) for the pair
        E_dip = gfac * (spin_dot / r^3 - 3 * (r_dot_i * r_dot_j) / r^5);

        % Add to total E_tot (J)
        E_tot = E_tot + E_dip;
    end
end

% Convert to meV
E_tot = E_tot * const.J2meV;
end