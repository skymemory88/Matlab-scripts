function E_dip = CntrDip(const, ion, pos, spin, idx, spin0, cutoff)
% Calculate the magnetic dipole-dipole interaction energy between a particular spin
% indexed by 'idx' with spin orientation 'spin0' and the rest of the spins in the ensemble.
% Input:
%   spin - n x 3 matrix containing the x, y, z components of n spins
%   pos - n x 3 matrix containing the x, y, z coordinates of n spins
%   idx - index of the spin in center
%   spin0 - the spin configuration of the spin in center
%   cutoff (optional) - dipolar interaction range cutoff
% Output:
%   E_dip - magnetic dipole-dipole interaction energy

% Initialize total interaction energy (SI units)
gfac = const.mu0/4/pi * (ion.gLande(ion.idx) * const.muB)^2; % prefactor, 10^30 to conver Angstron to m

% Calculate position vectors and distances
r_vec = pos - pos(idx, :); % [Angstrom]
r_vec = r_vec * 10^-10; % angstrom -> meter
r = vecnorm(r_vec, 2, 2); % [m]

% Apply isotropic cutoff
if nargin > 6
    mask = r <= cutoff;
else
    mask = true(size(r));
end

% Calculate dot products
spin_dot = dot(spin, repmat(spin0, size(spin, 1), 1), 2);
rm_dot1 = sum(r_vec .* repmat(spin0, size(r_vec, 1), 1), 2);
rm_dot2 = sum(r_vec .* spin, 2);

% Remove the self-interaction term
spin_dot(idx) = [];  
mask(idx) = [];
r(idx) = [];
rm_dot1(idx) = [];
rm_dot2(idx) = [];

% Calculate interaction energy (meV)
E_dip = gfac * (spin_dot ./ r.^3 - 3 * (rm_dot1 .* rm_dot2) ./ r.^5);
E_dip = sum(E_dip .* mask) * const.J2meV; % apply cutoff criteria

end