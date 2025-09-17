function E_dip = CntrDip(gfac, pos, spin, idx, spin0, cutoff)
% Calculate the magnetic dipole-dipole interaction energy between a particular spin
% indexed by 'idx' with spin orientation 'spin0' and the rest of the spins in the ensemble.
% Input:
%   spin - n x 3 matrix containing the x, y, z components of n spins
%   pos - n x 3 matrix containing the x, y, z coordinates of n spins
%   idx - index of the spin in center
%   spin0 - the spin configuration of the spin in center
%   cutoff (optional) - dipolar interaction range [Angstrom]
% Output:
%   E_dip - magnetic dipole-dipole interaction energy

% Calculate position vectors and distances
r_vec = pos - pos(idx, :); % [Angstrom]
r_vec = r_vec * 1e-10; % angstrom -> meter
r = vecnorm(r_vec, 2, 2); % [m]

% Apply isotropic cutoff
if nargin > 5
    mask = r <= cutoff * 1e-10; % angstrom -> meter
else
    mask = true(size(r));
end

% Calculate dot products
ss_dot = dot(spin, repmat(spin0, size(spin, 1), 1), 2);clc
sr_dot1 = sum(r_vec .* repmat(spin0, size(r_vec, 1), 1), 2);
sr_dot2 = sum(r_vec .* spin, 2);

% Remove the self-interaction term
ss_dot(idx) = [];  
mask(idx) = [];
r(idx) = [];
sr_dot1(idx) = [];
sr_dot2(idx) = [];

% Calculate interaction energy (meV)
E_dip = gfac * (ss_dot ./ r.^3 - 3 * (sr_dot1 .* sr_dot2) ./ r.^5);
E_dip = sum(E_dip .* mask); % apply cutoff & avoid double counting
end