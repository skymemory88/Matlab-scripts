function E_dip = CntrDip_a(params, gfac, spin, idx, spin0)
% Calculate magnetic dipole-dipole interaction energy using pre-computed neighbor lists
% Inputs:
%   gfac - dipole interaction prefactor
%   spin - N x 3 matrix containing spin orientations
%   idx - index of the central spin
%   spin0 - new spin orientation of the central spin
%   nList - cell array of neighbor indices
%   nVecs - cell array of displacement vectors [meters]
%   nDists - cell array of distances [meters]
% Output:
%   E_dip - magnetic dipole-dipole interaction energy [meV]

% Get neighbors for spin idx
nList = params.nList;
nVecs = params.nVecs;
nDists = params.nDists;

ngbr = nList{idx}; % retrive the neighbor list for spin(idx)
if isempty(ngbr)
    E_dip = 0;
    return;
end

% Get pre-computed vectors and distances
r_vec = nVecs{idx}; % Already in meters
r = nDists{idx};    % Already in meters

% Get spin orientations of neighbors
spin_ngbr = spin(ngbr, :);

% Vectorized dot products
ss_dot = spin_ngbr * spin0';
sr_dot1 = r_vec * spin0';
sr_dot2 = sum(r_vec .* spin_ngbr, 2);

% Calculate interaction energy
E_dip = gfac * sum(ss_dot ./ r.^3 - 3 * (sr_dot1 .* sr_dot2) ./ r.^5);
end