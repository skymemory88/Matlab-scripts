function [nList, nVecs, nDists] = neighborList_a(pos, dpRng)
% Build neighbor lists for dipole-dipole calculations
% Inputs:
%   pos - N x 3 array of positions [Angstrom]
%   dpRng - cutoff distance [Angstrom]
% Outputs:
%   nList{i} - indices of neighbors of atom i
%   nVecs{i} - displacement vectors from i to neighbors [meters]
%   nDists{i} - distances from i to neighbors [meters]

Nsite = size(pos, 1);
cutoff2 = dpRng^2; % Work with squared distances to avoid sqrt

% Initialize cell arrays
nList = cell(Nsite, 1);
nVecs = cell(Nsite, 1);
nDists = cell(Nsite, 1);

% Build lists for each atom
for ii = 1:Nsite
    % Calculate distances to all other atoms
    r_vec = pos(ii, :) - pos; % [Angstrom] - vectors FROM neighbors TO site ii (for consistency with fallback)
    r2 = sum(r_vec.^2, 2); % Squared distances
    
    % Find neighbors within cutoff (excluding self)
    mask = (r2 <= cutoff2) & (r2 > 0);
    ngbr_idx = find(mask);
    
    % Store neighbor data
    nList{ii} = ngbr_idx;
    nVecs{ii} = r_vec(ngbr_idx, :) * 1e-10; % Convert to meters
    nDists{ii} = sqrt(r2(ngbr_idx)) * 1e-10; % Convert to meters
end

% Report statistics
n_ngbr = cellfun(@length, nList);
fprintf('Neighbor list built: avg %.1f neighbors/atom (min %d, max %d)\n', ...
        mean(n_ngbr), min(n_ngbr), max(n_ngbr));
end