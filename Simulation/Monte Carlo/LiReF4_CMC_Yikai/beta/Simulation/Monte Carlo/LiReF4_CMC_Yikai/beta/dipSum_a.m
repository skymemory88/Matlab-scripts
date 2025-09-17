function E_dip = dipSum_a(params, gfac, spins)
% Calculate the total magnetic dipole-dipole interaction energy
% Can use either direct calculation or neighbor lists
%
% Inputs:
%   gfac - dipole interaction prefactor
%   pos - N x 3 positions [Angstrom] (for compatibility)
%   spins - N x 3 spin orientations
%   cutoff - cutoff distance [Angstrom]
%   nList - (optional) cell array of neighbor indices
%   nDists - (optional) cell array of distances [meters]

% Get neighbors for spin idx
nList = params.nList;
nVecs = params.nVecs;
nDists = params.nDists;
cutoff = params.dpRng;

spins = squeeze(spins);
N = size(spins, 1);

% Check if neighbor lists are provided
if nargin >= 6 && ~isempty(nList) && ~isempty(nDists)
    % Use neighbor list version
    E_dip = 0;
    
    % Need nVecs too!
    if ~exist('nVecs', 'var')
        error('dipSum_a: nVecs must be provided with nList and nDists');
    end
    
    % Loop over all spins
    for i = 1:N
        ngbr = nList{i};
        if isempty(ngbr)
            continue;
        end
        
        % Only count pairs where j > i to avoid double counting
        ngbr_upper = ngbr(ngbr > i);
        if isempty(ngbr_upper)
            continue;
        end
        
        % Get indices in the neighbor list for these upper neighbors
        idx_upper = ismember(ngbr, ngbr_upper);
        
        % Get distances and vectors for these neighbors
        r = nDists{i}(idx_upper);
        r_vec = nVecs{i}(idx_upper, :);
        
        % Get spin orientations
        spin_i = spins(i, :);
        spin_j = spins(ngbr_upper, :);
        
        % Dot products
        ss_dot = spin_j * spin_i';
        sr_dot_i = r_vec * spin_i';
        sr_dot_j = sum(r_vec .* spin_j, 2);
        
        % Add contribution
        E_dip = E_dip + gfac * sum(ss_dot ./ r.^3 - 3 * (sr_dot_i .* sr_dot_j) ./ r.^5);
    end
else
    % Original implementation (fallback)
    pos = params.pos * 1e-10; % Angstrom -> meter

    % Set cutoff if not provided
    if nargin < 4
        cutoff = 2 * max(vecnorm(pos, 2, 2));
    else
        cutoff = abs(cutoff) * 1e-10; % Convert to meters
    end

    % Create upper triangular indices
    [j_idx, i_idx] = meshgrid(1:N, 1:N);
    upper_mask = i_idx(:) < j_idx(:);
    i_pairs = i_idx(upper_mask);
    j_pairs = j_idx(upper_mask);

    % Vectorized calculation
    r_vec = pos(i_pairs, :) - pos(j_pairs, :);
    r = vecnorm(r_vec, 2, 2);

    % Apply cutoff
    mask = r <= cutoff & r > 0;

    if any(mask)
        r_valid = r(mask);
        r_vec_valid = r_vec(mask, :);
        spin_i = spins(i_pairs(mask), :);
        spin_j = spins(j_pairs(mask), :);

        % Dot products
        ss_dot = sum(spin_i .* spin_j, 2);
        sr_dot_i = sum(r_vec_valid .* spin_i, 2);
        sr_dot_j = sum(r_vec_valid .* spin_j, 2);

        % Total energy
        E_dip = gfac * sum(ss_dot ./ r_valid.^3 - 3 * (sr_dot_i .* sr_dot_j) ./ r_valid.^5);
    else
        E_dip = 0;
    end
end
end