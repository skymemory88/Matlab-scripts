function E_ex = exchange(params, ion, spin, cntr, spin0)

% Option 1: rank the distance from the center spin and 
origin = params.pos(cntr,:); % set the origin to the center spin
dVec = params.pos - origin; % displacement vector to the center spin
Nd = min( nonzeros( vecnorm(ion.tau * ion.abc{ion.idx},2,2) ) ); % distance to the nearest neighbour in the unit cell
dist = vecnorm(dVec,2,2);
idx = find(dist <= Nd);
idx = idx(idx ~= cntr); % exclude self-interaction

% % options 2: include the nearest N neighbours for exchange interaction
% nNeighbour = 4;
% dVec = params.pos - origin; % displacement vector to the center spin
% dist = vecnorm(dVec,2,2);
% [~, idx] = sort(dist);
% idx = idx(1:nNeighbour);

E_ex = sum(ion.ex(ion.idx) .* dot(repmat(spin0, length(idx), 1), spin(idx, :), 2)); % (meV) "/2" to avoid double counting
end