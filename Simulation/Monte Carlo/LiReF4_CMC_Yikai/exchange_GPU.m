function E_ex = exchange_GPU(pos, ion, nnD, spin, cntr, spin0)
% Compute distance vectors relative to the center spin
% Ensure inputs are on the GPU

if ~isa(pos, 'gpuArray')
    pos = gpuArray(pos);
end
if ~isa(spin, 'gpuArray')
    spin = gpuArray(spin);
end
if ~isa(spin0, 'gpuArray')
    spin0 = gpuArray(spin0);
end

origin = pos(cntr, :);  % The position of the center spin
dVec = pos - origin;   % Displacement vectors from the center spin

% Compute the Euclidean distance from the center spin
dist = vecnorm(dVec, 2, 2);

% Identify the nearest neighbours based on the unit cell distance
idx = find(dist > 0 & dist <= nnD);  % Include only those within the nearest neighbour distance and exclude the center itself

% Calculate exchange energy
if ~isempty(idx)
    spin_rep = repmat(spin0, length(idx), 1);  % Replicate spin0 to match the dimensions of spin(idx, :)
    int_sum = sum(spin_rep .* spin(idx, :), 2);  % Dot product of spins
    E_ex = ion.ex(ion.idx) * sum(int_sum);  % Sum the interactions and apply exchange strength
else
    E_ex = gpuArray(0);
end
end