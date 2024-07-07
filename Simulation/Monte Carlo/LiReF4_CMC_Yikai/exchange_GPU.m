function E_ex = exchange_GPU(pos, Jex, nnD, spin, cntr, spin0)
% Ensure all inputs are on the GPU
if ~isa(pos, 'gpuArray')
    pos = gpuArray(pos);
end
if ~isa(spin, 'gpuArray')
    spin = gpuArray(spin);
end
if ~isa(spin0, 'gpuArray')
    spin0 = gpuArray(spin0);
end

% Compute the position of the center spin
origin = pos(cntr, :);

% Calculate displacement vectors and distances
dVec = pos - origin; % Displacement vectors
dist = sqrt(sum(dVec .^ 2, 2)); % Euclidean distance

% Option 1:
% Identify the nearest neighbours based on the unit cell distance
idx = (dist > 0) & (dist <= nnD); % Logical indexing

% Calculate exchange energy
if any(idx)
    int_sum = sum(bsxfun(@times, spin(idx, :), spin0), 2); % Dot product
    E_ex = Jex * sum(int_sum); % Sum interactions
else
    E_ex = gpuArray(0);
end

% Option 2:
% Identify the nearest neighbours based on the unit cell distance
% idx = find(dist > 0 & dist <= nnD);  % Include only those within the nearest neighbour distance and exclude the center itself

% % Calculate exchange energy
% if ~isempty(idx)
%     % spin_rep = repmat(spin0, length(idx), 1);  % Replicate spin0 to match the dimensions of spin(idx, :)
%     % int_sum = sum(spin_rep .* spin(idx, :), 2);  % Dot product of spins
%     E_ex = ion.ex(ion.idx) * sum(int_sum);  % Sum the interactions and apply exchange strength
% else
%     E_ex = gpuArray(0);
% end
end