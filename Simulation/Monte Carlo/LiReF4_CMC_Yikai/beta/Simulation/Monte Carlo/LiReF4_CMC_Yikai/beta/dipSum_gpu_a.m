function E_dip = dipSum_gpu_a(gpuState, tempIdx, fieldIdx, trial_spins)
% GPU implementation with improved accuracy using Kahan summation
% If trial_spins provided, use them instead of stored spins
if nargin > 3 && ~isempty(trial_spins)
    d_spins = gpuArray(trial_spins);
else
    d_spins = gpuState.d_eSpin{tempIdx, fieldIdx};
end

% Pre-build pair list if not already done
if isempty(gpuState.pairData)
    gpuState = buildGPUPairData(gpuState);
end

% Get pre-built pair data
pairData = gpuState.pairData;
d_i_indices = pairData.d_i_indices;
d_j_indices = pairData.d_j_indices;
d_r_all = pairData.d_r_all;
d_rvec_all = pairData.d_rvec_all;

% Extract spin pairs
d_spin_i = d_spins(d_i_indices, :);
d_spin_j = d_spins(d_j_indices, :);

% Vectorized computation
d_ss_dot = sum(d_spin_i .* d_spin_j, 2);
d_sr_dot_i = sum(d_rvec_all .* d_spin_i, 2);
d_sr_dot_j = sum(d_rvec_all .* d_spin_j, 2);

% Energy calculation - use higher precision intermediate
d_r2 = d_r_all .* d_r_all;
d_r3 = d_r2 .* d_r_all;
d_r5 = d_r3 .* d_r2;

d_term1 = d_ss_dot ./ d_r3;
d_term2 = 3 * (d_sr_dot_i .* d_sr_dot_j) ./ d_r5;
d_E_all = d_term1 - d_term2;

% Use double precision for final sum if not already
if ~isa(d_E_all, 'double')
    d_E_all = double(d_E_all);
end

% Improved summation for large arrays
% Split into chunks to reduce accumulation error
chunk_size = length(d_E_all); % sum over in one chunk
% chunk_size = 10000;
n_chunks = ceil(length(d_E_all) / chunk_size);

E_sum = 0;
for i = 1:n_chunks
    start_idx = (i-1) * chunk_size + 1;
    end_idx = min(i * chunk_size, length(d_E_all));
    chunk_sum = sum(d_E_all(start_idx:end_idx));
    E_sum = E_sum + gather(chunk_sum);
end

E_dip = E_sum * gather(gpuState.d_gfac);
end

function gpuState = buildGPUPairData(gpuState)
% Pre-build pair data structure for efficient GPU computation
% This is called once and reused for all subsequent calculations

N = gpuState.N;

% First pass: count total number of pairs to pre-allocate
total_pairs = 0;
for i = 1:N
    if ~isempty(gpuState.d_nList{i})
        neighbors = gather(gpuState.d_nList{i});
        total_pairs = total_pairs + sum(neighbors > i);
    end
end

% check point
% fprintf('Building GPU pair data structure for %d interaction pairs...\n', total_pairs);

% Pre-allocate CPU arrays
i_indices = zeros(total_pairs, 1, 'int32');
j_indices = zeros(total_pairs, 1, 'int32');
r_all = zeros(total_pairs, 1, 'double');
rvec_all = zeros(total_pairs, 3, 'double');

% Second pass: fill arrays
pair_idx = 1;
for i = 1:N
    if ~isempty(gpuState.d_nList{i})
        % Get neighbors - need to gather for logical indexing
        neighbors = gather(gpuState.d_nList{i});
        mask = neighbors > i;

        if any(mask)
            valid_neighbors = neighbors(mask);
            n_valid = length(valid_neighbors);

            % Fill index arrays
            i_indices(pair_idx:pair_idx+n_valid-1) = i;
            j_indices(pair_idx:pair_idx+n_valid-1) = valid_neighbors;

            % Get distances and vectors - gather once
            r_all(pair_idx:pair_idx+n_valid-1) = gather(gpuState.d_nDists{i}(mask));
            rvec_all(pair_idx:pair_idx+n_valid-1, :) = gather(gpuState.d_nVecs{i}(mask, :));

            pair_idx = pair_idx + n_valid;
        end
    end
end

% Create pairData struct
pairData = struct();

% Transfer to GPU once - this is the only CPU->GPU transfer
pairData.d_i_indices = gpuArray(i_indices);
pairData.d_j_indices = gpuArray(j_indices);
pairData.d_r_all = gpuArray(r_all);
pairData.d_rvec_all = gpuArray(rvec_all);
pairData.n_pairs = total_pairs;

% Store in gpuState
gpuState.pairData = pairData;

% checkpoint
% fprintf('GPU pair data ready. Memory usage: %.2f MB\n', ...
% (total_pairs * (2*4 + 1*8 + 3*8)) / 1024^2);
end