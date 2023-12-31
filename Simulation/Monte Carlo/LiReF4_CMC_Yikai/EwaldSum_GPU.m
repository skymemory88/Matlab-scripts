function [totalEnergy] = EwaldSum_GPU(pos, spin, alpha, rcut, kcut)
    % Transfer data to GPU
    pos = gpuArray(pos);
    spin = gpuArray(spin);
    
    % Constants
    mu_0 = 4 * pi * 1e-7; % Permeability of vacuum in T*m/A
    const_factor = mu_0 / (4 * pi);
    
    % Initialize total energy
    totalEnergy = gpuArray(0);
    
    % Number of spins
    N = size(pos, 1);
    
    % Real-space sum
    [i_idx, j_idx] = ndgrid(1:N, 1:N);
    mask = i_idx ~= j_idx;
    
    r_ij = pos(i_idx(mask), :) - pos(j_idx(mask), :);
    r = vecnorm(r_ij, 2, 2);
    
    within_cutoff = r < rcut;
    r = r(within_cutoff);
    r_ij = r_ij(within_cutoff, :);
    
    r3 = r.^3;
    r5 = r.^5;
    erfc_term = erfc(alpha * r) ./ r3;
    exp_term = exp(-alpha^2 * r.^2) ./ r5;
    
    tensor_ij = const_factor * (2 * alpha / sqrt(pi) * exp_term + erfc_term) .* ...
                (3 * r_ij' * r_ij - r.^2 * eye(3));
    
    % Add to total energy
    energy_ij = sum(spin(i_idx(mask), :) .* tensor_ij .* spin(j_idx(mask), :), 'all');
    totalEnergy = totalEnergy + energy_ij;
    
    % Fourier-space sum
    [kx, ky, kz] = meshgrid(-kcut:kcut, -kcut:kcut, -kcut:kcut);
    k = gpuArray([kx(:), ky(:), kz(:)]);
    
    k_norm = vecnorm(k, 2, 2);
    valid_k = (k_norm > 0) & (k_norm < kcut);
    k = k(valid_k, :);
    k_norm = k_norm(valid_k);
    
    k2 = k_norm.^2;
    k4 = k2.^2;
    exp_term = exp(-k2 / (4 * alpha^2)) ./ k4;
    
    S_k = sum(spin .* exp(-1i * 2 * pi * pos * k'), 1);
    
    % Add to total energy
    totalEnergy = totalEnergy + const_factor * sum(exp_term .* abs(S_k).^2);
    
    % Transfer totalEnergy back to CPU
    totalEnergy = gather(totalEnergy);
end
