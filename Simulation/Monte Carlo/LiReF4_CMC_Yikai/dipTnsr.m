function [E_tot, Dij] = dipTnsr(pos, spin)
    % Constants
    mu_0 = 4 * pi * 1e-7; % Permeability of vacuum in T*m/A
    
    % Initialize total energy
    E_tot = 0;
    
    % Number of spins
    N = size(pos, 1);
    
    % Initialize interaction tensors
    Dij = zeros(3, 3, N, N);
    
    % Calculate pairwise distances and unit vectors
    r_vec = pos - permute(pos, [3, 2, 1]);
    r = sqrt(sum(r_vec.^2, 2));
    r = reshape(r, [N, N]);
    r_hat = r_vec ./ r;
    
    % Avoid division by zero for i=j
    r(r==0) = Inf;
    
    % Compute interaction energy and tensors
    for i = 1:N
        for j = 1:N
            if i ~= j
                r_ij = r(i, j);
                r_hat_ij = squeeze(r_hat(i, j, :));
                
                % Compute tensor for this pair
                tensor_ij = mu_0 / (4 * pi) * (3 * r_hat_ij * r_hat_ij' - eye(3)) / r_ij^3;
                
                % Add to total energy
                energy_ij = spin(i, :) * tensor_ij * spin(j, :)';
                E_tot = E_tot + energy_ij;
                
                % Store tensor
                Dij(:, :, i, j) = tensor_ij;
            end
        end
    end
end