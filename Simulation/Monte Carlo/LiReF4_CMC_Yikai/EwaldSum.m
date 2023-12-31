function E_tot = EwaldSum(pos, dip, boxSize, alpha, kMax)
    % Constants
    eps0 = 8.8542e-12; % Permittivity of vacuum in F/m
    
    % Initialize total energy
    E_tot = 0;
    
    % Real-space sum
    N = size(pos, 1);
    r_vec = pos - permute(pos, [3, 2, 1]);
    r_vec = r_vec - round(r_vec ./ boxSize) .* boxSize;
    r = sqrt(sum(r_vec.^2, 2));
    r = reshape(r, [N, N]);
    r(r==0) = Inf;  % Avoid division by zero for i=j
    energyMatrix = 1 / (4 * pi * eps0) * erfc(alpha * r) ./ r;
    E_tot = E_tot + sum(sum(energyMatrix));
    
    % Fourier-space sum
    [kx, ky, kz] = ndgrid(-kMax:kMax, -kMax:kMax, -kMax:kMax);
    k_vec = 2 * pi * [kx(:), ky(:), kz(:)] ./ boxSize;
    k = sqrt(sum(k_vec.^2, 2));
    validIdx = k > 0;
    k_vec = k_vec(validIdx, :);
    k = k(validIdx);
    
    dotProd = pos * k_vec';
    sum_cos = cos(dotProd)' * dip;
    sum_sin = sin(dotProd)' * dip;
    energyFourier = 1 / (4 * pi * eps0) * 4 * pi * exp(-k.^2 / (4 * alpha^2)) ./ (k.^2) .* (sum_cos.^2 + sum_sin.^2);
    E_tot = E_tot + sum(energyFourier);
    
    % Self-interaction term
    selfInteraction = -1 / (4 * pi * eps0) * alpha / sqrt(pi) * sum(sum(dip.^2));
    E_tot = E_tot + selfInteraction;
end