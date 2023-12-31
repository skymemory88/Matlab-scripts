function dip = dipDir(q, N, a)
    % Initialize tau matrix
    tau = [0, 0, 0;
           0, 1/2, 1/4;
           1/2, 1/2, 1/2;
           1/2, 0, 3/4];
    
    % Convert tau to a
    tau = tau * a;
    
    % Calculate unit cell volume
    vol = dot(a(1, :), cross(a(2, :), a(3, :)));
    
    % Calculate reciprocal lattice unit vectors
    b = 2 * pi * cross(a(2:3, :), a([3, 1], :), 2) / vol;
    
    % Convert q from Miller indices to reciprocal angstroms
    q = q * b;
    
    % Generate meshgrid for hkl
    [x, y, z] = meshgrid(-N:N, -N:N, -N:N);
    hkl = [x(:), y(:), z(:)];
    r0 = hkl * a;
    
    % Initialize dipole tensor
    dip = zeros(3, 3, size(tau, 1), size(tau, 1));
    
    % Loop over tau indices
    for nt = 1:size(tau, 1)
        for mt = 1:nt % Avoid double counting
            r = bsxfun(@minus, r0, tau(nt, :) - tau(mt, :));
            rr = sum(r .^ 2, 2);
            
            % Apply cutoff
            cut_off = find(rr > (N * min(vecnorm(a, 2, 2))) ^ 2 | rr < 0.01);
            r(cut_off, :) = [];
            rr(cut_off) = [];
            
            rr15 = rr .* sqrt(rr);
            rr25 = rr .* rr15;
            
            exp_qr = exp(-1i * q * r');
            
            % Calculate dipole tensor components
            for nn = 1:3
                for mm = 1:3
                    dip(nn, mm, nt, mt) = exp_qr * (3 * r(:, nn) .* r(:, mm) ./ rr25 - (nn == mm) ./ rr15);
                end
            end
            
            % Apply symmetry relation
            dip(:, :, mt, nt) = conj(dip(:, :, nt, mt));
        end
    end
end
