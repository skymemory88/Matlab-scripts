kB = 1.3806e-23; % [J/K]
J2meV = 6.24151e+21; % [mev/J]
temp = 0.1;

pos = rand(20,3);
spin = rand(20,3);
wavT = rand(17,2);
Jz = diag([-8:1:8]);
z = [1; 0];
% current dipole interaction energy
dipE0 = double.empty(size(pos,1),0);
for ii = 1:size(pos,1)
    dipE0(ii,1) = probFunc(pos,spin,ii,spin(ii,:)); % initial dipole-dipole interaction energy
end

spinT = spin;
    parfor jj = 1:size(pos,1)
        spin_i = [0 0 1];
        dEz = 0;
        dipE = probFunc(pos, spinT, jj, spin_i); % dipole interaction energy
        dEd = dipE - dipE0(jj,1); % change in dipole interaction energy
        dE = dEz + dEd; % total change of energy

        % Glauber algorithm to accept or reject rotation
        prob = 1 / ( 1 + exp(2*dE/kB/J2meV/temp) ); % probability critirion
        if prob >= rand || dE < 0; spin(jj,:) = spin_i; end
    end

function E_dip = probFunc(pos, spin, idx, spin0)
% Calculate the magnetic dipole-dipole interaction energy between a particular spin
% and the rest of the spins in the ensemble.
% Input:
%   spin - n x 3 matrix containing the x, y, z components of n spins
%   pos - n x 3 matrix containing the x, y, z coordinates of n spins
%   idx - index of the particular spin
% Output:
%   E_tot - magnetic dipole-dipole interaction energy

% Loop over all other spins
for j = 1:size(pos, 1)
    if j == idx
        spin(j,:) = spin0; % flip the spin at the center with supplied value
    end

    % Position vector between spins i and j
    r_vec = pos(idx, :) - pos(j, :);
    r = vecnorm(r_vec);

    % Calculate interaction energy between spins i and j (meV)
    E_dip = (dot(spin(idx, :), spin(j, :)) / r^3 - 3 * (dot(spin(idx, :), r_vec) * dot(spin(j, :), r_vec)) / r^5);
end
end