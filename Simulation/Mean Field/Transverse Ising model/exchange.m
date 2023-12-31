function [d]=exchange(q, Jex, a, tau)

% This function performs a brute force summation of
% the q-dependent exchange coupling fo a non-Bravais lattice.
% q=[h k l] is the q-vector given in Miller-indicies.
% nr is the number of unit cells that should be summed
% in each direction.

% For N moments in the unit cell, there will be N 
% coupling parameters J_ij. Many of these will be symmetry related
% (e.g. J_ij=J_ji), so we just calculate J_1j. 
% The result is a (3x3xN) matrix, where the first two dimensions
% are the x,y and z components. The last dimension holds the coupling
% between different ions in the unit cell.

% Convert tau to a
tau = tau*a;
% Unit cell volume
vol = sum(a(1,:).*cross(a(2,:),a(3,:)));
% Reciprocal lattice unit vectors
b = [2*pi*cross(a(2,:),a(3,:))/vol
     2*pi*cross(a(3,:),a(1,:))/vol
     2*pi*cross(a(1,:),a(2,:))/vol];
% Convert q from Miller indicies to reciprocal angstroms
q = q*b;

% % Kronecker delta in x,y,and z
% delta = [1 0 0;0 1 0;0 0 1]; % replaced by built-in eq(m,n) function, -Yikai 2021.03.08

[x,y,z] = meshgrid(-1:1,-1:1,-1:1); % expand one unit cell from the center
hkl = [x(:) y(:) z(:)];

d = zeros(3,3,size(tau,1),size(tau,1));
for nt = 1:size(tau,1)
    for mt = 1:nt % avoid double counting
        r = hkl*a;
        r(:,1) = r(:,1) - tau(nt,1) + tau(mt,1);
        r(:,2) = r(:,2) - tau(nt,2) + tau(mt,2);
        r(:,3) = r(:,3) - tau(nt,3) + tau(mt,3);
        rr = sum(r.^2,2);
        r(rr < 0.01 | rr > max(vecnorm(a,1).^2), :) = []; % remove singularity and cut off at one unit cell
        rr(rr < 0.01 | rr > max(vecnorm(a,1)).^2) = [];
        exp_qr = exp(-1i*q*r');
        for n = 1:3
            for m = 1:3
%                 d(n,m,nt,mt) = exp_qr*ones(length(rr),1)*Jex*eq(n,m); % exchange interaction
                d(n,m,nt,mt) = exp_qr*(rr<=14)*Jex*eq(n,m); % original code, 13.906 ~ 14 (NN distance^2)
            end
        end
        %  d(:,:,nt,mt) = d(:,:,nt,mt) + (4*pi/3)*0.01389*eye(3)/4; %Lorentz
        d(:,:,mt,nt) = conj(d(:,:,nt,mt)); % symmetrize the tensor
    end
end
return
end
