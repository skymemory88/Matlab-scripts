function dip = dipole_direct(q,N,a)
% This function performs a brute force summation of
% the q-dependent dipole coupling fo a non-Bravais lattice.
% q = [h k l] is the q-vector given in Miller-indicies.
% N is the number of unit cells that should be summed in each direction.

tau = [ 0    0    0
        0   1/2  1/4
       1/2  1/2  1/2
       1/2   0   3/4];

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
% % Length of q
% qq = sqrt(sum(q.*q));

% For NN moments in the unit cell, there will be NN 
% coupling parameters J_ij. Many of these will be symmetry related
% (e.g. J_ij=J_ji), so we just calculate J_1j. 
% The result is a (3x3xNNxNN) matrix, where the first two dimensions
% are the x,y and z components. The last dimension holds the coupling
% between different ions in the unit cell.

[x,y,z] = meshgrid(-N:N,-N:N,-N:N);
hkl = [x(:) y(:) z(:)];
r0 = hkl*a;

dip = zeros(3,3,size(tau,1),size(tau,1));
for nt = 1:size(tau,1)
    for mt = 1:nt % avoid double counting
        r = r0;
        r(:,1) = r(:,1) - tau(nt,1) + tau(mt,1);
        r(:,2) = r(:,2) - tau(nt,2) + tau(mt,2);
        r(:,3) = r(:,3) - tau(nt,3) + tau(mt,3);
        rr = sum(r.^2,2);
%         cut_off = find(rr>(N*5.162)^2 | rr<0.01); % original code
        cut_off = find(rr>(N * min(vecnorm(a,2,2)))^2 | rr<0.01); % points out of relevant range
        r(cut_off,:) = []; % clear Spins outside of the sphere
        rr(cut_off) = []; % and the central spin itself
        rr15 = rr.*sqrt(rr);
        rr25 = rr.*rr15;
        exp_qr = exp(-1i*q*r');
        for nn = 1:3
            for mm = 1:3
                dip(nn,mm,nt,mt) = exp_qr*(3*r(:,nn).*r(:,mm)./rr25-eq(nn,mm)./rr15);
            end
        end
        dip(:,:,mt,nt) = conj(dip(:,:,nt,mt)); % J_ij = J_ji*
        % d(:,:,nt,mt) = d(:,:,nt,mt)+(4*pi/3)*0.01389*eye(3)/4; % Lorentz
    end
end
return
end