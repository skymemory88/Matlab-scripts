function [dip, nN] = MF_dipole(q, N, a, tau)
% This function performs a brute force summation of
% the q-dependent dipole coupling fo a non-Bravais lattice.
% q=[h k l] is the q-vector given in Miller-indicies.
% N is the number of unit cells that should be summed in each direction.

tau = tau*a; % Convert tau to a
% vol = sum(a(1,:).*cross(a(2,:),a(3,:))); % Unit cell volume

% determine the dimensionality
if any(a(1,:))
    b1 = 2*pi./a(1,:); % reciprocal vector a*b = 2*pi
    b1(isinf(b1)) = 0;
    N1 = N; % coordination number in x
else
    b1 = [0 0 0];
    N1 = 0;
end

if any(a(2,:))
    b2 = 2*pi./a(2,:); % reciprocal vector a*b = 2*pi
    b2(isinf(b2)) = 0;
    N2 = N; % coordination number in y
else
    b2 = [0 0 0];
    N2 = 0;
end

if any(a(3,:))
    b3 = 2*pi./a(3,:); % reciprocal vector a*b = 2*pi
    b3(isinf(b3)) = 0;
    N3 = N; % coordination number in z
else
    b3 = [0 0 0];
    N3 = 0;
end

% Reciprocal lattice unit vectors
b = [b1; b2; b3];
% Convert q from Miller indicies to reciprocal angstroms
q = q*b;
% % Length of q
% qq=sqrt(sum(q.*q));

% For N moments in the unit cell, there will be N 
% coupling parameters J_ij. Many of these will be symmetry related
% (e.g. J_ij = J_ji), so we just calculate J_1j. 
% The result is a (3x3xN) matrix, where the first two dimensions
% are the x,y and z components. The last dimension holds the coupling
% between different ions in the unit cell.

[x,y,z] = meshgrid(-N1:N1,-N2:N2,-N3:N3);
hkl = [x(:) y(:) z(:)];
% hkl = [z(:) x(:) y(:)]; % use z x y to get nicer order in the list - but of no importance

r0 = hkl*a;
dip = zeros(3,3,size(tau,1),size(tau,1));
for nt = 1:size(tau,1)
    for mt = 1:nt
        r = r0;
        r(:,1) = r(:,1) - tau(nt,1) + tau(mt,1);
        r(:,2) = r(:,2) - tau(nt,2) + tau(mt,2);
        r(:,3) = r(:,3) - tau(nt,3) + tau(mt,3);
        r2 = vecnorm(r,2,2).^2;
        cut_off = find(r2 > (N * min(nonzeros(vecnorm(a,2,2))))^2 | r2<0.01); % define interaction range + remove singularity
        r(cut_off,:) = []; % clear Spins outside of the sphere
        r2(cut_off) = []; % and the central spin itself
        nN = length(r2);
        r3 = r2 .* vecnorm(r,2,2);
        r5 = r2 .* r3;
        exp_qr = exp(-1i*q*r');
        for nn = 1:3
            for mm = 1:3
                dip(nn,mm,nt,mt) = exp_qr*(3*r(:,nn).*r(:,mm)./r5 - eq(nn,mm)./r3);
                % dip(n,m,nt,mt) = dip(n,m,nt,mt) + exp_qr*(rr<14)*Jex*delta(n,m); % exchange
            end
        end
        %  dip(:,:,nt,mt) = dip(:,:,nt,mt) + (4*pi/3)*0.01389*eye(3)/4; %  term Lorentz (demagnetization)
        dip(:,:,mt,nt) = conj(dip(:,:,nt,mt)); % J_ij = J_ji*
    end
end
dip = squeeze(dip);
return
end