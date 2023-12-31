function [d, nN] = MF_exchange(q, Jex, a, tau)
% This function performs a brute force summation of
% the q-dependent isotropic exchange coupling fo a non-Bravais lattice.
% a: lattice constancs, tau: basis vectors
% q = [h k l] is the q-vector given in Miller-indicies.
% final result is a 3 x 3 diagonal matrix for the three Cartesan dimensions
% Jex can be arbitrarily long array to accommodate NN, NNN, etc interactions

% determine the dimensionality
if any(a(1,:))
    b1 = 2*pi./a(1,:); % reciprocal vector a*b = 2*pi
    b1(isinf(b1)) = 0;
    N1 = length(Jex); % coordination number in x
else
    b1 = [0 0 0];
    N1 = 0;
end

if any(a(2,:))
    b2 = 2*pi./a(2,:); % reciprocal vector a*b = 2*pi
    b2(isinf(b2)) = 0;
    N2 = length(Jex); % coordination number in y
else
    b2 = [0 0 0];
    N2 = 0;
end

if any(a(3,:))
    b3 = 2*pi./a(3,:); % reciprocal vector a*b = 2*pi
    b3(isinf(b3)) = 0;
    N3 = length(Jex); % coordination number in z
else
    b3 = [0 0 0];
    N3 = 0;
end

% Reciprocal lattice unit vectors
b = [b1; b2; b3];
% Convert q from Miller indicies to reciprocal angstroms
q = q*b;

[x,y,z] = meshgrid(-N1:N1, -N2:N2, -N3:N3); % location of neighbours
hkl = [x(:) y(:) z(:)]; % reciprocal lattice coordinates
r = hkl*a; % use real distance
tau = tau*a; % lattice basis

% place the neighbours in real space
d = zeros(3, 3);
if size(tau,1) > 1
    r = repmat(r',1,size(tau,1))'; % make copies for basis vectors
    for nt = 1:size(tau,1)
        r(1+nt*size(hkl,1):(nt+1)*size(hkl,1), 1) = r(1:size(hkl,1),1) + tau(nt,1);
        r(1+nt*size(hkl,1):(nt+1)*size(hkl,1), 2) = r(1:size(hkl,1),2) + tau(nt,2);
        r(1+nt*size(hkl,1):(nt+1)*size(hkl,1), 3) = r(1:size(hkl,1),3) + tau(nt,3);
    end
end
% rr = sum(r.^2, 2);
r2 = vecnorm(r, 2, 2).^2;
r(r2 == 0,:) = []; % avoid divergence
r2(r2 == 0) = [];
Jij = zeros(size(r2));
nDst2 = unique(r2); % remove duplicates of r^2 values
idx = cell(length(Jex),1);
for ii = 1:length(Jex) % iterate through nn, nnn, etc
    idx{ii} = find(r2 == nDst2(ii))';
    Jij(idx{ii}) = Jex(ii);
end
rridx = [idx{:}]';
r = r(rridx,:); % Remove points outside of cut-off range and singularities
r2 = r2(rridx,:);
Jij = nonzeros(Jij)';
nN = length(r2); % number of neighbours

exp_qr = exp(-1i*q*r')';
for n = 1:3 % m,n = x,y,z
    for m = 1:3
        d(n,m) = Jij * exp_qr * eq(n,m); % isotropic
    end
end
d(:,:) = conj(d(:,:)); % J_ij = J_ji*
return
end