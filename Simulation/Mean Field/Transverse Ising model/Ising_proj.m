
muB = 9.27401e-24; % Bohr magneton [J/T]
J2meV = 6.24151e+21; % [mev/J]

% % Crystal field parameter coefficients:
% ion.B = [-60.0   0.350   3.60   0.00   0.000400   0.0655   0.0098]; % Ho
ion.B = [-60.0   0.350   3.60   0.00   0.000400   0.0700   0.0098]; % Ho -- Phys. Rev. B 75, 054426 (2007)
ion.B = ion.B ./ 1e3; % ueV to meV
% ion.cfRot = 11;
ion.cfRot = 0;

ion.J = 8;
ion.gLande = 1.25;
Bx = linspace(0,9,101); % transverse field along x
By = zeros(size(Bx));
Bz = zeros(size(Bx));
fields = [Bx' By' Bz'];

% Pauli matrices
sigx = [0   1
        1   0];
sigy = [ 0   -1i
        1i  0];
sigz = [1   0
        0  -1];

Jz = diag(ion.J:-1:-ion.J); % Jz = -J, -J+1, ... ,J-1,J
Jp = diag(sqrt((ion.J-((ion.J-1):-1:-ion.J) ).*(ion.J+1+( (ion.J-1):-1:-ion.J) )),1); % electronic spin ladder operator
Jm = Jp'; % electronic spin ladder operator
Jx = (Jp+Jm)/2;
Jy = (Jp-Jm)/2i;

Hcf = CF(ion.J, ion.B, ion.cfRot);

Czz = double.empty(size(fields,1),0);
Cx = double.empty(size(fields,1),0);
Cxx = double.empty(size(fields,1),0);
Cxy = double.empty(size(fields,1),0);
Cxz = double.empty(size(fields,1),0);
Cy = double.empty(size(fields,1),0);
Cyy = double.empty(size(fields,1),0);
Cyx = double.empty(size(fields,1),0);
Cyz = double.empty(size(fields,1),0);
energies = double.empty(size(Hcf,1), size(fields,1), 0);
for ii = 1:size(fields,1)
    Hz = -ion.gLande * muB*J2meV * (fields(ii,1) * Jx + fields(ii,2) * Jy + fields(ii,3) * Jz);
    ham = Hcf + Hz;

    % Diagonalize
    [wav, En] = eig(ham);
    En = real(diag(En)); % Take only the real part of the eigen-energy to form a diaganol matrix
    [En, n] = sort(En); % sort the energy from lowest to the highest
    wav = wav(:,n); % sort the eigen-vectors in its basis accordingly
%     En = En - min(En); % Normalize the energy amplitude to the lowest eigen-energy
    energies(:,ii,1) = En;

    En = En(1:2); % truncate to Ising space
    wav = wav(:,1:2); % truncate to Ising subspace

    jz = wav'*Jz*wav;
    [uni, mz] = eig(jz); % further diagonalize Jz
    theta = angle(mz);
    theta = theta(1); % retain the phase angle of the ground state.
%     theta = 0; % for debugging
    Rm = exp(-1i*theta);
    mz = mz*Rm; % Rotate to keep ground state eigen value real
    Czz(ii) = mean(abs(real(diag(sigz\mz))));
    wav = wav*uni*Rm;

    jx = wav'*Jx*wav;
    on_coef = linsolve([diag(eye(2)) diag(sigz)], diag(jx)); % decompose the diagonal elements
    Cx(ii) = on_coef(1);
    Cxz(ii) = on_coef(2);
    off_coef = linsolve([1 -1i; 1 1i], [diag(jx,1); diag(jx,-1)]); % decompose the off diagonal elements
    Cxx(ii) = off_coef(1);
    Cxy(ii) = off_coef(2);

    jy = wav'*Jy*wav;
    on_coef = linsolve([diag(eye(2)) diag(sigz)], diag(jy)); % decompose the diagonal elements
    Cy(ii) = on_coef(1);
    Cyz(ii) = on_coef(2);
    off_coef = linsolve([1 -1i; 1 1i], [diag(jy,1); diag(jy,-1)]); % decompose the off diagonal elements
    Cyx(ii) = off_coef(1);
    Cyy(ii) = off_coef(2);
end