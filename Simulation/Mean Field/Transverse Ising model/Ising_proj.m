
muB = 9.27401e-24; % Bohr magneton [J/T]
J2meV = 6.24151e+21; % [mev/J]

% % Crystal field parameter coefficients:
ion.B = [-60.0   0.350   3.60   0.00   0.000400   0.0655   0.0098]; % Ho
% ion.B = [-60.0   0.350   3.60   0.00   0.000400   0.0700   0.0098]; % Ho -- Phys. Rev. B 75, 054426 (2007)
ion.B = ion.B ./ 1e3; % ueV to meV
% ion.cfRot = 11;
ion.cfRot = 0;

ion.J = 8;
ion.gLande = 1.25;
Bx = linspace(0,9,201); % transverse field along x
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

rePhs = [ 1  -1;
         -1   1]; % flip the relative phase between |0> and |1>
Czz = double.empty(size(fields,1),0);
Cx = double.empty(size(fields,1),0);
Cxx = double.empty(size(fields,1),0);
Cxy = double.empty(size(fields,1),0);
Cxz = double.empty(size(fields,1),0);
Cy = double.empty(size(fields,1),0);
Cyy = double.empty(size(fields,1),0);
Cyx = double.empty(size(fields,1),0);
Cyz = double.empty(size(fields,1),0);
debug = double.empty(size(fields,1),0);
energies = double.empty(size(Hcf,1), size(fields,1), 0);
for ii = 1:size(fields,1)
    Hz = -ion.gLande * muB*J2meV * (fields(ii,1) * Jx + fields(ii,2) * Jy + fields(ii,3) * Jz);
    ham = Hcf + Hz;

    % Diagonalize
    [uni, En] = eig(ham, 'nobalance'); % diagonalization by unitary transformation
%     En = real(diag(En)); % Take only the real part of the eigen-energy to form a diaganol matrix
    En = real(diag(En)); % Take only the real part of the eigen-energy to form a diaganol matrix
    [En, n] = sort(En); % sort the energy from lowest to the highest
%     En = En - min(En); % Normalize the energy amplitude to the lowest eigen-energy
    uni = uni(:,n);
    energies(:,ii,1) = En;
    En = En(1:2); % truncate to Ising space
    uni = uni(:,1:2); % truncate to Ising subspace

    jz = uni'*Jz*uni;    
    [ket, jz] = eig(jz); % second rotation to diagonalize Jz
    if jz(1) < 0; ket = flip(ket,2); jz = ket'*uni'*Jz*uni*ket; end % ensure Czz > 0

    % fix the relative phase between the two basis vectors
    thta = angle(ket(1,1));
    ket(:,1) = ket(:,1)*exp(-1i*thta);
    thta = angle(ket(1,2));
    ket(:,2) = ket(:,2)*exp(-1i*thta);

    jx = ket'*uni'*Jx*uni*ket;
    jy = ket'*uni'*Jy*uni*ket;
    
    % compute coefficients for the 2x2 spin operators
    Czz(ii) = mean(linsolve(sigz, diag(jz)));
    on_coef = linsolve([1  1; 1 -1], diag(jx)); % decompose the diagonal elements
    Cx(ii) = on_coef(1);
    Cxz(ii) = on_coef(2);
    debug(ii) = Cx(ii) - Cxz(ii);
    off_coef = linsolve([1 -1i; 1 1i], [diag(jx,1); diag(jx,-1)]); % decompose the off-diagonal elements
    Cxx(ii) = off_coef(1);
    Cxy(ii) = off_coef(2);

    on_coef = linsolve([1  1; 1 -1], diag(jy));
    Cy(ii) = on_coef(1);
    Cyz(ii) = on_coef(2);
    off_coef = linsolve([1 -1i; 1 1i], [diag(jy,1); diag(jy,-1)]);
    Cyx(ii) = off_coef(1);
    Cyy(ii) = off_coef(2);
end