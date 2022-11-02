Options.damp = 0.2; % update damping

muB = 9.27401e-24; % Bohr magneton [J/T]
muN = 3.15245e-5; % [meV/T]
mu0 = 4*pi*1e-7; % [H/m]
kB = 8.61733e-2; % Boltzmann constant [meV/K];
J2meV = 6.24151e+21; % [mev/J]

unitN = 4; % number of sites per unit cell
dpRng = 100; % dipole summation range

% Crystal field parameter coefficients:
% ion.B = [-60.0   0.350   3.60   0.00   0.000400   0.0655   0.0098]; % Ho
ion.B = [-60.0   0.350   3.60   0.00   0.000400   0.0700   0.0098]; % Ho -- Phys. Rev. B 75, 054426 (2007)
ion.B = ion.B ./ 1e3; % ueV to meV
% ion.cfRot = 11; % CF rotation
ion.cfRot = 0;
ion.J = 8; % Ho eletronic spin
ion.gLande = 1.25; % Lande factor
ion.I = 7/2; % Ho nuclear spin
ion.nLande = 1.192; % nuclear Lande factor
ion.A = -3.36e-3 .* [1 1 1]; % [meV] hyperfine interaction strength
% ion.A = [0 0 0]; % remove hyperfine interaction
ion.ex = 0.0001; % [meV] exchange interaction strength
% lattice constants
ion.abc = [5.175 0 0
           0 5.175 0
           0 0 10.75];
% basis vector
ion.tau = [ 0    0    0
            0   1/2  1/4
           1/2  1/2  1/2
           1/2   0   3/4];

Bx = linspace(0,9,201); % transverse field along x
By = zeros(size(Bx));
Bz = zeros(size(Bx));

qvec = [0 0 0];

temp = 0; % temperature
fields = [Bx' By' Bz'];

% Pauli matrices
sigx = [0   1
        1   0];
sigy = [ 0   -1i
        1i  0];
sigz = [1   0
        0  -1];

% Initiate J operators
Jz = diag(ion.J:-1:-ion.J); % Jz = -J, -J+1, ... ,J-1,J
Jp = diag(sqrt( (ion.J-((ion.J-1):-1:-ion.J) ).*(ion.J+1+( (ion.J-1):-1:-ion.J) )),1); % electronic spin ladder operator
Jm = Jp'; % electronic spin ladder operator
Jx = (Jp+Jm)/2;
Jy = (Jp-Jm)/2i;

if any(ion.A)
    % Initiate I operators
    Iz = diag(ion.I:-1:-ion.I);
    Ispan = size(Iz,1); % expansion factor of Hilbert space
    Ip = diag(sqrt( (ion.I-((ion.I-1):-1:-ion.I)).*(ion.I+1+((ion.I-1):-1:-ion.I)) ),1);
    Im = kron(eye(size(sigz)), Ip');
    Ip = kron(eye(size(sigz)), Ip);

    Ix = (Ip + Im)/2;
    Iy = (Ip - Im)/2i;
    Iz = kron(eye(size(sigz)), Iz);
else
    Ispan = 1;
end

Vc = sum( ion.abc(1,:) .* cross(ion.abc(2,:), ion.abc(3,:)) ); % Volume of unit cell [Ang^-3]
gfac = mu0/4/pi * (ion.gLande * muB)^2 * J2meV * 10^30; % mu0/(4pi).(gL.muB)^2 [meV.Ang^3]
eins = zeros(3,3,4,4);
eins(1,1,:,:) = 1; eins(2,2,:,:) = 1; eins(3,3,:,:) = 1;
D = zeros(3, 3, unitN, unitN, size(qvec,1)); % dipole interaction container
J = zeros(3, 3, unitN, unitN, size(qvec,1)); % dipole interaction container
parfor jj = 1:size(qvec,1)
% for jj = 1:size(qvec,1) % for debugging
    D(:,:,:,:,jj) = gfac * (MF_dipole(qvec(jj,:), dpRng, ion.abc, ion.tau) + 4*pi*eins/3/Vc); % dipole + Lorenz term
    J(:,:,:,:,jj) = exchange(qvec(jj,:), ion.ex, ion.abc, ion.tau); % [meV] AFM exchange interaction
end
Dip = sum(sum(D,4),3)/unitN; % average over the unit cell
Jex = sum(sum(J,4),3)/unitN;
Dip = Dip(3,3); % truncate to Ising interaction (Jzz)
Jex = Jex(3,3);
Jij = Dip-Jex; % effective interaction strength
% Jij = 6.3992e-3; % PRB 97, 214430 (2018)

Hcf = CF(ion.J, ion.B, ion.cfRot); % crystal field hamiltonian
Czz = double.empty(size(fields,1), 0);
Cx = double.empty(size(fields,1), 0);
Cxx = double.empty(size(fields,1), 0);
Cxy = double.empty(size(fields,1), 0);
Cxz = double.empty(size(fields,1), 0);
Cy = double.empty(size(fields,1), 0);
Cyy = double.empty(size(fields,1), 0);
Cyx = double.empty(size(fields,1), 0);
Cyz = double.empty(size(fields,1), 0);
Js = double.empty(size(fields,1), 3, 0);
Bx = double.empty(size(fields,1), 0);
if any(ion.A); Is = double.empty(size(fields,1), 3, 0); end
eigenE = double.empty(size(sigz,1)*Ispan, size(fields,1), 0);
eigenW = double.empty(size(sigz,1)*Ispan, size(sigz,1)*Ispan, size(fields,1), 0);
for ii = 1:size(fields,1)
    Hz = -ion.gLande * muB*J2meV * (fields(ii,1) * Jx + fields(ii,2) * Jy + fields(ii,3) * Jz);
    ham = Hcf + Hz;

    % Diagonalize
    [uni, ee] = eig(ham); % diagonalization by unitary transformation
    ee = real(diag(ee)); % Take only the real part of the eigen-energy to form a diaganol matrix
    [~, n] = sort(ee); % sort the energy from lowest to the highest
    uni = uni(:,n);
    uni = uni(:,1:2); % truncate to Ising subspace

    jz = uni'*Jz*uni; % truncate Jz operator
    [ket, jz] = eig(jz); % second rotation to diagonalize Jz
    if jz(1) < 0 % ensure Czz > 0
        ket = flip(ket,2);
        jz = ket'*uni'*Jz*uni*ket;
    end

    % fix the relative phase between the two basis vectors
    thta = angle(ket(1,1));
    ket(:,1) = ket(:,1)*exp(-1i*thta);
    thta = angle(ket(1,2));
    ket(:,2) = ket(:,2)*exp(-1i*thta);

    jx = ket'*uni'*Jx*uni*ket; % truncate Jx operator
    jy = ket'*uni'*Jy*uni*ket; % truncate Jy operator
    
    % compute coefficients for the 2x2 spin operators
    Czz(ii) = mean(linsolve(sigz, diag(jz)));
%     Jij = Czz(ii)^2 * Jij; % spin-spin interaction scaling
    on_coef = linsolve([1  1; 1 -1], diag(jx)); % decompose the diagonal elements
    Cx(ii) = on_coef(1);
    Cxz(ii) = on_coef(2);
    off_coef = linsolve([1 -1i; 1 1i], [diag(jx,1); diag(jx,-1)]); % decompose the off-diagonal elements
    Cxx(ii) = off_coef(1);
    Cxy(ii) = off_coef(2);

    on_coef = linsolve([1  1; 1 -1], diag(jy));
    Cy(ii) = on_coef(1);
    Cyz(ii) = on_coef(2);
    off_coef = linsolve([1 -1i; 1 1i], [diag(jy,1); diag(jy,-1)]);
    Cyx(ii) = off_coef(1);
    Cyy(ii) = off_coef(2);    
    
    Ht = ket'*uni'*ham*uni*ket; % (Hcf + Hx) hamiltonian
    Bx(ii,1) = mean([abs(Ht(2,1)) abs(Ht(1,2))]);
    S_mf = [0 0 jz(1)]; % initial guesses for electronic spin operators
    beta = 1/(kB*temp); % Boltzman constant
    En_old = 0; % update mark between iteractions
    if any(ion.A)
        jx = kron(jx, eye(2*ion.I+1));
        jy = kron(jy, eye(2*ion.I+1));
        jz = kron(jz, eye(2*ion.I+1));
        Ht = kron(Ht, eye(2*ion.I+1));
    end

    delta = 1e-6; % convergence criteria
    counter = 1; % while loop iterator
    while true
        h_int = S_mf(3) * Jij; % spin-spin interaction part 1
        if any(ion.A) % Case: finite hyperfine interaction
            Hint = -h_int * jz; % spin-spin interaction
            HzeeI = -ion.nLande * muN * (fields(ii,1)*Ix + fields(ii,2)*Iy + fields(ii,3)*Iz); % nuclear Zeeman interaction
%             HzeeI = zeros(size(Hint)); % Omit nuclear Zeeman interaction
            Hyp = ion.A(1)*jx*Ix + ion.A(2)*jy*Iy + ion.A(3)*jz*Iz; % hyperfine interaction
            ham = Ht + Hint + Hyp + HzeeI;
            [wav, En] = eig(ham);
            [En, n] = sort(real(diag(En))); % sort the energy in ascending order
            wav = wav(:,n); % sort the eigen-vector accordingly
            EnN = En - min(En); % normalize the eigenenergy to the ground state
            Z = exp(-EnN*beta)/sum(exp(-EnN*beta)); % Calculate the partition function weight
            if temp == 0
                % Calculate the expectation value of electronic spins
                newAvSx = real( wav(:,1)'*jx*wav(:,1) );
                newAvSy = real( wav(:,1)'*jy*wav(:,1) );
                newAvSz = real( wav(:,1)'*jz*wav(:,1) );
                % Calculate the expectation value of nuclear spins
                newAvIx = real(wav(:,1)'*Ix*wav(:,1))';
                newAvIy = real(wav(:,1)'*Iy*wav(:,1))';
                newAvIz = real(wav(:,1)'*Iz*wav(:,1))';
            else
                % Calculate the expectation value of electronic spins
                newAvSx = real( wav(:,1)'*jx*wav(:,1) );
                newAvSy = real( wav(:,1)'*jy*wav(:,1) );
                newAvSz = real( wav(:,1)'*jz*wav(:,1) );
                % Calculate the expectation value of nuclear spins
                newAvIx = real(wav(:,1)'*Ix*wav(:,1))';
                newAvIy = real(wav(:,1)'*Iy*wav(:,1))';
                newAvIz = real(wav(:,1)'*Iz*wav(:,1))';
            end
            newAvI = [newAvIx  newAvIy  newAvIz];
        else % Case: no hyperfine interaction
            Hint = -h_int * jz; % spin-spin interaction
            ham = Ht + Hint;
            [wav, En] = eig(ham);
            En = real(diag(En));
            EnN = En-min(En); % normalize the eigenenergy to the ground state
            Z = exp(-EnN*beta)/sum(exp(-EnN*beta)); % Calculate the partition function weight
            if temp == 0
                newAvSx = real( wav(:,1)'*jx*wav(:,1) );
                newAvSy = real( wav(:,1)'*jy*wav(:,1) );
                newAvSz = real( wav(:,1)'*jz*wav(:,1) );
            else
                newAvSx = real( diag(wav'*jx*wav) )'*Z;
                newAvSy = real( diag(wav'*jy*wav) )'*Z;
                newAvSz = real( diag(wav'*jz*wav) )'*Z;
            end
        end
        newAvS = [newAvSx  newAvSy  newAvSz];

        if ~any( abs(EnN - En_old) >= delta ) && ~any( abs(newAvS(3) - S_mf(3)) >= delta )
            eigenE(:,ii,1) = squeeze(En); % save energies before normalization
            eigenW(:,:,ii,1) = squeeze(wav); % save eigen-functions
            Js(ii,:,1) = S_mf;
            if exist('newAvI','var'); Is(ii,:,1) = newAvI; end
            break
        elseif counter > 1e5
            eigenE(:,ii,1) = squeeze(En); % save energies before normalization
            eigenW(:,:,ii,1) = squeeze(wav); % save eigen-functions
            Js(ii,:,1) = S_mf;
            if exist('newAvI','var'); Is(ii,:,1) = newAvI; end
            fprintf('Iteration exceeded limits, dE = %.2e \n', nonzeros(EnN - En_old));
            break
        else
            En_old = EnN; % store energy of the last round
            S_mf = Options.damp*S_mf + (1-Options.damp)*newAvS; % update the order parameter (internal field)
            counter = counter + 1;
            if mod(counter, 2e4) == 0
                fprintf('Iternation: %u, Magnetic field: %.2f, Temperature: %.2f\n', counter, fields(ii,1), temp)
            end
        end
    end
    if any(ismember(round(linspace(1,size(fields,1),20)), ii))
        fprintf('Calculating for |B/J| = %.2f, progress: %.1f%% \n', fields(ii,1), ii/size(fields,1)*100);
    end
end
