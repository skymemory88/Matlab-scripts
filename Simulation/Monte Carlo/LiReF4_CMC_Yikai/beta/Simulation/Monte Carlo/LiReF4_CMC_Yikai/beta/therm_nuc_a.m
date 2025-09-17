function nSpin = therm_nuc_a(const, beta, ion, params, ff, eSpin, nSpin)
for ii = params.isoIdx
    % electronic spin components
    spx = eSpin(ii,1);
    spy = eSpin(ii,2);
    spz = eSpin(ii,3);

    nspx = nSpin(ii,1);
    nspy = nSpin(ii,2);
    nspz = nSpin(ii,3);

    % Current energy
    EzI = const.Ngfac * (params.field(1,ff)*nspx + params.field(2,ff)*nspy + params.field(3,ff)*nspz); % Nuclear Zeeman interaction
    E_hyp = ion.A(ion.idx) * (spx*nspx + spy*nspy + spz*nspz); % Hyperfine interaction
    E_old = EzI + E_hyp;
    % new energy
    H_hyp = ion.A(ion.idx) * (spx*ion.Ix + spy*ion.Iy + spz*ion.Iz); % nuclear hyperfine interaction (meV)
    HzI = const.Ngfac * (params.field(1,ff)*ion.Ix + params.field(2,ff)*ion.Iy + params.field(3,ff)*ion.Iz); % nuclear Zeeman interaction (meV)
    ham_nuc = H_hyp + HzI; % total nuclear spin hamiltonian
    
    [eigen_n, ~] = eig(ham_nuc);
    coef_n = randSphN(1,size(ham_nuc,1)); % Generate a random point on the bloch sphere
    E_new = real(coef_n' * eigen_n' * ham_nuc * eigen_n * coef_n); % new energy

    % use Glauber algorithm to thermally excite the nuclear spin
    crit = rand; % update criteria
    % prob = 1 / ( 1 + exp( (E_new - E_old) * beta) ); % probability critirion
    prob = exp(-(E_new - E_old) * beta ); % metropolis algorithm
    if prob >= crit
        Ispx = real( coef_n' * eigen_n' * ion.Ix * eigen_n * coef_n)';
        Ispy = real( coef_n' * eigen_n' * ion.Iy * eigen_n * coef_n)';
        Ispz = real( coef_n' * eigen_n' * ion.Iz * eigen_n * coef_n)';
        nSpin(ii,:) = [Ispx Ispy Ispz];
    end
end

end