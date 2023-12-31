function [spin_nuc, Enn] = therm_nuc_a(const, ion, params, tt, ff, Espin, spin_nuc, Enn)
for ii = 1:size(params.pos,1)
    % electronic spin components
    spx = Espin(ii,1);
    spy = Espin(ii,2);
    spz = Espin(ii,3);

    H_hyp = ion.A(ion.idx) * (spx*ion.Ix + spy*ion.Iy + spz*ion.Iz); % nuclear hyperfine interaction (meV)
    HzI = -ion.nLande(ion.idx) * const.muN * const.J2meV * (params.field(1,ff)*ion.Ix +...,
        params.field(2,ff)*ion.Iy + params.field(3,ff)*ion.Iz); % nuclear Zeeman interaction (meV)
    ham_nuc = H_hyp + HzI; % total nuclear spin hamiltonian
    [eigen_n, ~] = eig(ham_nuc);
    coef_n = randSphN(1,size(ham_nuc,1)); % random point of the nuclear spin configuration
    E_new = real(coef_n' * eigen_n' * ham_nuc * eigen_n * coef_n); % new energy

    % use Glauber algorithm to thermally excite the nuclear spin
    crit = rand; % update criteria
    beta = 1/ (const.J2meV * const.kB * params.temp(tt)); % [meV]
    prob = 1 / ( 1 + exp( (E_new - Enn(ii)) * beta) ); % probability critirion
    if prob >= crit
        Ispx = real( coef_n' * eigen_n' * ion.Ix * eigen_n * coef_n)';
        Ispy = real( coef_n' * eigen_n' * ion.Iy * eigen_n * coef_n)';
        Ispz = real( coef_n' * eigen_n' * ion.Iz * eigen_n * coef_n)';
        Enn(ii) = E_new;
        spin_nuc(ii,:) = [Ispx Ispy Ispz];
    end
end

end