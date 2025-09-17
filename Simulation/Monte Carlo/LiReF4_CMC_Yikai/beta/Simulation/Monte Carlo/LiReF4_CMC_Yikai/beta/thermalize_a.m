function [dE, accpRate, Esi, coef, eSpin] = thermalize_a(ion, const, params, beta, Esi, hamI, coef, basis, eSpin, nSpin)
% Modified to use neighbor lists for dipole calculations

% Random update of the whole lattice
site = randperm(size(params.pos,1));

% Generate random rotations on the Bloch sphere
coords_new = randSph(length(site),'Cartesian');
[alp_new, bet_new, ~] = cart2sph(coords_new(:,1),coords_new(:,2),coords_new(:,3));
alp_new = alp_new + pi; % [-pi pi] -> [0 2pi]
bet_new = bet_new + pi/2; % [-pi/2 pi/2] -> [0 pi]
coef_new = [cos(bet_new'/2) ; sin(bet_new'/2).*exp(1i*alp_new')];

counter = 1; % iteration tracker
change = zeros(length(site),1); % Mark acceptances
dEi = zeros(length(site),1); % local energy change

for ii = site
    % Compute new effective spin moment
    wav = squeeze(coef_new(:, ii));
    spx = real(wav' * basis' * ion.Jx * basis * wav);
    spy = real(wav' * basis' * ion.Jy * basis * wav);
    spz = real(wav' * basis' * ion.Jz * basis * wav);
    spin_new = [spx spy spz];

    % Compute change of single-ion energy
    Esi_new = real(wav' * hamI * wav); % new single-ion energy
    if params.hyp && ismember(ii, params.isoIdx) % check for hyperfine interaction option
        Esi_new = Esi_new + ion.A(ion.idx) * dot(nSpin(ii,:), spin_new); % include hyperfine energy (meV)
    end
    dE_si = Esi_new - Esi(ii); % single-ion energy change

    % Use neighbor lists for dipole interaction
    Eint_old = CntrDip_a(params, const.gfac, eSpin, ii, eSpin(ii,:));
    Eint_new = CntrDip_a(params, const.gfac, eSpin, ii, spin_new);

    if ion.ex(ion.idx) % check exchange interaction strength
        Eint_old = Eint_old + exchange(params, ion, eSpin, ii, eSpin(ii,:)); % include exchange interaction energy
        Eint_new = Eint_new + exchange(params, ion, eSpin, ii, spin_new); % include exchange interaction energy
    end
    dE_int = Eint_new - Eint_old; % interaction energy change

    % Total energy difference
    dEi(counter) = dE_si + dE_int;

    % Metropolis algorithm
    crit = rand;
    prob = exp(-dEi(counter) * beta);
    if prob >= crit
        Esi(ii) = Esi_new;
        coef(:,ii) = wav;
        eSpin(ii, :) = spin_new;
        change(counter) = 1;
    end
    counter = counter + 1;
end

% Calculate global energy change and acceptance rate
dE = sum(dEi .* change);
accpRate = sum(change)/length(change);
end